"""Dependency-free center-star multiple sequence alignment for FASTA files."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence


@dataclass(frozen=True)
class Scoring:
    """Integer scores used by the global pairwise aligner."""

    match: int = 2
    mismatch: int = -1
    gap: int = -2


def read_fasta(path: str | Path) -> list[tuple[str, str]]:
    """Read unaligned FASTA records and return (identifier, sequence) pairs."""
    records: list[tuple[str, str]] = []
    identifier: str | None = None
    sequence_parts: list[str] = []

    def save_record() -> None:
        if identifier is None:
            return
        sequence = "".join(sequence_parts).upper()
        if not sequence:
            raise ValueError(f"FASTA record {identifier!r} has no sequence")
        if "-" in sequence or "." in sequence:
            raise ValueError(
                f"FASTA record {identifier!r} contains gaps; input must be unaligned"
            )
        records.append((identifier, sequence))

    try:
        with Path(path).open(encoding="utf-8") as fasta_file:
            for line_number, raw_line in enumerate(fasta_file, start=1):
                line = raw_line.strip()
                if not line or line.startswith(";"):
                    continue
                if line.startswith(">"):
                    save_record()
                    identifier = line[1:].strip()
                    if not identifier:
                        raise ValueError(
                            f"FASTA header on line {line_number} has no identifier"
                        )
                    sequence_parts = []
                else:
                    if identifier is None:
                        raise ValueError(
                            f"Sequence data before the first FASTA header "
                            f"on line {line_number}"
                        )
                    sequence_parts.append("".join(line.split()))
    except OSError as error:
        raise OSError(f"Could not read FASTA file {path!s}: {error}") from error

    save_record()
    if len(records) < 2:
        raise ValueError("At least two FASTA sequences are required")
    identifiers = [name for name, _ in records]
    if len(set(identifiers)) != len(identifiers):
        raise ValueError("FASTA identifiers must be unique")
    return records


def needleman_wunsch(
    first: str, second: str, scoring: Scoring = Scoring()
) -> tuple[str, str, int]:
    """Return a highest-scoring global alignment and its score."""
    rows, columns = len(first) + 1, len(second) + 1
    scores = [[0] * columns for _ in range(rows)]
    trace = [[""] * columns for _ in range(rows)]

    for i in range(1, rows):
        scores[i][0] = i * scoring.gap
        trace[i][0] = "up"
    for j in range(1, columns):
        scores[0][j] = j * scoring.gap
        trace[0][j] = "left"

    for i in range(1, rows):
        for j in range(1, columns):
            diagonal = scores[i - 1][j - 1] + (
                scoring.match if first[i - 1] == second[j - 1] else scoring.mismatch
            )
            up = scores[i - 1][j] + scoring.gap
            left = scores[i][j - 1] + scoring.gap
            best = max(diagonal, up, left)
            scores[i][j] = best
            # Fixed tie-breaking makes results reproducible.
            trace[i][j] = (
                "diagonal" if diagonal == best else "up" if up == best else "left"
            )

    aligned_first: list[str] = []
    aligned_second: list[str] = []
    i, j = len(first), len(second)
    while i or j:
        direction = trace[i][j]
        if direction == "diagonal":
            aligned_first.append(first[i - 1])
            aligned_second.append(second[j - 1])
            i -= 1
            j -= 1
        elif direction == "up":
            aligned_first.append(first[i - 1])
            aligned_second.append("-")
            i -= 1
        else:
            aligned_first.append("-")
            aligned_second.append(second[j - 1])
            j -= 1

    aligned_first.reverse()
    aligned_second.reverse()
    return "".join(aligned_first), "".join(aligned_second), scores[-1][-1]


def center_star_alignment(
    records: Sequence[tuple[str, str]], scoring: Scoring = Scoring()
) -> list[tuple[str, str]]:
    """Align sequences using the center-star progressive heuristic.

    The center is the sequence with the greatest sum of pairwise global-alignment
    scores. Insertions from other sequences are merged into shared slots around
    the center; within a slot, insertions are right-aligned.
    """
    if len(records) < 2:
        raise ValueError("At least two sequences are required")
    names = [name for name, _ in records]
    sequences = [sequence.upper() for _, sequence in records]
    if len(set(names)) != len(names):
        raise ValueError("Sequence identifiers must be unique")
    if any(not sequence for sequence in sequences):
        raise ValueError("Sequences must not be empty")
    if any("-" in sequence or "." in sequence for sequence in sequences):
        raise ValueError("Input sequences must be unaligned and contain no gap symbols")

    totals = [0] * len(sequences)
    for i in range(len(sequences)):
        for j in range(i + 1, len(sequences)):
            _, _, score = needleman_wunsch(sequences[i], sequences[j], scoring)
            totals[i] += score
            totals[j] += score
    center_index = max(range(len(sequences)), key=totals.__getitem__)
    center = sequences[center_index]
    aligned_by_index: dict[int, tuple[list[str], list[list[str]]]] = {}

    for index, sequence in enumerate(sequences):
        if index == center_index:
            continue
        aligned_center, aligned_sequence, _ = needleman_wunsch(
            center, sequence, scoring
        )
        insertions = [[] for _ in range(len(center) + 1)]
        residues = [""] * len(center)
        center_position = 0
        for center_char, sequence_char in zip(aligned_center, aligned_sequence):
            if center_char == "-":
                insertions[center_position].append(sequence_char)
            else:
                residues[center_position] = sequence_char
                center_position += 1
        aligned_by_index[index] = residues, insertions

    slot_widths = [
        max(
            (
                len(insertions[slot])
                for _, insertions in aligned_by_index.values()
            ),
            default=0,
        )
        for slot in range(len(center) + 1)
    ]

    result: list[tuple[str, str]] = []
    for index, name in enumerate(names):
        pieces: list[str] = []
        if index == center_index:
            for slot, residue in enumerate(center):
                pieces.extend(["-"] * slot_widths[slot])
                pieces.append(residue)
            pieces.extend(["-"] * slot_widths[-1])
        else:
            residues, insertions = aligned_by_index[index]
            for slot, residue in enumerate(residues):
                insertion = "".join(insertions[slot])
                pieces.extend(["-"] * (slot_widths[slot] - len(insertion)))
                pieces.append(insertion)
                pieces.append(residue)
            pieces.extend(["-"] * (slot_widths[-1] - len(insertions[-1])))
            pieces.append("".join(insertions[-1]))
        result.append((name, "".join(pieces)))

    lengths = {len(sequence) for _, sequence in result}
    if len(lengths) != 1:
        raise RuntimeError("Internal error: generated alignment rows have different lengths")
    return result


def write_fasta(records: Sequence[tuple[str, str]], path: str | Path) -> None:
    """Write aligned records as FASTA, wrapping sequence lines at 80 characters."""
    try:
        with Path(path).open("w", encoding="utf-8", newline="\n") as output_file:
            for name, sequence in records:
                output_file.write(f">{name}\n")
                for start in range(0, len(sequence), 80):
                    output_file.write(sequence[start : start + 80] + "\n")
    except OSError as error:
        raise OSError(f"Could not write FASTA file {path!s}: {error}") from error


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Create a multiple sequence alignment with the center-star method."
    )
    parser.add_argument("input", help="Input FASTA file containing unaligned sequences")
    parser.add_argument(
        "-o", "--output", default="aligned.fasta", help="Output FASTA path"
    )
    parser.add_argument("--match", type=int, default=2, help="Match score (default: 2)")
    parser.add_argument(
        "--mismatch", type=int, default=-1, help="Mismatch score (default: -1)"
    )
    parser.add_argument("--gap", type=int, default=-2, help="Gap score (default: -2)")
    args = parser.parse_args()

    try:
        records = read_fasta(args.input)
        aligned = center_star_alignment(
            records, Scoring(match=args.match, mismatch=args.mismatch, gap=args.gap)
        )
        write_fasta(aligned, args.output)
    except (OSError, ValueError, RuntimeError) as error:
        parser.error(str(error))

    print(f"Aligned {len(aligned)} sequences; saved to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
