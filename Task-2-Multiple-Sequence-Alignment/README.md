# Multiple Sequence Alignment in Python

A small, dependency-free multiple sequence alignment (MSA) program for unaligned
DNA or protein sequences in FASTA format. It uses the **center-star heuristic**:

1. Globally align every sequence pair with Needleman–Wunsch.
2. Select the sequence with the largest total pairwise alignment score as the
   center sequence.
3. Align every other sequence to the center and merge insertions into shared
   alignment columns.

This is an educational heuristic, not a replacement for established MSA
software such as MAFFT or MUSCLE. Pairwise alignments use a simple linear gap
penalty and match/mismatch scores; the score values can be configured on the
command line.

## Requirements

- Python 3.10 or newer
- No third-party packages

## Run

From this folder:

```text
python msa.py sample.fasta --output aligned.fasta
```

Set scoring values if needed:

```text
python msa.py sample.fasta -o aligned.fasta --match 2 --mismatch -1 --gap -2
```

The input must contain at least two uniquely named, unaligned FASTA records.
Sequence lines can span multiple lines. The output is aligned FASTA with `-`
characters marking gaps.

## Run tests

```text
python -m unittest discover -s tests -v
```

## Project files

- `msa.py` — FASTA parser/writer, global pairwise aligner, center-star MSA, CLI.
- `tests/test_msa.py` — tests for pairwise alignment, MSA invariants, and FASTA.
- `sample.fasta` — small example input.
