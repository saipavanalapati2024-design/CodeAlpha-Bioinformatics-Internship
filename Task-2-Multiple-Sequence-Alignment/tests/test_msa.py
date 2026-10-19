import tempfile
import unittest
from pathlib import Path

from msa import Scoring, center_star_alignment, needleman_wunsch, read_fasta, write_fasta


class NeedlemanWunschTests(unittest.TestCase):
    def test_global_alignment_preserves_both_sequences(self):
        first, second, _ = needleman_wunsch("ACGT", "AGT")
        self.assertEqual(first.replace("-", ""), "ACGT")
        self.assertEqual(second.replace("-", ""), "AGT")
        self.assertEqual(len(first), len(second))
        self.assertIn("-", second)

    def test_scoring_can_change_alignment_score(self):
        _, _, score = needleman_wunsch("AA", "A", Scoring(match=3, mismatch=-2, gap=-4))
        self.assertEqual(score, -1)


class CenterStarTests(unittest.TestCase):
    def test_rows_are_equal_length_and_ungap_to_original(self):
        source = [("a", "ACGT"), ("b", "AGT"), ("c", "ACGGT")]
        aligned = center_star_alignment(source)
        self.assertEqual([name for name, _ in aligned], ["a", "b", "c"])
        self.assertEqual(len({len(sequence) for _, sequence in aligned}), 1)
        for (name, original), (_, result) in zip(source, aligned):
            self.assertEqual(result.replace("-", ""), original)

    def test_rejects_duplicate_names(self):
        with self.assertRaisesRegex(ValueError, "unique"):
            center_star_alignment([("same", "AC"), ("same", "AG")])

    def test_rejects_already_gapped_sequences(self):
        with self.assertRaisesRegex(ValueError, "unaligned"):
            center_star_alignment([("a", "A-C"), ("b", "AGC")])


class FastaTests(unittest.TestCase):
    def test_reads_multiline_fasta_and_writes_aligned_fasta(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            source = Path(temp_dir) / "input.fasta"
            output = Path(temp_dir) / "output.fasta"
            source.write_text(">seq1 description\nac\ngt\n>seq2\nAGT\n", encoding="utf-8")
            records = read_fasta(source)
            self.assertEqual(records, [("seq1 description", "ACGT"), ("seq2", "AGT")])
            write_fasta(center_star_alignment(records), output)
            self.assertTrue(output.read_text(encoding="utf-8").startswith(">seq1 description\n"))

    def test_rejects_sequence_before_header(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            source = Path(temp_dir) / "input.fasta"
            source.write_text("ACGT\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "before the first FASTA header"):
                read_fasta(source)


if __name__ == "__main__":
    unittest.main()
