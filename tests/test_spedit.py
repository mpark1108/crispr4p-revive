import unittest

from crispr4p.spedit import (
    has_bsai,
    make_oligos,
)
from crispr4p.web_views import spedit_rows


# Published guide and oligo sequences from the SpEDIT paper's Table 2.
PUBLISHED_SPEDIT_OLIGOS = {
    "ade6": (
        "TTGATAGCAACAGTGGCGAC",
        "CTAGAGGTCTCGGACTTTGATAGCAACAGTGGCGACGTTTCGAGACCCTTCC",
        "GGAAGGGTCTCGAAACGTCGCCACTGTTGCTATCAAAGTCCGAGACCTCTAG",
    ),
    "ura4": (
        "CCTTGTATAATACCCTCGCC",
        "CTAGAGGTCTCGGACTCCTTGTATAATACCCTCGCCGTTTCGAGACCCTTCC",
        "GGAAGGGTCTCGAAACGGCGAGGGTATTATACAAGGAGTCCGAGACCTCTAG",
    ),
    "meu27": (
        "TATTAGCCTTTGAAGGATTT",
        "CTAGAGGTCTCGGACTTATTAGCCTTTGAAGGATTTGTTTCGAGACCCTTCC",
        "GGAAGGGTCTCGAAACAAATCCTTCAAAGGCTAATAAGTCCGAGACCTCTAG",
    ),
    "clr5": (
        "AGCTTGTGGCTGACCGTTAA",
        "CTAGAGGTCTCGGACTAGCTTGTGGCTGACCGTTAAGTTTCGAGACCCTTCC",
        "GGAAGGGTCTCGAAACTTAACGGTCAGCCACAAGCTAGTCCGAGACCTCTAG",
    ),
    "cup1:4xtetO": (
        "ATTTCTTTTGCTTTACGGTC",
        "CTAGAGGTCTCGGACTATTTCTTTTGCTTTACGGTCGTTTCGAGACCCTTCC",
        "GGAAGGGTCTCGAAACGACCGTAAAGCAAAAGAAATAGTCCGAGACCTCTAG",
    ),
    "cup1-GFP": (
        "GCTCAGGCTAAACGTCGGAA",
        "CTAGAGGTCTCGGACTGCTCAGGCTAAACGTCGGAAGTTTCGAGACCCTTCC",
        "GGAAGGGTCTCGAAACTTCCGACGTTTAGCCTGAGCAGTCCGAGACCTCTAG",
    ),
    "epe1": (
        "GGACTTTTAAGATGGATTCC",
        "CTAGAGGTCTCGGACTGGACTTTTAAGATGGATTCCGTTTCGAGACCCTTCC",
        "GGAAGGGTCTCGAAACGGAATCCATCTTAAAAGTCCAGTCCGAGACCTCTAG",
    ),
}


class TestSpeditOligos(unittest.TestCase):
    def test_published_spedit_examples(self) -> None:
        """Generated oligos must match all examples published by SpEDIT."""
        for name, (
            guide,
            expected_forward,
            expected_reverse,
        ) in PUBLISHED_SPEDIT_OLIGOS.items():
            with self.subTest(name=name):
                forward, reverse = make_oligos(guide)

                self.assertEqual(expected_forward, forward)
                self.assertEqual(expected_reverse, reverse)

    def test_normalizes_lowercase_and_whitespace(self) -> None:
        expected = make_oligos(
            "TTGATAGCAACAGTGGCGAC"
        )

        actual = make_oligos(
            "  ttgatagcaacagtggcgac\n"
        )

        self.assertEqual(expected, actual)

    def test_rejects_invalid_guides(self) -> None:
        cases = (
            ("ACGT", ValueError, "20-nt"),
            ("A" * 21, ValueError, "20-nt"),
            ("A" * 19 + "N", ValueError, "invalid nucleotide"),
            (12345, TypeError, "must be a string"),
        )
        for guide, error, message in cases:
            with self.subTest(guide=guide):
                with self.assertRaisesRegex(error, message):
                    make_oligos(guide)  # type: ignore[arg-type]

    def test_bsai_detection(self) -> None:
        for sequence in (
            "AAAAAGGTCTCAAAAAAAAA",
            "AAAAAGAGACCAAAAAAAAA",
            "  aaaaaggtctcaaaaaaaaa\n",
        ):
            with self.subTest(sequence=sequence):
                self.assertTrue(has_bsai(sequence))

        self.assertFalse(has_bsai("TTGATAGCAACAGTGGCGAC"))


def web_row(guide):
    primer = (
        guide,
        "LEGACY_FORWARD",
        "LEGACY_REVERSE",
        (100, 103),
        1,
        "TGG",
    )
    return [guide, primer]


class TestSpeditWebRows(unittest.TestCase):
    def test_preserves_candidate_order(self) -> None:
        first = "ACATTGGCTTACGACGGTCG"
        second = "TTGATAGCAACAGTGGCGAC"

        rows = spedit_rows([web_row(first), web_row(second)])

        self.assertEqual([first, second], [row["guide"] for row in rows])

    def test_generates_oligos_for_selected_candidate(self) -> None:
        guide = "TTGATAGCAACAGTGGCGAC"
        row = spedit_rows([web_row(guide)])[0]
        expected_forward, expected_reverse = make_oligos(guide)

        self.assertEqual(expected_forward, row["forward"])
        self.assertEqual(expected_reverse, row["reverse"])

    def test_internal_bsai_warning_stays_aligned(self) -> None:
        safe = "TTGATAGCAACAGTGGCGAC"
        unsafe = "TTTTGAATGGTCTCAGTTGT"
        rows = spedit_rows([web_row(safe), web_row(unsafe)])

        self.assertFalse(rows[0]["has_internal_bsai"])
        self.assertTrue(rows[1]["has_internal_bsai"])


if __name__ == "__main__":
    unittest.main()
