import unittest
from unittest.mock import patch

from crispr4p.crispr4p import NGG, PrimerDesign, chromosomeFasta
from crispr4p.guides import (
    find_guides,
    guide_primers,
)


PLUS_GUIDE = "ATACATACATACATACATAC"
MINUS_GUIDE = "TATGTATGTATGTATGTATG"


def synthetic_sequence():
    reverse_site = PrimerDesign.reverseComplement(MINUS_GUIDE + "AGG")
    return (
        "A" * 30
        + PLUS_GUIDE
        + "TGG"
        + "A" * 30
        + reverse_site
        + "A" * 30
    )


class TestGuideDiscovery(unittest.TestCase):
    def test_public_interval_keeps_guides_at_both_endpoints(self):
        spacer = "A" * 10
        interval = (
            PLUS_GUIDE
            + "TGG"
            + spacer
            + PrimerDesign.reverseComplement(MINUS_GUIDE + "AGG")
        )
        sequence = "A" * 20 + interval + "A" * 20
        chromosome = chromosomeFasta("synthetic description\n" + sequence)
        start = 21
        end = start + len(interval) - 1
        designer = PrimerDesign.__new__(PrimerDesign)
        designer.userNGGs = []

        designer._getUserNGGs(chromosome, start, end)
        primers = {
            guide.seed: designer.getPrimerGRNA(
                chromosome,
                start,
                end,
                guide,
            )
            for guide in designer.userNGGs
        }

        self.assertEqual({PLUS_GUIDE, MINUS_GUIDE}, set(primers))
        self.assertEqual((start + 20, start + 22), primers[PLUS_GUIDE][3])
        self.assertEqual((start + 33, start + 35), primers[MINUS_GUIDE][3])

    def test_discovers_both_strands_in_legacy_order_and_coordinates(self):
        sequence = synthetic_sequence()
        guides = find_guides(
            sequence,
            "synthetic",
            0,
            len(sequence) - 1,
            hit_factory=NGG,
            reverse_complement=PrimerDesign.reverseComplement,
        )

        self.assertEqual(
            [
                (51, 1, PLUS_GUIDE, "TGG"),
                (51, -1, MINUS_GUIDE, "AGG"),
            ],
            [
                (guide.pos, guide.strand, guide.seed, guide.pam)
                for guide in guides
            ],
        )
        self.assertEqual(
            [
                (
                    PLUS_GUIDE,
                    "ACATACATACgttttagagctagaaatagcaagttaaaataa",
                    "ATGTATGTATttcttcggtacaggttatgttttttggcaaca",
                    (51, 53),
                    1,
                    "TGG",
                ),
                (
                    MINUS_GUIDE,
                    "TGTATGTATGgttttagagctagaaatagcaagttaaaataa",
                    "TACATACATAttcttcggtacaggttatgttttttggcaaca",
                    (84, 86),
                    -1,
                    "AGG",
                ),
            ],
            [
                guide_primers(
                    sequence,
                    0,
                    len(sequence) - 1,
                    guide,
                    PrimerDesign.reverseComplement,
                )
                for guide in guides
            ],
        )

    def test_filters_a_seed_repeated_at_distinct_pam_sites(self):
        sequence = (
            "A" * 30
            + PLUS_GUIDE
            + "TGG"
            + "A" * 10
            + PLUS_GUIDE
            + "AGG"
            + "A" * 30
        )

        guides = find_guides(
            sequence,
            "synthetic",
            0,
            len(sequence) - 1,
            hit_factory=NGG,
            reverse_complement=PrimerDesign.reverseComplement,
        )

        self.assertEqual([], guides)

    def test_ignores_a_partial_site_at_the_interval_boundary(self):
        sequence = "GG" + "A" * 30 + PLUS_GUIDE + "TGG"

        guides = find_guides(
            sequence,
            "synthetic",
            0,
            len(sequence) - 1,
            hit_factory=NGG,
            reverse_complement=PrimerDesign.reverseComplement,
        )

        self.assertEqual(
            [(PLUS_GUIDE, "TGG")],
            [(guide.seed, guide.pam) for guide in guides],
        )

    def test_rejects_an_interval_without_ngg_sites(self):
        with self.assertRaisesRegex(AssertionError, "No nGG found in your input"):
            find_guides(
                "A" * 100,
                "synthetic",
                0,
                99,
                hit_factory=NGG,
                reverse_complement=PrimerDesign.reverseComplement,
            )

        chromosome = chromosomeFasta("synthetic description\n" + "A" * 100)
        designer = PrimerDesign.__new__(PrimerDesign)
        designer.userNGGs = [object()]
        with self.assertRaisesRegex(AssertionError, "No nGG found in your input"):
            designer._getUserNGGs(chromosome, 1, 100)
        self.assertEqual([], designer.userNGGs)

    def test_primer_design_methods_remain_legacy_adapters(self):
        sequence = synthetic_sequence()
        chromosome = chromosomeFasta("synthetic description\n" + sequence)
        designer = PrimerDesign.__new__(PrimerDesign)
        designer.userNGGs = []
        guide = NGG("synthetic", 51, 1, PLUS_GUIDE, "TGG")
        discovered = [guide]

        with patch(
            "crispr4p.crispr4p.find_guides",
            return_value=discovered,
        ) as discover:
            result = designer._getUserNGGs(
                chromosome,
                1,
                len(sequence),
            )

        self.assertIsNone(result)
        self.assertIs(discovered, designer.userNGGs)
        discover.assert_called_once_with(
            chromosome.sequence,
            chromosome.name,
            0,
            len(sequence) - 1,
            hit_factory=NGG,
            reverse_complement=designer.reverseComplement,
            seed_length=20,
        )

        sentinel = object()
        with patch(
            "crispr4p.crispr4p.guide_primers",
            return_value=sentinel,
        ) as build_primer:
            result = designer.getPrimerGRNA(
                chromosome,
                1,
                len(sequence),
                guide,
            )

        self.assertIs(sentinel, result)
        build_primer.assert_called_once_with(
            chromosome.sequence,
            0,
            len(sequence) - 1,
            guide,
            designer.reverseComplement,
        )


if __name__ == "__main__":
    unittest.main()
