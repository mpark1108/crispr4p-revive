import unittest
from unittest.mock import patch

from crispr4p.coordinates import (
    cut_from_pam,
    hit_coordinates,
    slice_bounds,
)
from crispr4p.crispr4p import NGG, PrimerDesign, chromosomeFasta


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


class TestCoordinateConversions(unittest.TestCase):
    def test_converts_biological_intervals_to_python_slices(self):
        self.assertEqual((0, 1), slice_bounds(1, 1))
        self.assertEqual((9, 20), slice_bounds(10, 20))

        with self.assertRaisesRegex(ValueError, "invalid 1-based interval"):
            slice_bounds(0, 20)

    def test_converts_plus_and_minus_pams_to_cut_boundaries(self):
        self.assertEqual(
            (47, 48),
            cut_from_pam(51, 53, 1),
        )
        self.assertEqual(
            (89, 90),
            cut_from_pam(84, 86, -1),
        )

    def test_normalizes_and_verifies_indexed_hits_on_both_strands(self):
        sequence = synthetic_sequence()

        self.assertEqual(
            ((51, 53), (47, 48)),
            hit_coordinates(
                51,
                1,
                "TGG",
                sequence,
                PrimerDesign.reverseComplement,
            ),
        )
        self.assertEqual(
            ((84, 86), (89, 90)),
            hit_coordinates(
                51,
                -1,
                "AGG",
                sequence,
                PrimerDesign.reverseComplement,
            ),
        )

    def test_preserves_invalid_strand_and_reference_mismatch_errors(self):
        sequence = synthetic_sequence()

        with self.assertRaisesRegex(
            ValueError,
            "Oligo hit strand must be 1 or -1",
        ):
            hit_coordinates(
                51,
                0,
                "TGG",
                sequence,
                PrimerDesign.reverseComplement,
            )
        with self.assertRaisesRegex(
            ValueError,
            "Normalized PAM coordinates do not match the FASTA",
        ):
            hit_coordinates(
                51,
                1,
                "AGG",
                sequence,
                PrimerDesign.reverseComplement,
            )

        with self.assertRaisesRegex(
            ValueError,
            "outside the FASTA",
        ):
            hit_coordinates(
                0,
                1,
                "",
                "GG" + "A" * 30,
                PrimerDesign.reverseComplement,
            )


class TestPrimerDesignCoordinateAdapter(unittest.TestCase):
    def test_coordinate_validation_accepts_the_last_reference_base(self):
        sequence = synthetic_sequence()
        chromosome = chromosomeFasta("synthetic description\n" + sequence)
        designer = PrimerDesign.__new__(PrimerDesign)
        designer.chromosomesData = {"synthetic": chromosome}

        self.assertTrue(
            designer.checkCoords_("synthetic", "1", str(len(sequence)))
        )
        with self.assertRaisesRegex(AssertionError, "Bad chromosomes"):
            designer.checkCoords_(
                "synthetic",
                "1",
                str(len(sequence) + 1),
            )

    def test_adapter_looks_up_chromosome_and_delegates(self):
        sequence = synthetic_sequence()
        chromosome = chromosomeFasta("synthetic description\n" + sequence)
        designer = PrimerDesign.__new__(PrimerDesign)
        designer.chromosomesData = {"synthetic": chromosome}
        hit = NGG("synthetic", 51, 1, PLUS_GUIDE, "TGG")
        sentinel = object()

        with patch(
            "crispr4p.crispr4p.hit_coordinates",
            return_value=sentinel,
        ) as normalize:
            result = designer.getOligoHitCoordinates(hit)

        self.assertIs(sentinel, result)
        normalize.assert_called_once_with(
            hit.pos,
            hit.strand,
            hit.pam,
            chromosome.sequence,
            designer.reverseComplement,
        )

    def test_adapter_preserves_missing_chromosome_error(self):
        designer = PrimerDesign.__new__(PrimerDesign)
        designer.chromosomesData = {}
        hit = NGG("missing", 51, 1, PLUS_GUIDE, "TGG")

        with self.assertRaisesRegex(
            ValueError,
            'Chromosome "missing" not found for oligo hit',
        ):
            designer.getOligoHitCoordinates(hit)


if __name__ == "__main__":
    unittest.main()
