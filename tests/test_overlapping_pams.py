import unittest
from types import SimpleNamespace

from crispr4p.crispr4p import PrimerDesign


GUIDE = "ATGATACATATGATACATAT"
ORIENTED_SITE = "A" + GUIDE + "GGG" + "A"


def make_designer(reference_sequence):
    """Build the small part of PrimerDesign needed for PAM scanning tests."""
    designer = PrimerDesign.__new__(PrimerDesign)
    designer.chromosomesData = {
        "synthetic": SimpleNamespace(sequence=reference_sequence)
    }
    designer.userNGGs = []
    return designer


class TestOverlappingPams(unittest.TestCase):
    def test_region_scan_finds_ggg_pam_on_both_strands(self):
        references = {
            1: ORIENTED_SITE,
            -1: PrimerDesign.reverseComplement(ORIENTED_SITE),
        }

        for expected_strand, reference in references.items():
            with self.subTest(strand=expected_strand):
                designer = make_designer(reference)
                chromosome = SimpleNamespace(
                    name="synthetic",
                    sequence=reference,
                )

                designer._getUserNGGs(chromosome, 1, len(reference))

                self.assertTrue(
                    any(
                        hit.seed == GUIDE
                        and hit.pam == "GGG"
                        and hit.strand == expected_strand
                        for hit in designer.userNGGs
                    )
                )

    def test_genome_index_finds_ggg_pam_on_both_strands(self):
        references = {
            1: ORIENTED_SITE,
            -1: PrimerDesign.reverseComplement(ORIENTED_SITE),
        }

        for expected_strand, reference in references.items():
            with self.subTest(strand=expected_strand):
                designer = make_designer(reference)

                designer.getNGGsFromGenome()
                matching_hits = [
                    hit
                    for hit in designer.NGGs.get(GUIDE[-8:], ())
                    if hit.seed == GUIDE and hit.pam == "GGG"
                ]

                self.assertEqual(1, len(matching_hits))
                self.assertEqual(expected_strand, matching_hits[0].strand)


if __name__ == "__main__":
    unittest.main()
