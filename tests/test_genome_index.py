import unittest

from crispr4p.crispr4p import NGG, PrimerDesign
from crispr4p.genome import GenomePamIndex


GUIDE = "ATGATACATATGATACATAT"
ORIENTED_SITE = "A" + GUIDE + "GGG" + "A"


class TestGenomePamIndex(unittest.TestCase):
    def test_ignores_sites_without_a_complete_protospacer(self) -> None:
        references = {
            "record_start": "GG" + "A" * 30,
            "near_start": "A" * 10 + "TGG" + "A" * 30,
        }

        for name, reference in references.items():
            with self.subTest(name=name):
                index = GenomePamIndex.build(
                    {"synthetic": reference},
                    hit_factory=NGG,
                )

                self.assertEqual(0, index.hit_count)
                self.assertEqual(0, len(index))

    def test_preserves_overlapping_pams_on_both_strands(self) -> None:
        references = {
            1: ORIENTED_SITE,
            -1: PrimerDesign.reverseComplement(ORIENTED_SITE),
        }

        for expected_strand, reference in references.items():
            with self.subTest(strand=expected_strand):
                index = GenomePamIndex.build(
                    {"synthetic": reference},
                    hit_factory=NGG,
                )
                matching_hits = [
                    hit
                    for hit in index.get(GUIDE[-8:], ())
                    if hit.seed == GUIDE and hit.pam == "GGG"
                ]

                self.assertEqual(1, len(matching_hits))
                self.assertEqual(expected_strand, matching_hits[0].strand)

    def test_mapping_and_buckets_are_read_only(self) -> None:
        index = GenomePamIndex.build(
            {"synthetic": ORIENTED_SITE},
            hit_factory=NGG,
        )
        suffix = GUIDE[-8:]

        with self.assertRaises(TypeError):
            index.by_suffix[suffix] = ()
        with self.assertRaises(AttributeError):
            index[suffix].append("hit")

    def test_index_build_is_ordered_and_idempotent_in_primer_design(self) -> None:
        designer = PrimerDesign.__new__(PrimerDesign)
        designer.chromosomesData = {
            "first": type("Record", (), {"sequence": ORIENTED_SITE})(),
            "second": type("Record", (), {"sequence": ORIENTED_SITE})(),
        }
        designer.genome_index = None

        first_index = designer.getNGGsFromGenome()
        second_index = designer.getNGGsFromGenome()

        self.assertIs(first_index, second_index)
        self.assertIs(first_index.by_suffix, designer.NGGs)
        chromosomes = [
            hit.chromosome
            for hit in first_index[GUIDE[-8:]]
            if hit.seed == GUIDE and hit.pam == "GGG"
        ]
        self.assertEqual(["first", "second"], chromosomes)


if __name__ == "__main__":
    unittest.main()
