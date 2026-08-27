import unittest
from pathlib import Path

from crispr4p.disruption import load_cassettes
from crispr4p.resources import read_fasta
from crispr4p.restoration import build_donor
from crispr4p.service import Crispr4pService
from crispr4p.spedit import has_bsai, make_oligos, reverse_complement
from crispr4p.web_views import cassette_data, render_design, restoration_data


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
FASTA = DATA / "Schizosaccharomyces_pombe.ASM294v2.26.dna.toplevel.fa"


class RestorationDesignTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.reference = read_fasta(FASTA)["III"].sequence
        cls.cassettes = load_cassettes(DATA / "stop_cassettes.json")

    def test_rescue_guide(self):
        cassette = self.cassettes[0]
        forward, reverse = make_oligos(cassette.guide)

        self.assertEqual("TGATGACGTGATGACCTAGT", cassette.guide)
        self.assertEqual("AGG", cassette.pam)
        self.assertEqual((17, 18), cassette.cut_bases)
        self.assertEqual(52, len(forward))
        self.assertEqual(reverse_complement(forward), reverse)
        self.assertFalse(has_bsai(cassette.guide))

    def test_restoration_donor(self):
        cut = (1316791, 1316792)
        donor = build_donor(self.reference, cut, 80)

        self.assertEqual(self.reference[1316711:1316791], donor.left_arm)
        self.assertEqual(self.reference[1316791:1316871], donor.right_arm)
        self.assertEqual(self.reference[1316711:1316871], donor.sequence)
        self.assertEqual(reverse_complement(donor.sequence), donor.reverse)
        self.assertEqual(160, donor.total_length)

        oligos = donor.oligos
        self.assertEqual(donor.left_arm + donor.right_arm[:20], oligos.forward)
        self.assertEqual(
            reverse_complement(donor.left_arm[-20:] + donor.right_arm),
            oligos.reverse,
        )
        self.assertEqual(
            donor.left_arm[-20:] + donor.right_arm[:20],
            oligos.overlap,
        )
        self.assertEqual(
            reverse_complement(oligos.overlap),
            oligos.reverse[-len(oligos.overlap):],
        )
        self.assertEqual(
            (100, 100, 40),
            (
                len(oligos.forward),
                len(oligos.reverse),
                len(oligos.overlap),
            ),
        )
        self.assertEqual(donor.sequence, oligos.product)
        self.assertEqual(160, oligos.product_length)

        edited = donor.left_arm + self.cassettes[0].sequence + donor.right_arm
        self.assertEqual(donor.sequence, edited[:80] + edited[103:])
        self.assertIn("ACATTGGCTTACGACGGTCGTGG", donor.sequence)
        self.assertNotIn(self.cassettes[0].sequence, donor.sequence)

        with self.assertRaisesRegex(ValueError, "adjacent bases"):
            build_donor(self.reference, (100, 102), 80)
        with self.assertRaisesRegex(ValueError, "positive integer"):
            build_donor(self.reference, cut, 0)
        with self.assertRaisesRegex(ValueError, "complete homology arms"):
            build_donor(self.reference, (40, 41), 80)

    def test_service_builds_aligned_donors(self):
        service = Crispr4pService.from_project_data()
        result = service.design_gene("ade6")
        annotations = service.annotate_guides(result.guides)
        donors = service.restoration_donors(
            result.guides,
            annotations,
            80,
            result.name,
        )

        self.assertEqual(len(result.guides), len(donors))
        self.assertEqual(
            self.reference[
                result.guides[0].cut_coordinates[0] - 80:
                result.guides[0].cut_coordinates[0] + 80
            ],
            donors[0].sequence,
        )

        bub1 = service.design_gene("bub1")
        bub1_annotations = service.annotate_guides(bub1.guides)
        bub1_donor = service.restoration_donors(
            bub1.guides,
            bub1_annotations,
            80,
            bub1.name,
        )[0]
        bub1_cut = bub1.guides[0].cut_coordinates[0]
        self.assertEqual(
            self.reference[bub1_cut - 80:bub1_cut + 80],
            bub1_donor.sequence,
        )

        noncoding = service.design_gene("SPNCRNA.1438")
        noncoding_annotations = service.annotate_guides(noncoding.guides)
        self.assertTrue(
            all(
                donor is None
                for donor in service.restoration_donors(
                    noncoding.guides,
                    noncoding_annotations,
                    80,
                    noncoding.name,
                )
            )
        )

        with self.assertRaisesRegex(ValueError, "counts must match"):
            service.restoration_donors(result.guides, (), 80, result.name)

    def test_web_section(self):
        service = Crispr4pService.from_project_data()
        result = service.design_gene("ade6")
        annotations = service.annotate_guides(result.guides)
        choices = service.cassette_choices(
            result.guides,
            annotations,
            result.name,
        )
        disruption = service.disruption_donors(
            result.guides,
            annotations,
            choices,
            80,
            result.name,
        )
        restoration = service.restoration_donors(
            result.guides,
            annotations,
            80,
            result.name,
        )
        template = (ROOT / "template/container_table.html").read_text(
            encoding="utf-8"
        )
        page = render_design(
            result,
            annotations,
            template,
            cassette_choices=choices,
            disruption_donors=disruption,
            restoration_donors=restoration,
        )

        cassette = cassette_data((choices[0],))["catalog"]["1"]
        donor = restoration_data((restoration[0],))[0]
        self.assertEqual(52, len(cassette["spedit_forward"]))
        self.assertEqual(52, len(cassette["spedit_reverse"]))
        self.assertFalse(cassette["has_internal_bsai"])
        self.assertEqual(80, donor["arm_length"])
        self.assertEqual(160, len(donor["left_arm"] + donor["right_arm"]))
        self.assertEqual(100, donor["hrfw_length"])
        self.assertEqual(100, donor["hrrv_length"])
        self.assertEqual(40, donor["overlap_length"])
        self.assertEqual(160, donor["hr_product_length"])

        self.assertIn("Wild-Type Restoration Design", page)
        self.assertIn('id="restoration_guide"', page)
        self.assertIn('id="restoration_sequence"', page)
        self.assertIn('id="restoration_hrfw"', page)
        self.assertIn('id="restoration_product_size"', page)
        self.assertIn(cassette["spedit_forward"], page)
        self.assertIn(donor["left_arm"], page)


if __name__ == "__main__":
    unittest.main()
