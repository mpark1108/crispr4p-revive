import contextlib
import io
import tempfile
import unittest
from pathlib import Path

from crispr4p.annotation_cli import cut_from_pam, format_cut, main
from crispr4p.annotations import GenomeAnnotations, VIABILITY_LABELS

from tests.test_annotations import GFF, PROJECT_GFF, PROJECT_VIABILITY, VIABILITY


class AnnotationCliTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.real_annotations = GenomeAnnotations.from_files(
            PROJECT_GFF,
            PROJECT_VIABILITY,
        )

        cls.temporary_directory = tempfile.TemporaryDirectory()
        directory = Path(cls.temporary_directory.name)
        cls.gff_path = directory / "annotations.gff3"
        cls.viability_path = directory / "gene_viability.tsv"
        cls.gff_path.write_text(GFF, encoding="utf-8")
        cls.viability_path.write_text(VIABILITY, encoding="utf-8")
        cls.synthetic_annotations = GenomeAnnotations.from_files(
            cls.gff_path,
            cls.viability_path,
        )

    @classmethod
    def tearDownClass(cls):
        cls.temporary_directory.cleanup()

    def test_crispr4p_pam_conversion_matches_existing_interface(self):
        self.assertEqual(
            (1316795, 1316797, 1316791, "+"),
            cut_from_pam("1316795 - 1316797", "1"),
        )
        self.assertEqual(
            (1001, 1003, 1006, "-"),
            cut_from_pam("1001 - 1003", "-1"),
        )

        with self.assertRaisesRegex(ValueError, "three bases"):
            cut_from_pam("1000 - 1003", "1")

    def test_only_viable_and_inviable_labels_add_essentiality(self):
        self.assertEqual(
            {
                "viable": "viable (non-essential)",
                "inviable": "inviable (essential)",
                "condition-dependent": "condition-dependent",
                "unknown": "unknown",
                None: "not available",
            },
            dict(VIABILITY_LABELS),
        )

    def test_real_ade6_report_includes_scientific_details(self):
        result = self.real_annotations.annotate_cut(
            "III",
            (1316791, 1316792),
        )
        report = format_cut(result)

        self.assertIn("III:1316791 | III:1316792", report)
        self.assertIn("SPCC1322.13 (ade6)", report)
        self.assertIn("exon (CDS) 1316337-1317995 (1659 bp)", report)
        self.assertIn("CDS: base 455/1659 (27.4%)", report)
        self.assertIn("SPNCRNA.7311", report)
        self.assertIn("non-coding exon 1316304-1317821 (1518 bp)", report)

    def test_feature_boundary_report_retains_both_regions(self):
        result = self.synthetic_annotations.annotate_cut("I", (219, 220))
        report = format_cut(result)

        self.assertIn(
            "cut crosses region boundary: exon (CDS) 120-219 (100 bp) "
            "-> intron 220-299 (80 bp)",
            report,
        )
        self.assertIn("gene viability (PomBase): viable (non-essential)", report)

    def test_intergenic_report_includes_nearest_gene_viability(self):
        result = self.synthetic_annotations.annotate_cut("I", (550, 551))
        report = format_cut(result)

        self.assertIn("I:550 | I:551", report)
        self.assertIn("intergenic on both sides", report)
        self.assertIn(
            "plus (plus_gene) 100-500; 50 bp from cut",
            report,
        )
        self.assertIn(
            "minus (minus_gene) 600-899; 49 bp from cut",
            report,
        )
        self.assertIn("viable (non-essential)", report)
        self.assertIn("inviable (essential)", report)

    def test_main_defaults_to_packaged_data(self):
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            main(
                [
                    "III",
                    "--pam",
                    "1316795 - 1316797",
                    "--strand",
                    "1",
                ]
            )

        report = output.getvalue()
        self.assertTrue(
            report.startswith(
                "true PAM: III:1316795-1316797  guide strand: +\n"
                "Cas9 cut boundary: III:1316791 | III:1316792\n"
            )
        )
        self.assertIn("gene: SPCC1322.13 (ade6)", report)
        self.assertIn(
            "gene viability (PomBase): unknown",
            report,
        )


if __name__ == "__main__":
    unittest.main()
