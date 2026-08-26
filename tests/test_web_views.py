import hashlib
import unittest
from pathlib import Path

from crispr4p.annotations import AmbiguousGeneNameError, GenomeAnnotations
from crispr4p.models import DesignResult, OligoAnalysisResult, OligoMatch
from crispr4p.web_views import (
    render_ambiguous_gene_error,
    render_design,
    render_error,
    render_gene_error,
    render_oligo,
    render_oligo_error,
    render_query_error,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]
GUIDE = "ACATTGGCTTACGACGGTCG"
SPEDIT_FORWARD = (
    "CTAGAGGTCTCGGACTACATTGGCTTACGACGGTCGGTTTCGAGACCCTTCC"
)
SPEDIT_REVERSE = (
    "GGAAGGGTCTCGAAACCGACCGTCGTAAGCCAATGTAGTCCGAGACCTCTAG"
)


class TestErrorViews(unittest.TestCase):
    def test_missing_query_error_is_unchanged(self) -> None:
        self.assertEqual(
            '<font color="red"><h3>Error: Please fill either Name, '
            'Coordinates, or Oligo Sequence</h3></font>',
            render_query_error(),
        )

    def test_execution_error_is_unchanged(self) -> None:
        self.assertEqual(
            '<font color="red"><h3>ERROR during execution: failure</h3></font>',
            render_error(ValueError("failure")),
        )

    def test_gene_name_not_found_is_clear_and_escapes_input(self) -> None:
        result = render_gene_error('<missing "gene">')

        self.assertEqual(
            '<font color="red"><h3>Gene &quot;&lt;missing &quot;gene&quot;&gt;'
            '&quot; was not found. Check the spelling or try a PomBase '
            'systematic gene ID, such as SPBPB2B2.01.</h3></font>',
            result,
        )
        self.assertNotIn("ERROR during execution", result)

    def test_ambiguous_gene_name_lists_ids_and_escapes_input(self) -> None:
        error = AmbiguousGeneNameError(
            '<shared "name">',
            ("SPAC1", "SPBC2"),
        )

        result = render_ambiguous_gene_error(error)

        self.assertIn("&lt;shared &quot;name&quot;&gt;", result)
        self.assertIn("SPAC1, SPBC2", result)
        self.assertNotIn("ERROR during execution", result)

    def test_oligo_length_error_is_unchanged(self) -> None:
        self.assertEqual(
            '<font color="red"><h3>Error: Oligo sequence must be 20 bp '
            '(seed only) or 23 bp (seed + PAM). Current length: 19</h3></font>',
            render_oligo_error(19),
        )


class TestResultViews(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        data_directory = PROJECT_ROOT / "data"
        annotations = GenomeAnnotations.from_files(
            data_directory / "Schizosaccharomyces_pombe_all_chromosomes.gff3",
            data_directory / "gene_viability.tsv",
        )
        cls.ade6_annotation = annotations.annotate_cut(
            "III",
            (1316791, 1316792),
        )

    def test_oligo_renderer_preserves_exact_html(self) -> None:
        analysis = OligoAnalysisResult(
            oligo_sequence=GUIDE,
            seed=GUIDE,
            pam="NGG",
            n_mismatch=0,
            spedit_forward=SPEDIT_FORWARD,
            spedit_reverse=SPEDIT_REVERSE,
            has_internal_bsai=False,
            match_counts={8: 5, 10: 1, 12: 1, 14: 1, 16: 1, 18: 1, 20: 1},
            full_matches=[
            OligoMatch(
                    chromosome="III",
                    pam_coordinates=(1316795, 1316797),
                    cut_coordinates=(1316791, 1316792),
                    strand=1,
                    seed=GUIDE,
                    pam="TGG",
                )
            ],
        )
        result = render_oligo(analysis)

        self.assertEqual(
            "dd389db45eaca6b3babe396638609b0752b865f5638c14c32aa6c0e7e4c4082e",
            hashlib.sha256(result.encode("utf-8")).hexdigest(),
        )

    def test_design_renderer_preserves_exact_html(self) -> None:
        primer_tuple = (
            GUIDE,
            "ACGACGGTCGgttttagagctagaaatagcaagttaaaataa",
            "AAGCCAATGTttcttcggtacaggttatgttttttggcaaca",
            (1316795, 1316797),
            1,
            "TGG",
        )
        legacy_result = (
            [[GUIDE, primer_tuple, 5, 1, 1, 1, 1, 1, 1]],
            ("HR_FORWARD", "HR_REVERSE", "DELETED_DNA"),
            [
                {
                    "PRIMER_LEFT_0_SEQUENCE": "LEFT",
                    "PRIMER_LEFT_0_TM": 59.4,
                    "PRIMER_RIGHT_0_SEQUENCE": "RIGHT",
                    "PRIMER_RIGHT_0_TM": 60.1,
                    "PRIMER_PAIR_0_PRODUCT_SIZE": 255,
                    "negative_result": 1913,
                }
            ],
            "ade6",
            "III",
            "1316337",
            "1317995",
        )
        template_text = (
            PROJECT_ROOT / "template" / "container_table.html"
        ).read_text(encoding="utf-8")

        result = render_design(
            DesignResult.from_legacy(legacy_result),
            (self.ade6_annotation,),
            template_text,
        )

        self.assertEqual(
            "b9021533ac8ba35e1b26b4e2f18d1f43085ee46aba086a5fef4c7d7ab0033a12",
            hashlib.sha256(result.encode("utf-8")).hexdigest(),
        )
        self.assertIn("Selected sgRNA and Cut Site", result)
        self.assertIn("Deletion Design &mdash; Original CRISPR4P", result)
        self.assertIn("SpEDIT/pLSB Golden Gate Cloning Oligos", result)
        self.assertIn('"gene_id":"SPCC1322.13"', result)
        self.assertIn('"viability":"unknown"', result)


if __name__ == "__main__":
    unittest.main()
