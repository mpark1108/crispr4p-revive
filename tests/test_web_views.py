import unittest
import urllib.parse
from io import BytesIO
from pathlib import Path
from unittest.mock import Mock, patch

import webapp
from crispr4p.annotations import AmbiguousGeneNameError, GenomeAnnotations
from crispr4p.models import DesignResult, OligoAnalysisResult, OligoMatch
from crispr4p.resources import GeneNameNotFoundError
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
PRIMER_TUPLE = (
    GUIDE,
    "ACGACGGTCGgttttagagctagaaatagcaagttaaaataa",
    "AAGCCAATGTttcttcggtacaggttatgttttttggcaaca",
    (1316795, 1316797),
    1,
    "TGG",
)
DESIGN_RESULT = DesignResult.from_legacy(
    (
        [[GUIDE, PRIMER_TUPLE, 5, 1, 1, 1, 1, 1, 1]],
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
)


def oligo_result():
    return OligoAnalysisResult(
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


def make_handler():
    return webapp.CRISPR4PHandler.__new__(webapp.CRISPR4PHandler)


class TestErrorViews(unittest.TestCase):
    def test_missing_query_error_names_the_available_inputs(self) -> None:
        result = render_query_error()

        self.assertIn("Name", result)
        self.assertIn("Coordinates", result)
        self.assertIn("Oligo Sequence", result)

    def test_execution_error_includes_the_failure(self) -> None:
        result = render_error(ValueError("failure"))

        self.assertIn("failure", result)

    def test_gene_name_not_found_is_clear_and_escapes_input(self) -> None:
        result = render_gene_error('<missing "gene">')

        self.assertIn("&lt;missing &quot;gene&quot;&gt;", result)
        self.assertIn("try a PomBase systematic gene ID", result)
        self.assertNotIn('<missing "gene">', result)

    def test_ambiguous_gene_name_lists_ids_and_escapes_input(self) -> None:
        error = AmbiguousGeneNameError(
            '<shared "name">',
            ("SPAC1", "SPBC2"),
        )

        result = render_ambiguous_gene_error(error)

        self.assertIn("&lt;shared &quot;name&quot;&gt;", result)
        self.assertIn("SPAC1, SPBC2", result)

    def test_oligo_length_error_reports_the_received_length(self) -> None:
        result = render_oligo_error(19)

        self.assertIn("20 bp", result)
        self.assertIn("23 bp", result)
        self.assertIn("Current length: 19", result)


class TestResultViews(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        data = PROJECT_ROOT / "data"
        annotations = GenomeAnnotations.from_files(
            data / "Schizosaccharomyces_pombe_all_chromosomes.gff3",
            data / "gene_viability.tsv",
        )
        cls.ade6_annotation = annotations.annotate_cut(
            "III",
            (1316791, 1316792),
        )

    def test_oligo_renderer_includes_scientific_results(self) -> None:
        result = render_oligo(oligo_result())

        self.assertIn("<tr><td>8 bp</td><td>5</td></tr>", result)
        self.assertIn("<tr><td>20 bp</td><td>1</td></tr>", result)
        self.assertIn("1316795 - 1316797", result)
        self.assertIn("1316791 | 1316792", result)
        self.assertIn(SPEDIT_FORWARD, result)

    def test_design_renderer_includes_scientific_results(self) -> None:
        template = (
            PROJECT_ROOT / "template" / "container_table.html"
        ).read_text(encoding="utf-8")
        result = render_design(
            DESIGN_RESULT,
            (self.ade6_annotation,),
            template,
        )

        self.assertIn(GUIDE, result)
        self.assertIn("5'-HR_FORWARD-3'", result)
        self.assertIn("LEFT", result)
        self.assertIn('"gene_id":"SPCC1322.13"', result)
        self.assertIn('"viability":"unknown"', result)


class TestWebHandler(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        data = PROJECT_ROOT / "data"
        annotations = GenomeAnnotations.from_files(
            data / "Schizosaccharomyces_pombe_all_chromosomes.gff3",
            data / "gene_viability.tsv",
        )
        cls.annotations = (
            annotations.annotate_cut("III", (1316791, 1316792)),
        )

    def service(self):
        service = Mock()
        service.design_gene.return_value = DESIGN_RESULT
        service.design_region.return_value = DESIGN_RESULT
        service.annotate_guides.return_value = self.annotations
        service.cassette_choices.return_value = ((),)
        service.disruption_donors.return_value = ((),)
        service.restoration_donors.return_value = (None,)
        return service

    def test_web_service_reuses_read_only_genome_data(self) -> None:
        webapp.create_service.cache_clear()
        try:
            first = webapp.create_service()
            second = webapp.create_service()
        finally:
            webapp.create_service.cache_clear()

        self.assertIs(first, second)

    def test_gene_design_calls_each_workflow_service(self) -> None:
        service = self.service()

        with patch.object(webapp, "create_service", return_value=service):
            result = make_handler().run_design_model(
                "ade6",
                None,
                None,
                None,
            )

        service.design_gene.assert_called_once_with("ade6", n_mismatch=0)
        service.annotate_guides.assert_called_once_with(DESIGN_RESULT.guides)
        service.cassette_choices.assert_called_once_with(
            DESIGN_RESULT.guides,
            self.annotations,
            "ade6",
        )
        service.disruption_donors.assert_called_once_with(
            DESIGN_RESULT.guides,
            self.annotations,
            ((),),
            80,
            "ade6",
        )
        service.restoration_donors.assert_called_once_with(
            DESIGN_RESULT.guides,
            self.annotations,
            80,
            "ade6",
        )
        self.assertIn(GUIDE, result)

    def test_coordinate_design_uses_region_service(self) -> None:
        service = self.service()

        with patch.object(webapp, "create_service", return_value=service):
            make_handler().run_design_model(
                None,
                "III",
                "1316337",
                "1317995",
            )

        service.design_region.assert_called_once_with(
            "III",
            "1316337",
            "1317995",
            n_mismatch=0,
        )
        service.design_gene.assert_not_called()


class TestHttpPost(unittest.TestCase):
    def post_handler(self, parameters):
        payload = urllib.parse.urlencode(parameters).encode("utf-8")
        handler = make_handler()
        handler.headers = {"Content-Length": str(len(payload))}
        handler.rfile = BytesIO(payload)
        handler.run_oligo_model = Mock(return_value="OLIGO RESULT")
        handler.run_design_model = Mock(return_value="DESIGN RESULT")
        handler.serve_form = Mock()
        return handler

    def test_oligo_takes_precedence_and_is_normalized(self) -> None:
        handler = self.post_handler(
            {
                "name": "ade6",
                "oligo_sequence": f"  {GUIDE.lower()}  ",
                "oligo_mismatch": "not-an-integer",
            }
        )

        handler.process_post()

        handler.run_oligo_model.assert_called_once_with(GUIDE, 0)
        handler.run_design_model.assert_not_called()
        handler.serve_form.assert_called_once_with("OLIGO RESULT")

    def test_gene_and_coordinate_queries_are_trimmed(self) -> None:
        gene_handler = self.post_handler({"name": "  ade6  "})
        coordinate_handler = self.post_handler(
            {
                "chromosome": " III ",
                "coor_lower": " 1316337 ",
                "coor_upper": " 1317995 ",
            }
        )

        gene_handler.process_post()
        coordinate_handler.process_post()

        gene_handler.run_design_model.assert_called_once_with(
            "ade6", None, None, None
        )
        coordinate_handler.run_design_model.assert_called_once_with(
            None, "III", "1316337", "1317995"
        )

    def test_gene_lookup_errors_give_specific_guidance(self) -> None:
        unknown = self.post_handler({"name": "missing"})
        unknown.run_design_model.side_effect = GeneNameNotFoundError("missing")
        ambiguous = self.post_handler({"name": "noc2"})
        ambiguous.run_design_model.side_effect = AmbiguousGeneNameError(
            "noc2",
            ("SPAC1142.04", "SPAC1B3.09c"),
        )

        unknown.process_post()
        ambiguous.process_post()

        unknown_result = unknown.serve_form.call_args.args[0]
        ambiguous_result = ambiguous.serve_form.call_args.args[0]
        self.assertIn("try a PomBase systematic gene ID", unknown_result)
        self.assertIn("SPAC1142.04, SPAC1B3.09c", ambiguous_result)

    def test_missing_query_returns_guidance_without_dispatch(self) -> None:
        handler = self.post_handler({})

        handler.process_post()

        handler.run_oligo_model.assert_not_called()
        handler.run_design_model.assert_not_called()
        self.assertIn(
            "Please fill either Name, Coordinates, or Oligo Sequence",
            handler.serve_form.call_args.args[0],
        )


if __name__ == "__main__":
    unittest.main()
