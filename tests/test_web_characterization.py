import hashlib
import unittest
import urllib.parse
from io import BytesIO
from pathlib import Path
from unittest.mock import Mock, patch

import webapp
from crispr4p.annotations import AmbiguousGeneNameError, GenomeAnnotations
from crispr4p.models import DesignResult, OligoAnalysisResult, OligoMatch
from crispr4p.resources import GeneNameNotFoundError
from crispr4p.service import OligoLengthError


GUIDE = "ACATTGGCTTACGACGGTCG"
PROJECT_ROOT = Path(__file__).resolve().parents[1]
SPEDIT_FORWARD = (
    "CTAGAGGTCTCGGACTACATTGGCTTACGACGGTCGGTTTCGAGACCCTTCC"
)


def make_handler():
    """Construct a handler without opening a network socket."""
    return webapp.CRISPR4PHandler.__new__(webapp.CRISPR4PHandler)


PRIMER_TUPLE = (
    GUIDE,
    "ACGACGGTCGgttttagagctagaaatagcaagttaaaataa",
    "AAGCCAATGTttcttcggtacaggttatgttttttggcaaca",
    (1316795, 1316797),
    1,
    "TGG",
)
DESIGN_TABLE = [[GUIDE, PRIMER_TUPLE, 5, 1, 1, 1, 1, 1, 1]]
HR_DNA = ("HR_FORWARD", "HR_REVERSE", "DELETED_DNA")
CHECKING_PRIMERS = [
    {
        "PRIMER_LEFT_0_SEQUENCE": "LEFT",
        "PRIMER_LEFT_0_TM": 59.4,
        "PRIMER_RIGHT_0_SEQUENCE": "RIGHT",
        "PRIMER_RIGHT_0_TM": 60.1,
        "PRIMER_PAIR_0_PRODUCT_SIZE": 255,
        "negative_result": 1913,
    }
]
DESIGN_RESULT = DesignResult.from_legacy(
    (
        DESIGN_TABLE,
        HR_DNA,
        CHECKING_PRIMERS,
        "ade6",
        "III",
        "1316337",
        "1317995",
    )
)
ANNOTATION_INDEX = GenomeAnnotations.from_files(
    PROJECT_ROOT / "data" / "Schizosaccharomyces_pombe_all_chromosomes.gff3",
    PROJECT_ROOT / "data" / "gene_viability.tsv",
)
DESIGN_ANNOTATIONS = (
    ANNOTATION_INDEX.annotate_cut("III", (1316791, 1316792)),
)
OLIGO_RESULT = OligoAnalysisResult(
    oligo_sequence=GUIDE,
    seed=GUIDE,
    pam="NGG",
    n_mismatch=0,
    spedit_forward=SPEDIT_FORWARD,
    spedit_reverse=(
        "GGAAGGGTCTCGAAACCGACCGTCGTAAGCCAATGTAGTCCGAGACCTCTAG"
    ),
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


class TestWebRenderingCharacterization(unittest.TestCase):
    def test_gene_query_label_accepts_name_or_id(self) -> None:
        template = (
            PROJECT_ROOT / "template" / "bahler_template.html"
        ).read_text(encoding="utf-8")

        self.assertIn("Gene name/ID:", template)
        self.assertNotIn("Gene name:", template)

    def test_web_service_is_reused_for_read_only_genome_data(self) -> None:
        webapp.create_service.cache_clear()
        try:
            first = webapp.create_service()
            second = webapp.create_service()
        finally:
            webapp.create_service.cache_clear()

        self.assertIs(first, second)

    def test_oligo_result_html_is_unchanged(self) -> None:
        handler = make_handler()
        service = Mock()
        service.analyze_oligo.return_value = OLIGO_RESULT

        with patch.object(webapp, "create_service", return_value=service):
            result = handler.run_oligo_model(GUIDE, 0)

        service.analyze_oligo.assert_called_once_with(GUIDE, n_mismatch=0)
        self.assertEqual(
            "dd389db45eaca6b3babe396638609b0752b865f5638c14c32aa6c0e7e4c4082e",
            hashlib.sha256(result.encode("utf-8")).hexdigest(),
        )
        self.assertIn("<tr><td>8 bp</td><td>5</td></tr>", result)
        self.assertIn("<tr><td>20 bp</td><td>1</td></tr>", result)
        self.assertIn("1316795 - 1316797", result)
        self.assertIn("1316791 | 1316792", result)
        self.assertIn(SPEDIT_FORWARD, result)

    def test_design_result_html_includes_annotation_workflows(self) -> None:
        handler = make_handler()
        service = Mock()
        service.design_gene.return_value = DESIGN_RESULT
        service.annotate_guides.return_value = DESIGN_ANNOTATIONS
        service.cassette_choices.return_value = ((),)
        service.disruption_donors.return_value = ((),)
        service.restoration_donors.return_value = (None,)

        with patch.object(webapp, "create_service", return_value=service):
            result = handler.run_design_model("ade6", None, None, None)

        service.design_gene.assert_called_once_with("ade6", n_mismatch=0)
        service.design_region.assert_not_called()
        service.annotate_guides.assert_called_once_with(DESIGN_RESULT.guides)
        service.disruption_donors.assert_called_once_with(
            DESIGN_RESULT.guides,
            DESIGN_ANNOTATIONS,
            ((),),
            80,
            "ade6",
        )
        service.restoration_donors.assert_called_once_with(
            DESIGN_RESULT.guides,
            DESIGN_ANNOTATIONS,
            80,
            "ade6",
        )
        self.assertIn("<b>Name</b>=ade6", result)
        self.assertIn(GUIDE, result)
        self.assertIn("5'-HR_FORWARD-3'", result)
        self.assertIn("LEFT", result)
        self.assertIn("59 &deg;C", result)
        self.assertIn("255 (bp)", result)
        self.assertIn(SPEDIT_FORWARD, result)
        self.assertIn("Selected sgRNA and Cut Site", result)
        self.assertIn("More cut-site details", result)
        self.assertLess(
            result.index("Deletion Design &mdash; Original CRISPR4P"),
            result.index("SpEDIT/pLSB Golden Gate Cloning Oligos"),
        )
        self.assertIn('"gene_id":"SPCC1322.13"', result)
        self.assertIn('"gene_id":"SPNCRNA.7311"', result)
        self.assertIn('"role":"Primary target"', result)
        self.assertIn('"region":"coding sequence (CDS)"', result)
        self.assertIn('"Gene viability"', result)
        self.assertNotIn("Gene viability (PomBase)", result)
        self.assertNotIn("cut also overlaps a protein-coding gene", result)
        self.assertNotIn("annotation_attention", result)

    def test_coordinate_design_uses_region_service_operation(self) -> None:
        handler = make_handler()
        service = Mock()
        service.design_region.return_value = DESIGN_RESULT
        service.annotate_guides.return_value = DESIGN_ANNOTATIONS
        service.cassette_choices.return_value = ((),)
        service.disruption_donors.return_value = ((),)
        service.restoration_donors.return_value = (None,)

        with patch.object(webapp, "create_service", return_value=service):
            handler.run_design_model(None, "III", "1316337", "1317995")

        service.design_region.assert_called_once_with(
            "III",
            "1316337",
            "1317995",
            n_mismatch=0,
        )
        service.design_gene.assert_not_called()
        service.annotate_guides.assert_called_once_with(DESIGN_RESULT.guides)

    def test_oligo_length_error_keeps_existing_html(self) -> None:
        handler = make_handler()
        service = Mock()
        service.analyze_oligo.side_effect = OligoLengthError(19)

        with patch.object(webapp, "create_service", return_value=service):
            result = handler.run_oligo_model("A" * 19, 0)

        self.assertEqual(
            '<font color="red"><h3>Error: Oligo sequence must be 20 bp '
            '(seed only) or 23 bp (seed + PAM). Current length: 19</h3></font>',
            result,
        )


class TestHttpPostCharacterization(unittest.TestCase):
    def make_post_handler(self, parameters):
        payload = urllib.parse.urlencode(parameters).encode("utf-8")
        handler = make_handler()
        handler.headers = {"Content-Length": str(len(payload))}
        handler.rfile = BytesIO(payload)
        handler.run_oligo_model = Mock(return_value="OLIGO RESULT")
        handler.run_design_model = Mock(return_value="DESIGN RESULT")
        handler.serve_form = Mock()
        return handler

    def test_oligo_takes_precedence_and_input_is_normalized(self) -> None:
        handler = self.make_post_handler(
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

    def test_gene_query_is_trimmed_and_forwarded(self) -> None:
        handler = self.make_post_handler({"name": "  ade6  "})

        handler.process_post()

        handler.run_design_model.assert_called_once_with(
            "ade6", None, None, None
        )
        handler.run_oligo_model.assert_not_called()
        handler.serve_form.assert_called_once_with("DESIGN RESULT")

    def test_unknown_gene_gets_specific_user_guidance(self) -> None:
        handler = self.make_post_handler({"name": "  missing  "})
        handler.run_design_model.side_effect = GeneNameNotFoundError("missing")

        handler.process_post()

        result = handler.serve_form.call_args.args[0]
        self.assertIn('Gene &quot;missing&quot; was not found.', result)
        self.assertIn("try a PomBase systematic gene ID", result)
        self.assertNotIn("ERROR during execution", result)

    def test_ambiguous_gene_gets_systematic_id_choices(self) -> None:
        handler = self.make_post_handler({"name": "noc2"})
        handler.run_design_model.side_effect = AmbiguousGeneNameError(
            "noc2",
            ("SPAC1142.04", "SPAC1B3.09c"),
        )

        handler.process_post()

        result = handler.serve_form.call_args.args[0]
        self.assertIn("Gene name &quot;noc2&quot; is ambiguous.", result)
        self.assertIn("SPAC1142.04, SPAC1B3.09c", result)
        self.assertNotIn("ERROR during execution", result)

    def test_coordinate_query_is_trimmed_and_forwarded(self) -> None:
        handler = self.make_post_handler(
            {
                "chromosome": " III ",
                "coor_lower": " 1316337 ",
                "coor_upper": " 1317995 ",
            }
        )

        handler.process_post()

        handler.run_design_model.assert_called_once_with(
            None, "III", "1316337", "1317995"
        )
        handler.run_oligo_model.assert_not_called()
        handler.serve_form.assert_called_once_with("DESIGN RESULT")

    def test_missing_query_keeps_existing_error_message(self) -> None:
        handler = self.make_post_handler({})

        handler.process_post()

        handler.run_oligo_model.assert_not_called()
        handler.run_design_model.assert_not_called()
        handler.serve_form.assert_called_once_with(
            '<font color="red"><h3>Error: Please fill either Name, '
            'Coordinates, or Oligo Sequence</h3></font>'
        )


if __name__ == "__main__":
    unittest.main()
