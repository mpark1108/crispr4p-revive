import subprocess
import sys
import unittest
from io import StringIO
from pathlib import Path
from unittest.mock import Mock

from crispr4p.cli import (
    format_design,
    format_oligo,
    run_cli,
)
from crispr4p.crispr4p import PrimerDesign
from crispr4p.models import DesignResult, OligoAnalysisResult, OligoMatch
from crispr4p.service import OligoLengthError


PROJECT_ROOT = Path(__file__).resolve().parents[1]
GUIDE = "ACATTGGCTTACGACGGTCG"
PRIMER_TUPLE = (
    GUIDE,
    "FORWARD",
    "REVERSE",
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
OLIGO_RESULT = OligoAnalysisResult(
    oligo_sequence=GUIDE,
    seed=GUIDE,
    pam="NGG",
    n_mismatch=0,
    spedit_forward="forward",
    spedit_reverse="reverse",
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


class TestCliRendering(unittest.TestCase):
    def test_design_renderer_keeps_legacy_fields(self) -> None:
        output = format_design(DESIGN_RESULT)

        self.assertIn(f"1 - {GUIDE} [5, 1, 1, 1, 1, 1, 1]", output)
        self.assertIn("PAM: 1316795 - 1316797 TGG 1", output)
        self.assertIn("HRfw:  HR_FORWARD", output)
        self.assertIn("Check primer left:  LEFT TM: 59.4", output)
        self.assertIn("Negative result product size:  1913", output)

    def test_oligo_renderer_keeps_legacy_fields(self) -> None:
        output = format_oligo(
            OLIGO_RESULT,
            include_loading_messages=False,
        )

        self.assertIn(f"oligo seed: {GUIDE} with PAM: NGG", output)
        self.assertIn("           8 | 5", output)
        self.assertIn("          20 | 1", output)
        self.assertIn("PAM coordinates: 1316795 - 1316797", output)
        self.assertIn("Cut: 1316791 | 1316792", output)


class TestCliDispatch(unittest.TestCase):
    def test_gene_query_uses_design_gene(self) -> None:
        service = Mock()
        service.design_gene.return_value = DESIGN_RESULT
        output = StringIO()

        status = run_cli(
            ["--name", "ade6", "--mismatch", "2"],
            service,
            output,
        )

        self.assertEqual(0, status)
        service.design_gene.assert_called_once_with("ade6", n_mismatch=2)
        service.design_region.assert_not_called()
        self.assertIn(GUIDE, output.getvalue())

    def test_coordinate_query_uses_design_region(self) -> None:
        service = Mock()
        service.design_region.return_value = DESIGN_RESULT
        output = StringIO()

        status = run_cli(
            [
                "--chromosome",
                "III",
                "--coords",
                "1316337...1317995",
            ],
            service,
            output,
        )

        self.assertEqual(0, status)
        service.design_region.assert_called_once_with(
            "III",
            "1316337",
            "1317995",
            n_mismatch=0,
        )
        service.design_gene.assert_not_called()

    def test_oligo_query_takes_precedence(self) -> None:
        service = Mock()
        output = StringIO()

        def analyze_oligo(*args, **kwargs):
            # Preserve the interactive legacy behavior: users see the loading
            # status before the potentially slow genome index is constructed.
            self.assertIn(
                "Loading genome data and indexing PAM sites...",
                output.getvalue(),
            )
            return OLIGO_RESULT

        service.analyze_oligo.side_effect = analyze_oligo

        status = run_cli(
            ["--name", "ade6", "--oligo", GUIDE],
            service,
            output,
        )

        self.assertEqual(0, status)
        service.analyze_oligo.assert_called_once_with(GUIDE, n_mismatch=0)
        service.design_gene.assert_not_called()
        self.assertIn("Loading genome data and indexing PAM sites", output.getvalue())

    def test_invalid_oligo_length_keeps_error_and_status(self) -> None:
        service = Mock()
        output = StringIO()

        status = run_cli(["--oligo", "A" * 19], service, output)

        self.assertEqual(1, status)
        self.assertIn("20 bp", output.getvalue())
        self.assertIn("23 bp", output.getvalue())
        self.assertIn("Received length: 19", output.getvalue())
        service.analyze_oligo.assert_not_called()

    def test_missing_query_keeps_help_and_exit(self) -> None:
        output = StringIO()

        with self.assertRaises(SystemExit) as raised:
            run_cli([], Mock(), output)

        self.assertIsNone(raised.exception.code)
        self.assertIn("crispr4p description", output.getvalue())
        self.assertIn("--oligo OLIGO", output.getvalue())


class TestCliEntrypoints(unittest.TestCase):
    def test_scientific_engine_does_not_build_cli_parser_eagerly(self) -> None:
        data_directory = PROJECT_ROOT / "data"
        designer = PrimerDesign(
            data_directory
            / "Schizosaccharomyces_pombe.ASM294v2.26.dna.toplevel.fa",
            data_directory / "COORDINATES.txt",
            data_directory / "SYNONIMS.txt",
        )

        self.assertFalse(hasattr(designer, "argp_"))

    def test_direct_legacy_script_invocation_reaches_new_cli(self) -> None:
        completed = subprocess.run(
            [
                sys.executable,
                "crispr4p/crispr4p.py",
                "--help",
            ],
            cwd=PROJECT_ROOT,
            text=True,
            capture_output=True,
            check=False,
        )

        self.assertEqual(0, completed.returncode, completed.stderr)
        self.assertIn("crispr4p description", completed.stdout)
        self.assertIn("--oligo OLIGO", completed.stdout)
        self.assertNotIn("run time", completed.stdout)

    def test_module_invocation_reaches_new_cli(self) -> None:
        completed = subprocess.run(
            [sys.executable, "-m", "crispr4p.cli", "--help"],
            cwd=PROJECT_ROOT,
            text=True,
            capture_output=True,
            check=False,
        )

        self.assertEqual(0, completed.returncode, completed.stderr)
        self.assertIn("crispr4p description", completed.stdout)
        self.assertIn("--oligo OLIGO", completed.stdout)


if __name__ == "__main__":
    unittest.main()
