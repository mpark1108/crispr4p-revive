import pickle
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock, patch

from crispr4p.cache import (
    cache_path,
    ensure_cache_dir,
    get_cached,
    save_cache,
)
from crispr4p.crispr4p import PRECOMPUTED_VERSION, PrimerDesign


CACHE_RESULT = (["table"], ("hr",), ["primers"], {"matches": []})


class TestResultCache(unittest.TestCase):
    def test_builds_legacy_gene_and_region_paths(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            cache_directory = Path(temporary_directory) / "cache"
            ensure_cache_dir(cache_directory)

            gene_path = cache_path(
                cache_directory,
                PRECOMPUTED_VERSION,
                0,
                "III",
                "1316337",
                "1317995",
                systematic_name="SPCC1322.13",
            )
            region_path = cache_path(
                cache_directory,
                PRECOMPUTED_VERSION,
                2,
                "III",
                "100",
                "200",
            )

            self.assertTrue(cache_directory.is_dir())
            self.assertEqual(
                "SPCC1322.13_v4_n0.pickle",
                Path(gene_path).name,
            )
            self.assertEqual(
                "III_100_200_v4_n2.pickle",
                Path(region_path).name,
            )

    def test_cold_write_and_warm_load_preserve_four_item_payload(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            cache_path = Path(temporary_directory) / "result.pickle"
            compute = Mock(return_value=CACHE_RESULT)

            cold_result = get_cached(cache_path, compute)
            with cache_path.open("rb") as cache_file:
                stored_result = pickle.load(cache_file)
            warm_result = get_cached(cache_path, compute)

        self.assertEqual(CACHE_RESULT, cold_result)
        self.assertEqual(CACHE_RESULT, stored_result)
        self.assertIsInstance(stored_result, tuple)
        self.assertEqual(CACHE_RESULT, warm_result)
        self.assertEqual(1, compute.call_count)

    def test_bad_warm_cache_is_removed_and_recomputed_without_rewrite(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            cache_path = Path(temporary_directory) / "result.pickle"
            cache_path.write_bytes(pickle.dumps(("only", "three", "items")))
            compute = Mock(return_value=CACHE_RESULT)

            result = get_cached(cache_path, compute)

            self.assertEqual(CACHE_RESULT, result)
            compute.assert_called_once_with()
            self.assertFalse(cache_path.exists())

    def test_failed_save_silently_removes_partial_or_existing_file(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            cache_path = Path(temporary_directory) / "result.pickle"
            cache_path.write_bytes(b"old cache")
            unpickleable_result = (lambda: None, "hr", "primers", "matches")

            save_cache(cache_path, unpickleable_result)

            self.assertFalse(cache_path.exists())

    def test_failed_cold_save_still_returns_computed_result(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            cache_path = Path(temporary_directory) / "result.pickle"
            unpickleable_result = (lambda: None, "hr", "primers", "matches")

            result = get_cached(
                cache_path,
                lambda: unpickleable_result,
            )

            self.assertIs(unpickleable_result[0], result[0])
            self.assertEqual(unpickleable_result[1:], result[1:])
            self.assertFalse(cache_path.exists())

class TestPrimerDesignCacheAdapters(unittest.TestCase):
    def test_legacy_path_adapter_resolves_systematic_name(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            cache_directory = Path(temporary_directory) / "cache"
            designer = PrimerDesign.__new__(PrimerDesign)
            designer.precomputed_folder = cache_directory
            designer.annotationParser_ = SimpleNamespace(
                synonims_=(
                    ("SPCC1322.13", "ade6"),
                ),
            )

            gene_path = designer._genPrecomputedName(
                "ade6",
                0,
                "III",
                "1316337",
                "1317995",
            )
            region_path = designer._genPrecomputedName(
                None,
                2,
                "III",
                "100",
                "200",
            )

            self.assertTrue(cache_directory.is_dir())
            self.assertEqual(
                "SPCC1322.13_v4_n0.pickle",
                Path(gene_path).name,
            )
            self.assertEqual(
                "III_100_200_v4_n2.pickle",
                Path(region_path).name,
            )

    def test_run_remains_validation_and_cache_orchestration_adapter(self):
        designer = PrimerDesign.__new__(PrimerDesign)
        designer.checkCoords_ = Mock()
        designer._genPrecomputedName = Mock(return_value="result.pickle")
        designer._isPrecomputed = Mock(return_value=None)
        designer.run_ = Mock(return_value=CACHE_RESULT)

        with patch(
            "crispr4p.crispr4p.get_cached",
            side_effect=lambda path, compute, is_cached: compute(),
        ) as cached:
            result = designer.run("III", "100", "200", 2, "ade6")

        self.assertEqual(CACHE_RESULT, result)
        designer.checkCoords_.assert_called_once_with("III", "100", "200")
        designer._genPrecomputedName.assert_called_once_with(
            "ade6",
            2,
            "III",
            "100",
            "200",
        )
        designer.run_.assert_called_once_with("III", 100, 200, 2, "ade6")
        cache_path, compute, is_cached = cached.call_args.args
        self.assertEqual("result.pickle", cache_path)
        self.assertTrue(callable(compute))
        self.assertIs(designer._isPrecomputed, is_cached)


if __name__ == "__main__":
    unittest.main()
