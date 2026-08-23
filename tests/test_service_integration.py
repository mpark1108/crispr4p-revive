import hashlib
import pickle
import tempfile
import unittest
from pathlib import Path

from crispr4p.cli import format_design, format_oligo
from crispr4p.service import Crispr4pService
from crispr4p.web_views import (
    annotation_rows,
    render_design,
    render_oligo,
)


GUIDE = "ACATTGGCTTACGACGGTCG"
PROJECT_ROOT = Path(__file__).resolve().parents[1]


class TestServiceIntegration(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        service = Crispr4pService.from_project_data(
            precomputed_folder="precomputed",
        )
        cls.service = service
        cls.result = service.analyze_oligo(GUIDE, n_mismatch=0)
        cls.design_result = service.design_gene("ade6", n_mismatch=0)
        cls.annotations_initially_unloaded = service.genome_annotations is None

    def test_real_ade6_oligo_analysis_is_unchanged(self) -> None:
        self.assertEqual(
            {8: 5, 10: 1, 12: 1, 14: 1, 16: 1, 18: 1, 20: 1},
            dict(self.result.match_counts),
        )
        self.assertEqual(1, len(self.result.full_matches))
        match = self.result.full_matches[0]
        self.assertEqual("III", match.chromosome)
        self.assertEqual((1316795, 1316797), match.pam_coordinates)
        self.assertEqual((1316791, 1316792), match.cut_coordinates)
        self.assertEqual(1, match.strand)
        self.assertEqual(GUIDE, match.seed)
        self.assertEqual("TGG", match.pam)
        self.assertEqual(
            "dd389db45eaca6b3babe396638609b0752b865f5638c14c32aa6c0e7e4c4082e",
            hashlib.sha256(
                render_oligo(self.result).encode("utf-8")
            ).hexdigest(),
        )

    def test_real_ade6_design_service_result_is_unchanged(self) -> None:
        self.assertEqual("ade6", self.design_result.name)
        self.assertEqual("III", self.design_result.chromosome)
        self.assertEqual("1316337", self.design_result.start)
        self.assertEqual("1317995", self.design_result.end)
        self.assertEqual(144, len(self.design_result.guide_table))
        self.assertEqual(GUIDE, self.design_result.guide_table[0][0])
        self.assertEqual(
            [5, 1, 1, 1, 1, 1, 1],
            self.design_result.guide_table[0][2:],
        )

    def test_real_ade6_design_exposes_typed_guide_candidates(self) -> None:
        self.assertEqual(144, len(self.design_result.guides))
        guide = self.design_result.guides[0]

        self.assertEqual("III", guide.chromosome)
        self.assertEqual(GUIDE, guide.seed)
        self.assertEqual((1316795, 1316797), guide.pam_coordinates)
        self.assertEqual((1316791, 1316792), guide.cut_coordinates)
        self.assertEqual(1, guide.strand)
        self.assertEqual("TGG", guide.pam)
        self.assertEqual(
            {8: 5, 10: 1, 12: 1, 14: 1, 16: 1, 18: 1, 20: 1},
            dict(guide.match_counts),
        )
        self.assertIs(
            self.design_result.guide_table[0],
            guide.to_legacy(),
        )

    def test_real_ade6_guide_annotation_uses_shared_packaged_index(self) -> None:
        first_guide = self.design_result.guides[0]

        self.assertTrue(self.annotations_initially_unloaded)
        result = self.service.annotate_guide(first_guide)
        shared_index = self.service.genome_annotations
        repeated = self.service.annotate_guides((first_guide,))[0]

        self.assertIsNotNone(shared_index)
        self.assertIs(shared_index, self.service.genome_annotations)
        self.assertEqual(result, repeated)
        self.assertEqual((1316791, 1316792), result.cut_coordinates)
        self.assertEqual(
            {"SPCC1322.13", "SPNCRNA.7311"},
            {gene.gene_id for gene in result.genes},
        )
        ade6 = next(
            context
            for context in result.transcripts
            if context.gene.gene_id == "SPCC1322.13"
        )
        self.assertEqual("CDS", ade6.block.feature_type)
        self.assertEqual((455, 1659), (
            ade6.cds_position.base,
            ade6.cds_position.total,
        ))
        self.assertEqual("unknown", ade6.gene.viability)

    def test_real_ade6_web_annotations_align_with_every_guide(self) -> None:
        annotations = self.service.annotate_guides(self.design_result.guides)
        candidates = annotation_rows(
            self.design_result.guides,
            annotations,
            target_name=self.design_result.name,
        )

        self.assertEqual(144, len(candidates))
        self.assertEqual(GUIDE, candidates[0]["seed"])
        self.assertEqual(
            self.design_result.guides[-1].seed,
            candidates[-1]["seed"],
        )
        first_genes = {
            gene["gene_id"]: gene for gene in candidates[0]["genes"]
        }
        self.assertEqual("Primary target", first_genes["SPCC1322.13"]["role"])
        self.assertEqual("unknown", first_genes["SPCC1322.13"]["viability"])
        self.assertEqual(
            "Additional overlap",
            first_genes["SPNCRNA.7311"]["role"],
        )

    def test_real_nrg1_name_fallback_preserves_legacy_design(self) -> None:
        by_name = self.service.design_gene("nrg1")
        by_id = self.service.design_gene("SPBPB2B2.01")

        self.assertEqual("nrg1", by_name.name)
        self.assertEqual(
            ("II", "4458010", "4459767"),
            (by_name.chromosome, by_name.start, by_name.end),
        )
        self.assertEqual(136, len(by_name.guides))
        self.assertEqual(by_id.guide_table, by_name.guide_table)
        self.assertEqual(by_id.hr_dna, by_name.hr_dna)
        self.assertEqual(by_id.checking_primers, by_name.checking_primers)

        html = render_design(
            by_name,
            self.service.annotate_guides(by_name.guides),
            (PROJECT_ROOT / "template" / "container_table.html").read_text(
                encoding="utf-8"
            ),
        )
        self.assertIn("<b>Name</b>=nrg1", html)
        self.assertIn('"gene_id":"SPBPB2B2.01"', html)
        self.assertIn('"role":"Primary target"', html)

    def test_gff_only_gene_design_uses_current_interval(self) -> None:
        result = self.service.design_gene("SPNCRNA.7311")
        annotations = self.service.annotate_guides(result.guides)
        choices = self.service.cassette_choices(
            result.guides,
            annotations,
            result.name,
        )
        donors = self.service.disruption_donors(
            result.guides,
            annotations,
            choices,
            80,
            result.name,
        )
        rows = annotation_rows(result.guides, annotations, result.name)

        self.assertEqual("SPNCRNA.7311", result.name)
        self.assertEqual(("III", "1316304", "1317821"), (
            result.chromosome,
            result.start,
            result.end,
        ))
        self.assertEqual(128, len(result.guides))
        self.assertTrue(all(row["coding_target"] is False for row in rows))
        self.assertTrue(all(not group for group in choices))
        self.assertTrue(all(not group for group in donors))

    def test_real_cli_result_bodies_are_byte_for_byte_compatible(self) -> None:
        self.assertEqual(
            "ff3d2cce276ae88fcd7596ba1986208bce0ca254e0ac7c4909be487cc6d9d4c9",
            hashlib.sha256(
                format_design(self.design_result).encode("utf-8")
            ).hexdigest(),
        )

    def test_version_4_cache_round_trip_excludes_shared_index(self) -> None:
        with tempfile.TemporaryDirectory() as cache_directory:
            cold_service = Crispr4pService.from_project_data(
                precomputed_folder=cache_directory,
            )
            cold_result = cold_service.design_gene("ade6")
            cache_path = (
                Path(cache_directory) / "SPCC1322.13_v4_n0.pickle"
            )

            with cache_path.open("rb") as cache_file:
                cached_result = pickle.load(cache_file)

            warm_service = Crispr4pService.from_project_data(
                precomputed_folder=cache_directory,
            )
            warm_result = warm_service.design_gene("ade6")

        self.assertEqual(4, len(cached_result))
        self.assertIsInstance(cached_result[3], dict)
        first_match_table = next(iter(cached_result[3].values()))
        self.assertIsInstance(first_match_table[8], list)
        self.assertEqual(
            cold_result.to_legacy()[:3],
            warm_result.to_legacy()[:3],
        )
        self.assertIsNotNone(cold_service.genome_index)
        self.assertIsNone(warm_service.genome_index)
        self.assertEqual(
            "74fa7e168045e969add6b88b2999025f60099aaf96911afac9a3b26c3968c53f",
            hashlib.sha256(
                format_oligo(self.result).encode("utf-8")
            ).hexdigest(),
        )


if __name__ == "__main__":
    unittest.main()
