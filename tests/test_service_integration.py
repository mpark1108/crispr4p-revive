import hashlib
import pickle
import tempfile
import unittest
from collections import Counter
from pathlib import Path

from crispr4p.crispr4p import PrimerDesign
from crispr4p.genome import GenomePamIndex
from crispr4p.service import Crispr4pService
from crispr4p.web_views import (
    annotation_rows,
    render_design,
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

        data = PROJECT_ROOT / "data"
        cls.core_designer = PrimerDesign(
            data / "Schizosaccharomyces_pombe.ASM294v2.26.dna.toplevel.fa",
            data / "COORDINATES.txt",
            data / "SYNONIMS.txt",
            regression=True,
        )
        (
            cls.core_table,
            cls.core_hr_dna,
            cls.core_checking_primers,
            cls.core_matches,
        ) = cls.core_designer.run_(
            "III",
            1316337,
            1317995,
            0,
            "ade6",
        )

    def test_real_ade6_oligo_values(self) -> None:
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

    def test_real_reverse_hit(self) -> None:
        result = self.service.analyze_oligo(
            "AGCTCTCCATAGCAGCCAAA",
            n_mismatch=0,
        )

        self.assertEqual(
            {8: 57, 10: 10, 12: 4, 14: 1, 16: 1, 18: 1, 20: 1},
            dict(result.match_counts),
        )
        self.assertEqual(1, len(result.full_matches))
        match = result.full_matches[0]
        self.assertEqual("III", match.chromosome)
        self.assertEqual((1317893, 1317895), match.pam_coordinates)
        self.assertEqual((1317898, 1317899), match.cut_coordinates)
        self.assertEqual(-1, match.strand)
        self.assertEqual("AGG", match.pam)

    def test_real_ade6_design_values(self) -> None:
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

    def test_uncached_ade6_outputs(self) -> None:
        expected_candidates = [
            ("ACATTGGCTTACGACGGTCG", [5, 1, 1, 1, 1, 1, 1]),
            ("GTGGCGACAGGGACACCTCG", [7, 1, 1, 1, 1, 1, 1]),
            ("AAACAGGTTGTAGGGATCCT", [10, 1, 1, 1, 1, 1, 1]),
            ("TGCCGCTCATTTGCCTGGTA", [12, 1, 1, 1, 1, 1, 1]),
            ("GCATAAGTAACCATGCGATC", [12, 1, 1, 1, 1, 1, 1]),
        ]
        self.assertEqual(144, len(self.core_table))
        self.assertEqual(
            expected_candidates,
            [(row[0], row[2:]) for row in self.core_table[:5]],
        )
        self.assertEqual(len(self.core_table), len(self.core_matches))
        self.assertEqual(
            (
                GUIDE,
                "ACGACGGTCGgttttagagctagaaatagcaagttaaaataa",
                "AAGCCAATGTttcttcggtacaggttatgttttttggcaaca",
                (1316795, 1316797),
                1,
                "TGG",
            ),
            self.core_table[0][1],
        )
        self.assertEqual(
            [
                "3948f29bd72bc5059eff88921869de90f9c99190232e4add15ecbd5eac7629dc",
                "62b0d6cbd581c1824d3437af1aeab7777ccc8127bc0642de72553f856ca0a78a",
                "66476b056e174bd5d37b0119377f2a829b9d0b01dc6295b7a6906faaeb4b5718",
            ],
            [
                hashlib.sha256(sequence.encode("ascii")).hexdigest()
                for sequence in self.core_hr_dna
            ],
        )

        first, second = self.core_checking_primers
        self.assertEqual(
            ("AGCCTGGTGCAGTATAAGGT", "CGTCGCAGCACATTATTCGG", 444, 2103),
            (
                first["PRIMER_LEFT_0_SEQUENCE"],
                first["PRIMER_RIGHT_0_SEQUENCE"],
                first["PRIMER_PAIR_0_PRODUCT_SIZE"],
                first["negative_result"],
            ),
        )
        self.assertEqual(
            ("ACTGCGCACTAACTCACTACA", "CGTCGCAGCACATTATTCGG", 254, 1913),
            (
                second["PRIMER_LEFT_1_SEQUENCE"],
                second["PRIMER_RIGHT_1_SEQUENCE"],
                second["PRIMER_PAIR_1_PRODUCT_SIZE"],
                second["negative_result"],
            ),
        )
        guides = {row[0]: row[1][3] for row in self.core_table}
        self.assertEqual((1316357, 1316359), guides["ATGAGCGAAAAACAGGTTGT"])
        self.assertNotIn("TTGGAAAAATTATTCTGCAT", guides)

    def test_real_index_matches_independent_scan(self) -> None:
        index = self.core_designer.genome_index
        self.assertIsInstance(index, GenomePamIndex)
        self.assertEqual(65366, len(index))
        self.assertEqual(2267530, index.hit_count)
        self.assertNotIn("", index)
        self.assertEqual(
            [
                ("III", 1316795, 1, GUIDE, "TGG"),
                ("III", 347120, 1, "GAGAGAAGTACGGACGGTCG", "GAG"),
                ("III", 1233769, 1, "TCTCAAGATTAAGACGGTCG", "TAG"),
                ("II", 397432, 1, "AGACATTCGCGGGACGGTCG", "TGG"),
                ("I", 3430748, 1, "GGCTACGCTAGGGACGGTCG", "AGG"),
            ],
            [
                (hit.chromosome, hit.pos, hit.strand, hit.seed, hit.pam)
                for hit in index[GUIDE[-8:]]
            ],
        )

        expected = Counter()
        for chromosome_name, chromosome in self.core_designer.chromosomesData.items():
            strands = {
                1: chromosome.sequence,
                -1: self.core_designer.reverseComplement(chromosome.sequence),
            }
            for strand, sequence in strands.items():
                for position in range(21, len(sequence) - 1):
                    if sequence[position] in "GA" and sequence[position + 1] == "G":
                        expected[
                            (
                                chromosome_name,
                                strand,
                                sequence[position - 1:position + 2],
                            )
                        ] += 1

        actual = Counter()
        for hits in index.values():
            for hit in hits:
                actual[(hit.chromosome, hit.strand, hit.pam)] += 1
        self.assertEqual(expected, actual)

    def test_orb6_without_deletion_primers_still_renders(self) -> None:
        with tempfile.TemporaryDirectory() as cache:
            service = Crispr4pService.from_project_data(
                precomputed_folder=cache,
            )
            result = service.design_gene("orb6")
            annotations = service.annotate_guides(result.guides)

        template = (
            PROJECT_ROOT / "template" / "container_table.html"
        ).read_text(encoding="utf-8")
        page = render_design(result, annotations, template)
        start = page.index("<h5>Deletion-checking primers</h5>")
        end = page.index("<h3", start)
        primer_section = page[start:end]

        self.assertEqual(128, len(result.guides))
        self.assertEqual([], result.checking_primers)
        self.assertIn("-", primer_section)
        self.assertIn("- &deg;C", primer_section)
        self.assertNotIn("5'---3'", primer_section)
        self.assertNotIn("0 &deg;C", primer_section)

    def test_terminal_gene_keeps_guides_without_invalid_deletion_output(self):
        with tempfile.TemporaryDirectory() as cache:
            service = Crispr4pService.from_project_data(
                precomputed_folder=cache,
                genome_index=self.service.genome_index,
                reference_resources=self.service.reference_resources,
            )
            result = service.design_gene("SPAC212.11")
            annotations = service.annotate_guides(result.guides)

        template = (
            PROJECT_ROOT / "template" / "container_table.html"
        ).read_text(encoding="utf-8")
        page = render_design(result, annotations, template)
        start = page.index("<h5>Deletion HR-template primers</h5>")
        end = page.index("<h3", start)
        deletion_section = page[start:end]

        self.assertTrue(result.guides)
        self.assertEqual(("", "", ""), result.hr_dna)
        self.assertEqual([], result.checking_primers)
        self.assertNotIn("5'---3'", deletion_section)
        self.assertIn("-", deletion_section)

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
        self.assertIn('"gene_id":"SPBPB2B2.01"', html)
        self.assertIn('"role":"Primary target"', html)

    def test_real_orb2_synonym_resolves_shk1_design(self) -> None:
        by_synonym = self.service.design_gene("orb2")
        by_name = self.service.design_gene("shk1")
        by_id = self.service.design_gene("SPBC1604.14c")
        annotations = self.service.annotate_guides(by_synonym.guides)
        choices = self.service.cassette_choices(
            by_synonym.guides,
            annotations,
            by_synonym.name,
        )
        rows = annotation_rows(
            by_synonym.guides,
            annotations,
            by_synonym.name,
        )

        self.assertEqual("shk1", by_synonym.name)
        self.assertEqual(
            ("II", "3905082", "3907058"),
            (by_synonym.chromosome, by_synonym.start, by_synonym.end),
        )
        self.assertEqual(by_id.guide_table, by_synonym.guide_table)
        self.assertEqual(by_name.guide_table, by_synonym.guide_table)
        self.assertEqual(by_name.hr_dna, by_synonym.hr_dna)
        self.assertEqual(by_name.checking_primers, by_synonym.checking_primers)
        first_genes = {
            gene["gene_id"]: gene for gene in rows[0]["genes"]
        }
        self.assertEqual(
            "Primary target",
            first_genes["SPBC1604.14c"]["role"],
        )
        self.assertTrue(rows[0]["coding_target"])
        self.assertTrue(choices[0])

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


if __name__ == "__main__":
    unittest.main()
