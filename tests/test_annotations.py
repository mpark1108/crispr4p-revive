import hashlib
import tempfile
import unittest
from dataclasses import FrozenInstanceError
from pathlib import Path

from crispr4p.annotations import (
    AmbiguousGeneNameError,
    GeneNames,
    GenomeAnnotations,
    read_viability,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]
PROJECT_GFF = (
    PROJECT_ROOT / "data" / "Schizosaccharomyces_pombe_all_chromosomes.gff3"
)
PROJECT_VIABILITY = PROJECT_ROOT / "data" / "gene_viability.tsv"
PROJECT_NAMES = PROJECT_ROOT / "data" / "gene_IDs_names.tsv"


GFF = """\
##gff-version 3
I\ttest\tgene\t100\t500\t.\t+\t.\tID=plus;Name=plus_gene;so_term_name=protein_coding_gene
I\ttest\tmRNA\t100\t500\t.\t+\t.\tID=plus.1;Parent=plus
I\ttest\tfive_prime_UTR\t100\t119\t.\t+\t.\tID=u1;Parent=plus.1
I\ttest\tCDS\t120\t219\t.\t+\t0\tID=c1;Parent=plus.1
I\ttest\tintron\t220\t299\t.\t+\t.\tID=i1;Parent=plus.1
I\ttest\tCDS\t300\t399\t.\t+\t0\tID=c2;Parent=plus.1
I\ttest\tthree_prime_UTR\t400\t500\t.\t+\t.\tID=u2;Parent=plus.1
I\ttest\tmRNA\t300\t500\t.\t+\t.\tID=plus.alt;Parent=plus
I\ttest\tCDS\t300\t450\t.\t+\t0\tID=ac1;Parent=plus.alt
I\ttest\tthree_prime_UTR\t451\t500\t.\t+\t.\tID=au1;Parent=plus.alt
I\ttest\tgene\t600\t899\t.\t-\t.\tID=minus;Name=minus_gene;so_term_name=protein_coding_gene
I\ttest\tmRNA\t600\t899\t.\t-\t.\tID=minus.1;Parent=minus
I\ttest\tCDS\t800\t899\t.\t-\t0\tID=mc1;Parent=minus.1
I\ttest\tintron\t700\t799\t.\t-\t.\tID=mi1;Parent=minus.1
I\ttest\tCDS\t600\t699\t.\t-\t0\tID=mc2;Parent=minus.1
I\ttest\tgene\t180\t260\t.\t-\t.\tID=antisense;Name=antisense_gene;so_term_name=lncRNA_gene
I\ttest\tlncRNA\t180\t260\t.\t-\t.\tID=antisense.1;Parent=antisense
I\ttest\texon\t180\t260\t.\t-\t.\tID=ae1;Parent=antisense.1
"""

VIABILITY = """\
plus\tviable
minus\tinviable
antisense\tcondition-dependent
"""

NAMES = """\
# Chado database date: test
gene_systematic_id\tgene_name\tsynonyms
plus\tplus_gene\tplus_alias,shared
minus\tminus_gene\tminus_alias,shared,plus,plus_gene
antisense\t\told_antisense
"""


class GeneNamesTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary_directory = tempfile.TemporaryDirectory()
        cls.names_path = Path(cls.temporary_directory.name) / "gene_names.tsv"
        cls.names_path.write_text(NAMES, encoding="utf-8")
        cls.names = GeneNames.from_file(cls.names_path)

    @classmethod
    def tearDownClass(cls):
        cls.temporary_directory.cleanup()

    def test_resolves_ids_names_and_comma_separated_synonyms(self):
        by_id = self.names.find(" PLUS ")
        by_name = self.names.find("PLUS_GENE")
        by_alias = self.names.find("plus_alias")

        self.assertIs(by_id, by_name)
        self.assertIs(by_name, by_alias)
        self.assertEqual("plus", by_alias.gene_id)
        self.assertEqual("plus_gene", by_alias.name)
        self.assertEqual(("plus_alias", "shared"), by_alias.synonyms)
        self.assertFalse(by_alias.is_alias("plus_gene"))
        self.assertTrue(by_alias.is_alias("plus_alias"))
        self.assertIsNone(self.names.find("missing"))

    def test_reports_ambiguous_names_with_systematic_ids(self):
        with self.assertRaises(AmbiguousGeneNameError) as raised:
            self.names.find("shared")

        self.assertEqual("shared", raised.exception.query)
        self.assertEqual(("plus", "minus"), raised.exception.gene_ids)
        self.assertIn("plus, minus", str(raised.exception))

    def test_rejects_malformed_and_duplicate_rows(self):
        cases = (
            ("wrong\theader\n", "header"),
            (
                "gene_systematic_id\tgene_name\tsynonyms\n"
                "plus\tplus_gene\n",
                "three columns",
            ),
            (
                "gene_systematic_id\tgene_name\tsynonyms\n"
                "plus\tplus_gene\t\nplus\tother\t\n",
                "duplicate",
            ),
        )
        for text, message in cases:
            with self.subTest(message=message):
                with tempfile.NamedTemporaryFile(
                    "w", encoding="utf-8"
                ) as handle:
                    handle.write(text)
                    handle.flush()
                    with self.assertRaisesRegex(ValueError, message):
                        GeneNames.from_file(handle.name)


class GenomeAnnotationsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary_directory = tempfile.TemporaryDirectory()
        directory = Path(cls.temporary_directory.name)
        cls.gff_path = directory / "annotations.gff3"
        cls.viability_path = directory / "gene_viability.tsv"
        cls.gff_path.write_text(GFF, encoding="utf-8")
        cls.viability_path.write_text(VIABILITY, encoding="utf-8")
        cls.annotations = GenomeAnnotations.from_files(
            cls.gff_path,
            cls.viability_path,
        )

    @classmethod
    def tearDownClass(cls):
        cls.temporary_directory.cleanup()

    def context(self, cut_left, transcript_id):
        result = self.annotations.annotate_cut("I", (cut_left, cut_left + 1))
        return next(
            context
            for context in result.transcripts
            if context.transcript_id == transcript_id
        )

    def test_classifies_cds_intron_utr_and_non_coding_exon(self):
        self.assertEqual(
            "five_prime_UTR",
            self.context(110, "plus.1").block.feature_type,
        )
        self.assertEqual("CDS", self.context(150, "plus.1").block.feature_type)
        self.assertEqual("intron", self.context(250, "plus.1").block.feature_type)
        self.assertEqual(
            "three_prime_UTR",
            self.context(450, "plus.1").block.feature_type,
        )
        self.assertEqual("exon", self.context(200, "antisense.1").block.feature_type)

    def test_plus_strand_spliced_cds_position(self):
        position = self.context(300, "plus.1").cds_position

        self.assertEqual((101, 200), (position.base, position.total))
        self.assertEqual(50.5, position.percent)

    def test_minus_strand_spliced_cds_position(self):
        position = self.context(699, "minus.1").cds_position

        self.assertEqual((101, 200), (position.base, position.total))
        self.assertEqual(50.5, position.percent)

    def test_block_length_cut_position_and_neighbors(self):
        context = self.context(150, "plus.1")

        self.assertEqual("within", context.relation)
        self.assertEqual((120, 219, 100), (
            context.block.start,
            context.block.end,
            context.block.length,
        ))
        self.assertEqual((31, 69), (context.lower_bases, context.higher_bases))
        self.assertEqual("five_prime_UTR", context.upstream.block.feature_type)
        self.assertEqual(31, context.upstream.distance)
        self.assertEqual("intron", context.downstream.block.feature_type)
        self.assertEqual(69, context.downstream.distance)

    def test_utr_without_an_upstream_region(self):
        context = self.context(110, "plus.1")

        self.assertEqual(20, context.block.length)
        self.assertEqual((11, 9), (context.lower_bases, context.higher_bases))
        self.assertIsNone(context.upstream)
        self.assertEqual("CDS", context.downstream.block.feature_type)
        self.assertEqual(9, context.downstream.distance)

    def test_parent_gene_overlap_does_not_imply_transcript_overlap(self):
        result = self.annotations.annotate_cut("I", (150, 151))

        self.assertEqual(
            {"plus.1"},
            {context.transcript_id for context in result.transcripts},
        )

    def test_cut_between_cds_and_intron_reports_both_regions(self):
        context = self.context(219, "plus.1")

        self.assertEqual("boundary", context.relation)
        self.assertIsNone(context.block)
        self.assertEqual("CDS", context.left.block.feature_type)
        self.assertEqual("intron", context.right.block.feature_type)
        self.assertEqual(0, context.upstream.distance)
        self.assertEqual(0, context.downstream.distance)

    def test_minus_strand_upstream_and_downstream_are_reversed(self):
        context = self.context(750, "minus.1")

        self.assertEqual("higher", context.upstream_direction)
        self.assertEqual("lower", context.downstream_direction)
        self.assertEqual((51, 49), (context.lower_bases, context.higher_bases))
        self.assertEqual((800, 899), (
            context.upstream.block.start,
            context.upstream.block.end,
        ))
        self.assertEqual(49, context.upstream.distance)
        self.assertEqual((600, 699), (
            context.downstream.block.start,
            context.downstream.block.end,
        ))
        self.assertEqual(51, context.downstream.distance)

    def test_overlapping_genes_retain_type_and_viability(self):
        result = self.annotations.annotate_cut("I", (200, 201))

        self.assertEqual({"plus", "antisense"}, {
            gene.gene_id for gene in result.genes
        })
        plus = next(gene for gene in result.genes if gene.gene_id == "plus")
        antisense = next(
            gene for gene in result.genes if gene.gene_id == "antisense"
        )
        self.assertTrue(plus.is_protein_coding)
        self.assertEqual("viable", plus.viability)
        self.assertFalse(antisense.is_protein_coding)
        self.assertEqual("condition-dependent", antisense.viability)

    def test_finds_current_gene_name_or_systematic_id_case_insensitively(self):
        by_name = self.annotations.find_gene("  PLUS_GENE ")
        by_id = self.annotations.find_gene("PLUS")

        self.assertEqual("plus", by_name.gene_id)
        self.assertEqual("plus_gene", by_name.name)
        self.assertEqual(by_name, by_id)
        self.assertIsNone(self.annotations.find_gene("missing"))

    def test_intergenic_cut_reports_nearest_genes_and_viability(self):
        result = self.annotations.annotate_cut("I", (550, 551))

        self.assertTrue(result.is_intergenic)
        self.assertEqual("plus", result.lower_gene.gene.gene_id)
        self.assertEqual(50, result.lower_gene.distance)
        self.assertEqual("viable", result.lower_gene.gene.viability)
        self.assertEqual("minus", result.higher_gene.gene.gene_id)
        self.assertEqual(49, result.higher_gene.distance)
        self.assertEqual("inviable", result.higher_gene.gene.viability)

    def test_missing_viability_is_distinct_from_explicit_unknown(self):
        with tempfile.NamedTemporaryFile("w", encoding="utf-8") as viability:
            viability.write("plus\tunknown\n")
            viability.flush()
            annotations = GenomeAnnotations.from_files(
                self.gff_path,
                viability.name,
            )

        result = annotations.annotate_cut("I", (650, 651))
        self.assertIsNone(result.genes[0].viability)
        plus = annotations.annotate_cut("I", (150, 151)).genes[0]
        self.assertEqual("unknown", plus.viability)

    def test_results_and_viability_mapping_are_read_only(self):
        result = self.annotations.annotate_cut("I", (150, 151))

        with self.assertRaises(FrozenInstanceError):
            result.chromosome = "II"
        with self.assertRaises(TypeError):
            self.annotations.viability["plus"] = "inviable"
        with self.assertRaises(FrozenInstanceError):
            self.annotations._viability = {}

    def test_rejects_non_adjacent_or_non_positive_cut_coordinates(self):
        with self.assertRaisesRegex(ValueError, "adjacent"):
            self.annotations.annotate_cut("I", (150, 152))
        with self.assertRaisesRegex(ValueError, "positive"):
            self.annotations.annotate_cut("I", (0, 1))


class GeneViabilityTests(unittest.TestCase):
    def test_accepts_all_pombase_statuses_and_comments(self):
        text = (
            "# downloaded from PomBase\n"
            "a\tviable\n"
            "b\tinviable\n"
            "c\tcondition-dependent\n"
            "d\tunknown\n"
        )
        with tempfile.NamedTemporaryFile("w", encoding="utf-8") as handle:
            handle.write(text)
            handle.flush()
            result = read_viability(handle.name)

        self.assertEqual(
            {
                "a": "viable",
                "b": "inviable",
                "c": "condition-dependent",
                "d": "unknown",
            },
            dict(result),
        )

    def test_rejects_malformed_unknown_and_duplicate_rows(self):
        cases = (
            ("a\tviable\textra\n", "two columns"),
            ("a\tmaybe\n", "unsupported"),
            ("a\tviable\na\tinviable\n", "duplicate"),
        )
        for text, message in cases:
            with self.subTest(text=text):
                with tempfile.NamedTemporaryFile("w", encoding="utf-8") as handle:
                    handle.write(text)
                    handle.flush()
                    with self.assertRaisesRegex(ValueError, message):
                        read_viability(handle.name)


class PackagedPomBaseDataTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.annotations = GenomeAnnotations.from_files(
            PROJECT_GFF,
            PROJECT_VIABILITY,
        )
        cls.names = GeneNames.from_file(PROJECT_NAMES)

    def test_reviewed_dataset_snapshots_are_packaged_unchanged(self):
        expected = {
            PROJECT_GFF: (
                "88e4a26c16762c6d97e7f6a1600cd7dac193bdaeff4753cc1dd2123a79f6a025"
            ),
            PROJECT_VIABILITY: (
                "e9399024327be0a2a6618c8fda1dfaef6ef72c493eb1050c7a6d86d41a3d3d09"
            ),
            PROJECT_NAMES: (
                "4c688312ebb5ab80356cf60c0b85ad9f409d049e5688c0915e477de60cde4899"
            ),
        }
        for path, expected_hash in expected.items():
            with self.subTest(path=path.name):
                self.assertTrue(path.is_file())
                self.assertEqual(
                    expected_hash,
                    hashlib.sha256(path.read_bytes()).hexdigest(),
                )

        self.assertEqual(12685, len(self.annotations.viability))
        self.assertEqual(12689, len(self.names))

    def test_current_synonym_resolves_orb2_to_shk1(self):
        gene = self.names.find(" ORB2 ")

        self.assertEqual("SPBC1604.14c", gene.gene_id)
        self.assertEqual("shk1", gene.name)
        self.assertEqual(("pak1", "orb2"), gene.synonyms)

    def test_current_names_precede_ambiguous_historical_aliases(self):
        self.assertEqual("SPAC9E9.12c", self.names.find("abc1").gene_id)

    def test_current_ambiguous_alias_lists_both_gene_ids(self):
        with self.assertRaises(AmbiguousGeneNameError) as raised:
            self.names.find("noc2")

        self.assertEqual(
            ("SPAC1142.04", "SPAC1B3.09c"),
            raised.exception.gene_ids,
        )
        self.assertEqual(
            "SPNCRNA.2001",
            self.names.find("SPNCRNA.2001").gene_id,
        )

    def test_real_ade6_cut_matches_reviewed_prototype(self):
        result = self.annotations.annotate_cut(
            "III",
            (1316791, 1316792),
        )

        self.assertEqual(
            {"SPCC1322.13", "SPNCRNA.7311"},
            {gene.gene_id for gene in result.genes},
        )
        ade6 = next(
            context
            for context in result.transcripts
            if context.gene.gene_id == "SPCC1322.13"
        )
        self.assertEqual("ade6", ade6.gene.name)
        self.assertEqual("unknown", ade6.gene.viability)
        self.assertEqual("CDS", ade6.block.feature_type)
        self.assertEqual((1316337, 1317995, 1659), (
            ade6.block.start,
            ade6.block.end,
            ade6.block.length,
        ))
        self.assertEqual((455, 1659), (
            ade6.cds_position.base,
            ade6.cds_position.total,
        ))
        self.assertAlmostEqual(27.4, ade6.cds_position.percent, places=1)

        overlap = next(
            context
            for context in result.transcripts
            if context.gene.gene_id == "SPNCRNA.7311"
        )
        self.assertEqual("lncRNA_gene", overlap.gene.gene_type)
        self.assertEqual("unknown", overlap.gene.viability)
        self.assertEqual("exon", overlap.block.feature_type)
        self.assertEqual(1518, overlap.block.length)

    def test_current_gff_resolves_nrg1_missing_from_legacy_synonyms(self):
        gene = self.annotations.find_gene("NRG1")

        self.assertEqual("SPBPB2B2.01", gene.gene_id)
        self.assertEqual("nrg1", gene.name)
        self.assertEqual("II", gene.chromosome)
        self.assertEqual((4457695, 4459842), (gene.start, gene.end))


if __name__ == "__main__":
    unittest.main()
