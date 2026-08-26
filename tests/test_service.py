import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from crispr4p.guides import match_guide
from crispr4p.models import DesignResult, OligoAnalysisResult
from crispr4p.resources import GeneNameNotFoundError
from crispr4p.service import Crispr4pService, OligoLengthError


LEGACY_RESULT = (
    [["guide"]],
    ("forward", "reverse", "deleted"),
    [{"primer": "checking"}],
    "ade6",
    "III",
    "1316337",
    "1317995",
)
GUIDE = "ACATTGGCTTACGACGGTCG"


class RecordingHit:
    def __init__(self, chromosome, strand, seed, pam):
        self.chromosome = chromosome
        self.strand = strand
        self.seed = seed
        self.pam = pam


class RecordingDesigner:
    instances = []

    def __init__(self, *args, **kwargs):
        self.constructor_args = args
        self.constructor_kwargs = kwargs
        self.run_web_calls = []
        self.index_calls = 0
        self.genome_index = kwargs.get("genome_index")
        self.reference_resources = kwargs.get("reference_resources")
        if self.reference_resources is None:
            self.reference_resources = object()
        self.__class__.instances.append(self)

    def runWeb(self, **kwargs):
        self.run_web_calls.append(kwargs)
        return LEGACY_RESULT

    def getNGGsFromGenome(self):
        self.index_calls += 1
        if self.genome_index is None:
            self.genome_index = object()
        hit = RecordingHit(
            chromosome="III",
            strand=1,
            seed=GUIDE,
            pam="TGG",
        )
        self.NGGs = {GUIDE[-8:]: [hit] * 5}
        return self.genome_index

    def getOligoHitCoordinates(self, match):
        return (1316795, 1316797), (1316791, 1316792)


class CurrentNameFallbackDesigner(RecordingDesigner):
    def runWeb(self, **kwargs):
        self.run_web_calls.append(kwargs)
        if kwargs["name"] == "nrg1":
            raise GeneNameNotFoundError("nrg1")
        return LEGACY_RESULT[:3] + (kwargs["name"],) + LEGACY_RESULT[4:]


class GffOnlyDesigner(RecordingDesigner):
    def runWeb(self, **kwargs):
        self.run_web_calls.append(kwargs)
        if kwargs["name"] is not None:
            raise GeneNameNotFoundError(kwargs["name"])
        return LEGACY_RESULT[:3] + (
            None,
            kwargs["cr"],
            kwargs["start"],
            kwargs["end"],
        )


class RecordingGenomeAnnotations:
    def __init__(self):
        self.calls = []

    def annotate_cut(self, chromosome, cut_coordinates):
        call = (chromosome, cut_coordinates)
        self.calls.append(call)
        return call


class TestCrispr4pService(unittest.TestCase):
    def setUp(self) -> None:
        RecordingDesigner.instances = []
        self.service = Crispr4pService(
            "genome.fa",
            "coordinates.txt",
            "synonyms.txt",
            precomputed_folder="cache",
            designer_factory=RecordingDesigner,
        )

    def test_gene_design_delegates_without_changing_legacy_result(self) -> None:
        result = self.service.design_gene("ade6", n_mismatch=2)

        self.assertIsInstance(result, DesignResult)
        self.assertIs(LEGACY_RESULT, result.to_legacy())
        self.assertIs(LEGACY_RESULT[0], result.guide_table)
        self.assertEqual("ade6", result.name)
        designer = RecordingDesigner.instances[0]
        self.assertEqual(
            ("genome.fa", "coordinates.txt", "synonyms.txt"),
            designer.constructor_args,
        )
        self.assertEqual(
            {
                "precomputed_folder": "cache",
                "verbose": False,
            },
            designer.constructor_kwargs,
        )
        self.assertEqual(
            [
                {
                    "name": "ade6",
                    "cr": None,
                    "start": None,
                    "end": None,
                    "strand": None,
                    "nMismatch": 2,
                }
            ],
            designer.run_web_calls,
        )

    def test_missing_legacy_name_falls_back_to_current_gff_name(self) -> None:
        annotation_index = SimpleNamespace(
            find_gene=lambda name: SimpleNamespace(
                gene_id="SPBPB2B2.01",
                name="nrg1",
            )
        )
        service = Crispr4pService(
            "genome.fa",
            "coordinates.txt",
            "synonyms.txt",
            designer_factory=CurrentNameFallbackDesigner,
            genome_annotations=annotation_index,
        )

        result = service.design_gene("nrg1", n_mismatch=2)

        self.assertEqual("nrg1", result.name)
        self.assertEqual(2, len(CurrentNameFallbackDesigner.instances))
        self.assertEqual(
            ["nrg1", "SPBPB2B2.01"],
            [
                designer.run_web_calls[0]["name"]
                for designer in CurrentNameFallbackDesigner.instances
            ],
        )
        self.assertIs(LEGACY_RESULT[0], result.guide_table)

    def test_current_synonym_uses_systematic_id_and_primary_name(self) -> None:
        gene = SimpleNamespace(
            gene_id="SPBC1604.14c",
            name="shk1",
            is_alias=lambda query: query.strip().casefold() == "orb2",
        )
        names = SimpleNamespace(find=lambda query: gene)
        service = Crispr4pService(
            "genome.fa",
            "coordinates.txt",
            "synonyms.txt",
            designer_factory=RecordingDesigner,
            gene_names=names,
        )

        result = service.design_gene("orb2", n_mismatch=2)

        self.assertEqual("shk1", result.name)
        self.assertEqual(1, len(RecordingDesigner.instances))
        self.assertEqual(
            "SPBC1604.14c",
            RecordingDesigner.instances[0].run_web_calls[0]["name"],
        )
        self.assertEqual(
            2,
            RecordingDesigner.instances[0].run_web_calls[0]["nMismatch"],
        )

    def test_current_name_table_loads_once(self) -> None:
        calls = []
        names = SimpleNamespace(find=lambda query: None)

        def load_names(path):
            calls.append(path)
            return names

        service = Crispr4pService(
            "genome.fa",
            "coordinates.txt",
            "synonyms.txt",
            designer_factory=RecordingDesigner,
            gene_names_file="gene_names.tsv",
            gene_names_factory=load_names,
        )

        service.design_gene("ade6")
        service.design_gene("ura4")

        self.assertEqual(["gene_names.tsv"], calls)
        self.assertIs(names, service.gene_names)

    def test_gff_only_gene_uses_its_interval(self) -> None:
        annotation_index = SimpleNamespace(
            find_gene=lambda name: SimpleNamespace(
                gene_id="SPNCRNA.7311",
                name=None,
                chromosome="III",
                start=1316304,
                end=1317821,
            )
        )
        service = Crispr4pService(
            "genome.fa",
            "coordinates.txt",
            "synonyms.txt",
            designer_factory=GffOnlyDesigner,
            genome_annotations=annotation_index,
        )

        result = service.design_gene("SPNCRNA.7311", n_mismatch=2)

        self.assertEqual("SPNCRNA.7311", result.name)
        self.assertEqual(("III", "1316304", "1317821"), (
            result.chromosome,
            result.start,
            result.end,
        ))
        self.assertEqual(3, len(GffOnlyDesigner.instances))
        self.assertEqual(
            [
                ("SPNCRNA.7311", None, None, None),
                ("SPNCRNA.7311", None, None, None),
                (None, "III", "1316304", "1317821"),
            ],
            [
                (
                    designer.run_web_calls[0]["name"],
                    designer.run_web_calls[0]["cr"],
                    designer.run_web_calls[0]["start"],
                    designer.run_web_calls[0]["end"],
                )
                for designer in GffOnlyDesigner.instances
            ],
        )
        self.assertTrue(all(
            designer.run_web_calls[0]["nMismatch"] == 2
            for designer in GffOnlyDesigner.instances
        ))

    def test_name_absent_from_legacy_and_gff_keeps_lookup_error(self) -> None:
        annotation_index = SimpleNamespace(find_gene=lambda name: None)
        service = Crispr4pService(
            "genome.fa",
            "coordinates.txt",
            "synonyms.txt",
            designer_factory=CurrentNameFallbackDesigner,
            genome_annotations=annotation_index,
        )

        with self.assertRaises(GeneNameNotFoundError) as raised:
            service.design_gene("nrg1")

        self.assertEqual("nrg1", raised.exception.query)

    def test_region_design_delegates_all_coordinates_and_options(self) -> None:
        result = self.service.design_region(
            "III",
            "100",
            "200",
            strand="-1",
            n_mismatch=1,
        )

        self.assertIsInstance(result, DesignResult)
        self.assertIs(LEGACY_RESULT, result.to_legacy())
        self.assertEqual("III", result.chromosome)
        self.assertEqual("1316337", result.start)
        self.assertEqual("1317995", result.end)
        self.assertEqual(
            [
                {
                    "name": None,
                    "cr": "III",
                    "start": "100",
                    "end": "200",
                    "strand": "-1",
                    "nMismatch": 1,
                }
            ],
            RecordingDesigner.instances[0].run_web_calls,
        )

    def test_each_request_gets_fresh_query_state(self) -> None:
        self.service.design_gene("ade6")
        self.service.design_region("III", "100", "200")
        self.service.analyze_oligo(GUIDE)
        self.service.analyze_oligo(GUIDE)

        self.assertEqual(4, len(RecordingDesigner.instances))
        self.assertIsNot(
            RecordingDesigner.instances[0],
            RecordingDesigner.instances[1],
        )
        self.assertIsNot(
            RecordingDesigner.instances[1],
            RecordingDesigner.instances[2],
        )
        self.assertIsNot(
            RecordingDesigner.instances[2],
            RecordingDesigner.instances[3],
        )
        shared_index = RecordingDesigner.instances[2].genome_index
        self.assertIs(shared_index, self.service.genome_index)
        self.assertIs(
            shared_index,
            RecordingDesigner.instances[3].constructor_kwargs["genome_index"],
        )
        shared_resources = RecordingDesigner.instances[0].reference_resources
        self.assertIs(shared_resources, self.service.reference_resources)
        for designer in RecordingDesigner.instances[1:]:
            self.assertIs(
                shared_resources,
                designer.constructor_kwargs["reference_resources"],
            )

    def test_injected_shared_components_are_forwarded(self) -> None:
        genome_index = object()
        reference_resources = object()
        service = Crispr4pService(
            "genome.fa",
            "coordinates.txt",
            "synonyms.txt",
            designer_factory=RecordingDesigner,
            genome_index=genome_index,
            reference_resources=reference_resources,
        )

        service.design_gene("ade6")

        self.assertIs(
            genome_index,
            RecordingDesigner.instances[0].constructor_kwargs["genome_index"],
        )
        self.assertIs(
            reference_resources,
            RecordingDesigner.instances[0].constructor_kwargs[
                "reference_resources"
            ],
        )

    def test_cut_annotations_load_lazily_and_are_reused(self) -> None:
        annotation_index = RecordingGenomeAnnotations()
        factory_calls = []

        def annotation_factory(gff_file, viability_file):
            factory_calls.append((gff_file, viability_file))
            return annotation_index

        service = Crispr4pService(
            "genome.fa",
            "coordinates.txt",
            "synonyms.txt",
            designer_factory=RecordingDesigner,
            annotation_gff_file="annotations.gff3",
            gene_viability_file="viability.tsv",
            annotation_factory=annotation_factory,
        )

        self.assertIsNone(service.genome_annotations)
        first = service.annotate_cut("III", (100, 101))
        second = service.annotate_cut("II", (200, 201))

        self.assertEqual(
            [("annotations.gff3", "viability.tsv")],
            factory_calls,
        )
        self.assertIs(annotation_index, service.genome_annotations)
        self.assertEqual(("III", (100, 101)), first)
        self.assertEqual(("II", (200, 201)), second)
        self.assertEqual(
            [("III", (100, 101)), ("II", (200, 201))],
            annotation_index.calls,
        )

    def test_guide_annotations_use_named_fields_and_preserve_order(self) -> None:
        annotation_index = RecordingGenomeAnnotations()
        service = Crispr4pService(
            "genome.fa",
            "coordinates.txt",
            "synonyms.txt",
            designer_factory=RecordingDesigner,
            genome_annotations=annotation_index,
        )
        guides = (
            SimpleNamespace(chromosome="III", cut_coordinates=(100, 101)),
            SimpleNamespace(chromosome="I", cut_coordinates=(300, 301)),
        )

        results = service.annotate_guides(guides)

        self.assertIsInstance(results, tuple)
        self.assertEqual(
            (("III", (100, 101)), ("I", (300, 301))),
            results,
        )
        self.assertEqual(
            [("III", (100, 101)), ("I", (300, 301))],
            annotation_index.calls,
        )

    def test_unconfigured_direct_service_reports_annotation_error(self) -> None:
        with self.assertRaisesRegex(RuntimeError, "GFF3 is not configured"):
            self.service.annotate_cut("III", (100, 101))

    def test_oligo_analysis_returns_structured_existing_results(self) -> None:
        with patch(
            "crispr4p.service.match_guide",
            wraps=match_guide,
        ) as matcher:
            result = self.service.analyze_oligo(
                f"  {GUIDE.lower()}  ",
                n_mismatch=1,
            )

        self.assertIsInstance(result, OligoAnalysisResult)
        self.assertEqual(GUIDE, result.oligo_sequence)
        self.assertEqual(GUIDE, result.seed)
        self.assertEqual("NGG", result.pam)
        self.assertEqual(1, result.n_mismatch)
        self.assertEqual(5, result.match_counts[8])
        self.assertEqual(1, result.match_counts[20])
        self.assertEqual(1, len(result.full_matches))
        self.assertEqual("III", result.full_matches[0].chromosome)
        self.assertEqual(
            (1316795, 1316797),
            result.full_matches[0].pam_coordinates,
        )
        self.assertEqual(
            (1316791, 1316792),
            result.full_matches[0].cut_coordinates,
        )
        designer = RecordingDesigner.instances[0]
        self.assertEqual(1, designer.index_calls)
        query, genome_hits, n_mismatch = matcher.call_args.args
        self.assertEqual(GUIDE, query.seed)
        self.assertEqual("NGG", query.pam)
        self.assertIs(designer.NGGs, genome_hits)
        self.assertEqual(1, n_mismatch)

    def test_oligo_analysis_preserves_supplied_pam(self) -> None:
        result = self.service.analyze_oligo(GUIDE + "TGG")

        self.assertEqual(GUIDE + "TGG", result.oligo_sequence)
        self.assertEqual(GUIDE, result.seed)
        self.assertEqual("TGG", result.pam)

    def test_oligo_analysis_rejects_invalid_length_before_loading_data(self) -> None:
        with self.assertRaises(OligoLengthError) as raised:
            self.service.analyze_oligo("A" * 19)

        self.assertEqual(19, raised.exception.sequence_length)
        self.assertEqual([], RecordingDesigner.instances)

    def test_project_data_factory_uses_packaged_reference_files(self) -> None:
        service = Crispr4pService.from_project_data(
            precomputed_folder="cache",
            designer_factory=RecordingDesigner,
        )
        service.design_gene("ade6")

        sequence_file, coordinates_file, synonyms_file = (
            Path(value)
            for value in RecordingDesigner.instances[0].constructor_args
        )
        self.assertEqual("data", sequence_file.parent.name)
        self.assertEqual(
            "Schizosaccharomyces_pombe.ASM294v2.26.dna.toplevel.fa",
            sequence_file.name,
        )
        self.assertEqual("COORDINATES.txt", coordinates_file.name)
        self.assertEqual("SYNONIMS.txt", synonyms_file.name)
        self.assertTrue(sequence_file.is_file())
        self.assertTrue(coordinates_file.is_file())
        self.assertTrue(synonyms_file.is_file())
        self.assertEqual(
            "Schizosaccharomyces_pombe_all_chromosomes.gff3",
            Path(service.annotation_gff_file).name,
        )
        self.assertEqual(
            "gene_viability.tsv",
            Path(service.gene_viability_file).name,
        )
        self.assertEqual(
            "gene_IDs_names.tsv",
            Path(service.gene_names_file).name,
        )
        self.assertTrue(Path(service.annotation_gff_file).is_file())
        self.assertTrue(Path(service.gene_viability_file).is_file())
        self.assertTrue(Path(service.gene_names_file).is_file())
        self.assertIsNone(service.genome_annotations)
        self.assertIsNotNone(service.gene_names)


if __name__ == "__main__":
    unittest.main()
