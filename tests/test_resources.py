import hashlib
import json
import tempfile
import unittest
from dataclasses import FrozenInstanceError
from pathlib import Path
from unittest.mock import patch

from crispr4p.crispr4p import AnnotationParser, PrimerDesign, chromosomeFasta
from crispr4p.resources import GeneNameNotFoundError, ReferenceResources


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DATA_DIRECTORY = PROJECT_ROOT / "data"
FASTA = DATA_DIRECTORY / "Schizosaccharomyces_pombe.ASM294v2.26.dna.toplevel.fa"
COORDINATES = DATA_DIRECTORY / "COORDINATES.txt"
SYNONIMS = DATA_DIRECTORY / "SYNONIMS.txt"


class TestReferenceResources(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.resources = ReferenceResources.from_files(
            FASTA,
            COORDINATES,
            SYNONIMS,
        )

    def test_real_fasta_records_are_unchanged(self):
        expected = [
            (
                "MT",
                19431,
                "24e17cf92f1cb680d2ecd1d660a7cafcbecd2e1238b3db068d4636866c8fe2b2",
            ),
            (
                "AB325691",
                20000,
                "73e22d3eb9358d3bd403e07f2204ffeae35f164ea98c5cf6ea4a49bfeac16b87",
            ),
            (
                "MTR",
                20128,
                "42e810d333a85e766b1ea4f2c1f1adbd09aac3ac470c52b69fead64b348e78d2",
            ),
            (
                "III",
                2452883,
                "82e0d0f3f9d8b71d75f8693f3a01849b85f8dc9bbaf57bdb08d35d082908aaf3",
            ),
            (
                "II",
                4539804,
                "7bd7223dad8e20998c15825e2eb336b26915fafccd1dcb03274c12366ce8fb58",
            ),
            (
                "I",
                5579133,
                "cdfcd06ba9cce85838d0e22dfe0fa33dba5bc2e4e07d1974e8b29643a34f5572",
            ),
        ]
        actual = [
            (
                name,
                len(chromosome.sequence),
                hashlib.sha256(chromosome.sequence.encode("ascii")).hexdigest(),
            )
            for name, chromosome in self.resources.chromosomes.items()
        ]

        self.assertEqual(expected, actual)

    def test_real_annotation_records_and_legacy_lookup_are_unchanged(self):
        annotations = self.resources.annotations

        self.assertEqual(7007, len(annotations.coordinates_))
        self.assertEqual(7015, len(annotations.synonims_))
        self.assertEqual(
            "25c73b8aab2a7084fc8ce7038f55863b15457964f0128682f0fbc36f91481bc9",
            hashlib.sha256(
                json.dumps(
                    annotations.coordinates_,
                    separators=(",", ":"),
                    ensure_ascii=False,
                ).encode()
            ).hexdigest(),
        )
        self.assertEqual(
            "983c64a037c4e0999c7547270164f732b75b9dd3622d35ac749e25d631e2a630",
            hashlib.sha256(
                json.dumps(
                    annotations.synonims_,
                    separators=(",", ":"),
                    ensure_ascii=False,
                ).encode()
            ).hexdigest(),
        )
        coordinates = annotations.getCoordsFromName("  ADE6  ")
        self.assertIsInstance(coordinates, list)
        self.assertEqual(["III", "1316337", "1317995", "1"], coordinates)
        self.assertEqual(
            ["II", "3905082", "3907058", "1"],
            annotations.getCoordsFromName("orb2"),
        )
        self.assertEqual(
            ["II", "3905082", "3907058", "1"],
            annotations.getCoordsFromName("spbc1604.14C"),
        )

        coordinate_rows = annotations.readCoordinates_(COORDINATES)
        synonym_rows = annotations.readSynonims_(SYNONIMS)
        self.assertIsInstance(coordinate_rows, list)
        self.assertIsInstance(coordinate_rows[0], list)
        self.assertIsInstance(synonym_rows, list)
        self.assertIsInstance(synonym_rows[0], list)

    def test_missing_legacy_name_uses_specific_lookup_error(self):
        with self.assertRaises(GeneNameNotFoundError) as raised:
            self.resources.annotations.getCoordsFromName("nrg1")

        self.assertEqual("nrg1", raised.exception.query)
        self.assertEqual('Gene "nrg1" was not found', str(raised.exception))

    def test_resource_containers_and_records_are_read_only(self):
        chromosome = self.resources.chromosomes["III"]
        annotations = self.resources.annotations

        with self.assertRaises(TypeError):
            self.resources.chromosomes["III"] = chromosome
        with self.assertRaises(FrozenInstanceError):
            chromosome.sequence = "A"
        with self.assertRaises(AttributeError):
            annotations.coordinates_.append(("new",))
        with self.assertRaises(TypeError):
            annotations.synonims_[0][0] = "new"
        with self.assertRaises(FrozenInstanceError):
            annotations.coordinates_ = ()

    def test_legacy_classes_remain_available_from_core_module(self):
        self.assertIs(AnnotationParser, type(self.resources.annotations))
        self.assertIs(
            chromosomeFasta,
            type(next(iter(self.resources.chromosomes.values()))),
        )

    def test_primer_design_accepts_shared_resources_without_copying(self):
        with patch("builtins.open", side_effect=AssertionError("unexpected read")):
            designer = PrimerDesign(
                FASTA,
                COORDINATES,
                SYNONIMS,
                reference_resources=self.resources,
            )

        self.assertIs(self.resources, designer.reference_resources)
        self.assertIs(self.resources.chromosomes, designer.chromosomesData)
        self.assertIs(self.resources.annotations, designer.annotationParser_)

        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8") as fasta:
            fasta.write(">tiny chromosome:test\nAC\nGT\n")
            fasta.flush()
            legacy_chromosomes = designer.readsequence(fasta.name)

        self.assertIsInstance(legacy_chromosomes, dict)
        self.assertEqual("ACGT", legacy_chromosomes["tiny"].sequence)


if __name__ == "__main__":
    unittest.main()
