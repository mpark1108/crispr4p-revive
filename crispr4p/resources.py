"""FASTA and legacy annotation data."""

from collections.abc import Mapping
from dataclasses import dataclass
from types import MappingProxyType


class GeneNameNotFoundError(LookupError):
    """Raised when a gene name or ID is not found."""

    def __init__(self, query):
        self.query = query
        super().__init__(f'Gene "{query}" was not found')


def _matches_name(row, name):
    # The legacy third column can contain several comma-separated synonyms.
    return any(
        name == synonym.strip()
        for value in row
        for synonym in value.split(",")
    )


@dataclass(frozen=True, slots=True, init=False, eq=False)
class chromosomeFasta:
    """One FASTA record."""

    header: str
    sequence: str
    name: str

    def __init__(self, data):
        lines = data.split("\n")
        header = lines[0]
        object.__setattr__(self, "header", header)
        object.__setattr__(self, "sequence", "".join(lines[1:]))
        object.__setattr__(self, "name", header[:header.index(" ")])

    def __str__(self):
        return " ".join(
            [
                "chromosome:",
                self.name,
                "Length:",
                str(len(self.sequence)),
                "Header:",
                self.header,
            ]
        )


def read_fasta(sequence_file):
    """Read FASTA records in file order."""
    with open(sequence_file, "r", encoding="utf-8") as fasta_file:
        data = fasta_file.read()

    chromosomes = {}
    for raw_record in data.split(">"):
        if raw_record:
            chromosome = chromosomeFasta(raw_record)
            chromosomes[chromosome.name] = chromosome
    return MappingProxyType(chromosomes)


@dataclass(frozen=True, slots=True, init=False, eq=False)
class AnnotationParser:
    """Gene coordinates and synonyms."""

    coordinates_: tuple
    synonims_: tuple

    def __init__(self, coordinates_txt, synonims_txt):
        coordinates = self.readCoordinates_(coordinates_txt)
        synonyms = self.readSynonims_(synonims_txt)
        object.__setattr__(
            self,
            "coordinates_",
            tuple(tuple(row) for row in coordinates),
        )
        object.__setattr__(
            self,
            "synonims_",
            tuple(tuple(row) for row in synonyms),
        )

    def readCoordinates_(self, coordinates_txt):
        with open(coordinates_txt, "r", encoding="utf-8") as coordinates_file:
            data = [line.rstrip() for line in coordinates_file.readlines()][1:]
        return [line.split("\t") for line in data]

    def readSynonims_(self, synonims_txt):
        with open(synonims_txt, "r", encoding="utf-8") as synonyms_file:
            data = [line.rstrip() for line in synonyms_file.readlines()][2:]
        rows = [line.split("\t") for line in data]
        return [[value for value in row if value] for row in rows]

    def normalize_name(self, name):
        name = name.strip()
        normalized = name.casefold()
        for row in self.coordinates_:
            if row[0].casefold() == normalized:
                return row[0]
        return name.lower()

    def getCoordsFromName(self, name):
        input_name = name
        name = self.normalize_name(name)

        try:
            found = next(
                row for row in self.synonims_ if _matches_name(row, name)
            )[0]
        except StopIteration:
            raise GeneNameNotFoundError(input_name)

        coordinates = next(
            row for row in self.coordinates_ if row[0] == found
        )[1:]
        # Callers expect a fresh list.
        return list(coordinates)


@dataclass(frozen=True, slots=True, eq=False)
class ReferenceResources:
    """Shared FASTA and gene-coordinate data."""

    chromosomes: Mapping
    annotations: AnnotationParser

    def __post_init__(self):
        # Copy the mapping so shared reference data cannot change underneath us.
        object.__setattr__(
            self,
            "chromosomes",
            MappingProxyType(dict(self.chromosomes)),
        )

    @classmethod
    def from_files(cls, sequence_file, coordinates_file, synonyms_file):
        return cls(
            chromosomes=read_fasta(sequence_file),
            annotations=AnnotationParser(coordinates_file, synonyms_file),
        )
