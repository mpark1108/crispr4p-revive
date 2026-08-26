"""PomBase cut-site annotations."""

from dataclasses import dataclass
from types import MappingProxyType


REGION_PRIORITY = (
    "CDS",
    "five_prime_UTR",
    "three_prime_UTR",
    "intron",
    "exon",
)

VALID_VIABILITY_STATUSES = frozenset(
    {"viable", "inviable", "condition-dependent", "unknown"}
)

VIABILITY_LABELS = MappingProxyType(
    {
        "viable": "viable (non-essential)",
        "inviable": "inviable (essential)",
        "condition-dependent": "condition-dependent",
        "unknown": "unknown",
        None: "not available",
    }
)


class AmbiguousGeneNameError(ValueError):
    """Raised when a name maps to more than one PomBase gene."""

    def __init__(self, query, gene_ids):
        self.query = query
        self.gene_ids = tuple(gene_ids)
        matches = ", ".join(self.gene_ids)
        super().__init__(
            f'Gene name "{query}" is ambiguous. Try one of these PomBase '
            f"systematic IDs: {matches}."
        )


@dataclass(frozen=True, slots=True)
class GeneName:
    """PomBase gene name and synonyms."""

    gene_id: str
    name: str | None
    synonyms: tuple[str, ...]

    def is_alias(self, query):
        normalized = str(query).strip().casefold()
        primary_names = {self.gene_id.casefold()}
        if self.name:
            primary_names.add(self.name.casefold())
        return normalized not in primary_names


def read_gene_names(path):
    """Read the PomBase gene names and identifiers table."""
    records = []
    seen_ids = set()
    header_found = False
    expected_header = ("gene_systematic_id", "gene_name", "synonyms")

    with open(path, encoding="utf-8-sig") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            columns = tuple(line.rstrip("\r\n").split("\t"))
            if not header_found:
                if columns != expected_header:
                    raise ValueError("unexpected PomBase gene names header")
                header_found = True
                continue
            if len(columns) != 3:
                raise ValueError(
                    f"gene names line {line_number} must contain three columns"
                )

            gene_id, name, synonym_text = columns
            if not gene_id:
                raise ValueError(f"gene names line {line_number} has no gene ID")
            if gene_id in seen_ids:
                raise ValueError(f"duplicate gene names ID {gene_id!r}")
            seen_ids.add(gene_id)
            synonyms = tuple(
                synonym.strip()
                for synonym in synonym_text.split(",")
                if synonym.strip()
            )
            records.append(GeneName(gene_id, name or None, synonyms))

    if not header_found:
        raise ValueError("PomBase gene names header was not found")
    return tuple(records)


@dataclass(frozen=True, slots=True, init=False, eq=False)
class GeneNames:
    """Read-only case-insensitive PomBase name index."""

    records: tuple[GeneName, ...]
    _by_id: object
    _by_primary: object
    _by_name: object

    def __init__(self, records):
        records = tuple(records)
        by_id = {
            record.gene_id.casefold(): record
            for record in records
        }
        by_primary = {}
        by_name = {}
        for record in records:
            if record.name:
                key = record.name.casefold()
                if key in by_primary:
                    raise ValueError(f"duplicate primary gene name {record.name!r}")
                by_primary[key] = record
            for identifier in (record.gene_id, record.name, *record.synonyms):
                if not identifier:
                    continue
                matches = by_name.setdefault(identifier.casefold(), {})
                matches[record.gene_id] = record

        object.__setattr__(self, "records", records)
        object.__setattr__(self, "_by_id", MappingProxyType(by_id))
        object.__setattr__(self, "_by_primary", MappingProxyType(by_primary))
        object.__setattr__(
            self,
            "_by_name",
            MappingProxyType(
                {
                    name: tuple(matches.values())
                    for name, matches in by_name.items()
                }
            ),
        )

    @classmethod
    def from_file(cls, path):
        return cls(read_gene_names(path))

    def find(self, name):
        normalized = str(name).strip().casefold()
        gene = self._by_id.get(normalized)
        if gene is not None:
            return gene
        gene = self._by_primary.get(normalized)
        if gene is not None:
            return gene
        matches = self._by_name.get(normalized, ())
        if not matches:
            return None
        if len(matches) > 1:
            raise AmbiguousGeneNameError(
                name,
                (match.gene_id for match in matches),
            )
        return matches[0]

    def __len__(self):
        return len(self.records)


def parse_attrs(text):
    """Parse a GFF3 attribute column."""
    attributes = {}
    for item in text.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            attributes[key] = value
    return MappingProxyType(attributes)


@dataclass(frozen=True, slots=True)
class RegionBlock:
    """Inclusive transcript interval."""

    feature_type: str
    start: int
    end: int

    @property
    def length(self):
        return self.end - self.start + 1


@dataclass(frozen=True, slots=True)
class CdsPosition:
    """Position in a spliced CDS."""

    base: int
    total: int

    @property
    def percent(self):
        return 100 * self.base / self.total


@dataclass(frozen=True, slots=True)
class BaseAnnotation:
    """Transcript feature covering one reference base."""

    block: RegionBlock
    cds_position: CdsPosition | None = None


@dataclass(frozen=True, slots=True)
class NeighborRegion:
    """Nearest different transcript block."""

    block: RegionBlock
    distance: int


@dataclass(frozen=True, slots=True)
class GeneAnnotation:
    """Gene identity and viability."""

    gene_id: str
    name: str | None
    gene_type: str | None
    start: int
    end: int
    strand: str
    viability: str | None
    chromosome: str | None = None

    @property
    def is_protein_coding(self):
        return self.gene_type == "protein_coding_gene"


@dataclass(frozen=True, slots=True)
class TranscriptCutAnnotation:
    """Transcript context at a Cas9 cut."""

    gene: GeneAnnotation
    transcript_id: str
    transcript_type: str
    strand: str
    left: BaseAnnotation | None
    right: BaseAnnotation | None
    block: RegionBlock | None
    lower_bases: int | None
    higher_bases: int | None
    upstream: NeighborRegion | None
    downstream: NeighborRegion | None
    upstream_direction: str
    downstream_direction: str

    @property
    def relation(self):
        return "within" if self.block is not None else "boundary"

    @property
    def cds_position(self):
        if self.left is not None and self.left.cds_position is not None:
            return self.left.cds_position
        if self.right is not None and self.right.cds_position is not None:
            return self.right.cds_position
        return None


@dataclass(frozen=True, slots=True)
class NearbyGene:
    """Nearest gene beside an intergenic cut."""

    gene: GeneAnnotation
    distance: int


@dataclass(frozen=True, slots=True)
class CutSiteAnnotation:
    """Annotation for a 1-based Cas9 cut boundary."""

    chromosome: str
    cut_coordinates: tuple[int, int]
    transcripts: tuple[TranscriptCutAnnotation, ...]
    lower_gene: NearbyGene | None = None
    higher_gene: NearbyGene | None = None

    def __post_init__(self):
        object.__setattr__(self, "cut_coordinates", tuple(self.cut_coordinates))
        object.__setattr__(self, "transcripts", tuple(self.transcripts))

    @property
    def is_intergenic(self):
        return not self.transcripts

    @property
    def genes(self):
        """Return unique overlapping genes in transcript order."""
        genes = []
        seen = set()
        for transcript in self.transcripts:
            if transcript.gene.gene_id not in seen:
                seen.add(transcript.gene.gene_id)
                genes.append(transcript.gene)
        return tuple(genes)


@dataclass(frozen=True, slots=True)
class _RawRecord:
    chromosome: str
    feature_type: str
    start: int
    end: int
    strand: str
    attributes: object


@dataclass(frozen=True, slots=True)
class _GeneRecord:
    gene_id: str
    chromosome: str
    start: int
    end: int
    strand: str
    name: str | None
    gene_type: str | None


@dataclass(frozen=True, slots=True)
class _TranscriptRecord:
    transcript_id: str
    gene_ids: tuple[str, ...]
    feature_type: str
    start: int
    end: int
    strand: str


def _parents(record):
    return tuple(
        parent
        for parent in record.attributes.get("Parent", "").split(",")
        if parent
    )


def _read_gff(path):
    records = []
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            columns = line.rstrip().split("\t")
            if len(columns) != 9:
                continue
            records.append(
                _RawRecord(
                    chromosome=columns[0],
                    feature_type=columns[2],
                    start=int(columns[3]),
                    end=int(columns[4]),
                    strand=columns[6],
                    attributes=parse_attrs(columns[8]),
                )
            )
    return records


def read_viability(path):
    """Read the PomBase gene viability table."""
    viability = {}
    with open(path, encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            columns = line.rstrip().split("\t")
            if len(columns) != 2:
                raise ValueError(
                    f"gene viability line {line_number} must contain two columns"
                )
            gene_id, status = columns
            if status not in VALID_VIABILITY_STATUSES:
                raise ValueError(
                    f"unsupported gene viability status {status!r} "
                    f"on line {line_number}"
                )
            if gene_id in viability:
                raise ValueError(f"duplicate gene viability ID {gene_id!r}")
            viability[gene_id] = status
    return MappingProxyType(viability)


@dataclass(frozen=True, slots=True, init=False, eq=False)
class GenomeAnnotations:
    """Read-only GFF3 and viability index."""

    _genes: object
    _transcripts: object
    _parts: object
    _viability: object
    _genes_by_chromosome: object
    _gene_ids_by_name: object
    _transcripts_by_gene: object

    def __init__(self, genes, transcripts, parts, viability=None):
        viability = {} if viability is None else viability
        object.__setattr__(self, "_genes", MappingProxyType(dict(genes)))
        object.__setattr__(
            self,
            "_transcripts",
            MappingProxyType(dict(transcripts)),
        )
        object.__setattr__(
            self,
            "_parts",
            MappingProxyType(
                {
                    transcript_id: tuple(items)
                    for transcript_id, items in parts.items()
                }
            ),
        )
        object.__setattr__(
            self,
            "_viability",
            MappingProxyType(dict(viability)),
        )

        genes_by_chromosome = {}
        for gene in self._genes.values():
            genes_by_chromosome.setdefault(gene.chromosome, []).append(gene)
        object.__setattr__(
            self,
            "_genes_by_chromosome",
            MappingProxyType(
                {
                    chromosome: tuple(
                        sorted(
                            items,
                            key=lambda gene: (
                                gene.start,
                                gene.end,
                                gene.gene_id,
                            ),
                        )
                    )
                    for chromosome, items in genes_by_chromosome.items()
                }
            ),
        )

        gene_ids_by_name = {}
        for gene in self._genes.values():
            for identifier in (gene.gene_id, gene.name):
                if identifier:
                    gene_ids_by_name.setdefault(
                        identifier.casefold(),
                        [],
                    ).append(gene.gene_id)
        object.__setattr__(
            self,
            "_gene_ids_by_name",
            MappingProxyType(
                {
                    name: tuple(dict.fromkeys(gene_ids))
                    for name, gene_ids in gene_ids_by_name.items()
                }
            ),
        )

        transcripts_by_gene = {}
        for transcript in self._transcripts.values():
            for gene_id in transcript.gene_ids:
                transcripts_by_gene.setdefault(gene_id, []).append(transcript)
        object.__setattr__(
            self,
            "_transcripts_by_gene",
            MappingProxyType(
                {
                    gene_id: tuple(
                        sorted(items, key=lambda item: item.transcript_id)
                    )
                    for gene_id, items in transcripts_by_gene.items()
                }
            ),
        )

    @classmethod
    def from_files(cls, gff_file, viability_file=None):
        """Load GFF3 and optional viability data."""
        records = _read_gff(gff_file)

        genes = {}
        for record in records:
            if record.feature_type != "gene":
                continue
            gene_id = record.attributes.get("ID")
            if gene_id is None:
                continue
            genes[gene_id] = _GeneRecord(
                gene_id=gene_id,
                chromosome=record.chromosome,
                start=record.start,
                end=record.end,
                strand=record.strand,
                name=record.attributes.get("Name"),
                gene_type=record.attributes.get("so_term_name"),
            )

        transcripts = {}
        for record in records:
            transcript_id = record.attributes.get("ID")
            gene_ids = tuple(
                parent for parent in _parents(record) if parent in genes
            )
            if transcript_id is None or not gene_ids:
                continue
            transcripts[transcript_id] = _TranscriptRecord(
                transcript_id=transcript_id,
                gene_ids=gene_ids,
                feature_type=record.feature_type,
                start=record.start,
                end=record.end,
                strand=record.strand,
            )

        parts = {transcript_id: [] for transcript_id in transcripts}
        for record in records:
            for parent in _parents(record):
                if parent in parts:
                    parts[parent].append(
                        RegionBlock(
                            feature_type=record.feature_type,
                            start=record.start,
                            end=record.end,
                        )
                    )

        viability = (
            read_viability(viability_file)
            if viability_file is not None
            else None
        )
        return cls(genes, transcripts, parts, viability)

    @property
    def viability(self):
        return self._viability

    def find_gene(self, name):
        """Find a gene by name or systematic ID."""
        normalized_name = str(name).strip().casefold()
        gene_ids = self._gene_ids_by_name.get(normalized_name, ())
        if not gene_ids:
            return None
        if len(gene_ids) > 1:
            matches = ", ".join(gene_ids)
            raise ValueError(
                f"gene name {name!r} is ambiguous; matches {matches}"
            )
        return self._gene(self._genes[gene_ids[0]])

    def annotate_cut(self, chromosome, cut_coordinates):
        """Annotate both bases flanking a Cas9 cut."""
        cut_left, cut_right = cut_coordinates
        if cut_right != cut_left + 1:
            raise ValueError("cut coordinates must identify adjacent reference bases")
        if cut_left < 1:
            raise ValueError("cut coordinates must be positive")

        left = self._base(chromosome, cut_left)
        right = self._base(chromosome, cut_right)
        contexts = []
        for key in sorted(set(left) | set(right)):
            left_annotation = left.get(key)
            right_annotation = right.get(key)
            representative = left_annotation or right_annotation
            gene, transcript = representative[0], representative[1]
            left_base = left_annotation[2] if left_annotation else None
            right_base = right_annotation[2] if right_annotation else None
            contexts.append(
                self._context(
                    cut_left,
                    gene,
                    transcript,
                    left_base,
                    right_base,
                )
            )

        lower_gene = None
        higher_gene = None
        if not contexts:
            lower_gene, higher_gene = self._near_genes(chromosome, cut_left)

        return CutSiteAnnotation(
            chromosome=chromosome,
            cut_coordinates=(cut_left, cut_right),
            transcripts=tuple(contexts),
            lower_gene=lower_gene,
            higher_gene=higher_gene,
        )

    def _gene(self, gene):
        return GeneAnnotation(
            gene_id=gene.gene_id,
            name=gene.name,
            gene_type=gene.gene_type,
            start=gene.start,
            end=gene.end,
            strand=gene.strand,
            viability=self._viability.get(gene.gene_id),
            chromosome=gene.chromosome,
        )

    def _base(self, chromosome, position):
        results = {}
        for gene in self._genes_by_chromosome.get(chromosome, ()):
            if gene.start > position:
                break
            if gene.end < position:
                continue
            for transcript in self._transcripts_by_gene.get(gene.gene_id, ()):
                if not (transcript.start <= position <= transcript.end):
                    continue
                block = self._feature(
                    position,
                    gene,
                    self._parts[transcript.transcript_id],
                )
                cds_position = None
                if block.feature_type == "CDS":
                    cds_position = self._cds(
                        position,
                        transcript,
                        self._parts[transcript.transcript_id],
                    )
                results[(gene.gene_id, transcript.transcript_id)] = (
                    gene,
                    transcript,
                    BaseAnnotation(block=block, cds_position=cds_position),
                )
        return results

    @staticmethod
    def _feature(position, gene, parts):
        overlapping = [
            part for part in parts if part.start <= position <= part.end
        ]
        for feature_type in REGION_PRIORITY:
            for part in overlapping:
                if part.feature_type == feature_type:
                    return part
        return RegionBlock("gene", gene.start, gene.end)

    @staticmethod
    def _cds(position, transcript, parts):
        segments = sorted(
            (part for part in parts if part.feature_type == "CDS"),
            key=lambda part: part.start,
            reverse=transcript.strand == "-",
        )
        total = sum(segment.length for segment in segments)
        passed = 0
        for segment in segments:
            if segment.start <= position <= segment.end:
                within = (
                    position - segment.start + 1
                    if transcript.strand == "+"
                    else segment.end - position + 1
                )
                return CdsPosition(base=passed + within, total=total)
            passed += segment.length
        return None

    def _context(self, cut_left, gene, transcript, left, right):
        regions = self._regions(
            self._parts[transcript.transcript_id]
        )
        same_block = (
            left is not None
            and right is not None
            and left.block == right.block
        )

        block = left.block if same_block else None
        lower_bases = cut_left - block.start + 1 if block else None
        higher_bases = block.end - cut_left if block else None

        if block is not None:
            lower = self._near_region(regions, cut_left, "lower", block)
            higher = self._near_region(regions, cut_left, "higher", block)
        else:
            lower = (
                NeighborRegion(left.block, 0)
                if left is not None
                else self._near_region(regions, cut_left, "lower")
            )
            higher = (
                NeighborRegion(right.block, 0)
                if right is not None
                else self._near_region(regions, cut_left, "higher")
            )

        if transcript.strand == "+":
            upstream, upstream_direction = lower, "lower"
            downstream, downstream_direction = higher, "higher"
        else:
            upstream, upstream_direction = higher, "higher"
            downstream, downstream_direction = lower, "lower"

        return TranscriptCutAnnotation(
            gene=self._gene(gene),
            transcript_id=transcript.transcript_id,
            transcript_type=transcript.feature_type,
            strand=transcript.strand,
            left=left,
            right=right,
            block=block,
            lower_bases=lower_bases,
            higher_bases=higher_bases,
            upstream=upstream,
            downstream=downstream,
            upstream_direction=upstream_direction,
            downstream_direction=downstream_direction,
        )

    @staticmethod
    def _regions(parts):
        regions = []
        seen = set()
        for part in parts:
            if part.feature_type not in REGION_PRIORITY:
                continue
            key = (part.feature_type, part.start, part.end)
            if key not in seen:
                seen.add(key)
                regions.append(part)
        return tuple(sorted(regions, key=lambda item: (item.start, item.end)))

    @staticmethod
    def _near_region(regions, cut_left, direction, current_block=None):
        if direction == "lower":
            if current_block is None:
                candidates = [region for region in regions if region.end <= cut_left]
            else:
                candidates = [
                    region for region in regions if region.end < current_block.start
                ]
            if not candidates:
                return None
            block = max(candidates, key=lambda item: (item.end, item.start))
            return NeighborRegion(block=block, distance=cut_left - block.end)

        if direction == "higher":
            if current_block is None:
                candidates = [
                    region for region in regions if region.start >= cut_left + 1
                ]
            else:
                candidates = [
                    region for region in regions if region.start > current_block.end
                ]
            if not candidates:
                return None
            block = min(candidates, key=lambda item: (item.start, item.end))
            return NeighborRegion(
                block=block,
                distance=block.start - cut_left - 1,
            )

        raise ValueError("direction must be lower or higher")

    def _near_genes(self, chromosome, cut_left):
        chromosome_genes = self._genes_by_chromosome.get(chromosome, ())
        lower = [gene for gene in chromosome_genes if gene.end <= cut_left]
        higher = [gene for gene in chromosome_genes if gene.start >= cut_left + 1]

        lower_result = None
        if lower:
            gene = max(lower, key=lambda item: item.end)
            lower_result = NearbyGene(
                gene=self._gene(gene),
                distance=cut_left - gene.end,
            )

        higher_result = None
        if higher:
            gene = min(higher, key=lambda item: item.start)
            higher_result = NearbyGene(
                gene=self._gene(gene),
                distance=gene.start - cut_left - 1,
            )

        return lower_result, higher_result
