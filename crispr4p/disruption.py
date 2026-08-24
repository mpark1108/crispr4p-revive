"""Stop-cassette disruption design."""

import json
from dataclasses import dataclass

from .primers import overlap_oligos
from .spedit import reverse_complement


STOP_CODONS = frozenset({"TAA", "TAG", "TGA"})
PAM_SUFFIXES = frozenset({"GG", "AG"})
SITE_LENGTH = 23


def _codons(sequence, frame):
    return tuple(
        sequence[start:start + 3]
        for start in range(frame, len(sequence) - 2, 3)
    )


def _has_tandem_stops(codons):
    stops = [codon in STOP_CODONS for codon in codons]
    return any(left and right for left, right in zip(stops, stops[1:]))


def has_junction_pam(sequence):
    """Detect a PAM within the first six bases in either orientation."""
    orientations = (sequence, reverse_complement(sequence))
    return any(
        oriented[start + 1:start + 3] in PAM_SUFFIXES
        for oriented in orientations
        for start in range(4)
    )


@dataclass(frozen=True, slots=True)
class RecutSite:
    """A target in the edited junction that may retain first-guide activity."""

    target: str
    pam: str
    target_strand: str
    cassette_strand: str
    mismatches: int


@dataclass(frozen=True, slots=True)
class StopCassette:
    """One selectable 23-nt stop cassette."""

    id: int
    sequence: str

    def __post_init__(self):
        sequence = self.sequence.upper()
        if len(sequence) != 23:
            raise ValueError("stop cassette must contain 23 nt")
        if set(sequence) - set("ACGT"):
            raise ValueError("stop cassette contains an invalid nucleotide")
        if sequence[-2:] != "GG":
            raise ValueError("stop cassette must end with an NGG PAM")
        if has_junction_pam(sequence):
            raise ValueError("stop cassette creates a PAM beside the junction")
        valid_frames = all(
            _has_tandem_stops(_codons(sequence, frame))
            for frame in range(3)
        )
        if not valid_frames:
            raise ValueError(
                "stop cassette must contain tandem stops in every frame"
            )
        object.__setattr__(self, "sequence", sequence)

    @property
    def guide(self):
        return self.sequence[:20]

    @property
    def pam(self):
        return self.sequence[20:]

    @property
    def cut_bases(self):
        """Return the cassette bases on either side of the Cas9 cut."""
        cut_left = len(self.guide) - 3
        return cut_left, cut_left + 1

    @property
    def gc_percent(self):
        return 100 * sum(base in "GC" for base in self.guide) / len(self.guide)

    @property
    def frames(self):
        return tuple(_codons(self.sequence, frame) for frame in range(3))

    def orient(self, strand):
        if strand in (None, "+"):
            return self.sequence
        if strand == "-":
            return reverse_complement(self.sequence)
        raise ValueError("coding strand must be + or -")


@dataclass(frozen=True, slots=True)
class CassetteFormat:
    """An optional diagnostic extension for a 23-nt stop cassette."""

    id: str
    label: str
    tail: str = ""
    enzyme: str | None = None
    site: str | None = None
    cut_offset: int | None = None

    def __post_init__(self):
        tail = self.tail.upper()
        site = self.site.upper() if self.site else None
        if not self.id or not self.label:
            raise ValueError("cassette format requires an ID and label")
        if set(tail) - set("ACGT"):
            raise ValueError("cassette format contains an invalid nucleotide")
        diagnostic = (self.enzyme, site, self.cut_offset)
        missing = tuple(value is None for value in diagnostic)
        if any(missing) and not all(missing):
            raise ValueError("diagnostic format is incomplete")
        if site is not None and not 0 < self.cut_offset < len(site):
            raise ValueError("restriction cut must fall within its site")
        if site is not None and site != reverse_complement(site):
            raise ValueError("diagnostic recognition site must be palindromic")
        object.__setattr__(self, "tail", tail)
        object.__setattr__(self, "site", site)

    @property
    def length(self):
        return SITE_LENGTH + len(self.tail)

    def sequence(self, cassette):
        sequence = cassette.sequence + self.tail
        if self.site and sequence.count(self.site) != 1:
            raise ValueError(
                f"{self.enzyme} format must contain one recognition site"
            )
        return sequence

    def orient(self, cassette, strand):
        sequence = self.sequence(cassette)
        if strand in (None, "+"):
            return sequence
        if strand == "-":
            return reverse_complement(sequence)
        raise ValueError("coding strand must be + or -")


CASSETTE_FORMATS = (
    CassetteFormat("none", "23 nt — no diagnostic site"),
    CassetteFormat(
        "asci",
        "29 nt — AscI",
        tail="CGCGCC",
        enzyme="AscI",
        site="GGCGCGCC",
        cut_offset=2,
    ),
    CassetteFormat(
        "paci",
        "31 nt — PacI",
        tail="TTAATTAA",
        enzyme="PacI",
        site="TTAATTAA",
        cut_offset=5,
    ),
    CassetteFormat(
        "swai",
        "31 nt — SwaI",
        tail="ATTTAAAT",
        enzyme="SwaI",
        site="ATTTAAAT",
        cut_offset=4,
    ),
)
NO_DIAGNOSTIC = CASSETTE_FORMATS[0]


def cassette_format(format_id):
    """Return a cassette format by its stable ID."""
    format_id = str(format_id).strip().lower()
    for item in CASSETTE_FORMATS:
        if item.id == format_id:
            return item
    raise ValueError(f"unknown cassette format: {format_id}")


@dataclass(frozen=True, slots=True)
class DisruptionDonor:
    """A stop cassette with its locus-specific homology arms."""

    cassette: StopCassette
    cassette_format: CassetteFormat
    coding_strand: str
    arm_length: int
    left_arm: str
    insert: str
    right_arm: str

    @property
    def sequence(self):
        return self.left_arm + self.insert + self.right_arm

    @property
    def reverse(self):
        return reverse_complement(self.sequence)

    @property
    def total_length(self):
        return len(self.sequence)

    @property
    def oligos(self):
        core = self.cassette.orient(self.coding_strand)
        start = len(self.left_arm) + self.insert.index(core)
        return overlap_oligos(
            self.sequence,
            start,
            start + len(core),
        )


def _cut_index(reference, cut):
    cut_left, cut_right = cut
    if cut_right != cut_left + 1:
        raise ValueError("cut coordinates must describe adjacent bases")
    if not 0 < cut_left < len(reference):
        raise ValueError("cut is outside the reference sequence")
    return cut_left


def build_donor(
    reference,
    cut,
    cassette,
    coding_strand,
    arm_length,
    cassette_format=NO_DIAGNOSTIC,
):
    """Build a disruption donor around a Cas9 cut boundary."""
    reference = reference.upper()
    cut_index = _cut_index(reference, cut)
    if coding_strand not in ("+", "-"):
        raise ValueError("coding strand must be + or -")
    if (
        not isinstance(arm_length, int)
        or isinstance(arm_length, bool)
        or arm_length < 1
    ):
        raise ValueError("arm length must be a positive integer")
    if cut_index < arm_length or len(reference) - cut_index < arm_length:
        raise ValueError("reference does not contain complete homology arms")

    return DisruptionDonor(
        cassette=cassette,
        cassette_format=cassette_format,
        coding_strand=coding_strand,
        arm_length=arm_length,
        left_arm=reference[cut_index - arm_length:cut_index],
        insert=cassette_format.orient(cassette, coding_strand),
        right_arm=reference[cut_index:cut_index + arm_length],
    )


def _recut_sites(window, insert_start, insert_end, guide, cassette_strand,
                 max_mismatches):
    sites = []
    window_length = len(window)
    orientations = (
        ("+", window, insert_start, insert_end),
        (
            "-",
            reverse_complement(window),
            window_length - insert_end,
            window_length - insert_start,
        ),
    )

    for target_strand, sequence, start, end in orientations:
        for pam_start in range(20, len(sequence) - 2):
            pam = sequence[pam_start:pam_start + 3]
            if pam[1:] not in PAM_SUFFIXES:
                continue

            site_start = pam_start - 20
            site_end = pam_start + 3
            if site_end <= start or site_start >= end:
                continue

            target = sequence[site_start:pam_start]
            mismatches = sum(
                left != right for left, right in zip(guide, target)
            )
            if mismatches <= max_mismatches:
                sites.append(
                    RecutSite(
                        target=target,
                        pam=pam,
                        target_strand=target_strand,
                        cassette_strand=cassette_strand,
                        mismatches=mismatches,
                    )
                )
    return sites


def recut_sites(reference, cut, guide, cassette, coding_strand=None,
                max_mismatches=4, cassette_format=NO_DIAGNOSTIC):
    """Find NGG/NAG junction targets similar to the first guide."""
    reference = reference.upper()
    guide = guide.upper()
    cut_left = _cut_index(reference, cut)
    if len(guide) != 20 or set(guide) - set("ACGT"):
        raise ValueError("guide must contain 20 DNA bases")
    if coding_strand not in (None, "+", "-"):
        raise ValueError("coding strand must be + or -")

    flank = SITE_LENGTH
    left = reference[max(0, cut_left - flank):cut_left]
    right = reference[cut_left:min(len(reference), cut_left + flank)]
    strands = ("+", "-") if coding_strand is None else (coding_strand,)
    sites = []

    for cassette_strand in strands:
        if cassette_format is NO_DIAGNOSTIC:
            insert = cassette.orient(cassette_strand)
        else:
            insert = cassette_format.orient(cassette, cassette_strand)
        window = left + insert + right
        insert_start = len(left)
        insert_end = insert_start + len(insert)
        sites.extend(
            _recut_sites(
                window,
                insert_start,
                insert_end,
                guide,
                cassette_strand,
                max_mismatches,
            )
        )
    return tuple(sites)


def load_cassettes(path):
    """Load the packaged stop-cassette catalog."""
    with open(path, encoding="utf-8") as handle:
        data = json.load(handle)

    if data.get("version") != 1:
        raise ValueError("unsupported stop-cassette data version")
    if data.get("assembly") != "ASM294v2.26":
        raise ValueError("stop-cassette assembly does not match CRISPR4P")

    cassettes = tuple(
        StopCassette(record["id"], record["sequence"])
        for record in data.get("candidates", ())
    )
    ids = tuple(cassette.id for cassette in cassettes)
    if ids != tuple(range(1, len(cassettes) + 1)):
        raise ValueError("stop-cassette IDs must be consecutive")
    if len({cassette.sequence for cassette in cassettes}) != len(cassettes):
        raise ValueError("stop-cassette sequences must be unique")
    return cassettes


def target_gene(annotation, target_name=None):
    """Return the named target or the only overlapping coding gene."""
    target = str(target_name).strip().casefold() if target_name else None
    if target:
        for gene in annotation.genes:
            identifiers = {gene.gene_id.casefold()}
            if gene.name:
                identifiers.add(gene.name.casefold())
            if target in identifiers:
                return gene

    coding_genes = [gene for gene in annotation.genes if gene.is_protein_coding]
    if len(coding_genes) == 1:
        return coding_genes[0]
    return None


def target_strand(annotation, target_name=None):
    """Return the coding strand when one target gene is identifiable."""
    gene = target_gene(annotation, target_name)
    return gene.strand if gene is not None and gene.is_protein_coding else None
