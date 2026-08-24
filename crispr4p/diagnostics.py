"""Restriction-digest checks for a cassette insertion."""

from dataclasses import dataclass

from .disruption import CassetteFormat
from .primers import InsertionChecks, InsertionPrimerPair


@dataclass(frozen=True, slots=True)
class DigestResult:
    enzyme: str
    site: str
    wt_sites: tuple[int, ...]
    edited_sites: tuple[int, ...]
    fragments: tuple[int, ...]

    @property
    def available(self):
        return len(self.wt_sites) == 0 and len(self.edited_sites) == 1


@dataclass(frozen=True, slots=True)
class CassetteOption:
    cassette_format: CassetteFormat
    coding_sequence: str
    insert: str
    available: bool
    spanning: InsertionPrimerPair | None
    checks: InsertionChecks | None
    digest: DigestResult | None


def find_sites(sequence, site):
    """Return overlapping recognition-site starts."""
    starts = []
    start = 0
    while True:
        start = sequence.find(site, start)
        if start < 0:
            return tuple(starts)
        starts.append(start)
        start += 1


def digest_fragments(length, sites, cut_offset):
    cuts = sorted(start + cut_offset for start in sites)
    boundaries = (0, *cuts, length)
    return tuple(
        sorted(
            (right - left for left, right in zip(boundaries, boundaries[1:])),
            reverse=True,
        )
    )


def pcr_rflp(reference, cut, insert, pair, cassette_format, window=300):
    """Evaluate one enzyme against the predicted PCR products."""
    if cassette_format.enzyme is None:
        raise ValueError("cassette format has no diagnostic enzyme")
    if pair.forward_start is None or pair.reverse_end is None:
        raise ValueError("primer positions are required for PCR-RFLP")

    cut_left, cut_right = cut
    if cut_right != cut_left + 1:
        raise ValueError("cut coordinates must describe adjacent bases")

    reference = reference.upper()
    left = reference[max(0, cut_left - window):cut_left]
    right = reference[cut_left:min(len(reference), cut_left + window)]
    wt_template = left + right
    edited_template = left + insert + right
    start = pair.forward_start
    wt_end = pair.reverse_end + 1
    edited_end = wt_end + len(insert)
    wt_product = wt_template[start:wt_end]
    edited_product = edited_template[start:edited_end]

    if len(wt_product) != pair.wt_product_size:
        raise ValueError("WT amplicon length does not match Primer3")
    if len(edited_product) != pair.disrupted_product_size:
        raise ValueError("edited amplicon length does not match Primer3")

    site = cassette_format.site
    wt_sites = find_sites(wt_product, site)
    edited_sites = find_sites(edited_product, site)
    fragments = digest_fragments(
        len(edited_product),
        edited_sites,
        cassette_format.cut_offset,
    )
    return DigestResult(
        enzyme=cassette_format.enzyme,
        site=site,
        wt_sites=wt_sites,
        edited_sites=edited_sites,
        fragments=fragments,
    )
