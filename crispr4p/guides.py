"""Guide discovery and matching."""

import re


SUFFIX_LENGTHS = (8, 10, 12, 14, 16, 18, 20)


def find_guides(
    chromosome_sequence,
    chromosome_name,
    start,
    end,
    hit_factory,
    reverse_complement,
    seed_length=20,
):
    """Find unique NGG guides in a genomic interval."""
    interval_sequence = chromosome_sequence[start:end + 1]
    hits_by_seed = {}
    strands = (
        (1, interval_sequence),
        (-1, reverse_complement(interval_sequence)),
    )

    for strand, sequence in strands:
        # Lookahead retains both PAMs when the interval contains GGG.
        for match in re.finditer(r"(?=GG)", sequence):
            position = match.start()
            if position < seed_length + 1:
                continue
            pam = sequence[position - 1:position + 2]
            seed = sequence[
                position - seed_length - 1:position - 1
            ]
            if seed:
                hits_by_seed.setdefault(seed, []).append(
                    hit_factory(
                        chromosome_name,
                        position,
                        strand,
                        seed,
                        pam,
                    )
                )

    if not hits_by_seed:
        raise AssertionError("No nGG found in your input")

    return [hits[0] for hits in hits_by_seed.values() if len(hits) == 1]


def guide_primers(
    chromosome_sequence,
    start,
    end,
    guide,
    reverse_complement,
):
    """Build cloning oligos and PAM fields for a guide."""
    if guide.strand == 1:
        start_index = start + 1 + guide.pos
        guide_sequence = chromosome_sequence[
            start_index - 22:start_index - 2
        ]
        pam = chromosome_sequence[start_index - 2:start_index + 1]
        pam_start = start_index - 2
    else:
        start_index = end - 1 - guide.pos
        pam = reverse_complement(
            chromosome_sequence[start_index:start_index + 3]
        )
        guide_sequence = reverse_complement(
            chromosome_sequence[start_index + 3:start_index + 23]
        )
        pam_start = start_index

    forward_oligo = (
        guide_sequence[-10:] + "gttttagagctagaaatagcaagttaaaataa"
    )
    reverse_oligo = (
        reverse_complement(guide_sequence[:10])
        + "ttcttcggtacaggttatgttttttggcaaca"
    )

    # Input slices use zero-based indexes; user-visible PAM coordinates are
    # one-based and inclusive.
    pam_coordinates = (pam_start + 1, pam_start + 3)
    return (
        guide_sequence,
        forward_oligo,
        reverse_oligo,
        pam_coordinates,
        guide.strand,
        pam,
    )


def is_match(first_sequence, second_sequence, allowed_mismatches):
    """Compare two sequences with the original mismatch rule."""
    if allowed_mismatches == 0:
        return first_sequence == second_sequence

    mismatch_count = len(
        [
            position
            for position in range(len(first_sequence))
            if first_sequence[position] != second_sequence[position]
        ]
    )
    return allowed_mismatches >= mismatch_count


def match_guide(guide, genome_hits_by_suffix, n_mismatch):
    """Build cumulative suffix matches from 8 through 20 nt."""
    suffix = guide.seed[-8:]
    genome_hits = list(genome_hits_by_suffix.get(suffix, ()))
    match_table = {8: genome_hits}

    for suffix_length in SUFFIX_LENGTHS[1:]:
        remaining_hits = []
        for genome_hit in genome_hits:
            if is_match(
                guide.seed[-suffix_length:],
                genome_hit.seed[-suffix_length:],
                n_mismatch,
            ):
                remaining_hits.append(genome_hit)
        genome_hits = list(set(remaining_hits))
        match_table[suffix_length] = genome_hits

    return guide, match_table
