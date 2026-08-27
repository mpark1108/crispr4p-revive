"""PAM and cut-site coordinates."""


def slice_bounds(start, end):
    """Convert a 1-based inclusive interval to a Python slice."""
    if start < 1 or end < start:
        raise ValueError("invalid 1-based interval")
    return start - 1, end


def cut_from_pam(pam_start, pam_end, strand):
    """Convert a 1-based inclusive PAM to the bases flanking the cut."""
    if strand == 1:
        cut_left = pam_start - 4
    elif strand == -1:
        cut_left = pam_end + 3
    else:
        raise ValueError("Oligo hit strand must be 1 or -1")
    return cut_left, cut_left + 1


def hit_coordinates(
    hit_position,
    hit_strand,
    hit_pam,
    chromosome_sequence,
    reverse_complement,
):
    """Convert an indexed PAM hit to checked reference coordinates."""
    chromosome_length = len(chromosome_sequence)
    if hit_strand == 1:
        pam_start = hit_position
        pam_end = hit_position + 2
    elif hit_strand == -1:
        pam_start = chromosome_length - hit_position - 1
        pam_end = chromosome_length - hit_position + 1
    else:
        raise ValueError("Oligo hit strand must be 1 or -1")

    if pam_start < 1 or pam_end > chromosome_length:
        raise ValueError("Normalized PAM coordinates fall outside the FASTA")

    reference_pam = chromosome_sequence[pam_start - 1:pam_end]
    if hit_strand == -1:
        reference_pam = reverse_complement(reference_pam)
    if reference_pam != hit_pam:
        raise ValueError("Normalized PAM coordinates do not match the FASTA")

    pam_coordinates = (pam_start, pam_end)
    cut_coordinates = cut_from_pam(
        pam_start,
        pam_end,
        hit_strand,
    )
    return pam_coordinates, cut_coordinates
