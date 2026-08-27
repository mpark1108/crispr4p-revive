"""Homology sequences and primer3 design."""

from dataclasses import dataclass, replace

from primer3 import bindings as primer3

if __package__:
    from .spedit import reverse_complement
else:
    from spedit import reverse_complement


def design_primers(sequence_args, global_args):
    """Call primer3-py, including its older API spelling."""
    design = getattr(primer3, "design_primers", None)
    if design is None:
        design = primer3.designPrimers
    return design(sequence_args, global_args)


class FlankError(ValueError):
    """Raised when a reference interval lacks the requested flanks."""


def _flanks(sequence, start, end, width, minimum):
    if not 0 <= start <= end <= len(sequence):
        raise ValueError("interval is outside the reference sequence")
    if width < 1:
        raise ValueError("flank width must be positive")
    left = sequence[max(0, start - width):start]
    right = sequence[end:min(len(sequence), end + width)]
    if len(left) < minimum or len(right) < minimum:
        raise FlankError("reference does not contain complete flanking sequence")
    return left, right


def build_hr_dna(
    chromosome_sequence,
    start,
    end,
    sequence_complement,
):
    """Build the original HR oligos and joined flanks."""
    previous_250, next_250 = _flanks(
        chromosome_sequence,
        start,
        end,
        250,
        80,
    )
    forward = previous_250[-80:] + next_250[:20]
    reverse = (
        "".join(reversed(next_250[:80]))
        + "".join(reversed(previous_250[-20:]))
    )
    reverse = sequence_complement(reverse)
    return forward, reverse, previous_250 + next_250


def _primer_settings(product_range):
    return {
        "PRIMER_OPT_SIZE": 20,
        "PRIMER_PICK_INTERNAL_OLIGO": 1,
        "PRIMER_INTERNAL_MAX_SELF_END": 8,
        "PRIMER_MIN_SIZE": 18,
        "PRIMER_MAX_SIZE": 25,
        "PRIMER_OPT_TM": 60.0,
        "PRIMER_MIN_TM": 57.0,
        "PRIMER_MAX_TM": 63.0,
        "PRIMER_MIN_GC": 20.0,
        "PRIMER_MAX_GC": 80.0,
        "PRIMER_MAX_POLY_X": 100,
        "PRIMER_INTERNAL_MAX_POLY_X": 100,
        "PRIMER_SALT_MONOVALENT": 50.0,
        "PRIMER_DNA_CONC": 50.0,
        "PRIMER_MAX_NS_ACCEPTED": 0,
        "PRIMER_MAX_SELF_ANY": 12,
        "PRIMER_MAX_SELF_END": 8,
        "PRIMER_PAIR_MAX_COMPL_ANY": 12,
        "PRIMER_PAIR_MAX_COMPL_END": 8,
        "PRIMER_PRODUCT_SIZE_RANGE": [list(product_range)],
    }


class PrimerNotFoundError(RuntimeError):
    """Raised when Primer3 cannot design a checking pair."""


@dataclass(frozen=True, slots=True)
class DonorOligos:
    """Two overlapping oligos used to build an HR template."""

    forward: str
    reverse: str
    overlap: str

    @property
    def product(self):
        reverse = reverse_complement(self.reverse)
        return self.forward + reverse[len(self.overlap):]

    @property
    def product_length(self):
        return len(self.product)


def overlap_oligos(sequence, overlap_start, overlap_end):
    """Split a donor around one shared sequence."""
    sequence = sequence.upper()
    if set(sequence) - set("ACGT"):
        raise ValueError("donor contains an invalid nucleotide")
    if not 0 <= overlap_start < overlap_end <= len(sequence):
        raise ValueError("overlap must be within the donor")

    return DonorOligos(
        forward=sequence[:overlap_end],
        reverse=reverse_complement(sequence[overlap_start:]),
        overlap=sequence[overlap_start:overlap_end],
    )


@dataclass(frozen=True, slots=True)
class InsertionPrimerPair:
    """One PCR pair spanning a cassette insertion."""

    forward: str
    reverse: str
    forward_tm: float
    reverse_tm: float
    wt_product_size: int
    insert_length: int
    forward_start: int | None = None
    reverse_end: int | None = None

    @property
    def disrupted_product_size(self):
        return self.wt_product_size + self.insert_length


@dataclass(frozen=True, slots=True)
class JunctionPrimerPair:
    forward: str
    reverse: str
    forward_tm: float
    reverse_tm: float
    product_size: int


@dataclass(frozen=True, slots=True)
class InsertionChecks:
    spanning: InsertionPrimerPair
    left: JunctionPrimerPair
    right: JunctionPrimerPair


def _check_pair(template, forward, reverse, product_size, primer_designer):
    sequence_args = {
        "SEQUENCE_ID": "junction_check",
        "SEQUENCE_TEMPLATE": template,
        "SEQUENCE_PRIMER": forward,
        "SEQUENCE_PRIMER_REVCOMP": reverse,
    }
    settings = _primer_settings((product_size, product_size))
    settings["PRIMER_TASK"] = "check_primers"
    settings["PRIMER_PICK_INTERNAL_OLIGO"] = 0
    settings["PRIMER_NUM_RETURN"] = 1

    answer = primer_designer(sequence_args, settings)
    if answer.get("PRIMER_PAIR_NUM_RETURNED", 0) < 1:
        raise PrimerNotFoundError("Primer3 could not check a junction pair")

    problem_keys = (
        "PRIMER_LEFT_0_PROBLEMS",
        "PRIMER_RIGHT_0_PROBLEMS",
        "PRIMER_PAIR_0_PROBLEMS",
    )
    problems = [str(answer[key]) for key in problem_keys if answer.get(key)]
    if problems:
        raise PrimerNotFoundError(
            "Primer3 rejected a junction primer: " + "; ".join(problems)
        )

    return JunctionPrimerPair(
        forward=answer["PRIMER_LEFT_0_SEQUENCE"],
        reverse=answer["PRIMER_RIGHT_0_SEQUENCE"],
        forward_tm=answer["PRIMER_LEFT_0_TM"],
        reverse_tm=answer["PRIMER_RIGHT_0_TM"],
        product_size=answer["PRIMER_PAIR_0_PRODUCT_SIZE"],
    )


def insertion_primers(
    reference,
    cut,
    arm_length=80,
    window=300,
    insert_length=23,
    primer_designer=design_primers,
):
    """Design one checking pair outside the insertion donor arms."""
    try:
        cut_left, cut_right = cut
    except (TypeError, ValueError):
        raise ValueError("cut coordinates must contain two bases") from None

    integers = (cut_left, cut_right, arm_length, window, insert_length)
    if any(not isinstance(value, int) or isinstance(value, bool)
           for value in integers):
        raise ValueError("primer coordinates and lengths must be integers")
    if cut_right != cut_left + 1:
        raise ValueError("cut coordinates must describe adjacent bases")
    if arm_length < 1 or window <= arm_length or insert_length < 1:
        raise ValueError("primer lengths and window are invalid")

    reference = reference.upper()
    if not 0 < cut_left < len(reference):
        raise ValueError("cut is outside the reference sequence")

    left = reference[max(0, cut_left - window):cut_left]
    right = reference[cut_left:min(len(reference), cut_left + window)]
    if len(left) < arm_length or len(right) < arm_length:
        raise ValueError("reference does not contain complete homology arms")

    template = left + right
    junction = len(left)
    right_start = junction + arm_length
    left_region = junction - arm_length
    right_region = len(template) - right_start
    sequence_args = {
        "SEQUENCE_ID": "insertion_check",
        "SEQUENCE_TEMPLATE": template,
        "SEQUENCE_INCLUDED_REGION": [0, len(template)],
        "SEQUENCE_EXCLUDED_REGION": [[left_region, 2 * arm_length]],
        "SEQUENCE_PRIMER_PAIR_OK_REGION_LIST": [
            [0, left_region, right_start, right_region]
        ],
    }
    answer = None
    for product_range in ((200, 300), (301, 350), (351, 400)):
        global_args = _primer_settings(product_range)
        global_args["PRIMER_PICK_INTERNAL_OLIGO"] = 0
        global_args["PRIMER_NUM_RETURN"] = 1
        answer = primer_designer(sequence_args, global_args)
        if answer.get("PRIMER_PAIR_NUM_RETURNED", 0) >= 1:
            break

    if answer.get("PRIMER_PAIR_NUM_RETURNED", 0) < 1:
        raise PrimerNotFoundError("Primer3 could not find a checking pair")

    left_position = answer.get("PRIMER_LEFT_0")
    right_position = answer.get("PRIMER_RIGHT_0")
    return InsertionPrimerPair(
        forward=answer["PRIMER_LEFT_0_SEQUENCE"],
        reverse=answer["PRIMER_RIGHT_0_SEQUENCE"],
        forward_tm=answer["PRIMER_LEFT_0_TM"],
        reverse_tm=answer["PRIMER_RIGHT_0_TM"],
        wt_product_size=answer["PRIMER_PAIR_0_PRODUCT_SIZE"],
        insert_length=insert_length,
        forward_start=left_position[0] if left_position else None,
        reverse_end=right_position[0] if right_position else None,
    )


def insertion_checks(
    reference,
    cut,
    insert,
    arm_length=80,
    window=300,
    primer_designer=design_primers,
    core=None,
    spanning=None,
):
    """Check the edit-spanning and cassette-junction PCR pairs."""
    insert = str(insert).upper()
    if not insert or set(insert) - set("ACGT"):
        raise ValueError("insert must contain only A, C, G, and T")
    core = insert if core is None else str(core).upper()
    if not core or set(core) - set("ACGT"):
        raise ValueError("cassette core must contain only A, C, G, and T")
    if insert.count(core) != 1:
        raise ValueError("insert must contain the cassette core once")

    if spanning is None:
        spanning = insertion_primers(
            reference,
            cut,
            arm_length=arm_length,
            window=window,
            insert_length=len(insert),
            primer_designer=primer_designer,
        )
    else:
        spanning = replace(spanning, insert_length=len(insert))
    if spanning.forward_start is None or spanning.reverse_end is None:
        raise PrimerNotFoundError("Primer3 did not return primer positions")

    cut_left = cut[0]
    reference = reference.upper()
    left = reference[max(0, cut_left - window):cut_left]
    right = reference[cut_left:min(len(reference), cut_left + window)]
    template = left + insert + right
    junction = len(left)
    core_start = insert.index(core)
    core_end = core_start + len(core)

    left_size = junction + core_end - spanning.forward_start
    right_size = (
        spanning.reverse_end
        + len(insert)
        - junction
        - core_start
        + 1
    )
    left_pair = _check_pair(
        template,
        spanning.forward,
        reverse_complement(core),
        left_size,
        primer_designer,
    )
    right_pair = _check_pair(
        template,
        core,
        spanning.reverse,
        right_size,
        primer_designer,
    )
    return InsertionChecks(spanning, left_pair, right_pair)


def checking_primers(
    chromosome_sequence,
    start,
    end,
    width,
    number_of_alternatives,
    primer_designer=design_primers,
):
    """Design checking primers around a region."""
    previous_sequence, next_sequence = _flanks(
        chromosome_sequence,
        start,
        end,
        width,
        80 + 18,
    )
    template_length = len(previous_sequence) + len(next_sequence)
    sequence_args = {
        "SEQUENCE_ID": "MH1000",
        "SEQUENCE_TEMPLATE": previous_sequence + next_sequence,
        "SEQUENCE_INCLUDED_REGION": [0, template_length],
        "SEQUENCE_EXCLUDED_REGION": [[len(previous_sequence) - 80, 160]],
    }
    global_args = _primer_settings(
        (min(width - 75, template_length), template_length)
    )

    primer3_result = primer_designer(sequence_args, global_args)
    result_keys = (
        "PRIMER_LEFT_%s_SEQUENCE",
        "PRIMER_LEFT_%s_SEQUENCE",
        "PRIMER_RIGHT_%s_SEQUENCE",
        "PRIMER_LEFT_%s_TM",
        "PRIMER_RIGHT_%s_TM",
        "PRIMER_LEFT_%s_GC_PERCENT",
        "PRIMER_RIGHT_%s_GC_PERCENT",
        "PRIMER_PAIR_%s_PRODUCT_SIZE",
        "PRIMER_LEFT_%s_TM",
        "PRIMER_RIGHT_%s_TM",
    )

    returned = primer3_result.get(
        "PRIMER_PAIR_NUM_RETURNED",
        number_of_alternatives,
    )
    alternatives = []
    for index in range(min(number_of_alternatives, returned)):
        keys = tuple(template % index for template in result_keys)
        if any(key not in primer3_result for key in keys):
            break
        alternative = {key: primer3_result[key] for key in keys}
        alternative["negative_result"] = (
            primer3_result["PRIMER_PAIR_%s_PRODUCT_SIZE" % index]
            + (end - start)
        )
        alternatives.append(alternative)
    return alternatives
