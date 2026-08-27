"""CRISPR4P command-line interface."""

import argparse
import contextlib
import io
import sys
import time

from .crispr4p import PRECOMPUTED
from .service import (
    Crispr4pService,
    OligoLengthError,
    parse_oligo,
)


def build_parser():
    """Build the argument parser."""
    parser = argparse.ArgumentParser(description="crispr4p description")
    parser.add_argument("--name", action="store", type=str, help="Name")
    parser.add_argument(
        "-cr",
        "--chromosome",
        action="store",
        type=str,
        help="Chromosome",
    )
    parser.add_argument(
        "-co",
        "--coords",
        action="store",
        type=str,
        help="Coordinates",
    )
    parser.add_argument(
        "--mismatch",
        action="store",
        type=int,
        default=0,
        help="Allowed amount of mismatches.",
    )
    parser.add_argument(
        "--oligo",
        action="store",
        type=str,
        help="Oligo/sgRNA sequence (20bp seed or 23bp seed+PAM) to analyze.",
    )
    return parser


def format_design(result):
    """Format a design result."""
    output = io.StringIO()

    for index, row in enumerate(result.guide_table):
        print(index + 1, "-", row[0], row[2:], file=output)

        grna = row[1]
        print(
            "gRNA: ",
            grna[0],
            "PAM: %d - %d" % grna[3],
            grna[5],
            grna[4],
            file=output,
        )
        print("gRNAfw: ", grna[1], file=output)
        print("gRNArv: ", grna[2], "\n", file=output)

    print("HRfw: ", result.hr_dna[0] or "-", file=output)
    print("HRrv: ", result.hr_dna[1] or "-", file=output)
    print("Deleted DNA: ", result.hr_dna[2] or "-", "\n", file=output)

    if result.checking_primers:
        primer = result.checking_primers[0]
        print(
            "Check primer left: ",
            primer["PRIMER_LEFT_0_SEQUENCE"],
            "TM:",
            primer["PRIMER_LEFT_0_TM"],
            file=output,
        )
        print(
            "Check primer right: ",
            primer["PRIMER_RIGHT_0_SEQUENCE"],
            "TM:",
            primer["PRIMER_RIGHT_0_TM"],
            file=output,
        )
        print(
            "Deleted DNA product size: ",
            primer["PRIMER_PAIR_0_PRODUCT_SIZE"],
            file=output,
        )
        print(
            "Negative result product size: ",
            primer["negative_result"],
            "\n",
            file=output,
        )

    return output.getvalue()


def format_oligo(
    result,
    include_query_header=True,
    include_loading_messages=True,
):
    """Format an oligo result."""
    output = io.StringIO()
    if include_query_header:
        print(
            f"Querying S. pombe genome for oligo seed: {result.seed} "
            f"with PAM: {result.pam} "
            f"(Allowed mismatches: {result.n_mismatch})",
            file=output,
        )
    if include_loading_messages:
        print("Loading genome data and indexing PAM sites...", file=output)
        print("Genome indexed successfully.", file=output)

    print(
        "\nSummary of genome occurrences matching the seed sequence:",
        file=output,
    )
    print("-" * 65, file=output)
    print(
        f"{'Seed Length':12s} | "
        f"{'Matching Sites (containing NGG or NAG PAM)':40s}",
        file=output,
    )
    print("-" * 65, file=output)
    for length in (8, 10, 12, 14, 16, 18, 20):
        print(
            f"{length:12d} | {result.match_counts.get(length, 0)}",
            file=output,
        )
    print("-" * 65, file=output)

    if result.full_matches:
        print(
            f"\nDetails of {len(result.full_matches)} genomic "
            "target/off-target sites (full 20bp matches):",
            file=output,
        )
        for index, match in enumerate(result.full_matches):
            strand = "+" if match.strand == 1 else "-"
            print(
                f"  {index+1:2d}. Chromosome: {match.chromosome:4s} | "
                f"PAM coordinates: {match.pam_coordinates[0]} - "
                f"{match.pam_coordinates[1]} | "
                f"Cut: {match.cut_coordinates[0]} | "
                f"{match.cut_coordinates[1]} | "
                f"Strand: {strand} | Sequence: {match.seed} | "
                f"PAM: {match.pam}",
                file=output,
            )
    else:
        print("\nNo full 20bp matches found in the genome.", file=output)

    return output.getvalue()


def run_cli(arguments, service, output=None):
    """Run one command-line request."""
    if output is None:
        output = sys.stdout

    parser = build_parser()
    args = parser.parse_args(arguments)

    if args.oligo:
        try:
            normalized_sequence, seed, pam = parse_oligo(
                args.oligo
            )
        except OligoLengthError as error:
            print(
                "Error: Oligo sequence must be 20 bp (seed only) or 23 bp "
                f"(seed + PAM). Received length: {error.sequence_length}",
                file=output,
            )
            return 1

        print(
            f"Querying S. pombe genome for oligo seed: {seed} "
            f"with PAM: {pam} (Allowed mismatches: {args.mismatch})",
            file=output,
        )
        print("Loading genome data and indexing PAM sites...", file=output)
        output.flush()

        result = service.analyze_oligo(
            normalized_sequence,
            n_mismatch=args.mismatch,
        )
        print("Genome indexed successfully.", file=output)
        output.write(
            format_oligo(
                result,
                include_query_header=False,
                include_loading_messages=False,
            )
        )
        return 0

    if args.name:
        with contextlib.redirect_stdout(output):
            result = service.design_gene(
                args.name,
                n_mismatch=args.mismatch,
            )
    elif args.coords and args.chromosome:
        assert "..." in args.coords, "Coordinates need 3 dots in the middle"
        start, end = [value.strip() for value in args.coords.split("...")]
        with contextlib.redirect_stdout(output):
            result = service.design_region(
                args.chromosome,
                start,
                end,
                n_mismatch=args.mismatch,
            )
    else:
        parser.print_help(file=output)
        raise SystemExit()

    output.write(format_design(result))
    return 0
def main(arguments=None):
    """Run the command-line program."""
    if arguments is None:
        arguments = sys.argv[1:]

    start_time = time.time()
    service = Crispr4pService.from_project_data(
        precomputed_folder=PRECOMPUTED,
        designer_verbose=True,
    )
    status = run_cli(arguments, service)
    if status == 0:
        print("run time", time.time() - start_time)
    return status


if __name__ == "__main__":
    raise SystemExit(main())
