import json
import unittest
from io import BytesIO
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock, patch

import webapp
from crispr4p.annotations import GeneAnnotation, GenomeAnnotations
from crispr4p.diagnostics import (
    CassetteOption,
    DigestResult,
    digest_fragments,
    find_sites,
    pcr_rflp,
)
from crispr4p.disruption import (
    CASSETTE_FORMATS,
    NO_DIAGNOSTIC,
    CassetteFormat,
    StopCassette,
    build_donor,
    cassette_format,
    has_junction_pam,
    load_cassettes,
    recut_sites,
    target_gene,
    target_strand,
    valid_rescue_site,
)
from crispr4p.models import DesignResult
from crispr4p.primers import (
    InsertionChecks,
    InsertionPrimerPair,
    JunctionPrimerPair,
    PrimerNotFoundError,
    insertion_checks,
    insertion_primers,
    overlap_oligos,
)
from crispr4p.resources import read_fasta
from crispr4p.service import Crispr4pService
from crispr4p.spedit import reverse_complement
from crispr4p.web_views import (
    annotation_rows,
    cassette_data,
    donor_data,
    render_design,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DATA = PROJECT_ROOT / "data"
CANDIDATES = DATA / "stop_cassettes.json"
FASTA = DATA / "Schizosaccharomyces_pombe.ASM294v2.26.dna.toplevel.fa"
FIRST = "TGATGAGCTGATGACCTAGTAGG"
OLD_FIRST = "TAGTAGCCTAGTGACCTAGTAGG"
ADE6_GUIDE = "ACATTGGCTTACGACGGTCG"
ADE6_REVERSE_GUIDE = "GTGGCGACAGGGACACCTCG"
BUB1_GUIDE = "TATCAGATTGCTCGGCCACA"


def gene(gene_id, name, strand, gene_type="protein_coding_gene"):
    return GeneAnnotation(
        gene_id=gene_id,
        name=name,
        gene_type=gene_type,
        start=100,
        end=300,
        strand=strand,
        viability="viable",
    )


def unsafe_cassette():
    return SimpleNamespace(
        orient=lambda strand: (
            reverse_complement(OLD_FIRST) if strand == "-" else OLD_FIRST
        )
    )


class DisruptionDesignTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cassettes = load_cassettes(CANDIDATES)
        records = read_fasta(FASTA)
        cls.reference = records["III"].sequence
        cls.references = tuple(record.sequence for record in records.values())

    def test_packaged_candidates(self):
        self.assertEqual(tuple(range(1, 11)), tuple(c.id for c in self.cassettes))
        self.assertFalse(any(has_junction_pam(c.sequence) for c in self.cassettes))
        self.assertEqual(FIRST, self.cassettes[0].sequence)
        self.assertEqual(
            "GTGATAGCTGATGACTAGTGAGG",
            self.cassettes[-1].sequence,
        )

        data = json.loads(CANDIDATES.read_text(encoding="utf-8"))
        self.assertEqual(1, data["version"])
        self.assertEqual("ASM294v2.26", data["assembly"])
        for cassette in self.cassettes:
            for reference in self.references:
                self.assertNotIn(cassette.sequence, reference)
                self.assertNotIn(
                    reverse_complement(cassette.sequence),
                    reference,
                )

    def test_candidate_details(self):
        cassette = self.cassettes[0]

        self.assertEqual("TGATGAGCTGATGACCTAGT", cassette.guide)
        self.assertEqual("AGG", cassette.pam)
        self.assertEqual(45.0, cassette.gc_percent)
        self.assertEqual(3, len(cassette.frames))
        self.assertEqual(
            "CCTACTAGGTCATCAGCTCATCA",
            cassette.orient("-"),
        )
        self.assertEqual(FIRST, cassette.orient("+"))

    def test_cassette_formats(self):
        expected = {
            "none": (23, None, None, None),
            "asci": (29, "AscI", "GGCGCGCC", 2),
            "paci": (31, "PacI", "TTAATTAA", 5),
            "swai": (31, "SwaI", "ATTTAAAT", 4),
        }

        self.assertIs(NO_DIAGNOSTIC, cassette_format("none"))
        with self.assertRaisesRegex(ValueError, "unknown cassette format"):
            cassette_format("other")

        for cassette in self.cassettes:
            for item in CASSETTE_FORMATS:
                sequence = item.sequence(cassette)
                length, enzyme, site, cut_offset = expected[item.id]
                self.assertEqual(length, len(sequence))
                self.assertEqual(length, item.length)
                self.assertEqual(enzyme, item.enzyme)
                self.assertEqual(cut_offset, item.cut_offset)
                self.assertEqual(sequence, item.orient(cassette, "+"))
                self.assertEqual(
                    reverse_complement(sequence),
                    item.orient(cassette, "-"),
                )
                if site:
                    self.assertEqual(1, sequence.count(site))

    def test_extended_donor_oligos_keep_core_overlap(self):
        cassette = self.cassettes[0]
        expected = {
            "none": (183, 103, 103),
            "asci": (189, 103, 109),
            "paci": (191, 103, 111),
            "swai": (191, 103, 111),
        }

        for item in CASSETTE_FORMATS:
            plus = build_donor(
                self.reference,
                (1316791, 1316792),
                cassette,
                "+",
                80,
                item,
            )
            minus = build_donor(
                self.reference,
                (1316791, 1316792),
                cassette,
                "-",
                80,
                item,
            )
            total, short, long = expected[item.id]

            self.assertEqual(item.sequence(cassette), plus.insert)
            self.assertEqual(total, plus.total_length)
            self.assertEqual((short, long), (
                len(plus.oligos.forward),
                len(plus.oligos.reverse),
            ))
            self.assertEqual((long, short), (
                len(minus.oligos.forward),
                len(minus.oligos.reverse),
            ))
            self.assertEqual(cassette.sequence, plus.oligos.overlap)
            self.assertEqual(
                reverse_complement(cassette.sequence),
                minus.oligos.overlap,
            )
            self.assertEqual(plus.sequence, plus.oligos.product)
            self.assertEqual(minus.sequence, minus.oligos.product)

    def test_invalid_candidate(self):
        with self.assertRaisesRegex(ValueError, "23 nt"):
            StopCassette(1, "TAG")
        with self.assertRaisesRegex(ValueError, "tandem stops"):
            StopCassette(1, "ACGTACGTACGTACGTACGTAGG")
        with self.assertRaisesRegex(ValueError, "beside the junction"):
            StopCassette(1, OLD_FIRST)

    def test_junction_recut(self):
        sites = recut_sites(
            self.reference,
            (1316791, 1316792),
            ADE6_GUIDE,
            unsafe_cassette(),
            "+",
        )

        self.assertTrue(has_junction_pam(OLD_FIRST))
        self.assertEqual(1, len(sites))
        self.assertEqual("ACATTGGCTTACGACGGTAG", sites[0].target)
        self.assertEqual("TAG", sites[0].pam)
        self.assertEqual(1, sites[0].mismatches)
        reverse_sites = recut_sites(
            self.reference,
            (1317806, 1317807),
            ADE6_REVERSE_GUIDE,
            unsafe_cassette(),
            "-",
        )
        self.assertEqual("-", reverse_sites[0].target_strand)
        self.assertEqual(1, reverse_sites[0].mismatches)
        self.assertFalse(
            recut_sites(
                self.reference,
                (1316791, 1316792),
                ADE6_GUIDE,
                self.cassettes[0],
            )
        )

    def test_rescue_site_is_unique(self):
        cassette = self.cassettes[0]
        cut = (1316791, 1316792)

        self.assertTrue(
            valid_rescue_site(self.reference, cut, cassette, "+")
        )
        self.assertTrue(
            valid_rescue_site(self.reference, cut, cassette, "-")
        )

        guide = list(cassette.guide)
        guide[0] = "A" if guide[0] != "A" else "C"
        extra_target = CassetteFormat(
            "extra",
            "extra rescue target",
            tail="".join(guide) + cassette.pam,
        )
        self.assertFalse(
            valid_rescue_site(
                self.reference,
                cut,
                cassette,
                "+",
                cassette_format=extra_target,
            )
        )

    def test_disruption_donor(self):
        donor = build_donor(
            self.reference,
            (1316791, 1316792),
            self.cassettes[0],
            "+",
            80,
        )

        self.assertEqual(self.reference[1316711:1316791], donor.left_arm)
        self.assertEqual(FIRST, donor.insert)
        self.assertEqual(self.reference[1316791:1316871], donor.right_arm)
        self.assertEqual(183, donor.total_length)
        self.assertEqual(reverse_complement(donor.sequence), donor.reverse)

        oligos = donor.oligos
        self.assertEqual(donor.left_arm + donor.insert, oligos.forward)
        self.assertEqual(
            reverse_complement(donor.insert + donor.right_arm),
            oligos.reverse,
        )
        self.assertEqual(donor.insert, oligos.overlap)
        self.assertEqual(
            reverse_complement(oligos.overlap),
            oligos.reverse[-len(oligos.overlap):],
        )
        self.assertEqual((103, 103), (len(oligos.forward), len(oligos.reverse)))
        self.assertEqual(donor.sequence, oligos.product)
        self.assertEqual(183, oligos.product_length)

        minus = build_donor(
            self.reference,
            (1316791, 1316792),
            self.cassettes[0],
            "-",
            80,
        )
        self.assertEqual(reverse_complement(FIRST), minus.insert)
        self.assertEqual(minus.insert, minus.oligos.overlap)
        self.assertEqual(minus.sequence, minus.oligos.product)
        with self.assertRaisesRegex(ValueError, "positive integer"):
            build_donor(
                self.reference,
                (1316791, 1316792),
                self.cassettes[0],
                "+",
                0,
            )
        with self.assertRaisesRegex(ValueError, "complete homology arms"):
            build_donor(
                self.reference,
                (40, 41),
                self.cassettes[0],
                "+",
                80,
            )

    def test_construction_oligos_cover_every_cassette_orientation(self):
        for cassette in self.cassettes:
            for strand in ("+", "-"):
                donor = build_donor(
                    self.reference,
                    (1316791, 1316792),
                    cassette,
                    strand,
                    80,
                )
                oligos = donor.oligos
                self.assertEqual(
                    (103, 103, 23),
                    (
                        len(oligos.forward),
                        len(oligos.reverse),
                        len(oligos.overlap),
                    ),
                )
                self.assertEqual(donor.insert, oligos.overlap)
                self.assertEqual(donor.sequence, oligos.product)

    def test_construction_oligo_validation(self):
        with self.assertRaisesRegex(ValueError, "invalid nucleotide"):
            overlap_oligos("AACN", 1, 3)
        with self.assertRaisesRegex(ValueError, "within the donor"):
            overlap_oligos("AACG", 2, 2)

    def test_insertion_primer_inputs(self):
        answer = {
            "PRIMER_PAIR_NUM_RETURNED": 1,
            "PRIMER_LEFT_0_SEQUENCE": "FORWARD",
            "PRIMER_RIGHT_0_SEQUENCE": "REVERSE",
            "PRIMER_LEFT_0_TM": 59.4,
            "PRIMER_RIGHT_0_TM": 60.2,
            "PRIMER_PAIR_0_PRODUCT_SIZE": 250,
            "PRIMER_LEFT_0": [130, 20],
            "PRIMER_RIGHT_0": [379, 20],
        }
        designer = Mock(return_value=answer)
        pair = insertion_primers(
            "A" * 300 + "C" * 300,
            (300, 301),
            primer_designer=designer,
        )
        sequence_args, global_args = designer.call_args.args

        self.assertEqual(
            "A" * 300 + "C" * 300,
            sequence_args["SEQUENCE_TEMPLATE"],
        )
        self.assertEqual(
            [[220, 160]],
            sequence_args["SEQUENCE_EXCLUDED_REGION"],
        )
        self.assertEqual(
            [[0, 220, 380, 220]],
            sequence_args["SEQUENCE_PRIMER_PAIR_OK_REGION_LIST"],
        )
        self.assertEqual(
            [[200, 300]],
            global_args["PRIMER_PRODUCT_SIZE_RANGE"],
        )
        self.assertEqual(0, global_args["PRIMER_PICK_INTERNAL_OLIGO"])
        self.assertEqual(1, global_args["PRIMER_NUM_RETURN"])
        self.assertEqual("FORWARD", pair.forward)
        self.assertEqual("REVERSE", pair.reverse)
        self.assertEqual(250, pair.wt_product_size)
        self.assertEqual(273, pair.disrupted_product_size)
        self.assertEqual(130, pair.forward_start)
        self.assertEqual(379, pair.reverse_end)

        no_pair = Mock(return_value={"PRIMER_PAIR_NUM_RETURNED": 0})
        with self.assertRaises(PrimerNotFoundError):
            insertion_primers(
                "A" * 300 + "C" * 300,
                (300, 301),
                primer_designer=no_pair,
            )
        self.assertEqual(3, no_pair.call_count)

        fallback_answer = dict(answer)
        fallback_answer["PRIMER_PAIR_0_PRODUCT_SIZE"] = 328
        fallback_answer["PRIMER_RIGHT_0"] = [457, 20]
        fallback = Mock(side_effect=[
            {"PRIMER_PAIR_NUM_RETURNED": 0},
            fallback_answer,
        ])
        pair = insertion_primers(
            "A" * 300 + "C" * 300,
            (300, 301),
            primer_designer=fallback,
        )
        self.assertEqual(328, pair.wt_product_size)
        self.assertEqual(
            [[200, 300]],
            fallback.call_args_list[0].args[1]["PRIMER_PRODUCT_SIZE_RANGE"],
        )
        self.assertEqual(
            [[301, 350]],
            fallback.call_args_list[1].args[1]["PRIMER_PRODUCT_SIZE_RANGE"],
        )

    def test_junction_primer_inputs(self):
        spanning = {
            "PRIMER_PAIR_NUM_RETURNED": 1,
            "PRIMER_LEFT_0_SEQUENCE": "FORWARD",
            "PRIMER_RIGHT_0_SEQUENCE": "REVERSE",
            "PRIMER_LEFT_0_TM": 59.4,
            "PRIMER_RIGHT_0_TM": 60.2,
            "PRIMER_PAIR_0_PRODUCT_SIZE": 298,
            "PRIMER_LEFT_0": [130, 20],
            "PRIMER_RIGHT_0": [427, 20],
        }
        left = {
            "PRIMER_PAIR_NUM_RETURNED": 1,
            "PRIMER_LEFT_0_SEQUENCE": "FORWARD",
            "PRIMER_RIGHT_0_SEQUENCE": reverse_complement(FIRST),
            "PRIMER_LEFT_0_TM": 59.4,
            "PRIMER_RIGHT_0_TM": 59.1,
            "PRIMER_PAIR_0_PRODUCT_SIZE": 193,
        }
        right = {
            "PRIMER_PAIR_NUM_RETURNED": 1,
            "PRIMER_LEFT_0_SEQUENCE": FIRST,
            "PRIMER_RIGHT_0_SEQUENCE": "REVERSE",
            "PRIMER_LEFT_0_TM": 59.1,
            "PRIMER_RIGHT_0_TM": 60.2,
            "PRIMER_PAIR_0_PRODUCT_SIZE": 151,
        }
        designer = Mock(side_effect=(spanning, left, right))

        checks = insertion_checks(
            "A" * 300 + "C" * 300,
            (300, 301),
            FIRST,
            primer_designer=designer,
        )

        self.assertEqual(3, designer.call_count)
        left_args, left_settings = designer.call_args_list[1].args
        right_args, right_settings = designer.call_args_list[2].args
        edited = "A" * 300 + FIRST + "C" * 300
        self.assertEqual(edited, left_args["SEQUENCE_TEMPLATE"])
        self.assertEqual("FORWARD", left_args["SEQUENCE_PRIMER"])
        self.assertEqual(
            reverse_complement(FIRST),
            left_args["SEQUENCE_PRIMER_REVCOMP"],
        )
        self.assertEqual(FIRST, right_args["SEQUENCE_PRIMER"])
        self.assertEqual("REVERSE", right_args["SEQUENCE_PRIMER_REVCOMP"])
        self.assertEqual("check_primers", left_settings["PRIMER_TASK"])
        self.assertEqual(57.0, left_settings["PRIMER_MIN_TM"])
        self.assertEqual(63.0, left_settings["PRIMER_MAX_TM"])
        self.assertEqual(
            [[193, 193]],
            left_settings["PRIMER_PRODUCT_SIZE_RANGE"],
        )
        self.assertEqual(
            [[151, 151]],
            right_settings["PRIMER_PRODUCT_SIZE_RANGE"],
        )
        self.assertEqual(193, checks.left.product_size)
        self.assertEqual(151, checks.right.product_size)

        failed = Mock(side_effect=(spanning, {"PRIMER_PAIR_NUM_RETURNED": 0}))
        with self.assertRaises(PrimerNotFoundError):
            insertion_checks(
                "A" * 300 + "C" * 300,
                (300, 301),
                FIRST,
                primer_designer=failed,
            )

        problem = dict(left)
        problem["PRIMER_RIGHT_0_PROBLEMS"] = "high hairpin stability"
        rejected = Mock(side_effect=(spanning, problem))
        with self.assertRaisesRegex(PrimerNotFoundError, "high hairpin"):
            insertion_checks(
                "A" * 300 + "C" * 300,
                (300, 301),
                FIRST,
                primer_designer=rejected,
            )

    def test_extended_insert_uses_core_junction_primer(self):
        spanning = {
            "PRIMER_PAIR_NUM_RETURNED": 1,
            "PRIMER_LEFT_0_SEQUENCE": "FORWARD",
            "PRIMER_RIGHT_0_SEQUENCE": "REVERSE",
            "PRIMER_LEFT_0_TM": 59.4,
            "PRIMER_RIGHT_0_TM": 60.2,
            "PRIMER_PAIR_0_PRODUCT_SIZE": 298,
            "PRIMER_LEFT_0": [130, 20],
            "PRIMER_RIGHT_0": [427, 20],
        }
        checked = {
            "PRIMER_PAIR_NUM_RETURNED": 1,
            "PRIMER_LEFT_0_SEQUENCE": "FORWARD",
            "PRIMER_RIGHT_0_SEQUENCE": "REVERSE",
            "PRIMER_LEFT_0_TM": 59.4,
            "PRIMER_RIGHT_0_TM": 59.1,
            "PRIMER_PAIR_0_PRODUCT_SIZE": 1,
        }
        asci = cassette_format("asci")
        full_insert = asci.sequence(self.cassettes[0])
        designer = Mock(side_effect=(spanning, checked, checked))

        checks = insertion_checks(
            "A" * 300 + "C" * 300,
            (300, 301),
            full_insert,
            core=FIRST,
            primer_designer=designer,
        )
        left_args, left_settings = designer.call_args_list[1].args
        right_args, right_settings = designer.call_args_list[2].args

        self.assertEqual(327, checks.spanning.disrupted_product_size)
        self.assertEqual(
            reverse_complement(FIRST),
            left_args["SEQUENCE_PRIMER_REVCOMP"],
        )
        self.assertEqual(FIRST, right_args["SEQUENCE_PRIMER"])
        self.assertEqual([[193, 193]], left_settings["PRIMER_PRODUCT_SIZE_RANGE"])
        self.assertEqual([[157, 157]], right_settings["PRIMER_PRODUCT_SIZE_RANGE"])

        reverse_insert = reverse_complement(full_insert)
        reverse_core = reverse_complement(FIRST)
        reverse_designer = Mock(side_effect=(spanning, checked, checked))
        insertion_checks(
            "A" * 300 + "C" * 300,
            (300, 301),
            reverse_insert,
            core=reverse_core,
            primer_designer=reverse_designer,
        )
        reverse_left = reverse_designer.call_args_list[1].args[1]
        reverse_right = reverse_designer.call_args_list[2].args[1]
        self.assertEqual(
            [[199, 199]],
            reverse_left["PRIMER_PRODUCT_SIZE_RANGE"],
        )
        self.assertEqual(
            [[151, 151]],
            reverse_right["PRIMER_PRODUCT_SIZE_RANGE"],
        )

    def test_diagnostic_digest(self):
        pair = InsertionPrimerPair(
            forward="FORWARD",
            reverse="REVERSE",
            forward_tm=59.4,
            reverse_tm=60.2,
            wt_product_size=298,
            insert_length=29,
            forward_start=130,
            reverse_end=427,
        )
        asci = cassette_format("asci")
        insert = asci.sequence(self.cassettes[0])
        result = pcr_rflp(
            "A" * 300 + "C" * 300,
            (300, 301),
            insert,
            pair,
            asci,
        )

        self.assertTrue(result.available)
        self.assertEqual((), result.wt_sites)
        self.assertEqual((191,), result.edited_sites)
        self.assertEqual((193, 134), result.fragments)
        self.assertEqual((0, 4), find_sites("TTAATTAATTAA", "TTAATTAA"))
        self.assertEqual((60, 40), digest_fragments(100, (38,), 2))

        reference = list("A" * 300 + "C" * 300)
        reference[200:208] = asci.site
        conflict = pcr_rflp(
            "".join(reference),
            (300, 301),
            insert,
            pair,
            asci,
        )
        self.assertFalse(conflict.available)
        self.assertEqual(1, len(conflict.wt_sites))
        self.assertEqual(2, len(conflict.edited_sites))

    def test_real_insertion_primers(self):
        service = Crispr4pService.from_project_data()
        guide = service.design_gene("ade6").guides[0]
        pair = service.insertion_primers(
            guide.chromosome,
            guide.cut_coordinates,
        )
        cut = guide.cut_coordinates[0]

        self.assertEqual("GCAACTCTGCGATGCATTCA", pair.forward)
        self.assertEqual("TGCGTACTACCATCACTGCA", pair.reverse)
        self.assertAlmostEqual(59.551009575121896, pair.forward_tm)
        self.assertAlmostEqual(59.10793666168439, pair.reverse_tm)
        self.assertEqual(298, pair.wt_product_size)
        self.assertEqual(321, pair.disrupted_product_size)
        self.assertNotEqual(
            -1,
            self.reference.find(pair.forward, cut - 300, cut - 80),
        )
        self.assertNotEqual(
            -1,
            self.reference.find(
                reverse_complement(pair.reverse),
                cut + 80,
                cut + 300,
            ),
        )

    def test_real_junction_primers(self):
        service = Crispr4pService.from_project_data()
        guide = service.design_gene("ade6").guides[0]

        for cassette in self.cassettes:
            checks = service.insertion_checks(
                guide.chromosome,
                guide.cut_coordinates,
                cassette.id,
                "+",
            )
            self.assertEqual(checks.spanning.forward, checks.left.forward)
            self.assertEqual(
                reverse_complement(cassette.sequence),
                checks.left.reverse,
            )
            self.assertEqual(cassette.sequence, checks.right.forward)
            self.assertEqual(checks.spanning.reverse, checks.right.reverse)
            self.assertEqual(193, checks.left.product_size)
            self.assertEqual(151, checks.right.product_size)

        reverse_guide = service.design_gene("bub1").guides[0]
        for cassette in self.cassettes:
            reverse_checks = service.insertion_checks(
                reverse_guide.chromosome,
                reverse_guide.cut_coordinates,
                cassette.id,
                "-",
            )
            self.assertEqual(cassette.sequence, reverse_checks.left.reverse)
            self.assertEqual(
                reverse_complement(cassette.sequence),
                reverse_checks.right.forward,
            )
            self.assertEqual(197, reverse_checks.left.product_size)
            self.assertEqual(123, reverse_checks.right.product_size)
        with self.assertRaisesRegex(ValueError, "coding strand"):
            service.insertion_checks(
                guide.chromosome,
                guide.cut_coordinates,
                1,
                None,
            )

    def test_real_cassette_options(self):
        service = Crispr4pService.from_project_data()
        guide = service.design_gene("ade6").guides[0]
        options = service.cassette_options(
            guide.chromosome,
            guide.cut_coordinates,
            guide.seed,
            1,
            "+",
        )

        self.assertEqual(("none", "asci", "paci", "swai"), tuple(
            option.cassette_format.id for option in options
        ))
        self.assertTrue(all(option.available for option in options))
        self.assertEqual((321, 327, 329, 329), tuple(
            option.spanning.disrupted_product_size for option in options
        ))
        self.assertEqual((151, 157, 159, 159), tuple(
            option.checks.right.product_size for option in options
        ))
        self.assertEqual((193, 134), options[1].digest.fragments)
        self.assertEqual((198, 131), options[2].digest.fragments)
        self.assertEqual((197, 132), options[3].digest.fragments)

        reverse_guide = service.design_gene("bub1").guides[0]
        reverse_options = service.cassette_options(
            reverse_guide.chromosome,
            reverse_guide.cut_coordinates,
            reverse_guide.seed,
            1,
            "-",
        )
        self.assertEqual((197, 203, 205, 205), tuple(
            option.checks.left.product_size for option in reverse_options
        ))
        self.assertEqual((123, 123, 123, 123), tuple(
            option.checks.right.product_size for option in reverse_options
        ))

        paci_conflict = service.cassette_options(
            "I",
            (1799341, 1799342),
            "CAGCTTCATTGAATTAATTA",
            1,
            "+",
        )
        self.assertEqual(
            (True, True, False, True),
            tuple(option.available for option in paci_conflict),
        )
        self.assertEqual(0, len(paci_conflict[2].digest.wt_sites))
        self.assertEqual(2, len(paci_conflict[2].digest.edited_sites))
        self.assertEqual(0, len(paci_conflict[3].digest.wt_sites))
        self.assertEqual(1, len(paci_conflict[3].digest.edited_sites))

        swai_conflict = service.cassette_options(
            "I",
            (251790, 251791),
            "TATATTATTGCATTTAACCT",
            1,
            "-",
        )
        self.assertEqual(
            (True, True, False, False),
            tuple(option.available for option in swai_conflict),
        )
        self.assertEqual(0, len(swai_conflict[3].digest.wt_sites))
        self.assertEqual(2, len(swai_conflict[3].digest.edited_sites))
        self.assertEqual((151, 143, 6), swai_conflict[3].digest.fragments)

    def test_atg7_primer_fallbacks(self):
        service = Crispr4pService.from_project_data()
        result = service.design_gene("atg7")
        guide = result.guides[8]
        annotation = service.annotate_guide(guide)
        choices = service.cassette_choices(
            (guide,),
            (annotation,),
            target_name="atg7",
        )[0]
        options = service.cassette_options(
            guide.chromosome,
            guide.cut_coordinates,
            guide.seed,
            choices[0].id,
            "-",
        )

        self.assertEqual("GGCATTTAACAGTGTACCCT", guide.seed)
        self.assertTrue(all(option.available for option in options))
        self.assertTrue(all(option.spanning is not None for option in options))
        self.assertTrue(all(option.checks is not None for option in options))
        self.assertEqual(328, options[0].spanning.wt_product_size)
        self.assertTrue(all(
            option.digest.available for option in options[1:]
        ))

        wider_pair = service.insertion_primers(
            result.guides[85].chromosome,
            result.guides[85].cut_coordinates,
        )
        self.assertEqual("CAAGCATTAGGTGAAATTGC", result.guides[85].seed)
        self.assertEqual(370, wider_pair.wt_product_size)

    def test_insertion_primer_endpoint(self):
        service = Mock()
        service.insertion_primers.return_value = InsertionPrimerPair(
            forward="FORWARD",
            reverse="REVERSE",
            forward_tm=59.4,
            reverse_tm=60.2,
            wt_product_size=250,
            insert_length=23,
        )
        handler = webapp.CRISPR4PHandler.__new__(webapp.CRISPR4PHandler)
        handler.send_response = Mock()
        handler.send_header = Mock()
        handler.end_headers = Mock()
        handler.wfile = BytesIO()

        with patch.object(webapp, "create_service", return_value=service):
            handler.serve_insertion_primers(
                {
                    "chromosome": ["III"],
                    "cut_left": ["1316791"],
                    "cut_right": ["1316792"],
                }
            )

        service.insertion_primers.assert_called_once_with(
            "III",
            (1316791, 1316792),
            arm_length=80,
            insert_length=23,
            window=300,
        )
        handler.send_response.assert_called_once_with(200)
        self.assertIn(
            ("Cache-Control", "no-store"),
            [call.args for call in handler.send_header.call_args_list],
        )
        self.assertEqual(
            {
                "forward": "FORWARD",
                "reverse": "REVERSE",
                "forward_tm": 59.4,
                "reverse_tm": 60.2,
                "wt_product_size": 250,
                "disrupted_product_size": 273,
            },
            json.loads(handler.wfile.getvalue()),
        )

    def test_junction_primer_endpoint(self):
        spanning = InsertionPrimerPair(
            forward="FORWARD",
            reverse="REVERSE",
            forward_tm=59.4,
            reverse_tm=60.2,
            wt_product_size=298,
            insert_length=23,
        )
        checks = InsertionChecks(
            spanning=spanning,
            left=JunctionPrimerPair(
                "FORWARD",
                "CASSETTE_REVERSE",
                59.4,
                59.1,
                193,
            ),
            right=JunctionPrimerPair(
                "CASSETTE_FORWARD",
                "REVERSE",
                59.1,
                60.2,
                151,
            ),
        )
        service = Mock()
        service.insertion_checks.return_value = checks
        handler = webapp.CRISPR4PHandler.__new__(webapp.CRISPR4PHandler)
        handler.send_response = Mock()
        handler.send_header = Mock()
        handler.end_headers = Mock()
        handler.wfile = BytesIO()

        with patch.object(webapp, "create_service", return_value=service):
            handler.serve_insertion_primers(
                {
                    "chromosome": ["III"],
                    "cut_left": ["1316791"],
                    "cut_right": ["1316792"],
                    "cassette_id": ["1"],
                    "coding_strand": ["+"],
                }
            )

        service.insertion_checks.assert_called_once_with(
            "III",
            (1316791, 1316792),
            1,
            "+",
            arm_length=80,
            window=300,
        )
        self.assertEqual(
            {
                "forward": "FORWARD",
                "reverse": "REVERSE",
                "forward_tm": 59.4,
                "reverse_tm": 60.2,
                "wt_product_size": 298,
                "disrupted_product_size": 321,
                "left_junction": {
                    "forward": "FORWARD",
                    "reverse": "CASSETTE_REVERSE",
                    "forward_tm": 59.4,
                    "reverse_tm": 59.1,
                    "product_size": 193,
                },
                "right_junction": {
                    "forward": "CASSETTE_FORWARD",
                    "reverse": "REVERSE",
                    "forward_tm": 59.1,
                    "reverse_tm": 60.2,
                    "product_size": 151,
                },
            },
            json.loads(handler.wfile.getvalue()),
        )

    def test_cassette_options_endpoint(self):
        pair = InsertionPrimerPair(
            forward="FORWARD",
            reverse="REVERSE",
            forward_tm=59.4,
            reverse_tm=60.2,
            wt_product_size=298,
            insert_length=29,
        )
        asci = cassette_format("asci")
        options = (
            CassetteOption(
                cassette_format=NO_DIAGNOSTIC,
                coding_sequence=FIRST,
                insert=FIRST,
                available=True,
                spanning=pair,
                checks=None,
                digest=None,
            ),
            CassetteOption(
                cassette_format=asci,
                coding_sequence=asci.sequence(self.cassettes[0]),
                insert=asci.sequence(self.cassettes[0]),
                available=True,
                spanning=pair,
                checks=None,
                digest=DigestResult(
                    enzyme="AscI",
                    site="GGCGCGCC",
                    wt_sites=(),
                    edited_sites=(191,),
                    fragments=(193, 134),
                ),
            ),
        )
        service = Mock()
        service.cassette_options.return_value = options
        handler = webapp.CRISPR4PHandler.__new__(webapp.CRISPR4PHandler)
        handler.send_response = Mock()
        handler.send_header = Mock()
        handler.end_headers = Mock()
        handler.wfile = BytesIO()

        with patch.object(webapp, "create_service", return_value=service):
            handler.serve_cassette_options({
                "chromosome": ["III"],
                "cut_left": ["1316791"],
                "cut_right": ["1316792"],
                "guide": [ADE6_GUIDE],
                "cassette_id": ["1"],
                "coding_strand": ["+"],
            })

        service.cassette_options.assert_called_once_with(
            "III",
            (1316791, 1316792),
            ADE6_GUIDE,
            1,
            "+",
            arm_length=80,
            window=300,
        )
        payload = json.loads(handler.wfile.getvalue())
        self.assertEqual(("none", "asci"), tuple(
            row["id"] for row in payload["formats"]
        ))
        self.assertEqual([193, 134], payload["formats"][1]["digest"]["fragments"])
        self.assertEqual(0, payload["formats"][1]["digest"]["wt_site_count"])
        self.assertTrue(payload["formats"][1]["available"])
        self.assertIsNone(payload["formats"][0]["product_url"])
        self.assertIsNone(payload["formats"][0]["protocol_url"])
        self.assertEqual(
            webapp.NEB_LINKS["AscI"],
            (
                payload["formats"][1]["product_url"],
                payload["formats"][1]["protocol_url"],
            ),
        )

    def test_minus_strand_donor(self):
        service = Crispr4pService.from_project_data()
        result = service.design_gene("bub1")
        guide = result.guides[0]
        annotation = service.annotate_guide(guide)
        choices = service.cassette_choices(
            (guide,),
            (annotation,),
            result.name,
        )
        donors = service.disruption_donors(
            (guide,),
            (annotation,),
            choices,
            80,
            result.name,
        )
        donor = donors[0][0]
        cut = guide.cut_coordinates[0]

        self.assertEqual(BUB1_GUIDE, guide.seed)
        self.assertEqual(1, guide.strand)
        self.assertEqual((1314935, 1314937), guide.pam_coordinates)
        self.assertEqual((1314931, 1314932), guide.cut_coordinates)
        self.assertEqual("-", donor.coding_strand)
        self.assertEqual(self.reference[cut - 80:cut], donor.left_arm)
        self.assertEqual(self.reference[cut:cut + 80], donor.right_arm)
        self.assertEqual(reverse_complement(FIRST), donor.insert)
        self.assertEqual(
            reverse_complement(donor.right_arm)
            + FIRST
            + reverse_complement(donor.left_arm),
            donor.reverse,
        )

    def test_target_strand(self):
        plus = gene("plus", "target", "+")
        minus = gene("minus", "other", "-")
        non_coding = gene("rna", None, "-", "lncRNA_gene")

        self.assertEqual(
            "+",
            target_strand(SimpleNamespace(genes=(plus, minus)), "target"),
        )
        self.assertEqual(
            "-",
            target_strand(SimpleNamespace(genes=(minus,))),
        )
        self.assertIsNone(
            target_strand(SimpleNamespace(genes=(plus, minus)))
        )
        annotation = SimpleNamespace(genes=(plus, non_coding))
        self.assertEqual(non_coding, target_gene(annotation, "rna"))
        self.assertIsNone(target_strand(annotation, "rna"))

    def test_service_filters_recut_sites(self):
        service = Crispr4pService.from_project_data(
            cassettes=(unsafe_cassette(),) + self.cassettes,
        )
        result = service.design_gene("ade6")
        guide = result.guides[0]
        annotation = service.annotate_guide(guide)

        choices = service.cassette_choices(
            (guide,),
            (annotation,),
            "ade6",
        )
        donors = service.disruption_donors(
            (guide,),
            (annotation,),
            choices,
            80,
            "ade6",
        )

        self.assertEqual((self.cassettes,), choices)
        self.assertEqual(10, len(donors[0]))
        self.assertEqual(183, donors[0][0].total_length)

    def test_format_filters_extra_rescue_sites(self):
        service = Crispr4pService.from_project_data()
        result = service.design_gene("ade6")
        guide = result.guides[0]

        with patch(
            "crispr4p.service.valid_rescue_site",
            return_value=False,
        ) as rescue_check:
            options = service.cassette_options(
                guide.chromosome,
                guide.cut_coordinates,
                guide.seed,
                1,
                "+",
            )

        self.assertTrue(all(not option.available for option in options))
        self.assertEqual(len(CASSETTE_FORMATS), rescue_check.call_count)

    def test_noncoding_target(self):
        service = Crispr4pService.from_project_data()
        result = service.design_gene("SPNCRNA.1438")
        annotations = service.annotate_guides(result.guides)
        choices = service.cassette_choices(
            result.guides,
            annotations,
            result.name,
        )
        donors = service.disruption_donors(
            result.guides,
            annotations,
            choices,
            80,
            result.name,
        )
        rows = annotation_rows(result.guides, annotations, result.name)
        template = (PROJECT_ROOT / "template/container_table.html").read_text(
            encoding="utf-8"
        )
        page = render_design(
            result,
            annotations,
            template,
            cassette_choices=choices,
            disruption_donors=donors,
        )

        self.assertTrue(all(not group for group in choices))
        self.assertTrue(all(not group for group in donors))
        self.assertTrue(all(row["coding_target"] is False for row in rows))
        self.assertTrue(all(row["coding_strand"] is None for row in rows))
        self.assertIn("Selected target is non-coding.", page)

    def test_compact_web_data(self):
        service = Crispr4pService.from_project_data()
        result = service.design_gene("ade6")
        annotations = service.annotate_guides(result.guides)
        choices = service.cassette_choices(
            result.guides,
            annotations,
            result.name,
        )
        donors = service.disruption_donors(
            result.guides,
            annotations,
            choices,
            80,
            result.name,
        )
        template = (PROJECT_ROOT / "template/container_table.html").read_text(
            encoding="utf-8"
        )
        page = render_design(
            result,
            annotations,
            template,
            cassette_choices=choices,
            disruption_donors=donors,
        )

        def browser_data(name):
            start = page.index(f"var {name} = ") + len(f"var {name} = ")
            end = page.index(";", start)
            return page[start:end], json.loads(page[start:end])

        cassette_json, cassettes = browser_data("cassette_data")
        donor_json, arms = browser_data("donor_arms")
        cassette = cassettes["catalog"]["1"]
        donor = donors[0][0]
        insert = cassette["sequence"]
        if arms[0]["coding_strand"] == "-":
            insert = reverse_complement(insert)
        sequence = arms[0]["left_arm"] + insert + arms[0]["right_arm"]
        overlap_start = len(sequence) - arms[0]["hrrv_length"]
        browser_forward = sequence[:arms[0]["hrfw_length"]]
        browser_reverse = reverse_complement(sequence[overlap_start:])
        browser_overlap = sequence[
            overlap_start:arms[0]["hrfw_length"]
        ]

        self.assertEqual(donor.sequence, sequence)
        self.assertEqual(donor.reverse, reverse_complement(sequence))
        self.assertEqual(donor.oligos.forward, browser_forward)
        self.assertEqual(donor.oligos.reverse, browser_reverse)
        self.assertEqual(donor.oligos.overlap, browser_overlap)
        self.assertEqual(103, arms[0]["hrfw_length"])
        self.assertEqual(103, arms[0]["hrrv_length"])
        self.assertEqual(23, arms[0]["overlap_length"])
        self.assertEqual(183, arms[0]["hr_product_length"])
        self.assertEqual(10, len(cassettes["catalog"]))
        self.assertEqual(144, len(cassettes["choices"]))
        self.assertEqual(144, len(arms))
        self.assertLess(len(cassette_json.encode("utf-8")), 10000)
        self.assertLess(len(donor_json.encode("utf-8")), 50000)

    def test_web_section(self):
        annotations = GenomeAnnotations.from_files(
            DATA / "Schizosaccharomyces_pombe_all_chromosomes.gff3",
            DATA / "gene_viability.tsv",
        )
        annotation = annotations.annotate_cut("III", (1316791, 1316792))
        guide = ADE6_GUIDE
        primer = (
            guide,
            "FORWARD",
            "REVERSE",
            (1316795, 1316797),
            1,
            "TGG",
        )
        result = DesignResult.from_legacy(
            (
                [[guide, primer, 5, 1, 1, 1, 1, 1, 1]],
                ("HR_FORWARD", "HR_REVERSE", "DELETED_DNA"),
                [],
                "ade6",
                "III",
                "1316337",
                "1317995",
            )
        )
        template = (PROJECT_ROOT / "template/container_table.html").read_text(
            encoding="utf-8"
        )

        cassettes = cassette_data((self.cassettes,))
        donors = tuple(
            build_donor(
                self.reference,
                (1316791, 1316792),
                cassette,
                "+",
                80,
            )
            for cassette in self.cassettes
        )
        arms = donor_data((donors,))
        guide_rows = annotation_rows(result.guides, (annotation,), "ade6")
        page = render_design(
            result,
            (annotation,),
            template,
            cassette_choices=(self.cassettes,),
            disruption_donors=(donors,),
        )

        self.assertEqual(10, len(cassettes["catalog"]))
        self.assertEqual(list(range(1, 11)), cassettes["choices"][0])
        self.assertEqual(
            "TGA* TGA* GCT GAT GAC CTA GTA",
            cassettes["catalog"]["1"]["frames"][0],
        )
        self.assertEqual(80, arms[0]["arm_length"])
        self.assertEqual(donors[0].left_arm, arms[0]["left_arm"])
        self.assertEqual(donors[0].right_arm, arms[0]["right_arm"])
        self.assertNotIn("sequence", arms[0])
        self.assertEqual("+", guide_rows[0]["coding_strand"])
        self.assertIn("Stop-Cassette Disruption Design", page)
        self.assertIn("Wild-Type Restoration Design", page)
        self.assertIn(FIRST, page)
        self.assertNotIn(OLD_FIRST, page)
        self.assertIn('id="stop_cassette_menu"', page)
        self.assertIn('id="cassette_format_menu"', page)
        report_start = page.index('id="design_report"')
        self.assertGreater(
            report_start,
            page.index('id="wild_type_restoration"'),
        )
        self.assertIn(
            '<textarea id="report_text" readonly '
            'aria-label="Copyable design report"></textarea>',
            page,
        )
        self.assertIn('onclick="copy_report();">Copy report</button>', page)


if __name__ == "__main__":
    unittest.main()
