import unittest
from unittest.mock import Mock, patch

import crispr4p.crispr4p as legacy_core
from crispr4p.crispr4p import PrimerDesign, chromosomeFasta
from crispr4p.primers import build_hr_dna, checking_primers


SEQUENCE = "A" * 300 + "G" * 200 + "C" * 300
START = 300
END = 500
PUBLIC_START = START + 1
PUBLIC_END = END
WIDTH = 300


def primer3_answer():
    answer = {}
    for index in range(2):
        answer.update(
            {
                f"PRIMER_LEFT_{index}_SEQUENCE": f"L{index}",
                f"PRIMER_RIGHT_{index}_SEQUENCE": f"R{index}",
                f"PRIMER_LEFT_{index}_TM": 60.0 + index,
                f"PRIMER_RIGHT_{index}_TM": 61.0 + index,
                f"PRIMER_LEFT_{index}_GC_PERCENT": 40.0 + index,
                f"PRIMER_RIGHT_{index}_GC_PERCENT": 41.0 + index,
                f"PRIMER_PAIR_{index}_PRODUCT_SIZE": 250 + index,
            }
        )
    return answer


class TestPrimerComputations(unittest.TestCase):
    def test_builds_legacy_hr_sequences(self):
        result = build_hr_dna(
            SEQUENCE,
            START,
            END,
            PrimerDesign.sequenceComplement_,
        )

        self.assertEqual(
            (
                "A" * 80 + "C" * 20,
                "G" * 80 + "T" * 20,
                "A" * 250 + "C" * 250,
            ),
            result,
        )

    def test_preserves_exact_primer3_inputs_and_result_shape(self):
        primer_designer = Mock(return_value=primer3_answer())

        result = checking_primers(
            SEQUENCE,
            START,
            END,
            WIDTH,
            number_of_alternatives=2,
            primer_designer=primer_designer,
        )

        sequence_args, global_args = primer_designer.call_args.args
        self.assertEqual(
            {
                "SEQUENCE_ID": "MH1000",
                "SEQUENCE_TEMPLATE": "A" * 300 + "C" * 300,
                "SEQUENCE_INCLUDED_REGION": [0, 600],
                "SEQUENCE_EXCLUDED_REGION": [[220, 160]],
            },
            sequence_args,
        )
        self.assertEqual(
            {
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
                "PRIMER_PRODUCT_SIZE_RANGE": [[225, 600]],
            },
            global_args,
        )
        self.assertEqual(
            [
                {
                    "PRIMER_LEFT_0_SEQUENCE": "L0",
                    "PRIMER_RIGHT_0_SEQUENCE": "R0",
                    "PRIMER_LEFT_0_TM": 60.0,
                    "PRIMER_RIGHT_0_TM": 61.0,
                    "PRIMER_LEFT_0_GC_PERCENT": 40.0,
                    "PRIMER_RIGHT_0_GC_PERCENT": 41.0,
                    "PRIMER_PAIR_0_PRODUCT_SIZE": 250,
                    "negative_result": 450,
                },
                {
                    "PRIMER_LEFT_1_SEQUENCE": "L1",
                    "PRIMER_RIGHT_1_SEQUENCE": "R1",
                    "PRIMER_LEFT_1_TM": 61.0,
                    "PRIMER_RIGHT_1_TM": 62.0,
                    "PRIMER_LEFT_1_GC_PERCENT": 41.0,
                    "PRIMER_RIGHT_1_GC_PERCENT": 42.0,
                    "PRIMER_PAIR_1_PRODUCT_SIZE": 251,
                    "negative_result": 451,
                },
            ],
            result,
        )

    def test_uses_only_returned_checking_pairs(self):
        one_pair = primer3_answer()
        one_pair["PRIMER_PAIR_NUM_RETURNED"] = 1
        no_pairs = {"PRIMER_PAIR_NUM_RETURNED": 0}

        first = checking_primers(
            SEQUENCE,
            START,
            END,
            WIDTH,
            number_of_alternatives=2,
            primer_designer=Mock(return_value=one_pair),
        )
        empty = checking_primers(
            SEQUENCE,
            START,
            END,
            WIDTH,
            number_of_alternatives=2,
            primer_designer=Mock(return_value=no_pairs),
        )

        self.assertEqual(1, len(first))
        self.assertEqual("L0", first[0]["PRIMER_LEFT_0_SEQUENCE"])
        self.assertEqual([], empty)

    def test_primer_design_methods_remain_legacy_adapters(self):
        chromosome = chromosomeFasta("synthetic description\n" + SEQUENCE)
        designer = PrimerDesign.__new__(PrimerDesign)
        designer._numAlternativeCheckings = 2
        hr_sentinel = object()
        primer_sentinel = object()

        with patch(
            "crispr4p.crispr4p.build_hr_dna",
            return_value=hr_sentinel,
        ) as build_hr:
            hr_result = designer.HR_DNA(
                chromosome,
                PUBLIC_START,
                PUBLIC_END,
            )

        self.assertIs(hr_sentinel, hr_result)
        build_hr.assert_called_once_with(
            chromosome.sequence,
            START,
            END,
            designer.sequenceComplement_,
        )

        with patch(
            "crispr4p.crispr4p.checking_primers",
            return_value=primer_sentinel,
        ) as design_checking:
            primer_result = designer.CheckingPrimersWidth_(
                chromosome,
                PUBLIC_START,
                PUBLIC_END,
                WIDTH,
            )

        self.assertIs(primer_sentinel, primer_result)
        design_checking.assert_called_once_with(
            chromosome.sequence,
            START,
            END,
            WIDTH,
            designer._numAlternativeCheckings,
            primer_designer=legacy_core._design_primers,
        )

    def test_default_checking_width_remains_300(self):
        chromosome = chromosomeFasta("synthetic description\n" + SEQUENCE)
        designer = PrimerDesign.__new__(PrimerDesign)
        sentinel = object()

        with patch.object(
            designer,
            "CheckingPrimersWidth_",
            return_value=sentinel,
        ) as design_width:
            result = designer.CheckingPrimers(
                chromosome,
                PUBLIC_START,
                PUBLIC_END,
            )

        self.assertIs(sentinel, result)
        design_width.assert_called_once_with(
            chromosome,
            PUBLIC_START,
            PUBLIC_END,
            300,
        )


if __name__ == "__main__":
    unittest.main()
