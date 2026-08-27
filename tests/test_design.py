import unittest
from unittest.mock import Mock, patch

import crispr4p.design as design
from crispr4p.crispr4p import NGG, PrimerDesign, chromosomeFasta


LENGTHS = (8, 10, 12, 14, 16, 18, 20)


def match_table(counts):
    return {
        length: [object()] * count
        for length, count in zip(LENGTHS, counts)
    }


class TestDesignTableConstruction(unittest.TestCase):
    def test_builds_rows_and_preserves_reverse_priority_sort(self):
        first = NGG("I", 1, 1, "A" * 20, "TGG")
        second = NGG("I", 2, 1, "C" * 20, "TGG")
        tied = NGG("I", 3, 1, "G" * 20, "TGG")
        first.primer = ("first",)
        second.primer = ("second",)
        tied.primer = ("tied",)
        guide_matches = {
            first: match_table((1, 1, 1, 1, 1, 2, 2)),
            second: match_table((5, 1, 1, 1, 1, 1, 1)),
            tied: match_table((5, 1, 1, 1, 1, 1, 1)),
        }

        table = design.guide_table(guide_matches)

        self.assertEqual(
            [
                ["C" * 20, ("second",), 5, 1, 1, 1, 1, 1, 1],
                ["G" * 20, ("tied",), 5, 1, 1, 1, 1, 1, 1],
                ["A" * 20, ("first",), 1, 1, 1, 1, 1, 2, 2],
            ],
            table,
        )

class TestDesignQueryOrchestration(unittest.TestCase):
    def test_runs_callbacks_in_order_and_returns_same_match_mapping(self):
        events = []
        chromosome = chromosomeFasta("I description\n" + "A" * 100)
        guide = NGG("I", 1, 1, "A" * 20, "TGG")
        guide_matches = {}

        def ensure_index():
            events.append("index")

        def discover(record, start, end):
            events.append(("discover", record, start, end))

        def match(n_mismatch):
            events.append(("match", n_mismatch))
            guide_matches[guide] = match_table((1,) * 7)

        def build_primer(record, start, end, matched_guide):
            events.append(("guide_primer", matched_guide))
            return ("primer",)

        def build_hr(record, start, end):
            events.append("hr")
            return ("hr",)

        def build_checking(record, start, end):
            events.append("checking")
            return ["checking"]

        result = design.run_design(
            "I",
            10,
            90,
            2,
            regression=False,
            chromosomes={"I": chromosome},
            guide_matches=guide_matches,
            ensure_genome_index=ensure_index,
            discover_guides=discover,
            match_guides=match,
            build_guide_primer=build_primer,
            build_hr_dna=build_hr,
            build_checking_primers=build_checking,
        )

        self.assertEqual(
            [
                "index",
                ("discover", chromosome, 10, 90),
                ("match", 2),
                ("guide_primer", guide),
                "hr",
                "checking",
            ],
            events,
        )
        self.assertEqual(("primer",), guide.primer)
        self.assertIs(guide_matches, result[3])
        self.assertEqual(("hr",), result[1])
        self.assertEqual(["checking"], result[2])

    def test_regression_mode_skips_index_build_and_keeps_accumulated_guides(self):
        chromosome = chromosomeFasta("I description\n" + "A" * 100)
        old_guide = NGG("I", 1, 1, "A" * 20, "TGG")
        new_guide = NGG("I", 2, 1, "C" * 20, "TGG")
        guide_matches = {old_guide: match_table((1,) * 7)}
        ensure_index = Mock()

        def add_new_match(n_mismatch):
            guide_matches[new_guide] = match_table((2,) * 7)

        primer_builder = Mock(side_effect=lambda *args: (args[3].seed,))
        result = design.run_design(
            "I",
            10,
            90,
            0,
            regression=True,
            chromosomes={"I": chromosome},
            guide_matches=guide_matches,
            ensure_genome_index=ensure_index,
            discover_guides=Mock(),
            match_guides=add_new_match,
            build_guide_primer=primer_builder,
            build_hr_dna=Mock(return_value=("hr",)),
            build_checking_primers=Mock(return_value=["checking"]),
        )

        ensure_index.assert_not_called()
        self.assertEqual(2, primer_builder.call_count)
        self.assertEqual((old_guide.seed,), old_guide.primer)
        self.assertEqual((new_guide.seed,), new_guide.primer)
        self.assertIs(guide_matches, result[3])
        self.assertEqual(2, len(result[0]))

    def test_primer_design_run_method_remains_thin_adapter(self):
        designer = PrimerDesign.__new__(PrimerDesign)
        designer.regression = False
        designer.chromosomesData = {"I": object()}
        designer.tableNGGs = {}
        designer.getNGGsFromGenome = Mock()
        designer._getUserNGGs = Mock()
        designer.gRNA_Table = Mock()
        designer.getPrimerGRNA = Mock()
        designer.HR_DNA = Mock()
        designer.CheckingPrimers = Mock()
        sentinel = object()

        with patch(
            "crispr4p.crispr4p.run_design",
            return_value=sentinel,
        ) as run_query:
            result = designer.run_("I", 10, 90, 2, "name")

        self.assertIs(sentinel, result)
        run_query.assert_called_once_with(
            "I",
            10,
            90,
            2,
            regression=designer.regression,
            chromosomes=designer.chromosomesData,
            guide_matches=designer.tableNGGs,
            ensure_genome_index=designer.getNGGsFromGenome,
            discover_guides=designer._getUserNGGs,
            match_guides=designer.gRNA_Table,
            build_guide_primer=designer.getPrimerGRNA,
            build_hr_dna=designer.HR_DNA,
            build_checking_primers=designer.CheckingPrimers,
        )


if __name__ == "__main__":
    unittest.main()
