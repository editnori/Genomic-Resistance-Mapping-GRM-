"""scripts/c4-screen/cmh_plain.py: the plain implementation the Rust adjusted screen is held to
(docs/screen-arms.md).

Every table here is SYNTHETIC: made-up labels, strata and rows from fixed seeds.
"""

import contextlib
import io
import json
import tempfile
import unittest
from fractions import Fraction
from pathlib import Path

from tests.c4_screen import common

plain = common.load("cmh_plain")
CASES = common.ROOT / "tests" / "c4_screen" / "cmh_cases.json"


class Plain(unittest.TestCase):
    def test_the_committed_cases_are_what_the_script_writes(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "cases.json"
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                self.assertEqual(plain.main([str(out)]), 0)
            self.assertEqual(out.read_text(), CASES.read_text())

    def test_with_one_stratum_the_statistic_is_the_chi_square_times_n_less_one_over_n(self):
        case = json.loads(CASES.read_text())["cases"][0]
        self.assertEqual(case["name"], "one-stratum")
        labels = [int(c) for c in case["labels"]]
        fit = [int(c) for c in case["fit"]]
        n = sum(fit)
        for row in case["rows"]:
            row = [int(c) for c in row]
            self.assertEqual(plain.statistic(row, labels, fit, case["strata"]),
                             plain.chi_square(row, labels, fit) * Fraction(n - 1, n))

    def test_a_stratum_of_one_class_or_one_sample_adds_nothing(self):
        labels = [1, 0, 1, 0, 1, 0, 1, 1, 0, 0]
        lineage = [0, 0, 0, 0, 1, 1, 1, 1, 1, 1]
        fit = [1] * 10
        row = [1, 0, 1, 1, 0, 1, 1, 0, 0, 1]
        base = plain.statistic(row, labels, fit, lineage)
        # Two resistant genomes in a stratum of their own, and one genome alone.
        more_labels = labels + [1, 1, 0]
        more_lineage = lineage + [7, 7, 9]
        self.assertEqual(plain.statistic(row + [1, 0, 1], more_labels, fit + [1, 1, 1], more_lineage), base)
        nr, slot, all_strata = plain.strata(more_labels, fit + [1, 1, 1], more_lineage)
        self.assertEqual((len(nr), all_strata, slot[-3:]), (2, 4, [None, None, None]))

    def test_a_genome_outside_the_fit_is_not_read(self):
        labels = [1, 0, 1, 0, 1, 0, 1, 1, 0, 0]
        lineage = [0, 0, 0, 0, 1, 1, 1, 1, 1, 1]
        fit = [1, 1, 1, 1, 1, 1, 1, 1, 0, 0]
        row = [1, 0, 1, 1, 0, 1, 1, 0, 0, 1]
        base = plain.statistic(row, labels, fit, lineage)
        self.assertEqual(plain.statistic(row[:8] + [1, 0], labels[:8] + [1, 1], fit, lineage[:8] + [5, 5]), base)

    def test_the_top_breaks_ties_by_the_lower_row(self):
        self.assertEqual(plain.top([Fraction(1), Fraction(3), Fraction(1), Fraction(3), Fraction(0)], 4), [1, 3, 0, 2])

    def test_every_case_is_synthetic_and_its_counts_follow_from_its_tables(self):
        doc = json.loads(CASES.read_text())
        self.assertEqual(doc["format"], plain.FORMAT)
        self.assertIn("SYNTHETIC", doc["note"])
        for case in doc["cases"]:
            labels = [int(c) for c in case["labels"]]
            fit = [int(c) for c in case["fit"]]
            nr, _, all_strata = plain.strata(labels, fit, case["strata"])
            self.assertEqual((len(nr), all_strata), (case["informative"], case["strata_in_fit"]), case["name"])
            self.assertEqual(len(case["top"]), min(case["k"], len(case["rows"])), case["name"])


if __name__ == "__main__":
    unittest.main()
