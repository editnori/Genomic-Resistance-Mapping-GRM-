"""The parts of the learner cross-check kit that need no reference library, on cases worked by hand.

panels.py decodes the panel file, scm_reference.py is the independent set covering machine, exact.py gives a tree's
split gains as fractions, check_env.py reads the pinned requirement files and table.py finds what is wrong with a
result table. compare.py and the three check_*.py files need scikit-learn, XGBoost and pyscm-ml and are not run here."""

import contextlib
import copy
import io
import json
import pathlib
import shutil
import tempfile
import unittest
from fractions import Fraction

from . import kit

table, panels = kit.table, kit.panels
scm, exact, check_env = kit.module("scm_reference"), kit.module("exact"), kit.module("check_env")


def bits(*samples):
    return sum(1 << sample for sample in samples)


class Panels(unittest.TestCase):
    def test_a_set_is_hex_nibbles_low_bit_first(self):
        # nibble k holds samples 4k to 4k + 3, and the low bit of a nibble is its first sample
        self.assertEqual(panels.decode_bits(8, "10"), bits(0))
        self.assertEqual(panels.decode_bits(8, "08"), bits(7))
        self.assertEqual(panels.decode_bits(8, "a3"), bits(1, 3, 4, 5))
        self.assertEqual(panels.decode_bits(5, "f1"), bits(0, 1, 2, 3, 4))
        self.assertEqual(panels.decode_bits(0, ""), 0)

    def test_encode_is_the_inverse_of_decode(self):
        for n, text in ((8, "a3"), (5, "f1"), (13, "0c21"), (4, "0"), (0, "")):
            self.assertEqual(panels.encode_bits(n, panels.decode_bits(n, text)), text)

    def test_a_set_of_the_wrong_length_or_past_the_last_sample_is_refused(self):
        with self.assertRaisesRegex(ValueError, "has 1 nibbles"):
            panels.decode_bits(8, "a")
        with self.assertRaisesRegex(ValueError, "has 3 nibbles"):
            panels.decode_bits(8, "a30")
        with self.assertRaisesRegex(ValueError, "past sample 5"):
            panels.decode_bits(5, "f2")

    def test_members_and_count(self):
        self.assertEqual(panels.members(bits(9, 0, 4)), [0, 4, 9])
        self.assertEqual(panels.members(0), [])
        self.assertEqual(panels.count(bits(9, 0, 4)), 3)

    def read(self, text):
        with tempfile.TemporaryDirectory() as folder:
            path = pathlib.Path(folder) / "panels.txt"
            path.write_text(text, encoding="ascii")
            return panels.read(path)

    def test_a_small_file(self):
        text = "# a comment\nfixture one 6 2 3 400\ny 30\nfit f3\nr 50\nr a2\npoint 0 0.5 0.25 3 1 0 0.0 1:0.5\n\nfixture two 4 1 3 400\ny 1\nfit f\nr 8\n"
        one, two = self.read(text)
        self.assertEqual(one, {"name": "one", "n": 6, "y": bits(0, 1), "fit": bits(0, 1, 2, 3, 4, 5), "rows": [bits(0, 2), bits(1, 3, 5)]})
        self.assertEqual(two, {"name": "two", "n": 4, "y": bits(0), "fit": bits(0, 1, 2, 3), "rows": [bits(3)]})
        self.assertEqual(
            panels.counts(one),
            {"samples": 6, "features": 2, "resistant": 2, "fit": 6, "fit_resistant": 2, "present": [2, 3]},
        )

    def test_a_file_that_is_not_a_panel_file_is_refused(self):
        good = "fixture one 6 2 3 400\ny 30\nfit f3\nr 50\nr a2\n"
        for text, message in (
            (good.replace("r a2\n", ""), "has 1 feature rows, not the count it declares"),
            (good + "y 30\n", "a second y line"),
            (good.replace("fit f3\n", ""), "lacks y or fit"),
            (good.replace("r a2", "row a2"), "unknown line type 'row'"),
            (good.replace("r a2", "r a2 0"), "a set line has 2 fields"),
            (good.replace("r a2", "r a4"), "line 5: a bit is set past sample 6"),
            (good.replace("fixture one 6 2 3 400", "fixture one 6 2 3"), "a fixture line has 6 fields"),
            ("y 30\n" + good, "no fixture yet"),
            (good + good, "a panel name is used twice"),
        ):
            with self.subTest(message=message), self.assertRaisesRegex(ValueError, message):
                self.read(text)


class SetCoveringMachine(unittest.TestCase):
    """Samples 0.. are bits; `y` is the resistant set and every sample is in training unless a case says otherwise."""

    def test_one_feature_that_marks_the_resistant_samples(self):
        # presence of feature 0 holds on both resistant samples and rules out both susceptible ones
        rules, steps = scm.fit([bits(0, 1)], bits(0, 1), bits(0, 1, 2, 3), 4, "conjunction", 1.0, 6)
        self.assertEqual((rules, steps), ([(0, True)], [(2, 0, 2)]))
        self.assertEqual(scm.predict([bits(0, 1)], 4, "conjunction", rules), bits(0, 1))

    def test_a_disjunction_swaps_the_classes_and_negates_the_rule(self):
        # the search rules out the resistant samples with "feature 0 absent"; the model stores "feature 0 present"
        rules, steps = scm.fit([bits(0, 1)], bits(0, 1), bits(0, 1, 2, 3), 4, "disjunction", 1.0, 6)
        self.assertEqual((rules, steps), ([(0, True)], [(2, 0, 2)]))
        self.assertEqual(scm.predict([bits(0, 1)], 4, "disjunction", rules), bits(0, 1))
        self.assertEqual(scm.predict([bits(0, 1), bits(3)], 4, "disjunction", [(0, True), (1, True)]), bits(0, 1, 3))
        self.assertEqual(scm.predict([bits(0, 1), bits(1, 3)], 4, "conjunction", [(0, True), (1, False)]), bits(0))

    def test_equal_rules_go_to_the_lower_feature(self):
        y, train = bits(0, 1), bits(0, 1, 2, 3)
        # two copies of a feature
        self.assertEqual(scm.fit([bits(0, 1), bits(0, 1)], y, train, 4, "conjunction", 1.0, 6)[0], [(0, True)])
        # a feature and its complement: presence of one is absence of the other, and the lower feature is taken
        self.assertEqual(scm.fit([bits(0, 1), bits(2, 3)], y, train, 4, "conjunction", 1.0, 6)[0], [(0, True)])
        self.assertEqual(scm.fit([bits(2, 3), bits(0, 1)], y, train, 4, "conjunction", 1.0, 6)[0], [(0, False)])

    def test_presence_is_taken_before_absence(self):
        # feature 0 is on samples 0 and 2: each of its two rules rules out one susceptible and one resistant sample
        rules, steps = scm.fit([bits(0, 2)], bits(0, 1), bits(0, 1, 2, 3), 4, "conjunction", 0.5, 6)
        self.assertEqual(rules, [(0, True), (0, False)])
        self.assertEqual(steps, [(1, 1, Fraction(1, 2)), (1, 1, Fraction(1, 2))])

    def test_fewer_lost_comes_before_the_lower_feature_and_a_rule_that_covers_nothing_stops_the_fit(self):
        # Resistant 0 1 2, susceptible 3 4 5, p = 1. Presence of feature 0 (on 0 1 3) rules out 4 and 5 and loses 2:
        # utility 1. Presence of feature 1 (on 0 1 2 3 4) rules out 5 and loses none: utility 1, and it is taken.
        # At the second step the largest utility is 0, held by feature 0 (covers 4, loses 2) and by feature 1 again
        # (covers nothing, loses nothing). Fewer lost wins, that rule covers nothing, and the fit stops.
        rows, y, train = [bits(0, 1, 3), bits(0, 1, 2, 3, 4)], bits(0, 1, 2), bits(0, 1, 2, 3, 4, 5)
        self.assertEqual(scm.fit(rows, y, train, 6, "conjunction", 1.0, 6), ([(1, True)], [(1, 0, 1)]))
        # at a smaller penalty feature 0's rule has the larger utility both times
        rules, steps = scm.fit(rows, y, train, 6, "conjunction", 0.5, 6)
        self.assertEqual(rules[0], (0, True))
        self.assertEqual(steps[0], (2, 1, Fraction(3, 2)))

    def test_the_fit_stops_at_max_rules(self):
        # three susceptible samples, each ruled out by its own feature
        rows, y, train = [bits(0, 2, 3), bits(0, 1, 3), bits(0, 1, 2)], bits(0), bits(0, 1, 2, 3)
        self.assertEqual(scm.fit(rows, y, train, 4, "conjunction", 1.0, 6)[0], [(0, True), (1, True), (2, True)])
        self.assertEqual(scm.fit(rows, y, train, 4, "conjunction", 1.0, 2)[0], [(0, True), (1, True)])
        self.assertEqual(scm.fit(rows, y, train, 4, "conjunction", 1.0, 0), ([], []))

    def test_a_rule_with_a_negative_utility_is_still_added_while_it_covers_a_sample(self):
        # One feature, on resistant 0 and 1 and susceptible 3. Its absence rules out 3 and loses 0 and 1: 1 - 2p.
        # Its presence rules out 4 and loses 2: 1 - p, the larger utility, negative at p = 3.
        rules, steps = scm.fit([bits(0, 1, 3)], bits(0, 1, 2), bits(0, 1, 2, 3, 4), 5, "conjunction", 3.0, 1)
        self.assertEqual((rules, steps), ([(0, True)], [(1, 1, -2)]))

    def test_one_class_in_training(self):
        rows, everyone = [bits(0, 1)], bits(0, 1, 2, 3)
        # nothing to rule out: no rule, and a conjunction of no rules calls every sample resistant
        self.assertEqual(scm.fit(rows, everyone, everyone, 4, "conjunction", 1.0, 6), ([], []))
        self.assertEqual(scm.predict(rows, 4, "conjunction", []), everyone)
        self.assertEqual(scm.fit(rows, 0, everyone, 4, "disjunction", 1.0, 6), ([], []))
        self.assertEqual(scm.predict(rows, 4, "disjunction", []), 0)

    def test_only_training_samples_count(self):
        # sample 3 is susceptible and has the feature, and it is outside training
        rules, steps = scm.fit([bits(0, 1, 3)], bits(0, 1), bits(0, 1, 2), 4, "conjunction", 1.0, 6)
        self.assertEqual((rules, steps), ([(0, True)], [(1, 0, 1)]))

    def test_the_utility_is_exact(self):
        # 0.1 is no exact binary fraction; the utility is taken with the float's exact value
        rules, steps = scm.fit([bits(0, 2)], bits(0, 1), bits(0, 1, 2, 3), 4, "conjunction", 0.1, 1)
        self.assertEqual(steps, [(1, 1, 1 - Fraction(0.1))])
        self.assertNotEqual(steps[0][2], Fraction(9, 10))

    def test_an_unknown_kind_is_refused(self):
        with self.assertRaisesRegex(ValueError, "unknown kind 'xor'"):
            scm.fit([bits(0)], bits(0), bits(0, 1), 2, "xor", 1.0, 6)


class ExactGains(unittest.TestCase):
    def test_balanced_weights(self):
        # 1 resistant and 3 susceptible training samples: a resistant sample weighs 4 / 2, a susceptible one 4 / 6
        weights = exact.Weights(1, 3)
        self.assertEqual((weights.resistant, weights.susceptible), (2, Fraction(2, 3)))
        # the root then has weighted counts 2 and 2: impurity times weight 2 * 2 * 2 / 4
        self.assertEqual(weights.mass(1, 3), 2)
        self.assertEqual(weights.value(1, 3), Fraction(1, 2))
        self.assertEqual(weights.value(1, 0), 1)
        self.assertEqual(weights.value(0, 2), 0)
        self.assertEqual(weights.value(1, 1), Fraction(3, 4))
        self.assertEqual((weights.mass(1, 0), weights.mass(0, 3), weights.mass(0, 0)), (0, 0, 0))
        self.assertEqual(weights.value(0, 0), Fraction(1, 2))

    def test_a_class_with_no_sample_counts_as_one(self):
        weights = exact.Weights(0, 4)
        self.assertEqual((weights.resistant, weights.susceptible), (2, Fraction(1, 2)))
        self.assertEqual(weights.value(0, 4), 0)

    def test_the_gain_of_a_split(self):
        weights = exact.Weights(2, 2)
        resistant, susceptible = bits(0, 1), bits(2, 3)
        # a feature on the two resistant samples leaves two pure children: the whole impurity, 2, is the gain
        self.assertEqual(exact.gain(weights, resistant, susceptible, bits(0, 1)), 2)
        # a feature on one sample of each class leaves two children as mixed as the node
        self.assertEqual(exact.gain(weights, resistant, susceptible, bits(0, 2)), 0)
        # on samples 0 1 2: that child has mass 2 * 2 * 1 / 3, the other is pure
        self.assertEqual(exact.gain(weights, resistant, susceptible, bits(0, 1, 2)), Fraction(2, 3))
        # a feature on every sample of the node, or on none, splits nothing
        self.assertIsNone(exact.gain(weights, resistant, susceptible, bits(0, 1, 2, 3)))
        self.assertIsNone(exact.gain(weights, resistant, susceptible, bits(5)))

    def test_the_best_set_holds_every_feature_of_the_largest_gain(self):
        weights = exact.Weights(2, 2)
        resistant, susceptible = bits(0, 1), bits(2, 3)
        rows = [bits(0, 1, 2), bits(2, 3), bits(0, 1), bits(0, 2), bits(0, 1, 2, 3)]
        # features 1 and 2 are present on opposite samples and split the node alike
        self.assertEqual(exact.best(weights, resistant, susceptible, rows, range(5)), (2, [1, 2]))
        self.assertEqual(exact.best(weights, resistant, susceptible, rows, [3, 0, 0]), (Fraction(2, 3), [0]))
        self.assertEqual(exact.best(weights, resistant, susceptible, rows, [3]), (0, [3]))
        self.assertEqual(exact.best(weights, resistant, susceptible, rows, [4]), (None, []))
        self.assertEqual(exact.best(weights, resistant, susceptible, rows, []), (None, []))

    def test_fractions_keep_a_tie_that_floats_can_break(self):
        # 3 resistant and 7 susceptible samples; features 0 and 1 are present on opposite samples, so they split the
        # node into the same two groups and their gains are equal. Computed in floats, one side directly and the other
        # by subtraction as src/learn/tree.rs does, the two gains need not be equal.
        weights = exact.Weights(3, 7)
        resistant, susceptible = bits(0, 1, 2), bits(3, 4, 5, 6, 7, 8, 9)
        row = bits(0, 3, 4, 5)
        other = (resistant | susceptible) & ~row
        first, second = exact.gain(weights, resistant, susceptible, row), exact.gain(weights, resistant, susceptible, other)
        self.assertEqual(first, second)
        self.assertGreater(first, 0)
        self.assertEqual(exact.best(weights, resistant, susceptible, [row, other], [0, 1])[1], [0, 1])


class Pins(unittest.TestCase):
    def test_the_committed_requirement_files(self):
        pins = check_env.pins()
        self.assertEqual(list(pins), ["requirements-build.txt", "requirements.txt", "requirements-scm.txt"])
        self.assertEqual([name for name, _, _ in pins["requirements-build.txt"]], ["setuptools", "wheel"])
        self.assertEqual(
            [name for name, _, _ in pins["requirements.txt"]],
            ["numpy", "scipy", "scikit-learn", "joblib", "threadpoolctl", "six", "xgboost"],
        )
        self.assertEqual([(name, version) for name, version, _ in pins["requirements-scm.txt"]], [("pyscm-ml", "1.1.2")])
        digests = [digest for file in pins.values() for _, _, digest in file]
        self.assertEqual(len(set(digests)), 10)

    def test_a_requirement_without_one_version_and_one_hash_is_refused(self):
        digest = "0" * 64
        good = f"# a comment\nnumpy==1.26.4 \\\n    --hash=sha256:{digest}\n"
        for bad in (
            "numpy==1.26.4\n",
            f"numpy>=1.26.4 \\\n    --hash=sha256:{digest}\n",
            f"numpy \\\n    --hash=sha256:{digest}\n",
            f"numpy==1.26.4 \\\n    --hash=sha256:{digest} \\\n    --hash=sha256:{digest}\n",
            f"numpy==1.26.4 \\\n    --hash=sha256:{digest[:-1]}\n",
            f"numpy==1.26.4 \\\n    --hash=md5:{digest}\n",
            good + "scipy\n",
            "--extra-index-url https://example.invalid/simple\n" + good,
            "# nothing pinned\n",
        ):
            with self.subTest(bad=bad), tempfile.TemporaryDirectory() as folder:
                for name in check_env.FILES:
                    (pathlib.Path(folder) / name).write_text(good, encoding="utf-8")
                self.assertEqual(check_env.pins(folder)["requirements.txt"], [("numpy", "1.26.4", digest)])
                (pathlib.Path(folder) / "requirements.txt").write_text(bad, encoding="utf-8")
                with self.assertRaisesRegex(ValueError, "requirements.txt: a line is no"):
                    check_env.pins(folder)


def small():
    """A table of two panels that has no problem: every outcome, and both kinds of check."""
    summary, differences = [], []
    outcomes = {
        ("a", "scm", "fits"): ("explained", 14, 12, 2, 0),
        ("a", "rf", "trees"): ("explained", 200, 199, 1, 0),
        ("a", "gbt", "end to end"): ("not applicable", 0, 0, 0, 0),
        ("b", "cart", "nodes"): ("unexplained", 9, 7, 1, 1),
        ("b", "gbt", "tree by tree"): ("unexplained", 400, 398, 1, 1),
    }
    for panel in ("a", "b"):
        for learner, check, units in table.CHECKS:
            outcome, compared, matched, explained, unexplained = outcomes.get((panel, learner, check), ("matched", 5, 5, 0, 0))
            summary.append(
                {"panel": panel, "samples": 10, "features": 3, "fit_samples": 8, "learner": learner, "check": check, "outcome": outcome,
                 "units": units, "compared": compared, "matched": matched, "explained": explained, "unexplained": unexplained}
            )  # fmt: skip

    def difference(panel, learner, check, unit, cause):
        differences.append({"panel": panel, "learner": learner, "check": check, "unit": unit, "what": "a rule", "cause": cause, "evidence": "numbers"})

    difference("a", "scm", "fits", "conjunction p=1", "the rule set")
    difference("a", "scm", "fits", "disjunction p=1", "the rule set")
    # a forest tree may have a row for each node
    difference("a", "rf", "trees", "tree 3", "rounding")
    difference("a", "rf", "trees", "tree 3", "rounding")
    difference("b", "cart", "nodes", "node 2", "rounding")
    difference("b", "cart", "nodes", "node 5", "")
    difference("b", "gbt", "tree by tree", "round 1", "equal split")
    difference("b", "gbt", "tree by tree", "round 2", "")
    return summary, differences


class Problems(unittest.TestCase):
    def test_a_table_with_no_problem(self):
        summary, differences = small()
        self.assertEqual(table.problems(summary, differences), [])
        self.assertEqual(table.problems(summary, differences, ["a", "b"]), [])

    def broken(self, change, panels_=None):
        summary, differences = small()
        change(summary, differences)
        return table.problems(summary, differences, panels_)

    def row(self, summary, panel, learner, check):
        return next(r for r in summary if (r["panel"], r["learner"], r["check"]) == (panel, learner, check))

    def test_each_broken_table_is_found(self):
        def set_row(panel, learner, check, **values):
            return lambda summary, differences: self.row(summary, panel, learner, check).update(values)

        cases = (
            ("the panels are not the fixture file's", lambda s, d: None, ["b", "a"]),
            ("the panels are not the fixture file's", lambda s, d: None, ["a", "b", "c"]),
            ("has not one row per panel and check", lambda s, d: s.pop(3), None),
            ("has not one row per panel and check", lambda s, d: s.insert(1, s.pop(0)), None),
            ("has not one row per panel and check", lambda s, d: s.append(dict(s[-1])), None),
            ("unknown outcome 'passed'", set_row("a", "cart", "nodes", outcome="passed"), None),
            ("the unit is not nodes", set_row("a", "cart", "nodes", units="trees"), None),
            ("sample or feature counts", set_row("a", "cart", "nodes", fit_samples=11), None),
            ("sample or feature counts", set_row("a", "cart", "nodes", features=0), None),
            ("do not add up to compared", set_row("a", "cart", "nodes", matched=4), None),
            ("do not add up to compared", set_row("a", "scm", "fits", compared=15), None),
            # a row said to match while a difference row names it
            ("a scm fits: matched, with a difference", set_row("a", "scm", "fits", outcome="matched"), None),
            ("a cart nodes: matched, with a difference", set_row("a", "cart", "nodes", matched=4, explained=1), None),
            # explained needs a difference row, and a cause on every one
            ("b rf trees: explained, without a cause", set_row("b", "rf", "trees", outcome="explained"), None),
            ("a scm fits: explained, without a cause", lambda s, d: d[0].update(cause=""), None),
            ("b cart nodes: explained, without a cause", set_row("b", "cart", "nodes", outcome="explained"), None),
            ("b cart nodes: unexplained, and no difference lacks a cause", lambda s, d: d[5].update(cause="found"), None),
            ("a gbt end to end: not applicable, with units", set_row("a", "gbt", "end to end", compared=1, matched=1), None),
            ("a scm fits: the unit counts are not the counts of its difference rows", lambda s, d: d.pop(0), None),
            ("b gbt tree by tree: the unit counts are not", set_row("b", "gbt", "tree by tree", matched=397, explained=2), None),
            ("a difference names ('c', 'scm', 'fits')", lambda s, d: d[0].update(panel="c"), None),
            ("names no unit or says not what differs", lambda s, d: d[2].update(unit=""), None),
            ("names no unit or says not what differs", lambda s, d: d[2].update(what=""), None),
        )
        for message, change, panels_ in cases:
            with self.subTest(message=message):
                found = self.broken(change, panels_)
                self.assertTrue(any(message in line for line in found), found)

    def test_totals(self):
        summary, _ = small()
        totals = {(row["learner"], row["check"]): row for row in table.totals(summary)}
        self.assertEqual(list(totals), [check[:2] for check in table.CHECKS])
        self.assertEqual(
            totals[("scm", "fits")],
            {"learner": "scm", "check": "fits", "units": "fits", "compared": 19, "units_matched": 17, "units_explained": 2, "units_unexplained": 0,
             "panels": {"matched": 1, "explained": 1, "unexplained": 0, "not applicable": 0}},
        )  # fmt: skip
        self.assertEqual(totals[("gbt", "end to end")]["panels"], {"matched": 1, "explained": 0, "unexplained": 0, "not applicable": 1})
        self.assertEqual(totals[("cart", "nodes")]["units_unexplained"], 1)


class Files(unittest.TestCase):
    def setUp(self):
        folder = tempfile.TemporaryDirectory()
        self.addCleanup(folder.cleanup)
        self.folder = pathlib.Path(folder.name)

    def write(self):
        summary, differences = small()
        table.write(self.folder / "summary.tsv", table.SUMMARY_COLUMNS, summary)
        table.write(self.folder / "differences.tsv", table.DIFFERENCE_COLUMNS, differences)
        (self.folder / "totals.json").write_text(json.dumps(table.totals(summary)), encoding="utf-8")
        return summary, differences

    def run_table(self):
        out = io.StringIO()
        with contextlib.redirect_stdout(out), self.assertRaises(SystemExit) as stop:
            table.main(["table.py", str(self.folder)])
        return stop.exception.code, out.getvalue()

    def test_what_is_written_is_what_is_read(self):
        summary, differences = self.write()
        self.assertEqual(table.read_summary(self.folder / "summary.tsv"), summary)
        self.assertEqual(table.read_differences(self.folder / "differences.tsv"), differences)

    def test_a_tab_or_a_line_break_in_a_field_is_refused(self):
        summary, differences = small()
        for text in ("a\tb", "a\nb", "a\rb"):
            broken = copy.deepcopy(differences)
            broken[0]["evidence"] = text
            with self.subTest(text=text), self.assertRaisesRegex(ValueError, "a tab or line break in evidence"):
                table.write(self.folder / "differences.tsv", table.DIFFERENCE_COLUMNS, broken)

    def test_a_file_of_another_shape_is_refused(self):
        self.write()
        path = self.folder / "summary.tsv"
        text = path.read_text(encoding="utf-8")
        for broken, message in (
            (text[:-1], "no line end at the end"),
            (text.replace("fit_samples", "training", 1), "the header is not"),
            (text.replace("\tscm\tfits\texplained\t", "\tscm\tfits\t", 1), "line 2: 11 fields"),
            (text.replace("\t14\t12\t2\t0\n", "\t14\t12\t2\tnone\n", 1), "unexplained of a is no count"),
            (text.replace("\t14\t12\t2\t0\n", "\t14\t12\t-2\t0\n", 1), "explained of a is no count"),
        ):
            path.write_text(broken, encoding="utf-8")
            with self.subTest(message=message), self.assertRaisesRegex(ValueError, message):
                table.read_summary(path)

    def test_the_script_exits_0_on_a_good_folder_and_1_on_a_bad_one(self):
        self.write()
        self.assertEqual(self.run_table(), (0, ""))
        totals = json.loads((self.folder / "totals.json").read_text(encoding="utf-8"))
        totals[0]["compared"] += 1
        (self.folder / "totals.json").write_text(json.dumps(totals), encoding="utf-8")
        self.assertEqual(self.run_table(), (1, "totals.json is not what summary.tsv adds up to\n"))

    def test_the_script_passes_the_committed_result(self):
        for name in ("summary.tsv", "differences.tsv", "totals.json"):
            shutil.copy(kit.REFERENCE / name, self.folder / name)
        self.assertEqual(self.run_table(), (0, ""))
        # one row of the committed summary said to match, though differences.tsv holds a difference for it
        path = self.folder / "summary.tsv"
        text = path.read_text(encoding="utf-8")
        self.assertEqual(text.count("\texplained\tnodes\t51\t50\t1\t0\n"), 1)
        path.write_text(text.replace("\texplained\tnodes\t51\t50\t1\t0\n", "\tmatched\tnodes\t51\t51\t0\t0\n"), encoding="utf-8")
        self.assertEqual(
            self.run_table(),
            (1, "campylobacter_jejuni__ciprofloxacin_lineage_fold0 cart nodes: matched, with a difference\ntotals.json is not what summary.tsv adds up to\n"),
        )
        # one cause removed from the committed differences
        path.write_text(text, encoding="utf-8")
        path = self.folder / "differences.tsv"
        text = path.read_text(encoding="utf-8")
        self.assertEqual(text.count("\tround 12\tat root.absent GRM splits on feature 11 and XGBoost on feature 0\tfloat32 rounding\t"), 1)
        path.write_text(text.replace("feature 11 and XGBoost on feature 0\tfloat32 rounding\t", "feature 11 and XGBoost on feature 0\t\t"), encoding="utf-8")
        code, printed = self.run_table()
        self.assertEqual(code, 1)
        self.assertIn("streptococcus_agalactiae__clindamycin_lineage_fold0 gbt tree by tree: explained, without a cause for every difference\n", printed)


if __name__ == "__main__":
    unittest.main()
