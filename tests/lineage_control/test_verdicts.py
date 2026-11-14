"""scripts/lineage-control/verdicts.py on synthetic reports.

Every report here is written in the test: a handful of isolates, their lab
calls, their lineage clusters and a few feature cards. No release is read.
"""

import base64
import contextlib
import importlib.util
import io
import json
import math
import pathlib
import tempfile
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]


def load_verdicts():
    spec = importlib.util.spec_from_file_location("lineage_control_verdicts", ROOT / "scripts" / "lineage-control" / "verdicts.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


V = load_verdicts()


def bits(flags):
    """Presence flags as a feature card stores them: base64 of little-endian 64-bit words."""
    value = sum(1 << i for i, flag in enumerate(flags) if flag)
    return base64.b64encode(value.to_bytes(8 * max(1, -(-len(flags) // 64)), "little")).decode()


def mask(flags):
    return sum(1 << i for i, flag in enumerate(flags) if flag)


def cluster(with_term_resistant, with_term_susceptible, without_resistant, without_susceptible):
    """One lineage cluster as (term, resistant) pairs."""
    return (
        [(True, True)] * with_term_resistant
        + [(True, False)] * with_term_susceptible
        + [(False, True)] * without_resistant
        + [(False, False)] * without_susceptible
    )


def report(clusters, kind="conjunction", rules=None, extra_cards=(), dataset="synthetic__drug"):
    """A report cut down to what verdicts.py reads. Feature 1 follows the clusters' term flags."""
    isolates = [(c, term, resistant) for c, members in enumerate(clusters) for term, resistant in members]
    cards = [{"feature": 1, "bits": bits([term for _, term, _ in isolates])}, *extra_cards]
    return {
        "id": dataset,
        "evaluation": {"revision": 14},
        "input_fingerprint": "0123456789abcdef",
        "samples": [{"resistant": resistant, "lineage": c} for c, _, resistant in isolates],
        "features": cards,
        "models": {"scm": {"kind": kind, "rules": rules or [{"feature": 1, "presence": True}]}},
    }


class Statistic(unittest.TestCase):
    def test_five_equal_clusters_give_the_hand_computed_test(self):
        # Each cluster: 3 resistant and 1 susceptible with the term, 1 and 3 without. Expected with the term
        # and resistant is 4 * 4 / 8 = 2, so the excess is 1 a cluster; the variance is 4*4*4*4 / (8*8*7).
        clusters = [cluster(3, 1, 1, 3)] * 5
        row = V.report_terms(report(clusters), "sha")[0]
        variance = 5 * 256 / 448
        self.assertEqual((row["isolates"], row["context_isolates"], row["informative_strata"]), (40, 40, 5))
        self.assertAlmostEqual(row["excess"], 5.0)
        self.assertAlmostEqual(row["p"], math.erfc(math.sqrt((5.0 - 0.5) ** 2 / variance / 2)))
        self.assertAlmostEqual(row["p"], 0.0077623456, places=9)

    def test_a_cluster_without_both_calls_and_both_states_is_not_informative(self):
        clusters = [
            cluster(3, 1, 1, 3),
            cluster(4, 0, 2, 0),  # every isolate resistant
            cluster(2, 2, 0, 0),  # every isolate has the term
            cluster(0, 0, 3, 3),  # no isolate has the term
            cluster(1, 0, 0, 0),  # one isolate
        ]
        row = V.report_terms(report(clusters), "sha")[0]
        self.assertEqual(row["informative_strata"], 1)
        self.assertAlmostEqual(row["excess"], 1.0)

    def test_no_informative_cluster_gives_no_p(self):
        row = V.report_terms(report([cluster(4, 0, 2, 0), cluster(0, 0, 3, 3)]), "sha")[0]
        self.assertEqual((row["informative_strata"], row["p"]), (0, None))

    def test_an_absence_rule_is_tested_on_the_isolates_without_the_feature(self):
        clusters = [cluster(3, 1, 1, 3)] * 5
        present = V.report_terms(report(clusters), "sha")[0]
        absent = V.report_terms(report(clusters, rules=[{"feature": 1, "presence": False}]), "sha")[0]
        self.assertAlmostEqual(absent["excess"], -present["excess"])
        self.assertAlmostEqual(absent["p"], present["p"])
        self.assertEqual(absent["presence"], 0)

    def test_the_test_agrees_with_counting_the_tables_one_isolate_at_a_time(self):
        # An uneven set of clusters, counted again here without bit sets.
        clusters = [cluster(5, 1, 2, 4), cluster(1, 3, 2, 2), cluster(7, 0, 1, 6), cluster(2, 2, 2, 2), cluster(0, 1, 4, 1), cluster(3, 0, 0, 9)]
        excess = variance = 0.0
        informative = 0
        for a, b, c, d in [(5, 1, 2, 4), (1, 3, 2, 2), (7, 0, 1, 6), (2, 2, 2, 2), (0, 1, 4, 1), (3, 0, 0, 9)]:
            n = a + b + c + d
            v = (a + b) * (c + d) * (a + c) * (b + d) / (n * n * (n - 1))
            if v > 0:
                informative += 1
                excess += a - (a + b) * (a + c) / n
                variance += v
        row = V.report_terms(report(clusters), "sha")[0]
        self.assertEqual(row["informative_strata"], informative)
        self.assertAlmostEqual(row["excess"], excess)
        self.assertAlmostEqual(row["p"], math.erfc(math.sqrt(max(abs(excess) - 0.5, 0.0) ** 2 / variance / 2)))


class Context(unittest.TestCase):
    def setUp(self):
        # Eight isolates in two clusters. Feature 1 holds for the first four, feature 2 for isolates 0, 1, 4 and 5.
        self.clusters = [cluster(2, 2, 0, 0), cluster(0, 0, 2, 2)]
        self.second = {"feature": 2, "bits": bits([True, True, False, False, True, True, False, False])}
        self.rules = [{"feature": 1, "presence": True}, {"feature": 2, "presence": True}]

    def test_a_conjunction_tests_a_later_term_on_the_isolates_the_earlier_terms_hold_for(self):
        rows = V.report_terms(report(self.clusters, "conjunction", self.rules, [self.second]), "sha")
        self.assertEqual([r["context_isolates"] for r in rows], [8, 4])
        self.assertEqual([r["position"] for r in rows], [1, 2])

    def test_a_disjunction_tests_a_later_term_on_the_isolates_no_earlier_term_holds_for(self):
        rules = [{"feature": 2, "presence": True}, {"feature": 1, "presence": False}]
        rows = V.report_terms(report(self.clusters, "disjunction", rules, [self.second]), "sha")
        # Feature 2 holds for four isolates, so the second term is tested on the other four.
        self.assertEqual([r["context_isolates"] for r in rows], [8, 4])

    def test_a_rule_whose_feature_has_no_card_holds_for_no_isolate(self):
        rows = V.report_terms(report(self.clusters, "conjunction", [{"feature": 9, "presence": True}, {"feature": 1, "presence": True}]), "sha")
        self.assertEqual((rows[0]["informative_strata"], rows[0]["p"]), (0, None))
        self.assertEqual(rows[1]["context_isolates"], 0)

    def test_another_kind_of_model_stops_the_run(self):
        with self.assertRaises(SystemExit):
            V.report_terms(report(self.clusters, "majority"), "sha")

    def test_presence_bits_are_read_as_the_report_stores_them(self):
        flags = [i % 3 == 0 for i in range(130)]
        self.assertEqual(V.bitset(bits(flags), 130), mask(flags))
        # Bits past the isolates are dropped.
        self.assertEqual(V.bitset(bits([True] * 130), 70), (1 << 70) - 1)


class Verdict(unittest.TestCase):
    def row(self, strata, excess, p):
        return {"informative_strata": strata, "excess": excess, "p": p}

    def test_benjamini_hochberg_adjusts_in_the_given_order(self):
        adjusted = V.benjamini_hochberg([0.01, 0.04, 0.03, 0.005])
        for got, want in zip(adjusted, [0.02, 0.04, 0.04, 0.02]):
            self.assertAlmostEqual(got, want)
        self.assertEqual(V.benjamini_hochberg([]), [])
        self.assertEqual(V.benjamini_hochberg([0.5, 0.5]), [0.5, 0.5])

    def test_each_verdict(self):
        rows = [
            self.row(5, 4.0, 0.001),  # assessable, under the level, more resistant than expected
            self.row(9, -4.0, 0.001),  # under the level, fewer
            self.row(6, 1.0, 0.6),  # assessable, not under the level
            self.row(4, 9.0, 1e-12),  # one cluster short, whatever its p
            self.row(0, 0.0, None),  # no informative cluster
        ]
        V.judge(rows)
        self.assertEqual([r["verdict"] for r in rows], ["seen", "reversed", "not_shown", "cannot_test", "cannot_test"])
        # Only the three assessable terms are adjusted, among themselves.
        self.assertEqual([r["q"] is None for r in rows], [False, False, False, True, True])
        self.assertAlmostEqual(rows[0]["q"], 0.0015)
        self.assertAlmostEqual(rows[2]["q"], 0.6)

    def test_a_term_is_adjusted_among_every_assessable_term_of_the_release(self):
        # Alone, 0.04 is under the level; among twenty assessable terms with larger p it is not.
        self.assertEqual(V.judge([self.row(5, 1.0, 0.04)])[0]["verdict"], "seen")
        many = [self.row(5, 1.0, 0.04)] + [self.row(5, 1.0, 0.5 + i / 100) for i in range(20)]
        self.assertEqual(V.judge(many)[0]["verdict"], "not_shown")

    def test_the_level_and_the_cluster_minimum_are_the_ones_fixed_before_the_results(self):
        self.assertEqual((V.MIN_STRATA, V.ALPHA), (5, 0.05))
        self.assertEqual(V.VERDICTS, ("seen", "not_shown", "cannot_test", "reversed"))


class Table(unittest.TestCase):
    def test_a_release_directory_gives_one_row_a_term_under_the_table_columns(self):
        strong = report([cluster(6, 0, 0, 6)] * 6, dataset="a__drug")
        weak = report([cluster(3, 1, 1, 3)] * 2, dataset="b__drug")
        with tempfile.TemporaryDirectory() as tmp:
            datasets = pathlib.Path(tmp) / "datasets"
            datasets.mkdir()
            for r in (strong, weak):
                (datasets / f"{r['id']}.json").write_text(json.dumps(r), encoding="utf-8")
            out = io.StringIO()
            with contextlib.redirect_stdout(out):
                V.main(["verdicts.py", tmp])
        lines = out.getvalue().splitlines()
        self.assertEqual(lines[0].split("\t"), list(V.COLUMNS))
        rows = [dict(zip(V.COLUMNS, line.split("\t"))) for line in lines[1:]]
        self.assertEqual([(r["dataset"], r["position"], r["verdict"]) for r in rows], [("a__drug", "1", "seen"), ("b__drug", "1", "cannot_test")])
        self.assertEqual([(r["revision"], r["input_fingerprint"], r["isolates"]) for r in rows], [("14", "0123456789abcdef", "72"), ("14", "0123456789abcdef", "16")])
        self.assertEqual(rows[1]["q"], "")
        self.assertTrue(all(len(r["report_sha256"]) == 64 for r in rows))

    def test_a_directory_without_reports_stops_the_run(self):
        with tempfile.TemporaryDirectory() as tmp, self.assertRaises(SystemExit):
            V.main(["verdicts.py", tmp])


if __name__ == "__main__":
    unittest.main()
