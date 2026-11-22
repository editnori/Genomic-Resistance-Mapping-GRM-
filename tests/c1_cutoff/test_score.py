"""scripts/c1-cutoff/score.py: the checks, the measures and the count (docs/cohort-draw.md).

Every result, report and line row is SYNTHETIC (tests/c1_cutoff/common.py
and the helpers below). Some tests plant a fault in a result and check that
the scorer's check fails; the others pin the measures and the count to
hand-worked figures, so a fault planted in score.py fails one of them.
"""

import csv
import json
import math
import tempfile
import unittest
from pathlib import Path

from tests.c1_cutoff import common

score = common.load("score")


def failed(e, doc, rep, sha="c" * 64):
    return [name for name, ok, _ in score.check(e, doc, rep, sha) if not ok]


class Intervals(unittest.TestCase):
    def test_wilson_and_clopper_pearson_match_published_values(self):
        lo, hi = score.wilson(3, 10)
        self.assertAlmostEqual(lo, 0.10779, places=4)
        self.assertAlmostEqual(hi, 0.60322, places=4)
        lo, hi = score.clopper_pearson(3, 10)
        self.assertAlmostEqual(lo, 0.066739, places=5)
        self.assertAlmostEqual(hi, 0.652453, places=5)
        # At the edges the exact interval has a closed form.
        self.assertEqual(score.clopper_pearson(0, 10)[0], 0.0)
        self.assertAlmostEqual(score.clopper_pearson(0, 10)[1], 1 - 0.025 ** (1 / 10), places=9)
        self.assertAlmostEqual(score.clopper_pearson(10, 10)[0], 0.025 ** (1 / 10), places=9)
        self.assertEqual(score.wilson(0, 0), (None, None))


class Checks(unittest.TestCase):
    def setUp(self):
        self.rep = common.report(common.samples())
        self.e = common.expect_row()

    def test_a_consistent_result_passes_every_check(self):
        self.assertEqual(failed(self.e, common.result(), self.rep), [])
        self.assertEqual(failed(common.expect_row(study="1"), common.result(study=True), self.rep), [])

    def test_a_moved_fingerprint_setting_seed_or_cohort_file_fails(self):
        self.assertEqual(failed(self.e, common.result(fingerprint="other"), self.rep), ["fingerprint is the release report's"])
        doc = common.result()
        doc["eval_config"]["folds"] = 10
        self.assertEqual(failed(self.e, doc, self.rep), ["settings are the report's but the seed"])
        self.assertEqual(failed(self.e, common.result(seed=45), self.rep), ["seed is the expected one"])
        self.assertEqual(failed(self.e, common.result(), self.rep, sha="d" * 64), ["cohort file sha256"])
        self.assertIn("each estimate scores the cohort", failed(self.e, common.result(study=True), self.rep))

    def test_the_study_estimate_carries_the_scheme_grm_writes(self):
        # grm learn --exp-cohort names the study estimate src/cutoff.rs's STUDY_SCHEME.
        src = (Path(__file__).resolve().parents[2] / "src" / "cutoff.rs").read_text()
        self.assertIn(f'pub const STUDY_SCHEME: &str = "{score.STUDY_SCHEME}";', src)
        doc = common.result(study=True)
        self.assertEqual([e["scheme"] for e in doc["estimates"]][2:], [score.STUDY_SCHEME])
        self.assertEqual(failed(common.expect_row(study="1"), doc, self.rep), [])
        doc["estimates"][2]["scheme"] = "study"
        self.assertIn("each estimate scores the cohort", failed(common.expect_row(study="1"), doc, self.rep))
        rows = score.score_rows(common.expect_row(study="1"), common.result(study=True))
        self.assertEqual(sorted({r["split"] for r in rows}), ["lineage", "random", "study"])

    def test_planted_an_estimate_that_is_not_its_calls_balanced_accuracy_fails(self):
        doc = common.result()
        self.assertAlmostEqual(score.balanced(doc["estimates"][1]["learners"][0]["calls"]), 0.82)
        doc["estimates"][1]["learners"][2]["balanced_accuracy"] = 0.9
        self.assertEqual(failed(self.e, doc, self.rep), ["each estimate is the balanced accuracy of its calls"])
        doc = common.result()
        doc["estimates"][0]["learners"][0]["balanced_accuracy"] = None
        self.assertEqual(failed(self.e, doc, self.rep), ["each estimate is the balanced accuracy of its calls"])

    def test_planted_a_cohort_genome_scored_as_left_out_fails_the_sum_check(self):
        doc = common.result()
        doc["scored"]["samples"] += 1
        doc["scored"]["resistant"] += 1
        for t in doc["truth"]["learners"]:
            t["all"]["n"] += 1
            t["all"]["tp"] += 1
        self.assertEqual(failed(self.e, doc, self.rep), ["scored, copies and cohort sum to the dataset"])

    def test_planted_a_cohort_projected_at_the_datasets_support_floor_fails(self):
        doc = common.result()
        self.assertEqual(doc["cohort"]["min_support"], 3)
        doc["cohort"]["min_support"] = 10  # support_floor of a 2,000-genome dataset
        self.assertEqual(failed(self.e, doc, self.rep), ["the cohort's support floor is its own"])

    def test_planted_a_class_short_of_r_fails_the_cohort_check(self):
        doc = common.result()
        doc["cohort"]["resistant"] = 24
        self.assertIn("not refused, r of each class", failed(self.e, doc, self.rep))

    def test_planted_a_truth_from_the_cohort_folds_fails_the_truth_check(self):
        doc = common.result(truth_from_folds=True)
        self.assertEqual(failed(self.e, doc, self.rep), ["each truth scores the scored genomes and new lineages"])


class Scores(unittest.TestCase):
    def test_rows_hold_the_error_both_intervals_and_rules_measures(self):
        e = common.expect_row()
        rows = score.score_rows(e, common.result())
        self.assertEqual(len(rows), 20)  # lineage and random split, two truth sets, five learners
        self.assertEqual([r["split"] for r in rows[::10]], ["lineage", "random"])
        row = rows[0]
        # Truth on all scored genomes: 140 - 50 - 2 = 88 genomes, 60 - 25 - 1 = 34 resistant.
        self.assertEqual((row["truth_recall_k"], row["truth_recall_n"]), (29, 34))
        self.assertEqual((row["truth_specificity_k"], row["truth_specificity_n"]), (44, 54))
        truth = (29 / 34 + 44 / 54) / 2
        self.assertAlmostEqual(row["truth"], truth)
        self.assertAlmostEqual(row["error"], abs(0.82 - truth))
        self.assertEqual((row["fold_range_lo"], row["fold_range_hi"]), (0.7, 0.9))
        self.assertAlmostEqual(row["fold_range_width"], 0.2)
        self.assertTrue(row["fold_range_covers"])
        # The levelled interval of 20 of 25 resistant and 21 of 25
        # susceptible called right, worked by hand from the Wilson 95%
        # intervals (0.608690, 0.911394) and (0.653464, 0.935965).
        self.assertAlmostEqual(row["estimate_from_calls"], 0.82)
        self.assertAlmostEqual(row["levelled_lo"], 0.686401, places=6)
        self.assertAlmostEqual(row["levelled_hi"], 0.893515, places=6)
        self.assertAlmostEqual(row["levelled_width"], 0.893515 - 0.686401, places=5)
        self.assertTrue(row["levelled_covers"])  # truth 0.834
        self.assertEqual((row["estimate_very_major_k"], row["estimate_very_major_n"]), (5, 25))
        self.assertEqual(row["truth_log_loss_n"], 88)
        self.assertAlmostEqual(row["truth_log_loss"], 12.5 / 88)
        self.assertEqual(rows[5]["truth_set"], "new_lineage")
        self.assertEqual(rows[5]["truth_recall_n"], 10)

    def test_the_levelled_interval_is_asymmetric_and_built_from_recall_and_specificity(self):
        # 10 of 25 resistant and 18 of 25 susceptible: Wilson (0.234033,
        # 0.592605) and (0.524234, 0.857161), so b 0.56 in (0.431675, 0.678227).
        b, lo, hi = score.levelled(common.calls(25, 25, 10, 18))
        self.assertAlmostEqual(b, 0.56)
        self.assertAlmostEqual(lo, 0.431675, places=6)
        self.assertAlmostEqual(hi, 0.678227, places=6)
        # Swapping the classes moves b to the same value and the interval with it.
        self.assertEqual(score.levelled(common.calls(25, 25, 18, 10))[0], b)
        self.assertEqual(score.levelled(common.calls(0, 25, 0, 18)), (None, None, None))

    def test_each_split_builds_its_levelled_interval_from_its_own_calls(self):
        doc = common.result(study=True)
        study = doc["estimates"][2]
        for x in study["learners"]:
            x["calls"] = common.calls(25, 25, 10, 18)
            x["balanced_accuracy"] = 0.56
        rows = score.score_rows(common.expect_row(study="1"), doc)
        by = {(r["split"], r["truth_set"], r["learner"]): r for r in rows}
        st, li = by[("study", "all", "rf")], by[("lineage", "all", "rf")]
        self.assertAlmostEqual(st["levelled_lo"], 0.431675, places=6)
        self.assertAlmostEqual(st["levelled_hi"], 0.678227, places=6)
        self.assertAlmostEqual(li["levelled_lo"], 0.686401, places=6)
        self.assertFalse(st["levelled_covers"])  # truth 0.834 above 0.678
        self.assertAlmostEqual(st["error"], abs(0.56 - st["truth"]))

    def test_an_edge_counts_as_covered(self):
        self.assertTrue(score.covers(0.5, 0.7, 0.5))
        self.assertTrue(score.covers(0.5, 0.7, 0.7))
        self.assertFalse(score.covers(0.5, 0.7, 0.70001))
        self.assertIsNone(score.covers(None, 0.7, 0.6))
        # Every genome called right in the cohort and in the scored set: the
        # levelled interval's upper edge is exactly 1, the truth's figure.
        doc = common.result(est_tp=25, est_tn=25, ba=1.0)
        for t in doc["truth"]["learners"]:
            a = t["all"]
            t["all"] = common.calls(a["tp"] + a["fn"], a["tn"] + a["fp"], a["tp"] + a["fn"], a["tn"] + a["fp"])
        row = score.score_rows(common.expect_row(), doc)[0]
        self.assertEqual((row["levelled_hi"], row["truth"], row["levelled_covers"]), (1.0, 1.0, True))

    def test_a_truth_on_the_fold_ranges_edge_is_covered(self):
        doc = common.result()
        truth = score.balanced(doc["truth"]["learners"][0]["all"])
        for folds in ((truth, 0.9, 0.95, 0.92, 0.93), (0.6, 0.65, 0.62, 0.61, truth)):
            row = score.score_rows(common.expect_row(), common.result(folds=folds))[0]
            self.assertEqual((row["truth"], row["fold_range_covers"]), (truth, True), folds)

    def test_a_fold_with_no_balanced_accuracy_leaves_the_draw_no_fold_range(self):
        row = score.score_rows(common.expect_row(), common.result(folds=(0.7, float("nan"), 0.75, 0.9, 0.85)))[0]
        self.assertTrue(row["available"])
        self.assertEqual((row["fold_range_lo"], row["fold_range_hi"], row["fold_range_covers"]), (None, None, None))
        self.assertTrue(row["levelled_covers"])

    def test_a_lines_figures_read_its_available_draws_only(self):
        ex = [common.expect_row(draw=j) for j in range(3)]
        rows = (score.score_rows(ex[0], common.result(draw=0)) + score.score_rows(ex[1], common.result(draw=1, ba=0.9))
                + score.score_rows(ex[2], common.result(draw=2, lineage_unavailable="too few lineages")))
        (x,) = [x for x in score.lines(rows, ex) if x["learner"] == "rf" and x["truth_set"] == "all" and x["split"] == "lineage"]
        truth = score.balanced(common.result()["truth"]["learners"][0]["all"])
        self.assertEqual((x["draws"], x["available"], x["available_draws"], x["fold_range_coverage_n"]), (3, 2, "0,1", 2))
        self.assertAlmostEqual(x["mean_error"], (abs(0.82 - truth) + abs(0.9 - truth)) / 2)
        self.assertEqual(x["errors"], ",".join(score.fmt(abs(b - truth)) for b in (0.82, 0.9)))
        self.assertEqual((x["levelled_covers"], x["levelled_covered"], x["levelled_coverage_n"]), ("1,1", 2, 2))
        self.assertEqual(sorted(score.per_draw(x)), [0, 1])

    def test_an_unavailable_lineage_split_leaves_the_draw_unavailable(self):
        rows = score.score_rows(common.expect_row(), common.result(lineage_unavailable="too few lineages"))
        lineage = [r for r in rows if r["split"] == "lineage"]
        self.assertTrue(lineage and all(not r["available"] and r["error"] is None for r in lineage))
        self.assertEqual(lineage[0]["estimate_unavailable"], "too few lineages")
        self.assertTrue(all(r["available"] for r in rows if r["split"] == "random"))


def line(d, line, errs, sp="synthia_coli", learner="rf", split="lineage", truth="all", draws=None,
         lcov=None, fcov=None):
    """A SYNTHETIC lines.tsv row: errs is one error a draw, or one error for every draw (0 to 4 by default)."""
    if not isinstance(errs, (list, tuple)):
        errs = [errs] * (len(draws) if draws is not None else 5)
    errs = list(errs)
    draws = list(draws) if draws is not None else list(range(len(errs)))
    lcov = list(lcov) if lcov is not None else [1] * len(draws)
    fcov = list(fcov) if fcov is not None else list(lcov)
    j = lambda xs: ",".join(score.fmt(x) for x in xs)
    return {"dataset": d, "species_slug": sp, "split": split, "truth_set": truth, "learner": learner, "line": line,
            "draws": 5, "available": len(draws), "available_draws": j(draws), "errors": j(errs),
            "levelled_covers": j(bool(c) for c in lcov), "levelled_widths": j(0.2 for _ in draws),
            "fold_range_covers": j(None if c is None else bool(c) for c in fcov), "fold_range_widths": j(0.3 for _ in draws)}


def get(rows, **key):
    """The one count row matching key; learner rows by default, levelled coverage by default."""
    key = {"row": "learner", "coverage": "levelled", **key}
    found = [r for r in rows if all(r.get(k) == v for k, v in key.items())]
    assert len(found) == 1, (key, len(found))
    return found[0]


class Count(unittest.TestCase):
    def test_better_worse_same_by_the_margin_and_three_paired_draws(self):
        rows = [
            line("a", "25", 0.10), line("a", "10", 0.13),          # worse by 0.03
            line("b", "25", 0.10), line("b", "10", 0.07),          # better by 0.03
            line("c", "25", 0.10), line("c", "10", 0.115),         # same
            line("d", "25", 0.10, draws=[0, 1, 2]), line("d", "10", 0.5, draws=[2, 3, 4]),  # one paired draw
        ]
        c = get(score.count(rows), line="10")
        self.assertEqual((c["datasets"], c["entered"], c["not_scored"], c["no_base_line"]), (4, 3, 1, 0))
        self.assertEqual((c["better"], c["worse"], c["same"]), (1, 1, 1))
        self.assertAlmostEqual(c["median_difference"], 0.015)
        # One species, so the interval is the median itself, below the margin,
        # and the coverage over the 15 shared results is 25's.
        self.assertAlmostEqual(c["interval_hi"], 0.015)
        self.assertEqual((c["covered"], c["coverage_n"], c["base_covered"], c["base_coverage_n"]), (15, 15, 15, 15))
        self.assertTrue(c["meets"])
        # A difference of exactly the margin is the same, not worse.
        self.assertEqual(get(score.count([line("a", "25", 0.10), line("a", "10", 0.12)]), line="10")["same"], 1)

    def test_the_difference_is_paired_on_the_shared_draws(self):
        # Line 10 has draws 0 to 2; 25 has 0 to 4 and its draws 3 and 4 are hard.
        rows = [line("a", "25", [0.10, 0.10, 0.10, 0.30, 0.30], lcov=[1, 1, 1, 0, 0]),
                line("a", "10", [0.11, 0.11, 0.11], lcov=[1, 1, 1])]
        c = get(score.count(rows), line="10")
        self.assertEqual(c["entered"], 1)
        self.assertAlmostEqual(c["median_difference"], 0.01)  # not 0.11 - 0.18
        self.assertEqual((c["same"], c["better"]), (1, 0))
        self.assertEqual(c["shared_draws"], 3)
        # Coverage over the same three results at both lines; each line's own beside.
        self.assertEqual((c["covered"], c["coverage_n"], c["base_covered"], c["base_coverage_n"]), (3, 3, 3, 3))
        self.assertEqual((c["line_available"], c["line_covered"], c["line_coverage_n"]), (3, 3, 3))
        self.assertEqual((c["base_available"], c["base_covered_own"], c["base_coverage_n_own"]), (5, 3, 5))
        self.assertAlmostEqual(c["median_width"], 0.2)
        self.assertTrue(c["meets"])

    def test_a_dataset_with_no_line_at_25_is_counted_as_not_scored(self):
        rows = [line("a", "25", 0.10), line("a", "10", 0.10), line("b", "10", 0.10)]
        c = get(score.count(rows), line="10")
        self.assertEqual((c["datasets"], c["entered"], c["not_scored"], c["no_base_line"]), (2, 1, 1, 1))
        self.assertEqual((c["line_available"], c["base_available"]), (10, 5))

    def test_an_interval_end_of_exactly_002_is_not_wholly_below_it(self):
        # One species, two datasets: differences -0.108 and 0.148, whose median
        # is exactly 0.02; in floating point -0.108 + 0.148 is 0.039999999999999994.
        rows = [line("a", "25", 0.2), line("a", "10", 0.092), line("b", "25", 0.2), line("b", "10", 0.348)]
        c = get(score.count(rows), line="10")
        self.assertEqual(c["entered"], 2)
        self.assertEqual((c["median_difference"], c["interval_lo"], c["interval_hi"]), (0.02, 0.02, 0.02))
        self.assertEqual((c["covered"], c["coverage_n"], c["base_covered"], c["base_coverage_n"]), (10, 10, 10, 10))
        self.assertFalse(c["meets"])
        # The written row says the same: an upper end of 0.02 is not below 0.02.
        self.assertEqual(score.fmt(c["interval_hi"]), "0.02")

    def test_coverage_within_005_of_25s_meets_at_the_edge(self):
        # 20 shared results: 19 covered at 10, 20 at 25, a difference of exactly 0.05.
        rows = []
        for k in range(4):
            rows += [line(f"d{k}", "25", 0.10, lcov=[1] * 5), line(f"d{k}", "10", 0.10, lcov=[0 if k == 0 else 1] + [1] * 4)]
        c = get(score.count(rows), line="10")
        self.assertEqual((c["covered"], c["coverage_n"], c["base_covered"], c["base_coverage_n"]), (19, 20, 20, 20))
        self.assertTrue(c["meets"])
        rows[1] = line("d0", "10", 0.10, lcov=[0, 0, 1, 1, 1])  # 18 of 20: 0.10 apart
        self.assertFalse(get(score.count(rows), line="10")["meets"])

    def test_the_species_interval_reads_the_101st_and_3900th_of_4000_resampled_medians(self):
        rows = []
        for k in range(6):
            d = f"sp{k}__drug"
            rows += [line(d, "25", 0.10, sp=f"sp{k}"), line(d, "15", 0.10 - 0.01 * k, sp=f"sp{k}")]
        c = get(score.count(rows), line="15")
        self.assertEqual(c["species"], 6)
        self.assertLessEqual(c["interval_lo"], c["median_difference"])
        self.assertLessEqual(c["median_difference"], c["interval_hi"])
        self.assertTrue(c["interval_hi"] < score.MARGIN)
        self.assertTrue(c["meets"])
        self.assertEqual(score.count(rows), score.count(rows))  # its resamples are seeded

    def test_a_rows_interval_does_not_move_when_rows_are_added(self):
        # 61 distinct differences in 20 species, so the 101st and 3,900th
        # resampled medians move with the resample stream.
        rows = []
        for k in range(61):
            sp, d = f"sp{k % 20}", f"sp{k % 20}__drug{k}"
            rows += [line(d, "25", 0.10, sp=sp), line(d, "15", round(0.10 + 0.03 * math.sin(1.7 * k + 0.3), 9), sp=sp)]
        alone = get(score.count(rows), line="15", learner="rf", split="lineage")
        self.assertNotEqual(alone["interval_lo"], alone["interval_hi"])
        more = list(rows)
        for k in range(61):
            sp, d = f"sp{k % 20}", f"sp{k % 20}__drug{k}"
            more += [line(d, "25", 0.10, sp=sp, learner="cart"), line(d, "15", 0.10 + (k % 7) * 0.01, sp=sp, learner="cart"),
                     line(d, "10", 0.10 + (k % 5) * 0.01, sp=sp), line(d, "25-s45", 0.10, sp=sp),
                     line(d, "25", 0.10, sp=sp, split="random"), line(d, "15", 0.11, sp=sp, split="random")]
        after = get(score.count(more), line="15", learner="rf", split="lineage")
        self.assertEqual((after["interval_lo"], after["interval_hi"]), (alone["interval_lo"], alone["interval_hi"]))
        self.assertEqual(score.row_seed("lineage", "all", "rf", "15"),
                         int(__import__("hashlib").sha256(b"lineage|all|rf|15").hexdigest(), 16))

    def learners(self, good, split="lineage", truth="all", line_="10"):
        """Rows of the five learners on 3 species: those in `good` gain nothing, the rest lose 0.05."""
        rows = []
        for k in range(3):
            d, sp = f"sp{k}__drug", f"sp{k}"
            for l in common.LEARNERS:
                rows += [line(d, "25", 0.10, sp=sp, learner=l, split=split, truth=truth),
                         line(d, line_, 0.10 if l in good else 0.15, sp=sp, learner=l, split=split, truth=truth)]
        return rows

    def test_a_line_is_a_candidate_when_3_of_5_learners_meet_it(self):
        out = score.count(self.learners({"scm", "rf", "lasso"}))
        c = get(out, row="line", line="10")
        self.assertEqual((c["learner"], c["learners_meeting"], c["meeting"], c["candidate"]), ("3 of 5", 3, "scm,rf,lasso", True))
        self.assertTrue(get(out, line="10", learner="rf")["meets"])
        self.assertFalse(get(out, line="10", learner="gbt")["meets"])
        c = get(score.count(self.learners({"scm", "rf"})), row="line", line="10")
        self.assertEqual((c["learners_meeting"], c["candidate"]), (2, False))

    def test_only_the_lineage_all_levelled_line_row_decides(self):
        rows = (self.learners({"scm", "rf"}) + self.learners(set(common.LEARNERS), split="random")
                + self.learners(set(common.LEARNERS), split="study") + self.learners(set(common.LEARNERS), truth="new_lineage")
                + self.learners(set(common.LEARNERS), line_="25-s45"))
        out = score.count(rows)
        deciding = [(r["row"], r["split"], r["truth_set"], r["learner"], r["line"], r["coverage"]) for r in out if r["decides"]]
        self.assertEqual(deciding, [("line", "lineage", "all", "3 of 5", "10", "levelled")])
        # It reads its own learners: 2 of 5 there, though 5 of 5 meet on the
        # random and study splits, on new lineages and in the twin.
        self.assertFalse(get(out, row="line", line="10", split="lineage", truth_set="all")["candidate"])
        for key in ({"split": "random"}, {"split": "study"}, {"truth_set": "new_lineage"}):
            self.assertTrue(get(out, row="line", line="10", **{"split": "lineage", "truth_set": "all", **key})["candidate"])
        twin = get(out, row="line", line="25-s45", split="lineage", truth_set="all")
        self.assertEqual((twin["decides"], twin["learners_meeting"], twin["candidate"]), (False, 5, True))

    def test_the_fold_ranges_coverage_is_written_and_decides_nothing(self):
        rows = []
        for k in range(3):
            d, sp = f"sp{k}__drug", f"sp{k}"
            for l in common.LEARNERS:
                rows += [line(d, "25", 0.10, sp=sp, learner=l, fcov=[0] * 5),
                         line(d, "10", 0.10, sp=sp, learner=l, fcov=[1] * 5)]  # the range widens at 10
        out = score.count(rows)
        lv, fr = get(out, line="10", learner="rf"), get(out, line="10", learner="rf", coverage="fold_range")
        self.assertEqual((lv["covered"], lv["base_covered"], lv["meets"]), (15, 15, True))
        self.assertEqual((fr["covered"], fr["base_covered"], fr["meets"], fr["decides"]), (15, 0, False, False))
        self.assertEqual((lv["interval_lo"], lv["interval_hi"]), (fr["interval_lo"], fr["interval_hi"]))
        self.assertTrue(get(out, row="line", line="10")["decides"])
        self.assertTrue(get(out, row="line", line="10")["candidate"])
        f = get(out, row="line", line="10", coverage="fold_range")
        self.assertEqual((f["candidate"], f["decides"]), (False, False))

    def test_the_twin_counts_against_25(self):
        rows = [line("a", "25", 0.10), line("a", "25-s45", 0.10)]
        c = get(score.count(rows), line="25-s45")
        self.assertEqual((c["same"], c["entered"]), (1, 1))

    def test_a_lines_row_whose_lists_disagree_is_refused(self):
        bad = line("a", "10", 0.10)
        bad["errors"] = "0.1,0.1"
        with self.assertRaises(ValueError):
            score.count([line("a", "25", 0.10), bad])


class Main(unittest.TestCase):
    def test_main_checks_scores_and_exits_1_on_a_failed_check(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            smp = common.samples()
            common.write_release(tmp / "rel", [common.report(smp)])
            (tmp / "coh").mkdir()
            (tmp / "res").mkdir()
            cohort = b'{"synthetic": true}\n'
            import hashlib
            sha = hashlib.sha256(cohort).hexdigest()
            rows = [common.expect_row(r=25, draw=j, cohort_sha=sha) for j in range(3)]
            for e in rows:
                (tmp / "coh" / e["cohort"]).write_bytes(cohort)
                doc = common.result(draw=int(e["draw"]), cohort_sha=sha)
                (tmp / "res" / e["file"]).write_text(json.dumps(doc))
            with open(tmp / "expect.tsv", "w") as f:
                cols = list(rows[0])
                f.write("\t".join(cols) + "\n")
                for e in rows:
                    f.write("\t".join(e[c] for c in cols) + "\n")
            args = [str(tmp / "expect.tsv"), str(tmp / "coh"), str(tmp / "res"), str(tmp / "rel"), str(tmp / "out")]
            import contextlib, io
            with contextlib.redirect_stdout(io.StringIO()) as out:
                self.assertEqual(score.main(args), 0)
            self.assertIn("3 of 3 results pass", out.getvalue())
            with open(tmp / "out" / "lines.tsv") as f:
                lines = list(csv.DictReader(f, delimiter="\t"))
            self.assertEqual(len(lines), 20)  # lineage and random split
            self.assertEqual({l["available"] for l in lines}, {"3"})
            # The count read back from the file equals the count from memory.
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(score.main(["--count", str(tmp / "again"), str(tmp / "out" / "lines.tsv")]), 0)
            self.assertEqual((tmp / "again" / "count.tsv").read_text(), (tmp / "out" / "count.tsv").read_text())
            with open(tmp / "out" / "count.tsv") as f:
                head = next(csv.reader(f, delimiter="\t"))
            self.assertEqual(head, score.COUNT_COLUMNS)
            self.assertNotIn("9000.", (tmp / "out" / "scores.tsv").read_text())
            (tmp / "res" / rows[0]["file"]).unlink()
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(score.main(args), 1)


if __name__ == "__main__":
    unittest.main()
