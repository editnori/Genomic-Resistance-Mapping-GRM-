"""scripts/c1-cutoff/score.py: the checks, the measures and the count (docs/cohort-draw.md).

Every result and report is SYNTHETIC (tests/c1_cutoff/common.py). Three
tests plant a fault in a result and check that the scorer's check fails.
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

    def test_planted_a_cohort_genome_scored_as_left_out_fails_the_sum_check(self):
        doc = common.result()
        doc["scored"]["samples"] += 1
        doc["scored"]["resistant"] += 1
        for t in doc["truth"]["learners"]:
            t["all"]["n"] += 1
            t["all"]["tp"] += 1
        self.assertEqual(failed(self.e, doc, self.rep), ["scored, copies and cohort sum to the dataset"])

    def test_planted_a_class_short_of_r_fails_the_cohort_check(self):
        doc = common.result()
        doc["cohort"]["resistant"] = 24
        self.assertIn("not refused, r of each class", failed(self.e, doc, self.rep))

    def test_planted_a_truth_from_the_cohort_folds_fails_the_truth_check(self):
        doc = common.result(truth_from_folds=True)
        self.assertEqual(failed(self.e, doc, self.rep), ["each truth scores the scored genomes and new lineages"])


class Scores(unittest.TestCase):
    def test_rows_hold_the_error_coverage_and_rules_measures(self):
        e = common.expect_row()
        rows = score.score_rows(e, common.result())
        self.assertEqual(len(rows), 10)  # lineage split, two truth sets, five learners
        row = rows[0]
        # Truth on all scored genomes: 140 - 50 - 2 = 88 genomes, 60 - 25 - 1 = 34 resistant.
        self.assertEqual((row["truth_recall_k"], row["truth_recall_n"]), (29, 34))
        self.assertEqual((row["truth_specificity_k"], row["truth_specificity_n"]), (44, 54))
        truth = (29 / 34 + 44 / 54) / 2
        self.assertAlmostEqual(row["truth"], truth)
        self.assertAlmostEqual(row["error"], abs(0.82 - truth))
        self.assertEqual((row["interval_lo"], row["interval_hi"]), (0.7, 0.9))
        self.assertTrue(row["covers"])
        self.assertEqual((row["estimate_very_major_k"], row["estimate_very_major_n"]), (5, 25))
        self.assertEqual(row["truth_log_loss_n"], 88)
        self.assertAlmostEqual(row["truth_log_loss"], 12.5 / 88)
        self.assertEqual(rows[5]["truth_set"], "new_lineage")
        self.assertEqual(rows[5]["truth_recall_n"], 10)

    def test_an_unavailable_lineage_split_leaves_the_draw_unavailable(self):
        rows = score.score_rows(common.expect_row(), common.result(lineage_unavailable="too few lineages"))
        self.assertTrue(all(not r["available"] and r["error"] is None for r in rows))
        self.assertEqual(rows[0]["estimate_unavailable"], "too few lineages")


class Count(unittest.TestCase):
    def line(self, d, line, err, draws="0,1,2,3,4", sp="synthia_coli", covered=3, n=5):
        return {"dataset": d, "species_slug": sp, "split": "lineage", "truth_set": "all", "learner": "rf", "line": line,
                "draws": 5, "available": len(draws.split(",")), "available_draws": draws, "mean_error": err,
                "covered": covered, "coverage_n": n}

    def test_better_worse_same_by_the_margin_and_three_paired_draws(self):
        rows = [
            self.line("a", "25", 0.10), self.line("a", "10", 0.13),          # worse by 0.03
            self.line("b", "25", 0.10), self.line("b", "10", 0.07),          # better by 0.03
            self.line("c", "25", 0.10), self.line("c", "10", 0.115),         # same
            self.line("d", "25", 0.10, draws="0,1,2"), self.line("d", "10", 0.5, draws="2,3,4"),  # one paired draw
        ]
        (c,) = score.count(rows)
        self.assertEqual((c["line"], c["datasets"], c["entered"], c["not_scored"]), ("10", 4, 3, 1))
        self.assertEqual((c["better"], c["worse"], c["same"]), (1, 1, 1))
        self.assertAlmostEqual(c["median_difference"], 0.015)
        # The release row reads the median, its interval and the coverage,
        # not the better and worse counts: one species, so the interval is
        # the median itself, below the margin, and the coverage is 25's.
        self.assertAlmostEqual(c["interval_hi"], 0.015)
        self.assertTrue(c["candidate"])
        # Coverage at 10 of 3 in 15 against 9 in 15 at 25: no candidate.
        for row in rows:
            if row["line"] == "10":
                row["covered"] = 1
        self.assertFalse(score.count(rows)[0]["candidate"])

    def test_the_species_interval_reads_the_101st_and_3900th_of_4000_resampled_medians(self):
        rows = []
        for k in range(6):
            d = f"sp{k}__drug"
            rows += [self.line(d, "25", 0.10, sp=f"sp{k}"), self.line(d, "15", 0.10 - 0.01 * k, sp=f"sp{k}")]
        (c,) = score.count(rows)
        self.assertEqual(c["species"], 6)
        self.assertLessEqual(c["interval_lo"], c["median_difference"])
        self.assertLessEqual(c["median_difference"], c["interval_hi"])
        self.assertTrue(c["interval_hi"] < score.MARGIN)
        self.assertTrue(c["candidate"])
        self.assertEqual(score.count(rows), [c])  # its resamples are seeded

    def test_the_twin_counts_against_25(self):
        rows = [self.line("a", "25", 0.10), self.line("a", "25-s45", 0.10)]
        (c,) = score.count(rows)
        self.assertEqual((c["line"], c["same"]), ("25-s45", 1))


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
            self.assertEqual(len(lines), 10)
            self.assertEqual({l["available"] for l in lines}, {"3"})
            self.assertNotIn("9000.", (tmp / "out" / "scores.tsv").read_text())
            (tmp / "res" / rows[0]["file"]).unlink()
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(score.main(args), 1)


if __name__ == "__main__":
    unittest.main()
