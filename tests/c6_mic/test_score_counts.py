"""scripts/c6-mic/score.py: denominators and checks that a balanced SYNTHETIC
cohort cannot tell apart (docs/mic-plan.md).

Cohort B has five folds of two made-up genomes, 3 resistant and 7
susceptible, so a rate over the wrong class shows. r_above is 2 in every
fold, so a step above 2 is called R:

  fold  MIC        call  released  ordered
  0     exact 1    S      1.0       1
  0     exact 3    R      1.0       3
  1     exact 0    S      3.0       0
  1     exact 1    S      0.0       1
  2     exact 4    R      4.0       4
  2     exact 0    S      0.0       0
  3     exact 1    S      1.0       1
  3     exact 2    S      2.0       2
  4     above 5    R      5.0       5
  4     exact -1   S     -1.0      -1
"""

import unittest

from tests.c6_mic import common
from tests.c6_mic.test_score import by_arm, cohort_a, run

R = common.row


def cohort_b():
    return [
        R(0, 1, "Exact", False, 1.0, 1), R(0, 3, "Exact", True, 1.0, 3),
        R(1, 0, "Exact", False, 3.0, 0), R(1, 1, "Exact", False, 0.0, 1),
        R(2, 4, "Exact", True, 4.0, 4), R(2, 0, "Exact", False, 0.0, 0),
        R(3, 1, "Exact", False, 1.0, 1), R(3, 2, "Exact", False, 2.0, 2),
        R(4, 5, "Above", True, 5.0, 5), R(4, -1, "Exact", False, -1.0, -1),
    ]


class Denominators(unittest.TestCase):
    def test_each_error_rate_is_taken_over_its_own_class(self):
        # Released calls S S | R S | R S | S S | R S against S R | S S | R S | S S | R S:
        # 8 of 10 agree, 1 of 3 resistant called S (fold 0), 1 of 7
        # susceptible called R (fold 1). Ordered calls every genome right:
        # 10 of 10, 0 of 3, 0 of 7.
        code, _, scores = run([(common.expect_row(), common.result(cohort_b()))])
        self.assertEqual(code, 0)
        arm = by_arm(scores)
        got = {a: tuple(arm[a][k] for k in ["category_count", "category_n", "very_major_count",
                                             "very_major_n", "major_count", "major_n"])
               for a in ("released", "ordered")}
        self.assertEqual(got["released"], ("8", "10", "1", "3", "1", "7"))
        self.assertEqual(got["ordered"], ("10", "10", "0", "3", "0", "7"))
        self.assertEqual(arm["released"]["very_major_rate"], f"{1 / 3:.6f}")
        self.assertEqual(arm["released"]["major_rate"], f"{1 / 7:.6f}")


class FoldChecks(unittest.TestCase):
    def test_each_fold_count_is_recounted_from_the_rows(self):
        for key in ["test_mic", "test_exact", "test_resistant"]:
            doc = common.result(cohort_a())
            doc["folds"][1][key] += 1
            code, checks, scores = run([(common.expect_row(), doc)])
            self.assertEqual((code, checks[0]["fold_counts"]), (1, "differs"), key)
            self.assertIsNone(scores, key)

    def test_a_row_median_unlike_its_folds_stops_the_scorer(self):
        # Fold 0's median is 3 and its first row reads 4: the scorer would
        # count the row's 4, not the fold's 3.
        doc = common.result(cohort_a())
        doc["tested"][0]["median"] = doc["folds"][0]["median"] + 1
        code, checks, scores = run([(common.expect_row(), doc)])
        self.assertEqual((code, checks[0]["median_check"]), (1, "differs"))
        self.assertIsNone(scores)


class SeedChecks(unittest.TestCase):
    def test_a_seed_42_result_with_the_report_key_but_other_settings_anchors_nothing(self):
        # The report's key with 100 trees cannot come from grm learn, since the
        # key hashes the settings; the scorer still takes it as no anchor, so
        # the seed-43 result has nothing to be checked against.
        s42 = common.result(cohort_a(), seed=42)
        s42["eval_config"]["forest"]["n_trees"] = 100
        exp43 = common.expect_row(seed=43, plan="lineage-s43")
        code, checks, scores = run([(common.expect_row(), s42), (exp43, common.result(cohort_a(), seed=43))])
        self.assertEqual(code, 1)
        self.assertEqual([c["fingerprint"] for c in checks], ["differs", "differs"])
        self.assertIsNone(scores)

    def test_a_result_whose_settings_seed_is_not_its_plans_differs(self):
        # The plan says seed 43 and the settings say 44.
        s42 = (common.expect_row(seed=42, plan="lineage-s42"), common.result(cohort_a(), seed=42))
        doc = common.result(cohort_a(), seed=43)
        doc["eval_config"]["seed"] = 44
        doc["input_fingerprint"] = "f00d44"
        code, checks, scores = run([s42, (common.expect_row(seed=43, plan="lineage-s43"), doc)])
        self.assertEqual((code, checks[1]["fingerprint"]), (1, "differs"))
        self.assertIsNone(scores)


if __name__ == "__main__":
    unittest.main()
