"""scripts/c6-mic/score.py on SYNTHETIC results (docs/mic-plan.md).

Each cohort below is made up, and each expected figure is counted by hand in
the comment beside it. Cohort A has five folds of two genomes:

  fold  MIC        call  released  ordered  training median (other folds' exact MICs)
  0     exact 1    S      1.4       1       3 of 0 2 3 4 5
  0     at most -3 S     -2.6      -3
  1     exact 2    S      0.5       4       3 of 0 1 3 4 5
  1     above 6    R      6.2       6
  2     exact 4    R      2.5       4       1 of 0 1 2 5
  2     exact 3    R      1.49      3
  3     exact 0    S     -0.5       0       3 of 1 2 3 4 5
  3     at most -3 S     -1.5      -3
  4     exact 5    R      5.5       5       2 of 0 1 2 3 4
  4     above 6    R      5.0       6

A censored value sits at each end: at most -3 and above 6.
"""

import contextlib
import io
import json
import os
import tempfile
import unittest

from tests.c6_mic import common

score = common.load("score")
R = common.row


def cohort_a():
    return [
        R(0, 1, "Exact", False, 1.4, 1), R(0, -3, "AtMost", False, -2.6, -3),
        R(1, 2, "Exact", False, 0.5, 4), R(1, 6, "Above", True, 6.2, 6),
        R(2, 4, "Exact", True, 2.5, 4), R(2, 3, "Exact", True, 1.49, 3),
        R(3, 0, "Exact", False, -0.5, 0), R(3, -3, "AtMost", False, -1.5, -3),
        R(4, 5, "Exact", True, 5.5, 5), R(4, 6, "Above", True, 5.0, 6),
    ]


def run(entries, probe=None):
    """score.main on a SYNTHETIC run: (exit code, checks, scores or None).
    `probe` maps a file name to a probe result document."""
    with tempfile.TemporaryDirectory() as tmp:
        expect, results, release = common.write_run(tmp, entries)
        out = os.path.join(tmp, "out")
        args = [str(expect), str(results), str(release), out]
        if probe is not None:
            os.makedirs(os.path.join(tmp, "probe"))
            for name, doc in probe.items():
                with open(os.path.join(tmp, "probe", name), "w") as f:
                    json.dump(doc, f)
            args.append(os.path.join(tmp, "probe"))
        with contextlib.redirect_stdout(io.StringIO()):
            code = score.main(args)
        _, checks = common.read_tsv(os.path.join(out, "checks.tsv"))
        path = os.path.join(out, "scores.tsv")
        scores = common.read_tsv(path) if os.path.exists(path) else None
        return code, checks, scores


def by_arm(scores):
    head, rows = scores
    return {r["arm"]: r for r in rows}


class EaExact(unittest.TestCase):
    def test_ea_exact_reads_exact_mics_only(self):
        # Released steps 1, -3 | 1, 6 | 3, 1 | -1, -2 | 6, 5. Exact MICs 1 2 4 3 0 5
        # against steps 1 1 3 1 -1 6: within one for all but exact 3 (step 1), so
        # 5 of 6. The four censored values would add 4 more if read as exact.
        code, _, scores = run([(common.expect_row(), common.result(cohort_a()))])
        self.assertEqual(code, 0)
        arm = by_arm(scores)
        self.assertEqual((arm["released"]["ea_exact_count"], arm["released"]["ea_exact_n"]), ("5", "6"))
        self.assertEqual(arm["released"]["ea_exact"], f"{5 / 6:.6f}")
        # Ordered steps on the exact MICs: 1 4 4 3 0 5, all within one but 4 for exact 2.
        self.assertEqual((arm["ordered"]["ea_exact_count"], arm["ordered"]["ea_exact_n"]), ("5", "6"))
        # Training medians on the exact MICs: 3 3 1 1 3 2, within one only for exact 2.
        self.assertEqual((arm["median"]["ea_exact_count"], arm["median"]["ea_exact_n"]), ("1", "6"))

    def test_the_released_margin_rounds_half_away_from_zero(self):
        # Exact 4 at margin 2.5 is step 3, within one; rounded half to even it
        # would be step 2 and miss. Exact 0 at -0.5 is step -1, still within one.
        self.assertEqual([score.rround(v) for v in (2.5, -2.5, 0.5, -0.5, 1.49, -1.5)], [3, -3, 1, -1, 1, -2])
        code, _, scores = run([(common.expect_row(), common.result(cohort_a()))])
        self.assertEqual(by_arm(scores)["released"]["ea_exact_count"], "5")

    def test_the_training_median_can_beat_the_model_and_the_scorer_only_counts(self):
        # SYNTHETIC: ten genomes, all exact 2, two a fold. The released margin is
        # 5.0 everywhere (0 of 10 within one), the ordered step 2 and every fold's
        # median 2 (10 of 10). The scorer writes both rows and names no winner.
        rows = [R(i % 5, 2, "Exact", i % 2 == 0, 5.0, 2) for i in range(10)]
        code, _, scores = run([(common.expect_row(), common.result(rows))])
        self.assertEqual(code, 0)
        arm = by_arm(scores)
        self.assertEqual(arm["released"]["ea_exact_count"], "0")
        self.assertEqual(arm["median"]["ea_exact_count"], "10")
        self.assertEqual(arm["ordered"]["ea_exact_count"], "10")
        head, _ = scores
        self.assertFalse([c for c in head if "win" in c or "better" in c or "diff" in c])


class Category(unittest.TestCase):
    def test_category_agreement_and_both_error_rates_count_every_tested_genome(self):
        # r_above 2 in every fold, so a step above 2 is called R.
        # Released calls S S | S R | R S | S S | R R against S S | S R | R R | S S | R R:
        # 9 of 10 agree, 1 of 5 resistant called S, 0 of 5 susceptible called R.
        # Ordered calls S S | R R | R R | S S | R R: 9 of 10, 0 of 5, 1 of 5.
        # Median calls R R | R R | S S | R R | S S: 1 of 10, 4 of 5, 5 of 5.
        # The censored values count here: they have a call.
        code, _, scores = run([(common.expect_row(), common.result(cohort_a()))])
        arm = by_arm(scores)
        got = {a: tuple(arm[a][k] for k in ["category_count", "category_n", "very_major_count",
                                             "very_major_n", "major_count", "major_n"])
               for a in arm}
        self.assertEqual(got["released"], ("9", "10", "1", "5", "0", "5"))
        self.assertEqual(got["ordered"], ("9", "10", "0", "5", "1", "5"))
        self.assertEqual(got["median"], ("1", "10", "4", "5", "5", "5"))
        self.assertEqual(arm["released"]["very_major_rate"], f"{1 / 5:.6f}")
        self.assertEqual(arm["median"]["major_rate"], f"{1:.6f}")

    def test_a_fold_without_a_breakpoint_leaves_the_category_figures_empty(self):
        # Fold 2 has MIC test genomes and no breakpoint: no category figure is
        # taken over the other four folds. EA-exact still reads all folds.
        doc = common.result(cohort_a(), r_above=(2, 2, None, 2, 2))
        code, _, scores = run([(common.expect_row(), doc)])
        self.assertEqual(code, 0)
        for arm in by_arm(scores).values():
            self.assertEqual(arm["folds_without_breakpoint"], "1")
            for key in ["category_agreement", "category_count", "very_major_rate", "very_major_count",
                        "major_rate", "major_count"]:
                self.assertEqual(arm[key], "", (arm["arm"], key))
            self.assertEqual((arm["category_n"], arm["very_major_n"], arm["major_n"]), ("10", "5", "5"))
            self.assertNotEqual(arm["ea_exact_count"], "")


class Folds(unittest.TestCase):
    def test_a_fold_with_no_mic_test_genome_is_not_scored_and_not_hidden(self):
        # Cohort A without fold 3's two genomes: four folds score, fold 3 holds
        # no median and no breakpoint, and the medians are recounted without it:
        # fold 0 from 2 4 3 5 is 3, fold 1 from 1 4 3 5 is 3, fold 2 from 1 2 5
        # is 2, fold 4 from 1 2 4 3 is 2.
        rows = [r for r in cohort_a() if r["fold"] != 3]
        doc = common.result(rows)
        self.assertEqual([f["median"] for f in doc["folds"]], [3, 3, 2, None, 2])
        code, checks, scores = run([(common.expect_row(), doc)])
        self.assertEqual(code, 0)
        self.assertEqual(checks[0]["median_check"], "equal")
        for arm in by_arm(scores).values():
            self.assertEqual((arm["folds"], arm["folds_scored"], arm["folds_undefined"]), ("5", "4", "0"))
            self.assertEqual(arm["folds_without_breakpoint"], "0")
            self.assertEqual((arm["ea_exact_n"], arm["category_n"]), ("5", "8"))

    def test_an_unavailable_arm_writes_no_figure_and_counts_its_folds(self):
        # The ordered arm fails in fold 2, so grm learn fits it in folds 0 and 1
        # only and wipes its predictions. Its row keeps the denominators and
        # counts three folds undefined; the other arms are scored.
        doc = common.result(cohort_a(), ordered_from_fold=2)
        code, checks, scores = run([(common.expect_row(), doc)])
        self.assertEqual(code, 0)
        self.assertEqual(checks[0]["predictions"], "whole")
        arm = by_arm(scores)
        self.assertIn("fold 2", arm["ordered"]["unavailable"])
        self.assertEqual(arm["ordered"]["folds_undefined"], "3")
        for key in ["ea_exact", "ea_exact_count", "category_agreement", "category_count"]:
            self.assertEqual(arm["ordered"][key], "", key)
        self.assertEqual((arm["ordered"]["ea_exact_n"], arm["ordered"]["category_n"]), ("6", "10"))
        self.assertEqual(arm["released"]["ea_exact_count"], "5")

    def test_a_result_under_the_gate_writes_rows_that_say_why(self):
        doc = common.result(cohort_a(), unavailable="40 genomes with an MIC and 20 exact, under the gate of 100 and 50")
        code, checks, scores = run([(common.expect_row(), doc)])
        self.assertEqual(code, 0)
        arm = by_arm(scores)
        self.assertEqual(sorted(arm), ["median", "ordered", "released"])
        for row in arm.values():
            self.assertIn("under the gate", row["unavailable"])
            self.assertEqual((row["ea_exact"], row["folds"]), ("", "0"))

    def test_an_available_arm_missing_a_prediction_stops_the_scorer(self):
        # SYNTHETIC fault in a result: the released arm is available but one row
        # has no prediction. Scoring the rest would shrink the denominator
        # without saying so, so nothing is scored.
        doc = common.result(cohort_a())
        doc["tested"][4]["released"] = None
        code, checks, scores = run([(common.expect_row(), doc)])
        self.assertEqual(code, 1)
        self.assertEqual((checks[0]["predictions"], checks[0]["ok"]), ("partial", "no"))
        self.assertIsNone(scores)

    def test_an_unavailable_arm_that_still_predicts_some_folds_stops_the_scorer(self):
        # SYNTHETIC fault in a result: the ordered arm is marked unavailable from
        # fold 2 on, yet folds 0 and 1 keep their predictions. grm learn wipes
        # them; a file that keeps them is not scored over the folds it fitted.
        doc = common.result(cohort_a(), ordered_from_fold=2)
        for r, steps in zip(doc["tested"], cohort_a()):
            if r["fold"] < 2:
                r["ordered"] = steps["ordered"]
        code, checks, scores = run([(common.expect_row(), doc)])
        self.assertEqual(code, 1)
        self.assertEqual(checks[0]["predictions"], "partial")
        self.assertIsNone(scores)


class Checks(unittest.TestCase):
    def test_the_training_median_is_cross_checked_from_the_other_folds(self):
        # The medians 3 3 1 3 2 come from the other folds only. A result whose
        # medians read the test fold too (2 in every fold, from all six exact
        # MICs) differs, and the run stops.
        good = common.result(cohort_a())
        code, checks, _ = run([(common.expect_row(), good)])
        self.assertEqual((code, checks[0]["median_check"]), (0, "equal"))
        bad = common.result(cohort_a())
        for f in bad["folds"]:
            f["median"] = 2
        for r in bad["tested"]:
            r["median"] = 2
        code, checks, scores = run([(common.expect_row(), bad)])
        self.assertEqual((code, checks[0]["median_check"]), (1, "differs"))
        self.assertIsNone(scores)

    def test_a_missing_result_file_stops_the_scorer(self):
        # grm learn writes no file for a dataset it skips and exits 0.
        entries = [
            (common.expect_row(seed=42, plan="lineage-s42"), common.result(cohort_a(), seed=42)),
            (common.expect_row(seed=43, plan="lineage-s43"), None),
        ]
        code, checks, scores = run(entries)
        self.assertEqual(code, 2)
        self.assertEqual([c["exists"] for c in checks], ["yes", "no"])
        self.assertIsNone(scores)

    def test_the_fingerprint_plan_and_arms_must_be_the_ones_expected(self):
        cases = [
            ("fingerprint", common.result(cohort_a(), fingerprint="beef"), common.expect_row()),
            ("plan_and_arms", common.result(cohort_a(), seed=43), common.expect_row()),
            ("plan_and_arms", common.result(cohort_a(), arms=("released",)), common.expect_row()),
            ("format", common.result(cohort_a(), dataset="synthia_coli__drug_b"), common.expect_row()),
        ]
        for column, doc, exp in cases:
            code, checks, scores = run([(exp, doc)])
            self.assertEqual(code, 1, column)
            self.assertNotIn(checks[0][column], ("ok", "equal"), column)
            self.assertIsNone(scores)

    def test_a_lineage_result_at_another_seed_is_checked_against_the_seed_42_result(self):
        # The fingerprint hashes the settings, seed included, and the release
        # holds reports at seed 42 only: at seed 43 the result's key is
        # "f00d43", not the report's "f00d". It passes when its settings are
        # the report's but the seed and its cohort is the seed-42 result's.
        s42 = (common.expect_row(seed=42, plan="lineage-s42"), common.result(cohort_a(), seed=42))
        exp43 = common.expect_row(seed=43, plan="lineage-s43")
        code, checks, _ = run([s42, (exp43, common.result(cohort_a(), seed=43))])
        self.assertEqual(code, 0)
        self.assertEqual([c["fingerprint"] for c in checks], ["equal", "seed 43, cohort as at 42"])
        # Each of these stops the run: one MIC changed (fold 3's exact 0 read
        # as 1), a second setting changed (100 trees), the seed-42 key kept at
        # seed 43, and a seed-42 result that is not at the release's settings.
        other = cohort_a()
        other[6] = R(3, 1, "Exact", False, -0.5, 0)
        moved = common.result(cohort_a(), seed=43)
        moved["eval_config"]["forest"]["n_trees"] = 100
        cases = [
            [s42, (exp43, common.result(other, seed=43))],
            [s42, (exp43, moved)],
            [s42, (exp43, common.result(cohort_a(), seed=43, fingerprint="f00d"))],
            [(s42[0], common.result(cohort_a(), seed=42, fingerprint="beef")),
             (exp43, common.result(cohort_a(), seed=43))],
        ]
        for i, entries in enumerate(cases):
            code, checks, scores = run(entries)
            self.assertEqual((code, checks[1]["fingerprint"], checks[1]["ok"]), (1, "differs", "no"), i)
            self.assertIsNone(scores)

    def test_a_seed_43_cohort_must_hold_the_seed_42_rows_in_order(self):
        # grm learn writes the tested rows in the cohort's sample order at
        # every seed. Swapping two genomes' MIC, censor and call (fold 0's
        # exact 1, S and fold 4's above 6, R) keeps the counts and the rows
        # as a set, and is a different cohort.
        s42 = (common.expect_row(seed=42, plan="lineage-s42"), common.result(cohort_a(), seed=42))
        rows = cohort_a()
        rows[0], rows[9] = (R(0, 6, "Above", True, 1.4, 1), R(4, 1, "Exact", False, 5.0, 6))
        code, checks, scores = run([s42, (common.expect_row(seed=43, plan="lineage-s43"),
                                          common.result(rows, seed=43))])
        self.assertEqual((code, checks[1]["fingerprint"]), (1, "differs"))
        self.assertIsNone(scores)

    def test_a_rebuilt_breakpoint_passes_only_through_a_probe_equal_in_all_else(self):
        # The rebuilt table fits the dataset breakpoint on one more MIC pair,
        # and the fingerprint hashes the breakpoint: the run's key at the
        # release's settings is "beef", the report's "f00d". The probe, the
        # seed-42 fit with the release's breakpoint, keys "f00d" and is the
        # same result but for its wall time.
        exp42 = common.expect_row(seed=42, plan="lineage-s42")
        exp43 = common.expect_row(seed=43, plan="lineage-s43")
        exp47 = common.expect_row(split="study", plan="c6mic-s47", seed=47, arms="released")
        s42 = common.result(cohort_a(), seed=42, fingerprint="beef")
        s43 = common.result(cohort_a(), seed=43)
        plan = {"kind": "fold-file", "scheme": "c6mic-s47", "fold_file_sha256": "ab"}

        def study(key):
            return common.result(cohort_a(), arms=("released",), plan=plan, fingerprint=key)

        def probe(rows=None, key="f00d"):
            doc = common.result(rows or cohort_a(), seed=42, fingerprint=key)
            doc["seconds"] = 9.5
            return {"synthia_coli__drug_a.lineage-s42.mic.json": doc}

        entries = [(exp42, s42), (exp43, s43), (exp47, study("beef"))]
        code, checks, scores = run(entries, probe())
        self.assertEqual(code, 0)
        self.assertEqual([c["fingerprint"] for c in checks], [
            "equal with the release's breakpoint", "seed 43, cohort as at 42",
            "equal with the release's breakpoint"])
        self.assertEqual(len(scores[1]), 3 + 3 + 2)
        # Each of these stops the run: no probe; a probe that does not give
        # the report's key; a probe that differs in one MIC (fold 3's exact 0
        # read as 1); a study result keyed apart from the seed-42 result.
        other = cohort_a()
        other[6] = R(3, 1, "Exact", False, -0.5, 0)
        for i, (ents, prb) in enumerate([
            (entries, None),
            (entries, probe(key="cafe")),
            (entries, probe(rows=other)),
            ([(exp42, s42), (exp43, s43), (exp47, study("babe"))], probe()),
        ]):
            code, checks, scores = run(ents, prb)
            self.assertEqual(code, 1, i)
            self.assertIn("differs", [c["fingerprint"] for c in checks], i)
            self.assertIsNone(scores, i)

    def test_a_study_result_is_matched_by_its_scheme(self):
        exp = common.expect_row(split="study", plan="c6mic-s47", seed=47, arms="released")
        doc = common.result(cohort_a(), arms=("released",),
                            plan={"kind": "fold-file", "scheme": "c6mic-s47", "fold_file_sha256": "ab"})
        code, checks, scores = run([(exp, doc)])
        self.assertEqual((code, checks[0]["ok"]), (0, "yes"))
        rows = scores[1]
        self.assertEqual(sorted(r["arm"] for r in rows), ["median", "released"])
        self.assertEqual({r["draw"] for r in rows}, {"twin"})


class Output(unittest.TestCase):
    def test_one_aggregate_row_per_dataset_split_seed_and_arm(self):
        entries = []
        for seed in (42, 43, 44, 45, 46, 47):
            arms = ("released", "ordered") if seed < 45 else ("released",)
            exp = common.expect_row(seed=seed, plan=f"lineage-s{seed}", arms=",".join(arms))
            entries.append((exp, common.result(cohort_a(), seed=seed, arms=arms)))
        code, _, scores = run(entries)
        self.assertEqual(code, 0)
        head, rows = scores
        self.assertEqual(head, score.SCORE_COLUMNS)
        self.assertEqual(head.index("ea_exact"), head.index("arm") + 6)
        self.assertLess(head.index("ea_exact"), head.index("category_agreement"))
        keys = [(r["dataset"], r["split"], r["seed"], r["arm"]) for r in rows]
        self.assertEqual(len(keys), len(set(keys)))
        self.assertEqual(len(rows), 3 * 3 + 3 * 2)
        self.assertEqual({r["draw"] for r in rows if r["seed"] in ("42", "43", "44")}, {"draw"})
        self.assertEqual({r["ea_exact_reads"] for r in rows}, {"MIC values only"})
        self.assertEqual({r["category_reads"] for r in rows}, {"R and S labels and the training breakpoint"})
        for r in rows:
            self.assertFalse({"id", "genome", "genome_id", "tested", "samples"} & set(r))


if __name__ == "__main__":
    unittest.main()
