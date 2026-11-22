"""scripts/c9-pick/score.py: the checks, the served figure, the default, ties, the count (docs/served-pick.md).

Every result, report and dataset name is SYNTHETIC.
"""

import contextlib
import csv
import io
import math
import tempfile
import unittest
from pathlib import Path

from tests.c9_pick import common

score = common.load("score")


def run(rows, docs, releases, *options):
    """score.py's first form on a SYNTHETIC run: (exit code, stdout, checks, draws)."""
    with tempfile.TemporaryDirectory() as tmp:
        exp, picks, rel, out = common.write_run(tmp, docs, rows, releases)
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            rc = score.main([str(exp), str(picks), str(rel), str(out), *options])
        return rc, buf.getvalue(), table(out / "checks.tsv"), table(out / "draws.tsv"), table(out / "datasets.tsv")


def table(path):
    if not path.exists():
        return []
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def failed(checks):
    return sorted({c["check"] for c in checks if c["passed"] != "1"})


class Checks(unittest.TestCase):
    def test_a_whole_run_passes_and_prints_its_check_counts_only(self):
        rows, docs, rels = common.full_run(datasets=("synthia_coli__drug_a", "synthia_coli__drug_b"),
                                           study=("synthia_coli__drug_b",))
        rc, out, checks, draws, datasets = run(rows, docs, rels)
        self.assertEqual(rc, 0, failed(checks))
        self.assertEqual(out.count("\n"), 1)
        self.assertTrue(out.startswith(f"checks: {len(checks)} of {len(checks)} pass; results 13 of 13;"))
        self.assertNotIn("drug", out)
        # Two datasets: 3 draws with both arms and 3 twins each; one study run with both arms.
        self.assertEqual(len(draws), 2 * (3 * 2 + 3) + 2)
        self.assertEqual({(r["split"], r["arm"]) for r in datasets},
                         {("lineage", "default"), ("lineage", "inner"), ("lineage", "twin"), ("study", "default"),
                          ("study", "inner")})
        self.assertEqual({r["species"] for r in datasets}, {"Synthia coli"})
        self.assertEqual({r["antibiotic"] for r in datasets}, {"drug a", "drug b"})

    def assert_caught(self, change, check, **kw):
        rows, docs, rels = common.full_run(change=change, **kw)
        rc, _, checks, draws, _ = run(rows, docs, rels)
        self.assertEqual(rc, 1)
        self.assertIn(check, " | ".join(failed(checks)))
        self.assertEqual(draws, [], "nothing is scored when a check fails")

    def test_an_inner_fold_holding_an_outer_test_genome_fails(self):
        # Planted fault 1 as the file shows it: the inner folds test more genomes than the outer fold trains on.
        def change(r, d):
            if r[6] == "draw":
                d["folds"][1]["inner"]["inner_folds"][0]["test"] += 5
        self.assert_caught(change, "fold 1: the inner folds cover the training genomes exactly")

    def test_a_pick_its_inner_figures_do_not_give_fails(self):
        def change(r, d):
            if r[4] == 43 and r[6] == "draw":
                d["folds"][2]["inner"]["pick"] = "gbt"
        self.assert_caught(change, "fold 2: the pick and tied set follow from the inner figures")

    def test_a_tie_sent_against_the_files_rule_fails(self):
        def tie(f, li):
            return [(11, 16)] * 3 if li in (1, 2) else [(8, 14)] * 3

        def change(r, d):
            if r[6] == "draw" and r[4] == 44:
                fresh = common.doc(dataset=r[0], seed=44, inner_calls=tie)
                d["folds"] = fresh["folds"]
                for fold in d["folds"]:
                    fold["inner"]["pick"] = "cart"
        self.assert_caught(change, "the pick and tied set follow from the inner figures")

    def test_outer_calls_that_do_not_give_the_split_figure_fail(self):
        def change(r, d):
            if r[6] == "twin":
                d["folds"][0]["outer"][2]["tp"] += 1
                d["folds"][0]["outer"][2]["fn"] -= 1
        self.assert_caught(change, "rf: the outer calls give the split's balanced accuracy")

    def test_a_result_of_another_seed_or_inner_count_fails(self):
        def change(r, d):
            if r[4] == 43:
                d["seed"] = 42
        self.assert_caught(change, "dataset, scheme, seed and inner folds")

    def test_a_twin_with_an_inner_fit_fails(self):
        def change(r, d):
            if r[6] == "twin":
                d["folds"][0]["inner"]["status"] = "picked"
        self.assert_caught(change, "fold 0: no inner fit")

    def test_a_seed_42_figure_other_than_the_release_reports_fails(self):
        rows, docs, rels = common.full_run()
        rels[common.DATASET]["splits"][1]["learners"][3]["metrics"]["balanced_accuracy"] += 1e-9
        rc, _, checks, _, _ = run(rows, docs, rels)
        self.assertEqual(rc, 1)
        self.assertEqual(failed(checks), ["seed 42: each learner's lineage figure is the release report's"])

    def test_a_default_read_from_the_lineage_split_fails_at_seed_42(self):
        # Planted fault 5 as the file shows it: the default's figure is not the release's random-split one.
        def change(r, d):
            if r[6] == "draw" and r[4] == 42:
                d["default"]["random_balanced_accuracy"] = d["learners"][2]["balanced_accuracy"]
        rows, docs, rels = common.full_run()
        rows2, docs2, _ = common.full_run(change=change)
        rc, _, checks, _, _ = run(rows2, docs2, rels)
        self.assertEqual(rc, 1)
        self.assertEqual(failed(checks), ["seed 42: the default's random-split figure is the release report's"])

    def test_a_default_other_than_the_learner_the_release_serves_fails_at_seed_42(self):
        # A default naming another learner with that learner's own random figure and tied set passes the figure
        # checks; only the release's index, which names the learner it serves, shows it.
        def change(r, d):
            if r[6] == "draw":
                d["default"] = {"learner": "gbt", "random_balanced_accuracy": 0.7, "tied": ["gbt"]}
        rows, docs, rels = common.full_run()
        rows2, docs2, _ = common.full_run(change=change)
        rc, _, checks, _, _ = run(rows2, docs2, rels)
        self.assertEqual(rc, 1)
        self.assertEqual(failed(checks), ["seed 42: the default's learner is the one the release serves"])

    def test_a_default_that_breaks_a_random_split_tie_to_the_earlier_learner_fails(self):
        def change(r, d):
            if r[6] == "draw":
                d["default"] = {"learner": "rf", "random_balanced_accuracy": 0.8, "tied": ["rf", "gbt"]}
        rows, docs, rels = common.full_run(change=change)
        rc, _, checks, _, _ = run(rows, docs, rels)
        self.assertEqual(rc, 1)
        self.assertEqual(failed(checks), ["the default's learner is the later of its tied set, as the release breaks a tie"])

    def test_a_missing_result_fails(self):
        rows, docs, rels = common.full_run()
        del docs[rows[1][7]]
        rc, out, checks, _, _ = run(rows, docs, rels)
        self.assertEqual(rc, 1)
        self.assertEqual(failed(checks), ["exists"])
        self.assertIn("results 6 of 7", out)

    def test_an_unavailable_split_is_counted_not_scored(self):
        def change(r, d):
            if r[0] == common.DATASET and r[2] == "lineage" and r[4] != 42:
                d.update(unavailable_reason="too few lineages", folds=[], tested=0, tested_resistant=0)
        rows, docs, rels = common.full_run(change=change)
        rc, out, checks, draws, datasets = run(rows, docs, rels)
        self.assertEqual(rc, 0, failed(checks))
        self.assertIn("unavailable in 5", out)
        self.assertEqual(sum(r["why"] == "split unavailable" for r in draws), 2 * 2 + 3)
        lineage = {r["arm"]: r for r in datasets if r["split"] == "lineage"}
        self.assertEqual({a: r["scored"] for a, r in lineage.items()}, {"default": "0", "inner": "0", "twin": "0"})


def picked(f_pick):
    """SYNTHETIC inner calls whose best learner on fold f is f_pick[f]."""
    def inner(f, li):
        good = common.LEARNERS.index(f_pick[f]) == li
        return [(12 if good else 9, 18 if good else 15)] * 3
    return inner


class Served(unittest.TestCase):
    def test_the_served_figure_is_each_folds_pick_on_its_outer_test_fold(self):
        # Planted fault 6 would score the pick on its inner folds: the inner calls here differ from every outer row.
        d = common.doc(inner_calls=picked(["scm", "gbt", "lasso"]))
        c, counts = score.served(d, "rf", "later", "default")
        want = {k: 0 for k in ("n", "tp", "fp", "tn", "fn", "loss_sum")}
        for f, l in enumerate(["scm", "gbt", "lasso"]):
            tp, tn = common.outer_default(f, common.LEARNERS.index(l))
            for k, v in common.calls(tp, tn, common.OUTER_POS, common.OUTER_NEG).items():
                want[k] += v
        self.assertEqual({k: c[k] for k in ("n", "tp", "fp", "tn", "fn")}, {k: want[k] for k in ("n", "tp", "fp", "tn", "fn")})
        self.assertEqual(c["n"], 150)
        self.assertEqual((counts["picks_scm"], counts["picks_gbt"], counts["picks_lasso"], counts["picks_rf"]), (1, 1, 1, 0))
        self.assertEqual(counts["same_as_default"], 0)
        m = score.measures(c)
        self.assertEqual(m["balanced_accuracy"], common.balanced(want))

    def test_a_tie_is_picked_again_by_the_scorers_rule(self):
        def tie(f, li):
            return [(11, 16)] * 3 if li in (1, 2) else [(8, 14)] * 3
        d = common.doc(inner_calls=tie)
        rows = {}
        for rule, default in (("later", "gbt"), ("earlier", "gbt"), ("default", "cart"), ("default", "gbt")):
            _, counts = score.served(d, default, rule, "default")
            rows[(rule, default)] = [l for l in common.LEARNERS if counts["picks_" + l]]
            self.assertEqual(counts["ties"], 3)
        self.assertEqual(rows, {("later", "gbt"): ["rf"], ("earlier", "gbt"): ["cart"], ("default", "cart"): ["cart"],
                                ("default", "gbt"): ["rf"]})

    def test_scores_equal_at_12_decimals_tie_and_others_do_not(self):
        self.assertEqual(score.pick({"scm": 0.8, "cart": 0.8 + 1e-15, "rf": 0.7}, "later"), ("cart", ["scm", "cart"]))
        self.assertEqual(score.pick({"scm": 0.8, "cart": 0.8 + 1e-15, "rf": 0.7}, "earlier"), ("scm", ["scm", "cart"]))
        self.assertEqual(score.pick({"scm": 0.8, "cart": 0.8000001}, "earlier"), ("cart", ["cart"]))
        self.assertEqual(score.pick({"scm": None, "cart": None}, "later"), (None, []))

    def test_a_fold_with_no_pick_serves_the_default_or_drops_the_draw(self):
        d = common.doc(no_pick=(1,))
        c, counts = score.served(d, "rf", "later", "default")
        self.assertEqual((counts["no_pick"], counts["picks_rf"], counts["same_as_default"]), (1, 3, 2))
        self.assertEqual(c["n"], 150)
        c, counts = score.served(d, "rf", "later", "drop")
        self.assertIsNone(c)

    def test_same_as_default_counts_folds_whose_inner_pick_is_the_defaults_learner(self):
        # Design section 3: a fold with no pick serves the default's learner, counted under no_pick and not here.
        d = common.doc(no_pick=(0, 1))
        _, counts = score.served(d, "rf", "later", "default")
        self.assertEqual((counts["no_pick"], counts["picks_rf"], counts["same_as_default"]), (2, 3, 1))
        _, counts = score.served(d, "gbt", "later", "default")
        self.assertEqual((counts["no_pick"], counts["picks_gbt"], counts["same_as_default"]), (2, 2, 0))

    def test_a_failed_inner_fit_does_not_compete(self):
        def inner(f, li):
            return None if li == 2 else common.inner_default(f, li)
        d = common.doc(inner_calls=inner)
        _, counts = score.served(d, "rf", "later", "default")
        # rf failed; cart and gbt tie behind it, and the tie goes to the later.
        self.assertEqual((counts["picks_rf"], counts["picks_gbt"], counts["ties"]), (0, 3, 3))

    def test_the_default_is_each_draws_own_pick_or_seed_42s(self):
        def change(r, d):
            if r[4] == 43 and r[6] == "draw":
                d["default"] = {"learner": "gbt", "random_balanced_accuracy": 0.81, "tied": ["gbt"]}
        rows, docs, rels = common.full_run(change=change)
        for option, at43 in (("per-draw", "gbt"), ("seed42", "rf")):
            rc, _, checks, draws, _ = run(rows, docs, rels, "--default", option)
            self.assertEqual(rc, 0, failed(checks))
            got = {(r["split"], r["seed"], r["arm"]): r["learner"] for r in draws if r["arm"] != "inner"}
            self.assertEqual(got[("lineage", "43", "default")], at43)
            self.assertEqual(got[("lineage", "42", "default")], "rf")
            # The study split's default is always the seed-42 lineage pick.
            self.assertEqual(got[("study", "42", "default")], "rf")

    def test_the_default_serves_its_learner_on_every_outer_fold(self):
        d = common.doc()
        c = score.fixed(d, "gbt")
        tp = sum(common.outer_default(f, 3)[0] for f in range(3))
        self.assertEqual((c["n"], c["tp"]), (150, tp))

    def test_an_option_outside_its_choices_is_refused(self):
        with self.assertRaises(SystemExit):
            score.options(["--tie", "first"])
        self.assertEqual(score.options(["a", "--no-pick", "drop"]), ({"--tie": "later", "--no-pick": "drop",
                                                                      "--default": "per-draw"}, ["a"]))


class Measures(unittest.TestCase):
    def test_rates_with_their_wilson_intervals_and_log_loss(self):
        m = score.measures({"n": 50, "tp": 15, "fp": 3, "tn": 27, "fn": 5, "loss_sum": 10.0})
        self.assertEqual((m["recall"], m["specificity"], m["vme"], m["me"]), (0.75, 0.9, 0.25, 0.1))
        self.assertAlmostEqual(m["balanced_accuracy"], 0.825, places=15)
        self.assertEqual(m["log_loss"], 0.2)
        lo, hi = score.wilson(5, 10)
        self.assertAlmostEqual(lo, 0.2365930, places=6)
        self.assertAlmostEqual(hi, 0.7634070, places=6)
        self.assertEqual(score.wilson(0, 10)[0], 0.0)
        self.assertAlmostEqual(score.wilson(10, 10)[1], 1.0, places=12)
        self.assertEqual(score.wilson(0, 0), (None, None))

    def test_a_dataset_is_the_mean_of_its_draws_and_unscored_if_any_draw_is(self):
        base = {"dataset": "synthia_coli__drug_a", "species": "Synthia coli", "antibiotic": "drug a", "split": "lineage",
                "arm": "default", "n": 150, "tp": 1, "fp": 1, "tn": 1, "fn": 1}
        draws = [dict(base, seed=s, scored=1, balanced_accuracy=v, recall=v, specificity=v, vme=v, me=v, log_loss=v)
                 for s, v in ((42, 0.7), (43, 0.8), (44, 0.75))]
        row = score.dataset_rows(draws)[0]
        self.assertEqual((row["draws"], row["scored"], row["balanced_accuracy"]), (3, 1, 0.75))
        draws[1] = dict(base, seed=43, scored=0)
        row = score.dataset_rows(draws)[0]
        self.assertEqual((row["scored"], row["balanced_accuracy"], row["n"]), (0, None, None))


def cell(d, split, arm, sp, **m):
    row = {"dataset": d, "species_slug": sp, "split": split, "arm": arm, "scored": "1"}
    for k in score.MEASURES:
        row[k] = repr(m.get(k, 0.5))
    return row


def pairs(values, split="lineage", arm="inner", measure="balanced_accuracy", species=None, default=True):
    """SYNTHETIC datasets.tsv rows: per dataset, the default at 0.70 (unless default is False) and the arm at 0.70 + value."""
    rows = []
    for i, v in enumerate(values):
        sp = species[i] if species else f"sp{i % 4}"
        if default:
            rows.append(cell(f"{sp}__d{i}", split, "default", sp, **{measure: 0.7}))
        rows.append(cell(f"{sp}__d{i}", split, arm, sp, **{measure: 0.7 + v}))
    return rows


def at(rows, split, pair, measure):
    return [r for r in score.count(rows) if (r["split"], r["pair"], r["measure"]) == (split, pair, measure)][0]


class Count(unittest.TestCase):
    def test_the_margin_is_counted_on_differences_rounded_to_12_decimals(self):
        # 0.72 - 0.70 is 0.020000000000000018 in floating point; rounded it is the margin, so the same.
        rows = pairs([0.02, -0.02, 0.0200000001, -0.0200000001, 0.0])
        r = at(rows, "lineage", "inner-default", "balanced_accuracy")
        self.assertEqual((r["better"], r["worse"], r["same"], r["entered"]), (1, 1, 3, 5))

    def test_lower_is_better_for_log_loss_and_both_error_rates(self):
        for m in ("log_loss", "vme", "me"):
            r = at(pairs([-0.05, -0.05, 0.05], measure=m), "lineage", "inner-default", m)
            self.assertEqual((r["better"], r["worse"], r["orientation"]), (2, 1, "lower"))
        r = at(pairs([-0.05, -0.05, 0.05], measure="recall"), "lineage", "inner-default", "recall")
        self.assertEqual((r["better"], r["worse"]), (1, 2))

    def test_the_median_and_its_interval_are_the_arm_less_the_default_for_every_measure(self):
        # A lower log loss or error rate is better, but its difference stays the arm's value less the default's.
        for m in score.MEASURES:
            r = at(pairs([-0.05] * 6, measure=m), "lineage", "inner-default", m)
            self.assertEqual((r["median_difference"], r["interval_lo"], r["interval_hi"]), (-0.05, -0.05, -0.05), m)
            self.assertEqual((r["better"], r["worse"]), (6, 0) if m in score.LOWER_IS_BETTER else (0, 6), m)

    def test_the_sign_test_is_one_sided_and_holm_runs_over_the_inner_pick_and_the_twin(self):
        self.assertEqual(score.sign_p(5, 0), 1 / 32)
        self.assertEqual(score.sign_p(3, 2), 0.5)
        self.assertEqual(score.sign_p(0, 0), 1.0)
        self.assertEqual(score.holm([0.01, 0.04]), [0.02, 0.04])
        self.assertEqual(score.holm([0.04, 0.01]), [0.04, 0.02])
        self.assertEqual(score.holm([0.03, 0.02]), [0.04, 0.04])
        rows = pairs([0.05] * 6) + pairs([0.05] * 5 + [-0.05], arm="twin", default=False)
        inner = at(rows, "lineage", "inner-default", "balanced_accuracy")
        twin = at(rows, "lineage", "twin-default", "balanced_accuracy")
        self.assertEqual((inner["sign_p"], twin["sign_p"]), (1 / 64, 7 / 64))
        self.assertEqual((inner["holm_p"], twin["holm_p"]), (2 / 64, 7 / 64))
        self.assertIsNone(at(rows, "lineage", "inner-default", "recall")["holm_p"])

    def test_the_species_interval_is_seeded_by_its_row_and_read_at_the_101st_and_3900th_values(self):
        self.assertEqual((score.RESAMPLES, score.LO_INDEX, score.HI_INDEX), (4000, 100, 3899))
        self.assertNotEqual(score.row_seed("lineage", "inner-default", "recall"),
                            score.row_seed("lineage", "inner-default", "balanced_accuracy"))
        flat = at(pairs([0.05] * 8), "lineage", "inner-default", "balanced_accuracy")
        self.assertEqual((flat["median_difference"], flat["interval_lo"], flat["interval_hi"]), (0.05, 0.05, 0.05))
        spread = pairs([0.01, 0.03, -0.02, 0.06, 0.0, 0.04, 0.02, -0.01], species=["a", "a", "b", "b", "c", "c", "d", "d"])
        a, b = at(spread, "lineage", "inner-default", "balanced_accuracy"), at(spread, "lineage", "inner-default", "balanced_accuracy")
        self.assertEqual((a["interval_lo"], a["interval_hi"]), (b["interval_lo"], b["interval_hi"]))
        self.assertLessEqual(a["interval_lo"], a["median_difference"])
        self.assertLessEqual(a["median_difference"], a["interval_hi"])
        self.assertEqual(a["species"], 4)

    def test_an_unscored_dataset_is_counted_apart(self):
        rows = pairs([0.05, 0.05])
        rows[1]["scored"] = "0"  # the arm row of the first dataset
        r = at(rows, "lineage", "inner-default", "balanced_accuracy")
        self.assertEqual((r["datasets"], r["entered"], r["not_scored"]), (2, 1, 1))

    def test_the_conditions_read_the_count(self):
        good = pairs([0.05] * 8) + pairs([0.0] * 8, arm="twin", default=False) + pairs([0.05] * 3, split="study")
        rows = {r["condition"]: r["met"] for r in score.conditions(score.count(good))}
        self.assertEqual(rows, {1: 1, 2: 1, 3: 1, 4: 1, 5: 1, 6: ""})
        # The twin moving as much as the arm: condition 3 fails.
        bad_twin = pairs([0.05] * 8) + pairs([0.05] * 8, arm="twin", default=False) + pairs([0.05] * 3, split="study")
        self.assertEqual(score.conditions(score.count(bad_twin))[2]["met"], 0)
        # A study median of the other sign: condition 5 fails.
        other = pairs([0.05] * 8) + pairs([0.0] * 8, arm="twin", default=False) + pairs([-0.05] * 3, split="study")
        self.assertEqual(score.conditions(score.count(other))[4]["met"], 0)

    def test_the_count_form_writes_its_tables_and_refuses_a_row_listed_twice(self):
        rows = pairs([0.05] * 4) + pairs([0.0] * 4, arm="twin", default=False) + pairs([0.01], split="study")
        cols = ["dataset", "species_slug", "split", "arm", "scored"] + score.MEASURES
        with tempfile.TemporaryDirectory() as tmp:
            path = f"{tmp}/datasets.tsv"
            with open(path, "w") as f:
                f.write("\t".join(cols) + "\n")
                for r in rows:
                    f.write("\t".join(r[c] for c in cols) + "\n")
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                self.assertEqual(score.main(["--count", f"{tmp}/out", path]), 0)
            self.assertEqual(buf.getvalue(), f"count: {len(rows)} dataset rows from 1 files\n")
            counted = table(Path(tmp) / "out" / "count.tsv")
            self.assertEqual(len(counted), 3 * len(score.MEASURES))
            self.assertTrue(all(not math.isnan(float(r["sign_p"])) for r in counted))
            with self.assertRaises(SystemExit):
                with contextlib.redirect_stdout(io.StringIO()):
                    score.main(["--count", f"{tmp}/out", path, path])


if __name__ == "__main__":
    unittest.main()
