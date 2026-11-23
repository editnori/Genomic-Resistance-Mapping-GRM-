"""SYNTHETIC tests of scripts/c5-settings/drive.py against a fake grm (fake_grm.py): no fit runs, no data is read."""

import copy
import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from tests.c5_settings import common

D = common.load("drive")
F = common.fake()
LEARNERS = common.LEARNERS
# Distinct settings among lasso's 30 draws at the sampler seed 20261030.
LASSO_DISTINCT = 28


def run(top, doc=None, fold=0, trials=D.TRIALS, **extra):
    grm, root, folds, pick, work, out = common.setup(top, doc, **extra)
    path = D.drive(grm, root, common.SPECIES, folds, pick, fold, work, out, trials=trials)
    return json.loads(path.read_text()), root, work


def outer_test(doc, fold):
    return {g for g, o in zip(doc["genomes"], doc["outer"]) if o == fold}


def stem(doc, fold):
    return f"{doc['dataset']}.{doc['scheme']}-s{doc['seed']}.f{fold}"


class Drive(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.top = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def test_a_fold_runs_and_its_choices_hold_counts_only(self):
        doc = common.folds_doc()
        got, root, work = run(self.top, doc, fold=2)
        self.assertEqual(got["format"], "grm-c5-choices-1")
        self.assertEqual([r["learner"] for r in got["learners"]], LEARNERS)
        self.assertEqual(got["counts"], {"train": 48, "train_resistant": 24, "test": 12, "test_resistant": 6})
        for r in got["learners"]:
            self.assertEqual(r["trials"], 31)
            for arm in ("default", "search", "cutoff", "plain-best"):
                self.assertEqual(r["outer"][arm]["n"], 12)
        # Ten grm runs, one inner and one outer per learner, each timed.
        self.assertEqual([t["step"] for t in got["timing"]], [f"{s}-{l}" for s in ("inner", "outer") for l in LEARNERS])
        text = json.dumps(got)
        self.assertFalse(any(g in text for g in doc["genomes"]))

    def test_no_outer_test_genome_is_read_before_the_pick(self):
        doc = common.folds_doc()
        for fold in (0, 3):
            top = self.top / str(fold)
            _, root, _ = run(top, doc, fold=fold)
            log = common.calls_log(root)
            test = outer_test(doc, fold)
            before = [c for c in log if not c["pick_written"]]
            self.assertEqual(len(before), 5)
            for c in before:
                self.assertTrue(c["scheme"].startswith("c5i-"))
                self.assertFalse(test & set(c["tested"] + c["trained"]), c["scheme"])
                self.assertEqual(len(c["tested"]) + len(c["trained"]), 48)
            after = [c for c in log if c["pick_written"]]
            self.assertEqual(len(after), 5)
            for c in after:
                self.assertEqual(set(c["tested"]), test)

    def test_a_folds_file_that_splits_a_lineage_is_refused(self):
        doc = common.folds_doc()
        inner = doc["inner"][1]
        # Move one genome of a lineage to another inner fold.
        i = next(i for i, j in enumerate(inner) if j == 0)
        inner[i] = 1
        with self.assertRaisesRegex(D.DriveError, "splits a lineage"):
            run(self.top, doc, fold=1)

    def test_a_folds_file_whose_inner_folds_leave_out_a_training_genome_is_refused(self):
        doc = common.folds_doc()
        i = next(i for i, j in enumerate(doc["inner"][0]) if j is not None)
        doc["inner"][0][i] = None
        with self.assertRaisesRegex(D.DriveError, "cover exactly"):
            run(self.top, doc, fold=0)

    def test_the_cutoff_is_chosen_on_inner_scores_alone(self):
        # Inner scores sit 0.25 below the outer ones, so the inner best cutoff is well under 0.5 and the outer
        # fold would choose 0.5.
        doc = common.folds_doc()
        got, root, work = run(self.top, doc, fold=1, inner_shift=-0.25)
        here = work / stem(doc, 1)
        pick = json.loads((here / "pick.json").read_text())
        w = json.loads((root / "fake.json").read_text())
        test = sorted(outer_test(doc, 1))
        moved = 0
        for r in got["learners"]:
            l = r["learner"]
            inner0 = D.read_trials(here / f"inner-{l}", doc["dataset"], f"c5i-lineage-s42-f1", [0])[0]
            cut, _, _ = D.choose_cutoff(D.grid_calls(inner0, l))
            self.assertEqual(pick[l]["cutoff"], cut)
            self.assertEqual(r["cutoff"], cut)
            self.assertEqual(r["outer"]["cutoff"], F.calls(w, test, l, 0, cut, True), l)
            moved += cut < 0.5
        self.assertGreaterEqual(moved, 3)

    def test_a_tie_goes_to_the_lowest_trial_by_inner_scores_alone(self):
        # rf's trials 3 and 8 score alike on every inner fold and beat the release; on the outer fold 8 is far
        # better. The pick is 3, and the search arm serves 3's outer calls.
        doc = common.folds_doc()
        got, root, _ = run(self.top, doc, fold=0, quality={"rf:3": 0.97, "rf:8": 0.97},
                           alias={"rf:8": 3}, outer_quality={"rf:3": 0.6, "rf:8": 0.99})
        w = json.loads((root / "fake.json").read_text())
        rf = next(r for r in got["learners"] if r["learner"] == "rf")
        self.assertEqual(rf["tied_at_best"], 2)
        self.assertEqual((rf["picked_trial"], rf["plain_best_trial"]), (3, 3))
        self.assertEqual(rf["outer"]["search"], F.calls(w, sorted(outer_test(doc, 0)), "rf", 3, 0.5, True))
        trained = [c for c in common.calls_log(root) if c["scheme"].startswith("c5o-")]
        self.assertTrue(all(8 not in c["trials"] for c in trained))

    def test_the_default_arm_is_the_defaults_own_fit(self):
        # gbt's trial 5 beats the release on the inner folds, so the search arm differs from the default; the
        # outer refit of trial 0 is the release's setting, and its calls are the default's.
        doc = common.folds_doc()
        got, root, work = run(self.top, doc, fold=4, quality={"gbt:5": 0.97}, outer_quality={"gbt:5": 0.97})
        pick = json.loads((self.top / "pick.json").read_text())
        here = work / stem(doc, 4)
        gbt = next(r for r in got["learners"] if r["learner"] == "gbt")
        self.assertEqual(gbt["picked_trial"], 5)
        self.assertNotEqual(gbt["outer"]["search"], gbt["outer"]["default"])
        for li, r in enumerate(got["learners"]):
            want = {k: pick["folds"][4]["outer"][li][k] for k in D.ZERO}
            self.assertEqual(r["outer"]["default"], want, r["learner"])
            s = json.loads((here / f"outer-{r['learner']}.settings.json").read_text())
            self.assertEqual(s["trials"][0], {"id": 0, "learners": {r["learner"]: "release"},
                                              "cutoffs": {r["learner"]: r["cutoff"]}})

    def test_a_default_that_is_not_the_release_refit_stops_the_fold(self):
        grm, root, folds, pick, work, out = common.setup(self.top)
        p = json.loads(pick.read_text())
        p["folds"][2]["outer"][1]["tp"] -= 1
        p["folds"][2]["outer"][1]["fn"] += 1
        pick.write_text(json.dumps(p))
        with self.assertRaisesRegex(D.DriveError, "cart: the release's refit on fold 2 is not the default's fit"):
            D.drive(grm, root, common.SPECIES, folds, pick, 2, work, out)

    def test_a_restart_runs_only_the_missing_steps(self):
        doc = common.folds_doc()
        got, root, work = run(self.top, doc, fold=3)
        self.assertEqual(len(common.calls_log(root)), 10)
        out = self.top / "out" / f"{stem(doc, 3)}.choices.json"
        out.unlink()
        (work / stem(doc, 3) / "outer-lasso.done").unlink()
        again = D.drive(self.top / "grm", root, common.SPECIES, self.top / "folds.json", self.top / "pick.json", 3,
                        work, self.top / "out")
        log = common.calls_log(root)
        self.assertEqual(len(log), 11)
        self.assertTrue(log[-1]["scheme"].startswith("c5o-"))
        self.assertEqual(json.loads(again.read_text())["learners"], got["learners"])

    def test_a_choices_file_that_would_name_a_genome_is_not_written(self):
        doc = common.folds_doc(dataset=f"synthia_coli__drug-{common.genome(0)}")
        with self.assertRaisesRegex(D.DriveError, "would name a genome"):
            run(self.top, doc, fold=0, trials=2)
        self.assertFalse(any((self.top / "out").glob("*")) if (self.top / "out").exists() else False)


    def test_the_plain_best_is_refitted_beside_only_when_it_is_neither_the_release_nor_the_pick(self):
        # SYNTHETIC: the rule's answers are set per learner, in drive order. scm's plain best is its pick (4);
        # cart's is trial 0; rf's is trial 6 with the pick at the release (its lead within gamma); gbt's release
        # failed inside and its plain best is 9; every lasso trial failed.
        forced = {"scm": (4, 4, "picked"), "cart": (0, 0, "picked"), "rf": (0, 6, "picked"),
                  "gbt": (0, 9, "release failed inside"), "lasso": (0, None, "release failed inside")}
        order, real = iter(LEARNERS), D.pick

        def rule(means, gamma=D.GAMMA):
            picked, plain, status = forced[next(order)]
            return dict(real(means, gamma), picked=picked, plain_best=plain, status=status)

        doc = common.folds_doc()
        with mock.patch.object(D, "pick", rule):
            got, root, _ = run(self.top, doc, fold=1)
        w = json.loads((root / "fake.json").read_text())
        test = sorted(outer_test(doc, 1))
        by = {r["learner"]: r for r in got["learners"]}
        outer = [c["trials"] for c in common.calls_log(root) if c["scheme"].startswith("c5o-")]
        self.assertEqual(outer, [[0, 4], [0], [0, 6], [0, 9], [0]])
        self.assertEqual([by[l]["plain_best_refitted"] for l in LEARNERS], [False, False, True, True, False])
        self.assertEqual(by["scm"]["outer"]["plain-best"], by["scm"]["outer"]["search"])
        self.assertEqual(by["cart"]["outer"]["plain-best"], by["cart"]["outer"]["default"])
        self.assertEqual(by["rf"]["outer"]["plain-best"], F.calls(w, test, "rf", 6, 0.5, True))
        self.assertEqual(by["rf"]["outer"]["search"], by["rf"]["outer"]["default"])
        self.assertEqual(by["gbt"]["outer"]["plain-best"], F.calls(w, test, "gbt", 9, 0.5, True))
        self.assertEqual(by["lasso"]["outer"]["plain-best"], by["lasso"]["outer"]["default"])
        self.assertIsNotNone(by["lasso"]["plain_best_note"])
        self.assertTrue(all(by[l]["plain_best_note"] is None for l in LEARNERS[:4]))


class Rule(unittest.TestCase):
    def test_an_inner_fold_with_no_figure_drops_the_trial_and_is_never_read_as_zero(self):
        # SYNTHETIC: an inner test fold that lacks a class has no balanced accuracy. Its trial has no mean, as a
        # failed fit has none; read as 0, trial 1 would keep a mean of 0.633..., and a trial 0 with such a fold
        # would anchor the rule at 0.6 in place of the release failed inside.
        def doc(bas, status="ok"):
            return {"learners": [{"learner": "rf", "status": status, "fold_balanced_accuracy": bas}]}

        docs = {0: doc([0.7, 0.7, 0.7]), 1: doc([0.95, 0.95, None]), 2: doc([0.8, 0.8, 0.8], "failed")}
        self.assertEqual(D.history(docs, "rf"), {0: 0.7, 1: None, 2: None})
        release = D.history({0: doc([0.9, 0.9, None]), 1: doc([0.75, 0.75, 0.75])}, "rf")
        self.assertEqual(release, {0: None, 1: 0.75})
        got = D.pick(release)
        self.assertEqual((got["picked"], got["plain_best"], got["status"]), (0, 1, "release failed inside"))

    def test_the_plain_best_is_the_best_inner_mean_ties_lowest_also_when_the_release_failed(self):
        self.assertEqual(D.pick({0: 0.70, 3: 0.73, 5: 0.73})["plain_best"], 3)
        self.assertEqual(D.pick({0: 0.70, 3: 0.73, 5: 0.73})["picked"], 0)
        got = D.pick({0: None, 2: 0.9, 1: 0.9, 4: 0.5})
        self.assertEqual((got["picked"], got["plain_best"], got["status"]), (0, 1, "release failed inside"))
        self.assertIsNone(D.pick({0: None, 1: None})["plain_best"])

    def test_the_release_is_kept_unless_the_best_beats_it_by_more_than_gamma(self):
        self.assertEqual(D.pick({0: 0.70, 1: 0.74, 2: 0.73})["picked"], 0)
        self.assertEqual(D.pick({0: 0.70, 1: 0.740000000001, 2: 0.73})["picked"], 1)
        # 12 decimals: a lead of gamma plus 1e-13 is gamma.
        self.assertEqual(D.pick({0: 0.70, 1: 0.7400000000001})["picked"], 0)
        got = D.pick({0: 0.70, 1: 0.80, 2: 0.80, 3: None})
        self.assertEqual((got["picked"], got["tied"], got["failed"], got["better"]), (1, [1, 2], 1, 2))

    def test_a_failed_release_keeps_the_release(self):
        got = D.pick({0: None, 1: 0.9})
        self.assertEqual((got["picked"], got["status"]), (0, "release failed inside"))

    def test_the_cutoff_tie_goes_nearest_one_half_then_lower(self):
        def c(tp, tn):
            return {"n": 20, "tp": tp, "fn": 10 - tp, "tn": tn, "fp": 10 - tn, "loss_sum": 1.0}
        flat = {g: c(5, 5) for g in D.GRID}
        self.assertEqual(D.choose_cutoff(flat)[0], 0.5)
        two = {**flat, 0.45: c(8, 8), 0.55: c(8, 8)}
        self.assertEqual(D.choose_cutoff(two)[0], 0.45)
        far = {**flat, 0.2: c(9, 9), 0.7: c(9, 9)}
        self.assertEqual(D.choose_cutoff(far)[0], 0.7)
        self.assertEqual(D.choose_cutoff(None), (0.5, None, None))

    def test_the_thirty_settings_are_the_protocols_draws_the_same_everywhere(self):
        from hyperopt import protocol, search
        p = {"families": protocol.FAMILIES, "seeds": {"sampler": 20261030}}
        for l in LEARNERS:
            self.assertEqual(D.setting(l, 0), "release")
            drawn = [D.setting(l, t) for t in range(1, 31)]
            self.assertEqual(drawn, [search.sample_config(p, l, t) for t in range(1, 31)])
            # lasso's space is one whole number, n_lambdas in 10 to 100, so its 30 draws repeat some values.
            distinct = len({json.dumps(s, sort_keys=True) for s in drawn})
            self.assertEqual(distinct, 30 if l != "lasso" else LASSO_DISTINCT, l)
            for s in drawn:
                self.assertEqual(sorted(s), sorted(protocol.FAMILIES[l]["space"]))

    def test_the_inner_seed_is_fixed_by_the_draw_and_fold(self):
        self.assertEqual(D.inner_seed(42, 0), D.inner_seed(42, 0))
        self.assertNotEqual(D.inner_seed(42, 0), D.inner_seed(42, 1))
        self.assertLess(D.inner_seed(42, 3), 2**63)
        self.assertEqual(D.fold_seed(42, 3), 42 + 303)


if __name__ == "__main__":
    unittest.main()
