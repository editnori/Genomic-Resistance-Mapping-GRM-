"""The outside-learner runner (scripts/c7-outside/outside.py) on SYNTHETIC exports, unitig builds and fold files.

The arm tests import XGBoost 3.0.5, LightGBM 4.6.0 and scikit-learn 1.7.2, the pins of
scripts/c7-outside/requirements.txt; without them they fail, and none is skipped.
"""

import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np

from tests.c7_outside import common as C
from tests.c7_outside import plain

O = C.O


class Temp(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.top = Path(self._tmp.name)

    def tearDown(self):
        self._tmp.cleanup()

    def export(self, name="a", **kw):
        d, x, y = C.write_export(self.top / name, **kw)
        return d


class Rules(unittest.TestCase):
    def test_grouped_folds_keep_each_group_whole_and_balance_the_classes(self):
        rng = np.random.default_rng(3)
        groups = rng.integers(0, 20, 200).tolist()
        y = rng.random(200) < 0.4
        fold = O.grouped_folds(y, groups, 5, np.random.default_rng(1))
        for g in set(groups):
            self.assertEqual(len({int(fold[i]) for i in range(200) if groups[i] == g}), 1)
        self.assertEqual(sorted(set(fold.tolist())), [0, 1, 2, 3, 4])

    def test_stratified_folds_deal_each_class_round_the_folds(self):
        y = np.array([True] * 23 + [False] * 37)
        fold = O.stratified_folds(y, 3, np.random.default_rng(2))
        for cls in (True, False):
            counts = np.bincount(fold[y == cls], minlength=3)
            self.assertLessEqual(int(counts.max() - counts.min()), 1)

    def test_a_grouped_holdout_is_the_fold_of_seven_closest_to_15_percent(self):
        rng = np.random.default_rng(4)
        groups = rng.integers(0, 30, 300).tolist()
        y = rng.random(300) < 0.5
        train = np.ones(300, dtype=bool)
        train[:40] = False
        fit, hold = O.tuning_holdout(train, y, groups, 0.15, np.random.default_rng(9))
        self.assertFalse(np.any(hold & ~train))
        self.assertTrue(np.array_equal(fit | hold, train))
        folds = O.training_folds(y, train, groups, 7, np.random.default_rng(9))
        sizes = [int(np.sum(train & (folds == f))) for f in range(7)]
        self.assertIn(int(hold.sum()), sizes)
        for g in set(groups):
            members = [i for i in range(300) if groups[i] == g and train[i]]
            self.assertEqual(len({bool(hold[i]) for i in members}), min(1, len(members)))

    def test_an_ungrouped_holdout_takes_15_percent_of_each_class(self):
        y = np.array([True] * 40 + [False] * 60)
        train = np.ones(100, dtype=bool)
        _, hold = O.tuning_holdout(train, y, None, 0.15, np.random.default_rng(1))
        self.assertEqual(int(np.sum(hold & y)), 6)
        self.assertEqual(int(np.sum(hold & ~y)), 9)

    def test_rust_round_takes_halves_away_from_zero(self):
        self.assertEqual([O.rust_round(v) for v in (0.5, 1.5, 2.5, 6.666)], [1, 2, 3, 7])


class Screen(unittest.TestCase):
    def test_the_screen_is_the_exports_rule_written_plainly(self):
        rng = np.random.default_rng(12)
        n = 90
        x = rng.random((300, n)) < rng.uniform(0.02, 0.6, (300, 1))
        x[200:230] = x[100:130]  # repeated patterns
        x[230:240] = x[0]  # one pattern many times
        y = rng.random(n) < 0.45
        train = rng.random(n) < 0.8
        fit = train & (rng.random(n) < 0.85)
        words = C.packed(x)
        for k in (5, 40, 4000):
            n_cand, tops = O.screen_fold(words, train, y, 3, {"train": train, "fit": fit}, k)
            for name, mask in (("train", train), ("fit", fit)):
                want_n, want = plain.screen(x, y, train, 3, mask, k)
                self.assertEqual(n_cand, want_n)
                self.assertEqual([int(r) for r in tops[name]], want, f"{name} at {k}")

    def test_the_fixture_screens_are_the_plain_screens(self):
        with tempfile.TemporaryDirectory() as t:
            d, x, _ = C.write_export(t)
            meta = json.loads((d / "dataset.json").read_text())
            y = np.array(meta["resistant"])
            for s in meta["splits"]:
                outer = np.array(s["plan"]["outer"])
                for sc in s["screens"]:
                    train = outer != sc["fold"]
                    want_n, want = plain.screen(x, y, train, 3, train, 4000)
                    self.assertEqual((sc["n_candidates"], sc["screened_rows"]), (want_n, want))


class Plans(Temp):
    def test_the_exports_plan_is_read_with_its_holdout_inside_the_training_genomes(self):
        ex = O.Export(self.export())
        plan, folds = O.export_plan(ex, "lineage", 42)
        self.assertEqual(len(folds), 5)
        O.check_folds(folds, ex.y, ex.n)
        for fo in folds:
            self.assertFalse(np.any(fo["hold"] & fo["test"]))
            self.assertTrue(np.array_equal(fo["fit"] | fo["hold"], fo["train"]))
            self.assertEqual(fo["seed"], 42 + 101 * fo["fold"])

    def test_a_holdout_that_holds_an_outer_test_genome_is_refused(self):
        def change(meta):
            p = meta["splits"][1]["plan"]
            test = [i for i, f in enumerate(p["outer"]) if f == 0]
            p["holdouts"][0]["gbt"] = sorted(set(p["holdouts"][0]["gbt"]) | {test[0]})
            p["sha256"] = O.canonical_sha256(p)
        d = self.export(change=change)
        with self.assertRaisesRegex(O.Refuse, "holdout holds a genome outside"):
            O.run("xgboost", d, "lineage", 42, fitter=lambda *a: (None, {}))

    def test_a_plan_whose_hash_is_not_its_content_is_refused(self):
        def change(meta):
            meta["splits"][1]["plan"]["outer"][0] = (meta["splits"][1]["plan"]["outer"][0] + 1) % 5
        d = self.export(change=change)
        with self.assertRaisesRegex(O.Refuse, "sha256"):
            O.run("xgboost", d, "lineage", 42)

    def test_a_plan_of_another_seed_is_refused(self):
        d = self.export(seed=43)
        with self.assertRaisesRegex(O.Refuse, "seed 43, not 42"):
            O.run("xgboost", d, "lineage", 42)

    def test_a_fold_file_is_not_read_on_a_drawn_split(self):
        d = self.export()
        f = C.write_fold_file(self.top)
        with self.assertRaisesRegex(O.Refuse, "not read on the lineage split"):
            O.run("xgboost", d, "lineage", 42, fold_file=f)

    def test_the_study_split_needs_a_study_whole_fold_file(self):
        d = self.export()
        f = C.write_fold_file(self.top)
        doc = json.loads(f.read_text())
        doc["scheme"] = "first-public-past"
        f.write_text(json.dumps(doc))
        with self.assertRaisesRegex(O.Refuse, "not study-whole"):
            O.run("xgboost", d, "study-whole", 42, fold_file=f)

    def test_study_holdouts_and_inner_folds_lie_inside_each_folds_training_genomes(self):
        ex = O.Export(self.export())
        plan, folds, sha = O.study_plan(ex, C.write_fold_file(self.top), 42)
        self.assertEqual(len(folds), 3)
        self.assertEqual(plan["floor"], 3)
        O.check_folds(folds, ex.y, ex.n)
        for fo in folds:
            self.assertFalse(np.any(fo["hold"] & fo["test"]))
            self.assertFalse(np.any(fo["hold"] & ~fo["train"]))
            self.assertTrue(np.array_equal(fo["inner"] >= 0, fo["train"]))
            self.assertFalse(np.any(fo["train"][:2]) or np.any(fo["test"][:2]), "an unused genome reached a fold")
            self.assertEqual(bool(fo["train"][3]), bool(fo["train"][4]))

    def test_the_study_draws_are_the_same_on_a_second_read(self):
        ex = O.Export(self.export())
        f = C.write_fold_file(self.top)
        _, a, _ = O.study_plan(ex, f, 42)
        _, b, _ = O.study_plan(ex, f, 42)
        for fa, fb in zip(a, b):
            self.assertTrue(np.array_equal(fa["hold"], fb["hold"]) and np.array_equal(fa["inner"], fb["inner"]))

    def test_a_fold_file_that_leaves_a_genome_out_is_refused(self):
        d = self.export()
        f = C.write_fold_file(self.top)
        doc = json.loads(f.read_text())
        doc["train"] = doc["train"][1:]
        f.write_text(json.dumps(doc))
        with self.assertRaisesRegex(O.Refuse, "leaves a cohort genome out"):
            O.run("xgboost", d, "study-whole", 42, fold_file=f)


def flipped_scores(case, arm, scheme="lineage", fold=2, **kw):
    """The scores of one outer fold's test genomes with the true labels and with those genomes' labels turned round."""
    a = case.export("true")
    b = case.export("flip", flip=(scheme, fold))
    ra, sa = O.run(arm, a, scheme, 42, **kw)
    rb, sb = O.run(arm, b, scheme, 42, **kw)
    meta = json.loads((a / "dataset.json").read_text())
    outer = np.array(next(s["plan"]["outer"] for s in meta["splits"] if s["plan"]["scheme"] == scheme))
    test = outer == fold
    return sa[test], sb[test], ra, rb


class Arms(Temp):
    def unitigs(self):
        return C.write_unitigs(self.top)

    def test_no_outer_test_label_moves_an_xgboost_score(self):
        a, b, ra, _ = flipped_scores(self, "xgboost")
        self.assertTrue(np.array_equal(a, b))
        self.assertEqual(ra["constant_folds"], 0)

    def test_no_outer_test_label_moves_a_lightgbm_score(self):
        a, b, ra, _ = flipped_scores(self, "lightgbm")
        self.assertTrue(np.array_equal(a, b))
        self.assertEqual(ra["constant_folds"], 0)

    def test_no_outer_test_label_moves_an_elastic_net_score(self):
        # On the random split: the SYNTHETIC labels follow rows that follow lineages, so on the lineage split the
        # pick rightly takes the null fit on most folds, and a pick read on the test genomes would go unseen.
        a, b, ra, _ = flipped_scores(self, "elasticnet", scheme="random", unitigs=self.unitigs())
        self.assertTrue(np.array_equal(a, b))
        self.assertEqual(ra["null_fit_folds"], 0)

    def test_two_runs_write_the_same_bytes(self):
        d = self.export()
        u = self.unitigs()
        for arm, kw in (("xgboost", {}), ("lightgbm", {}), ("elasticnet", {"unitigs": u})):
            self.assertEqual(O.dump(O.run(arm, d, "lineage", 42, **kw)[0]), O.dump(O.run(arm, d, "lineage", 42, **kw)[0]), arm)

    def test_the_study_split_runs_on_its_own_folds(self):
        d = self.export()
        r, _ = O.run("xgboost", d, "study-whole", 42, fold_file=C.write_fold_file(self.top))
        self.assertEqual(len(r["folds"]), 3)
        self.assertEqual(r["tested"] + 2 + sum(1 for _ in json.loads(C.write_fold_file(self.top).read_text())["train"]), 150)
        self.assertIsNotNone(r["fold_file_sha256"])

    def test_each_folds_calls_are_its_test_genomes(self):
        d = self.export()
        r, s = O.run("lightgbm", d, "lineage", 42)
        for row in r["folds"]:
            c = row["calls"]
            self.assertEqual(c["n"], row["test"])
            self.assertEqual(c["tp"] + c["fn"], row["test_resistant"])
            self.assertEqual(c["tp"] + c["fp"] + c["tn"] + c["fn"], c["n"])
        self.assertEqual(int(np.sum(~np.isnan(s))), r["tested"])

    def test_a_fold_with_no_screened_row_scores_its_training_share(self):
        def change(meta):
            meta["splits"][1]["screens"][1]["screened_rows"] = []
        d = self.export(change=change)
        r, s = O.run("xgboost", d, "lineage", 42)
        row = r["folds"][1]
        self.assertEqual((row["constant"], row["screened"], r["constant_folds"]), (1, 0, 1))
        meta = json.loads((d / "dataset.json").read_text())
        outer = np.array(meta["splits"][1]["plan"]["outer"])
        share = row["train_resistant"] / row["train"]
        self.assertTrue(np.all(s[outer == 1] == share))

    def test_the_report_holds_only_its_listed_keys_and_no_genome(self):
        d = self.export()
        for arm, kw in (("xgboost", {}), ("elasticnet", {"unitigs": self.unitigs()})):
            r, _ = O.run(arm, d, "lineage", 42, **kw)
            self.assertLessEqual(set(r), O.REPORT_KEYS)
            for row in r["folds"]:
                self.assertLessEqual(set(row), O.FOLD_KEYS)
                self.assertTrue(all(isinstance(v, (int, float)) for k, v in row.items() if k != "calls"))
            text = O.dump(r)
            self.assertNotIn("s0001", text)
        bad = dict(r, features=["row 1"])
        with self.assertRaisesRegex(O.Refuse, "unlisted keys"):
            O.check_report(bad)

    def test_another_library_version_is_refused(self):
        d = self.export()
        with mock.patch.dict(O.VERSIONS, {"xgboost": "3.0.4"}):
            with self.assertRaisesRegex(O.Refuse, "not 3.0.4"):
                O.run("xgboost", d, "lineage", 42)


class Plain(Temp):
    """The arms against the plain paths of plain.py, on SYNTHETIC tables."""

    def test_the_xgboost_arm_is_the_plain_booster_with_xgboosts_rules(self):
        d = self.export()
        with mock.patch.dict(O.GRM_GBT, {"colsample": 1.0, "max_depth": 2}):
            got = {}

            def keep(X, y, fo, w):
                s, info = O.fit_xgboost(X, y, fo, w)
                fit, hold, test = O.sides(fo)
                want, run, kept = plain.boost(X, y, w, fit, hold, test, 400, 0.1, 2, 1.0, 1.0, 25)
                got[fo["fold"]] = (s, want, info, run, kept)
                return s, info
            O.run("xgboost", d, "lineage", 42, fitter=keep)
        self.assertEqual(len(got), 5)
        for f, (s, want, info, run, kept) in got.items():
            self.assertEqual((info["rounds_run"], info["rounds_kept"]), (run, kept), f"fold {f}")
            self.assertLess(float(np.max(np.abs(s - want))), 1e-4, f"fold {f}")

    def test_the_lightgbm_arm_is_the_plain_booster_with_lightgbms_rules(self):
        d = self.export()
        patch = {"feature_fraction": 1.0, "max_depth": 2, "num_leaves": 4}
        with mock.patch.dict(O.LIGHTGBM, patch):
            got = {}

            def keep(X, y, fo, w):
                s, info = O.fit_lightgbm(X, y, fo, w)
                fit, hold, test = O.sides(fo)
                want, run, kept = plain.boost(X, y, w, fit, hold, test, 400, 0.1, 2, 1.0, 1.0, 25, min_count=20,
                                              split_eps=0.0, eps=1e-15, hessian_count=True)
                got[fo["fold"]] = (s, want, info, run, kept)
                return s, info
            O.run("lightgbm", d, "lineage", 42, fitter=keep)
        self.assertEqual(len(got), 5)
        for f, (s, want, info, run, kept) in got.items():
            self.assertEqual((info["rounds_run"], info["rounds_kept"]), (run, kept), f"fold {f}")
            self.assertLess(float(np.max(np.abs(s - want))), 1e-4, f"fold {f}")

    def saga_against_fista(self, rel):
        x, y, _ = C.table(120, 60, seed=21)
        X = x.T.astype(np.float64)
        fit = np.arange(120)
        w = O.class_weights(y, np.ones(120, dtype=bool))
        scale = O.scaled(X, fit)
        Xs = X * scale
        for mix in (0.5, 1.0):
            top = O.lambda_max(X, y, w, scale, mix)
            lams = list(np.geomspace(top, top * 0.05, 4))
            path = O.enet_path(X, y, w, scale, mix, lams, 42)
            for lam, (coef, b, ok) in zip(lams, path):
                self.assertTrue(ok)
                pc, pb = plain.enet_fista(Xs, y, w, lam, mix)
                mine = plain.enet_objective(Xs, y, w, coef, b, lam, mix)
                theirs = plain.enet_objective(Xs, y, w, pc, pb, lam, mix)
                self.assertLess(mine - theirs, rel * abs(theirs), f"mix {mix} penalty {lam}")
                self.assertGreater(mine - theirs, -1e-6)
            self.assertEqual(int(np.sum(path[0][0] != 0)), 0, "lambda_max leaves every weight at 0")
            self.assertGreater(int(np.sum(path[-1][0] != 0)), 0)

    def test_saga_solves_the_stated_objective(self):
        # Run to a tight tolerance, saga and the plain solver reach the same minimum: the objective, the class
        # weights, the column scale and C = 1 / (lambda n) are as the runner states them.
        with mock.patch.multiple(O, ENET_TOL=1e-9, ENET_MAX_ITER=200000):
            self.saga_against_fista(1e-6)

    def test_saga_at_the_arms_tolerance_is_near_the_minimum(self):
        # At the arm's own stopping rule (tol 1e-4) saga stops a little short of the minimum: within 1e-5 of it here.
        self.saga_against_fista(1e-5)

    def test_the_elastic_net_projects_the_unitig_build_on_the_cohort(self):
        d = self.export()
        u = C.write_unitigs(self.top)
        r, _ = O.run("elasticnet", d, "lineage", 42, unitigs=u)
        x, _, _ = C.table(150, 160)
        meta = json.loads((u / "matrix.json").read_text())
        pats = np.fromfile(u / "patterns.bin", dtype="<u8").reshape(meta["n_patterns"], -1)
        bits = np.unpackbits(pats.view(np.uint8), axis=1, bitorder="little")[:, :150].astype(bool)
        keep = (bits.sum(axis=1) >= 3) & (bits.sum(axis=1) <= 147)
        self.assertEqual(r["input_rows"], len({tuple(b) for b in bits[keep]}))
        for row in r["folds"]:
            self.assertEqual(row["grid_points"], 40)
            self.assertLessEqual(row["screened"], row["candidates"])
            self.assertIn(row["picked_mix"], O.MIXES)


if __name__ == "__main__":
    unittest.main()
