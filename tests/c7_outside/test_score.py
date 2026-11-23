"""scripts/c7-outside/score.py and expect.py: the checks, the pairing, the count and its conditions
(docs/outside-learners.md). Every result, report and dataset name is SYNTHETIC; the default's and the twin's
results are C9's writer's at no inner fold (tests/c9_pick/common.py), as C3's tests make them.
"""

import contextlib
import csv
import io
import json
import tempfile
import unittest
from pathlib import Path

from tests.c3_inputs import common as c3c
from tests.c7_outside import common as C
from tests.c9_pick import common as c9

score = C.load("score")
expect = C.load("expect")
O = C.O


def arm_report(base, arm, scheme, seed):
    """A SYNTHETIC arm report on the folds of the default's result `base`."""
    shape = [(f["train"], f["train_resistant"], f["test"], f["test_resistant"]) for f in base["folds"]]
    return C.outside_report(dataset=base["dataset"], scheme=scheme, seed=seed, arm=arm, shape=shape)


def full_run(datasets=(C.DATASET, C.OTHER), study=(C.DATASET,), seeds=(42, 43, 44), twin=(45, 46, 47), arm_seeds=(42,),
             change=None):
    """SYNTHETIC expect rows and results of one species' whole run. change(row, doc) may alter a result."""
    rows = expect.rows(C.SPECIES, list(datasets), list(study), list(seeds), list(twin), list(arm_seeds))
    docs, releases = {}, {}
    for r in rows:
        dataset, _, arm, split, scheme, seed, role, name = r
        base = c3c.pick(dataset=dataset, scheme=scheme, seed=seed)
        d = base if arm in ("k31", "twin") else arm_report(base, arm, scheme, seed)
        if change:
            change(r, d)
        docs[name] = d
        if arm == "k31" and split == "lineage" and seed == 42:
            releases[dataset] = c9.release(base, dataset)
    return rows, docs, releases


def run(rows, docs, releases, drop=()):
    """score.py's first form on a SYNTHETIC run: (exit code, stdout, checks, draws, datasets)."""
    with tempfile.TemporaryDirectory() as tmp:
        _, top, rel, out = c3c.write_run(tmp, [], {k: v for k, v in docs.items() if k not in drop}, {}, releases)
        exp = Path(tmp) / "c7-expect.tsv"
        with open(exp, "w") as f:
            f.write("\t".join(expect.COLUMNS) + "\n")
            for r in rows:
                f.write("\t".join(str(x) for x in r) + "\n")
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            rc = score.main([str(exp), str(top), str(rel), str(out)])
        return rc, buf.getvalue(), table(out / "checks.tsv"), table(out / "draws.tsv"), table(out / "datasets.tsv")


def table(path):
    if not path.exists():
        return []
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def failed(checks):
    return sorted({c["check"] for c in checks if c["passed"] != "1"})


def at(arm, scheme="lineage", seed=42, dataset=C.DATASET):
    return lambda r, d: r[0] == dataset and r[2] == arm and r[4] == scheme and r[5] == seed


class Expect(unittest.TestCase):
    def test_the_default_runs_every_seed_and_the_arms_their_own(self):
        rows = expect.rows(C.SPECIES, [C.DATASET, C.OTHER], [C.DATASET], [42, 43, 44], [45, 46, 47], [42])
        n = {}
        for r in rows:
            n[(r[2], r[3])] = n.get((r[2], r[3]), 0) + 1
        self.assertEqual(n[("k31", "lineage")], 6)
        self.assertEqual(n[("twin", "lineage")], 6)
        for arm in ("xgboost", "lightgbm", "elasticnet"):
            self.assertEqual(n[(arm, "lineage")], 2)
            self.assertEqual(n[(arm, "study")], 1)
        self.assertEqual(n[("k31", "study")], 1)
        self.assertNotIn(("twin", "study"), n)
        self.assertIn(f"results/outside/elasticnet/{C.DATASET}.study-whole-s42.outside.json", [r[7] for r in rows])
        self.assertIn(f"c3/results/picks/twin/{C.OTHER}.lineage-s47.pick.json", [r[7] for r in rows])

    def test_the_elastic_net_left_out_has_no_row(self):
        rows = expect.rows(C.SPECIES, [C.DATASET], [C.DATASET], [42, 43, 44], [45, 46, 47], [42, 43, 44],
                           ("xgboost", "lightgbm"))
        self.assertEqual({r[2] for r in rows}, {"k31", "twin", "xgboost", "lightgbm"})
        self.assertEqual(sum(r[2] == "xgboost" for r in rows), 4)
        with self.assertRaisesRegex(ValueError, "two boosters"):
            expect.rows(C.SPECIES, [C.DATASET], [], [42], [45], [42], ("xgboost", "elasticnet"))

    def test_bad_lists_are_refused(self):
        for args, msg in (
            (([C.DATASET], [], [43], [45], [43]), "release's seed"),
            (([C.DATASET], [], [42], [42], [42]), "both"),
            (([C.DATASET], [], [42], [45], [42, 43]), "not a seed of the default"),
            (([C.DATASET], [C.OTHER], [42], [45], [42]), "not a listed dataset"),
            ((["other_species__drug_a"], [], [42], [45], [42]), "not a dataset of"),
        ):
            with self.assertRaisesRegex(ValueError, msg):
                expect.rows(C.SPECIES, *args)


class Checks(unittest.TestCase):
    def test_a_whole_run_passes_and_scores_each_arm_against_its_grm_learner(self):
        rows, docs, releases = full_run()
        rc, out, checks, draws, datasets = run(rows, docs, releases)
        self.assertEqual(rc, 0, failed(checks))
        self.assertEqual({r["arm"] for r in draws}, {"k31", "twin", "xgboost", "lightgbm", "elasticnet"})
        learners = {}
        for r in draws:
            learners.setdefault(r["arm"], set()).add(r["learner"])
        self.assertEqual(learners, {"k31": {"gbt", "lasso"}, "twin": {"gbt", "lasso"}, "xgboost": {"gbt"},
                                    "lightgbm": {"gbt"}, "elasticnet": {"lasso"}})
        self.assertIn("checks:", out)
        self.assertNotIn("balanced", out)
        self.assertEqual(datasets[0]["species"], "Synthia coli")

    def test_an_unconverged_refit_is_carried_to_the_draw_and_the_dataset(self):
        def change(r, d):
            if at("elasticnet", seed=43)(r, d):
                d["unconverged_refits"] = 1
        rows, docs, releases = full_run(arm_seeds=(42, 43, 44), change=change)
        rc, out, checks, draws, datasets = run(rows, docs, releases)
        self.assertEqual(rc, 0, failed(checks))
        hit = [(r["arm"], r["split"], r["seed"]) for r in draws if r["unconverged"] == "1"]
        self.assertEqual(hit, [("elasticnet", "lineage", "43")])
        hit = [(r["arm"], r["split"], r["unconverged_draws"]) for r in datasets
               if r["dataset"] == C.DATASET and r["unconverged_draws"] not in ("0", "")]
        self.assertEqual(hit, [("elasticnet", "lineage", "1")])

    def test_each_arm_pairs_with_its_own_grm_learner(self):
        self.assertEqual(O.ARMS, {"xgboost": "gbt", "lightgbm": "gbt", "elasticnet": "lasso"})
        self.assertEqual(score.FAMILY, [("xgboost", "gbt"), ("lightgbm", "gbt"), ("elasticnet", "lasso"),
                                        ("twin", "gbt"), ("twin", "lasso")])

    def test_an_arm_reported_against_another_learner_fails(self):
        def change(r, d):
            if at("elasticnet")(r, d):
                d["pairs_with"] = "gbt"

        rc, _, checks, _, _ = run(*full_run(change=change))
        self.assertEqual(rc, 1)
        self.assertEqual(failed(checks), ["pairs with its GRM learner"])

    def test_an_arm_on_other_folds_fails(self):
        def change(r, d):
            if at("lightgbm")(r, d):
                d["folds"][0]["test"] += 1
                d["folds"][0]["calls"]["n"] += 1
                d["folds"][0]["calls"]["tn"] += 1

        rc, _, checks, _, _ = run(*full_run(change=change))
        self.assertEqual(rc, 1)
        self.assertEqual(failed(checks), ["folds are the default's"])

    def test_a_study_report_read_as_a_lineage_one_fails(self):
        # The study default's folds have the same shape here, so only the split names it.
        def change(r, d):
            if at("xgboost")(r, d):
                d["scheme"] = "study-whole"

        rc, _, checks, _, _ = run(*full_run(change=change))
        self.assertEqual(rc, 1)
        self.assertEqual(failed(checks), ["dataset, split and seed"])

    def test_a_report_with_a_field_outside_its_list_fails(self):
        def change(r, d):
            if at("xgboost")(r, d):
                d["folds"][1]["features"] = ["row 17"]

        rc, _, checks, _, _ = run(*full_run(change=change))
        self.assertEqual(rc, 1)
        self.assertIn("listed keys only", failed(checks))

    def test_another_library_version_fails(self):
        def change(r, d):
            if at("lightgbm")(r, d):
                d["version"] = "4.5.0"

        rc, _, checks, _, _ = run(*full_run(change=change))
        self.assertEqual(rc, 1)
        self.assertEqual(failed(checks), ["arm and library"])

    def test_calls_that_are_not_the_test_genomes_fail(self):
        def change(r, d):
            if at("elasticnet")(r, d):
                d["folds"][2]["calls"]["tp"] += 1
                d["folds"][2]["calls"]["fn"] -= 1
                d["folds"][2]["calls"]["fp"] += 1
                d["folds"][2]["calls"]["n"] += 1

        rc, _, checks, _, _ = run(*full_run(change=change))
        self.assertEqual(rc, 1)
        self.assertEqual(failed(checks), ["calls are the test genomes'"])

    def test_a_missing_draw_fails_and_is_never_read_as_a_zero(self):
        rows, docs, releases = full_run()
        gone = f"results/outside/xgboost/{C.DATASET}.lineage-s42.outside.json"
        rc, _, checks, draws, _ = run(rows, docs, releases, drop=(gone,))
        self.assertEqual(rc, 1)
        self.assertEqual(failed(checks), ["exists"])
        self.assertEqual(draws, [])
        # Scored directly, the missing draw leaves its dataset unscored, with no value at all.
        names = {d: ("Synthia coli", "drug a") for d in (C.DATASET, C.OTHER)}
        rows_ = [dict(zip(expect.COLUMNS, r)) for r in rows]
        kept = {k: v for k, v in docs.items() if k != gone}
        got = score.score_rows(rows_, kept, names)
        miss = [r for r in got if r["arm"] == "xgboost" and r["dataset"] == C.DATASET and r["split"] == "lineage"]
        self.assertEqual([(r["scored"], r["why"]) for r in miss], [(0, "draw missing")])
        self.assertTrue(all(m not in miss[0] for m in score.MEASURES))
        ds = [r for r in score.dataset_rows(got) if r["arm"] == "xgboost" and r["dataset"] == C.DATASET
              and r["split"] == "lineage"]
        self.assertEqual(ds[0]["scored"], 0)
        self.assertTrue(all(ds[0][m] is None for m in score.MEASURES))

    def test_an_arm_meets_the_default_at_its_own_seeds(self):
        rows, docs, releases = full_run(arm_seeds=(42,))
        rc, _, checks, draws, datasets = run(rows, docs, releases)
        self.assertEqual(rc, 0, failed(checks))
        n = {(r["arm"], r["split"], r["learner"]): int(r["draws"]) for r in datasets if r["dataset"] == C.DATASET}
        self.assertEqual(n[("k31", "lineage", "gbt")], 3)
        self.assertEqual(n[(score.ARM_BASE, "lineage", "gbt")], 1)
        self.assertEqual(n[(score.ARM_BASE, "study", "lasso")], 1)
        self.assertEqual(n[("xgboost", "lineage", "gbt")], 1)
        self.assertEqual(n[("twin", "lineage", "lasso")], 3)

    def test_an_arm_at_a_seed_the_default_did_not_run_fails(self):
        rows, docs, releases = full_run()
        extra = [C.DATASET, C.SPECIES, "lightgbm", "lineage", "lineage", 49, "draw",
                 f"results/outside/lightgbm/{C.DATASET}.lineage-s49.outside.json"]
        docs[extra[7]] = arm_report(c3c.pick(dataset=C.DATASET, seed=49), "lightgbm", "lineage", 49)
        rc, _, checks, _, _ = run(rows + [extra], docs, releases)
        self.assertEqual(rc, 1)
        self.assertEqual(failed(checks), ["folds are the default's"])

    def test_the_scores_hold_no_genome_or_feature_field(self):
        rows, docs, releases = full_run()
        _, _, _, draws, datasets = run(rows, docs, releases)
        for t in (draws, datasets):
            for name in t[0]:
                self.assertNotIn(name, {"genome", "genomes", "id", "ids", "feature", "features", "rows", "kmer"})


def count_rows(shift, n=12, species=4, study=True, base_shift=0.0, constant=(), unconverged=()):
    """SYNTHETIC datasets.tsv rows: n datasets in `species` species; the default's values, the default at the arms'
    seeds (plus base_shift), and each pair's arm the default's plus shift(arm, learner, measure)."""
    out = []
    base = {"balanced_accuracy": 0.8, "recall": 0.7, "specificity": 0.9, "vme": 0.3, "me": 0.1, "log_loss": 0.4}
    for k in range(n):
        sp = f"synthia_{k % species}"
        d = f"{sp}__drug_{k}"
        for split in ("lineage", "study") if study else ("lineage",):
            cells = [("k31", l) for l in score.PAIRED] + [(score.ARM_BASE, l) for l in score.PAIRED] + list(score.FAMILY)
            for arm, learner in cells:
                if split == "study" and arm == "twin":
                    continue
                row = {"dataset": d, "species_slug": sp, "species": sp, "antibiotic": f"drug {k}", "split": split,
                       "arm": arm, "learner": learner, "draws": 3, "scored_draws": 3, "scored": 1,
                       "constant_draws": int(arm in constant and k == 0),
                       "unconverged_draws": int(arm in unconverged and k in (0, 1))}
                for m, v in base.items():
                    add = 0.0 if arm == "k31" else base_shift if arm == score.ARM_BASE else shift(arm, learner, m)
                    row[m] = v + add
                out.append(row)
    return out


def count(rows, withdrawn=()):
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "datasets.tsv"
        score.write_tsv(path, rows, score.DATASET_COLUMNS)
        out = Path(tmp) / "count"
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            rc = score.main(["--count", str(out)] + [a for w in withdrawn for a in ("--withdrawn", w)] + [str(path)])
        return rc, {p.name: table(p) for p in sorted(out.glob("*.tsv"))}


def by(rows):
    return {(r["split"], r["pair"], r["learner"], r["measure"]): r for r in rows}


class Count(unittest.TestCase):
    def test_an_arm_is_read_against_the_default_at_its_seeds_and_the_twin_against_every_draw(self):
        # The default at the arms' seeds is 0.1 above the default; each arm equals it, the twin equals the default.
        rc, t = count(count_rows(lambda a, l, m: 0.0 if a == "twin" else 0.1, base_shift=0.1))
        self.assertEqual(rc, 0)
        at_ = by(t["count.tsv"])
        for arm, learner in score.FAMILY:
            r = at_[("lineage", f"{arm}-k31", learner, "balanced_accuracy")]
            self.assertEqual(float(r["median_difference"]), 0.0, arm)
            self.assertEqual(int(r["entered"]), 12)

    def test_holm_runs_over_the_five_lineage_pairs_of_the_primary_measure(self):
        rc, t = count(count_rows(lambda a, l, m: 0.05 if a == "lightgbm" and m == "balanced_accuracy" else 0.0))
        self.assertEqual(rc, 0)
        family = [r for r in t["count.tsv"] if r["split"] == "lineage" and r["measure"] == "balanced_accuracy"]
        self.assertEqual(sorted((r["pair"], r["learner"]) for r in family),
                         [("elasticnet-k31", "lasso"), ("lightgbm-k31", "gbt"), ("twin-k31", "gbt"), ("twin-k31", "lasso"),
                          ("xgboost-k31", "gbt")])
        p = 0.5 ** 12
        lg = [r for r in family if r["pair"] == "lightgbm-k31"][0]
        self.assertEqual(float(lg["sign_p"]), p)
        self.assertEqual(float(lg["holm_p"]), 5 * p)
        self.assertTrue(all(r["holm_p"] == "" for r in t["count.tsv"] if r not in family))
        self.assertFalse(any(r["pair"] == "twin-k31" for r in t["count.tsv"] if r["split"] == "study"))

    def test_the_count_refuses_part_of_the_family(self):
        rows = [r for r in count_rows(lambda a, l, m: 0.0) if not (r["arm"] == "twin" and r["learner"] == "lasso")]
        with self.assertRaises(SystemExit) as e:
            count(rows)
        self.assertIn("twin (lasso)", str(e.exception))

    def test_the_count_refuses_a_dataset_in_m_and_counts_the_first_run_datasets_it_holds(self):
        with open(C.KIT / "lists.tsv", encoding="utf-8", newline="") as f:
            lists = list(csv.DictReader(f, delimiter="\t"))
        in_m = next(r for r in lists if r["first_run"] == "0")
        first = next(r for r in lists if r["first_run"] == "1")
        rows = count_rows(lambda a, l, m: 0.0, n=4)
        for r in rows:
            if r["dataset"] == "synthia_1__drug_1":
                r.update(dataset=in_m["dataset"], species_slug=in_m["species_slug"])
        with self.assertRaises(SystemExit) as e:
            count(rows)
        self.assertIn(f"1 datasets are in M, outside C7's first run, and never counted: {in_m['dataset']}", str(e.exception))
        for r in rows:
            if r["dataset"] == in_m["dataset"]:
                r.update(dataset=first["dataset"], species_slug=first["species_slug"])
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "datasets.tsv"
            score.write_tsv(path, rows, score.DATASET_COLUMNS)
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                rc = score.main(["--count", str(Path(tmp) / "count"), str(path)])
        self.assertEqual(rc, 0)
        n_first = sum(r["first_run"] == "1" for r in lists)
        self.assertEqual(buf.getvalue(), f"count: {len(rows)} dataset rows from 1 files, 1 of {n_first} first-run datasets\n")

    def test_conditions_are_read_for_the_three_arms_and_never_for_the_twin(self):
        rc, t = count(count_rows(lambda a, l, m: 0.0))
        self.assertEqual(rc, 0)
        cond = t["conditions.tsv"]
        self.assertEqual(sorted({r["arm"] for r in cond}), ["elasticnet", "lightgbm", "xgboost"])
        self.assertEqual(len(cond), 21)
        self.assertEqual(sum(r["decides"] == "0" for r in cond), 3)

    def test_an_arm_that_gains_everywhere_meets_one_to_five(self):
        def shift(a, l, m):
            return {"balanced_accuracy": 0.05, "recall": 0.03}.get(m, 0.0) if a == "elasticnet" else 0.0

        rc, t = count(count_rows(shift))
        self.assertEqual(rc, 0)
        met = {(r["arm"], int(r["condition"])): r["met"] for r in t["conditions.tsv"] if r["decides"] == "1"}
        self.assertEqual([met[("elasticnet", c)] for c in range(1, 7)], ["1", "1", "1", "1", "1", ""])
        self.assertEqual([met[("xgboost", c)] for c in range(1, 6)], ["0", "0", "1", "1", "0"])

    def test_condition_three_fails_when_the_twin_meets_one_in_either_pair(self):
        for learner in ("gbt", "lasso"):
            def shift(a, l, m, learner=learner):
                if m != "balanced_accuracy":
                    return 0.0
                return 0.05 if a == "elasticnet" or (a == "twin" and l == learner) else 0.0

            rc, t = count(count_rows(shift))
            self.assertEqual(rc, 0)
            met = {(r["arm"], int(r["condition"])): r["met"] for r in t["conditions.tsv"] if r["decides"] == "1"}
            self.assertEqual(met[("elasticnet", 1)], "1")
            self.assertEqual(met[("elasticnet", 3)], "0", learner)

    def test_a_recall_interval_below_the_margin_fails_condition_four(self):
        def shift(a, l, m):
            return {"balanced_accuracy": 0.05, "recall": -0.05}.get(m, 0.0) if a == "xgboost" else 0.0

        rc, t = count(count_rows(shift))
        met = {(r["arm"], int(r["condition"])): r["met"] for r in t["conditions.tsv"] if r["decides"] == "1"}
        self.assertEqual((met[("xgboost", 1)], met[("xgboost", 4)]), ("1", "0"))

    def test_a_constant_fold_dataset_enters_the_count_and_is_left_out_beside_it(self):
        rc, t = count(count_rows(lambda a, l, m: 0.0, constant=("lightgbm",)))
        self.assertEqual(rc, 0)
        a = by(t["count.tsv"])[("lineage", "lightgbm-k31", "gbt", "balanced_accuracy")]
        b = by(t["count-no-constant.tsv"])[("lineage", "lightgbm-k31", "gbt", "balanced_accuracy")]
        self.assertEqual((int(a["entered"]), int(a["constant"])), (12, 1))
        self.assertEqual((int(b["entered"]), int(b["not_scored"])), (11, 1))

    def test_each_arm_is_read_beside_against_its_own_learners_twin_pair(self):
        # The twin meets 1 on gbt only: condition 3 fails for every arm, and the reading beside, which decides
        # nothing, fails for the two boosters and holds for the elastic net, whose learner is the lasso.
        def shift(a, l, m):
            return 0.05 if m == "balanced_accuracy" and a == "twin" and l == "gbt" else 0.0

        rc, t = count(count_rows(shift))
        self.assertEqual(rc, 0)
        met = {(r["arm"], int(r["condition"]), r["decides"]): r["met"] for r in t["conditions.tsv"]}
        self.assertEqual([met[(a, 3, "1")] for a in ("xgboost", "lightgbm", "elasticnet")], ["0", "0", "0"])
        self.assertEqual([met[(a, 3, "0")] for a in ("xgboost", "lightgbm", "elasticnet")], ["0", "0", "1"])

    def test_a_withdrawn_arm_keeps_its_pair_at_p_one_and_is_named(self):
        rows = [r for r in count_rows(lambda a, l, m: 0.05 if a == "lightgbm" and m == "balanced_accuracy" else 0.0)
                if r["arm"] != "elasticnet"]
        rc, t = count(rows, withdrawn=["elasticnet"])
        self.assertEqual(rc, 0)
        family = {r["pair"] + " " + r["learner"]: r for r in t["count.tsv"]
                  if r["split"] == "lineage" and r["measure"] == "balanced_accuracy"}
        self.assertEqual(len(family), 5)
        e = family["elasticnet-k31 lasso"]
        self.assertEqual((e["withdrawn"], e["entered"], float(e["sign_p"]), float(e["holm_p"])), ("1", "0", 1.0, 1.0))
        self.assertEqual(float(family["lightgbm-k31 gbt"]["holm_p"]), 5 * 0.5 ** 12)
        self.assertEqual({r["withdrawn"] for r in t["count.tsv"] if not r["pair"].startswith("elasticnet")}, {"0"})

    def test_a_withdrawn_arm_with_results_is_refused(self):
        with self.assertRaises(SystemExit) as e:
            count(count_rows(lambda a, l, m: 0.0), withdrawn=["elasticnet"])
        self.assertIn("a withdrawn arm has results: elasticnet", str(e.exception))

    def test_the_twin_or_the_default_cannot_be_withdrawn(self):
        rows = count_rows(lambda a, l, m: 0.0)
        for arm in ("twin", "k31"):
            with self.assertRaises(SystemExit) as e:
                count([r for r in rows if r["arm"] != arm], withdrawn=[arm])
            self.assertIn("not the twin or the default", str(e.exception), arm)

    def test_an_unconverged_refit_dataset_enters_the_count_and_is_left_out_beside_it(self):
        rc, t = count(count_rows(lambda a, l, m: 0.0, unconverged=("elasticnet",)))
        self.assertEqual(rc, 0)
        key = ("lineage", "elasticnet-k31", "lasso", "balanced_accuracy")
        a, b = by(t["count.tsv"])[key], by(t["count-no-unconverged.tsv"])[key]
        self.assertEqual((int(a["entered"]), int(a["unconverged"])), (12, 2))
        self.assertEqual((int(b["entered"]), int(b["not_scored"])), (10, 2))
        self.assertEqual(by(t["count-no-unconverged.tsv"])[("lineage", "xgboost-k31", "gbt", "balanced_accuracy")]["entered"],
                         "12")

    def test_a_dataset_listed_twice_is_refused(self):
        rows = count_rows(lambda a, l, m: 0.0)
        with self.assertRaises(SystemExit):
            count(rows + rows[:1])


if __name__ == "__main__":
    unittest.main()
