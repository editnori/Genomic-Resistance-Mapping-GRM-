"""SYNTHETIC tests of scripts/c5-settings/score.py: made-up results, choices files, reports and dataset rows."""

import contextlib
import copy
import csv
import io
import json
import tempfile
import unittest
from pathlib import Path

from tests.c5_settings import common
from tests.c9_pick import common as c9

S = common.load("score")
E = common.load("expect")
LEARNERS = common.LEARNERS


def choices(base, k, dataset, scheme, seed, change=None):
    """One SYNTHETIC choices file for fold k of the default's result `base`, as drive.py writes it."""
    f = base["folds"][k]
    rows = []
    for li, l in enumerate(LEARNERS):
        d = {key: f["outer"][li][key] for key in S.CALL_KEYS}
        moved = dict(d)
        if l != "scm" and d["fn"] > 0:
            moved.update(tp=d["tp"] + 1, fn=d["fn"] - 1)
        rows.append({"learner": l, "status": "picked", "trials": 31, "failed_trials": 0, "picked_trial": 0,
                     "picked_setting": "release", "plain_best_trial": 4, "tied_at_best": 1,
                     "release_inner_mean": 0.7, "best_inner_mean": 0.72, "trials_above_gamma": 0,
                     "trials_better": 3, "trials_same": 0, "trials_worse": 27, "cutoff": 0.45 if l != "scm" else 0.5,
                     "cutoff_inner_balanced_accuracy": 0.71, "half_inner_balanced_accuracy": 0.7, "search_note": None,
                     "outer": {"default": d, "search": dict(d), "cutoff": moved},
                     "seconds_inner": 1.0, "seconds_outer": 0.1})
    doc = {"format": "grm-c5-choices-1", "evaluation_revision": 14, "dataset": dataset, "scheme": scheme,
           "seed": seed, "fold": k, "tuning": "lineage" if scheme == "lineage" else "copies", "inner_seed": 1,
           "outer_seed": seed + 101 * k, "sampler_seed": 20261030, "gamma": 0.04, "grid": S.GRID,
           "random_trials": 30, "counts": {"train": f["train"], "train_resistant": f["train_resistant"],
                                           "test": f["test"], "test_resistant": f["test_resistant"]},
           "inner_folds": [], "min_support_inner": 3, "min_support_outer": 3, "learners": rows, "timing": []}
    if change:
        change(doc)
    return doc


def write_run(top, seeds=(42, 43, 44), twin=(45, 46, 47), change=None, drop=()):
    """expect.tsv, the default's and the twin's results, the choices files and the release reports of a
    SYNTHETIC run of two datasets, the first with a study plan. change(name, doc) may alter a choices file;
    drop names choices files left unwritten."""
    top = Path(top)
    rows = E.rows(common.SPECIES, [common.DATASET, common.OTHER], [common.DATASET], list(seeds), list(twin))
    docs, files, releases = {}, {}, {}
    for d, sp, arm, split, scheme, seed, role, name in rows:
        if arm in ("default", "twin"):
            doc = c9.doc(dataset=d, scheme=scheme, seed=seed, inner=0)
            docs[name] = doc
            if arm == "default":
                for k in range(len(doc["folds"])):
                    n = E.choices_name(d, scheme, seed, k)
                    if n in drop:
                        continue
                    files[n] = choices(doc, k, d, scheme, seed, (lambda doc, n=n: change(n, doc)) if change else None)
                if split == "lineage" and seed == 42:
                    releases[d] = c9.release(doc, d)
    with open(top / "expect.tsv", "w") as f:
        f.write("\t".join(E.COLUMNS) + "\n")
        for r in rows:
            f.write("\t".join(str(x) for x in r) + "\n")
    for name, doc in list(docs.items()) + list(files.items()):
        path = top / "picks" / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(doc))
    (top / "release" / "datasets").mkdir(parents=True, exist_ok=True)
    index = []
    for d, rep in releases.items():
        (top / "release" / "datasets" / f"{d}.json").write_text(json.dumps(rep))
        rnd = [(l["metrics"]["balanced_accuracy"], i, l["learner"]) for s in rep["splits"] if s["scheme"] == "random"
               for i, l in enumerate(s["learners"])]
        index.append({"id": d, "best_learner": max(rnd)[2]})
    (top / "release" / "index.json").write_text(json.dumps({"datasets": index}))
    return [str(top / "expect.tsv"), str(top / "picks"), str(top / "release"), str(top / "score")]


def score(args):
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        rc = S.main(args)
    return rc, buf.getvalue()


def tsv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


class Score(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.top = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def test_a_whole_run_passes_and_scores_every_arm(self):
        rc, text = score(write_run(self.top))
        self.assertEqual(rc, 0, text)
        self.assertIn("missing draws 0", text)
        ds = tsv(self.top / "score" / "datasets.tsv")
        arms = {(r["dataset"], r["split"], r["arm"]) for r in ds}
        self.assertIn((common.DATASET, "study", "search"), arms)
        self.assertIn((common.OTHER, "lineage", "twin"), arms)
        self.assertNotIn((common.OTHER, "study", "search"), arms)
        self.assertTrue(all(r["scored"] == "1" for r in ds))
        # The cutoff arm moved one call per fold and learner from fn to tp, scm aside.
        cut = {r["learner"]: r for r in ds if r["arm"] == "cutoff" and r["split"] == "lineage" and r["dataset"] == common.OTHER}
        base = {r["learner"]: r for r in ds if r["arm"] == "default" and r["split"] == "lineage" and r["dataset"] == common.OTHER}
        self.assertEqual(int(cut["scm"]["tp"]), int(base["scm"]["tp"]))
        folds = len(c9.doc()["folds"])
        self.assertEqual(int(cut["rf"]["tp"]), int(base["rf"]["tp"]) + 3 * folds)

    def test_a_missing_draw_is_not_scored_and_never_read_as_zero(self):
        gone = E.choices_name(common.OTHER, "lineage", 43, 2)
        rc, text = score(write_run(self.top, drop={gone}))
        self.assertEqual(rc, 0, text)
        self.assertIn("missing draws 2", text)
        draws = [r for r in tsv(self.top / "score" / "draws.tsv")
                 if r["dataset"] == common.OTHER and r["seed"] == "43" and r["arm"] in ("search", "cutoff")]
        self.assertEqual(len(draws), 10)
        for r in draws:
            self.assertEqual((r["scored"], r["why"], r["balanced_accuracy"], r["n"]), ("0", "missing draw", "", ""))
        for r in tsv(self.top / "score" / "datasets.tsv"):
            if r["dataset"] == common.OTHER and r["split"] == "lineage" and r["arm"] in ("search", "cutoff"):
                self.assertEqual((r["scored"], r["balanced_accuracy"]), ("0", ""), r["arm"])

    def test_a_choices_file_whose_default_is_not_the_defaults_fit_fails(self):
        def change(name, doc):
            if name == E.choices_name(common.DATASET, "lineage", 44, 1):
                d = doc["learners"][2]["outer"]["default"]
                d["tp"], d["fn"] = d["tp"] - 1, d["fn"] + 1
        rc, _ = score(write_run(self.top, change=change))
        self.assertEqual(rc, 1)
        bad = [r for r in tsv(self.top / "score" / "checks.tsv") if r["passed"] != "1"]
        self.assertEqual([r["check"] for r in bad], ["fold 1 rf: the default's calls are the default's own fit"])

    def test_a_pick_against_the_rule_fails(self):
        def change(name, doc):
            if name == E.choices_name(common.OTHER, "lineage", 42, 0):
                doc["learners"][3]["picked_trial"] = 4
        rc, _ = score(write_run(self.top, change=change))
        self.assertEqual(rc, 1)
        bad = [r["check"] for r in tsv(self.top / "score" / "checks.tsv") if r["passed"] != "1"]
        self.assertEqual(bad, ["fold 0 gbt: the pick follows the default-anchored rule"])

    def test_scms_cutoff_arm_must_be_the_defaults_calls(self):
        def change(name, doc):
            if name == E.choices_name(common.DATASET, "study-whole", 42, 2):
                c = doc["learners"][0]["outer"]["cutoff"]
                c["tp"], c["fn"] = c["tp"] + 1, c["fn"] - 1
        rc, _ = score(write_run(self.top, change=change))
        self.assertEqual(rc, 1)
        bad = [r["check"] for r in tsv(self.top / "score" / "checks.tsv") if r["passed"] != "1"]
        self.assertEqual(bad, ["fold 2 scm: the cutoff arm's calls are the default's"])

    def test_an_arm_at_a_twin_seed_fails(self):
        args = write_run(self.top)
        with open(args[0], "a") as f:
            f.write("\t".join([common.DATASET, common.SPECIES, "search", "lineage", "lineage", "45", "draw",
                               E.choices_name(common.DATASET, "lineage", 45)]) + "\n")
        rc, _ = score(args)
        self.assertEqual(rc, 1)
        # The row has no choices files either; the test holds only when the twin-seed check is the one that fails.
        failed = {(r["file"], r["check"]) for r in tsv(self.top / "score" / "checks.tsv") if r["passed"] != "1"}
        self.assertIn((E.choices_name(common.DATASET, "lineage", 45),
                       "an arm runs at a draw's seed or the study's, never the twin's"), failed)

    def test_a_default_that_differs_from_the_release_fails(self):
        args = write_run(self.top)
        rel = Path(args[2]) / "datasets" / f"{common.OTHER}.json"
        doc = json.loads(rel.read_text())
        lineage = next(x for x in doc["splits"] if x["scheme"] == "lineage")
        lineage["learners"][1]["metrics"]["balanced_accuracy"] += 0.01
        rel.write_text(json.dumps(doc))
        rc, _ = score(args)
        self.assertEqual(rc, 1)


def dataset_rows(arms=("default", "search", "cutoff", "twin"), species=8, better=("search",)):
    """SYNTHETIC datasets.tsv rows: 3 datasets in each of `species` species; the arms in `better` gain 0.05
    on every dataset, the others match the default."""
    rows = []
    for s in range(species):
        sp = f"synthia_{s}"
        for k in range(3):
            d = f"{sp}__drug_{k}"
            for split in ("lineage", "study"):
                for arm in arms:
                    if split == "study" and arm == "twin":
                        continue
                    for li, l in enumerate(LEARNERS):
                        v = 0.7 + 0.01 * k + 0.002 * li + (0.05 if arm in better else 0.0)
                        row = {"dataset": d, "species_slug": sp, "species": sp, "antibiotic": f"drug {k}",
                               "split": split, "arm": arm, "learner": l, "draws": 3, "scored_draws": 3,
                               "constant_draws": 0, "scored": 1}
                        for m in S.C3.MEASURES:
                            row[m] = v if m in ("balanced_accuracy", "recall", "specificity") else 0.3 - (v - 0.7)
                        rows.append(row)
    return rows


class Count(unittest.TestCase):
    def test_holm_runs_over_the_fifteen_lineage_rows_of_the_primary_measure(self):
        rows = S.C3.count(dataset_rows())
        held = [r for r in rows if r["holm_p"] is not None]
        self.assertEqual(len(held), 15)
        self.assertEqual({r["pair"] for r in held}, {"search-default", "cutoff-default", "twin-default"})
        self.assertTrue(all(r["split"] == "lineage" and r["measure"] == "balanced_accuracy" for r in held))
        self.assertEqual({r["pair"] for r in rows if r["split"] == "study"}, {"search-default", "cutoff-default"})
        search = [r for r in held if r["pair"] == "search-default"]
        self.assertTrue(all(r["better"] == 24 and r["worse"] == 0 for r in search))

    def test_conditions_read_the_two_arms_with_the_twin_for_condition_3(self):
        cond = S.C3.conditions(S.C3.count(dataset_rows()))
        self.assertEqual({r["arm"] for r in cond}, {"search", "cutoff"})
        first = next(r for r in cond if r["arm"] == "search" and r["reading"] == "ruling 15" and r["condition"] == "1 and 2")
        self.assertIn("after Holm over 15;", first["reads"])
        self.assertEqual((first["met"], first["learners"]), (1, 5))
        three = next(r for r in cond if r["arm"] == "search" and r["reading"] == "ruling 15" and r["condition"] == "3")
        self.assertEqual((three["met"], three["learners"]), (1, 0))
        cut = next(r for r in cond if r["arm"] == "cutoff" and r["reading"] == "ruling 15" and r["condition"] == "1 and 2")
        self.assertEqual(cut["met"], 0)

    def test_the_count_refuses_a_family_without_the_twin(self):
        with self.assertRaises(SystemExit):
            S.C3.count(dataset_rows(arms=("default", "search", "cutoff")))

    def test_c3s_own_scorer_keeps_its_family_of_thirty(self):
        c3 = common.load("score").C3
        self.assertEqual(c3.FAMILY, ["search", "cutoff", "twin"])
        import importlib.util
        spec = importlib.util.spec_from_file_location("c3_alone", common.ROOT / "scripts" / "c3-inputs" / "score.py")
        alone = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(alone)
        self.assertEqual((len(alone.FAMILY), alone.DEFAULT), (6, "k31"))


if __name__ == "__main__":
    unittest.main()
