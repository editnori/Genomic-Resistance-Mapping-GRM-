"""scripts/c2-cap/score.py: the checks before any score, the scores per learner, and the count (docs/cap-arms.md).

Every genome, plan, prediction row and pick file is SYNTHETIC. The dataset
names are the committed first-run list's (scripts/c2-cap/lists.tsv), which
the scorer reads; no call or score of theirs is here.
"""

import csv
import json
import shutil
import tempfile
import unittest
from pathlib import Path

from tests.c2_cap import common

score = common.load("score")

D = "enterococcus_faecium__vancomycin"
D2 = "enterococcus_faecium__ampicillin"
M_DATASET = "enterococcus_faecium__linezolid"
DRAWS, TWIN = 2, 2
COLUMNS = ["dataset", "species_slug", "arm", "split", "slot", "role", "fit"]


def table(path):
    if not path.is_file():
        return None
    with open(path, encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def lineage_rows(d, twin_fit="cap-2000", twin_slots=None):
    rows = [[d, "enterococcus_faecium", "cap-2000", "lineage", str(s), "draw", "cap-2000"] for s in range(DRAWS)]
    for s in (twin_slots if twin_slots is not None else range(DRAWS, DRAWS + TWIN)):
        rows.append([d, "enterococcus_faecium", "twin", "lineage", str(s), "twin", twin_fit])
    rows += [[d, "enterococcus_faecium", "cap-all", "lineage", str(s), "draw", "cap-all"] for s in range(DRAWS)]
    rows += [[d, "enterococcus_faecium", "cap-4000", "lineage", str(s), "draw", "cap-all"] for s in range(DRAWS)]
    return rows


def pick(d, folds, learners=None):
    """A SYNTHETIC grm learn --exp-pick file at seed 42 with no inner fold."""
    out = []
    for i, (train, test, test_r) in enumerate(folds):
        outer = [{"n": test, "tp": test_r - 1, "fn": 1, "tn": test - test_r - j, "fp": j, "loss_sum": 0.3 * test}
                 for j in range(5)]
        out.append({"fold": i, "train": train, "train_resistant": train // 3, "test": test, "test_resistant": test_r,
                    "outer": outer})
    return {"format": "grm-c9-pick-1", "evaluation_revision": 14, "dataset": d, "scheme": "study-whole", "seed": 42,
            "inner_folds": 0, "learners": [{"learner": l} for l in (learners or common.LEARNERS)], "folds": out,
            "tested": sum(f[1] for f in folds), "tested_resistant": sum(f[2] for f in folds), "unavailable_reason": None}


BASE_FOLDS = [(40, 10, 4), (42, 8, 3), (41, 9, 4)]
WIDE_FOLDS = [(52, 10, 4), (55, 8, 3), (50, 9, 4)]


class Run:
    """A SYNTHETIC C2 run of one or two datasets: the default's and cap-all's lineage fits and study picks."""

    def __init__(self, tmp, datasets=(D,), study=True, keyed=None, refused=()):
        self.t = Path(tmp)
        self.sp = common.Species()
        self.expect = []
        for d in datasets:
            common.lineage_fits(self.sp, self.t / "fits", d, DRAWS + TWIN, DRAWS, keyed=keyed, refused=refused)
            # With a refused default plan grm arm writes no lineage row for the default, nor for an arm that
            # reaches a refused slot among its draws.
            common.predictions(self.t / "fits", "cap-2000", d, [] if refused else range(DRAWS + TWIN))
            common.predictions(self.t / "fits", "cap-all", d, [] if any(s < DRAWS for s in refused) else range(DRAWS))
            self.expect += lineage_rows(d)
            if study:
                common.write_json(self.t / "picks" / "cap-2000" / f"{d}.study-whole-s42.pick.json", pick(d, BASE_FOLDS))
                common.write_json(self.t / "picks" / "cap-all" / f"{d}.study-whole-s42.pick.json", pick(d, WIDE_FOLDS))
                self.expect += [[d, "enterococcus_faecium", "cap-2000", "study", "-", "study", "cap-2000"],
                                [d, "enterococcus_faecium", "cap-all", "study", "-", "study", "cap-all"],
                                [d, "enterococcus_faecium", "cap-4000", "study", "-", "study", "cap-all"]]

    def run(self):
        exp = self.t / "expect.tsv"
        exp.write_text("\t".join(COLUMNS) + "\n" + "".join("\t".join(r) + "\n" for r in self.expect), encoding="utf-8")
        out = self.t / "score"
        rc, text = common.quiet(score.main, [str(exp), str(self.t / "fits"), str(self.t / "picks"), str(out)])
        return rc, text, table(out / "checks.tsv"), table(out / "draws.tsv"), table(out / "datasets.tsv")


def failed(checks):
    return [c["check"] for c in checks if c["passed"] != "1"]


class Checks(unittest.TestCase):
    def test_a_correct_run_scores_every_row_and_learner(self):
        with tempfile.TemporaryDirectory() as tmp:
            rc, text, checks, draws, datasets = Run(tmp).run()
        self.assertEqual(rc, 0, text)
        self.assertEqual(failed(checks), [])
        lineage = [r for r in draws if r["split"] == "lineage"]
        self.assertEqual(len(lineage), (DRAWS + TWIN + 2 * DRAWS) * 5)
        self.assertEqual({r["fit"] for r in lineage if r["arm"] == "cap-4000"}, {"cap-all"})
        self.assertEqual({(r["arm"], r["slot"]) for r in lineage if r["role"] == "twin"}, {("twin", "2"), ("twin", "3")})
        # Each learner's lineage calls are the 60 default genomes, once each.
        self.assertTrue(all(r["n"] == "60" for r in lineage))
        study = [r for r in draws if r["split"] == "study"]
        self.assertEqual(len(study), 3 * 5)
        self.assertEqual({r["n"] for r in study}, {"27"})
        self.assertEqual(sorted({(r["arm"], r["draws"]) for r in datasets if r["split"] == "lineage"}),
                         [("cap-2000", "2"), ("cap-4000", "2"), ("cap-all", "2"), ("twin", "2")])
        self.assertIn(f"checks: {len(checks)} of {len(checks)} pass", text)

    def test_a_replaced_draw_is_scored_at_its_slot_not_its_keyed_repeat(self):
        # Slot 1's refused draw replaced by keyed repeat 4, slot 3's by 5: plans hold 4 and 5, prediction rows the
        # slots 1 and 3, as grm arm writes them (src/ablate/arm.rs).
        with tempfile.TemporaryDirectory() as tmp:
            rc, text, checks, draws, _ = Run(tmp, keyed={1: DRAWS + TWIN, 3: DRAWS + TWIN + 1}).run()
        self.assertEqual(rc, 0, failed(checks))
        lineage = [r for r in draws if r["split"] == "lineage"]
        self.assertEqual(len(lineage), (DRAWS + TWIN + 2 * DRAWS) * 5)
        self.assertEqual(sorted({r["slot"] for r in lineage}), ["0", "1", "2", "3"])
        self.assertTrue(all(r["n"] == "60" for r in lineage))

    def test_an_unavailable_lineage_split_is_scored_zero_on_every_arm_and_the_run_goes_on(self):
        # A draw slot (1) and a twin slot (3) refused past grm arm's replacements: in the second, cap-all's own
        # draws ran and wrote rows, yet no pair with the default exists, so none is read.
        for refused in ({1}, {3}):
            with self.subTest(refused=refused), tempfile.TemporaryDirectory() as tmp:
                rc, text, checks, draws, datasets = Run(tmp, refused=refused).run()
                self.assertEqual(rc, 0, failed(checks))
                lineage = [r for r in draws if r["split"] == "lineage"]
                self.assertEqual(len(lineage), (DRAWS + TWIN + 2 * DRAWS) * 5)
                self.assertEqual({(r["scored"], r["why"]) for r in lineage}, {("0", "lineage split unavailable")})
                self.assertEqual({r["scored"] for r in draws if r["split"] == "study"}, {"1"})
                self.assertEqual({r["scored"] for r in datasets if r["split"] == "lineage"}, {"0"})
                self.assertEqual({r["scored"] for r in datasets if r["split"] == "study"}, {"1"})

    def test_an_arm_scored_on_other_genomes_than_the_defaults_outer_test_folds_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            r = Run(tmp, study=False)
            p = r.t / "fits" / "cap-all" / D / "predictions.jsonl"
            lines = p.read_text(encoding="utf-8").splitlines()
            # One scored genome of slot 0 moves to another outer fold, for every learner.
            first = json.loads(lines[0])
            moved = []
            for line in lines:
                row = json.loads(line)
                if row["sample_id"] == first["sample_id"] and row["fold_plan_sha256"] == "plan0":
                    row["fold"] = (row["fold"] + 1) % common.K
                moved.append(json.dumps(row))
            p.write_text("\n".join(moved) + "\n", encoding="utf-8")
            plan = r.t / "fits" / "cap-all" / D / "plans" / "lineage-0.cohort.json"
            doc = json.loads(plan.read_text(encoding="utf-8"))
            ids = json.loads((r.t / "fits" / "cap-all" / D / "training.jsonl").read_text(encoding="utf-8")
                             .splitlines()[0])["sample_ids"]
            i = ids.index(first["sample_id"])
            doc["outer"][i] = (doc["outer"][i] + 1) % common.K
            plan.write_text(json.dumps(doc), encoding="utf-8")
            rc, text, checks, draws, _ = r.run()
        self.assertEqual(rc, 1)
        self.assertIn("the scored genomes of each outer fold are the default's", failed(checks))
        self.assertIsNone(draws)

    def test_a_twin_at_a_draw_slot_or_from_a_wider_fit_stops_the_scorer(self):
        for fit, slots in (("cap-2000", [0, 3]), ("cap-all", [2, 3])):
            with self.subTest(fit=fit, slots=slots), tempfile.TemporaryDirectory() as tmp:
                r = Run(tmp, study=False)
                r.expect = [e for e in r.expect if e[5] != "twin"]
                r.expect += [[D, "enterococcus_faecium", "twin", "lineage", str(s), "twin", fit] for s in slots]
                if fit == "cap-all":
                    for s in slots:
                        shutil.copy(r.t / "fits" / "cap-all" / D / "plans" / "lineage-0.cohort.json",
                                    r.t / "fits" / "cap-all" / D / "plans" / f"lineage-{s}.cohort.json")
                rc, text, checks, draws, _ = r.run()
                self.assertEqual(rc, 1)
                self.assertIn("the twin is the default's fit at a slot no draw has", failed(checks))
                self.assertIsNone(draws)

    def test_a_missing_draw_stops_the_scorer_and_is_never_a_zero(self):
        with tempfile.TemporaryDirectory() as tmp:
            r = Run(tmp, datasets=(D, D2), study=False)
            shutil.rmtree(r.t / "fits" / "cap-2000" / D2)
            shutil.rmtree(r.t / "fits" / "cap-all" / D2)
            rc, text, checks, draws, datasets = r.run()
        self.assertEqual(rc, 1)
        self.assertIn("the plan and the prediction rows exist", failed(checks))
        self.assertIsNone(draws)
        self.assertIsNone(datasets)

    def test_a_learner_or_slot_with_no_rows_stops_the_scorer(self):
        with tempfile.TemporaryDirectory() as tmp:
            r = Run(tmp, study=False)
            common.predictions(r.t / "fits", "cap-2000", D, range(DRAWS + TWIN), drop_slot=3)
            rc, _, checks, draws, _ = r.run()
        self.assertEqual(rc, 1)
        self.assertIn("scm has rows", failed(checks))
        self.assertIsNone(draws)

    def test_a_scored_genome_missing_from_a_fold_stops_the_scorer(self):
        with tempfile.TemporaryDirectory() as tmp:
            r = Run(tmp, study=False)
            p = r.t / "fits" / "cap-2000" / D / "predictions.jsonl"
            lines = p.read_text(encoding="utf-8").splitlines()
            p.write_text("\n".join(lines[1:]) + "\n", encoding="utf-8")
            rc, _, checks, _, _ = r.run()
        self.assertEqual(rc, 1)
        self.assertIn(f"scm: one row per scored genome in each of {common.K} folds", failed(checks))

    def test_a_dataset_of_m_or_outside_the_first_run_is_refused(self):
        for d in (M_DATASET, "synthia_coli__drug_a"):
            with self.subTest(d=d), tempfile.TemporaryDirectory() as tmp:
                rc, _, checks, draws, _ = Run(tmp, datasets=(d,), study=False).run()
                self.assertEqual(rc, 1)
                self.assertIn("the dataset is in the first run (CAP, not M)", failed(checks))
                self.assertIsNone(draws)

    def test_a_study_arm_on_other_test_folds_or_with_fewer_training_genomes_fails(self):
        for folds in ([(52, 10, 4), (55, 8, 2), (50, 9, 4)], [(52, 10, 4), (41, 8, 3), (50, 9, 4)]):
            with self.subTest(folds=folds), tempfile.TemporaryDirectory() as tmp:
                r = Run(tmp)
                common.write_json(r.t / "picks" / "cap-all" / f"{D}.study-whole-s42.pick.json", pick(D, folds))
                rc, _, checks, _, _ = r.run()
                self.assertEqual(rc, 1)
                self.assertIn("the test folds are the default's, with at least its training genomes", failed(checks))

    def test_a_pick_file_of_another_fit_fails(self):
        for change in ({"seed": 7}, {"inner_folds": 5}, {"scheme": "study"}, {"evaluation_revision": 13}):
            with self.subTest(change=change), tempfile.TemporaryDirectory() as tmp:
                r = Run(tmp)
                common.write_json(r.t / "picks" / "cap-all" / f"{D}.study-whole-s42.pick.json",
                                  dict(pick(D, WIDE_FOLDS), **change))
                rc, _, checks, _, _ = r.run()
                self.assertEqual(rc, 1)
                self.assertIn("the pick file is the dataset's study-whole fit at seed 42 with no inner fold",
                              failed(checks))

    def test_an_arm_that_adds_no_study_genome_is_not_scored(self):
        with tempfile.TemporaryDirectory() as tmp:
            r = Run(tmp)
            (r.t / "picks" / "cap-all" / f"{D}.study-whole-s42.pick.json").unlink()
            r.expect = [e[:5] + (["none-added"] if e[3] == "study" and e[2] != "cap-2000" else [e[5]]) + e[6:]
                        for e in r.expect]
            rc, text, checks, draws, _ = r.run()
        self.assertEqual(rc, 0, text)
        self.assertEqual({r["arm"] for r in draws if r["split"] == "study"}, {"cap-2000"})

    def test_a_none_added_or_same_cohort_lineage_row_is_not_scored(self):
        with tempfile.TemporaryDirectory() as tmp:
            r = Run(tmp, study=False)
            shutil.rmtree(r.t / "fits" / "cap-all")
            r.expect = [e[:5] + (["none-added" if e[2] == "cap-all" else "same-cohort"] if e[2] in ("cap-all", "cap-4000")
                                 else [e[5]]) + e[6:] for e in r.expect]
            rc, text, checks, draws, _ = r.run()
        self.assertEqual(rc, 0, text)
        self.assertEqual({x["arm"] for x in draws}, {"cap-2000", "twin"})

    def test_the_log_loss_is_clipped(self):
        self.assertAlmostEqual(score.loss(1, 0.0), -__import__("math").log(1e-4))
        self.assertAlmostEqual(score.loss(0, 1.0), -__import__("math").log(1e-4))
        self.assertAlmostEqual(score.loss(1, 0.5), -__import__("math").log(0.5))


class Count(unittest.TestCase):
    def dataset_rows(self, arms=("cap-2000", "cap-4000", "cap-all", "twin")):
        rows = []
        for i, d in enumerate([D, D2, "enterococcus_faecium__gentamicin", "staphylococcus_aureus__oxacillin"]):
            for split in ("lineage", "study"):
                for arm in arms:
                    if split == "study" and arm == "twin":
                        continue
                    for j, learner in enumerate(common.LEARNERS):
                        v = 0.7 + 0.01 * i + (0.05 if arm == "cap-all" else 0.0) + 0.001 * j
                        row = {"dataset": d, "species_slug": d.split("__")[0], "species": "x", "antibiotic": "y",
                               "split": split, "arm": arm, "learner": learner, "draws": 3, "scored_draws": 3,
                               "constant_draws": 0, "scored": 1, "n": 100, "tp": 1, "fp": 1, "tn": 1, "fn": 1}
                        for m in score.MEASURES:
                            row[m] = v
                        rows.append(row)
        return rows

    def run_count(self, rows, withdrawn=()):
        with tempfile.TemporaryDirectory() as tmp:
            t = Path(tmp)
            cols = list(rows[0])
            with open(t / "datasets.tsv", "w", encoding="utf-8", newline="") as f:
                w = csv.DictWriter(f, cols, delimiter="\t")
                w.writeheader()
                w.writerows(rows)
            flags = [x for a in withdrawn for x in ("--withdrawn", a)]
            rc, text = common.quiet(score.main, ["--count", str(t / "out")] + flags + [str(t / "datasets.tsv")])
            return rc, text, table(t / "out" / "count.tsv"), table(t / "out" / "conditions.tsv")

    def test_holm_runs_over_the_fifteen_lineage_rows_of_the_primary_measure(self):
        rc, text, count, conditions = self.run_count(self.dataset_rows())
        self.assertEqual(rc, 0, text)
        primary = [r for r in count if r["split"] == "lineage" and r["measure"] == score.C9.PRIMARY]
        self.assertEqual(len(primary), 15)
        self.assertEqual({r["pair"] for r in primary}, {"cap-4000-cap-2000", "cap-all-cap-2000", "twin-cap-2000"})
        self.assertTrue(all(r["holm_p"] not in ("", None) for r in primary))
        self.assertFalse([r for r in count if r["split"] == "study" and r["pair"] == "twin-cap-2000"])
        text_all = " ".join(" ".join(r.values()) for r in conditions)
        self.assertIn("Holm over 15", text_all)
        self.assertNotIn("Holm over 30", text_all)

    def test_every_difference_is_the_arms_value_less_the_defaults(self):
        _, _, count, _ = self.run_count(self.dataset_rows())
        wide = [r for r in count if r["pair"] == "cap-all-cap-2000" and r["measure"] == "balanced_accuracy"]
        self.assertTrue(all(abs(float(r["median_difference"]) - 0.05) < 1e-9 for r in wide))

    def test_a_withdrawn_arm_keeps_its_five_pairs_at_p_one_and_the_tables_name_it_withdrawn(self):
        rc, text, count, conditions = self.run_count(self.dataset_rows(arms=("cap-2000", "cap-all", "twin")),
                                                     withdrawn=["cap-4000"])
        self.assertEqual(rc, 0, text)
        primary = [r for r in count if r["split"] == "lineage" and r["measure"] == score.C9.PRIMARY]
        self.assertEqual(len(primary), 15)
        gone = [r for r in primary if r["pair"] == "cap-4000-cap-2000"]
        self.assertEqual([(r["withdrawn"], float(r["sign_p"]), r["datasets"]) for r in gone], [("1", 1.0, "0")] * 5)
        self.assertEqual({r["withdrawn"] for r in primary if r["pair"] != "cap-4000-cap-2000"}, {"0"})
        mid = [r for r in conditions if r["arm"] == "cap-4000"]
        self.assertTrue(mid)
        self.assertEqual({(r["withdrawn"], r["met"], r["learners"]) for r in mid}, {("1", "", "")})
        self.assertEqual({r["withdrawn"] for r in conditions if r["arm"] == "cap-all"}, {"0"})
        self.assertTrue(any(r["met"] != "" for r in conditions if r["arm"] == "cap-all"))

    def test_a_withdrawn_arm_with_rows_or_the_twin_withdrawn_is_refused(self):
        with self.assertRaises(SystemExit):
            self.run_count(self.dataset_rows(), withdrawn=["cap-4000"])
        with self.assertRaises(SystemExit):
            self.run_count(self.dataset_rows(arms=("cap-2000", "cap-4000", "cap-all")), withdrawn=["twin"])

    def test_a_none_added_dataset_is_out_of_the_arms_denominator(self):
        rows = [r for r in self.dataset_rows(arms=("cap-2000", "cap-all", "twin"))
                if not (r["arm"] == "cap-all" and r["dataset"] == D2)]
        rc, text, count, _ = self.run_count(rows, withdrawn=["cap-4000"])
        self.assertEqual(rc, 0, text)
        lin = {r["pair"]: r for r in count if r["split"] == "lineage" and r["measure"] == score.C9.PRIMARY
               and r["learner"] == common.LEARNERS[0]}
        self.assertEqual((lin["cap-all-cap-2000"]["datasets"], lin["twin-cap-2000"]["datasets"]), ("3", "4"))

    def test_the_count_refuses_part_of_the_family_or_a_row_twice(self):
        with self.assertRaises(SystemExit):
            self.run_count(self.dataset_rows(arms=("cap-2000", "cap-all", "twin")))
        rows = self.dataset_rows()
        with self.assertRaises(SystemExit):
            self.run_count(rows + rows[:1])


if __name__ == "__main__":
    unittest.main()
