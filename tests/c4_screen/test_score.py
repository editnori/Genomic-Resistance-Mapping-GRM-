"""scripts/c4-screen/score.py: the checks, the scores, the count over the 21 pairs and the conditions
(docs/screen-arms.md).

Every result, inputs file, report and dataset name is SYNTHETIC.
"""

import contextlib
import io
import re
import tempfile
import unittest
from pathlib import Path

from tests.c3_inputs import test_score as c3_test
from tests.c4_screen import common

score = common.load("score")
table, failed = c3_test.table, c3_test.failed
D = common.DATASET
ARM = f"cmh-4000/{D}.lineage-s42.pick.json"
ARM_INPUTS = f"cmh-4000/{D}.drawn-s42.inputs.json"
FLOOR = f"floor-3/{D}.lineage-s43.pick.json"
TWIN = f"twin/{D}.lineage-s45.pick.json"
STUDY_ARM = f"cmh-16000/{D}.study-whole-s42.pick.json"


def run(rows, docs, files, releases, extra=None):
    """score.py's first form on a SYNTHETIC run: (exit code, stdout, stderr, checks, draws, datasets)."""
    with tempfile.TemporaryDirectory() as tmp:
        exp, picks, rel, out = common.write_run(tmp, rows, docs, dict(files, **(extra or {})), releases)
        buf, err = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(err):
            rc = score.main([str(exp), str(picks), str(rel), str(out)])
        return (rc, buf.getvalue(), err.getvalue(), table(out / "checks.tsv"), table(out / "draws.tsv"),
                table(out / "datasets.tsv"))


def changed(files_or_docs, name, edit):
    edit(files_or_docs[name])
    return files_or_docs


class Checks(unittest.TestCase):
    def test_a_whole_run_passes_and_scores_each_learner_of_each_result(self):
        rows, docs, files, releases = common.full_run()
        rc, out, err, checks, draws, datasets = run(rows, docs, files, releases)
        self.assertEqual((rc, failed(checks), err), (0, [], ""))
        self.assertEqual(len(draws), len(rows) * 5)
        self.assertTrue(out.startswith("checks: "))
        # Per dataset, split, arm and learner; with the dataset's genome count.
        self.assertEqual(len(datasets), (2 * 4 + 2) * 5 + 4 * 5)
        self.assertEqual({r["genomes"] for r in datasets}, {"150"})
        self.assertEqual({r["scored"] for r in datasets}, {"1"})

    def test_each_arm_result_meets_the_inputs_checks(self):
        rows, docs, files, releases = common.full_run()
        checks = run(rows, docs, files, releases)[3]
        names = {c["check"] for c in checks if c["file"] == ARM}
        for want in ("inputs: every screen by the arm's rule", "inputs: rf, gbt and lasso screen once in each outer fold",
                     "scm and cart, which read no screen, make the default's calls in every fold",
                     "the result's scheme is its split's"):
            self.assertIn(want, names)
        floor = {c["check"] for c in checks if c["file"] == FLOOR}
        self.assertIn("below 800 genomes the floor arm's result is the default's, timing apart", floor)
        self.assertIn("inputs: the chi-square screens name no strata", floor)

    def test_an_arm_screened_by_the_defaults_rule_fails(self):
        rows, docs, files, releases = common.full_run()

        def default_screen(i):
            for s in i["screens"]:
                s["rule"] = "chi-square"
        rc, _, err, checks, draws, _ = run(rows, docs, changed(files, ARM_INPUTS, default_screen), releases)
        self.assertEqual((rc, failed(checks), draws), (1, ["inputs: every screen by the arm's rule"], []))
        self.assertIn(f"failed: {ARM}: inputs: every screen by the arm's rule", err)

    def test_an_arm_with_no_inputs_file_or_another_arms_fails(self):
        rows, docs, files, releases = common.full_run()
        del files[ARM_INPUTS]
        self.assertEqual(failed(run(rows, docs, files, releases)[3]), ["inputs: the arm's file exists"])
        rows, docs, files, releases = common.full_run()
        files[ARM_INPUTS] = common.inputs(docs[ARM], arm="cmh-16000")
        self.assertEqual(failed(run(rows, docs, files, releases)[3]), ["inputs: the arm's rule and budget"])
        rows, docs, files, releases = common.full_run()
        files[ARM_INPUTS] = common.inputs(docs[ARM], arm="floor-3")
        self.assertIn("inputs: format, revision, dataset, tag, seed and arm", failed(run(rows, docs, files, releases)[3]))

    def test_a_floor_arm_at_another_budget_fails(self):
        rows, docs, files, releases = common.full_run()
        name = f"floor-3/{D}.drawn-s43.inputs.json"
        files[name]["budget"] = 16000
        self.assertEqual(failed(run(rows, docs, files, releases)[3]), ["inputs: the arm's rule and budget"])

    def test_a_screen_that_reads_more_genomes_than_its_fold_trains_fails(self):
        rows, docs, files, releases = common.full_run()

        def wider(i):
            for s in i["screens"]:
                if s["split"] == "lineage" and s["learner"] == "gbt":
                    s["fitting"] = 101
        self.assertEqual(failed(run(rows, docs, changed(files, ARM_INPUTS, wider), releases)[3]),
                         ["inputs: each outer screen reads its fold's candidates, rf its training genomes, gbt and lasso no more"])

    def test_a_missing_or_doubled_screen_or_fit_time_fails(self):
        rows, docs, files, releases = common.full_run()
        files[ARM_INPUTS]["screens"].pop(0)
        self.assertEqual(failed(run(rows, docs, files, releases)[3]),
                         ["inputs: rf, gbt and lasso screen once in each outer fold"])
        rows, docs, files, releases = common.full_run()
        i = files[ARM_INPUTS]
        i["fit_times"].append(dict(i["fit_times"][0]))
        self.assertEqual(failed(run(rows, docs, files, releases)[3]),
                         ["inputs: a fit time for each learner of each outer fold"])

    def test_a_screen_over_the_budget_or_strata_below_the_informative_fail(self):
        rows, docs, files, releases = common.full_run()
        files[ARM_INPUTS]["screens"][0]["kept"] = 499
        self.assertEqual(failed(run(rows, docs, files, releases)[3]),
                         ["inputs: each screen keeps the budget, or every candidate if fewer"])
        rows, docs, files, releases = common.full_run()
        files[ARM_INPUTS]["screens"][0]["informative_strata"] = 13
        self.assertEqual(failed(run(rows, docs, files, releases)[3]),
                         ["inputs: each adjusted screen names its strata, no fewer than the informative ones"])

    def test_an_arm_at_another_floor_fails(self):
        rows, docs, files, releases = common.full_run()
        files[ARM_INPUTS]["min_support"] = 4
        self.assertEqual(failed(run(rows, docs, files, releases)[3]), ["inputs: the arm's floor, and the result's"])

    def test_a_screen_arm_whose_scm_or_cart_differs_from_the_default_fails(self):
        rows, docs, files, releases = common.full_run()
        fold = docs[ARM]["folds"][0]["outer"][1]
        fold["loss_sum"] += 0.5
        self.assertEqual(failed(run(rows, docs, files, releases)[3]),
                         ["scm and cart, which read no screen, make the default's calls in every fold"])

    def test_a_floor_arm_below_800_genomes_unlike_the_default_fails_and_at_800_may_differ(self):
        def lower(row, d, i):
            if row[2] == "floor-3":
                d["folds"][0]["outer"][2]["loss_sum"] += 0.5
        rows, docs, files, releases = common.full_run(change=lower)
        self.assertEqual(failed(run(rows, docs, files, releases)[3]),
                         ["below 800 genomes the floor arm's result is the default's, timing apart"])
        rows, docs, files, releases = common.full_run(samples=800, change=lower)
        rc, _, _, checks, _, datasets = run(rows, docs, files, releases)
        self.assertEqual((rc, failed(checks)), (0, []))
        self.assertEqual({r["genomes"] for r in datasets}, {"800"})

    def test_a_floor_arm_differing_only_in_its_timing_passes(self):
        def slower(row, d, i):
            if row[2] == "floor-3":
                d["seconds_inner"] = 9.0
        rows, docs, files, releases = common.full_run(change=slower)
        self.assertEqual(run(rows, docs, files, releases)[0], 0)

    def test_the_twin_passes_as_the_default_at_its_own_seeds_and_fails_with_an_inputs_file(self):
        rows, docs, files, releases = common.full_run()
        checks = run(rows, docs, files, releases)[3]
        twin = {c["check"] for c in checks if c["file"] == TWIN}
        self.assertIn("the default and the twin write no inputs file", twin)
        self.assertIn("a twin seed is the twin's and no draw's", twin)
        self.assertFalse(any(c.startswith("inputs:") or c.startswith("the folds are") for c in twin))
        extra = {f"twin/{D}.drawn-s45.inputs.json": common.inputs(docs[TWIN], arm="cmh-4000")}
        self.assertEqual(failed(run(rows, docs, files, releases, extra)[3]),
                         ["the default and the twin write no inputs file"])

    def test_a_study_result_read_as_a_lineage_result_fails(self):
        rows, docs, files, releases = common.full_run()
        for r in rows:
            if r[7] == STUDY_ARM:
                r[3] = "lineage"
        self.assertIn("the result's scheme is its split's", failed(run(rows, docs, files, releases)[3]))

    def test_an_arm_on_other_folds_fails(self):
        def shift(row, d, i):
            if row[7] == ARM:
                c3_test.shift_folds(d)
                i.update(common.inputs(d, arm="cmh-4000"))
        rows, docs, files, releases = common.full_run(change=shift)
        self.assertIn("the folds are the default's at the same dataset, split and seed, by class",
                      failed(run(rows, docs, files, releases)[3]))

    def test_a_missing_result_stops_the_scores(self):
        rows, docs, files, releases = common.full_run()
        del docs[FLOOR]
        rc, _, err, checks, draws, datasets = run(rows, docs, files, releases)
        self.assertEqual((rc, failed(checks), draws, datasets), (1, ["exists"], [], []))
        self.assertIn(f"failed: {FLOOR}: exists", err)


class Scores(unittest.TestCase):
    def draws(self):
        rows, docs, files, releases = common.full_run()
        expect = [dict(zip(common.c3.COLUMNS, map(str, r))) for r in rows]
        names = {d: ("Synthia coli", d.split("__")[1]) for d in (common.DATASET, common.OTHER)}
        return expect, score.C3.score_rows(expect, docs, names)

    def test_a_dataset_with_a_missing_draw_is_not_scored(self):
        expect, draws = self.draws()
        kept = [r for r in draws if not (r["arm"] == "floor-3" and r["seed"] == 43 and r["dataset"] == D)]
        rows = score.dataset_rows(kept, expect, {D: 150, common.OTHER: 150})
        by = {(r["dataset"], r["split"], r["arm"], r["learner"]): r for r in rows}
        for learner in common.LEARNERS:
            r = by[(D, "lineage", "floor-3", learner)]
            self.assertEqual((r["draws"], r["scored_draws"], r["scored"]), (3, 2, 0))
            self.assertTrue(all(r[m] is None for m in score.MEASURES))
            self.assertEqual(by[(D, "lineage", "cmh-4000", learner)]["scored"], 1)

    def test_every_draw_present_scores_each_dataset_as_c3_does(self):
        expect, draws = self.draws()
        mine = score.dataset_rows(draws, expect, {D: 150, common.OTHER: 150})
        theirs = score.C3.dataset_rows(draws)
        self.assertEqual([{k: v for k, v in r.items() if k != "genomes"} for r in mine], theirs)


def big(i):
    """SYNTHETIC genome counts: every fourth dataset has 900 genomes, the rest 150."""
    return 900 if i % 4 == 0 else 150


def count_rows(arms=("cmh-4000", "cmh-16000", "floor-3", "twin"), better=None, genomes=big, datasets=12):
    """SYNTHETIC datasets.tsv rows over three species: better(arm, learner, i) gives the arm's lift on dataset i;
    genomes(i) its genome count."""
    rows = []
    for i in range(datasets):
        sp = f"synthia_{'abc'[i % 3]}"
        d = f"{sp}__drug_{i}"
        for split in ("lineage", "study"):
            for arm in ("k31",) + tuple(arms):
                if arm == "twin" and split == "study":
                    continue
                for learner in common.LEARNERS:
                    lift = better(arm, learner, i) if better and arm != "k31" else 0.0
                    base = 0.7 + 0.001 * i
                    row = {"dataset": d, "species_slug": sp, "species": "Synthia", "antibiotic": f"drug {i}",
                           "genomes": genomes(i) if callable(genomes) else genomes, "split": split, "arm": arm, "learner": learner, "draws": 3,
                           "scored_draws": 3, "constant_draws": 0, "scored": 1}
                    for m in score.MEASURES:
                        row[m] = base + (lift if m in ("balanced_accuracy", "recall", "specificity") else -lift)
                    rows.append(row)
    return rows


class Count(unittest.TestCase):
    def test_the_family_is_twenty_one_pairs_and_holm_runs_over_them(self):
        rows = score.count(count_rows())
        fam = [r for r in rows if r["split"] == "lineage" and r["measure"] == score.PRIMARY and r["holm_p"] is not None]
        self.assertEqual(len(fam), 21)
        self.assertEqual(sorted({(r["pair"], r["set"]) for r in fam}),
                         sorted([("cmh-4000-k31", "all"), ("cmh-16000-k31", "all"), ("floor-3-k31", "800-or-more"),
                                 ("twin-k31", "all"), ("twin-k31", "800-or-more")]))
        self.assertEqual(sorted({(r["pair"], r["learner"]) for r in fam if r["pair"].startswith("cmh")}),
                         sorted((f"{a}-k31", l) for a in ("cmh-4000", "cmh-16000") for l in ("rf", "gbt", "lasso")))
        self.assertTrue(all(r["holm_p"] is None for r in rows if r not in fam))

    def test_holm_over_twenty_one_is_held_to_the_plain_step_down(self):
        lift = lambda arm, learner, i: 0.05 if arm == "floor-3" and i % 2 == 0 else 0.0
        rows = score.count(count_rows(better=lift, datasets=24))
        fam = [r for part in score.FAMILY_PARTS for r in rows if (score.arm_of(r["pair"]), r["set"]) == part
               and r["split"] == "lineage" and r["measure"] == score.PRIMARY]
        ps = [r["sign_p"] for r in fam]
        order = sorted(range(21), key=lambda k: ps[k])
        want, run_max = [None] * 21, 0.0
        for rank, k in enumerate(order):
            run_max = max(run_max, min(1.0, (21 - rank) * ps[k]))
            want[k] = run_max
        self.assertEqual([round(r["holm_p"], 12) for r in fam], [round(x, 12) for x in want])
        self.assertLess(min(r["sign_p"] for r in fam), 0.05)

    def test_the_floor_arm_is_counted_over_the_datasets_of_800_genomes_or_more(self):
        # SYNTHETIC: the floor arm lifts every dataset of 900 genomes and none of 150, as a floor that moves only
        # at 800 or more would; over the 6 of 24 it can move it is better in all 6, over all 24 its median is 0.
        lift = lambda arm, learner, i: 0.05 if arm == "floor-3" and big(i) >= 800 else 0.0
        rows = count_rows(better=lift, datasets=24)
        counted = score.count(rows)
        floor = [r for r in counted if r["pair"] == "floor-3-k31" and r["split"] == "lineage"
                 and r["measure"] == score.PRIMARY]
        self.assertEqual({(r["set"], r["datasets"], r["better"], r["median_difference"]) for r in floor},
                         {("800-or-more", 6, 6, 0.05)})
        beside = [r for r in score.floor_all_count(rows) if r["split"] == "lineage" and r["measure"] == score.PRIMARY]
        self.assertEqual({(r["set"], r["datasets"], r["better"], r["median_difference"], r["holm_p"]) for r in beside},
                         {("all", 24, 6, 0.0, None)})
        got = {(r["arm"], r["reading"], r["condition"]): r for r in score.conditions(counted)}
        self.assertEqual(got[("floor-3", "stricter, 800 genomes or more", "1")]["learners"], 5)

    def test_a_dataset_of_exactly_800_genomes_is_in_the_floor_arms_set_and_one_of_799_is_not(self):
        # SYNTHETIC: datasets alternate between 800 and 799 genomes. max(3, 800 / 200) is 4, so the floor arm
        # can move a dataset of 800; at 799 the default's floor is 3 already. The set holds the 800s alone.
        edge = lambda i: 800 if i % 2 == 0 else 799
        self.assertEqual([score.floor_moves(g) for g in (799, 800, "800")], [False, True, True])
        rows = count_rows(genomes=edge, datasets=10)
        self.assertEqual(sorted({r["dataset"] for r in score.moved_rows(rows)}),
                         sorted({r["dataset"] for r in rows if r["genomes"] == 800}))
        floor = [r for r in score.count(rows) if r["pair"] == "floor-3-k31" and r["split"] == "lineage"
                 and r["measure"] == score.PRIMARY]
        self.assertEqual({(r["set"], r["datasets"]) for r in floor}, {("800-or-more", 5)})

    def test_the_twin_is_read_on_each_arms_set(self):
        # SYNTHETIC: the twin moves on the datasets of 800 genomes or more alone. It voids condition 3 for the
        # floor arm, read over those datasets, and not for the screen arms, read over every dataset.
        lift = lambda arm, learner, i: 0.05 if arm == "twin" and big(i) >= 800 else 0.0
        got = {(r["arm"], r["reading"], r["condition"]): r
               for r in score.conditions(score.count(count_rows(better=lift, datasets=24)))}
        self.assertEqual(got[("floor-3", "ruling 15, 800 genomes or more", "3")]["met"], 0)
        self.assertEqual(got[("cmh-4000", "the learners it touches, every dataset", "3")]["met"], 1)
        self.assertEqual(got[("cmh-16000", "the learners it touches, every dataset", "3")]["met"], 1)

    def test_the_count_refuses_a_run_with_no_dataset_of_800_genomes(self):
        with self.assertRaises(SystemExit):
            score.count(count_rows(genomes=150))

    def test_a_screen_arms_scm_and_cart_are_not_counted(self):
        rows = score.count(count_rows())
        self.assertFalse([r for r in rows if r["pair"].startswith("cmh") and r["learner"] in ("scm", "cart")])
        self.assertEqual(len([r for r in rows if r["pair"] == "floor-3-k31" and r["split"] == "lineage"]),
                         5 * len(score.MEASURES))

    def test_the_count_refuses_part_of_the_family_unless_withdrawn(self):
        with self.assertRaises(SystemExit):
            score.count(count_rows(arms=("cmh-4000", "floor-3", "twin")))
        rows = score.count(count_rows(arms=("cmh-4000", "floor-3", "twin")), withdrawn=["cmh-16000"])
        gone = [r for r in rows if r["pair"] == "cmh-16000-k31" and r["split"] == "lineage"]
        self.assertEqual({(r["datasets"], r["sign_p"]) for r in gone}, {(0, 1.0)})
        with self.assertRaises(SystemExit):
            score.count(count_rows(arms=("cmh-4000", "cmh-16000", "floor-3")), withdrawn=["twin"])

    def test_the_main_form_writes_the_four_tables(self):
        rows = count_rows()
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "datasets.tsv"
            score.write_tsv(path, rows, score.DATASET_COLUMNS)
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                self.assertEqual(score.main(["--count", str(Path(tmp) / "out"), str(path)]), 0)
            for name in ("count", "count-no-constant", "conditions", "count-floor-all"):
                self.assertTrue((Path(tmp) / "out" / f"{name}.tsv").exists(), name)
        self.assertEqual(buf.getvalue().strip(), f"count: {len(rows)} dataset rows from 1 files")


class Conditions(unittest.TestCase):
    def read(self, better):
        return {(r["arm"], r["reading"], r["condition"]): r for r in score.conditions(score.count(count_rows(better=better)))}

    def test_a_screen_arm_is_read_on_the_three_learners_it_touches(self):
        lift = lambda arm, learner, i: 0.05 if arm == "cmh-4000" and learner in ("rf", "gbt", "lasso") else 0.0
        got = self.read(lift)
        reading = "the learners it touches, every dataset"
        self.assertEqual(sorted(k[2] for k in got if k[:2] == ("cmh-4000", reading)), ["1", "2", "3", "4", "5", "6"])
        self.assertEqual((got[("cmh-4000", reading, "1")]["met"], got[("cmh-4000", reading, "1")]["learners"]), (1, 3))
        self.assertEqual(got[("cmh-4000", reading, "5")]["met"], 1)
        self.assertEqual(got[("cmh-4000", reading, "3")]["met"], 1)
        self.assertEqual(got[("cmh-16000", reading, "1")]["met"], 0)
        self.assertIn("Holm over 21", got[("cmh-4000", reading, "2")]["reads"])

    def test_one_screened_learner_short_fails_condition_1(self):
        lift = lambda arm, learner, i: 0.05 if arm == "cmh-4000" and learner in ("rf", "gbt") else 0.0
        got = self.read(lift)
        self.assertEqual(got[("cmh-4000", "the learners it touches, every dataset", "1")]["learners"], 2)
        self.assertEqual(got[("cmh-4000", "the learners it touches, every dataset", "1")]["met"], 0)

    def test_the_floor_arm_takes_c3s_two_readings_at_holm_over_twenty_one(self):
        got = self.read(None)
        self.assertEqual(sorted({k[1] for k in got if k[0] == "floor-3"}),
                         ["ruling 15, 800 genomes or more", "stricter, 800 genomes or more"])
        # Every Holm figure of every row is 21, whatever family size C3's own text writes.
        self.assertEqual({m for r in got.values() for m in re.findall(r"Holm over (\d+)", r["reads"])}, {"21"})
        self.assertIn("Holm over 21", got[("floor-3", "ruling 15, 800 genomes or more", "1 and 2")]["reads"])
        self.assertNotIn("twin", {k[0] for k in got})

    def test_a_twin_that_moves_voids_condition_3_for_every_arm(self):
        lift = lambda arm, learner, i: 0.05 if arm in ("twin", "cmh-4000") and learner in ("scm", "cart", "rf") else 0.0
        got = self.read(lift)
        self.assertEqual(got[("cmh-4000", "the learners it touches, every dataset", "3")]["met"], 0)
        self.assertEqual(got[("floor-3", "ruling 15, 800 genomes or more", "3")]["met"], 0)


if __name__ == "__main__":
    unittest.main()
