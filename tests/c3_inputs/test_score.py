"""scripts/c3-inputs/score.py: the checks, the scores per learner, the count and its family (docs/input-arms.md).

Every result, inputs file, report and dataset name is SYNTHETIC.
"""

import contextlib
import csv
import io
import json
import tempfile
import unittest
from pathlib import Path

from tests.c3_inputs import common

score = common.load("score")


def run(rows, docs, files, releases, extra=None):
    """score.py's first form on a SYNTHETIC run: (exit code, stdout, checks, draws, datasets).

    extra maps a path under picks/ to a document written beside the run's.
    """
    with tempfile.TemporaryDirectory() as tmp:
        exp, picks, rel, out = common.write_run(tmp, rows, docs, dict(files, **(extra or {})), releases)
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            rc = score.main([str(exp), str(picks), str(rel), str(out)])
        return rc, buf.getvalue(), table(out / "checks.tsv"), table(out / "draws.tsv"), table(out / "datasets.tsv")


def table(path):
    if not path.exists():
        return []
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def failed(checks):
    return sorted({c["check"] for c in checks if c["passed"] != "1"})


FOLDS = "the folds are the default's at the same dataset, split and seed, by class"


def shift_folds(d):
    """Move one genome of each class from fold 0's test genomes to fold 1's, keeping the totals."""
    for f, dt in ((0, -1), (1, 1)):
        fold = d["folds"][f]
        fold["test"] += 2 * dt
        fold["test_resistant"] += dt
        fold["train"] -= 2 * dt
        fold["train_resistant"] -= dt
        for o in fold["outer"]:
            o["n"] += 2 * dt
            o["tp"] += dt
            o["tn"] += dt


class Checks(unittest.TestCase):
    def test_a_whole_run_passes_and_scores_each_learner_of_each_result(self):
        rows, docs, files, releases = common.full_run()
        rc, out, checks, draws, datasets = run(rows, docs, files, releases)
        self.assertEqual(rc, 0, failed(checks))
        self.assertEqual(len(draws), 5 * len(rows))
        self.assertEqual({r["arm"] for r in draws}, {"k31", "twin", "program", "k31-program"})
        # The program arms run on the A dataset only; the default and the twin on both.
        self.assertEqual({r["dataset"] for r in draws if r["arm"] == "program"}, {common.DATASET})
        self.assertEqual({r["dataset"] for r in draws if r["arm"] == "twin"}, {common.DATASET, common.OTHER})
        self.assertEqual(datasets[0]["species"], "Synthia coli")
        self.assertIn("checks:", out)
        self.assertNotIn("balanced", out)

    def test_an_arm_on_other_folds_fails(self):
        def change(r, d, i):
            if r[2] == "program" and r[5] == 43:
                shift_folds(d)

        rows, docs, files, releases = common.full_run(change=change)
        rc, _, checks, draws, _ = run(rows, docs, files, releases)
        self.assertEqual(rc, 1)
        self.assertEqual(failed(checks), [FOLDS])
        self.assertEqual(draws, [])

    def test_an_arm_is_held_to_the_default_at_its_own_seed(self):
        # The default's folds at seed 43 differ from seed 42's, and the arm at
        # 43 has seed 42's: a pairing with the wrong seed would pass it.
        def change(r, d, i):
            if r[5] == 43 and r[2] == "k31":
                shift_folds(d)

        rows, docs, files, releases = common.full_run(change=change)
        rc, _, checks, _, _ = run(rows, docs, files, releases)
        self.assertEqual(rc, 1)
        self.assertEqual(failed(checks), [FOLDS])

    def test_a_twin_at_a_draw_seed_fails(self):
        rows, docs, files, releases = common.full_run(twin=(45,))
        twin = [r for r in rows if r[2] == "twin"][0]
        moved = list(twin)
        moved[5], moved[7] = 43, f"twin/{twin[0]}.lineage-s43.pick.json"
        docs[moved[7]] = common.pick(dataset=twin[0], seed=43)
        rows = [r for r in rows if r is not twin] + [moved]
        rc, _, checks, _, _ = run(rows, docs, files, releases)
        self.assertEqual(rc, 1)
        self.assertEqual(failed(checks), ["a twin seed is the twin's and no draw's"])

    def test_a_twin_with_an_inputs_file_fails(self):
        rows, docs, files, releases = common.full_run()
        extra = {f"twin/{common.DATASET}.drawn.inputs.json": common.inputs(arm="program")}
        rc, _, checks, _, _ = run(rows, docs, files, releases, extra)
        self.assertEqual(rc, 1)
        self.assertEqual(failed(checks), ["the default and the twin write no inputs file"])

    def test_a_default_with_an_inputs_file_fails(self):
        rows, docs, files, releases = common.full_run()
        extra = {f"k31/{common.OTHER}.drawn.inputs.json": common.inputs(dataset=common.OTHER)}
        rc, _, checks, _, _ = run(rows, docs, files, releases, extra)
        self.assertEqual(rc, 1)
        self.assertEqual(failed(checks), ["the default and the twin write no inputs file"])

    def test_a_default_that_differs_from_the_release_fails(self):
        rows, docs, files, releases = common.full_run()
        releases[common.DATASET]["splits"][1]["learners"][0]["metrics"]["balanced_accuracy"] = 0.5
        rc, _, checks, _, _ = run(rows, docs, files, releases)
        self.assertEqual(rc, 1)
        self.assertEqual(failed(checks), ["seed 42: each learner's lineage figure is the release report's"])

    def test_a_default_with_no_release_report_fails(self):
        rows, docs, files, releases = common.full_run()
        del releases[common.OTHER]
        rc, _, checks, _, _ = run(rows, docs, files, releases)
        self.assertEqual(rc, 1)
        self.assertIn("seed 42: the default is held to the release report", failed(checks))

    def test_an_arm_result_is_not_held_to_the_release(self):
        # An arm's figures differ from the release's by design; only the default is held to it.
        def change(r, d, i):
            if r[2] == "program":
                for f in d["folds"]:
                    for o in f["outer"]:
                        o["tp"], o["fn"] = 0, o["tp"] + o["fn"]
                        o["tn"], o["fp"] = o["tn"] + o["fp"], 0
                d["learners"] = [{"learner": l, "balanced_accuracy": 0.5} for l in common.LEARNERS]

        rows, docs, files, releases = common.full_run(change=change)
        rc, _, checks, _, _ = run(rows, docs, files, releases)
        self.assertEqual(rc, 0, failed(checks))

    def test_an_arm_with_no_inputs_file_fails(self):
        rows, docs, files, releases = common.full_run()
        files = {k: v for k, v in files.items() if not k.startswith("k31-program/")}
        rc, _, checks, _, _ = run(rows, docs, files, releases)
        self.assertEqual(rc, 1)
        self.assertEqual(failed(checks), ["inputs: the arm's file exists"])

    def test_an_inputs_file_of_another_arm_or_tag_fails(self):
        for key, value in (("arm", "k31-program"), ("tag", "lineage"), ("dataset", common.OTHER),
                           ("evaluation_revision", 13)):
            with self.subTest(key=key):
                def change(r, d, i, key=key, value=value):
                    if i is not None and r[2] == "program":
                        i[key] = value

                rows, docs, files, releases = common.full_run(change=change)
                rc, _, checks, _, _ = run(rows, docs, files, releases)
                self.assertEqual(rc, 1)
                self.assertEqual(failed(checks), ["inputs: format, revision, dataset, arm and tag"])

    def test_a_joined_arm_that_drops_or_doubles_a_block_fails(self):
        for rows_held in (4, 100, 2 * 4 + 100):
            with self.subTest(arm_rows=rows_held):
                def change(r, d, i, n=rows_held):
                    if i is not None and r[2] == "k31-program":
                        i["arm_rows"] = n

                rows, docs, files, releases = common.full_run(change=change)
                rc, _, checks, _, _ = run(rows, docs, files, releases)
                self.assertEqual(rc, 1)
                self.assertEqual(failed(checks), ["inputs: the arm's rows are its blocks', each once"])

    def test_a_cohort_genome_with_no_table_read_fails(self):
        def change(r, d, i):
            if i is not None:
                i["tables_read"] = i["samples"] - 1

        rows, docs, files, releases = common.full_run(change=change)
        rc, _, checks, _, _ = run(rows, docs, files, releases)
        self.assertEqual(rc, 1)
        self.assertEqual(failed(checks), ["inputs: a table read for every cohort genome"])

    def test_program_names_that_do_not_add_up_fail(self):
        def change(r, d, i):
            if i is not None:
                i["program_names_seen"] += 1

        rows, docs, files, releases = common.full_run(change=change)
        rc, _, checks, _, _ = run(rows, docs, files, releases)
        self.assertEqual(rc, 1)
        self.assertEqual(failed(checks), ["inputs: the program names seen are those out of range, collapsed and kept"])

    def test_a_missing_result_fails(self):
        rows, docs, files, releases = common.full_run()
        del docs[f"program/{common.DATASET}.lineage-s44.pick.json"]
        rc, _, checks, _, _ = run(rows, docs, files, releases)
        self.assertEqual(rc, 1)
        self.assertEqual(failed(checks), ["exists"])

    def test_the_scores_hold_no_genome_or_sequence_field(self):
        rows, docs, files, releases = common.full_run()
        _, _, _, draws, datasets = run(rows, docs, files, releases)
        for t in (draws, datasets):
            self.assertFalse({k for k in t[0] if "genome" in k or "kmer" in k or "sequence" in k})


def count_rows(arm_shift, n=12, species=4, arms=score.FAMILY, study=True):
    """SYNTHETIC datasets.tsv rows: n datasets in `species` species, each arm's value the default's plus
    arm_shift(arm, learner, measure)."""
    out = []
    base = {"balanced_accuracy": 0.8, "recall": 0.7, "specificity": 0.9, "vme": 0.3, "me": 0.1, "log_loss": 0.4}
    for k in range(n):
        sp = f"synthia_{k % species}"
        d = f"{sp}__drug_{k}"
        for split in ("lineage", "study") if study else ("lineage",):
            for arm in ["k31"] + list(arms):
                if split == "study" and arm == "twin":
                    continue
                for learner in common.LEARNERS:
                    row = {"dataset": d, "species_slug": sp, "species": sp, "antibiotic": f"drug {k}", "split": split,
                           "arm": arm, "learner": learner, "draws": 3, "scored_draws": 3, "scored": 1}
                    for m, v in base.items():
                        row[m] = v + (0.0 if arm == "k31" else arm_shift(arm, learner, m))
                    out.append(row)
    return out


def count(rows):
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "datasets.tsv"
        score.write_tsv(path, rows, score.DATASET_COLUMNS)
        out = Path(tmp) / "count"
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            rc = score.main(["--count", str(out), str(path)])
        return rc, table(out / "count.tsv"), table(out / "conditions.tsv")


class Count(unittest.TestCase):
    def test_every_difference_is_the_arms_value_less_the_defaults(self):
        def shift(arm, learner, m):
            return {"balanced_accuracy": 0.05, "vme": -0.05}.get(m, 0.0) if arm == "program" else 0.0

        rc, rows, _ = count(count_rows(shift))
        self.assertEqual(rc, 0)
        at = {(r["split"], r["pair"], r["learner"], r["measure"]): r for r in rows}
        ba = at[("lineage", "program-k31", "rf", "balanced_accuracy")]
        self.assertEqual(float(ba["median_difference"]), 0.05)
        self.assertEqual((int(ba["better"]), int(ba["worse"])), (12, 0))
        # Lower is better for the very major error rate: the median keeps the arm's
        # value less the default's, and only better and worse turn round.
        vme = at[("lineage", "program-k31", "rf", "vme")]
        self.assertEqual(float(vme["median_difference"]), -0.05)
        self.assertEqual((int(vme["better"]), int(vme["worse"])), (12, 0))
        same = at[("lineage", "unitigs-k31", "rf", "balanced_accuracy")]
        self.assertEqual((int(same["better"]), int(same["worse"]), int(same["same"])), (0, 0, 12))

    def test_holm_runs_over_the_thirty_lineage_pairs_of_the_primary_measure(self):
        def shift(arm, learner, m):
            return 0.05 if arm == "program" and m == "balanced_accuracy" else 0.0

        rc, rows, _ = count(count_rows(shift))
        self.assertEqual(rc, 0)
        family = [r for r in rows if r["split"] == "lineage" and r["measure"] == "balanced_accuracy"]
        self.assertEqual(len(family), 30)
        self.assertEqual({r["pair"] for r in family}, {f"{a}-k31" for a in score.FAMILY})
        p = 0.5 ** 12
        prog = [r for r in family if r["pair"] == "program-k31"]
        self.assertTrue(all(float(r["sign_p"]) == p for r in prog))
        self.assertTrue(all(float(r["holm_p"]) == 30 * p for r in prog))
        self.assertTrue(all(r["holm_p"] == "" for r in rows if r not in family))

    def test_the_count_refuses_part_of_the_family(self):
        rows = count_rows(lambda a, l, m: 0.0, arms=[a for a in score.FAMILY if a != "unitigs"])
        with self.assertRaises(SystemExit) as e:
            count(rows)
        self.assertIn("unitigs", str(e.exception))

    def test_conditions_read_three_of_five_learners_and_the_twin(self):
        def shift(arm, learner, m):
            if arm == "k31-program" and m == "balanced_accuracy" and learner in ("rf", "gbt", "lasso"):
                return 0.05
            return 0.0

        rc, _, conditions = count(count_rows(shift))
        self.assertEqual(rc, 0)
        at = {(r["arm"], r["condition"]): r for r in conditions}
        self.assertEqual(at[("k31-program", "1")]["met"], "1")
        self.assertEqual(at[("k31-program", "1")]["learners"], "3")
        self.assertEqual(at[("k31-program", "2")]["met"], "1")
        self.assertEqual(at[("k31-program", "3")]["met"], "1")
        self.assertEqual(at[("k31-program", "5")]["met"], "1")
        self.assertEqual(at[("program", "1")]["met"], "0")
        self.assertEqual(at[("k31-program", "6")]["met"], "")

    def test_a_twin_that_moves_voids_condition_three(self):
        def shift(arm, learner, m):
            return 0.05 if arm in ("twin", "program") and m == "balanced_accuracy" else 0.0

        _, _, conditions = count(count_rows(shift))
        at = {(r["arm"], r["condition"]): r for r in conditions}
        self.assertEqual(at[("program", "1")]["met"], "1")
        self.assertEqual(at[("program", "3")]["met"], "0")

    def test_a_recall_interval_below_the_margin_fails_condition_four(self):
        def shift(arm, learner, m):
            return -0.05 if arm == "program" and m == "recall" and learner == "scm" else 0.0

        _, _, conditions = count(count_rows(shift))
        at = {(r["arm"], r["condition"]): r for r in conditions}
        self.assertEqual(at[("program", "4")]["met"], "0")
        self.assertEqual(at[("program", "4")]["learners"], "4")

    def test_the_species_interval_resamples_species(self):
        # Two species: one where the arm gains 0.1 on every dataset, one where it gains nothing.
        rows = count_rows(lambda arm, learner, m: 0.0, n=8, species=2)
        for r in rows:
            if r["arm"] == "program" and r["species_slug"] == "synthia_0":
                r["balanced_accuracy"] += 0.1
        _, counted, _ = count(rows)
        at = {(r["split"], r["pair"], r["learner"], r["measure"]): r for r in counted}
        r = at[("lineage", "program-k31", "scm", "balanced_accuracy")]
        # A species resample can draw synthia_0 twice or not at all, so the
        # interval spans 0 to 0.1; a dataset resample of 8 would rarely reach either end.
        self.assertEqual((float(r["interval_lo"]), float(r["interval_hi"])), (0.0, 0.1))
        self.assertEqual(int(r["species"]), 2)


if __name__ == "__main__":
    unittest.main()
