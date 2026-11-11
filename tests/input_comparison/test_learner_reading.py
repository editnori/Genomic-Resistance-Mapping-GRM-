"""Tests for scripts/input-comparison/learner_reading.py: the paired call counts and their rates, the pairing rule on
genomes only one arm scored, the split kept apart, the two weightings of a learner's share, the run's call rule, the
cohort lines, the sensitivity reading's recomputation, the command, and what it refuses.

Every folder here is written the way run.py writes a WORK: prediction rows per arm and dataset, and J's own records
from stats.pair_units, stats.compare_units and stats.finalize on those rows.

usage: python3 -B -m unittest discover -s tests/input_comparison -t .
"""
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

from tests.input_comparison import test_run as base

reading = base.load("learner_reading")
stats = reading.stats
DS = "species_a__drug"
N = 20  # genomes 0..7 resistant, 8..19 susceptible
TWIN = "twin"
TWINS = {TWIN: ["scm", "rf", "gbt"]}


def y_of(i):
    return int(i < 8)


def fold_of(repeat, i):
    return (i + repeat) % 2


def score_row(arm, learner, split, repeat, i, score, called=None, dataset=DS):
    return {"dataset": dataset, "arm": arm, "learner": learner, "split": split, "repeat": repeat,
            "sample_id": f"g{i:02d}", "y": y_of(i), "fold": fold_of(repeat, i), "score": score,
            "called": score >= 0.5 if called is None else called, "cluster": i,
            "fold_plan_sha256": f"{split}-plan-{repeat}"}


def flat(arm, learner, split, r_score, s_score, *, repeats=3, keep=lambda repeat, i: True, dataset=DS):
    """Rows with one score for every resistant genome and one for every susceptible one, nudged per genome so no
    two units are equal; `keep` picks the genomes the arm scores in each repeat."""
    return [score_row(arm, learner, split, repeat, i,
                      (r_score if y_of(i) else s_score) + 0.001 * ((i * 7 + repeat) % 5), dataset=dataset)
            for repeat in range(repeats) for i in range(N) if keep(repeat, i)]


def write_work(root, arms, *, species="Species a", commit="c" * 40, default=None, dataset=DS):
    """A WORK folder: `arms` maps an arm to (pairing, rows); `default` is the default's rows (a calm scorer on
    every learner and split by default). J's records are made from the rows as run.pending_contrasts makes them,
    the twin's family C judged on its own p."""
    root = Path(root)
    default = default if default is not None else [
        r for learner in reading.LEARNERS for split in reading.SPLITS
        for r in flat("default", learner, split, 0.7, 0.3, dataset=dataset)]
    rows = {"default": default, **{arm: r for arm, (_, r) in arms.items()}}
    for arm, arm_rows in rows.items():
        out = root / "arms" / arm / dataset
        out.mkdir(parents=True)
        (out / "predictions.jsonl").write_text("".join(json.dumps(r) + "\n" for r in arm_rows))
        (out / "arm.json").write_text(json.dumps({"samples": len({r["sample_id"] for r in arm_rows}) + 5}))
    (root / "arms" / "default" / dataset / "plans").mkdir()
    for name in ("lineage-0.json", "lineage-1.json", "refused-lineage-2.json"):
        (root / "arms" / "default" / dataset / "plans" / name).write_text(json.dumps(
            {"sample_ids": [f"g{i:02d}" for i in range(N)], "y": [y_of(i) for i in range(N)]}
            if not name.startswith("refused") else {"refused": True}))
    pending, planned, missing, unadjusted = [], [], {}, []
    every = default + [r for _, r in arms.values() for r in r]
    for arm, (pairing, _) in arms.items():
        for learner in reading.LEARNERS:
            for split in reading.SPLITS:
                family = "C" if split == "lineage" else "D"
                pid = f"{dataset}/{stats.contrast_key(family, learner, split, 'default', arm, 'log_loss')}"
                planned.append(pid)
                if arm == TWIN and family == "C":
                    unadjusted.append(pid)
                try:
                    paired = stats.pair_units(every, dataset=dataset, learner=learner, split=split,
                                              baseline="default", treatment=arm, pairing=pairing)
                except stats.PairingError as e:
                    missing[pid] = f"pairing refused: {e}"
                    continue
                pending += stats.compare_units(paired, family=family, margins={"log_loss": 0.02}, development=[],
                                               test=stats.LEARNER_TESTS[learner])
    records = stats.finalize(pending, planned, missing, unadjusted)
    (root / "j").mkdir()
    for name, value in (("contrasts", records), ("planned", planned), ("missing", missing)):
        (root / "j" / f"{name}.json").write_text(json.dumps(value))
    (root / "k").mkdir()
    manifest = {"datasets": [{"id": dataset, "species": species, "partition": "confirmatory"}],
                "sentinel_arms": [TWIN], "sentinel_learners": TWINS, "sentinel_adjust": "none",
                "source_commit": commit, "binary_sha256": "b" * 64, "synthetic": True,
                "runs": [{"arm": "default", "dataset": dataset, "learner": learner, "split": "lineage", "status": "ok"}
                         for learner in reading.LEARNERS]}
    (root / "k" / "manifest.json").write_text(json.dumps(manifest))
    # The report as second_reading reads it, cut down to the keys it takes (test_second_reading.report).
    (root / "report.json").write_text(json.dumps({
        "contrasts": records, "planned_contrasts": planned, "undefined_contrasts": [], "unreliable": [],
        "not_produced": [{"id": pid, "reason": reason} for pid, reason in missing.items()], "manifest": manifest}))
    return root


def all_learners(arm, r_score, s_score, **kw):
    return [r for learner in reading.LEARNERS for split in reading.SPLITS
            for r in flat(arm, learner, split, r_score, s_score, **kw)]


def quiet(tmp, name="w", dataset=DS, **kw):
    """A folder with one same_plan arm that sharpens every learner and a twin that changes nothing."""
    return write_work(Path(tmp) / name, {"budget": ("same_plan", all_learners("budget", 0.8, 0.2, dataset=dataset)),
                                         TWIN: ("same_plan", all_learners(TWIN, 0.7, 0.3, dataset=dataset))},
                      dataset=dataset, **kw)


def run_reading(*folders, sensitivity=False):
    read = [reading.read_folder(f, sensitivity) for f in folders]
    reports = [json.loads((Path(f) / "report.json").read_text()) for f in folders] if sensitivity else None
    return reading.read(read, reports)


class CallCountCases(unittest.TestCase):
    def test_paired_counts_and_rates_of_known_calls(self):
        # scm on the lineage split: the default calls every resistant genome (0.7) and every susceptible one too
        # (0.6); the arm misses resistant genomes 0..2 (0.4) and calls susceptible 8..10 (0.9). Per repeat: default
        # TP 8 FN 0 FP 12 TN 0, arm TP 5 FN 3 FP 3 TN 9; three repeats.
        def arm_score(repeat, i):
            return 0.4 if i < 3 else 0.7 if i < 8 else 0.9 if i < 11 else 0.1
        default = [score_row("default", "scm", "lineage", r, i, 0.7 if y_of(i) else 0.6)
                   for r in range(3) for i in range(N)]
        arm = [score_row("budget", "scm", "lineage", r, i, arm_score(r, i)) for r in range(3) for i in range(N)]
        with tempfile.TemporaryDirectory() as tmp:
            got = run_reading(write_work(tmp, {"budget": ("same_plan", arm)}, default=default))
        t = got["tables"]["lineage"]["budget"]["scm"]
        self.assertEqual(t["default"], {"tp": 24, "fn": 0, "fp": 36, "tn": 0})
        self.assertEqual(t["arm"], {"tp": 15, "fn": 9, "fp": 9, "tn": 27})
        self.assertEqual(t["rates_default"], {"recall": 1.0, "specificity": 0.0, "balanced_accuracy": 0.5,
                                              "precision": 0.4, "f1": 48 / 84})
        self.assertEqual(t["rates_arm"], {"recall": 15 / 24, "specificity": 27 / 36,
                                          "balanced_accuracy": (15 / 24 + 27 / 36) / 2, "precision": 15 / 24,
                                          "f1": 30 / 48})
        self.assertAlmostEqual(t["rates_change"]["recall"], 15 / 24 - 1.0, places=15)
        # Log loss improves (the susceptible genomes leave 0.6) while recall falls: the disagreement count.
        self.assertGreater(t["mean"], 0)
        self.assertEqual((t["improves_recall_falls"], t["worse_recall_rises"], t["recall_defined"]), (1, 0, 1))
        (row,) = [r for r in got["rows"] if r["split"] == "lineage"]
        self.assertEqual((row["repeats"], row["paired_calls"]), (3, 60))

    def test_worse_while_recall_rises_and_an_undefined_rate(self):
        # The arm calls every genome resistant at about 0.95: recall rises from 0 to 1 and log loss gets worse.
        default = flat("default", "rf", "lineage", 0.3, 0.2)
        arm = flat("budget", "rf", "lineage", 0.95, 0.95)
        with tempfile.TemporaryDirectory() as tmp:
            got = run_reading(write_work(tmp, {"budget": ("same_plan", arm)}, default=default))
        t = got["tables"]["lineage"]["budget"]["rf"]
        self.assertLess(t["mean"], 0)
        self.assertEqual((t["improves_recall_falls"], t["worse_recall_rises"], t["recall_defined"]), (0, 1, 1))
        self.assertEqual(reading.rates({"tp": 0, "fn": 0, "fp": 2, "tn": 1})["recall"], None)
        self.assertEqual(reading.rates({"tp": 0, "fn": 0, "fp": 0, "tn": 0}),
                         dict.fromkeys(reading.RATES, None))

    def test_the_call_is_the_rows_called_and_not_the_score(self):
        # A cohort arm's no-call row: no model, scored at its own training prevalence 0.7, called false. The run
        # counts it as a susceptible call; a rule of score >= 0.5 would call it resistant.
        arm = [score_row("cap", "scm", "lineage", r, i, 0.7 if i == 0 else 0.9 if y_of(i) else 0.1,
                         called=False if i == 0 else None) for r in range(3) for i in range(N)]
        default = [score_row("default", "scm", "lineage", r, i, 0.9 if y_of(i) else 0.1)
                   for r in range(3) for i in range(N)]
        with tempfile.TemporaryDirectory() as tmp:
            got = run_reading(write_work(tmp, {"cap": ("common_outer_test_set", arm)}, default=default))
        t = got["tables"]["lineage"]["cap"]["scm"]
        self.assertEqual(t["arm"], {"tp": 21, "fn": 3, "fp": 0, "tn": 36})
        self.assertEqual(t["default"], {"tp": 24, "fn": 0, "fp": 0, "tn": 36})


class PairingCases(unittest.TestCase):
    def test_a_genome_only_the_default_scored_is_left_out(self):
        # The cohort arm scores genomes 0..5 and 10..19 in repeat 0 and every genome after: the default's 6..9 in
        # repeat 0 aren't paired, so its counts there come from the 16 genomes both scored.
        keep = lambda repeat, i: repeat > 0 or not 6 <= i < 10  # noqa: E731
        arm = flat("cap", "rf", "lineage", 0.9, 0.1, keep=keep)
        default = flat("default", "rf", "lineage", 0.9, 0.6)
        with tempfile.TemporaryDirectory() as tmp:
            got = run_reading(write_work(tmp, {"cap": ("common_outer_test_set", arm)}, default=default))
        t = got["tables"]["lineage"]["cap"]["rf"]
        # Repeat 0: resistant 0..5, susceptible 10..19; repeats 1 and 2: 8 and 12.
        self.assertEqual(t["default"], {"tp": 6 + 16, "fn": 0, "fp": 10 + 24, "tn": 0})
        self.assertEqual(t["arm"], {"tp": 22, "fn": 0, "fp": 0, "tn": 34})
        (row,) = got["rows"]
        self.assertEqual(row["paired_calls"], 16 + 40)

    def test_a_genome_only_the_arm_scored_refuses_the_folder(self):
        arm = flat("cap", "rf", "lineage", 0.9, 0.1, keep=lambda repeat, i: i != 3)
        default = flat("default", "rf", "lineage", 0.9, 0.6, keep=lambda repeat, i: i != 3)
        with tempfile.TemporaryDirectory() as tmp:
            work = write_work(tmp, {"cap": ("common_outer_test_set", arm)}, default=default)
            path = work / "arms" / "cap" / DS / "predictions.jsonl"
            extra = score_row("cap", "rf", "lineage", 0, 3, 0.9)
            path.write_text(path.read_text() + json.dumps(extra) + "\n")
            with self.assertRaisesRegex(SystemExit, r"rows don't pair: .*samples only in 'cap'"):
                run_reading(work)

    def test_a_same_plan_arm_missing_a_genome_refuses_the_folder(self):
        with tempfile.TemporaryDirectory() as tmp:
            work = quiet(tmp)
            path = work / "arms" / "budget" / DS / "predictions.jsonl"
            lines = path.read_text().splitlines(keepends=True)
            path.write_text("".join(line for line in lines if '"sample_id": "g04"' not in line
                                    or '"learner": "gbt"' not in line))
            with self.assertRaisesRegex(SystemExit, r"rows don't pair: .*samples only in 'default'"):
                run_reading(work)

    def test_a_contrast_j_did_not_make_is_left_out(self):
        # The arm has no gbt rows: J's pairing refuses, the contrast is in missing.json, and no row reads it.
        arm = [r for r in all_learners("budget", 0.8, 0.2) if r["learner"] != "gbt"]
        with tempfile.TemporaryDirectory() as tmp:
            work = write_work(tmp, {"budget": ("same_plan", arm)})
            got = run_reading(work)
            missing = json.loads((work / "j" / "missing.json").read_text())
        self.assertEqual(len(missing), 2)
        self.assertNotIn("gbt", got["tables"]["lineage"]["budget"])
        self.assertEqual(sorted(got["tables"]["lineage"]["budget"]), ["cart", "lasso", "rf", "scm"])


class SplitAndWeightCases(unittest.TestCase):
    def test_the_random_split_is_kept_apart(self):
        # The arm is sharp on lineage and blunt on random: the two tables hold the two.
        arm = flat("budget", "cart", "lineage", 0.9, 0.1) + flat("budget", "cart", "random", 0.45, 0.55)
        default = flat("default", "cart", "lineage", 0.7, 0.3) + flat("default", "cart", "random", 0.7, 0.3)
        with tempfile.TemporaryDirectory() as tmp:
            got = run_reading(write_work(tmp, {"budget": ("same_plan", arm)}, default=default))
        lineage, rand = got["tables"]["lineage"]["budget"]["cart"], got["tables"]["random"]["budget"]["cart"]
        self.assertEqual(lineage["arm"], {"tp": 24, "fn": 0, "fp": 0, "tn": 36})
        self.assertEqual(rand["arm"], {"tp": 0, "fn": 24, "fp": 36, "tn": 0})
        self.assertEqual((lineage["better"], lineage["worse"], rand["better"], rand["worse"]), (1, 0, 0, 1))
        self.assertEqual(list(got["tables"]), ["lineage", "random"])

    def test_the_two_weightings_of_a_share(self):
        # A cohort arm that scores 5 genomes of repeat 0 (fold 0: 0 and 8; fold 1: 1, 3 and 9) and all 20 after.
        # Units: every (repeat, fold) weighs the same, so repeat 0's two small, unequal units count as much as the
        # others. Genomes: repeat 0's mean over its 5 genomes, then the mean of the three repeats.
        keep = lambda repeat, i: repeat > 0 or i in (0, 1, 3, 8, 9)  # noqa: E731
        arms = {"cap": ("common_outer_test_set", [
            *flat("cap", "scm", "lineage", 0.9, 0.1, keep=keep),
            *flat("cap", "lasso", "lineage", 0.6, 0.4, keep=keep)])}
        default = flat("default", "scm", "lineage", 0.6, 0.4) + flat("default", "lasso", "lineage", 0.6, 0.4)
        with tempfile.TemporaryDirectory() as tmp:
            work = write_work(tmp, arms, default=default)
            got = run_reading(work)
            recorded = {x["learner"]: x["estimate"] for x in json.loads((work / "j" / "contrasts.json").read_text())
                        if x["estimate"] is not None}
        by = {r["learner"]: r for r in got["rows"]}
        rows = {learner: [r for r in arms["cap"][1] if r["learner"] == learner] for learner in ("scm", "lasso")}

        def delta(r):
            d = next(x for x in default if (x["learner"], x["repeat"], x["sample_id"])
                     == (r["learner"], r["repeat"], r["sample_id"]))
            y = np.array([r["y"]])
            return float(stats.sample_log_loss(y, np.array([d["score"]]))[0]
                         - stats.sample_log_loss(y, np.array([r["score"]]))[0])

        def mean(xs):
            return sum(xs) / len(xs)

        for learner in ("scm", "lasso"):
            units = [mean([delta(r) for r in rows[learner] if (r["repeat"], r["fold"]) == u])
                     for u in sorted({(r["repeat"], r["fold"]) for r in rows[learner]})]
            genomes = [mean([delta(r) for r in rows[learner] if r["repeat"] == repeat]) for repeat in range(3)]
            self.assertAlmostEqual(by[learner]["estimate"], mean(units), places=12)
            self.assertEqual(by[learner]["estimate"], recorded[learner])  # the record's, not a recomputation
            self.assertAlmostEqual(by[learner]["genome_weighted"], mean(genomes), places=12)
        self.assertNotAlmostEqual(by["scm"]["estimate"], by["scm"]["genome_weighted"], places=6)
        share = {learner: got["tables"]["lineage"]["cap"][learner]["share"] for learner in ("scm", "lasso")}
        units_total = by["scm"]["estimate"] + by["lasso"]["estimate"]
        genomes_total = by["scm"]["genome_weighted"] + by["lasso"]["genome_weighted"]
        self.assertAlmostEqual(share["scm"]["units"], by["scm"]["estimate"] / units_total, places=12)
        self.assertAlmostEqual(share["scm"]["genomes"], by["scm"]["genome_weighted"] / genomes_total, places=12)
        self.assertAlmostEqual(share["scm"]["units"] + share["lasso"]["units"], 1.0, places=12)
        self.assertEqual(got["pooled"]["lineage"]["cap"]["contrasts"], 2)
        self.assertAlmostEqual(got["pooled"]["lineage"]["cap"]["genomes"], genomes_total, places=12)

    def test_sign_counts_and_twins(self):
        with tempfile.TemporaryDirectory() as tmp:
            got = run_reading(quiet(tmp))
        budget, twin = got["tables"]["lineage"]["budget"], got["tables"]["lineage"][TWIN]
        self.assertEqual({learner: (t["better"], t["worse"], t["datasets"]) for learner, t in budget.items()},
                         dict.fromkeys(reading.LEARNERS, (1, 0, 1)))
        # Log loss improves and recall stays at 1: a tie, which is no fall.
        self.assertEqual({(t["rates_change"]["recall"], t["improves_recall_falls"], t["recall_defined"])
                          for t in budget.values()}, {(0.0, 0, 1)})
        # The twin's rows are the default's again, so every estimate is exactly zero.
        self.assertEqual({(t["unchanged"], t["improves_recall_falls"], t["twin"]) for t in twin.values()},
                         {(1, 0, True)})
        self.assertTrue(all(t["share"] == {"units": None, "genomes": None} for t in twin.values()))


class CohortCases(unittest.TestCase):
    def test_a_cohort_arm_has_its_default_class_counts_and_no_mix_of_its_own(self):
        with tempfile.TemporaryDirectory() as tmp:
            arms = {"cap": ("common_outer_test_set", flat("cap", "rf", "lineage", 0.9, 0.1)),
                    "budget": ("same_plan", flat("budget", "rf", "lineage", 0.9, 0.1))}
            got = run_reading(write_work(tmp, arms, default=flat("default", "rf", "lineage", 0.7, 0.3)))
        self.assertEqual(got["cohort"], [{"species": "Species a", "dataset": DS, "arm": "cap", "default_genomes": 25,
                                          "default_resistant": 8, "default_susceptible": 12, "arm_genomes": 25,
                                          "arm_class_mix": "not recorded"}])

    def test_plan_files_that_disagree_on_the_labels_refuse(self):
        with tempfile.TemporaryDirectory() as tmp:
            work = write_work(tmp, {"cap": ("common_outer_test_set", flat("cap", "rf", "lineage", 0.9, 0.1))},
                              default=flat("default", "rf", "lineage", 0.7, 0.3))
            (work / "arms" / "default" / DS / "plans" / "lineage-1.json").write_text(json.dumps(
                {"sample_ids": [f"g{i:02d}" for i in range(N)], "y": [1] * N}))
            with self.assertRaisesRegex(SystemExit, "plan files disagree"):
                run_reading(work)


class SensitivityCases(unittest.TestCase):
    def test_the_registered_tests_give_the_second_reading_and_one_test_moves_it(self):
        # The arm improves every learner on every unit by a similar amount. The twin judges scm, rf and gbt and
        # changes nothing, so none of them is flagged.
        with tempfile.TemporaryDirectory() as tmp:
            work = quiet(tmp)
            got = run_reading(work, sensitivity=True)["sensitivity"]
            report = json.loads((work / "report.json").read_text())
        second = reading.second_reading.read([report])["panel"]
        self.assertEqual(got["registered"]["firm"], second["firm"])
        self.assertEqual(got["registered"]["decisions"], second["decisions"])
        self.assertIn("decides nothing", got["caveat"])
        self.assertEqual(got["tests"]["corrected_t"], dict.fromkeys(reading.LEARNERS, stats.UNIT_CORRECTED_T))
        # scm, rf and gbt are judged and improve under every test; cart and the lasso have no twin, so no
        # decision of theirs is firm.
        for name in ("registered", "paired_t", "corrected_t"):
            self.assertEqual(got[name]["firm_by_learner"], {"scm": 1, "cart": 0, "rf": 1, "gbt": 1, "lasso": 0}, name)
        # The uniform test reaches the records: the corrected t's intervals on scm are wider than the paired t's.
        width = {}
        for name in ("paired_t", "corrected_t"):
            (d,) = [d for d in got[name]["decisions"] if ":scm:" in d["id"]]
            width[name] = d["ci"][1] - d["ci"][0]
        self.assertGreater(width["corrected_t"], width["paired_t"])

    def test_a_report_that_isnt_the_folders_refuses(self):
        with tempfile.TemporaryDirectory() as tmp:
            work = quiet(tmp)
            report = json.loads((work / "report.json").read_text())
            next(x for x in report["contrasts"] if x["p"] is not None)["p_adj"] = 0.999
            (work / "report.json").write_text(json.dumps(report))
            with self.assertRaisesRegex(SystemExit, "aren't its report's"):
                run_reading(work, sensitivity=True)


class RefusalCases(unittest.TestCase):
    def test_a_record_the_rows_dont_give_back_refuses(self):
        with tempfile.TemporaryDirectory() as tmp:
            work = quiet(tmp)
            path = work / "j" / "contrasts.json"
            records = json.loads(path.read_text())
            x = next(x for x in records if x["treatment"] == "budget" and x["learner"] == "rf")
            x["estimate"] += 1e-6
            path.write_text(json.dumps(records))
            with self.assertRaisesRegex(SystemExit, "the rows give estimate"):
                run_reading(work)

    def test_a_planned_contrast_neither_made_nor_missing_refuses(self):
        with tempfile.TemporaryDirectory() as tmp:
            work = quiet(tmp)
            path = work / "j" / "contrasts.json"
            path.write_text(json.dumps(json.loads(path.read_text())[1:]))
            with self.assertRaisesRegex(SystemExit, "doesn't account for every planned contrast"):
                run_reading(work)

    def test_a_missing_file_refuses(self):
        with tempfile.TemporaryDirectory() as tmp:
            work = quiet(tmp)
            (work / "arms" / "budget" / DS / "predictions.jsonl").unlink()
            with self.assertRaisesRegex(SystemExit, "predictions.jsonl is missing"):
                run_reading(work)
            (work / "j" / "planned.json").write_text("{not json")
            with self.assertRaisesRegex(SystemExit, "planned.json is not JSON"):
                run_reading(work)

    def test_folders_that_arent_one_comparison_refuse(self):
        with tempfile.TemporaryDirectory() as tmp:
            a = quiet(tmp, "a")
            with self.assertRaisesRegex(SystemExit, "these have more: Species a"):
                run_reading(a, quiet(tmp, "b"))
            with self.assertRaisesRegex(SystemExit, "differ in source_commit"):
                run_reading(a, quiet(tmp, "c", species="Species b", commit="d" * 40))
            mixed = quiet(tmp, "d")
            manifest = json.loads((mixed / "k" / "manifest.json").read_text())
            manifest["datasets"].append({"id": "species_b__drug", "species": "Species b", "partition": "confirmatory"})
            (mixed / "k" / "manifest.json").write_text(json.dumps(manifest))
            with self.assertRaisesRegex(SystemExit, "holds Species a, Species b"):
                run_reading(mixed)


class CommandCases(unittest.TestCase):
    def test_the_command_writes_the_reading_and_prints_its_text(self):
        with tempfile.TemporaryDirectory() as tmp:
            a = quiet(tmp, "a")
            b = quiet(tmp, "b", species="Species b", dataset="species_b__drug")
            out = Path(tmp) / "reading.json"
            p = subprocess.run([sys.executable, "-B", str(base.HERE / "learner_reading.py"), str(a), str(b),
                                "--json", str(out)], capture_output=True, text=True, check=True)
            want = run_reading(a, b)
            self.assertEqual(json.loads(out.read_text()), json.loads(json.dumps(want)))
            self.assertEqual(p.stdout, reading.render(want) + "\n")
        lines = p.stdout.splitlines()
        self.assertIn("decides nothing", lines[0])
        self.assertIn("== lineage split (decides) ==", lines)
        self.assertIn("== random split (kept apart; decides nothing) ==", lines)
        budget = [line for line in lines if line.startswith("budget") and " scm " in line]
        self.assertTrue(budget[0].split()[2:6] == ["2", "2", "0", "0"], budget[0])  # n, better, worse, zero
        self.assertTrue(any(line.startswith("budget") and "pooled" in line and "10 contrasts" in line
                            for line in lines))

    def test_the_text_labels_the_sensitivity_reading(self):
        with tempfile.TemporaryDirectory() as tmp:
            work = quiet(tmp)
            text = reading.render(run_reading(work, sensitivity=True))
            self.assertIsNone(run_reading(work)["sensitivity"])
        self.assertIn(f"SENSITIVITY READING ({reading.CAVEAT})", text)
        self.assertIn("corrected_t  firm", text)


if __name__ == "__main__":
    unittest.main()
