"""scripts/c6-mic/compare.py, the reading of the MIC comparison, on SYNTHETIC score rows (docs/mic-plan.md).

Every row here is SYNTHETIC: built in the test, in the columns score.py
writes, for made-up species and drugs. No count comes from a real run.
"""

import contextlib
import importlib.util
import io
import os
import tempfile
import unittest
from fractions import Fraction
from unittest import mock

from tests.c6_mic import common

compare = common.load("compare")
score = common.load("score")

STUDY_PLANS = {"42": "study-whole", "43": "study-d1", "44": "study-d2",
               "45": "study-d3", "46": "study-d4", "47": "c6mic-s47"}


def load_record():
    path = common.ROOT / "scripts" / "amrfinderplus-comparison" / "record.py"
    spec = importlib.util.spec_from_file_location("c6_mic_test_record", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def plan(split, seed):
    return f"lineage-s{seed}" if split == "lineage" else STUDY_PLANS[seed]


def score_row(slug, dataset, split, seed, arm, ea=(50, 100), unavailable="", folds_undefined=0,
              category=(80, 120), very_major=(5, 40), major=(4, 80)):
    """One SYNTHETIC scores.tsv row. `ea` and the three beside it are (count, denominator)."""
    blank = bool(unavailable)
    row = {
        "dataset": dataset, "species_slug": slug, "antibiotic": "drug " + dataset.split("__")[1],
        "split": split, "plan": plan(split, seed), "seed": seed,
        "draw": "draw" if seed in ("42", "43", "44") else "twin", "arm": arm,
        "unavailable": unavailable, "folds": 5, "folds_scored": 5, "folds_undefined": folds_undefined,
        "folds_without_breakpoint": 0,
        "ea_exact_reads": "MIC values only", "category_reads": "R and S labels and the training breakpoint",
    }
    for name, (count, n) in [("ea_exact", ea), ("category", category),
                             ("very_major", very_major), ("major", major)]:
        rate = {"ea_exact": "ea_exact", "category": "category_agreement",
                "very_major": "very_major_rate", "major": "major_rate"}[name]
        row[f"{name}_count"] = "" if blank else count
        row[f"{name}_n"] = n
        row[rate] = "" if blank or n == 0 else f"{count / n:.6f}"
    return row


def write_species(root, slug, sets, *, extra=None, drop=None, check_not_ok=None):
    """A SYNTHETIC species run under root/slug: expect.tsv beside score/.
    `sets` maps (dataset, split) to a function (seed, arm) -> score_row keyword
    overrides. `extra` adds a row once more, `drop` leaves one out, both as
    (dataset, split, seed, arm). `check_not_ok` marks a check row (dataset,
    split, seed) not ok. Returns the score folder."""
    top = os.path.join(root, slug)
    os.makedirs(os.path.join(top, "score"))
    expect, checks, scores = [], [], []
    for (dataset, split), over in sets.items():
        for seed in compare.SEEDS:
            arms = compare.EXPECT_ARMS[seed]
            file = f"{split}/{dataset}.{plan(split, seed)}.mic.json"
            expect.append([dataset, slug, split, plan(split, seed), seed, arms, file])
            ok = "no" if check_not_ok == (dataset, split, seed) else "yes"
            checks.append([dataset, split, plan(split, seed), seed, file, "yes", "ok", "equal", "ok",
                           "equal", "equal", "whole", ok])
            for arm in ["released", "ordered", "median"] if seed in compare.DRAW_SEEDS else ["released", "median"]:
                key = (dataset, split, seed, arm)
                if key == drop:
                    continue
                row = score_row(slug, dataset, split, seed, arm, **over(seed, arm))
                scores.append([row[c] for c in compare.SCORE_COLUMNS])
                if key == extra:
                    scores.append([row[c] for c in compare.SCORE_COLUMNS])
    for path, cols, rows in [(os.path.join(top, "expect.tsv"), compare.EXPECT_COLUMNS, expect),
                             (os.path.join(top, "score", "checks.tsv"), compare.CHECK_COLUMNS, checks),
                             (os.path.join(top, "score", "scores.tsv"), compare.SCORE_COLUMNS, scores)]:
        with open(path, "w") as f:
            f.write("\t".join(cols) + "\n")
            for r in rows:
                f.write("\t".join(str(c) for c in r) + "\n")
    return os.path.join(top, "score")


def ea_by_arm(released, median, ordered=None, twin=None):
    """Overrides giving each arm an EA-exact count over 100 at every draw.
    Each argument is one count or a tuple of three, one per seed."""
    def per_seed(value, seed, seeds):
        return value[seeds.index(seed)] if isinstance(value, tuple) else value

    def over(seed, arm):
        if seed in compare.DRAW_SEEDS:
            value = {"released": released, "median": median,
                     "ordered": released if ordered is None else ordered}[arm]
            return {"ea": (per_seed(value, seed, compare.DRAW_SEEDS), 100)}
        value = released if twin is None else twin
        return {"ea": (per_seed(value, seed, compare.TWIN_SEEDS) if arm == "released" else 50, 100)}
    return over


def run(species, *flags):
    """compare.main on SYNTHETIC species: (exit code, stdout, stderr, out dir files)."""
    with tempfile.TemporaryDirectory() as tmp:
        dirs = [write_species(os.path.join(tmp, "in"), slug, sets, **kw) for slug, sets, kw in species]
        out = os.path.join(tmp, "out")
        so, se = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(so), contextlib.redirect_stderr(se):
            code = compare.main(["--out", out, *flags, *dirs])
        files = {}
        if os.path.isdir(out):
            for name in sorted(os.listdir(out)):
                with open(os.path.join(out, name)) as f:
                    files[name] = f.read()
        return code, so.getvalue(), se.getvalue(), files


def tsv(text):
    lines = text.rstrip("\n").split("\n")
    head = lines[0].split("\t")
    return [dict(zip(head, line.split("\t"))) for line in lines[1:]]


def pair(files, split, name):
    return next(r for r in tsv(files["pairs.tsv"]) if r["split"] == split and r["pair"] == name)


def edge_sets():
    """SYNTHETIC lineage datasets around the 0.02 edge of P1, released minus median.
    The median arm reads 49/100 at every draw."""
    return {
        ("synthia_coli__drug_a", "lineage"): ea_by_arm(51, 49),          # +0.02 exactly: same
        ("synthia_coli__drug_b", "lineage"): ea_by_arm(47, 49),          # -0.02 exactly: same
        ("synthia_coli__drug_c", "lineage"): ea_by_arm(52, 49),          # +0.03: better
        ("synthia_coli__drug_d", "lineage"): ea_by_arm(46, 49),          # -0.03: worse
        ("synthia_coli__drug_e", "lineage"): ea_by_arm((51, 51, 52), 49),  # +7/300: better
    }


class Columns(unittest.TestCase):
    def test_the_columns_are_the_scorers(self):
        self.assertEqual(compare.SCORE_COLUMNS, score.SCORE_COLUMNS)
        self.assertEqual(compare.CHECK_COLUMNS, score.CHECK_COLUMNS)
        self.assertEqual(compare.EXPECT_COLUMNS, score.EXPECT_COLUMNS)


class Counts(unittest.TestCase):
    def test_the_classes_at_the_002_edge_are_exact(self):
        # 0.51 - 0.49 is 0.020000000000000018 in floats, so a float reading would call drug_a better.
        self.assertGreater(0.51 - 0.49, 0.02)
        code, _, _, files = run([("synthia_coli", edge_sets(), {})])
        self.assertEqual(code, 0)
        p1 = pair(files, "lineage", "P1")
        self.assertEqual((p1["compared"], p1["better"], p1["same"], p1["worse"]), ("5", "2", "2", "1"))
        cls = {r["dataset"]: r["P1_class"] for r in tsv(files["datasets.tsv"])}
        self.assertEqual(cls, {"synthia_coli__drug_a": "same", "synthia_coli__drug_b": "same",
                               "synthia_coli__drug_c": "better", "synthia_coli__drug_d": "worse",
                               "synthia_coli__drug_e": "better"})

    def test_a_value_is_the_mean_of_its_three_draws_and_the_twin_reads_45_to_47(self):
        sets = {("synthia_coli__drug_a", "lineage"): ea_by_arm((60, 61, 65), 49, twin=(70, 70, 71))}
        code, _, _, files = run([("synthia_coli", sets, {})])
        self.assertEqual(code, 0)
        (row,) = tsv(files["datasets.tsv"])
        self.assertEqual(row["released_ea_exact"], "%.6f" % (Fraction(186, 300)))
        self.assertEqual((row["released_ea_count"], row["released_ea_n"]), ("186", "300"))
        self.assertEqual(row["twin_ea_exact"], "%.6f" % (Fraction(211, 300)))
        self.assertEqual(row["P3_difference"], "%.6f" % (Fraction(211 - 186, 300)))
        self.assertEqual(row["mic_genomes"], "120")


class LeftOut(unittest.TestCase):
    def test_an_undefined_draw_leaves_the_dataset_out_of_its_pairs_only(self):
        base = ea_by_arm(60, 49)

        def over(seed, arm):
            o = base(seed, arm)
            if (seed, arm) == ("43", "ordered"):
                o["unavailable"] = "fold 2's training genomes give an empty fitting or holdout part"
            return o

        sets = {("synthia_coli__drug_a", "lineage"): over,
                ("synthia_coli__drug_b", "lineage"): ea_by_arm(60, 49)}
        code, out, _, files = run([("synthia_coli", sets, {})])
        self.assertEqual(code, 0)
        p2 = pair(files, "lineage", "P2")
        self.assertEqual((p2["datasets"], p2["compared"], p2["left_out"]), ("2", "1", "1"))
        self.assertEqual(pair(files, "lineage", "P1")["compared"], "2")
        self.assertEqual(pair(files, "lineage", "P3")["compared"], "2")
        row = next(r for r in tsv(files["datasets.tsv"]) if r["dataset"] == "synthia_coli__drug_a")
        self.assertEqual((row["P2_difference"], row["P2_class"]), ("", "left out"))
        self.assertEqual(row["P2_left_out"], "ordered: unavailable at seed 43")
        self.assertEqual(row["ordered_ea_exact"], "")
        self.assertIn("lineage P2 synthia_coli__drug_a: ordered: unavailable at seed 43", out)

    def test_an_undefined_fold_or_no_exact_mic_leaves_it_out(self):
        base = ea_by_arm(60, 49)

        def over(seed, arm):
            o = base(seed, arm)
            if (seed, arm) == ("46", "released"):
                o["folds_undefined"] = 1
            if (seed, arm) == ("44", "median"):
                o["ea"] = (0, 0)
            return o

        code, _, _, files = run([("synthia_coli", {("synthia_coli__drug_a", "lineage"): over}, {})])
        self.assertEqual(code, 0)
        (row,) = tsv(files["datasets.tsv"])
        self.assertEqual(row["P3_left_out"], "twin: 1 fold(s) undefined at seed 46")
        self.assertEqual(row["P1_left_out"], "median: no exact MIC tested at seed 44")
        self.assertEqual(row["P2_class"], "same")


class Stops(unittest.TestCase):
    def check_stop(self, kw, words, *flags):
        code, out, err, files = run([("synthia_coli", edge_sets(), kw)], *flags)
        self.assertEqual(code, 2)
        self.assertIn(words, err)
        self.assertEqual((out, files), ("", {}))

    def test_a_missing_row_stops_the_reading(self):
        self.check_stop({"drop": ("synthia_coli__drug_c", "lineage", "45", "median")}, "lacks the row")

    def test_a_doubled_row_stops_it(self):
        self.check_stop({"extra": ("synthia_coli__drug_c", "lineage", "42", "ordered")}, "twice")

    def test_a_check_that_is_not_ok_stops_it(self):
        self.check_stop({"check_not_ok": ("synthia_coli__drug_b", "lineage", "47")}, "is not ok")

    def test_a_total_that_differs_stops_it(self):
        self.check_stop({}, "results: 30, want 768", "--results", "768")
        code, _, _, _ = run([("synthia_coli", edge_sets(), {})], "--results", "30", "--lineage", "5",
                            "--study", "0", "--species", "1")
        self.assertEqual(code, 0)


class Interval(unittest.TestCase):
    def test_it_reads_the_101st_and_the_3900th_sorted_values(self):
        # The resample medians come out as 3999, 3998, ..., 0; sorted, the 101st is 100.
        values = iter(range(compare.RESAMPLES - 1, -1, -1))
        with mock.patch.object(compare.statistics, "median", side_effect=lambda xs: next(values)):
            self.assertEqual(compare.species_interval({"a": [0.1], "b": [0.2]}, "seed"), (100, 3899))

    def test_it_gives_the_ends_of_record_py(self):
        record = load_record()
        by_species = {f"synthia_{k}": [Fraction(k * 7 % 11 - 5, 100), Fraction(k, 300)] for k in range(1, 8)}
        seed = compare.seed_string("lineage", "P1")
        ends = record.redraw(by_species, seed)
        self.assertEqual(compare.species_interval(by_species, seed), (ends[0], ends[1]))


class SignAndHolm(unittest.TestCase):
    def test_the_sign_test_is_one_sided_and_exact(self):
        self.assertEqual(compare.sign_p(8, 1), Fraction(10, 512))
        self.assertEqual(compare.sign_p(1, 8), Fraction(511, 512))
        self.assertEqual(compare.sign_p(0, 0), 1)

    def test_holm_steps_down_and_keeps_order(self):
        p = [Fraction(1, 100), Fraction(4, 100), Fraction(3, 100)]
        self.assertEqual(compare.holm(p), [Fraction(3, 100), Fraction(6, 100), Fraction(6, 100)])
        p = [Fraction(1, 1000), Fraction(20, 1000), Fraction(19, 1000)]
        self.assertEqual(compare.holm(p), [Fraction(3, 1000), Fraction(38, 1000), Fraction(38, 1000)])

    def test_holm_reads_the_lineage_split_alone(self):
        sets = dict(edge_sets())
        sets[("synthia_coli__drug_a", "study")] = ea_by_arm(60, 49)
        code, _, _, files = run([("synthia_coli", sets, {})])
        self.assertEqual(code, 0)
        self.assertEqual(pair(files, "lineage", "P1")["holm_p"], "%.6g" % min(1, 3 * float(compare.sign_p(2, 1))))
        study = pair(files, "study", "P1")
        self.assertEqual((study["holm_p"], study["decides"], study["median_sign"]), ("", "no, shown beside", "above 0"))


class ReleaseRule(unittest.TestCase):
    def test_it_counts_differences_under_002_and_exactly_002_is_not_under(self):
        # Differences +0.02, -0.02, +0.03, -0.03 and +7/300: under 0.02 are -0.02 and -0.03.
        code, out, _, files = run([("synthia_coli", edge_sets(), {})])
        self.assertEqual(code, 0)
        p1 = pair(files, "lineage", "P1")
        self.assertEqual((p1["under_0.02"], p1["compared"], p1["under_0.02_over_half"]), ("2", "5", "no"))
        self.assertIn("release rule (lineage P1 under 0.02): 2 of 5; more than half: no", out)
        self.assertEqual(pair(files, "lineage", "P2")["under_0.02"], "")

    def test_half_is_not_more_than_half(self):
        sets = dict(edge_sets())
        sets[("synthia_coli__drug_f", "lineage")] = ea_by_arm(50, 49)  # +0.01: under
        code, _, _, files = run([("synthia_coli", sets, {})])
        p1 = pair(files, "lineage", "P1")
        self.assertEqual((code, p1["under_0.02"], p1["compared"], p1["under_0.02_over_half"]), (0, "3", "6", "no"))
        sets[("synthia_coli__drug_g", "lineage")] = ea_by_arm(49, 49)  # 0: under
        code, _, _, files = run([("synthia_coli", sets, {})])
        p1 = pair(files, "lineage", "P1")
        self.assertEqual((code, p1["under_0.02"], p1["compared"], p1["under_0.02_over_half"]), (0, "4", "7", "yes"))


if __name__ == "__main__":
    unittest.main()
