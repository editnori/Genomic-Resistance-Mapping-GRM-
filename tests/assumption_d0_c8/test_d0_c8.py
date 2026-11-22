"""Synthetic tests of the D0 and C8 readers (scripts/assumption-d0-c8).

Run as CI runs the python suites:
    python3 -s -P -B -W error -m unittest discover -s tests/assumption_d0_c8 -t . -v
"""

import contextlib
import csv
import io
import json
import os
import sys
import tempfile
import unittest
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "scripts", "assumption-d0-c8"))
sys.path.insert(0, HERE)

import numpy as np  # noqa: E402

import c8  # noqa: E402
import d0  # noqa: E402
import foldread as fr  # noqa: E402
import synth as S  # noqa: E402

ALL = list(range(S.N))
# Resistance in the first three samples of every cluster: both classes in
# every fold, spread over every cluster.
SPREAD = [int(i % 6 < 3) for i in ALL]
# Resistance in clusters 0 and 1 only.
NARROW = [int(S.cluster(i) < 2) for i in ALL]
EMPTY5 = [set() for _ in range(S.FOLDS)]


def only(fold, rows):
    """A learner's five fold sets: `rows` on `fold`, nothing elsewhere."""
    return [set(rows) if f == fold else set() for f in range(S.FOLDS)]


class Case:
    """Build one synthetic dataset and read it with d0."""

    def __init__(self, tmp, name, **kw):
        self.paths = S.make(os.path.join(tmp, name), **kw)

    def read(self, reference=True, learners=None):
        report = fr.load_report(self.paths["report"])
        export = d0.Export(self.paths["export"])
        calls = d0.read_calls(self.paths["calls"], report["id"], "tool_default")
        ref = d0.read_reference(self.paths["reference"], report["id"]) if reference else None
        folds, data = d0.read_dataset(report, export, calls, ref, learners)
        return folds, {r["learner"]: r for r in data}


def fold_row(folds, learner, f):
    return next(r for r in folds if r["learner"] == learner and r["fold"] == f)


class D0Causes(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)

    def test_absent_when_no_row_meets_the_program(self):
        rows = S.noise(40, seed=1)
        c = Case(self.tmp.name, "absent", truth=SPREAD, program=SPREAD, rows=rows,
                 scores={"scm": S.scores_for(SPREAD, wrong_folds=(0,))})
        folds, data = c.read()
        self.assertTrue(data["scm"]["in_miss"])
        self.assertEqual(data["scm"]["cause"], "absent")
        self.assertEqual(data["scm"]["missed_folds"], 1)
        r = fold_row(folds, "scm", 0)
        self.assertEqual((r["cause"], r["matching_rows"]), ("absent", 0))
        self.assertIsNone(fold_row(folds, "scm", 1).get("cause"))
        for learner in ("cart", "rf", "gbt", "lasso"):
            self.assertFalse(data[learner]["in_miss"])
            self.assertIsNone(data[learner]["cause"])

    def test_removed_by_the_floor_and_flagged_one_class(self):
        program = [int(i in (0, 1)) for i in ALL]
        rows = S.noise(30, seed=2) + [[0] * S.N, program]
        c = Case(self.tmp.name, "floor", truth=SPREAD, program=program, rows=rows,
                 scores={"scm": S.scores_for(SPREAD, wrong_folds=range(5))})
        folds, data = c.read()
        self.assertEqual(data["scm"]["cause"], "removed")
        self.assertEqual(data["scm"]["folds_removed"], 5)
        for f in range(S.FOLDS):
            r = fold_row(folds, "scm", f)
            self.assertEqual((r["matching_rows"], r["inside_floor"], r["offered"]), (2, 0, 0))
            self.assertTrue(r["program_one_class"])

    def screen_case(self, chosen):
        flips = (12, 18, 24, 40)
        program = S.flip(SPREAD, flips)
        rows = [SPREAD, S.flip(SPREAD, (13,)), S.flip(SPREAD, (41,)), program] + S.noise(20, seed=3)
        miss = S.scores_for(SPREAD, wrong_folds=(0,))
        return Case(self.tmp.name, "screen", truth=SPREAD, program=program, rows=rows, k=3,
                    scores={l: miss for l in S.LEARNERS}, chosen=chosen)

    def test_past_the_screen_is_removed_for_rf_gbt_lasso_only(self):
        noise_row = 10
        c = self.screen_case({"scm": only(0, [noise_row]), "rf": only(0, [0]), "gbt": only(0, [1]),
                              "lasso": only(0, [2])})
        folds, data = c.read()
        for learner in ("rf", "gbt", "lasso"):
            self.assertEqual(data[learner]["cause"], "removed", learner)
            r = fold_row(folds, learner, 0)
            self.assertEqual((r["inside_floor"], r["offered"], r["past_screen"]), (1, 0, 1), learner)
        # The set covering machine and the tree read every candidate: the
        # same row is offered to them, and they did not choose it.
        for learner in ("scm", "cart"):
            self.assertEqual(data[learner]["cause"], "unused", learner)
            self.assertEqual(fold_row(folds, learner, 0)["offered"], 1, learner)

    def test_lineage_when_the_chosen_proxy_lies_in_two_clusters(self):
        rows = S.noise(20, seed=4) + [NARROW]
        m = len(rows) - 1
        c = Case(self.tmp.name, "lineage", truth=NARROW, program=NARROW, rows=rows,
                 scores={"scm": S.scores_for(NARROW, wrong_folds=(0,)),
                         "rf": S.scores_for(NARROW, wrong_folds=(0,))},
                 chosen={"scm": only(0, [m]), "rf": only(0, [m, 3])})
        folds, data = c.read()
        self.assertEqual(data["scm"]["cause"], "lineage")
        self.assertEqual(data["rf"]["cause"], "lineage")
        r = fold_row(folds, "rf", 0)
        self.assertEqual((r["chosen_agreeing"], r["chosen_agreeing_in_2"]), (1, 1))

    def test_lineage_counts_two_training_clusters_and_no_test_cluster(self):
        # Resistance in clusters 0, 1 and 2. Fold 0 tests clusters 0 and 5,
        # so the chosen proxy lies in exactly two training clusters (1 and 2)
        # and is Lineage at the boundary. Counting the tested cluster 0 too,
        # or asking for fewer than two clusters, makes it Chosen.
        three = [int(S.cluster(i) < 3) for i in ALL]
        rows = S.noise(20, seed=13) + [three]
        m = len(rows) - 1
        c = Case(self.tmp.name, "two", truth=three, program=three, rows=rows,
                 scores={"scm": S.scores_for(three, wrong_folds=(0,))}, chosen={"scm": only(0, [m])})
        folds, data = c.read()
        r = fold_row(folds, "scm", 0)
        self.assertEqual((r["chosen_agreeing"], r["chosen_agreeing_in_2"]), (1, 1))
        self.assertEqual(data["scm"]["cause"], "lineage")

    def test_spread_choice_is_chosen_and_absence_polarity_counts(self):
        rows = S.noise(20, seed=5) + [SPREAD, [1 - v for v in SPREAD]]
        m, neg = len(rows) - 2, len(rows) - 1
        miss = S.scores_for(SPREAD, wrong_folds=(0,))
        c = Case(self.tmp.name, "chosen", truth=SPREAD, program=SPREAD, rows=rows,
                 scores={"scm": miss, "rf": miss, "cart": miss},
                 chosen={"scm": only(0, [m]), "rf": only(0, [neg]), "cart": only(0, [2])})
        folds, data = c.read()
        self.assertEqual(data["scm"]["cause"], "chosen")
        self.assertEqual(data["rf"]["cause"], "chosen")
        self.assertEqual(data["cart"]["cause"], "unused")
        # Both polarities match, so each fold has two matching rows.
        self.assertEqual(fold_row(folds, "cart", 0)["matching_rows"], 2)

    def test_tie_counts_in_both_causes(self):
        row = S.flip(SPREAD, [i for i in ALL if S.TEST_FOLD[i] == 1])
        rows = S.noise(20, seed=6) + [row]
        c = Case(self.tmp.name, "tie", truth=SPREAD, program=SPREAD, rows=rows,
                 scores={"scm": S.scores_for(SPREAD, wrong_folds=(0, 1))})
        folds, data = c.read()
        self.assertEqual(fold_row(folds, "scm", 0)["cause"], "absent")
        self.assertEqual(fold_row(folds, "scm", 1)["cause"], "unused")
        self.assertEqual(data["scm"]["cause"], "absent+unused")
        table, every, _ = d0.combine(list(data.values()))
        scm = next(t for t in table if t[1] == "scm")
        col = {h: i for i, h in enumerate(d0.SUMMARY_HEADER)}
        self.assertEqual((scm[col["datasets_absent"]], scm[col["datasets_unused"]]), (1, 1))
        self.assertAlmostEqual(scm[col["gap_share_absent"]], 0.5)
        self.assertAlmostEqual(scm[col["gap_share_unused"]], 0.5)
        self.assertEqual(scm[col["ties"]], 1)

    def test_no_missed_fold_falls_back_to_the_largest_gap(self):
        truth = [0 if S.cluster(i) in (2, 7) else SPREAD[i] for i in ALL]
        rows = S.noise(20, seed=7) + [truth]
        # Wrong on fold 2, whose test genomes are all susceptible, and on
        # one resistant genome of fold 0: fold 0's gap is 1/12, above 0 and
        # under 0.10, the largest of the folds that test both classes.
        c = Case(self.tmp.name, "fallback", truth=truth, program=truth, rows=rows,
                 scores={"scm": S.scores_for(truth, wrong_folds=(2,), wrong_samples=(0,))})
        folds, data = c.read()
        self.assertTrue(data["scm"]["in_miss"])
        self.assertTrue(data["scm"]["fallback"])
        self.assertEqual(data["scm"]["missed_folds"], 1)
        self.assertFalse(fold_row(folds, "scm", 2)["both_classes"])
        self.assertFalse(fold_row(folds, "scm", 2)["missed"])
        self.assertAlmostEqual(fold_row(folds, "scm", 0)["fold_gap"], 1 / 12)
        self.assertTrue(fold_row(folds, "scm", 0)["missed"])
        self.assertEqual(data["scm"]["cause"], "unused")

    def test_no_fallback_to_a_fold_where_the_learner_did_as_well(self):
        truth = [0 if S.cluster(i) in (2, 7) else SPREAD[i] for i in ALL]
        rows = S.noise(20, seed=7) + [truth]
        # Wrong only on fold 2, whose test genomes are all susceptible: the
        # pooled gap is 0.10 or more, and every fold that tests both classes
        # has a gap of 0. The pair has no missed fold and no cause.
        c = Case(self.tmp.name, "nofallback", truth=truth, program=truth, rows=rows,
                 scores={"scm": S.scores_for(truth, wrong_folds=(2,))})
        folds, data = c.read()
        self.assertTrue(data["scm"]["in_miss"])
        self.assertFalse(data["scm"]["fallback"])
        self.assertEqual(data["scm"]["missed_folds"], 0)
        self.assertIsNone(data["scm"]["cause"])
        self.assertFalse(any(r["missed"] for r in folds))
        self.assertTrue(all(r["fold_gap"] <= 0 for r in folds if r["learner"] == "scm" and r["both_classes"]))
        table, every, pairs = d0.combine(list(data.values()))
        col = {h: i for i, h in enumerate(d0.SUMMARY_HEADER)}
        scm = next(t for t in table if t[1] == "scm")
        self.assertEqual((scm[col["miss_datasets"]], scm[col["no_missed_fold"]]), (1, 1))
        self.assertEqual(sum(every[k] for k in d0.CAUSES), 0)
        self.assertEqual(pairs, 1)

    def test_a_fold_without_a_miss_gets_no_cause(self):
        rows = S.noise(20, seed=8) + [SPREAD]
        c = Case(self.tmp.name, "nomiss", truth=SPREAD, program=SPREAD, rows=rows,
                 scores={"scm": S.scores_for(SPREAD, wrong_samples=(0,))})
        folds, data = c.read()
        self.assertFalse(data["scm"]["in_miss"])
        self.assertEqual(data["scm"]["missed_folds"], 0)
        self.assertTrue(all(r.get("cause") is None for r in folds))

    def test_a_null_lasso_chooses_nothing(self):
        c = self.screen_case({"lasso": [None] * S.FOLDS})
        report = fr.load_report(c.paths["report"])
        self.assertEqual(fr.fold_sets(report, "lineage", "lasso"), [frozenset()] * S.FOLDS)

    def test_no_call_genomes_are_left_out(self):
        rows = S.noise(20, seed=9) + [SPREAD]
        c = Case(self.tmp.name, "nocall", truth=SPREAD, program=SPREAD, rows=rows, no_call=(40, 41),
                 scores={"scm": S.scores_for(SPREAD, wrong_folds=(0,))})
        folds, data = c.read()
        self.assertEqual(data["scm"]["cause"], "unused")
        self.assertEqual(fold_row(folds, "scm", 0)["matching_rows"], 1)


class D0Screen(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)

    def test_candidates_scores_and_ties_follow_learn_screen(self):
        # Row 0 equals row 5 on fold 0's training samples: only row 0 is a
        # candidate there. Rows 1 and 2 tie; the lower row comes first.
        base = S.noise(12, seed=11)
        rows = [S.flip(SPREAD, [i for i in ALL if S.TEST_FOLD[i] == 0])] + base[:4] + [SPREAD] + base[4:]
        rows[2] = list(rows[1])
        rows[2][0] = 1 - rows[2][0]  # differs from row 1 only on a fold-0 test sample
        p = S.make(os.path.join(self.tmp.name, "s"), truth=SPREAD, program=SPREAD, rows=rows, k=4)
        export = d0.Export(p["export"])
        n = S.N
        f = 0
        t = S.training(f)
        train = d0.to_words([i in set(t) for i in ALL], n)
        called = d0.to_words([True] * n, n)
        program = d0.to_words(SPREAD, n)
        y = program
        fits = {"rf": train}
        for name, held in S.holdouts(f).items():
            fits[name] = train & ~d0.to_words([i in set(held) for i in ALL], n)
        scan = d0.scan_fold(export.x, train, called, program, fits, y, 3)
        self.assertEqual(list(scan["candidates"]), S.candidates(rows, t, 3))
        self.assertNotIn(5, list(scan["candidates"]))
        self.assertNotIn(2, list(scan["candidates"]))
        for name, held in S.holdouts(f).items():
            fit = [i for i in t if i not in held]
            want = [S.chi(rows[r], SPREAD, fit) for r in scan["candidates"]]
            self.assertEqual(list(scan["scores"][name]), want, name)
        self.assertEqual([int(r) for r in d0.top(scan["scores"]["rf"], scan["candidates"], 4)],
                         S.screen(rows, SPREAD, t, t, 3, 4))
        # The matching row 5 maps to its pattern's candidate, row 0.
        self.assertEqual(sorted(int(r) for r in scan["match_rep"]), [0, 0])

    def test_read_dataset_screens_boosting_and_lasso_without_their_holdouts(self):
        # The export holds the forest's screen only, so nothing gates the
        # other two: the fits handed to scan_fold are checked here.
        rows = S.noise(20, seed=14) + [SPREAD]
        miss = S.scores_for(SPREAD, wrong_folds=(0,))
        c = Case(self.tmp.name, "holdout", truth=SPREAD, program=SPREAD, rows=rows,
                 scores={"gbt": miss, "lasso": miss}, chosen={"gbt": only(0, [3]), "lasso": only(0, [3])})
        seen = []
        real = d0.scan_fold

        def spy(x, train, called, program, fits, y, floor):
            seen.append((train.copy(), {k: v.copy() for k, v in fits.items()}))
            return real(x, train, called, program, fits, y, floor)

        with mock.patch.object(d0, "scan_fold", spy):
            c.read()
        self.assertEqual(len(seen), 1)
        train, fits = seen[0]
        self.assertTrue(np.array_equal(fits["rf"], train))
        for name, held in S.holdouts(0).items():
            want = d0.to_words([S.TEST_FOLD[i] != 0 and i not in held for i in ALL], S.N)
            self.assertTrue(np.array_equal(fits[name], want), name)

    def test_ties_at_the_cut_keep_the_lower_rows(self):
        scores = np.array([1.0, 2.0, 2.0, 2.0, 0.5])
        rows = np.array([0, 1, 2, 3, 4])
        self.assertEqual(list(d0.top(scores, rows, 2)), [1, 2])

    def test_popcount_and_words(self):
        w = d0.to_words([1, 0, 1] + [0] * 60 + [1], 64)
        self.assertEqual(int(d0.popcount(w[None, :])[0]), 3)
        self.assertEqual(int(w[0]), 1 | 4 | (1 << 63))


class D0Gates(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        rows = S.noise(20, seed=12) + [SPREAD]
        self.c = Case(self.tmp.name, "g", truth=SPREAD, program=SPREAD, rows=rows,
                      scores={"scm": S.scores_for(SPREAD, wrong_folds=(0,))}, chosen={"scm": only(0, [3])})

    def edit_json(self, path, change):
        with open(path, encoding="utf-8") as f:
            doc = json.load(f)
        change(doc)
        with open(path, "w", encoding="utf-8") as f:
            json.dump(doc, f)

    def test_reads_clean(self):
        self.c.read()

    def test_forest_screen_must_equal_the_exports(self):
        meta = os.path.join(self.c.paths["export"], "dataset.json")
        self.edit_json(meta, lambda d: d["splits"][0]["screens"][0]["screened_rows"].reverse())
        with self.assertRaisesRegex(d0.Gate, "forest screen"):
            self.c.read()

    def test_matrix_sha_must_match(self):
        path = os.path.join(self.c.paths["export"], "matrix.bin")
        with open(path, "r+b") as f:
            b = f.read(1)
            f.seek(0)
            f.write(bytes([b[0] ^ 1]))
        with self.assertRaisesRegex(d0.Gate, "sha256"):
            d0.Export(self.c.paths["export"])

    def test_fold_model_pattern_must_be_its_row(self):
        self.edit_json(self.c.paths["report"],
                       lambda d: d["fold_models"]["features"][0].update(bits=S.bits_text([1] * S.N)))
        with self.assertRaisesRegex(d0.Gate, "pattern"):
            self.c.read()

    def test_plan_must_be_the_reports_folds(self):
        meta = os.path.join(self.c.paths["export"], "dataset.json")

        def swap(d):
            o = d["splits"][0]["plan"]["outer"]
            o[0], o[6] = o[6], o[0]
        self.edit_json(meta, swap)
        with self.assertRaisesRegex(d0.Gate, "outer folds"):
            self.c.read()

    def test_reference_counts_must_match(self):
        with open(self.c.paths["reference"], encoding="utf-8") as f:
            text = f.read()
        lines = text.splitlines()
        cells = lines[1].split("\t")
        cells[2] = str(int(cells[2]) + 1)
        lines[1] = "\t".join(cells)
        with open(self.c.paths["reference"], "w", encoding="utf-8") as f:
            f.write("\n".join(lines) + "\n")
        with self.assertRaisesRegex(d0.Gate, "reference"):
            self.c.read()

    def test_an_export_without_fold_models_reads_every_screen(self):
        flips = (12, 18, 24, 40)
        program = S.flip(SPREAD, flips)
        rows = [SPREAD, S.flip(SPREAD, (13,)), S.flip(SPREAD, (41,)), program] + S.noise(20, seed=37)
        miss = S.scores_for(SPREAD, wrong_folds=(0,))
        p = S.make(os.path.join(self.tmp.name, "nofm"), truth=SPREAD, program=program, rows=rows, k=3,
                   scores={l: miss for l in S.LEARNERS}, release=True)
        report = fr.load_report(p["report"], need_models=False)
        calls = d0.read_calls(p["calls"], report["id"], "tool_default")
        _, data = d0.read_dataset(report, d0.Export(p["export"]), calls)
        got = {l: data_l["cause"] for l, data_l in ((r["learner"], r) for r in data)}
        self.assertEqual(got, {"scm": "offered", "cart": "offered", "rf": "removed", "gbt": "removed",
                               "lasso": "removed"})

    def test_no_fold_models_is_refused(self):
        self.edit_json(self.c.paths["report"], lambda d: d.pop("fold_models"))
        with self.assertRaisesRegex(fr.ReportError, "no fold_models"):
            fr.load_report(self.c.paths["report"])


class D0Release(unittest.TestCase):
    """A released report has no fold models: the dataset matrix comes from
    the species matrix and the causes after Removed stay open."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)

    def build(self, name, truth, program, rows, scores, cards=()):
        p = S.make(os.path.join(self.tmp.name, name), truth=truth, program=program, rows=rows, k=3,
                   scores=scores, cards=cards, release=True)
        p["species"] = S.make_species(os.path.join(self.tmp.name, name, "species"), rows)
        return p

    def read(self, p):
        report = fr.load_report(p["report"], need_models=False)
        calls = d0.read_calls(p["calls"], report["id"], "tool_default")
        ref = d0.read_reference(p["reference"], report["id"])
        folds, data = d0.read_release(report, p["species"], calls, ref)
        return folds, {r["learner"]: r for r in data}

    def test_projection_gives_the_dataset_rows_back(self):
        rows = S.noise(30, seed=31)
        p = self.build("proj", SPREAD, SPREAD, rows, {})
        x, kept = d0.project(p["species"], [f"g{i:03d}" for i in ALL], 3)
        want = np.array([[sum(1 << i for i, v in enumerate(r) if v)] for r in rows], dtype=np.uint64)
        self.assertEqual(x.tolist(), want.tolist())
        self.assertEqual(kept, 2 * len(rows))

    def test_projection_keeps_rows_on_the_floor(self):
        # Present in exactly 3 report genomes, and absent from exactly 3:
        # both bounds of [floor, n - floor] are inclusive.
        on_floor = [int(i in (5, 17, 44)) for i in ALL]
        under_top = [int(i not in (2, 30, 59)) for i in ALL]
        rows = [on_floor, under_top] + S.noise(10, seed=33)
        p = self.build("floor", SPREAD, SPREAD, rows, {})
        x, kept = d0.project(p["species"], [f"g{i:03d}" for i in ALL], 3)
        self.assertEqual(x.shape[0], len(rows))
        self.assertEqual(int(x[0, 0]), sum(1 << i for i in (5, 17, 44)))
        self.assertEqual(kept, 2 * len(rows))

    def test_open_causes_without_fold_models(self):
        flips = (12, 18, 24, 40)
        program = S.flip(SPREAD, flips)
        rows = [SPREAD, S.flip(SPREAD, (13,)), S.flip(SPREAD, (41,)), program] + S.noise(20, seed=32)
        miss = S.scores_for(SPREAD, wrong_folds=(0,))
        p = self.build("open", SPREAD, program, rows, {l: miss for l in S.LEARNERS}, cards=[(3, ["blaFAM"])])
        folds, data = self.read(p)
        self.assertEqual(data["rf"]["cause"], "removed")
        self.assertEqual(data["scm"]["cause"], "offered")
        self.assertEqual(data["cart"]["cause"], "offered")
        self.assertEqual(data["gbt"]["cause"], "inside_floor")
        self.assertEqual(data["lasso"]["cause"], "inside_floor")
        self.assertEqual(data["rf"]["source"], "release report, no fold models")
        _, every, pairs = d0.combine(list(data.values()))
        self.assertIn("not settled", d0.ordering(every, pairs))

    def test_absent_and_removed_are_read_without_fold_models(self):
        rows = S.noise(30, seed=33)
        p = self.build("abs", SPREAD, SPREAD, rows, {"scm": S.scores_for(SPREAD, wrong_folds=(0,))})
        _, data = self.read(p)
        self.assertEqual(data["scm"]["cause"], "absent")
        # A dataset row must be inside the floor on all 60 genomes, so a row
        # is removed only where a fold's training genomes leave it under.
        four = [int(i in (0, 1, 2, 6)) for i in ALL]
        p2 = self.build("rem", SPREAD, four, S.noise(30, seed=34) + [four],
                        {"lasso": S.scores_for(SPREAD, wrong_folds=range(5))})
        folds2, data2 = self.read(p2)
        self.assertEqual(fold_row(folds2, "lasso", 0)["cause"], "removed")
        self.assertTrue(fold_row(folds2, "lasso", 0)["program_one_class"])
        self.assertEqual(fold_row(folds2, "lasso", 1)["cause"], "inside_floor")
        self.assertEqual((data2["lasso"]["folds_removed"], data2["lasso"]["cause"]), (1, "inside_floor"))

    def test_a_card_must_be_its_projected_row(self):
        rows = S.noise(30, seed=35)
        p = self.build("card", SPREAD, SPREAD, rows, {}, cards=[(4, ["blaFAM"])])
        with open(p["report"], encoding="utf-8") as f:
            doc = json.load(f)
        doc["features"][0]["bits"] = S.bits_text(rows[5])
        with open(p["report"], "w", encoding="utf-8") as f:
            json.dump(doc, f)
        with self.assertRaisesRegex(d0.Gate, "feature card"):
            self.read(p)

    def test_the_row_count_must_be_the_reports(self):
        rows = S.noise(30, seed=36)
        p = self.build("count", SPREAD, SPREAD, rows, {})
        with open(p["report"], encoding="utf-8") as f:
            doc = json.load(f)
        doc["counts"]["features"] += 1
        with open(p["report"], "w", encoding="utf-8") as f:
            json.dump(doc, f)
        with self.assertRaisesRegex(d0.Gate, "projected rows"):
            self.read(p)


class D0Combine(unittest.TestCase):
    def row(self, learner, cause, gap):
        return {"learner": learner, "in_miss": True, "gap": gap, "cause": cause}

    def test_shares_and_ordering(self):
        rows = [self.row("scm", "absent", 0.3), self.row("scm", "unused", 0.1),
                self.row("rf", "absent", 0.2), self.row("rf", "removed+unused", 0.2),
                {"learner": "cart", "in_miss": False, "gap": 0.0, "cause": None}]
        table, every, pairs = d0.combine(rows)
        self.assertEqual(pairs, 4)
        col = {h: i for i, h in enumerate(d0.SUMMARY_HEADER)}
        scm = next(t for t in table if t[1] == "scm")
        self.assertAlmostEqual(scm[col["gap_share_absent"]], 0.75)
        self.assertEqual(scm[col["miss_datasets"]], 2)
        rf = next(t for t in table if t[1] == "rf")
        self.assertAlmostEqual(rf[col["gap_share_removed"]], 0.25)
        self.assertEqual(every, {"absent": 2, "removed": 1, "lineage": 0, "chosen": 0, "unused": 2,
                                 "offered": 0, "inside_floor": 0})
        self.assertIn("no single cause", d0.ordering(every, 4))
        # 40% of the miss pairs, not of the placements, and a tie never orders.
        self.assertIn("C3 first", d0.ordering({"absent": 2, "removed": 1, "lineage": 0, "chosen": 1, "unused": 1}, 5))
        self.assertIn("no single cause", d0.ordering({"absent": 1, "removed": 0, "lineage": 0, "chosen": 3, "unused": 0}, 4))

    def test_combine_refuses_a_twice_given_or_missing_directory(self):
        rows = [dict(self.row("scm", "absent", 0.3), dataset="a__x"),
                dict(self.row("rf", "unused", 0.2), dataset="a__x"),
                dict(self.row("scm", "unused", 0.1), dataset="b__y")]
        d0.complete(rows, 2, 3)
        with self.assertRaises(d0.Gate):
            d0.complete(rows + rows[:1], 2, 3)
        with self.assertRaises(d0.Gate):
            d0.complete(rows[:2], 2, 3)
        with self.assertRaises(d0.Gate):
            d0.complete(rows, 2, 4)

    def test_combine_reads_written_tables(self):
        with tempfile.TemporaryDirectory() as tmp:
            rows = [{"dataset": "a__b", "species": "a", "learner": "scm", "program_ba": 0.9,
                     "learner_lineage_ba": 0.7, "gap": 0.2, "in_miss": True, "missed_folds": 1,
                     "fallback": False, "folds_absent": 1, "folds_removed": 0, "folds_lineage": 0,
                     "folds_chosen": 0, "folds_unused": 0, "cause": "absent"}]
            d0.write_dataset(tmp, [], rows)
            back = d0.read_rows(os.path.join(tmp, "d0-datasets.tsv"))
            self.assertEqual((back[0]["cause"], back[0]["in_miss"], back[0]["gap"]), ("absent", True, 0.2))
            out = os.path.join(tmp, "sum")
            with contextlib.redirect_stdout(io.StringIO()):
                d0.main(["combine", "--out", out, tmp])
            with open(os.path.join(out, "d0-summary.tsv"), encoding="utf-8") as f:
                self.assertIn("tests nothing", f.read())
            # With intcount's marks the table adds the M and not-M rows.
            marks = os.path.join(tmp, "marks.tsv")
            with open(marks, "w", encoding="utf-8") as f:
                f.write("dataset\tpairs_S+I\tpairs_I+R\na__b\t0\t2\n")
            with contextlib.redirect_stdout(io.StringIO()):
                d0.main(["combine", "--marks", marks, "--out", out, tmp])
            with open(os.path.join(out, "d0-summary.tsv"), encoding="utf-8", newline="") as f:
                got = {(r["subset"], r["learner"]): r["miss_datasets"] for r in csv.DictReader(f, delimiter="\t")}
            self.assertEqual((got[("all", "scm")], got[("M", "scm")], got[("not_M", "scm")]), ("1", "1", "0"))
            with open(marks, "w", encoding="utf-8") as f:
                f.write("dataset\tpairs_S+I\tpairs_I+R\nc__d\t0\t0\n")
            with self.assertRaises(d0.Gate):
                d0.main(["combine", "--marks", marks, "--out", out, tmp])


class C8(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)

    def test_jaccard_pairs(self):
        self.assertEqual(fr.pairwise([{1, 2}, {2, 3}, set()]), [1 / 3, 0.0, 0.0])
        self.assertEqual(fr.pairwise([set(), set(), {1}]), [0.0, 0.0])
        self.assertIsNone(fr.median(fr.pairwise([set(), set()])))

    def make(self, chosen, cards=()):
        rows = S.noise(10, seed=13)
        rows[1] = S.flip(rows[0], [i for i in ALL if S.TEST_FOLD[i] == 0])
        p = S.make(os.path.join(self.tmp.name, "c8"), truth=SPREAD, program=SPREAD, rows=rows,
                   chosen=chosen, cards=cards)
        return fr.load_report(p["report"])

    def test_dataset_rows(self):
        chosen = {
            "scm": [{1}, {0}, {0}, {0}, {0}],  # rows 0 and 1 differ only on fold 0's test samples
            "rf": [{2, 3}, {2, 3}, {2, 3}, {2, 3}, {2, 3}],
            "lasso": [None, {4}, {4, 5}, {4}, {4}],
        }
        report = self.make(chosen, cards=[(2, ["blaFAM"]), (3, []), (4, ["aacFAM", "aphFAM"])])
        rows = c8.dataset_rows(report, True)
        col = {h: i for i, h in enumerate(c8.DATASET_HEADER)}
        by = {r[col["learner"]]: r for r in rows}
        scm = by["scm"]
        # By row number: fold 0 shares nothing with the others (4 pairs at 0),
        # folds 1 to 4 agree (6 pairs at 1).
        self.assertEqual(scm[col["lineage_pairs"]], 10)
        self.assertEqual(scm[col["lineage_median_jaccard"]], 1.0)
        self.assertEqual(scm[col["lineage_median_jaccard_pattern"]], 1.0)
        self.assertEqual(scm[col["lineage_distinct"]], 2)
        rf = by["rf"]
        self.assertEqual(rf[col["lineage_median_jaccard"]], 1.0)
        self.assertIsNone(rf[col["lineage_median_jaccard_pattern"]])
        self.assertEqual((rf[col["lineage_distinct"]], rf[col["lineage_with_card"]],
                          rf[col["lineage_with_family"]]), (2, 2, 1))
        lasso = by["lasso"]
        # A null fold against four with {4}: 4 pairs at 0; {4} vs {4,5}: 3 at 0.5; 3 at 1.
        self.assertEqual(lasso[col["lineage_pairs"]], 10)
        self.assertEqual(lasso[col["lineage_median_jaccard"]], 0.5)
        self.assertEqual((lasso[col["lineage_with_card"]], lasso[col["lineage_with_family"]]), (1, 1))
        cart = by["cart"]
        self.assertEqual((cart[col["lineage_pairs"]], cart[col["lineage_median_jaccard"]]), (0, None))
        # No random split ran.
        self.assertEqual((scm[col["random_folds"]], scm[col["random_median_jaccard"]]), (0, None))

    def test_pattern_reading_differs_from_numbers(self):
        report = self.make({"scm": [{1}, {0}, set(), set(), set()]})
        sets = fr.fold_sets(report, "lineage", "scm")
        # Seven pairs hold a chosen feature; the three of empty folds are left out.
        self.assertEqual(fr.pairwise(sets), [0.0] * 7)
        pats = fr.patterns(report)
        by_pattern = fr.pairwise_by_pattern(sets, pats, S.TEST_FOLD, S.FOLDS)
        self.assertEqual(by_pattern[0], 1.0)

    def test_pattern_reading_leaves_out_both_folds_test_samples(self):
        def as_int(p):
            return sum(1 << i for i, v in enumerate(p) if v)

        sets = [{0}, {1}, set(), set(), set()]
        for tested_by, same in ((0, True), (1, True), (2, False)):
            i = next(i for i in ALL if S.TEST_FOLD[i] == tested_by)
            pats = {0: as_int(SPREAD), 1: as_int(S.flip(SPREAD, (i,)))}
            got = fr.pairwise_by_pattern(sets, pats, S.TEST_FOLD, S.FOLDS)
            # Pair (0, 1) comes first; the two features differ only on a
            # sample that fold `tested_by` tests.
            self.assertEqual(got[0], 1.0 if same else 0.0, tested_by)

    def test_summary_counts(self):
        rows = []
        col = {h: i for i, h in enumerate(c8.DATASET_HEADER)}
        for j, (v, mark) in enumerate([(0.1, True), (0.5, False), (0.9, False), (0.85, None)]):
            r = [None] * len(c8.DATASET_HEADER)
            r[col["dataset"]], r[col["species"]], r[col["m_mark"]], r[col["learner"]] = f"d{j}__x", "d", mark, "rf"
            for s in fr.SCHEMES:
                r[col[f"{s}_folds"]] = 5 if s == "lineage" else 0
                r[col[f"{s}_pairs"]] = 10
                r[col[f"{s}_median_jaccard"]] = v if s == "lineage" else None
                r[col[f"{s}_distinct"]] = 10
                r[col[f"{s}_with_card"]] = 4
                r[col[f"{s}_with_family"]] = 1
            rows.append(r)
        out = c8.summarize(rows)
        h = {k: i for i, k in enumerate(c8.SUMMARY_HEADER)}
        allrf = next(o for o in out if o[h["subset"]] == "all" and o[h["learner"]] == "rf"
                     and o[h["scheme"]] == "lineage")
        self.assertEqual((allrf[h["datasets"]], allrf[h["under_0.2"]], allrf[h["over_0.8"]]), (4, 1, 2))
        self.assertAlmostEqual(allrf[h["median"]], 0.675)
        self.assertAlmostEqual(allrf[h["q1"]], 0.4)
        self.assertAlmostEqual(allrf[h["family_share_of_carded"]], 0.25)
        m = next(o for o in out if o[h["subset"]] == "M" and o[h["learner"]] == "rf" and o[h["scheme"]] == "lineage")
        self.assertEqual(m[h["datasets"]], 1)
        self.assertTrue(any(o[h["subset"]] == "unmarked" for o in out))

    def test_final_models(self):
        report = self.make({}, cards=[(2, ["blaFAM"]), (3, []), (4, ["aacFAM"])])
        report["models"] = {"scm": {"kind": "conjunction", "p": 1.0, "rules": [{"feature": 2, "presence": True}]},
                            "cart": {"nodes": [{"split": {"feature": 3}}, {"value": 0.0}]},
                            "rf_importance": [{"feature": 2, "weight": 0.5}, {"feature": 9, "weight": 0.1}],
                            "gbt_importance": [], "lasso": None}
        rows = c8.final_rows(report, False)
        col = {h: i for i, h in enumerate(c8.FINAL_HEADER)}
        got = {r[col["learner"]]: (r[col["features"]], r[col["with_card"]], r[col["with_family"]]) for r in rows}
        self.assertEqual(got, {"scm": (1, 1, 1), "cart": (1, 1, 0), "rf": (2, 1, 1), "gbt": (0, 0, 0),
                               "lasso": (0, 0, 0)})
        out = c8.final_summary(rows)
        h = {k: i for i, k in enumerate(c8.FINAL_SUMMARY_HEADER)}
        # Feature 2 is chosen by scm and rf and counts twice in the summed row.
        anyrow = next(o for o in out if o[h["subset"]] == "all" and o[h["learner"]] == "any_summed")
        self.assertEqual((anyrow[h["datasets"]], anyrow[h["with_a_family_feature"]], anyrow[h["features"]],
                          anyrow[h["with_card"]], anyrow[h["with_family"]]), (1, 1, 4, 3, 2))

    def test_marks(self):
        path = os.path.join(self.tmp.name, "marks.tsv")
        with open(path, "w", encoding="utf-8") as f:
            f.write("dataset\tpairs_S\tpairs_S+I\tpairs_I+R\n")
            f.write("a__x\t3\t0\t0\na__y\t0\t2\t0\na__z\t0\t0\t1\n")
        self.assertEqual(fr.read_marks(path), {"a__x": False, "a__y": True, "a__z": True})

    def test_main_writes_only_counts(self):
        self.make({"scm": [{1}, {0}, {0}, {0}, {0}]}, cards=[(0, ["blaFAM"])])
        out = os.path.join(self.tmp.name, "out")
        rdir = os.path.join(self.tmp.name, "reports")
        os.makedirs(rdir)
        os.replace(os.path.join(self.tmp.name, "c8", "report.json"),
                   os.path.join(rdir, "synthia_example__drugamycin.json"))
        with contextlib.redirect_stdout(io.StringIO()):
            c8.main(["--reports", rdir, "--out", out])
        with open(os.path.join(out, "c8-datasets.tsv"), encoding="utf-8") as f:
            text = f.read()
        self.assertNotIn("blaFAM", text)
        self.assertNotIn("g00", text)
        self.assertEqual(len(text.splitlines()), 1 + len(fr.LEARNERS))


if __name__ == "__main__":
    unittest.main()
