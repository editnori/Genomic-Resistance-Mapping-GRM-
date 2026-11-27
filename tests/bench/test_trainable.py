"""Synthetic unit tests for trainable mode: fold parity via plans.py, strict
inside-outer-training-fold fitting (anti-leakage), grm-fold-file-v1 (study-whole
and time), per-fold no_call/failure accounting, and scratch-only OOF vectors."""

from __future__ import annotations

import copy
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

from scripts.bench.protocol import BenchError, SELECTION_ROLE_EVAL_ONLY, load_adapter
from scripts.bench.trainable import (
    FoldFitContext,
    FoldPredictions,
    fold_call_counts,
    run_trainable_dataset,
    verify_fold_parity,
    write_oof_scratch,
)
from tests.bench.common import (
    build_matching_reference_docs,
    build_synthetic_export,
    build_synthetic_fold_file,
)


def _synthetic_linear_fitter(ctx: FoldFitContext) -> FoldPredictions:
    """Deterministic training-only prototype fitter (requires numpy only).

    Computes class-balanced feature prevalence difference strictly on `ctx.fit_mask`,
    picks threshold on `ctx.hold_mask`, and predicts probabilities on `ctx.X_test`.
    """
    X_fit = ctx.X_train[ctx.fit_mask].astype(np.float64)
    y_fit = ctx.y_train[ctx.fit_mask]
    X_hold = ctx.X_train[ctx.hold_mask].astype(np.float64)
    y_hold = ctx.y_train[ctx.hold_mask]

    p_r = X_fit[y_fit].mean(axis=0)
    p_s = X_fit[~y_fit].mean(axis=0)
    w = p_r - p_s
    # Tune bias on hold_mask only
    z_hold = X_hold @ w
    bias = float(-0.5 * (z_hold[y_hold].mean() + z_hold[~y_hold].mean()))
    z_test = ctx.X_test.astype(np.float64) @ w + bias
    probs = 1.0 / (1.0 + np.exp(-np.clip(z_test, -20.0, 20.0)))
    scores = tuple(float(p) for p in probs)
    calls = tuple("R" if p >= 0.5 else "S" for p in scores)
    return FoldPredictions(calls=calls, scores=scores, meta={"bias": bias, "w_sum": float(np.sum(np.abs(w)))})


class TrainableFoldParityTests(unittest.TestCase):
    def test_plans_py_fold_parity_passes_and_catches_single_swapped_genome(self):
        with tempfile.TemporaryDirectory() as tmp:
            ds_dir, meta, _ = build_synthetic_export(Path(tmp), n_samples=100, seed=42)
            ref_pick, ref_report, release_report = build_matching_reference_docs(meta, "lineage")

            ok = verify_fold_parity(
                ds_dir,
                42,
                scheme="lineage",
                reference_result_doc=ref_pick,
                reference_report_doc=ref_report,
                release_report_doc=release_report,
            )
            self.assertTrue(all(ok["checks"].values()))

            # Swap two SAME-CLASS genomes between two outer folds: class counts stay identical,
            # but plans.py genome-membership SHA-256 MUST catch it!
            swapped_report = copy.deepcopy(ref_report)
            tf = swapped_report["splits"][0]["test_fold"]
            y = meta["resistant"]
            i0 = next(i for i in range(len(tf)) if tf[i] == 0 and y[i])
            i1 = next(i for i in range(len(tf)) if tf[i] == 1 and y[i])
            tf[i0], tf[i1] = tf[i1], tf[i0]
            with self.assertRaisesRegex(BenchError, "membership differs"):
                verify_fold_parity(ds_dir, 42, scheme="lineage", reference_result_doc=ref_pick, reference_report_doc=swapped_report)

            # Corrupt fold class count -> caught by result_shape
            bad_pick = copy.deepcopy(ref_pick)
            bad_pick["folds"][0]["test_resistant"] += 1
            with self.assertRaisesRegex(BenchError, "counts by class differ"):
                verify_fold_parity(ds_dir, 42, scheme="lineage", reference_result_doc=bad_pick)

            # Corrupt release input_fingerprint -> caught
            bad_rel = dict(release_report, input_fingerprint="0" * 32)
            with self.assertRaisesRegex(BenchError, "input_fingerprint differs"):
                verify_fold_parity(ds_dir, 42, scheme="lineage", release_report_doc=bad_rel)


class TrainableExecutionAndLeakageTests(unittest.TestCase):
    def test_trainable_runs_on_lineage_study_whole_and_time_splits(self):
        with tempfile.TemporaryDirectory() as tmp:
            tdir = Path(tmp)
            ds_dir, meta, _ = build_synthetic_export(tdir, n_samples=100, seed=42)
            ref_pick, ref_report, release_report = build_matching_reference_docs(meta, "lineage")
            ff_study = build_synthetic_fold_file(tdir / "study.json", meta, scheme="study-whole", n_test_folds=2, n_unused=4)
            ff_time = build_synthetic_fold_file(tdir / "time.json", meta, scheme="time", n_test_folds=1, n_unused=4)

            for tool in ("kover", "xgboost", "lightgbm", "elasticnet"):
                adapter = load_adapter(tool)
                res_lin = run_trainable_dataset(
                    adapter,
                    ds_dir,
                    "lineage",
                    42,
                    fitter=_synthetic_linear_fitter,
                    reference_result_doc=ref_pick,
                    reference_report_doc=ref_report,
                    release_report_doc=release_report,
                )
                self.assertEqual(len(res_lin["folds"]), 5)
                self.assertEqual(sum(f["n"] for f in res_lin["folds"]), 100)
                self.assertTrue(all("loss_sum" in f for f in res_lin["folds"]))

            # study-whole and time fold-file plans
            xgb = load_adapter("xgboost")
            res_study = run_trainable_dataset(xgb, ds_dir, "study-whole", 42, fold_file=ff_study, fitter=_synthetic_linear_fitter)
            self.assertEqual(len(res_study["folds"]), 2)
            self.assertEqual(sum(f["n"] for f in res_study["folds"]), 96)
            self.assertIsNotNone(res_study["fold_file_sha256"])

            res_time = run_trainable_dataset(xgb, ds_dir, "time", 42, fold_file=ff_time, fitter=_synthetic_linear_fitter)
            self.assertEqual(len(res_time["folds"]), 1)
            self.assertIsNotNone(res_time["fold_file_sha256"])

    def test_strict_anti_leakage_flipping_outer_test_labels_changes_no_training_fit(self):
        """Flipping every label in outer test fold 0 must leave fold 0's training
        screen, fitted weights, and predicted test scores bit-for-bit identical."""
        with tempfile.TemporaryDirectory() as tmp:
            tdir = Path(tmp)
            ds_dir_a, meta_a, _ = build_synthetic_export(tdir / "a", n_samples=100, seed=42)
            ff_a = build_synthetic_fold_file(tdir / "ff_a.json", meta_a, scheme="time", n_test_folds=1, n_unused=4)
            ff_doc = json.loads(ff_a.read_text(encoding="utf-8"))
            test_set = set(ff_doc["test"][0])

            # Copy export and flip EVERY test genome's label in dataset.json
            ds_dir_b = tdir / "b" / meta_a["id"]
            ds_dir_b.mkdir(parents=True, exist_ok=True)
            (ds_dir_b / "matrix.bin").write_bytes((ds_dir_a / "matrix.bin").read_bytes())
            meta_b = copy.deepcopy(meta_a)
            for i, g in enumerate(meta_b["samples"]):
                if g in test_set:
                    meta_b["resistant"][i] = not meta_b["resistant"][i]
            (ds_dir_b / "dataset.json").write_text(json.dumps(meta_b), encoding="utf-8")

            captured_a: list[tuple[np.ndarray, FoldPredictions]] = []
            captured_b: list[tuple[np.ndarray, FoldPredictions]] = []

            def spy_fitter(bucket):
                def _fit(ctx: FoldFitContext) -> FoldPredictions:
                    preds = _synthetic_linear_fitter(ctx)
                    bucket.append((ctx.X_train.copy(), preds))
                    return preds
                return _fit

            adapter = load_adapter("lightgbm")
            run_trainable_dataset(adapter, ds_dir_a, "time", 42, fold_file=ff_a, fitter=spy_fitter(captured_a))
            run_trainable_dataset(adapter, ds_dir_b, "time", 42, fold_file=ff_a, fitter=spy_fitter(captured_b))

            self.assertEqual(len(captured_a), 1)
            self.assertEqual(len(captured_b), 1)
            np.testing.assert_array_equal(captured_a[0][0], captured_b[0][0])
            self.assertEqual(captured_a[0][1].scores, captured_b[0][1].scores)
            self.assertEqual(captured_a[0][1].calls, captured_b[0][1].calls)
            self.assertEqual(captured_a[0][1].meta, captured_b[0][1].meta)

    def test_refuses_training_on_evaluation_only_cohort(self):
        with tempfile.TemporaryDirectory() as tmp:
            ds_dir, _, _ = build_synthetic_export(Path(tmp), n_samples=100, seed=42)
            adapter = load_adapter("xgboost")
            with self.assertRaisesRegex(BenchError, "evaluation_only"):
                run_trainable_dataset(
                    adapter,
                    ds_dir,
                    "lineage",
                    42,
                    selection_role=SELECTION_ROLE_EVAL_ONLY,
                    fitter=_synthetic_linear_fitter,
                )

    def test_no_call_i_and_failure_accounting_never_map_to_susceptible(self):
        truth = np.array([True, True, True, False, False, False], dtype=bool)
        calls = ("R", "no_call", "failure", "S", "I", "no_call")
        scores = (0.9, None, None, 0.1, None, None)
        c = fold_call_counts(truth, calls, scores)
        self.assertEqual(c["n"], 6)
        self.assertEqual((c["tp"], c["fp"], c["tn"], c["fn"]), (1, 0, 1, 0))
        self.assertEqual((c["no_call_r"], c["no_call_s"]), (1, 1))
        self.assertEqual((c["i_r"], c["i_s"]), (0, 1))
        self.assertEqual((c["failure_r"], c["failure_s"]), (1, 0))
        self.assertIn("loss_sum", c)

        # Fitter exception -> all test genomes counted as failure_r / failure_s, NOT fn/tn!
        with tempfile.TemporaryDirectory() as tmp:
            ds_dir, _, _ = build_synthetic_export(Path(tmp), n_samples=100, seed=42)
            adapter = load_adapter("elasticnet")

            def crashing_fitter(_ctx):
                raise ValueError("numerical crash")

            res = run_trainable_dataset(adapter, ds_dir, "lineage", 42, fitter=crashing_fitter)
            for f in res["folds"]:
                self.assertEqual(f["tp"] + f["fp"] + f["tn"] + f["fn"], 0)
                self.assertEqual(f["failure_r"] + f["failure_s"], f["n"])

    def test_oof_scratch_isolation_refuses_reference_and_result_overlap(self):
        with tempfile.TemporaryDirectory() as tmp:
            tdir = Path(tmp)
            ds_dir, meta, _ = build_synthetic_export(tdir / "exp", n_samples=100, seed=42)
            res_dir = tdir / "results"
            scratch_dir = tdir / "scratch_oof"
            adapter = load_adapter("xgboost")

            res = run_trainable_dataset(
                adapter,
                ds_dir,
                "lineage",
                42,
                fitter=_synthetic_linear_fitter,
                oof_scratch_dir=scratch_dir,
                result_dir=res_dir,
            )
            oof_file = Path(res["oof_scratch_path"])
            self.assertTrue(oof_file.is_file())
            oof_doc = json.loads(oof_file.read_text(encoding="utf-8"))
            self.assertTrue(oof_doc["local_only"])
            self.assertTrue(oof_doc["never_commit"])
            self.assertEqual(len(oof_doc["scores"]), 100)
            self.assertNotIn("samples", oof_doc)

            # Overlap with result_dir refused
            with self.assertRaisesRegex(BenchError, "strictly separate"):
                write_oof_scratch(
                    res_dir / "oof",
                    result_dir=res_dir,
                    tool="xgboost",
                    dataset=meta["id"],
                    split_scheme="lineage",
                    seed=42,
                    input_fingerprint=meta["input_fingerprint"],
                    fold_of=[0],
                    calls=["R"],
                    scores=[0.9],
                )

            # Inside repo reference/ refused
            from scripts.bench.vocab import REPO_ROOT
            with self.assertRaisesRegex(BenchError, "reference/"):
                write_oof_scratch(
                    REPO_ROOT / "reference" / "bad_oof",
                    result_dir=None,
                    tool="xgboost",
                    dataset=meta["id"],
                    split_scheme="lineage",
                    seed=42,
                    input_fingerprint=meta["input_fingerprint"],
                    fold_of=[0],
                    calls=["R"],
                    scores=[0.9],
                )


if __name__ == "__main__":
    unittest.main()
