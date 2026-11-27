"""End-to-end unit tests for fixed-learner comparator (`compare.py`)."""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

from tests.point_mutations.common import C7_SCAN, CMP, FT, SCRIPTS_DIR, build_synthetic_dataset


class L2LogisticSolverTests(unittest.TestCase):
    def _assert_kkt(self, X: np.ndarray, y: np.ndarray, lam: float, res: dict):
        self.assertTrue(res["converged"], f"solver did not converge: grad_norm={res['grad_norm']}")
        w_norm = CMP.class_balanced_weights(y) / float(X.shape[0])
        p = CMP.predict_l2_logistic(X, res["coef"], res["intercept"])
        resid = w_norm * (p - y.astype(np.float64))
        g_w = X.T @ resid + lam * res["coef"]
        g_b = float(np.sum(resid))
        self.assertLess(float(np.max(np.abs(g_w))), 1e-7)
        self.assertLess(abs(g_b), 1e-7)

    def test_primal_and_woodbury_solvers_satisfy_kkt_and_agree(self):
        rng = np.random.default_rng(12345)
        # Case 1: p <= n (primal Newton path)
        X_small = rng.normal(0.0, 1.0, size=(60, 12))
        y_small = (X_small[:, 0] - 0.5 * X_small[:, 1] + rng.normal(0.0, 0.4, 60)) > 0.0
        for lam in (1e-3, 1e-1, 10.0):
            res = CMP.fit_l2_logistic(X_small, y_small, lam)
            self._assert_kkt(X_small, y_small, lam, res)

        # Case 2: p > n (Woodbury Newton path)
        X_wide = rng.normal(0.0, 1.0, size=(28, 120))
        y_wide = (X_wide[:, 0] + X_wide[:, 5] > 0.0)
        for lam in (1e-2, 1.0, 50.0):
            res_wide = CMP.fit_l2_logistic(X_wide, y_wide, lam)
            self._assert_kkt(X_wide, y_wide, lam, res_wide)

        # Case 3: direct numerical cross-check between primal and Woodbury at boundary
        X_eq = rng.normal(0.0, 1.0, size=(20, 21))
        y_eq = (X_eq[:, 0] > 0.0)
        res_wb = CMP.fit_l2_logistic(X_eq, y_eq, 0.5)
        self._assert_kkt(X_eq, y_eq, 0.5, res_wb)
        # Pad 1 dummy row? Instead check exact optimality against padded duplicate-free primal:
        # Compare against 20-col matrix augmented with a zero column (which forces Woodbury p=21 > n=20).
        X_20 = X_eq[:, :20]
        res_primal = CMP.fit_l2_logistic(X_20, y_eq, 0.5)
        X_20_pad = np.hstack([X_20, np.zeros((20, 1), dtype=np.float64)])
        res_woodbury = CMP.fit_l2_logistic(X_20_pad, y_eq, 0.5)
        np.testing.assert_allclose(res_primal["coef"], res_woodbury["coef"][:20], atol=1e-10, rtol=1e-10)
        self.assertAlmostEqual(res_primal["intercept"], res_woodbury["intercept"], places=10)
        self.assertEqual(res_woodbury["coef"][20], 0.0)

    def test_zero_feature_matrix_returns_null_fit(self):
        X = np.zeros((16, 0), dtype=np.float64)
        y = np.array([i % 2 == 0 for i in range(16)], dtype=bool)
        res = CMP.fit_l2_logistic(X, y, 1.0)
        self.assertTrue(res["null_fit"])
        self.assertTrue(res["converged"])


class CompareEndToEndTests(unittest.TestCase):
    def test_multi_arm_comparison_on_lineage_and_study_whole_with_hygiene_scan(self):
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            fx = build_synthetic_dataset(td_path)
            export = CMP.ExportDataset(fx["export_dir"])

            for scheme, fold_file in (("lineage", None), ("study-whole", fx["fold_file"])):
                with self.subTest(scheme=scheme):
                    _, _, floor, _ = CMP.load_folds(export, scheme, 42, fold_file=fold_file)
                    arm_kmer = CMP.make_kmer31_arm(export, floor=floor)
                    arm_var = CMP.make_protein_variant_arm(
                        fx["isolates_variants"],
                        target_genes=fx["target_genes"],
                    )
                    arm_burden = CMP.make_burden_arm(
                        fx["isolates_variants"],
                        target_genes=fx["target_genes"],
                        regions=fx["regions"],
                    )
                    adapter = FT.FoundationScoreAdapter(
                        fx["foundation_file"],
                        export_meta=export.meta,
                        fold_doc=fx["fold_doc"] if fold_file is not None else None,
                        pca_components=4,
                    )
                    arm_esm = CMP.make_foundation_arm(adapter)
                    arm_hybrid = CMP.make_concat_arm("hybrid_kmer31_esm2", [arm_kmer, arm_esm])

                    report, _ = CMP.compare_representations(
                        export,
                        scheme,
                        42,
                        [arm_kmer, arm_var, arm_burden, arm_esm, arm_hybrid],
                        fold_file=fold_file,
                    )
                    self.assertEqual(report["format"], "grm-rep-compare-v1")
                    self.assertEqual(report["evaluation_revision"], 14)
                    self.assertEqual(report["baseline_arm"], "kmer31")
                    self.assertEqual(len(report["arms"]), 5)
                    self.assertEqual(len(report["pairwise_vs_baseline"]), 4)

                    by_name = {a["arm"]: a for a in report["arms"]}
                    self.assertEqual(by_name["kmer31"]["claim_ceiling"], "retrospective_internal_split")
                    self.assertEqual(by_name["foundation_esm2"]["claim_ceiling"], "retrospective_only")
                    self.assertEqual(by_name["hybrid_kmer31_esm2"]["claim_ceiling"], "retrospective_only")
                    self.assertEqual(len(by_name["hybrid_kmer31_esm2"]["exposure_audits"]), 1)

                    # Every arm scores exact same tested denominator.
                    for arm in report["arms"]:
                        self.assertEqual(arm["pooled_metrics"]["n"], report["tested"])
                        # Recompute pooled BA from fold calls.
                        tp = sum(f["calls"]["tp"] for f in arm["folds"])
                        fn = sum(f["calls"]["fn"] for f in arm["folds"])
                        tn = sum(f["calls"]["tn"] for f in arm["folds"])
                        fp = sum(f["calls"]["fp"] for f in arm["folds"])
                        ba = 0.5 * (tp / (tp + fn) + tn / (tn + fp))
                        self.assertAlmostEqual(arm["pooled_metrics"]["balanced_accuracy"], ba, places=12)

                    # Verify report passes scripts/c7-outside/scan.py with 0 genome IDs.
                    out_file = td_path / f"report_{scheme}.json"
                    CMP.write_atomic(out_file, CMP.dump_report(report))
                    ids = C7_SCAN.ids_of(str(fx["species_file"]))
                    hits, n_files = C7_SCAN.scan(ids, [str(out_file)])
                    self.assertEqual(hits, 0)
                    self.assertEqual(n_files, 1)

    def test_cli_determinism_and_refusal_on_corrupted_plan_sha(self):
        script = SCRIPTS_DIR / "compare.py"
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            fx = build_synthetic_dataset(td_path)
            out1 = td_path / "run1.json"
            out2 = td_path / "run2.json"
            for out_path in (out1, out2):
                done = subprocess.run(
                    [
                        sys.executable,
                        "-I",
                        "-B",
                        str(script),
                        "run",
                        "--export",
                        str(fx["export_dir"]),
                        "--scheme",
                        "lineage",
                        "--seed",
                        "42",
                        "--variants",
                        str(fx["variants_file"]),
                        "--foundation",
                        str(fx["foundation_file"]),
                        "--pca-components",
                        "4",
                        "--out",
                        str(out_path),
                    ],
                    capture_output=True,
                    text=True,
                    check=False,
                )
                self.assertEqual(done.returncode, 0, done.stderr)
            self.assertEqual(out1.read_bytes(), out2.read_bytes())

            # Run scripts/c7-outside/scan.py CLI on the generated report.
            scan_run = subprocess.run(
                [sys.executable, "-I", "-B", str(ROOT_SCAN_PATH()), str(fx["species_file"]), str(out1)],
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(scan_run.returncode, 0, scan_run.stderr)
            self.assertIn("genome ids: 0 in 1 files", scan_run.stdout)

            # Corrupt plan sha256 and confirm refusal.
            meta_path = fx["export_dir"] / "dataset.json"
            bad_meta = json.loads(meta_path.read_text(encoding="utf-8"))
            bad_meta["splits"][0]["plan"]["sha256"] = "0" * 64
            meta_path.write_text(json.dumps(bad_meta), encoding="utf-8")
            done_bad = subprocess.run(
                [
                    sys.executable,
                    "-I",
                    "-B",
                    str(script),
                    "run",
                    "--export",
                    str(fx["export_dir"]),
                    "--scheme",
                    "lineage",
                    "--seed",
                    "42",
                    "--out",
                    str(td_path / "should_fail.json"),
                ],
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(done_bad.returncode, 1)
            self.assertIn("sha256 mismatch", done_bad.stderr)


def ROOT_SCAN_PATH() -> Path:
    from tests.point_mutations.common import SCAN_SCRIPT
    return SCAN_SCRIPT


if __name__ == "__main__":
    unittest.main()
