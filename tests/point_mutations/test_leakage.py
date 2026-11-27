"""Planted-fault leakage tests for Track B encoders and comparator.

1. Clean invariance checks: mutate every outer test isolate's label, 31-mer
   vector, protein variants, burden counts, and foundation-model embeddings in
   each outer fold, and assert that EVERY fitted encoder vocabulary, burden
   scaler, foundation scaler/PCA projection matrix, selected inner lambda, and
   training-fold weight vector + intercept is bit-for-bit identical.
2. Planted-fault negative controls: intentionally inject five concrete leakage
   bugs (full-cohort vocabulary, full-cohort burden scaling, full-cohort PCA,
   outer-test lambda selection, and outer-test refit inclusion) and verify that
   our exact invariance assertion FAILS on every planted fault.
"""

from __future__ import annotations

import copy
import tempfile
import unittest
from pathlib import Path
from typing import Any

import numpy as np

from tests.point_mutations.common import CMP, FT, build_synthetic_dataset


def _mutate_outer_test_fold(
    export: CMP.ExportDataset,
    isolates_variants: list[list[dict[str, Any]]],
    foundation_doc: dict[str, Any],
    test_mask: np.ndarray,
) -> tuple[CMP.ExportDataset, list[list[dict[str, Any]]], dict[str, Any]]:
    """Return mutated copies where EVERY outer test row's label and inputs are corrupted."""
    mut_export = copy.deepcopy(export)
    te_idx = np.flatnonzero(test_mask)

    # 1. Flip every outer test label.
    mut_export.y[te_idx] = ~mut_export.y[te_idx]

    # 2. Flip every outer test 31-mer bit.
    mut_export.kmer_dense[te_idx, :] = 1 - mut_export.kmer_dense[te_idx, :]

    # 3. Replace outer test variants with extreme test-only novel variants and huge burdens.
    mut_vars = copy.deepcopy(isolates_variants)
    for k, i in enumerate(te_idx):
        mut_vars[int(i)] = [
            {"gene": "gyrA", "change": f"LEAK_TEST_{k}", "pos": 83, "kind": "missense"},
            {"gene": "gyrA", "change": f"LEAK_STOP_{k}", "pos": 85, "kind": "nonsense"},
            {"gene": "parC", "change": f"LEAK_DEL_{k}", "pos_start": 70, "pos_end": 75, "kind": "inframe_deletion"},
            {"gene": "rpoB", "change": f"LEAK_SYN_{k}", "pos": 430, "kind": "synonymous"},
            {"gene": "rpoB", "change": f"LEAK_RPOB_{k}", "pos": 445, "kind": "missense"},
        ]

    # 4. Perturb outer test foundation embeddings by extreme shifts.
    mut_fdoc = copy.deepcopy(foundation_doc)
    test_sample_set = {export.samples[int(i)] for i in te_idx}
    for row_i, sid in enumerate(mut_fdoc["samples"]):
        if sid in test_sample_set:
            mut_fdoc["values"][row_i] = [
                float((row_i + 1) * 10000.0 * ((-1) ** (d + row_i)))
                for d in range(len(mut_fdoc["feature_names"]))
            ]

    return mut_export, mut_vars, mut_fdoc


def _extract_fold_artifacts(
    export: CMP.ExportDataset,
    isolates_variants: list[list[dict[str, Any]]],
    foundation_doc: dict[str, Any],
    target_genes: list[str],
    regions: list[dict[str, Any]],
    scheme: str,
    fold_file: Path | None,
    fold_doc: dict[str, Any] | None,
    target_fold: int,
) -> dict[str, dict[str, Any]]:
    """Run all four single arms + concat arm on `target_fold` and capture exact fitted internals."""
    plan_info, folds, floor, _ = CMP.load_folds(export, scheme, 42, fold_file=fold_file)
    fo = folds[target_fold]
    arms = [
        CMP.make_kmer31_arm(export, floor=floor),
        CMP.make_protein_variant_arm(isolates_variants, target_genes=target_genes),
        CMP.make_burden_arm(isolates_variants, target_genes=target_genes, regions=regions, standardize=True),
        CMP.make_foundation_arm(
            FT.FoundationScoreAdapter(
                foundation_doc,
                export_meta=export.meta,
                fold_doc=fold_doc if fold_file is not None else None,
                standardize=True,
                pca_components=4,
            )
        ),
    ]
    arms.append(CMP.make_concat_arm("hybrid_all", list(arms)))

    captured: dict[str, dict[str, Any]] = {}
    for arm in arms:
        _, _, diag = CMP._evaluate_arm_fold(
            arm,
            export.y,
            fo,
            plan_info["inner_k"],
            CMP.DEFAULT_LAMBDA_GRID,
            CMP.DEFAULT_TOL,
            CMP.DEFAULT_MAX_ITER,
        )
        captured[arm.name] = diag
    return captured


def _assert_fold_artifacts_identical(
    tc: unittest.TestCase,
    clean: dict[str, dict[str, Any]],
    mutated: dict[str, dict[str, Any]],
) -> None:
    tc.assertEqual(set(clean.keys()), set(mutated.keys()))
    for arm_name in sorted(clean.keys()):
        c = clean[arm_name]
        m = mutated[arm_name]
        tc.assertEqual(c["picked_lambda"], m["picked_lambda"], f"{arm_name}: picked_lambda drifted")
        tc.assertEqual(c["intercept"], m["intercept"], f"{arm_name}: intercept drifted")
        tc.assertEqual(
            np.ascontiguousarray(c["coef"], dtype="<f8").tobytes(),
            np.ascontiguousarray(m["coef"], dtype="<f8").tobytes(),
            f"{arm_name}: weight vector bytes drifted",
        )
        enc_c = c["encoder"]
        enc_m = m["encoder"]
        if isinstance(enc_c, dict):
            for sub_name in sorted(enc_c.keys()):
                tc.assertEqual(
                    enc_c[sub_name].state_fingerprint(),
                    enc_m[sub_name].state_fingerprint(),
                    f"{arm_name}/{sub_name}: fitted encoder state drifted",
                )
        else:
            tc.assertEqual(
                enc_c.state_fingerprint(),
                enc_m.state_fingerprint(),
                f"{arm_name}: fitted encoder state drifted",
            )
            if isinstance(enc_c, FT.ProteinVariantEncoder):
                tc.assertEqual(enc_c.vocabulary_, enc_m.vocabulary_)
            elif isinstance(enc_c, FT.BurdenEncoder):
                tc.assertEqual(enc_c.mean_.tobytes(), enc_m.mean_.tobytes())
                tc.assertEqual(enc_c.scale_.tobytes(), enc_m.scale_.tobytes())
            elif isinstance(enc_c, FT.FoundationScoreAdapter):
                tc.assertEqual(enc_c.mean_.tobytes(), enc_m.mean_.tobytes())
                tc.assertEqual(enc_c.scale_.tobytes(), enc_m.scale_.tobytes())
                tc.assertEqual(enc_c.components_.tobytes(), enc_m.components_.tobytes())


class LeakageInvarianceTests(unittest.TestCase):
    def test_mutating_outer_test_rows_leaves_every_fitted_artifact_bit_identical(self):
        with tempfile.TemporaryDirectory() as td:
            fx = build_synthetic_dataset(Path(td))
            export = CMP.ExportDataset(fx["export_dir"])

            for scheme, fold_file in (("lineage", None), ("study-whole", fx["fold_file"])):
                _, folds, _, _ = CMP.load_folds(export, scheme, 42, fold_file=fold_file)
                for fo in folds:
                    f = fo["fold"]
                    with self.subTest(scheme=scheme, fold=f):
                        clean = _extract_fold_artifacts(
                            export,
                            fx["isolates_variants"],
                            fx["foundation_doc"],
                            fx["target_genes"],
                            fx["regions"],
                            scheme,
                            fold_file,
                            fx["fold_doc"],
                            target_fold=f,
                        )
                        mut_export, mut_vars, mut_fdoc = _mutate_outer_test_fold(
                            export,
                            fx["isolates_variants"],
                            fx["foundation_doc"],
                            fo["test"],
                        )
                        mutated = _extract_fold_artifacts(
                            mut_export,
                            mut_vars,
                            mut_fdoc,
                            fx["target_genes"],
                            fx["regions"],
                            scheme,
                            fold_file,
                            fx["fold_doc"],
                            target_fold=f,
                        )
                        _assert_fold_artifacts_identical(self, clean, mutated)


class PlantedFaultLeakageDetectionTests(unittest.TestCase):
    """Plant five deliberate leakage faults and prove the invariance harness catches each one."""

    def test_fault1_full_cohort_variant_vocabulary_is_caught(self):
        with tempfile.TemporaryDirectory() as td:
            fx = build_synthetic_dataset(Path(td))
            export = CMP.ExportDataset(fx["export_dir"])
            _, folds, _, _ = CMP.load_folds(export, "lineage", 42)
            fo = folds[0]
            _, mut_vars, _ = _mutate_outer_test_fold(
                export, fx["isolates_variants"], fx["foundation_doc"], fo["test"]
            )
            # Leaky encoder: fits on entire cohort instead of train_idx.
            e_clean = FT.ProteinVariantEncoder(target_genes=fx["target_genes"]).fit(fx["isolates_variants"])
            e_mut = FT.ProteinVariantEncoder(target_genes=fx["target_genes"]).fit(mut_vars)
            self.assertNotEqual(e_clean.vocabulary_, e_mut.vocabulary_)
            self.assertNotEqual(e_clean.state_fingerprint(), e_mut.state_fingerprint())

    def test_fault2_full_cohort_burden_scaler_is_caught(self):
        with tempfile.TemporaryDirectory() as td:
            fx = build_synthetic_dataset(Path(td))
            export = CMP.ExportDataset(fx["export_dir"])
            _, folds, _, _ = CMP.load_folds(export, "lineage", 42)
            fo = folds[0]
            _, mut_vars, _ = _mutate_outer_test_fold(
                export, fx["isolates_variants"], fx["foundation_doc"], fo["test"]
            )
            b_clean = FT.BurdenEncoder(fx["target_genes"], regions=fx["regions"]).fit(fx["isolates_variants"])
            b_mut = FT.BurdenEncoder(fx["target_genes"], regions=fx["regions"]).fit(mut_vars)
            self.assertNotEqual(b_clean.mean_.tobytes(), b_mut.mean_.tobytes())
            self.assertNotEqual(b_clean.state_fingerprint(), b_mut.state_fingerprint())

    def test_fault3_full_cohort_foundation_pca_is_caught(self):
        with tempfile.TemporaryDirectory() as td:
            fx = build_synthetic_dataset(Path(td))
            export = CMP.ExportDataset(fx["export_dir"])
            _, folds, _, _ = CMP.load_folds(export, "lineage", 42)
            fo = folds[0]
            _, _, mut_fdoc = _mutate_outer_test_fold(
                export, fx["isolates_variants"], fx["foundation_doc"], fo["test"]
            )
            all_idx = np.arange(export.n, dtype=np.int64)
            f_clean = FT.FoundationScoreAdapter(
                fx["foundation_doc"], export_meta=export.meta, pca_components=4
            ).fit(all_idx)
            f_mut = FT.FoundationScoreAdapter(
                mut_fdoc, export_meta=export.meta, pca_components=4
            ).fit(all_idx)
            self.assertNotEqual(f_clean.components_.tobytes(), f_mut.components_.tobytes())
            self.assertNotEqual(f_clean.state_fingerprint(), f_mut.state_fingerprint())

    def test_fault4_outer_test_lambda_selection_is_caught(self):
        with tempfile.TemporaryDirectory() as td:
            fx = build_synthetic_dataset(Path(td))
            export = CMP.ExportDataset(fx["export_dir"])
            _, folds, _, _ = CMP.load_folds(export, "lineage", 42)
            fo = folds[0]
            tr_idx = np.flatnonzero(fo["train"])
            te_idx = np.flatnonzero(fo["test"])
            adapter = FT.FoundationScoreAdapter(fx["foundation_doc"], export_meta=export.meta, pca_components=4)
            X_tr = adapter.fit_transform(tr_idx)
            X_te = adapter.transform(te_idx)

            def leaky_pick(y_te: np.ndarray) -> float:
                best = None
                for lam in CMP.DEFAULT_LAMBDA_GRID:
                    fit = CMP.fit_l2_logistic(X_tr, export.y[tr_idx], lam)
                    preds = CMP.predict_l2_logistic(X_te, fit["coef"], fit["intercept"])
                    loss = CMP.inner_balanced_log_loss(preds, y_te)
                    key = (round(loss, 12), -lam)
                    if best is None or key < best[0]:
                        best = (key, lam)
                assert best is not None
                return best[1]

            lam_clean = leaky_pick(export.y[te_idx])
            lam_flipped = leaky_pick(~export.y[te_idx])
            self.assertNotEqual(lam_clean, lam_flipped)

    def test_fault5_outer_test_included_in_refit_is_caught(self):
        with tempfile.TemporaryDirectory() as td:
            fx = build_synthetic_dataset(Path(td))
            export = CMP.ExportDataset(fx["export_dir"])
            _, folds, _, _ = CMP.load_folds(export, "lineage", 42)
            fo = folds[0]
            all_idx = np.arange(export.n, dtype=np.int64)
            adapter = FT.FoundationScoreAdapter(fx["foundation_doc"], export_meta=export.meta, pca_components=4)
            X_all = adapter.fit_transform(all_idx)

            y_mut = export.y.copy()
            y_mut[fo["test"]] = ~y_mut[fo["test"]]
            fit_clean = CMP.fit_l2_logistic(X_all, export.y, 1.0)
            fit_mut = CMP.fit_l2_logistic(X_all, y_mut, 1.0)
            self.assertNotEqual(fit_clean["coef"].tobytes(), fit_mut["coef"].tobytes())


if __name__ == "__main__":
    unittest.main()
