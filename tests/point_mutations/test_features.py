"""Unit tests for fold-leak-free feature builders (`features.py`)."""

from __future__ import annotations

import copy
import hashlib
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

from tests.point_mutations.common import EA, FT, build_synthetic_dataset, sample_exposure_audit


class ProteinVariantEncoderTests(unittest.TestCase):
    def test_train_only_vocabulary_and_unseen_variant_count(self):
        train = [
            [
                {"gene": "gyrA", "change": "S83L", "pos": 83, "kind": "missense"},
                {"gene": "gyrA", "change": "S83L", "pos": 83, "kind": "missense"},  # duplicate ignored
                {"gene": "gyrA", "change": "A90A", "pos": 90, "kind": "synonymous"},  # excluded by default
                {"gene": "offTarget", "change": "K1R", "pos": 1, "kind": "missense"},  # outside target_genes
            ],
            [
                {"gene": "parC", "change": "S80I", "pos": 80, "kind": "missense"},
            ],
        ]
        test = [
            [
                {"gene": "gyrA", "change": "S83L", "pos": 83, "kind": "missense"},
                {"gene": "gyrA", "change": "D87N", "pos": 87, "kind": "missense"},  # unseen 1
                {"gene": "parC", "change": "E84K", "pos": 84, "kind": "missense"},  # unseen 2
                {"gene": "parC", "change": "E84K", "pos": 84, "kind": "missense"},  # dup unseen counted once
                {"gene": "gyrA", "change": "V99V", "pos": 99, "kind": "synonymous"},  # ignored
            ],
            [],
        ]
        enc = FT.ProteinVariantEncoder(target_genes=["gyrA", "parC"])
        X_tr = enc.fit_transform(train)
        self.assertEqual(enc.vocabulary_, (("gyrA", "S83L"), ("parC", "S80I")))
        self.assertEqual(enc.feature_names_, ("gyrA:S83L", "parC:S80I", "unseen_variant_count"))
        np.testing.assert_array_equal(X_tr, np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]))

        vocab_before = enc.vocabulary_
        X_te = enc.transform(test)
        self.assertEqual(enc.vocabulary_, vocab_before)
        np.testing.assert_array_equal(X_te, np.array([[1.0, 0.0, 2.0], [0.0, 0.0, 0.0]]))

    def test_catalogue_only_and_union_modes_and_support_bounds(self):
        train = [
            [{"gene": "rpoB", "change": "S450L", "kind": "missense"}],
            [{"gene": "rpoB", "change": "S450L", "kind": "missense"}],
            [{"gene": "rpoB", "change": "H445Y", "kind": "missense"}],
        ]
        cat_only = FT.ProteinVariantEncoder(
            target_genes=["rpoB"],
            frozen_catalogue_variants=[("rpoB", "D435V")],
            vocabulary_mode="catalogue_only",
        )
        X_cat = cat_only.fit_transform(train)
        self.assertEqual(cat_only.vocabulary_, (("rpoB", "D435V"),))
        # Every isolate carries 1 unseen non-catalogue variant.
        np.testing.assert_array_equal(X_cat, np.array([[0.0, 1.0], [0.0, 1.0], [0.0, 1.0]]))

        union = FT.ProteinVariantEncoder(
            target_genes=["rpoB"],
            frozen_catalogue_variants=[("rpoB", "D435V")],
            vocabulary_mode="train_union_catalogue",
            min_support=2,
        )
        union.fit(train)
        self.assertEqual(union.vocabulary_, (("rpoB", "D435V"), ("rpoB", "S450L")))

    def test_transform_before_fit_raises(self):
        enc = FT.ProteinVariantEncoder()
        with self.assertRaises(FT.FeatureError):
            enc.transform([[]])


class BurdenEncoderTests(unittest.TestCase):
    def test_gene_and_region_burden_counts_and_train_only_standardization(self):
        regions = [
            {"gene": "rpoB", "region": "RRDR", "start": 426, "end": 452},
        ]
        train = [
            [
                {"gene": "rpoB", "change": "S450L", "pos": 450, "kind": "missense"},
                {"gene": "rpoB", "change": "A100A", "pos": 100, "kind": "synonymous"},
            ],
            [
                {"gene": "rpoB", "change": "del_451_455", "pos_start": 451, "pos_end": 455, "kind": "inframe_deletion"},
                {"gene": "rpoB", "change": "Q432*", "pos": 432, "kind": "nonsense"},
            ],
        ]
        enc_raw = FT.BurdenEncoder(target_genes=["rpoB"], regions=regions, standardize=False)
        X_raw = enc_raw.fit_transform(train)
        self.assertEqual(
            enc_raw.feature_names_,
            (
                "gene:rpoB:synonymous",
                "gene:rpoB:nonsynonymous",
                "gene:rpoB:truncating",
                "gene:rpoB:indel",
                "region:rpoB:RRDR:synonymous",
                "region:rpoB:RRDR:nonsynonymous",
                "region:rpoB:RRDR:truncating",
                "region:rpoB:RRDR:indel",
            ),
        )
        # Row 0: gene syn=1, nonsyn=1; region nonsyn=1 (A100A outside RRDR)
        # Row 1: gene trunc=1, indel=1; region trunc=1, indel=1 (span 451..455 overlaps 426..452)
        expected = np.array([
            [1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0],
        ])
        np.testing.assert_array_equal(X_raw, expected)

        enc_std = FT.BurdenEncoder(target_genes=["rpoB"], regions=regions, standardize=True)
        X_std = enc_std.fit_transform(train)
        # Constant column (region:rpoB:RRDR:synonymous) gets scale 0 -> 0.0 even on test with huge count.
        self.assertEqual(enc_std.scale_[4], 0.0)
        self.assertTrue(np.all(np.isfinite(X_std)))
        test_mut = [[{"gene": "rpoB", "change": "L440L", "pos": 440, "kind": "synonymous"}]]
        X_te = enc_std.transform(test_mut)
        self.assertEqual(X_te[0, 4], 0.0)

    def test_invalid_region_or_unknown_kind_raises(self):
        with self.assertRaises(FT.FeatureError):
            FT.BurdenEncoder(target_genes=["gyrA"], regions=[{"gene": "rpoB", "region": "RRDR", "start": 426, "end": 452}])
        with self.assertRaises(FT.FeatureError):
            FT.BurdenEncoder(target_genes=["gyrA"], regions=[{"gene": "gyrA", "region": "QRDR", "start": 106, "end": 67}])
        enc = FT.BurdenEncoder(target_genes=["gyrA"])
        with self.assertRaisesRegex(FT.FeatureError, "known variant kind"):
            enc.fit([[{"gene": "gyrA", "change": "S83L", "pos": 83, "kind": "unknown_kind"}]])


class FoundationScoreAdapterTests(unittest.TestCase):
    def test_missing_or_invalid_exposure_audit_is_refused(self):
        with tempfile.TemporaryDirectory() as td:
            fx = build_synthetic_dataset(Path(td))
            bad = copy.deepcopy(fx["foundation_doc"])
            del bad["exposure_audit"]
            with self.assertRaisesRegex(EA.AuditError, "missing required 'exposure_audit'"):
                FT.FoundationScoreAdapter(bad, export_meta=fx["export_meta"])

            ag2 = copy.deepcopy(fx["foundation_doc"])
            ag2["exposure_audit"]["model_family"] = "alphagenome2"
            with self.assertRaisesRegex(EA.AuditError, "AlphaGenome"):
                FT.FoundationScoreAdapter(ag2, export_meta=fx["export_meta"])

    def test_alignment_to_export_and_fold_file_and_pca_determinism(self):
        with tempfile.TemporaryDirectory() as td:
            fx = build_synthetic_dataset(Path(td))
            adapter = FT.FoundationScoreAdapter(
                fx["foundation_file"],
                export_meta=fx["export_meta"],
                fold_doc=fx["fold_doc"],
                standardize=True,
                pca_components=4,
            )
            self.assertEqual(adapter.samples_, tuple(fx["samples"]))
            self.assertEqual(adapter.claim_ceiling_, "retrospective_only")

            tr_idx = np.arange(0, 30, dtype=np.int64)
            te_idx = np.arange(30, 45, dtype=np.int64)
            X_tr = adapter.fit_transform(tr_idx)
            X_te = adapter.transform(te_idx)
            self.assertEqual(X_tr.shape, (30, 4))
            self.assertEqual(X_te.shape, (15, 4))
            self.assertEqual(adapter.feature_names_, ("pc_0", "pc_1", "pc_2", "pc_3"))

            # Aligning strictly via fold_doc alone works identically.
            adapter_ff = FT.FoundationScoreAdapter(
                fx["foundation_doc"],
                fold_doc=fx["fold_doc"],
                sample_order=fx["samples"],
                standardize=True,
                pca_components=4,
            )
            X_tr_2 = adapter_ff.fit_transform(tr_idx)
            X_te_2 = adapter_ff.transform(te_idx)
            np.testing.assert_array_equal(X_tr, X_tr_2)
            np.testing.assert_array_equal(X_te, X_te_2)

    def test_npy_loading_and_sha256_verification(self):
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            fx = build_synthetic_dataset(td_path)
            arr = np.asarray(fx["foundation_doc"]["values"], dtype=np.float64)
            npy_path = td_path / "emb.npy"
            np.save(npy_path, arr)
            sha = hashlib.sha256(npy_path.read_bytes()).hexdigest()

            manifest = copy.deepcopy(fx["foundation_doc"])
            del manifest["values"]
            manifest["values_npy"] = "emb.npy"
            manifest["values_sha256"] = sha
            mpath = td_path / "manifest_npy.json"
            mpath.write_text(json.dumps(manifest), encoding="utf-8")

            ad_inline = FT.FoundationScoreAdapter(fx["foundation_doc"], export_meta=fx["export_meta"])
            ad_npy = FT.FoundationScoreAdapter(mpath, export_meta=fx["export_meta"])
            idx = np.arange(len(fx["samples"]), dtype=np.int64)
            np.testing.assert_array_equal(ad_inline.fit_transform(idx), ad_npy.fit_transform(idx))

            # Tampered hash must be refused.
            manifest["values_sha256"] = "0" * 64
            with self.assertRaisesRegex(FT.FeatureError, "sha256 mismatch"):
                FT.FoundationScoreAdapter(manifest, export_meta=fx["export_meta"], base_dir=td_path)

    def test_sample_mismatch_and_missing_value_gating(self):
        with tempfile.TemporaryDirectory() as td:
            fx = build_synthetic_dataset(Path(td))
            bad_samples = copy.deepcopy(fx["foundation_doc"])
            bad_samples["samples"][0] = "999.9999"
            with self.assertRaisesRegex(FT.FeatureError, "sample mismatch"):
                FT.FoundationScoreAdapter(bad_samples, export_meta=fx["export_meta"])

            with_nan = copy.deepcopy(fx["foundation_doc"])
            with_nan["values"][0][0] = None
            with self.assertRaisesRegex(FT.FeatureError, "missing/NaN"):
                FT.FoundationScoreAdapter(with_nan, export_meta=fx["export_meta"], impute_missing=False)

            ad_imp = FT.FoundationScoreAdapter(with_nan, export_meta=fx["export_meta"], impute_missing=True)
            X_imp = ad_imp.fit_transform(np.arange(30, dtype=np.int64))
            self.assertTrue(np.all(np.isfinite(X_imp)))


if __name__ == "__main__":
    unittest.main()
