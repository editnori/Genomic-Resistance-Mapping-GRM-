"""Unit tests for Track B three-exposure pretraining audit (`exposure_audit.py`)."""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from tests.point_mutations.common import EA, SCRIPTS_DIR, sample_exposure_audit


class ExposureAuditTests(unittest.TestCase):
    def test_allowed_prokaryotic_model_families_pass(self):
        for family in ("esm2", "evo2", "prokbert", "bacformer"):
            with self.subTest(family=family):
                rec = sample_exposure_audit(
                    model_family=family,
                    sequence_overlap_status="verified_disjoint",
                    phenotype_annotation_status="none_verified",
                    eval_tuning_status="frozen_preregistered",
                    checked_at_utc="2026-11-27T00:00:00Z",
                )
                out = EA.validate_exposure_audit(rec)
                self.assertEqual(out["model_family"], family)
                self.assertEqual(out["derived_claim_ceiling"], "certified_unseen_sequence")
                self.assertEqual(out["effective_claim_ceiling"], "certified_unseen_sequence")
                self.assertEqual(out["claim_ceiling_reasons"], [])

    def test_unknown_sequence_overlap_caps_at_retrospective_only(self):
        rec = sample_exposure_audit(
            sequence_overlap_status="unknown",
            phenotype_annotation_status="none_verified",
            eval_tuning_status="frozen_preregistered",
        )
        out = EA.validate_exposure_audit(rec)
        self.assertEqual(out["derived_claim_ceiling"], "retrospective_only")
        self.assertNotEqual(out["derived_claim_ceiling"], "certified_unseen_sequence")
        self.assertTrue(any("unknown" in r for r in out["claim_ceiling_reasons"]))

    def test_partial_and_confirmed_sequence_overlap_cap_at_retrospective_only(self):
        for seq_status in ("partial_overlap", "confirmed_overlap"):
            with self.subTest(seq_status=seq_status):
                rec = sample_exposure_audit(
                    sequence_overlap_status=seq_status,
                    phenotype_annotation_status="none_verified",
                    eval_tuning_status="outer_train_only",
                    checked_at_utc="2026-11-27T00:00:00Z",
                )
                out = EA.validate_exposure_audit(rec)
                self.assertEqual(out["derived_claim_ceiling"], "retrospective_only")

    def test_phenotype_or_catalogue_exposure_caps_at_non_independent_retrospective(self):
        for pheno in ("indirect_catalogue", "supervised_phenotype", "unknown"):
            with self.subTest(pheno=pheno):
                rec = sample_exposure_audit(
                    sequence_overlap_status="verified_disjoint",
                    phenotype_annotation_status=pheno,
                    eval_tuning_status="frozen_preregistered",
                    checked_at_utc="2026-11-27T00:00:00Z",
                )
                out = EA.validate_exposure_audit(rec)
                self.assertEqual(out["derived_claim_ceiling"], "non_independent_retrospective")

    def test_eval_touched_or_unknown_tuning_caps_at_contaminated_exploratory(self):
        for tune in ("eval_touched", "unknown"):
            with self.subTest(tune=tune):
                rec = sample_exposure_audit(
                    sequence_overlap_status="verified_disjoint",
                    phenotype_annotation_status="none_verified",
                    eval_tuning_status=tune,
                    checked_at_utc="2026-11-27T00:00:00Z",
                )
                out = EA.validate_exposure_audit(rec)
                self.assertEqual(out["derived_claim_ceiling"], "contaminated_exploratory")

    def test_overclaiming_ceiling_is_rejected_and_conservative_claim_accepted(self):
        over = sample_exposure_audit(
            sequence_overlap_status="unknown",
            phenotype_annotation_status="none_verified",
            eval_tuning_status="frozen_preregistered",
            claimed_ceiling="certified_unseen_sequence",
        )
        with self.assertRaisesRegex(EA.AuditError, "exceeds derived claim ceiling"):
            EA.validate_exposure_audit(over)

        conservative = sample_exposure_audit(
            sequence_overlap_status="verified_disjoint",
            phenotype_annotation_status="none_verified",
            eval_tuning_status="frozen_preregistered",
            checked_at_utc="2026-11-27T00:00:00Z",
            claimed_ceiling="retrospective_only",
        )
        out = EA.validate_exposure_audit(conservative)
        self.assertEqual(out["derived_claim_ceiling"], "certified_unseen_sequence")
        self.assertEqual(out["effective_claim_ceiling"], "retrospective_only")

    def test_alphagenome2_is_explicitly_refused(self):
        for name in ("alphagenome2", "AlphaGenome-2", "AlphaGenome_2", "alphagenome"):
            with self.subTest(name=name):
                rec = sample_exposure_audit(model_family=name)
                with self.assertRaisesRegex(EA.AuditError, "AlphaGenome.*excluded for bacterial"):
                    EA.validate_exposure_audit(rec)

    def test_eukaryotic_taxonomic_scope_is_refused(self):
        for scope in ("eukaryotic_only", "mammalian_only", "human_only"):
            with self.subTest(scope=scope):
                rec = sample_exposure_audit()
                rec["training_corpus"]["taxonomic_scope"] = scope
                with self.assertRaisesRegex(EA.AuditError, "incompatible with bacterial"):
                    EA.validate_exposure_audit(rec)

    def test_missing_extra_or_malformed_fields_are_refused(self):
        with self.assertRaises(EA.AuditError):
            EA.validate_exposure_audit(None)

        missing_key = sample_exposure_audit()
        del missing_key["checkpoint_sha256"]
        with self.assertRaisesRegex(EA.AuditError, "missing required keys"):
            EA.validate_exposure_audit(missing_key)

        extra_key = sample_exposure_audit()
        extra_key["unexpected"] = True
        with self.assertRaisesRegex(EA.AuditError, "unlisted keys"):
            EA.validate_exposure_audit(extra_key)

        bad_sha = sample_exposure_audit()
        bad_sha["checkpoint_sha256"] = "not-a-sha256"
        with self.assertRaisesRegex(EA.AuditError, "checkpoint_sha256"):
            EA.validate_exposure_audit(bad_sha)

        unverified_disjoint = sample_exposure_audit(
            sequence_overlap_status="verified_disjoint",
            checked_at_utc="unverified",
        )
        with self.assertRaisesRegex(EA.AuditError, "requires a checked_at_utc timestamp"):
            EA.validate_exposure_audit(unverified_disjoint)

    def test_cli_validate_passes_valid_and_refuses_invalid(self):
        script = SCRIPTS_DIR / "exposure_audit.py"
        with tempfile.TemporaryDirectory() as td:
            valid_path = Path(td) / "valid.json"
            valid_path.write_text(json.dumps(sample_exposure_audit()), encoding="utf-8")
            done_ok = subprocess.run(
                [sys.executable, "-I", "-B", str(script), "validate", str(valid_path)],
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(done_ok.returncode, 0, done_ok.stderr)
            parsed = json.loads(done_ok.stdout)
            self.assertEqual(parsed["derived_claim_ceiling"], "retrospective_only")

            bad_path = Path(td) / "bad.json"
            bad_path.write_text(json.dumps(sample_exposure_audit(model_family="alphagenome2")), encoding="utf-8")
            done_bad = subprocess.run(
                [sys.executable, "-I", "-B", str(script), "validate", str(bad_path)],
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(done_bad.returncode, 1)
            self.assertIn("AlphaGenome", done_bad.stderr)


if __name__ == "__main__":
    unittest.main()
