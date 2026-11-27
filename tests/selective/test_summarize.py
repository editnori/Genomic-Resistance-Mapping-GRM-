"""Unit tests for scripts/selective/summarize.py on SYNTHETIC grm-selective-v1 inputs."""

import importlib.util
import json
import math
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "selective" / "summarize.py"
C7_SCAN = ROOT / "scripts" / "c7-outside" / "scan.py"


def _load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


summarize = _load_module("selective_summarize", SCRIPT)
c7_scan = _load_module("c7_outside_scan", C7_SCAN)


def _rate_ci(x, n):
    sparse = n < summarize.SPARSE_DENOMINATOR_MIN
    if n == 0:
        return {"x": x, "n": n, "rate": None, "ci_lo": None, "ci_hi": None, "sparse_denominator": sparse}
    p = x / n
    return {"x": x, "n": n, "rate": p, "ci_lo": max(0.0, p - 0.05), "ci_hi": min(1.0, p + 0.05), "sparse_denominator": sparse}


def _make_accounting(n_r, n_s, called_r, called_s, fn, fp, reason_key="uncertain_score", conformal=False):
    tp = called_r - fn
    tn = called_s - fp
    abstain_r = n_r - called_r
    abstain_s = n_s - called_s
    abstain_total = abstain_r + abstain_s
    reasons = {
        "null_fit": 0,
        "ood_mash_distance": 0,
        "uncertain_score": 0,
        "conformal_both": 0,
        "conformal_empty": 0,
        "random_thinning": 0,
    }
    reasons[reason_key] = abstain_total
    n = n_r + n_s
    called = called_r + called_s
    sparse_r = called_r < summarize.SPARSE_DENOMINATOR_MIN
    sparse_s = called_s < summarize.SPARSE_DENOMINATOR_MIN
    acc = {
        "n": n,
        "n_r": n_r,
        "n_s": n_s,
        "called": called,
        "called_r": called_r,
        "called_s": called_s,
        "abstain_r": abstain_r,
        "abstain_s": abstain_s,
        "tp": tp,
        "fp": fp,
        "tn": tn,
        "fn": fn,
        "abstain_reasons": reasons,
        "coverage": called / n if n else None,
        "coverage_r": called_r / n_r if n_r else None,
        "coverage_s": called_s / n_s if n_s else None,
        "vme_called": _rate_ci(fn, called_r),
        "me_called": _rate_ci(fp, called_s),
        "npv": _rate_ci(tn, tn + fn),
        "for_s": _rate_ci(fn, tn + fn),
        "ppv": _rate_ci(tp, tp + fp),
        "balanced_accuracy_called": 0.5 * (tp / called_r + tn / called_s) if (called_r and called_s) else None,
        "vme_all": _rate_ci(fn, n_r),
        "me_all": _rate_ci(fp, n_s),
        "per_100_r": {"caught": 100.0 * tp / n_r, "missed": 100.0 * fn / n_r, "deferred": 100.0 * abstain_r / n_r},
        "per_100_s": {"cleared": 100.0 * tn / n_s, "overcalled": 100.0 * fp / n_s, "deferred": 100.0 * abstain_s / n_s},
        "sparse_denominator": sparse_r or sparse_s,
        "sparse_denominator_r": sparse_r,
        "sparse_denominator_s": sparse_s,
    }
    if conformal:
        acc["set_true_label_coverage"] = _rate_ci(tp + tn + abstain_total, n)
    return acc


def _make_selective_doc(
    species_slug,
    antibiotic,
    scheme,
    vme_by_policy,
    *,
    n_r=40,
    n_s=60,
    called_r_at_85=34,
    seed=42,
    distance_available=True,
):
    """Build a valid SYNTHETIC grm-selective-v1 document with exact fold counts."""
    dataset = f"{species_slug}__{antibiotic}"
    n_total = n_r + n_s
    # 2 equal outer folds of size n_total // 2
    assert n_total % 2 == 0 and n_r % 2 == 0 and n_s % 2 == 0
    f_n = n_total // 2
    f_nr = n_r // 2
    f_ns = n_s // 2

    def build_policy(pname):
        if pname in ("mash_distance", "margin_plus_distance") and not distance_available:
            return {
                "policy": pname,
                "available": False,
                "unavailable_reason": "synthetic missing sketches",
                "calibrated_operating_point": [],
                "matched_coverage": [],
            }
        reason_key = {
            "outcome_blind": "random_thinning",
            "score_margin": "uncertain_score",
            "mash_distance": "ood_mash_distance",
            "margin_plus_distance": "uncertain_score",
        }[pname]
        cal_pts = []
        match_pts = []
        target_vme = vme_by_policy[pname]
        for c in summarize.COVERAGE_LADDER:
            keep_fold = summarize.retention_count(c, f_n)
            keep_total = 2 * keep_fold
            if math.isclose(c, summarize.PRIMARY_COVERAGE, abs_tol=1e-12):
                cr = min(called_r_at_85, keep_total, n_r)
            else:
                cr = min(round(n_r * c), keep_total, n_r)
            cs = keep_total - cr
            if cs > n_s:
                cs = n_s
                cr = keep_total - cs
            fn = round(target_vme * cr)
            fp = min(2, cs)
            acc = _make_accounting(n_r, n_s, cr, cs, fn, fp, reason_key=reason_key)
            # Split across 2 folds cleanly so sum(fold) == pooled and fold.called == keep_fold.
            cr0 = min(cr // 2, f_nr, keep_fold)
            cr1 = cr - cr0
            if cr1 > f_nr:
                cr0 += cr1 - f_nr
                cr1 = f_nr
            cs0 = keep_fold - cr0
            cs1 = keep_fold - cr1
            fn0 = min(fn // 2, cr0)
            fn1 = fn - fn0
            fp0 = min(fp // 2, cs0)
            fp1 = fp - fp0
            f0 = _make_accounting(f_nr, f_ns, cr0, cs0, fn0, fp0, reason_key=reason_key)
            f1 = _make_accounting(f_nr, f_ns, cr1, cs1, fn1, fp1, reason_key=reason_key)
            fold_recs = [
                {
                    "fold": idx,
                    "n": fd["n"],
                    "n_r": fd["n_r"],
                    "n_s": fd["n_s"],
                    "called": fd["called"],
                    "called_r": fd["called_r"],
                    "called_s": fd["called_s"],
                    "abstain_r": fd["abstain_r"],
                    "abstain_s": fd["abstain_s"],
                    "tp": fd["tp"],
                    "fp": fd["fp"],
                    "tn": fd["tn"],
                    "fn": fd["fn"],
                    "abstain_reasons": fd["abstain_reasons"],
                }
                for idx, fd in enumerate((f0, f1))
            ]
            pt = {"target_coverage": c, **acc, "folds": fold_recs}
            cal_pts.append(pt)
            match_pts.append(pt)
        return {
            "policy": pname,
            "available": True,
            "unavailable_reason": None,
            "calibrated_operating_point": cal_pts,
            "matched_coverage": match_pts,
        }

    learners = []
    for lname in summarize.LEARNERS:
        learners.append(
            {
                "learner": lname,
                "outer_null_folds": 0,
                "calibrations": [],
                "policies": [build_policy(p) for p in summarize.POLICIES],
                "mondrian_conformal_available": False,
                "mondrian_conformal_unavailable_reason": "synthetic omitted",
                "mondrian_conformal": [],
            }
        )

    return {
        "format": summarize.FORMAT,
        "evaluation_revision": summarize.EVALUATION_REVISION,
        "dataset": dataset,
        "species_slug": species_slug,
        "antibiotic": antibiotic,
        "scheme": scheme,
        "seed": seed,
        "folds": 2,
        "inner_folds": 2,
        "min_support": 3,
        "samples": n_total,
        "tested": n_total,
        "tested_resistant": n_r,
        "tested_susceptible": n_s,
        "fold_file_sha256": None if scheme in ("random", "lineage") else "0" * 64,
        "unavailable_reason": None,
        "constants": {
            "primary_coverage": summarize.PRIMARY_COVERAGE,
            "coverage_ladder": list(summarize.COVERAGE_LADDER),
            "conformal_alphas": [0.01, 0.02, 0.05, 0.1, 0.15, 0.2],
            "sparse_denominator_min": summarize.SPARSE_DENOMINATOR_MIN,
            "ci_alpha": 0.05,
            "joint_l2_lambda": 0.05,
            "tie_rule": "fnv1a64(selective-tie/<seed>/<scheme>/<fold>/<sample_pos>)",
        },
        "distances": {
            "available": distance_available,
            "source": "species_check_cache" if distance_available else None,
            "k": 21 if distance_available else None,
            "sketch_size": 1000 if distance_available else None,
            "unavailable_reason": None if distance_available else "synthetic missing sketches",
        },
        "learners": learners,
        "seconds": 0.1,
    }


class TestSelectiveSummarize(unittest.TestCase):
    def test_exact_one_sided_binomial_sign_test(self):
        # 5 wins, 0 losses -> (1/2)^5 = 0.03125
        self.assertAlmostEqual(summarize.exact_one_sided_sign_test(5, 0), 1 / 32, places=12)
        # 6 wins, 0 losses -> 1/64 = 0.015625
        self.assertAlmostEqual(summarize.exact_one_sided_sign_test(6, 0), 1 / 64, places=12)
        # 5 wins, 1 loss -> (C(6,5)+C(6,6))/64 = 7/64 = 0.109375
        self.assertAlmostEqual(summarize.exact_one_sided_sign_test(5, 1), 7 / 64, places=12)
        # 0 non-ties -> 1.0
        self.assertEqual(summarize.exact_one_sided_sign_test(0, 0), 1.0)

    def test_species_aggregation_sign_test_sparse_filter_and_joint_criterion(self):
        # 6 synthetic species:
        # - sp1..sp5 all show score_margin reducing vme_called vs outcome_blind on BOTH splits
        #   (5/5 wins -> p = 1/32 = 0.03125 < 0.05 -> BOTH pass -> joint_claim_supported = 1).
        # - sp6 is a tie (0 reduction) so it must NOT enter the non-tie denominator (still 5/5 -> p = 0.03125).
        # - sp7 has called_r = 24 (< 25) -> sparse_denominator excluded from primary median/sign test,
        #   while raw counts stay recorded in datasets.tsv!
        # - For margin_plus_distance vs score_margin:
        #   lineage has 5/5 wins (passes), but study-whole has only 2 wins and 3 losses (fails)
        #   -> joint_claim_supported MUST be 0!
        docs = []
        for idx in range(1, 6):
            sp = f"species_{idx}"
            # Two drugs in sp1 to test within-species median aggregation (still counts as 1 species!)
            drugs = ["drug_a", "drug_b"] if idx == 1 else ["drug_a"]
            for drug in drugs:
                docs.append(
                    _make_selective_doc(
                        sp,
                        drug,
                        "lineage",
                        {
                            "outcome_blind": 0.30,
                            "score_margin": 0.15,
                            "mash_distance": 0.25,
                            "margin_plus_distance": 0.05,
                        },
                        n_r=40,
                        n_s=60,
                        called_r_at_85=34,
                    )
                )
                mpd_study = 0.05 if idx <= 2 else 0.25  # 2 wins, 3 losses on study-whole
                docs.append(
                    _make_selective_doc(
                        sp,
                        drug,
                        "study-whole",
                        {
                            "outcome_blind": 0.30,
                            "score_margin": 0.15,
                            "mash_distance": 0.25,
                            "margin_plus_distance": mpd_study,
                        },
                        n_r=40,
                        n_s=60,
                        called_r_at_85=34,
                    )
                )

        # sp6: exact tie on both contrasts and both splits.
        for sch in ("lineage", "study-whole"):
            docs.append(
                _make_selective_doc(
                    "species_6",
                    "drug_a",
                    sch,
                    {
                        "outcome_blind": 0.20,
                        "score_margin": 0.20,
                        "mash_distance": 0.20,
                        "margin_plus_distance": 0.20,
                    },
                    n_r=40,
                    n_s=60,
                    called_r_at_85=34,
                )
            )

        # sp7: sparse denominator (called_r = 24 < 25) with huge adversarial loss; MUST be excluded!
        for sch in ("lineage", "study-whole"):
            docs.append(
                _make_selective_doc(
                    "species_7",
                    "drug_a",
                    sch,
                    {
                        "outcome_blind": 0.00,
                        "score_margin": 0.50,
                        "mash_distance": 0.50,
                        "margin_plus_distance": 0.75,
                    },
                    n_r=30,
                    n_s=70,
                    called_r_at_85=24,
                )
            )

        for d in docs:
            summarize.validate_doc(d)

        res = summarize.summarize(docs, alpha=0.05)

        # Check explicit sign-test direction-not-magnitude annotation.
        self.assertIn("direction", res["sign_test_note"])
        self.assertIn("not mean magnitude", res["sign_test_note"])

        # Check sparse dataset sp7 stays recorded in datasets with raw counts and excluded_sparse=1.
        sp7_rows = [
            r
            for r in res["datasets"]
            if r["species_slug"] == "species_7"
            and r["layer"] == "matched_coverage"
            and r["learner"] == "rf"
            and r["contrast"] == "score_margin_vs_outcome_blind"
        ]
        self.assertEqual(len(sp7_rows), 2)
        for row in sp7_rows:
            self.assertEqual(row["excluded_sparse"], 1)
            self.assertEqual(row["eligible"], 0)
            self.assertEqual(row["called_r_arm"], 24)
            self.assertEqual(row["n_r"], 30)

        # Check sign test rows on matched_coverage for rf.
        st_map = {
            (r["scheme"], r["contrast"]): r
            for r in res["sign_tests"]
            if r["layer"] == "matched_coverage" and r["learner"] == "rf"
        }
        sm_lin = st_map[("lineage", "score_margin_vs_outcome_blind")]
        self.assertEqual((sm_lin["wins"], sm_lin["ties"], sm_lin["losses"], sm_lin["non_ties"]), (5, 1, 0, 5))
        self.assertEqual(sm_lin["species_eligible"], 6)
        self.assertEqual(sm_lin["datasets_sparse_excluded"], 1)
        self.assertAlmostEqual(sm_lin["sign_test_p"], 1 / 32, places=12)
        self.assertEqual(sm_lin["split_pass"], 1)
        self.assertIn("not mean magnitude", sm_lin["sign_test_note"])

        mpd_lin = st_map[("lineage", "margin_plus_distance_vs_score_margin")]
        mpd_stu = st_map[("study-whole", "margin_plus_distance_vs_score_margin")]
        self.assertEqual(mpd_lin["split_pass"], 1)
        self.assertEqual((mpd_stu["wins"], mpd_stu["ties"], mpd_stu["losses"]), (2, 1, 3))
        self.assertEqual(mpd_stu["split_pass"], 0)

        # Check joint cross-split claims:
        # score_margin_vs_outcome_blind passes BOTH lineage AND study-whole -> supported = 1.
        # margin_plus_distance_vs_score_margin passes lineage only -> supported = 0!
        jc_map = {
            r["contrast"]: r
            for r in res["joint_claims"]
            if r["layer"] == "matched_coverage" and r["learner"] == "rf"
        }
        self.assertEqual(jc_map["score_margin_vs_outcome_blind"]["joint_claim_supported"], 1)
        self.assertEqual(jc_map["margin_plus_distance_vs_score_margin"]["joint_claim_supported"], 0)
        self.assertIn("study-whole did not pass", jc_map["margin_plus_distance_vs_score_margin"]["reason"])

    def test_missing_study_whole_blocks_joint_claim(self):
        docs = [
            _make_selective_doc(
                f"sp_{i}",
                "drug",
                "lineage",
                {"outcome_blind": 0.30, "score_margin": 0.10, "mash_distance": 0.25, "margin_plus_distance": 0.05},
            )
            for i in range(6)
        ]
        res = summarize.summarize(docs, alpha=0.05)
        jc = next(
            r
            for r in res["joint_claims"]
            if r["layer"] == "matched_coverage"
            and r["learner"] == "gbt"
            and r["contrast"] == "score_margin_vs_outcome_blind"
        )
        self.assertEqual(jc["lineage_pass"], 1)
        self.assertEqual(jc["study_whole_present"], 0)
        self.assertEqual(jc["joint_claim_supported"], 0)
        self.assertIn("study-whole split missing", jc["reason"])

    def test_validation_rejects_tampered_counts_and_wrong_matched_retention(self):
        doc = _make_selective_doc(
            "sp_1",
            "drug",
            "lineage",
            {"outcome_blind": 0.30, "score_margin": 0.15, "mash_distance": 0.25, "margin_plus_distance": 0.10},
        )
        # Tamper with matched_coverage fold retention count.
        bad = json.loads(json.dumps(doc))
        bad["learners"][0]["policies"][0]["matched_coverage"][3]["folds"][0]["called"] += 1
        with self.assertRaises(ValueError):
            summarize.validate_doc(bad)

        # Tamper with sparse_denominator flag.
        bad2 = json.loads(json.dumps(doc))
        bad2["learners"][0]["policies"][0]["matched_coverage"][3]["vme_called"]["sparse_denominator"] = True
        with self.assertRaises(ValueError):
            summarize.validate_doc(bad2)

    def test_cli_determinism_and_zero_id_scan_matches_c7_outside_scan(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            in_dir = root / "in"
            in_dir.mkdir()
            out_a = root / "out_a"
            out_b = root / "out_b"

            for i in range(1, 7):
                sp = f"escherichia_coli_{i}"
                for sch in ("lineage", "study-whole"):
                    doc = _make_selective_doc(
                        sp,
                        "meropenem",
                        sch,
                        {"outcome_blind": 0.30, "score_margin": 0.15, "mash_distance": 0.20, "margin_plus_distance": 0.10},
                    )
                    (in_dir / f"{sp}__meropenem.{sch}-s42.selective.json").write_text(
                        json.dumps(doc), encoding="utf-8"
                    )

            # Synthetic species JSON holding realistic BV-BRC-shaped genome IDs.
            species_json = root / "species.json"
            fake_ids = [f"562.{80000 + i}" for i in range(50)]
            species_json.write_text(
                json.dumps({"genomes": [{"id": gid} for gid in fake_ids]}), encoding="utf-8"
            )

            # Run CLI twice and verify byte-for-byte determinism and zero genome ID hits.
            for out_dir in (out_a, out_b):
                proc = subprocess.run(
                    [
                        sys.executable,
                        "-I",
                        "-B",
                        "-W",
                        "error",
                        str(SCRIPT),
                        str(in_dir),
                        "--out",
                        str(out_dir),
                        "--scan-species-json",
                        str(species_json),
                    ],
                    capture_output=True,
                    text=True,
                    check=False,
                )
                self.assertEqual(proc.returncode, 0, proc.stderr + proc.stdout)
                self.assertIn("genome ids: 0 in ", proc.stdout)

            # Confirm scripts/c7-outside/scan.py directly agrees: 0 hits on all outputs.
            hits, n_files = c7_scan.scan(set(fake_ids), [str(in_dir), str(out_a)])
            self.assertEqual(hits, 0)
            self.assertGreater(n_files, 0)

            # Byte-identical outputs across runs.
            for name in ("datasets.tsv", "species_differences.tsv", "sign_tests.tsv", "joint_claims.tsv", "summary.json"):
                self.assertEqual((out_a / name).read_bytes(), (out_b / name).read_bytes(), name)

            # Planted-fault test: planting a genome ID into an output file MUST fail the scan!
            leaked = root / "leaked"
            leaked.mkdir()
            (leaked / "bad.tsv").write_text(f"note\t{fake_ids[0]}\n", encoding="utf-8")
            hits_bad, _ = summarize.scan_ids_against_species_json([species_json], [leaked])
            self.assertEqual(hits_bad, 1)


if __name__ == "__main__":
    unittest.main()
