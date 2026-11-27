"""Synthetic unit tests for unified scoring core (`scripts/bench/score.py`):
- Exact Clopper-Pearson and Wilson 95% CIs
- Overflow-safe exact two-sided McNemar test
- Worst-case / best-case bounds under no_call + failure
- Pairwise common-support (>= 25 R and >= 25 S) + full-denominator coverage
- Strict 4-tier separation (ungraded & proxy never pooled into primary medians)
- 4,000-draw deterministic SHA-256 species & lineage cluster bootstraps
- Time-split provenance enforcement (historical_prospective vs retrospective_modern_db)
- Zero-ID hygiene validation via scripts/c7-outside/scan.py
"""

from __future__ import annotations

import copy
import math
import tempfile
import unittest
from pathlib import Path

from scripts.bench.pretrained import ProjectedEndpoint
from scripts.bench.protocol import (
    BenchError,
    RESULT_FORMAT,
    TIER_PROXY_CLASS,
    TIER_SUPPORTED_PRIMARY,
    TIER_SUPPORTED_UNGRADED,
    TIER_UNSUPPORTED,
    TIME_PROVENANCE_HISTORICAL,
    TIME_PROVENANCE_RETRO_MODERN_DB,
    TIME_PROVENANCE_RETRO_SAME_DB,
    canonical_json_bytes,
    load_adapter,
    validate_time_provenance,
)
from scripts.bench.score import (
    BOOTSTRAP_RESAMPLES,
    bounds_under_uncalled,
    check_zero_id_hygiene,
    clopper_pearson_ci,
    lineage_cluster_bootstrap,
    mcnemar_exact_p,
    score_pretrained_cohort,
    score_trainable_cohort,
    species_cluster_bootstrap,
    wilson_ci,
    write_result_atomic,
)
from scripts.bench.trainable import run_trainable_dataset
from tests.bench.common import build_synthetic_export
from tests.bench.test_trainable import _synthetic_linear_fitter


def _binom_upper_tail(n: int, k: int, p: float) -> float:
    return sum(math.comb(n, i) * (p ** i) * ((1.0 - p) ** (n - i)) for i in range(k, n + 1))


def _binom_lower_tail(n: int, k: int, p: float) -> float:
    return sum(math.comb(n, i) * (p ** i) * ((1.0 - p) ** (n - i)) for i in range(0, k + 1))


def _make_endpoint(
    dataset: str,
    species_slug: str,
    drug: str,
    tier: str,
    *,
    n_r: int = 30,
    n_s: int = 30,
    tp: int = 27,
    tn: int = 28,
    no_call_r: int = 0,
    no_call_s: int = 0,
    i_r: int = 0,
    i_s: int = 0,
    failure_r: int = 0,
    failure_s: int = 0,
    learner_tp: int | None = None,
    learner_tn: int | None = None,
    id_offset: int = 0,
    unsupported_reason: str = "-",
) -> ProjectedEndpoint:
    gids: list[str] = []
    truth: list[bool] = []
    lineages: list[int] = []
    calls: list[str] = []
    l_calls: list[str] = []

    fn = n_r - tp - no_call_r - i_r - failure_r
    fp = n_s - tn - no_call_s - i_s - failure_s
    assert fn >= 0 and fp >= 0

    ltp = tp - 2 if learner_tp is None else learner_tp
    ltn = tn - 1 if learner_tn is None else learner_tn

    r_calls = (
        ["R"] * tp
        + ["S"] * fn
        + ["no_call"] * no_call_r
        + ["I"] * i_r
        + ["failure"] * failure_r
    )
    s_calls = (
        ["S"] * tn
        + ["R"] * fp
        + ["no_call"] * no_call_s
        + ["I"] * i_s
        + ["failure"] * failure_s
    )
    if tier == TIER_UNSUPPORTED:
        r_calls = ["unsupported"] * n_r
        s_calls = ["unsupported"] * n_s

    for j, c in enumerate(r_calls):
        gids.append(f"888000.{id_offset + j}")
        truth.append(True)
        lineages.append(j % 12)
        calls.append(c)
        l_calls.append("R" if j < ltp else "S")
    for j, c in enumerate(s_calls):
        gids.append(f"888000.{id_offset + n_r + j}")
        truth.append(False)
        lineages.append(j % 12)
        calls.append(c)
        l_calls.append("S" if j < ltn else "R")

    return ProjectedEndpoint(
        dataset=dataset,
        species_slug=species_slug,
        canonical_drug=drug,
        split_scheme="lineage",
        selection_role="internal_cv",
        support_tier=tier,
        unsupported_reason=unsupported_reason,
        reading="who_assoc",
        genome_ids=tuple(gids),
        truth=tuple(truth),
        lineages=tuple(lineages),
        calls=tuple(calls),
        fold_of=tuple(0 for _ in gids),
        comparator_calls={"scm:lineage": tuple(l_calls)},
    )


class StatisticalPrimitivesTests(unittest.TestCase):
    def test_clopper_pearson_exact_binomial_tails(self):
        self.assertEqual(clopper_pearson_ci(0, 0), (None, None))
        lo0, hi0 = clopper_pearson_ci(0, 20)
        self.assertEqual(lo0, 0.0)
        self.assertAlmostEqual(hi0, 1.0 - 0.025 ** (1.0 / 20.0), places=11)

        lo1, hi1 = clopper_pearson_ci(20, 20)
        self.assertEqual(hi1, 1.0)
        self.assertAlmostEqual(lo1, 0.025 ** (1.0 / 20.0), places=11)

        for k, n in ((1, 20), (5, 10), (18, 25), (45, 50), (80, 100)):
            lo, hi = clopper_pearson_ci(k, n)
            self.assertAlmostEqual(_binom_upper_tail(n, k, lo), 0.025, places=10)
            self.assertAlmostEqual(_binom_lower_tail(n, k, hi), 0.025, places=10)

        # Large n stability check (n=1000, k=950)
        lo_l, hi_l = clopper_pearson_ci(950, 1000)
        self.assertTrue(0.934 < lo_l < 0.936)
        self.assertTrue(0.962 < hi_l < 0.963)

    def test_wilson_ci_boundaries_and_values(self):
        import importlib.util
        from scripts.bench.vocab import REPO_ROOT

        spec = importlib.util.spec_from_file_location("c9_score_ref", REPO_ROOT / "scripts" / "c9-pick" / "score.py")
        c9 = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(c9)

        self.assertEqual(wilson_ci(0, 0), (None, None))
        for k, n in ((0, 50), (1, 20), (25, 50), (80, 100), (950, 1000), (50, 50)):
            lo, hi = wilson_ci(k, n)
            ref_lo, ref_hi = c9.wilson(k, n)
            self.assertAlmostEqual(lo, ref_lo, places=11)
            self.assertAlmostEqual(hi, ref_hi, places=11)
        self.assertEqual(wilson_ci(0, 50)[0], 0.0)
        self.assertEqual(wilson_ci(50, 50)[1], 1.0)

    def test_mcnemar_exact_small_and_large_overflow_safe(self):
        self.assertEqual(mcnemar_exact_p(0, 0), 1.0)
        self.assertEqual(mcnemar_exact_p(5, 5), 1.0)
        self.assertAlmostEqual(mcnemar_exact_p(10, 0), 2.0 / 1024.0, places=12)
        self.assertAlmostEqual(mcnemar_exact_p(0, 10), 2.0 / 1024.0, places=12)
        # Large discordant sum > 1024 must not overflow int->float!
        p_large_eq = mcnemar_exact_p(600, 600)
        self.assertAlmostEqual(p_large_eq, 1.0, places=9)
        p_large_diff = mcnemar_exact_p(650, 550)
        self.assertTrue(0.0 < p_large_diff < 0.01)

    def test_worst_and_best_bounds_under_no_call_and_failure(self):
        # 40 R called (36 TP, 4 FN) + 10 R uncalled; 40 S called (38 TN, 2 FP) + 10 S uncalled
        b = bounds_under_uncalled(tp=36, fp=2, tn=38, fn=4, uncalled_r=10, uncalled_s=10)
        self.assertAlmostEqual(b["sensitivity_worst"], 36 / 50, places=12)
        self.assertAlmostEqual(b["sensitivity_best"], 46 / 50, places=12)
        self.assertAlmostEqual(b["specificity_worst"], 38 / 50, places=12)
        self.assertAlmostEqual(b["specificity_best"], 48 / 50, places=12)
        self.assertAlmostEqual(b["vme_best"], 4 / 50, places=12)
        self.assertAlmostEqual(b["vme_worst"], 14 / 50, places=12)
        self.assertAlmostEqual(b["me_best"], 2 / 50, places=12)
        self.assertAlmostEqual(b["me_worst"], 12 / 50, places=12)
        # Invariants: worst <= called <= best for sens/spec; best <= called <= worst for VME/ME
        self.assertLessEqual(b["sensitivity_worst"], 36 / 40)
        self.assertLessEqual(36 / 40, b["sensitivity_best"])
        self.assertLessEqual(b["vme_best"], 4 / 40)
        self.assertLessEqual(4 / 40, b["vme_worst"])


class TierSeparationBootstrapAndProvenanceTests(unittest.TestCase):
    def test_min_class_25_and_tier_separation_and_lineage_bootstrap(self):
        adapter = load_adapter("tbprofiler")
        eps = [
            # Eligible primary (28 R called, 28 S called; plus 2 no_call R, 2 failure S)
            _make_endpoint("mycobacterium_tuberculosis__rifampin", "mycobacterium_tuberculosis", "rifampin", TIER_SUPPORTED_PRIMARY, n_r=30, n_s=30, tp=26, tn=27, no_call_r=2, failure_s=2, id_offset=0),
            # Eligible primary (25 R called, 25 S called - exact boundary!)
            _make_endpoint("mycobacterium_tuberculosis__isoniazid", "mycobacterium_tuberculosis", "isoniazid", TIER_SUPPORTED_PRIMARY, n_r=25, n_s=25, tp=23, tn=24, id_offset=100),
            # Below min_class (24 R called, 25 S called -> ineligible!)
            _make_endpoint("mycobacterium_tuberculosis__linezolid", "mycobacterium_tuberculosis", "linezolid", TIER_SUPPORTED_PRIMARY, n_r=26, n_s=25, tp=20, tn=24, no_call_r=2, id_offset=200),
            # Supported ungraded (must NEVER pool into summary_primary!)
            _make_endpoint("mycobacterium_tuberculosis__cycloserine", "mycobacterium_tuberculosis", "cycloserine", TIER_SUPPORTED_UNGRADED, n_r=30, n_s=30, tp=0, tn=30, id_offset=300),
            # Proxy class (must NEVER pool into summary_primary!)
            _make_endpoint("mycobacterium_tuberculosis__ciprofloxacin", "mycobacterium_tuberculosis", "ciprofloxacin", TIER_PROXY_CLASS, n_r=30, n_s=30, tp=27, tn=25, id_offset=400),
            # Unsupported (counted in full-denominator coverage only, never mapped to S!)
            _make_endpoint("mycobacterium_tuberculosis__amoxicillin", "mycobacterium_tuberculosis", "amoxicillin", TIER_UNSUPPORTED, n_r=30, n_s=30, tp=0, tn=0, id_offset=500, unsupported_reason="drug_not_in_tbdb"),
        ]

        prov = {"label": TIME_PROVENANCE_RETRO_SAME_DB}
        known_ids = {g for ep in eps for g in ep.genome_ids}
        doc1 = score_pretrained_cohort(adapter, eps, split_scheme="lineage", time_provenance=prov, known_genome_ids=known_ids)
        doc2 = score_pretrained_cohort(adapter, list(reversed(eps)), split_scheme="lineage", time_provenance=prov, known_genome_ids=known_ids)

        # Byte-for-byte determinism regardless of input list order!
        self.assertEqual(canonical_json_bytes(doc1), canonical_json_bytes(doc2))
        self.assertEqual(doc1["format"], RESULT_FORMAT)
        self.assertEqual(doc1["bootstrap_resamples"], BOOTSTRAP_RESAMPLES)
        self.assertEqual(doc1["bootstrap_cluster_unit"], "lineage")

        # Full-denominator coverage includes ALL 6 datasets and ALL 351 rows
        cov = doc1["full_denominator_coverage"]
        self.assertEqual(cov["total_datasets"], 6)
        self.assertEqual(cov["supported_primary_datasets"], 3)
        self.assertEqual(cov["eligible_primary_datasets"], 2)  # linezolid has 24 R called < 25!
        self.assertEqual(cov["supported_ungraded_datasets"], 1)
        self.assertEqual(cov["proxy_class_datasets"], 1)
        self.assertEqual(cov["unsupported_datasets"], 1)
        self.assertEqual(cov["unsupported_reasons"], {"drug_not_in_tbdb": 1})
        self.assertEqual(cov["unsupported_rows"], 60)
        self.assertEqual(cov["no_call_rows"], 4)
        self.assertEqual(cov["failure_rows"], 2)

        # Primary medians pool ONLY the 2 eligible supported_primary datasets (NOT cycloserine 0.50 BA or ciprofloxacin!)
        self.assertEqual(doc1["summary_primary"]["eligible_endpoints"], 2)
        self.assertGreater(doc1["summary_primary"]["metrics"]["balanced_accuracy"]["median"], 0.90)
        self.assertEqual(doc1["summary_ungraded"]["eligible_endpoints"], 1)
        self.assertAlmostEqual(doc1["summary_ungraded"]["metrics"]["balanced_accuracy"]["median"], 0.50, places=6)
        self.assertEqual(doc1["summary_proxy_class"]["eligible_endpoints"], 1)

        # Lineage bootstrap ran 4,000 draws and produced valid CIs
        lb = doc1["summary_primary"]["lineage_bootstrap"]
        self.assertIsNotNone(lb)
        self.assertEqual(lb["resamples"], 4000)
        self.assertIsNotNone(lb["median_metrics"]["balanced_accuracy"]["median_ci_lo"])

        # Pairwise discordance conservation
        ep_by_id = {e["dataset"]: e for e in doc1["endpoints"]}
        pw = ep_by_id["mycobacterium_tuberculosis__rifampin"]["pairwise_common_support"]["scm:lineage"]
        disc = pw["discordance"]
        self.assertEqual(
            disc["tool_only_right"] + disc["learner_only_right"] + disc["calls_agree"],
            pw["common_called_n"],
        )
        self.assertEqual(ep_by_id["mycobacterium_tuberculosis__linezolid"]["eligible_common_support"], False)
        self.assertEqual(ep_by_id["mycobacterium_tuberculosis__linezolid"]["ineligible_reason"], "below_min_class_25")

    def test_species_cluster_bootstrap_4000_draws_deterministic_and_coarse_flag(self):
        vals_4 = {f"sp_{i}": [0.7 + 0.05 * i, 0.72 + 0.05 * i] for i in range(4)}
        b1 = species_cluster_bootstrap(vals_4, "bench|amrfinderplus|lineage|ba")
        b2 = species_cluster_bootstrap(vals_4, "bench|amrfinderplus|lineage|ba")
        self.assertEqual(b1, b2)
        self.assertEqual(b1["resamples"], 4000)
        self.assertTrue(b1["coarse_cluster_warning"])
        self.assertLessEqual(b1["median_ci_lo"], b1["median"])
        self.assertLessEqual(b1["median"], b1["median_ci_hi"])

        vals_6 = {f"sp_{i}": [0.7 + 0.03 * i] for i in range(6)}
        b6 = species_cluster_bootstrap(vals_6, "bench|amrfinderplus|lineage|ba")
        self.assertFalse(b6["coarse_cluster_warning"])

    def test_time_split_provenance_strict_enforcement(self):
        # Valid historical_prospective: ALL 4 dates <= cutoff
        ok = validate_time_provenance(
            split_scheme="time",
            declared_label=TIME_PROVENANCE_HISTORICAL,
            fold_cutoff_date="2023-01-01",
            software_release_date="2022-06-01",
            database_release_date="2022-12-15",
            weights_release_date="2022-12-15",
            training_labels_max_date="2022-12-31",
            mode="pretrained",
            is_bvbrc_insertion_time_split=True,
        )
        self.assertEqual(ok["label"], TIME_PROVENANCE_HISTORICAL)

        # Post-cutoff DB on BV-BRC insertion-date time split MUST be retrospective_modern_db
        with self.assertRaisesRegex(BenchError, "retrospective_modern_db"):
            validate_time_provenance(
                split_scheme="time",
                declared_label=TIME_PROVENANCE_HISTORICAL,
                fold_cutoff_date="2023-01-01",
                software_release_date="2022-06-01",
                database_release_date="2026-08-07",
                training_labels_max_date="2022-12-31",
                mode="pretrained",
                is_bvbrc_insertion_time_split=True,
            )
        with self.assertRaisesRegex(BenchError, "retrospective_modern_db"):
            validate_time_provenance(
                split_scheme="time",
                declared_label=TIME_PROVENANCE_RETRO_SAME_DB,
                fold_cutoff_date="2023-01-01",
                software_release_date="2026-08-07",
                database_release_date="2026-08-07",
                mode="pretrained",
                is_bvbrc_insertion_time_split=True,
            )

        # Properly labeled retrospective_modern_db passes
        retro = validate_time_provenance(
            split_scheme="time",
            declared_label=TIME_PROVENANCE_RETRO_MODERN_DB,
            fold_cutoff_date="2023-01-01",
            software_release_date="2026-08-07",
            database_release_date="2026-08-07",
            mode="pretrained",
            is_bvbrc_insertion_time_split=True,
        )
        self.assertEqual(retro["label"], TIME_PROVENANCE_RETRO_MODERN_DB)

        # Post-cutoff training labels refused for historical_prospective
        with self.assertRaisesRegex(BenchError, "training_labels_max_date"):
            validate_time_provenance(
                split_scheme="time",
                declared_label=TIME_PROVENANCE_HISTORICAL,
                fold_cutoff_date="2023-01-01",
                software_release_date="2022-01-01",
                training_labels_max_date="2024-05-01",
                mode="trainable",
            )

    def test_zero_id_hygiene_scan_catches_leaks_sequences_accessions_and_timing(self):
        adapter = load_adapter("tbprofiler")
        ep = _make_endpoint("mycobacterium_tuberculosis__rifampin", "mycobacterium_tuberculosis", "rifampin", TIER_SUPPORTED_PRIMARY)
        known_ids = set(ep.genome_ids)
        doc = score_pretrained_cohort(
            adapter,
            [ep],
            split_scheme="lineage",
            time_provenance={"label": TIME_PROVENANCE_RETRO_SAME_DB},
            known_genome_ids=known_ids,
        )

        with tempfile.TemporaryDirectory() as tmp:
            out_path = Path(tmp) / "result.json"
            digest = write_result_atomic(out_path, doc, known_genome_ids=known_ids)
            self.assertEqual(len(digest), 64)
            self.assertEqual(check_zero_id_hygiene(out_path, known_genome_ids=known_ids)["genome_id_hits"], 0)

        # 1. Inject a known genome ID into a string field
        leaked_id = copy.deepcopy(doc)
        leaked_id["genotype_to_drug_rule"] += " saw 888000.0 here"
        with self.assertRaisesRegex(BenchError, "genome id hit"):
            check_zero_id_hygiene(leaked_id, known_genome_ids=known_ids)

        # 2. Inject a 25-bp nucleotide run
        leaked_seq = copy.deepcopy(doc)
        leaked_seq["genotype_to_drug_rule"] += " ACGTACGTACGTACGTACGTACGTA"
        with self.assertRaisesRegex(BenchError, "nucleotide_run"):
            check_zero_id_hygiene(leaked_seq)

        # 3. Inject a BioProject / BioSample / SRA accession
        for acc in ("PRJNA123456", "SAMN12345678", "GCA_000195955", "SRR1234567"):
            leaked_acc = copy.deepcopy(doc)
            leaked_acc["genotype_to_drug_rule"] += f" {acc}"
            with self.assertRaisesRegex(BenchError, "banned pattern"):
                check_zero_id_hygiene(leaked_acc)

        # 4. Inject forbidden timing or per-genome keys
        leaked_timing = copy.deepcopy(doc)
        leaked_timing["endpoints"][0]["wall_s"] = 1.23
        with self.assertRaisesRegex(BenchError, "forbidden"):
            check_zero_id_hygiene(leaked_timing)

        # 5. Inject unlisted top-level key
        extra_top = copy.deepcopy(doc)
        extra_top["unexpected_field"] = 1
        with self.assertRaisesRegex(BenchError, "unlisted top-level keys"):
            check_zero_id_hygiene(extra_top)

    def test_mykrobe_intermediate_calls_and_trainable_cohort_scoring(self):
        myk = load_adapter("mykrobe")
        ep_i = _make_endpoint(
            "mycobacterium_tuberculosis__isoniazid",
            "mycobacterium_tuberculosis",
            "isoniazid",
            TIER_SUPPORTED_PRIMARY,
            n_r=35,
            n_s=35,
            tp=26,
            tn=27,
            no_call_r=2,
            i_r=3,
            i_s=4,
            failure_s=1,
            id_offset=700,
        )
        ep_i = ProjectedEndpoint(**{**ep_i.__dict__, "reading": "tool_default"})
        doc_myk = score_pretrained_cohort(
            myk,
            [ep_i],
            split_scheme="lineage",
            time_provenance={"label": TIME_PROVENANCE_RETRO_SAME_DB},
            known_genome_ids=set(ep_i.genome_ids),
        )
        cnt = doc_myk["endpoints"][0]["counts"]
        self.assertEqual((cnt["i_r"], cnt["i_s"]), (3, 4))
        self.assertEqual((cnt["no_call_r"], cnt["failure_s"]), (2, 1))
        self.assertEqual(doc_myk["full_denominator_coverage"]["i_rows"], 7)
        # bounds_all_uncalled includes I; bounds_no_call_failure excludes I
        self.assertEqual(doc_myk["endpoints"][0]["bounds_no_call_failure"]["denom_r"], 32)
        self.assertEqual(doc_myk["endpoints"][0]["bounds_all_uncalled"]["denom_r"], 35)

        # Trainable cohort scoring across 2 species
        xgb = load_adapter("xgboost")
        with tempfile.TemporaryDirectory() as tmp:
            tdir = Path(tmp)
            d1, m1, _ = build_synthetic_export(tdir, species_slug="escherichia_coli", antibiotic="ciprofloxacin", n_samples=100, seed=42)
            d2, m2, _ = build_synthetic_export(tdir, species_slug="klebsiella_pneumoniae", antibiotic="meropenem", n_samples=100, seed=42)
            r1 = run_trainable_dataset(xgb, d1, "lineage", 42, fitter=_synthetic_linear_fitter)
            r2 = run_trainable_dataset(xgb, d2, "lineage", 42, fitter=_synthetic_linear_fitter)
            known = set(m1["samples"]) | set(m2["samples"])
            tdoc = score_trainable_cohort(
                xgb,
                [r1, r2],
                split_scheme="lineage",
                time_provenance={"label": TIME_PROVENANCE_RETRO_SAME_DB},
                known_genome_ids=known,
            )
            self.assertEqual(tdoc["comparison_kind"], "fixed_learner_representation")
            self.assertEqual(tdoc["summary_primary"]["eligible_endpoints"], 2)
            self.assertIn("log_loss", tdoc["endpoints"][0]["metrics"])
            self.assertEqual(check_zero_id_hygiene(tdoc, known_genome_ids=known)["genome_id_hits"], 0)


if __name__ == "__main__":
    unittest.main()
