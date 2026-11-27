"""Unified scoring core (`grm-bench-result-v1`) and zero-ID hygiene validator.

Implements (using Python stdlib + numpy only; zero scipy dependency):
- Exact Clopper-Pearson 95% CIs (Lentz continued-fraction regularized incomplete
  beta + 64-step bisection, accurate to < 1e-13).
- Wilson 95% CIs (matching scripts/c9-pick/score.py).
- Exact two-sided McNemar p-values on paired discordance counts
  (tool_only_right vs learner_only_right), overflow-safe for large m.
- Worst-case and best-case sensitivity, specificity, VME, ME, and balanced
  accuracy bounds under `no_call + failure` (and all uncalled categories).
- Pairwise common-support comparison on identical called isolates requiring
  >= 25 R and >= 25 S, alongside full-denominator coverage across all 4 tiers.
- Strict separation of support tiers (`supported_primary` only in primary
  medians; `supported_ungraded` and `proxy_class` reported separately;
  `unsupported` counted in coverage denominator only and NEVER mapped to S).
- 4,000-draw deterministic SHA-256 species-clustered bootstrap (multi-species)
  and lineage-clustered bootstrap (single-species tools: tbprofiler, mykrobe).
- Zero-ID hygiene validation reusing `scripts/c7-outside/scan.py` rules plus
  strict recursive schema allow-lists and banned sequence/accession patterns.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import math
import re
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np

from scripts.bench.pretrained import ProjectedEndpoint
from scripts.bench.protocol import (
    AdapterSpec,
    BenchError,
    EVALUATION_REVISION,
    RESULT_FORMAT,
    TIER_PROXY_CLASS,
    TIER_SUPPORTED_PRIMARY,
    TIER_SUPPORTED_UNGRADED,
    TIER_UNSUPPORTED,
    canonical_json_bytes,
    require,
    validate_time_provenance,
)

SCRIPTS_DIR = Path(__file__).resolve().parent.parent
SCAN_PY = SCRIPTS_DIR / "c7-outside" / "scan.py"

MIN_CLASS = 25
Z_95 = 1.959963984540054
BOOTSTRAP_RESAMPLES = 4000
BOOTSTRAP_LO_INDEX = 100   # 101st of 4,000 (0-based)
BOOTSTRAP_HI_INDEX = 3899  # 3,900th of 4,000 (0-based)
MIN_VALID_RESAMPLES = 40
MIN_SPECIES_UNWARNED = 6

# Banned patterns for zero-ID / zero-sequence hygiene (mirroring scan.py + tbprofiler/amrfinderplus)
BANNED_PATTERNS = (
    ("nucleotide_run", re.compile(r"(?i)[acgt]{20,}")),
    ("bioproject_accession", re.compile(r"\bPRJ[A-Z]{1,2}\d+")),
    ("biosample_accession", re.compile(r"\bSAM[NED][A-Z]?\d+")),
    ("assembly_accession", re.compile(r"\bGC[AF]_\d+")),
    ("sra_accession", re.compile(r"\b[SE]R[RSX]\d{5,}")),
)

# Forbidden keys anywhere inside grm-bench-result-v1
FORBIDDEN_KEYS = frozenset({
    "genome_id",
    "genome_ids",
    "samples",
    "sample_ids",
    "sid",
    "predictions",
    "scores",
    "raw_features",
    "features",
    "kmers",
    "weights",
    "wall_s",
    "cpu_s",
    "peak_rss_kb",
    "peak_kb",
    "oof_scratch_path",
    "fasta_path",
})

RESULT_TOP_KEYS = frozenset({
    "format",
    "evaluation_revision",
    "tool",
    "display_name",
    "mode",
    "comparison_kind",
    "reading",
    "split_scheme",
    "selection_role",
    "time_provenance",
    "bootstrap_cluster_unit",
    "bootstrap_resamples",
    "min_class",
    "software_pin",
    "database_pin",
    "i_handling",
    "genotype_to_drug_rule",
    "adapter_hashes",
    "full_denominator_coverage",
    "summary_primary",
    "summary_ungraded",
    "summary_proxy_class",
    "endpoints",
})


def _load_scan() -> Any:
    spec = importlib.util.spec_from_file_location("bench_c7_scan", SCAN_PY)
    if spec is None or spec.loader is None:
        raise BenchError(f"cannot load {SCAN_PY}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _r12(x: float | None) -> float | None:
    if x is None:
        return None
    val = round(float(x), 12)
    return 0.0 if val == 0.0 else val


# ------------------------------------------------------------------ exact Clopper-Pearson & Wilson


def wilson_ci(k: int, n: int, z: float = Z_95) -> tuple[float | None, float | None]:
    """Wilson score 95% CI for binomial proportion k / n."""
    if n <= 0:
        return None, None
    require(0 <= k <= n, f"wilson_ci: invalid k={k}, n={n}")
    p = k / n
    z2 = z * z
    d = 1.0 + z2 / n
    c = (p + z2 / (2.0 * n)) / d
    h = z * math.sqrt(p * (1.0 - p) / n + z2 / (4.0 * n * n)) / d
    return _r12(max(0.0, c - h)), _r12(min(1.0, c + h))


def _betacf(a: float, b: float, x: float) -> float:
    """Modified Lentz continued fraction for regularized incomplete beta I_x(a, b)."""
    max_iter = 400
    eps = 3.0e-15
    fpmin = 1.0e-300
    qab = a + b
    qap = a + 1.0
    qam = a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < fpmin:
        d = fpmin
    d = 1.0 / d
    h = d
    for m in range(1, max_iter + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < fpmin:
            d = fpmin
        c = 1.0 + aa / c
        if abs(c) < fpmin:
            c = fpmin
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < fpmin:
            d = fpmin
        c = 1.0 + aa / c
        if abs(c) < fpmin:
            c = fpmin
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < eps:
            break
    return h


def regularized_incomplete_beta(x: float, a: float, b: float) -> float:
    """Compute regularized incomplete beta function I_x(a, b) in [0, 1]."""
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    ln_bt = (
        math.lgamma(a + b)
        - math.lgamma(a)
        - math.lgamma(b)
        + a * math.log(x)
        + b * math.log(1.0 - x)
    )
    bt = math.exp(ln_bt) if ln_bt > -745.0 else 0.0
    if x < (a + 1.0) / (a + b + 2.0):
        return bt * _betacf(a, b, x) / a
    return 1.0 - bt * _betacf(b, a, 1.0 - x) / b


def _beta_ppf(q: float, a: float, b: float) -> float:
    """Inverse of regularized incomplete beta I_x(a, b) = q via 64-step bisection."""
    if q <= 0.0:
        return 0.0
    if q >= 1.0:
        return 1.0
    lo, hi = 0.0, 1.0
    for _ in range(64):
        mid = 0.5 * (lo + hi)
        if regularized_incomplete_beta(mid, a, b) < q:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def clopper_pearson_ci(k: int, n: int, alpha: float = 0.05) -> tuple[float | None, float | None]:
    """Exact Clopper-Pearson (1 - alpha) CI for binomial proportion k / n."""
    if n <= 0:
        return None, None
    require(0 <= k <= n, f"clopper_pearson_ci: invalid k={k}, n={n}")
    half = alpha / 2.0
    lo = 0.0 if k == 0 else _beta_ppf(half, float(k), float(n - k + 1))
    hi = 1.0 if k == n else _beta_ppf(1.0 - half, float(k + 1), float(n - k))
    return _r12(max(0.0, lo)), _r12(min(1.0, hi))


# ------------------------------------------------------------------ exact McNemar test


def mcnemar_exact_p(tool_only_right: int, learner_only_right: int) -> float:
    """Exact two-sided McNemar test p-value from discordant pair counts (b, c).

    Overflow-safe for arbitrary m = b + c (exact integer comb for m <= 1000,
    log-space binomial probabilities via lgamma for m > 1000).
    """
    require(tool_only_right >= 0 and learner_only_right >= 0, "McNemar counts must be non-negative")
    m = tool_only_right + learner_only_right
    if m == 0:
        return 1.0
    k = min(tool_only_right, learner_only_right)
    if m <= 1000:
        numer = sum(math.comb(m, i) for i in range(k + 1))
        p = math.ldexp(float(numer), 1 - m)
        return _r12(min(1.0, p)) or 0.0
    # Log-space sum for large m > 1000
    ln_half_m = -m * math.log(2.0)
    lg_m1 = math.lgamma(m + 1)
    log_terms = [
        lg_m1 - math.lgamma(i + 1) - math.lgamma(m - i + 1) + ln_half_m
        for i in range(k + 1)
    ]
    max_lt = max(log_terms)
    tail = math.exp(max_lt) * math.fsum(math.exp(lt - max_lt) for lt in log_terms)
    return _r12(min(1.0, 2.0 * tail)) or 0.0


# ------------------------------------------------------------------ rates and bounds


def rate_with_cis(k: int, n: int) -> dict[str, Any]:
    val = _r12(k / n) if n > 0 else None
    cp_lo, cp_hi = clopper_pearson_ci(k, n)
    w_lo, w_hi = wilson_ci(k, n)
    return {
        "count": int(k),
        "denom": int(n),
        "value": val,
        "clopper_pearson_lo": cp_lo,
        "clopper_pearson_hi": cp_hi,
        "wilson_lo": w_lo,
        "wilson_hi": w_hi,
    }


def bounds_under_uncalled(
    tp: int,
    fp: int,
    tn: int,
    fn: int,
    uncalled_r: int,
    uncalled_s: int,
) -> dict[str, float | None]:
    """Worst-case and best-case sensitivity, specificity, VME, ME, and balanced accuracy bounds."""
    denom_r = tp + fn + uncalled_r
    denom_s = tn + fp + uncalled_s
    sens_worst = _r12(tp / denom_r) if denom_r > 0 else None
    sens_best = _r12((tp + uncalled_r) / denom_r) if denom_r > 0 else None
    spec_worst = _r12(tn / denom_s) if denom_s > 0 else None
    spec_best = _r12((tn + uncalled_s) / denom_s) if denom_s > 0 else None
    vme_best = _r12(fn / denom_r) if denom_r > 0 else None
    vme_worst = _r12((fn + uncalled_r) / denom_r) if denom_r > 0 else None
    me_best = _r12(fp / denom_s) if denom_s > 0 else None
    me_worst = _r12((fp + uncalled_s) / denom_s) if denom_s > 0 else None
    ba_worst = _r12(0.5 * (sens_worst + spec_worst)) if (sens_worst is not None and spec_worst is not None) else None
    ba_best = _r12(0.5 * (sens_best + spec_best)) if (sens_best is not None and spec_best is not None) else None
    return {
        "denom_r": int(denom_r),
        "denom_s": int(denom_s),
        "sensitivity_worst": sens_worst,
        "sensitivity_best": sens_best,
        "specificity_worst": spec_worst,
        "specificity_best": spec_best,
        "vme_best": vme_best,
        "vme_worst": vme_worst,
        "me_best": me_best,
        "me_worst": me_worst,
        "balanced_accuracy_worst": ba_worst,
        "balanced_accuracy_best": ba_best,
    }


def compute_binary_metrics(tp: int, fp: int, tn: int, fn: int, loss_sum: float | None = None) -> dict[str, Any]:
    called_r = tp + fn
    called_s = tn + fp
    called_n = called_r + called_s
    sens = rate_with_cis(tp, called_r)
    spec = rate_with_cis(tn, called_s)
    vme = rate_with_cis(fn, called_r)
    me = rate_with_cis(fp, called_s)
    ba = (
        _r12(0.5 * (sens["value"] + spec["value"]))
        if (sens["value"] is not None and spec["value"] is not None)
        else None
    )
    out: dict[str, Any] = {
        "called_n": int(called_n),
        "called_r": int(called_r),
        "called_s": int(called_s),
        "tp": int(tp),
        "fp": int(fp),
        "tn": int(tn),
        "fn": int(fn),
        "sensitivity": sens,
        "specificity": spec,
        "balanced_accuracy": ba,
        "vme": vme,
        "me": me,
    }
    if loss_sum is not None and called_n > 0:
        out["log_loss"] = _r12(loss_sum / called_n)
    return out


# ------------------------------------------------------------------ deterministic SHA-256 bootstraps


def sha256_seed_int(seed_key: str) -> int:
    digest = hashlib.sha256(seed_key.encode("utf-8")).digest()
    return int.from_bytes(digest[:8], "little")


def species_cluster_bootstrap(
    values_by_species: dict[str, list[float]],
    seed_key: str,
    n_resamples: int = BOOTSTRAP_RESAMPLES,
) -> dict[str, Any]:
    """Species-clustered bootstrap (default 4,000 draws, deterministic SHA-256 seed).

    Resamples sorted species with replacement, pools all eligible endpoint values
    of each drawn species, and returns 95% percentile intervals for median and mean.
    """
    sp_keys = sorted(k for k, vals in values_by_species.items() if vals)
    n_species = len(sp_keys)
    all_vals = [float(x) for k in sp_keys for x in values_by_species[k]]
    if not all_vals:
        return {
            "cluster_unit": "species",
            "n_clusters": 0,
            "n_endpoints": 0,
            "resamples": n_resamples,
            "seed_sha256_prefix": hashlib.sha256(seed_key.encode("utf-8")).hexdigest()[:16],
            "median": None,
            "median_ci_lo": None,
            "median_ci_hi": None,
            "mean": None,
            "mean_ci_lo": None,
            "mean_ci_hi": None,
            "coarse_cluster_warning": True,
        }

    pt_med = _r12(float(np.median(all_vals)))
    pt_mean = _r12(float(np.mean(all_vals)))
    if n_species < 2:
        return {
            "cluster_unit": "species",
            "n_clusters": n_species,
            "n_endpoints": len(all_vals),
            "resamples": n_resamples,
            "seed_sha256_prefix": hashlib.sha256(seed_key.encode("utf-8")).hexdigest()[:16],
            "median": pt_med,
            "median_ci_lo": None,
            "median_ci_hi": None,
            "mean": pt_mean,
            "mean_ci_lo": None,
            "mean_ci_hi": None,
            "coarse_cluster_warning": True,
        }

    arrays = [np.asarray(values_by_species[k], dtype=np.float64) for k in sp_keys]
    rng = np.random.default_rng(sha256_seed_int(seed_key))
    meds = np.empty(n_resamples, dtype=np.float64)
    means = np.empty(n_resamples, dtype=np.float64)
    for b in range(n_resamples):
        drawn = rng.integers(0, n_species, size=n_species)
        sample = np.concatenate([arrays[int(i)] for i in drawn])
        meds[b] = float(np.median(sample))
        means[b] = float(np.mean(sample))
    meds.sort()
    means.sort()
    lo_idx = round(0.025 * (n_resamples - 1))
    hi_idx = n_resamples - 1 - lo_idx
    return {
        "cluster_unit": "species",
        "n_clusters": n_species,
        "n_endpoints": len(all_vals),
        "resamples": n_resamples,
        "seed_sha256_prefix": hashlib.sha256(seed_key.encode("utf-8")).hexdigest()[:16],
        "median": pt_med,
        "median_ci_lo": _r12(float(meds[lo_idx])),
        "median_ci_hi": _r12(float(meds[hi_idx])),
        "mean": pt_mean,
        "mean_ci_lo": _r12(float(means[lo_idx])),
        "mean_ci_hi": _r12(float(means[hi_idx])),
        "coarse_cluster_warning": n_species < MIN_SPECIES_UNWARNED,
    }


def join_lineage_groups(endpoints: list[ProjectedEndpoint]) -> tuple[dict[str, int], int]:
    """Union-find lineage clusters across endpoints, numbered by smallest genome_id.

    Mirrors scripts/tbprofiler-comparison/record.py::lineage_groups.
    Genomes with None lineage form singleton groups.
    """
    parent: dict[str, str] = {}

    def find(a: str) -> str:
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for ep in endpoints:
        first_in_lin: dict[tuple[str, str], str] = {}
        for gid, lin in zip(ep.genome_ids, ep.lineages):
            parent.setdefault(gid, gid)
            if lin is None:
                continue
            lkey = (ep.species_slug, str(lin))
            if lkey in first_in_lin:
                ra, rb = find(first_in_lin[lkey]), find(gid)
                if ra != rb:
                    parent[max(ra, rb)] = min(ra, rb)
            else:
                first_in_lin[lkey] = gid
    roots = sorted({find(g) for g in parent})
    idx = {r: i for i, r in enumerate(roots)}
    return {g: idx[find(g)] for g in parent}, len(roots)


def _interval_from_valid(values: list[float | None]) -> tuple[float | None, float | None, int]:
    valid = sorted(x for x in values if x is not None)
    nv = len(valid)
    if nv < MIN_VALID_RESAMPLES:
        return None, None, nv
    lo_idx = round(0.025 * (nv - 1))
    hi_idx = nv - 1 - lo_idx
    return _r12(valid[lo_idx]), _r12(valid[hi_idx]), nv


def lineage_cluster_bootstrap(
    endpoints: list[ProjectedEndpoint],
    seed_key: str,
    n_resamples: int = BOOTSTRAP_RESAMPLES,
) -> dict[str, Any]:
    """Lineage-cluster bootstrap (4,000 draws, deterministic SHA-256 seed) for single-species tools.

    Resamples joined lineage clusters with replacement across all passed endpoints
    and computes per-endpoint and across-endpoint median 95% CIs for sensitivity,
    specificity, balanced accuracy, VME, ME, and comparator BA differences.
    """
    if not endpoints:
        return {
            "cluster_unit": "lineage",
            "n_clusters": 0,
            "n_endpoints": 0,
            "resamples": n_resamples,
            "seed_sha256_prefix": hashlib.sha256(seed_key.encode("utf-8")).hexdigest()[:16],
            "per_endpoint": {},
            "median_metrics": {},
            "median_comparator_diffs": {},
        }

    group_of, n_groups = join_lineage_groups(endpoints)
    if n_groups < 2:
        return {
            "cluster_unit": "lineage",
            "n_clusters": n_groups,
            "n_endpoints": len(endpoints),
            "resamples": n_resamples,
            "seed_sha256_prefix": hashlib.sha256(seed_key.encode("utf-8")).hexdigest()[:16],
            "per_endpoint": {},
            "median_metrics": {},
            "median_comparator_diffs": {},
        }

    rng = np.random.default_rng(sha256_seed_int(seed_key))
    # Pre-draw (n_resamples, n_groups) multinomial weight matrix W where each row sums to n_groups.
    draws = rng.integers(0, n_groups, size=(n_resamples, n_groups))
    weights = np.zeros((n_resamples, n_groups), dtype=np.float64)
    row_ids = np.repeat(np.arange(n_resamples), n_groups)
    np.add.at(weights, (row_ids, draws.ravel()), 1.0)

    comp_names = sorted({
        cname
        for ep in endpoints
        if ep.comparator_calls
        for cname in ep.comparator_calls
    })

    # Per-endpoint draw vectors: shape (n_resamples,) with NaN when class absent in draw
    ep_metric_draws: dict[str, dict[str, np.ndarray]] = {}
    ep_diff_draws: dict[tuple[str, str], np.ndarray] = {}
    per_endpoint_ci: dict[str, dict[str, Any]] = {}

    for ep in endpoints:
        # Pack tool (tp, fp, tn, fn) per lineage group over called isolates
        tool_mat = np.zeros((n_groups, 4), dtype=np.float64)
        for gid, y, c in zip(ep.genome_ids, ep.truth, ep.calls):
            if c not in ("R", "S"):
                continue
            g = group_of[gid]
            if y and c == "R":
                tool_mat[g, 0] += 1.0
            elif (not y) and c == "R":
                tool_mat[g, 1] += 1.0
            elif (not y) and c == "S":
                tool_mat[g, 2] += 1.0
            else:
                tool_mat[g, 3] += 1.0

        agg = weights @ tool_mat  # (n_resamples, 4)
        tp, fp, tn, fn = agg[:, 0], agg[:, 1], agg[:, 2], agg[:, 3]
        pos = tp + fn
        neg = tn + fp
        valid_both = (pos > 0) & (neg > 0)

        sens = np.where(pos > 0, tp / np.where(pos > 0, pos, 1.0), np.nan)
        spec = np.where(neg > 0, tn / np.where(neg > 0, neg, 1.0), np.nan)
        ba = np.where(valid_both, 0.5 * (sens + spec), np.nan)
        vme = np.where(pos > 0, fn / np.where(pos > 0, pos, 1.0), np.nan)
        me = np.where(neg > 0, fp / np.where(neg > 0, neg, 1.0), np.nan)

        mdict = {
            "sensitivity": sens,
            "specificity": spec,
            "balanced_accuracy": ba,
            "vme": vme,
            "me": me,
        }
        ep_metric_draws[ep.dataset] = mdict
        ep_out: dict[str, Any] = {}
        for mname, arr in mdict.items():
            vals = [None if math.isnan(float(v)) else float(v) for v in arr]
            lo, hi, nv = _interval_from_valid(vals)
            ep_out[mname] = {"ci_lo": lo, "ci_hi": hi, "valid_resamples": nv}

        # Paired comparator common-support bootstrap
        comp_out: dict[str, Any] = {}
        if ep.comparator_calls:
            for cname, ccalls in sorted(ep.comparator_calls.items()):
                t_pair = np.zeros((n_groups, 4), dtype=np.float64)
                l_pair = np.zeros((n_groups, 4), dtype=np.float64)
                for gid, y, tc, lc in zip(ep.genome_ids, ep.truth, ep.calls, ccalls):
                    if tc not in ("R", "S") or lc not in ("R", "S"):
                        continue
                    g = group_of[gid]
                    t_pair[g, 0 if (y and tc == "R") else 1 if ((not y) and tc == "R") else 2 if ((not y) and tc == "S") else 3] += 1.0
                    l_pair[g, 0 if (y and lc == "R") else 1 if ((not y) and lc == "R") else 2 if ((not y) and lc == "S") else 3] += 1.0
                ta = weights @ t_pair
                la = weights @ l_pair
                p_pos = ta[:, 0] + ta[:, 3]
                p_neg = ta[:, 2] + ta[:, 1]
                ok = (p_pos > 0) & (p_neg > 0)
                t_ba = np.where(ok, 0.5 * (ta[:, 0] / np.where(ok, p_pos, 1.0) + ta[:, 2] / np.where(ok, p_neg, 1.0)), np.nan)
                l_ba = np.where(ok, 0.5 * (la[:, 0] / np.where(ok, p_pos, 1.0) + la[:, 2] / np.where(ok, p_neg, 1.0)), np.nan)
                diff = t_ba - l_ba
                ep_diff_draws[(ep.dataset, cname)] = diff
                dvals = [None if math.isnan(float(v)) else float(v) for v in diff]
                lo, hi, nv = _interval_from_valid(dvals)
                comp_out[cname] = {"ba_diff_ci_lo": lo, "ba_diff_ci_hi": hi, "valid_resamples": nv}
        ep_out["comparator_diffs"] = comp_out
        per_endpoint_ci[ep.dataset] = ep_out

    # Across-endpoint medians per resample
    median_metrics: dict[str, Any] = {}
    for mname in ("sensitivity", "specificity", "balanced_accuracy", "vme", "me"):
        mat = np.vstack([ep_metric_draws[ep.dataset][mname] for ep in endpoints])  # (E, B)
        med_vals: list[float | None] = []
        for b in range(n_resamples):
            col = mat[:, b]
            finite = col[~np.isnan(col)]
            med_vals.append(float(np.median(finite)) if finite.size else None)
        lo, hi, nv = _interval_from_valid(med_vals)
        median_metrics[mname] = {"median_ci_lo": lo, "median_ci_hi": hi, "valid_resamples": nv}

    median_comp_diffs: dict[str, Any] = {}
    for cname in comp_names:
        cols = [ep_diff_draws[(ep.dataset, cname)] for ep in endpoints if (ep.dataset, cname) in ep_diff_draws]
        if not cols:
            continue
        mat = np.vstack(cols)
        med_vals = []
        for b in range(n_resamples):
            col = mat[:, b]
            finite = col[~np.isnan(col)]
            med_vals.append(float(np.median(finite)) if finite.size else None)
        lo, hi, nv = _interval_from_valid(med_vals)
        median_comp_diffs[cname] = {"median_ba_diff_ci_lo": lo, "median_ba_diff_ci_hi": hi, "valid_resamples": nv}

    return {
        "cluster_unit": "lineage",
        "n_clusters": n_groups,
        "n_endpoints": len(endpoints),
        "resamples": n_resamples,
        "seed_sha256_prefix": hashlib.sha256(seed_key.encode("utf-8")).hexdigest()[:16],
        "per_endpoint": per_endpoint_ci,
        "median_metrics": median_metrics,
        "median_comparator_diffs": median_comp_diffs,
    }


# ------------------------------------------------------------------ endpoint & cohort scoring


def _score_projected_endpoint(ep: ProjectedEndpoint, min_class: int = MIN_CLASS) -> dict[str, Any]:
    tp = fp = tn = fn = 0
    no_call_r = no_call_s = 0
    i_r = i_s = 0
    failure_r = failure_s = 0
    unsupported_r = unsupported_s = 0

    for y, c in zip(ep.truth, ep.calls):
        if c == "R":
            tp += int(y)
            fp += int(not y)
        elif c == "S":
            fn += int(y)
            tn += int(not y)
        elif c == "no_call":
            no_call_r += int(y)
            no_call_s += int(not y)
        elif c == "I":
            i_r += int(y)
            i_s += int(not y)
        elif c == "failure":
            failure_r += int(y)
            failure_s += int(not y)
        elif c == "unsupported":
            unsupported_r += int(y)
            unsupported_s += int(not y)
        else:
            raise BenchError(f"{ep.dataset}: invalid projected call {c!r}")

    total_r = tp + fn + no_call_r + i_r + failure_r + unsupported_r
    total_s = tn + fp + no_call_s + i_s + failure_s + unsupported_s
    total_n = total_r + total_s
    called_r = tp + fn
    called_s = tn + fp
    called_n = called_r + called_s

    eligible = (
        ep.support_tier != TIER_UNSUPPORTED
        and called_r >= min_class
        and called_s >= min_class
    )
    if ep.support_tier == TIER_UNSUPPORTED:
        ineligible_reason: str | None = f"unsupported:{ep.unsupported_reason}"
    elif not eligible:
        ineligible_reason = f"below_min_class_{min_class}"
    else:
        ineligible_reason = None

    # Pairwise common-support comparisons on identical called isolates
    pairwise: dict[str, Any] = {}
    if ep.support_tier != TIER_UNSUPPORTED and ep.comparator_calls:
        for cname, ccalls in sorted(ep.comparator_calls.items()):
            p_t_tp = p_t_fp = p_t_tn = p_t_fn = 0
            p_l_tp = p_l_fp = p_l_tn = p_l_fn = 0
            tool_only = learner_only = both_right = both_wrong = 0
            for y, tc, lc in zip(ep.truth, ep.calls, ccalls):
                if tc not in ("R", "S") or lc not in ("R", "S"):
                    continue
                t_r = tc == "R"
                l_r = lc == "R"
                if y:
                    p_t_tp += int(t_r)
                    p_t_fn += int(not t_r)
                    p_l_tp += int(l_r)
                    p_l_fn += int(not l_r)
                else:
                    p_t_fp += int(t_r)
                    p_t_tn += int(not t_r)
                    p_l_fp += int(l_r)
                    p_l_tn += int(not l_r)
                t_ok = t_r == y
                l_ok = l_r == y
                if t_ok and not l_ok:
                    tool_only += 1
                elif l_ok and not t_ok:
                    learner_only += 1
                elif t_ok and l_ok:
                    both_right += 1
                else:
                    both_wrong += 1

            p_r = p_t_tp + p_t_fn
            p_s = p_t_tn + p_t_fp
            p_n = p_r + p_s
            p_elig = p_r >= min_class and p_s >= min_class
            t_m = compute_binary_metrics(p_t_tp, p_t_fp, p_t_tn, p_t_fn)
            l_m = compute_binary_metrics(p_l_tp, p_l_fp, p_l_tn, p_l_fn)

            def _diff(a: float | None, b: float | None) -> float | None:
                return None if (a is None or b is None) else _r12(a - b)

            pairwise[cname] = {
                "eligible_common_support": p_elig,
                "common_called_n": p_n,
                "common_called_r": p_r,
                "common_called_s": p_s,
                "tool_metrics": t_m,
                "learner_metrics": l_m,
                "diff_tool_minus_learner": {
                    "balanced_accuracy": _diff(t_m["balanced_accuracy"], l_m["balanced_accuracy"]),
                    "sensitivity": _diff(t_m["sensitivity"]["value"], l_m["sensitivity"]["value"]),
                    "specificity": _diff(t_m["specificity"]["value"], l_m["specificity"]["value"]),
                    "vme": _diff(t_m["vme"]["value"], l_m["vme"]["value"]),
                    "me": _diff(t_m["me"]["value"], l_m["me"]["value"]),
                },
                "discordance": {
                    "tool_only_right": tool_only,
                    "learner_only_right": learner_only,
                    "both_right": both_right,
                    "both_wrong": both_wrong,
                    "calls_agree": both_right + both_wrong,
                    "mcnemar_p": mcnemar_exact_p(tool_only, learner_only),
                },
            }

    return {
        "dataset": ep.dataset,
        "species_slug": ep.species_slug,
        "canonical_drug": ep.canonical_drug,
        "support_tier": ep.support_tier,
        "unsupported_reason": ep.unsupported_reason,
        "eligible_common_support": eligible,
        "ineligible_reason": ineligible_reason,
        "fold_file_sha256": ep.fold_file_sha256,
        "counts": {
            "total_n": total_n,
            "total_r": total_r,
            "total_s": total_s,
            "called_n": called_n,
            "called_r": called_r,
            "called_s": called_s,
            "tp": tp,
            "fp": fp,
            "tn": tn,
            "fn": fn,
            "no_call_r": no_call_r,
            "no_call_s": no_call_s,
            "i_r": i_r,
            "i_s": i_s,
            "failure_r": failure_r,
            "failure_s": failure_s,
            "unsupported_r": unsupported_r,
            "unsupported_s": unsupported_s,
        },
        "coverage": {
            "all": _r12(called_n / total_n) if total_n > 0 else None,
            "resistant": _r12(called_r / total_r) if total_r > 0 else None,
            "susceptible": _r12(called_s / total_s) if total_s > 0 else None,
        },
        "metrics": compute_binary_metrics(tp, fp, tn, fn),
        "bounds_no_call_failure": bounds_under_uncalled(
            tp, fp, tn, fn, no_call_r + failure_r, no_call_s + failure_s
        ),
        "bounds_all_uncalled": bounds_under_uncalled(
            tp,
            fp,
            tn,
            fn,
            no_call_r + failure_r + i_r + unsupported_r,
            no_call_s + failure_s + i_s + unsupported_s,
        ),
        "pairwise_common_support": pairwise,
    }


def _summarize_tier(
    tier_name: str,
    endpoint_rows: list[dict[str, Any]],
    projected_map: dict[str, ProjectedEndpoint],
    cluster_unit: str,
    seed_prefix: str,
    n_resamples: int,
) -> dict[str, Any]:
    """Summarize one support tier (never pooling other tiers into it)."""
    eligible_rows = [r for r in endpoint_rows if r["eligible_common_support"]]
    eligible_eps = [projected_map[r["dataset"]] for r in eligible_rows if r["dataset"] in projected_map]

    metrics_summary: dict[str, Any] = {}
    for mname in ("balanced_accuracy", "sensitivity", "specificity", "vme", "me"):
        by_sp: dict[str, list[float]] = {}
        for r in eligible_rows:
            val = r["metrics"]["balanced_accuracy"] if mname == "balanced_accuracy" else r["metrics"][mname]["value"]
            if val is not None:
                by_sp.setdefault(r["species_slug"], []).append(float(val))
        if cluster_unit == "species":
            metrics_summary[mname] = species_cluster_bootstrap(
                by_sp, f"{seed_prefix}|{tier_name}|metric|{mname}", n_resamples=n_resamples
            )
        else:
            vals = [x for xs in by_sp.values() for x in xs]
            metrics_summary[mname] = {
                "cluster_unit": "lineage",
                "n_endpoints": len(vals),
                "median": _r12(float(np.median(vals))) if vals else None,
                "mean": _r12(float(np.mean(vals))) if vals else None,
            }

    comp_names = sorted({
        cname
        for r in eligible_rows
        for cname in r["pairwise_common_support"]
    })
    pairwise_summary: dict[str, Any] = {}
    for cname in comp_names:
        by_sp_diff: dict[str, list[float]] = {}
        tool_ahead = learner_ahead = equal = 0
        tot_t_only = tot_l_only = 0
        for r in eligible_rows:
            pw = r["pairwise_common_support"].get(cname)
            if not pw or not pw["eligible_common_support"]:
                continue
            d = pw["diff_tool_minus_learner"]["balanced_accuracy"]
            if d is None:
                continue
            by_sp_diff.setdefault(r["species_slug"], []).append(float(d))
            if d > 0:
                tool_ahead += 1
            elif d < 0:
                learner_ahead += 1
            else:
                equal += 1
            tot_t_only += pw["discordance"]["tool_only_right"]
            tot_l_only += pw["discordance"]["learner_only_right"]
        if cluster_unit == "species":
            boot = species_cluster_bootstrap(
                by_sp_diff, f"{seed_prefix}|{tier_name}|pairwise|{cname}|ba", n_resamples=n_resamples
            )
        else:
            vals = [x for xs in by_sp_diff.values() for x in xs]
            boot = {
                "cluster_unit": "lineage",
                "n_endpoints": len(vals),
                "median": _r12(float(np.median(vals))) if vals else None,
                "mean": _r12(float(np.mean(vals))) if vals else None,
            }
        pairwise_summary[cname] = {
            "tool_ahead": tool_ahead,
            "learner_ahead": learner_ahead,
            "equal": equal,
            "pooled_tool_only_right": tot_t_only,
            "pooled_learner_only_right": tot_l_only,
            "pooled_mcnemar_p": mcnemar_exact_p(tot_t_only, tot_l_only),
            "ba_diff_tool_minus_learner": boot,
        }

    lineage_boot: dict[str, Any] | None = None
    if cluster_unit == "lineage" and eligible_eps:
        lineage_boot = lineage_cluster_bootstrap(
            eligible_eps, f"{seed_prefix}|{tier_name}|lineage", n_resamples=n_resamples
        )
        # Attach lineage bootstrap CIs into metrics_summary and pairwise_summary
        for mname, ci in lineage_boot.get("median_metrics", {}).items():
            if mname in metrics_summary:
                metrics_summary[mname].update(ci)
        for cname, ci in lineage_boot.get("median_comparator_diffs", {}).items():
            if cname in pairwise_summary:
                pairwise_summary[cname]["ba_diff_tool_minus_learner"].update(ci)

    return {
        "support_tier": tier_name,
        "total_endpoints": len(endpoint_rows),
        "eligible_endpoints": len(eligible_rows),
        "species_count": len({r["species_slug"] for r in eligible_rows}),
        "metrics": metrics_summary,
        "pairwise_common_support": pairwise_summary,
        "lineage_bootstrap": lineage_boot,
    }


def score_pretrained_cohort(
    adapter: AdapterSpec,
    endpoints: list[ProjectedEndpoint],
    *,
    split_scheme: str,
    time_provenance: dict[str, Any],
    min_class: int = MIN_CLASS,
    n_resamples: int = BOOTSTRAP_RESAMPLES,
    known_genome_ids: set[str] | None = None,
) -> dict[str, Any]:
    """Build a deterministic, zero-ID `grm-bench-result-v1` document for a pretrained tool."""
    require(len(endpoints) >= 1, "score_pretrained_cohort requires at least one endpoint")
    sorted_eps = sorted(endpoints, key=lambda e: e.dataset)
    require(len({e.dataset for e in sorted_eps}) == len(sorted_eps), "duplicate dataset in endpoints")
    roles = {e.selection_role for e in sorted_eps}
    require(len(roles) == 1, f"cannot mix selection_roles in one result: {sorted(roles)}")
    selection_role = next(iter(roles))
    readings = {e.reading for e in sorted_eps}
    require(len(readings) == 1, f"cannot mix readings in one result: {sorted(readings)}")
    reading = next(iter(readings))

    # Validate time provenance
    prov = validate_time_provenance(
        split_scheme=split_scheme,
        declared_label=time_provenance["label"],
        fold_cutoff_date=time_provenance.get("fold_cutoff_date"),
        software_release_date=time_provenance.get("software_release_date", adapter.manifest["software_pin"].get("release_date")),
        database_release_date=time_provenance.get("database_release_date", (adapter.manifest["database_pin"] or {}).get("release_date")),
        weights_release_date=time_provenance.get("weights_release_date"),
        training_labels_max_date=time_provenance.get("training_labels_max_date"),
        mode=adapter.mode,
        is_bvbrc_insertion_time_split=bool(time_provenance.get("is_bvbrc_insertion_time_split", False)),
    )

    projected_map = {e.dataset: e for e in sorted_eps}
    ep_rows = [_score_projected_endpoint(e, min_class=min_class) for e in sorted_eps]

    # Full-denominator coverage across ALL tiers
    tier_counts = Counter(r["support_tier"] for r in ep_rows)
    unsup_reasons = Counter(r["unsupported_reason"] for r in ep_rows if r["support_tier"] == TIER_UNSUPPORTED)
    tot_rows = sum(r["counts"]["total_n"] for r in ep_rows)
    called_rows = sum(r["counts"]["called_n"] for r in ep_rows)
    nc_rows = sum(r["counts"]["no_call_r"] + r["counts"]["no_call_s"] for r in ep_rows)
    i_rows = sum(r["counts"]["i_r"] + r["counts"]["i_s"] for r in ep_rows)
    fail_rows = sum(r["counts"]["failure_r"] + r["counts"]["failure_s"] for r in ep_rows)
    unsup_rows = sum(r["counts"]["unsupported_r"] + r["counts"]["unsupported_s"] for r in ep_rows)

    primary_rows = [r for r in ep_rows if r["support_tier"] == TIER_SUPPORTED_PRIMARY]
    ungraded_rows = [r for r in ep_rows if r["support_tier"] == TIER_SUPPORTED_UNGRADED]
    proxy_rows = [r for r in ep_rows if r["support_tier"] == TIER_PROXY_CLASS]

    prim_tot = sum(r["counts"]["total_n"] for r in primary_rows)
    prim_called = sum(r["counts"]["called_n"] for r in primary_rows)

    coverage = {
        "total_datasets": len(ep_rows),
        "supported_primary_datasets": tier_counts[TIER_SUPPORTED_PRIMARY],
        "eligible_primary_datasets": sum(1 for r in primary_rows if r["eligible_common_support"]),
        "supported_ungraded_datasets": tier_counts[TIER_SUPPORTED_UNGRADED],
        "proxy_class_datasets": tier_counts[TIER_PROXY_CLASS],
        "unsupported_datasets": tier_counts[TIER_UNSUPPORTED],
        "unsupported_reasons": dict(sorted(unsup_reasons.items())),
        "total_genome_drug_rows": tot_rows,
        "called_rows": called_rows,
        "no_call_rows": nc_rows,
        "i_rows": i_rows,
        "failure_rows": fail_rows,
        "unsupported_rows": unsup_rows,
        "full_cohort_call_rate": _r12(called_rows / tot_rows) if tot_rows > 0 else None,
        "supported_primary_call_rate": _r12(prim_called / prim_tot) if prim_tot > 0 else None,
    }

    seed_prefix = f"grm-bench-v1|{adapter.tool}|{reading}|{split_scheme}"
    cluster_unit = adapter.bootstrap_cluster_unit

    doc = {
        "format": RESULT_FORMAT,
        "evaluation_revision": EVALUATION_REVISION,
        "tool": adapter.tool,
        "display_name": adapter.manifest["display_name"],
        "mode": adapter.mode,
        "comparison_kind": adapter.comparison_kind,
        "reading": reading,
        "split_scheme": split_scheme,
        "selection_role": selection_role,
        "time_provenance": prov,
        "bootstrap_cluster_unit": cluster_unit,
        "bootstrap_resamples": n_resamples,
        "min_class": min_class,
        "software_pin": adapter.manifest["software_pin"],
        "database_pin": adapter.manifest["database_pin"],
        "i_handling": adapter.manifest["i_handling"],
        "genotype_to_drug_rule": adapter.manifest["genotype_to_drug_rule"],
        "adapter_hashes": {
            "manifest_sha256": adapter.manifest_sha256,
            "drug_map_sha256": adapter.drug_map_sha256,
            "species_support_sha256": adapter.species_support_sha256,
        },
        "full_denominator_coverage": coverage,
        "summary_primary": _summarize_tier(
            TIER_SUPPORTED_PRIMARY, primary_rows, projected_map, cluster_unit, seed_prefix, n_resamples
        ),
        "summary_ungraded": _summarize_tier(
            TIER_SUPPORTED_UNGRADED, ungraded_rows, projected_map, cluster_unit, seed_prefix, n_resamples
        ),
        "summary_proxy_class": _summarize_tier(
            TIER_PROXY_CLASS, proxy_rows, projected_map, cluster_unit, seed_prefix, n_resamples
        ),
        "endpoints": ep_rows,
    }

    check_zero_id_hygiene(doc, known_genome_ids=known_genome_ids)
    return doc


def score_trainable_cohort(
    adapter: AdapterSpec,
    dataset_runs: list[dict[str, Any]],
    *,
    split_scheme: str,
    time_provenance: dict[str, Any],
    min_class: int = MIN_CLASS,
    n_resamples: int = BOOTSTRAP_RESAMPLES,
    known_genome_ids: set[str] | None = None,
) -> dict[str, Any]:
    """Aggregate trainable per-fold runs into a deterministic, zero-ID `grm-bench-result-v1`."""
    require(len(dataset_runs) >= 1, "score_trainable_cohort requires at least one dataset run")
    sorted_runs = sorted(dataset_runs, key=lambda r: r["dataset"])
    require(len({r["dataset"] for r in sorted_runs}) == len(sorted_runs), "duplicate dataset in trainable runs")
    roles = {r["selection_role"] for r in sorted_runs}
    require(len(roles) == 1, f"cannot mix selection_roles: {sorted(roles)}")
    selection_role = next(iter(roles))

    prov = validate_time_provenance(
        split_scheme=split_scheme,
        declared_label=time_provenance["label"],
        fold_cutoff_date=time_provenance.get("fold_cutoff_date"),
        software_release_date=time_provenance.get("software_release_date", adapter.manifest["software_pin"].get("release_date")),
        database_release_date=time_provenance.get("database_release_date"),
        weights_release_date=time_provenance.get("weights_release_date"),
        training_labels_max_date=time_provenance.get("training_labels_max_date"),
        mode=adapter.mode,
        is_bvbrc_insertion_time_split=bool(time_provenance.get("is_bvbrc_insertion_time_split", False)),
    )

    ep_rows: list[dict[str, Any]] = []
    for run in sorted_runs:
        folds = run.get("folds") or []
        tp = sum(int(f["tp"]) for f in folds)
        fp = sum(int(f["fp"]) for f in folds)
        tn = sum(int(f["tn"]) for f in folds)
        fn = sum(int(f["fn"]) for f in folds)
        nc_r = sum(int(f["no_call_r"]) for f in folds)
        nc_s = sum(int(f["no_call_s"]) for f in folds)
        i_r = sum(int(f.get("i_r", 0)) for f in folds)
        i_s = sum(int(f.get("i_s", 0)) for f in folds)
        fail_r = sum(int(f["failure_r"]) for f in folds)
        fail_s = sum(int(f["failure_s"]) for f in folds)
        loss_vals = [f["loss_sum"] for f in folds if "loss_sum" in f and f["loss_sum"] is not None]
        loss_sum = float(math.fsum(loss_vals)) if (loss_vals and len(loss_vals) == len(folds)) else None

        called_r = tp + fn
        called_s = tn + fp
        called_n = called_r + called_s
        total_r = called_r + nc_r + i_r + fail_r
        total_s = called_s + nc_s + i_s + fail_s
        total_n = total_r + total_s

        tier = run["support_tier"]
        unavail = run.get("unavailable_reason")
        eligible = unavail is None and tier != TIER_UNSUPPORTED and called_r >= min_class and called_s >= min_class
        inelig_reason = (
            f"unavailable:{unavail}"
            if unavail is not None
            else (f"below_min_class_{min_class}" if not eligible else None)
        )

        ep_rows.append({
            "dataset": run["dataset"],
            "species_slug": run["species_slug"],
            "canonical_drug": run["canonical_drug"],
            "support_tier": tier,
            "unsupported_reason": run.get("unsupported_reason", "-"),
            "eligible_common_support": eligible,
            "ineligible_reason": inelig_reason,
            "input_fingerprint": run.get("input_fingerprint"),
            "matrix_sha256": run.get("matrix_sha256"),
            "plan_sha256": run.get("plan_sha256"),
            "fold_file_sha256": run.get("fold_file_sha256"),
            "folds": folds,
            "counts": {
                "total_n": total_n,
                "total_r": total_r,
                "total_s": total_s,
                "called_n": called_n,
                "called_r": called_r,
                "called_s": called_s,
                "tp": tp,
                "fp": fp,
                "tn": tn,
                "fn": fn,
                "no_call_r": nc_r,
                "no_call_s": nc_s,
                "i_r": i_r,
                "i_s": i_s,
                "failure_r": fail_r,
                "failure_s": fail_s,
                "unsupported_r": 0,
                "unsupported_s": 0,
            },
            "coverage": {
                "all": _r12(called_n / total_n) if total_n > 0 else None,
                "resistant": _r12(called_r / total_r) if total_r > 0 else None,
                "susceptible": _r12(called_s / total_s) if total_s > 0 else None,
            },
            "metrics": compute_binary_metrics(tp, fp, tn, fn, loss_sum=loss_sum),
            "bounds_no_call_failure": bounds_under_uncalled(tp, fp, tn, fn, nc_r + fail_r, nc_s + fail_s),
            "bounds_all_uncalled": bounds_under_uncalled(tp, fp, tn, fn, nc_r + fail_r + i_r, nc_s + fail_s + i_s),
            "pairwise_common_support": {},
        })

    tot_rows = sum(r["counts"]["total_n"] for r in ep_rows)
    called_rows = sum(r["counts"]["called_n"] for r in ep_rows)
    primary_rows = [r for r in ep_rows if r["support_tier"] == TIER_SUPPORTED_PRIMARY]
    coverage = {
        "total_datasets": len(ep_rows),
        "supported_primary_datasets": len(primary_rows),
        "eligible_primary_datasets": sum(1 for r in primary_rows if r["eligible_common_support"]),
        "supported_ungraded_datasets": 0,
        "proxy_class_datasets": 0,
        "unsupported_datasets": 0,
        "unsupported_reasons": {},
        "total_genome_drug_rows": tot_rows,
        "called_rows": called_rows,
        "no_call_rows": sum(r["counts"]["no_call_r"] + r["counts"]["no_call_s"] for r in ep_rows),
        "i_rows": sum(r["counts"]["i_r"] + r["counts"]["i_s"] for r in ep_rows),
        "failure_rows": sum(r["counts"]["failure_r"] + r["counts"]["failure_s"] for r in ep_rows),
        "unsupported_rows": 0,
        "full_cohort_call_rate": _r12(called_rows / tot_rows) if tot_rows > 0 else None,
        "supported_primary_call_rate": _r12(called_rows / tot_rows) if tot_rows > 0 else None,
    }

    seed_prefix = f"grm-bench-v1|{adapter.tool}|{adapter.primary_reading}|{split_scheme}"
    doc = {
        "format": RESULT_FORMAT,
        "evaluation_revision": EVALUATION_REVISION,
        "tool": adapter.tool,
        "display_name": adapter.manifest["display_name"],
        "mode": adapter.mode,
        "comparison_kind": adapter.comparison_kind,
        "reading": adapter.primary_reading,
        "split_scheme": split_scheme,
        "selection_role": selection_role,
        "time_provenance": prov,
        "bootstrap_cluster_unit": adapter.bootstrap_cluster_unit,
        "bootstrap_resamples": n_resamples,
        "min_class": min_class,
        "software_pin": adapter.manifest["software_pin"],
        "database_pin": adapter.manifest["database_pin"],
        "i_handling": adapter.manifest["i_handling"],
        "genotype_to_drug_rule": adapter.manifest["genotype_to_drug_rule"],
        "adapter_hashes": {
            "manifest_sha256": adapter.manifest_sha256,
            "drug_map_sha256": adapter.drug_map_sha256,
            "species_support_sha256": adapter.species_support_sha256,
        },
        "full_denominator_coverage": coverage,
        "summary_primary": _summarize_tier(
            TIER_SUPPORTED_PRIMARY, primary_rows, {}, adapter.bootstrap_cluster_unit, seed_prefix, n_resamples
        ),
        "summary_ungraded": _summarize_tier(
            TIER_SUPPORTED_UNGRADED, [], {}, adapter.bootstrap_cluster_unit, seed_prefix, n_resamples
        ),
        "summary_proxy_class": _summarize_tier(
            TIER_PROXY_CLASS, [], {}, adapter.bootstrap_cluster_unit, seed_prefix, n_resamples
        ),
        "endpoints": ep_rows,
    }
    check_zero_id_hygiene(doc, known_genome_ids=known_genome_ids)
    return doc


# ------------------------------------------------------------------ zero-ID hygiene (scan.py rules)


def _walk_keys(obj: Any, path: str = "$") -> list[str]:
    bad: list[str] = []
    if isinstance(obj, dict):
        for k, v in obj.items():
            if k in FORBIDDEN_KEYS:
                bad.append(f"{path}.{k}")
            bad.extend(_walk_keys(v, f"{path}.{k}"))
    elif isinstance(obj, list):
        for idx, item in enumerate(obj):
            bad.extend(_walk_keys(item, f"{path}[{idx}]"))
    return bad


def check_zero_id_hygiene(
    doc_or_path: dict[str, Any] | Path | str,
    *,
    known_genome_ids: set[str] | None = None,
    species_json_paths: list[Path | str] | None = None,
) -> dict[str, int]:
    """Validate `grm-bench-result-v1` schema keys and zero-ID hygiene using scripts/c7-outside/scan.py.

    Refuses any genome ID hit, banned sequence/accession shape, unlisted top-level
    key, or forbidden nested key (e.g., genome_id, predictions, scores, wall_s, peak_rss_kb).
    Never prints matched genome IDs.
    """
    if isinstance(doc_or_path, (Path, str)):
        p = Path(doc_or_path)
        raw_text = p.read_text(encoding="utf-8")
        doc = json.loads(raw_text)
    else:
        doc = doc_or_path
        raw_text = canonical_json_bytes(doc).decode("utf-8")

    require(isinstance(doc, dict), "result document must be a JSON object")
    require(doc.get("format") == RESULT_FORMAT, f"result format {doc.get('format')!r} != {RESULT_FORMAT!r}")
    extra_top = set(doc) - RESULT_TOP_KEYS
    missing_top = RESULT_TOP_KEYS - set(doc)
    require(not extra_top, f"result holds unlisted top-level keys: {sorted(extra_top)}")
    require(not missing_top, f"result missing top-level keys: {sorted(missing_top)}")

    forbidden = _walk_keys(doc)
    require(not forbidden, f"result holds forbidden zero-ID/timing keys: {forbidden[:10]}")

    for label, pat in BANNED_PATTERNS:
        hits = len(pat.findall(raw_text))
        require(hits == 0, f"result hygiene violation: {hits} hit(s) matching banned pattern {label}")

    scan = _load_scan()
    ids: set[str] = set(known_genome_ids or ())
    for sp_path in species_json_paths or ():
        ids |= scan.ids_of(str(sp_path))

    genome_hits = 0
    if ids:
        genome_hits = sum(1 for m in scan.TOKEN.finditer(raw_text) if m.group(0) in ids)
        require(genome_hits == 0, f"zero-ID scan violation: {genome_hits} genome id hit(s) in result")

    return {"genome_id_hits": genome_hits, "banned_pattern_hits": 0, "forbidden_key_hits": 0}


def write_result_atomic(path: Path | str, doc: dict[str, Any], *, known_genome_ids: set[str] | None = None) -> str:
    """Validate hygiene and atomically write deterministic `grm-bench-result-v1` JSON bytes."""
    check_zero_id_hygiene(doc, known_genome_ids=known_genome_ids)
    raw = canonical_json_bytes(doc)
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_name(p.name + ".tmp")
    tmp.write_bytes(raw)
    tmp.replace(p)
    return hashlib.sha256(raw).hexdigest()
