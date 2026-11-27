"""Species-level cross-split aggregator for nested selective calling (`grm-selective-v1`).

Usage:
  python3 -I scripts/selective/summarize.py <input...> --out <out_dir> [--alpha 0.05]
      [--tie-tol 1e-12] [--scan-species-json <species.json>...]

Reads `grm-selective-v1` JSON outputs written by `grm learn --exp-selective <DIR>`
across `lineage` and `study-whole` splits (as well as any other splits present),
validates schema, pre-committed constants, and raw count accounting, and computes:

1. Pre-committed primary contrasts on `vme_called` at target coverage c = 0.85:
   - `score_margin_vs_outcome_blind`: vme_called(outcome_blind) - vme_called(score_margin)
   - `margin_plus_distance_vs_score_margin`: vme_called(score_margin) - vme_called(margin_plus_distance)
   Positive values indicate an incremental reduction in `vme_called` by the arm.
2. Strict sparse-denominator filtering: datasets with called_r < 25 on either side of
   a paired contrast are flagged `excluded_sparse = 1` and excluded from primary
   species medians and sign tests while raw counts stay recorded in `datasets.tsv`.
3. Species-level paired aggregation: each eligible species contributes one paired
   difference (median over its eligible datasets), yielding absolute species
   differences, wins/ties/losses, and the one-sided exact binomial sign test on
   non-tied species (explicitly annotated as testing direction across species, not
   mean magnitude).
4. Cross-split joint criterion: a joint claim requires success on BOTH `lineage`
   and `study-whole` splits (direction positive, wins > losses, sign test p < alpha).
5. Zero-ID hygiene verification using the exact token rule of `scripts/c7-outside/scan.py`.
"""

import argparse
import collections
import json
import math
import os
import re
import statistics
import sys
from pathlib import Path

FORMAT = "grm-selective-v1"
EVALUATION_REVISION = 14
PRIMARY_COVERAGE = 0.85
COVERAGE_LADDER = (1.0, 0.95, 0.9, 0.85, 0.8, 0.7, 0.6, 0.5)
SPARSE_DENOMINATOR_MIN = 25
LEARNERS = ("scm", "cart", "rf", "gbt", "lasso")
POLICIES = ("outcome_blind", "score_margin", "mash_distance", "margin_plus_distance")
LAYERS = ("matched_coverage", "calibrated_operating_point")
REQUIRED_SPLITS = ("lineage", "study-whole")
PRIMARY_CONTRASTS = (
    ("score_margin_vs_outcome_blind", "score_margin", "outcome_blind"),
    ("margin_plus_distance_vs_score_margin", "margin_plus_distance", "score_margin"),
)
REQUIRED_ABSTAIN_KEYS = (
    "null_fit",
    "ood_mash_distance",
    "uncertain_score",
    "conformal_both",
    "conformal_empty",
)
SIGN_TEST_ANNOTATION = (
    "One-sided exact binomial sign test on non-tied species (H0: P(win) <= 0.5); "
    "tests direction of paired differences across species, not mean magnitude."
)
# Exact genome-id token rule from scripts/c7-outside/scan.py.
GENOME_ID_TOKEN = re.compile(r"(?<![\w.])\d+\.\d+(?!\d)")


def r12(value):
    return None if value is None else round(float(value), 12)


def retention_count(coverage, n):
    if coverage >= 1.0:
        return n
    if coverage <= 0.0 or n == 0:
        return 0
    pct = round(coverage * 100.0)
    if abs(coverage * 100.0 - pct) < 1e-12 and 0 <= pct <= 100:
        return (pct * n) // 100
    return int(math.floor(coverage * n))


def exact_one_sided_sign_test(wins, losses):
    """Exact one-sided binomial sign test P(X >= wins | n = wins + losses, p = 0.5)."""
    n = wins + losses
    if n == 0:
        return 1.0
    return sum(math.comb(n, k) for k in range(wins, n + 1)) / (2 ** n)


def _is_nonneg_int(val):
    return isinstance(val, int) and not isinstance(val, bool) and val >= 0


def validate_accounting(acc, context, require_conformal=False):
    for key in ("n", "n_r", "n_s", "called", "called_r", "called_s", "abstain_r", "abstain_s", "tp", "fp", "tn", "fn"):
        if not _is_nonneg_int(acc.get(key)):
            raise ValueError(f"{context}: invalid count field {key!r}")
    if acc["n_r"] + acc["n_s"] != acc["n"]:
        raise ValueError(f"{context}: n_r + n_s != n")
    if acc["tp"] + acc["fn"] != acc["called_r"]:
        raise ValueError(f"{context}: tp + fn != called_r")
    if acc["tn"] + acc["fp"] != acc["called_s"]:
        raise ValueError(f"{context}: tn + fp != called_s")
    if acc["called_r"] + acc["called_s"] != acc["called"]:
        raise ValueError(f"{context}: called_r + called_s != called")
    if acc["called_r"] + acc["abstain_r"] != acc["n_r"]:
        raise ValueError(f"{context}: called_r + abstain_r != n_r")
    if acc["called_s"] + acc["abstain_s"] != acc["n_s"]:
        raise ValueError(f"{context}: called_s + abstain_s != n_s")

    reasons = acc.get("abstain_reasons")
    if not isinstance(reasons, dict):
        raise ValueError(f"{context}: missing abstain_reasons object")
    for k in REQUIRED_ABSTAIN_KEYS:
        if not _is_nonneg_int(reasons.get(k)):
            raise ValueError(f"{context}: invalid abstain_reasons.{k}")
    if any(not _is_nonneg_int(v) for v in reasons.values()):
        raise ValueError(f"{context}: non-integer value in abstain_reasons")
    if sum(reasons.values()) != acc["abstain_r"] + acc["abstain_s"]:
        raise ValueError(f"{context}: abstain_reasons sum != abstain_r + abstain_s")

    def check_rate(name, want_x, want_n):
        obj = acc.get(name)
        if not isinstance(obj, dict):
            raise ValueError(f"{context}: missing rate object {name}")
        if obj.get("x") != want_x or obj.get("n") != want_n:
            raise ValueError(f"{context}: {name} counts ({obj.get('x')}/{obj.get('n')}) != ({want_x}/{want_n})")
        want_sparse = want_n < SPARSE_DENOMINATOR_MIN
        if obj.get("sparse_denominator") is not want_sparse:
            raise ValueError(f"{context}: {name}.sparse_denominator mismatch")
        if want_n == 0:
            if obj.get("rate") is not None:
                raise ValueError(f"{context}: {name}.rate must be null when n == 0")
        else:
            if not math.isclose(obj["rate"], want_x / want_n, abs_tol=1e-12):
                raise ValueError(f"{context}: {name}.rate differs from x / n")

    check_rate("vme_called", acc["fn"], acc["called_r"])
    check_rate("me_called", acc["fp"], acc["called_s"])
    check_rate("npv", acc["tn"], acc["tn"] + acc["fn"])
    check_rate("for_s", acc["fn"], acc["tn"] + acc["fn"])
    check_rate("ppv", acc["tp"], acc["tp"] + acc["fp"])
    check_rate("vme_all", acc["fn"], acc["n_r"])
    check_rate("me_all", acc["fp"], acc["n_s"])

    want_sparse_r = acc["called_r"] < SPARSE_DENOMINATOR_MIN
    want_sparse_s = acc["called_s"] < SPARSE_DENOMINATOR_MIN
    if acc.get("sparse_denominator_r") is not want_sparse_r:
        raise ValueError(f"{context}: sparse_denominator_r mismatch")
    if acc.get("sparse_denominator_s") is not want_sparse_s:
        raise ValueError(f"{context}: sparse_denominator_s mismatch")
    if acc.get("sparse_denominator") is not (want_sparse_r or want_sparse_s):
        raise ValueError(f"{context}: sparse_denominator mismatch")

    if require_conformal and not isinstance(acc.get("set_true_label_coverage"), dict):
        raise ValueError(f"{context}: conformal point missing set_true_label_coverage")


def validate_doc(doc, source="<memory>"):
    if not isinstance(doc, dict):
        raise ValueError(f"{source}: top-level JSON must be an object")
    if doc.get("format") != FORMAT:
        raise ValueError(f"{source}: format {doc.get('format')!r} != {FORMAT!r}")
    if doc.get("evaluation_revision") != EVALUATION_REVISION:
        raise ValueError(f"{source}: evaluation_revision {doc.get('evaluation_revision')!r} != {EVALUATION_REVISION}")
    for key in ("dataset", "species_slug", "antibiotic", "scheme"):
        if not isinstance(doc.get(key), str) or not doc[key]:
            raise ValueError(f"{source}: missing string field {key!r}")
    if not _is_nonneg_int(doc.get("samples")):
        raise ValueError(f"{source}: samples must be a non-negative integer count, never a list")
    for key in ("seed", "folds", "inner_folds", "tested", "tested_resistant", "tested_susceptible"):
        if not _is_nonneg_int(doc.get(key)):
            raise ValueError(f"{source}: invalid integer field {key!r}")

    consts = doc.get("constants")
    if not isinstance(consts, dict):
        raise ValueError(f"{source}: missing constants metadata")
    if not math.isclose(float(consts.get("primary_coverage", -1.0)), PRIMARY_COVERAGE, abs_tol=1e-12):
        raise ValueError(f"{source}: primary_coverage != {PRIMARY_COVERAGE}")
    if consts.get("sparse_denominator_min") != SPARSE_DENOMINATOR_MIN:
        raise ValueError(f"{source}: sparse_denominator_min != {SPARSE_DENOMINATOR_MIN}")
    ladder = consts.get("coverage_ladder")
    if not isinstance(ladder, list) or len(ladder) != len(COVERAGE_LADDER) or not all(
        math.isclose(float(a), b, abs_tol=1e-12) for a, b in zip(ladder, COVERAGE_LADDER)
    ):
        raise ValueError(f"{source}: coverage_ladder mismatch")

    if doc.get("unavailable_reason") is not None:
        if doc.get("learners"):
            raise ValueError(f"{source}: unavailable split must have empty learners")
        return

    if doc["tested_resistant"] + doc["tested_susceptible"] != doc["tested"]:
        raise ValueError(f"{source}: tested_resistant + tested_susceptible != tested")

    learners = doc.get("learners")
    if not isinstance(learners, list) or [lr.get("learner") for lr in learners] != list(LEARNERS):
        raise ValueError(f"{source}: learners must match {list(LEARNERS)} in pipeline order")

    for lr in learners:
        lname = lr["learner"]
        policies = lr.get("policies")
        if not isinstance(policies, list) or [p.get("policy") for p in policies] != list(POLICIES):
            raise ValueError(f"{source} {lname}: policies must match {list(POLICIES)} in order")
        for pol in policies:
            pname = pol["policy"]
            if not pol.get("available"):
                if not pol.get("unavailable_reason"):
                    raise ValueError(f"{source} {lname} {pname}: unavailable policy must state unavailable_reason")
                continue
            for layer in LAYERS:
                pts = pol.get(layer)
                if not isinstance(pts, list) or len(pts) != len(COVERAGE_LADDER):
                    raise ValueError(f"{source} {lname} {pname} {layer}: expected {len(COVERAGE_LADDER)} rungs")
                for idx, (pt, want_c) in enumerate(zip(pts, COVERAGE_LADDER)):
                    ctx = f"{source} {lname} {pname} {layer}[{idx}]"
                    if not math.isclose(float(pt.get("target_coverage", -1.0)), want_c, abs_tol=1e-12):
                        raise ValueError(f"{ctx}: target_coverage mismatch")
                    validate_accounting(pt, ctx, require_conformal=False)
                    if pt["n"] != doc["tested"] or pt["n_r"] != doc["tested_resistant"]:
                        raise ValueError(f"{ctx}: accounting totals differ from split tested counts")
                    folds = pt.get("folds", [])
                    if not isinstance(folds, list) or len(folds) != doc["folds"]:
                        raise ValueError(f"{ctx}: expected {doc['folds']} fold records")
                    if sum(f["called"] for f in folds) != pt["called"] or sum(f["n"] for f in folds) != pt["n"]:
                        raise ValueError(f"{ctx}: fold counts do not sum to pooled counts")
                    if layer == "matched_coverage":
                        for f in folds:
                            want_keep = retention_count(want_c, f["n"])
                            if f["called"] != want_keep:
                                raise ValueError(
                                    f"{ctx} fold {f['fold']}: matched_coverage called={f['called']} "
                                    f"!= floor({want_c} * {f['n']})={want_keep}"
                                )
        if lr.get("mondrian_conformal_available"):
            for ci, pt in enumerate(lr.get("mondrian_conformal", [])):
                validate_accounting(pt, f"{source} {lname} mondrian_conformal[{ci}]", require_conformal=True)


def discover_result_files(paths):
    found = []
    for raw in paths:
        p = Path(raw)
        if p.is_file():
            found.append(p)
        elif p.is_dir():
            for child in sorted(p.rglob("*.selective.json")):
                if child.is_file():
                    found.append(child)
        else:
            raise FileNotFoundError(f"input path does not exist: {p}")
    return sorted(set(found), key=lambda x: x.as_posix())


def load_results(paths):
    files = discover_result_files(paths)
    docs = []
    seen = set()
    for path in files:
        doc = json.loads(path.read_text(encoding="utf-8"))
        validate_doc(doc, source=str(path))
        key = (doc["dataset"], doc["scheme"], doc["seed"])
        if key in seen:
            raise ValueError(f"duplicate result for dataset={key[0]}, scheme={key[1]}, seed={key[2]}")
        seen.add(key)
        docs.append(doc)
    return docs


def _find_coverage_point(points, target_c):
    for pt in points:
        if math.isclose(float(pt["target_coverage"]), target_c, abs_tol=1e-12):
            return pt
    return None


def summarize(docs, alpha=0.05, tie_tol=1e-12):
    """Compute dataset paired rows, species differences, exact sign tests, and joint claims."""
    if not (0.0 < alpha < 1.0):
        raise ValueError(f"alpha must lie in (0, 1), not {alpha}")
    if tie_tol < 0.0:
        raise ValueError(f"tie_tol must be non-negative, not {tie_tol}")

    # Group docs by (dataset, scheme) to handle single-seed and multi-seed cleanly.
    grouped = collections.defaultdict(list)
    for doc in docs:
        grouped[(doc["dataset"], doc["scheme"])].append(doc)
    for k in grouped:
        grouped[k].sort(key=lambda d: d["seed"])

    dataset_rows = []
    schemes_seen = sorted({doc["scheme"] for doc in docs})

    for layer in LAYERS:
        for (dataset, scheme), seed_docs in sorted(grouped.items()):
            species_slug = seed_docs[0]["species_slug"]
            antibiotic = seed_docs[0]["antibiotic"]
            for learner in LEARNERS:
                for contrast_name, arm_policy, base_policy in PRIMARY_CONTRASTS:
                    deltas = []
                    arm_vmes = []
                    base_vmes = []
                    arm_covs = []
                    base_covs = []
                    arm_cov_rs = []
                    base_cov_rs = []
                    min_called_r_arm = None
                    min_called_r_base = None
                    sum_n = sum_n_r = sum_called_r_arm = sum_called_r_base = 0
                    sum_fn_arm = sum_fn_base = 0
                    excluded_reason = ""
                    excluded_sparse = False
                    unavailable = False

                    for doc in seed_docs:
                        if doc.get("unavailable_reason") is not None:
                            unavailable = True
                            excluded_reason = "split_unavailable"
                            break
                        lr = next(x for x in doc["learners"] if x["learner"] == learner)
                        pol_map = {p["policy"]: p for p in lr["policies"]}
                        p_arm = pol_map[arm_policy]
                        p_base = pol_map[base_policy]
                        if not p_arm["available"] or not p_base["available"]:
                            unavailable = True
                            excluded_reason = "policy_unavailable"
                            break
                        pt_arm = _find_coverage_point(p_arm[layer], PRIMARY_COVERAGE)
                        pt_base = _find_coverage_point(p_base[layer], PRIMARY_COVERAGE)
                        if pt_arm is None or pt_base is None:
                            unavailable = True
                            excluded_reason = "missing_primary_coverage"
                            break

                        sum_n += pt_arm["n"]
                        sum_n_r += pt_arm["n_r"]
                        cr_arm = pt_arm["called_r"]
                        cr_base = pt_base["called_r"]
                        sum_called_r_arm += cr_arm
                        sum_called_r_base += cr_base
                        sum_fn_arm += pt_arm["fn"]
                        sum_fn_base += pt_base["fn"]
                        min_called_r_arm = cr_arm if min_called_r_arm is None else min(min_called_r_arm, cr_arm)
                        min_called_r_base = cr_base if min_called_r_base is None else min(min_called_r_base, cr_base)

                        if pt_arm["vme_called"]["sparse_denominator"] or pt_base["vme_called"]["sparse_denominator"]:
                            excluded_sparse = True
                            excluded_reason = "sparse_denominator"

                        r_arm = pt_arm["vme_called"]["rate"]
                        r_base = pt_base["vme_called"]["rate"]
                        if r_arm is not None and r_base is not None:
                            arm_vmes.append(r_arm)
                            base_vmes.append(r_base)
                            deltas.append(r_base - r_arm)
                        if pt_arm["coverage"] is not None and pt_base["coverage"] is not None:
                            arm_covs.append(pt_arm["coverage"])
                            base_covs.append(pt_base["coverage"])
                        if pt_arm["coverage_r"] is not None and pt_base["coverage_r"] is not None:
                            arm_cov_rs.append(pt_arm["coverage_r"])
                            base_cov_rs.append(pt_base["coverage_r"])

                    eligible = (not unavailable) and (not excluded_sparse) and len(deltas) == len(seed_docs)
                    mean_delta = r12(sum(deltas) / len(deltas)) if len(deltas) == len(seed_docs) and deltas else None
                    dataset_rows.append(
                        {
                            "layer": layer,
                            "scheme": scheme,
                            "species_slug": species_slug,
                            "dataset": dataset,
                            "antibiotic": antibiotic,
                            "learner": learner,
                            "contrast": contrast_name,
                            "arm_policy": arm_policy,
                            "base_policy": base_policy,
                            "target_coverage": PRIMARY_COVERAGE,
                            "seeds": len(seed_docs),
                            "n": sum_n,
                            "n_r": sum_n_r,
                            "called_r_arm": sum_called_r_arm,
                            "called_r_base": sum_called_r_base,
                            "min_called_r_arm": min_called_r_arm,
                            "min_called_r_base": min_called_r_base,
                            "fn_arm": sum_fn_arm,
                            "fn_base": sum_fn_base,
                            "arm_vme_called": r12(sum(arm_vmes) / len(arm_vmes)) if arm_vmes else None,
                            "base_vme_called": r12(sum(base_vmes) / len(base_vmes)) if base_vmes else None,
                            "vme_reduction": mean_delta,
                            "arm_minus_base": r12(-mean_delta) if mean_delta is not None else None,
                            "arm_coverage": r12(sum(arm_covs) / len(arm_covs)) if arm_covs else None,
                            "base_coverage": r12(sum(base_covs) / len(base_covs)) if base_covs else None,
                            "arm_coverage_r": r12(sum(arm_cov_rs) / len(arm_cov_rs)) if arm_cov_rs else None,
                            "base_coverage_r": r12(sum(base_cov_rs) / len(base_cov_rs)) if base_cov_rs else None,
                            "excluded_sparse": int(excluded_sparse),
                            "unavailable": int(unavailable),
                            "eligible": int(eligible),
                            "exclusion_reason": excluded_reason,
                        }
                    )

    # Aggregate by species.
    species_rows = []
    by_species_bucket = collections.defaultdict(list)
    for row in dataset_rows:
        by_species_bucket[
            (row["layer"], row["scheme"], row["learner"], row["contrast"], row["species_slug"])
        ].append(row)

    for (layer, scheme, learner, contrast, sp), rows in sorted(by_species_bucket.items()):
        elig = [r for r in rows if r["eligible"]]
        sparse_excl = sum(1 for r in rows if r["excluded_sparse"] and not r["unavailable"])
        unavail = sum(1 for r in rows if r["unavailable"])
        if elig:
            med_red = r12(statistics.median(r["vme_reduction"] for r in elig))
            mean_red = r12(sum(r["vme_reduction"] for r in elig) / len(elig))
            med_arm = r12(statistics.median(r["arm_vme_called"] for r in elig))
            med_base = r12(statistics.median(r["base_vme_called"] for r in elig))
            med_cov_arm = r12(statistics.median(r["arm_coverage"] for r in elig))
            med_cov_base = r12(statistics.median(r["base_coverage"] for r in elig))
            med_cov_r_arm = r12(statistics.median(r["arm_coverage_r"] for r in elig))
            med_cov_r_base = r12(statistics.median(r["base_coverage_r"] for r in elig))
            if med_red > tie_tol:
                direction = "win"
            elif med_red < -tie_tol:
                direction = "loss"
            else:
                direction = "tie"
        else:
            med_red = mean_red = med_arm = med_base = None
            med_cov_arm = med_cov_base = med_cov_r_arm = med_cov_r_base = None
            direction = "ineligible"

        species_rows.append(
            {
                "layer": layer,
                "scheme": scheme,
                "learner": learner,
                "contrast": contrast,
                "species_slug": sp,
                "datasets_total": len(rows),
                "datasets_eligible": len(elig),
                "datasets_sparse_excluded": sparse_excl,
                "datasets_unavailable": unavail,
                "median_vme_reduction": med_red,
                "mean_vme_reduction": mean_red,
                "median_arm_vme_called": med_arm,
                "median_base_vme_called": med_base,
                "median_arm_coverage": med_cov_arm,
                "median_base_coverage": med_cov_base,
                "median_arm_coverage_r": med_cov_r_arm,
                "median_base_coverage_r": med_cov_r_base,
                "direction": direction,
            }
        )

    # Exact one-sided binomial sign test across species per (layer, scheme, learner, contrast).
    sign_test_rows = []
    sign_index = {}
    for layer in LAYERS:
        for scheme in schemes_seen:
            for learner in LEARNERS:
                for contrast_name, arm_policy, base_policy in PRIMARY_CONTRASTS:
                    sp_rows = [
                        r
                        for r in species_rows
                        if (r["layer"], r["scheme"], r["learner"], r["contrast"])
                        == (layer, scheme, learner, contrast_name)
                    ]
                    elig_sp = [r for r in sp_rows if r["datasets_eligible"] > 0]
                    wins = sum(1 for r in elig_sp if r["direction"] == "win")
                    losses = sum(1 for r in elig_sp if r["direction"] == "loss")
                    ties = sum(1 for r in elig_sp if r["direction"] == "tie")
                    non_ties = wins + losses
                    med_species = (
                        r12(statistics.median(r["median_vme_reduction"] for r in elig_sp)) if elig_sp else None
                    )
                    p_val = r12(exact_one_sided_sign_test(wins, losses))
                    ds_total = sum(r["datasets_total"] for r in sp_rows)
                    ds_elig = sum(r["datasets_eligible"] for r in sp_rows)
                    ds_sparse = sum(r["datasets_sparse_excluded"] for r in sp_rows)
                    ds_unavail = sum(r["datasets_unavailable"] for r in sp_rows)
                    split_pass = (
                        non_ties > 0
                        and wins > losses
                        and med_species is not None
                        and med_species > tie_tol
                        and p_val < alpha
                    )
                    row = {
                        "layer": layer,
                        "scheme": scheme,
                        "learner": learner,
                        "contrast": contrast_name,
                        "arm_policy": arm_policy,
                        "base_policy": base_policy,
                        "target_coverage": PRIMARY_COVERAGE,
                        "datasets_total": ds_total,
                        "datasets_eligible": ds_elig,
                        "datasets_sparse_excluded": ds_sparse,
                        "datasets_unavailable": ds_unavail,
                        "species_total": len(sp_rows),
                        "species_eligible": len(elig_sp),
                        "wins": wins,
                        "ties": ties,
                        "losses": losses,
                        "non_ties": non_ties,
                        "median_species_vme_reduction": med_species,
                        "sign_test_p": p_val,
                        "alpha": alpha,
                        "split_pass": int(split_pass),
                        "sign_test_note": SIGN_TEST_ANNOTATION,
                    }
                    sign_test_rows.append(row)
                    sign_index[(layer, scheme, learner, contrast_name)] = row

    # Joint cross-split verdict requiring BOTH lineage AND study-whole to pass.
    joint_rows = []
    for layer in LAYERS:
        for learner in LEARNERS:
            for contrast_name, arm_policy, base_policy in PRIMARY_CONTRASTS:
                lin = sign_index.get((layer, "lineage", learner, contrast_name))
                stu = sign_index.get((layer, "study-whole", learner, contrast_name))
                lin_pass = bool(lin and lin["split_pass"])
                stu_pass = bool(stu and stu["split_pass"])
                joint_ok = lin_pass and stu_pass
                reasons = []
                if lin is None:
                    reasons.append("lineage split missing")
                elif not lin_pass:
                    reasons.append("lineage did not pass")
                if stu is None:
                    reasons.append("study-whole split missing")
                elif not stu_pass:
                    reasons.append("study-whole did not pass")
                reason = "both lineage and study-whole passed" if joint_ok else "; ".join(reasons)
                joint_rows.append(
                    {
                        "layer": layer,
                        "learner": learner,
                        "contrast": contrast_name,
                        "arm_policy": arm_policy,
                        "base_policy": base_policy,
                        "target_coverage": PRIMARY_COVERAGE,
                        "alpha": alpha,
                        "lineage_present": int(lin is not None),
                        "lineage_species_eligible": lin["species_eligible"] if lin else 0,
                        "lineage_wins": lin["wins"] if lin else 0,
                        "lineage_ties": lin["ties"] if lin else 0,
                        "lineage_losses": lin["losses"] if lin else 0,
                        "lineage_median_vme_reduction": lin["median_species_vme_reduction"] if lin else None,
                        "lineage_sign_test_p": lin["sign_test_p"] if lin else None,
                        "lineage_pass": int(lin_pass),
                        "study_whole_present": int(stu is not None),
                        "study_whole_species_eligible": stu["species_eligible"] if stu else 0,
                        "study_whole_wins": stu["wins"] if stu else 0,
                        "study_whole_ties": stu["ties"] if stu else 0,
                        "study_whole_losses": stu["losses"] if stu else 0,
                        "study_whole_median_vme_reduction": stu["median_species_vme_reduction"] if stu else None,
                        "study_whole_sign_test_p": stu["sign_test_p"] if stu else None,
                        "study_whole_pass": int(stu_pass),
                        "joint_claim_supported": int(joint_ok),
                        "reason": reason,
                        "sign_test_note": SIGN_TEST_ANNOTATION,
                    }
                )

    return {
        "format": "grm-selective-summary-v1",
        "source_format": FORMAT,
        "evaluation_revision": EVALUATION_REVISION,
        "primary_coverage": PRIMARY_COVERAGE,
        "sparse_denominator_min": SPARSE_DENOMINATOR_MIN,
        "alpha": alpha,
        "tie_tol": tie_tol,
        "required_joint_splits": list(REQUIRED_SPLITS),
        "sign_test_note": SIGN_TEST_ANNOTATION,
        "datasets": dataset_rows,
        "species_differences": species_rows,
        "sign_tests": sign_test_rows,
        "joint_claims": joint_rows,
    }


def _fmt(value):
    if value is None:
        return ""
    if isinstance(value, bool):
        return str(int(value))
    if isinstance(value, float):
        return "" if math.isnan(value) else repr(round(value, 12))
    return str(value)


def write_tsv(path, rows, columns):
    with open(path, "w", encoding="utf-8", newline="") as f:
        f.write("\t".join(columns) + "\n")
        for row in rows:
            f.write("\t".join(_fmt(row.get(c)) for c in columns) + "\n")


DATASET_COLUMNS = [
    "layer", "scheme", "species_slug", "dataset", "antibiotic", "learner", "contrast",
    "arm_policy", "base_policy", "target_coverage", "seeds", "n", "n_r",
    "called_r_arm", "called_r_base", "min_called_r_arm", "min_called_r_base",
    "fn_arm", "fn_base", "arm_vme_called", "base_vme_called", "vme_reduction",
    "arm_minus_base", "arm_coverage", "base_coverage", "arm_coverage_r",
    "base_coverage_r", "excluded_sparse", "unavailable", "eligible", "exclusion_reason",
]
SPECIES_COLUMNS = [
    "layer", "scheme", "learner", "contrast", "species_slug", "datasets_total",
    "datasets_eligible", "datasets_sparse_excluded", "datasets_unavailable",
    "median_vme_reduction", "mean_vme_reduction", "median_arm_vme_called",
    "median_base_vme_called", "median_arm_coverage", "median_base_coverage",
    "median_arm_coverage_r", "median_base_coverage_r", "direction",
]
SIGN_TEST_COLUMNS = [
    "layer", "scheme", "learner", "contrast", "arm_policy", "base_policy",
    "target_coverage", "datasets_total", "datasets_eligible",
    "datasets_sparse_excluded", "datasets_unavailable", "species_total",
    "species_eligible", "wins", "ties", "losses", "non_ties",
    "median_species_vme_reduction", "sign_test_p", "alpha", "split_pass",
    "sign_test_note",
]
JOINT_COLUMNS = [
    "layer", "learner", "contrast", "arm_policy", "base_policy", "target_coverage",
    "alpha", "lineage_present", "lineage_species_eligible", "lineage_wins",
    "lineage_ties", "lineage_losses", "lineage_median_vme_reduction",
    "lineage_sign_test_p", "lineage_pass", "study_whole_present",
    "study_whole_species_eligible", "study_whole_wins", "study_whole_ties",
    "study_whole_losses", "study_whole_median_vme_reduction",
    "study_whole_sign_test_p", "study_whole_pass", "joint_claim_supported",
    "reason", "sign_test_note",
]


def write_summary(out_dir, summary):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    write_tsv(out / "datasets.tsv", summary["datasets"], DATASET_COLUMNS)
    write_tsv(out / "species_differences.tsv", summary["species_differences"], SPECIES_COLUMNS)
    write_tsv(out / "sign_tests.tsv", summary["sign_tests"], SIGN_TEST_COLUMNS)
    write_tsv(out / "joint_claims.tsv", summary["joint_claims"], JOINT_COLUMNS)
    (out / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def scan_ids_against_species_json(species_json_paths, target_paths):
    """Verify zero-ID hygiene using the exact token rule of scripts/c7-outside/scan.py."""
    ids = set()
    for sp in species_json_paths:
        with open(sp, encoding="utf-8") as f:
            ids.update(g["id"] for g in json.load(f)["genomes"])
    if not ids:
        raise ValueError("species json files list no genome ids")

    files_to_scan = []
    for raw in target_paths:
        p = Path(raw)
        if p.is_file():
            files_to_scan.append(p)
        elif p.is_dir():
            for root, _, names in os.walk(p):
                for name in sorted(names):
                    files_to_scan.append(Path(root) / name)

    hits = 0
    for fp in sorted(files_to_scan, key=lambda x: x.as_posix()):
        text = fp.read_text(encoding="utf-8", errors="replace")
        hits += sum(1 for m in GENOME_ID_TOKEN.finditer(text) if m.group(0) in ids)
    return hits, len(files_to_scan)


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Species-level cross-split summary and exact sign test for grm-selective-v1 results."
    )
    parser.add_argument("inputs", nargs="+", help="Input .selective.json files or directories.")
    parser.add_argument("--out", required=True, type=Path, help="Output directory for TSV and JSON summaries.")
    parser.add_argument("--alpha", type=float, default=0.05, help="One-sided sign-test threshold (default: 0.05).")
    parser.add_argument("--tie-tol", type=float, default=1e-12, help="Absolute tie tolerance on paired differences.")
    parser.add_argument(
        "--scan-species-json",
        action="append",
        default=[],
        help="Optional release species JSON file(s) to run scripts/c7-outside/scan.py zero-ID verification.",
    )
    args = parser.parse_args(argv)

    docs = load_results(args.inputs)
    if not docs:
        print("no grm-selective-v1 files found", file=sys.stderr)
        return 1

    summary = summarize(docs, alpha=args.alpha, tie_tol=args.tie_tol)
    write_summary(args.out, summary)

    if args.scan_species_json:
        hits, n_files = scan_ids_against_species_json(args.scan_species_json, [*args.inputs, args.out])
        print(f"genome ids: {hits} in {n_files} files")
        if hits != 0:
            return 1

    supported = sum(r["joint_claim_supported"] for r in summary["joint_claims"])
    print(
        f"selective summary: {len(docs)} input files, {len(summary['datasets'])} dataset rows, "
        f"{len(summary['species_differences'])} species rows, {supported} of {len(summary['joint_claims'])} "
        f"joint claims supported"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
