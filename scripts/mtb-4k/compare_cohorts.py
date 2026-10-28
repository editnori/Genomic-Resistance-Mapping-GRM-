#!/usr/bin/env python3
"""Matched-cohort comparator for the MTB 2k / added-2k / 4k design (stdlib only).

usage: python3 -I compare_cohorts.py
           --cohort prefix=ROOT --cohort added=ROOT --cohort full=ROOT
           --cohort matched2k=ROOT
           --checker-receipt prefix=PATH ... (one per cohort)
           --duplicates SCREEN.json [--accept-reviewed-duplicates SHA256]
           --out OUT.json [--slug mycobacterium_tuberculosis]

ROOT is a cohort's data root: reports are read from
species/<slug>/results/*.json and the cohort from species/<slug>/genomes.txt.
Each --checker-receipt is the output of the reviewed revision-8 checker for
that cohort (see the execution plan); this script records its SHA-256 and
refuses to run without one, but does not replace it.

What it does, all from the saved reports:

1. Per report: requires schema 2, evaluation revision 8 and
   evaluation.valid; recomputes sample and class counts, every learner's
   confusion matrix and no-calls from its saved held-out scores at 0.5, and
   per-fold denominators, and requires them to equal the stored rows; checks
   that every revision-8 lasso row has one solver record per fold, none
   unconverged. Anything else marks the report invalid (with reasons).
2. Not-interpretable drugs: from the full cohort's labels, the added
   genomes' resistant and susceptible counts per drug; a drug with fewer than
   MIN_CLASS of either is not interpretable for before/after claims. The
   derived set must equal the frozen expectation (EXPECTED_SPARSE_ADDED), or
   the run stops.
3. Eligibility and the invalid-in-one-cohort rule: a drug is eligible in a
   cohort when the cohort's genomes hold at least MIN_CLASS resistant and
   MIN_CLASS susceptible labels (from the full cohort's report, the same
   labels `grm learn` uses); the added cohort's eligible set must equal the
   frozen 14 (EXPECTED_ADDED_ELIGIBLE). An ineligible drug is listed, not an
   error. A drug whose report is missing or invalid in a cohort it is
   eligible in is dropped from every comparison that involves that cohort
   and from the size-matched decomposition; it is listed under `excluded`,
   never imputed, and the run exits 3 so a parent reviews it before any
   before/after claim is written.
4. Comparisons, on the primary split (lineage) and the secondary (random):
   - full vs prefix, paired by genome ID on the prefix isolates: per learner,
     correct/wrong/no-call transitions at 0.5 (fixed, broken), score deltas,
     and an exact McNemar p-value on fixed vs broken, Holm-adjusted across
     the interpretable drugs of each learner and split; descriptive only;
   - added vs prefix, side by side (disjoint cohorts, not paired);
   - the size-matched control: matched2k (a frozen random 2,000 of the 4,000)
     vs prefix (composition at equal n) and full vs matched2k (sample size at
     fixed composition, paired on matched2k's isolates).
   Each row has per-learner metrics with Clopper-Pearson 95% intervals for
   sensitivity and specificity, and the median of the five learners'
   balanced accuracy (only when all five have metrics).
5. Duplicates: the screen_duplicates.py output must report no exact
   cross-cohort duplicate, unless --accept-reviewed-duplicates gives that
   file's SHA-256 (a recorded parent decision).

Exit status: 0 complete; 2 input or integrity error (no output); 3 output
written but a drug was excluded or a report invalid (stop for review).
"""

import argparse
import hashlib
import json
import math
import statistics
import sys
from pathlib import Path

LEARNERS = ["scm", "cart", "rf", "gbt", "lasso"]
SCHEMES = ["lineage", "random"]  # the first is the primary split
COHORTS = ["prefix", "added", "full", "matched2k"]
MIN_CLASS = 25
THRESHOLD = 0.5
# Added-cohort class counts from the frozen cohort evidence (full 4k minus
# the 2k prefix): before/after differences for these drugs mostly reflect
# prevalence and composition, not added information.
EXPECTED_ADDED_ELIGIBLE = {
    "amikacin", "capreomycin", "delamanid", "ethambutol", "ethionamide", "isoniazid", "kanamycin",
    "moxifloxacin", "ofloxacin", "prothionamide", "pyrazinamide", "rifabutin", "rifampin", "streptomycin",
}
EXPECTED_SPARSE_ADDED = {
    "amoxicillin": "no resistant isolates among the added 2,000",
    "ciprofloxacin": "no resistant isolates among the added 2,000",
    "clofazimine": "no resistant isolates among the added 2,000",
    "cycloserine": "no resistant isolates among the added 2,000",
    "linezolid": "no resistant isolates among the added 2,000",
    "nicotinamide": "no resistant isolates among the added 2,000",
    "para-aminosalicylic acid": "5 resistant isolates among the added 2,000",
    "levofloxacin": "774 resistant and 1 susceptible isolate among the added 2,000",
}


class Stop(Exception):
    pass


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# ---------------------------------------------------------------- statistics

def _betacf(a, b, x):
    """Continued fraction for the regularized incomplete beta (Lentz)."""
    tiny, eps = 1e-300, 3e-16
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c, d = 1.0, 1.0 - qab * x / qap
    d = 1.0 / (d if abs(d) > tiny else tiny)
    h = d
    for m in range(1, 10000):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        d = 1.0 / (d if abs(d) > tiny else tiny)
        c = 1.0 + aa / c
        c = c if abs(c) > tiny else tiny
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        d = 1.0 / (d if abs(d) > tiny else tiny)
        c = 1.0 + aa / c
        c = c if abs(c) > tiny else tiny
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < eps:
            return h
    raise Stop("incomplete beta did not converge")


def betainc(a, b, x):
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    lbeta = math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
    front = math.exp(lbeta + a * math.log(x) + b * math.log1p(-x))
    if x < (a + 1.0) / (a + b + 2.0):
        return front * _betacf(a, b, x) / a
    return 1.0 - front * _betacf(b, a, 1.0 - x) / b


def beta_quantile(q, a, b):
    lo, hi = 0.0, 1.0
    for _ in range(200):
        mid = (lo + hi) / 2
        if betainc(a, b, mid) < q:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def clopper_pearson(k, n, alpha=0.05):
    """Exact two-sided interval for k successes in n; None when n == 0."""
    if n == 0:
        return None
    lower = 0.0 if k == 0 else beta_quantile(alpha / 2, k, n - k + 1)
    upper = 1.0 if k == n else beta_quantile(1 - alpha / 2, k + 1, n - k)
    return [lower, upper]


def mcnemar_exact(b, c):
    """Two-sided exact McNemar p-value from discordant counts b and c."""
    n = b + c
    if n == 0:
        return 1.0
    k = min(b, c)
    tail = sum(math.comb(n, i) for i in range(k + 1)) / 2 ** n
    return min(1.0, 2 * tail)


def holm(pvalues):
    """Holm step-down adjustment of {key: p}."""
    order = sorted(pvalues, key=lambda key: pvalues[key])
    m, running, out = len(order), 0.0, {}
    for i, key in enumerate(order):
        running = max(running, min(1.0, (m - i) * pvalues[key]))
        out[key] = running
    return out


# ------------------------------------------------------------------ reports

def ratio(a, b):
    return a / b if b else 0.0


def recompute(truth, scores, folds):
    """The confusion matrix, no-calls and per-fold denominators of saved
    held-out scores (None = no call) at THRESHOLD."""
    tp = fp = tn = fn = nc_r = nc_s = 0
    per_fold = {}
    for t, s, f in zip(truth, scores, folds):
        d = per_fold.setdefault(f, {"fold": f, "resistant": 0, "susceptible": 0, "called_resistant": 0, "called_susceptible": 0})
        d["resistant" if t else "susceptible"] += 1
        if s is None:
            if t:
                nc_r += 1
            else:
                nc_s += 1
            continue
        d["called_resistant" if t else "called_susceptible"] += 1
        call = s >= THRESHOLD
        if t and call:
            tp += 1
        elif t:
            fn += 1
        elif call:
            fp += 1
        else:
            tn += 1
    sens, spec = ratio(tp, tp + fn), ratio(tn, tn + fp)
    return {
        "n_called": tp + fp + tn + fn, "tp": tp, "fp": fp, "tn": tn, "fn": fn,
        "no_call_resistant": nc_r, "no_call_susceptible": nc_s,
        "sensitivity": sens, "specificity": spec, "balanced_accuracy": (sens + spec) / 2,
        "vme": ratio(fn, tp + fn), "me": ratio(fp, fp + tn),
        "sensitivity_ci": clopper_pearson(tp, tp + fn), "specificity_ci": clopper_pearson(tn, tn + fp),
        "folds": [per_fold[f] for f in sorted(per_fold)],
    }


def lasso_unconverged(fit):
    if fit is None:
        return True
    conv = [fit.get("selected_convergence"), fit.get("refit_convergence")]
    return fit.get("path_unconverged", 1) > 0 or any(c is not None and not c.get("converged", False) for c in conv)


def load_report(path: Path):
    """(report summary, reasons it is invalid)."""
    rep = json.loads(path.read_text())
    reasons = []
    if rep.get("schema") != 2:
        reasons.append(f"schema {rep.get('schema')} is not 2")
    ev = rep.get("evaluation") or {}
    if ev.get("revision") != 8:
        reasons.append(f"evaluation revision {ev.get('revision')} is not 8")
    if ev.get("valid") is not True:
        reasons.append("evaluation.valid is not true: " + "; ".join(ev.get("reasons") or []))
    samples = rep["samples"]
    ids = [s["id"] for s in samples]
    truth = [bool(s["resistant"]) for s in samples]
    counts = rep["counts"]
    if len(set(ids)) != len(ids):
        reasons.append("repeated sample ids")
    if (counts["samples"], counts["resistant"], counts["susceptible"]) != (len(ids), sum(truth), len(truth) - sum(truth)):
        reasons.append("stored sample or class counts differ from the samples")
    splits = {}
    for split in rep["splits"]:
        scheme = split["scheme"]
        if split.get("unavailable_reason"):
            splits[scheme] = {"unavailable": split["unavailable_reason"]}
            continue
        folds = split["test_fold"]
        if len(folds) != len(ids):
            reasons.append(f"{scheme}: {len(folds)} test folds for {len(ids)} samples")
            continue
        rows = {}
        for row in split["learners"]:
            name = row["learner"]
            if name not in LEARNERS:
                continue
            scores = row["scores"]
            if len(scores) != len(ids):
                reasons.append(f"{scheme}/{name}: {len(scores)} scores for {len(ids)} samples")
                continue
            m = recompute(truth, scores, folds)
            stored = row.get("metrics")
            if m["n_called"] == 0:
                if stored is not None:
                    reasons.append(f"{scheme}/{name}: metrics stored for a row without calls")
            elif stored is None or [stored[x] for x in ("tp", "fp", "tn", "fn")] != [m[x] for x in ("tp", "fp", "tn", "fn")]:
                reasons.append(f"{scheme}/{name}: stored confusion matrix differs from its saved scores")
            nc = row.get("no_call")
            if nc is not None and (nc["resistant"], nc["susceptible"]) != (m["no_call_resistant"], m["no_call_susceptible"]):
                reasons.append(f"{scheme}/{name}: stored no-calls differ from its saved scores")
            if nc is None and (m["no_call_resistant"] or m["no_call_susceptible"]):
                reasons.append(f"{scheme}/{name}: uncalled samples without a no_call record")
            if name == "lasso":
                fits = row.get("lasso_fits")
                if fits is None or len(fits) != split["folds"]:
                    reasons.append(f"{scheme}/lasso: no per-fold solver record for each of {split['folds']} folds")
                elif any(lasso_unconverged(f) for f in fits):
                    reasons.append(f"{scheme}/lasso: an unconverged two-class lasso fit")
            m["scores"] = scores
            rows[name] = m
        missing = [x for x in LEARNERS if x not in rows]
        if missing:
            reasons.append(f"{scheme}: missing learner rows {missing}")
        splits[scheme] = {"folds_of": dict(zip(ids, folds)), "rows": rows}
    return {"id": rep["id"], "antibiotic": rep["antibiotic"], "ids": ids, "truth": dict(zip(ids, truth)),
            "eval_config": rep.get("eval_config"), "fingerprint": rep.get("input_fingerprint"), "splits": splits}, reasons


def load_cohort(root: Path, slug: str):
    sp = root / "species" / slug
    genomes = [x.strip() for x in (sp / "genomes.txt").read_text().splitlines() if x.strip()]
    reports, invalid, parsed = {}, {}, {}
    results = sp / "results"
    for path in sorted(results.glob("*.json")):
        rep, reasons = load_report(path)
        drug = rep["antibiotic"]
        if drug in parsed:
            raise Stop(f"{root}: two reports for {drug}")
        parsed[drug] = rep
        if reasons:
            invalid[drug] = reasons
        else:
            reports[drug] = rep
        if not set(rep["ids"]) <= set(genomes):
            raise Stop(f"{path}: samples outside the cohort's genomes.txt")
    return {"root": str(root), "genomes": genomes, "genomes_sha256": hashlib.sha256((sp / "genomes.txt").read_bytes()).hexdigest(),
            "reports": reports, "invalid": invalid, "parsed": parsed}


# -------------------------------------------------------------- comparisons

def median_ba(rows):
    vals = [rows[x]["balanced_accuracy"] for x in LEARNERS if x in rows and rows[x]["n_called"]]
    return statistics.median(vals) if len(vals) == len(LEARNERS) else None


def public(m):
    return {k: v for k, v in m.items() if k != "scores"}


def outcome(truth, score):
    if score is None:
        return "no_call"
    return "correct" if (score >= THRESHOLD) == truth else "wrong"


def paired(a, b, scheme):
    """b against a on the isolates both hold (genome IDs), per learner."""
    sa, sb = a["splits"].get(scheme), b["splits"].get(scheme)
    if not sa or not sb or "rows" not in sa or "rows" not in sb:
        return None
    common = [g for g in a["ids"] if g in b["truth"]]
    for g in common:
        if a["truth"][g] != b["truth"][g]:
            raise Stop(f"{a['antibiotic']}: {g} has different labels in the two cohorts")
    ia = {g: i for i, g in enumerate(a["ids"])}
    ib = {g: i for i, g in enumerate(b["ids"])}
    out = {}
    for learner in LEARNERS:
        ra, rb = sa["rows"][learner], sb["rows"][learner]
        moves, deltas = {}, []
        for g in common:
            x, y = ra["scores"][ia[g]], rb["scores"][ib[g]]
            key = f"{outcome(a['truth'][g], x)}->{outcome(b['truth'][g], y)}"
            moves[key] = moves.get(key, 0) + 1
            if x is not None and y is not None:
                deltas.append(y - x)
        fixed, broken = moves.get("wrong->correct", 0), moves.get("correct->wrong", 0)
        out[learner] = {"isolates": len(common), "transitions": dict(sorted(moves.items())), "fixed": fixed, "broken": broken,
                        "mcnemar_p": mcnemar_exact(fixed, broken),
                        "score_delta_mean": statistics.fmean(deltas) if deltas else None,
                        "score_delta_median": statistics.median(deltas) if deltas else None,
                        "fold_changes_recorded": True}
    return out


def side_by_side(a, b, scheme):
    sa, sb = a["splits"].get(scheme), b["splits"].get(scheme)
    if not sa or not sb or "rows" not in sa or "rows" not in sb:
        return None
    ma, mb = median_ba(sa["rows"]), median_ba(sb["rows"])
    return {
        "learners": {x: {"a": public(sa["rows"][x]), "b": public(sb["rows"][x])} for x in LEARNERS},
        "median_ba": {"a": ma, "b": mb, "delta": None if ma is None or mb is None else mb - ma},
    }


def compare(cohorts, a_name, b_name, kind, eligible, not_interpretable, excluded):
    a, b = cohorts[a_name], cohorts[b_name]
    rows = []
    for drug in sorted(eligible):
        if {a_name, b_name} & {e["cohort"] for e in excluded.get(drug, [])}:
            continue
        ra, rb = a["reports"].get(drug), b["reports"].get(drug)
        if ra is None or rb is None:
            continue
        if ra["eval_config"] != rb["eval_config"]:
            raise Stop(f"{drug}: {a_name} and {b_name} used different eval_config")
        row = {"drug": drug, "before_after_interpretable": drug not in not_interpretable}
        if drug in not_interpretable:
            row["not_interpretable_reason"] = not_interpretable[drug]
        for scheme in SCHEMES:
            row[scheme] = {"side_by_side": side_by_side(ra, rb, scheme)}
            if kind == "paired":
                row[scheme]["paired"] = paired(ra, rb, scheme)
        rows.append(row)
    if kind == "paired":
        for scheme in SCHEMES:
            for learner in LEARNERS:
                ps = {r["drug"]: r[scheme]["paired"][learner]["mcnemar_p"] for r in rows
                      if r["before_after_interpretable"] and r[scheme].get("paired")}
                for drug, p in holm(ps).items():
                    next(r for r in rows if r["drug"] == drug)[scheme]["paired"][learner]["mcnemar_p_holm"] = p
    return {"a": a_name, "b": b_name, "kind": kind, "rows": rows}


def class_counts(full_report, genomes):
    labels = [t for g, t in full_report["truth"].items() if g in genomes]
    return sum(labels), len(labels) - sum(labels)


def derive_not_interpretable(cohorts):
    added = set(cohorts["added"]["genomes"])
    out = {}
    for drug, rep in cohorts["full"]["parsed"].items():
        r, s = class_counts(rep, added)
        if r < MIN_CLASS or s < MIN_CLASS:
            out[drug] = f"{r} resistant and {s} susceptible isolates among the added genomes (fewer than {MIN_CLASS} of a class)"
    return out


def derive_eligibility(cohorts):
    """{cohort: set of drugs} that pass the class gate on its genomes."""
    out = {}
    for name, c in cohorts.items():
        genomes = set(c["genomes"])
        out[name] = {drug for drug, rep in cohorts["full"]["parsed"].items()
                     if min(class_counts(rep, genomes)) >= MIN_CLASS}
    return out


def run(args):
    cohorts = {}
    for spec in args.cohort:
        name, _, root = spec.partition("=")
        if name not in COHORTS or name in cohorts:
            raise Stop(f"--cohort {spec}: expected one each of {COHORTS}")
        cohorts[name] = load_cohort(Path(root), args.slug)
    if sorted(cohorts) != sorted(COHORTS):
        raise Stop(f"need --cohort for each of {COHORTS}")
    receipts = {}
    for spec in args.checker_receipt:
        name, _, path = spec.partition("=")
        p = Path(path)
        if name not in COHORTS or not p.is_file():
            raise Stop(f"--checker-receipt {spec}: missing")
        receipts[name] = {"path": str(p), "sha256": sha256(p)}
    if sorted(receipts) != sorted(COHORTS):
        raise Stop("the reviewed revision-8 checker must have passed on every cohort (--checker-receipt per cohort)")
    screen_path = Path(args.duplicates)
    screen = json.loads(screen_path.read_text())
    screen_sha = sha256(screen_path)
    if screen.get("cross_cohort_exact", 1) and args.accept_reviewed_duplicates != screen_sha:
        raise Stop(f"the duplicate screen reports {screen.get('cross_cohort_exact')} exact cross-cohort duplicates; stop for a parent decision")

    prefix, full, added, matched = (cohorts[x] for x in ("prefix", "full", "added", "matched2k"))
    if prefix["genomes"] != full["genomes"][:len(prefix["genomes"])]:
        raise Stop("the prefix cohort is not the full cohort's prefix")
    if added["genomes"] != full["genomes"][len(prefix["genomes"]):]:
        raise Stop("the added cohort is not the rest of the full cohort")
    if not set(matched["genomes"]) <= set(full["genomes"]) or len(matched["genomes"]) != len(prefix["genomes"]):
        raise Stop("matched2k is not a same-size subset of the full cohort")

    not_interpretable = derive_not_interpretable(cohorts)
    if set(not_interpretable) != set(EXPECTED_SPARSE_ADDED):
        raise Stop(f"added-cohort class counts disagree with the frozen evidence: derived {sorted(not_interpretable)}, expected {sorted(EXPECTED_SPARSE_ADDED)}")
    eligible = derive_eligibility(cohorts)
    if eligible["added"] != EXPECTED_ADDED_ELIGIBLE:
        raise Stop(f"added-cohort eligibility disagrees with the frozen 14: {sorted(eligible['added'] ^ EXPECTED_ADDED_ELIGIBLE)}")
    excluded, not_eligible = {}, {}
    for name, c in cohorts.items():
        for drug in sorted(set(c["parsed"]) - eligible[name]):
            raise Stop(f"{name} has a report for {drug}, which fails the {MIN_CLASS}-per-class gate on its genomes")
        not_eligible[name] = sorted(set(cohorts["full"]["parsed"]) - eligible[name])
        for drug in sorted(eligible[name]):
            if drug in c["invalid"]:
                excluded.setdefault(drug, []).append({"cohort": name, "reasons": c["invalid"][drug]})
            elif drug not in c["reports"]:
                excluded.setdefault(drug, []).append({"cohort": name, "reasons": ["no report although eligible"]})

    both = lambda x, y: eligible[x] & eligible[y]
    comparisons = {
        "full_vs_prefix": compare(cohorts, "prefix", "full", "paired", both("prefix", "full"), not_interpretable, excluded),
        "added_vs_prefix": compare(cohorts, "prefix", "added", "side_by_side", both("prefix", "added"), not_interpretable, excluded),
        "matched2k_vs_prefix": compare(cohorts, "prefix", "matched2k", "side_by_side", both("prefix", "matched2k"), not_interpretable, excluded),
        "full_vs_matched2k": compare(cohorts, "matched2k", "full", "paired", both("matched2k", "full"), not_interpretable, excluded),
    }
    # Size-matched decomposition of the primary split's median BA change:
    # total (full - prefix) = composition (matched2k - prefix) + size (full - matched2k).
    decomposition = []
    for drug in sorted(both("prefix", "full") & eligible["matched2k"] - set(excluded)):
        meds = {}
        for name in ("prefix", "matched2k", "full"):
            split = cohorts[name]["reports"][drug]["splits"].get(SCHEMES[0], {})
            meds[name] = median_ba(split["rows"]) if "rows" in split else None
        if None in meds.values():
            continue
        decomposition.append({"drug": drug, "before_after_interpretable": drug not in not_interpretable,
                              "total": meds["full"] - meds["prefix"], "composition": meds["matched2k"] - meds["prefix"],
                              "size": meds["full"] - meds["matched2k"], "median_ba": meds})
    return {
        "inputs": {name: {k: c[k] for k in ("root", "genomes_sha256")} | {"valid_reports": sorted(c["reports"]), "invalid_reports": c["invalid"]}
                   for name, c in cohorts.items()},
        "checker_receipts": receipts,
        "duplicate_screen": {"path": str(screen_path), "sha256": screen_sha, "cross_cohort_exact": screen.get("cross_cohort_exact"),
                             "cross_cohort_near": screen.get("cross_cohort_near")},
        "primary_split": SCHEMES[0],
        "threshold": THRESHOLD,
        "not_interpretable_for_before_after": not_interpretable,
        "not_eligible": not_eligible,
        "excluded": excluded,
        "comparisons": comparisons,
        "size_matched_decomposition": decomposition,
        "scope": "Descriptive internal cross-validation comparison of matched cohorts; no external, prospective, clinical or causal validation. McNemar p-values are descriptive, Holm-adjusted within learner and split across interpretable drugs.",
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--cohort", action="append", default=[], required=True)
    ap.add_argument("--checker-receipt", action="append", default=[])
    ap.add_argument("--duplicates", required=True)
    ap.add_argument("--accept-reviewed-duplicates")
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--slug", default="mycobacterium_tuberculosis")
    args = ap.parse_args()
    try:
        out = run(args)
    except (Stop, KeyError, ValueError, OSError) as err:
        print(f"compare_cohorts: STOP: {err}", file=sys.stderr)
        return 2
    with args.out.open("x") as f:
        json.dump(out, f, indent=1)
        f.write("\n")
    if out["excluded"]:
        print(f"compare_cohorts: {len(out['excluded'])} drugs excluded; stop for review before any before/after claim", file=sys.stderr)
        return 3
    return 0


if __name__ == "__main__":
    sys.exit(main())
