#!/usr/bin/env python3
"""Step two of PLAN.md (sections 5 to 7): score the arms of one dataset.

Inputs, for one dataset d:
  --released   the released report (rev14-274-20261112/datasets/<d>.json)
  --with       this lane's refit of the unchanged cohort (arm "with")
  --without    the refit with the chosen study's genomes removed (arm "without")
  --removal    the removed genome ids (results/removal/<d>.ids)
  --placebo    zero or more placebo refits, each as <removal-ids>=<report.json>
  --out        output directory (results/scores/<d>/)

Writes:
  acceptance.txt   arm "with" against the release; arm "without" against arm "with" (PLAN.md section 5)
  views.tsv        V1, V2, V3 and the study alone, per split and row (five learners and the known-gene rule)
  effects.tsv      total = V3 - V1, composition = V2 - V1, training = V3 - V2, with bootstrap intervals
  placebo.tsv      the total and training effects of each placebo arm (point values)
  caveat.tsv       the rule of PLAN.md section 7, row by row
  summary.json     the same, for the report

Run: python3 -I scripts/score_arms.py --released ... --with ... --without ... --removal ... --out ...
"""
import argparse
import collections
import json
import math
import operator
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from stats_util import clopper_pearson  # noqa: E402

EPS = 1e-4  # scripts/ablation/stats.py EPS: scores clipped to [1e-4, 1 - 1e-4] for log loss
CUTOFF = 0.5  # a score of 0.5 or more is a resistant call (checked against the released counts)
LEARNERS = ["scm", "cart", "rf", "gbt", "lasso", "catalog"]
MEASURES = ["recall", "specificity", "balanced_accuracy", "log_loss"]
RULE_MEASURES = ["balanced_accuracy", "recall", "specificity"]  # PLAN.md section 7; log loss does not enter
N_BOOT = 2000
SEED_CLUSTER = 20261114  # PLAN.md section 6
SEED_SAMPLE = 20261115   # PLAN-ADDENDUM-20261115.md A7
THRESHOLD = 0.05         # PLAN.md section 7: an adapted research screen


def load(path):
    with open(path) as fh:
        return json.load(fh)


def splits_of(rep):
    """scheme -> {learner: scores}, or scheme -> None when the split is unavailable."""
    out = {}
    for s in rep["splits"]:
        if s.get("unavailable_reason"):
            out[s["scheme"]] = None
            continue
        out[s["scheme"]] = {l["learner"]: l["scores"] for l in s["learners"]}
    return out


def metrics_of(rep):
    out = {}
    for s in rep["splits"]:
        if s.get("unavailable_reason"):
            continue
        for l in s["learners"]:
            out[(s["scheme"], l["learner"])] = l.get("metrics")
    return out


def sample_ll(y, p):
    p = min(max(p, EPS), 1 - EPS)
    return -math.log(p) if y else -math.log1p(-p)


def measures(nr, ns, tp, tn, ll, n):
    """recall, specificity, balanced accuracy, log loss from summed counts; None where undefined."""
    rec = tp / nr if nr > 0 else None
    spe = tn / ns if ns > 0 else None
    ba = (rec + spe) / 2 if rec is not None and spe is not None else None
    lls = ll / n if n > 0 and nr > 0 and ns > 0 else None
    return {"recall": rec, "specificity": spe, "balanced_accuracy": ba, "log_loss": lls}


def quantile(sorted_vals, q):
    """Linear interpolation between order statistics (numpy's default)."""
    if not sorted_vals:
        return None
    pos = q * (len(sorted_vals) - 1)
    lo = math.floor(pos)
    hi = min(lo + 1, len(sorted_vals) - 1)
    return sorted_vals[lo] + (sorted_vals[hi] - sorted_vals[lo]) * (pos - lo)


def fmt(x, nd=4):
    if x is None:
        return "NA"
    if isinstance(x, float) and x != x:
        return "NA"
    return f"{x:.{nd}f}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--released", required=True)
    ap.add_argument("--with", dest="with_", required=True)
    ap.add_argument("--without", required=True)
    ap.add_argument("--removal", required=True)
    ap.add_argument("--placebo", nargs="*", default=[])
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    acc = []

    rel, wit, wo = load(a.released), load(a.with_), load(a.without)
    removed = set(open(a.removal).read().split())

    # ---- Acceptance: arm "with" against the release (PLAN.md section 5).
    rel_s, wit_s = rel["samples"], wit["samples"]
    same_ids = [s["id"] for s in rel_s] == [s["id"] for s in wit_s]
    acc.append(f"with vs release: sample ids equal and in the same order: {same_ids} "
               f"(release {len(rel_s)}, refit {len(wit_s)})")
    for field in ("resistant", "catalog_resistant", "lineage"):
        eq = same_ids and all(x[field] == y[field] for x, y in zip(rel_s, wit_s))
        acc.append(f"with vs release: samples[].{field} equal: {eq}")
    acc.append(f"with vs release: input_fingerprint equal: {rel['input_fingerprint'] == wit['input_fingerprint']} "
               f"({rel['input_fingerprint']} / {wit['input_fingerprint']})")
    acc.append(f"with vs release: eval_config equal: {rel['eval_config'] == wit['eval_config']}")
    rel_sp, wit_sp, wo_sp = splits_of(rel), splits_of(wit), splits_of(wo)
    scores_equal = True
    for scheme in sorted(rel_sp):
        for lr in LEARNERS:
            r = (rel_sp[scheme] or {}).get(lr)
            w = (wit_sp.get(scheme) or {}).get(lr)
            if r is None or w is None or len(r) != len(w):
                acc.append(f"with vs release: {scheme}/{lr}: scores missing or of another length")
                scores_equal = False
                continue
            diff = max((abs(x - y) for x, y in zip(r, w)), default=0.0)
            n_diff = sum(1 for x, y in zip(r, w) if x != y)
            if n_diff:
                scores_equal = False
            acc.append(f"with vs release: {scheme}/{lr}: scores differing {n_diff} of {len(r)}, max abs diff {diff:.3g}")
    acc.append(f"with vs release: every out-of-fold score equal: {scores_equal}")

    # The release's cutoff: recompute the released counts from the released scores.
    y_rel = [bool(s["resistant"]) for s in rel_s]
    cut_ok = True
    for (scheme, lr), m in sorted(metrics_of(rel).items()):
        sc = rel_sp[scheme][lr]
        tp = sum(1 for y, p in zip(y_rel, sc) if y and p >= CUTOFF)
        fn = sum(1 for y, p in zip(y_rel, sc) if y and p < CUTOFF)
        tn = sum(1 for y, p in zip(y_rel, sc) if not y and p < CUTOFF)
        fp = sum(1 for y, p in zip(y_rel, sc) if not y and p >= CUTOFF)
        if (tp, fn, tn, fp) != (m["tp"], m["fn"], m["tn"], m["fp"]):
            cut_ok = False
            acc.append(f"cutoff check FAILED {scheme}/{lr}: recomputed {(tp, fn, tn, fp)} released "
                       f"{(m['tp'], m['fn'], m['tn'], m['fp'])}")
    acc.append(f"cutoff check: released tp, fn, tn, fp recomputed from released scores at >= {CUTOFF}: {cut_ok}")
    cat_ok = all(float(s["catalog_resistant"]) == rel_sp[sch]["catalog"][i]
                 for sch in rel_sp if rel_sp[sch] for i, s in enumerate(rel_s))
    acc.append(f"catalog row scores equal samples[].catalog_resistant: {cat_ok}")

    # The arm whose scores stand for "with": the lane's own refit (identical to the release when accepted).
    base = wit
    base_s = base["samples"]
    base_sp = wit_sp

    # ---- Arm "without" against arm "with".
    wo_s = wo["samples"]
    expect = [s["id"] for s in base_s if s["id"] not in removed]
    acc.append(f"without vs with: sample ids = with minus the removed genomes, same order: "
               f"{[s['id'] for s in wo_s] == expect} (without {len(wo_s)}, expected {len(expect)}, "
               f"removed listed {len(removed)}, removed present in with {sum(1 for s in base_s if s['id'] in removed)})")
    by_id = {s["id"]: s for s in base_s}
    for field in ("resistant", "catalog_resistant", "lineage"):
        eq = all(by_id[s["id"]][field] == s[field] for s in wo_s if s["id"] in by_id)
        acc.append(f"without vs with: samples[].{field} equal on every remaining sample: {eq}")
    acc.append(f"without vs with: eval_config equal: {wo['eval_config'] == base['eval_config']}")
    for k in ("samples", "resistant", "susceptible", "min_support", "lineages", "k", "features", "patterns"):
        acc.append(f"counts.{k}: with {base['counts'].get(k)} without {wo['counts'].get(k)}")

    # ---- Per-sample table.
    n = len(base_s)
    y = [bool(s["resistant"]) for s in base_s]
    ids = [s["id"] for s in base_s]
    lineage = [s["lineage"] for s in base_s]
    inside = [g in removed for g in ids]
    wo_pos = {s["id"]: i for i, s in enumerate(wo_s)}
    schemes = [s["scheme"] for s in base["splits"]]
    avail = [sch for sch in schemes if base_sp.get(sch) and wo_sp.get(sch)]
    for sch in schemes:
        if sch not in avail:
            acc.append(f"split {sch}: unavailable in an arm (with: {base_sp.get(sch) is not None}, "
                       f"without: {wo_sp.get(sch) is not None}); enters no rule")

    # Cells: (scheme, learner, view) -> per-sample (in_view, score); views V1, V2, V3, VT.
    def cell_vectors(sch, lr):
        sw = base_sp[sch][lr]
        so = wo_sp[sch][lr]
        v = {}
        v["V1"] = [(True, sw[i]) for i in range(n)]
        v["V2"] = [(not inside[i], sw[i]) for i in range(n)]
        v["V3"] = [(not inside[i], so[wo_pos[ids[i]]] if not inside[i] else 0.0) for i in range(n)]
        v["VT"] = [(inside[i], sw[i]) for i in range(n)]
        return v

    # Units for the bootstrap: lineage clusters, or single samples.
    def unit_index(kind):
        if kind == "cluster":
            keys = sorted(set(lineage))
            pos = {k: j for j, k in enumerate(keys)}
            return [pos[lineage[i]] for i in range(n)], len(keys)
        return list(range(n)), n

    # Per unit sums for each cell: nr, ns, tp, tn, ll, nn.
    def unit_sums(units, n_units, vec):
        nr = [0] * n_units
        ns = [0] * n_units
        tp = [0] * n_units
        tn = [0] * n_units
        ll = [0.0] * n_units
        nn = [0] * n_units
        for i, (keep, p) in enumerate(vec):
            if not keep:
                continue
            u = units[i]
            nn[u] += 1
            ll[u] += sample_ll(y[i], p)
            if y[i]:
                nr[u] += 1
                tp[u] += p >= CUTOFF
            else:
                ns[u] += 1
                tn[u] += p < CUTOFF
        return nr, ns, tp, tn, ll, nn

    def point(vec):
        nr = sum(1 for i, (k, _) in enumerate(vec) if k and y[i])
        ns = sum(1 for i, (k, _) in enumerate(vec) if k and not y[i])
        tp = sum(1 for i, (k, p) in enumerate(vec) if k and y[i] and p >= CUTOFF)
        tn = sum(1 for i, (k, p) in enumerate(vec) if k and not y[i] and p < CUTOFF)
        ll = math.fsum(sample_ll(y[i], p) for i, (k, p) in enumerate(vec) if k)
        return nr, ns, tp, tn, ll, nr + ns

    cells = {(sch, lr): cell_vectors(sch, lr) for sch in avail for lr in LEARNERS}
    pts = {key: {v: point(vec) for v, vec in views.items()} for key, views in cells.items()}

    boots = {}
    meta = {}
    for kind, seed in (("cluster", SEED_CLUSTER), ("sample", SEED_SAMPLE)):
        units, n_units = unit_index(kind)
        meta[kind] = {"units": n_units}
        sums = {key: {v: unit_sums(units, n_units, vec) for v, vec in views.items()} for key, views in cells.items()}
        rng = random.Random(seed)
        draws = []
        for _ in range(N_BOOT):
            w = [0] * n_units
            for _ in range(n_units):
                w[rng.randrange(n_units)] += 1
            draws.append(w)
        res = collections.defaultdict(list)  # (sch, lr, view or effect, measure) -> values
        dropped = collections.Counter()
        for w in draws:
            for key, views in sums.items():
                ms = {}
                for v, arrs in views.items():
                    tot = [sum(map(operator.mul, w, arr)) for arr in arrs]
                    ms[v] = measures(*tot)
                for m in MEASURES:
                    for v in ("V1", "V2", "V3", "VT"):
                        if ms[v][m] is not None:
                            res[key + (v, m)].append(ms[v][m])
                    for eff, (hi_v, lo_v) in (("total", ("V3", "V1")), ("composition", ("V2", "V1")),
                                              ("training", ("V3", "V2"))):
                        a_, b_ = ms[hi_v][m], ms[lo_v][m]
                        if a_ is None or b_ is None:
                            dropped[key + (eff, m)] += 1
                        else:
                            res[key + (eff, m)].append(a_ - b_)
        boots[kind] = (res, dropped)

    def ci(kind, key):
        res, dropped = boots[kind]
        vals = sorted(res.get(key, []))
        return quantile(vals, 0.025), quantile(vals, 0.975), dropped.get(key, 0), len(vals)

    few = meta["cluster"]["units"] < 20

    # ---- views.tsv
    with open(os.path.join(a.out, "views.tsv"), "w") as fh:
        fh.write("\t".join(["dataset", "split", "row", "view", "n", "reported_R", "reported_S", "tp", "fn", "tn",
                            "fp", "recall", "recall_cp_low", "recall_cp_high", "specificity", "spec_cp_low",
                            "spec_cp_high", "balanced_accuracy", "ba_cluster_low", "ba_cluster_high", "log_loss",
                            "ll_cluster_low", "ll_cluster_high", "scores_binary"]) + "\n")
        for (sch, lr), views in pts.items():
            sc = base_sp[sch][lr]
            binary = all(p in (0.0, 1.0) for p in sc)
            for v in ("V1", "V2", "V3", "VT"):
                nr, ns, tp, tn, ll, nn = views[v]
                m = measures(nr, ns, tp, tn, ll, nn)
                rlo, rhi = clopper_pearson(tp, nr) if nr else (None, None)
                slo, shi = clopper_pearson(tn, ns) if ns else (None, None)
                blo, bhi, _, _ = ci("cluster", (sch, lr, v, "balanced_accuracy"))
                llo, lhi, _, _ = ci("cluster", (sch, lr, v, "log_loss"))
                fh.write("\t".join([a.dataset, sch, lr, v, str(nn), str(nr), str(ns), str(tp), str(nr - tp),
                                    str(tn), str(ns - tn), fmt(m["recall"]), fmt(rlo), fmt(rhi),
                                    fmt(m["specificity"]), fmt(slo), fmt(shi), fmt(m["balanced_accuracy"]),
                                    fmt(blo), fmt(bhi), fmt(m["log_loss"]), fmt(llo), fmt(lhi),
                                    "yes" if binary else "no"]) + "\n")

    # ---- effects.tsv and the caveat rule.
    effects = []
    with open(os.path.join(a.out, "effects.tsv"), "w") as fh:
        fh.write("\t".join(["dataset", "split", "row", "measure", "effect", "estimate", "cluster_low",
                            "cluster_high", "cluster_dropped", "sample_low", "sample_high", "sample_dropped",
                            "cluster_excludes_0", "passes_rule"]) + "\n")
        for (sch, lr), views in pts.items():
            ms = {v: measures(*views[v]) for v in views}
            for m in MEASURES:
                for eff, (hi_v, lo_v) in (("total", ("V3", "V1")), ("composition", ("V2", "V1")),
                                          ("training", ("V3", "V2"))):
                    a_, b_ = ms[hi_v][m], ms[lo_v][m]
                    est = a_ - b_ if a_ is not None and b_ is not None else None
                    clo, chi, cdrop, _ = ci("cluster", (sch, lr, eff, m))
                    slo, shi, sdrop, _ = ci("sample", (sch, lr, eff, m))
                    excl = clo is not None and (clo > 0 or chi < 0)
                    passes = (eff == "total" and m in RULE_MEASURES and est is not None
                              and abs(est) >= THRESHOLD and excl)
                    row = {"split": sch, "row": lr, "measure": m, "effect": eff, "estimate": est,
                           "cluster_ci": [clo, chi], "cluster_dropped": cdrop, "sample_ci": [slo, shi],
                           "sample_dropped": sdrop, "cluster_excludes_0": excl, "passes": passes}
                    effects.append(row)
                    fh.write("\t".join([a.dataset, sch, lr, m, eff, fmt(est), fmt(clo), fmt(chi), str(cdrop),
                                        fmt(slo), fmt(shi), str(sdrop), "yes" if excl else "no",
                                        "yes" if passes else "no"]) + "\n")

    # ---- Placebo arms: point values of the total and training effects.
    placebo_rows = []
    for spec in a.placebo:
        ids_path, rep_path = spec.split("=", 1)
        prem = set(open(ids_path).read().split())
        prep = load(rep_path)
        psp = splits_of(prep)
        ppos = {s["id"]: i for i, s in enumerate(prep["samples"])}
        p_expect = [s["id"] for s in base_s if s["id"] not in prem]
        ok = [s["id"] for s in prep["samples"]] == p_expect
        tag = os.path.basename(ids_path)
        acc.append(f"placebo {tag}: sample ids = with minus its removal: {ok}; eval_config equal: "
                   f"{prep['eval_config'] == base['eval_config']}")
        for sch in avail:
            if not psp.get(sch):
                continue
            for lr in LEARNERS:
                sw, sp_ = base_sp[sch][lr], psp[sch][lr]
                v1 = point([(True, sw[i]) for i in range(n)])
                v2 = point([(ids[i] not in prem, sw[i]) for i in range(n)])
                v3 = point([(ids[i] not in prem, sp_[ppos[ids[i]]] if ids[i] in ppos else 0.0) for i in range(n)])
                m1, m2, m3 = measures(*v1), measures(*v2), measures(*v3)
                for m in MEASURES:
                    tot = m3[m] - m1[m] if m3[m] is not None and m1[m] is not None else None
                    trn = m3[m] - m2[m] if m3[m] is not None and m2[m] is not None else None
                    placebo_rows.append({"placebo": tag, "split": sch, "row": lr, "measure": m, "total": tot,
                                         "training": trn})
    with open(os.path.join(a.out, "placebo.tsv"), "w") as fh:
        fh.write("dataset\tplacebo\tsplit\trow\tmeasure\ttotal\ttraining\n")
        for r in placebo_rows:
            fh.write("\t".join([a.dataset, r["placebo"], r["split"], r["row"], r["measure"], fmt(r["total"]),
                                fmt(r["training"])]) + "\n")

    def spread(sch, lr, m):
        vals = [r["total"] for r in placebo_rows if (r["split"], r["row"], r["measure"]) == (sch, lr, m)
                and r["total"] is not None]
        if not vals:
            return None
        return {"n": len(vals), "min": min(vals), "max": max(vals),
                "max_abs": max(abs(v) for v in vals)}

    tested = [e for e in effects if e["effect"] == "total" and e["measure"] in RULE_MEASURES]
    passing = [e for e in tested if e["passes"]]
    with open(os.path.join(a.out, "caveat.tsv"), "w") as fh:
        fh.write("dataset\tsplit\trow\tmeasure\ttotal\tcluster_low\tcluster_high\tpasses\tplacebo_n\t"
                 "placebo_min\tplacebo_max\tplacebo_max_abs\n")
        for e in tested:
            sp_ = spread(e["split"], e["row"], e["measure"])
            fh.write("\t".join([a.dataset, e["split"], e["row"], e["measure"], fmt(e["estimate"]),
                                fmt(e["cluster_ci"][0]), fmt(e["cluster_ci"][1]), "yes" if e["passes"] else "no",
                                str(sp_["n"]) if sp_ else "0", fmt(sp_["min"]) if sp_ else "NA",
                                fmt(sp_["max"]) if sp_ else "NA", fmt(sp_["max_abs"]) if sp_ else "NA"]) + "\n")

    summary = {
        "dataset": a.dataset, "lineage_clusters": meta["cluster"]["units"], "few_clusters": few,
        "removed": len(removed), "samples_with": n, "samples_without": len(wo_s),
        "splits_scored": avail, "rows_tested": len(tested), "rows_passing": len(passing),
        "passing": [{k: e[k] for k in ("split", "row", "measure", "estimate", "cluster_ci")} |
                    {"placebo": spread(e["split"], e["row"], e["measure"])} for e in passing],
        "dependence_caveat": bool(passing), "placebo_arms": len(a.placebo),
        "bootstrap": {"replicates": N_BOOT, "cluster_seed": SEED_CLUSTER, "sample_seed": SEED_SAMPLE},
        "with_equals_release": scores_equal and same_ids,
    }
    with open(os.path.join(a.out, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=1)
    with open(os.path.join(a.out, "acceptance.txt"), "w") as fh:
        fh.write("\n".join(acc) + "\n")
    print("\n".join(acc))
    print(json.dumps({k: summary[k] for k in ("dataset", "lineage_clusters", "few_clusters", "rows_tested",
                                              "rows_passing", "with_equals_release")}))


if __name__ == "__main__":
    main()
