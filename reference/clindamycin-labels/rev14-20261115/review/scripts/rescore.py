#!/usr/bin/env python3
"""Brief point 3.4: the reviewer's own recompute of lane clinda's refit table from saved scores.

Written without reuse of the lane's code. For one dataset:

  python3 -I scripts/rescore.py --dataset D --with W.json --without WO.json --removal IDS --seed S [--seed S2 ...]
         [--placebo IDS=REPORT ...] --out DIR

Views, per split and row (five learners and the known-gene rule, "catalog"):
  V1  arm "with", all samples;  V2  arm "with", samples outside the study;  V3  arm "without", the same samples.
Measures at the 0.5 cutoff: recall, specificity, balanced accuracy, and log loss with scores clipped to
[1e-4, 1 - 1e-4]. Changes: total V3-V1, composition V2-V1, training V3-V2.
Intervals: percentile bootstrap over lineage clusters (samples[].lineage of arm "with"), 2,000 replicates, one run per
seed given. Rule of PLAN.md section 7: |total| >= 0.05 and the 95% interval excludes 0, for balanced accuracy, recall
and specificity, per split and row.

Writes DIR/table.tsv (every row, measure and change with one interval per seed), DIR/rule.tsv and DIR/placebo.tsv.
"""
import argparse
import json
import math
import os
import random

CUT = 0.5
EPS = 1e-4
ROWS = ["scm", "cart", "rf", "gbt", "lasso", "catalog"]
RULE = ["balanced_accuracy", "recall", "specificity"]
ALLM = ["recall", "specificity", "balanced_accuracy", "log_loss"]
B = 2000


def ll1(y, p):
    p = min(max(p, EPS), 1 - EPS)
    return -math.log(p) if y else -math.log(1 - p)


def stats(nr, ns, tp, tn, ll, n):
    if nr <= 0 or ns <= 0:
        return None
    rec, spe = tp / nr, tn / ns
    return {"recall": rec, "specificity": spe, "balanced_accuracy": (rec + spe) / 2, "log_loss": ll / n}


def pct(v, q):
    v = sorted(v)
    if not v:
        return None
    h = (len(v) - 1) * q
    lo = int(math.floor(h))
    hi = min(lo + 1, len(v) - 1)
    return v[lo] + (v[hi] - v[lo]) * (h - lo)


def scores(rep):
    return {s["scheme"]: ({l["learner"]: l["scores"] for l in s["learners"]} if not s.get("unavailable_reason")
                          else None) for s in rep["splits"]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--with", dest="w", required=True)
    ap.add_argument("--without", required=True)
    ap.add_argument("--removal", required=True)
    ap.add_argument("--seed", type=int, action="append", required=True)
    ap.add_argument("--placebo", nargs="*", default=[])
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    W, WO = json.load(open(a.w)), json.load(open(a.without))
    rem = set(open(a.removal).read().split())
    S = W["samples"]
    n = len(S)
    y = [bool(s["resistant"]) for s in S]
    ids = [s["id"] for s in S]
    out_ = [g not in rem for g in ids]
    pos_wo = {s["id"]: i for i, s in enumerate(WO["samples"])}
    assert all((g in pos_wo) == o for g, o in zip(ids, out_)), "without arm is not with minus the removal"
    lin = sorted({s["lineage"] for s in S})
    lpos = {k: j for j, k in enumerate(lin)}
    unit = [lpos[s["lineage"]] for s in S]
    K = len(lin)
    sw, so = scores(W), scores(WO)
    splits = [s["scheme"] for s in W["splits"] if sw.get(s["scheme"]) and so.get(s["scheme"])]

    def vec(view, sch, r):
        if view == "V1":
            return [(i, sw[sch][r][i]) for i in range(n)]
        if view == "V2":
            return [(i, sw[sch][r][i]) for i in range(n) if out_[i]]
        return [(i, so[sch][r][pos_wo[ids[i]]]) for i in range(n) if out_[i]]

    def per_unit(v):
        arr = [[0, 0, 0, 0, 0.0, 0] for _ in range(K)]
        for i, p in v:
            u = arr[unit[i]]
            u[5] += 1
            u[4] += ll1(y[i], p)
            if y[i]:
                u[0] += 1
                u[2] += p >= CUT
            else:
                u[1] += 1
                u[3] += p < CUT
        return arr

    def total(arr, w=None):
        t = [0, 0, 0, 0, 0.0, 0]
        for j, u in enumerate(arr):
            c = 1 if w is None else w[j]
            if c:
                for k in range(6):
                    t[k] += c * u[k]
        return stats(*t)

    cells = {}
    for sch in splits:
        for r in ROWS:
            cells[(sch, r)] = {v: per_unit(vec(v, sch, r)) for v in ("V1", "V2", "V3")}
    point = {k: {v: total(arr) for v, arr in c.items()} for k, c in cells.items()}
    ci = {}
    for seed in a.seed:
        rng = random.Random(seed)
        acc = {}
        for _ in range(B):
            w = [0] * K
            for _ in range(K):
                w[rng.randrange(K)] += 1
            for k, c in cells.items():
                st = {v: total(arr, w) for v, arr in c.items()}
                for m in ALLM:
                    for eff, (h, l) in (("total", ("V3", "V1")), ("composition", ("V2", "V1")),
                                        ("training", ("V3", "V2"))):
                        if st[h] is not None and st[l] is not None:
                            acc.setdefault(k + (m, eff), []).append(st[h][m] - st[l][m])
        ci[seed] = {key: (pct(v, 0.025), pct(v, 0.975), len(v)) for key, v in acc.items()}

    rows = []
    rule = []
    for (sch, r), pt in point.items():
        binary = all(p in (0.0, 1.0) for p in sw[sch][r])
        for m in ALLM:
            for eff, (h, l) in (("total", ("V3", "V1")), ("composition", ("V2", "V1")), ("training", ("V3", "V2"))):
                est = pt[h][m] - pt[l][m]
                line = [a.dataset, sch, r, m, eff, f"{pt['V1'][m]:.4f}", f"{pt['V2'][m]:.4f}", f"{pt['V3'][m]:.4f}",
                        f"{est:+.4f}", "binary" if binary else "graded"]
                for seed in a.seed:
                    lo, hi, nb = ci[seed][(sch, r, m, eff)]
                    line += [f"{lo:+.4f}", f"{hi:+.4f}", str(nb)]
                rows.append(line)
                if eff == "total" and m in RULE:
                    passes = []
                    for seed in a.seed:
                        lo, hi, _ = ci[seed][(sch, r, m, eff)]
                        passes.append(abs(est) >= 0.05 and (lo > 0 or hi < 0))
                    rule.append([a.dataset, sch, r, m, f"{est:+.4f}"] + ["yes" if p else "no" for p in passes])
    head = ["dataset", "split", "row", "measure", "change", "V1", "V2", "V3", "estimate", "scores"]
    for seed in a.seed:
        head += [f"low_s{seed}", f"high_s{seed}", f"reps_s{seed}"]
    with open(os.path.join(a.out, "table.tsv"), "w") as fh:
        fh.write("\t".join(head) + "\n")
        for line in rows:
            fh.write("\t".join(line) + "\n")
    with open(os.path.join(a.out, "rule.tsv"), "w") as fh:
        fh.write("\t".join(["dataset", "split", "row", "measure", "total"] + [f"passes_s{s}" for s in a.seed]) + "\n")
        for line in rule:
            fh.write("\t".join(line) + "\n")
    # Placebo arms: point values of total change for every row.
    prow = []
    for spec in a.placebo:
        pid, prep = spec.split("=", 1)
        prem = set(open(pid).read().split())
        P = json.load(open(prep))
        sp = scores(P)
        ppos = {s["id"]: i for i, s in enumerate(P["samples"])}
        assert [s["id"] for s in P["samples"]] == [g for g in ids if g not in prem]
        for sch in splits:
            for r in ROWS:
                v1 = total(per_unit([(i, sw[sch][r][i]) for i in range(n)]))
                v2 = total(per_unit([(i, sw[sch][r][i]) for i in range(n) if ids[i] not in prem]))
                v3 = total(per_unit([(i, sp[sch][r][ppos[ids[i]]]) for i in range(n) if ids[i] not in prem]))
                for m in ALLM:
                    prow.append([a.dataset, os.path.basename(pid), sch, r, m, f"{v3[m] - v1[m]:+.4f}",
                                 f"{v3[m] - v2[m]:+.4f}"])
    with open(os.path.join(a.out, "placebo.tsv"), "w") as fh:
        fh.write("dataset\tplacebo\tsplit\trow\tmeasure\ttotal\ttraining\n")
        for line in prow:
            fh.write("\t".join(line) + "\n")
    npass = {s: sum(1 for x in rule if x[5 + a.seed.index(s)] == "yes") for s in a.seed}
    print(f"{a.dataset}: clusters={K} samples with={n} without={len(WO['samples'])} rows tested={len(rule)} "
          f"passing per seed={npass}")


if __name__ == "__main__":
    main()
