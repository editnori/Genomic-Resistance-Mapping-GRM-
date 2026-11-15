#!/usr/bin/env python3
"""Brief point 4, as fixed in ADDENDUM.md (exploratory): fold plan spread for staphylococcus_aureus__clindamycin.

Reads /home/user/scratch/revclinda/runs/c1f/{with-sa,without-sa}.seed-<s>.json for s in 42, 1..10, checks the default
seed against lane clinda's reports, and writes results/foldplans/{check.txt,spread.tsv,tests.tsv,summary.txt}.

Run: python3 -I scripts/foldplans.py
"""
import json
import os

RUNS = "/home/user/scratch/revclinda/runs/c1f"
LANE = "/home/user/scratch/clinda/reports"
REM = "/home/user/grm-ops/program-20261114/clinda/results/removal/staphylococcus_aureus__clindamycin.ids"
RULE = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results", "refit",
                    "rescore-staphylococcus_aureus", "rule.tsv")
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results", "foldplans")
D = "staphylococcus_aureus__clindamycin"
SEEDS = [42] + list(range(1, 11))
ROWS = ["scm", "cart", "rf", "gbt", "lasso", "catalog"]
MEAS = ["balanced_accuracy", "recall", "specificity"]


def scores(rep):
    return {s["scheme"]: {l["learner"]: l["scores"] for l in s["learners"]} for s in rep["splits"]}


def meas(pairs):
    tp = sum(1 for t, p in pairs if t and p >= .5)
    fn = sum(1 for t, p in pairs if t and p < .5)
    tn = sum(1 for t, p in pairs if not t and p < .5)
    fp = sum(1 for t, p in pairs if not t and p >= .5)
    r, s = tp / (tp + fn), tn / (tn + fp)
    return {"recall": r, "specificity": s, "balanced_accuracy": (r + s) / 2}


def main():
    os.makedirs(OUT, exist_ok=True)
    rem = set(open(REM).read().split())
    W = {s: json.load(open(os.path.join(RUNS, f"with-sa.seed-{s}.json"))) for s in SEEDS}
    WO = {s: json.load(open(os.path.join(RUNS, f"without-sa.seed-{s}.json"))) for s in SEEDS}
    check = []
    lw, lwo = json.load(open(os.path.join(LANE, D + ".with.json"))), json.load(open(os.path.join(LANE, D + ".without.json")))
    for name, mine, lane in (("with", W[42], lw), ("without", WO[42], lwo)):
        a, b = scores(mine), scores(lane)
        n = d = 0
        for sch in a:
            for r in a[sch]:
                n += len(a[sch][r])
                d += sum(1 for x, y in zip(a[sch][r], b[sch][r]) if x != y) + abs(len(a[sch][r]) - len(b[sch][r]))
        same_ids = [s["id"] for s in mine["samples"]] == [s["id"] for s in lane["samples"]]
        check.append(f"default seed, arm {name}: c1f36fc scores vs lane clinda's report: compared {n}, differing {d}; "
                     f"sample ids equal {same_ids}; eval_config.seed {mine['eval_config'].get('seed')}")
    for s in SEEDS:
        check.append(f"seed {s}: eval_config.seed with {W[s]['eval_config'].get('seed')} without "
                     f"{WO[s]['eval_config'].get('seed')}; n with {len(W[s]['samples'])} without {len(WO[s]['samples'])}")
    # Distinct fold plans: compare test_fold vectors across seeds.
    for arm, R in (("with", W), ("without", WO)):
        for sch in ("random", "lineage"):
            plans = {tuple(x["test_fold"]) for s in SEEDS for x in R[s]["splits"] if x["scheme"] == sch}
            check.append(f"arm {arm} split {sch}: distinct fold plans over {len(SEEDS)} seeds: {len(plans)}")
    ids = [s["id"] for s in W[42]["samples"]]
    y = [s["resistant"] for s in W[42]["samples"]]
    for s in SEEDS:
        assert [x["id"] for x in W[s]["samples"]] == ids
    vals = {}  # (split,row,view,measure) -> list over seeds
    for s in SEEDS:
        sw, so = scores(W[s]), scores(WO[s])
        pos = {x["id"]: i for i, x in enumerate(WO[s]["samples"])}
        for sch in ("random", "lineage"):
            for r in ROWS:
                v1 = meas(list(zip(y, sw[sch][r])))
                v2 = meas([(y[i], sw[sch][r][i]) for i in range(len(y)) if ids[i] not in rem])
                v3 = meas([(y[i], so[sch][r][pos[ids[i]]]) for i in range(len(y)) if ids[i] not in rem])
                for v, m_ in (("V1", v1), ("V2", v2), ("V3", v3)):
                    for m in MEAS:
                        vals.setdefault((sch, r, v, m), []).append(m_[m])
    with open(os.path.join(OUT, "spread.tsv"), "w") as fh:
        fh.write("split\trow\tmeasure\tview\tseed42\tmin\tmax\trange\tvalues_seed42_then_1_to_10\n")
        for (sch, r, v, m), xs in sorted(vals.items()):
            fh.write(f"{sch}\t{r}\t{m}\t{v}\t{xs[0]:.4f}\t{min(xs):.4f}\t{max(xs):.4f}\t{max(xs) - min(xs):.4f}\t"
                     + ",".join(f"{x:.4f}" for x in xs) + "\n")
    passing = [l.rstrip("\n").split("\t") for l in open(RULE)][1:]
    passing = [p for p in passing if p[5] == "yes" and p[3] in ("balanced_accuracy", "recall")]
    tests = []
    for p in passing:
        sch, r, m, est = p[1], p[2], p[3], float(p[4])
        V1, V2, V3 = vals[(sch, r, "V1", m)], vals[(sch, r, "V2", m)], vals[(sch, r, "V3", m)]
        up = est > 0
        t1 = (min(V3) > max(V1)) if up else (max(V3) < min(V1))
        t2 = (min(V3) > max(V2)) if up else (max(V3) < min(V2))
        per_seed = [b - a for a, b in zip(V1, V3)]
        tests.append([sch, r, m, f"{est:+.4f}", f"{max(V1) - min(V1):.4f}", f"{max(V2) - min(V2):.4f}",
                      f"{max(V3) - min(V3):.4f}", f"{min(V1):.4f}..{max(V1):.4f}", f"{min(V2):.4f}..{max(V2):.4f}",
                      f"{min(V3):.4f}..{max(V3):.4f}", f"{min(per_seed):+.4f}..{max(per_seed):+.4f}",
                      "yes" if t1 else "no", "yes" if t2 else "no"])
    with open(os.path.join(OUT, "tests.tsv"), "w") as fh:
        fh.write("split\trow\tmeasure\tlane_change\trange_V1\trange_V2\trange_V3\tV1_min..max\tV2_min..max\t"
                 "V3_min..max\tper_seed_change_min..max\ttest1_total\ttest2_training\n")
        for t in tests:
            fh.write("\t".join(t) + "\n")
    summ = []
    for sch in ("random", "lineage"):
        ts = [t for t in tests if t[0] == sch]
        summ.append(f"{sch}: lane's passing rows (balanced accuracy, recall) {len(ts)}; meet test 1 (total) "
                    f"{sum(t[11] == 'yes' for t in ts)}; meet test 2 (training) {sum(t[12] == 'yes' for t in ts)}")
    with open(os.path.join(OUT, "check.txt"), "w") as fh:
        fh.write("\n".join(check) + "\n")
    with open(os.path.join(OUT, "summary.txt"), "w") as fh:
        fh.write("\n".join(summ) + "\n")
    print("\n".join(check + summ))


if __name__ == "__main__":
    main()
