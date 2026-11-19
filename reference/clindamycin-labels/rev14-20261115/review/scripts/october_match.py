#!/usr/bin/env python3
"""Exploratory (decided after results were seen): search revision 14 views of staphylococcus_aureus__clindamycin for
the October record's two figures, 0.772809 and 0.870555, over rows, splits, views V1/V2/V3 and measures, within 5e-7.
Also prints the lineage-split set covering machine balanced accuracy of V1, V2, V3 for each pair at six decimals.

  python3 -I scripts/october_match.py <removal ids> <with.json>=<without.json> [more pairs ...]
"""
import json, math, sys
rem = set(open(sys.argv[1]).read().split())
def meas(pairs):
    tp = sum(1 for t, p in pairs if t and p >= .5); fn = sum(1 for t, p in pairs if t and p < .5)
    tn = sum(1 for t, p in pairs if not t and p < .5); fp = sum(1 for t, p in pairs if not t and p >= .5)
    d = (tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)
    return {"accuracy": (tp + tn) / len(pairs), "recall": tp / (tp + fn), "specificity": tn / (tn + fp),
            "balanced_accuracy": (tp / (tp + fn) + tn / (tn + fp)) / 2,
            "mcc": (tp * tn - fp * fn) / math.sqrt(d) if d else float("nan"), "f1": 2 * tp / (2 * tp + fp + fn)}
targets = (0.772809, 0.870555)
for pair in sys.argv[2:]:
    wp, wop = pair.split("=")
    W = json.load(open(wp)); WO = json.load(open(wop))
    ids = [s["id"] for s in W["samples"]]; y = [s["resistant"] for s in W["samples"]]
    pos = {s["id"]: i for i, s in enumerate(WO["samples"])}
    tag = wop.split("/")[-1]
    for sp in W["splits"]:
        sch = sp["scheme"]; sc = {l["learner"]: l["scores"] for l in sp["learners"]}
        so = {l["learner"]: l["scores"] for l in [x for x in WO["splits"] if x["scheme"] == sch][0]["learners"]}
        for lr in sc:
            views = {"V1": list(zip(y, sc[lr])), "V2": [(y[i], sc[lr][i]) for i in range(len(y)) if ids[i] not in rem],
                     "V3": [(y[i], so[lr][pos[ids[i]]]) for i in range(len(y)) if ids[i] not in rem]}
            ms = {v: meas(p) for v, p in views.items()}
            if sch == "lineage" and lr == "scm":
                print(f"{tag}: lineage scm balanced accuracy V1 {ms['V1']['balanced_accuracy']:.6f} "
                      f"V2 {ms['V2']['balanced_accuracy']:.6f} V3 {ms['V3']['balanced_accuracy']:.6f}")
            for v, m_ in ms.items():
                for m, x in m_.items():
                    for t in targets:
                        if abs(x - t) < 5e-7:
                            print(f"  MATCH {t}: {tag} {sch} {lr} {v} {m}")
print("pairs searched:", len(sys.argv) - 2)
