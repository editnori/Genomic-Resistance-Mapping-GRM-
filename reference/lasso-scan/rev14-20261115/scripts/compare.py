"""Compares the revision 14 scan's lasso fits with the released reports (PLAN.md,
"The release's fit against the scan's").

For each dataset: pass 1 (timed, unmodified binary) against pass 2 (instrumented
binary) fit by fit; the scan's audits, weight counts and final model against
the report; and each split's out-of-fold scores, assembled from pass 2's dumps,
against the report's lasso scores, with the released measures recomputed on
both sides.

usage: python3 -I compare.py --scan DIR --scan2 DIR --dump DIR --release DIR --out FILE
"""
import argparse
import json
import math
import sys
from pathlib import Path


def events(path):
    return [json.loads(l) for l in Path(path).read_text().splitlines() if l.strip()]


def fits(evs):
    return [e for e in evs if e.get("event") == "fit"]


def metrics(truth, scores):
    tp = sum(1 for t, s in zip(truth, scores) if t and s >= 0.5)
    fn = sum(1 for t, s in zip(truth, scores) if t and s < 0.5)
    fp = sum(1 for t, s in zip(truth, scores) if not t and s >= 0.5)
    tn = sum(1 for t, s in zip(truth, scores) if not t and s < 0.5)
    pos, neg = tp + fn, tn + fp
    rec = tp / pos if pos else 0.0
    spec = tn / neg if neg else 0.0
    eps = 1e-15
    ll = -math.fsum(math.log(min(max(s, eps), 1 - eps)) if t else math.log(1 - min(max(s, eps), 1 - eps))
                    for t, s in zip(truth, scores)) / len(truth)
    # AUROC by ranks, ties shared (Mann-Whitney).
    order = sorted(range(len(scores)), key=lambda i: scores[i])
    ranks = [0.0] * len(scores)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and scores[order[j + 1]] == scores[order[i]]:
            j += 1
        r = (i + j) / 2 + 1
        for k in range(i, j + 1):
            ranks[order[k]] = r
        i = j + 1
    rs = math.fsum(r for r, t in zip(ranks, truth) if t)
    auroc = (rs - pos * (pos + 1) / 2) / (pos * neg) if pos and neg else None
    return {"tp": tp, "fp": fp, "tn": tn, "fn": fn, "recall": rec, "log_loss": ll,
            "balanced_accuracy": (rec + spec) / 2, "vme": fn / pos if pos else 0.0, "me": fp / neg if neg else 0.0,
            "auroc": auroc}


def strip(ev):
    return {k: v for k, v in ev.items() if k != "seconds"}


def dataset(args, did):
    rel = json.loads(Path(args.release, "datasets", did + ".json").read_text())
    truth = [bool(s["resistant"]) for s in rel["samples"]]
    p1 = fits(events(Path(args.scan, did + ".events.jsonl")))
    p2 = fits(events(Path(args.scan2, did + ".events.jsonl")))
    out = {"id": did, "pass1_fits": len(p1), "pass2_fits": len(p2)}
    out["pass2_equals_pass1"] = len(p1) == len(p2) and all(strip(a) == strip(b) for a, b in zip(p1, p2))
    diffs = []
    splits = [s for s in rel["splits"] if not s.get("unavailable_reason")]
    split_rows = []
    for sp in splits:
        L = next(l for l in sp["learners"] if l["learner"] == "lasso")
        mine = [f for f in p1 if f["scheme"] == sp["scheme"]]
        if len(mine) != sp["folds"]:
            diffs.append(f"{sp['scheme']}: {len(mine)} scan fits for {sp['folds']} folds")
            continue
        oof = [None] * len(truth)
        part_ok = True
        for fi, f in enumerate(mine):
            if f["outer_fold"] != fi:
                diffs.append(f"{sp['scheme']} fold order")
            if f["audit"] != L["lasso_fits"][fi]:
                diffs.append(f"{sp['scheme']} fold {fi} audit differs: scan {f['audit']} release {L['lasso_fits'][fi]}")
            w = len(f["model"]["coef"]) if f.get("model") else 0
            if float(w) != L["fold_complexity"][fi]:
                diffs.append(f"{sp['scheme']} fold {fi} weights {w} release {L['fold_complexity'][fi]}")
            d = json.loads(Path(args.dump, f"{did}__{sp['scheme']}__{fi}.json").read_text())
            train = set(d["train"])
            test = [i for i in range(d["n"]) if i not in train]
            if d["n"] != len(truth) or sorted(d["y"]) != [i for i, t in enumerate(truth) if t]:
                diffs.append(f"{sp['scheme']} fold {fi}: dump labels are not the report's sample order")
                part_ok = False
            if test != [i for i, t in enumerate(sp["test_fold"]) if t == fi]:
                diffs.append(f"{sp['scheme']} fold {fi}: test samples are not the report's test_fold")
                part_ok = False
            if d["scores"] is None:
                prev = sum(1 for i in train if truth[i]) / len(train)
                for i in test:
                    oof[i] = prev
            else:
                for i in test:
                    oof[i] = d["scores"][i]
        if not part_ok or any(s is None for s in oof):
            diffs.append(f"{sp['scheme']}: out-of-fold scores incomplete")
            continue
        rs = L["scores"]
        m_scan, m_rel = metrics(truth, oof), metrics(truth, rs)
        # Check my metric code against the report's own fields.
        rm = L["metrics"]
        own_check = (m_rel["tp"], m_rel["fp"], m_rel["tn"], m_rel["fn"]) == (rm["tp"], rm["fp"], rm["tn"], rm["fn"]) and \
            abs(m_rel["recall"] - rm["sensitivity"]) < 1e-12 and abs(m_rel["balanced_accuracy"] - rm["balanced_accuracy"]) < 1e-12 and \
            (rm["auroc"] is None or abs(m_rel["auroc"] - rm["auroc"]) < 1e-9)
        calls_changed = sum(1 for a, b in zip(oof, rs) if (a >= 0.5) != (b >= 0.5))
        split_rows.append({
            "scheme": sp["scheme"], "n": len(truth),
            "scores_identical": sum(1 for a, b in zip(oof, rs) if a == b),
            "max_abs_score_diff": max(abs(a - b) for a, b in zip(oof, rs)),
            "calls_changed": calls_changed,
            "release": {k: m_rel[k] for k in ("recall", "log_loss", "balanced_accuracy", "vme", "me", "auroc")},
            "scan": {k: m_scan[k] for k in ("recall", "log_loss", "balanced_accuracy", "vme", "me", "auroc")},
            "delta": {k: (m_scan[k] - m_rel[k]) if m_scan[k] is not None and m_rel[k] is not None else None
                      for k in ("recall", "log_loss", "balanced_accuracy", "vme", "me", "auroc")},
            "metric_code_matches_report": own_check,
        })
    fin = [f for f in p1 if f["scheme"] == "final"]
    if len(fin) != 1:
        diffs.append("no single final fit")
    else:
        if fin[0]["audit"] != rel["models"]["lasso_fit"]:
            diffs.append(f"final audit differs: scan {fin[0]['audit']} release {rel['models']['lasso_fit']}")
        rm = rel["models"]["lasso"]
        sm = fin[0]["model"]
        if (rm is None) != (sm is None) or (sm is not None and (sm["intercept"] != rm["intercept"] or sm["lambda"] != rm["lambda"]
                                                              or [list(c) for c in sm["coef"]] != [list(c) for c in rm["coef"]])):
            diffs.append("final model differs from models.lasso")
        out["final_weights"] = len(sm["coef"]) if sm else 0
    out["splits"] = split_rows
    out["differences"] = diffs
    out["same_as_release"] = not diffs and all(r["calls_changed"] == 0 and r["max_abs_score_diff"] == 0 for r in split_rows)
    return out


def main(argv):
    ap = argparse.ArgumentParser()
    for a in ("--scan", "--scan2", "--dump", "--release", "--out"):
        ap.add_argument(a, type=Path, required=True)
    args = ap.parse_args(argv)
    summary = json.loads(Path(args.scan, "summary.json").read_text())
    rows = [dataset(args, d["id"]) for d in summary["datasets"]]
    args.out.write_text(json.dumps(rows, indent=1) + "\n")
    for r in rows:
        print(r["id"], "pass2=pass1" if r["pass2_equals_pass1"] else "PASS2 DIFFERS", "same_as_release" if r["same_as_release"] else "DIFFERS",
              len(r["differences"]), [(s["scheme"], s["calls_changed"], s["max_abs_score_diff"]) for s in r["splits"]])


if __name__ == "__main__":
    main(sys.argv[1:])
