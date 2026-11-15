"""revlasso: recount the lassoscan scan from its saved events, and compare its
out-of-fold lasso scores with the revision 14 release, written independently
of the lane's scan14.py, compare.py and tables.py.

usage: python3 -I recount.py EVENTS_DIR DUMP_DIR RELEASE_DIR LIST OUT_JSON

EVENTS_DIR holds <id>.events.jsonl(.gz) of pass 1; DUMP_DIR the pass 2 dumps;
RELEASE_DIR the unpacked release; LIST the registered scan-list.json.
Standard library only.
"""
import gzip
import json
import math
import sys
from pathlib import Path


def read_events(d, did):
    p = Path(d, did + ".events.jsonl")
    text = p.read_text() if p.exists() else gzip.open(str(p) + ".gz", "rt").read()
    return [json.loads(x) for x in text.splitlines() if x.strip()]


def audit_bad(a):
    # LassoFitAudit::has_unconverged, as scan.py's unconverged() states it.
    if a is None:
        return True
    if a["path_unconverged"] != 0:
        return True
    for k in ("selected_convergence", "refit_convergence"):
        if a[k] is not None and a[k]["converged"] is not True:
            return True
    return False


def auroc(y, s):
    # Probability a random positive outranks a random negative, ties half,
    # by grouping equal scores.
    pairs = sorted(zip(s, y))
    neg_below = 0
    total = 0.0
    i = 0
    npos = sum(y)
    nneg = len(y) - npos
    if npos == 0 or nneg == 0:
        return None
    while i < len(pairs):
        j = i
        while j < len(pairs) and pairs[j][0] == pairs[i][0]:
            j += 1
        p = sum(1 for k in range(i, j) if pairs[k][1])
        n = (j - i) - p
        total += p * neg_below + 0.5 * p * n
        neg_below += n
        i = j
    return total / (npos * nneg)


def measures(y, s):
    tp = fn = fp = tn = 0
    for t, v in zip(y, s):
        c = v >= 0.5
        if t and c:
            tp += 1
        elif t:
            fn += 1
        elif c:
            fp += 1
        else:
            tn += 1
    rec = tp / (tp + fn)
    spec = tn / (tn + fp)
    lo, hi = 1e-15, 1 - 1e-15
    ll = 0.0
    terms = []
    for t, v in zip(y, s):
        q = min(max(v, lo), hi)
        terms.append(-math.log(q) if t else -math.log(1 - q))
    ll = math.fsum(terms) / len(y)
    return {"tp": tp, "fn": fn, "fp": fp, "tn": tn, "recall": rec, "log_loss": ll,
            "balanced_accuracy": 0.5 * (rec + spec), "auroc": auroc(y, s)}


def main(evdir, dumpdir, reldir, listpath, outpath):
    entries = json.loads(Path(listpath).read_text())["datasets"]
    out = {"datasets": []}
    tot = {"datasets": 0, "converged": 0, "fits": 0, "points": 0, "path_points": 0, "refit_points": 0,
           "unconverged_points": 0, "short_solve_points": 0, "max_irls": 0, "irls_total": 0,
           "scores_compared": 0, "scores_bit_identical": 0, "scores_from_binary": 0,
           "scores_from_binary_identical": 0, "null_fold_scores": 0, "calls_changed": 0,
           "audit_mismatch": 0, "weights_mismatch": 0, "final_mismatch": 0, "trace_audit_mismatch": 0}
    null_folds = []
    for e in entries:
        did = e["id"]
        rel = json.loads(Path(reldir, "datasets", did + ".json").read_text())
        evs = read_events(evdir, did)
        inp = [x for x in evs if x["event"] == "inputs"]
        fits = [x for x in evs if x["event"] == "fit"]
        starts = [x for x in evs if x["event"] == "start"]
        y = [bool(s["resistant"]) for s in rel["samples"]]
        row = {"id": did, "fits": len(fits), "starts": len(starts),
               "fingerprint": inp[0]["input_fingerprint"] if inp else None,
               "release_fingerprint": rel["input_fingerprint"],
               "samples": inp[0]["samples"] if inp else None, "release_samples": rel["counts"]["samples"],
               "features": inp[0]["features"] if inp else None, "release_features": rel["counts"]["features"],
               "resistant": inp[0]["resistant"] if inp else None, "release_resistant": rel["counts"]["resistant"]}
        pts = unc = cap = mx = irls = pp = rp = 0
        bad_fits = []
        for f in fits:
            tr = f["trace"]
            allp = tr["path"] + tr["refit"]
            pp += len(tr["path"])
            rp += len(tr["refit"])
            pts += len(allp)
            unc += sum(1 for p in allp if not p["converged"])
            cap += sum(1 for p in allp if p["sweep_cap_hits"] > 0)
            mx = max([mx] + [p["irls_steps"] for p in allp])
            irls += sum(p["irls_steps"] for p in allp)
            a = f["audit"]
            if f.get("error") is not None or audit_bad(a):
                bad_fits.append(f"{f['scheme']} {f['outer_fold']}")
            if a is not None and a["path_unconverged"] != sum(1 for p in tr["path"] if not p["converged"]):
                tot["trace_audit_mismatch"] += 1
        avail = [s for s in rel["splits"] if not s.get("unavailable_reason")]
        want_outer = sum(s["folds"] for s in avail)
        outer = [f for f in fits if f["scheme"] != "final"]
        has_final = bool(fits) and fits[-1]["scheme"] == "final"
        conv = (not bad_fits and row["fingerprint"] == row["release_fingerprint"] and len(outer) == want_outer
                and has_final)
        row.update(points=pts, path_points=pp, refit_points=rp, unconverged_points=unc, short_solve_points=cap,
                   max_irls=mx, irls_total=irls, outer_fits=len(outer), release_outer_fits=want_outer,
                   has_final=has_final, bad_fits=bad_fits, converged=conv,
                   lasso_seconds_outer=sum(f["seconds"] for f in outer),
                   final_seconds=fits[-1]["seconds"] if has_final else None)
        # Comparison with the release, split by split.
        splits = []
        for sp in avail:
            L = [l for l in sp["learners"] if l["learner"] == "lasso"]
            assert len(L) == 1
            L = L[0]
            mine = [f for f in fits if f["scheme"] == sp["scheme"]]
            n = len(y)
            score = [None] * n
            source = [None] * n
            for fi in range(sp["folds"]):
                f = [g for g in mine if g["outer_fold"] == fi]
                assert len(f) == 1, (did, sp["scheme"], fi)
                f = f[0]
                if f["audit"] != L["lasso_fits"][fi]:
                    tot["audit_mismatch"] += 1
                w = len(f["model"]["coef"]) if f["model"] is not None else 0
                if w != L["fold_complexity"][fi]:
                    tot["weights_mismatch"] += 1
                d = json.loads(Path(dumpdir, f"{did}__{sp['scheme']}__{fi}.json").read_text())
                assert d["n"] == n and set(d["y"]) == {i for i in range(n) if y[i]}, "labels"
                tr = set(d["train"])
                test = [i for i in range(n) if sp["test_fold"][i] == fi]
                assert test == [i for i in range(n) if i not in tr], "test fold"
                if d["scores"] is None:
                    prev = sum(1 for i in d["train"] if y[i]) / len(d["train"])
                    null_folds.append({"id": did, "scheme": sp["scheme"], "fold": fi, "weights": w,
                                       "complexity": d["complexity"], "train": len(d["train"]),
                                       "train_resistant": sum(1 for i in d["train"] if y[i]),
                                       "prevalence": prev, "test": len(test),
                                       "release_scores_all_equal_prevalence": all(L["scores"][i] == prev for i in test),
                                       "release_called_metrics_present": "called_metrics" in L,
                                       "audit": f["audit"]})
                    for i in test:
                        score[i] = prev
                        source[i] = "prevalence"
                else:
                    for i in test:
                        score[i] = d["scores"][i]
                        source[i] = "binary"
            assert all(s is not None for s in score)
            rs = L["scores"]
            same = sum(1 for a, b in zip(score, rs) if a == b and math.copysign(1, a) == math.copysign(1, b))
            fb = sum(1 for s in source if s == "binary")
            fbs = sum(1 for a, b, s in zip(score, rs, source) if s == "binary" and a == b)
            calls = sum(1 for a, b in zip(score, rs) if (a >= 0.5) != (b >= 0.5))
            ms, mr = measures(y, score), measures(y, rs)
            rm = L["metrics"]
            own = {"tp_fn_fp_tn": (mr["tp"], mr["fn"], mr["fp"], mr["tn"]) == (rm["tp"], rm["fn"], rm["fp"], rm["tn"]),
                   "recall": mr["recall"] - rm["sensitivity"],
                   "balanced_accuracy": mr["balanced_accuracy"] - rm["balanced_accuracy"],
                   "auroc": None if rm["auroc"] is None else mr["auroc"] - rm["auroc"]}
            splits.append({"scheme": sp["scheme"], "n": n, "bit_identical": same, "from_binary": fb,
                           "from_binary_identical": fbs, "calls_changed": calls,
                           "scan": {k: ms[k] for k in ("recall", "log_loss", "balanced_accuracy", "auroc")},
                           "release": {k: mr[k] for k in ("recall", "log_loss", "balanced_accuracy", "auroc")},
                           "release_field": {"recall": rm["sensitivity"], "balanced_accuracy": rm["balanced_accuracy"],
                                             "auroc": rm["auroc"], "log_loss": rm.get("log_loss")},
                           "my_code_vs_release_fields": own})
            tot["scores_compared"] += n
            tot["scores_bit_identical"] += same
            tot["scores_from_binary"] += fb
            tot["scores_from_binary_identical"] += fbs
            tot["null_fold_scores"] += n - fb
            tot["calls_changed"] += calls
        fin = [f for f in fits if f["scheme"] == "final"]
        rmod = rel["models"]["lasso"]
        smod = fin[0]["model"] if fin else None
        fin_ok = bool(fin) and fin[0]["audit"] == rel["models"]["lasso_fit"] and smod == rmod
        if not fin_ok:
            tot["final_mismatch"] += 1
        row["final_model_equal"] = fin_ok
        row["splits"] = splits
        out["datasets"].append(row)
        tot["datasets"] += 1
        tot["converged"] += int(conv)
        tot["fits"] += len(fits)
        tot["points"] += pts
        tot["path_points"] += pp
        tot["refit_points"] += rp
        tot["unconverged_points"] += unc
        tot["short_solve_points"] += cap
        tot["max_irls"] = max(tot["max_irls"], mx)
        tot["irls_total"] += irls
    out["totals"] = tot
    out["null_folds"] = null_folds
    Path(outpath).write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(tot, indent=1))
    for nf in null_folds:
        print("null fold", nf["id"], nf["scheme"], nf["fold"], "weights", nf["weights"], "prev", nf["prevalence"],
              "release all equal", nf["release_scores_all_equal_prevalence"], "test", nf["test"], "audit", nf["audit"])
    for r in out["datasets"]:
        print(r["id"], r["converged"], r["fits"], r["points"], r["irls_total"], r["max_irls"], r["path_points"],
              r["short_solve_points"], round(r["lasso_seconds_outer"], 2), round(r["final_seconds"], 2),
              [(s["scheme"], s["bit_identical"], s["n"], s["calls_changed"], round(s["scan"]["recall"], 4),
                round(s["scan"]["log_loss"], 4)) for s in r["splits"]])


if __name__ == "__main__":
    main(*sys.argv[1:6])
