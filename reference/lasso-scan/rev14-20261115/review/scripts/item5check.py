"""revlasso: independent checks of lassoscan's one-fit measurement (brief point 4)
on mycobacterium_tuberculosis__delamanid, random split, outer fold 0.

usage: VENV/bin/python -I -B item5check.py DUMP.json LANE_RESULTS_ITEM5_DIR LANECODE_DIR OUT_JSON SCRATCH_DIR

1. Twins, recounted with Python integers as bit patterns over the fitting
   samples (no X'X): exact copies by pattern, distance one by flipping each bit.
2. Certificate. The lane saved no dual point, only its value. Its point is
   rebuilt here by running the lane's own polish() unchanged (a copy in
   LANECODE_DIR) and the construction of its dual2(), then checked by my own
   code: q in [0, 1], the intercept condition sum a(q - y) = 0, and
   |sum over the samples holding each column of a(q - y)| <= lambda sd for every column,
   with exactly rounded sums (math.fsum). The bound D = sum a H(q) is
   recomputed, and P_grm (my own objective code) minus D bounds grm's excess.
   A second certificate starts from zero, not from grm's solution: a path of
   my own, each point warm-started from my previous point, so the optimum's
   bound does not lean on grm's active set.
3. Held-out calls: grm's selection from the path by grm's lasso::log_loss on
   the tuning-holdout samples, recomputed from the path models; grm's refit
   test scores recomputed from the refit model and checked against the dump.
"""
import json
import math
import sys
import time
from pathlib import Path

import numpy as np


def sig(e):
    return 1.0 / (1.0 + math.exp(-e)) if e >= 0 else math.exp(e) / (1.0 + math.exp(e))


def main(argv):
    dump, lane5, lanecode, outp, scratch = (Path(a) for a in argv[:5])
    d = json.loads(dump.read_text())
    n_all = d["n"]
    feats = d["features"]
    cols = d["columns"]
    fit = sorted(d["fit"])
    val = sorted(d["val"])
    train = sorted(d["train"])
    yset = set(d["y"])
    out = {"fit": len(fit), "val": len(val), "train": len(train), "screened": len(feats)}

    # ---- 1. twins, with Python integers ----
    loc = {s: j for j, s in enumerate(fit)}
    pats = []
    for c in cols:
        v = 0
        for s in c:
            j = loc.get(s)
            if j is not None:
                v |= 1 << j
        pats.append(v)
    nfit = len(fit)
    full = (1 << nfit) - 1
    groups = {}
    for k, v in enumerate(pats):
        groups.setdefault(v, []).append(k)
    exact = [len(groups[v]) > 1 for v in pats]
    one = []
    for k, v in enumerate(pats):
        hit = False
        for j in range(nfit):
            if (v ^ (1 << j)) in groups:
                hit = True
                break
        one.append(hit)
    out["twins"] = {
        "with_exact_copy": sum(exact),
        "min_distance_exactly_one": sum(1 for e, o in zip(exact, one) if o and not e),
        "within_one": sum(1 for e, o in zip(exact, one) if e or o),
        "redundant_copies": sum(len(g) - 1 for g in groups.values()),
        "distinct_patterns": len(groups),
        "present_on_every_fitting_sample": sum(1 for v in pats if v == full),
        "absent_from_every_fitting_sample": sum(1 for v in pats if v == 0),
    }
    print("twins", out["twins"], flush=True)

    # ---- the problem, as Problem::new builds it ----
    n = nfit
    yv = [1.0 if s in yset else 0.0 for s in fit]
    pos = sum(yv)
    a = [(n / (2.0 * max(pos, 1.0)) if t == 1.0 else n / (2.0 * max(n - pos, 1.0))) / n for t in yv]
    supp = []
    for v in pats:
        if v == full:
            supp.append([])  # cleared: present on every fitting sample
        else:
            supp.append([j for j in range(n) if (v >> j) & 1])
    sd = [max(math.sqrt((len(s) / n) * (1 - len(s) / n)), 1e-6) for s in supp]
    pos_of = {f: k for k, f in enumerate(feats)}

    def model_vec(m):
        b = [0.0] * len(feats)
        for f, w in m["coef"]:
            b[pos_of[f]] = w
        return m["intercept"], b

    def etas(b0, b):
        e = [b0] * n
        for k, w in enumerate(b):
            if w != 0.0:
                for j in supp[k]:
                    e[j] += w
        return e

    def objective(lam, b0, b):
        e = etas(b0, b)
        terms = []
        for ej, t, aj in zip(e, yv, a):
            z = -ej if t == 1.0 else ej
            terms.append(aj * (max(z, 0.0) + math.log1p(math.exp(-abs(z)))))
        return math.fsum(terms) + lam * math.fsum(s * abs(w) for s, w in zip(sd, b))

    def kkt(lam, b0, b):
        e = etas(b0, b)
        r = [aj * (sig(ej) - t) for ej, t, aj in zip(e, yv, a)]
        worst = abs(math.fsum(r))
        for k in range(len(b)):
            if not supp[k]:
                continue
            g = math.fsum(r[j] for j in supp[k]) / sd[k]
            v = abs(g + lam * math.copysign(1.0, b[k])) if b[k] != 0.0 else abs(g) - lam
            worst = max(worst, v)
        return worst

    def check_dual(lam, q):
        """My feasibility check of the dual point q and its value."""
        rng_bad = sum(1 for x in q if not (0.0 <= x <= 1.0))
        u = [aj * (qj - t) for aj, qj, t in zip(a, q, yv)]
        eps = math.fsum(u)
        ratio = 0.0
        for k in range(len(supp)):
            if not supp[k]:
                continue
            s = abs(math.fsum(u[j] for j in supp[k]))
            ratio = max(ratio, s / (lam * sd[k]))
        shrink = 1.0
        if ratio > 1.0:
            shrink = 1.0 / ratio
            q = [t + shrink * (qj - t) for qj, t in zip(q, yv)]
        h = []
        for aj, x in zip(a, q):
            if 0.0 < x < 1.0:
                h.append(aj * (-(x * math.log(x)) - (1.0 - x) * math.log1p(-x)))
        return {"q_out_of_range": rng_bad, "intercept_sum": eps, "max_ratio": ratio, "shrink": shrink,
                "D": math.fsum(h)}

    # ---- 2. certificate ----
    sys.path.insert(0, str(lanecode))
    import item5 as L5  # noqa: E402  (the lane's code, copied, unchanged)
    import item5_polish as LP  # noqa: E402
    Xall, yall = L5.dense(d)
    P = L5.Prob(Xall, yall, fit)
    lane_pol = json.loads((lane5 / "item5-polish.json").read_text())["points"]
    lane_main = json.loads((lane5 / "item5.json").read_text())["points"]

    def dual2_point(lam, b0, b):
        # The construction of the lane's dual2(), returning the point q.
        aa = P.cw / P.n
        p = L5.sigmoid(P.eta(b0, b))
        w = p * (1.0 - p)
        sw = math.fsum((aa * w).tolist())
        c = math.fsum((aa * (p - P.y)).tolist()) / sw if sw > 0 else 0.0
        q0 = p - c * w
        if (q0 < 0).any() or (q0 > 1).any():
            return None
        v = aa * (q0 - P.y)
        corr = np.abs(P.X.T @ v)[~P.empty] / P.sd[~P.empty]
        s = min(1.0, lam / corr.max()) if corr.max() > 0 else 1.0
        return (P.y + s * (q0 - P.y)).tolist()

    rows = []
    zb0 = math.log(1.0)  # balanced weights: the null intercept is 0
    zb = np.zeros(len(feats))
    qsave = {}
    for k, f in enumerate(d["path"]):
        t0 = time.monotonic()
        lam = f["model"]["lambda"]
        gb0, gb = model_vec(f["model"])
        Pg = objective(lam, gb0, gb)
        Kg = kkt(lam, gb0, gb)
        # the lane's certified point, rebuilt with its own code
        pb0, pb, _ = LP.polish(P, lam, gb0, np.array(gb))
        q = dual2_point(lam, pb0, pb)
        cl = check_dual(lam, q) if q is not None else None
        # my certificate from zero: my own path of polish() from the null model
        zb0, zb, zlog = LP.polish(P, lam, zb0, zb, max_outer=2000)
        qz = dual2_point(lam, zb0, zb)
        cz = check_dual(lam, qz) if qz is not None else None
        Pz = objective(lam, zb0, list(zb))
        if q is not None:
            qsave[f"q{k:02d}"] = np.array(q)
        lp = lane_pol[k]
        row = {"point": k, "lambda": lam, "P_grm": Pg, "kkt_grm": Kg,
               "lane_P_grm": lp["P_grm"], "lane_D2_from_polished": lp["D2_from_polished"],
               "lane_bound": lp["certified_bound_grm"],
               "rebuilt_check": cl, "rebuilt_bound": None if cl is None else Pg - cl["D"],
               "rebuilt_D_minus_lane_D2": None if cl is None or lp["D2_from_polished"] is None else cl["D"] - lp["D2_from_polished"],
               "zero_start_P": Pz, "zero_start_check": cz, "zero_start_bound": None if cz is None else Pg - cz["D"],
               "zero_start_active": int((zb != 0).sum()), "zero_start_kkt": kkt(lam, zb0, list(zb)),
               "zero_start_P_minus_D": None if cz is None else Pz - cz["D"],
               "grm_active": sum(1 for w in gb if w != 0.0),
               "lane_P_ind": lane_main[k]["P_ind"], "lane_clarabel_status": lane_main[k]["clarabel"]["status"],
               "seconds": time.monotonic() - t0}
        rows.append(row)
        print(f"point {k} lam {lam:.3e} Pg {Pg:.15f} lanePg-Pg {lp['P_grm'] - Pg:+.1e} kkt {Kg:.1e} "
              f"rebuilt: rng {cl and cl['q_out_of_range']} eps {cl and cl['intercept_sum']:.1e} ratio {cl and cl['max_ratio']!r} "
              f"bound {row['rebuilt_bound']:.2e} dD {row['rebuilt_D_minus_lane_D2']} | zero-start: active {row['zero_start_active']} "
              f"bound {row['zero_start_bound']} P-D {row['zero_start_P_minus_D']} ratio {cz and cz['max_ratio']!r} {row['seconds']:.1f}s",
              flush=True)
    out["certificate"] = rows
    np.savez_compressed(scratch / "item5-dual-points.npz", **qsave)

    # ---- 3. held-out calls on grm's side ----
    def predict(m, samples):
        b0, b = model_vec(m)
        e = {s: b0 for s in samples}
        sset = set(samples)
        for k, w in enumerate(b):
            if w != 0.0:
                for s in cols[k]:
                    if s in sset:
                        e[s] += w
        return {s: sig(v) for s, v in e.items()}

    def grm_log_loss(sc, mask):
        p_ = sum(1 for s in mask if s in yset)
        nn = len(mask)
        tot = 0.0
        for s in mask:  # grm sums in sample order, without fsum
            p = min(max(sc[s], 1e-9), 1.0 - 1e-9)
            tot += -math.log(p) * nn / (2.0 * max(p_, 1)) if s in yset else -math.log(1.0 - p) * nn / (2.0 * max(nn - p_, 1))
        return tot / max(nn, 1)

    vll = [grm_log_loss(predict(f["model"], val), val) for f in d["path"]]
    sel = min(range(len(vll)), key=lambda i: (vll[i], i))
    trs = set(train)
    test = [s for s in range(n_all) if s not in trs]
    refit = d["refit"][-1]["model"]
    ts = predict(refit, test)
    dump_scores = d["scores"]
    out["calls"] = {"selected_recomputed": sel, "selected_in_audit": d["audit"]["selected_point"],
                    "refit_points": len(d["refit"]), "refit_weights": len(refit["coef"]),
                    "refit_lambda_equals_selected_lambda": refit["lambda"] == d["path"][sel]["model"]["lambda"],
                    "test_samples": len(test),
                    "test_score_max_abs_diff_vs_dump": max(abs(ts[s] - dump_scores[s]) for s in test),
                    "test_resistant_calls": sum(1 for s in test if ts[s] >= 0.5),
                    "val_log_loss_grm_point3": vll[3]}
    print("calls", out["calls"], flush=True)
    Path(outp).write_text(json.dumps(out, indent=1) + "\n")


if __name__ == "__main__":
    main(sys.argv[1:])
