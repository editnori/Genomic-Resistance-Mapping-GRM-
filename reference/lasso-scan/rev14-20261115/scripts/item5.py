"""Item 5 of the lassoscan brief (PLAN.md, "Item 5"): on
mycobacterium_tuberculosis__delamanid, random split, outer fold 0, at revision
14: near-twin columns, the gap between grm's objective and an independent
solve at every tuning-path point with KKT residuals and a duality-gap bound,
and the held-out calls of the two solves.

Reads the instrumented binary's full dump of that fit (GRM_LASSO_DUMP_FULL).
The objective, KKT residual and log loss follow src/learn/lasso.rs.

usage: venv/bin/python -I item5.py DUMP.json OUTDIR [--saga-seconds 1800]

Run 2 (2026-11-15T14:58Z) changes two things from run 1, neither a method:
the bound is None, not a crash, when no dual point can be built at a point;
and each independent path solution is saved in OUTDIR/ind/ so a stopped run
resumes without solving again.
"""
import json
import math
import sys
import time
from pathlib import Path

import numpy as np
import cvxpy as cp

THRESH = 1e-8  # |b| x sd below this is set to zero before the independent KKT residual


def sigmoid(e):
    return np.where(e >= 0, 1.0 / (1.0 + np.exp(-np.abs(e))), np.exp(-np.abs(e)) / (1.0 + np.exp(-np.abs(e))))


class Prob:
    """grm's Problem::new on `rows` (sorted sample indices) and the screened columns."""

    def __init__(self, Xall, yall, rows):
        self.rows = np.asarray(rows)
        X = Xall[self.rows].astype(np.float64)
        n = len(rows)
        cnt = X.sum(axis=0)
        self.cleared = cnt == n  # grm clears a column present on every sample
        X[:, self.cleared] = 0.0
        cnt = X.sum(axis=0)
        self.empty = cnt == 0
        p = cnt / n
        self.sd = np.maximum(np.sqrt(p * (1 - p)), 1e-6)
        self.X = X
        self.y = yall[self.rows].astype(np.float64)
        pos = self.y.sum()
        self.cw = np.where(self.y == 1.0, n / (2.0 * max(pos, 1.0)), n / (2.0 * max(n - pos, 1.0)))
        self.n = n

    def eta(self, b0, b):
        return b0 + self.X @ b

    def objective(self, lam, b0, b):
        e = self.eta(b0, b)
        z = np.where(self.y == 1.0, -e, e)
        loss = self.cw * (np.maximum(z, 0.0) + np.log1p(np.exp(-np.abs(z))))
        return math.fsum(loss.tolist()) / self.n + lam * math.fsum((self.sd * np.abs(b)).tolist())

    def kkt(self, lam, b0, b):
        e = self.eta(b0, b)
        res = self.cw * (sigmoid(e) - self.y) / self.n
        worst = abs(math.fsum(res.tolist()))
        g = (self.X.T @ res) / self.sd
        v = np.where(b != 0.0, np.abs(g + lam * np.sign(b)), np.abs(g) - lam)
        v = v[~self.empty]
        return max(worst, float(v.max()) if len(v) else 0.0)

    def dual(self, lam, b0, b):
        """A dual-feasible point built from (b0, b) and its value, a lower bound on the optimum."""
        a = self.cw / self.n
        p = sigmoid(self.eta(b0, b))
        c = math.fsum((a * (p - self.y)).tolist()) / math.fsum(a.tolist())
        v = a * (p - c - self.y)
        corr = np.abs(self.X.T @ v)[~self.empty] / self.sd[~self.empty]
        s = min(1.0, lam / corr.max()) if corr.max() > 0 else 1.0
        q = self.y + s * (p - c - self.y)
        if (q < 0).any() or (q > 1).any():
            return None
        with np.errstate(divide="ignore", invalid="ignore"):
            h = np.where((q > 0) & (q < 1), q * np.log(q) + (1 - q) * np.log1p(-q), 0.0)
        return -math.fsum((a * h).tolist())


def dense(d):
    n, cols = d["n"], d["columns"]
    X = np.zeros((n, len(cols)), dtype=np.uint8)
    for k, c in enumerate(cols):
        X[c, k] = 1
    y = np.zeros(n, dtype=np.uint8)
    y[d["y"]] = 1
    return X, y


def grm_coef(model, feats):
    pos = {f: k for k, f in enumerate(feats)}
    b = np.zeros(len(feats))
    for f, w in model["coef"]:
        b[pos[f]] = w
    return model["intercept"], b


def clarabel_solve(P, lam):
    keep = np.where(~P.empty)[0]
    Xk = P.X[:, keep]
    b0 = cp.Variable()
    b = cp.Variable(len(keep))
    z = b0 + Xk @ b
    # n times grm's objective: class weights of order 1 keep the conic problem well scaled.
    obj = cp.sum(cp.multiply(P.cw, cp.logistic(z) - cp.multiply(P.y, z))) + P.n * lam * cp.norm1(cp.multiply(P.sd[keep], b))
    prob = cp.Problem(cp.Minimize(obj))
    t0 = time.monotonic()
    prob.solve(solver=cp.CLARABEL, tol_gap_abs=1e-10, tol_gap_rel=1e-10, tol_feas=1e-10, max_iter=400)
    secs = time.monotonic() - t0
    full = np.zeros(P.X.shape[1])
    full[keep] = b.value
    st = prob.solver_stats
    return float(b0.value), full, {"status": prob.status, "seconds": secs, "iterations": st.num_iters,
                                   "solver": st.solver_name}


def thresholded(P, b):
    t = b.copy()
    t[np.abs(t) * P.sd < THRESH] = 0.0
    return t


def balanced_log_loss(scores, yall, mask):
    """grm's lasso::log_loss: mean class-balanced log loss on `mask`, clamp 1e-9."""
    m = np.asarray(mask)
    yy = yall[m].astype(bool)
    p = np.clip(scores[m], 1e-9, 1 - 1e-9)
    n, pos = len(m), yy.sum()
    t = np.where(yy, -np.log(p) * n / (2.0 * max(pos, 1)), -np.log(1 - p) * n / (2.0 * max(n - pos, 1)))
    return math.fsum(t.tolist()) / max(n, 1)


def plain_log_loss(scores, yy):
    p = np.clip(scores, 1e-15, 1 - 1e-15)
    return -math.fsum((yy * np.log(p) + (1 - yy) * np.log(1 - p)).tolist()) / len(yy)


def twins(P_fit_X_raw):
    """Columns with another column differing in at most one fitting sample."""
    X = P_fit_X_raw.astype(np.float64)
    c = X.sum(axis=0)
    G = X.T @ X
    H = c[:, None] + c[None, :] - 2.0 * G
    np.fill_diagonal(H, np.inf)
    mn = H.min(axis=1)
    p = X.shape[1]
    # lowest-index identical column (grm's representatives on raw columns)
    rep = np.array([int(np.where(H[k, :k] == 0)[0][0]) if (H[k, :k] == 0).any() else k for k in range(p)])
    return {
        "screened_columns": int(p),
        "fitting_samples": int(X.shape[0]),
        "with_twin_at_most_one": int((mn <= 1).sum()),
        "with_exact_copy": int((mn == 0).sum()),
        "with_twin_exactly_one": int(((mn == 1)).sum()),
        "redundant_copies": int((rep != np.arange(p)).sum()),
        "present_on_every_fitting_sample": int((c == X.shape[0]).sum()),
        "absent_from_every_fitting_sample": int((c == 0).sum()),
        "min_hamming_histogram_0_to_5": [int((mn == h).sum()) for h in range(6)],
    }


def main(argv):
    dump, outdir = Path(argv[0]), Path(argv[1])
    saga_seconds = float(argv[argv.index("--saga-seconds") + 1]) if "--saga-seconds" in argv else 0.0
    outdir.mkdir(parents=True, exist_ok=True)
    d = json.loads(dump.read_text())
    feats = d["features"]
    Xall, yall = dense(d)
    fit, val, train = sorted(d["fit"]), sorted(d["val"]), sorted(d["train"])
    test = [i for i in range(d["n"]) if i not in set(train)]
    out = {"dataset": d["dataset"], "scheme": d["scheme"], "outer_fold": d["outer_fold"], "seed": d["seed"],
           "samples": d["n"], "train": len(train), "fit": len(fit), "val": len(val), "test": len(test),
           "audit": d["audit"]}
    out["twins"] = twins(Xall[fit])
    print("twins", out["twins"], flush=True)
    (outdir / "item5.json").write_text(json.dumps(out, indent=1) + "\n")

    P = Prob(Xall, yall, fit)
    rows = []
    ind_models = []
    (outdir / "ind").mkdir(exist_ok=True)
    for k, f in enumerate(d["path"]):
        lam = f["model"]["lambda"]
        gb0, gb = grm_coef(f["model"], feats)
        saved = outdir / "ind" / f"point{k:02d}.npz"
        if saved.exists():
            z = np.load(saved)
            ib0, ib, info = float(z["b0"]), z["b"], json.loads(str(z["info"]))
            info["resumed"] = True
        else:
            ib0, ib, info = clarabel_solve(P, lam)
            np.savez(saved, b0=ib0, b=ib, info=json.dumps(info))
        ibt = thresholded(P, ib)
        Pg = P.objective(lam, gb0, gb)
        Pi = P.objective(lam, ib0, ib)
        Pit = P.objective(lam, ib0, ibt)
        Dg, Di = P.dual(lam, gb0, gb), P.dual(lam, ib0, ib)
        Dbest = max((x for x in (Dg, Di) if x is not None), default=None)
        row = {"point": k, "lambda": lam, "grm_active": int((gb != 0).sum()), "grm_converged": f["convergence"]["converged"],
               "grm_irls_steps": f["convergence"]["irls_steps"],
               "P_grm": Pg, "P_ind": Pi, "P_ind_thresholded": Pit, "P_grm_minus_P_ind": Pg - Pi,
               "kkt_grm": P.kkt(lam, gb0, gb), "kkt_ind_thresholded": P.kkt(lam, ib0, ibt),
               "ind_active_thresholded": int((ibt != 0).sum()),
               "D_from_grm": Dg, "D_from_ind": Di, "certified_bound_grm": None if Dbest is None else Pg - Dbest,
               "certified_bound_ind": None if Dbest is None else Pi - Dbest,
               "clarabel": info}
        rows.append(row)
        ind_models.append((lam, ib0, ib))
        print(f"point {k}: lam {lam:.3e} grm-ind {Pg - Pi:+.2e} bound_grm {row['certified_bound_grm']} kkt_ind {row['kkt_ind_thresholded']:.1e} "
              f"{info['status']} {info['seconds']:.1f}s", flush=True)
        out["points"] = rows
        (outdir / "item5.json").write_text(json.dumps(out, indent=1) + "\n")

    # Held-out calls.
    vmask = np.array(val)
    sel_grm = d["audit"]["selected_point"]
    vll = []
    for lam, ib0, ib in ind_models:
        sc = sigmoid(ib0 + Xall.astype(np.float64)[:, :] @ np.where(P.empty, 0.0, ib))
        vll.append(balanced_log_loss(sc, yall, vmask))
    sel_ind = int(min(range(len(vll)), key=lambda i: (vll[i], i)))
    grm_vll = []
    for f in d["path"]:
        gb0, gb = grm_coef(f["model"], feats)
        grm_vll.append(balanced_log_loss(sigmoid(gb0 + Xall.astype(np.float64) @ gb), yall, vmask))
    R = Prob(Xall, yall, train)
    tmask = np.array(test)
    gscores = np.array(d["scores"])
    yt = yall[tmask].astype(np.float64)
    calls = {"selected_grm": sel_grm, "selected_ind": sel_ind, "val_log_loss_ind": vll, "val_log_loss_grm": grm_vll,
             "grm_refit_weights": len(d["refit"][-1]["model"]["coef"]) if d["refit"] else 0}
    for name, sel in (("at_ind_selection", sel_ind), ("at_grm_selection", sel_grm)):
        if name == "at_grm_selection" and sel_grm == sel_ind:
            continue
        lam = d["path"][sel]["model"]["lambda"]
        rb0, rb, info = clarabel_solve(R, lam)
        rbt = thresholded(R, rb)
        isc = sigmoid(rb0 + Xall.astype(np.float64) @ rb)
        gb0, gb = grm_coef(d["refit"][-1]["model"], feats)
        calls[name] = {
            "lambda": lam, "clarabel": info,
            "ind_refit_active_thresholded": int((rbt != 0).sum()),
            "refit_P_grm": R.objective(lam, gb0, gb) if sel == sel_grm else None,
            "refit_P_ind": R.objective(lam, rb0, rb),
            "refit_kkt_ind_thresholded": R.kkt(lam, rb0, rbt),
            "test_calls_differ": int(((gscores[tmask] >= 0.5) != (isc[tmask] >= 0.5)).sum()),
            "test_calls_differ_thresholded_ind": int(((gscores[tmask] >= 0.5) != (sigmoid(rb0 + Xall.astype(np.float64) @ rbt)[tmask] >= 0.5)).sum()),
            "max_abs_test_score_change": float(np.abs(gscores[tmask] - isc[tmask]).max()),
            "test_log_loss_grm": plain_log_loss(gscores[tmask], yt),
            "test_log_loss_ind": plain_log_loss(isc[tmask], yt),
            "test_balanced_log_loss_grm": balanced_log_loss(gscores, yall, tmask),
            "test_balanced_log_loss_ind": balanced_log_loss(isc, yall, tmask),
            "test_samples": len(test),
        }
    out["calls"] = calls
    (outdir / "item5.json").write_text(json.dumps(out, indent=1) + "\n")
    print("calls", json.dumps({k: v for k, v in calls.items() if not k.startswith("val_")}, indent=1), flush=True)

    if saga_seconds > 0:
        from sklearn.linear_model import LogisticRegression
        keep = np.where(~P.empty)[0]
        Xs = P.X[:, keep] / P.sd[keep]
        saga = []
        t_all = time.monotonic()
        for k, (lam, ib0, ib) in enumerate(ind_models):
            if time.monotonic() - t_all > saga_seconds:
                saga.append({"point": k, "status": "not run: time box"})
                continue
            t0 = time.monotonic()
            m = LogisticRegression(penalty="l1", solver="saga", C=1.0 / (lam * P.n), tol=1e-10, max_iter=200000)
            m.fit(Xs, P.y, sample_weight=P.cw)
            b = np.zeros(P.X.shape[1])
            b[keep] = m.coef_[0] / P.sd[keep]
            b0 = float(m.intercept_[0])
            saga.append({"point": k, "lambda": lam, "seconds": time.monotonic() - t0, "n_iter": int(m.n_iter_[0]),
                         "P_saga": P.objective(lam, b0, b), "P_saga_minus_P_ind": P.objective(lam, b0, b) - P.objective(lam, ib0, ib),
                         "kkt_saga": P.kkt(lam, b0, b), "status": "ok" if m.n_iter_[0] < 200000 else "max_iter"})
            print("saga", saga[-1], flush=True)
            out["saga"] = saga
            (outdir / "item5.json").write_text(json.dumps(out, indent=1) + "\n")
        out["saga"] = saga
    (outdir / "item5.json").write_text(json.dumps(out, indent=1) + "\n")
    (outdir / "item5.done").write_text(time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()) + "\n")


if __name__ == "__main__":
    main(sys.argv[1:])
