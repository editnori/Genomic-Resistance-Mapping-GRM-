"""Item 5, exploratory addition (decided at 2026-11-15T14:55Z, after reading
point 11 of item5.py's output, where the pre-registered duality-gap bound for
grm's point was 1.10e-6, above PLAN.md's 1e-6 reading).

The pre-registered dual point scales grm's residual into the dual box, so its
bound loosens roughly as (KKT residual / lambda) times the loss scale, and lambda
falls along the path. This script builds a tighter dual-feasible point: it
solves each path point's lasso to a KKT residual near machine precision with a
sign-constrained Newton active-set method (started from grm's active set,
adding violating columns and dropping columns that reach zero), then builds the
same dual point (item5.Prob.dual) from that solution. Any dual-feasible point
gives a valid lower bound D on the optimum, whatever produced it, so
P_grm - D bounds how far grm's point is above the optimum.

usage: venv/bin/python -I item5-polish.py DUMP.json OUTDIR
"""
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from item5 import Prob, dense, grm_coef, sigmoid  # noqa: E402


def smooth_parts(P, lam, b0, b, A, s):
    """Objective, gradient and Hessian of the sign-fixed smooth problem on the active set A."""
    e = P.eta(b0, b)
    p = sigmoid(e)
    r = P.cw * (p - P.y) / P.n
    w = P.cw * p * (1.0 - p) / P.n
    XA = P.X[:, A]
    g = np.concatenate(([r.sum()], XA.T @ r + lam * P.sd[A] * s))
    Z = np.hstack((np.ones((P.n, 1)), XA))
    H = Z.T @ (Z * w[:, None])
    return g, H


def polish(P, lam, b0, b, max_outer=200):
    b = b.copy()
    elig = ~P.empty
    log = {"newton_steps": 0, "added": 0, "dropped": 0}
    for _ in range(max_outer):
        for _ in range(200):
            A = np.where((b != 0.0) & elig)[0]
            s = np.sign(b[A])
            g, H = smooth_parts(P, lam, b0, b, A, s)
            try:
                d = -np.linalg.solve(H, g)
            except np.linalg.LinAlgError:
                d = -np.linalg.lstsq(H, g, rcond=None)[0]
            dec = float(-(g @ d))
            if not np.isfinite(dec) or dec <= 1e-32:
                break
            dA = d[1:]
            ratios = np.full(len(A), np.inf)
            m = (np.sign(b[A] + dA) != s) & (dA != 0.0)
            ratios[m] = -b[A][m] / dA[m]
            tmax = min(1.0, float(ratios.min())) if len(A) else 1.0
            t = tmax
            if dec > 1e-20:
                # Armijo backtracking; below 1e-20 the decrease is under round-off, so take the Newton step.
                f0 = P.objective(lam, b0, b)
                while t > 1e-12:
                    nb = b.copy()
                    nb[A] = b[A] + t * dA
                    if P.objective(lam, b0 + t * d[0], nb) <= f0 - 1e-4 * t * dec:
                        break
                    t *= 0.5
            b0 = b0 + t * d[0]
            b[A] = b[A] + t * dA
            log["newton_steps"] += 1
            if t == tmax and tmax < 1.0:
                # The first coefficient to reach zero leaves the active set.
                b[A[int(np.argmin(ratios))]] = 0.0
                log["dropped"] += 1
                continue
            if dec < 1e-28:
                break
        # Add the worst violating inactive column, with the sign that lowers the objective.
        r = P.cw * (sigmoid(P.eta(b0, b)) - P.y) / P.n
        gfull = P.X.T @ r
        viol = np.where((b == 0.0) & elig, np.abs(gfull) - lam * P.sd, -np.inf)
        k = int(np.argmax(viol))
        if viol[k] <= 1e-14:
            break
        b[k] = -np.sign(gfull[k]) * 1e-12 / P.sd[k]
        log["added"] += 1
    return b0, b, log


def dual2(P, lam, b0, b):
    """Like item5.Prob.dual, but the intercept condition is met by moving each
    probability in proportion to p(1 - p), which keeps every q in [0, 1] (the
    shift of Prob.dual can push a probability near 0 or 1 out of range)."""
    a = P.cw / P.n
    p = sigmoid(P.eta(b0, b))
    w = p * (1.0 - p)
    sw = math.fsum((a * w).tolist())
    c = math.fsum((a * (p - P.y)).tolist()) / sw if sw > 0 else 0.0
    q0 = p - c * w
    if (q0 < 0).any() or (q0 > 1).any():
        return None
    v = a * (q0 - P.y)
    corr = np.abs(P.X.T @ v)[~P.empty] / P.sd[~P.empty]
    s = min(1.0, lam / corr.max()) if corr.max() > 0 else 1.0
    q = P.y + s * (q0 - P.y)
    with np.errstate(divide="ignore", invalid="ignore"):
        h = np.where((q > 0) & (q < 1), q * np.log(q) + (1 - q) * np.log1p(-q), 0.0)
    return -math.fsum((a * h).tolist())


def main(argv):
    dump, outdir = Path(argv[0]), Path(argv[1])
    d = json.loads(dump.read_text())
    feats = d["features"]
    Xall, yall = dense(d)
    P = Prob(Xall, yall, sorted(d["fit"]))
    rows = []
    for k, f in enumerate(d["path"]):
        lam = f["model"]["lambda"]
        gb0, gb = grm_coef(f["model"], feats)
        t0 = time.monotonic()
        pb0, pb, log = polish(P, lam, gb0, gb)
        Pg, Pp = P.objective(lam, gb0, gb), P.objective(lam, pb0, pb)
        Dg, Dp = P.dual(lam, gb0, gb), P.dual(lam, pb0, pb)
        Dg2, Dp2 = dual2(P, lam, gb0, gb), dual2(P, lam, pb0, pb)
        Ds = [x for x in (Dg, Dp, Dg2, Dp2) if x is not None]
        row = {"point": k, "lambda": lam, "grm_active": int((gb != 0).sum()), "polished_active": int((pb != 0).sum()),
               "active_sets_equal": bool(np.array_equal(gb != 0, pb != 0)),
               "P_grm": Pg, "P_polished": Pp, "P_grm_minus_P_polished": Pg - Pp,
               "kkt_grm": P.kkt(lam, gb0, gb), "kkt_polished": P.kkt(lam, pb0, pb),
               "D_from_grm": Dg, "D_from_polished": Dp, "D2_from_grm": Dg2, "D2_from_polished": Dp2,
               "certified_bound_grm": Pg - max(Ds) if Ds else None,
               "polished_gap": Pp - max(Ds) if Ds else None,
               "max_abs_coef_change_x_sd": float(np.max(np.abs((pb - gb) * P.sd))),
               "seconds": time.monotonic() - t0, **log}
        rows.append(row)
        print(f"point {k}: lam {lam:.3e} active {row['grm_active']}/{row['polished_active']} same {row['active_sets_equal']} "
              f"grm-pol {Pg - Pp:+.2e} bound_grm {row['certified_bound_grm']} kkt_pol {row['kkt_polished']:.1e} "
              f"gap_pol {row['polished_gap']} {row['seconds']:.1f}s", flush=True)
    out = {"dataset": d["dataset"], "scheme": d["scheme"], "outer_fold": d["outer_fold"],
           "label": "exploratory, decided after reading item5.py point 11", "points": rows}
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "item5-polish.json").write_text(json.dumps(out, indent=1) + "\n")


if __name__ == "__main__":
    main(sys.argv[1:])
