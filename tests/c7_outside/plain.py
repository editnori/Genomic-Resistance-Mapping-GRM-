"""The second, plain paths the outside-learner runner is checked against, on SYNTHETIC tables only.

- `boost`: second-order boosting of depth-limited trees on 0/1 features, grown level by level with an exact search,
  with XGBoost's or LightGBM's split floor, leaf value and early stopping written out in numpy. LightGBM's floor
  counts a side's genomes from its hessian (the side's hessian times the node's genomes over the node's hessian,
  rounded half up) and drops, before the first tree, a feature that cannot leave `min_count` genomes on both sides
  of the fitting genomes (its feature pre-filter).
- `enet_objective` and `enet_fista`: the elastic net's objective as the runner states it, and an accelerated
  proximal gradient solver for it.
- `screen`: the candidates and the chi-square top k by the export's rule (docs/dataset-export.md), with Python loops.

None of them is fast; each is short enough to read against its rule.
"""

import math

import numpy as np


def _sigmoid(z):
    return 1.0 / (1.0 + np.exp(-z))


def _tree(X, g, h, rows, depth, lam, min_hess, min_count, split_eps, hessian_count, usable):
    """One tree as nested dicts, grown on `rows` with the weighted gradients g and hessians h."""
    G, H = float(np.sum(g[rows])), float(np.sum(h[rows]))
    leaf = {"leaf": -G / (H + lam)}
    if depth == 0:
        return leaf
    best = None
    root = G * G / (H + lam)
    for j in usable:
        on = rows[X[rows, j] == 1]
        off = rows[X[rows, j] == 0]
        GR, HR = float(np.sum(g[on])), float(np.sum(h[on]))
        GL, HL = G - GR, H - HR
        if hessian_count:
            n_on = int(math.floor(HR * rows.shape[0] / H + 0.5)) if H > 0 else 0
            n_off = rows.shape[0] - n_on
        else:
            n_on, n_off = on.shape[0], off.shape[0]
        if n_on < min_count or n_off < min_count:
            continue
        if HL < min_hess or HR < min_hess:
            continue
        gain = GL * GL / (HL + lam) + GR * GR / (HR + lam) - root
        if gain > split_eps and (best is None or gain > best[0]):
            best = (gain, j, off, on)
    if best is None:
        return leaf
    _, j, off, on = best
    return {"feature": j,
            "absent": _tree(X, g, h, off, depth - 1, lam, min_hess, min_count, split_eps, hessian_count, usable),
            "present": _tree(X, g, h, on, depth - 1, lam, min_hess, min_count, split_eps, hessian_count, usable)}


def _apply(tree, X):
    out = np.zeros(X.shape[0])
    stack = [(tree, np.arange(X.shape[0]))]
    while stack:
        node, idx = stack.pop()
        if "leaf" in node:
            out[idx] = node["leaf"]
            continue
        bit = X[idx, node["feature"]] == 1
        stack.append((node["present"], idx[bit]))
        stack.append((node["absent"], idx[~bit]))
    return out


def weighted_log_loss(p, y, w, eps):
    p = np.clip(p, eps, 1.0 - eps)
    return float(np.sum(w * np.where(y == 1, -np.log(p), -np.log(1.0 - p))) / np.sum(w))


def boost(X, y, w, fit, hold, test, rounds, eta, depth, lam, min_hess, patience, min_count=0, split_eps=1e-6,
          eps=1e-16, hessian_count=False):
    """Boosting from margin 0 on the fitting rows, the weighted log loss read on the holdout rows after each
    round, stopped after `patience` rounds without a strictly lower loss, the rounds up to the lowest kept.
    hessian_count: LightGBM's count of a side's genomes and its feature pre-filter (see above).
    Returns (scores of the test rows, rounds run, rounds kept)."""
    ones = X[fit].sum(axis=0)
    if hessian_count:
        usable = [j for j in range(X.shape[1]) if ones[j] >= min_count and fit.shape[0] - ones[j] >= min_count]
    else:
        usable = list(range(X.shape[1]))
    y = y.astype(np.float64)
    margin = np.zeros(X.shape[0])
    trees, losses = [], []
    best, best_round = math.inf, -1
    for r in range(rounds):
        p = _sigmoid(margin)
        g, h = w * (p - y), w * p * (1.0 - p)
        tree = _tree(X, g, h, fit, depth, lam, min_hess, min_count, split_eps, hessian_count, usable)
        step = _apply(tree, X) * eta
        margin = margin + step
        trees.append(step)
        loss = weighted_log_loss(_sigmoid(margin[hold]), y[hold], w[hold], eps)
        losses.append(loss)
        if loss < best:
            best, best_round = loss, r
        elif r - best_round >= patience:
            break
    kept = best_round + 1
    final = np.sum(trees[:kept], axis=0)
    return _sigmoid(final[test]), len(losses), kept


def enet_objective(X, y, w, coef, intercept, lam, mix):
    """Mean weighted log loss plus lam * ((1 - mix) / 2 * |coef|^2 + mix * |coef|_1), the intercept unpenalized."""
    z = X @ coef + intercept
    loss = np.mean(w * (np.logaddexp(0.0, z) - y * z))
    return float(loss + lam * ((1.0 - mix) / 2.0 * np.sum(coef * coef) + mix * np.sum(np.abs(coef))))


def enet_fista(X, y, w, lam, mix, iters=20000, tol=1e-12):
    """The same objective by accelerated proximal gradient with backtracking. Returns (coef, intercept)."""
    n, k = X.shape
    y = y.astype(np.float64)

    def smooth(c, b):
        z = X @ c + b
        p = _sigmoid(z)
        r = w * (p - y) / n
        val = float(np.mean(w * (np.logaddexp(0.0, z) - y * z)) + lam * (1.0 - mix) / 2.0 * np.sum(c * c))
        return val, X.T @ r + lam * (1.0 - mix) * c, float(np.sum(r))

    c, b = np.zeros(k), 0.0
    yc, yb, t, step = c.copy(), b, 1.0, 1.0
    last = math.inf
    for _ in range(iters):
        f, gc, gb = smooth(yc, yb)
        while True:
            zc = yc - step * gc
            nc = np.sign(zc) * np.maximum(np.abs(zc) - step * lam * mix, 0.0)
            nb = yb - step * gb
            fn, _, _ = smooth(nc, nb)
            dc, db = nc - yc, nb - yb
            if fn <= f + float(gc @ dc) + gb * db + (float(dc @ dc) + db * db) / (2.0 * step) + 1e-15:
                break
            step *= 0.5
        t2 = (1.0 + math.sqrt(1.0 + 4.0 * t * t)) / 2.0
        yc = nc + (t - 1.0) / t2 * (nc - c)
        yb = nb + (t - 1.0) / t2 * (nb - b)
        c, b, t = nc, nb, t2
        obj = enet_objective(X, y, w, c, b, lam, mix)
        if abs(last - obj) < tol:
            break
        last = obj
    return c, b


def screen(x, y, train, floor, fit, k):
    """Candidates on `train` (support in [floor, n_train - floor], the lowest row of each training pattern) and the
    chi-square top k on the genomes of `fit`, highest first, equal scores by lower row."""
    n_train = int(np.sum(train))
    seen, cands = set(), []
    for r in range(x.shape[0]):
        s = int(np.sum(x[r] & train))
        if s < floor or s > n_train - floor:
            continue
        key = tuple(np.flatnonzero(x[r] & train))
        if key in seen:
            continue
        seen.add(key)
        cands.append(r)
    pos = float(np.sum(fit & y))
    neg = float(np.sum(fit & ~y))
    scored = []
    for r in cands:
        a = float(np.sum(x[r] & fit & y))
        b = float(np.sum(x[r] & fit & ~y))
        c, d = pos - a, neg - b
        denom = (a + b) * (c + d) * pos * neg
        chi = 0.0 if denom == 0 else (pos + neg) * (a * d - b * c) ** 2 / denom
        scored.append((-chi, r))
    scored.sort()
    return len(cands), [r for _, r in scored[:k]]
