"""The runner of plan block C7 (docs/outside-learners.md): one outside learner on GRM's own folds.

  python3 -I outside.py fit --arm ARM --export DIR --scheme SCHEME --seed S --out FILE
                            [--fold-file F] [--unitigs MATRIX] [--timing FILE]

ARM is xgboost, lightgbm or elasticnet. DIR is one dataset's `grm export`
directory (matrix.bin, dataset.json) at seed S. SCHEME is lineage or random,
read from the export's plan of that name, or study-whole with --fold-file F,
the reviewed fold file of that dataset, whose outer folds are read and whose
inner folds and boosting holdout are drawn here by GRM's rule (question 2).
The elastic net reads the unitig build MATRIX (`grm exp-input-build
unitigs`, its matrix folder) projected on the export's genomes.

Per outer fold it fits the arm on that fold's training genomes only and
calls the fold's test genomes at 0.5, and it writes FILE, format
grm-c7-outside-1: per fold the genomes by class, the fitting and holdout
genomes by class, the counts of rows offered and kept, the rounds run and
kept or the grid point picked, and the calls (n, tp, fp, tn, fn, loss_sum).
FILE holds no genome id, no row per genome, no row index, no feature, no
model and no weight, gain or importance. --timing writes the wall and CPU
seconds and the peak memory beside it; FILE holds no time, so two runs write
the same bytes.

The libraries are imported only by the arm that uses them: XGBoost 3.0.5,
LightGBM 4.6.0 and scikit-learn 1.7.2, pinned in requirements.txt. Every
other part needs numpy alone.
"""

import argparse
import hashlib
import importlib.util
import json
import math
import os
import resource
import sys
import time
import warnings
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
SCRIPTS = HERE.parent


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# The export reader, the chi-square screen and the cohort projection are D0's (scripts/assumption-d0-c8/d0.py),
# which recomputes `learn::screen` and `Dataset::project` and is held to the export by its own tests.
D0 = _load("d0_for_c7", SCRIPTS / "assumption-d0-c8" / "d0.py")

FORMAT = "grm-c7-outside-1"
REVISION = 14
ARMS = {"xgboost": "gbt", "lightgbm": "gbt", "elasticnet": "lasso"}
DRAWN = ("lineage", "random")
STUDY = "study-whole"
SCREEN = 4000
INNER_K = 3
GBT_HOLDOUT = 0.15
CUTOFF = 0.5
LOSS_CLIP = 1e-4  # src/cutoff.rs
LASSO_CLIP = 1e-9  # learn::lasso::log_loss
FOLD_SEED_STRIDE = 101  # ablate::plan::fold_seed
VERSIONS = {"xgboost": "3.0.5", "lightgbm": "4.6.0", "sklearn": "1.7.2", "numpy": "1.26.4"}
LIBRARY = {"xgboost": "xgboost", "lightgbm": "lightgbm", "elasticnet": "sklearn"}

# GRM's boosted trees (src/learn/gbt.rs, GbtParams::default).
GRM_GBT = {"rounds": 400, "learning_rate": 0.1, "max_depth": 4, "lambda": 1.0, "min_child_weight": 1.0,
           "colsample": 0.5, "patience": 25}
GBT_PARAMS_DEFAULT = dict(GRM_GBT)  # what the cross-check's map is held to
LIGHTGBM = {"objective": "binary", "boost_from_average": False, "learning_rate": GRM_GBT["learning_rate"],
            "max_depth": GRM_GBT["max_depth"], "num_leaves": 2 ** GRM_GBT["max_depth"], "lambda_l2": GRM_GBT["lambda"],
            "min_sum_hessian_in_leaf": GRM_GBT["min_child_weight"], "feature_fraction": GRM_GBT["colsample"],
            "metric": "binary_logloss", "deterministic": True, "force_col_wise": True, "num_threads": 1,
            "verbosity": -1}
# The elastic net (question 7 and section 4 of the design).
MIXES = (0.25, 0.5, 0.75, 1.0)
PENALTIES = 10
PENALTY_SPAN = 1e-3
ENET_TOL = 1e-4
ENET_MAX_ITER = 2000
ZERO = 1e-10  # a saga weight this close to 0 is 0 (a fit at lambda_max leaves ~1e-15)

REPORT_KEYS = {"format", "evaluation_revision", "dataset", "scheme", "seed", "arm", "pairs_with", "library", "version",
               "settings", "input", "plan_sha256", "input_fingerprint", "matrix_sha256", "fold_file_sha256", "samples",
               "resistant", "tested", "tested_resistant", "min_support", "input_rows", "unavailable_reason", "folds",
               "constant_folds", "one_tree_folds", "null_fit_folds", "unconverged_refits"}
FOLD_KEYS = {"fold", "train", "train_resistant", "test", "test_resistant", "fitting", "fitting_resistant", "holdout",
             "holdout_resistant", "candidates", "screened", "constant", "rounds_run", "rounds_kept", "grid_points",
             "grid_unconverged", "picked_mix", "picked_penalty", "refit_converged", "nonzero", "null_fit", "calls"}
CALL_KEYS = {"n", "tp", "fp", "tn", "fn", "loss_sum"}


class Refuse(Exception):
    """An input this runner will not fit on."""


def need(ok, message):
    if not ok:
        raise Refuse(message)


# ---------------------------------------------------------------- bits

def words_of(mask):
    """A boolean vector as little-endian 64-bit words."""
    n = mask.shape[0]
    packed = np.packbits(mask.astype(np.uint8), bitorder="little")
    out = np.zeros(((n + 63) // 64) * 8, dtype=np.uint8)
    out[: packed.shape[0]] = packed
    return out.view("<u8")


def dense(rows, n):
    """Pattern rows (k x words) as a genomes x k 0/1 array."""
    if rows.shape[0] == 0:
        return np.zeros((n, 0), dtype=np.uint8)
    bits = np.unpackbits(np.ascontiguousarray(rows).view(np.uint8), axis=1, bitorder="little")[:, :n]
    return np.ascontiguousarray(bits.T)


def fold_seed(repeat_seed, fold):
    return (repeat_seed + FOLD_SEED_STRIDE * fold) % (1 << 64)


def rust_round(x):
    """f64::round for x >= 0: halves away from zero."""
    return int(math.floor(x + 0.5))


def canonical_sha256(plan):
    body = {k: v for k, v in plan.items() if k != "sha256"}
    return hashlib.sha256(json.dumps(body, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode()).hexdigest()


def draw_rng(*parts):
    """numpy's generator seeded from the sha256 of the parts, for the draws GRM makes in Rust (question 2)."""
    digest = hashlib.sha256("|".join(str(p) for p in parts).encode()).digest()
    return np.random.default_rng(int.from_bytes(digest[:8], "little"))


# ---------------------------------------------------------------- GRM's partition rules, for the study split

def stratified_folds(y, k, rng):
    """learn::stratified_folds: each class shuffled and dealt round the k folds."""
    fold = np.zeros(y.shape[0], dtype=np.int64)
    for cls in (True, False):
        members = np.flatnonzero(y == cls)
        members = members[rng.permutation(members.shape[0])]
        fold[members] = np.arange(members.shape[0]) % k
    return fold


def grouped_folds(y, groups, k, rng):
    """learn::grouped_folds: groups shuffled, placed largest first where the squared class imbalance rises least."""
    members = {}
    for i, g in enumerate(groups):
        members.setdefault(g, []).append(i)
    order = [members[g] for g in sorted(members)]
    order = [order[j] for j in rng.permutation(len(order))]
    order.sort(key=lambda m: -len(m))  # stable, as Rust's sort_by_key
    total_r = float(np.sum(y))
    total_s = float(y.shape[0]) - total_r
    load = [[0.0, 0.0] for _ in range(k)]
    fold = np.zeros(y.shape[0], dtype=np.int64)
    for m in order:
        r = float(sum(1 for i in m if y[i]))
        s = len(m) - r

        def cost(f):
            lr, ls = load[f]
            return ((lr + r) ** 2 - lr ** 2) / max(total_r, 1.0) ** 2 + ((ls + s) ** 2 - ls ** 2) / max(total_s, 1.0) ** 2

        best = min(range(k), key=cost)  # the first of equal costs, as Iterator::min_by
        load[best][0] += r
        load[best][1] += s
        fold[m] = best
    return fold


def training_folds(y, train, groups, k, rng):
    """learn::training_folds: folds of the training genomes only; -1 elsewhere."""
    idx = np.flatnonzero(train)
    local = grouped_folds(y[idx], [groups[i] for i in idx], k, rng) if groups is not None else stratified_folds(y[idx], k, rng)
    out = np.full(y.shape[0], -1, dtype=np.int64)
    out[idx] = local
    return out


def both(mask, y):
    pos = int(np.sum(mask & y))
    return 0 < pos < int(np.sum(mask))


def tuning_holdout(train, y, groups, frac, rng):
    """learn::tuning_holdout: (fit, holdout) masks."""
    if groups is not None:
        n_groups = len({groups[i] for i in np.flatnonzero(train)})
        need(n_groups >= 2, "tuning needs at least two training groups")
        k = min(max(rust_round(1.0 / frac), 2), n_groups)
        folds = training_folds(y, train, groups, k, rng)
        target = float(np.sum(train)) * frac
        best = None
        for f in range(k):
            val = train & (folds == f)
            fit = train & ~val
            if both(fit, y) and both(val, y):
                gap = abs(float(np.sum(val)) - target)
                if best is None or gap < best[0]:
                    best = (gap, fit, val)
        need(best is not None, "no group tuning partition contains both classes")
        fit, val = best[1], best[2]
    else:
        val = np.zeros(y.shape[0], dtype=bool)
        for cls in (True, False):
            members = np.flatnonzero(train & (y == cls))
            members = members[rng.permutation(members.shape[0])]
            take = min(max(rust_round(members.shape[0] * frac), 1), max(members.shape[0] - 1, 0))
            val[members[:take]] = True
        fit = train & ~val
    need(both(fit, y) and both(val, y), "tuning partition lacks one or both classes")
    return fit, val


def tuning_groups(scheme, lineage, duplicates, n):
    """eval::tuning_groups: lineages on the lineage split, else duplicate groups whole, else None."""
    if scheme == "lineage":
        return list(lineage)
    if not duplicates:
        return None
    ids = list(range(n))
    for members in duplicates:
        for i in members:
            ids[i] = members[0]
    return ids


# ---------------------------------------------------------------- plans

class Export:
    """One dataset's grm export: labels, lineages, duplicate groups and the matrix (D0's reader)."""

    def __init__(self, path):
        self.e = D0.Export(str(path))
        m = self.e.meta
        need(m.get("evaluation_revision") == REVISION, f"the export is of revision {m.get('evaluation_revision')}, not {REVISION}")
        self.meta = m
        self.n = len(m["samples"])
        self.y = np.array([bool(v) for v in m["resistant"]], dtype=bool)
        self.x = self.e.x
        need(m["cols"] == self.n and self.y.shape[0] == self.n, "the export's columns are not its samples")
        need(m["screened_features"] == SCREEN, f"the export screened {m['screened_features']} rows, not {SCREEN}")
        need(m["inner_folds"] == INNER_K, f"the export has {m['inner_folds']} inner folds, not {INNER_K}")


def export_plan(ex, scheme, seed):
    """The export's own plan of a drawn split, checked: (plan dict, folds)."""
    need(scheme in DRAWN, f"{scheme} is not a split the export draws")
    split = ex.e.split(scheme)
    need(split is not None, f"the export has no {scheme} plan")
    plan = split["plan"]
    need(plan["scheme"] == scheme, f"the {scheme} entry holds a {plan['scheme']} plan")
    need(plan["seed"] == seed and plan["repeat"] == 0, f"the export's plan is seed {plan['seed']}, not {seed}")
    need(canonical_sha256(plan) == plan["sha256"], f"the {scheme} plan's sha256 is not its content's")
    if plan.get("unavailable_reason"):
        return plan, None
    screens = {s["fold"]: s for s in split["screens"]}
    outer = np.array(plan["outer"], dtype=np.int64)
    folds = []
    for f in range(plan["k"]):
        test = outer == f
        train = ~test
        hold = np.zeros(ex.n, dtype=bool)
        hold[np.array(plan["holdouts"][f]["gbt"], dtype=np.int64)] = True
        inner = np.array([-1 if v is None else v for v in plan["inner"][f]], dtype=np.int64)
        need(f in screens and screens[f]["n_train"] == int(train.sum()), f"fold {f}'s screen is not of its training genomes")
        folds.append({"fold": f, "train": train, "test": test, "hold": hold, "fit": train & ~hold, "inner": inner,
                      "seed": fold_seed(plan["repeat_seed"], f), "screened_rows": np.array(screens[f]["screened_rows"], dtype=np.int64),
                      "candidates": screens[f]["n_candidates"]})
    return plan, folds


def read_fold_file(path, ex):
    """A reviewed fold file's roles on the export's genomes: test fold, -1 train, -2 unused."""
    with open(path, encoding="utf-8") as f:
        doc = json.load(f)
    need(doc.get("format") == "grm-fold-file-v1", "not a grm-fold-file-v1 file")
    need(doc["dataset"] == ex.meta["id"], "the fold file names another dataset")
    need(doc["scheme"] not in DRAWN, f"a fold file may not take the name {doc['scheme']}")
    col = {g: i for i, g in enumerate(ex.meta["samples"])}
    roles = np.full(ex.n, -3, dtype=np.int64)
    listed = [(g, f) for f, ids in enumerate(doc["test"]) for g in ids] + [(g, -1) for g in doc["train"]] + [
        (g, -2) for g in doc["unused"]]
    need(len({g for g, _ in listed}) == len(listed), "the fold file lists a genome twice")
    need(all(g in col for g, _ in listed), "the fold file lists a genome outside the cohort")
    for g, r in listed:
        roles[col[g]] = r
    need(bool(np.all(roles != -3)), "the fold file leaves a cohort genome out")
    with open(path, "rb") as f:
        digest = hashlib.sha256(f.read()).hexdigest()
    return doc, roles, digest


def study_plan(ex, fold_file, seed):
    """Outer folds from the fold file; inner folds, boosting's holdout and the forest's screen drawn here by
    GRM's rule (question 2): the same rule as GRM's fold-file fit, not the same genomes."""
    doc, roles, digest = read_fold_file(fold_file, ex)
    dataset = ex.meta["id"]
    need(doc["scheme"] == STUDY, f"the fold file is {doc['scheme']}, not {STUDY}")
    used = roles != -2
    floor = max(3, int(np.sum(used)) // 200)  # eval::supplied_support_floor
    dups = [[i for i in g if used[i]] for g in ex.meta["duplicates"]]
    dups = [g for g in dups if len(g) > 1]
    groups = tuning_groups(doc["scheme"], None, dups, ex.n)
    folds = []
    for f in range(len(doc["test"])):
        test = roles == f
        train = (roles == -1) | ((roles >= 0) & (roles != f))
        need(both(test, ex.y) and both(train, ex.y), f"study outer fold {f} lacks a class")
        for g in dups:
            need(not (any(test[i] for i in g) and any(train[i] for i in g)), f"study outer fold {f} splits a duplicate group")
        inner = training_folds(ex.y, train, groups, INNER_K, draw_rng(dataset, doc["scheme"], f, "inner"))
        fit, hold = tuning_holdout(train, ex.y, groups, GBT_HOLDOUT, draw_rng(dataset, doc["scheme"], f, "gbt"))
        folds.append({"fold": f, "train": train, "test": test, "hold": hold, "fit": fit, "inner": inner,
                      "seed": fold_seed(seed, f), "screened_rows": None, "candidates": None})
    return {"scheme": doc["scheme"], "sha256": None, "unavailable_reason": None, "floor": floor}, folds, digest


def check_folds(folds, y, n):
    """No outer test genome on any side a fit reads; both classes wherever a fit tunes or scores."""
    for fo in folds:
        f, train, test = fo["fold"], fo["train"], fo["test"]
        need(not bool(np.any(train & test)), f"fold {f} trains on a test genome")
        need(not bool(np.any(fo["hold"] & ~train)), f"fold {f}'s holdout holds a genome outside its training genomes")
        need(not bool(np.any(fo["fit"] & ~train)), f"fold {f}'s fitting side holds a genome outside its training genomes")
        need(not bool(np.any(fo["fit"] & fo["hold"])), f"fold {f}'s fitting side and holdout overlap")
        need(bool(np.all((fo["fit"] | fo["hold"]) == train)), f"fold {f}'s fitting side and holdout are not its training genomes")
        need(bool(np.all((fo["inner"] >= 0) == train)), f"fold {f}'s inner folds are not its training genomes")
        need(bool(np.all(fo["inner"] < INNER_K)), f"fold {f} has an inner fold outside 0..{INNER_K}")
        for side, mask in (("test", test), ("training", train), ("fitting", fo["fit"]), ("holdout", fo["hold"])):
            need(both(mask, y), f"fold {f}'s {side} genomes lack a class")
        for j in range(INNER_K):
            need(both(train & (fo["inner"] != j), y) and both(train & (fo["inner"] == j), y),
                 f"fold {f}'s inner fold {j} lacks a class")


# ---------------------------------------------------------------- screens

def screen_fold(x, train, y, floor, fits, k=SCREEN):
    """learn::training_candidates on `train` and learn::screen's top k on each fitting set of `fits`."""
    zero = np.zeros(x.shape[1], dtype=np.uint64)
    scan = D0.scan_fold(x, words_of(train), zero, zero, {name: words_of(m) for name, m in fits.items()}, words_of(y), floor)
    tops = {name: D0.top(scan["scores"][name], scan["candidates"], k) for name in fits}
    return scan["candidates"].shape[0], tops


# ---------------------------------------------------------------- the arms

def class_weights(y, fit):
    """GRM's sample weights: balanced on the fitting genomes, given to every genome."""
    pos = float(np.sum(y & fit))
    neg = float(np.sum(fit)) - pos
    total = pos + neg
    w = np.where(y, total / (2.0 * max(pos, 1.0)), total / (2.0 * max(neg, 1.0)))
    return w.astype(np.float64)


def import_xgboost():
    """XGBoost 3.0.5's manylinux2014 wheel, the CPU build the cross-check pins, warns at import that the glibc it was
    built for is old. That one warning is expected and let pass; any other stays a warning (an error under -W error)."""
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Your system has an old version of glibc", category=FutureWarning)
        import xgboost
    return xgboost


def xgboost_params(seed):
    import_xgboost()
    sys.path.insert(0, str(SCRIPTS / "learner-crosscheck"))
    try:
        import check_gbt  # the cross-check's map of GRM's settings to XGBoost's names
    finally:
        sys.path.pop(0)
    p = dict(check_gbt.PARAMS)
    names = (("max_depth", "max_depth"), ("learning_rate", "eta"), ("lambda", "lambda"), ("min_child_weight", "min_child_weight"))
    for ours, theirs in names:
        need(p[theirs] == GBT_PARAMS_DEFAULT[ours], f"the cross-check maps {ours} to {p[theirs]}, not GRM's {GBT_PARAMS_DEFAULT[ours]}")
    p.update({theirs: GRM_GBT[ours] for ours, theirs in names})
    p.update(colsample_bytree=GRM_GBT["colsample"], seed=seed % (1 << 31), nthread=1, eval_metric="logloss")
    return p


def sides(fo):
    """One outer fold's fitting genomes, its boosting holdout and its test genomes, as indices."""
    fit = np.flatnonzero(fo["fit"])
    hold = np.flatnonzero(fo["hold"])
    return fit, hold, np.flatnonzero(fo["test"])


def fit_xgboost(X, y, fo, w):
    xgboost = import_xgboost()
    need(xgboost.__version__ == VERSIONS["xgboost"], f"xgboost {xgboost.__version__}, not {VERSIONS['xgboost']}")
    fit, hold, test = sides(fo)
    yf = y.astype(np.float64)
    dfit = xgboost.DMatrix(X[fit], label=yf[fit], weight=w[fit], missing=np.nan, nthread=1)
    dhold = xgboost.DMatrix(X[hold], label=yf[hold], weight=w[hold], missing=np.nan, nthread=1)
    history = {}
    model = xgboost.train(xgboost_params(fo["seed"]), dfit, num_boost_round=GRM_GBT["rounds"], evals=[(dhold, "validation")],
                          early_stopping_rounds=GRM_GBT["patience"], evals_result=history, verbose_eval=False)
    kept = model.best_iteration + 1
    scores = model.predict(xgboost.DMatrix(X[test], missing=np.nan, nthread=1), iteration_range=(0, kept))
    return scores.astype(np.float64), {"rounds_run": len(history["validation"]["logloss"]), "rounds_kept": kept}


def fit_lightgbm(X, y, fo, w):
    import lightgbm
    need(lightgbm.__version__ == VERSIONS["lightgbm"], f"lightgbm {lightgbm.__version__}, not {VERSIONS['lightgbm']}")
    fit, hold, test = sides(fo)
    yf = y.astype(np.float64)
    params = dict(LIGHTGBM, seed=fo["seed"] % (1 << 31))
    dfit = lightgbm.Dataset(X[fit], label=yf[fit], weight=w[fit], params={"verbosity": -1})
    dhold = lightgbm.Dataset(X[hold], label=yf[hold], weight=w[hold], reference=dfit, params={"verbosity": -1})
    history = {}
    model = lightgbm.train(params, dfit, num_boost_round=GRM_GBT["rounds"], valid_sets=[dhold], valid_names=["validation"],
                           callbacks=[lightgbm.early_stopping(GRM_GBT["patience"], first_metric_only=True, verbose=False),
                                      lightgbm.record_evaluation(history)])
    run = len(history["validation"]["binary_logloss"])
    kept = model.best_iteration if model.best_iteration > 0 else run
    return model.predict(X[test], num_iteration=kept).astype(np.float64), {"rounds_run": run, "rounds_kept": kept}


def lasso_log_loss(scores, y):
    """learn::lasso::log_loss: mean class-balanced log loss."""
    n = y.shape[0]
    pos = float(np.sum(y))
    p = np.clip(scores, LASSO_CLIP, 1.0 - LASSO_CLIP)
    terms = np.where(y, -np.log(p) * n / (2.0 * max(pos, 1.0)), -np.log(1.0 - p) * n / (2.0 * max(n - pos, 1.0)))
    return float(np.sum(terms)) / max(n, 1)


def scaled(X, fit):
    """Columns divided by their standard deviation on the fitting genomes (no centring); 0 where constant."""
    p = X[fit].mean(axis=0, dtype=np.float64)
    sd = np.sqrt(p * (1.0 - p))
    scale = np.divide(1.0, sd, out=np.zeros_like(sd), where=sd > 0)
    return scale


def lambda_max(X, y, w, scale, mix):
    g = (w * (y.astype(np.float64) - 0.5)) @ X.astype(np.float64)
    return float(np.max(np.abs(g * scale))) / y.shape[0] / mix if X.shape[1] else 0.0


def enet_path(X, y, w, scale, mix, penalties, seed):
    """scikit-learn's saga down one mix's penalties, warm-started: [(coef, intercept, converged)].

    The columns, divided by their sd, are also centred on the genomes fitted and passed dense. With the intercept
    unpenalized that changes no weight, only the intercept, which is moved back (b - mean . coef), so the model is
    the one stated on the scaled columns; it is done because saga on a sparse matrix damps its intercept steps
    (0.01) and on uncentred 0/1 columns needs several times more epochs. saga visits the genomes in a random order
    drawn from `seed` (the fold seed), so a fit is the same bytes every time. A weight within ZERO of 0 is 0."""
    import sklearn
    from sklearn.exceptions import ConvergenceWarning
    from sklearn.linear_model import LogisticRegression
    need(sklearn.__version__ == VERSIONS["sklearn"], f"scikit-learn {sklearn.__version__}, not {VERSIONS['sklearn']}")
    Xs = X.astype(np.float64) * scale
    centre = Xs.mean(axis=0)
    Xc = np.ascontiguousarray(Xs - centre)
    n = y.shape[0]
    model = LogisticRegression(penalty="elasticnet", solver="saga", l1_ratio=mix, C=1.0, tol=ENET_TOL,
                               max_iter=ENET_MAX_ITER, warm_start=True, fit_intercept=True,
                               random_state=seed % (1 << 31))
    out = []
    for lam in penalties:
        model.set_params(C=1.0 / (lam * n))
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            model.fit(Xc, y.astype(np.int64), sample_weight=w)
        unconverged = any(issubclass(c.category, ConvergenceWarning) for c in caught) or int(model.n_iter_[0]) >= ENET_MAX_ITER
        others = [c for c in caught if not issubclass(c.category, ConvergenceWarning)]
        need(not others, f"scikit-learn warned: {others[0].message if others else ''}")
        coef = model.coef_[0].copy()
        coef[np.abs(coef) <= ZERO] = 0.0
        out.append((coef, float(model.intercept_[0]) - float(centre @ coef), not unconverged))
    return out


def predict_linear(X, scale, coef, intercept):
    z = (X.astype(np.float64) * scale) @ coef + intercept
    return 1.0 / (1.0 + np.exp(-z))


def grid(Xtrain, ytrain, wtrain, scale):
    """The 40 points: per mix, PENALTIES log-spaced from lambda_max(mix) down to PENALTY_SPAN of it."""
    points = []
    for mix in MIXES:
        top = lambda_max(Xtrain, ytrain, wtrain, scale, mix)
        lams = list(np.geomspace(top, top * PENALTY_SPAN, PENALTIES)) if top > 0 else []
        points.append((mix, lams))
    return points


def fit_elasticnet(U, y, fo, floor):
    """One outer fold: grid on the plan's inner folds, each with its own screen, then the refit on every training genome."""
    train, test = fo["train"], fo["test"]
    fits = {"refit": train}
    for j in range(INNER_K):
        fits[f"inner{j}"] = train & (fo["inner"] != j)
    n_cand, tops = screen_fold(U, train, y, floor, fits)
    info = {"candidates": n_cand, "screened": int(tops["refit"].shape[0]), "grid_points": 0, "grid_unconverged": 0}
    if n_cand == 0:
        return None, info
    Xr = dense(np.asarray(U[np.sort(tops["refit"])]), y.shape[0])
    tr = np.flatnonzero(train)
    wr = class_weights(y, train)
    scale_r = scaled(Xr, tr)
    points = grid(Xr[tr], y[tr], wr[tr], scale_r)
    n_points = sum(len(l) for _, l in points)
    info["grid_points"] = n_points
    if n_points == 0:
        return None, info
    pooled = {(m, i): np.full(y.shape[0], np.nan) for m, lams in points for i in range(len(lams))}
    bad = set()
    for j in range(INNER_K):
        part = fits[f"inner{j}"]
        # Every genome outside this inner fit: its held genomes, and the outer test genomes, which no pick reads.
        held = np.flatnonzero(~part)
        Xj = dense(np.asarray(U[np.sort(tops[f"inner{j}"])]), y.shape[0])
        pj = np.flatnonzero(part)
        wj = class_weights(y, part)
        sj = scaled(Xj, pj)
        for mix, lams in points:
            for i, (coef, b, ok) in enumerate(enet_path(Xj[pj], y[pj], wj[pj], sj, mix, lams, fo["seed"])):
                if not ok:
                    bad.add((mix, i))
                pooled[(mix, i)][held] = predict_linear(Xj[held], sj, coef, b)
    info["grid_unconverged"] = len(bad)
    pick_on = np.flatnonzero(train)  # each training genome scored by the inner fit that held it out
    best = None
    for mix, lams in points:
        for i, lam in enumerate(lams):
            if (mix, i) in bad:
                continue
            loss = round(lasso_log_loss(pooled[(mix, i)][pick_on], y[pick_on]), 12)
            # the lowest loss; a tie to the larger penalty, then the larger mix
            key = (loss, -round(float(lam), 12), -mix)
            if best is None or key < best[0]:
                best = (key, mix, i)
    if best is None:
        info.update(refit_converged=0, nonzero=0, null_fit=1)
        return None, info
    _, mix, i = best
    lams = dict(points)[mix][: i + 1]
    coef, b, ok = enet_path(Xr[tr], y[tr], wr[tr], scale_r, mix, lams, fo["seed"])[-1]
    nonzero = int(np.sum(coef != 0))
    info.update(picked_mix=mix, picked_penalty=i, refit_converged=int(ok), nonzero=nonzero, null_fit=int(nonzero == 0))
    if nonzero == 0:
        return None, info
    te = np.flatnonzero(test)
    return predict_linear(Xr[te], scale_r, coef, b), info


# ---------------------------------------------------------------- calls and the report

def calls(scores, y):
    """The 0.5 cutoff's calls and the clipped log loss, summed."""
    p = np.clip(scores, LOSS_CLIP, 1.0 - LOSS_CLIP)
    r = scores >= CUTOFF
    loss = float(np.sum(np.where(y, -np.log(p), -np.log(1.0 - p))))
    return {"n": int(y.shape[0]), "tp": int(np.sum(r & y)), "fp": int(np.sum(r & ~y)), "tn": int(np.sum(~r & ~y)),
            "fn": int(np.sum(~r & y)), "loss_sum": loss}


def settings(arm):
    if arm == "xgboost":
        p = xgboost_params(0)
        p.pop("seed")
        return dict(sorted(p.items()), num_boost_round=GRM_GBT["rounds"], early_stopping_rounds=GRM_GBT["patience"],
                    holdout="the plan's gbt holdout", seed="fold seed mod 2^31")
    if arm == "lightgbm":
        return dict(sorted(LIGHTGBM.items()), num_boost_round=GRM_GBT["rounds"], early_stopping_rounds=GRM_GBT["patience"],
                    holdout="the plan's gbt holdout", seed="fold seed mod 2^31")
    return {"solver": "saga", "penalty": "elasticnet", "class_weight": "balanced on the fitting genomes",
            "columns": "divided by their sd on the fitting genomes; centred for saga, the intercept moved back", "zero": ZERO, "mixes": list(MIXES), "penalties": PENALTIES,
            "penalty_span": PENALTY_SPAN, "tol": ENET_TOL, "max_iter": ENET_MAX_ITER, "inner_folds": INNER_K,
            "screen": SCREEN, "pick": "lowest class-balanced log loss pooled over the inner folds",
            "seed": "fold seed mod 2^31"}


def run(arm, export_dir, scheme, seed, fold_file=None, unitigs=None, fitter=None):
    """Fit `arm` on every outer fold. Returns (report, scores): scores holds each tested genome's score
    in memory for the tests and is never written."""
    need(arm in ARMS, f"no such arm: {arm}")
    ex = Export(export_dir)
    dataset = ex.meta["id"]
    need(Path(export_dir).name == dataset, "the export's folder is not named for its dataset")
    if scheme == STUDY:
        need(fold_file is not None, "the study split needs its fold file")
        plan, folds, ff_sha = study_plan(ex, fold_file, seed)
        floor = plan["floor"]
    else:
        need(fold_file is None, f"a fold file is not read on the {scheme} split")
        plan, folds = export_plan(ex, scheme, seed)
        ff_sha = None
        floor = ex.meta["min_support"]
    report = {"format": FORMAT, "evaluation_revision": REVISION, "dataset": dataset, "scheme": scheme, "seed": seed,
              "arm": arm, "pairs_with": ARMS[arm], "library": LIBRARY[arm], "version": VERSIONS[LIBRARY[arm]],
              "settings": settings(arm), "input": "unitigs" if arm == "elasticnet" else "export",
              "plan_sha256": plan.get("sha256"), "input_fingerprint": ex.meta["input_fingerprint"],
              "matrix_sha256": ex.meta["matrix_sha256"], "fold_file_sha256": ff_sha, "samples": ex.n,
              "resistant": int(ex.y.sum()), "min_support": floor, "unavailable_reason": plan.get("unavailable_reason")}
    scores = np.full(ex.n, np.nan)
    if folds is None:
        report.update(tested=0, tested_resistant=0, input_rows=None, folds=[], constant_folds=0, one_tree_folds=0,
                      null_fit_folds=0, unconverged_refits=0)
        return report, scores
    check_folds(folds, ex.y, ex.n)
    y = ex.y
    if arm == "elasticnet":
        need(unitigs is not None, "the elastic net reads the unitig build")
        U, _ = D0.project(str(unitigs), ex.meta["samples"], floor)
        report["input_rows"] = int(U.shape[0])
    else:
        need(unitigs is None, "a library arm reads the export's rows, not the unitig build")
        report["input_rows"] = int(ex.meta["rows"])
    out_folds = []
    tested = np.zeros(ex.n, dtype=bool)
    for fo in folds:
        train, test = fo["train"], fo["test"]
        row = {"fold": fo["fold"], "train": int(train.sum()), "train_resistant": int((train & y).sum()),
               "test": int(test.sum()), "test_resistant": int((test & y).sum()), "fitting": int(fo["fit"].sum()),
               "fitting_resistant": int((fo["fit"] & y).sum()), "holdout": int(fo["hold"].sum()),
               "holdout_resistant": int((fo["hold"] & y).sum())}
        if arm == "elasticnet":
            s, info = (fitter or fit_elasticnet)(U, y, fo, floor)
        else:
            rows = fo["screened_rows"]
            if rows is None:
                n_cand, tops = screen_fold(ex.x, train, y, floor, {"rf": train})
                rows, fo["candidates"] = tops["rf"], n_cand
            info = {"candidates": int(fo["candidates"]), "screened": int(rows.shape[0])}
            s = None
            if rows.shape[0]:
                # The rows in matrix order: a booster's column order is no rank of the screen.
                X = dense(np.asarray(ex.x[np.sort(rows)]), ex.n).astype(np.float32)
                w = class_weights(y, fo["fit"])
                s, more = (fitter or {"xgboost": fit_xgboost, "lightgbm": fit_lightgbm}[arm])(X, y, fo, w)
                info.update(more)
        te = np.flatnonzero(test)
        if s is None:
            # No feature, or a null fit: the fold's training share of resistant genomes (GRM's null fit,
            # ruling 16's constant model).
            s = np.full(te.shape[0], float((train & y).sum()) / float(train.sum()))
            info["constant"] = 1
        else:
            info["constant"] = 0
        scores[te] = s
        tested |= test
        row.update(info)
        row["calls"] = calls(np.asarray(s, dtype=np.float64), y[te])
        out_folds.append(row)
    report.update(tested=int(tested.sum()), tested_resistant=int((tested & y).sum()), folds=out_folds,
                  constant_folds=sum(r["constant"] for r in out_folds),
                  one_tree_folds=sum(1 for r in out_folds if r.get("rounds_kept") == 1),
                  null_fit_folds=sum(r.get("null_fit", 0) for r in out_folds),
                  unconverged_refits=sum(1 for r in out_folds if r.get("refit_converged") == 0))
    check_report(report)
    return report, scores


def check_report(report):
    """The report's keys are the listed ones: no genome, row, feature or model field can enter it."""
    extra = set(report) - REPORT_KEYS
    need(not extra, f"the report holds unlisted keys: {sorted(extra)}")
    for row in report["folds"]:
        extra = set(row) - FOLD_KEYS
        need(not extra, f"a fold row holds unlisted keys: {sorted(extra)}")
        need(set(row["calls"]) == CALL_KEYS, "a fold's calls hold other keys than n, tp, fp, tn, fn, loss_sum")
        for k, v in row.items():
            need(k == "calls" or isinstance(v, (int, float)), f"a fold row's {k} is not a count")


def dump(report):
    return json.dumps(report, sort_keys=True, indent=1) + "\n"


def write_atomic(path, text):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(text)
    os.replace(tmp, path)


def main(argv):
    ap = argparse.ArgumentParser(prog="outside.py")
    sub = ap.add_subparsers(dest="cmd", required=True)
    f = sub.add_parser("fit")
    f.add_argument("--arm", required=True, choices=sorted(ARMS))
    f.add_argument("--export", required=True)
    f.add_argument("--scheme", required=True, choices=DRAWN + (STUDY,))
    f.add_argument("--seed", required=True, type=int)
    f.add_argument("--out", required=True)
    f.add_argument("--fold-file")
    f.add_argument("--unitigs")
    f.add_argument("--timing")
    a = ap.parse_args(argv)
    started, cpu0 = time.time(), os.times()
    try:
        report, _ = run(a.arm, a.export, a.scheme, a.seed, a.fold_file, a.unitigs)
    except Refuse as e:
        print(f"outside.py: refused: {e}", file=sys.stderr)
        return 1
    write_atomic(a.out, dump(report))
    if a.timing:
        cpu1 = os.times()
        write_atomic(a.timing, json.dumps({"arm": a.arm, "dataset": report["dataset"], "scheme": a.scheme, "seed": a.seed,
                                           "wall_s": round(time.time() - started, 3),
                                           "cpu_s": round((cpu1.user - cpu0.user) + (cpu1.system - cpu0.system), 3),
                                           "peak_kb": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}, sort_keys=True) + "\n")
    print(f"{a.arm} {report['scheme']} s{a.seed}: {len(report['folds'])} folds, tested {report['tested']}, "
          f"constant folds {report['constant_folds']}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
