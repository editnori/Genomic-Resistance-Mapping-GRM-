"""Fixed-learner representation comparator for Track B (`grm-rep-compare-v1`).

Compares representation arms (31-mer baseline from `grm export`, protein-aligned
known-target variants, gene/active-site burden, pretrained foundation-model
scores/embeddings, and concatenated hybrid arms) under the EXACT SAME:

- outer folds (`lineage` or `random` from `grm-dataset-export-v1`, or
  `study-whole` from `grm-fold-file-v1`),
- inner folds (fold-safe group-preserving inner partition inside each outer
  training fold),
- cohort row support,
- downstream learner (deterministic pure-NumPy L2-regularized class-balanced
  logistic regression solved via exact primal/Woodbury Newton-IRLS),
- regularization grid and tie-breaking rule (larger lambda on exact ties).

All encoders, screens, scalers, and PCA projections are fit strictly inside
each training partition (inner training folds during hyperparameter selection,
and outer training fold for the final refit). Outer test rows never inform
feature extraction, scaling, projection, lambda selection, or model weights.

Output JSON (`grm-rep-compare-v1`) holds zero genome IDs, zero per-genome rows,
zero sequences/variants/positions, zero weights, and no wall-clock timestamps,
so repeated runs emit byte-identical JSON passing `scripts/c7-outside/scan.py`.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import os
import re
import sys
from pathlib import Path
from typing import Any, Callable, Sequence

import numpy as np


def _load_sibling(name: str, filename: str):
    key = f"grm_point_mutations.{name}"
    if key in sys.modules:
        return sys.modules[key]
    path = Path(__file__).resolve().parent / filename
    spec = importlib.util.spec_from_file_location(key, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[key] = module
    spec.loader.exec_module(module)
    return module


EA = _load_sibling("exposure_audit", "exposure_audit.py")
FT = _load_sibling("features", "features.py")

FORMAT = "grm-rep-compare-v1"
REVISION = 14
DRAWN_SCHEMES = ("lineage", "random")
STUDY_SCHEME = "study-whole"
INNER_K = 3
CUTOFF = 0.5
LOSS_CLIP = 1e-4  # matches src/cutoff.rs and scripts/c7-outside/outside.py
INNER_LOSS_CLIP = 1e-9  # matches learn::lasso::log_loss
DEFAULT_LAMBDA_GRID = (1e-3, 1e-2, 1e-1, 1.0, 10.0, 100.0, 1000.0)
DEFAULT_TOL = 1e-9
DEFAULT_MAX_ITER = 60
ZERO_COEF_TOL = 1e-12
KMER_SCREEN_DEFAULT = 4000

REPORT_KEYS = frozenset({
    "format",
    "evaluation_revision",
    "dataset",
    "scheme",
    "seed",
    "input_fingerprint",
    "plan_sha256",
    "matrix_sha256",
    "fold_file_sha256",
    "samples",
    "resistant",
    "tested",
    "tested_resistant",
    "min_support",
    "baseline_arm",
    "learner_control",
    "scope_note",
    "arms",
    "pairwise_vs_baseline",
})

ARM_KEYS = frozenset({
    "arm",
    "kind",
    "claim_ceiling",
    "exposure_audits",
    "pooled_metrics",
    "constant_folds",
    "null_fit_folds",
    "unconverged_folds",
    "folds",
})

FOLD_KEYS = frozenset({
    "fold",
    "train",
    "train_resistant",
    "test",
    "test_resistant",
    "n_features",
    "picked_lambda",
    "converged",
    "null_fit",
    "constant",
    "calls",
})

CALL_KEYS = frozenset({"n", "tp", "fp", "tn", "fn", "loss_sum"})

POOLED_KEYS = frozenset({
    "n",
    "tp",
    "fp",
    "tn",
    "fn",
    "balanced_accuracy",
    "recall",
    "specificity",
    "vme",
    "me",
    "log_loss",
})

PAIRWISE_KEYS = frozenset({
    "arm",
    "baseline_arm",
    "common_tested_n",
    "delta_balanced_accuracy",
    "delta_recall",
    "delta_specificity",
    "delta_vme",
    "delta_me",
    "delta_log_loss",
})

GENOME_ID_TOKEN_RE = re.compile(r"(?<![\w.])\d+\.\d+(?!\d)")


class Refuse(Exception):
    """Raised when an input fails scientific, alignment, audit, or hygiene gates."""


def _need(cond: bool, message: str) -> None:
    if not cond:
        raise Refuse(message)


# ---------------------------------------------------------------- export & fold loading

def canonical_plan_sha256(plan: dict[str, Any]) -> str:
    body = {k: v for k, v in plan.items() if k != "sha256"}
    raw = json.dumps(body, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()


def _draw_rng(*parts: Any) -> np.random.Generator:
    digest = hashlib.sha256("|".join(str(p) for p in parts).encode("utf-8")).digest()
    return np.random.default_rng(int.from_bytes(digest[:8], "little"))


def _both_classes(mask: np.ndarray, y: np.ndarray) -> bool:
    n = int(np.sum(mask))
    pos = int(np.sum(mask & y))
    return 0 < pos < n


def _stratified_folds(y: np.ndarray, k: int, rng: np.random.Generator) -> np.ndarray:
    fold = np.zeros(y.shape[0], dtype=np.int64)
    for cls in (True, False):
        members = np.flatnonzero(y == cls)
        members = members[rng.permutation(members.shape[0])]
        fold[members] = np.arange(members.shape[0]) % k
    return fold


def _grouped_folds(y: np.ndarray, groups: Sequence[Any], k: int, rng: np.random.Generator) -> np.ndarray:
    members: dict[Any, list[int]] = {}
    for i, g in enumerate(groups):
        members.setdefault(g, []).append(i)
    order = [members[g] for g in sorted(members)]
    order = [order[j] for j in rng.permutation(len(order))]
    order.sort(key=lambda m: -len(m))
    total_r = float(np.sum(y))
    total_s = float(y.shape[0]) - total_r
    load = [[0.0, 0.0] for _ in range(k)]
    fold = np.zeros(y.shape[0], dtype=np.int64)
    for m in order:
        r = float(sum(1 for i in m if y[i]))
        s = len(m) - r

        def cost(f: int) -> float:
            lr, ls = load[f]
            return (
                ((lr + r) ** 2 - lr ** 2) / max(total_r, 1.0) ** 2
                + ((ls + s) ** 2 - ls ** 2) / max(total_s, 1.0) ** 2
            )

        best = min(range(k), key=cost)
        load[best][0] += r
        load[best][1] += s
        fold[m] = best
    return fold


def _training_folds(
    y: np.ndarray,
    train: np.ndarray,
    groups: Sequence[Any] | None,
    k: int,
    rng: np.random.Generator,
) -> np.ndarray:
    idx = np.flatnonzero(train)
    local = (
        _grouped_folds(y[idx], [groups[i] for i in idx], k, rng)
        if groups is not None
        else _stratified_folds(y[idx], k, rng)
    )
    out = np.full(y.shape[0], -1, dtype=np.int64)
    out[idx] = local
    return out


class ExportDataset:
    """Validated reader for a `grm export` directory (`dataset.json` + `matrix.bin`)."""

    def __init__(self, export_dir: str | Path) -> None:
        self.path = Path(export_dir)
        meta_path = self.path / "dataset.json"
        _need(meta_path.is_file(), f"missing dataset.json in {self.path}")
        self.meta = json.loads(meta_path.read_text(encoding="utf-8"))
        m = self.meta
        _need(m.get("format") == "grm-dataset-export-v1", "not a grm-dataset-export-v1 file")
        _need(m.get("evaluation_revision") == REVISION,
              f"evaluation_revision is {m.get('evaluation_revision')}, expected {REVISION}")
        self.samples: tuple[str, ...] = tuple(m["samples"])
        self.n: int = len(self.samples)
        _need(len(set(self.samples)) == self.n and self.n > 0, "invalid samples list in export")
        self.y: np.ndarray = np.array([bool(v) for v in m["resistant"]], dtype=bool)
        _need(m["cols"] == self.n and self.y.shape[0] == self.n, "export cols/labels mismatch")
        words = (self.n + 63) // 64
        _need(m["words"] == words, "export words count does not match samples")
        mpath = self.path / "matrix.bin"
        _need(mpath.is_file(), f"missing matrix.bin in {self.path}")
        expected_bytes = m["rows"] * m["words"] * 8
        _need(mpath.stat().st_size == expected_bytes, "matrix.bin size mismatch")
        actual_sha = hashlib.sha256(mpath.read_bytes()).hexdigest()
        _need(actual_sha == m["matrix_sha256"], "matrix.bin sha256 mismatch")
        if m["rows"] > 0:
            raw = np.frombuffer(mpath.read_bytes(), dtype="<u8").reshape((m["rows"], m["words"]))
            bits = np.unpackbits(raw.view(np.uint8), axis=1, bitorder="little")[:, : self.n]
            self.kmer_dense: np.ndarray = np.ascontiguousarray(bits.T, dtype=np.uint8)
        else:
            self.kmer_dense = np.zeros((self.n, 0), dtype=np.uint8)

    def split_entry(self, scheme: str) -> dict[str, Any] | None:
        for s in self.meta.get("splits", ()):
            if s.get("plan", {}).get("scheme") == scheme:
                return s
        return None


def load_folds(
    export: ExportDataset,
    scheme: str,
    seed: int,
    fold_file: str | Path | None = None,
) -> tuple[dict[str, Any], list[dict[str, Any]], int, str | None]:
    """Load and verify outer + inner fold definitions without leaking test rows."""
    if scheme in DRAWN_SCHEMES:
        _need(fold_file is None, f"fold_file cannot be provided for drawn scheme {scheme!r}")
        split = export.split_entry(scheme)
        _need(split is not None, f"export lacks split {scheme!r}")
        plan = split["plan"]
        _need(plan["scheme"] == scheme, f"plan scheme mismatch: {plan['scheme']!r}")
        _need(plan["seed"] == seed and plan.get("repeat", 0) == 0,
              f"plan seed {plan['seed']} != {seed}")
        _need(canonical_plan_sha256(plan) == plan["sha256"], "split plan sha256 mismatch")
        _need(not plan.get("unavailable_reason"),
              f"split {scheme!r} is unavailable: {plan.get('unavailable_reason')}")
        outer = np.array(plan["outer"], dtype=np.int64)
        _need(outer.shape[0] == export.n, "plan outer length mismatch")
        folds: list[dict[str, Any]] = []
        for f in range(plan["k"]):
            test = outer == f
            train = ~test
            inner = np.array([-1 if v is None else int(v) for v in plan["inner"][f]], dtype=np.int64)
            folds.append({"fold": f, "train": train, "test": test, "inner": inner})
        floor = int(export.meta["min_support"])
        ff_sha = None
        plan_info = {"sha256": plan["sha256"], "inner_k": int(export.meta["inner_folds"])}
    elif scheme == STUDY_SCHEME:
        _need(fold_file is not None, "study-whole scheme requires a grm-fold-file-v1 file")
        fpath = Path(fold_file)
        fbytes = fpath.read_bytes()
        ff_sha = hashlib.sha256(fbytes).hexdigest()
        doc = json.loads(fbytes.decode("utf-8"))
        _need(doc.get("format") == "grm-fold-file-v1", "not a grm-fold-file-v1 file")
        _need(doc.get("dataset") == export.meta["id"], "fold file dataset differs from export")
        _need(doc.get("scheme") == STUDY_SCHEME, f"fold file scheme is {doc.get('scheme')!r}, not {STUDY_SCHEME!r}")

        col = {g: i for i, g in enumerate(export.samples)}
        roles = np.full(export.n, -3, dtype=np.int64)
        listed = (
            [(g, f) for f, ids in enumerate(doc["test"]) for g in ids]
            + [(g, -1) for g in doc["train"]]
            + [(g, -2) for g in doc["unused"]]
        )
        _need(len({g for g, _ in listed}) == len(listed), "fold file lists a sample twice")
        _need(all(g in col for g, _ in listed), "fold file lists a sample outside export cohort")
        for g, r in listed:
            roles[col[g]] = r
        _need(bool(np.all(roles != -3)), "fold file omits a cohort sample")

        used = roles != -2
        floor = max(3, int(np.sum(used)) // 200)
        dups = [[i for i in g if used[i]] for g in export.meta.get("duplicates", ())]
        dups = [g for g in dups if len(g) > 1]

        # Group-preserve duplicate groups (and lineages when available) inside study inner folds.
        groups: list[int] | None = None
        if dups:
            groups = list(range(export.n))
            for members in dups:
                for i in members:
                    groups[i] = members[0]

        folds = []
        for f in range(len(doc["test"])):
            test = roles == f
            train = (roles == -1) | ((roles >= 0) & (roles != f))
            for g in dups:
                _need(not (any(test[i] for i in g) and any(train[i] for i in g)),
                      f"study outer fold {f} splits a k-mer-identical duplicate group")
            inner = _training_folds(
                export.y,
                train,
                groups,
                INNER_K,
                _draw_rng(export.meta["id"], STUDY_SCHEME, seed, f, "inner"),
            )
            folds.append({"fold": f, "train": train, "test": test, "inner": inner})
        plan_info = {"sha256": None, "inner_k": INNER_K}
    else:
        raise Refuse(f"unsupported scheme: {scheme!r}")

    # Strict fold hygiene invariants.
    inner_k = plan_info["inner_k"]
    for fo in folds:
        f, train, test, inner = fo["fold"], fo["train"], fo["test"], fo["inner"]
        _need(not bool(np.any(train & test)), f"fold {f} overlaps train and test")
        _need(bool(np.all((inner >= 0) == train)), f"fold {f} inner folds do not match training rows")
        _need(bool(np.all(inner < inner_k)), f"fold {f} has inner fold index >= {inner_k}")
        _need(_both_classes(train, export.y), f"fold {f} training rows lack a class")
        _need(_both_classes(test, export.y), f"fold {f} test rows lack a class")
        for j in range(inner_k):
            _need(_both_classes(train & (inner != j), export.y), f"fold {f} inner train {j} lacks a class")
            _need(_both_classes(train & (inner == j), export.y), f"fold {f} inner val {j} lacks a class")

    return plan_info, folds, floor, ff_sha


# ---------------------------------------------------------------- 31-mer fold-safe encoder

def _chi_square_screen(
    X_bin: np.ndarray,
    y: np.ndarray,
    floor: int,
    max_features: int,
) -> np.ndarray:
    """Fold-safe 31-mer candidate filter + chi-square top-k (matching GRM rules)."""
    n_train, n_cols = X_bin.shape
    if n_cols == 0:
        return np.zeros(0, dtype=np.int64)
    lo, hi = floor, n_train - floor
    support = X_bin.sum(axis=0, dtype=np.int64)
    keep_mask = (support >= lo) & (support <= hi)
    cand = np.flatnonzero(keep_mask)
    if cand.size == 0:
        return np.zeros(0, dtype=np.int64)

    # Deduplicate identical training patterns keeping lowest feature index.
    sub = X_bin[:, cand]
    seen: dict[bytes, int] = {}
    unique_cols: list[int] = []
    for pos_i in range(cand.size):
        col_bytes = np.ascontiguousarray(sub[:, pos_i]).tobytes()
        if col_bytes not in seen:
            col_idx = int(cand[pos_i])
            seen[col_bytes] = col_idx
            unique_cols.append(col_idx)
    uniq = np.array(unique_cols, dtype=np.int64)

    pos_mask = y.astype(bool)
    np_ = float(np.sum(pos_mask))
    nn = float(n_train) - np_
    sub_u = X_bin[:, uniq]
    a = sub_u[pos_mask].sum(axis=0, dtype=np.float64)
    b = sub_u[~pos_mask].sum(axis=0, dtype=np.float64)
    c = np_ - a
    d = nn - b
    denom = (a + b) * (c + d) * np_ * nn
    diff = a * d - b * c
    scores = np.zeros(uniq.shape[0], dtype=np.float64)
    nz = denom != 0.0
    scores[nz] = float(n_train) * (diff[nz] * diff[nz]) / denom[nz]

    order = np.lexsort((uniq, -scores))
    chosen = uniq[order[:max_features]]
    chosen.sort()
    return chosen


class Kmer31FoldEncoder:
    """Fold-safe 31-mer screener and standardizer built on `grm export` binary matrix."""

    def __init__(
        self,
        kmer_dense: np.ndarray,
        floor: int,
        max_features: int = KMER_SCREEN_DEFAULT,
        standardize: bool = True,
    ) -> None:
        self._dense = kmer_dense
        self.floor = int(floor)
        self.max_features = int(max_features)
        self.standardize = bool(standardize)
        self.selected_cols_: np.ndarray = np.zeros(0, dtype=np.int64)
        self.mean_: np.ndarray = np.zeros(0, dtype=np.float64)
        self.scale_: np.ndarray = np.zeros(0, dtype=np.float64)
        self.fitted_: bool = False

    def fit(self, train_idx: np.ndarray, y_train: np.ndarray) -> "Kmer31FoldEncoder":
        X_tr_bin = self._dense[train_idx]
        self.selected_cols_ = _chi_square_screen(X_tr_bin, y_train, self.floor, self.max_features)
        if self.selected_cols_.size == 0:
            self.mean_ = np.zeros(0, dtype=np.float64)
            self.scale_ = np.zeros(0, dtype=np.float64)
        else:
            X_sel = X_tr_bin[:, self.selected_cols_].astype(np.float64)
            if self.standardize:
                self.mean_, self.scale_ = FT._fit_standardizer(X_sel)
            else:
                self.mean_ = np.zeros(X_sel.shape[1], dtype=np.float64)
                self.scale_ = np.ones(X_sel.shape[1], dtype=np.float64)
        self.fitted_ = True
        return self

    def transform(self, idx: np.ndarray) -> np.ndarray:
        _need(self.fitted_, "Kmer31FoldEncoder must be fitted before transform")
        if self.selected_cols_.size == 0:
            return np.zeros((idx.shape[0], 0), dtype=np.float64)
        X = self._dense[idx][:, self.selected_cols_].astype(np.float64)
        if not self.standardize:
            return X
        out = (X - self.mean_) * self.scale_
        return np.where(out == 0.0, 0.0, out)

    def state_fingerprint(self) -> str:
        _need(self.fitted_, "Kmer31FoldEncoder must be fitted before fingerprinting")
        h = hashlib.sha256()
        h.update(np.ascontiguousarray(self.selected_cols_, dtype="<i8").tobytes())
        h.update(np.ascontiguousarray(self.mean_, dtype="<f8").tobytes())
        h.update(np.ascontiguousarray(self.scale_, dtype="<f8").tobytes())
        return h.hexdigest()


# ---------------------------------------------------------------- deterministic L2 logistic regression

def class_balanced_weights(y: np.ndarray) -> np.ndarray:
    """Balanced sample weights summing to len(y), matching GRM's class weights."""
    yb = y.astype(bool)
    n = float(yb.shape[0])
    pos = float(np.sum(yb))
    neg = n - pos
    return np.where(yb, n / (2.0 * max(pos, 1.0)), n / (2.0 * max(neg, 1.0))).astype(np.float64)


def _sigmoid(z: np.ndarray) -> np.ndarray:
    zc = np.clip(z, -500.0, 500.0)
    return np.where(zc >= 0.0, 1.0 / (1.0 + np.exp(-zc)), np.exp(zc) / (1.0 + np.exp(zc)))


def _logistic_objective(
    X: np.ndarray,
    y_f: np.ndarray,
    w_norm: np.ndarray,
    coef: np.ndarray,
    intercept: float,
    lam: float,
) -> float:
    z = X @ coef + intercept
    # Stable log(1 + exp(z)) - y * z = softplus(z) - y * z
    log1p_exp = np.logaddexp(0.0, z)
    nll = float(np.sum(w_norm * (log1p_exp - y_f * z)))
    return nll + 0.5 * lam * float(np.dot(coef, coef))


def fit_l2_logistic(
    X: np.ndarray,
    y: np.ndarray,
    lam: float,
    *,
    tol: float = DEFAULT_TOL,
    max_iter: int = DEFAULT_MAX_ITER,
    gram: np.ndarray | None = None,
) -> dict[str, Any]:
    """Deterministic class-balanced L2-regularized logistic regression (pure NumPy).

    Solves min_{w, b} (1/n) sum_i omega_i log(1 + exp(-y'_i (x_i . w + b))) + (lam/2) ||w||^2
    via exact Newton-IRLS with Armijo line search. Uses a direct (p+1)x(p+1)
    primal solve when p <= n, and an exact Woodbury n x n SPD solve when p > n.
    """
    _need(lam > 0.0 and math.isfinite(lam), f"lambda must be finite and > 0, got {lam!r}")
    X = np.ascontiguousarray(X, dtype=np.float64)
    n, p = X.shape
    yb = y.astype(bool)
    y_f = yb.astype(np.float64)
    w_sample = class_balanced_weights(yb)
    w_norm = w_sample / float(n)

    if p == 0:
        return {
            "coef": np.zeros(0, dtype=np.float64),
            "intercept": 0.0,
            "converged": True,
            "n_iter": 0,
            "null_fit": True,
            "grad_norm": 0.0,
        }

    # Check if every column has zero gradient at (0, 0) (e.g. all-zero / constant matrix).
    p0 = np.full(n, 0.5, dtype=np.float64)
    g0_w = X.T @ (w_norm * (p0 - y_f))
    if float(np.max(np.abs(g0_w))) <= ZERO_COEF_TOL:
        return {
            "coef": np.zeros(p, dtype=np.float64),
            "intercept": 0.0,
            "converged": True,
            "n_iter": 0,
            "null_fit": True,
            "grad_norm": 0.0,
        }

    coef = np.zeros(p, dtype=np.float64)
    intercept = 0.0
    use_woodbury = p > n
    K = gram if (use_woodbury and gram is not None) else (X @ X.T if use_woodbury else None)

    converged = False
    n_iter = 0
    grad_norm = float("inf")
    obj = _logistic_objective(X, y_f, w_norm, coef, intercept, lam)

    for it in range(1, max_iter + 1):
        n_iter = it
        z = X @ coef + intercept
        prob = _sigmoid(z)
        resid = w_norm * (prob - y_f)
        g_w = X.T @ resid + lam * coef
        g_b = float(np.sum(resid))
        grad_norm = float(max(np.max(np.abs(g_w)), abs(g_b)))
        if grad_norm <= tol:
            converged = True
            break

        s = np.maximum(w_norm * prob * (1.0 - prob), 1e-16)
        if not use_woodbury:
            Xs = X * s[:, None]
            H = np.empty((p + 1, p + 1), dtype=np.float64)
            H[:p, :p] = X.T @ Xs
            H[:p, :p].flat[:: p + 1] += lam
            col_b = X.T @ s
            H[:p, p] = col_b
            H[p, :p] = col_b
            H[p, p] = float(np.sum(s))
            rhs = np.empty(p + 1, dtype=np.float64)
            rhs[:p] = -g_w
            rhs[p] = -g_b
            step = np.linalg.solve(H, rhs)
            dw = step[:p]
            db = float(step[p])
        else:
            # Exact Woodbury solve for H_ww = X^T diag(s) X + lam I_p.
            sqrt_s = np.sqrt(s)
            assert K is not None
            M = (K * (sqrt_s[:, None] * sqrt_s[None, :])) / lam
            M.flat[:: n + 1] += 1.0
            h_wb = X.T @ s
            h_bb = float(np.sum(s))
            # Solve H_ww [u1, u2] = [-g_w, h_wb] simultaneously.
            R = np.column_stack((-g_w, h_wb))  # p x 2
            BR = sqrt_s[:, None] * (X @ R)  # n x 2
            Z_sol = np.linalg.solve(M, BR)  # n x 2
            U = (R - (X.T @ (sqrt_s[:, None] * Z_sol)) / lam) / lam  # p x 2
            u1 = U[:, 0]
            u2 = U[:, 1]
            schur = max(h_bb - float(np.dot(h_wb, u2)), 1e-16)
            db = (-g_b - float(np.dot(h_wb, u1))) / schur
            dw = u1 - db * u2

        directional = float(np.dot(g_w, dw)) + g_b * db
        if directional >= 0.0:
            # Fallback to steepest descent if numerical noise flips curvature.
            dw = -g_w
            db = -g_b
            directional = float(np.dot(g_w, dw)) + g_b * db

        alpha = 1.0
        accepted = False
        for _ in range(30):
            cand_w = coef + alpha * dw
            cand_b = intercept + alpha * db
            cand_obj = _logistic_objective(X, y_f, w_norm, cand_w, cand_b, lam)
            if cand_obj <= obj + 1e-4 * alpha * directional + 1e-15:
                coef = cand_w
                intercept = cand_b
                obj = cand_obj
                accepted = True
                break
            alpha *= 0.5
        if not accepted:
            break

    # Final gradient check after last step.
    z = X @ coef + intercept
    prob = _sigmoid(z)
    resid = w_norm * (prob - y_f)
    g_w = X.T @ resid + lam * coef
    g_b = float(np.sum(resid))
    grad_norm = float(max(np.max(np.abs(g_w)), abs(g_b)))
    if grad_norm <= tol:
        converged = True

    coef = np.where(np.abs(coef) <= ZERO_COEF_TOL, 0.0, coef)
    null_fit = bool(np.all(coef == 0.0))
    return {
        "coef": coef,
        "intercept": float(intercept),
        "converged": converged,
        "n_iter": n_iter,
        "null_fit": null_fit,
        "grad_norm": grad_norm,
    }


def predict_l2_logistic(X: np.ndarray, coef: np.ndarray, intercept: float) -> np.ndarray:
    return _sigmoid(np.ascontiguousarray(X, dtype=np.float64) @ coef + float(intercept))


def inner_balanced_log_loss(scores: np.ndarray, y: np.ndarray) -> float:
    """Class-balanced mean log loss for inner-fold hyperparameter selection."""
    yb = y.astype(bool)
    n = yb.shape[0]
    pos = float(np.sum(yb))
    p = np.clip(scores, INNER_LOSS_CLIP, 1.0 - INNER_LOSS_CLIP)
    terms = np.where(
        yb,
        -np.log(p) * n / (2.0 * max(pos, 1.0)),
        -np.log(1.0 - p) * n / (2.0 * max(n - pos, 1.0)),
    )
    return float(np.sum(terms)) / max(n, 1)


def fold_calls(scores: np.ndarray, y: np.ndarray) -> dict[str, Any]:
    """Confusion counts at 0.5 cutoff and clipped log loss sum (matches src/cutoff.rs)."""
    yb = y.astype(bool)
    p = np.clip(scores, LOSS_CLIP, 1.0 - LOSS_CLIP)
    pred = scores >= CUTOFF
    loss_sum = float(np.sum(np.where(yb, -np.log(p), -np.log(1.0 - p))))
    return {
        "n": int(yb.shape[0]),
        "tp": int(np.sum(pred & yb)),
        "fp": int(np.sum(pred & ~yb)),
        "tn": int(np.sum(~pred & ~yb)),
        "fn": int(np.sum(~pred & yb)),
        "loss_sum": loss_sum,
    }


def pooled_metrics(fold_rows: Sequence[dict[str, Any]]) -> dict[str, Any]:
    """Pool outer-fold calls into dataset-level classification and calibration metrics."""
    n = sum(r["calls"]["n"] for r in fold_rows)
    tp = sum(r["calls"]["tp"] for r in fold_rows)
    fp = sum(r["calls"]["fp"] for r in fold_rows)
    tn = sum(r["calls"]["tn"] for r in fold_rows)
    fn = sum(r["calls"]["fn"] for r in fold_rows)
    loss_sum = sum(r["calls"]["loss_sum"] for r in fold_rows)
    n_pos = tp + fn
    n_neg = tn + fp
    recall = tp / n_pos if n_pos > 0 else 0.0
    specificity = tn / n_neg if n_neg > 0 else 0.0
    ba = 0.5 * (recall + specificity)
    vme = fn / n_pos if n_pos > 0 else 0.0
    me = fp / n_neg if n_neg > 0 else 0.0
    log_loss = loss_sum / n if n > 0 else 0.0
    return {
        "n": n,
        "tp": tp,
        "fp": fp,
        "tn": tn,
        "fn": fn,
        "balanced_accuracy": ba,
        "recall": recall,
        "specificity": specificity,
        "vme": vme,
        "me": me,
        "log_loss": log_loss,
    }


# ---------------------------------------------------------------- representation arm specs

class ArmSpec:
    """One representation arm in a fixed-learner comparison.

    `build_fold_pair(train_idx, eval_idx, y_train)` must fit strictly on
    `train_idx` (never reading `eval_idx` labels or values) and return
    `(X_train, X_eval, fitted_state)`.
    """

    def __init__(
        self,
        name: str,
        kind: str,
        build_fn: Callable[[np.ndarray, np.ndarray, np.ndarray], tuple[np.ndarray, np.ndarray, Any]],
        *,
        exposure_audits: list[dict[str, Any]] | None = None,
    ) -> None:
        _need(isinstance(name, str) and bool(name.strip()), "arm name must be non-empty")
        self.name = name.strip()
        self.kind = kind
        self._build_fn = build_fn
        self.exposure_audits: list[dict[str, Any]] = list(exposure_audits or ())

    @property
    def claim_ceiling(self) -> str:
        if not self.exposure_audits:
            return "retrospective_internal_split"
        return EA.min_claim_ceiling([a["effective_claim_ceiling"] for a in self.exposure_audits])

    def fit_transform_pair(
        self,
        train_idx: np.ndarray,
        eval_idx: np.ndarray,
        y_train: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray, Any]:
        X_tr, X_ev, state = self._build_fn(train_idx, eval_idx, y_train)
        _need(X_tr.ndim == 2 and X_ev.ndim == 2, f"arm {self.name}: feature matrices must be 2D")
        _need(X_tr.shape[0] == train_idx.shape[0], f"arm {self.name}: X_train row mismatch")
        _need(X_ev.shape[0] == eval_idx.shape[0], f"arm {self.name}: X_eval row mismatch")
        _need(X_tr.shape[1] == X_ev.shape[1], f"arm {self.name}: X_train/X_eval column mismatch")
        _need(bool(np.all(np.isfinite(X_tr))) and bool(np.all(np.isfinite(X_ev))),
              f"arm {self.name}: non-finite feature values")
        return X_tr, X_ev, state


def make_kmer31_arm(
    export: ExportDataset,
    floor: int,
    *,
    name: str = "kmer31",
    max_features: int = KMER_SCREEN_DEFAULT,
    standardize: bool = True,
) -> ArmSpec:
    def _build(train_idx: np.ndarray, eval_idx: np.ndarray, y_train: np.ndarray):
        enc = Kmer31FoldEncoder(export.kmer_dense, floor=floor, max_features=max_features, standardize=standardize)
        enc.fit(train_idx, y_train)
        return enc.transform(train_idx), enc.transform(eval_idx), enc

    return ArmSpec(name=name, kind="kmer31", build_fn=_build)


def make_protein_variant_arm(
    isolates_aligned: Sequence[Sequence[dict[str, Any]]],
    *,
    name: str = "protein_variants",
    target_genes: Sequence[str] | None = None,
    frozen_catalogue_variants: Sequence[tuple[str, str]] | None = None,
    vocabulary_mode: str = "train_only",
    min_support: int = 1,
    max_support_margin: int = 0,
    include_synonymous: bool = False,
    include_unseen_side_feature: bool = True,
) -> ArmSpec:
    def _build(train_idx: np.ndarray, eval_idx: np.ndarray, y_train: np.ndarray):
        del y_train
        enc = FT.ProteinVariantEncoder(
            target_genes=target_genes,
            frozen_catalogue_variants=frozen_catalogue_variants,
            vocabulary_mode=vocabulary_mode,
            min_support=min_support,
            max_support_margin=max_support_margin,
            include_synonymous=include_synonymous,
            include_unseen_side_feature=include_unseen_side_feature,
        )
        tr_isos = [isolates_aligned[int(i)] for i in train_idx]
        ev_isos = [isolates_aligned[int(i)] for i in eval_idx]
        X_tr = enc.fit_transform(tr_isos)
        X_ev = enc.transform(ev_isos)
        return X_tr, X_ev, enc

    return ArmSpec(name=name, kind="protein_variants", build_fn=_build)


def make_burden_arm(
    isolates_aligned: Sequence[Sequence[dict[str, Any]]],
    target_genes: Sequence[str],
    *,
    regions: Sequence[dict[str, Any]] | None = None,
    standardize: bool = True,
    name: str = "gene_burden",
) -> ArmSpec:
    def _build(train_idx: np.ndarray, eval_idx: np.ndarray, y_train: np.ndarray):
        del y_train
        enc = FT.BurdenEncoder(target_genes=target_genes, regions=regions, standardize=standardize)
        tr_isos = [isolates_aligned[int(i)] for i in train_idx]
        ev_isos = [isolates_aligned[int(i)] for i in eval_idx]
        X_tr = enc.fit_transform(tr_isos)
        X_ev = enc.transform(ev_isos)
        return X_tr, X_ev, enc

    return ArmSpec(name=name, kind="gene_burden", build_fn=_build)


def make_foundation_arm(
    adapter: FT.FoundationScoreAdapter,
    *,
    name: str | None = None,
) -> ArmSpec:
    _need(isinstance(adapter, FT.FoundationScoreAdapter), "foundation arm requires a FoundationScoreAdapter")
    arm_name = name or f"foundation_{adapter.exposure_audit_['model_family']}"

    def _build(train_idx: np.ndarray, eval_idx: np.ndarray, y_train: np.ndarray):
        del y_train
        enc = adapter.clone_unfitted()
        X_tr = enc.fit_transform(train_idx)
        X_ev = enc.transform(eval_idx)
        return X_tr, X_ev, enc

    return ArmSpec(
        name=arm_name,
        kind="foundation",
        build_fn=_build,
        exposure_audits=[adapter.exposure_audit_],
    )


def make_concat_arm(name: str, arms: Sequence[ArmSpec]) -> ArmSpec:
    """Concatenate two or more representation arms fold-safely under one learner."""
    _need(len(arms) >= 2, "concat arm requires at least two component arms")
    audits: list[dict[str, Any]] = []
    for a in arms:
        audits.extend(a.exposure_audits)

    def _build(train_idx: np.ndarray, eval_idx: np.ndarray, y_train: np.ndarray):
        tr_parts: list[np.ndarray] = []
        ev_parts: list[np.ndarray] = []
        states: dict[str, Any] = {}
        for a in arms:
            X_tr, X_ev, st = a.fit_transform_pair(train_idx, eval_idx, y_train)
            tr_parts.append(X_tr)
            ev_parts.append(X_ev)
            states[a.name] = st
        X_tr_all = np.concatenate(tr_parts, axis=1) if tr_parts else np.zeros((train_idx.shape[0], 0), dtype=np.float64)
        X_ev_all = np.concatenate(ev_parts, axis=1) if ev_parts else np.zeros((eval_idx.shape[0], 0), dtype=np.float64)
        return X_tr_all, X_ev_all, states

    return ArmSpec(name=name, kind="concat", build_fn=_build, exposure_audits=audits)


# ---------------------------------------------------------------- comparator core

def _evaluate_arm_fold(
    arm: ArmSpec,
    y: np.ndarray,
    fo: dict[str, Any],
    inner_k: int,
    lambda_grid: Sequence[float],
    tol: float,
    max_iter: int,
) -> tuple[dict[str, Any], np.ndarray, dict[str, Any]]:
    """Run inner-fold lambda selection + outer-train refit strictly inside outer fold `fo`."""
    train_mask = fo["train"]
    test_mask = fo["test"]
    inner = fo["inner"]
    tr_idx = np.flatnonzero(train_mask)
    te_idx = np.flatnonzero(test_mask)

    # Inner-fold CV strictly on outer training rows; every encoder refits inside each inner train split.
    n_lam = len(lambda_grid)
    pooled_inner_preds = np.full((n_lam, tr_idx.shape[0]), np.nan, dtype=np.float64)
    local_pos = {int(g_i): l_i for l_i, g_i in enumerate(tr_idx)}
    inner_failed = [False] * n_lam

    for j in range(inner_k):
        in_tr = np.flatnonzero(train_mask & (inner != j))
        in_va = np.flatnonzero(train_mask & (inner == j))
        X_in_tr, X_in_va, _ = arm.fit_transform_pair(in_tr, in_va, y[in_tr])
        gram = (X_in_tr @ X_in_tr.T) if X_in_tr.shape[1] > X_in_tr.shape[0] else None
        va_local = np.array([local_pos[int(g_i)] for g_i in in_va], dtype=np.int64)
        prev_rate = float(np.sum(y[in_tr])) / float(in_tr.shape[0])
        for li, lam in enumerate(lambda_grid):
            fit_j = fit_l2_logistic(X_in_tr, y[in_tr], lam, tol=tol, max_iter=max_iter, gram=gram)
            if not fit_j["converged"]:
                inner_failed[li] = True
            if fit_j["null_fit"]:
                preds = np.full(in_va.shape[0], prev_rate, dtype=np.float64)
            else:
                preds = predict_l2_logistic(X_in_va, fit_j["coef"], fit_j["intercept"])
            pooled_inner_preds[li, va_local] = preds

    # Pick lambda minimizing pooled inner class-balanced log loss; ties go to larger lambda.
    best_choice: tuple[tuple[float, float], int] | None = None
    for li, lam in enumerate(lambda_grid):
        if inner_failed[li]:
            continue
        loss = round(inner_balanced_log_loss(pooled_inner_preds[li], y[tr_idx]), 12)
        key = (loss, -round(float(lam), 12))
        if best_choice is None or key < best_choice[0]:
            best_choice = (key, li)

    # Fallback if every inner grid point failed convergence: still pick lowest loss deterministically.
    if best_choice is None:
        for li, lam in enumerate(lambda_grid):
            loss = round(inner_balanced_log_loss(pooled_inner_preds[li], y[tr_idx]), 12)
            key = (loss, -round(float(lam), 12))
            if best_choice is None or key < best_choice[0]:
                best_choice = (key, li)

    assert best_choice is not None
    picked_lam = float(lambda_grid[best_choice[1]])

    # Outer-fold refit strictly on outer training rows, then predict outer test rows.
    X_tr, X_te, fitted_encoder = arm.fit_transform_pair(tr_idx, te_idx, y[tr_idx])
    refit = fit_l2_logistic(X_tr, y[tr_idx], picked_lam, tol=tol, max_iter=max_iter)
    if refit["null_fit"]:
        te_scores = np.full(te_idx.shape[0], float(np.sum(y[tr_idx])) / float(tr_idx.shape[0]), dtype=np.float64)
        constant = 1
    else:
        te_scores = predict_l2_logistic(X_te, refit["coef"], refit["intercept"])
        constant = 0

    fold_row = {
        "fold": int(fo["fold"]),
        "train": int(tr_idx.shape[0]),
        "train_resistant": int(np.sum(y[tr_idx])),
        "test": int(te_idx.shape[0]),
        "test_resistant": int(np.sum(y[te_idx])),
        "n_features": int(X_tr.shape[1]),
        "picked_lambda": picked_lam,
        "converged": int(refit["converged"]),
        "null_fit": int(refit["null_fit"]),
        "constant": constant,
        "calls": fold_calls(te_scores, y[te_idx]),
    }
    diagnostics = {
        "fold": int(fo["fold"]),
        "picked_lambda": picked_lam,
        "encoder": fitted_encoder,
        "coef": refit["coef"].copy(),
        "intercept": float(refit["intercept"]),
        "converged": bool(refit["converged"]),
        "null_fit": bool(refit["null_fit"]),
    }
    return fold_row, te_scores, diagnostics


def compare_representations(
    export: ExportDataset,
    scheme: str,
    seed: int,
    arms: Sequence[ArmSpec],
    *,
    fold_file: str | Path | None = None,
    baseline_arm: str | None = None,
    lambda_grid: Sequence[float] = DEFAULT_LAMBDA_GRID,
    tol: float = DEFAULT_TOL,
    max_iter: int = DEFAULT_MAX_ITER,
) -> tuple[dict[str, Any], dict[str, list[dict[str, Any]]]]:
    """Run fixed-learner comparison across `arms` on identical folds and support.

    Returns `(report, diagnostics)`. `report` is the zero-ID `grm-rep-compare-v1`
    JSON-ready dict; `diagnostics` holds in-memory fitted objects per arm/fold
    for unit and leakage verification and is never written to the report.
    """
    _need(len(arms) > 0, "at least one representation arm is required")
    arm_names = [a.name for a in arms]
    _need(len(set(arm_names)) == len(arm_names), "duplicate arm names")
    lam_grid = tuple(float(x) for x in lambda_grid)
    _need(len(lam_grid) > 0 and all(x > 0.0 and math.isfinite(x) for x in lam_grid),
          "lambda_grid must contain finite values > 0")

    plan_info, folds, floor, ff_sha = load_folds(export, scheme, seed, fold_file=fold_file)

    if baseline_arm is None:
        base_name = "kmer31" if "kmer31" in arm_names else arm_names[0]
    else:
        _need(baseline_arm in arm_names, f"baseline_arm {baseline_arm!r} not in arms {arm_names}")
        base_name = baseline_arm

    tested_mask = np.zeros(export.n, dtype=bool)
    for fo in folds:
        tested_mask |= fo["test"]

    arm_reports: list[dict[str, Any]] = []
    arm_diagnostics: dict[str, list[dict[str, Any]]] = {}
    arm_pooled: dict[str, dict[str, Any]] = {}

    for arm in arms:
        fold_rows: list[dict[str, Any]] = []
        fold_diags: list[dict[str, Any]] = []
        for fo in folds:
            row, _, diag = _evaluate_arm_fold(
                arm, export.y, fo, plan_info["inner_k"], lam_grid, tol, max_iter
            )
            fold_rows.append(row)
            fold_diags.append(diag)
        pooled = pooled_metrics(fold_rows)
        arm_pooled[arm.name] = pooled
        arm_diagnostics[arm.name] = fold_diags
        arm_reports.append({
            "arm": arm.name,
            "kind": arm.kind,
            "claim_ceiling": arm.claim_ceiling,
            "exposure_audits": [EA.summarize_audit_for_report(a) for a in arm.exposure_audits],
            "pooled_metrics": pooled,
            "constant_folds": sum(r["constant"] for r in fold_rows),
            "null_fit_folds": sum(r["null_fit"] for r in fold_rows),
            "unconverged_folds": sum(1 - r["converged"] for r in fold_rows),
            "folds": fold_rows,
        })

    base_metrics = arm_pooled[base_name]
    pairwise: list[dict[str, Any]] = []
    for arm in arms:
        if arm.name == base_name:
            continue
        m = arm_pooled[arm.name]
        _need(m["n"] == base_metrics["n"], f"arm {arm.name} tested support differs from baseline {base_name}")
        pairwise.append({
            "arm": arm.name,
            "baseline_arm": base_name,
            "common_tested_n": m["n"],
            "delta_balanced_accuracy": round(m["balanced_accuracy"] - base_metrics["balanced_accuracy"], 12),
            "delta_recall": round(m["recall"] - base_metrics["recall"], 12),
            "delta_specificity": round(m["specificity"] - base_metrics["specificity"], 12),
            "delta_vme": round(m["vme"] - base_metrics["vme"], 12),
            "delta_me": round(m["me"] - base_metrics["me"], 12),
            "delta_log_loss": round(m["log_loss"] - base_metrics["log_loss"], 12),
        })

    report = {
        "format": FORMAT,
        "evaluation_revision": REVISION,
        "dataset": export.meta["id"],
        "scheme": scheme,
        "seed": int(seed),
        "input_fingerprint": export.meta["input_fingerprint"],
        "plan_sha256": plan_info["sha256"],
        "matrix_sha256": export.meta["matrix_sha256"],
        "fold_file_sha256": ff_sha,
        "samples": export.n,
        "resistant": int(np.sum(export.y)),
        "tested": int(np.sum(tested_mask)),
        "tested_resistant": int(np.sum(tested_mask & export.y)),
        "min_support": floor,
        "baseline_arm": base_name,
        "learner_control": {
            "learner": "l2_logistic_newton",
            "class_weight": "balanced_on_training_fold",
            "lambda_grid": list(lam_grid),
            "inner_folds": plan_info["inner_k"],
            "cutoff": CUTOFF,
            "loss_clip": LOSS_CLIP,
            "tol": tol,
            "max_iter": max_iter,
            "tie_rule": "larger_lambda",
        },
        "scope_note": "Retrospective comparison inside one database; not external or clinical validation.",
        "arms": arm_reports,
        "pairwise_vs_baseline": pairwise,
    }
    validate_report(report, forbidden_ids=export.samples)
    return report, arm_diagnostics


def validate_report(report: dict[str, Any], forbidden_ids: Sequence[str] | None = None) -> None:
    """Enforce exact `grm-rep-compare-v1` schema keys and zero-ID hygiene."""
    _need(set(report.keys()) == REPORT_KEYS, f"report keys mismatch: {sorted(set(report) ^ REPORT_KEYS)}")
    _need(report["format"] == FORMAT, f"unexpected format: {report['format']!r}")
    _need(report["evaluation_revision"] == REVISION, f"unexpected revision: {report['evaluation_revision']!r}")
    for arm in report["arms"]:
        _need(set(arm.keys()) == ARM_KEYS, f"arm keys mismatch: {sorted(set(arm) ^ ARM_KEYS)}")
        _need(set(arm["pooled_metrics"].keys()) == POOLED_KEYS, "pooled_metrics keys mismatch")
        for row in arm["folds"]:
            _need(set(row.keys()) == FOLD_KEYS, f"fold keys mismatch: {sorted(set(row) ^ FOLD_KEYS)}")
            _need(set(row["calls"].keys()) == CALL_KEYS, "fold calls keys mismatch")
            for k, v in row.items():
                if k != "calls":
                    _need(isinstance(v, (int, float)) and not isinstance(v, bool),
                          f"fold field {k!r} must be numeric")
    for pw in report["pairwise_vs_baseline"]:
        _need(set(pw.keys()) == PAIRWISE_KEYS, f"pairwise keys mismatch: {sorted(set(pw) ^ PAIRWISE_KEYS)}")

    if forbidden_ids:
        id_set = set(forbidden_ids)
        text = dump_report(report)
        hits = [m.group(0) for m in GENOME_ID_TOKEN_RE.finditer(text) if m.group(0) in id_set]
        _need(not hits, f"report leaked {len(hits)} genome id token(s)")


def dump_report(report: dict[str, Any]) -> str:
    """Deterministic canonical JSON serialization (no timestamps, sorted keys)."""
    return json.dumps(report, sort_keys=True, indent=2, ensure_ascii=True) + "\n"


def write_atomic(path: str | Path, text: str) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_name(p.name + ".tmp")
    tmp.write_text(text, encoding="utf-8")
    os.replace(tmp, p)


# ---------------------------------------------------------------- CLI

def _load_variants_manifest(path: str | Path, export: ExportDataset) -> tuple[list[list[dict[str, Any]]], dict[str, Any]]:
    doc = json.loads(Path(path).read_text(encoding="utf-8"))
    _need(isinstance(doc, dict) and doc.get("dataset") == export.meta["id"],
          "variants manifest dataset differs from export")
    by_sample = doc.get("isolates")
    _need(isinstance(by_sample, dict), "variants manifest requires 'isolates' mapping")
    _need(set(by_sample.keys()) == set(export.samples),
          "variants manifest sample set does not match export cohort")
    aligned = [list(by_sample[s]) for s in export.samples]
    return aligned, doc


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="compare.py", description="Track B fixed-learner representation comparator.")
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run", help="Compare representation arms on one exported dataset.")
    r.add_argument("--export", required=True, help="grm export directory (dataset.json + matrix.bin)")
    r.add_argument("--scheme", required=True, choices=DRAWN_SCHEMES + (STUDY_SCHEME,))
    r.add_argument("--seed", required=True, type=int)
    r.add_argument("--fold-file", help="grm-fold-file-v1 JSON (required for study-whole)")
    r.add_argument("--variants", help="JSON manifest of per-isolate protein-aligned target variants")
    r.add_argument("--foundation", action="append", default=[],
                   help="Path to grm-rep-foundation-v1 manifest (may be repeated)")
    r.add_argument("--pca-components", type=int, default=None,
                   help="Optional fold-safe PCA dimension for foundation arms")
    r.add_argument("--out", required=True, help="Output path for grm-rep-compare-v1 JSON")
    args = ap.parse_args(argv)

    try:
        export = ExportDataset(args.export)
        _, _, floor, _ = load_folds(export, args.scheme, args.seed, fold_file=args.fold_file)
        fold_doc = json.loads(Path(args.fold_file).read_text(encoding="utf-8")) if args.fold_file else None
        arms: list[ArmSpec] = [make_kmer31_arm(export, floor=floor)]
        if args.variants:
            aligned_vars, vdoc = _load_variants_manifest(args.variants, export)
            target_genes = vdoc.get("target_genes")
            regions = vdoc.get("regions")
            arms.append(make_protein_variant_arm(aligned_vars, target_genes=target_genes))
            if target_genes:
                arms.append(make_burden_arm(aligned_vars, target_genes=target_genes, regions=regions))
        for fpath in args.foundation:
            adapter = FT.FoundationScoreAdapter(
                fpath,
                export_meta=export.meta,
                fold_doc=fold_doc,
                pca_components=args.pca_components,
            )
            arms.append(make_foundation_arm(adapter))
        report, _ = compare_representations(
            export,
            args.scheme,
            args.seed,
            arms,
            fold_file=args.fold_file,
        )
    except (Refuse, EA.AuditError, FT.FeatureError, OSError, json.JSONDecodeError) as exc:
        print(f"compare.py: refused: {exc}", file=sys.stderr)
        return 1

    write_atomic(args.out, dump_report(report))
    print(f"compare.py: wrote {len(report['arms'])} arms on {report['dataset']} ({report['scheme']} s{report['seed']})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
