"""Trainable AMR adapter execution on GRM exports and grm-fold-file-v1 splits.

Runs per (dataset, split, outer_fold) consuming:
- `grm export` directories (`dataset.json` + `matrix.bin`) on `random` or `lineage`, or
- `grm-fold-file-v1` files (`study-whole`, `time`, or other reviewed fold-file schemes).

Guarantees:
- Fold parity enforcement via `scripts/c7-outside/plans.py` (shape by class AND
  membership by genome SHA-256, plus release `input_fingerprint` at seed 42).
- Every learned preprocessing step (support filtering, chi-square feature
  screening, column scaling, hyperparameter tuning over `inner` / `holdouts`)
  fits strictly inside outer training folds; outer test labels are physically
  withheld from every fitter.
- External cohorts (`selection_role='evaluation_only'`) are refused for tuning
  or model selection.
- Emits per-fold counts `{fold, n, tp, fp, tn, fn, no_call_r, no_call_s, i_r, i_s,
  failure_r, failure_s, loss_sum}` with zero genome IDs, plus optional local-only
  out-of-fold score vectors strictly in a scratch directory (never committed).
"""

from __future__ import annotations

import importlib.util
import json
import math
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable

import numpy as np

from scripts.bench.protocol import (
    AdapterSpec,
    BenchError,
    DRAWN_SCHEMES,
    EVALUATION_REVISION,
    MODE_TRAINABLE,
    SELECTION_ROLE_INTERNAL_CV,
    assert_selection_allowed,
    canonical_json_bytes,
    load_fold_file,
    require,
)

SCRIPTS_DIR = Path(__file__).resolve().parent.parent
C7_DIR = SCRIPTS_DIR / "c7-outside"
REPO_ROOT = SCRIPTS_DIR.parent

LOSS_CLIP = 1e-4  # src/cutoff.rs
CUTOFF = 0.5
INNER_K = 3
GBT_HOLDOUT = 0.15


def _load_c7_module(name: str, filename: str) -> Any:
    path = C7_DIR / filename
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise BenchError(f"cannot load {path}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


_OUTSIDE = None
_PLANS = None


def _get_outside() -> Any:
    global _OUTSIDE
    if _OUTSIDE is None:
        _OUTSIDE = _load_c7_module("bench_c7_outside", "outside.py")
    return _OUTSIDE


def _get_plans() -> Any:
    global _PLANS
    if _PLANS is None:
        _PLANS = _load_c7_module("bench_c7_plans", "plans.py")
    return _PLANS


# ------------------------------------------------------------------ fold parity (plans.py)


def verify_fold_parity(
    export_dir: Path | str,
    seed: int,
    *,
    scheme: str = "lineage",
    reference_result_doc: dict[str, Any] | None = None,
    reference_report_doc: dict[str, Any] | None = None,
    release_report_doc: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Enforce exact fold parity via scripts/c7-outside/plans.py.

    Verifies:
    - The export's plan at `seed` is valid and its canonical SHA-256 holds.
    - If `reference_result_doc` is provided: outer fold shape by class matches
      `plans.result_shape(reference_result_doc)`.
    - If `reference_report_doc` is provided: outer fold membership by genome
      SHA-256 matches `plans.report_members(reference_report_doc, seed)` (so
      swapping even two same-class genomes between folds fails).
    - If `release_report_doc` is provided (at seed 42): `input_fingerprint` matches.

    Never emits genome IDs or member hashes in outputs or errors.
    """
    outside = _get_outside()
    plans = _get_plans()
    try:
        ex = outside.Export(str(export_dir))
        plan, folds = outside.export_plan(ex, scheme, seed)
    except (outside.Refuse, OSError, KeyError, ValueError) as exc:
        raise BenchError(f"fold parity export load failed: {exc}") from exc

    checks_run: dict[str, bool] = {"export_plan_valid": True}

    if scheme == "lineage":
        exp_shape = plans.export_shape(ex, seed)
        exp_members = plans.export_members(ex, seed)
    else:
        if folds is None:
            exp_shape = ("unavailable",)
            exp_members = ("unavailable",)
        else:
            y = ex.y
            shape = tuple(
                (int(f["train"].sum()), int((f["train"] & y).sum()), int(f["test"].sum()), int((f["test"] & y).sum()))
                for f in folds
            )
            tested = sum(int(f["test"].sum()) for f in folds)
            exp_shape = shape + ((tested, sum(int((f["test"] & y).sum()) for f in folds)),)
            exp_members = plans.members(ex.meta["samples"], plan["outer"])

    if reference_result_doc is not None:
        want_shape = plans.result_shape(reference_result_doc)
        same_shape = exp_shape == want_shape
        checks_run["folds_by_class_match"] = same_shape
        require(same_shape, f"fold parity failed for {ex.meta['id']} ({scheme} seed {seed}): fold counts by class differ")

    if reference_report_doc is not None:
        try:
            if scheme == "lineage":
                want_members = plans.report_members(reference_report_doc, seed)
            else:
                require(
                    reference_report_doc["eval_config"]["seed"] == seed,
                    f"reference report seed {reference_report_doc['eval_config']['seed']} != {seed}",
                )
                sp_list = [x for x in reference_report_doc["splits"] if x["scheme"] == scheme]
                require(len(sp_list) == 1, f"reference report has {len(sp_list)} {scheme} splits")
                if sp_list[0].get("unavailable_reason") is not None:
                    want_members = ("unavailable",)
                else:
                    ids = [x["id"] for x in reference_report_doc["samples"]]
                    want_members = plans.members(ids, sp_list[0]["test_fold"])
        except (KeyError, TypeError, ValueError) as exc:
            raise BenchError(f"fold parity reference report invalid: {exc}") from exc
        same_members = exp_members == want_members
        checks_run["folds_by_genome_match"] = same_members
        require(
            same_members,
            f"fold parity failed for {ex.meta['id']} ({scheme} seed {seed}): outer fold genome membership differs",
        )

    if release_report_doc is not None:
        same_fp = ex.meta["input_fingerprint"] == release_report_doc.get("input_fingerprint")
        checks_run["input_fingerprint_match"] = same_fp
        require(
            same_fp,
            f"fold parity failed for {ex.meta['id']}: export input_fingerprint differs from release report",
        )

    return {
        "dataset": ex.meta["id"],
        "scheme": scheme,
        "seed": seed,
        "plan_sha256": plan.get("sha256"),
        "input_fingerprint": ex.meta["input_fingerprint"],
        "checks": checks_run,
    }


# ------------------------------------------------------------------ fold-file plans (study-whole, time, etc.)


def build_fold_file_plan(ex: Any, fold_file_path: Path | str, seed: int) -> tuple[dict[str, Any], list[dict[str, Any]], str]:
    """Build outer folds from any valid grm-fold-file-v1 (study-whole, time, etc.)
    with inner folds and boosting holdouts drawn strictly inside each outer training fold."""
    outside = _get_outside()
    dataset = ex.meta["id"]
    ff = load_fold_file(fold_file_path, expected_dataset=dataset, sample_ids=ex.meta["samples"])

    col = {g: i for i, g in enumerate(ex.meta["samples"])}
    roles = np.full(ex.n, -3, dtype=np.int64)
    for f_idx, fold_ids in enumerate(ff.test):
        for g in fold_ids:
            roles[col[g]] = f_idx
    for g in ff.train:
        roles[col[g]] = -1
    for g in ff.unused:
        roles[col[g]] = -2
    require(bool(np.all(roles != -3)), f"{dataset}: fold file left a sample unassigned")

    used = roles != -2
    floor = max(3, int(np.sum(used)) // 200)
    dups = [[i for i in grp if used[i]] for grp in ex.meta.get("duplicates", [])]
    dups = [grp for grp in dups if len(grp) > 1]
    groups = outside.tuning_groups(ff.scheme, None, dups, ex.n)

    folds: list[dict[str, Any]] = []
    for f_idx in range(len(ff.test)):
        test = roles == f_idx
        train = (roles == -1) | ((roles >= 0) & (roles != f_idx))
        require(
            outside.both(test, ex.y) and outside.both(train, ex.y),
            f"{dataset} {ff.scheme} outer fold {f_idx} lacks one of the two classes",
        )
        for grp in dups:
            require(
                not (any(test[i] for i in grp) and any(train[i] for i in grp)),
                f"{dataset} {ff.scheme} outer fold {f_idx} splits a duplicate group",
            )
        inner = outside.training_folds(ex.y, train, groups, INNER_K, outside.draw_rng(dataset, ff.scheme, f_idx, "inner"))
        fit, hold = outside.tuning_holdout(train, ex.y, groups, GBT_HOLDOUT, outside.draw_rng(dataset, ff.scheme, f_idx, "gbt"))
        folds.append({
            "fold": f_idx,
            "train": train,
            "test": test,
            "hold": hold,
            "fit": fit,
            "inner": inner,
            "seed": outside.fold_seed(seed, f_idx),
            "screened_rows": None,
            "candidates": None,
        })
    try:
        outside.check_folds(folds, ex.y, ex.n)
    except outside.Refuse as exc:
        raise BenchError(str(exc)) from exc

    return (
        {"scheme": ff.scheme, "sha256": None, "unavailable_reason": None, "floor": floor},
        folds,
        ff.sha256,
    )


# ------------------------------------------------------------------ strict training-only fold context


@dataclass(frozen=True)
class FoldFitContext:
    """Strictly training-side context passed to a trainable fitter.

    Contains ZERO outer test labels and ZERO genome IDs, making outer-test
    label leakage structurally impossible.
    """

    tool: str
    fold: int
    seed: int
    X_train: np.ndarray  # shape (n_train, n_features)
    y_train: np.ndarray  # shape (n_train,) bool
    X_test: np.ndarray   # shape (n_test, n_features)
    fit_mask: np.ndarray   # shape (n_train,) bool (inside X_train)
    hold_mask: np.ndarray  # shape (n_train,) bool (inside X_train)
    inner_folds: np.ndarray  # shape (n_train,) int in 0..INNER_K-1
    sample_weights_fit: np.ndarray  # shape (n_train,) float64 balanced on fit_mask
    sample_weights_train: np.ndarray  # shape (n_train,) float64 balanced on train


@dataclass(frozen=True)
class FoldPredictions:
    """Per-test-genome predictions returned by a trainable fitter.

    Either continuous probabilities `scores` (float in [0, 1] or None when
    abstained/failed) and/or explicit `calls` in ('R', 'S', 'I', 'no_call', 'failure').
    """

    calls: tuple[str, ...]
    scores: tuple[float | None, ...]
    meta: dict[str, Any]


def _default_library_fitter(ctx: FoldFitContext) -> FoldPredictions:
    """Default trainable fitter delegating to pinned libraries (lazy-imported)."""
    outside = _get_outside()
    n_train = ctx.X_train.shape[0]
    n_test = ctx.X_test.shape[0]
    if ctx.X_train.shape[1] == 0:
        prev = float(np.mean(ctx.y_train))
        scores = tuple(prev for _ in range(n_test))
        calls = tuple("R" if prev >= CUTOFF else "S" for _ in range(n_test))
        return FoldPredictions(calls=calls, scores=scores, meta={"constant": 1})

    # Reconstruct local combined arrays ONLY over train + test, with y_test masked to False
    # (never read by fitters; only indices inside fit/hold/train are passed to library fits).
    X_comb = np.vstack([ctx.X_train, ctx.X_test]).astype(np.float32)
    y_comb = np.concatenate([ctx.y_train, np.zeros(n_test, dtype=bool)])
    w_comb = np.concatenate([ctx.sample_weights_fit, np.ones(n_test, dtype=np.float64)])
    fit_full = np.concatenate([ctx.fit_mask, np.zeros(n_test, dtype=bool)])
    hold_full = np.concatenate([ctx.hold_mask, np.zeros(n_test, dtype=bool)])
    test_full = np.concatenate([np.zeros(n_train, dtype=bool), np.ones(n_test, dtype=bool)])
    fo_local = {"fit": fit_full, "hold": hold_full, "test": test_full, "seed": ctx.seed}

    if ctx.tool == "xgboost":
        sc, info = outside.fit_xgboost(X_comb, y_comb, fo_local, w_comb)
        scores = tuple(float(v) for v in sc)
        calls = tuple("R" if v >= CUTOFF else "S" for v in scores)
        return FoldPredictions(calls=calls, scores=scores, meta=info)
    if ctx.tool == "lightgbm":
        sc, info = outside.fit_lightgbm(X_comb, y_comb, fo_local, w_comb)
        scores = tuple(float(v) for v in sc)
        calls = tuple("R" if v >= CUTOFF else "S" for v in scores)
        return FoldPredictions(calls=calls, scores=scores, meta=info)
    if ctx.tool == "elasticnet":
        # Train-only scaled grid search over inner folds using outside.enet_path
        Xtr_u8 = ctx.X_train.astype(np.uint8)
        Xte_u8 = ctx.X_test.astype(np.uint8)
        tr_idx = np.arange(n_train)
        scale_r = outside.scaled(Xtr_u8, tr_idx)
        points = outside.grid(Xtr_u8, ctx.y_train, ctx.sample_weights_train, scale_r)
        n_points = sum(len(l) for _, l in points)
        if n_points == 0:
            prev = float(np.mean(ctx.y_train))
            return FoldPredictions(
                calls=tuple("R" if prev >= CUTOFF else "S" for _ in range(n_test)),
                scores=tuple(prev for _ in range(n_test)),
                meta={"constant": 1, "null_fit": 1},
            )
        pooled = {(m, i): np.full(n_train, np.nan) for m, lams in points for i in range(len(lams))}
        bad: set[tuple[float, int]] = set()
        for j in range(INNER_K):
            part = ctx.inner_folds != j
            held = np.flatnonzero(~part)
            pj = np.flatnonzero(part)
            wj = outside.class_weights(ctx.y_train, part)
            sj = outside.scaled(Xtr_u8, pj)
            for mix, lams in points:
                for i, (coef, b, ok) in enumerate(outside.enet_path(Xtr_u8[pj], ctx.y_train[pj], wj[pj], sj, mix, lams, ctx.seed)):
                    if not ok:
                        bad.add((mix, i))
                    pooled[(mix, i)][held] = outside.predict_linear(Xtr_u8[held], sj, coef, b)
        best = None
        for mix, lams in points:
            for i, lam in enumerate(lams):
                if (mix, i) in bad:
                    continue
                loss = round(outside.lasso_log_loss(pooled[(mix, i)], ctx.y_train), 12)
                key = (loss, -round(float(lam), 12), -mix)
                if best is None or key < best[0]:
                    best = (key, mix, i)
        if best is None:
            prev = float(np.mean(ctx.y_train))
            return FoldPredictions(
                calls=tuple("R" if prev >= CUTOFF else "S" for _ in range(n_test)),
                scores=tuple(prev for _ in range(n_test)),
                meta={"constant": 1, "null_fit": 1, "refit_converged": 0},
            )
        _, mix, i = best
        lams = dict(points)[mix][: i + 1]
        coef, b, ok = outside.enet_path(Xtr_u8, ctx.y_train, ctx.sample_weights_train, scale_r, mix, lams, ctx.seed)[-1]
        nonzero = int(np.sum(coef != 0))
        if nonzero == 0:
            prev = float(np.mean(ctx.y_train))
            return FoldPredictions(
                calls=tuple("R" if prev >= CUTOFF else "S" for _ in range(n_test)),
                scores=tuple(prev for _ in range(n_test)),
                meta={"constant": 1, "null_fit": 1, "nonzero": 0, "refit_converged": int(ok)},
            )
        sc = outside.predict_linear(Xte_u8, scale_r, coef, b)
        scores = tuple(float(v) for v in sc)
        calls = tuple("R" if v >= CUTOFF else "S" for v in scores)
        return FoldPredictions(
            calls=calls,
            scores=scores,
            meta={"constant": 0, "null_fit": 0, "nonzero": nonzero, "picked_mix": mix, "picked_penalty": i, "refit_converged": int(ok)},
        )
    if ctx.tool == "kover":
        raise BenchError(
            "kover real CLI requires FASTA contigs via run_lambda executor; supply a custom fitter for matrix-mode tests"
        )
    raise BenchError(f"unsupported trainable tool: {ctx.tool!r}")


# ------------------------------------------------------------------ fold accounting


def fold_call_counts(
    truth: np.ndarray,
    calls: tuple[str, ...],
    scores: tuple[float | None, ...] | None = None,
) -> dict[str, Any]:
    """Compute per-fold counts with strict separation of no_call, I, and failure.

    NEVER maps no_call, I, or failure to S.
    """
    n = int(truth.shape[0])
    require(len(calls) == n, f"calls length {len(calls)} != test fold size {n}")
    if scores is not None:
        require(len(scores) == n, f"scores length {len(scores)} != test fold size {n}")

    tp = fp = tn = fn = 0
    no_call_r = no_call_s = 0
    i_r = i_s = 0
    failure_r = failure_s = 0
    loss_terms: list[float] = []
    has_all_called_scores = scores is not None

    for idx in range(n):
        y = bool(truth[idx])
        c = calls[idx]
        if c == "R":
            if y:
                tp += 1
            else:
                fp += 1
        elif c == "S":
            if y:
                fn += 1
            else:
                tn += 1
        elif c == "no_call":
            if y:
                no_call_r += 1
            else:
                no_call_s += 1
        elif c == "I":
            if y:
                i_r += 1
            else:
                i_s += 1
        elif c == "failure":
            if y:
                failure_r += 1
            else:
                failure_s += 1
        else:
            raise BenchError(f"invalid fold prediction call {c!r}")

        if c in ("R", "S") and scores is not None:
            sc = scores[idx]
            if sc is None:
                has_all_called_scores = False
            else:
                require(math.isfinite(float(sc)) and 0.0 <= float(sc) <= 1.0, f"score out of [0, 1]: {sc!r}")
                p = min(max(float(sc), LOSS_CLIP), 1.0 - LOSS_CLIP)
                loss_terms.append(-math.log(p) if y else -math.log(1.0 - p))

    out: dict[str, Any] = {
        "n": n,
        "tp": tp,
        "fp": fp,
        "tn": tn,
        "fn": fn,
        "no_call_r": no_call_r,
        "no_call_s": no_call_s,
        "i_r": i_r,
        "i_s": i_s,
        "failure_r": failure_r,
        "failure_s": failure_s,
    }
    if has_all_called_scores and (tp + fp + tn + fn) > 0:
        out["loss_sum"] = round(float(math.fsum(loss_terms)), 12)
    return out


# ------------------------------------------------------------------ local-only OOF scratch writer


def write_oof_scratch(
    scratch_dir: Path | str,
    *,
    result_dir: Path | str | None,
    tool: str,
    dataset: str,
    split_scheme: str,
    seed: int,
    input_fingerprint: str,
    fold_of: list[int],
    calls: list[str | None],
    scores: list[float | None],
) -> Path:
    """Write local-only out-of-fold score vector strictly in scratch (never committed).

    Refuses if `scratch_dir` is inside `reference/` or overlaps `result_dir`.
    """
    sdir = Path(scratch_dir).resolve()
    ref_dir = (REPO_ROOT / "reference").resolve()
    require(
        sdir != ref_dir and ref_dir not in sdir.parents,
        "refused: out-of-fold score vectors may never be written inside reference/",
    )
    if result_dir is not None:
        rdir = Path(result_dir).resolve()
        require(
            sdir != rdir and rdir not in sdir.parents and sdir not in rdir.parents,
            "refused: OOF scratch directory must be strictly separate from the result output directory",
        )
    sdir.mkdir(parents=True, exist_ok=True)
    out_path = sdir / f"{tool}__{dataset}__{split_scheme}__s{seed}.oof.json"
    doc = {
        "format": "grm-bench-oof-scratch-v1",
        "local_only": True,
        "never_commit": True,
        "tool": tool,
        "dataset": dataset,
        "split_scheme": split_scheme,
        "seed": seed,
        "input_fingerprint": input_fingerprint,
        "n_samples": len(calls),
        "fold_of": fold_of,
        "calls": calls,
        "scores": [None if s is None else round(float(s), 12) for s in scores],
    }
    tmp = out_path.with_name(out_path.name + ".tmp")
    tmp.write_bytes(canonical_json_bytes(doc))
    os.replace(tmp, out_path)
    return out_path


# ------------------------------------------------------------------ main trainable runner


def run_trainable_dataset(
    adapter: AdapterSpec,
    export_dir: Path | str,
    scheme: str,
    seed: int,
    *,
    fold_file: Path | str | None = None,
    selection_role: str = SELECTION_ROLE_INTERNAL_CV,
    fitter: Callable[[FoldFitContext], FoldPredictions] | None = None,
    reference_result_doc: dict[str, Any] | None = None,
    reference_report_doc: dict[str, Any] | None = None,
    release_report_doc: dict[str, Any] | None = None,
    oof_scratch_dir: Path | str | None = None,
    result_dir: Path | str | None = None,
) -> dict[str, Any]:
    """Run one trainable adapter across all outer folds of a dataset export.

    Returns a zero-ID per-dataset summary dictionary containing per-fold counts
    and fold-parity metadata, ready for aggregation by `score.py`.
    """
    require(adapter.mode == MODE_TRAINABLE, f"{adapter.tool}: expected trainable adapter, got {adapter.mode!r}")
    assert_selection_allowed(selection_role, f"training/tuning trainable adapter {adapter.tool}")

    outside = _get_outside()
    try:
        ex = outside.Export(str(export_dir))
    except (outside.Refuse, OSError, KeyError, ValueError) as exc:
        raise BenchError(f"cannot load export {export_dir}: {exc}") from exc

    dataset = ex.meta["id"]
    sp = ex.meta["species_slug"]
    ab = ex.meta["antibiotic"]
    require(ex.meta["evaluation_revision"] == EVALUATION_REVISION, f"{dataset}: revision != {EVALUATION_REVISION}")

    parity_info: dict[str, Any] | None = None
    if scheme in DRAWN_SCHEMES:
        require(fold_file is None, f"fold_file must not be supplied for drawn split {scheme!r}")
        if reference_result_doc is not None or reference_report_doc is not None or release_report_doc is not None:
            parity_info = verify_fold_parity(
                export_dir,
                seed,
                scheme=scheme,
                reference_result_doc=reference_result_doc,
                reference_report_doc=reference_report_doc,
                release_report_doc=release_report_doc,
            )
        plan, folds = outside.export_plan(ex, scheme, seed)
        ff_sha = None
        floor = ex.meta["min_support"]
    else:
        require(fold_file is not None, f"split {scheme!r} requires a grm-fold-file-v1 path")
        plan, folds, ff_sha = build_fold_file_plan(ex, fold_file, seed)
        require(plan["scheme"] == scheme, f"fold file scheme {plan['scheme']!r} != requested {scheme!r}")
        floor = plan["floor"]

    support = adapter.support_for(sp, ab)
    use_fitter = fitter or _default_library_fitter

    if folds is None:
        return {
            "tool": adapter.tool,
            "mode": adapter.mode,
            "comparison_kind": adapter.comparison_kind,
            "dataset": dataset,
            "species_slug": sp,
            "canonical_drug": ab,
            "split_scheme": scheme,
            "seed": seed,
            "selection_role": selection_role,
            "support_tier": support.support_tier,
            "unsupported_reason": support.unsupported_reason,
            "evaluation_revision": EVALUATION_REVISION,
            "input_fingerprint": ex.meta["input_fingerprint"],
            "matrix_sha256": ex.meta["matrix_sha256"],
            "plan_sha256": plan.get("sha256"),
            "fold_file_sha256": ff_sha,
            "fold_parity": parity_info,
            "unavailable_reason": plan.get("unavailable_reason"),
            "folds": [],
            "oof_scratch_path": None,
        }

    outside.check_folds(folds, ex.y, ex.n)
    y = ex.y
    oof_scores: list[float | None] = [None] * ex.n
    oof_calls: list[str | None] = [None] * ex.n
    oof_fold_of: list[int] = [-1] * ex.n
    fold_rows: list[dict[str, Any]] = []

    for fo in folds:
        train = fo["train"]
        test = fo["test"]
        tr_idx = np.flatnonzero(train)
        te_idx = np.flatnonzero(test)

        # Screening strictly on outer training genomes ONLY (never reading test genomes!)
        rows = fo["screened_rows"]
        if rows is None:
            n_cand, tops = outside.screen_fold(ex.x, train, y, floor, {"train": train})
            rows = tops["train"]
            n_candidates = int(n_cand)
        else:
            n_candidates = int(fo["candidates"])

        if rows.shape[0] > 0:
            X_all = outside.dense(np.asarray(ex.x[np.sort(rows)]), ex.n)
            X_train = np.ascontiguousarray(X_all[tr_idx])
            X_test = np.ascontiguousarray(X_all[te_idx])
        else:
            X_train = np.zeros((tr_idx.shape[0], 0), dtype=np.uint8)
            X_test = np.zeros((te_idx.shape[0], 0), dtype=np.uint8)

        # Map fit, hold, inner masks strictly into local 0..n_train-1 coordinates
        fit_local = fo["fit"][tr_idx]
        hold_local = fo["hold"][tr_idx]
        inner_local = fo["inner"][tr_idx]
        y_train = y[tr_idx].copy()
        w_fit = outside.class_weights(y, fo["fit"])[tr_idx]
        w_train = outside.class_weights(y, train)[tr_idx]

        ctx = FoldFitContext(
            tool=adapter.tool,
            fold=int(fo["fold"]),
            seed=int(fo["seed"]),
            X_train=X_train,
            y_train=y_train,
            X_test=X_test,
            fit_mask=fit_local,
            hold_mask=hold_local,
            inner_folds=inner_local,
            sample_weights_fit=w_fit,
            sample_weights_train=w_train,
        )

        try:
            preds = use_fitter(ctx)
            require(isinstance(preds, FoldPredictions), "fitter must return FoldPredictions")
            require(len(preds.calls) == te_idx.shape[0], "fitter returned wrong number of test calls")
            require(len(preds.scores) == te_idx.shape[0], "fitter returned wrong number of test scores")
        except Exception as exc:  # noqa: BLE001
            # Runtime failure inside one outer fold: explicit failure counts, NEVER mapped to S!
            preds = FoldPredictions(
                calls=tuple("failure" for _ in range(te_idx.shape[0])),
                scores=tuple(None for _ in range(te_idx.shape[0])),
                meta={"failure_reason": f"fitter_exception:{type(exc).__name__}"},
            )

        counts = fold_call_counts(y[te_idx], preds.calls, preds.scores)
        fold_row: dict[str, Any] = {
            "fold": int(fo["fold"]),
            "n_train": int(tr_idx.shape[0]),
            "n_train_resistant": int(np.sum(y_train)),
            "n_candidates": n_candidates,
            "n_screened": int(rows.shape[0]),
            **counts,
        }
        fold_rows.append(fold_row)

        for local_i, global_i in enumerate(te_idx):
            oof_fold_of[int(global_i)] = int(fo["fold"])
            oof_calls[int(global_i)] = preds.calls[local_i]
            oof_scores[int(global_i)] = preds.scores[local_i]

    oof_path_str: str | None = None
    if oof_scratch_dir is not None:
        oof_path = write_oof_scratch(
            oof_scratch_dir,
            result_dir=result_dir,
            tool=adapter.tool,
            dataset=dataset,
            split_scheme=scheme,
            seed=seed,
            input_fingerprint=ex.meta["input_fingerprint"],
            fold_of=oof_fold_of,
            calls=oof_calls,
            scores=oof_scores,
        )
        oof_path_str = str(oof_path)

    return {
        "tool": adapter.tool,
        "mode": adapter.mode,
        "comparison_kind": adapter.comparison_kind,
        "dataset": dataset,
        "species_slug": sp,
        "canonical_drug": ab,
        "split_scheme": scheme,
        "seed": seed,
        "selection_role": selection_role,
        "support_tier": support.support_tier,
        "unsupported_reason": support.unsupported_reason,
        "evaluation_revision": EVALUATION_REVISION,
        "input_fingerprint": ex.meta["input_fingerprint"],
        "matrix_sha256": ex.meta["matrix_sha256"],
        "plan_sha256": plan.get("sha256"),
        "fold_file_sha256": ff_sha,
        "fold_parity": parity_info,
        "unavailable_reason": None,
        "folds": fold_rows,
        "oof_scratch_path": oof_path_str,
    }
