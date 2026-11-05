"""Masked multi-task sparse-group logistic regression, a development-only oracle.

One logistic task per drug shares the genome features. The fit minimises

    sum_t sum_{i observed in t} w_it * logloss(y_it, x_i' beta_t + b_t)
        + lambda1 * ||B||_1 + lambda_group * sum_j ||B_j.||_2
        + (lambda2 / 2) * ||B||_F^2

where B is features x tasks, row B_j. groups feature j across tasks, and the
intercepts b are unpenalised. lambda2 > 0 is required: it makes the objective
strongly convex in the coefficients B. The intercepts get no ridge, so the
objective has a finite minimiser only when every task observes both classes.
fit, objective and kkt_violation refuse a task with one observed class under
every normalisation; callers mask such tasks out first. Then the minimiser is
unique, so solvers can be compared on it. lambda_group = 0 is the per-task
elastic net through the same solver.

The loss weight w_it is 0 for an unobserved label, and otherwise set by the
normalisation:

  per_task (default)   1 / n_t, with n_t observed labels in task t. Each task
                       contributes its mean log loss, so every task carries the
                       same weight however many labels it has.
  unweighted           T / N, with N observed labels across all T tasks. Every
                       label weighs the same, so tasks with more labels count
                       more. The total weight is T, as for per_task, and the two
                       are identical when every task has the same n_t.
  inverse_prevalence   1 / (2 n_t+) for a positive label and 1 / (2 n_t-) for a
                       negative one, so within each task the two classes share
                       the task's unit weight equally. Every task needs both.

The solver is FISTA with the constant step 1/L, where L bounds the Lipschitz
constant of the smooth part, and a monotone safeguard: a step that would raise
the objective is discarded and momentum restarts from the current point, so the
recorded objective never increases beyond rounding (ROUNDING_SLACK). A fit
converges only when the relative objective change and the proximal-gradient
norm are both within tolerance; otherwise its status is unconverged and
select_candidate rejects it.

Known limits:
  - converged means a small proximal-gradient norm and objective change, not
    accurate coefficients. Coefficient error grows roughly like
    kkt_tol / lambda2: at the default kkt_tol = 1e-7 on separable data, an
    interior-point solver measured coefficient gaps of 3e-6, 3e-4 and 2e-2 at
    lambda2 = 1e-2, 1e-4 and 1e-6, with objective gaps of 1e-9 or less. Pass
    much tighter tolerances to use the coefficients as an oracle.
  - The step 1/L includes lambda2, but the intercepts have no ridge curvature,
    so a large lambda2 slows them: at lambda2 = 1e6 an intercept-only solution
    took about 15,600 of the default 20,000 iterations.
  - X must be dense, and L costs one SVD per task.
  - X is used as given. Standardise its columns first if the penalty should
    treat features equally.
"""
import math
from dataclasses import dataclass

import numpy as np

CONVERGED = "converged"
UNCONVERGED = "unconverged"
NORMALISATIONS = ("per_task", "unweighted", "inverse_prevalence")
# Relative objective rise a plain proximal-gradient step may show from rounding alone.
ROUNDING_SLACK = 64 * np.finfo(float).eps


@dataclass(frozen=True, eq=False)
class MultitaskFit:
    coef: np.ndarray  # features x tasks
    intercept: np.ndarray  # tasks
    status: str
    diagnostics: dict
    lambda1: float
    lambda_group: float
    lambda2: float
    normalisation: str

    @property
    def converged(self):
        return self.status == CONVERGED

    def predict_proba(self, X):
        """Positive-label probability for every genome and task, genomes x tasks."""
        X = _features(X)
        if X.shape[1] != self.coef.shape[0]:
            raise ValueError(f"X has {X.shape[1]} features; the fit has {self.coef.shape[0]}")
        return _sigmoid(X @ self.coef + self.intercept)


@dataclass(frozen=True, eq=False)
class _Problem:
    X: np.ndarray
    Y: np.ndarray  # observed labels, 0 where unobserved
    W: np.ndarray  # loss weights, 0 where unobserved
    lambda1: float
    lambda_group: float
    lambda2: float
    normalisation: str


def _features(X):
    X = np.asarray(X, dtype=float)
    if X.ndim != 2 or not np.isfinite(X).all():
        raise ValueError("X must be a finite genomes x features matrix")
    return X


def _labels(Y, mask):
    """Observed labels with unobserved entries zeroed, and the boolean mask.

    mask=None treats exactly the non-NaN entries of Y as observed.
    """
    Y = np.asarray(Y, dtype=float)
    if Y.ndim != 2:
        raise ValueError("Y must be a genomes x tasks matrix")
    mask = ~np.isnan(Y) if mask is None else np.asarray(mask)
    if mask.dtype != bool or mask.shape != Y.shape:
        raise ValueError("mask must be a boolean array shaped like Y")
    observed = Y[mask]
    if not np.isin(observed, (0.0, 1.0)).all():
        raise ValueError("observed labels must be 0 or 1")
    return np.where(mask, Y, 0.0), mask


def loss_weights(Y, mask=None, normalisation="per_task"):
    """Loss weight of every label, genomes x tasks; see the module docstring.

    A task with one observed class is allowed here, as when scoring a
    validation fold, except under inverse_prevalence, whose weights need both.
    """
    return _weights(*_labels(Y, mask), normalisation)


def _weights(Y, mask, normalisation):
    if normalisation not in NORMALISATIONS:
        raise ValueError(f"normalisation must be one of {NORMALISATIONS}, not {normalisation!r}")
    counts = mask.sum(axis=0)
    if (counts == 0).any():
        raise ValueError(f"tasks {np.flatnonzero(counts == 0).tolist()} have no observed labels")
    if normalisation == "per_task":
        return mask / counts
    if normalisation == "unweighted":
        return mask * (Y.shape[1] / counts.sum())
    positive = _require_both_classes(Y, mask)
    return np.where(Y == 1, 0.5 / positive, 0.5 / (counts - positive)) * mask


def _require_both_classes(Y, mask):
    """Positive counts per task; raises naming every task with one observed class."""
    positive = (mask & (Y == 1)).sum(axis=0)
    one_class = np.flatnonzero((positive == 0) | (positive == mask.sum(axis=0))).tolist()
    if one_class:
        raise ValueError(f"tasks {one_class} observe only one class, so their unpenalised "
                         "intercepts have no finite minimiser; mask those tasks out first")
    return positive


def _checked(name, value, positive=False):
    value = float(value)
    if not (value > 0 if positive else value >= 0) or not math.isfinite(value):
        raise ValueError(f"{name} must be finite and {'> 0' if positive else '>= 0'}, not {value}")
    return value


def _problem(X, Y, mask, lambda1, lambda_group, lambda2, normalisation):
    X = _features(X)
    Y, mask = _labels(Y, mask)
    if Y.shape[0] != X.shape[0]:
        raise ValueError(f"X has {X.shape[0]} genomes; Y has {Y.shape[0]}")
    W = _weights(Y, mask, normalisation)
    _require_both_classes(Y, mask)
    return _Problem(
        X=X, Y=Y, W=W,
        lambda1=_checked("lambda1", lambda1),
        lambda_group=_checked("lambda_group", lambda_group),
        lambda2=_checked("lambda2", lambda2, positive=True),
        normalisation=normalisation,
    )


def _sigmoid(Z):
    return np.exp(-np.logaddexp(0.0, -Z))


def _smooth(problem, coef, intercept):
    """Weighted log loss plus the ridge term, and its gradients in coef and intercept."""
    Z = problem.X @ coef + intercept
    value = (np.sum(problem.W * (np.logaddexp(0.0, Z) - problem.Y * Z))
             + 0.5 * problem.lambda2 * np.sum(coef * coef))
    R = problem.W * (_sigmoid(Z) - problem.Y)
    return value, problem.X.T @ R + problem.lambda2 * coef, R.sum(axis=0)


def _nonsmooth(problem, coef):
    return (problem.lambda1 * np.abs(coef).sum()
            + problem.lambda_group * np.linalg.norm(coef, axis=1).sum())


def objective(fit, X, Y, mask=None):
    """The fit's objective value on (X, Y, mask) under the fit's penalties."""
    problem = _problem(X, Y, mask, fit.lambda1, fit.lambda_group, fit.lambda2, fit.normalisation)
    return _smooth(problem, fit.coef, fit.intercept)[0] + _nonsmooth(problem, fit.coef)


def prox_sparse_group(V, step, lambda1, lambda_group):
    """Proximal operator of step * (lambda1 ||.||_1 + lambda_group sum_j ||V_j.||_2).

    Soft-threshold every entry by step * lambda1, then shrink each row (one
    feature's coefficients across tasks) toward zero by step * lambda_group in
    Euclidean norm. A row whose thresholded norm is at most that becomes zero.
    """
    V = np.asarray(V, dtype=float)
    S = np.sign(V) * np.maximum(np.abs(V) - step * lambda1, 0.0)
    norms = np.linalg.norm(S, axis=1, keepdims=True)
    threshold = step * lambda_group
    keep = norms > threshold
    return np.where(keep, 1.0 - threshold / np.where(keep, norms, 1.0), 0.0) * S


def _lipschitz(problem):
    """Bound on the Lipschitz constant of the smooth part's gradient.

    Each task's log-loss Hessian is at most [X 1]' diag(w_t) [X 1] / 4, tasks
    form independent blocks, and the ridge adds lambda2.
    """
    X1 = np.hstack([problem.X, np.ones((problem.X.shape[0], 1))])
    largest = max(np.linalg.norm(np.sqrt(w)[:, None] * X1, 2) ** 2 for w in problem.W.T)
    return 0.25 * largest + problem.lambda2


def _prox_grad_norm(problem, step, coef, grad_coef, grad_intercept):
    """Largest entry of the proximal-gradient map (x - prox(x - step * grad)) / step.

    It is zero exactly at the minimiser.
    """
    moved = prox_sparse_group(coef - step * grad_coef, step, problem.lambda1, problem.lambda_group)
    return max(np.abs(coef - moved).max(initial=0.0) / step, np.abs(grad_intercept).max())


def _kkt_violation(problem, coef, intercept):
    """Largest violation of the subgradient optimality conditions, checked directly.

    Intercept gradients must be 0. In a nonzero row, a nonzero entry needs
    grad + lambda1 sign(beta) + lambda_group beta / ||row|| = 0 and a zero entry
    needs |grad| <= lambda1. A zero row needs ||soft(grad, lambda1)|| <= lambda_group.
    """
    _, grad, grad_intercept = _smooth(problem, coef, intercept)
    norms = np.linalg.norm(coef, axis=1, keepdims=True)
    zero_row = norms[:, 0] == 0
    active = coef != 0
    stationarity = np.abs(grad + problem.lambda1 * np.sign(coef)
                          + problem.lambda_group * coef / np.where(zero_row[:, None], 1.0, norms))
    excess = np.maximum(np.abs(grad) - problem.lambda1, 0.0)
    return max(
        np.abs(grad_intercept).max(),
        stationarity[active].max(initial=0.0),
        excess[~active & ~zero_row[:, None]].max(initial=0.0),
        (np.linalg.norm(excess[zero_row], axis=1) - problem.lambda_group).max(initial=0.0),
    )


def kkt_violation(fit, X, Y, mask=None):
    """Largest optimality-condition violation of the fit on (X, Y, mask)."""
    problem = _problem(X, Y, mask, fit.lambda1, fit.lambda_group, fit.lambda2, fit.normalisation)
    return _kkt_violation(problem, fit.coef, fit.intercept)


def fit(X, Y, mask=None, *, lambda1, lambda_group, lambda2, normalisation="per_task",
        max_iter=20000, objective_tol=1e-10, kkt_tol=1e-7):
    """Fit the sparse-group model from zero coefficients and intercepts.

    Converges when an accepted step changes the objective by at most
    objective_tol relative to max(1, |objective|) and the proximal-gradient norm
    is at most kkt_tol. The status is unconverged when max_iter runs out, or
    when a plain proximal-gradient step raises the objective beyond rounding
    (stalled), which a valid step 1/L cannot do.
    """
    problem = _problem(X, Y, mask, lambda1, lambda_group, lambda2, normalisation)
    if int(max_iter) != max_iter or max_iter < 1:
        raise ValueError(f"max_iter must be a positive integer, not {max_iter}")
    objective_tol = _checked("objective_tol", objective_tol, positive=True)
    kkt_tol = _checked("kkt_tol", kkt_tol, positive=True)
    step = 1.0 / _lipschitz(problem)
    coef = np.zeros((problem.X.shape[1], problem.Y.shape[1]))
    intercept = np.zeros(problem.Y.shape[1])
    smooth, grad_coef, grad_intercept = _smooth(problem, coef, intercept)
    value = smooth + _nonsmooth(problem, coef)
    prox_grad = _prox_grad_norm(problem, step, coef, grad_coef, grad_intercept)
    history = [value]
    relative_change = math.inf
    momentum, extrapolated, at_current = 1.0, (coef, intercept), True
    restarts = 0
    reason = "max_iter"
    for iteration in range(1, int(max_iter) + 1):
        if at_current:
            step_grad_coef, step_grad_intercept = grad_coef, grad_intercept
        else:
            step_grad_coef, step_grad_intercept = _smooth(problem, *extrapolated)[1:]
        new_coef = prox_sparse_group(extrapolated[0] - step * step_grad_coef, step,
                                     problem.lambda1, problem.lambda_group)
        new_intercept = extrapolated[1] - step * step_grad_intercept
        new_smooth, new_grad_coef, new_grad_intercept = _smooth(problem, new_coef, new_intercept)
        new_value = new_smooth + _nonsmooth(problem, new_coef)
        # A plain step from the current point lowers the objective in exact
        # arithmetic because step <= 1/L. Near the minimiser that decrease
        # falls below rounding, so such a step may tie within ROUNDING_SLACK;
        # an extrapolated step must not raise the objective at all.
        slack = ROUNDING_SLACK * max(1.0, abs(value)) if at_current else 0.0
        if new_value <= value + slack:
            relative_change = abs(value - new_value) / max(1.0, abs(new_value))
            previous = coef, intercept
            coef, intercept, value = new_coef, new_intercept, new_value
            grad_coef, grad_intercept = new_grad_coef, new_grad_intercept
            history.append(value)
            prox_grad = _prox_grad_norm(problem, step, coef, grad_coef, grad_intercept)
            if relative_change <= objective_tol and prox_grad <= kkt_tol:
                reason = "tolerances met"
                break
            next_momentum = (1.0 + math.sqrt(1.0 + 4.0 * momentum * momentum)) / 2.0
            weight = (momentum - 1.0) / next_momentum
            extrapolated = (coef + weight * (coef - previous[0]),
                            intercept + weight * (intercept - previous[1]))
            momentum, at_current = next_momentum, weight == 0.0
            continue
        history.append(value)
        if at_current:
            # Even a plain step raised the objective beyond rounding.
            reason = "stalled"
            break
        restarts += 1
        momentum, extrapolated, at_current = 1.0, (coef, intercept), True
    status = CONVERGED if reason == "tolerances met" else UNCONVERGED
    diagnostics = {
        "reason": reason,
        "iterations": iteration,
        "objective": value,
        "relative_change": relative_change,
        "prox_grad_norm": prox_grad,
        "kkt_violation": _kkt_violation(problem, coef, intercept),
        "step": step,
        "restarts": restarts,
        "objective_history": tuple(history),
    }
    return MultitaskFit(coef=coef, intercept=intercept, status=status, diagnostics=diagnostics,
                        lambda1=problem.lambda1, lambda_group=problem.lambda_group,
                        lambda2=problem.lambda2, normalisation=normalisation)


def fit_elastic_net(X, Y, mask=None, *, lambda1, lambda2, **options):
    """The per-task elastic-net control: fit() with lambda_group = 0."""
    return fit(X, Y, mask, lambda1=lambda1, lambda_group=0.0, lambda2=lambda2, **options)


def select_candidate(candidates, scores):
    """Index of the lowest-scoring candidate whose fits all converged.

    A candidate is one MultitaskFit or a non-empty sequence of them, such as one
    per cross-validation fold, and is eligible only if every fit converged.
    Lower scores are better, such as a validation log loss; ties go to the
    first index. Raises ValueError when no candidate is eligible or an eligible
    candidate's score is not finite (NaN or +-inf).
    """
    candidates, scores = list(candidates), [float(score) for score in scores]
    if len(candidates) != len(scores):
        raise ValueError(f"{len(candidates)} candidates but {len(scores)} scores")
    best = None
    for index, (candidate, score) in enumerate(zip(candidates, scores)):
        fits = (candidate,) if isinstance(candidate, MultitaskFit) else tuple(candidate)
        if not fits:
            raise ValueError(f"candidate {index} has no fits")
        if not all(member.converged for member in fits):
            continue
        if not math.isfinite(score):
            raise ValueError(f"candidate {index} converged but its score is {score}")
        if best is None or score < scores[best]:
            best = index
    if best is None:
        raise ValueError("no candidate converged")
    return best
