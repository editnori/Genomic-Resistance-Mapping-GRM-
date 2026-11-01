"""scripts/multitask/model.py on tiny synthetic problems, against solvers written
here independently of its proximal-gradient code.

usage: python3 -I -B -m unittest discover -s tests/multitask -t .
"""
import itertools
import math
import unittest
from unittest import mock

import numpy as np

from scripts.multitask import model


def sigmoid(z):
    return 1.0 / (1.0 + np.exp(-z))


def intercept_gradient(fit, X, Y, mask=None):
    """Largest |d loss / d b_t|, written out from the loss weights."""
    residual = model.loss_weights(Y, mask) * (sigmoid(X @ fit.coef + fit.intercept) - np.nan_to_num(Y))
    return np.abs(residual.sum(axis=0)).max()


def tiny_problem():
    """40 genomes, 3 features and 2 tasks: feature 0 drives both tasks, feature
    1 neither, feature 2 only task 0. Task 1 misses four labels."""
    rng = np.random.default_rng(20)
    X = rng.normal(size=(40, 3))
    truth = np.array([[2.0, 1.5], [0.0, 0.0], [1.2, 0.0]])
    Y = (rng.random((40, 2)) < sigmoid(X @ truth - 0.3)).astype(float)
    mask = np.ones_like(Y, dtype=bool)
    mask[[5, 11, 12, 30], 1] = False
    Y[~mask] = np.nan
    return X, Y, mask


def enumeration_oracle(X, Y, mask, lambda1, lambda_group, lambda2):
    """Minimiser of the per_task objective found by trying every sign pattern.

    For a pattern s in {-1, 0, +1} per coefficient, F_s replaces |beta| with
    s * beta and fixes the s = 0 coefficients at zero. F_s is smooth and convex
    there, so damped Newton finds its minimiser. The true minimiser's pattern
    reproduces it, and any pattern whose minimiser keeps its own signs has
    F_s equal to the true objective there, so the lowest such value wins.
    A pattern is dropped once a coefficient leaves its sign, which only ever
    discards wrong patterns: dropping the true one would make the comparison
    with the solver fail, never pass. Parameters are a (1 + features) x tasks
    matrix whose row 0 is the intercepts.
    """
    tasks = Y.shape[1]
    X1 = np.hstack([np.ones((len(X), 1)), X])
    weight = mask / mask.sum(axis=0)
    labels = np.where(mask, Y, 0.0)

    def value(signs, theta):
        Z, B = X1 @ theta, theta[1:]
        return (np.sum(weight * (np.logaddexp(0.0, Z) - labels * Z)) + lambda1 * np.sum(signs * B)
                + lambda_group * np.linalg.norm(B, axis=1).sum() + 0.5 * lambda2 * np.sum(B * B))

    def gradient_hessian(signs, theta):
        B = theta[1:]
        P = sigmoid(X1 @ theta)
        norms = np.linalg.norm(B, axis=1, keepdims=True)
        grad = X1.T @ (weight * (P - labels))
        grad[1:] += lambda1 * signs + lambda2 * B + lambda_group * B / np.where(norms > 0, norms, 1.0)
        size = theta.size
        hess = np.zeros((size, size))
        for t in range(tasks):
            hess[t::tasks, t::tasks] = X1.T @ (X1 * (weight[:, t] * P[:, t] * (1 - P[:, t]))[:, None])
        hess[tasks:, tasks:] += lambda2 * np.eye(size - tasks)
        for j, (row, norm) in enumerate(zip(B, norms[:, 0]), start=1):
            if norm > 0:
                block = slice(j * tasks, (j + 1) * tasks)
                hess[block, block] += lambda_group * (np.eye(tasks) - np.outer(row, row) / norm ** 2) / norm
        return grad.ravel(), hess

    best = None
    for pattern in itertools.product((-1.0, 0.0, 1.0), repeat=X.shape[1] * tasks):
        signs = np.reshape(pattern, (X.shape[1], tasks))
        free = np.concatenate([np.ones(tasks, dtype=bool), signs.ravel() != 0])
        theta = np.vstack([np.zeros((1, tasks)), 0.1 * signs])
        support = signs != 0
        for _ in range(50):
            if np.any(np.sign(theta[1:][support]) != signs[support]):
                break
            grad, hess = gradient_hessian(signs, theta)
            grad, hess = grad[free], hess[np.ix_(free, free)]
            direction = np.zeros(theta.size)
            direction[free] = np.linalg.solve(hess, -grad)
            if -grad @ direction[free] < 1e-26:
                break  # Newton decrement at rounding level
            current, size = value(signs, theta), 1.0
            while size > 1e-10:
                trial = theta + size * direction.reshape(theta.shape)
                if value(signs, trial) <= current + 1e-4 * size * grad @ direction[free]:
                    break
                size /= 2
            else:
                break  # no further decrease to rounding
            theta = trial
        coef = theta[1:]
        if np.any(np.sign(coef[support]) != signs[support]):
            continue
        if np.abs(gradient_hessian(signs, theta)[0][free]).max() > 1e-9:
            continue
        objective = value(np.sign(coef), theta)
        if best is None or objective < best[0]:
            best = objective, coef, theta[0]
    return best


def coordinate_descent_elastic_net(x, y, lambda1, lambda2, sweeps=2000):
    """One task's (1/n) logloss + lambda1 |beta|_1 + (lambda2/2) |beta|^2.

    Cyclic exact coordinate minimisation, each coordinate found by bisection on
    its monotone derivative; no proximal operator or step size is involved.
    """
    beta, intercept = np.zeros(x.shape[1]), 0.0

    def slope(j, value):
        if j is None:
            return np.mean(sigmoid(x @ beta + value) - y)
        trial = beta.copy()
        trial[j] = value
        return np.mean(x[:, j] * (sigmoid(x @ trial + intercept) - y)) + lambda2 * value

    def root(f, low, high):
        while f(low) > 0:
            low *= 2
        while f(high) < 0:
            high *= 2
        for _ in range(200):
            middle = 0.5 * (low + high)
            if middle in (low, high):
                break
            low, high = (middle, high) if f(middle) < 0 else (low, middle)
        return 0.5 * (low + high)

    for _ in range(sweeps):
        before = np.append(beta, intercept)
        intercept = root(lambda b: slope(None, b), -1.0, 1.0)
        for j in range(len(beta)):
            at_zero = slope(j, 0.0)
            if abs(at_zero) <= lambda1:
                beta[j] = 0.0
            elif at_zero < -lambda1:
                beta[j] = root(lambda v: slope(j, v) + lambda1, 1e-300, 1.0)
            else:
                beta[j] = root(lambda v: slope(j, v) - lambda1, -1.0, -1e-300)
        if np.abs(np.append(beta, intercept) - before).max() < 1e-14:
            break
    return beta, intercept


TIGHT = {"objective_tol": 1e-15, "kkt_tol": 1e-11, "max_iter": 200000}


class ProxTests(unittest.TestCase):
    def test_soft_threshold_alone_can_zero_a_group(self):
        V = np.array([[0.3, -0.5], [2.0, 0.0]])
        out = model.prox_sparse_group(V, step=1.0, lambda1=0.5, lambda_group=0.0)
        np.testing.assert_array_equal(out, [[0.0, 0.0], [1.5, 0.0]])

    def test_one_task_reduces_to_scalar_soft_threshold_at_lambda1_plus_lambda_group(self):
        V = np.array([[3.0], [-2.0], [0.7], [-0.75]])
        out = model.prox_sparse_group(V, step=0.5, lambda1=1.0, lambda_group=0.5)
        np.testing.assert_allclose(out, [[2.25], [-1.25], [0.0], [0.0]], rtol=0, atol=1e-15)

    def test_group_threshold_edge(self):
        # Soft-thresholding by 2 * 0.25 leaves (3, -4), norm 5; the group
        # threshold is 2 * lambda_group.
        V = np.array([[3.5, -4.5]])
        at_edge = model.prox_sparse_group(V, step=2.0, lambda1=0.25, lambda_group=2.5)
        np.testing.assert_array_equal(at_edge, [[0.0, 0.0]])
        above = model.prox_sparse_group(V, step=2.0, lambda1=0.25, lambda_group=2.495)
        np.testing.assert_allclose(above, [[0.006, -0.008]], rtol=1e-9, atol=0)

    def test_group_shrink_keeps_the_thresholded_direction(self):
        V = np.array([[3.0, -1.0, 0.2]])
        out = model.prox_sparse_group(V, step=1.0, lambda1=0.5, lambda_group=1.0)
        thresholded = np.array([2.5, -0.5, 0.0])
        np.testing.assert_allclose(out[0], (1 - 1 / np.linalg.norm(thresholded)) * thresholded, rtol=1e-15)


class SolverTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.X, cls.Y, cls.mask = tiny_problem()
        cls.penalties = {"lambda1": 0.02, "lambda_group": 0.04, "lambda2": 0.05}
        cls.fit = model.fit(cls.X, cls.Y, cls.mask, **cls.penalties, **TIGHT)

    def test_matches_sign_pattern_enumeration_oracle(self):
        value, coef, intercept = enumeration_oracle(self.X, self.Y, self.mask, **self.penalties)
        self.assertEqual(self.fit.status, model.CONVERGED)
        # The case must reach a zero group, a zero inside a nonzero group, and
        # a dense group, or it would not exercise the sparse-group prox.
        zero_rows = np.all(coef == 0, axis=1)
        partial_rows = np.any(coef == 0, axis=1) & ~zero_rows
        self.assertTrue(zero_rows.any() and partial_rows.any() and np.all(coef != 0, axis=1).any(), coef)
        # The group penalty, not soft-thresholding alone, zeroes that row.
        control = model.fit_elastic_net(self.X, self.Y, self.mask, lambda1=0.02, lambda2=0.05, **TIGHT)
        self.assertTrue(np.any(control.coef[zero_rows] != 0), control.coef)
        np.testing.assert_array_equal(self.fit.coef == 0, coef == 0)
        np.testing.assert_allclose(self.fit.coef, coef, rtol=0, atol=1e-8)
        np.testing.assert_allclose(self.fit.intercept, intercept, rtol=0, atol=1e-8)
        self.assertAlmostEqual(self.fit.diagnostics["objective"], value, delta=1e-12)
        self.assertAlmostEqual(model.objective(self.fit, self.X, self.Y, self.mask), value, delta=1e-12)

    def test_elastic_net_control_matches_coordinate_descent(self):
        control = model.fit_elastic_net(self.X, self.Y, self.mask, lambda1=0.02, lambda2=0.05, **TIGHT)
        self.assertEqual(control.status, model.CONVERGED)
        self.assertTrue(np.any(control.coef == 0), control.coef)
        for t in range(2):
            rows = self.mask[:, t]
            beta, intercept = coordinate_descent_elastic_net(self.X[rows], self.Y[rows, t], 0.02, 0.05)
            np.testing.assert_array_equal(control.coef[:, t] == 0, beta == 0)
            np.testing.assert_allclose(control.coef[:, t], beta, rtol=0, atol=1e-8)
            self.assertAlmostEqual(control.intercept[t], intercept, delta=1e-8)

    def test_zero_group_penalty_equals_separate_task_fits(self):
        for normalisation in model.NORMALISATIONS:
            joint = model.fit(self.X, self.Y, self.mask, lambda1=0.02, lambda_group=0.0, lambda2=0.1,
                              normalisation=normalisation, **TIGHT)
            weights = model.loss_weights(self.Y, self.mask, normalisation)
            for t in range(2):
                rows = self.mask[:, t]
                # A lone task's weights total 1, so divide the lambdas by the
                # joint fit's total weight on this task to keep the same minimiser.
                scale = 1.0 / weights[rows, t].sum()
                alone = "inverse_prevalence" if normalisation == "inverse_prevalence" else "per_task"
                single = model.fit(self.X[rows], self.Y[rows, t:t + 1], lambda1=0.02 * scale,
                                   lambda_group=0.0, lambda2=0.1 * scale, normalisation=alone, **TIGHT)
                np.testing.assert_allclose(joint.coef[:, t], single.coef[:, 0], rtol=0, atol=1e-8)
                np.testing.assert_allclose(joint.intercept[t], single.intercept[0], rtol=0, atol=1e-8)

    def test_objective_never_increases(self):
        # Correlated features and a weak ridge make FISTA overshoot, so the
        # monotone safeguard has to discard steps.
        rng = np.random.default_rng(3)
        base = rng.normal(size=(60, 1))
        X = np.hstack([base + 0.05 * rng.normal(size=(60, 3)), rng.normal(size=(60, 1))])
        Y = (rng.random((60, 2)) < sigmoid(np.hstack([2 * base, -base]))).astype(float)
        result = model.fit(X, Y, lambda1=1e-3, lambda_group=1e-3, lambda2=1e-4)
        history = np.array(result.diagnostics["objective_history"])
        self.assertEqual(result.status, model.CONVERGED)
        self.assertGreater(result.diagnostics["restarts"], 0)
        slack = model.ROUNDING_SLACK * np.maximum(1.0, np.abs(history[:-1]))
        self.assertTrue(np.all(np.diff(history) <= slack), np.diff(history).max())
        self.assertLess(history[-1], history[0])

    def test_kkt_conditions_hold_at_the_solution_only(self):
        self.assertLess(model.kkt_violation(self.fit, self.X, self.Y, self.mask), 1e-9)
        self.assertEqual(self.fit.diagnostics["kkt_violation"],
                         model.kkt_violation(self.fit, self.X, self.Y, self.mask))
        moved = model.MultitaskFit(**{**vars(self.fit), "coef": self.fit.coef * 0.9})
        self.assertGreater(model.kkt_violation(moved, self.X, self.Y, self.mask), 1e-3)
        # Shifting the intercepts also moves the coefficient gradients, so
        # require the violation to reach the intercept gradient itself.
        shifted = model.MultitaskFit(**{**vars(self.fit), "intercept": self.fit.intercept + 1e-3})
        gradient = intercept_gradient(shifted, self.X, self.Y, self.mask)
        self.assertGreater(gradient, 1e-5)
        self.assertGreaterEqual(model.kkt_violation(shifted, self.X, self.Y, self.mask), gradient * (1 - 1e-9))

    def test_masked_labels_do_not_affect_the_fit(self):
        for replacement in (0.0, 1.0, 7.0):
            Y = np.where(self.mask, self.Y, replacement)
            other = model.fit(self.X, Y, self.mask, **self.penalties, **TIGHT)
            np.testing.assert_array_equal(other.coef, self.fit.coef)
            np.testing.assert_array_equal(other.intercept, self.fit.intercept)
            self.assertEqual(other.diagnostics["objective_history"], self.fit.diagnostics["objective_history"])
        # Without a mask, NaN marks the unobserved labels.
        implicit = model.fit(self.X, self.Y, **self.penalties, **TIGHT)
        np.testing.assert_array_equal(implicit.coef, self.fit.coef)

    def test_unconverged_fit_is_flagged_and_rejected(self):
        short = model.fit(self.X, self.Y, self.mask, **self.penalties, max_iter=3)
        self.assertEqual(short.status, model.UNCONVERGED)
        self.assertEqual(short.diagnostics["reason"], "max_iter")
        self.assertEqual(short.diagnostics["iterations"], 3)
        self.assertGreater(short.diagnostics["prox_grad_norm"], 1e-7)
        self.assertEqual(model.select_candidate([short, self.fit], [0.1, 0.5]), 1)
        self.assertEqual(model.select_candidate([[self.fit, short], [self.fit]], [0.1, 0.5]), 1)
        with self.assertRaisesRegex(ValueError, "no candidate converged"):
            model.select_candidate([short], [0.1])

    def test_select_candidate_prefers_the_lowest_score_then_the_first(self):
        self.assertEqual(model.select_candidate([self.fit] * 3, [0.4, 0.2, 0.2]), 1)
        for score in (math.nan, math.inf, -math.inf):
            with self.assertRaisesRegex(ValueError, "score is"):
                model.select_candidate([self.fit], [score])
        short = model.fit(self.X, self.Y, self.mask, **self.penalties, max_iter=3)
        self.assertEqual(model.select_candidate([short, self.fit], [-math.inf, 0.3]), 1)

    def test_one_class_task_is_refused_under_every_normalisation(self):
        Y = self.Y.copy()
        Y[self.mask[:, 1], 1] = 1.0
        for normalisation in model.NORMALISATIONS:
            with self.subTest(normalisation=normalisation):
                with self.assertRaisesRegex(ValueError, r"tasks \[1\] observe only one class"):
                    model.fit(self.X, Y, self.mask, **self.penalties, normalisation=normalisation)
        Y[self.mask[:, 0], 0] = 0.0
        with self.assertRaisesRegex(ValueError, r"tasks \[0, 1\] observe only one class"):
            model.fit(self.X, Y, self.mask, **self.penalties)
        # objective and kkt_violation build the same problem, so they refuse too.
        for check in (model.objective, model.kkt_violation):
            with self.subTest(check=check.__name__):
                with self.assertRaisesRegex(ValueError, r"tasks \[0, 1\] observe only one class"):
                    check(self.fit, self.X, Y, self.mask)
        # Scoring a validation fold may still meet a one-class task.
        self.assertEqual(model.loss_weights(Y, self.mask).shape, Y.shape)

    def test_intercept_gradient_is_within_kkt_tol_on_binary_features(self):
        # Uncentred 0/1 features and a rare positive class tie the intercepts
        # to the coefficients, and lambda2 = 0.5 slows them, so the objective
        # settles before the intercept gradient does.
        rng = np.random.default_rng(1)
        X = (rng.random((60, 4)) < 0.7).astype(float)
        Y = (rng.random((60, 3)) < 0.1).astype(float)
        result = model.fit(X, Y, lambda1=0.01, lambda_group=0.02, lambda2=0.5, kkt_tol=1e-7)
        self.assertEqual(result.status, model.CONVERGED)
        self.assertLessEqual(intercept_gradient(result, X, Y), 1e-7)

    def test_too_large_a_step_stalls_instead_of_converging(self):
        lipschitz = model._lipschitz
        with mock.patch.object(model, "_lipschitz", lambda problem: lipschitz(problem) / 10):
            result = model.fit(self.X, self.Y, self.mask, **self.penalties)
        self.assertEqual(result.status, model.UNCONVERGED)
        self.assertEqual(result.diagnostics["reason"], "stalled")
        with self.assertRaisesRegex(ValueError, "no candidate converged"):
            model.select_candidate([result], [0.1])

    def test_lambda2_must_be_positive(self):
        for lambda2 in (0.0, -0.1, math.nan, math.inf):
            with self.assertRaisesRegex(ValueError, "lambda2"):
                model.fit(self.X, self.Y, self.mask, lambda1=0.1, lambda_group=0.1, lambda2=lambda2)
        with self.assertRaisesRegex(ValueError, "lambda_group"):
            model.fit(self.X, self.Y, self.mask, lambda1=0.1, lambda_group=-1.0, lambda2=0.1)

    def test_fit_is_deterministic(self):
        again = model.fit(self.X, self.Y, self.mask, **self.penalties, **TIGHT)
        np.testing.assert_array_equal(again.coef, self.fit.coef)
        np.testing.assert_array_equal(again.intercept, self.fit.intercept)
        self.assertEqual(again.diagnostics, self.fit.diagnostics)

    def test_predict_proba_is_the_logistic_of_each_task(self):
        probability = self.fit.predict_proba(self.X[:4])
        self.assertEqual(probability.shape, (4, 2))
        np.testing.assert_allclose(probability, sigmoid(self.X[:4] @ self.fit.coef + self.fit.intercept),
                                   rtol=1e-14)


class WeightTests(unittest.TestCase):
    def setUp(self):
        self.Y = np.array([[1.0, 0.0], [0.0, np.nan], [1.0, 1.0], [1.0, np.nan]])

    def test_per_task_weights_are_one_over_observed_count(self):
        np.testing.assert_array_equal(model.loss_weights(self.Y),
                                      [[0.25, 0.5], [0.25, 0.0], [0.25, 0.5], [0.25, 0.0]])

    def test_unweighted_gives_every_label_the_same_weight_totalling_the_task_count(self):
        np.testing.assert_allclose(model.loss_weights(self.Y, normalisation="unweighted"),
                                   [[1 / 3, 1 / 3], [1 / 3, 0], [1 / 3, 1 / 3], [1 / 3, 0]], rtol=1e-15)
        full = np.nan_to_num(self.Y)
        np.testing.assert_array_equal(model.loss_weights(full, normalisation="unweighted"),
                                      model.loss_weights(full))

    def test_inverse_prevalence_splits_each_task_between_its_classes(self):
        np.testing.assert_allclose(model.loss_weights(self.Y, normalisation="inverse_prevalence"),
                                   [[1 / 6, 0.5], [0.5, 0], [1 / 6, 0.5], [1 / 6, 0]], rtol=1e-15)
        with self.assertRaisesRegex(ValueError, r"tasks \[0\] observe only one class"):
            model.loss_weights(np.ones((3, 1)), normalisation="inverse_prevalence")

    def test_rejects_bad_labels_and_empty_tasks(self):
        with self.assertRaisesRegex(ValueError, "0 or 1"):
            model.loss_weights(np.array([[0.0], [0.5]]))
        with self.assertRaisesRegex(ValueError, "no observed labels"):
            model.loss_weights(np.array([[0.0, np.nan], [1.0, np.nan]]))
        with self.assertRaisesRegex(ValueError, "normalisation"):
            model.loss_weights(np.zeros((2, 1)), normalisation="balanced")


if __name__ == "__main__":
    unittest.main()
