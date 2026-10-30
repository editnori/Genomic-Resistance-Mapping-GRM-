//! L1-regularized logistic regression on screened presence features.
//!
//! IRLS outer loop with warm starts down a lambda path. Each IRLS step's
//! weighted least-squares problem is solved exactly on its active set by
//! feature-sign search (Lee et al. 2007), and the step is halved while it
//! raises the objective. Each feature's penalty is scaled by its standard
//! deviation on the fitting samples, which is the raw-scale form of fitting
//! standardized features: rare and common features compete on equal terms.
//! Class-balancing weights match the tree learners.
//!
//! A point converges when the largest violation of its optimality conditions,
//! the KKT residual, is at most 1e-6. The IRLS and active-set loops have
//! guards; a point that stops at one is not a converged fit, and its
//! `Convergence` says so. See docs/lasso-stopping.md.

use crate::bits::{self, BitMatrix, Bits};
use rustc_hash::FxHashMap;
use serde::{Deserialize, Serialize};
use std::time::Instant;

/// Guards for one penalty point. Neither is the stopping rule: a point stops
/// when its KKT residual reaches `KKT_TOL`, and one that reaches a guard
/// first is unconverged.
#[derive(Clone, Copy, Debug)]
struct Limits {
    irls_steps: usize,
    /// Active-set iterations in one IRLS step's least-squares solve.
    inner_iterations: usize,
}

/// On the saved fixture panels and MTB delamanid fold 0 a point takes at
/// most 7 IRLS steps (C. jejuni lineage fold 0, point 4) and a solve at most
/// 44 iterations.
const LIMITS: Limits = Limits {
    irls_steps: 100,
    inner_iterations: 5_000,
};

/// A point converges when no optimality condition of the penalized objective
/// is violated by more than this on the standardized scale. The saved
/// fixture panels are accepted at the same tolerance.
const KKT_TOL: f64 = 1e-6;

/// One IRLS step's least-squares solve stops when no zero coefficient
/// violates its condition by more than this, 1,000 times below `KKT_TOL`.
const INNER_TOL: f64 = 1e-9;

/// Lower bound on the IRLS weight p(1 - p).
const WEIGHT_FLOOR: f64 = 1e-5;

/// A column whose Cholesky pivot is below this share of its own weighted sum
/// of squares counts as a linear combination of the active columns.
const PIVOT_TOL: f64 = 1e-10;

/// Halvings of an IRLS step that raises the objective before the point is
/// given up.
const HALVINGS: u32 = 30;

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct LassoModel {
    pub intercept: f64,
    /// Non-zero coefficients as (feature, weight).
    pub coef: Vec<(u32, f64)>,
    pub lambda: f64,
}

/// How the solve at one penalty point stopped.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Convergence {
    /// IRLS steps taken; 0 when the warm start already met the tolerance.
    pub irls_steps: u32,
    /// IRLS steps whose active-set solve stopped short of its own optimality
    /// test, at its iteration guard or where rounding left it no descent.
    /// The name dates from the coordinate-descent solver's sweep cap.
    pub sweep_cap_hits: u32,
    /// The KKT residual met its tolerance on a two-class fit with finite
    /// parameters. False when a guard or a failed step search ended the
    /// solve.
    pub converged: bool,
}

#[derive(Clone, Debug)]
pub struct LassoFit {
    pub model: LassoModel,
    pub convergence: Convergence,
}

/// One penalty point, reported as soon as it is solved.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PointProgress {
    pub point: usize,
    pub lambda: f64,
    pub irls_steps: u32,
    /// Non-zero coefficients.
    pub active: usize,
    /// KKT residual at the point's final parameters.
    pub kkt: f64,
    /// Wall-clock seconds spent on this point.
    pub seconds: f64,
    pub converged: bool,
}

impl LassoModel {
    pub fn predict_proba(&self, x: &BitMatrix) -> Vec<f64> {
        let mut eta = vec![self.intercept; x.cols];
        for &(f, w) in &self.coef {
            let row = x.row(f as usize);
            for (wi, &word) in row.iter().enumerate() {
                let mut word = word;
                while word != 0 {
                    eta[wi * 64 + word.trailing_zeros() as usize] += w;
                    word &= word - 1;
                }
            }
        }
        eta.into_iter().map(super::gbt::sigmoid).collect()
    }
}

struct Problem {
    /// Fitting-sample indices (into the full sample space) per feature.
    columns: Vec<Vec<u32>>,
    features: Vec<u32>,
    sd: Vec<f64>,
    samples: Vec<usize>,
    y: Vec<f64>,
    class_weight: Vec<f64>,
}

impl Problem {
    fn new(x: &BitMatrix, y: &Bits, fit: &Bits, features: &[u32]) -> Problem {
        let samples: Vec<usize> = fit.iter_ones().collect();
        let mut local = vec![u32::MAX; x.cols];
        for (j, &i) in samples.iter().enumerate() {
            local[i] = j as u32;
        }
        let n = samples.len() as f64;
        let columns: Vec<Vec<u32>> = features
            .iter()
            .map(|&f| {
                let mut col = Vec::new();
                let row = x.row(f as usize);
                for (wi, (&r, &m)) in row.iter().zip(fit.words()).enumerate() {
                    let mut w = r & m;
                    while w != 0 {
                        col.push(local[wi * 64 + w.trailing_zeros() as usize]);
                        w &= w - 1;
                    }
                }
                // A feature on every fitting sample has no variance to
                // standardize. Fitted, it would act as a second intercept
                // and shift every sample that lacks it.
                if col.len() == samples.len() {
                    col.clear();
                }
                col
            })
            .collect();
        let sd = columns
            .iter()
            .map(|c| {
                let p = c.len() as f64 / n;
                (p * (1.0 - p)).sqrt().max(1e-6)
            })
            .collect();
        let yv: Vec<f64> = samples
            .iter()
            .map(|&i| y.contains(i) as u8 as f64)
            .collect();
        let pos = yv.iter().sum::<f64>();
        let class_weight = yv
            .iter()
            .map(|&t| {
                if t == 1.0 {
                    n / (2.0 * pos.max(1.0))
                } else {
                    n / (2.0 * (n - pos).max(1.0))
                }
            })
            .collect();
        Problem {
            columns,
            features: features.to_vec(),
            sd,
            samples,
            y: yv,
            class_weight,
        }
    }

    fn lambda_max(&self) -> f64 {
        // At beta = 0 with balanced weights the fitted probability is 0.5.
        let n = self.samples.len() as f64;
        self.columns
            .iter()
            .zip(&self.sd)
            .map(|(c, sd)| {
                let grad: f64 = c
                    .iter()
                    .map(|&j| self.class_weight[j as usize] * (self.y[j as usize] - 0.5))
                    .sum();
                grad.abs() / n / sd
            })
            .fold(0.0, f64::max)
    }
}

/// Fits the path at each lambda (descending), warm-starting each from the
/// last. Stops early once more than `max_nonzero` coefficients are active.
/// Points that stopped unconverged are reported on stderr; use
/// `fit_path_with_status` to act on them.
pub fn fit_path(
    x: &BitMatrix,
    y: &Bits,
    fit: &Bits,
    features: &[u32],
    n_lambdas: usize,
    max_nonzero: usize,
) -> Vec<LassoModel> {
    report_unconverged(
        "fit_path",
        fit_path_with_status(x, y, fit, features, n_lambdas, max_nonzero, &mut |_| {}),
    )
}

/// `fit_path`, with how each point stopped. `observe` sees each point as
/// soon as it is solved.
pub fn fit_path_with_status(
    x: &BitMatrix,
    y: &Bits,
    fit: &Bits,
    features: &[u32],
    n_lambdas: usize,
    max_nonzero: usize,
    observe: &mut dyn FnMut(&PointProgress),
) -> Vec<LassoFit> {
    if !has_both_classes(y, fit) {
        return vec![];
    }
    path_problem(x, y, fit, features, n_lambdas).map_or_else(Vec::new, |(prob, lambdas)| {
        solve_path(&prob, &lambdas, max_nonzero, LIMITS, observe)
    })
}

/// One-class logistic regression has no finite intercept optimum. Do not run a
/// long divergent path whose final status could depend on floating-point
/// overflow rather than convergence. The evaluator also requires both classes
/// in every fitting and tuning partition.
fn has_both_classes(y: &Bits, fit: &Bits) -> bool {
    let pos = fit.and(y).count();
    pos > 0 && pos < fit.count()
}

/// The problem on `fit` and its penalty grid: `n_lambdas` log-spaced points
/// from the smallest penalty that keeps every coefficient at zero down to
/// 1e-3 of it. None when there is nothing to fit.
fn path_problem(
    x: &BitMatrix,
    y: &Bits,
    fit: &Bits,
    features: &[u32],
    n_lambdas: usize,
) -> Option<(Problem, Vec<f64>)> {
    if n_lambdas == 0 || fit.count() == 0 {
        return None;
    }
    let prob = Problem::new(x, y, fit, features);
    let lmax = prob.lambda_max();
    if lmax == 0.0 {
        return None;
    }
    let ratio: f64 = 1e-3;
    let lambdas: Vec<f64> = (0..n_lambdas)
        .map(|k| lmax * ratio.powf(k as f64 / (n_lambdas - 1).max(1) as f64))
        .collect();
    Some((prob, lambdas))
}

/// `fit_path_with_status` under other limits, and without its one-class
/// guard, so the fixture accounting keeps the raw one-class trace.
#[cfg(test)]
fn fit_path_limited(
    x: &BitMatrix,
    y: &Bits,
    fit: &Bits,
    features: &[u32],
    n_lambdas: usize,
    max_nonzero: usize,
    limits: Limits,
) -> Vec<LassoFit> {
    path_problem(x, y, fit, features, n_lambdas).map_or_else(Vec::new, |(prob, lambdas)| {
        solve_path(&prob, &lambdas, max_nonzero, limits, &mut |_| {})
    })
}

/// Refit the selected penalty path without recomputing its scale or spacing.
/// Points that stopped unconverged are reported on stderr; use
/// `refit_path_with_status` to act on them.
pub fn refit_path(
    x: &BitMatrix,
    y: &Bits,
    fit: &Bits,
    features: &[u32],
    lambdas: &[f64],
) -> Vec<LassoModel> {
    report_unconverged(
        "refit_path",
        refit_path_with_status(x, y, fit, features, lambdas, &mut |_| {}),
    )
}

/// `refit_path`, with how each point stopped. `observe` sees each point as
/// soon as it is solved.
pub fn refit_path_with_status(
    x: &BitMatrix,
    y: &Bits,
    fit: &Bits,
    features: &[u32],
    lambdas: &[f64],
    observe: &mut dyn FnMut(&PointProgress),
) -> Vec<LassoFit> {
    assert!(lambdas.iter().all(|v| v.is_finite() && *v > 0.0));
    if fit.count() == 0 || !has_both_classes(y, fit) {
        return vec![];
    }
    solve_path(
        &Problem::new(x, y, fit, features),
        lambdas,
        usize::MAX,
        LIMITS,
        observe,
    )
}

fn report_unconverged(caller: &str, path: Vec<LassoFit>) -> Vec<LassoModel> {
    for (i, f) in path.iter().enumerate() {
        let c = f.convergence;
        if !c.converged {
            eprintln!(
                "lasso {caller}: point {i} (lambda {:e}) unconverged: {} IRLS steps (guard {}), \
                 {} active-set solves stopped short (guard {} iterations)",
                f.model.lambda,
                c.irls_steps,
                LIMITS.irls_steps,
                c.sweep_cap_hits,
                LIMITS.inner_iterations
            );
        }
    }
    path.into_iter().map(|f| f.model).collect()
}

fn solve_path(
    prob: &Problem,
    lambdas: &[f64],
    max_nonzero: usize,
    limits: Limits,
    observe: &mut dyn FnMut(&PointProgress),
) -> Vec<LassoFit> {
    let rep = representatives(&prob.columns);
    // One-class logistic regression has no finite optimum: its gradient
    // vanishes only as the intercept diverges, so a small KKT residual
    // certifies nothing there.
    let two_class = prob.y.contains(&1.0) && prob.y.contains(&0.0);

    let mut beta = vec![0.0; prob.columns.len()];
    let mut b0 = 0.0;
    let mut eta = vec![0.0; prob.samples.len()];
    let mut path = Vec::new();

    for (point, &lambda) in lambdas.iter().enumerate() {
        let started = Instant::now();
        let mut convergence = Convergence {
            irls_steps: 0,
            sweep_cap_hits: 0,
            converged: false,
        };
        let mut kkt = kkt_residual(prob, &rep, lambda, &beta, &eta);
        while (kkt.is_nan() || kkt > KKT_TOL)
            && (convergence.irls_steps as usize) < limits.irls_steps
        {
            convergence.irls_steps += 1;
            let Some(mut wls) = Wls::new(prob, &rep, lambda, &eta) else {
                convergence.sweep_cap_hits += 1;
                break;
            };
            let (mut next_b0, mut next_beta) = (b0, beta.clone());
            if !wls.solve(&mut next_b0, &mut next_beta, limits.inner_iterations) {
                convergence.sweep_cap_hits += 1;
            }
            if !take_step(
                prob, lambda, &mut b0, &mut beta, &mut eta, next_b0, &next_beta,
            ) {
                break;
            }
            kkt = kkt_residual(prob, &rep, lambda, &beta, &eta);
        }
        convergence.converged = two_class && kkt <= KKT_TOL;
        let convergence = require_finite(convergence, b0, &beta);
        let coef: Vec<(u32, f64)> = beta
            .iter()
            .enumerate()
            .filter(|(_, b)| **b != 0.0)
            .map(|(k, b)| (prob.features[k], *b))
            .collect();
        // Exact copies share one coefficient (`representatives`), so this
        // counts distinct columns.
        let active = coef.len();
        observe(&PointProgress {
            point,
            lambda,
            irls_steps: convergence.irls_steps,
            active,
            kkt,
            seconds: started.elapsed().as_secs_f64(),
            converged: convergence.converged,
        });
        path.push(LassoFit {
            model: LassoModel {
                intercept: b0,
                coef,
                lambda,
            },
            convergence,
        });
        if active > max_nonzero {
            break;
        }
    }
    path
}

fn require_finite(
    mut convergence: Convergence,
    intercept: f64,
    coefficients: &[f64],
) -> Convergence {
    convergence.converged &= intercept.is_finite() && coefficients.iter().all(|b| b.is_finite());
    convergence
}

/// Each column's representative: the lowest-index column identical to it on
/// the fitting samples. Only representatives enter the active set, so exact
/// copies share one coefficient, kept on the first. Copies have the same
/// penalty factor, so every same-sign split of that weight among them has
/// the same loss and penalty, and the optimal objective does not change.
fn representatives(columns: &[Vec<u32>]) -> Vec<usize> {
    let mut first: FxHashMap<&[u32], usize> = FxHashMap::default();
    columns
        .iter()
        .enumerate()
        .map(|(k, c)| *first.entry(c.as_slice()).or_insert(k))
        .collect()
}

/// The largest violation of the penalized objective's optimality conditions
/// at `beta` with linear predictor `eta`, on the standardized scale: the
/// intercept's gradient, and for each column with a carrier, with g its
/// gradient over its standard deviation, |g + lambda sign(b)| when its
/// coefficient b is non-zero and max(|g| - lambda, 0) when it is zero.
/// NaN when any term is NaN.
fn kkt_residual(prob: &Problem, rep: &[usize], lambda: f64, beta: &[f64], eta: &[f64]) -> f64 {
    let n = eta.len() as f64;
    let res: Vec<f64> = eta
        .iter()
        .zip(&prob.y)
        .zip(&prob.class_weight)
        .map(|((&e, &t), &cw)| cw * (super::gbt::sigmoid(e) - t) / n)
        .collect();
    let mut worst = res.iter().sum::<f64>().abs();
    // f64::max would drop a NaN, here and below.
    if worst.is_nan() {
        return f64::NAN;
    }
    let mut grad = vec![0.0; beta.len()];
    for (k, col) in prob.columns.iter().enumerate() {
        if col.is_empty() {
            continue;
        }
        // A copy's representative comes first, so its gradient is ready.
        if rep[k] == k {
            grad[k] = col.iter().map(|&j| res[j as usize]).sum();
        }
        let g = grad[rep[k]] / prob.sd[k];
        let v = if beta[k] != 0.0 {
            (g + lambda * beta[k].signum()).abs()
        } else {
            g.abs() - lambda
        };
        if v.is_nan() {
            return f64::NAN;
        }
        worst = worst.max(v);
    }
    worst
}

/// Mean class-weighted log loss on the fitting samples plus lambda times the
/// standard-deviation-weighted L1 norm.
fn objective(prob: &Problem, lambda: f64, beta: &[f64], eta: &[f64]) -> f64 {
    let loss: f64 = eta
        .iter()
        .zip(&prob.y)
        .zip(&prob.class_weight)
        .map(|((&e, &t), &cw)| {
            let z = if t == 1.0 { -e } else { e };
            cw * (z.max(0.0) + (-z.abs()).exp().ln_1p())
        })
        .sum();
    let penalty: f64 = beta.iter().zip(&prob.sd).map(|(b, sd)| sd * b.abs()).sum();
    loss / eta.len() as f64 + lambda * penalty
}

fn linear_predictor(prob: &Problem, b0: f64, beta: &[f64], eta: &mut [f64]) {
    eta.fill(b0);
    for (col, &b) in prob.columns.iter().zip(beta) {
        if b != 0.0 {
            for &j in col {
                eta[j as usize] += b;
            }
        }
    }
}

/// Moves (b0, beta) toward the IRLS solution (next_b0, next_beta), halving
/// the step while it raises the penalized objective: near separation a full
/// Newton step can overshoot. Keeps `eta` the linear predictor. False, with
/// nothing moved, when the solution is the current point or no halving
/// lowers the objective.
fn take_step(
    prob: &Problem,
    lambda: f64,
    b0: &mut f64,
    beta: &mut Vec<f64>,
    eta: &mut Vec<f64>,
    next_b0: f64,
    next_beta: &[f64],
) -> bool {
    if next_b0 == *b0 && next_beta == beta.as_slice() {
        return false;
    }
    let current = objective(prob, lambda, beta, eta);
    // Rounding in the sum over samples; near the optimum a step that lowers
    // the KKT residual can change the objective by less than this.
    let slack = 1e-12 * current.abs().max(1.0);
    let (mut trial_b0, mut trial) = (next_b0, next_beta.to_vec());
    let mut trial_eta = vec![0.0; eta.len()];
    let mut t = 1.0;
    for _ in 0..=HALVINGS {
        linear_predictor(prob, trial_b0, &trial, &mut trial_eta);
        if objective(prob, lambda, &trial, &trial_eta) <= current + slack {
            *b0 = trial_b0;
            *beta = trial;
            *eta = trial_eta;
            return true;
        }
        t *= 0.5;
        trial_b0 = *b0 + t * (next_b0 - *b0);
        for ((v, &a), &b) in trial.iter_mut().zip(beta.iter()).zip(next_beta) {
            *v = a + t * (b - a);
        }
    }
    false
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    // Four running sums let the compiler vectorize the loop.
    let mut acc = [0.0; 4];
    let ((ca, ta), (cb, tb)) = (a.as_chunks::<4>(), b.as_chunks::<4>());
    let tail: f64 = ta.iter().zip(tb).map(|(x, y)| x * y).sum();
    for (x, y) in ca.iter().zip(cb) {
        for ((s, u), v) in acc.iter_mut().zip(x).zip(y) {
            *s += u * v;
        }
    }
    acc[0] + acc[1] + acc[2] + acc[3] + tail
}

/// Where coefficient c reaches zero on c + t d, 0 < t <= 1.
fn crossing(c: f64, d: f64) -> Option<f64> {
    (c != 0.0 && c * (c + d) <= 0.0).then(|| -c / d)
}

/// The first point t < 1 where a non-zero coefficient of c + t d, past the
/// unpenalized intercept at index 0, reaches zero.
fn first_crossing(c: &[f64], d: &[f64]) -> Option<f64> {
    (1..c.len())
        .filter_map(|i| crossing(c[i], d[i]))
        .filter(|&t| t < 1.0)
        .min_by(f64::total_cmp)
}

/// The feature-sign line search on c + t d, with index 0 the unpenalized
/// intercept and `pen[i - 1]` the penalty of coefficient i. It evaluates the
/// least-squares objective's change at each point where a non-zero
/// coefficient reaches zero (that coefficient set to exactly zero) and at the
/// full step; the smooth part changes by t slope + t^2 curv / 2. Returns the
/// step with the lowest change, and that change.
fn line_search(c: &[f64], d: &[f64], slope: f64, curv: f64, pen: &[f64]) -> (f64, f64) {
    let change = |t: f64| -> f64 {
        let penalty: f64 = (1..c.len())
            .map(|i| {
                let v = if crossing(c[i], d[i]) == Some(t) {
                    0.0
                } else {
                    c[i] + t * d[i]
                };
                pen[i - 1] * (v.abs() - c[i].abs())
            })
            .sum();
        t * slope + 0.5 * t * t * curv + penalty
    };
    let mut best = (1.0, change(1.0));
    for t in (1..c.len())
        .filter_map(|i| crossing(c[i], d[i]))
        .filter(|&t| t < 1.0)
    {
        let dq = change(t);
        if dq < best.1 {
            best = (t, dq);
        }
    }
    best
}

/// One IRLS step's weighted least-squares problem,
///
///   min over (b0, beta) of 1/2 sum_j w_j (z_j - b0 - x_j beta)^2 + sum_k pen_k |beta_k|,
///
/// with pen_k = lambda sd_k and the intercept unpenalized, solved exactly by
/// feature-sign search on a Cholesky-factored active set. The weights are
/// w_j = c_j max(p_j (1 - p_j), 1e-5) / n for class weight c_j, so at the
/// current point the gradient equals the penalized objective's and the solve's
/// violations are on the KKT scale.
struct Wls<'a> {
    prob: &'a Problem,
    rep: &'a [usize],
    lambda: f64,
    w: Vec<f64>,
    sw: Vec<f64>,
    z: Vec<f64>,
    /// Active columns. Position i is index i + 1 of the Gram matrix and its
    /// factor; index 0 is the intercept.
    active: Vec<usize>,
    /// The sign fixed for each active coefficient: its own when non-zero, its
    /// gradient's while it enters at zero.
    sign: Vec<f64>,
    /// sqrt(w) times the intercept and each active column, dense.
    rows: Vec<Vec<f64>>,
    gram: Vec<Vec<f64>>,
    /// Lower-triangular Cholesky factor of `gram`, row by row.
    chol: Vec<Vec<f64>>,
    /// X'Wz for the intercept and each active column.
    xwz: Vec<f64>,
}

/// A column that cannot join the factor: L^-1 times its Gram products with
/// the intercept and the active columns, and its pivot.
struct Dependent {
    l: Vec<f64>,
    pivot: f64,
}

impl<'a> Wls<'a> {
    /// The problem linearized at `eta`, with only the intercept factored.
    /// None when the weights do not sum to a positive finite number.
    fn new(prob: &'a Problem, rep: &'a [usize], lambda: f64, eta: &[f64]) -> Option<Self> {
        let n = eta.len() as f64;
        let (mut w, mut z) = (Vec::with_capacity(eta.len()), Vec::with_capacity(eta.len()));
        for (j, &e) in eta.iter().enumerate() {
            let p = super::gbt::sigmoid(e);
            let v = (p * (1.0 - p)).max(WEIGHT_FLOOR);
            w.push(prob.class_weight[j] * v / n);
            z.push(e + (prob.y[j] - p) / v);
        }
        let sw: Vec<f64> = w.iter().map(|v: &f64| v.sqrt()).collect();
        let total = dot(&sw, &sw);
        if !(total > 0.0 && total.is_finite()) {
            return None;
        }
        let xwz = vec![dot(&w, &z)];
        Some(Wls {
            prob,
            rep,
            lambda,
            rows: vec![sw.clone()],
            w,
            sw,
            z,
            active: Vec::new(),
            sign: Vec::new(),
            gram: vec![vec![total]],
            chol: vec![vec![total.sqrt()]],
            xwz,
        })
    }

    fn pen(&self, k: usize) -> f64 {
        self.lambda * self.prob.sd[k]
    }

    /// Solves L y = b.
    fn forward(&self, b: &[f64]) -> Vec<f64> {
        let mut y = Vec::with_capacity(b.len());
        for (i, row) in self.chol.iter().enumerate() {
            let s = b[i] - dot(&row[..i], &y);
            y.push(s / row[i]);
        }
        y
    }

    /// Solves L' x = y.
    fn backward(&self, y: &[f64]) -> Vec<f64> {
        let m = y.len();
        let mut x = vec![0.0; m];
        for i in (0..m).rev() {
            let mut s = y[i];
            for (row, xj) in self.chol[i + 1..].iter().zip(&x[i + 1..]) {
                s -= row[i] * xj;
            }
            x[i] = s / self.chol[i][i];
        }
        x
    }

    fn gram_times(&self, v: &[f64]) -> Vec<f64> {
        self.gram.iter().map(|g| dot(g, v)).collect()
    }

    /// The intercept, then each active coefficient.
    fn coefficients(&self, b0: f64, beta: &[f64]) -> Vec<f64> {
        std::iter::once(b0)
            .chain(self.active.iter().map(|&k| beta[k]))
            .collect()
    }

    /// Adds column k to the active set with `sign` fixed, extending the Gram
    /// matrix and its factor by one row. Refuses, leaving everything as it
    /// was, when the new pivot is at most `tol` times the column's weighted
    /// sum of squares.
    fn push(&mut self, k: usize, sign: f64, tol: f64) -> Result<(), Dependent> {
        let col = &self.prob.columns[k];
        let mut row = vec![0.0; self.sw.len()];
        for &j in col {
            row[j as usize] = self.sw[j as usize];
        }
        let mut cross: Vec<f64> = self.rows.iter().map(|r| dot(r, &row)).collect();
        let diag = dot(&row, &row);
        let mut l = self.forward(&cross);
        let pivot = diag - dot(&l, &l);
        if pivot.is_nan() || pivot <= tol * diag {
            return Err(Dependent { l, pivot });
        }
        for (g, &c) in self.gram.iter_mut().zip(&cross) {
            g.push(c);
        }
        cross.push(diag);
        self.gram.push(cross);
        l.push(pivot.sqrt());
        self.chol.push(l);
        self.rows.push(row);
        self.xwz.push(
            col.iter()
                .map(|&j| self.w[j as usize] * self.z[j as usize])
                .sum(),
        );
        self.active.push(k);
        self.sign.push(sign);
        Ok(())
    }

    /// Drops the active column at position `pos`. Removing its row from the
    /// factor leaves one entry past the diagonal in each later row; Givens
    /// rotations on adjacent columns clear them and keep L L' = G.
    fn remove(&mut self, pos: usize) {
        let q = pos + 1;
        self.active.remove(pos);
        self.sign.remove(pos);
        self.rows.remove(q);
        self.xwz.remove(q);
        self.gram.remove(q);
        for g in &mut self.gram {
            g.remove(q);
        }
        self.chol.remove(q);
        for i in q..self.chol.len() {
            let (a, b) = (self.chol[i][i], self.chol[i][i + 1]);
            let r = a.hypot(b);
            let (c, s) = (a / r, b / r);
            for row in &mut self.chol[i..] {
                let (u, v) = (row[i], row[i + 1]);
                row[i] = c * u + s * v;
                row[i + 1] = c * v - s * u;
            }
            self.chol[i].pop();
        }
    }

    fn position(&self, k: usize) -> usize {
        self.active
            .iter()
            .position(|&j| j == k)
            .expect("entering columns are active")
    }

    /// Back to the intercept alone.
    fn reset(&mut self) {
        self.active.clear();
        self.sign.clear();
        self.rows.truncate(1);
        self.xwz.truncate(1);
        self.gram.truncate(1);
        self.gram[0].truncate(1);
        self.chol.truncate(1);
    }

    /// The minimizer over the active set with every sign held fixed:
    /// G c = X'Wz - pen sign, the intercept unpenalized.
    fn sign_fixed(&self) -> Vec<f64> {
        let mut rhs = self.xwz.clone();
        for (i, (&k, &s)) in self.active.iter().zip(&self.sign).enumerate() {
            rhs[i + 1] -= self.pen(k) * s;
        }
        self.backward(&self.forward(&rhs))
    }

    /// Feature-sign search from (b0, beta), in place. True when it stops at
    /// its optimality test; false at the iteration guard or where rounding
    /// leaves it no descent.
    fn solve(&mut self, b0: &mut f64, beta: &mut [f64], limit: usize) -> bool {
        // Warm start on the non-zero coefficients. One that these weights
        // make numerically dependent on the others starts at zero instead;
        // the IRLS step search compares the result with the current point,
        // not with this start.
        for (k, b) in beta.iter_mut().enumerate() {
            if *b != 0.0 && self.push(k, b.signum(), PIVOT_TOL).is_err() {
                *b = 0.0;
            }
        }
        let mut blocked = vec![false; beta.len()];
        let mut exact = true;
        // Active columns still at zero, most violating first.
        let mut entering: Vec<usize> = Vec::new();
        for _ in 0..limit {
            let c = self.coefficients(*b0, beta);
            let d: Vec<f64> = self
                .sign_fixed()
                .iter()
                .zip(&c)
                .map(|(a, b)| a - b)
                .collect();
            let slope: f64 = self
                .gram_times(&c)
                .iter()
                .zip(&self.xwz)
                .zip(&d)
                .map(|((g, b), di)| (g - b) * di)
                .sum();
            let curv = dot(&d, &self.gram_times(&d));
            if !(slope.is_finite() && curv.is_finite()) {
                return false;
            }
            let pen: Vec<f64> = self.active.iter().map(|&k| self.pen(k)).collect();
            let (t, change) = line_search(&c, &d, slope, curv, &pen);
            // Every non-zero active coefficient has its fixed sign. When each
            // entering one moves in the direction of its own, the objective
            // equals the sign-fixed quadratic along d up to the first point
            // where a non-zero coefficient reaches zero, so it falls all the
            // way there. A coefficient within rounding of zero crosses so
            // early that the fall is below rounding and no step shows a
            // change below zero. Step to that crossing anyway, so that it
            // leaves; otherwise it blocks every step, and each entering
            // column is blocked in turn or the zero columns are checked at a
            // point that is not stationary.
            let with_signs = self
                .active
                .iter()
                .zip(&self.sign)
                .zip(&d[1..])
                .all(|((&k, &s), &dk)| beta[k] != 0.0 || dk * s > 0.0);
            let step = if change < 0.0 {
                Some(t)
            } else if with_signs {
                first_crossing(&c, &d)
            } else {
                None
            };
            if let Some(t) = step {
                *b0 = c[0] + t * d[0];
                for (i, &k) in self.active.iter().enumerate() {
                    let (ci, di) = (c[i + 1], d[i + 1]);
                    beta[k] = if crossing(ci, di) == Some(t) {
                        0.0
                    } else {
                        ci + t * di
                    };
                }
                let flipped = self
                    .active
                    .iter()
                    .zip(&self.sign)
                    .any(|(&k, &s)| beta[k] != 0.0 && beta[k].signum() != s);
                for i in (0..self.active.len()).rev() {
                    let k = self.active[i];
                    if beta[k] == 0.0 {
                        self.remove(i);
                    } else {
                        self.sign[i] = beta[k].signum();
                    }
                }
                entering.clear();
                if t < 1.0 || flipped {
                    continue;
                }
                // The full step with every sign as fixed: each active
                // coefficient is stationary.
            } else if entering.len() > 1 {
                // A batch whose step does not lower the objective. From a
                // stationary point, one entering column's step always does,
                // so keep only the most violating.
                for &k in &entering[1..] {
                    let i = self.position(k);
                    self.remove(i);
                }
                entering.truncate(1);
                continue;
            } else if let Some(k) = entering.pop() {
                // Even that step does not, so rounding has the last word on
                // this column for the rest of this solve.
                let i = self.position(k);
                self.remove(i);
                blocked[k] = true;
                exact = false;
                continue;
            }
            // Every active coefficient is stationary; check the zero ones.
            let Some(violators) = self.violators(*b0, beta, &blocked) else {
                return false;
            };
            let Some(&(k, top, g)) = violators.first() else {
                return exact;
            };
            let room = self.active.len() / 10 + 1;
            match self.push(k, g.signum(), PIVOT_TOL) {
                Ok(()) => {
                    entering.push(k);
                    for &(k, _, g) in violators[1..]
                        .iter()
                        .take_while(|v| v.1 > top / 2.0)
                        .take(room - 1)
                    {
                        if self.push(k, g.signum(), PIVOT_TOL).is_ok() {
                            entering.push(k);
                        }
                    }
                }
                Err(dep) => {
                    if !self.exchange(k, g, dep, b0, beta) {
                        blocked[k] = true;
                        exact = false;
                    }
                }
            }
        }
        false
    }

    /// Zero, unblocked representative columns whose condition
    /// |g_k| <= pen_k fails by more than `INNER_TOL` on the standardized
    /// scale, as (column, violation, g_k), most violating first; g_k is the
    /// least-squares objective's descent rate along column k. None when a
    /// gradient is NaN.
    fn violators(&self, b0: f64, beta: &[f64], blocked: &[bool]) -> Option<Vec<(usize, f64, f64)>> {
        let mut fit = vec![b0; self.w.len()];
        for &k in &self.active {
            for &j in &self.prob.columns[k] {
                fit[j as usize] += beta[k];
            }
        }
        let q: Vec<f64> = self
            .w
            .iter()
            .zip(&self.z)
            .zip(&fit)
            .map(|((w, z), f)| w * (z - f))
            .collect();
        let mut out = Vec::new();
        for (k, col) in self.prob.columns.iter().enumerate() {
            if self.rep[k] != k || col.is_empty() || beta[k] != 0.0 || blocked[k] {
                continue;
            }
            let g: f64 = col.iter().map(|&j| q[j as usize]).sum();
            let v = (g.abs() - self.pen(k)) / self.prob.sd[k];
            if v.is_nan() {
                return None;
            }
            if v > INNER_TOL {
                out.push((k, v, g));
            }
        }
        // A tie goes to the lowest column, so the order is a property of the
        // violations alone, not of how the list was built.
        out.sort_unstable_by(|a, b| b.1.total_cmp(&a.1).then(a.0.cmp(&b.0)));
        Some(out)
    }

    /// Column k violates its condition but is, to rounding, a linear
    /// combination of the intercept and the active columns, so it cannot
    /// join the factor. Moving beta_k by theta t and the others by
    /// -theta t u, where G_A u is k's Gram column and theta the sign of g,
    /// leaves the fit unchanged up to the pivot and, from a stationary point,
    /// lowers the objective at the rate pen_k - |g| < 0. Follow that
    /// direction to the first active coefficient that reaches zero, which
    /// leaves so that k can join, or to the minimum when the pivot is
    /// positive. False, with nothing changed, unless the exchange strictly
    /// improves: its rate must be below zero by more than `INNER_TOL` on the
    /// standardized scale, the margin a zero column needs to enter at all.
    fn exchange(
        &mut self,
        k: usize,
        g: f64,
        dep: Dependent,
        b0: &mut f64,
        beta: &mut [f64],
    ) -> bool {
        let theta = g.signum();
        let v: Vec<f64> = self.backward(&dep.l).iter().map(|u| -theta * u).collect();
        let c = self.coefficients(*b0, beta);
        let grad: Vec<f64> = self
            .gram_times(&c)
            .iter()
            .zip(&self.xwz)
            .map(|(a, b)| a - b)
            .collect();
        let penalty_rate: f64 = self
            .active
            .iter()
            .zip(&self.sign)
            .zip(&v[1..])
            .map(|((&j, &s), &vi)| self.pen(j) * s * vi)
            .sum();
        let rate = dot(&grad, &v) - theta * g + penalty_rate + self.pen(k);
        // The smooth part of the rate is zero up to rounding, since the fit
        // does not change, so the rate is pen_k less the active penalties it
        // displaces. An exact complement of an active column (x_k = 1 - x_j,
        // the same standard deviation) displaces exactly its own penalty, as
        // a copy would without the collapse: its rate is zero, and acting
        // on the sign rounding gives it would swap the two back and forth.
        if rate.is_nan() || rate >= -INNER_TOL * self.prob.sd[k] {
            return false;
        }
        let curv = dep.pivot.max(0.0);
        let mut t = if curv > 0.0 {
            -rate / curv
        } else {
            f64::INFINITY
        };
        for (&ci, &vi) in c[1..].iter().zip(&v[1..]) {
            if ci * vi < 0.0 && -ci / vi < t {
                t = -ci / vi;
            }
        }
        if !t.is_finite() {
            return false;
        }
        let saved: Vec<(usize, f64)> = self.active.iter().map(|&j| (j, beta[j])).collect();
        *b0 = c[0] + t * v[0];
        for (i, &j) in self.active.iter().enumerate() {
            let (ci, next) = (c[i + 1], c[i + 1] + t * v[i + 1]);
            // Every coefficient that reaches zero at t, to rounding, leaves:
            // the first one, any that ties with it, and any that rounding
            // carries past zero, so each one still active keeps its fixed
            // sign.
            let leaves = ci * next <= 0.0 || next.abs() <= 8.0 * f64::EPSILON * ci.abs();
            beta[j] = if leaves { 0.0 } else { next };
        }
        for i in (0..self.active.len()).rev() {
            if beta[self.active[i]] == 0.0 {
                self.remove(i);
            }
        }
        beta[k] = theta * t;
        if self.push(k, theta, 0.0).is_ok() {
            return true;
        }
        // Rounding leaves k dependent even so: undo.
        beta[k] = 0.0;
        *b0 = c[0];
        self.reset();
        for (j, b) in saved {
            beta[j] = b;
            if self.push(j, b.signum(), 0.0).is_err() {
                beta[j] = 0.0;
            }
        }
        false
    }
}

/// Mean class-balanced log loss on `mask`.
pub fn log_loss(scores: &[f64], y: &Bits, mask: &Bits) -> f64 {
    let pos = mask.and(y).count() as f64;
    let n = mask.count() as f64;
    let mut total = 0.0;
    for i in mask.iter_ones() {
        let p = scores[i].clamp(1e-9, 1.0 - 1e-9);
        total += if y.contains(i) {
            -p.ln() * n / (2.0 * pos.max(1.0))
        } else {
            -(1.0 - p).ln() * n / (2.0 * (n - pos).max(1.0))
        };
    }
    total / n.max(1.0)
}

/// Non-zero count of a feature in `mask` (used to report support).
pub fn support(x: &BitMatrix, f: u32, mask: &Bits) -> u32 {
    bits::and_count(x.row(f as usize), mask.words())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn refit_preserves_the_selected_lambda() {
        let n = 80;
        let y = Bits::from_indices(n, (0..n).filter(|i| i % 2 == 0));
        let mut x = BitMatrix::with_cols(n);
        x.push_row(y.words());
        let fit = Bits::from_indices(n, 0..60);
        let path = fit_path(&x, &y, &fit, &[0], 7, 10);
        let lambdas: Vec<f64> = path[..3].iter().map(|m| m.lambda).collect();
        let refit = refit_path(&x, &y, &Bits::ones(n), &[0], &lambdas);
        assert_eq!(refit.last().unwrap().lambda, path[2].lambda);
        // Rebuilding a three-point path instead would choose its 0.001
        // endpoint, not the selected point on the original seven-point path.
        let wrong = fit_path(&x, &y, &Bits::ones(n), &[0], 3, 10);
        assert_ne!(refit.last().unwrap().lambda, wrong.last().unwrap().lambda);
        assert!(fit_path(&x, &y, &fit, &[0], 0, 10).is_empty());
    }

    #[test]
    fn selects_the_informative_feature_and_drops_noise() {
        let n = 300;
        let mut x = BitMatrix::with_cols(n);
        let marker = Bits::from_indices(n, (0..n).filter(|i| i % 3 == 0));
        x.push_row(Bits::from_indices(n, (0..n).filter(|i| i % 5 == 0)).words()); // noise
        x.push_row(marker.words());
        x.push_row(Bits::from_indices(n, (0..n).filter(|i| i % 7 == 1)).words()); // noise
        // 90% of marker carriers resistant, 5% of the rest.
        let y = Bits::from_indices(
            n,
            (0..n).filter(|i| (i % 3 == 0 && i % 10 != 0) || i % 20 == 1),
        );
        let path = fit_path(&x, &y, &Bits::ones(n), &[0, 1, 2], 20, 10);
        // The first feature to enter the path is the marker.
        let first = path.iter().find(|m| !m.coef.is_empty()).unwrap();
        assert_eq!(first.coef.iter().map(|c| c.0).collect::<Vec<_>>(), vec![1]);
        assert!(first.coef[0].1 > 0.0);
        let last = path.last().unwrap();
        let p = last.predict_proba(&x);
        assert!(p[3] > 0.5 && p[4] < 0.5);
    }

    #[test]
    fn features_constant_on_fitting_samples_get_no_weight() {
        // Every fitting sample carries f1 and no held-out sample does. Any
        // weight on f1 duplicates the intercept and shifts held-out scores.
        let n = 100;
        let y = Bits::from_indices(
            n,
            (0..n).filter(|i| (i % 3 == 0 && i % 10 != 0) || i % 20 == 1),
        );
        let mut x = BitMatrix::with_cols(n);
        x.push_row(Bits::from_indices(n, (0..n).filter(|i| i % 3 == 0)).words());
        x.push_row(Bits::from_indices(n, 0..75).words());
        let path = fit_path(&x, &y, &Bits::from_indices(n, 0..75), &[0, 1], 10, 10);
        assert!(path.iter().any(|m| !m.coef.is_empty()));
        for m in &path {
            assert!(m.coef.iter().all(|&(f, _)| f == 0), "{:?}", m.coef);
        }
    }

    #[test]
    fn a_guard_that_ends_the_solve_is_flagged() {
        let n = 300;
        let marker = Bits::from_indices(n, (0..n).filter(|i| i % 3 == 0));
        let y = Bits::from_indices(
            n,
            (0..n).filter(|i| (i % 3 == 0 && i % 10 != 0) || i % 20 == 1),
        );
        let mut x = BitMatrix::with_cols(n);
        x.push_row(marker.words());
        // One active-set iteration cannot finish a solve, and one IRLS step
        // cannot reach the tolerance from the warm start.
        let tight = Limits {
            irls_steps: 1,
            inner_iterations: 1,
        };
        let capped = fit_path_limited(&x, &y, &Bits::ones(n), &[0], 3, 10, tight);
        let last = capped.last().unwrap().convergence;
        assert!(!last.converged);
        assert_eq!((last.irls_steps, last.sweep_cap_hits), (1, 1));
        // The IRLS guard alone leaves the point unconverged too.
        let irls_only = Limits {
            irls_steps: 1,
            ..LIMITS
        };
        let capped = fit_path_limited(&x, &y, &Bits::ones(n), &[0], 3, 10, irls_only);
        let last = capped.last().unwrap().convergence;
        assert!(!last.converged);
        assert_eq!((last.irls_steps, last.sweep_cap_hits), (1, 0));
        let fits = fit_path_with_status(&x, &y, &Bits::ones(n), &[0], 3, 10, &mut |_| {});
        assert!(fits.iter().all(|f| f.convergence.converged));
        assert_eq!(fits[2].convergence.sweep_cap_hits, 0);
        let direct = fit_path_limited(&x, &y, &Bits::ones(n), &[0], 3, 10, LIMITS);
        assert_eq!(fits.len(), direct.len());
        for (f, d) in fits.iter().zip(&direct) {
            assert!(same(&f.model, &d.model));
            assert_eq!(f.convergence, d.convergence);
        }
        // The stderr-reporting wrappers return the same complete paths.
        let path = fit_path(&x, &y, &Bits::ones(n), &[0], 3, 10);
        assert_eq!(path.len(), fits.len());
        for (m, f) in path.iter().zip(&fits) {
            assert!(same(m, &f.model));
        }
        let lambdas: Vec<f64> = fits.iter().map(|f| f.model.lambda).collect();
        let refits = refit_path_with_status(&x, &y, &Bits::ones(n), &[0], &lambdas, &mut |_| {});
        let refit = refit_path(&x, &y, &Bits::ones(n), &[0], &lambdas);
        assert_eq!((refits.len(), refit.len()), (fits.len(), fits.len()));
        for ((m, r), f) in refit.iter().zip(&refits).zip(&fits) {
            assert!(same(m, &r.model) && same(m, &f.model));
            assert_eq!(r.convergence, f.convergence);
        }
    }

    #[test]
    fn public_wrapper_certifies_a_saved_panel_the_old_cap_left_unconverged() {
        let f = fixture("streptococcus_agalactiae__clindamycin_lineage_fold0");
        // The 300-sweep solver stopped at its cap in all 16 IRLS steps of
        // point 6, so revision 7 left it unconverged.
        let old = &f.points[6];
        assert_eq!((old.old_irls_steps, old.old_sweep_cap_hits), (16, 16));
        // This default test goes through the public wrapper: a guard that
        // stopped the solve early would leave a point unconverged or above
        // the tolerance.
        let feats: Vec<u32> = (0..f.x.rows as u32).collect();
        let actual = fit_path_with_status(
            &f.x,
            &f.y,
            &f.fit,
            &feats,
            f.n_lambdas,
            f.max_nonzero,
            &mut |_| {},
        );
        assert_eq!(actual.len(), f.points.len());
        for (i, (a, saved)) in actual.iter().zip(&f.points).enumerate() {
            assert!(a.convergence.converged, "point {i} {:?}", a.convergence);
            let (obj, kkt) = objective_and_kkt(&f.x, &f.y, &f.fit, &a.model);
            assert!(kkt <= KKT_TOL, "point {i}: KKT {kkt:e}");
            assert!(obj - saved.best_objective <= 1e-6, "point {i}");
        }
        for (a, r) in actual.iter().zip(f.run()) {
            assert!(same(&a.model, &r.model));
            assert_eq!(a.convergence, r.convergence);
        }
    }

    #[test]
    fn public_refit_and_model_only_wrappers_repeat_the_converged_path() {
        let f = fixture("streptococcus_agalactiae__clindamycin_lineage_fold0");
        let feats: Vec<u32> = (0..f.x.rows as u32).collect();
        let new = f.run();
        // The evaluator scores and saves the refit, so a default test pins
        // it too: the same penalties on the same mask repeat the converged
        // path.
        let lambdas: Vec<f64> = new.iter().map(|n| n.model.lambda).collect();
        let refit = refit_path_with_status(&f.x, &f.y, &f.fit, &feats, &lambdas, &mut |_| {});
        assert_eq!(refit.len(), new.len());
        for (r, n) in refit.iter().zip(new) {
            assert!(r.convergence.converged);
            assert!(same(&r.model, &n.model));
            assert_eq!(r.convergence, n.convergence);
        }
        // Nothing in production calls the model-only wrappers, so nothing
        // else would notice them drifting from the status wrappers.
        let path = fit_path(&f.x, &f.y, &f.fit, &feats, f.n_lambdas, f.max_nonzero);
        let refit = refit_path(&f.x, &f.y, &f.fit, &feats, &lambdas);
        assert_eq!((path.len(), refit.len()), (new.len(), new.len()));
        for ((p, r), n) in path.iter().zip(&refit).zip(new) {
            assert!(same(p, &n.model) && same(r, &n.model));
        }
    }

    #[test]
    fn one_class_long_paths_are_empty_by_design() {
        let n = 80;
        let mut x = BitMatrix::with_cols(n);
        x.push_row(Bits::from_indices(n, (0..n).filter(|i| i % 2 == 0)).words());
        let fit = Bits::ones(n);
        let lambdas = vec![0.1; 4001];
        for y in [Bits::zeros(n), Bits::ones(n)] {
            let path = fit_path_with_status(&x, &y, &fit, &[0], 4001, usize::MAX, &mut |_| {});
            assert!(path.is_empty());
            assert!(refit_path_with_status(&x, &y, &fit, &[0], &lambdas, &mut |_| {}).is_empty());
        }
    }

    #[test]
    fn non_finite_parameters_never_report_convergence() {
        let converged = Convergence {
            irls_steps: 1,
            sweep_cap_hits: 0,
            converged: true,
        };
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(!require_finite(converged, bad, &[0.0]).converged);
            assert!(!require_finite(converged, 0.0, &[0.0, bad]).converged);
        }
        assert!(require_finite(converged, 0.0, &[1.0]).converged);

        // f64::max drops NaN, so the residual must not.
        let prob = |class_weight| Problem {
            columns: vec![vec![0]],
            features: vec![0],
            sd: vec![1.0],
            samples: vec![0, 1],
            y: vec![1.0, 0.0],
            class_weight,
        };
        let finite = prob(vec![1.0, 1.0]);
        let rep = representatives(&finite.columns);
        assert!(kkt_residual(&finite, &rep, 0.3, &[0.0], &[0.0, f64::NAN]).is_nan());
        assert!(kkt_residual(&finite, &rep, 0.3, &[0.0], &[0.0, 0.0]) <= KKT_TOL);

        // A NaN or infinite working weight leaves no least-squares problem to
        // solve: the point stops at once, unconverged, and its residual is
        // never reported as small.
        for bad in [f64::NAN, f64::INFINITY] {
            let mut seen = Vec::new();
            let fits = solve_path(&prob(vec![bad, 1.0]), &[0.1], 10, LIMITS, &mut |p| {
                seen.push(*p)
            });
            let c = fits[0].convergence;
            assert!(!c.converged);
            assert_eq!((c.irls_steps, c.sweep_cap_hits), (1, 1));
            assert!(
                seen[0].kkt.is_nan() || seen[0].kkt > KKT_TOL,
                "{:?}",
                seen[0]
            );
        }
    }

    #[test]
    fn the_kkt_residual_covers_the_intercept_and_every_condition() {
        // Two samples, one resistant, class weights 1: at eta = 0 the
        // gradient of the column carried by the resistant sample is
        // (0.5 - 1) / 2 = -0.25, and the intercept's is 0.
        let prob = Problem {
            columns: vec![vec![0], vec![0]],
            features: vec![0, 1],
            sd: vec![1.0, 1.0],
            samples: vec![0, 1],
            y: vec![1.0, 0.0],
            class_weight: vec![1.0, 1.0],
        };
        let rep = representatives(&prob.columns);
        assert_eq!(rep, [0, 0]);
        let kkt = |beta: &[f64], eta: &[f64]| kkt_residual(&prob, &rep, 0.1, beta, eta);
        // A zero coefficient: |g| - lambda; its copy shares that gradient.
        assert!((kkt(&[0.0, 0.0], &[0.0, 0.0]) - 0.15).abs() < 1e-15);
        // A non-zero coefficient: |g + lambda sign(b)|, so its sign counts.
        // The copy at zero then violates by |g| - lambda.
        assert!((kkt(&[0.5, 0.0], &[0.0, 0.0]) - 0.15).abs() < 1e-15);
        assert!((kkt(&[-0.5, 0.0], &[0.0, 0.0]) - 0.35).abs() < 1e-15);
        // A NaN coefficient makes its own condition NaN; f64::max over the
        // conditions would drop it.
        assert!(kkt(&[f64::NAN, 0.0], &[0.0, 0.0]).is_nan());
        // The intercept: with every coefficient's condition met under a large
        // penalty, an intercept of 1 still leaves (2 sigmoid(1) - 1) / 2.
        let off = kkt_residual(&prob, &rep, 10.0, &[0.0, 0.0], &[1.0, 1.0]);
        let expected = (2.0 * crate::learn::gbt::sigmoid(1.0) - 1.0) / 2.0;
        assert!((off - expected).abs() < 1e-15, "{off} {expected}");
    }

    #[test]
    fn every_point_is_reported_as_it_is_solved() {
        let f = fixture("planted_logistic");
        let feats: Vec<u32> = (0..f.x.rows as u32).collect();
        let mut seen = Vec::new();
        let fits = fit_path_with_status(
            &f.x,
            &f.y,
            &f.fit,
            &feats,
            f.n_lambdas,
            f.max_nonzero,
            &mut |p| seen.push(*p),
        );
        assert_eq!(seen.len(), fits.len());
        for (i, (p, fit)) in seen.iter().zip(&fits).enumerate() {
            let c = fit.convergence;
            assert_eq!(
                (p.point, p.lambda.to_bits(), p.active),
                (i, fit.model.lambda.to_bits(), fit.model.coef.len())
            );
            assert_eq!((p.irls_steps, p.converged), (c.irls_steps, c.converged));
            // The solver's residual is its convergence test; it matches the
            // one computed here independently.
            let (_, kkt) = objective_and_kkt(&f.x, &f.y, &f.fit, &fit.model);
            assert!((p.kkt - kkt).abs() <= 1e-12, "point {i}: {} {kkt}", p.kkt);
        }
    }

    #[test]
    fn the_line_search_stops_where_a_coefficient_reaches_zero() {
        // One coefficient at 1 whose full step goes to -2, crossing zero at
        // t = 1/3. The smooth part alone is lowest at t = 1 (-1); with the
        // penalty, t = 1/3 gives -2/3 + 1/9 - 1 and t = 1 gives -1 + 1.
        let (t, dq) = line_search(&[0.0, 1.0], &[0.0, -3.0], -2.0, 2.0, &[1.0]);
        assert_eq!(t, 1.0 / 3.0);
        assert!((dq - (-14.0 / 9.0)).abs() < 1e-15);
        // With a small penalty the full step is lower.
        let (t, dq) = line_search(&[0.0, 1.0], &[0.0, -3.0], -2.0, 2.0, &[0.1]);
        assert_eq!(t, 1.0);
        assert!((dq - (-1.0 + 0.1)).abs() < 1e-15);
        // Of two crossings it takes the lower: coefficient 2 crosses at
        // t = 0.5, coefficient 1 at t = 0.25, and the heavy penalty on
        // coefficient 2 favors its crossing.
        let (t, _) = line_search(&[0.0, 1.0, 1.0], &[0.0, -4.0, -2.0], -1.0, 1.0, &[0.0, 5.0]);
        assert_eq!(t, 0.5);
        // A coefficient entering at zero has no crossing.
        let (t, dq) = line_search(&[0.3, 0.0], &[0.1, 2.0], -1.0, 1.0, &[0.2]);
        assert_eq!(t, 1.0);
        assert!((dq - (-1.0 + 0.5 + 0.4)).abs() < 1e-15);
    }

    #[test]
    fn a_coefficient_warm_started_with_the_wrong_sign_crosses_zero() {
        let (x, y) = union_panel();
        let prob = Problem::new(&x, &y, &Bits::ones(x.cols), &[0, 3]);
        let rep = representatives(&prob.columns);
        let lambda = 0.02;
        let eta = vec![0.0; x.cols];
        let mut wls = Wls::new(&prob, &rep, lambda, &eta).unwrap();
        // Column 0 is carried mostly by resistant samples; start it negative.
        let (mut b0, mut beta) = (0.0, vec![-1.0, 0.0]);
        assert!(wls.solve(&mut b0, &mut beta, 1000));
        assert!(beta[0] > 0.0);
        let found = wls_objective(&wls, b0, &beta);
        let best = brute_force(&wls);
        assert!(
            found - best <= 1e-12 * best.abs().max(1.0),
            "{found} {best}"
        );
    }

    #[test]
    fn a_column_dependent_on_the_active_set_is_exchanged_in() {
        // Column 2 is the union of the disjoint columns 0 and 1, so
        // x2 = x0 + x1 and it can never join a factor that holds both. Its
        // penalty factor, sqrt(p (1 - p)) at twice the frequency, is below
        // the sum of theirs, so while both carry positive weight it violates
        // its condition.
        let (x, y) = union_panel();
        let prob = Problem::new(&x, &y, &Bits::ones(x.cols), &[0, 1, 2, 3]);
        let rep = representatives(&prob.columns);
        let lambda = 0.01;
        let eta = vec![0.0; x.cols];
        let mut wls = Wls::new(&prob, &rep, lambda, &eta).unwrap();
        let (mut b0, mut beta) = (-0.5, vec![0.5, 0.5, 0.0, 0.0]);
        assert!(wls.solve(&mut b0, &mut beta, 1000));
        assert!(beta[2] > 0.0, "{beta:?}");
        assert!(beta[0] == 0.0 || beta[1] == 0.0, "{beta:?}");
        let found = wls_objective(&wls, b0, &beta);
        let best = brute_force(&wls);
        assert!(
            found - best <= 1e-12 * best.abs().max(1.0),
            "{found} {best}"
        );

        // The full path certifies every point on the same panel.
        let feats = [0, 1, 2, 3];
        let fits = fit_path_with_status(&x, &y, &Bits::ones(x.cols), &feats, 15, 400, &mut |_| {});
        assert!(fits.iter().any(|f| f.model.coef.iter().any(|c| c.0 == 2)));
        for (i, f) in fits.iter().enumerate() {
            let c = f.convergence;
            assert!(c.converged && c.sweep_cap_hits == 0, "point {i} {c:?}");
            let (_, kkt) = objective_and_kkt(&x, &y, &Bits::ones(x.cols), &f.model);
            assert!(kkt <= KKT_TOL, "point {i}: KKT {kkt:e}");
        }
    }

    #[test]
    fn exact_copies_share_one_coefficient_and_leave_the_objective_unchanged() {
        // Each column maps to the first column identical to it. The pivot
        // test would also keep a copy out of the factor, but only through
        // rounding; this keeps the weight on the first copy by rule.
        let columns = [vec![0, 1], vec![2], vec![0, 1], vec![], vec![2], vec![]];
        assert_eq!(representatives(&columns), [0, 1, 0, 3, 1, 3]);

        let n = 200;
        let fit = Bits::from_indices(n, 0..160);
        let mut rng = Lcg(11);
        let base: Vec<Bits> = (0..6).map(|_| random_bits(&mut rng, n, 30)).collect();
        let y = Bits::from_indices(
            n,
            (0..n).filter(|&i| (base[0].contains(i) || base[1].contains(i)) != (i % 9 == 0)),
        );
        // Columns 6 to 11 copy columns 0 to 5 on the fitting samples. The
        // copy of column 0 also covers held-out samples 160 to 179.
        let mut lone = BitMatrix::with_cols(n);
        let mut full = BitMatrix::with_cols(n);
        for b in &base {
            lone.push_row(b.words());
            full.push_row(b.words());
        }
        full.push_row(base[0].or(&Bits::from_indices(n, 160..180)).words());
        for b in &base[1..] {
            full.push_row(b.words());
        }
        let feats: Vec<u32> = (0..12).collect();
        let a = fit_path_with_status(&lone, &y, &fit, &feats[..6], 20, 400, &mut |_| {});
        let b = fit_path_with_status(&full, &y, &fit, &feats, 20, 400, &mut |_| {});
        assert_eq!(a.len(), b.len());
        for (i, (a, b)) in a.iter().zip(&b).enumerate() {
            assert_eq!(a.model.lambda.to_bits(), b.model.lambda.to_bits());
            assert!(b.convergence.converged, "point {i}");
            // Only the first of each pair carries weight, the same weight as
            // without the copies.
            assert_eq!(a.model.coef.len(), b.model.coef.len(), "point {i}");
            for (p, q) in a.model.coef.iter().zip(&b.model.coef) {
                assert!(p.0 == q.0 && (p.1 - q.1).abs() <= 1e-12, "point {i}");
            }
            let (oa, _) = objective_and_kkt(&lone, &y, &fit, &a.model);
            let (ob, kkt) = objective_and_kkt(&full, &y, &fit, &b.model);
            assert!((oa - ob).abs() <= 1e-12 && kkt <= KKT_TOL, "point {i}");
        }
        // max_nonzero counts distinct columns: a path capped at c stops right
        // after its first point with more than c weights, and no copy counts.
        let counts: Vec<usize> = b.iter().map(|f| f.model.coef.len()).collect();
        for cap in 0..=*counts.iter().max().unwrap() {
            let capped = fit_path_with_status(&full, &y, &fit, &feats, 20, cap, &mut |_| {});
            let stop = counts
                .iter()
                .position(|&c| c > cap)
                .map_or(counts.len(), |i| i + 1);
            assert_eq!(capped.len(), stop, "cap {cap}: {counts:?}");
        }
    }

    #[test]
    fn only_an_exactly_dependent_column_is_refused_by_the_factor() {
        // Column 2 is the union of columns 0 and 1. Column 4 is that union
        // plus one more sample, so it is independent of them, but barely.
        let (mut x, y) = union_panel();
        x.push_row(Bits::from_indices(x.cols, (0..20).chain(25..26)).words());
        let prob = Problem::new(&x, &y, &Bits::ones(x.cols), &[0, 1, 2, 3, 4]);
        let rep = representatives(&prob.columns);
        // A large linear predictor on the extra sample puts its IRLS weight
        // at the floor, the smallest weight any sample can have.
        let mut eta = vec![0.0; x.cols];
        eta[25] = 30.0;
        let mut wls = Wls::new(&prob, &rep, 0.01, &eta).unwrap();
        assert!(wls.push(0, 1.0, PIVOT_TOL).is_ok() && wls.push(1, 1.0, PIVOT_TOL).is_ok());
        let Err(dependent) = wls.push(2, 1.0, PIVOT_TOL) else {
            panic!("the union of two active columns joined the factor");
        };
        let diag: f64 = prob.columns[2].iter().map(|&j| wls.w[j as usize]).sum();
        assert!(dependent.pivot.abs() <= 1e-14 * diag, "{}", dependent.pivot);
        assert!(wls.push(4, 1.0, PIVOT_TOL).is_ok());
        assert_eq!(wls.active, [0, 1, 4]);
    }

    #[test]
    fn near_twin_columns_are_certified_in_few_irls_steps() {
        // Twelve base columns, each with three twins that differ from it in
        // one sample, as 2,820 of the 4,000 screened MTB delamanid columns
        // have a twin within one sample. Coordinate descent crawls between
        // such columns.
        let n = 240;
        let mut rng = Lcg(7);
        let base: Vec<Bits> = (0..12).map(|_| random_bits(&mut rng, n, 35)).collect();
        let mut x = BitMatrix::with_cols(n);
        for (b, col) in base.iter().enumerate() {
            x.push_row(col.words());
            for t in 0..3 {
                let mut twin = col.clone();
                let j = (b * 17 + t * 53 + 5) % n;
                if twin.contains(j) {
                    twin.remove(j);
                } else {
                    twin.insert(j);
                }
                x.push_row(twin.words());
            }
        }
        let noise = random_bits(&mut rng, n, 12);
        let y = Bits::from_indices(
            n,
            (0..n).filter(|&i| {
                let votes = (0..3).filter(|&b| base[b].contains(i)).count();
                (votes >= 2) != noise.contains(i)
            }),
        );
        let feats: Vec<u32> = (0..x.rows as u32).collect();
        let fits = fit_path_with_status(&x, &y, &Bits::ones(n), &feats, 20, 400, &mut |_| {});
        assert_eq!(fits.len(), 20);
        for (i, f) in fits.iter().enumerate() {
            let c = f.convergence;
            assert!(c.converged && c.sweep_cap_hits == 0, "point {i} {c:?}");
            assert!(c.irls_steps <= 10, "point {i} {c:?}");
            let (_, kkt) = objective_and_kkt(&x, &y, &Bits::ones(n), &f.model);
            assert!(kkt <= KKT_TOL, "point {i}: KKT {kkt:e}");
        }
        // Twins carry weight side by side at the small-penalty end, so the
        // active set holds nearly dependent columns.
        let last = &fits.last().unwrap().model;
        assert!(
            (0..12).any(|b| last.coef.iter().filter(|c| c.0 / 4 == b).count() >= 2),
            "{:?}",
            last.coef
        );
    }

    #[test]
    fn an_irls_step_that_raises_the_objective_is_halved() {
        // A marker on 18 of 20 resistant samples and 2 of 20 susceptible.
        let n = 40;
        let y = Bits::from_indices(n, 0..20);
        let mut x = BitMatrix::with_cols(n);
        x.push_row(Bits::from_indices(n, (2..20).chain(20..22)).words());
        let prob = Problem::new(&x, &y, &Bits::ones(n), &[0]);
        let lambda = 0.01;
        let start = objective(&prob, lambda, &[0.0], &[0.0; 40]);
        // A step to beta = 50, far past the optimum near 4.
        let (mut b0, mut beta, mut eta) = (0.0, vec![0.0], vec![0.0; n]);
        assert!(take_step(
            &prob,
            lambda,
            &mut b0,
            &mut beta,
            &mut eta,
            0.0,
            &[50.0]
        ));
        assert!(beta[0] > 0.0 && beta[0] < 50.0, "{beta:?}");
        assert!(objective(&prob, lambda, &beta, &eta) < start);
        let mut expected = vec![0.0; n];
        linear_predictor(&prob, b0, &beta, &mut expected);
        assert_eq!(eta, expected);
        // No halving rescues a step the wrong way, and the solution being
        // the current point is no step; neither moves anything.
        let (mut b0, mut beta, mut eta) = (0.0, vec![0.0], vec![0.0; n]);
        assert!(!take_step(
            &prob,
            lambda,
            &mut b0,
            &mut beta,
            &mut eta,
            0.0,
            &[-50.0]
        ));
        assert!(!take_step(
            &prob,
            lambda,
            &mut b0,
            &mut beta,
            &mut eta,
            0.0,
            &[0.0]
        ));
        assert_eq!((b0, beta, eta), (0.0, vec![0.0], vec![0.0; n]));
    }

    /// 60 samples. Columns 0 and 1 are disjoint (samples 0..10 and 10..20),
    /// column 2 is their union, column 3 is noise. 16 of the 20 samples in
    /// the union are resistant, and 8 of the other 40.
    fn union_panel() -> (BitMatrix, Bits) {
        let n = 60;
        let mut x = BitMatrix::with_cols(n);
        x.push_row(Bits::from_indices(n, 0..10).words());
        x.push_row(Bits::from_indices(n, 10..20).words());
        x.push_row(Bits::from_indices(n, 0..20).words());
        x.push_row(Bits::from_indices(n, (0..n).filter(|i| i % 4 == 1)).words());
        let y = Bits::from_indices(
            n,
            (0..n).filter(|&i| (i < 20 && i % 5 != 0) || (i >= 20 && i % 5 == 0)),
        );
        (x, y)
    }

    /// The least-squares objective of `wls` at (b0, beta).
    fn wls_objective(wls: &Wls, b0: f64, beta: &[f64]) -> f64 {
        let mut fit = vec![0.0; wls.w.len()];
        linear_predictor(wls.prob, b0, beta, &mut fit);
        let loss: f64 = (0..fit.len())
            .map(|j| 0.5 * wls.w[j] * (wls.z[j] - fit[j]).powi(2))
            .sum();
        loss + (0..beta.len())
            .map(|k| wls.pen(k) * beta[k].abs())
            .sum::<f64>()
    }

    /// The minimum of a small least-squares problem over every sign pattern
    /// of its coefficients. Each pattern's sign-fixed normal equations are
    /// solved by Gaussian elimination and kept only when the solution has the
    /// pattern's signs; a singular pattern is skipped, since some optimum
    /// has linearly independent non-zero columns.
    fn brute_force(wls: &Wls) -> f64 {
        let all: Vec<usize> = (0..wls.prob.columns.len()).collect();
        brute_force_over(wls, &all).0
    }

    /// `brute_force` with only the columns in `free` allowed off zero: the
    /// minimum and its intercept and coefficients.
    fn brute_force_over(wls: &Wls, free: &[usize]) -> (f64, f64, Vec<f64>) {
        let p = wls.prob.columns.len();
        let n = wls.w.len();
        let mut best = (f64::INFINITY, 0.0, vec![0.0; p]);
        for code in 0..3usize.pow(free.len() as u32) {
            let mut signs = vec![0.0; p];
            for (i, &k) in free.iter().enumerate() {
                signs[k] = [0.0, 1.0, -1.0][code / 3usize.pow(i as u32) % 3];
            }
            let on: Vec<usize> = (0..p).filter(|&k| signs[k] != 0.0).collect();
            let column = |i: usize, j: usize| -> f64 {
                if i == 0 {
                    1.0
                } else {
                    wls.prob.columns[on[i - 1]].contains(&(j as u32)) as u8 as f64
                }
            };
            let m = on.len() + 1;
            let mut a = vec![vec![0.0; m + 1]; m];
            for (r, row) in a.iter_mut().enumerate() {
                for (c, v) in row[..m].iter_mut().enumerate() {
                    *v = (0..n).map(|j| wls.w[j] * column(r, j) * column(c, j)).sum();
                }
                row[m] = (0..n)
                    .map(|j| wls.w[j] * column(r, j) * wls.z[j])
                    .sum::<f64>()
                    - if r == 0 {
                        0.0
                    } else {
                        wls.pen(on[r - 1]) * signs[on[r - 1]]
                    };
            }
            let Some(sol) = gauss(a) else { continue };
            let mut beta = vec![0.0; p];
            for (i, &k) in on.iter().enumerate() {
                beta[k] = sol[i + 1];
            }
            if on.iter().all(|&k| beta[k] * signs[k] >= 0.0) {
                let found = wls_objective(wls, sol[0], &beta);
                if found < best.0 {
                    best = (found, sol[0], beta);
                }
            }
        }
        best
    }

    fn gauss(mut a: Vec<Vec<f64>>) -> Option<Vec<f64>> {
        let m = a.len();
        for c in 0..m {
            let r = (c..m).max_by(|&i, &j| a[i][c].abs().total_cmp(&a[j][c].abs()))?;
            if a[r][c].abs() < 1e-9 {
                return None;
            }
            a.swap(c, r);
            let pivot = a[c].clone();
            for (i, row) in a.iter_mut().enumerate() {
                if i != c {
                    let f = row[c] / pivot[c];
                    for (v, p) in row[c..].iter_mut().zip(&pivot[c..]) {
                        *v -= f * p;
                    }
                }
            }
        }
        Some((0..m).map(|i| a[i][m] / a[i][i]).collect())
    }

    /// Deterministic pseudo-random draws for synthetic panels.
    struct Lcg(u64);

    impl Lcg {
        fn below(&mut self, percent: u64) -> bool {
            self.0 = self
                .0
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (self.0 >> 33) % 100 < percent
        }
    }

    fn random_bits(rng: &mut Lcg, n: usize, percent: u64) -> Bits {
        Bits::from_indices(n, (0..n).filter(|_| rng.below(percent)))
    }

    /// The saved cross-check panels with each saved candidate path; see the
    /// header of the fixture file and docs/lasso-stopping.md.
    const FIXTURES: &str = include_str!("../../tests/data/lasso-stopping-fixtures.txt");

    /// One saved point: the best objective any solver reached there, and the
    /// revision-7 solver's (300 coordinate-descent sweeps per IRLS step)
    /// trace and model.
    struct SavedPoint {
        lambda: f64,
        best_objective: f64,
        old_irls_steps: u32,
        old_converged: bool,
        old_sweep_cap_hits: u32,
        old: LassoModel,
    }

    struct Fixture {
        name: String,
        x: BitMatrix,
        y: Bits,
        fit: Bits,
        n_lambdas: usize,
        max_nonzero: usize,
        points: Vec<SavedPoint>,
        /// The path under the shipped limits, solved on first use.
        run: std::sync::OnceLock<Vec<LassoFit>>,
    }

    impl Fixture {
        /// The path under the shipped limits. It skips the public one-class
        /// guard, so `one_class_train` keeps its raw trace for the
        /// accounting.
        fn run(&self) -> &[LassoFit] {
            self.run.get_or_init(|| {
                let feats: Vec<u32> = (0..self.x.rows as u32).collect();
                fit_path_limited(
                    &self.x,
                    &self.y,
                    &self.fit,
                    &feats,
                    self.n_lambdas,
                    self.max_nonzero,
                    LIMITS,
                )
            })
        }

        fn one_class(&self) -> bool {
            let pos = self.fit.and(&self.y).count();
            pos == 0 || pos == self.fit.count()
        }
    }

    fn hex_bits(n: usize, s: &str) -> Bits {
        Bits::from_indices(
            n,
            s.bytes().enumerate().flat_map(|(i, c)| {
                let v = (c as char).to_digit(16).unwrap();
                (0..4)
                    .filter(move |b| (v >> b) & 1 == 1)
                    .map(move |b| 4 * i + b)
            }),
        )
    }

    fn fixtures() -> &'static [Fixture] {
        static PARSED: std::sync::OnceLock<Vec<Fixture>> = std::sync::OnceLock::new();
        PARSED.get_or_init(|| {
            let mut out: Vec<Fixture> = Vec::new();
            let mut n = 0;
            for line in FIXTURES.lines().filter(|l| !l.starts_with('#')) {
                let t: Vec<&str> = line.split(' ').collect();
                match t[0] {
                    "fixture" => {
                        n = t[2].parse().unwrap();
                        out.push(Fixture {
                            name: t[1].to_string(),
                            x: BitMatrix::with_cols(n),
                            y: Bits::zeros(n),
                            fit: Bits::zeros(n),
                            n_lambdas: t[4].parse().unwrap(),
                            max_nonzero: t[5].parse().unwrap(),
                            points: vec![],
                            run: std::sync::OnceLock::new(),
                        });
                    }
                    "y" => out.last_mut().unwrap().y = hex_bits(n, t[1]),
                    "fit" => out.last_mut().unwrap().fit = hex_bits(n, t[1]),
                    "r" => out
                        .last_mut()
                        .unwrap()
                        .x
                        .push_row(hex_bits(n, t[1]).words()),
                    "point" => {
                        let lambda = t[2].parse().unwrap();
                        out.last_mut().unwrap().points.push(SavedPoint {
                            lambda,
                            best_objective: t[3].parse().unwrap(),
                            old_irls_steps: t[4].parse().unwrap(),
                            old_converged: t[5] == "1",
                            old_sweep_cap_hits: t[6].parse().unwrap(),
                            old: LassoModel {
                                intercept: t[7].parse().unwrap(),
                                coef: t[8..]
                                    .iter()
                                    .map(|c| {
                                        let (f, w) = c.split_once(':').unwrap();
                                        (f.parse().unwrap(), w.parse().unwrap())
                                    })
                                    .collect(),
                                lambda,
                            },
                        });
                    }
                    other => panic!("unknown fixture line {other}"),
                }
            }
            assert_eq!(out.len(), 16);
            out
        })
    }

    fn fixture(name: &str) -> &'static Fixture {
        fixtures().iter().find(|f| f.name == name).unwrap()
    }

    /// Penalized objective and the largest stationarity violation on the
    /// standardized scale, over the intercept and the coefficients the solver
    /// may fit, computed on `fit` independently of `Problem`.
    fn objective_and_kkt(x: &BitMatrix, y: &Bits, fit: &Bits, m: &LassoModel) -> (f64, f64) {
        let samples: Vec<usize> = fit.iter_ones().collect();
        let nt = samples.len() as f64;
        let pos = samples.iter().filter(|&&i| y.contains(i)).count() as f64;
        let mut beta = vec![0.0; x.rows];
        for &(k, w) in &m.coef {
            beta[k as usize] = w;
        }
        let mut res = vec![0.0; samples.len()];
        let mut loss = 0.0;
        for (s, &i) in samples.iter().enumerate() {
            let eta = m.intercept
                + beta
                    .iter()
                    .enumerate()
                    .filter(|&(k, b)| *b != 0.0 && x.get(k, i))
                    .map(|(_, b)| b)
                    .sum::<f64>();
            let yi = y.contains(i);
            let cw = if yi {
                nt / (2.0 * pos.max(1.0))
            } else {
                nt / (2.0 * (nt - pos).max(1.0))
            };
            let z = if yi { -eta } else { eta };
            loss += cw * (z.max(0.0) + (-z.abs()).exp().ln_1p());
            res[s] = cw * (crate::learn::gbt::sigmoid(eta) - yi as u8 as f64) / nt;
        }
        let mut penalty = 0.0;
        let mut kkt = res.iter().sum::<f64>().abs();
        for (k, &b) in beta.iter().enumerate() {
            let carriers: Vec<usize> = (0..samples.len())
                .filter(|&s| x.get(k, samples[s]))
                .collect();
            let p = carriers.len() as f64 / nt;
            let sd = (p * (1.0 - p)).sqrt().max(1e-6);
            penalty += sd * b.abs();
            if carriers.is_empty() || carriers.len() == samples.len() {
                continue;
            }
            let g = carriers.iter().map(|&s| res[s]).sum::<f64>() / sd;
            let v = if b != 0.0 {
                (g + m.lambda * b.signum()).abs()
            } else {
                (g.abs() - m.lambda).max(0.0)
            };
            kkt = kkt.max(v);
        }
        (loss / nt + m.lambda * penalty, kkt)
    }

    fn same(a: &LassoModel, b: &LassoModel) -> bool {
        a.lambda.to_bits() == b.lambda.to_bits()
            && a.intercept.to_bits() == b.intercept.to_bits()
            && a.coef.len() == b.coef.len()
            && a.coef
                .iter()
                .zip(&b.coef)
                .all(|(p, q)| p.0 == q.0 && p.1.to_bits() == q.1.to_bits())
    }

    /// Fitted probabilities on the fitting samples.
    fn fitted(x: &BitMatrix, fit: &Bits, m: &LassoModel) -> Vec<f64> {
        let p = m.predict_proba(x);
        fit.iter_ones().map(|i| p[i]).collect()
    }

    #[derive(Debug, Default, PartialEq)]
    struct Tally {
        points: usize,
        /// Saved 300-sweep trace: points with a sweep-cap hit in any IRLS
        /// step.
        old_sweep_capped: usize,
        /// Saved 300-sweep trace: points that ran out of IRLS steps.
        old_irls_capped: usize,
        /// Points with an active-set solve that stopped short.
        short_solves: usize,
        /// Points flagged unconverged.
        flagged: usize,
        /// Points the 300-sweep solver converged at without a cap, where the
        /// two solutions agree.
        agree_with_saved: usize,
    }

    /// Checks every point of `f` and counts how each stopped.
    fn check(f: &Fixture, t: &mut Tally) {
        let new = f.run();
        assert_eq!(new.len(), f.points.len(), "{}", f.name);
        for (i, saved) in f.points.iter().enumerate() {
            let c = new[i].convergence;
            assert_eq!(new[i].model.lambda.to_bits(), saved.lambda.to_bits());
            if f.one_class() {
                // The public entry points return an empty one-class path;
                // the raw trace must not claim convergence either.
                assert!(!c.converged, "{} point {i}", f.name);
            } else {
                assert!(c.converged, "{} point {i} {c:?}", f.name);
                let (obj, kkt) = objective_and_kkt(&f.x, &f.y, &f.fit, &new[i].model);
                assert!(kkt <= KKT_TOL, "{} point {i}: KKT {kkt:e}", f.name);
                assert!(obj - saved.best_objective <= 1e-6, "{} point {i}", f.name);
            }
            // Where the old solver settled without a cap, both solved the
            // same problem, to a coefficient change of 1e-6 and a KKT
            // residual of 1e-6. Neither is exact, so they agree within
            // tolerance, not to the bit: at most 5.4e-10 apart in objective
            // and 7.5e-5 in a fitted probability (planted_logistic point 6).
            if saved.old_converged && saved.old_sweep_cap_hits == 0 {
                let (obj, _) = objective_and_kkt(&f.x, &f.y, &f.fit, &new[i].model);
                let (old, _) = objective_and_kkt(&f.x, &f.y, &f.fit, &saved.old);
                assert!((obj - old).abs() <= 1e-9, "{} point {i}", f.name);
                let gap = fitted(&f.x, &f.fit, &new[i].model)
                    .iter()
                    .zip(fitted(&f.x, &f.fit, &saved.old))
                    .map(|(a, b)| (a - b).abs())
                    .fold(0.0, f64::max);
                assert!(gap <= 1e-4, "{} point {i}: {gap:e}", f.name);
                t.agree_with_saved += 1;
            }
            t.points += 1;
            t.old_sweep_capped += (saved.old_sweep_cap_hits > 0) as usize;
            t.old_irls_capped += !saved.old_converged as usize;
            t.short_solves += (c.sweep_cap_hits > 0) as usize;
            t.flagged += !c.converged as usize;
        }
    }

    #[test]
    fn saved_panels_converge_and_agree_with_the_old_solver_where_it_settled() {
        let mut t = Tally::default();
        for f in fixtures() {
            check(f, &mut t);
        }
        // The 7 flagged points are one_class_train.
        let expected = Tally {
            points: 91,
            old_sweep_capped: 46,
            old_irls_capped: 15,
            short_solves: 0,
            flagged: 7,
            agree_with_saved: 38,
        };
        assert_eq!(t, expected);
    }

    #[test]
    fn campylobacter_failures_reach_the_reference_optimum() {
        let random = "campylobacter_jejuni__ciprofloxacin_random_fold0";
        let lineage = "campylobacter_jejuni__ciprofloxacin_lineage_fold0";
        // (panel, point, IRLS steps) under the shipped limits; coordinate
        // descent with 100,000 sweeps per step took 7, 9 and 18.
        for (name, i, irls_steps) in [(random, 5, 6), (random, 6, 6), (lineage, 6, 5)] {
            let f = fixture(name);
            let saved = &f.points[i];
            // The 300-sweep solver ran out of IRLS steps there, each at its
            // sweep cap, far from the optimum.
            assert_eq!(
                (
                    saved.old_irls_steps,
                    saved.old_sweep_cap_hits,
                    saved.old_converged
                ),
                (25, 25, false),
                "{name} point {i}"
            );
            let (obj, kkt) = objective_and_kkt(&f.x, &f.y, &f.fit, &saved.old);
            assert!(
                kkt > 1e-4 && obj - saved.best_objective > 1e-4,
                "{name} {i}"
            );

            let new = f.run();
            let c = new[i].convergence;
            assert!(c.converged, "{name} {i} {c:?}");
            assert_eq!(
                (c.irls_steps, c.sweep_cap_hits),
                (irls_steps, 0),
                "{name} {i}"
            );
            let (obj, kkt) = objective_and_kkt(&f.x, &f.y, &f.fit, &new[i].model);
            assert!(kkt <= KKT_TOL, "{name} point {i}: KKT {kkt:e}");
            assert!(obj - saved.best_objective <= 1e-6, "{name} point {i}");

            // Refit semantics: the same penalties on the same mask give the
            // same model and status.
            let feats: Vec<u32> = (0..f.x.rows as u32).collect();
            let lambdas: Vec<f64> = new[..=i].iter().map(|m| m.model.lambda).collect();
            let refit = refit_path_with_status(&f.x, &f.y, &f.fit, &feats, &lambdas, &mut |_| {});
            assert!(same(&refit[i].model, &new[i].model));
            assert_eq!(refit[i].convergence, c);
        }
    }

    /// Random fold 0 of M. tuberculosis delamanid, the fit that stopped the
    /// revision-8 scan, from the slowdown diagnosis
    /// (data/research/grm-lasso-r8-slowdown-20261028, gitignored, so local
    /// only): the 4,000 screened columns on all 1,248 samples with the fold's
    /// fitting, tuning-holdout and training masks (runs/x-random0), the
    /// prototype's optimum at each of the 30 tuning-path points, certified to
    /// a KKT residual of 1e-9 (runs/proto-asnewton/path-full-v3.json), and
    /// the revision-7 report.
    #[test]
    #[ignore = "MTB delamanid fold 0 from the local diagnosis folder: run with --release -- --include-ignored"]
    fn mtb_delamanid_fold0_reaches_the_certified_optima() {
        let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("data/research/grm-lasso-r8-slowdown-20261028");
        let read = |path: &str| -> serde_json::Value {
            let text = std::fs::read_to_string(dir.join(path)).unwrap();
            // Python wrote the prototype's undefined duality gaps as NaN.
            serde_json::from_str(&text.replace("NaN", "null")).unwrap()
        };
        let d = read("runs/x-random0/x-random0.json");
        let n = d["n"].as_u64().unwrap() as usize;
        let ids = |v: &serde_json::Value| -> Vec<usize> {
            v.as_array()
                .unwrap()
                .iter()
                .map(|i| i.as_u64().unwrap() as usize)
                .collect()
        };
        let mask = |key: &str| Bits::from_indices(n, ids(&d[key]));
        let (y, fit, val, train) = (mask("y"), mask("fit"), mask("val"), mask("train"));
        let mut x = BitMatrix::with_cols(n);
        for col in d["columns"].as_array().unwrap() {
            x.push_row(Bits::from_indices(n, ids(col)).words());
        }
        let feats: Vec<u32> = (0..x.rows as u32).collect();
        let certified = read("runs/proto-asnewton/path-full-v3.json");
        let certified = certified["points"].as_array().unwrap();

        let started = Instant::now();
        let path = fit_path_with_status(&x, &y, &fit, &feats, 30, 400, &mut |p| {
            eprintln!("path {p:?}")
        });
        let path_seconds = started.elapsed().as_secs_f64();
        assert_eq!((path.len(), certified.len()), (30, 30));
        let (mut most_above, mut worst_kkt) = (f64::NEG_INFINITY, 0.0_f64);
        for (i, (f, cert)) in path.iter().zip(certified).enumerate() {
            // The same problem: the diagnosis's penalty grid, to the bit.
            let lambda = cert["lambda"].as_f64().unwrap();
            assert_eq!(f.model.lambda.to_bits(), lambda.to_bits(), "point {i}");
            let c = f.convergence;
            assert!(c.converged && c.sweep_cap_hits == 0, "point {i} {c:?}");
            let (obj, kkt) = objective_and_kkt(&x, &y, &fit, &f.model);
            let above = obj - cert["obj"].as_f64().unwrap();
            eprintln!("point {i}: objective - certified {above:e}, KKT {kkt:e}");
            assert!(kkt <= KKT_TOL, "point {i}: KKT {kkt:e}");
            // At most 1e-7 above the certified optimum, and below it by no
            // more than that optimum's own residual allows. The shipped worst
            // is 8.4e-9; an INNER_TOL of 1e-6 lands at 6.1e-7 and fails here.
            assert!((-1e-9..=1e-7).contains(&above), "point {i}: {above:e}");
            most_above = most_above.max(above);
            worst_kkt = worst_kkt.max(kkt);
        }

        // Selection and refit as eval::fit_lasso_from_path does them: the
        // converged point with the lowest class-balanced log loss on the
        // tuning holdout, then the training fold refit down to its penalty.
        let loss = |f: &LassoFit| log_loss(&f.model.predict_proba(&x), &y, &val);
        let best = (0..path.len())
            .min_by(|&a, &b| loss(&path[a]).total_cmp(&loss(&path[b])))
            .unwrap();
        assert_eq!(best, 3);
        let lambdas: Vec<f64> = path[..=best].iter().map(|f| f.model.lambda).collect();
        let started = Instant::now();
        let refit = refit_path_with_status(&x, &y, &train, &feats, &lambdas, &mut |p| {
            eprintln!("refit {p:?}")
        });
        let refit_seconds = started.elapsed().as_secs_f64();
        assert!(refit.iter().all(|f| f.convergence.converged));
        let model = &refit.last().unwrap().model;
        assert_eq!(model.coef.len(), 53);

        // Calls on the fold's test samples against revision 7's refit.
        let report = read("inputs/revision7-report/mycobacterium_tuberculosis__delamanid.json");
        let r7: Vec<f64> = report["splits"]
            .as_array()
            .unwrap()
            .iter()
            .find(|s| s["scheme"] == "random")
            .unwrap()["learners"]
            .as_array()
            .unwrap()
            .iter()
            .find(|l| l["learner"] == "lasso")
            .unwrap()["scores"]
            .as_array()
            .unwrap()
            .iter()
            .map(|s| s.as_f64().unwrap())
            .collect();
        let scores = model.predict_proba(&x);
        let test = Bits::ones(n).and_not(&train);
        let changed = test
            .iter_ones()
            .filter(|&i| (scores[i] >= 0.5) != (r7[i] >= 0.5))
            .count();
        let largest = test
            .iter_ones()
            .map(|i| (scores[i] - r7[i]).abs())
            .fold(0.0, f64::max);
        eprintln!(
            "MTB delamanid random fold 0: path {path_seconds:.2} s, refit {refit_seconds:.2} s; \
             largest objective above the certified optimum {most_above:e}, largest KKT \
             {worst_kkt:e}; selected point {best}, {} refit weights; {} test samples, {changed} \
             calls changed against revision 7, largest score change {largest:e}; test log \
             loss {:.5} against revision 7's {:.5}",
            model.coef.len(),
            test.count(),
            log_loss(&scores, &y, &test),
            log_loss(&r7, &y, &test)
        );
        assert_eq!(changed, 0);
    }

    /// Under unequal IRLS weights an exactly dependent column's computed
    /// pivot can round to a small positive number, so `PIVOT_TOL`, not exact
    /// arithmetic, keeps it out of the factor. Checked for a disjoint union
    /// and for the complement of an active column.
    /// A uniform draw in [lo, hi) from a 64-bit LCG state.
    fn uniform(state: &mut u64, lo: f64, hi: f64) -> f64 {
        *state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        lo + (hi - lo) * ((*state >> 11) as f64 / (1u64 << 53) as f64)
    }

    #[test]
    fn push_refuses_exactly_dependent_columns_under_unequal_weights() {
        let n = 60;
        let (x, y) = union_panel();
        let prob = Problem::new(&x, &y, &Bits::ones(n), &[0, 1, 2, 3]);
        let rep = representatives(&prob.columns);
        let mut comp = BitMatrix::with_cols(n);
        comp.push_row(Bits::from_indices(n, 0..23).words());
        comp.push_row(Bits::from_indices(n, 23..n).words());
        let cprob = Problem::new(&comp, &y, &Bits::ones(n), &[0, 1]);
        let crep = representatives(&cprob.columns);
        let mut state = 0x2545_f491_4f6c_dd1d_u64;
        for draw in 0..2_000 {
            let spread = uniform(&mut state, 0.0, 12.0);
            let eta: Vec<f64> = (0..n)
                .map(|_| uniform(&mut state, -spread, spread))
                .collect();
            let mut wls = Wls::new(&prob, &rep, 0.01, &eta).unwrap();
            assert!(wls.push(0, 1.0, PIVOT_TOL).is_ok() && wls.push(1, 1.0, PIVOT_TOL).is_ok());
            assert!(
                wls.push(2, 1.0, PIVOT_TOL).is_err(),
                "union joined the factor, draw {draw}"
            );
            let mut wls = Wls::new(&cprob, &crep, 0.01, &eta).unwrap();
            assert!(wls.push(0, 1.0, PIVOT_TOL).is_ok());
            assert!(
                wls.push(1, 1.0, PIVOT_TOL).is_err(),
                "complement joined the factor, draw {draw}"
            );
        }
    }

    /// A tuning path saved from real data, in the format of the header of
    /// `tests/data/lasso-vibrio-fold0-fixture.txt`: only the fitting samples,
    /// all of them fitted.
    struct PathFixture {
        x: BitMatrix,
        y: Bits,
        fit: Bits,
        n_lambdas: usize,
        max_nonzero: usize,
    }

    impl PathFixture {
        fn parse(text: &str) -> PathFixture {
            let mut lines = text.lines().filter(|l| !l.starts_with('#'));
            let head: Vec<&str> = lines.next().unwrap().split(' ').collect();
            assert_eq!(head[0], "fixture");
            let n: usize = head[1].parse().unwrap();
            let mut f = PathFixture {
                x: BitMatrix::with_cols(n),
                y: Bits::zeros(n),
                fit: Bits::ones(n),
                n_lambdas: head[3].parse().unwrap(),
                max_nonzero: head[4].parse().unwrap(),
            };
            for line in lines {
                let (key, bits) = line.split_once(' ').unwrap();
                match key {
                    "y" => f.y = hex_bits(n, bits),
                    "r" => f.x.push_row(hex_bits(n, bits).words()),
                    other => panic!("unknown fixture line {other}"),
                }
            }
            assert_eq!(f.x.rows, head[2].parse::<usize>().unwrap());
            f
        }

        fn column(&self, k: usize) -> Bits {
            Bits::from_words(self.x.cols, self.x.row(k).to_vec())
        }

        fn path(&self) -> Vec<LassoFit> {
            let feats: Vec<u32> = (0..self.x.rows as u32).collect();
            fit_path_with_status(
                &self.x,
                &self.y,
                &self.fit,
                &feats,
                self.n_lambdas,
                self.max_nonzero,
                &mut |_| {},
            )
        }
    }

    /// V. cholerae nalidixic acid, random split, outer fold 0: the tuning
    /// path that stopped the revision-8 run at point 13.
    fn vibrio() -> &'static PathFixture {
        static PARSED: std::sync::OnceLock<PathFixture> = std::sync::OnceLock::new();
        PARSED.get_or_init(|| {
            PathFixture::parse(include_str!(
                "../../tests/data/lasso-vibrio-fold0-fixture.txt"
            ))
        })
    }

    #[test]
    fn vibrio_cholerae_fold0_converges_at_every_point() {
        let f = vibrio();
        let n = f.x.cols;
        assert_eq!((n, f.x.rows, f.y.count()), (57, 4000, 31));
        // Two exact complement pairs on the fitting samples: disjoint, and
        // together on every sample.
        for (a, b) in [(2209, 2985), (2389, 3174)] {
            let (a, b) = (f.column(a), f.column(b));
            assert_eq!((a.and(&b).count(), a.or(&b).count()), (0, n as u32));
        }
        let path = f.path();
        assert_eq!(path.len(), 30);
        // Revision 8 stopped at point 13 after 100 IRLS steps, 99 of them
        // with a solve at its 5,000-iteration guard, at KKT 6.1e-4: column
        // 2209 sat at 2.4e-15, its crossing blocked every step, and
        // exchanges swapped it with its complement 2985 on rounding.
        assert_eq!(path[13].model.lambda, 2.1903036636112934e-2);
        for (i, p) in path.iter().enumerate() {
            let c = p.convergence;
            assert!(c.converged && c.sweep_cap_hits == 0, "point {i} {c:?}");
            let (_, kkt) = objective_and_kkt(&f.x, &f.y, &f.fit, &p.model);
            assert!(kkt <= KKT_TOL, "point {i}: KKT {kkt:e}");
        }
        assert!(path[13].convergence.irls_steps <= 5, "{:?}", path[13]);
    }

    /// S. enterica azithromycin, random split, outer fold 3: the tuning path
    /// that stopped the revision-8 run at point 19.
    fn salmonella() -> &'static PathFixture {
        static PARSED: std::sync::OnceLock<PathFixture> = std::sync::OnceLock::new();
        PARSED.get_or_init(|| {
            use std::io::Read;
            let gz = include_bytes!("../../tests/data/lasso-salmonella-fold3-fixture.txt.gz");
            let mut text = String::new();
            flate2::read::GzDecoder::new(&gz[..])
                .read_to_string(&mut text)
                .unwrap();
            PathFixture::parse(&text)
        })
    }

    #[test]
    fn salmonella_enterica_fold3_converges_at_every_point() {
        let f = salmonella();
        assert_eq!((f.x.cols, f.x.rows, f.y.count()), (867, 4000, 118));
        // Near-copies of column 0 on the fitting samples, none exact.
        for (k, differ) in [(1, 1), (3, 2), (4, 2)] {
            let (a, b) = (f.column(0), f.column(k));
            assert_eq!(a.or(&b).count() - a.and(&b).count(), differ, "column {k}");
        }
        let path = f.path();
        assert_eq!(path.len(), 30);
        // Revision 8 stopped at point 19 after 3 IRLS steps at KKT 2.8e-3:
        // a coefficient within rounding of zero blocked every step, so 31
        // entering columns, near-copies of active ones among them, were
        // blocked in turn, the solves stopped short, and the third IRLS
        // step found no halving that lowered the objective.
        assert_eq!(path[19].model.lambda, 6.6846962608593284e-3);
        for (i, p) in path.iter().enumerate() {
            let c = p.convergence;
            assert!(c.converged && c.sweep_cap_hits == 0, "point {i} {c:?}");
            let (_, kkt) = objective_and_kkt(&f.x, &f.y, &f.fit, &p.model);
            assert!(kkt <= KKT_TOL, "point {i}: KKT {kkt:e}");
        }
    }

    #[test]
    fn the_saved_real_paths_are_the_same_at_any_thread_count() {
        for f in [vibrio(), salmonella()] {
            let run = |threads| {
                rayon::ThreadPoolBuilder::new()
                    .num_threads(threads)
                    .build()
                    .unwrap()
                    .install(|| f.path())
            };
            let (one, many, again) = (run(1), run(64), run(64));
            assert_eq!((one.len(), many.len(), again.len()), (30, 30, 30));
            for ((a, b), c) in one.iter().zip(&many).zip(&again) {
                assert!(same(&a.model, &b.model) && same(&b.model, &c.model));
                assert_eq!(
                    (a.convergence, b.convergence),
                    (c.convergence, c.convergence)
                );
            }
        }
    }

    /// How the column that starts within rounding of zero, column 0, depends
    /// on its partners.
    #[derive(Clone, Copy, Debug, PartialEq)]
    enum Dependence {
        /// Its exact complement, the V. cholerae pair.
        Complement,
        /// An exact copy, with the collapse into one coefficient turned off.
        Copy,
        /// It is the union of two disjoint partners.
        Union,
        /// Its partners differ from it in one sample and in two, the
        /// S. enterica near-copies; nothing is exactly dependent.
        NearCopy,
    }

    /// One IRLS step's least-squares problem in the state that stalled
    /// V. cholerae fold 0: a random panel whose columns 1 to m - 1 sit at
    /// their optimum scaled up by 1 + kappa, so the solve has to shrink
    /// them, column 0 at plus or minus 1e-200, and its partners at zero. The
    /// draws come from `seed`. Returns the problem, its representatives,
    /// lambda, the linear predictor and the warm start (b0, beta).
    #[allow(clippy::type_complexity)]
    fn stalled(seed: u64, dep: Dependence) -> (Problem, Vec<usize>, f64, Vec<f64>, f64, Vec<f64>) {
        let mut st = seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) ^ 0xA5A5;
        let n = 30 + uniform(&mut st, 0.0, 30.0) as usize;
        let m = 3 + uniform(&mut st, 0.0, 4.0) as usize;
        let mut columns: Vec<Vec<u32>> = (0..m)
            .map(|_| {
                (0..n as u32)
                    .filter(|_| uniform(&mut st, 0.0, 1.0) < 0.4)
                    .collect()
            })
            .collect();
        if uniform(&mut st, 0.0, 1.0) < 0.6 {
            // Column 1 is column 0 with one sample switched.
            let mut c = columns[0].clone();
            let flip = uniform(&mut st, 0.0, n as f64) as u32;
            match c.iter().position(|&v| v == flip) {
                Some(i) => {
                    c.remove(i);
                }
                None => {
                    c.push(flip);
                    c.sort();
                }
            }
            columns[1] = c;
        }
        assert!(columns.iter().all(|c| c.len() >= 2 && c.len() + 2 <= n));
        match dep {
            Dependence::Complement => {
                let c = (0..n as u32).filter(|j| !columns[0].contains(j)).collect();
                columns.push(c);
            }
            Dependence::Copy => columns.push(columns[0].clone()),
            Dependence::Union => {
                let (even, odd) = columns[0].iter().partition(|&&j| j % 2 == 0);
                columns.push(even);
                columns.push(odd);
            }
            Dependence::NearCopy => {
                let outside = (0..n as u32).find(|j| !columns[0].contains(j)).unwrap();
                let mut one = columns[0].clone();
                one.push(outside);
                one.sort();
                let two = one
                    .iter()
                    .copied()
                    .filter(|&j| j != columns[0][0])
                    .collect();
                columns.push(one);
                columns.push(two);
            }
        }
        let p = columns.len();
        let weights = [1.5, -1.0, 0.8, 0.5, -0.5, 0.3, 0.7];
        let y: Vec<f64> = (0..n)
            .map(|j| {
                let s: f64 = columns[..m]
                    .iter()
                    .zip(weights)
                    .filter(|(c, _)| c.contains(&(j as u32)))
                    .map(|(_, w)| w)
                    .sum();
                let p = crate::learn::gbt::sigmoid(s - 0.3);
                (uniform(&mut st, 0.0, 1.0) < p) as u8 as f64
            })
            .collect();
        let (nf, pos) = (n as f64, y.iter().sum::<f64>());
        assert!(pos >= 2.0 && pos <= nf - 2.0);
        let prob = Problem {
            sd: columns
                .iter()
                .map(|c| {
                    let q = c.len() as f64 / nf;
                    (q * (1.0 - q)).sqrt().max(1e-6)
                })
                .collect(),
            columns,
            features: (0..p as u32).collect(),
            samples: (0..n).collect(),
            class_weight: y
                .iter()
                .map(|&t| nf / (2.0 * if t == 1.0 { pos } else { nf - pos }))
                .collect(),
            y,
        };
        let rep: Vec<usize> = (0..p).collect();
        if dep != Dependence::Copy {
            assert_eq!(representatives(&prob.columns), rep);
        }
        let lambda = uniform(&mut st, 0.002, 0.06);
        let spread = uniform(&mut st, 0.0, 3.0);
        let eta: Vec<f64> = (0..n).map(|_| uniform(&mut st, -spread, spread)).collect();
        let wls = Wls::new(&prob, &rep, lambda, &eta).unwrap();
        let others: Vec<usize> = (1..m).collect();
        let (_, b0, optimum) = brute_force_over(&wls, &others);
        assert!(others.iter().any(|&k| optimum[k] != 0.0));
        let kappa = uniform(&mut st, 0.0, 0.3);
        let mut beta: Vec<f64> = optimum.iter().map(|b| b * (1.0 + kappa)).collect();
        beta[0] = if uniform(&mut st, 0.0, 1.0) < 0.5 {
            1e-200
        } else {
            -1e-200
        };
        (prob, rep, lambda, eta, b0, beta)
    }

    #[test]
    fn a_coefficient_within_rounding_of_zero_does_not_stall_the_solve() {
        for dep in [
            Dependence::Complement,
            Dependence::Copy,
            Dependence::Union,
            Dependence::NearCopy,
        ] {
            let (prob, rep, lambda, eta, b0, beta) = stalled(53, dep);
            // The stalled state. Column 0 crosses zero so early that no step
            // shows the objective falling, though the rest of the active set
            // still has far to go.
            let mut probe = Wls::new(&prob, &rep, lambda, &eta).unwrap();
            for (k, &b) in beta.iter().enumerate() {
                if b != 0.0 {
                    assert!(probe.push(k, b.signum(), PIVOT_TOL).is_ok());
                }
            }
            assert_eq!(probe.active[0], 0);
            let c = probe.coefficients(b0, &beta);
            let d: Vec<f64> = probe
                .sign_fixed()
                .iter()
                .zip(&c)
                .map(|(a, b)| a - b)
                .collect();
            let slope: f64 = probe
                .gram_times(&c)
                .iter()
                .zip(&probe.xwz)
                .zip(&d)
                .map(|((g, b), di)| (g - b) * di)
                .sum();
            let curv = dot(&d, &probe.gram_times(&d));
            let pen: Vec<f64> = probe.active.iter().map(|&k| probe.pen(k)).collect();
            let (_, change) = line_search(&c, &d, slope, curv, &pen);
            assert!(change >= 0.0, "{dep:?}: {change:e}");
            assert!(crossing(c[1], d[1]).is_some_and(|t| t < 1e-100), "{dep:?}");
            assert!(d[2..].iter().any(|v| v.abs() > 1e-3), "{dep:?}: {d:?}");
            if dep == Dependence::Complement {
                // Revision 8 checked the zero columns there. As in
                // V. cholerae, the complement violated its condition and
                // could not join the factor, so it was exchanged with column
                // 0 on rounding. With the other partners, each entering
                // column was blocked in turn instead.
                let violators = probe
                    .violators(b0, &beta, &vec![false; beta.len()])
                    .unwrap();
                let partner = prob.columns.len() - 1;
                assert!(violators.iter().any(|v| v.0 == partner), "{dep:?}");
                assert!(probe.push(partner, 1.0, PIVOT_TOL).is_err(), "{dep:?}");
            }

            let mut wls = Wls::new(&prob, &rep, lambda, &eta).unwrap();
            let (mut b0, mut beta) = (b0, beta);
            assert!(wls.solve(&mut b0, &mut beta, 1000), "{dep:?}");
            let found = wls_objective(&wls, b0, &beta);
            let best = brute_force(&wls);
            assert!(
                found - best <= 1e-12 * best.abs().max(1.0),
                "{dep:?}: {found} {best}"
            );
        }
    }

    #[test]
    fn an_exchange_must_lower_the_objective_by_more_than_rounding() {
        // Column 2 is the union of the disjoint columns 0 and 1. Moving the
        // weight of both onto it leaves the fit as it is and changes the
        // objective at the rate lambda (sd_2 - sd_0 - sd_1); each standard
        // deviation is set here to fix that rate.
        let (x, y) = union_panel();
        let lambda = 0.01;
        let eta = vec![0.0; x.cols];
        for (gain, taken) in [(1e-12, false), (0.05, true)] {
            let mut prob = Problem::new(&x, &y, &Bits::ones(x.cols), &[0, 1, 2, 3]);
            prob.sd[2] = prob.sd[0] + prob.sd[1] - gain;
            let rep = representatives(&prob.columns);
            let mut wls = Wls::new(&prob, &rep, lambda, &eta).unwrap();
            let (mut b0, mut beta) = (-0.5, vec![0.5, 0.5, 0.0, 0.0]);
            assert!(wls.push(0, 1.0, PIVOT_TOL).is_ok() && wls.push(1, 1.0, PIVOT_TOL).is_ok());
            let Err(dep) = wls.push(2, 1.0, PIVOT_TOL) else {
                panic!("the union joined the factor");
            };
            // Column 2's descent rate: 16 of its 20 samples are resistant.
            let mut fit = vec![0.0; x.cols];
            linear_predictor(&prob, b0, &beta, &mut fit);
            let g: f64 = prob.columns[2]
                .iter()
                .map(|&j| wls.w[j as usize] * (wls.z[j as usize] - fit[j as usize]))
                .sum();
            assert!(g > 0.0);
            assert_eq!(wls.exchange(2, g, dep, &mut b0, &mut beta), taken, "{gain}");
            if taken {
                // Both columns reach zero together, and both leave.
                assert_eq!(wls.active, [2]);
                assert_eq!((beta[0], beta[1]), (0.0, 0.0));
                assert!((beta[2] - 0.5).abs() < 1e-12, "{beta:?}");
            } else {
                assert_eq!(wls.active, [0, 1]);
                assert_eq!((b0, beta), (-0.5, vec![0.5, 0.5, 0.0, 0.0]));
            }
        }
    }

    #[test]
    fn violators_tie_to_the_lowest_column() {
        // 64 disjoint columns of 4 samples. The even ones hold 3 resistant
        // samples and the odd ones 4, so at zero every even column has the
        // same violation, bit for bit, and every odd one a larger one.
        let blocks = 64;
        let n = 4 * blocks;
        let columns: Vec<Vec<u32>> = (0..blocks as u32)
            .map(|b| (4 * b..4 * b + 4).collect())
            .collect();
        let y: Vec<f64> = (0..n)
            .map(|j| ((j / 4) % 2 == 1 || j % 4 != 3) as u8 as f64)
            .collect();
        let prob = Problem {
            sd: vec![0.5; blocks],
            features: (0..blocks as u32).collect(),
            samples: (0..n).collect(),
            class_weight: vec![1.0; n],
            columns,
            y,
        };
        let rep = representatives(&prob.columns);
        let wls = Wls::new(&prob, &rep, 0.001, &vec![0.0; n]).unwrap();
        let found: Vec<usize> = wls
            .violators(0.0, &vec![0.0; blocks], &vec![false; blocks])
            .unwrap()
            .iter()
            .map(|v| v.0)
            .collect();
        let expected: Vec<usize> = (1..blocks)
            .step_by(2)
            .chain((0..blocks).step_by(2))
            .collect();
        assert_eq!(found, expected);
    }
}
