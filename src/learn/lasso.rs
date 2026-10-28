//! L1-regularized logistic regression on screened presence features.
//!
//! IRLS outer loop, cyclic coordinate descent inner loop (glmnet-style), with
//! warm starts down a lambda path. Each feature's penalty is scaled by its
//! standard deviation on the fitting samples, which is the raw-scale form of
//! fitting standardized features: rare and common features compete on equal
//! terms. Class-balancing weights match the tree learners.
//!
//! Both loops are capped. A point that stops at a cap is not a converged fit,
//! and its `Convergence` says so; see docs/lasso-stopping.md.

use crate::bits::{self, BitMatrix, Bits};
use serde::{Deserialize, Serialize};

/// Stopping rule for one penalty point.
#[derive(Clone, Copy, Debug)]
struct Limits {
    irls_steps: usize,
    sweeps: usize,
}

/// 25 IRLS steps; 100,000 coordinate-descent sweeps per step, glmnet's
/// default `maxit`. A 300-sweep cap stopped ill-conditioned panels well short
/// of the optimum (C. jejuni ciprofloxacin, excess objective up to 2.7e-4).
const LIMITS: Limits = Limits {
    irls_steps: 25,
    sweeps: 100_000,
};

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
    /// IRLS steps taken.
    pub irls_steps: u32,
    /// IRLS steps whose coordinate descent stopped at the sweep cap unsettled.
    pub sweep_cap_hits: u32,
    /// The IRLS change fell below tolerance and that step's coordinate
    /// descent settled. False when either cap ended the solve.
    pub converged: bool,
}

#[derive(Clone, Debug)]
pub struct LassoFit {
    pub model: LassoModel,
    pub convergence: Convergence,
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
/// Points that stopped at a cap are reported on stderr; use
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
        fit_path_with_status(x, y, fit, features, n_lambdas, max_nonzero),
    )
}

/// `fit_path`, with how each point stopped.
pub fn fit_path_with_status(
    x: &BitMatrix,
    y: &Bits,
    fit: &Bits,
    features: &[u32],
    n_lambdas: usize,
    max_nonzero: usize,
) -> Vec<LassoFit> {
    if !has_both_classes(y, fit) {
        return vec![];
    }
    fit_path_limited(x, y, fit, features, n_lambdas, max_nonzero, LIMITS)
}

/// One-class logistic regression has no finite intercept optimum. Do not run a
/// long divergent path whose final status could depend on floating-point
/// overflow rather than convergence. The evaluator also requires both classes
/// in every fitting and tuning partition.
fn has_both_classes(y: &Bits, fit: &Bits) -> bool {
    let pos = fit.and(y).count();
    pos > 0 && pos < fit.count()
}

fn fit_path_limited(
    x: &BitMatrix,
    y: &Bits,
    fit: &Bits,
    features: &[u32],
    n_lambdas: usize,
    max_nonzero: usize,
    limits: Limits,
) -> Vec<LassoFit> {
    if n_lambdas == 0 || fit.count() == 0 {
        return vec![];
    }
    let prob = Problem::new(x, y, fit, features);
    let lmax = prob.lambda_max();
    if lmax == 0.0 {
        return vec![];
    }
    let ratio: f64 = 1e-3;
    let lambdas: Vec<f64> = (0..n_lambdas)
        .map(|k| lmax * ratio.powf(k as f64 / (n_lambdas - 1).max(1) as f64))
        .collect();
    solve_path(&prob, &lambdas, max_nonzero, limits)
}

/// Refit the selected penalty path without recomputing its scale or spacing.
/// Points that stopped at a cap are reported on stderr; use
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
        refit_path_with_status(x, y, fit, features, lambdas),
    )
}

/// `refit_path`, with how each point stopped.
pub fn refit_path_with_status(
    x: &BitMatrix,
    y: &Bits,
    fit: &Bits,
    features: &[u32],
    lambdas: &[f64],
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
    )
}

fn report_unconverged(caller: &str, path: Vec<LassoFit>) -> Vec<LassoModel> {
    for (i, f) in path.iter().enumerate() {
        let c = f.convergence;
        if !c.converged {
            eprintln!(
                "lasso {caller}: point {i} (lambda {:e}) stopped at a cap unconverged: \
                 {} IRLS steps (cap {}), {} sweep-cap hits (cap {} sweeps)",
                f.model.lambda, c.irls_steps, LIMITS.irls_steps, c.sweep_cap_hits, LIMITS.sweeps
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
) -> Vec<LassoFit> {
    let n = prob.samples.len();
    let nf = n as f64;

    let mut beta = vec![0.0; prob.columns.len()];
    let mut b0 = 0.0;
    let mut eta = vec![0.0; n];
    let mut path = Vec::new();

    for &lambda in lambdas {
        let mut convergence = Convergence {
            irls_steps: 0,
            sweep_cap_hits: 0,
            converged: false,
        };
        for _irls in 0..limits.irls_steps {
            convergence.irls_steps += 1;
            let mut w = vec![0.0; n];
            let mut r = vec![0.0; n];
            for j in 0..n {
                let p = super::gbt::sigmoid(eta[j]);
                let wj = (p * (1.0 - p)).max(1e-5);
                w[j] = wj * prob.class_weight[j];
                r[j] = (prob.y[j] - p) / wj; // working residual z - eta
            }
            let mut max_change: f64 = 0.0;
            // Active-set iteration: after a full sweep, sweep only non-zero
            // coefficients until they settle, then confirm with a full sweep
            // that nothing new enters.
            let mut active_only = false;
            let mut inner_settled = false;
            for _cd in 0..limits.sweeps {
                let mut sweep_change: f64 = 0.0;
                for (k, col) in prob.columns.iter().enumerate() {
                    if active_only && beta[k] == 0.0 {
                        continue;
                    }
                    let (mut num, mut den) = (0.0, 0.0);
                    for &j in col {
                        let j = j as usize;
                        num += w[j] * (r[j] + beta[k]);
                        den += w[j];
                    }
                    if den == 0.0 {
                        continue;
                    }
                    let thresh = lambda * prob.sd[k] * nf;
                    let new = soft(num, thresh) / den;
                    let delta = new - beta[k];
                    if delta != 0.0 {
                        for &j in col {
                            r[j as usize] -= delta;
                        }
                        beta[k] = new;
                        sweep_change = sweep_change.max(delta.abs() * den.sqrt());
                    }
                }
                let wsum: f64 = w.iter().sum();
                let d0 = w.iter().zip(&r).map(|(a, b)| a * b).sum::<f64>() / wsum;
                b0 += d0;
                r.iter_mut().for_each(|v| *v -= d0);
                sweep_change = sweep_change.max(d0.abs());
                max_change = max_change.max(sweep_change);
                let settled = sweep_change < 1e-6;
                if active_only {
                    active_only = !settled;
                } else if settled {
                    inner_settled = true;
                    break;
                } else {
                    active_only = true;
                }
            }
            if !inner_settled {
                convergence.sweep_cap_hits += 1;
            }
            // Recompute the linear predictor from the coefficients.
            eta.iter_mut().for_each(|e| *e = b0);
            for (k, col) in prob.columns.iter().enumerate() {
                if beta[k] != 0.0 {
                    for &j in col {
                        eta[j as usize] += beta[k];
                    }
                }
            }
            if max_change < 1e-5 {
                convergence.converged = inner_settled;
                break;
            }
        }
        let convergence = require_finite(convergence, b0, &beta);
        let coef: Vec<(u32, f64)> = beta
            .iter()
            .enumerate()
            .filter(|(_, b)| **b != 0.0)
            .map(|(k, b)| (prob.features[k], *b))
            .collect();
        let active = coef.len();
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

fn soft(z: f64, t: f64) -> f64 {
    if z > t {
        z - t
    } else if z < -t {
        z + t
    } else {
        0.0
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
    fn a_cap_that_ends_the_solve_is_flagged() {
        let n = 300;
        let marker = Bits::from_indices(n, (0..n).filter(|i| i % 3 == 0));
        let y = Bits::from_indices(
            n,
            (0..n).filter(|i| (i % 3 == 0 && i % 10 != 0) || i % 20 == 1),
        );
        let mut x = BitMatrix::with_cols(n);
        x.push_row(marker.words());
        let tight = Limits {
            irls_steps: 1,
            sweeps: 1,
        };
        let capped = fit_path_limited(&x, &y, &Bits::ones(n), &[0], 3, 10, tight);
        let last = capped.last().unwrap().convergence;
        assert!(!last.converged);
        assert_eq!((last.irls_steps, last.sweep_cap_hits), (1, 1));
        let fits = fit_path_with_status(&x, &y, &Bits::ones(n), &[0], 3, 10);
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
        let refits = refit_path_with_status(&x, &y, &Bits::ones(n), &[0], &lambdas);
        let refit = refit_path(&x, &y, &Bits::ones(n), &[0], &lambdas);
        assert_eq!((refits.len(), refit.len()), (fits.len(), fits.len()));
        for ((m, r), f) in refit.iter().zip(&refits).zip(&fits) {
            assert!(same(m, &r.model) && same(m, &f.model));
            assert_eq!(r.convergence, f.convergence);
        }
    }

    #[test]
    fn public_wrapper_uses_the_raised_cap_on_a_saved_panel() {
        let f = fixtures()
            .iter()
            .find(|f| f.name == "streptococcus_agalactiae__clindamycin_lineage_fold0")
            .unwrap();
        let feats: Vec<u32> = (0..f.x.rows as u32).collect();
        // This default (not ignored) regression goes through the public
        // wrapper: reverting LIMITS.sweeps to 300 leaves point 6 unconverged.
        let actual = fit_path_with_status(&f.x, &f.y, &f.fit, &feats, f.n_lambdas, f.max_nonzero);
        assert!(actual.iter().all(|f| f.convergence.converged));
        let (old, new) = f.runs();
        assert!(old.iter().any(|f| !f.convergence.converged));
        assert_eq!(actual.len(), new.len());
        for (a, n) in actual.iter().zip(new) {
            assert!(same(&a.model, &n.model));
            assert_eq!(a.convergence, n.convergence);
        }
    }

    #[test]
    fn public_refit_and_model_only_wrappers_use_the_raised_cap_on_a_saved_panel() {
        let f = fixtures()
            .iter()
            .find(|f| f.name == "streptococcus_agalactiae__clindamycin_lineage_fold0")
            .unwrap();
        let feats: Vec<u32> = (0..f.x.rows as u32).collect();
        let (old, new) = f.runs();
        assert!(old.iter().any(|f| !f.convergence.converged));
        // The evaluator scores and saves the refit, so its cap is pinned by a
        // default test too: the same penalties on the same mask repeat the
        // converged path, and a 300-sweep refit leaves point 6 unconverged.
        let lambdas: Vec<f64> = new.iter().map(|n| n.model.lambda).collect();
        let refit = refit_path_with_status(&f.x, &f.y, &f.fit, &feats, &lambdas);
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
            assert!(fit_path_with_status(&x, &y, &fit, &[0], 4001, usize::MAX).is_empty());
            assert!(refit_path_with_status(&x, &y, &fit, &[0], &lambdas).is_empty());
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

        // Exercise the solver itself with a NaN working weight. f64::max must
        // not let the resulting NaN intercept and coefficient look settled.
        let prob = Problem {
            columns: vec![vec![0]],
            features: vec![0],
            sd: vec![1.0],
            samples: vec![0, 1],
            y: vec![1.0, 0.0],
            class_weight: vec![f64::NAN, 1.0],
        };
        let fits = solve_path(
            &prob,
            &[0.1],
            10,
            Limits {
                irls_steps: 1,
                sweeps: 10,
            },
        );
        assert!(!fits[0].model.intercept.is_finite());
        assert!(fits[0].model.coef.iter().any(|(_, b)| !b.is_finite()));
        assert_eq!(fits[0].convergence.sweep_cap_hits, 0);
        assert!(!fits[0].convergence.converged);
    }

    /// The saved cross-check panels with each saved candidate path; see the
    /// header of the fixture file and docs/lasso-stopping.md.
    const FIXTURES: &str = include_str!("../../tests/data/lasso-stopping-fixtures.txt");

    /// Limits before the fix.
    const OLD: Limits = Limits {
        irls_steps: 25,
        sweeps: 300,
    };

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
        /// Paths under the old and the current limits, solved on first use.
        runs: std::sync::OnceLock<(Vec<LassoFit>, Vec<LassoFit>)>,
    }

    impl Fixture {
        fn runs(&self) -> (&[LassoFit], &[LassoFit]) {
            let (old, new) = self.runs.get_or_init(|| {
                let feats: Vec<u32> = (0..self.x.rows as u32).collect();
                let run = |l| {
                    fit_path_limited(
                        &self.x,
                        &self.y,
                        &self.fit,
                        &feats,
                        self.n_lambdas,
                        self.max_nonzero,
                        l,
                    )
                };
                (run(OLD), run(LIMITS))
            });
            (old, new)
        }

        /// The two C. jejuni panels take about 2.3 million sweeps under the
        /// current limits: seconds optimized, many minutes unoptimized.
        fn heavy(&self) -> bool {
            self.name.starts_with("campylobacter_jejuni")
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
                            runs: std::sync::OnceLock::new(),
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

    /// Penalized objective and the largest stationarity violation on the
    /// standardized scale, over the intercept and the coefficients the solver
    /// may fit, computed independently of `Problem`.
    fn objective_and_kkt(f: &Fixture, m: &LassoModel) -> (f64, f64) {
        let samples: Vec<usize> = f.fit.iter_ones().collect();
        let nt = samples.len() as f64;
        let pos = samples.iter().filter(|&&i| f.y.contains(i)).count() as f64;
        let mut beta = vec![0.0; f.x.rows];
        for &(k, w) in &m.coef {
            beta[k as usize] = w;
        }
        let mut res = vec![0.0; samples.len()];
        let mut loss = 0.0;
        for (s, &i) in samples.iter().enumerate() {
            let eta = m.intercept
                + (0..f.x.rows)
                    .filter(|&k| f.x.get(k, i))
                    .map(|k| beta[k])
                    .sum::<f64>();
            let yi = f.y.contains(i);
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
                .filter(|&s| f.x.get(k, samples[s]))
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

    #[derive(Debug, Default, PartialEq)]
    struct Tally {
        points: usize,
        /// Old limits: points with a sweep-cap hit in any IRLS step.
        old_sweep_capped: usize,
        /// Old limits: points that ran out of IRLS steps.
        old_irls_capped: usize,
        /// Old limits: points flagged unconverged.
        old_flagged: usize,
        /// Current limits: points that still hit the sweep cap in an
        /// intermediate IRLS step.
        new_sweep_capped: usize,
        /// Current limits: points flagged unconverged.
        new_flagged: usize,
        /// Points bit-identical to the saved candidate path.
        unchanged: usize,
    }

    /// Checks every saved point of `f` and counts how each stopped.
    fn check(f: &Fixture, t: &mut Tally) {
        let (old, new) = f.runs();
        assert_eq!(old.len(), f.points.len(), "{}", f.name);
        assert_eq!(new.len(), f.points.len(), "{}", f.name);
        let mut prefix_clean = true;
        for (i, saved) in f.points.iter().enumerate() {
            // The old limits reproduce the saved candidate path and its trace.
            let o = old[i].convergence;
            assert!(same(&old[i].model, &saved.old), "{} point {i}", f.name);
            assert_eq!(o.irls_steps, saved.old_irls_steps, "{} point {i}", f.name);
            assert_eq!(o.sweep_cap_hits, saved.old_sweep_cap_hits, "{} {i}", f.name);
            // Every cap stop is flagged; a flag needs a cap stop.
            if !saved.old_converged {
                assert!(!o.converged, "{} point {i}", f.name);
            }
            if !o.converged {
                assert!(o.sweep_cap_hits > 0 || o.irls_steps as usize == OLD.irls_steps);
            }

            let c = new[i].convergence;
            assert_eq!(new[i].model.lambda.to_bits(), saved.lambda.to_bits());
            if f.one_class() {
                // Keep the raw seven-point trace for fixture accounting. The
                // public entry points return an empty one-class path instead;
                // these saved points must not claim convergence either.
                assert!(!c.converged, "{} point {i}", f.name);
                assert_eq!(c.irls_steps as usize, LIMITS.irls_steps);
            } else {
                assert!(c.converged, "{} point {i} {c:?}", f.name);
                let (obj, kkt) = objective_and_kkt(f, &new[i].model);
                assert!(kkt <= 1e-6, "{} point {i}: KKT {kkt:e}", f.name);
                assert!(obj - saved.best_objective <= 1e-6, "{} point {i}", f.name);
            }

            // Warm starts carry each point into the next, so a point is
            // unchanged only if no earlier point on its path stopped at a cap.
            prefix_clean &= o.converged && o.sweep_cap_hits == 0;
            if prefix_clean {
                assert!(same(&new[i].model, &saved.old), "{} point {i}", f.name);
                assert_eq!(c, o, "{} point {i}", f.name);
                t.unchanged += 1;
            }
            t.points += 1;
            t.old_sweep_capped += (o.sweep_cap_hits > 0) as usize;
            t.old_irls_capped += !saved.old_converged as usize;
            t.old_flagged += !o.converged as usize;
            t.new_sweep_capped += (c.sweep_cap_hits > 0) as usize;
            t.new_flagged += !c.converged as usize;
        }
    }

    #[test]
    fn saved_panels_converge_flag_their_caps_and_keep_clean_paths() {
        let mut t = Tally::default();
        for f in fixtures().iter().filter(|f| !f.heavy()) {
            check(f, &mut t);
        }
        // The 7 flagged points under the current limits are one_class_train.
        let expected = Tally {
            points: 77,
            old_sweep_capped: 34,
            old_irls_capped: 9,
            // 9 points ran out of IRLS steps; at 7 more the last inner solve
            // stopped at 300 sweeps.
            old_flagged: 16,
            new_sweep_capped: 0,
            new_flagged: 7,
            unchanged: 34,
        };
        assert_eq!(t, expected);
    }

    #[test]
    #[ignore = "C. jejuni panels: run with --release -- --include-ignored"]
    fn all_91_saved_points() {
        let mut t = Tally::default();
        for f in fixtures() {
            check(f, &mut t);
        }
        let expected = Tally {
            points: 91,
            old_sweep_capped: 46,
            old_irls_capped: 15,
            // 15 points ran out of IRLS steps; at 11 more the IRLS change met
            // its tolerance while the last inner solve stopped at 300 sweeps.
            old_flagged: 26,
            // Four C. jejuni points still reach 100,000 sweeps in an early
            // IRLS step, but a later step settles.
            new_sweep_capped: 4,
            new_flagged: 7,
            unchanged: 36,
        };
        assert_eq!(t, expected);
    }

    #[test]
    #[ignore = "C. jejuni panels: run with --release -- --include-ignored"]
    fn campylobacter_failures_reach_the_reference_optimum() {
        let all = fixtures();
        let random = "campylobacter_jejuni__ciprofloxacin_random_fold0";
        let lineage = "campylobacter_jejuni__ciprofloxacin_lineage_fold0";
        // (panel, point, IRLS steps, sweep-cap hits) under the new limits.
        for (name, i, irls_steps, sweep_cap_hits) in
            [(random, 5, 7, 6), (random, 6, 9, 6), (lineage, 6, 18, 3)]
        {
            let f = all.iter().find(|f| f.name == name).unwrap();
            let (old, new) = f.runs();
            let saved = &f.points[i];
            assert_eq!(
                old[i].convergence,
                Convergence {
                    irls_steps: 25,
                    sweep_cap_hits: 25,
                    converged: false
                },
                "{name} point {i}"
            );
            let (obj, kkt) = objective_and_kkt(f, &old[i].model);
            assert!(
                kkt > 1e-4 && obj - saved.best_objective > 1e-4,
                "{name} {i}"
            );

            let c = new[i].convergence;
            assert!(c.converged, "{name} {i} {c:?}");
            assert_eq!(c.irls_steps, irls_steps, "{name} {i}");
            assert_eq!(c.sweep_cap_hits, sweep_cap_hits, "{name} {i}");
            let (obj, kkt) = objective_and_kkt(f, &new[i].model);
            assert!(kkt <= 1e-6, "{name} point {i}: KKT {kkt:e}");
            assert!(obj - saved.best_objective <= 1e-6, "{name} point {i}");

            // Refit semantics: the same penalties on the same mask give the
            // same model and status.
            let feats: Vec<u32> = (0..f.x.rows as u32).collect();
            let lambdas: Vec<f64> = new[..=i].iter().map(|m| m.model.lambda).collect();
            let refit = refit_path_with_status(&f.x, &f.y, &f.fit, &feats, &lambdas);
            assert!(same(&refit[i].model, &new[i].model));
            assert_eq!(refit[i].convergence, c);
        }
    }
}
