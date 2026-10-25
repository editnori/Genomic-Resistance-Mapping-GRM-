//! L1-regularized logistic regression on screened presence features.
//!
//! IRLS outer loop, cyclic coordinate descent inner loop (glmnet-style), with
//! warm starts down a lambda path. Each feature's penalty is scaled by its
//! standard deviation on the fitting samples, which is the raw-scale form of
//! fitting standardized features: rare and common features compete on equal
//! terms. Class-balancing weights match the tree learners.

use crate::bits::{self, BitMatrix, Bits};
use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct LassoModel {
    pub intercept: f64,
    /// Non-zero coefficients as (feature, weight).
    pub coef: Vec<(u32, f64)>,
    pub lambda: f64,
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
pub fn fit_path(
    x: &BitMatrix,
    y: &Bits,
    fit: &Bits,
    features: &[u32],
    n_lambdas: usize,
    max_nonzero: usize,
) -> Vec<LassoModel> {
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
    solve_path(&prob, &lambdas, max_nonzero)
}

/// Refit the selected penalty path without recomputing its scale or spacing.
pub fn refit_path(
    x: &BitMatrix,
    y: &Bits,
    fit: &Bits,
    features: &[u32],
    lambdas: &[f64],
) -> Vec<LassoModel> {
    assert!(lambdas.iter().all(|v| v.is_finite() && *v > 0.0));
    if fit.count() == 0 {
        return vec![];
    }
    solve_path(&Problem::new(x, y, fit, features), lambdas, usize::MAX)
}

fn solve_path(prob: &Problem, lambdas: &[f64], max_nonzero: usize) -> Vec<LassoModel> {
    let n = prob.samples.len();
    let nf = n as f64;

    let mut beta = vec![0.0; prob.columns.len()];
    let mut b0 = 0.0;
    let mut eta = vec![0.0; n];
    let mut path = Vec::new();

    for &lambda in lambdas {
        for _irls in 0..25 {
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
            for _cd in 0..300 {
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
                    break;
                } else {
                    active_only = true;
                }
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
                break;
            }
        }
        let coef: Vec<(u32, f64)> = beta
            .iter()
            .enumerate()
            .filter(|(_, b)| **b != 0.0)
            .map(|(k, b)| (prob.features[k], *b))
            .collect();
        let active = coef.len();
        path.push(LassoModel {
            intercept: b0,
            coef,
            lambda,
        });
        if active > max_nonzero {
            break;
        }
    }
    path
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
}
