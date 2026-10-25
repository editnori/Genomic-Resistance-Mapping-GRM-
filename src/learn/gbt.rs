//! Gradient-boosted trees on presence/absence features.
//!
//! Second-order boosting (XGBoost-style leaf weights and split gain) with two
//! objectives:
//! - `Logistic`: resistant vs susceptible, class-balanced;
//! - `Censored`: log2 MIC regression where "<= v" and "> v" values only
//!   cost something when the prediction crosses their bound (squared hinge),
//!   so off-scale MICs inform the fit without being treated as exact.
//!
//! Splits only consider features pre-screened on the fitting samples. Gradient
//! sums over a split's present side iterate the set bits of node & feature.
//!
//! A node splits only on gain above 1e-9, never on zero gain. Every logistic
//! fit starts at margin 0 with class-balanced weights, so for an exact
//! interaction such as resistance = f0 XOR f1, each combination on equally
//! many samples, the gradients on each side of either feature sum to zero:
//! no tree lowers validation loss, early stopping keeps none, and the model
//! is the 0.5 prior. Evaluation counts a model with no tree as a null fit.

use super::sample_indices;
use crate::bits::{BitMatrix, Bits};
use crate::curate::Censor;
use rand::SeedableRng;
use rand_chacha::ChaCha8Rng;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

pub enum Objective<'a> {
    Logistic {
        y: &'a Bits,
    },
    Censored {
        value: &'a [f64],
        censor: &'a [Censor],
    },
}

#[derive(Clone, Copy, Debug, Serialize)]
pub struct GbtParams {
    pub rounds: usize,
    pub learning_rate: f64,
    pub max_depth: u32,
    pub lambda: f64,
    pub min_child_weight: f64,
    pub colsample: f64,
    pub patience: usize,
    /// Evaluation replaces this per fit with a seed from `EvalConfig::seed`
    /// and the outer fold, so a recorded config leaves it out.
    #[serde(skip)]
    pub seed: u64,
}

impl Default for GbtParams {
    fn default() -> Self {
        GbtParams {
            rounds: 400,
            learning_rate: 0.1,
            max_depth: 4,
            lambda: 1.0,
            min_child_weight: 1.0,
            colsample: 0.5,
            patience: 25,
            seed: 23,
        }
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct GNode {
    /// None for a leaf.
    pub feature: Option<u32>,
    pub absent: u32,
    pub present: u32,
    pub value: f64,
    pub gain: f64,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct GTree {
    pub nodes: Vec<GNode>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Gbt {
    pub base: f64,
    pub trees: Vec<GTree>,
    pub logistic: bool,
}

impl GTree {
    fn add_to(&self, x: &BitMatrix, out: &mut [f64]) {
        let mut stack = vec![(0usize, Bits::ones(x.cols))];
        while let Some((id, mask)) = stack.pop() {
            let node = &self.nodes[id];
            match node.feature {
                Some(f) => {
                    let row = x.row(f as usize);
                    stack.push((node.absent as usize, mask.and_not_slice(row)));
                    stack.push((node.present as usize, mask.and_slice(row)));
                }
                None => {
                    for i in mask.iter_ones() {
                        out[i] += node.value;
                    }
                }
            }
        }
    }
}

impl Gbt {
    /// Raw margins (logistic) or predicted log2 MIC (censored).
    pub fn margin(&self, x: &BitMatrix) -> Vec<f64> {
        let mut out = vec![self.base; x.cols];
        for t in &self.trees {
            t.add_to(x, &mut out);
        }
        out
    }

    /// Resistant probability (logistic objective).
    pub fn predict_proba(&self, x: &BitMatrix) -> Vec<f64> {
        self.margin(x).into_iter().map(sigmoid).collect()
    }

    pub fn importance(&self) -> std::collections::HashMap<u32, f64> {
        let mut imp = std::collections::HashMap::new();
        for t in &self.trees {
            for n in &t.nodes {
                if let Some(f) = n.feature {
                    *imp.entry(f).or_default() += n.gain;
                }
            }
        }
        imp
    }
}

pub fn sigmoid(z: f64) -> f64 {
    1.0 / (1.0 + (-z).exp())
}

struct Grad {
    g: Vec<f64>,
    h: Vec<f64>,
}

fn gradients(obj: &Objective, f: &[f64], weight: &[f64]) -> Grad {
    let n = f.len();
    let mut g = vec![0.0; n];
    let mut h = vec![0.0; n];
    for i in 0..n {
        let (gi, hi) = match obj {
            Objective::Logistic { y } => {
                let p = sigmoid(f[i]);
                (p - y.contains(i) as u8 as f64, (p * (1.0 - p)).max(1e-6))
            }
            Objective::Censored { value, censor } => {
                let v = value[i];
                match censor[i] {
                    Censor::Exact => (f[i] - v, 1.0),
                    Censor::AtMost if f[i] > v => (f[i] - v, 1.0),
                    // A "> v" MIC is at least one dilution above v.
                    Censor::Above if f[i] < v + 1.0 => (f[i] - (v + 1.0), 1.0),
                    _ => (0.0, 0.0),
                }
            }
        };
        g[i] = gi * weight[i];
        h[i] = hi * weight[i];
    }
    Grad { g, h }
}

fn loss(obj: &Objective, f: &[f64], mask: &Bits, weight: &[f64]) -> f64 {
    let mut total = 0.0;
    let mut wsum = 0.0;
    for i in mask.iter_ones() {
        let l = match obj {
            Objective::Logistic { y } => {
                let p = sigmoid(f[i]).clamp(1e-12, 1.0 - 1e-12);
                if y.contains(i) {
                    -p.ln()
                } else {
                    -(1.0 - p).ln()
                }
            }
            Objective::Censored { value, censor } => {
                let v = value[i];
                let d = match censor[i] {
                    Censor::Exact => f[i] - v,
                    Censor::AtMost => (f[i] - v).max(0.0),
                    Censor::Above => ((v + 1.0) - f[i]).max(0.0),
                };
                d * d / 2.0
            }
        };
        total += l * weight[i];
        wsum += weight[i];
    }
    total / wsum.max(1e-12)
}

/// Sum of g and h over samples in `mask & row`.
fn masked_sums(mask: &[u64], row: &[u64], grad: &Grad) -> (f64, f64) {
    let (mut gs, mut hs) = (0.0, 0.0);
    for (wi, (&m, &r)) in mask.iter().zip(row).enumerate() {
        let mut w = m & r;
        while w != 0 {
            let i = wi * 64 + w.trailing_zeros() as usize;
            gs += grad.g[i];
            hs += grad.h[i];
            w &= w - 1;
        }
    }
    (gs, hs)
}

fn grow(x: &BitMatrix, fit: &Bits, grad: &Grad, candidates: &[u32], p: &GbtParams) -> GTree {
    let mut nodes = Vec::new();
    let mut stack = vec![(0usize, fit.clone(), 0u32)];
    nodes.push(GNode {
        feature: None,
        absent: 0,
        present: 0,
        value: 0.0,
        gain: 0.0,
    });
    while let Some((id, mask, depth)) = stack.pop() {
        let (gt, ht) = mask
            .iter_ones()
            .fold((0.0, 0.0), |(a, b), i| (a + grad.g[i], b + grad.h[i]));
        nodes[id].value = -gt / (ht + p.lambda) * p.learning_rate;
        if depth >= p.max_depth || ht < 2.0 * p.min_child_weight {
            continue;
        }
        let parent_score = gt * gt / (ht + p.lambda);
        let best = candidates
            .par_iter()
            .filter_map(|&f| {
                let (gl, hl) = masked_sums(mask.words(), x.row(f as usize), grad);
                let (gr, hr) = (gt - gl, ht - hl);
                if hl < p.min_child_weight || hr < p.min_child_weight {
                    return None;
                }
                let gain =
                    0.5 * (gl * gl / (hl + p.lambda) + gr * gr / (hr + p.lambda) - parent_score);
                (gain > 1e-9).then_some((gain, f))
            })
            .reduce_with(|a, b| {
                if a.0 > b.0 || (a.0 == b.0 && a.1 < b.1) {
                    a
                } else {
                    b
                }
            });
        let Some((gain, f)) = best else { continue };
        let row = x.row(f as usize);
        let absent = nodes.len();
        nodes.push(GNode {
            feature: None,
            absent: 0,
            present: 0,
            value: 0.0,
            gain: 0.0,
        });
        let present = nodes.len();
        nodes.push(GNode {
            feature: None,
            absent: 0,
            present: 0,
            value: 0.0,
            gain: 0.0,
        });
        nodes[id].feature = Some(f);
        nodes[id].absent = absent as u32;
        nodes[id].present = present as u32;
        nodes[id].gain = gain;
        stack.push((absent, mask.and_not_slice(row), depth + 1));
        stack.push((present, mask.and_slice(row), depth + 1));
    }
    GTree { nodes }
}

/// Class-balancing weights for the logistic objective; 1 otherwise.
pub fn sample_weights(obj: &Objective, fit: &Bits) -> Vec<f64> {
    let n = fit.len();
    match obj {
        Objective::Logistic { y } => {
            let pos = fit.and(y).count() as f64;
            let neg = fit.count() as f64 - pos;
            let total = pos + neg;
            (0..n)
                .map(|i| {
                    if y.contains(i) {
                        total / (2.0 * pos.max(1.0))
                    } else {
                        total / (2.0 * neg.max(1.0))
                    }
                })
                .collect()
        }
        Objective::Censored { .. } => vec![1.0; n],
    }
}

/// Boosts on `fit`, early-stopping on `valid`. `features` are the screened
/// candidate features.
pub fn fit(
    x: &BitMatrix,
    obj: &Objective,
    fit_mask: &Bits,
    valid: &Bits,
    features: &[u32],
    params: &GbtParams,
) -> Gbt {
    let weight = sample_weights(obj, fit_mask);
    let base = match obj {
        Objective::Logistic { .. } => 0.0, // balanced weights put the prior at 1:1
        Objective::Censored { value, censor } => {
            let mut v: Vec<f64> = fit_mask
                .iter_ones()
                .filter(|&i| censor[i] == Censor::Exact)
                .map(|i| value[i])
                .collect();
            v.sort_by(f64::total_cmp);
            v.get(v.len() / 2).copied().unwrap_or(0.0)
        }
    };
    let mut f = vec![base; x.cols];
    let mut rng = ChaCha8Rng::seed_from_u64(params.seed);
    let mut trees = Vec::new();
    let (mut best_loss, mut best_len, mut since_best) = (loss(obj, &f, valid, &weight), 0, 0);
    let per_tree = ((features.len() as f64 * params.colsample).ceil() as usize).max(1);

    for _ in 0..params.rounds {
        let grad = gradients(obj, &f, &weight);
        let pick = sample_indices(features.len(), per_tree, &mut rng);
        let candidates: Vec<u32> = pick.into_iter().map(|i| features[i as usize]).collect();
        let tree = grow(x, fit_mask, &grad, &candidates, params);
        tree.add_to(x, &mut f);
        trees.push(tree);
        let l = loss(obj, &f, valid, &weight);
        if l < best_loss - 1e-9 {
            best_loss = l;
            best_len = trees.len();
            since_best = 0;
        } else {
            since_best += 1;
            if since_best >= params.patience {
                break;
            }
        }
    }
    trees.truncate(best_len);
    Gbt {
        base,
        trees,
        logistic: matches!(obj, Objective::Logistic { .. }),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn early_stopping_can_keep_the_unboosted_baseline() {
        let n = 80;
        let x = matrix(&[(0..n).filter(|i| i % 2 == 0).collect()], n);
        // The only feature has the opposite association in validation.
        let y = Bits::from_indices(n, (0..n).filter(|i| (i % 2 == 0) == (*i < 40)));
        let fit_mask = Bits::from_indices(n, 0..40);
        let valid = fit_mask.not();
        let model = fit(
            &x,
            &Objective::Logistic { y: &y },
            &fit_mask,
            &valid,
            &[0],
            &GbtParams {
                rounds: 5,
                patience: 2,
                ..Default::default()
            },
        );
        assert!(model.trees.is_empty());
        assert!(model.predict_proba(&x).iter().all(|&p| p == 0.5));
    }

    #[test]
    fn pure_xor_keeps_no_tree() {
        // Resistant iff exactly one of f0, f1 is present, each combination
        // on 20 of 80 samples. At margin 0 the class-balanced gradients on
        // each side of either feature sum to zero, so no split has gain and
        // no tree lowers validation loss.
        let f0: Vec<usize> = (40..80).collect();
        let f1: Vec<usize> = (0..80).filter(|i| i / 20 % 2 == 1).collect();
        let x = matrix(&[f0, f1], 80);
        let y = Bits::from_indices(80, (0..80).filter(|i| (*i >= 40) != (i / 20 % 2 == 1)));
        let all = Bits::ones(80);
        let model = fit(
            &x,
            &Objective::Logistic { y: &y },
            &all,
            &all,
            &[0, 1],
            &GbtParams {
                colsample: 1.0,
                ..Default::default()
            },
        );
        assert!(model.trees.is_empty());
        assert!(model.predict_proba(&x).iter().all(|&p| p == 0.5));
    }

    fn matrix(rows: &[Vec<usize>], cols: usize) -> BitMatrix {
        let mut x = BitMatrix::with_cols(cols);
        for r in rows {
            x.push_row(Bits::from_indices(cols, r.iter().copied()).words());
        }
        x
    }

    #[test]
    fn logistic_boosting_learns_an_interaction() {
        // Resistant iff f0 and f1 both present; f2 is noise.
        let n = 200;
        let f0: Vec<usize> = (0..n).filter(|i| i % 2 == 0).collect();
        let f1: Vec<usize> = (0..n).filter(|i| i % 3 == 0).collect();
        let f2: Vec<usize> = (0..n).filter(|i| i % 7 == 0).collect();
        let x = matrix(&[f0, f1, f2], n);
        let y = Bits::from_indices(n, (0..n).filter(|i| i % 6 == 0));
        let all = Bits::ones(n);
        let (fit_m, val) = super::super::holdout(&all, &y, 0.2, 3);
        let model = fit(
            &x,
            &Objective::Logistic { y: &y },
            &fit_m,
            &val,
            &[0, 1, 2],
            &GbtParams::default(),
        );
        let p = model.predict_proba(&x);
        for (i, &p) in p[..n].iter().enumerate() {
            assert_eq!(p >= 0.5, y.contains(i), "sample {i} p={}", p);
        }
    }

    #[test]
    fn censored_regression_respects_bounds() {
        // Carriers of f0 have MIC 8 (log2 3); others report "<= 0.25".
        let n = 100;
        let carriers: Vec<usize> = (0..n).filter(|i| i % 4 == 0).collect();
        let x = matrix(std::slice::from_ref(&carriers), n);
        let mut value = vec![-2.0; n];
        let mut censor = vec![Censor::AtMost; n];
        for &i in &carriers {
            value[i] = 3.0;
            censor[i] = Censor::Exact;
        }
        let all = Bits::ones(n);
        let obj = Objective::Censored {
            value: &value,
            censor: &censor,
        };
        let model = fit(
            &x,
            &obj,
            &all,
            &all,
            &[0],
            &GbtParams {
                patience: 400,
                ..Default::default()
            },
        );
        let pred = model.margin(&x);
        assert!(
            (pred[0] - 3.0).abs() < 0.2,
            "carrier prediction {}",
            pred[0]
        );
        assert!(
            pred[1] <= -2.0 + 0.2,
            "non-carrier prediction {} should sit at or under its bound",
            pred[1]
        );
    }
}
