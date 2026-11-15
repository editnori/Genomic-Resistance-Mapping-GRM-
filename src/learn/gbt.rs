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
//!
//! A model with no tree is a no-call, by design, and stays one: it makes no
//! call for any sample rather than reporting the constant 0.5 prior as a
//! prediction. Abstaining is the honest output of a fit that found nothing,
//! and it matches selective calling. XGBoost always keeps at least one tree,
//! so on such data it calls where GRM abstains; that difference is kept, not
//! a defect (confirmed 5 Nov 2026, DECISIONS X1). The eval.rs tests
//! `a_boosted_model_that_keeps_no_tree_is_a_null_fit_without_a_call` and
//! `null_fits_make_no_call_while_fitted_half_scores_are_calls` pin it.
//!
//! Among candidates of equal 64-bit gain a node splits on the lowest-numbered
//! feature. `TieRule::Exact`, an experiment that is off by default, moves a
//! split only to a lower feature with the same two groups: the
//! lowest-numbered candidate that divides the node's samples as the pick
//! does. Candidates with equal gains and different groups still follow
//! rounding (`tie`).

use super::sample_indices;
use super::tie::{self, Tally, TieRule};
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
    /// Evaluation sets this per fit from `EvalConfig::tie_rule`, so a
    /// recorded config leaves it out.
    #[serde(skip)]
    pub tie: TieRule,
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
            tie: TieRule::Float,
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
    /// FNV-1a digest of every node's split, children, value and gain.
    pub fn digest(&self) -> String {
        tie::digest(self.nodes.iter().flat_map(|n| {
            [
                n.feature.map_or(u64::MAX, u64::from),
                u64::from(n.absent),
                u64::from(n.present),
                n.value.to_bits(),
                n.gain.to_bits(),
            ]
        }))
    }

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

/// One tree. When `traced`, every searched node also finds the exact rule's
/// pick and counts it in `tally`; the split follows `p.tie`.
fn grow(
    x: &BitMatrix,
    fit: &Bits,
    grad: &Grad,
    candidates: &[u32],
    p: &GbtParams,
    traced: bool,
    tally: &mut Tally,
) -> GTree {
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
        let best = if p.tie.is_float() && !traced {
            best
        } else {
            let lowest = best.map(|(_, f)| {
                tie::lowest_same_partition(&mask, |g| x.row(g as usize), candidates, f)
            });
            if traced {
                tally.searched += 1;
                tally.differ_partition += (lowest != best.map(|b| b.1)) as u64;
            }
            match (p.tie, best, lowest) {
                (TieRule::Exact, Some((_, f)), Some(g)) if g != f => {
                    let (gl, hl) = masked_sums(mask.words(), x.row(g as usize), grad);
                    let (gr, hr) = (gt - gl, ht - hl);
                    let gain = 0.5
                        * (gl * gl / (hl + p.lambda) + gr * gr / (hr + p.lambda) - parent_score);
                    Some((gain, g))
                }
                _ => best,
            }
        };
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
    let traced = tie::tracing();
    let mut tally = Tally::default();
    let mut grown = 0;

    for _ in 0..params.rounds {
        let grad = gradients(obj, &f, &weight);
        let pick = sample_indices(features.len(), per_tree, &mut rng);
        let candidates: Vec<u32> = pick.into_iter().map(|i| features[i as usize]).collect();
        let tree = grow(x, fit_mask, &grad, &candidates, params, traced, &mut tally);
        grown += 1;
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
    if traced {
        tie::record(&tie::FitRecord {
            learner: "gbt",
            objective: Some(match obj {
                Objective::Logistic { .. } => "logistic",
                Objective::Censored { .. } => "censored",
            }),
            rule: params.tie,
            seed: params.seed,
            fit_samples: fit_mask.count() as usize,
            max_depth: params.max_depth,
            rounds: Some(grown),
            trees: trees.iter().map(GTree::digest).collect(),
            tally,
        });
    }
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

    /// Synthetic gradients (0.1, 0.1, 0.1, -0.1) and unit hessians on four
    /// fitting samples. Feature 0 is on samples 0 and 1, feature 1 on samples
    /// 2 and 3: the same two groups. Their 64-bit gains are
    /// 0.002666666666666666 and 0.002666666666666667 (`scripts/find_gbt_ties.py`
    /// of lane tierule). Feature 2, on sample 0 alone, divides them otherwise
    /// and has no gain. Sample 4 carries features 0 and 1 and is not fitted.
    fn mirrored_gradients() -> (BitMatrix, Grad) {
        let x = matrix(&[vec![0, 1, 4], vec![2, 3, 4], vec![0]], 5);
        let grad = Grad {
            g: vec![0.1, 0.1, 0.1, -0.1, 0.0],
            h: vec![1.0, 1.0, 1.0, 1.0, 0.0],
        };
        (x, grad)
    }

    fn stump(candidates: &[u32], tie: TieRule, traced: bool) -> (GTree, Tally) {
        let (x, grad) = mirrored_gradients();
        let params = GbtParams {
            max_depth: 1,
            tie,
            ..Default::default()
        };
        let mut tally = Tally::default();
        let fit = Bits::from_indices(5, 0..4);
        let tree = grow(&x, &fit, &grad, candidates, &params, traced, &mut tally);
        (tree, tally)
    }

    #[test]
    fn rounding_puts_the_mirrored_feature_first_unless_the_rule_is_exact() {
        let (float, _) = stump(&[1, 0], TieRule::Float, false);
        let (exact, _) = stump(&[1, 0], TieRule::Exact, false);
        let root = |t: &GTree| (t.nodes[0].feature, t.nodes[0].gain);
        assert_eq!(root(&float), (Some(1), 0.002666666666666667));
        assert_eq!(root(&exact), (Some(0), 0.002666666666666666));
        // The fitted samples get the same leaves either way; sample 4, which
        // carries both features, falls on the other side.
        let (x, _) = mirrored_gradients();
        let (mut a, mut b) = (vec![0.0; 5], vec![0.0; 5]);
        float.add_to(&x, &mut a);
        exact.add_to(&x, &mut b);
        assert_eq!(a[..4], b[..4]);
        assert_ne!(a[4], b[4]);
    }

    #[test]
    fn the_exact_rule_moves_a_split_only_to_the_same_two_groups() {
        // Without feature 0 among the tree's candidates the pick stays
        // feature 1: feature 2 divides the samples otherwise.
        let (exact, _) = stump(&[2, 1], TieRule::Exact, false);
        assert_eq!(exact.nodes[0].feature, Some(1));
        // With feature 0 drawn as well, the split moves to it.
        let (exact, _) = stump(&[2, 1, 0], TieRule::Exact, false);
        assert_eq!(exact.nodes[0].feature, Some(0));
    }

    /// The exact rule moves a boosted split only to a lower feature with the
    /// same two groups; two features with equal counts but other samples are
    /// left to the 64-bit gains, under both rules. Synthetic, from the review
    /// of lane revtie: in the first round of a logistic fit the gradients
    /// are equal within a class, so two features present on 3 of 7 resistant
    /// and 5 of 13 susceptible fit samples have equal gains in exact
    /// arithmetic. Summed over other samples their 64-bit gains differ in the
    /// last digits, and both rules split on the larger, feature 1.
    #[test]
    fn equal_counts_on_other_samples_follow_the_rounded_gain_under_both_rules() {
        use rand::Rng;
        let (np, nn) = (7usize, 13usize);
        let n = np + nn;
        // Resistant samples interleaved with susceptible ones, so subsets
        // with equal counts sum their gradients in different orders.
        let res_ix: Vec<usize> = (0..np).map(|i| 3 * i + 1).collect();
        let sus_ix: Vec<usize> = (0..n).filter(|i| !res_ix.contains(i)).collect();
        let y = Bits::from_indices(n, res_ix.iter().copied());
        let fit = Bits::ones(n);
        let obj = Objective::Logistic { y: &y };
        let w = sample_weights(&obj, &fit);
        let grad = gradients(&obj, &vec![0.0; n], &w);
        let p = GbtParams {
            max_depth: 1,
            ..Default::default()
        };
        let (gt, ht) = (0..n).fold((0.0, 0.0), |(a, b), i| (a + grad.g[i], b + grad.h[i]));
        let parent = gt * gt / (ht + p.lambda);
        let gain = |rows: &[usize]| {
            let b = Bits::from_indices(n, rows.iter().copied());
            let (gl, hl) = masked_sums(fit.words(), b.words(), &grad);
            let (gr, hr) = (gt - gl, ht - hl);
            0.5 * (gl * gl / (hl + p.lambda) + gr * gr / (hr + p.lambda) - parent)
        };
        let mut rng = ChaCha8Rng::seed_from_u64(5);
        let pick = |rng: &mut ChaCha8Rng| -> Vec<usize> {
            let (mut res, mut sus) = (res_ix.clone(), sus_ix.clone());
            let mut r: Vec<usize> = (0..3)
                .map(|_| res.swap_remove(rng.random_range(0..res.len())))
                .collect();
            r.extend((0..5).map(|_| sus.swap_remove(rng.random_range(0..sus.len()))));
            r.sort();
            r
        };
        let (low, high) = (0..10_000)
            .find_map(|_| {
                let (a, b) = (pick(&mut rng), pick(&mut rng));
                let (ga, gb) = (gain(&a), gain(&b));
                (ga != gb && ga > 1e-9 && gb > 1e-9).then_some(if ga < gb {
                    (a, b)
                } else {
                    (b, a)
                })
            })
            .expect("a pair with unequal 64-bit gains");
        assert_eq!(low, [3, 4, 8, 11, 12, 14, 16, 19]);
        assert_eq!(high, [0, 4, 5, 7, 9, 12, 15, 16]);
        assert_eq!(
            (gain(&low), gain(&high)),
            (0.014051146172066311, 0.014051146172066318)
        );
        let x = matrix(&[low, high], n);
        for tie in [TieRule::Float, TieRule::Exact] {
            let params = GbtParams {
                max_depth: 1,
                tie,
                ..Default::default()
            };
            let t = grow(
                &x,
                &fit,
                &grad,
                &[0, 1],
                &params,
                false,
                &mut Tally::default(),
            );
            assert_eq!(
                (t.nodes[0].feature, t.nodes[0].gain),
                (Some(1), 0.014051146172066318),
                "{tie:?}"
            );
        }
    }

    /// The trace's tree digest moves with each field of a node, the gain
    /// among them, so a traced fit whose gains differ is seen to differ.
    #[test]
    fn a_boosted_tree_digest_changes_with_any_node_field() {
        let (tree, _) = stump(&[1, 0, 2], TieRule::Float, false);
        let base = tree.digest();
        let edits: [fn(&mut GTree); 5] = [
            |t| t.nodes[0].feature = Some(2),
            |t| t.nodes[0].absent += 1,
            |t| t.nodes[0].present += 1,
            |t| t.nodes[1].value += 0.25,
            |t| t.nodes[0].gain *= 2.0,
        ];
        for (i, edit) in edits.iter().enumerate() {
            let mut u = tree.clone();
            edit(&mut u);
            assert_ne!(u.digest(), base, "edit {i}");
        }
    }

    #[test]
    fn a_traced_tree_is_the_untraced_tree_and_counts_the_moved_split() {
        for tie in [TieRule::Float, TieRule::Exact] {
            let (plain, _) = stump(&[1, 0, 2], tie, false);
            let (traced, tally) = stump(&[1, 0, 2], tie, true);
            assert_eq!(plain.digest(), traced.digest(), "{tie:?}");
            assert_eq!((tally.searched, tally.differ_partition), (1, 1), "{tie:?}");
        }
    }
}
