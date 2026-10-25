//! Set Covering Machine (Marchand and Shawe-Taylor 2002), as used by Kover.
//!
//! A conjunction predicts resistant when every rule holds; a disjunction
//! predicts resistant when any rule holds. A rule is the presence or absence
//! of one feature. Rules are added greedily: each step picks the rule that
//! maximizes `covered - p * lost`, where for a conjunction `covered` counts
//! susceptible samples the rule newly rules out and `lost` counts resistant
//! samples it wrongly rules out. A disjunction is the same search with the
//! classes swapped and each chosen rule negated.
//!
//! The greedy path is prefix-nested: the model with m rules is the first m
//! rules of the model with m + 1. Cross-validation therefore scores every
//! rule count from one fit per configuration.

use super::fold_masks;
use crate::bits::{self, BitMatrix, Bits};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Kind {
    Conjunction,
    Disjunction,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Rule {
    pub feature: u32,
    pub presence: bool,
}

impl Rule {
    /// Samples on which the rule holds.
    pub fn truth(&self, x: &BitMatrix) -> Bits {
        let row = Bits::from_words(x.cols, x.row(self.feature as usize).to_vec());
        if self.presence { row } else { row.not() }
    }

    fn negate(self) -> Rule {
        Rule {
            presence: !self.presence,
            ..self
        }
    }
}

/// What one greedy step did on the training samples.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Step {
    pub covered: u32,
    pub lost: u32,
    pub utility: f64,
    /// Other rules with exactly the same covered and lost counts at this
    /// step (first few), and how many there were in total.
    pub equivalent: Vec<Rule>,
    pub n_equivalent: usize,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ScmModel {
    pub kind: Kind,
    pub p: f64,
    pub rules: Vec<Rule>,
    pub steps: Vec<Step>,
}

impl ScmModel {
    pub fn truncated(&self, n_rules: usize) -> ScmModel {
        let n = n_rules.min(self.rules.len());
        ScmModel {
            kind: self.kind,
            p: self.p,
            rules: self.rules[..n].to_vec(),
            steps: self.steps[..n].to_vec(),
        }
    }

    /// Resistant predictions over all samples.
    pub fn predict(&self, x: &BitMatrix) -> Bits {
        predict_prefix(self.kind, &self.rules, x)
    }
}

fn predict_prefix(kind: Kind, rules: &[Rule], x: &BitMatrix) -> Bits {
    match kind {
        Kind::Conjunction => rules
            .iter()
            .fold(Bits::ones(x.cols), |acc, r| acc.and(&r.truth(x))),
        Kind::Disjunction => rules
            .iter()
            .fold(Bits::zeros(x.cols), |acc, r| acc.or(&r.truth(x))),
    }
}

const MAX_EQUIVALENT: usize = 25;

/// Greedy fit on the `train` samples, up to `max_rules` rules. Listing
/// equivalent rules costs a scan over every feature per step, so
/// cross-validation fits skip it.
pub fn fit(
    x: &BitMatrix,
    y: &Bits,
    train: &Bits,
    kind: Kind,
    p: f64,
    max_rules: usize,
    track_equivalent: bool,
) -> ScmModel {
    let (mut positives, mut negatives) = match kind {
        Kind::Conjunction => (train.and(y), train.and_not(y)),
        Kind::Disjunction => (train.and_not(y), train.and(y)),
    };
    let mut rules = Vec::new();
    let mut steps = Vec::new();
    let mut counts: Vec<(u32, u32)> = vec![(0, 0); x.rows];

    while rules.len() < max_rules && negatives.count() > 0 {
        let (n_neg, n_pos) = (negatives.count(), positives.count());
        counts.par_iter_mut().enumerate().for_each(|(f, c)| {
            let row = x.row(f);
            *c = (
                bits::and_count(row, negatives.words()),
                bits::and_count(row, positives.words()),
            );
        });
        // Presence of f rules out the samples lacking f; absence rules out
        // the samples carrying it.
        let score = |f: usize, presence: bool| -> (u32, u32) {
            let (cn, cp) = counts[f];
            if presence {
                (n_neg - cn, n_pos - cp)
            } else {
                (cn, cp)
            }
        };
        let best = (0..x.rows)
            .into_par_iter()
            .flat_map_iter(|f| [(f, true), (f, false)])
            .map(|(f, presence)| {
                let (covered, lost) = score(f, presence);
                (covered as f64 - p * lost as f64, covered, lost, f, presence)
            })
            .reduce_with(|a, b| {
                // Higher utility, then fewer errors, then lower index,
                // presence before absence: a fixed order keeps fits
                // reproducible across thread counts.
                let ord =
                    a.0.total_cmp(&b.0)
                        .then(b.2.cmp(&a.2))
                        .then(b.3.cmp(&a.3))
                        .then(a.4.cmp(&b.4));
                if ord.is_ge() { a } else { b }
            });
        let Some((utility, covered, lost, f, presence)) = best else {
            break;
        };
        if covered == 0 {
            break;
        }
        let chosen = Rule {
            feature: f as u32,
            presence,
        };
        let mut equivalent = Vec::new();
        let mut n_equivalent = 0;
        if track_equivalent {
            for g in 0..x.rows {
                for pres in [true, false] {
                    if (g, pres) != (f, presence) && score(g, pres) == (covered, lost) {
                        n_equivalent += 1;
                        if equivalent.len() < MAX_EQUIVALENT {
                            equivalent.push(Rule {
                                feature: g as u32,
                                presence: pres,
                            });
                        }
                    }
                }
            }
        }
        let truth = chosen.truth(x);
        negatives = negatives.and(&truth);
        positives = positives.and(&truth);
        let (stored, equivalent) = match kind {
            Kind::Conjunction => (chosen, equivalent),
            Kind::Disjunction => (
                chosen.negate(),
                equivalent.into_iter().map(Rule::negate).collect(),
            ),
        };
        rules.push(stored);
        steps.push(Step {
            covered,
            lost,
            utility,
            equivalent,
            n_equivalent,
        });
    }
    ScmModel {
        kind,
        p,
        rules,
        steps,
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Choice {
    pub kind: Kind,
    pub p: f64,
    pub n_rules: usize,
    /// Mean balanced error across inner folds.
    pub cv_error: f64,
}

pub const P_GRID: [f64; 7] = [0.1, 0.316, 0.562, 1.0, 1.778, 3.162, 10.0];

/// Chooses kind, p, and rule count by inner cross-validation on `train`.
pub fn select(
    x: &BitMatrix,
    y: &Bits,
    train: &Bits,
    max_rules: usize,
    inner_folds: &[usize],
    k: usize,
) -> Choice {
    let masks: Vec<(Bits, Bits)> = fold_masks(inner_folds, k)
        .into_iter()
        .map(|(fit_m, val)| (fit_m.and(train), val.and(train)))
        .collect();
    let mut best: Option<Choice> = None;
    for kind in [Kind::Conjunction, Kind::Disjunction] {
        for &p in &P_GRID {
            let mut err = vec![0.0; max_rules];
            for (fit_mask, val) in &masks {
                let model = fit(x, y, fit_mask, kind, p, max_rules, false);
                let val_pos = val.and(y);
                let val_neg = val.and_not(y);
                for (m, e) in err.iter_mut().enumerate() {
                    let n = (m + 1).min(model.rules.len());
                    let pred = predict_prefix(kind, &model.rules[..n], x);
                    let fnr = val_pos.and_not(&pred).count() as f64 / val_pos.count().max(1) as f64;
                    let fpr = val_neg.and(&pred).count() as f64 / val_neg.count().max(1) as f64;
                    *e += (fnr + fpr) / 2.0 / masks.len() as f64;
                }
            }
            for (m, &e) in err.iter().enumerate() {
                // Strictly better only, so ties keep the earlier (simpler)
                // rule count and the earlier grid point.
                if best.as_ref().is_none_or(|b| e < b.cv_error - 1e-12) {
                    best = Some(Choice {
                        kind,
                        p,
                        n_rules: m + 1,
                        cv_error: e,
                    });
                }
            }
        }
    }
    best.expect("grid is not empty")
}

/// Increase in training balanced error when each rule is dropped.
pub fn rule_importance(model: &ScmModel, x: &BitMatrix, y: &Bits, train: &Bits) -> Vec<f64> {
    let error = |rules: &[Rule]| {
        let pred = predict_prefix(model.kind, rules, x);
        let pos = train.and(y);
        let neg = train.and_not(y);
        let fnr = pos.and_not(&pred).count() as f64 / pos.count().max(1) as f64;
        let fpr = neg.and(&pred).count() as f64 / neg.count().max(1) as f64;
        (fnr + fpr) / 2.0
    };
    let base = error(&model.rules);
    (0..model.rules.len())
        .map(|i| {
            let without: Vec<Rule> = model
                .rules
                .iter()
                .enumerate()
                .filter(|(j, _)| *j != i)
                .map(|(_, r)| *r)
                .collect();
            error(&without) - base
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn matrix(rows: &[&[usize]], cols: usize) -> BitMatrix {
        let mut x = BitMatrix::with_cols(cols);
        for r in rows {
            x.push_row(Bits::from_indices(cols, r.iter().copied()).words());
        }
        x
    }

    #[test]
    fn conjunction_finds_two_required_features() {
        // Resistant = samples 0..4, which carry both f0 and f1. Each feature
        // alone also appears in some susceptible samples.
        let x = matrix(&[&[0, 1, 2, 3, 4, 5], &[0, 1, 2, 3, 6, 7], &[0, 4, 6]], 10);
        let y = Bits::from_indices(10, [0, 1, 2, 3]);
        let train = Bits::ones(10);
        let m = fit(&x, &y, &train, Kind::Conjunction, 1.0, 5, true);
        let mut feats: Vec<u32> = m.rules.iter().map(|r| r.feature).collect();
        feats.sort();
        assert_eq!(feats, vec![0, 1]);
        assert!(m.rules.iter().all(|r| r.presence));
        assert_eq!(m.predict(&x), y);
    }

    #[test]
    fn disjunction_finds_alternative_mechanisms() {
        // Two different genes each explain part of the resistant set.
        let x = matrix(&[&[0, 1], &[2, 3], &[0, 5, 6]], 8);
        let y = Bits::from_indices(8, [0, 1, 2, 3]);
        let train = Bits::ones(8);
        let m = fit(&x, &y, &train, Kind::Disjunction, 1.0, 5, true);
        let feats: std::collections::BTreeSet<u32> = m.rules.iter().map(|r| r.feature).collect();
        assert_eq!(feats, [0u32, 1].into_iter().collect());
        assert!(m.rules.iter().all(|r| r.presence));
        assert_eq!(m.predict(&x), y);
    }

    #[test]
    fn records_equivalent_rules() {
        // f0 and f1 have identical presence on every sample.
        let x = matrix(&[&[0, 1], &[0, 1], &[0, 3]], 4);
        let y = Bits::from_indices(4, [0, 1]);
        let m = fit(&x, &y, &Bits::ones(4), Kind::Conjunction, 1.0, 1, true);
        assert_eq!(m.rules[0].feature, 0);
        assert_eq!(m.steps[0].n_equivalent, 1);
        assert_eq!(m.steps[0].equivalent[0].feature, 1);
    }

    #[test]
    fn fit_ignores_samples_outside_train() {
        // On samples 0..4 f0 separates perfectly; sample 5 would contradict
        // it but is held out.
        let x = matrix(&[&[0, 1, 5], &[0, 2]], 6);
        let y = Bits::from_indices(6, [0, 1]);
        let train = Bits::from_indices(6, [0, 1, 2, 3, 4]);
        let m = fit(&x, &y, &train, Kind::Conjunction, 1.0, 3, true);
        assert_eq!(
            m.rules,
            vec![Rule {
                feature: 0,
                presence: true
            }]
        );
        assert_eq!(m.steps[0].lost, 0);
    }
}
