//! Learners and the evaluation that compares them.
//!
//! Every learner reads the same inputs: a feature x sample bit matrix, a
//! label bitset (resistant = 1), and a training mask. Nothing a learner does
//! may look at samples outside its training mask; that includes feature
//! screening, which is why `screen` takes the mask.

pub mod eval;
pub mod gbt;
pub mod lasso;
pub mod scm;
pub mod tree;

use crate::bits::{self, BitMatrix, Bits};
use rand::seq::SliceRandom;
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

/// Binary-classification agreement, in the clinical vocabulary.
/// VME (very major error) = resistant called susceptible, over all resistant.
/// ME (major error) = susceptible called resistant, over all susceptible.
#[derive(Clone, Debug, Default, Serialize, Deserialize, PartialEq)]
pub struct BinaryMetrics {
    pub n: usize,
    pub tp: usize,
    pub fp: usize,
    pub tn: usize,
    #[serde(rename = "fn")]
    pub fn_: usize,
    pub accuracy: f64,
    pub balanced_accuracy: f64,
    pub sensitivity: f64,
    pub specificity: f64,
    pub precision: f64,
    pub f1: f64,
    pub mcc: f64,
    pub vme: f64,
    pub me: f64,
    pub auroc: Option<f64>,
    pub auprc: Option<f64>,
}

/// Samples that received no call, by reported class. Metrics beside it
/// cover only the called samples, so called plus uncalled is every sample.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct NoCall {
    pub resistant: usize,
    pub susceptible: usize,
}

fn ratio(a: usize, b: usize) -> f64 {
    if b == 0 { 0.0 } else { a as f64 / b as f64 }
}

/// Metrics for `scores` >= `threshold` predicted resistant. Scores that
/// carry ranking information also get AUROC and AUPRC.
pub fn binary_metrics(
    truth: &[bool],
    scores: &[f64],
    threshold: f64,
    ranked: bool,
) -> BinaryMetrics {
    assert_eq!(
        truth.len(),
        scores.len(),
        "one score is required per sample"
    );
    assert!(
        scores.iter().all(|s| s.is_finite()),
        "scores must be finite"
    );
    assert!(threshold.is_finite(), "threshold must be finite");
    let (mut tp, mut fp, mut tn, mut fn_) = (0, 0, 0, 0);
    for (&t, &s) in truth.iter().zip(scores) {
        match (t, s >= threshold) {
            (true, true) => tp += 1,
            (false, true) => fp += 1,
            (false, false) => tn += 1,
            (true, false) => fn_ += 1,
        }
    }
    let sensitivity = ratio(tp, tp + fn_);
    let specificity = ratio(tn, tn + fp);
    let precision = ratio(tp, tp + fp);
    let denom =
        ((tp + fp) as f64 * (tp + fn_) as f64 * (tn + fp) as f64 * (tn + fn_) as f64).sqrt();
    BinaryMetrics {
        n: truth.len(),
        tp,
        fp,
        tn,
        fn_,
        accuracy: ratio(tp + tn, truth.len()),
        balanced_accuracy: (sensitivity + specificity) / 2.0,
        sensitivity,
        specificity,
        precision,
        f1: if precision + sensitivity == 0.0 {
            0.0
        } else {
            2.0 * precision * sensitivity / (precision + sensitivity)
        },
        mcc: if denom == 0.0 {
            0.0
        } else {
            (tp as f64 * tn as f64 - fp as f64 * fn_ as f64) / denom
        },
        vme: ratio(fn_, tp + fn_),
        me: ratio(fp, fp + tn),
        auroc: ranked.then(|| auroc(truth, scores)).flatten(),
        auprc: ranked.then(|| auprc(truth, scores)).flatten(),
    }
}

/// Area under the ROC curve via the rank-sum statistic (ties get half credit).
pub fn auroc(truth: &[bool], scores: &[f64]) -> Option<f64> {
    if truth.len() != scores.len() || scores.iter().any(|s| !s.is_finite()) {
        return None;
    }
    let mut idx: Vec<usize> = (0..scores.len()).collect();
    idx.sort_by(|&a, &b| scores[a].total_cmp(&scores[b]));
    let n_pos = truth.iter().filter(|&&t| t).count();
    let n_neg = truth.len() - n_pos;
    if n_pos == 0 || n_neg == 0 {
        return None;
    }
    let mut rank_sum = 0.0;
    let mut i = 0;
    while i < idx.len() {
        let mut j = i;
        while j + 1 < idx.len() && scores[idx[j + 1]] == scores[idx[i]] {
            j += 1;
        }
        let avg_rank = (i + j) as f64 / 2.0 + 1.0;
        rank_sum += idx[i..=j].iter().filter(|&&k| truth[k]).count() as f64 * avg_rank;
        i = j + 1;
    }
    Some((rank_sum - (n_pos * (n_pos + 1)) as f64 / 2.0) / (n_pos as f64 * n_neg as f64))
}

/// Average precision (area under the precision-recall step curve).
pub fn auprc(truth: &[bool], scores: &[f64]) -> Option<f64> {
    if truth.len() != scores.len() || scores.iter().any(|s| !s.is_finite()) {
        return None;
    }
    let n_pos = truth.iter().filter(|&&t| t).count();
    if n_pos == 0 {
        return None;
    }
    let mut idx: Vec<usize> = (0..scores.len()).collect();
    idx.sort_by(|&a, &b| scores[b].total_cmp(&scores[a]));
    let (mut tp, mut seen, mut ap) = (0usize, 0usize, 0.0);
    let mut i = 0;
    while i < idx.len() {
        let mut j = i;
        while j + 1 < idx.len() && scores[idx[j + 1]] == scores[idx[i]] {
            j += 1;
        }
        let block_tp = idx[i..=j].iter().filter(|&&k| truth[k]).count();
        seen += j - i + 1;
        if block_tp > 0 {
            tp += block_tp;
            ap += block_tp as f64 / n_pos as f64 * (tp as f64 / seen as f64);
        }
        i = j + 1;
    }
    Some(ap)
}

/// ROC curve points (fpr, tpr), thinned to at most `max_points`.
pub fn roc_curve(truth: &[bool], scores: &[f64], max_points: usize) -> Vec<(f64, f64)> {
    if max_points < 2
        || truth.len() != scores.len()
        || scores.iter().any(|s| !s.is_finite())
        || !truth.iter().any(|&v| v)
        || !truth.iter().any(|&v| !v)
    {
        return vec![];
    }
    let n_pos = truth.iter().filter(|&&t| t).count().max(1) as f64;
    let n_neg = truth.iter().filter(|&&t| !t).count().max(1) as f64;
    let mut idx: Vec<usize> = (0..scores.len()).collect();
    idx.sort_by(|&a, &b| scores[b].total_cmp(&scores[a]));
    let mut pts = vec![(0.0, 0.0)];
    let (mut tp, mut fp) = (0.0, 0.0);
    let mut i = 0;
    while i < idx.len() {
        let mut j = i;
        while j + 1 < idx.len() && scores[idx[j + 1]] == scores[idx[i]] {
            j += 1;
        }
        for &k in &idx[i..=j] {
            if truth[k] { tp += 1.0 } else { fp += 1.0 }
        }
        pts.push((fp / n_neg, tp / n_pos));
        i = j + 1;
    }
    if pts.len() > max_points {
        pts = (0..max_points)
            .map(|i| pts[i * (pts.len() - 1) / (max_points - 1)])
            .collect();
    }
    pts
}

/// Assigns each sample to one of `k` folds, stratified by class.
pub fn stratified_folds(y: &Bits, k: usize, seed: u64) -> Vec<usize> {
    assert!(k >= 2, "cross-validation needs at least two folds");
    let mut rng = ChaCha8Rng::seed_from_u64(seed);
    let mut fold = vec![0; y.len()];
    for class in [true, false] {
        let mut members: Vec<usize> = (0..y.len()).filter(|&i| y.contains(i) == class).collect();
        members.shuffle(&mut rng);
        for (j, i) in members.into_iter().enumerate() {
            fold[i] = j % k;
        }
    }
    fold
}

/// Assigns whole groups (lineage clusters) to folds so no group spans
/// train and test. Groups are placed largest first into the fold that is
/// minimizes the increase in squared class imbalance across all folds.
pub fn grouped_folds(y: &Bits, groups: &[u32], k: usize, seed: u64) -> Vec<usize> {
    assert!(k >= 2, "cross-validation needs at least two folds");
    assert_eq!(y.len(), groups.len());
    let mut members: std::collections::BTreeMap<u32, Vec<usize>> = Default::default();
    for (i, &g) in groups.iter().enumerate() {
        members.entry(g).or_default().push(i);
    }
    let mut order: Vec<(u32, Vec<usize>)> = members.into_iter().collect();
    let mut rng = ChaCha8Rng::seed_from_u64(seed);
    order.shuffle(&mut rng);
    order.sort_by_key(|(_, m)| std::cmp::Reverse(m.len()));
    let total_r = y.count() as f64;
    let total_s = y.len() as f64 - total_r;
    let mut load = vec![(0.0f64, 0.0f64); k];
    let mut fold = vec![0; y.len()];
    for (_, m) in order {
        let r = m.iter().filter(|&&i| y.contains(i)).count() as f64;
        let s = m.len() as f64 - r;
        let best = (0..k)
            .min_by(|&a, &b| {
                let cost = |f: usize| {
                    let (lr, ls) = load[f];
                    // Compare the change in the global objective, not this
                    // fold's absolute load. An existing large R-only group
                    // must not deter adding S-only groups to the same fold.
                    ((lr + r).powi(2) - lr.powi(2)) / total_r.max(1.0).powi(2)
                        + ((ls + s).powi(2) - ls.powi(2)) / total_s.max(1.0).powi(2)
                };
                cost(a).total_cmp(&cost(b))
            })
            .unwrap();
        load[best].0 += r;
        load[best].1 += s;
        for i in m {
            fold[i] = best;
        }
    }
    fold
}

/// Moves each of `groups` whole into one of `k` folds: the fold holding most
/// of its members, a tie going to the fold of its lowest-index member among
/// the tied folds. Every other sample keeps its fold. Groups are disjoint and
/// list their members in ascending order. A tie broken by the lowest fold
/// index would pile split pairs into fold 0; a member's own fold is spread
/// evenly by the shuffle that assigned it.
pub fn isolate_groups(folds: &[usize], groups: &[Vec<usize>], k: usize) -> Vec<usize> {
    let mut isolated = folds.to_vec();
    for members in groups {
        let mut count = vec![0usize; k];
        for &i in members {
            count[folds[i]] += 1;
        }
        let top = count.iter().copied().max().unwrap_or(0);
        if let Some(anchor) = members.iter().map(|&i| folds[i]).find(|&f| count[f] == top) {
            for &i in members {
                isolated[i] = anchor;
            }
        }
    }
    isolated
}

pub fn fold_masks(folds: &[usize], k: usize) -> Vec<(Bits, Bits)> {
    (0..k)
        .map(|f| {
            let test = Bits::from_indices(folds.len(), (0..folds.len()).filter(|&i| folds[i] == f));
            (test.not(), test)
        })
        .collect()
}

/// Build inner folds from training samples only. Outside samples receive a
/// sentinel and cannot change the shuffle, class balance, or group allocation.
pub fn training_folds(
    y: &Bits,
    train: &Bits,
    groups: Option<&[u32]>,
    k: usize,
    seed: u64,
) -> Vec<usize> {
    assert_eq!(y.len(), train.len());
    let members: Vec<usize> = train.iter_ones().collect();
    let local_y = Bits::from_indices(
        members.len(),
        members
            .iter()
            .enumerate()
            .filter_map(|(j, &i)| y.contains(i).then_some(j)),
    );
    let local_folds = if let Some(groups) = groups {
        assert_eq!(y.len(), groups.len());
        let local_groups: Vec<u32> = members.iter().map(|&i| groups[i]).collect();
        grouped_folds(&local_y, &local_groups, k, seed)
    } else {
        stratified_folds(&local_y, k, seed)
    };
    let mut folds = vec![usize::MAX; y.len()];
    for (&i, f) in members.iter().zip(local_folds) {
        folds[i] = f;
    }
    folds
}

/// Reject partitions whose balanced classification scores are undefined.
pub fn validate_folds(y: &Bits, train: &Bits, folds: &[usize], k: usize) -> anyhow::Result<()> {
    anyhow::ensure!(k >= 2, "cross-validation needs at least two folds");
    anyhow::ensure!(
        folds.len() == y.len() && train.len() == y.len(),
        "fold dimensions differ"
    );
    anyhow::ensure!(
        train.iter_ones().all(|i| folds[i] < k),
        "a sample has no test fold"
    );
    for (f, (fit, val)) in fold_masks(folds, k).into_iter().enumerate() {
        for (name, mask) in [("training", fit.and(train)), ("validation", val.and(train))] {
            let pos = mask.and(y).count();
            anyhow::ensure!(
                pos > 0 && pos < mask.count(),
                "fold {f} {name} lacks one or both classes"
            );
        }
    }
    Ok(())
}

/// Hold out whole lineages for tuning when the outer evaluation is grouped.
pub fn tuning_holdout(
    train: &Bits,
    y: &Bits,
    groups: Option<&[u32]>,
    frac: f64,
    seed: u64,
) -> anyhow::Result<(Bits, Bits)> {
    let (fit, val) = if let Some(groups) = groups {
        let n_groups = train
            .iter_ones()
            .map(|i| groups[i])
            .collect::<std::collections::BTreeSet<_>>()
            .len();
        anyhow::ensure!(n_groups >= 2, "tuning needs at least two training lineages");
        let k = ((1.0 / frac).round() as usize).max(2).min(n_groups);
        let folds = training_folds(y, train, Some(groups), k, seed);
        fold_masks(&folds, k)
            .into_iter()
            .map(|(fit, val)| (fit.and(train), val.and(train)))
            .filter(|(fit, val)| {
                [fit, val].iter().all(|m| {
                    let pos = m.and(y).count();
                    pos > 0 && pos < m.count()
                })
            })
            .min_by(|a, b| {
                let target = train.count() as f64 * frac;
                (a.1.count() as f64 - target)
                    .abs()
                    .total_cmp(&(b.1.count() as f64 - target).abs())
            })
            .ok_or_else(|| anyhow::anyhow!("no lineage tuning partition contains both classes"))?
    } else {
        holdout(train, y, frac, seed)
    };
    for mask in [&fit, &val] {
        let pos = mask.and(y).count();
        anyhow::ensure!(
            pos > 0 && pos < mask.count(),
            "tuning partition lacks one or both classes"
        );
    }
    Ok((fit, val))
}

/// Splits `train` into (fit, validation) with `frac` of each class held out.
pub fn holdout(train: &Bits, y: &Bits, frac: f64, seed: u64) -> (Bits, Bits) {
    let mut rng = ChaCha8Rng::seed_from_u64(seed);
    let mut val = Bits::zeros(train.len());
    for class in [true, false] {
        let mut members: Vec<usize> = train
            .iter_ones()
            .filter(|&i| y.contains(i) == class)
            .collect();
        members.shuffle(&mut rng);
        let take = ((members.len() as f64 * frac).round() as usize)
            .max(1)
            .min(members.len().saturating_sub(1));
        for &i in &members[..take] {
            val.insert(i);
        }
    }
    (train.and_not(&val), val)
}

/// The features a learner may use when fitting on `train`: rows of `x` whose
/// support among training samples lies in [min_support, |train| - min_support],
/// one row per distinct training pattern (the lowest index). Returns the
/// candidate matrix over every column, so held-out samples can be scored, and
/// the source row of each candidate.
///
/// Held-out columns cannot add, drop, or duplicate a candidate. The dataset
/// matrix is a safe superset: it applies the same `min_support` to all
/// samples, which every training subset implies, and its rows are ordered by
/// smallest k-mer, so the lowest-index representative does not depend on
/// held-out presence either.
pub fn training_candidates(x: &BitMatrix, train: &Bits, min_support: u32) -> (BitMatrix, Vec<u32>) {
    assert_eq!(x.cols, train.len());
    let max = train.count().saturating_sub(min_support);
    let mut seen: rustc_hash::FxHashSet<Box<[u64]>> = Default::default();
    let mut sub = BitMatrix::with_cols(x.cols);
    let mut rows = Vec::new();
    let mut masked = vec![0u64; train.words().len()];
    for f in 0..x.rows {
        let row = x.row(f);
        for ((m, &r), &t) in masked.iter_mut().zip(row).zip(train.words()) {
            *m = r & t;
        }
        let s = bits::count(&masked);
        if s < min_support || s > max || !seen.insert(masked.as_slice().into()) {
            continue;
        }
        sub.push_row(row);
        rows.push(f as u32);
    }
    (sub, rows)
}

/// Ranks features by the chi-square statistic of the 2x2 presence-by-class
/// table on the training samples, and returns the top `k` feature indices.
pub fn screen(x: &BitMatrix, y: &Bits, train: &Bits, k: usize) -> Vec<u32> {
    let pos = train.and(y);
    let neg = train.and_not(y);
    let (np, nn) = (pos.count() as f64, neg.count() as f64);
    let n = np + nn;
    let mut scored: Vec<(f64, u32)> = (0..x.rows)
        .into_par_iter()
        .map(|f| {
            let row = x.row(f);
            let a = bits::and_count(row, pos.words()) as f64; // present, resistant
            let b = bits::and_count(row, neg.words()) as f64; // present, susceptible
            let (c, d) = (np - a, nn - b);
            let denom = (a + b) * (c + d) * np * nn;
            let chi = if denom == 0.0 {
                0.0
            } else {
                n * (a * d - b * c).powi(2) / denom
            };
            (chi, f as u32)
        })
        .collect();
    let k = k.min(scored.len());
    if k < scored.len() {
        scored.select_nth_unstable_by(k, |a, b| b.0.total_cmp(&a.0).then(a.1.cmp(&b.1)));
        scored.truncate(k);
    }
    scored.sort_by(|a, b| b.0.total_cmp(&a.0).then(a.1.cmp(&b.1)));
    scored.into_iter().map(|(_, f)| f).collect()
}

/// Ranks features for MIC regression by the between-group sum of squares of
/// log2 MIC (carriers vs non-carriers) over fitting samples. Censored values
/// enter at their bound, which is enough for ranking.
pub fn screen_regression(x: &BitMatrix, value: &[f64], fit: &Bits, k: usize) -> Vec<u32> {
    let n = fit.count() as f64;
    let total: f64 = fit.iter_ones().map(|i| value[i]).sum();
    let mut scored: Vec<(f64, u32)> = (0..x.rows)
        .into_par_iter()
        .map(|f| {
            let (mut sum, mut cnt) = (0.0, 0.0);
            for (wi, (&r, &m)) in x.row(f).iter().zip(fit.words()).enumerate() {
                let mut w = r & m;
                while w != 0 {
                    sum += value[wi * 64 + w.trailing_zeros() as usize];
                    cnt += 1.0;
                    w &= w - 1;
                }
            }
            let rest = n - cnt;
            let score = if cnt == 0.0 || rest == 0.0 {
                0.0
            } else {
                let diff = sum / cnt - (total - sum) / rest;
                diff * diff * cnt * rest / n
            };
            (score, f as u32)
        })
        .collect();
    let k = k.min(scored.len());
    if k < scored.len() {
        scored.select_nth_unstable_by(k, |a, b| b.0.total_cmp(&a.0).then(a.1.cmp(&b.1)));
        scored.truncate(k);
    }
    scored.sort_by(|a, b| b.0.total_cmp(&a.0).then(a.1.cmp(&b.1)));
    scored.into_iter().map(|(_, f)| f).collect()
}

/// A random subset of `k` of `n` indices.
pub fn sample_indices(n: usize, k: usize, rng: &mut impl Rng) -> Vec<u32> {
    if k >= n {
        return (0..n as u32).collect();
    }
    rand::seq::index::sample(rng, n, k)
        .into_iter()
        .map(|i| i as u32)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn isolated_groups_move_whole_to_their_majority_fold_and_tie_on_a_member_fold() {
        let folds = [3, 1, 1, 2, 4, 0, 0, 2, 4];
        let groups = [
            // Two of three members sit in fold 1.
            vec![0, 1, 2],
            // Folds 2 and 0 tie, and the lowest-index member (3) is in fold
            // 2. A lowest-fold-index rule would have picked fold 0.
            vec![3, 5],
            // Folds 4 and 2 tie and sample 4 is in fold 4.
            vec![4, 7],
        ];
        // Samples 6 and 8 are in no group and keep their folds.
        assert_eq!(
            isolate_groups(&folds, &groups, 5),
            [1, 1, 1, 2, 4, 2, 0, 4, 4]
        );
        assert_eq!(isolate_groups(&folds, &[], 5), folds);
    }

    #[test]
    fn a_large_single_class_group_does_not_exclude_the_other_class() {
        let y = Bits::from_indices(200, 0..100);
        let groups: Vec<u32> = (0..200)
            .map(|i| if i < 80 { 0 } else { i as u32 })
            .collect();
        let folds = grouped_folds(&y, &groups, 5, 42);
        validate_folds(&y, &Bits::ones(200), &folds, 5).unwrap();
        assert!(folds[..80].iter().all(|&f| f == folds[0]));
        for (_, test) in fold_masks(&folds, 5) {
            assert_eq!(test.and_not(&y).count(), 20);
        }
    }

    #[test]
    fn inner_partitions_ignore_outer_test_labels_and_groups() {
        let y = Bits::from_indices(24, (0..24).filter(|i| i % 2 == 0));
        let train = Bits::from_indices(24, 0..16);
        let changed_y = Bits::from_indices(24, (0..16).filter(|i| i % 2 == 0));
        let groups: Vec<u32> = (0..24).map(|i| i / 2).collect();
        let mut changed_groups = groups.clone();
        changed_groups[16..].fill(0);
        for grouped in [false, true] {
            let a = training_folds(&y, &train, grouped.then_some(groups.as_slice()), 2, 42);
            let b = training_folds(
                &changed_y,
                &train,
                grouped.then_some(changed_groups.as_slice()),
                2,
                42,
            );
            assert_eq!(a, b);
            assert!(a[16..].iter().all(|&f| f == usize::MAX));
            validate_folds(&y, &train, &a, 2).unwrap();
        }
    }

    #[test]
    fn lineage_tuning_keeps_whole_groups_inside_outer_training() {
        let y = Bits::from_indices(40, (0..40).filter(|i| i % 2 == 0));
        let groups: Vec<u32> = (0..40).map(|i| i / 4).collect();
        let train = Bits::from_indices(40, 0..32);
        let (fit, val) = tuning_holdout(&train, &y, Some(&groups), 0.2, 42).unwrap();
        assert_eq!(fit.and(&val).count(), 0);
        assert_eq!(fit.count() + val.count(), train.count());
        assert_eq!(fit.and_not(&train).count() + val.and_not(&train).count(), 0);
        let fit_groups: std::collections::BTreeSet<_> =
            fit.iter_ones().map(|i| groups[i]).collect();
        assert!(val.iter_ones().all(|i| !fit_groups.contains(&groups[i])));
    }

    #[test]
    fn invalid_folds_cannot_produce_balanced_scores() {
        let y = Bits::from_indices(4, [0, 1]);
        let all = Bits::ones(4);
        assert!(validate_folds(&y, &all, &[0, 0, 1, 1], 2).is_err());
        assert!(validate_folds(&y, &all, &[0, 0, 0, 0], 2).is_err());
        assert!(validate_folds(&y, &all, &[0, 1, 0, usize::MAX], 2).is_err());
        assert!(validate_folds(&y, &all, &[0, 1, 0, 1], 2).is_ok());
    }

    #[test]
    #[should_panic(expected = "one score is required")]
    fn missing_scores_are_not_silently_dropped() {
        binary_metrics(&[true, false], &[0.9], 0.5, true);
    }

    #[test]
    #[should_panic(expected = "scores must be finite")]
    fn nonfinite_scores_are_not_susceptible_calls() {
        binary_metrics(&[true, false], &[f64::NAN, 0.2], 0.5, true);
    }

    #[test]
    fn ranking_metrics_handle_ties_and_reject_invalid_inputs() {
        assert_eq!(auroc(&[true, false], &[0.5, 0.5]), Some(0.5));
        assert_eq!(auprc(&[true, false], &[0.5, 0.5]), Some(0.5));
        assert_eq!(auroc(&[true, false], &[0.5]), None);
        assert_eq!(auprc(&[true, false], &[0.5, f64::NAN]), None);
        let truth = [true, false, true, false];
        let scores = [0.9, 0.7, 0.6, 0.1];
        let points = roc_curve(&truth, &scores, 3);
        assert_eq!(points.len(), 3);
        assert_eq!(points[0], (0.0, 0.0));
        assert_eq!(points[2], (1.0, 1.0));
        assert!(roc_curve(&truth, &scores, 0).is_empty());
        assert!(roc_curve(&[true, true], &[0.9, 0.6], 3).is_empty());
    }

    #[test]
    fn metrics_follow_clinical_definitions() {
        let truth = [
            true, true, true, true, false, false, false, false, false, false,
        ];
        let scores = [0.9, 0.8, 0.7, 0.1, 0.6, 0.2, 0.2, 0.1, 0.05, 0.0];
        let m = binary_metrics(&truth, &scores, 0.5, true);
        assert_eq!((m.tp, m.fn_, m.fp, m.tn), (3, 1, 1, 5));
        assert!((m.vme - 0.25).abs() < 1e-12); // 1 of 4 resistant missed
        assert!((m.me - 1.0 / 6.0).abs() < 1e-12); // 1 of 6 susceptible over-called
        assert!((m.balanced_accuracy - (0.75 + 5.0 / 6.0) / 2.0).abs() < 1e-12);
        // 20 of 24 positive-negative pairs ordered correctly, one tie at 0.1.
        assert!((m.auroc.unwrap() - 20.5 / 24.0).abs() < 1e-12);
    }

    #[test]
    fn auprc_is_one_for_perfect_ranking() {
        let truth = [true, true, false, false];
        assert_eq!(auprc(&truth, &[0.9, 0.8, 0.2, 0.1]), Some(1.0));
        assert_eq!(auroc(&truth, &[0.9, 0.8, 0.2, 0.1]), Some(1.0));
    }

    #[test]
    fn grouped_folds_never_split_a_group() {
        let y = Bits::from_indices(12, [0, 1, 2, 6, 7, 8]);
        let groups = [0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5];
        let folds = grouped_folds(&y, &groups, 3, 1);
        for pair in groups.chunks(2).enumerate() {
            let i = pair.0 * 2;
            assert_eq!(folds[i], folds[i + 1]);
        }
        let used: std::collections::BTreeSet<_> = folds.iter().collect();
        assert_eq!(used.len(), 3);
    }

    #[test]
    fn training_candidates_filter_and_collapse_on_training_samples_only() {
        // Samples 0..6 train, 6..8 held out.
        let train = Bits::from_indices(8, 0..6);
        let mut x = BitMatrix::with_cols(8);
        x.push_row(Bits::from_indices(8, [0, 1, 6]).words()); // kept
        x.push_row(Bits::from_indices(8, [0, 1, 7]).words()); // same training pattern
        x.push_row(Bits::from_indices(8, [0, 6, 7]).words()); // training support 1
        x.push_row(Bits::from_indices(8, [0, 1, 2, 3, 4]).words()); // training absence 1
        x.push_row(Bits::from_indices(8, [2, 3, 7]).words()); // kept
        let (sub, rows) = training_candidates(&x, &train, 2);
        assert_eq!(rows, vec![0, 4]);
        assert_eq!(sub.cols, 8);
        assert_eq!(sub.row(0), x.row(0));
        assert_eq!(sub.row(1), x.row(4));
    }

    #[test]
    fn screening_ranks_the_separating_feature_first() {
        let y = Bits::from_indices(8, [0, 1, 2, 3]);
        let mut x = BitMatrix::with_cols(8);
        x.push_row(Bits::from_indices(8, [0, 4]).words()); // uninformative
        x.push_row(Bits::from_indices(8, [0, 1, 2, 3]).words()); // perfect
        x.push_row(Bits::from_indices(8, [0, 1, 4]).words()); // weak
        let train = Bits::ones(8);
        assert_eq!(screen(&x, &y, &train, 2), vec![1, 2]);
    }
}
