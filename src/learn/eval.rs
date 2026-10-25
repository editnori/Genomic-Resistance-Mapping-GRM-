//! Cross-validated comparison of every learner on one dataset.
//!
//! Two split schemes run side by side:
//! - `random`: stratified k-fold, the scheme most published numbers use;
//! - `lineage`: whole lineage clusters held out together, so a test genome
//!   never has a near-identical relative in its training fold.
//!
//! Their difference measures sensitivity to this split scheme; it does not
//! identify a causal resistance mechanism.
//!
//! Hyperparameters are chosen inside each training fold (inner folds follow
//! the outer scheme), and feature screening sees only fitting samples.
//! Each outer fit also chooses its candidate features from its own training
//! samples (`training_candidates`); inner tuning reuses that outer-training
//! set, so the whole procedure stays a function of outer training data.
//! Every random draw in an outer fold (inner folds, tuning holdouts, boosting
//! column samples, forest subsamples and node features) is seeded from
//! `EvalConfig::seed` and the fold, so `fit_fold` rebuilds any held-out fit.
//!
//! A boosted model that keeps no tree, or a lasso that keeps no weight, is a
//! null fit: it scores the class-balanced 0.5 prior, so it makes no call for
//! its test fold. Metrics cover the called samples and `no_call` counts the
//! rest; a learner with no called sample has no metrics. A fitted model's
//! score is a call even when it is exactly 0.5.

use super::gbt::{self, GbtParams, Objective};
use super::lasso;
use super::scm::{self, Kind, Rule};
use super::tree::{self, ForestParams, TreeParams};
use super::{
    BinaryMetrics, NoCall, binary_metrics, fold_masks, grouped_folds, holdout, roc_curve, screen,
    stratified_folds, training_candidates, training_folds, tuning_holdout, validate_folds,
};
use crate::bits::{BitMatrix, Bits};
use crate::curate::{Breakpoint, Call, Censor, Mic, infer_breakpoint};
use anyhow::{Context, Result, ensure};
use serde::{Deserialize, Serialize};
use std::time::Instant;

#[derive(Clone, Debug, Serialize)]
pub struct EvalConfig {
    pub folds: usize,
    pub inner_folds: usize,
    pub scm_max_rules: usize,
    pub cart_max_depth: u32,
    pub forest: ForestParams,
    pub gbt: GbtParams,
    /// Chi-square-screened features offered to boosting, the lasso, and the
    /// forest, screened on each fit's own fitting samples.
    pub screened_features: usize,
    pub lasso_lambdas: usize,
    pub seed: u64,
}

impl Default for EvalConfig {
    fn default() -> Self {
        EvalConfig {
            folds: 5,
            inner_folds: 3,
            scm_max_rules: 6,
            cart_max_depth: 8,
            forest: ForestParams::default(),
            gbt: GbtParams::default(),
            screened_features: 4000,
            lasso_lambdas: 30,
            seed: 42,
        }
    }
}

pub const LEARNERS: [&str; 5] = ["scm", "cart", "rf", "gbt", "lasso"];

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct LearnerResult {
    pub learner: String,
    /// Metrics on the pooled out-of-fold calls; samples without a call are
    /// left out and counted in `no_call`. None (null in JSON) when no sample
    /// has a call: an empty call set has no error rates, not 0% error.
    pub metrics: Option<BinaryMetrics>,
    /// Balanced accuracy of each outer fold's calls, by fold id; NaN (null
    /// in JSON) for a fold whose fit made no call.
    #[serde(deserialize_with = "deserialize_nullable_scores")]
    pub fold_balanced_accuracy: Vec<f64>,
    pub roc: Vec<(f64, f64)>,
    /// Out-of-fold resistant score per sample, in dataset sample order; NaN
    /// (null in JSON) for a sample without a call.
    #[serde(deserialize_with = "deserialize_nullable_scores")]
    pub scores: Vec<f64>,
    pub conformal: Vec<ConformalPoint>,
    /// Per-fold model size: SCM rules, CART leaves, forest trees, boosted
    /// trees kept by early stopping, or lasso non-zero weights. Zero for
    /// boosting or the lasso is a null fit. Rows without `no_call` scored
    /// that fold's samples 0.5, which the 0.5 cutoff calls resistant; rows
    /// with it leave them uncalled.
    pub fold_complexity: Vec<f64>,
    pub seconds: f64,
    /// Samples in folds whose fit was null, by reported class. Absent on
    /// catalog rows, which call every sample, and on rows written before
    /// null fits became no-calls.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub no_call: Option<NoCall>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SplitResult {
    pub scheme: String,
    pub folds: usize,
    pub n_groups: usize,
    /// Outer fold that held out each sample, in dataset sample order: the
    /// fold whose fit gave its score. Empty when the split is unavailable.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub test_fold: Vec<usize>,
    pub learners: Vec<LearnerResult>,
    /// The planned partitions could not support evaluation; no scores exist.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub unavailable_reason: Option<String>,
}

/// Class-conditional (Mondrian) split-conformal calls at one error level.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ConformalPoint {
    pub alpha: f64,
    /// Share of samples that receive a single call (R or S).
    pub coverage: f64,
    pub vme: f64,
    pub me: f64,
    /// Share receiving both labels ("cannot call").
    pub ambiguous: f64,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct MicResult {
    pub n: usize,
    pub n_exact: usize,
    /// Within one doubling dilution, exact MICs only.
    pub essential_agreement: f64,
    /// Censored MICs counted as agreeing when the prediction is consistent
    /// with the bound to within one dilution.
    pub essential_agreement_all: f64,
    pub categorical: Option<BinaryMetrics>,
    pub breakpoint: Option<Breakpoint>,
    /// Categorical predictions use only their training fold's empirical
    /// threshold and are scored against independently reported R/S labels.
    #[serde(default)]
    pub fold_breakpoints: Vec<Option<Breakpoint>>,
    /// Observed (exact) by predicted dilution counts: (observed, predicted, n).
    pub confusion: Vec<(i32, i32, u32)>,
    /// Out-of-fold predicted log2 MIC per sample (NaN without an MIC).
    #[serde(deserialize_with = "deserialize_nullable_scores")]
    pub predicted: Vec<f64>,
}

/// Reads values that serde_json wrote as null from NaN back as NaN.
fn deserialize_nullable_scores<'de, D: serde::Deserializer<'de>>(
    d: D,
) -> Result<Vec<f64>, D::Error> {
    Ok(Vec::<Option<f64>>::deserialize(d)?
        .into_iter()
        .map(|v| v.unwrap_or(f64::NAN))
        .collect())
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct FinalScm {
    pub kind: Kind,
    pub p: f64,
    pub cv_error: f64,
    pub rules: Vec<Rule>,
    pub steps: Vec<scm::Step>,
    pub importance: Vec<f64>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct WeightedFeature {
    pub feature: u32,
    pub weight: f64,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct FinalModels {
    pub scm: FinalScm,
    pub cart: tree::Tree,
    pub cart_depth: u32,
    pub lasso: Option<lasso::LassoModel>,
    pub rf_importance: Vec<WeightedFeature>,
    pub gbt_importance: Vec<WeightedFeature>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct EvalOutput {
    pub splits: Vec<SplitResult>,
    pub mic: Option<MicResult>,
    pub finals: FinalModels,
    pub lineage_threshold: Option<f64>,
    pub seconds: f64,
}

pub struct MicInput<'a> {
    /// The `y` passed to evaluate must contain independently reported labels.
    pub mic: &'a [Option<Mic>],
    pub breakpoint: Option<Breakpoint>,
}

struct Scored {
    /// A score for every sample, or None for a null fit, which makes no call.
    scores: Option<Vec<f64>>,
    complexity: f64,
}

fn fit_scm(x: &BitMatrix, y: &Bits, train: &Bits, inner: &[usize], cfg: &EvalConfig) -> Scored {
    let c = scm::select(x, y, train, cfg.scm_max_rules, inner, cfg.inner_folds);
    let m = scm::fit(x, y, train, c.kind, c.p, c.n_rules, false);
    let pred = m.predict(x);
    Scored {
        scores: Some((0..x.cols).map(|i| pred.contains(i) as u8 as f64).collect()),
        complexity: m.rules.len() as f64,
    }
}

fn fit_cart(x: &BitMatrix, y: &Bits, train: &Bits, inner: &[usize], cfg: &EvalConfig) -> Scored {
    let (depth, _) = tree::select_depth(x, y, train, inner, cfg.inner_folds, cfg.cart_max_depth);
    let t = tree::fit_tree(
        x,
        y,
        train,
        &TreeParams {
            max_depth: depth,
            min_leaf: 1,
            max_features: None,
            seed: 0,
        },
    );
    let leaves = t.nodes.iter().filter(|n| n.split.is_none()).count();
    Scored {
        scores: Some(t.predict(x, depth)),
        complexity: leaves as f64,
    }
}

/// The forest fits on the same number of chi-square-screened features as
/// boosting and the lasso, screened on its training samples. Its trees name
/// positions in that screen, so importance is returned keyed by rows of `x`.
fn fit_rf(
    x: &BitMatrix,
    y: &Bits,
    train: &Bits,
    cfg: &EvalConfig,
    seed: u64,
) -> (Scored, std::collections::HashMap<u32, f64>) {
    let feats = screen(x, y, train, cfg.screened_features);
    let mut screened = BitMatrix::with_cols(x.cols);
    for &f in &feats {
        screened.push_row(x.row(f as usize));
    }
    let forest = tree::fit_forest(&screened, y, train, &ForestParams { seed, ..cfg.forest });
    let importance = forest
        .importance()
        .into_iter()
        .map(|(f, w)| (feats[f as usize], w))
        .collect();
    (
        Scored {
            scores: Some(forest.predict(&screened)),
            complexity: forest.trees.len() as f64,
        },
        importance,
    )
}

/// Scores from a boosted model, or None when early stopping kept no tree:
/// that model is the 0.5 prior, a null fit. A model with trees is scored
/// even where it gives exactly 0.5.
fn gbt_scores(m: &gbt::Gbt, x: &BitMatrix) -> Option<Vec<f64>> {
    (!m.trees.is_empty()).then(|| m.predict_proba(x))
}

/// Scores from a lasso refit, or None when there is no refit or it keeps no
/// weight: the intercept alone is the 0.5 prior, a null fit. A refit with
/// weights is scored even where it gives exactly 0.5.
fn lasso_scores(m: Option<&lasso::LassoModel>, x: &BitMatrix) -> Option<Vec<f64>> {
    m.filter(|m| !m.coef.is_empty()).map(|m| m.predict_proba(x))
}

fn fit_gbt(
    x: &BitMatrix,
    y: &Bits,
    train: &Bits,
    groups: Option<&[u32]>,
    cfg: &EvalConfig,
    seed: u64,
) -> Result<(Scored, gbt::Gbt)> {
    let (fit_m, val) = tuning_holdout(train, y, groups, 0.15, seed)?;
    let feats = screen(x, y, &fit_m, cfg.screened_features);
    let m = gbt::fit(
        x,
        &Objective::Logistic { y },
        &fit_m,
        &val,
        &feats,
        &GbtParams { seed, ..cfg.gbt },
    );
    Ok((
        Scored {
            scores: gbt_scores(&m, x),
            complexity: m.trees.len() as f64,
        },
        m,
    ))
}

fn fit_lasso(
    x: &BitMatrix,
    y: &Bits,
    train: &Bits,
    groups: Option<&[u32]>,
    cfg: &EvalConfig,
    seed: u64,
) -> Result<(Scored, Option<lasso::LassoModel>)> {
    let (fit_m, val) = tuning_holdout(train, y, groups, 0.2, seed)?;
    let feats = screen(x, y, &fit_m, cfg.screened_features);
    let path = lasso::fit_path(x, y, &fit_m, &feats, cfg.lasso_lambdas, 400);
    let best = path
        .iter()
        .enumerate()
        .map(|(i, m)| (i, lasso::log_loss(&m.predict_proba(x), y, &val)))
        .min_by(|a, b| a.1.total_cmp(&b.1))
        .map(|(i, _)| i);
    // An empty path means no penalty admits a feature: no refit.
    let model = best.and_then(|best| {
        // Refit on the whole training fold down to the chosen lambda.
        let lambdas: Vec<f64> = path[..=best].iter().map(|m| m.lambda).collect();
        lasso::refit_path(x, y, train, &feats, &lambdas)
            .into_iter()
            .last()
    });
    Ok((
        Scored {
            scores: lasso_scores(model.as_ref(), x),
            complexity: model.as_ref().map_or(0, |m| m.coef.len()) as f64,
        },
        model,
    ))
}

/// The seed of every random draw in outer fold `fold`.
fn fold_seed(cfg: &EvalConfig, fold: usize) -> u64 {
    cfg.seed + fold as u64 * 101
}

/// One learner's fit on outer fold `fold`'s training samples, scoring every
/// sample. Evaluation and any rebuild of a saved score call this alone.
fn fit_fold(
    learner: &str,
    x: &BitMatrix,
    y: &Bits,
    train: &Bits,
    groups: Option<&[u32]>,
    cfg: &EvalConfig,
    fold: usize,
) -> Result<Scored> {
    let seed = fold_seed(cfg, fold);
    let inner = || training_folds(y, train, groups, cfg.inner_folds, seed);
    Ok(match learner {
        "scm" => fit_scm(x, y, train, &inner(), cfg),
        "cart" => fit_cart(x, y, train, &inner(), cfg),
        "rf" => fit_rf(x, y, train, cfg, seed).0,
        "gbt" => fit_gbt(x, y, train, groups, cfg, seed)?.0,
        "lasso" => fit_lasso(x, y, train, groups, cfg, seed)?.0,
        _ => unreachable!(),
    })
}

/// Pools one learner's outer-fold fits, in fold order. Sample i takes its
/// score from the fit of fold `test_fold[i]`, the one fit that did not train
/// on it. A null fit leaves its fold uncalled: those scores stay NaN, the
/// fold's balanced accuracy is NaN, and every metric covers called samples.
/// When every fit is null there are no metrics and no ROC points.
fn pool_folds(
    learner: &str,
    truth: &[bool],
    test_fold: &[usize],
    fits: &[Scored],
    seconds: f64,
) -> Result<LearnerResult> {
    let n = truth.len();
    ensure!(
        test_fold.len() == n && test_fold.iter().all(|&f| f < fits.len()),
        "a sample has no outer test fold"
    );
    let mut scores = vec![f64::NAN; n];
    for (i, &f) in test_fold.iter().enumerate() {
        if let Some(s) = &fits[f].scores {
            ensure!(
                s[i].is_finite(),
                "{learner} produced a non-finite test score"
            );
            scores[i] = s[i];
        }
    }
    let calls = |keep: &dyn Fn(usize) -> bool| -> (Vec<bool>, Vec<f64>) {
        (0..n)
            .filter(|&i| keep(i) && !scores[i].is_nan())
            .map(|i| (truth[i], scores[i]))
            .unzip()
    };
    let fold_balanced_accuracy = fits
        .iter()
        .enumerate()
        .map(|(f, fit)| {
            if fit.scores.is_none() {
                return f64::NAN;
            }
            let (t, s) = calls(&|i| test_fold[i] == f);
            binary_metrics(&t, &s, 0.5, false).balanced_accuracy
        })
        .collect();
    let uncalled = |class: bool| {
        (0..n)
            .filter(|&i| scores[i].is_nan() && truth[i] == class)
            .count()
    };
    let no_call = NoCall {
        resistant: uncalled(true),
        susceptible: uncalled(false),
    };
    let ranked = !matches!(learner, "scm");
    let (t, s) = calls(&|_| true);
    Ok(LearnerResult {
        learner: learner.to_string(),
        metrics: (!t.is_empty()).then(|| binary_metrics(&t, &s, 0.5, ranked)),
        fold_balanced_accuracy,
        roc: if ranked {
            roc_curve(&t, &s, 60)
        } else {
            vec![]
        },
        // Preserve the precision used for thresholding and ranks.
        scores,
        // A separate calibration set is not currently fitted.
        conformal: vec![],
        fold_complexity: fits.iter().map(|f| f.complexity).collect(),
        seconds,
        no_call: Some(no_call),
    })
}

/// Exploratory threshold diagnostic on pooled OOF scores. These predictions
/// share training data across halves, so this is NOT split-conformal
/// validation and must not be included as such in evaluation reports.
pub fn conformal(truth: &[bool], scores: &[f64], alphas: &[f64], seed: u64) -> Vec<ConformalPoint> {
    let folds = stratified_folds(
        &Bits::from_indices(truth.len(), (0..truth.len()).filter(|&i| truth[i])),
        2,
        seed,
    );
    alphas
        .iter()
        .map(|&alpha| {
            let (mut single, mut amb, mut fnr_n, mut fnr_d, mut fpr_n, mut fpr_d) =
                (0.0f64, 0.0f64, 0.0f64, 0.0f64, 0.0f64, 0.0f64);
            for half in 0..2 {
                let cal = |class: bool| -> f64 {
                    // Nonconformity: 1 - score for R, score for S.
                    let mut nc: Vec<f64> = (0..truth.len())
                        .filter(|&i| folds[i] == half && truth[i] == class)
                        .map(|i| if class { 1.0 - scores[i] } else { scores[i] })
                        .collect();
                    nc.sort_by(f64::total_cmp);
                    let n = nc.len();
                    if n == 0 {
                        return 1.0;
                    }
                    let rank = (((n + 1) as f64) * (1.0 - alpha)).ceil() as usize;
                    if rank > n {
                        f64::INFINITY
                    } else {
                        nc[rank - 1]
                    }
                };
                let (q_r, q_s) = (cal(true), cal(false));
                for i in (0..truth.len()).filter(|&i| folds[i] != half) {
                    let has_r = 1.0 - scores[i] <= q_r;
                    let has_s = scores[i] <= q_s;
                    match (has_r, has_s) {
                        (true, true) | (false, false) => amb += 1.0,
                        (true, false) => {
                            single += 1.0;
                            if !truth[i] {
                                fpr_n += 1.0;
                            }
                        }
                        (false, true) => {
                            single += 1.0;
                            if truth[i] {
                                fnr_n += 1.0;
                            }
                        }
                    }
                    if truth[i] { fnr_d += 1.0 } else { fpr_d += 1.0 }
                }
            }
            let n = truth.len() as f64;
            ConformalPoint {
                alpha,
                coverage: single / n,
                ambiguous: amb / n,
                vme: fnr_n / fnr_d.max(1.0),
                me: fpr_n / fpr_d.max(1.0),
            }
        })
        .collect()
}

/// `min_support` is the per-fit candidate policy: a feature must be present
/// in, and absent from, at least that many training samples.
pub fn evaluate(
    x: &BitMatrix,
    y: &Bits,
    groups: &[u32],
    min_support: u32,
    mic: Option<MicInput>,
    cfg: &EvalConfig,
) -> Result<EvalOutput> {
    let started = Instant::now();
    let n = x.cols;
    ensure!(
        n == y.len() && n == groups.len(),
        "evaluation dimensions differ"
    );
    ensure!(min_support >= 1, "candidate features need min_support >= 1");
    ensure!(
        cfg.folds >= 2 && cfg.inner_folds >= 2,
        "cross-validation needs at least two folds"
    );
    ensure!(
        cfg.scm_max_rules > 0
            && cfg.cart_max_depth > 0
            && cfg.lasso_lambdas > 0
            && cfg.forest.n_trees > 0,
        "learner search sizes must be positive"
    );
    if let Some(m) = &mic {
        ensure!(m.mic.len() == n, "MIC dimensions differ");
    }
    let all = Bits::ones(n);
    let truth: Vec<bool> = (0..n).map(|i| y.contains(i)).collect();
    let n_groups = groups
        .iter()
        .collect::<std::collections::BTreeSet<_>>()
        .len();
    let mut splits = Vec::new();

    let schemes = [
        ("random", stratified_folds(y, cfg.folds, cfg.seed)),
        ("lineage", grouped_folds(y, groups, cfg.folds, cfg.seed)),
    ];

    for (scheme, folds) in &schemes {
        let tuning_groups = (*scheme == "lineage").then_some(groups);
        // Check the complete nested partition plan before fitting any model.
        // A missing lineage estimate must not erase a valid random estimate.
        let plan = (|| -> Result<()> {
            if *scheme == "lineage" {
                ensure!(
                    n_groups >= cfg.folds * 2,
                    "fewer than twice as many lineages as requested folds"
                );
            }
            validate_folds(y, &all, folds, cfg.folds)
                .with_context(|| format!("invalid {scheme} outer evaluation"))?;
            for (fi, (train, _)) in fold_masks(folds, cfg.folds).into_iter().enumerate() {
                let seed = fold_seed(cfg, fi);
                let inner = training_folds(y, &train, tuning_groups, cfg.inner_folds, seed);
                validate_folds(y, &train, &inner, cfg.inner_folds).with_context(|| {
                    format!("invalid {scheme} inner evaluation in outer fold {fi}")
                })?;
                for fraction in [0.15, 0.2] {
                    tuning_holdout(&train, y, tuning_groups, fraction, seed).with_context(
                        || format!("invalid {scheme} tuning holdout in outer fold {fi}"),
                    )?;
                }
            }
            Ok(())
        })();
        match plan {
            Ok(()) => {}
            Err(error) if *scheme == "lineage" => {
                splits.push(SplitResult {
                    scheme: scheme.to_string(),
                    folds: cfg.folds,
                    n_groups,
                    test_fold: vec![],
                    learners: vec![],
                    unavailable_reason: Some(format!("{error:#}")),
                });
                continue;
            }
            Err(error) => return Err(error),
        }
        let mut fits: Vec<Vec<Scored>> = LEARNERS.iter().map(|_| Vec::new()).collect();
        let mut seconds = vec![0.0; LEARNERS.len()];
        for (fi, (train, _)) in fold_masks(folds, cfg.folds).into_iter().enumerate() {
            let (cx, _) = training_candidates(x, &train, min_support);
            for (li, learner) in LEARNERS.iter().enumerate() {
                let t0 = Instant::now();
                fits[li].push(fit_fold(learner, &cx, y, &train, tuning_groups, cfg, fi)?);
                seconds[li] += t0.elapsed().as_secs_f64();
            }
        }
        let learners = LEARNERS
            .iter()
            .zip(&fits)
            .zip(seconds)
            .map(|((name, fits), seconds)| {
                pool_folds(name, &truth, folds, fits, seconds)
                    .with_context(|| format!("{scheme} split"))
            })
            .collect::<Result<_>>()?;
        splits.push(SplitResult {
            scheme: scheme.to_string(),
            folds: cfg.folds,
            n_groups,
            test_fold: folds.clone(),
            learners,
            unavailable_reason: None,
        });
    }

    let mic_result = mic.and_then(|m| evaluate_mic(x, y, min_support, &m, cfg));
    let finals = final_models(x, y, min_support, cfg)?;
    Ok(EvalOutput {
        splits,
        mic: mic_result,
        finals,
        lineage_threshold: None,
        seconds: started.elapsed().as_secs_f64(),
    })
}

fn top_weighted(imp: std::collections::HashMap<u32, f64>, k: usize) -> Vec<WeightedFeature> {
    let mut v: Vec<WeightedFeature> = imp
        .into_iter()
        .map(|(feature, weight)| WeightedFeature { feature, weight })
        .collect();
    v.sort_by(|a, b| {
        b.weight
            .total_cmp(&a.weight)
            .then(a.feature.cmp(&b.feature))
    });
    v.truncate(k);
    v
}

fn final_models(
    x: &BitMatrix,
    y: &Bits,
    min_support: u32,
    cfg: &EvalConfig,
) -> Result<FinalModels> {
    let all = Bits::ones(x.cols);
    // Final models train on every sample, so their candidates are the dataset
    // matrix itself; that keeps their feature ids valid for annotation.
    ensure!(
        training_candidates(x, &all, min_support).1.len() == x.rows,
        "final models need a matrix filtered and collapsed on all samples at min_support {min_support}"
    );
    let inner = stratified_folds(y, cfg.inner_folds, cfg.seed + 7);
    validate_folds(y, &all, &inner, cfg.inner_folds)?;
    let choice = scm::select(x, y, &all, cfg.scm_max_rules, &inner, cfg.inner_folds);
    let model = scm::fit(x, y, &all, choice.kind, choice.p, choice.n_rules, true);
    let importance = scm::rule_importance(&model, x, y, &all);
    let (depth, _) = tree::select_depth(x, y, &all, &inner, cfg.inner_folds, cfg.cart_max_depth);
    let cart = tree::fit_tree(
        x,
        y,
        &all,
        &TreeParams {
            max_depth: depth,
            min_leaf: 1,
            max_features: None,
            seed: 0,
        },
    );
    let (_, lasso_model) = fit_lasso(x, y, &all, None, cfg, cfg.seed + 11)?;
    let (_, rf_importance) = fit_rf(x, y, &all, cfg, cfg.seed + 17);
    let (_, gbt_model) = fit_gbt(x, y, &all, None, cfg, cfg.seed + 13)?;
    Ok(FinalModels {
        scm: FinalScm {
            kind: model.kind,
            p: model.p,
            cv_error: choice.cv_error,
            rules: model.rules,
            steps: model.steps,
            importance,
        },
        cart,
        cart_depth: depth,
        lasso: lasso_model,
        rf_importance: top_weighted(rf_importance, 50),
        gbt_importance: top_weighted(gbt_model.importance(), 50),
    })
}

fn evaluate_mic(
    x: &BitMatrix,
    y: &Bits,
    min_support: u32,
    m: &MicInput,
    cfg: &EvalConfig,
) -> Option<MicResult> {
    let n = x.cols;
    let has: Vec<usize> = (0..n).filter(|&i| m.mic[i].is_some()).collect();
    let n_exact = has
        .iter()
        .filter(|&&i| m.mic[i].unwrap().censor == Censor::Exact)
        .count();
    if has.len() < 100 || n_exact < 50 {
        return None;
    }
    let value: Vec<f64> = (0..n)
        .map(|i| m.mic[i].map_or(0.0, |v| v.log2 as f64))
        .collect();
    let censor: Vec<Censor> = (0..n)
        .map(|i| m.mic[i].map_or(Censor::Exact, |v| v.censor))
        .collect();
    let with_mic = Bits::from_indices(n, has.iter().copied());
    let folds = training_folds(y, &with_mic, None, cfg.folds, cfg.seed);
    let mut predicted = vec![f64::NAN; n];
    let mut fold_breakpoints = Vec::new();
    for (fi, (train, test)) in fold_masks(&folds, cfg.folds).into_iter().enumerate() {
        let train = train.and(&with_mic);
        let test = test.and(&with_mic);
        if train.count() == 0 || test.count() == 0 {
            return None;
        }
        let pairs: Vec<(Mic, Call)> = train
            .iter_ones()
            .map(|i| {
                (
                    m.mic[i].unwrap(),
                    if y.contains(i) { Call::R } else { Call::S },
                )
            })
            .collect();
        fold_breakpoints.push(infer_breakpoint(&pairs));
        let (fit_m, val) = holdout(&train, y, 0.15, cfg.seed + fi as u64);
        if fit_m.count() == 0 || val.count() == 0 {
            return None;
        }
        let (cx, _) = training_candidates(x, &train, min_support);
        let x = &cx;
        let feats = super::screen_regression(x, &value, &fit_m, cfg.screened_features);
        let obj = Objective::Censored {
            value: &value,
            censor: &censor,
        };
        let model = gbt::fit(
            x,
            &obj,
            &fit_m,
            &val,
            &feats,
            &GbtParams {
                seed: cfg.seed + fi as u64,
                ..cfg.gbt
            },
        );
        let pred = model.margin(x);
        for i in test.iter_ones() {
            predicted[i] = pred[i];
        }
    }
    if has.iter().any(|&i| !predicted[i].is_finite()) {
        return None;
    }
    let (mut ea, mut ea_all) = (0usize, 0usize);
    let mut confusion: std::collections::BTreeMap<(i32, i32), u32> = Default::default();
    for &i in &has {
        let p = predicted[i].round() as i32;
        let mic = m.mic[i].unwrap();
        let ok = match mic.censor {
            Censor::Exact => {
                *confusion.entry((mic.log2, p)).or_default() += 1;
                let ok = (p - mic.log2).abs() <= 1;
                ea += ok as usize;
                ok
            }
            Censor::AtMost => p <= mic.log2 + 1,
            Censor::Above => p >= mic.log2,
        };
        ea_all += ok as usize;
    }
    let categorical = mic_categorical(y, &has, &predicted, &folds, &fold_breakpoints);
    Some(MicResult {
        n: has.len(),
        n_exact,
        essential_agreement: ea as f64 / n_exact.max(1) as f64,
        essential_agreement_all: ea_all as f64 / has.len() as f64,
        categorical,
        breakpoint: m.breakpoint,
        fold_breakpoints,
        confusion: confusion.into_iter().map(|((o, p), c)| (o, p, c)).collect(),
        predicted,
    })
}

fn mic_categorical(
    y: &Bits,
    has: &[usize],
    predicted: &[f64],
    folds: &[usize],
    breakpoints: &[Option<Breakpoint>],
) -> Option<BinaryMetrics> {
    // Do not silently score a favorable subset if one fold cannot estimate
    // a threshold. The reference remains the reported phenotype, even when
    // its measured MIC disagrees with an empirical threshold.
    let truth: Vec<bool> = has.iter().map(|&i| y.contains(i)).collect();
    let calls: Option<Vec<f64>> = has
        .iter()
        .map(|&i| {
            breakpoints[folds[i]].map(|bp| (predicted[i].round() as i32 > bp.r_above) as u8 as f64)
        })
        .collect();
    Some(binary_metrics(&truth, &calls?, 0.5, false))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn small_config() -> EvalConfig {
        EvalConfig {
            folds: 2,
            inner_folds: 2,
            scm_max_rules: 1,
            cart_max_depth: 1,
            forest: ForestParams {
                n_trees: 2,
                max_depth: 2,
                ..Default::default()
            },
            gbt: GbtParams {
                rounds: 3,
                max_depth: 1,
                ..Default::default()
            },
            screened_features: 2,
            lasso_lambdas: 3,
            ..Default::default()
        }
    }

    #[test]
    fn evaluation_rejects_invalid_folds_and_emits_only_held_out_scores() {
        let n = 80;
        let y = Bits::from_indices(n, (0..n).filter(|i| i % 2 == 0));
        let mut x = BitMatrix::with_cols(n);
        x.push_row(y.words());
        let groups: Vec<u32> = (0..n).map(|i| (i / 4) as u32).collect();
        let mut cfg = small_config();
        cfg.folds = 1;
        assert!(evaluate(&x, &y, &groups, 3, None, &cfg).is_err());
        cfg.folds = 2;
        let mut duplicated = x.clone();
        duplicated.push_row(y.words());
        assert!(evaluate(&duplicated, &y, &groups, 3, None, &cfg).is_err());
        let out = evaluate(&x, &y, &groups, 3, None, &cfg).unwrap();
        assert_eq!(out.splits.len(), 2);
        for split in out.splits {
            assert!(split.unavailable_reason.is_none());
            for learner in split.learners {
                assert_eq!(learner.scores.len(), n);
                assert!(learner.scores.iter().all(|s| s.is_finite()));
                let m = learner.metrics.unwrap();
                assert_eq!(m.n, n);
                assert_eq!(m.tp + m.fp + m.tn + m.fn_, n);
                assert_eq!(learner.fold_balanced_accuracy.len(), cfg.folds);
                assert!(learner.conformal.is_empty());
            }
        }
    }

    #[test]
    fn unsupported_lineage_plan_preserves_the_random_evaluation() {
        let n = 80;
        let y = Bits::from_indices(n, 0..40);
        let mut x = BitMatrix::with_cols(n);
        x.push_row(y.words());
        // Every R sample belongs to one lineage, so a held-out lineage fold
        // necessarily loses R from either training or validation.
        let groups: Vec<u32> = (0..n).map(|i| if i < 40 { 0 } else { i as u32 }).collect();
        let out = evaluate(&x, &y, &groups, 3, None, &small_config()).unwrap();
        let random = &out.splits[0];
        assert_eq!(random.scheme, "random");
        assert!(random.unavailable_reason.is_none());
        assert_eq!(random.learners.len(), LEARNERS.len());
        let lineage = &out.splits[1];
        assert_eq!(lineage.scheme, "lineage");
        assert!(lineage.learners.is_empty());
        assert!(
            lineage
                .unavailable_reason
                .as_ref()
                .unwrap()
                .contains("lacks one or both classes")
        );
        let encoded = serde_json::to_string(lineage).unwrap();
        let decoded: SplitResult = serde_json::from_str(&encoded).unwrap();
        assert_eq!(decoded.unavailable_reason, lineage.unavailable_reason);
    }

    #[test]
    fn lineage_splits_without_enough_whole_lineages_are_unavailable() {
        let n = 40;
        let y = Bits::from_indices(n, 0..20);
        let mut x = BitMatrix::with_cols(n);
        x.push_row(y.words());
        let cases: [(Vec<u32>, &str); 2] = [
            // Three mixed lineages cannot fill two folds twice over.
            ((0..n).map(|i| (i % 3) as u32).collect(), "fewer than twice"),
            // Two resistant lineages and twenty susceptible singletons. Each
            // outer training fold keeps one resistant lineage, so lineage
            // inner folds cannot put resistant samples on both sides. Random
            // inner folds could, and must not stand in for them.
            (
                (0..n)
                    .map(|i| if i < 20 { (i / 10) as u32 } else { i as u32 })
                    .collect(),
                "inner evaluation",
            ),
        ];
        for (groups, expected) in cases {
            let out = evaluate(&x, &y, &groups, 3, None, &small_config()).unwrap();
            assert!(out.splits[0].unavailable_reason.is_none());
            let lineage = &out.splits[1];
            assert!(lineage.learners.is_empty());
            let reason = lineage.unavailable_reason.as_deref().unwrap();
            assert!(reason.contains(expected), "{reason}");
        }
    }

    #[test]
    fn saved_scores_are_held_out_fits_that_reproduce_every_reported_metric() {
        use rand::{Rng, SeedableRng};
        let n = 60;
        let y = Bits::from_indices(n, (0..n).filter(|i| i % 2 == 0));
        let mut rng = rand_chacha::ChaCha8Rng::seed_from_u64(5);
        let mut raw = BitMatrix::with_cols(n);
        for k in 0..12 {
            let rate = |r: bool| match (k < 3, r) {
                (true, true) => 0.7,
                (true, false) => 0.2,
                _ => 0.4,
            };
            let present = (0..n).filter(|&i| rng.random_bool(rate(y.contains(i))));
            raw.push_row(Bits::from_indices(n, present.collect::<Vec<_>>()).words());
        }
        // At exactly min_support, this is a candidate for all samples but not
        // for every training fold.
        raw.push_row(Bits::from_indices(n, [1, 2, 3]).words());
        let (x, _) = training_candidates(&raw, &Bits::ones(n), 3);
        let groups: Vec<u32> = (0..n).map(|i| (i % 15) as u32).collect();
        let cfg = small_config();
        let out = evaluate(&x, &y, &groups, 3, None, &cfg).unwrap();
        let truth: Vec<bool> = (0..n).map(|i| y.contains(i)).collect();
        let metrics_on = |idx: &[usize], scores: &[f64], ranked: bool| {
            binary_metrics(
                &idx.iter().map(|&i| truth[i]).collect::<Vec<_>>(),
                &idx.iter().map(|&i| scores[i]).collect::<Vec<_>>(),
                0.5,
                ranked,
            )
        };
        for split in &out.splits {
            assert!(split.unavailable_reason.is_none(), "{split:?}");
            // The saved fold ids are the partition the split was scored on.
            let folds = match split.scheme.as_str() {
                "random" => stratified_folds(&y, cfg.folds, cfg.seed),
                _ => grouped_folds(&y, &groups, cfg.folds, cfg.seed),
            };
            assert_eq!(split.test_fold, folds);
            let tuning_groups = (split.scheme == "lineage").then_some(groups.as_slice());
            for learner in &split.learners {
                let called: Vec<usize> = (0..n).filter(|&i| !learner.scores[i].is_nan()).collect();
                let ranked = learner.learner != "scm";
                assert_eq!(
                    learner.metrics,
                    (!called.is_empty()).then(|| metrics_on(&called, &learner.scores, ranked))
                );
                let uncalled = |class: bool| {
                    (0..n)
                        .filter(|&i| learner.scores[i].is_nan() && truth[i] == class)
                        .count()
                };
                assert_eq!(
                    learner.no_call,
                    Some(NoCall {
                        resistant: uncalled(true),
                        susceptible: uncalled(false)
                    })
                );
                for (f, (train, test)) in fold_masks(&split.test_fold, cfg.folds)
                    .into_iter()
                    .enumerate()
                {
                    // Rebuild the fold's fit from its training samples alone.
                    let (cx, _) = training_candidates(&x, &train, 3);
                    let refit = fit_fold(&learner.learner, &cx, &y, &train, tuning_groups, &cfg, f)
                        .unwrap();
                    assert_eq!(learner.fold_complexity[f], refit.complexity);
                    let idx: Vec<usize> = test.iter_ones().collect();
                    match &refit.scores {
                        None => {
                            assert!(idx.iter().all(|&i| learner.scores[i].is_nan()));
                            assert!(learner.fold_balanced_accuracy[f].is_nan());
                        }
                        Some(s) => {
                            assert!(idx.iter().all(|&i| learner.scores[i] == s[i]));
                            assert_eq!(
                                learner.fold_balanced_accuracy[f],
                                metrics_on(&idx, &learner.scores, false).balanced_accuracy
                            );
                        }
                    }
                    if learner.learner == "rf" {
                        // The forest scores its own training samples apart
                        // from the held-out fits whose scores were saved.
                        let s = refit.scores.as_ref().unwrap();
                        assert!(train.iter_ones().any(|i| s[i] != learner.scores[i]));
                    }
                }
            }
        }
    }

    #[test]
    fn pooled_metrics_cover_called_samples_and_count_null_folds_apart() {
        // Six samples in three folds; fold 1's fit is null. Each fit also
        // scores samples outside its test fold, which must not be read.
        // Sample 4 scores exactly 0.5 from a fitted model: a resistant call.
        let truth = [true, false, true, false, true, false];
        let test_fold = [0, 0, 1, 1, 2, 2];
        let fitted = |s: [f64; 6]| Scored {
            scores: Some(s.to_vec()),
            complexity: 3.0,
        };
        let fits = [
            fitted([0.9, 0.2, 0.0, 0.0, 0.0, 0.0]),
            Scored {
                scores: None,
                complexity: 0.0,
            },
            fitted([0.0, 0.0, 0.0, 0.0, 0.5, 0.6]),
        ];
        let l = pool_folds("gbt", &truth, &test_fold, &fits, 1.0).unwrap();
        assert_eq!(l.scores[..2], [0.9, 0.2]);
        assert!(l.scores[2].is_nan() && l.scores[3].is_nan());
        assert_eq!(l.scores[4..], [0.5, 0.6]);
        let m = l.metrics.as_ref().unwrap();
        assert_eq!((m.n, m.tp, m.fn_, m.tn, m.fp), (4, 2, 0, 1, 1));
        assert_eq!(
            l.no_call,
            Some(NoCall {
                resistant: 1,
                susceptible: 1
            })
        );
        let (called, called_scores) = ([true, false, true, false], [0.9, 0.2, 0.5, 0.6]);
        assert_eq!(m.auroc, crate::learn::auroc(&called, &called_scores));
        assert_eq!(l.roc, roc_curve(&called, &called_scores, 60));
        assert!(!l.roc.is_empty());
        assert_eq!(l.fold_balanced_accuracy[0], 1.0);
        assert!(l.fold_balanced_accuracy[1].is_nan());
        assert_eq!(l.fold_balanced_accuracy[2], 0.5);
        assert_eq!(l.fold_complexity, [3.0, 0.0, 3.0]);
    }

    #[test]
    fn null_fits_make_no_call_while_fitted_half_scores_are_calls() {
        // With no candidate feature, boosting keeps no tree and the lasso
        // path is empty in every fold: null fits. CART's one leaf and the
        // forest's leaf-only trees score exactly 0.5 too, but they are
        // fitted models, so the 0.5 cutoff calls every sample resistant.
        let n = 80;
        let y = Bits::from_indices(n, (0..n).filter(|i| i % 4 == 0));
        let groups: Vec<u32> = (0..n).map(|i| (i / 4) as u32).collect();
        let cfg = small_config();
        let out = evaluate(&BitMatrix::with_cols(n), &y, &groups, 3, None, &cfg).unwrap();
        let none = NoCall {
            resistant: 0,
            susceptible: 0,
        };
        for split in &out.splits {
            assert!(split.unavailable_reason.is_none());
            for l in &split.learners {
                match l.learner.as_str() {
                    "gbt" | "lasso" => {
                        assert!(l.scores.iter().all(|s| s.is_nan()));
                        assert_eq!(l.metrics, None);
                        assert!(l.roc.is_empty());
                        assert_eq!(
                            l.no_call,
                            Some(NoCall {
                                resistant: 20,
                                susceptible: 60
                            })
                        );
                        assert!(l.fold_balanced_accuracy.iter().all(|b| b.is_nan()));
                        assert_eq!(l.fold_complexity, vec![0.0; cfg.folds]);
                    }
                    "cart" | "rf" => {
                        assert!(l.scores.iter().all(|&s| s == 0.5));
                        let m = l.metrics.as_ref().unwrap();
                        assert_eq!((m.n, m.tp, m.fp), (n, 20, 60));
                        assert_eq!(l.no_call, Some(none));
                    }
                    _ => assert_eq!(l.no_call, Some(none)),
                }
            }
        }
        // A sample without a call is saved as null and read back as NaN.
        let gbt = &out.splits[0].learners[3];
        let saved = serde_json::to_value(gbt).unwrap();
        assert!(
            saved["scores"]
                .as_array()
                .unwrap()
                .iter()
                .all(|s| s.is_null())
        );
        let back: LearnerResult = serde_json::from_value(saved).unwrap();
        assert!(back.scores.iter().all(|s| s.is_nan()));
        assert!(back.fold_balanced_accuracy.iter().all(|b| b.is_nan()));
    }

    #[test]
    fn a_learner_with_every_fold_null_has_no_metrics_roc_or_fold_scores() {
        // E. coli amikacin in revision 5: 1,012 isolates, 30 resistant, and
        // lineage boosting kept no tree in any of its 5 folds. Revision 5
        // scored those isolates 0.5 and reported 0% VME; without a call
        // there is no rate to report.
        let truth: Vec<bool> = (0..1012).map(|i| i < 30).collect();
        let test_fold: Vec<usize> = (0..1012).map(|i| i % 5).collect();
        let fits: Vec<Scored> = (0..5)
            .map(|_| Scored {
                scores: None,
                complexity: 0.0,
            })
            .collect();
        let l = pool_folds("gbt", &truth, &test_fold, &fits, 1.0).unwrap();
        assert_eq!(l.metrics, None);
        assert!(l.roc.is_empty());
        assert!(l.fold_balanced_accuracy.iter().all(|b| b.is_nan()));
        assert_eq!(
            l.no_call,
            Some(NoCall {
                resistant: 30,
                susceptible: 982
            })
        );
        let saved = serde_json::to_value(&l).unwrap();
        assert!(saved["metrics"].is_null());
        assert_eq!(saved["roc"], serde_json::json!([]));
        for key in ["scores", "fold_balanced_accuracy"] {
            assert!(saved[key].as_array().unwrap().iter().all(|v| v.is_null()));
        }
        let back: LearnerResult = serde_json::from_value(saved).unwrap();
        assert_eq!(back.metrics, None);
    }

    #[test]
    fn held_out_labels_cannot_change_any_training_fit() {
        // Samples 0..60 train and 60..120 are held out; y2 changes only
        // held-out labels. Row 0 is absent from every training sample and
        // marks the held-out resistant samples under y. Row 2 marks 24 of 30
        // resistant and 4 of 30 susceptible training samples. Over every
        // sample, y ranks row 0 first (chi-square 40 against 18.6) and y2
        // ranks row 2 first, so a one-feature screen that read held-out
        // labels would give the two fits different features.
        let n = 120;
        let train = Bits::from_indices(n, 0..60);
        let y = Bits::from_indices(n, (0..n).filter(|i| i % 2 == 0));
        let y2 = Bits::from_indices(
            n,
            (0..n).filter(|&i| if i < 60 { i % 2 == 0 } else { i < 90 }),
        );
        let mut x = BitMatrix::with_cols(n);
        x.push_row(Bits::from_indices(n, (60..n).filter(|i| i % 2 == 0)).words());
        x.push_row(Bits::from_indices(n, (0..n).filter(|i| i % 3 == 0)).words());
        let marker = (0..60).filter(|&i| if i % 2 == 0 { i % 10 != 0 } else { i < 8 });
        x.push_row(Bits::from_indices(n, marker).words());
        assert_eq!(screen(&x, &y, &train, 1), vec![2]);
        assert_eq!(screen(&x, &y, &Bits::ones(n), 1), vec![0]);
        assert_eq!(screen(&x, &y2, &Bits::ones(n), 1), vec![2]);

        let cfg = EvalConfig {
            screened_features: 1,
            ..small_config()
        };
        for learner in LEARNERS {
            let fit = |labels: &Bits| {
                let s = fit_fold(learner, &x, labels, &train, None, &cfg, 1).unwrap();
                (s.scores, s.complexity)
            };
            assert_eq!(fit(&y), fit(&y2), "{learner}");
        }
        // The forest fits on screen position 0; its importance names row 2.
        let (_, importance) = fit_rf(&x, &y, &train, &cfg, 7);
        assert!(
            !importance.is_empty() && importance.keys().all(|&f| f == 2),
            "{importance:?}"
        );
    }

    #[test]
    fn a_boosted_model_that_keeps_no_tree_is_a_null_fit_without_a_call() {
        // The only feature marks resistant samples among boosting's fitting
        // samples and susceptible ones in its validation slice, so every tree
        // raises validation loss and early stopping keeps none.
        let n = 80;
        let y = Bits::from_indices(n, (0..n).filter(|i| i % 2 == 0));
        let all = Bits::ones(n);
        let (fit, val) = tuning_holdout(&all, &y, None, 0.15, 42).unwrap();
        let mut x = BitMatrix::with_cols(n);
        x.push_row(fit.and(&y).or(&val.and_not(&y)).words());
        let (scored, model) = fit_gbt(&x, &y, &all, None, &small_config(), 42).unwrap();
        assert!(model.trees.is_empty());
        assert!(scored.scores.is_none());
        assert_eq!(scored.complexity, 0.0);
    }

    #[test]
    fn fitted_models_are_called_even_where_they_score_exactly_one_half() {
        let mut x = BitMatrix::with_cols(4);
        x.push_row(Bits::from_indices(4, [0, 1]).words());
        // One kept tree whose leaves add nothing: every margin stays 0.
        let leaf = gbt::GNode {
            feature: None,
            absent: 0,
            present: 0,
            value: 0.0,
            gain: 0.0,
        };
        let split = gbt::GNode {
            feature: Some(0),
            absent: 1,
            present: 2,
            gain: 1.0,
            ..leaf.clone()
        };
        let kept = gbt::Gbt {
            base: 0.0,
            trees: vec![gbt::GTree {
                nodes: vec![split, leaf.clone(), leaf],
            }],
            logistic: true,
        };
        assert_eq!(gbt_scores(&kept, &x), Some(vec![0.5; 4]));
        let no_tree = gbt::Gbt {
            trees: vec![],
            ..kept
        };
        assert_eq!(gbt_scores(&no_tree, &x), None);
        // The intercept cancels the weight on carriers of feature 1, which
        // every sample carries.
        x.push_row(Bits::ones(4).words());
        let model = lasso::LassoModel {
            intercept: -0.75,
            coef: vec![(1, 0.75)],
            lambda: 0.1,
        };
        assert_eq!(lasso_scores(Some(&model), &x), Some(vec![0.5; 4]));
        let no_weight = lasso::LassoModel {
            coef: vec![],
            ..model
        };
        assert_eq!(lasso_scores(Some(&no_weight), &x), None);
        assert_eq!(lasso_scores(None, &x), None);
    }

    #[test]
    fn final_forest_importance_names_dataset_rows() {
        // Row 3 separates the classes. The forest fits on the two screened
        // rows, where row 3 is position 0, so its importance has to be
        // mapped back to row 3 for annotation to name the right feature.
        let n = 60;
        let y = Bits::from_indices(n, (0..n).filter(|i| i % 2 == 0));
        let mut x = BitMatrix::with_cols(n);
        for rule in [|i: usize| i.is_multiple_of(3), |i| i % 5 < 2, |i| i % 7 < 3] {
            x.push_row(Bits::from_indices(n, (0..n).filter(|&i| rule(i))).words());
        }
        x.push_row(y.words());
        let cfg = small_config();
        let screened = screen(&x, &y, &Bits::ones(n), cfg.screened_features);
        assert_eq!(screened[0], 3);
        let finals = final_models(&x, &y, 3, &cfg).unwrap();
        assert_eq!(finals.rf_importance[0].feature, 3);
        assert!(
            finals
                .rf_importance
                .iter()
                .all(|w| screened.contains(&w.feature))
        );
    }

    #[test]
    fn forest_draws_follow_the_evaluation_seed_and_fold() {
        use rand::{Rng, SeedableRng};
        let n = 60;
        let y = Bits::from_indices(n, (0..n).filter(|i| i % 2 == 0));
        let mut rng = rand_chacha::ChaCha8Rng::seed_from_u64(9);
        let mut x = BitMatrix::with_cols(n);
        for _ in 0..6 {
            let present: Vec<usize> = (0..n).filter(|_| rng.random_bool(0.4)).collect();
            x.push_row(Bits::from_indices(n, present).words());
        }
        let train = Bits::from_indices(n, 0..45);
        let forest = |cfg: &EvalConfig, fold| {
            fit_fold("rf", &x, &y, &train, None, cfg, fold)
                .unwrap()
                .scores
                .unwrap()
        };
        let cfg = small_config();
        assert_eq!(forest(&cfg, 0), forest(&cfg, 0));
        assert_ne!(forest(&cfg, 0), forest(&cfg, 1));
        assert_ne!(
            forest(&cfg, 0),
            forest(
                &EvalConfig {
                    seed: 7,
                    ..cfg.clone()
                },
                0
            )
        );
    }

    #[test]
    fn the_recorded_config_leaves_out_seeds_each_fit_replaces() {
        let saved = serde_json::to_value(EvalConfig::default()).unwrap();
        assert_eq!(saved["folds"], 5);
        assert_eq!(saved["seed"], 42);
        assert_eq!(saved["screened_features"], 4000);
        assert_eq!(saved["forest"]["n_trees"], 200);
        assert_eq!(saved["gbt"]["rounds"], 400);
        assert!(saved["forest"].get("seed").is_none());
        assert!(saved["gbt"].get("seed").is_none());
    }

    #[test]
    fn an_empty_lasso_path_is_a_null_fit_without_a_call() {
        // A quarter of the samples are resistant. With no candidate feature,
        // or only one present on every fitting sample, no penalty admits a
        // feature and the path is empty.
        let n = 80;
        let y = Bits::from_indices(n, (0..n).filter(|i| i % 4 == 0));
        let all = Bits::ones(n);
        let (fit, _) = tuning_holdout(&all, &y, None, 0.2, 42).unwrap();
        let mut constant_on_fit = BitMatrix::with_cols(n);
        constant_on_fit.push_row(fit.words());
        for x in [BitMatrix::with_cols(n), constant_on_fit] {
            let (scored, model) = fit_lasso(&x, &y, &all, None, &small_config(), 42).unwrap();
            assert!(model.is_none());
            assert!(scored.scores.is_none());
            assert_eq!(scored.complexity, 0.0);
        }
    }

    #[test]
    fn a_lasso_refit_that_keeps_no_feature_is_a_null_fit_without_a_call() {
        let n = 260;
        let y = Bits::from_indices(n, (0..n).filter(|i| i % 9 == 0));
        let all = Bits::ones(n);
        let cfg = small_config();
        let (fit, val) = tuning_holdout(&all, &y, None, 0.2, 42).unwrap();
        let mut x = BitMatrix::with_cols(n);
        x.push_row(fit.and(&y).or(&val.and_not(&y)).words());
        assert!(!lasso::fit_path(&x, &y, &fit, &[0], cfg.lasso_lambdas, 400).is_empty());
        let (scored, model) = fit_lasso(&x, &y, &all, None, &cfg, 42).unwrap();
        assert!(model.unwrap().coef.is_empty());
        assert!(scored.scores.is_none());
        assert_eq!(scored.complexity, 0.0);
    }

    #[test]
    fn lasso_refit_uses_the_validation_selected_penalty() {
        let n = 120;
        let y = Bits::from_indices(n, (0..n).filter(|i| i % 3 == 0));
        let mut x = BitMatrix::with_cols(n);
        x.push_row(Bits::from_indices(n, (0..n).filter(|i| i % 3 == 0 && i % 7 != 0)).words());
        let all = Bits::ones(n);
        let mut cfg = small_config();
        cfg.lasso_lambdas = 9;
        let (fit, val) = tuning_holdout(&all, &y, None, 0.2, 42).unwrap();
        let feats = screen(&x, &y, &fit, cfg.screened_features);
        let path = lasso::fit_path(&x, &y, &fit, &feats, cfg.lasso_lambdas, 400);
        let selected = path
            .iter()
            .min_by(|a, b| {
                lasso::log_loss(&a.predict_proba(&x), &y, &val).total_cmp(&lasso::log_loss(
                    &b.predict_proba(&x),
                    &y,
                    &val,
                ))
            })
            .unwrap();
        let (_, model) = fit_lasso(&x, &y, &all, None, &cfg, 42).unwrap();
        assert_eq!(model.unwrap().lambda, selected.lambda);
    }

    #[test]
    fn serialized_scores_preserve_threshold_calls_and_ranks() {
        let truth = [false, true, false, true];
        let scores = [0.5 - 1e-12, 0.5 + 1e-12, 0.5, 0.5 - 2e-12];
        let result = LearnerResult {
            learner: "lasso".into(),
            metrics: Some(binary_metrics(&truth, &scores, 0.5, true)),
            fold_balanced_accuracy: vec![],
            roc: roc_curve(&truth, &scores, 60),
            // Infer the stored field's type so this regression also runs
            // when that type is mutated back to the old narrower format.
            scores: scores.iter().map(|&s| s as _).collect(),
            conformal: vec![],
            fold_complexity: vec![],
            seconds: 0.0,
            no_call: None,
        };
        let saved: serde_json::Value =
            serde_json::from_slice(&serde_json::to_vec(&result).unwrap()).unwrap();
        let restored: Vec<f64> = saved["scores"]
            .as_array()
            .unwrap()
            .iter()
            .map(|s| s.as_f64().unwrap())
            .collect();
        assert_eq!(restored, scores);
        assert_eq!(
            Some(binary_metrics(&truth, &restored, 0.5, true)),
            result.metrics
        );
    }

    #[test]
    fn mic_categorical_uses_reported_truth_and_each_training_threshold() {
        let y = Bits::from_indices(4, [0, 2]);
        let bp = |r_above| {
            Some(Breakpoint {
                r_above,
                agreement: 1.0,
                n: 50,
            })
        };
        let thresholds = [bp(0), bp(2)];
        let m = mic_categorical(&y, &[0, 1, 2, 3], &[1.0; 4], &[0, 0, 1, 1], &thresholds).unwrap();
        assert_eq!((m.tp, m.fp, m.tn, m.fn_), (1, 1, 1, 1));
        assert!(
            mic_categorical(&y, &[0, 1, 2, 3], &[1.0; 4], &[0, 0, 1, 1], &[bp(0), None]).is_none()
        );
    }

    #[test]
    fn mic_scores_round_trip_missing_measurements_and_fit_fold_thresholds() {
        let n = 160;
        let y = Bits::from_indices(n, (0..n).filter(|i| i % 2 == 0));
        let mut x = BitMatrix::with_cols(n);
        x.push_row(y.words());
        let mic: Vec<Option<Mic>> = (0..n)
            .map(|i| {
                (i < 120).then_some(Mic {
                    log2: if y.contains(i) { 2 } else { -2 },
                    censor: Censor::Exact,
                })
            })
            .collect();
        let input = MicInput {
            mic: &mic,
            // Deliberately wrong cohort-level breakpoint: it must not decide
            // any held-out categorical call.
            breakpoint: Some(Breakpoint {
                r_above: 99,
                agreement: 1.0,
                n: 999,
            }),
        };
        let mut report = evaluate_mic(&x, &y, 3, &input, &small_config()).unwrap();
        assert_eq!(report.n, 120);
        assert!(report.predicted[..120].iter().all(|v| v.is_finite()));
        assert!(report.predicted[120..].iter().all(|v| v.is_nan()));
        assert_eq!(report.fold_breakpoints.len(), 2);
        assert!(
            report
                .fold_breakpoints
                .iter()
                .all(|b| b.unwrap().n == 60 && b.unwrap().r_above == -2)
        );
        assert_eq!(report.categorical.as_ref().unwrap().n, 120);
        // Rounding to f32 here would cross the MIC dilution boundary.
        report.predicted[0] = (1.5 - 1e-12) as _;
        let encoded = serde_json::to_string(&report).unwrap();
        let decoded: MicResult = serde_json::from_str(&encoded).unwrap();
        assert_eq!(decoded.predicted[0], 1.5 - 1e-12);
        assert_eq!(decoded.predicted[..120], report.predicted[..120]);
        assert!(decoded.predicted[120..].iter().all(|v| v.is_nan()));
    }

    #[test]
    fn held_out_presence_cannot_change_any_training_fit() {
        // Two dataset matrices that agree on the 40 training samples. The
        // second inverts every held-out bit, splits one training pattern into
        // two global features, and adds one present in a single training
        // sample, so the raw matrices differ in row count.
        let n = 60;
        let y = Bits::from_indices(n, (0..n).filter(|i| i % 2 == 0));
        let train = Bits::from_indices(n, 0..40);
        let rules: [fn(usize) -> bool; 4] = [
            |i| i % 2 == 0 && i % 5 != 0,
            |i| i % 3 == 0,
            |i| i % 4 < 2,
            |i| i % 2 == 0 && i % 7 != 1,
        ];
        let row = |rule: fn(usize) -> bool, flip: bool| {
            Bits::from_indices(n, (0..n).filter(|&i| rule(i) != (flip && i >= 40)))
        };
        let mut x = BitMatrix::with_cols(n);
        let mut x2 = BitMatrix::with_cols(n);
        for (k, &rule) in rules.iter().enumerate() {
            x.push_row(row(rule, false).words());
            x2.push_row(row(rule, true).words());
            if k == 1 {
                x2.push_row(row(rule, false).words());
            }
        }
        x2.push_row(Bits::from_indices(n, std::iter::once(0).chain(40..n)).words());
        assert_ne!(x.rows, x2.rows);

        let cfg = small_config();
        let fits = |m: &BitMatrix| -> Vec<Vec<f64>> {
            let (cx, _) = training_candidates(m, &train, 2);
            LEARNERS
                .iter()
                .map(|learner| {
                    let s = fit_fold(learner, &cx, &y, &train, None, &cfg, 1).unwrap();
                    let s = s.scores.expect("every learner fits this matrix");
                    assert!(s.iter().all(|v| v.is_finite()));
                    s[..40].to_vec()
                })
                .collect()
        };
        assert_eq!(fits(&x), fits(&x2));
    }

    #[test]
    fn conformal_trades_coverage_for_error() {
        // Well-separated scores with a band of overlap.
        let n = 400;
        let truth: Vec<bool> = (0..n).map(|i| i % 2 == 0).collect();
        let scores: Vec<f64> = (0..n)
            .map(|i| {
                if truth[i] {
                    0.55 + (i % 9) as f64 * 0.05
                } else {
                    0.45 - (i % 9) as f64 * 0.05 + if i % 10 == 1 { 0.5 } else { 0.0 }
                }
            })
            .collect();
        let pts = conformal(&truth, &scores, &[0.01, 0.2], 3);
        assert!(pts[0].coverage <= pts[1].coverage);
        assert!(pts[1].me <= 0.25);
    }
}
