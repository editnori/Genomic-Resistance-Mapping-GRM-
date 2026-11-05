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
//! Both outer splits keep each k-mer-identical duplicate group (`duplicates`)
//! whole in one test fold, so no held-out score in `splits` comes from a fit
//! that trained on the genome's own copy. The random folds move each group
//! by `isolate_groups`; the lineage folds hold them already, since identical
//! genomes share a lineage. Inner tuning keeps each group whole too: the
//! random split tunes on duplicate groups (`tuning_groups`), the lineage
//! split on lineages. From revision 12 the final models tune on the
//! duplicate groups as well (`final_tuning`). Only the MIC regression's own
//! folds (`evaluate_mic`) are drawn as before and may split a group. A move
//! can leave a random test fold without one class, and then the dataset
//! fails like any invalid random split (docs/duplicate-isolation.md).
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
    BinaryMetrics, NoCall, binary_metrics, fold_masks, grouped_folds, holdout, isolate_groups,
    roc_curve, screen, stratified_folds, training_candidates, training_folds, tuning_holdout,
    validate_folds,
};
use crate::bits::{BitMatrix, Bits};
use crate::curate::{Breakpoint, Call, Censor, Mic, infer_breakpoint};
use anyhow::{Context, Result, ensure};
use serde::{Deserialize, Serialize};
use std::borrow::Cow;
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

/// Solver evidence for one two-class lasso fit. `selected_point` and both
/// convergence records are present whenever the tuning path is non-empty;
/// an empty path is an intentional null fit, not a solver failure.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct LassoFitAudit {
    pub path_points: usize,
    pub path_unconverged: usize,
    pub selected_point: Option<usize>,
    pub selected_convergence: Option<lasso::Convergence>,
    pub refit_convergence: Option<lasso::Convergence>,
}

impl LassoFitAudit {
    fn empty() -> Self {
        Self {
            path_points: 0,
            path_unconverged: 0,
            selected_point: None,
            selected_convergence: None,
            refit_convergence: None,
        }
    }

    pub fn has_unconverged(&self) -> bool {
        self.path_unconverged > 0
            || self.selected_convergence.is_some_and(|c| !c.converged)
            || self.refit_convergence.is_some_and(|c| !c.converged)
    }
}

/// How every penalty point of one lasso fit stopped: the tuning path, then
/// the refit down to the selected penalty, whose last point is scored.
/// Reports keep only `LassoFitAudit`; the large-species scan keeps this.
#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize)]
pub struct LassoTrace {
    pub path: Vec<lasso::Convergence>,
    pub refit: Vec<lasso::Convergence>,
}

/// One lasso fit run by `lasso_only`.
#[derive(Clone, Debug, Serialize)]
pub struct LassoOnlyFit {
    /// "random", "lineage", or "final".
    pub scheme: String,
    /// None for the final model.
    pub outer_fold: Option<usize>,
    /// Wall-clock seconds of the fit alone, as `evaluate` times a learner.
    pub seconds: f64,
    /// The record a report would keep; None when the fit failed.
    pub audit: Option<LassoFitAudit>,
    /// The lasso error, naming the dataset, split, and fold.
    pub error: Option<String>,
    pub trace: LassoTrace,
    /// The scored refit; None for an empty path or a failed fit.
    pub model: Option<lasso::LassoModel>,
}

/// What `lasso_only`, and `analysis::lasso_scan` around it, report as they go.
#[derive(Clone, Debug, Serialize)]
#[serde(tag = "event", rename_all = "snake_case")]
pub enum LassoOnlyEvent {
    /// A dataset's evaluation inputs, which `analysis::lasso_scan` reports
    /// before the dataset's fits.
    Inputs {
        samples: usize,
        resistant: usize,
        features: usize,
        lineages: usize,
        min_support: u32,
        input_fingerprint: String,
        seconds: f64,
    },
    /// A split whose partition plan cannot support evaluation. `evaluate`
    /// records the same reason and fits nothing on it.
    Unavailable { scheme: String, reason: String },
    /// A fit is about to start.
    Start {
        scheme: String,
        outer_fold: Option<usize>,
    },
    /// A fit finished or failed.
    Fit(LassoOnlyFit),
}

/// One penalty point of a `lasso_only` fit, reported as soon as the solver
/// finishes it, so a fit that is stopped part-way shows how far it got.
#[derive(Clone, Copy, Debug)]
pub struct LassoProgress<'a> {
    /// "random", "lineage", or "final".
    pub scheme: &'a str,
    /// None for the final model.
    pub outer_fold: Option<usize>,
    /// "path" for the tuning path, "refit" for the refit.
    pub stage: &'a str,
    pub point: lasso::PointProgress,
}

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
    /// One solver-audit record per outer fold on revision-8 lasso rows.
    /// Missing on other learners and when reading older reports.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub lasso_fits: Option<Vec<LassoFitAudit>>,
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
    /// Solver evidence for the final lasso fit, including an empty-path null
    /// fit. Missing only when deserializing a pre-revision-8 report.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub lasso_fit: Option<LassoFitAudit>,
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
    lasso_fit: Option<LassoFitAudit>,
}

fn fit_scm(x: &BitMatrix, y: &Bits, train: &Bits, inner: &[usize], cfg: &EvalConfig) -> Scored {
    let c = scm::select(x, y, train, cfg.scm_max_rules, inner, cfg.inner_folds);
    let m = scm::fit(x, y, train, c.kind, c.p, c.n_rules, false);
    let pred = m.predict(x);
    Scored {
        scores: Some((0..x.cols).map(|i| pred.contains(i) as u8 as f64).collect()),
        complexity: m.rules.len() as f64,
        lasso_fit: None,
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
        lasso_fit: None,
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
            lasso_fit: None,
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

/// The share of training samples boosting holds out for early stopping.
pub(crate) const GBT_HOLDOUT: f64 = 0.15;
/// The share of training samples the lasso holds out to choose its penalty.
pub(crate) const LASSO_HOLDOUT: f64 = 0.2;

fn fit_gbt(
    x: &BitMatrix,
    y: &Bits,
    train: &Bits,
    groups: Option<&[u32]>,
    unit: &str,
    cfg: &EvalConfig,
    seed: u64,
) -> Result<(Scored, gbt::Gbt)> {
    let (fit_m, val) = tuning_holdout(train, y, groups, unit, GBT_HOLDOUT, seed)?;
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
            lasso_fit: None,
        },
        m,
    ))
}

#[derive(Clone, Copy, Debug)]
struct LassoContext<'a> {
    dataset: &'a str,
    scheme: &'a str,
    outer_fold: Option<usize>,
}

impl LassoContext<'_> {
    fn location(self) -> String {
        match self.outer_fold {
            Some(fold) => format!(
                "dataset {} split {} outer fold {fold}",
                self.dataset, self.scheme
            ),
            None => format!(
                "dataset {} final model on all samples (no outer fold)",
                self.dataset
            ),
        }
    }
}

fn lasso_error(context: LassoContext<'_>, stage: &str, detail: String) -> anyhow::Error {
    anyhow::anyhow!("lasso {stage} failure at {}: {detail}", context.location())
}

fn fit_lasso(
    x: &BitMatrix,
    y: &Bits,
    train: &Bits,
    groups: Option<&[u32]>,
    cfg: &EvalConfig,
    seed: u64,
    context: LassoContext<'_>,
) -> Result<(Scored, Option<lasso::LassoModel>)> {
    fit_lasso_traced(
        x,
        y,
        train,
        groups,
        cfg,
        seed,
        context,
        &mut LassoTrace::default(),
        &mut |_, _| {},
    )
}

/// `fit_lasso`, recording in `trace` how each point stopped as it is solved,
/// so a fit that fails still shows its points. `progress` sees each point of
/// the tuning path ("path") and the refit ("refit") as it is solved.
#[allow(clippy::too_many_arguments)]
fn fit_lasso_traced(
    x: &BitMatrix,
    y: &Bits,
    train: &Bits,
    groups: Option<&[u32]>,
    cfg: &EvalConfig,
    seed: u64,
    context: LassoContext<'_>,
    trace: &mut LassoTrace,
    progress: &mut dyn FnMut(&str, &lasso::PointProgress),
) -> Result<(Scored, Option<lasso::LassoModel>)> {
    let unit = tuning_unit(context.scheme);
    let (fit_m, val) = tuning_holdout(train, y, groups, unit, LASSO_HOLDOUT, seed)?;
    let feats = screen(x, y, &fit_m, cfg.screened_features);
    let path =
        lasso::fit_path_with_status(x, y, &fit_m, &feats, cfg.lasso_lambdas, 400, &mut |p| {
            progress("path", p)
        });
    trace.path = path.iter().map(|f| f.convergence).collect();
    fit_lasso_from_path(x, y, &val, path, context, |lambdas| {
        let refit = lasso::refit_path_with_status(x, y, train, &feats, lambdas, &mut |p| {
            progress("refit", p)
        });
        trace.refit = refit.iter().map(|f| f.convergence).collect();
        refit
    })
}

/// Select and refit one already-solved path. Keeping this separate lets the
/// negative tests inject capped statuses without changing solver limits.
fn fit_lasso_from_path(
    x: &BitMatrix,
    y: &Bits,
    val: &Bits,
    path: Vec<lasso::LassoFit>,
    context: LassoContext<'_>,
    refit: impl FnOnce(&[f64]) -> Vec<lasso::LassoFit>,
) -> Result<(Scored, Option<lasso::LassoModel>)> {
    let path_unconverged = path.iter().filter(|f| !f.convergence.converged).count();
    if path.is_empty() {
        return Ok((
            Scored {
                scores: None,
                complexity: 0.0,
                lasso_fit: Some(LassoFitAudit::empty()),
            },
            None,
        ));
    }
    let best = path
        .iter()
        .enumerate()
        .filter(|(_, f)| f.convergence.converged)
        .map(|(i, f)| (i, lasso::log_loss(&f.model.predict_proba(x), y, val)))
        .min_by(|a, b| a.1.total_cmp(&b.1))
        .map(|(i, _)| i)
        .ok_or_else(|| {
            lasso_error(
                context,
                "path",
                format!("no converged candidate among {} points", path.len()),
            )
        })?;
    let selected = path[best].convergence;
    if !selected.converged {
        return Err(lasso_error(
            context,
            "path selection",
            format!("selected point {best} is unconverged: {selected:?}"),
        ));
    }

    // Refit on the whole training fold down to the chosen lambda. Only the
    // last refit is scored or saved; its own convergence is mandatory.
    let lambdas: Vec<f64> = path[..=best].iter().map(|f| f.model.lambda).collect();
    let refit_path = refit(&lambdas);
    let final_fit = refit_path.last().ok_or_else(|| {
        lasso_error(
            context,
            "refit",
            format!("selected point {best} returned an empty refit path"),
        )
    })?;
    let final_convergence = final_fit.convergence;
    if !final_convergence.converged {
        return Err(lasso_error(
            context,
            "refit",
            format!("selected point {best} is unconverged: {final_convergence:?}"),
        ));
    }
    let model = final_fit.model.clone();
    let audit = LassoFitAudit {
        path_points: path.len(),
        path_unconverged,
        selected_point: Some(best),
        selected_convergence: Some(selected),
        refit_convergence: Some(final_convergence),
    };
    Ok((
        Scored {
            scores: lasso_scores(Some(&model), x),
            complexity: model.coef.len() as f64,
            lasso_fit: Some(audit),
        },
        Some(model),
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
    context: LassoContext<'_>,
) -> Result<Scored> {
    let fold = context.outer_fold.expect("an outer fit needs its fold");
    let seed = fold_seed(cfg, fold);
    let inner = || training_folds(y, train, groups, cfg.inner_folds, seed);
    Ok(match learner {
        "scm" => fit_scm(x, y, train, &inner(), cfg),
        "cart" => fit_cart(x, y, train, &inner(), cfg),
        "rf" => fit_rf(x, y, train, cfg, seed).0,
        "gbt" => fit_gbt(x, y, train, groups, tuning_unit(context.scheme), cfg, seed)?.0,
        "lasso" => fit_lasso(x, y, train, groups, cfg, seed, context)?.0,
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
    let lasso_fits = if learner == "lasso" {
        let audit = fits
            .iter()
            .map(|f| f.lasso_fit.clone())
            .collect::<Option<Vec<_>>>();
        ensure!(audit.is_some(), "a lasso fold is missing solver status");
        audit
    } else {
        ensure!(
            fits.iter().all(|f| f.lasso_fit.is_none()),
            "{learner} unexpectedly has lasso solver status"
        );
        None
    };
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
        lasso_fits,
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
/// in, and absent from, at least that many training samples. `duplicates`
/// lists the k-mer-identical groups `check_inputs` describes.
#[allow(clippy::too_many_arguments)]
pub fn evaluate(
    dataset: &str,
    x: &BitMatrix,
    y: &Bits,
    groups: &[u32],
    duplicates: &[Vec<usize>],
    min_support: u32,
    mic: Option<MicInput>,
    cfg: &EvalConfig,
) -> Result<EvalOutput> {
    let started = Instant::now();
    let n = x.cols;
    check_inputs(x, y, groups, duplicates, min_support, cfg)?;
    if let Some(m) = &mic {
        ensure!(m.mic.len() == n, "MIC dimensions differ");
    }
    let truth: Vec<bool> = (0..n).map(|i| y.contains(i)).collect();
    let n_groups = count_groups(groups);
    let mut splits = Vec::new();

    for (scheme, folds) in &outer_splits(y, groups, duplicates, cfg) {
        let tuning = tuning_groups(scheme, groups, duplicates);
        // A missing lineage estimate must not erase a valid random estimate.
        match check_split_plan(scheme, folds, y, tuning.as_deref(), n_groups, cfg) {
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
                fits[li].push(fit_fold(
                    learner,
                    &cx,
                    y,
                    &train,
                    tuning.as_deref(),
                    cfg,
                    LassoContext {
                        dataset,
                        scheme,
                        outer_fold: Some(fi),
                    },
                )?);
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
    let tuning = final_tuning(groups, duplicates);
    let finals = final_models(dataset, x, y, tuning.as_deref(), min_support, cfg)?;
    Ok(EvalOutput {
        splits,
        mic: mic_result,
        finals,
        lineage_threshold: None,
        seconds: started.elapsed().as_secs_f64(),
    })
}

/// What `evaluate` requires of its inputs and settings. Each duplicate
/// group lists two or more samples in ascending order, no sample is in two
/// groups, and a group's samples share a lineage, so no lineage fold splits
/// one.
fn check_inputs(
    x: &BitMatrix,
    y: &Bits,
    groups: &[u32],
    duplicates: &[Vec<usize>],
    min_support: u32,
    cfg: &EvalConfig,
) -> Result<()> {
    let n = x.cols;
    ensure!(
        n == y.len() && n == groups.len(),
        "evaluation dimensions differ"
    );
    let mut grouped = Bits::zeros(n);
    for group in duplicates {
        ensure!(
            group.len() >= 2 && group.windows(2).all(|w| w[0] < w[1]) && group[group.len() - 1] < n,
            "a duplicate group needs two or more ascending sample indices below {n}: {group:?}"
        );
        for &i in group {
            ensure!(
                !grouped.contains(i),
                "sample {i} is in two duplicate groups"
            );
            grouped.insert(i);
            ensure!(
                groups[i] == groups[group[0]],
                "duplicate samples {} and {i} are in different lineages",
                group[0]
            );
        }
    }
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
    Ok(())
}

fn count_groups(groups: &[u32]) -> usize {
    groups
        .iter()
        .collect::<std::collections::BTreeSet<_>>()
        .len()
}

/// The outer partitions `evaluate` scores, in order. The random folds move
/// each duplicate group whole into one fold; a dataset without duplicates
/// keeps its stratified folds unchanged.
fn outer_splits(
    y: &Bits,
    groups: &[u32],
    duplicates: &[Vec<usize>],
    cfg: &EvalConfig,
) -> [(&'static str, Vec<usize>); 2] {
    let random = stratified_folds(y, cfg.folds, cfg.seed);
    [
        ("random", isolate_groups(&random, duplicates, cfg.folds)),
        ("lineage", grouped_folds(y, groups, cfg.folds, cfg.seed)),
    ]
}

/// The groups that inner tuning (scm and cart's inner folds, boosting's and
/// the lasso's holdouts) keeps whole inside each outer training set:
/// - `lineage`: the lineage groups;
/// - `random` without duplicate groups: none, so tuning stays stratified;
/// - `random` with duplicate groups: each group takes its first member's
///   index as its id and every other sample its own index, so a genome and
///   its copy never sit on both sides of a tuning partition.
///
/// The first member, the lowest because `check_inputs` requires ascending
/// groups, is part of the contract, not one choice among equals. Any
/// member's index would keep groups whole, but ids set the order the groups
/// are shuffled from, so another member would redraw every random tuning
/// partition and change reported results. This function doesn't check the
/// order itself.
///
/// The outer random folds hold each group whole, so a training set holds a
/// group entirely or not at all, and ids are sample indices: no outer test
/// sample changes a training sample's id or the order of training ids.
///
/// Groups move random tuning onto the grouped path, which redraws every
/// sample's tuning partition, not only a duplicate's. Inner folds come from
/// one shuffle of the groups instead of one per class; groups are placed
/// largest first, so duplicate groups go before singletons, each into the
/// fold that least raises the squared class imbalance. The singletons then
/// even out each class across folds as stratification does. A holdout is
/// the valid one of round(1 / fraction) grouped folds closest in size to the
/// fraction, not that fraction of each class: boosting holds out about 1/7
/// of training (not 15%), the lasso about 1/5.
pub(crate) fn tuning_groups<'a>(
    scheme: &str,
    groups: &'a [u32],
    duplicates: &[Vec<usize>],
) -> Option<Cow<'a, [u32]>> {
    if scheme == "lineage" {
        return Some(Cow::Borrowed(groups));
    }
    if duplicates.is_empty() {
        return None;
    }
    let mut ids: Vec<u32> = (0..groups.len() as u32).collect();
    for members in duplicates {
        for &i in members {
            ids[i] = members[0] as u32;
        }
    }
    Some(Cow::Owned(ids))
}

/// What tuning errors call one of `scheme`'s tuning groups. The lineage
/// split keeps its revision-10 error text.
pub(crate) fn tuning_unit(scheme: &str) -> &'static str {
    if scheme == "lineage" {
        "lineage"
    } else {
        "group"
    }
}

/// Checks one split's complete nested partition plan before any model is
/// fitted: its outer folds, each outer fold's inner folds, and both tuning
/// holdouts, drawn on `tuning` groups (`tuning_groups`).
fn check_split_plan(
    scheme: &str,
    folds: &[usize],
    y: &Bits,
    tuning: Option<&[u32]>,
    n_groups: usize,
    cfg: &EvalConfig,
) -> Result<()> {
    if scheme == "lineage" {
        ensure!(
            n_groups >= cfg.folds * 2,
            "fewer than twice as many lineages as requested folds"
        );
    }
    validate_folds(y, &Bits::ones(y.len()), folds, cfg.folds)
        .with_context(|| format!("invalid {scheme} outer evaluation"))?;
    for (fi, (train, _)) in fold_masks(folds, cfg.folds).into_iter().enumerate() {
        let seed = fold_seed(cfg, fi);
        let inner = training_folds(y, &train, tuning, cfg.inner_folds, seed);
        validate_folds(y, &train, &inner, cfg.inner_folds)
            .with_context(|| format!("invalid {scheme} inner evaluation in outer fold {fi}"))?;
        for fraction in [GBT_HOLDOUT, LASSO_HOLDOUT] {
            tuning_holdout(&train, y, tuning, tuning_unit(scheme), fraction, seed)
                .with_context(|| format!("invalid {scheme} tuning holdout in outer fold {fi}"))?;
        }
    }
    Ok(())
}

/// Runs only the lasso fits of `evaluate`, in its order and with its
/// partitions, seeds, candidates, and solver: each available split's outer
/// folds, then the final model. `observe` sees each fit start and finish, so
/// a caller can time the shipped fits one at a time, as the large-species
/// scan (`scripts/lasso-scan`) does. A lasso error goes to `observe` and is
/// then returned, as `evaluate` would return it. `progress` sees each penalty
/// point as it is solved.
#[allow(clippy::too_many_arguments)]
pub fn lasso_only(
    dataset: &str,
    x: &BitMatrix,
    y: &Bits,
    groups: &[u32],
    duplicates: &[Vec<usize>],
    min_support: u32,
    cfg: &EvalConfig,
    observe: &mut dyn FnMut(&LassoOnlyEvent) -> Result<()>,
    progress: &mut dyn FnMut(&LassoProgress),
) -> Result<()> {
    check_inputs(x, y, groups, duplicates, min_support, cfg)?;
    let n_groups = count_groups(groups);
    for (scheme, folds) in &outer_splits(y, groups, duplicates, cfg) {
        let tuning = tuning_groups(scheme, groups, duplicates);
        match check_split_plan(scheme, folds, y, tuning.as_deref(), n_groups, cfg) {
            Ok(()) => {}
            Err(error) if *scheme == "lineage" => {
                observe(&LassoOnlyEvent::Unavailable {
                    scheme: scheme.to_string(),
                    reason: format!("{error:#}"),
                })?;
                continue;
            }
            Err(error) => return Err(error),
        }
        for (fi, (train, _)) in fold_masks(folds, cfg.folds).into_iter().enumerate() {
            let (cx, _) = training_candidates(x, &train, min_support);
            let context = LassoContext {
                dataset,
                scheme,
                outer_fold: Some(fi),
            };
            observed_lasso(scheme, Some(fi), observe, |trace| {
                let seed = fold_seed(cfg, fi);
                fit_lasso_traced(
                    &cx,
                    y,
                    &train,
                    tuning.as_deref(),
                    cfg,
                    seed,
                    context,
                    trace,
                    &mut |stage, point| {
                        progress(&LassoProgress {
                            scheme,
                            outer_fold: Some(fi),
                            stage,
                            point: *point,
                        })
                    },
                )
            })?;
        }
    }
    let tuning = final_tuning(groups, duplicates);
    final_plan(x, y, tuning.as_deref(), min_support, cfg)?;
    observed_lasso("final", None, observe, |trace| {
        final_lasso(
            dataset,
            x,
            y,
            tuning.as_deref(),
            cfg,
            trace,
            &mut |stage, point| {
                progress(&LassoProgress {
                    scheme: "final",
                    outer_fold: None,
                    stage,
                    point: *point,
                })
            },
        )
    })
}

/// Runs one lasso fit between its `Start` and `Fit` events.
fn observed_lasso(
    scheme: &str,
    outer_fold: Option<usize>,
    observe: &mut dyn FnMut(&LassoOnlyEvent) -> Result<()>,
    fit: impl FnOnce(&mut LassoTrace) -> Result<(Scored, Option<lasso::LassoModel>)>,
) -> Result<()> {
    observe(&LassoOnlyEvent::Start {
        scheme: scheme.to_string(),
        outer_fold,
    })?;
    let mut trace = LassoTrace::default();
    let started = Instant::now();
    let result = fit(&mut trace);
    let seconds = started.elapsed().as_secs_f64();
    let (audit, model, error) = match result {
        Ok((scored, model)) => (scored.lasso_fit, model, None),
        Err(error) => (None, None, Some(error)),
    };
    observe(&LassoOnlyEvent::Fit(LassoOnlyFit {
        scheme: scheme.to_string(),
        outer_fold,
        seconds,
        audit,
        error: error.as_ref().map(|e| format!("{e:#}")),
        trace,
        model,
    }))?;
    error.map_or(Ok(()), Err)
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

/// The groups the final models tune on (scm and cart's inner folds,
/// boosting's and the lasso's holdouts): the duplicate groups, as the random
/// split's inner tuning draws them (`tuning_groups`), so a genome and its
/// copy never sit on both sides of a final tuning partition. Lineages play
/// no part. Without duplicate groups it is None, and the final models tune
/// by class as revision 11's did.
fn final_tuning<'a>(groups: &'a [u32], duplicates: &[Vec<usize>]) -> Option<Cow<'a, [u32]>> {
    tuning_groups("random", groups, duplicates)
}

fn final_models(
    dataset: &str,
    x: &BitMatrix,
    y: &Bits,
    tuning: Option<&[u32]>,
    min_support: u32,
    cfg: &EvalConfig,
) -> Result<FinalModels> {
    let all = Bits::ones(x.cols);
    let inner = final_plan(x, y, tuning, min_support, cfg)?;
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
    let (lasso_scored, lasso_model) = final_lasso(
        dataset,
        x,
        y,
        tuning,
        cfg,
        &mut LassoTrace::default(),
        &mut |_, _| {},
    )?;
    let (_, rf_importance) = fit_rf(x, y, &all, cfg, cfg.seed + 17);
    let (_, gbt_model) = fit_gbt(x, y, &all, tuning, "group", cfg, cfg.seed + 13)?;
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
        lasso_fit: lasso_scored.lasso_fit,
        rf_importance: top_weighted(rf_importance, 50),
        gbt_importance: top_weighted(gbt_model.importance(), 50),
    })
}

/// The final models' preconditions, and their inner folds, drawn on
/// `tuning` (`final_tuning`). They train on every sample, so their
/// candidates are the dataset matrix itself; that keeps their feature ids
/// valid for annotation. Without groups the folds are revision 11's
/// stratified ones: `training_folds` over every sample is
/// `stratified_folds`.
fn final_plan(
    x: &BitMatrix,
    y: &Bits,
    tuning: Option<&[u32]>,
    min_support: u32,
    cfg: &EvalConfig,
) -> Result<Vec<usize>> {
    let all = Bits::ones(x.cols);
    ensure!(
        training_candidates(x, &all, min_support).1.len() == x.rows,
        "final models need a matrix filtered and collapsed on all samples at min_support {min_support}"
    );
    let inner = training_folds(y, &all, tuning, cfg.inner_folds, cfg.seed + 7);
    validate_folds(y, &all, &inner, cfg.inner_folds)
        .context("invalid final-model inner evaluation")?;
    Ok(inner)
}

/// The final lasso: every sample, its holdout drawn on `tuning`
/// (`final_tuning`), and its own seed.
fn final_lasso(
    dataset: &str,
    x: &BitMatrix,
    y: &Bits,
    tuning: Option<&[u32]>,
    cfg: &EvalConfig,
    trace: &mut LassoTrace,
    progress: &mut dyn FnMut(&str, &lasso::PointProgress),
) -> Result<(Scored, Option<lasso::LassoModel>)> {
    fit_lasso_traced(
        x,
        y,
        &Bits::ones(x.cols),
        tuning,
        cfg,
        cfg.seed + 11,
        LassoContext {
            dataset,
            scheme: "final",
            outer_fold: None,
        },
        trace,
        progress,
    )
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

    fn test_context<'a>(scheme: &'a str, outer_fold: usize) -> LassoContext<'a> {
        LassoContext {
            dataset: "fixture__fixturemycin",
            scheme,
            outer_fold: Some(outer_fold),
        }
    }

    fn canned_lasso(lambda: f64, intercept: f64, weight: f64, converged: bool) -> lasso::LassoFit {
        lasso::LassoFit {
            model: lasso::LassoModel {
                intercept,
                coef: vec![(0, weight)],
                lambda,
            },
            convergence: lasso::Convergence {
                irls_steps: if converged { 3 } else { 25 },
                sweep_cap_hits: u32::from(!converged),
                converged,
            },
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
        assert!(evaluate("fixture__fixturemycin", &x, &y, &groups, &[], 3, None, &cfg).is_err());
        cfg.folds = 2;
        let mut duplicated = x.clone();
        duplicated.push_row(y.words());
        assert!(
            evaluate(
                "fixture__fixturemycin",
                &duplicated,
                &y,
                &groups,
                &[],
                3,
                None,
                &cfg
            )
            .is_err()
        );
        let out = evaluate("fixture__fixturemycin", &x, &y, &groups, &[], 3, None, &cfg).unwrap();
        assert_eq!(out.splits.len(), 2);
        assert!(!out.finals.lasso_fit.as_ref().unwrap().has_unconverged());
        for split in out.splits {
            assert!(split.unavailable_reason.is_none());
            for learner in split.learners {
                assert_eq!(learner.scores.len(), n);
                assert!(learner.scores.iter().all(|s| s.is_finite()));
                if learner.learner == "lasso" {
                    let fits = learner.lasso_fits.as_ref().unwrap();
                    assert_eq!(fits.len(), cfg.folds);
                    assert!(fits.iter().all(|f| !f.has_unconverged()));
                } else {
                    assert!(learner.lasso_fits.is_none());
                }
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
        let out = evaluate(
            "fixture__fixturemycin",
            &x,
            &y,
            &groups,
            &[],
            3,
            None,
            &small_config(),
        )
        .unwrap();
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

    /// A split's results with fit timings removed, as JSON bytes.
    fn timeless(split: &SplitResult) -> Vec<u8> {
        let mut v = serde_json::to_value(split).unwrap();
        for learner in v["learners"].as_array_mut().unwrap() {
            learner.as_object_mut().unwrap().remove("seconds");
        }
        serde_json::to_vec(&v).unwrap()
    }

    #[test]
    fn duplicate_groups_sit_whole_in_one_outer_fold_of_every_split() {
        // Lineages of four; samples 4j and 4j + 1 are copies, one R and
        // one S. Every split copy pair is a tie, so it joins the fold of its
        // lower-index member.
        let n = 80;
        let y = Bits::from_indices(n, (0..n).filter(|i| i % 2 == 0));
        let mut x = BitMatrix::with_cols(n);
        x.push_row(y.words());
        let groups: Vec<u32> = (0..n).map(|i| (i / 4) as u32).collect();
        let duplicates: Vec<Vec<usize>> = (0..n / 4).map(|j| vec![4 * j, 4 * j + 1]).collect();
        let mut cfg = small_config();
        cfg.folds = 4;
        let out = evaluate(
            "fixture__fixturemycin",
            &x,
            &y,
            &groups,
            &duplicates,
            3,
            None,
            &cfg,
        )
        .unwrap();
        let base = stratified_folds(&y, cfg.folds, cfg.seed);
        let random = &out.splits[0].test_fold;
        assert_eq!(*random, isolate_groups(&base, &duplicates, cfg.folds));
        let split_pairs: Vec<&Vec<usize>> = duplicates
            .iter()
            .filter(|g| base[g[0]] != base[g[1]])
            .collect();
        assert!(!split_pairs.is_empty());
        for g in &split_pairs {
            assert_eq!((random[g[0]], random[g[1]]), (base[g[0]], base[g[0]]));
        }
        assert!(split_pairs.iter().any(|g| base[g[0]] > base[g[1]]));
        assert_eq!(
            out.splits[1].test_fold,
            grouped_folds(&y, &groups, cfg.folds, cfg.seed)
        );
        for split in &out.splits {
            assert!(split.unavailable_reason.is_none());
            assert!(
                duplicates
                    .iter()
                    .all(|g| split.test_fold[g[0]] == split.test_fold[g[1]])
            );
        }
    }

    #[test]
    fn duplicate_groups_already_whole_in_their_folds_change_only_random_tuning() {
        let n = 80;
        let y = Bits::from_indices(n, (0..n).filter(|i| i % 2 == 0));
        let mut x = BitMatrix::with_cols(n);
        x.push_row(y.words());
        x.push_row(Bits::from_indices(n, (0..n).filter(|i| i % 3 == 0)).words());
        let groups: Vec<u32> = (0..n).map(|i| (i / 4) as u32).collect();
        let cfg = small_config();
        let base = stratified_folds(&y, cfg.folds, cfg.seed);
        // Same-lineage pairs that stratified folds already keep together.
        let duplicates: Vec<Vec<usize>> = (0..n / 4)
            .flat_map(|j| [(4 * j, 4 * j + 2), (4 * j + 1, 4 * j + 3)])
            .filter(|&(a, b)| base[a] == base[b])
            .map(|(a, b)| vec![a, b])
            .collect();
        assert!(!duplicates.is_empty());
        let plain = evaluate("fixture__fixturemycin", &x, &y, &groups, &[], 3, None, &cfg).unwrap();
        let grouped = evaluate(
            "fixture__fixturemycin",
            &x,
            &y,
            &groups,
            &duplicates,
            3,
            None,
            &cfg,
        )
        .unwrap();
        assert_eq!(plain.splits[0].test_fold, base);
        // No sample moves between outer folds. Random tuning now keeps the
        // pairs whole, so only the random split's fits may change.
        for (a, b) in plain.splits.iter().zip(&grouped.splits) {
            assert_eq!(a.test_fold, b.test_fold);
        }
        assert_eq!(timeless(&plain.splits[1]), timeless(&grouped.splits[1]));
        assert_eq!(
            serde_json::to_vec(&plain.finals).unwrap(),
            serde_json::to_vec(&grouped.finals).unwrap()
        );
    }

    #[test]
    fn a_move_that_empties_a_class_from_a_random_fold_fails_the_dataset() {
        // Ten resistant genomes in 80, two per stratified fold. Each resistant
        // genome of one fold is paired with a lower-index resistant genome in
        // another fold, so every tie sends it there and the fold keeps none.
        // A known limit: min_class >= folds no longer guarantees a valid
        // random split once groups move.
        let n = 80;
        let y = Bits::from_indices(n, 0..10);
        let mut x = BitMatrix::with_cols(n);
        x.push_row(y.words());
        x.push_row(Bits::from_indices(n, (0..n).filter(|i| i % 3 == 0)).words());
        let mut groups: Vec<u32> = (0..n).map(|i| (i / 4) as u32).collect();
        let mut cfg = small_config();
        cfg.folds = 5;
        let base = stratified_folds(&y, cfg.folds, cfg.seed);
        let duplicates = (0..cfg.folds)
            .find_map(|t| {
                let mut used = vec![];
                let mut pairs = vec![];
                for r in (0..10).filter(|&i| base[i] == t) {
                    let p = (0..r).find(|&p| base[p] != t && !used.contains(&p))?;
                    used.push(p);
                    pairs.push(vec![p, r]);
                }
                Some(pairs)
            })
            .unwrap();
        for g in &duplicates {
            groups[g[1]] = groups[g[0]];
        }
        let plain = evaluate("fixture__fixturemycin", &x, &y, &groups, &[], 3, None, &cfg);
        assert!(plain.is_ok(), "{:#}", plain.err().unwrap());
        let error = evaluate(
            "fixture__fixturemycin",
            &x,
            &y,
            &groups,
            &duplicates,
            3,
            None,
            &cfg,
        )
        .err()
        .unwrap();
        assert!(
            format!("{error:#}").contains("invalid random outer evaluation"),
            "{error:#}"
        );
    }

    #[test]
    fn malformed_duplicate_groups_are_rejected_before_any_fit() {
        let n = 40;
        let y = Bits::from_indices(n, (0..n).filter(|i| i % 2 == 0));
        let mut x = BitMatrix::with_cols(n);
        x.push_row(y.words());
        let groups: Vec<u32> = (0..n).map(|i| (i / 4) as u32).collect();
        let cases: [(Vec<Vec<usize>>, &str); 5] = [
            (vec![vec![3]], "two or more"),
            (vec![vec![2, 1]], "ascending"),
            (vec![vec![38, 40]], "below 40"),
            (
                vec![vec![0, 1], vec![1, 2]],
                "sample 1 is in two duplicate groups",
            ),
            (
                vec![vec![2, 4]],
                "duplicate samples 2 and 4 are in different lineages",
            ),
        ];
        for (duplicates, expected) in cases {
            let error = evaluate(
                "fixture__fixturemycin",
                &x,
                &y,
                &groups,
                &duplicates,
                3,
                None,
                &small_config(),
            )
            .err()
            .unwrap();
            assert!(format!("{error:#}").contains(expected), "{error:#}");
        }
    }

    /// 120 genomes in lineages of four, 48 resistant. Each feature row is a
    /// bit of a genome id, so samples have identical columns exactly when
    /// they share an id: ten same-lineage pairs and ten triples, some of
    /// them with both labels.
    fn duplicate_fixture() -> (BitMatrix, Bits, Vec<u32>) {
        let n = 120;
        let genome = |i: usize| match ((i / 4) % 3, i % 4) {
            (0, 1) => i - 1,
            (1, 2 | 3) => i / 4 * 4 + 1,
            _ => i,
        };
        let mut x = BitMatrix::with_cols(n);
        for bit in 0..7 {
            x.push_row(
                Bits::from_indices(n, (0..n).filter(|&i| (genome(i) >> bit) & 1 == 1)).words(),
            );
        }
        let y = Bits::from_indices(n, (0..n).filter(|i| i % 5 < 2));
        (x, y, (0..n).map(|i| (i / 4) as u32).collect())
    }

    /// The groups of samples whose columns of `x` are identical, as the
    /// report's `kset_duplicates` lists them.
    fn identical_columns(x: &BitMatrix) -> Vec<Vec<usize>> {
        let mut by_column: std::collections::BTreeMap<Vec<bool>, Vec<usize>> = Default::default();
        for i in 0..x.cols {
            let column = (0..x.rows).map(|r| x.get(r, i)).collect();
            by_column.entry(column).or_default().push(i);
        }
        let mut groups: Vec<Vec<usize>> = by_column.into_values().filter(|g| g.len() > 1).collect();
        groups.sort();
        groups
    }

    type InnerPlan = (Vec<usize>, Vec<(Bits, Bits)>);

    /// One outer fold's tuning partitions, drawn as `check_split_plan` and
    /// the fits draw them: scm and cart's inner folds, then the boosting and
    /// lasso holdouts as (fit, validation).
    fn inner_plan(
        y: &Bits,
        train: &Bits,
        tuning: Option<&[u32]>,
        cfg: &EvalConfig,
        fold: usize,
    ) -> InnerPlan {
        let seed = fold_seed(cfg, fold);
        let holdouts = [GBT_HOLDOUT, LASSO_HOLDOUT]
            .map(|frac| tuning_holdout(train, y, tuning, "group", frac, seed).unwrap())
            .to_vec();
        (
            training_folds(y, train, tuning, cfg.inner_folds, seed),
            holdouts,
        )
    }

    #[test]
    fn random_tuning_keeps_each_duplicate_group_whole_inside_outer_training() {
        // A group takes its lowest member's index, every other sample its own.
        assert_eq!(
            tuning_groups("random", &[7; 6], &[vec![1, 3], vec![2, 4, 5]]).as_deref(),
            Some([0, 1, 2, 1, 2, 2].as_slice())
        );
        let (x, y, groups) = duplicate_fixture();
        let duplicates = identical_columns(&x);
        assert_eq!(duplicates.len(), 20);
        let cfg = EvalConfig::default();
        let [(_, random), (_, lineage)] = outer_splits(&y, &groups, &duplicates, &cfg);
        // The outer folds are revision 10's.
        let stratified = stratified_folds(&y, cfg.folds, cfg.seed);
        assert_eq!(random, isolate_groups(&stratified, &duplicates, cfg.folds));
        assert_eq!(lineage, grouped_folds(&y, &groups, cfg.folds, cfg.seed));
        let tuning = tuning_groups("random", &groups, &duplicates);
        check_split_plan("random", &random, &y, tuning.as_deref(), 30, &cfg).unwrap();
        let mut split_without_groups = 0;
        for (fi, (train, _)) in fold_masks(&random, cfg.folds).into_iter().enumerate() {
            assert!(
                duplicates
                    .iter()
                    .all(|g| g.iter().all(|&i| train.contains(i) == train.contains(g[0])))
            );
            let inside: Vec<&Vec<usize>> =
                duplicates.iter().filter(|g| train.contains(g[0])).collect();
            assert!(!inside.is_empty());
            // Groups that straddle an inner fold or either holdout.
            let straddling = |(inner, holdouts): &InnerPlan| {
                inside
                    .iter()
                    .filter(|g| {
                        g.iter().any(|&i| {
                            inner[i] != inner[g[0]]
                                || holdouts
                                    .iter()
                                    .any(|(_, val)| val.contains(i) != val.contains(g[0]))
                        })
                    })
                    .count()
            };
            let plan = inner_plan(&y, &train, tuning.as_deref(), &cfg, fi);
            assert_eq!(straddling(&plan), 0, "outer fold {fi}");
            split_without_groups += straddling(&inner_plan(&y, &train, None, &cfg, fi));
            // The grouped path stays as balanced as stratification: each
            // class is within one sample of an even share of every inner
            // fold and of each holdout's round(1 / frac) folds.
            let (inner, holdouts) = &plan;
            for class in [true, false] {
                let members: Vec<usize> = train
                    .iter_ones()
                    .filter(|&i| y.contains(i) == class)
                    .collect();
                let count =
                    |m: &dyn Fn(usize) -> bool| members.iter().filter(|&&i| m(i)).count() as f64;
                let even = members.len() as f64 / cfg.inner_folds as f64;
                for f in 0..cfg.inner_folds {
                    assert!(
                        (count(&|i| inner[i] == f) - even).abs() < 1.0,
                        "outer fold {fi}"
                    );
                }
                for ((_, val), frac) in holdouts.iter().zip([GBT_HOLDOUT, LASSO_HOLDOUT]) {
                    let even = members.len() as f64 / (1.0 / frac).round();
                    assert!(
                        (count(&|i| val.contains(i)) - even).abs() < 1.0,
                        "outer fold {fi}"
                    );
                }
            }
        }
        // Revision 10's tuning split copies in this fixture.
        assert!(split_without_groups > 0);
    }

    #[test]
    fn outer_test_samples_cannot_change_inner_tuning() {
        // Inner plans read labels, lineages and duplicate groups, which come
        // from features. Reverse the outer test samples' labels among them
        // and, separately, rotate their features and lineages, which moves
        // their duplicate groups. No training sample's plan may change.
        let (x, y, groups) = duplicate_fixture();
        let n = x.cols;
        let duplicates = identical_columns(&x);
        let cfg = EvalConfig::default();
        for (scheme, folds) in outer_splits(&y, &groups, &duplicates, &cfg) {
            let tuning = tuning_groups(scheme, &groups, &duplicates);
            for (fi, (train, test)) in fold_masks(&folds, cfg.folds).into_iter().enumerate() {
                let t: Vec<usize> = test.iter_ones().collect();
                let m = t.len();
                let mut label_from: Vec<usize> = (0..n).collect();
                let mut genome_from = label_from.clone();
                for (j, &i) in t.iter().enumerate() {
                    label_from[i] = t[m - 1 - j];
                    genome_from[t[(j + 1) % m]] = i;
                }
                let y2 = Bits::from_indices(n, (0..n).filter(|&i| y.contains(label_from[i])));
                let groups2: Vec<u32> = genome_from.iter().map(|&s| groups[s]).collect();
                let mut x2 = BitMatrix::with_cols(n);
                for r in 0..x.rows {
                    x2.push_row(
                        Bits::from_indices(n, (0..n).filter(|&i| x.get(r, genome_from[i]))).words(),
                    );
                }
                let duplicates2 = identical_columns(&x2);
                assert!(
                    y2 != y && duplicates2 != duplicates,
                    "{scheme} outer fold {fi}"
                );
                check_inputs(&x2, &y2, &groups2, &duplicates2, 1, &cfg).unwrap();
                let tuning2 = tuning_groups(scheme, &groups2, &duplicates2);
                assert_eq!(
                    inner_plan(&y2, &train, tuning2.as_deref(), &cfg, fi),
                    inner_plan(&y, &train, tuning.as_deref(), &cfg, fi),
                    "{scheme} outer fold {fi}"
                );
            }
        }
    }

    #[test]
    fn random_tuning_without_duplicates_is_revision_10s_stratified_tuning() {
        let (_, y, groups) = duplicate_fixture();
        let cfg = EvalConfig::default();
        let tuning = tuning_groups("random", &groups, &[]);
        assert!(tuning.is_none());
        let [(_, random), _] = outer_splits(&y, &groups, &[], &cfg);
        assert_eq!(random, stratified_folds(&y, cfg.folds, cfg.seed));
        for (fi, (train, _)) in fold_masks(&random, cfg.folds).into_iter().enumerate() {
            // Stratified inner folds over the training samples, and `holdout`.
            let seed = fold_seed(&cfg, fi);
            let members: Vec<usize> = train.iter_ones().collect();
            let local_y = Bits::from_indices(
                members.len(),
                (0..members.len()).filter(|&j| y.contains(members[j])),
            );
            let mut inner = vec![usize::MAX; y.len()];
            for (&i, f) in members
                .iter()
                .zip(stratified_folds(&local_y, cfg.inner_folds, seed))
            {
                inner[i] = f;
            }
            let holdouts = [GBT_HOLDOUT, LASSO_HOLDOUT].map(|frac| holdout(&train, &y, frac, seed));
            assert_eq!(
                inner_plan(&y, &train, tuning.as_deref(), &cfg, fi),
                (inner, holdouts.to_vec())
            );
        }
    }

    /// The final models' tuning partitions over every sample, drawn as
    /// `final_models` draws them: scm and cart's inner folds, then the
    /// boosting and lasso holdouts with their own seeds.
    fn final_tuning_plan(x: &BitMatrix, y: &Bits, tuning: Option<&[u32]>) -> InnerPlan {
        let cfg = EvalConfig::default();
        let all = Bits::ones(y.len());
        let holdouts = [(GBT_HOLDOUT, cfg.seed + 13), (LASSO_HOLDOUT, cfg.seed + 11)]
            .map(|(frac, seed)| tuning_holdout(&all, y, tuning, "group", frac, seed).unwrap())
            .to_vec();
        (final_plan(x, y, tuning, 1, &cfg).unwrap(), holdouts)
    }

    #[test]
    fn final_tuning_keeps_each_duplicate_group_whole() {
        let (x, y, groups) = duplicate_fixture();
        let duplicates = identical_columns(&x);
        assert_eq!(duplicates.len(), 20);
        // Groups that straddle a final inner fold or either final holdout.
        let straddling = |(inner, holdouts): &InnerPlan| {
            duplicates
                .iter()
                .filter(|g| {
                    g.iter().any(|&i| inner[i] != inner[g[0]])
                        || holdouts.iter().any(|(_, val)| {
                            g.iter().any(|&i| val.contains(i) != val.contains(g[0]))
                        })
                })
                .count()
        };
        let tuning = final_tuning(&groups, &duplicates);
        assert_eq!(tuning, tuning_groups("random", &groups, &duplicates));
        assert_eq!(straddling(&final_tuning_plan(&x, &y, tuning.as_deref())), 0);
        // Revision 11's final tuning, by class, split some of them.
        assert!(straddling(&final_tuning_plan(&x, &y, None)) > 0);
        // Without duplicate groups the final folds are revision 11's.
        assert_eq!(final_tuning(&groups, &[]), None);
        let cfg = EvalConfig::default();
        assert_eq!(
            final_tuning_plan(&x, &y, None).0,
            stratified_folds(&y, cfg.inner_folds, cfg.seed + 7)
        );
    }

    #[test]
    fn final_models_tune_on_the_duplicate_groups() {
        // evaluate's final models are the ones fitted on final_tuning's
        // partitions: scm's choice and cart's depth come from the grouped
        // inner folds, and boosting from the grouped holdout.
        let (x, y, groups) = duplicate_fixture();
        let duplicates = identical_columns(&x);
        let cfg = EvalConfig::default();
        let all = Bits::ones(y.len());
        let tuning = final_tuning(&groups, &duplicates);
        let (inner, _) = final_tuning_plan(&x, &y, tuning.as_deref());
        let out = evaluate(
            "fixture__fixturemycin",
            &x,
            &y,
            &groups,
            &duplicates,
            1,
            None,
            &cfg,
        )
        .unwrap();
        let choice = scm::select(&x, &y, &all, cfg.scm_max_rules, &inner, cfg.inner_folds);
        assert_eq!(
            (
                out.finals.scm.kind,
                out.finals.scm.p,
                out.finals.scm.cv_error
            ),
            (choice.kind, choice.p, choice.cv_error)
        );
        let (depth, _) =
            tree::select_depth(&x, &y, &all, &inner, cfg.inner_folds, cfg.cart_max_depth);
        assert_eq!(out.finals.cart_depth, depth);
        let (_, gbt) = fit_gbt(
            &x,
            &y,
            &all,
            tuning.as_deref(),
            "group",
            &cfg,
            cfg.seed + 13,
        )
        .unwrap();
        assert_eq!(
            serde_json::to_value(&out.finals.gbt_importance).unwrap(),
            serde_json::to_value(top_weighted(gbt.importance(), 50)).unwrap()
        );
        let direct =
            final_models("fixture__fixturemycin", &x, &y, tuning.as_deref(), 1, &cfg).unwrap();
        assert_eq!(
            serde_json::to_value(&out.finals).unwrap(),
            serde_json::to_value(&direct).unwrap()
        );
    }

    #[test]
    fn lineage_tuning_ignores_duplicate_groups() {
        let (x, y, groups) = duplicate_fixture();
        let duplicates = identical_columns(&x);
        let cfg = EvalConfig::default();
        for listed in [duplicates, vec![]] {
            let tuning = tuning_groups("lineage", &groups, &listed);
            assert_eq!(tuning.as_deref(), Some(groups.as_slice()));
            let [_, (_, lineage)] = outer_splits(&y, &groups, &listed, &cfg);
            assert_eq!(lineage, grouped_folds(&y, &groups, cfg.folds, cfg.seed));
            check_split_plan("lineage", &lineage, &y, tuning.as_deref(), 30, &cfg).unwrap();
        }
    }

    #[test]
    fn random_tuning_that_cannot_keep_duplicates_whole_fails_before_any_fit() {
        // Four resistant genomes in 80, two per stratified fold, and each
        // fold's pair is one duplicate group. The outer folds stay valid,
        // but every training set's resistant genomes are then one group, so
        // no grouped inner fold split holds them on both sides. Stratified
        // tuning, with no duplicates listed, splits the pair and passes.
        let n = 80;
        let y = Bits::from_indices(n, 0..4);
        let mut x = BitMatrix::with_cols(n);
        x.push_row(y.words());
        x.push_row(Bits::from_indices(n, (0..n).filter(|i| i % 3 == 0)).words());
        let mut groups: Vec<u32> = (0..n).map(|i| (i / 4) as u32).collect();
        let cfg = small_config();
        let base = stratified_folds(&y, cfg.folds, cfg.seed);
        let duplicates: Vec<Vec<usize>> = (0..cfg.folds)
            .map(|f| (0..4).filter(|&i| base[i] == f).collect())
            .collect();
        for g in &duplicates {
            groups[g[1]] = groups[g[0]];
        }
        assert_eq!(outer_splits(&y, &groups, &duplicates, &cfg)[0].1, base);
        let dataset = "fixture__fixturemycin";
        let plain = evaluate(dataset, &x, &y, &groups, &[], 3, None, &cfg);
        assert!(plain.is_ok(), "{:#}", plain.err().unwrap());
        let expected = "invalid random inner evaluation in outer fold 0";
        let error = evaluate(dataset, &x, &y, &groups, &duplicates, 3, None, &cfg).unwrap_err();
        assert!(format!("{error:#}").contains(expected), "{error:#}");
        let mut events = Vec::new();
        let error = lasso_only(
            dataset,
            &x,
            &y,
            &groups,
            &duplicates,
            3,
            &cfg,
            &mut |e| {
                events.push(e.clone());
                Ok(())
            },
            &mut |_| {},
        )
        .unwrap_err();
        assert!(format!("{error:#}").contains(expected), "{error:#}");
        assert!(events.is_empty(), "{events:?}");
    }

    #[test]
    fn a_grouped_holdout_that_fails_after_its_inner_folds_pass_stops_the_plan() {
        // Outer fold 0 tests two resistant and two susceptible genomes. It
        // trains on susceptible groups of 100, a group of 98 resistant
        // genomes and one susceptible, a group of one resistant and 97
        // susceptible genomes, and three susceptible singletons. Two inner
        // folds part the two resistant groups. A grouped holdout with one
        // fold more than there are groups of 100 gives each of them a fold
        // and puts both resistant groups in the last, where they add the
        // least susceptible imbalance, so no fold leaves resistant genomes on
        // both sides; with any other fold count here they part. Six groups
        // of 100 fail boosting's seven folds alone, and four fail the
        // lasso's five alone. A stratified holdout has no such problem. The
        // groups serve as lineages and as duplicate groups alike, so each
        // scheme fails in its own words.
        for (big, failing) in [(6, GBT_HOLDOUT), (4, LASSO_HOLDOUT)] {
            let mut units: Vec<(usize, usize)> = vec![(98, 1), (1, 97)];
            units.extend(std::iter::repeat_n((0, 100), big));
            let mut labels = vec![true, true, false, false];
            let mut duplicates = Vec::new();
            for (r, s) in units {
                let start = labels.len();
                labels.extend((0..r + s).map(|j| j < r));
                duplicates.push((start..labels.len()).collect::<Vec<usize>>());
            }
            labels.extend([false; 3]);
            let n = labels.len();
            let y = Bits::from_indices(n, (0..n).filter(|&i| labels[i]));
            let mut groups: Vec<u32> = (0..n as u32).collect();
            for g in &duplicates {
                for &i in g {
                    groups[i] = g[0] as u32;
                }
            }
            let folds: Vec<usize> = (0..n).map(|i| usize::from(i >= 4)).collect();
            let cfg = small_config();
            let (train, _) = fold_masks(&folds, cfg.folds).remove(0);
            let seed = fold_seed(&cfg, 0);
            for frac in [GBT_HOLDOUT, LASSO_HOLDOUT] {
                assert!(tuning_holdout(&train, &y, None, "group", frac, seed).is_ok());
            }
            for (scheme, listed, reason) in [
                ("random", duplicates, "no group tuning partition"),
                // Revision 10's text, byte for byte.
                ("lineage", vec![], "no lineage tuning partition"),
            ] {
                let tuning = tuning_groups(scheme, &groups, &listed);
                let inner = training_folds(&y, &train, tuning.as_deref(), cfg.inner_folds, seed);
                validate_folds(&y, &train, &inner, cfg.inner_folds).unwrap();
                for frac in [GBT_HOLDOUT, LASSO_HOLDOUT] {
                    let holdout =
                        tuning_holdout(&train, &y, tuning.as_deref(), "group", frac, seed);
                    assert_eq!(
                        holdout.is_err(),
                        frac == failing,
                        "{big} groups of 100, {scheme}, holdout {frac}"
                    );
                }
                let error = check_split_plan(
                    scheme,
                    &folds,
                    &y,
                    tuning.as_deref(),
                    count_groups(&groups),
                    &cfg,
                )
                .unwrap_err();
                assert_eq!(
                    format!("{error:#}"),
                    format!(
                        "invalid {scheme} tuning holdout in outer fold 0: \
                         {reason} contains both classes"
                    )
                );
            }
        }
    }

    #[test]
    fn two_training_groups_are_enough_for_a_grouped_holdout() {
        // I1: k is 2 here, and each lineage is one side with both classes.
        let y = Bits::from_indices(8, [0, 4]);
        let groups = [1u32, 1, 1, 1, 2, 2, 2, 2];
        let (fit, val) = tuning_holdout(
            &Bits::ones(8),
            &y,
            Some(&groups),
            "lineage",
            LASSO_HOLDOUT,
            42,
        )
        .unwrap();
        assert_eq!(fit.count() + val.count(), 8);
        for side in [&fit, &val] {
            let members: Vec<usize> = side.iter_ones().collect();
            assert!(
                members == [0, 1, 2, 3] || members == [4, 5, 6, 7],
                "{members:?}"
            );
        }
    }

    #[test]
    fn tuning_holdout_names_lineages_in_revision_10s_words() {
        // One training group. The lineage text is revision 10's, byte for
        // byte; duplicate groups get their own word.
        let y = Bits::from_indices(8, [0, 1]);
        let one_group = [3u32; 8];
        for (unit, text) in [
            ("lineage", "tuning needs at least two training lineages"),
            ("group", "tuning needs at least two training groups"),
        ] {
            let ids = Some(one_group.as_slice());
            let error =
                tuning_holdout(&Bits::ones(8), &y, ids, unit, LASSO_HOLDOUT, 42).unwrap_err();
            assert_eq!(format!("{error:#}"), text);
        }
    }

    #[test]
    fn scm_cart_boosting_and_lasso_tune_on_the_splits_tuning_groups() {
        // 160 genomes in lineages of four; every other lineage is four
        // copies of one genome. A marker calls resistance, but seven copy
        // groups carry it and are susceptible, each with a feature of its
        // own, and eight rows spell each genome's id. A fit can learn a
        // copy's label from its id or private feature only when tuning
        // splits the copies. Each fit must equal one built here from inner
        // folds and holdouts drawn on the split's tuning groups, and some fit
        // of each learner must differ from one built without them, so a fit
        // that drops the groups cannot pass.
        let n = 160;
        let copies = |i: usize| (i / 4).is_multiple_of(2);
        let genome = |i: usize| if copies(i) { i / 4 * 4 } else { i };
        let marker = |i: usize| copies(i) || i.is_multiple_of(2);
        let mislabeled = |i: usize| copies(i) && (i / 4).is_multiple_of(6);
        let y = Bits::from_indices(n, (0..n).filter(|&i| marker(i) && !mislabeled(i)));
        let mut x = BitMatrix::with_cols(n);
        x.push_row(Bits::from_indices(n, (0..n).filter(|&i| marker(i))).words());
        for bit in 0..8 {
            x.push_row(
                Bits::from_indices(n, (0..n).filter(|&i| (genome(i) >> bit) & 1 == 1)).words(),
            );
        }
        for lineage in (0..n / 4).step_by(6) {
            x.push_row(Bits::from_indices(n, (0..n).filter(|&i| i / 4 == lineage)).words());
        }
        let groups: Vec<u32> = (0..n).map(|i| (i / 4) as u32).collect();
        let duplicates = identical_columns(&x);
        assert_eq!(duplicates.len(), 20);
        let cfg = EvalConfig::default();
        let mut differs = [false; 4];
        // The lasso must also differ on the random split itself, so a fit
        // that drops the groups only there cannot pass.
        let mut random_lasso_differs = false;
        for (scheme, folds) in outer_splits(&y, &groups, &duplicates, &cfg) {
            let tuning = tuning_groups(scheme, &groups, &duplicates);
            check_split_plan(scheme, &folds, &y, tuning.as_deref(), 40, &cfg).unwrap();
            for (fi, (train, _)) in fold_masks(&folds, cfg.folds).into_iter().enumerate() {
                let (cx, _) = training_candidates(&x, &train, 3);
                let seed = fold_seed(&cfg, fi);
                let built = |ids: Option<&[u32]>| {
                    let inner = training_folds(&y, &train, ids, cfg.inner_folds, seed);
                    // The lasso's penalty is chosen on a holdout drawn here.
                    let (fit, val) =
                        tuning_holdout(&train, &y, ids, "group", LASSO_HOLDOUT, seed).unwrap();
                    let feats = screen(&cx, &y, &fit, cfg.screened_features);
                    let path = lasso::fit_path_with_status(
                        &cx,
                        &y,
                        &fit,
                        &feats,
                        cfg.lasso_lambdas,
                        400,
                        &mut |_| {},
                    );
                    let penalized = fit_lasso_from_path(
                        &cx,
                        &y,
                        &val,
                        path,
                        test_context(scheme, fi),
                        |lambdas| {
                            lasso::refit_path_with_status(
                                &cx,
                                &y,
                                &train,
                                &feats,
                                lambdas,
                                &mut |_| {},
                            )
                        },
                    )
                    .unwrap()
                    .0;
                    [
                        fit_scm(&cx, &y, &train, &inner, &cfg),
                        fit_cart(&cx, &y, &train, &inner, &cfg),
                        fit_gbt(&cx, &y, &train, ids, "group", &cfg, seed)
                            .unwrap()
                            .0,
                        penalized,
                    ]
                    .map(|s| (s.scores, s.complexity, s.lasso_fit))
                };
                let (grouped, plain) = (built(tuning.as_deref()), built(None));
                for (li, learner) in ["scm", "cart", "gbt", "lasso"].into_iter().enumerate() {
                    let context = test_context(scheme, fi);
                    let fit = fit_fold(learner, &cx, &y, &train, tuning.as_deref(), &cfg, context)
                        .unwrap();
                    assert_eq!(
                        (fit.scores, fit.complexity, fit.lasso_fit),
                        grouped[li],
                        "{learner} {scheme} outer fold {fi}"
                    );
                    differs[li] |= plain[li] != grouped[li];
                }
                random_lasso_differs |= scheme == "random" && plain[3] != grouped[3];
            }
        }
        assert_eq!(differs, [true; 4]);
        assert!(random_lasso_differs);
    }

    #[test]
    fn isolating_duplicates_takes_away_the_calls_a_genome_borrows_from_its_copy() {
        // The S. flexneri ciprofloxacin shape at small scale: 70 pairs of
        // identical genomes, each pair its own lineage. A marker carries 25
        // resistant pairs; 20 more are resistant without it and differ from
        // the 25 susceptible pairs only by a 24-bit barcode they share with
        // their copy. A copy in training lets the forest recall those 40
        // genomes. Isolated, it has only the marker, so VME rises and
        // balanced accuracy falls.
        use rand::{Rng, SeedableRng};
        let pairs = 70;
        let n = 2 * pairs;
        let resistant = |p: usize| p < 45;
        let y = Bits::from_indices(n, (0..n).filter(|&i| resistant(i / 2)));
        let mut x = BitMatrix::with_cols(n);
        x.push_row(Bits::from_indices(n, (0..n).filter(|&i| i / 2 < 25)).words());
        let mut rng = rand_chacha::ChaCha8Rng::seed_from_u64(9);
        for _ in 0..24 {
            let on: Vec<usize> = (0..pairs).filter(|_| rng.random_bool(0.5)).collect();
            x.push_row(Bits::from_indices(n, on.iter().flat_map(|&p| [2 * p, 2 * p + 1])).words());
        }
        let groups: Vec<u32> = (0..n).map(|i| (i / 2) as u32).collect();
        let duplicates: Vec<Vec<usize>> = (0..pairs).map(|p| vec![2 * p, 2 * p + 1]).collect();
        let cfg = EvalConfig {
            folds: 5,
            inner_folds: 3,
            forest: ForestParams {
                n_trees: 60,
                ..Default::default()
            },
            gbt: GbtParams {
                rounds: 30,
                ..Default::default()
            },
            lasso_lambdas: 5,
            ..Default::default()
        };
        let rf = |duplicates: &[Vec<usize>]| {
            let out = evaluate(
                "fixture__fixturemycin",
                &x,
                &y,
                &groups,
                duplicates,
                3,
                None,
                &cfg,
            )
            .unwrap();
            let random = &out.splits[0];
            let rf = random.learners.iter().find(|l| l.learner == "rf").unwrap();
            rf.metrics.clone().unwrap()
        };
        let (leaky, isolated) = (rf(&[]), rf(&duplicates));
        assert!(
            isolated.vme > leaky.vme + 0.1
                && isolated.balanced_accuracy < leaky.balanced_accuracy - 0.05,
            "leaky {leaky:?}, isolated {isolated:?}"
        );
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
            let out = evaluate(
                "fixture__fixturemycin",
                &x,
                &y,
                &groups,
                &[],
                3,
                None,
                &small_config(),
            )
            .unwrap();
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
        // Same-lineage pairs, so random tuning also runs on duplicate groups.
        let pairs: Vec<Vec<usize>> = (0..15).step_by(3).map(|i| vec![i, i + 15]).collect();
        let cfg = small_config();
        let truth: Vec<bool> = (0..n).map(|i| y.contains(i)).collect();
        let metrics_on = |idx: &[usize], scores: &[f64], ranked: bool| {
            binary_metrics(
                &idx.iter().map(|&i| truth[i]).collect::<Vec<_>>(),
                &idx.iter().map(|&i| scores[i]).collect::<Vec<_>>(),
                0.5,
                ranked,
            )
        };
        for duplicates in [vec![], pairs] {
            let out = evaluate(
                "fixture__fixturemycin",
                &x,
                &y,
                &groups,
                &duplicates,
                3,
                None,
                &cfg,
            )
            .unwrap();
            for split in &out.splits {
                assert!(split.unavailable_reason.is_none(), "{split:?}");
                // The saved fold ids are the partition the split was scored on.
                let folds = match split.scheme.as_str() {
                    "random" => isolate_groups(
                        &stratified_folds(&y, cfg.folds, cfg.seed),
                        &duplicates,
                        cfg.folds,
                    ),
                    _ => grouped_folds(&y, &groups, cfg.folds, cfg.seed),
                };
                assert_eq!(split.test_fold, folds);
                let tuning = tuning_groups(&split.scheme, &groups, &duplicates);
                for learner in &split.learners {
                    let called: Vec<usize> =
                        (0..n).filter(|&i| !learner.scores[i].is_nan()).collect();
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
                        let refit = fit_fold(
                            &learner.learner,
                            &cx,
                            &y,
                            &train,
                            tuning.as_deref(),
                            &cfg,
                            test_context(&split.scheme, f),
                        )
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
            lasso_fit: None,
        };
        let fits = [
            fitted([0.9, 0.2, 0.0, 0.0, 0.0, 0.0]),
            Scored {
                scores: None,
                complexity: 0.0,
                lasso_fit: None,
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
        let out = evaluate(
            "fixture__fixturemycin",
            &BitMatrix::with_cols(n),
            &y,
            &groups,
            &[],
            3,
            None,
            &cfg,
        )
        .unwrap();
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
                lasso_fit: None,
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
                let s = fit_fold(
                    learner,
                    &x,
                    labels,
                    &train,
                    None,
                    &cfg,
                    test_context("random", 1),
                )
                .unwrap();
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
        let (fit, val) = tuning_holdout(&all, &y, None, "group", 0.15, 42).unwrap();
        let mut x = BitMatrix::with_cols(n);
        x.push_row(fit.and(&y).or(&val.and_not(&y)).words());
        let (scored, model) = fit_gbt(&x, &y, &all, None, "group", &small_config(), 42).unwrap();
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
        let finals = final_models("fixture__fixturemycin", &x, &y, None, 3, &cfg).unwrap();
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
            fit_fold(
                "rf",
                &x,
                &y,
                &train,
                None,
                cfg,
                test_context("random", fold),
            )
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
        let (fit, _) = tuning_holdout(&all, &y, None, "group", 0.2, 42).unwrap();
        let mut constant_on_fit = BitMatrix::with_cols(n);
        constant_on_fit.push_row(fit.words());
        for x in [BitMatrix::with_cols(n), constant_on_fit] {
            let (scored, model) = fit_lasso(
                &x,
                &y,
                &all,
                None,
                &small_config(),
                42,
                test_context("random", 0),
            )
            .unwrap();
            assert!(model.is_none());
            assert!(scored.scores.is_none());
            assert_eq!(scored.complexity, 0.0);
            assert_eq!(scored.lasso_fit, Some(LassoFitAudit::empty()));
        }
    }

    #[test]
    fn a_lasso_refit_that_keeps_no_feature_is_a_null_fit_without_a_call() {
        let n = 260;
        let y = Bits::from_indices(n, (0..n).filter(|i| i % 9 == 0));
        let all = Bits::ones(n);
        let cfg = small_config();
        let (fit, val) = tuning_holdout(&all, &y, None, "group", 0.2, 42).unwrap();
        let mut x = BitMatrix::with_cols(n);
        x.push_row(fit.and(&y).or(&val.and_not(&y)).words());
        assert!(!lasso::fit_path(&x, &y, &fit, &[0], cfg.lasso_lambdas, 400).is_empty());
        let (scored, model) =
            fit_lasso(&x, &y, &all, None, &cfg, 42, test_context("random", 0)).unwrap();
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
        let (fit, val) = tuning_holdout(&all, &y, None, "group", 0.2, 42).unwrap();
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
        let (_, model) =
            fit_lasso(&x, &y, &all, None, &cfg, 42, test_context("random", 0)).unwrap();
        assert_eq!(model.unwrap().lambda, selected.lambda);
    }

    #[test]
    fn unconverged_lasso_candidates_cannot_win_validation_selection() {
        let n = 6;
        let y = Bits::from_indices(n, [0, 2, 4]);
        let mut x = BitMatrix::with_cols(n);
        x.push_row(y.words());
        let bad = canned_lasso(0.5, -10.0, 20.0, false);
        let good = canned_lasso(0.25, -1.0, 2.0, true);
        let (scored, model) = fit_lasso_from_path(
            &x,
            &y,
            &Bits::ones(n),
            vec![bad, good.clone()],
            test_context("random", 2),
            move |lambdas| {
                assert_eq!(lambdas, [0.5, 0.25]);
                vec![good]
            },
        )
        .unwrap();
        assert_eq!(model.unwrap().lambda, 0.25);
        let audit = scored.lasso_fit.unwrap();
        assert_eq!((audit.path_points, audit.path_unconverged), (2, 1));
        assert_eq!(audit.selected_point, Some(1));
        assert!(audit.selected_convergence.unwrap().converged);
        assert!(audit.refit_convergence.unwrap().converged);
    }

    #[test]
    fn a_path_without_a_converged_candidate_fails_with_full_context() {
        let n = 6;
        let y = Bits::from_indices(n, [0, 2, 4]);
        let mut x = BitMatrix::with_cols(n);
        x.push_row(y.words());
        let error = fit_lasso_from_path(
            &x,
            &y,
            &Bits::ones(n),
            vec![canned_lasso(0.5, -10.0, 20.0, false)],
            test_context("random", 2),
            |_| unreachable!("an unconverged path must not be refit"),
        )
        .err()
        .unwrap()
        .to_string();
        for expected in [
            "lasso path failure",
            "no converged candidate",
            "dataset fixture__fixturemycin",
            "split random",
            "outer fold 2",
        ] {
            assert!(error.contains(expected), "{error}");
        }
    }

    #[test]
    fn an_unconverged_final_lasso_refit_fails_with_full_context() {
        let n = 6;
        let y = Bits::from_indices(n, [0, 2, 4]);
        let mut x = BitMatrix::with_cols(n);
        x.push_row(y.words());
        let error = fit_lasso_from_path(
            &x,
            &y,
            &Bits::ones(n),
            vec![canned_lasso(0.25, -1.0, 2.0, true)],
            test_context("lineage", 3),
            |_| vec![canned_lasso(0.25, -1.0, 2.0, false)],
        )
        .err()
        .unwrap()
        .to_string();
        for expected in [
            "lasso refit failure",
            "selected point 0 is unconverged",
            "dataset fixture__fixturemycin",
            "split lineage",
            "outer fold 3",
        ] {
            assert!(error.contains(expected), "{error}");
        }
    }

    #[test]
    fn final_model_lasso_failures_name_the_final_fit() {
        let n = 6;
        let y = Bits::from_indices(n, [0, 2, 4]);
        let mut x = BitMatrix::with_cols(n);
        x.push_row(y.words());
        let finals = LassoContext {
            dataset: "fixture__fixturemycin",
            scheme: "final",
            outer_fold: None,
        };
        let refits: [Vec<lasso::LassoFit>; 2] =
            [vec![canned_lasso(0.25, -1.0, 2.0, false)], vec![]];
        for (refit, detail) in refits.into_iter().zip([
            "selected point 0 is unconverged",
            "selected point 0 returned an empty refit path",
        ]) {
            let error = fit_lasso_from_path(
                &x,
                &y,
                &Bits::ones(n),
                vec![canned_lasso(0.25, -1.0, 2.0, true)],
                finals,
                |_| refit,
            )
            .err()
            .unwrap()
            .to_string();
            assert!(error.contains("lasso refit failure"), "{error}");
            assert!(error.contains(detail), "{error}");
            assert!(
                error.contains("dataset fixture__fixturemycin final model on all samples"),
                "{error}"
            );
            assert!(!error.contains("split"), "{error}");
        }
    }

    #[test]
    fn lasso_only_repeats_the_evaluator_lasso_fits() {
        // Without duplicates, and with same-lineage pairs that put random
        // tuning on duplicate groups.
        let pairs: Vec<Vec<usize>> = (0..200).step_by(10).map(|i| vec![i, i + 1]).collect();
        for duplicates in [vec![], pairs] {
            lasso_only_repeats_the_evaluator_lasso_fits_with(&duplicates);
        }
    }

    fn lasso_only_repeats_the_evaluator_lasso_fits_with(duplicates: &[Vec<usize>]) {
        // 200 samples in 40 lineages: two markers of resistance, noise, and
        // a resistant-only feature on three samples. At min_support 3 it is
        // a candidate on all samples but not in folds that hold one out.
        let n = 200;
        let y = Bits::from_indices(n, (0..n).filter(|i| i % 3 == 0));
        let rules: [fn(usize) -> bool; 6] = [
            |i| (i % 3 == 0 && i % 10 != 0) || i % 20 == 1,
            |i| i % 3 == 0 && i % 7 != 0,
            |i| i % 5 == 0,
            |i| i % 7 == 1,
            |i| i % 4 < 2,
            |i| [0, 3, 6].contains(&i),
        ];
        let mut x = BitMatrix::with_cols(n);
        for rule in rules {
            x.push_row(Bits::from_indices(n, (0..n).filter(|&i| rule(i))).words());
        }
        let groups: Vec<u32> = (0..n).map(|i| (i / 5) as u32).collect();
        let mut cfg = small_config();
        cfg.folds = 3;
        cfg.screened_features = 6;
        cfg.lasso_lambdas = 8;
        let dataset = "fixture__fixturemycin";
        let out = evaluate(dataset, &x, &y, &groups, duplicates, 3, None, &cfg).unwrap();
        let mut events = Vec::new();
        let mut points = Vec::new();
        lasso_only(
            dataset,
            &x,
            &y,
            &groups,
            duplicates,
            3,
            &cfg,
            &mut |e| {
                events.push(e.clone());
                Ok(())
            },
            &mut |p| {
                points.push((
                    p.scheme.to_string(),
                    p.outer_fold,
                    p.stage.to_string(),
                    p.point,
                ))
            },
        )
        .unwrap();

        // Each fit is announced, then reported, in evaluation order, and its
        // penalty points were reported one by one as they were solved.
        let mut fits = Vec::new();
        let mut points = points.into_iter();
        for pair in events.chunks(2) {
            let (LassoOnlyEvent::Start { scheme, outer_fold }, LassoOnlyEvent::Fit(fit)) =
                (&pair[0], &pair[1])
            else {
                panic!("{pair:?}");
            };
            assert_eq!((scheme, outer_fold), (&fit.scheme, &fit.outer_fold));
            assert!(fit.error.is_none() && fit.seconds >= 0.0);
            for (stage, trace) in [("path", &fit.trace.path), ("refit", &fit.trace.refit)] {
                for (i, c) in trace.iter().enumerate() {
                    let (s, f, st, p) = points.next().unwrap();
                    assert_eq!(
                        (&s, f, st.as_str(), p.point),
                        (scheme, *outer_fold, stage, i)
                    );
                    assert_eq!((p.irls_steps, p.converged), (c.irls_steps, c.converged));
                    assert!(p.kkt <= 1e-6 && p.seconds >= 0.0, "{p:?}");
                }
            }
            fits.push(fit);
        }
        assert!(points.next().is_none());
        let order: Vec<(&str, Option<usize>)> = fits
            .iter()
            .map(|f| (f.scheme.as_str(), f.outer_fold))
            .collect();
        assert_eq!(
            order,
            [
                ("random", Some(0)),
                ("random", Some(1)),
                ("random", Some(2)),
                ("lineage", Some(0)),
                ("lineage", Some(1)),
                ("lineage", Some(2)),
                ("final", None),
            ]
        );
        for split in &out.splits {
            let lasso = split
                .learners
                .iter()
                .find(|l| l.learner == "lasso")
                .unwrap();
            for (fold, audit) in lasso.lasso_fits.as_ref().unwrap().iter().enumerate() {
                let fit = fits
                    .iter()
                    .find(|f| f.scheme == split.scheme && f.outer_fold == Some(fold))
                    .unwrap();
                assert_eq!(fit.audit.as_ref(), Some(audit));
                let weights = fit.model.as_ref().map_or(0, |m| m.coef.len());
                assert_eq!(weights as f64, lasso.fold_complexity[fold]);
                assert_eq!(fit.trace.path.len(), audit.path_points);
                assert_eq!(fit.trace.refit.last().copied(), audit.refit_convergence);
            }
        }
        let last = fits.last().unwrap();
        assert_eq!(last.audit, out.finals.lasso_fit);
        let (a, b) = (
            last.model.as_ref().unwrap(),
            out.finals.lasso.as_ref().unwrap(),
        );
        assert!(!a.coef.is_empty());
        assert_eq!(
            (a.lambda.to_bits(), a.intercept.to_bits(), &a.coef),
            (b.lambda.to_bits(), b.intercept.to_bits(), &b.coef)
        );
    }

    #[test]
    fn a_failed_lasso_fit_is_reported_with_its_trace_before_it_is_returned() {
        let capped = lasso::Convergence {
            irls_steps: 25,
            sweep_cap_hits: 25,
            converged: false,
        };
        let mut events = Vec::new();
        let error = observed_lasso(
            "lineage",
            Some(3),
            &mut |e| {
                events.push(e.clone());
                Ok(())
            },
            |trace| {
                trace.path = vec![capped];
                Err(anyhow::anyhow!("lasso path failure at outer fold 3"))
            },
        )
        .unwrap_err();
        assert!(error.to_string().contains("lasso path failure"));
        let [LassoOnlyEvent::Start { .. }, LassoOnlyEvent::Fit(fit)] = &events[..] else {
            panic!("{events:?}");
        };
        assert_eq!((fit.scheme.as_str(), fit.outer_fold), ("lineage", Some(3)));
        assert!(fit.error.as_ref().unwrap().contains("lasso path failure"));
        assert!(fit.audit.is_none() && fit.model.is_none());
        assert_eq!(fit.trace.path, [capped]);

        // An observer that refuses the start stops the scan before the fit.
        let refused = observed_lasso(
            "random",
            Some(0),
            &mut |_| anyhow::bail!("stop rule"),
            |_| unreachable!("a refused fit must not run"),
        );
        assert_eq!(refused.unwrap_err().to_string(), "stop rule");
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
            lasso_fits: None,
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
                    let s = fit_fold(
                        learner,
                        &cx,
                        &y,
                        &train,
                        None,
                        &cfg,
                        test_context("random", 1),
                    )
                    .unwrap();
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
