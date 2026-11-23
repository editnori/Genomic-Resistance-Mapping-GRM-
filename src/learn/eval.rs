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
//! `EvalConfig::learner_seed`, when set, seeds the forest's and boosting's
//! own draws in its place and leaves every partition as `seed` draws it.
//!
//! A boosted model that keeps no tree, or a lasso that keeps no weight, is a
//! null fit. A null fit is still a model: it predicts its training
//! prevalence (`train_prevalence`, resistant training samples over training
//! samples of the outer fold) for every sample, and the 0.5 cutoff calls
//! that score like any other, so a prevalence of exactly 0.5 is a resistant
//! call. From revision 14 every learner is therefore scored on every sample
//! of a split: `metrics`, the ROC points and each fold's balanced accuracy
//! cover them all. `no_call` still counts, by class, the samples of null-fit
//! folds, and a row with any such sample keeps the figure over the other
//! samples in `called_metrics` (null when every fold's fit was null), the
//! figure revisions 7 to 13 published as `metrics`. A fitted model's score
//! is a call even when it is exactly 0.5.
//!
//! `EvalOutput::fold_models` keeps what each outer fold's model reads
//! (`FoldSplit`), its features numbered as rows of the dataset matrix like
//! the final models'. The final models are fitted on every sample, so the
//! fold models are what shows whether the folds chose the same features.
//! They record fits `evaluate` already makes and change no score.

use super::gbt::{self, GbtParams, Objective};
use super::inner_pick;
use super::lasso;
use super::null_refit;
use super::scm::{self, Kind, Rule};
use super::screen_arm;
use super::tie::TieRule;
use super::tree::{self, ForestParams, TreeParams};
use super::{
    BinaryMetrics, NoCall, any_training_candidate, binary_metrics, fold_masks, grouped_folds,
    holdout, isolate_groups, roc_curve, screen, stratified_folds, train_prevalence,
    training_candidates, training_folds, tuning_holdout, validate_folds,
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
    /// Seeds every partition (both splits' outer folds, each fit's inner
    /// folds and tuning holdouts, the MIC regression's folds) and, unless
    /// `learner_seed` is set, the learners' own draws.
    pub seed: u64,
    /// Seeds the learners' own draws in place of `seed`: the forest's
    /// subsamples and node features and boosting's column samples, in every
    /// outer fold, the final models and the MIC regression. No partition
    /// reads it. None, the default and every release's setting, is left out
    /// of a report's `eval_config`, so such a report is unchanged.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub learner_seed: Option<u64>,
    /// The tree learners' tie rule. `TieRule::Exact` is an experiment, off by
    /// default; a report records the field only when it is on.
    #[serde(skip_serializing_if = "TieRule::is_float")]
    pub tie_rule: TieRule,
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
            learner_seed: None,
            tie_rule: TieRule::Float,
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
    /// Metrics on the pooled out-of-fold scores of every sample, a null-fit
    /// fold's samples at that fold's training prevalence. Rows written before
    /// revision 14 left a null fit's samples out (they are counted in
    /// `no_call`) and are None (null in JSON) when no sample had a call.
    pub metrics: Option<BinaryMetrics>,
    /// Revision 14 rows with a null-fit fold only: the metrics over the
    /// samples of fitted folds, which revisions 7 to 13 published as
    /// `metrics`. `Some(None)` (null in JSON) when every fold's fit was null.
    /// Absent on a row without a null fit, whose `metrics` already are that
    /// figure, and on every row written before revision 14.
    #[serde(
        default,
        skip_serializing_if = "Option::is_none",
        deserialize_with = "deserialize_present"
    )]
    pub called_metrics: Option<Option<BinaryMetrics>>,
    /// Balanced accuracy of each outer fold's calls, by fold id; a null-fit
    /// fold's is that of its constant prevalence call. Rows written before
    /// revision 14 hold NaN (null in JSON) for a null-fit fold.
    #[serde(deserialize_with = "deserialize_nullable_scores")]
    pub fold_balanced_accuracy: Vec<f64>,
    pub roc: Vec<(f64, f64)>,
    /// Out-of-fold resistant score per sample, in dataset sample order; a
    /// null-fit fold's samples carry its training prevalence. Rows written
    /// before revision 14 hold NaN (null in JSON) for those samples.
    #[serde(deserialize_with = "deserialize_nullable_scores")]
    pub scores: Vec<f64>,
    pub conformal: Vec<ConformalPoint>,
    /// Per-fold model size: SCM rules, CART leaves, forest trees, boosted
    /// trees kept by early stopping, or lasso non-zero weights. Zero for
    /// boosting or the lasso is a null fit. Rows without `no_call` scored
    /// that fold's samples 0.5, which the 0.5 cutoff calls resistant; rows
    /// of revisions 7 to 13 leave them uncalled; revision 14 rows score them
    /// at the fold's training prevalence.
    pub fold_complexity: Vec<f64>,
    pub seconds: f64,
    /// Samples in folds whose fit was null, by reported class: those without
    /// a fitted model's call. Absent on catalog rows, which call every
    /// sample, and on rows written before null fits became no-calls.
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

/// Reads a field that is present, null included, as Some: `called_metrics`
/// is absent (None), null (Some(None)) or an object.
fn deserialize_present<'de, D, T>(d: D) -> Result<Option<T>, D::Error>
where
    D: serde::Deserializer<'de>,
    T: Deserialize<'de>,
{
    T::deserialize(d).map(Some)
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

/// How many features a forest's or a boosted model's importance list
/// keeps: the final models' lists and each fold model's.
pub const TOP_WEIGHTED: usize = 50;

/// One outer fold's set covering machine: the kind and trade-off its inner
/// folds chose, and its rules in the order the greedy fit added them. The
/// fold fits list no tied rules (`scm::fit` without `track_equivalent`).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct FoldScm {
    pub kind: Kind,
    pub p: f64,
    pub rules: Vec<Rule>,
}

/// The outer-fold models of one split, each learner's by fold id: the model
/// fitted for the samples `test_fold` gives that fold, whose scores the
/// split's row holds for them unless the fit was null. Each has the shape
/// of the final model it stands beside:
///
/// - `scm`: as `FinalModels::scm`'s kind, p and rules;
/// - `cart`: the tree, as `FinalModels::cart`, grown to its selected depth;
/// - `rf` and `gbt`: the `TOP_WEIGHTED` features by importance, heaviest
///   first and equal weights by lower feature first, as `rf_importance` and
///   `gbt_importance`; a boosted model that kept no tree has none;
/// - `lasso`: the scored refit, as `FinalModels::lasso`; None when the
///   tuning path was empty, and without weights when the refit kept none.
///   A lasso fold model that is None or has no weights, like a boosting
///   fold with an empty list, is a null fit: its fold scored each test
///   sample at the fold's training prevalence (`pool_folds`,
///   docs/equal-scoring.md), not at the sigmoid of its intercept.
///
/// A feature is a row of the dataset matrix, the numbering of the final
/// models and of a report's feature cards. A fold chooses among the rows
/// that differ on its own training samples and names the lowest row of each
/// set that agrees there (`training_candidates`), so two folds, or a fold
/// and the final model, can name different rows for patterns that differ
/// only on samples one of them did not train on. Equal ids are therefore
/// the strict reading of "the same feature". A report also carries the
/// presence pattern of every feature in `pattern_features`, so those can be
/// compared sample by sample against `SplitResult::test_fold`.
///
/// A feature here is a marker a model reads and an association with the
/// reported label. It is not a cause.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct FoldSplit {
    /// The split's `SplitResult::scheme`.
    pub scheme: String,
    pub scm: Vec<FoldScm>,
    pub cart: Vec<tree::Tree>,
    pub rf: Vec<Vec<WeightedFeature>>,
    pub gbt: Vec<Vec<WeightedFeature>>,
    pub lasso: Vec<Option<lasso::LassoModel>>,
}

impl FoldSplit {
    fn new(scheme: &str) -> Self {
        FoldSplit {
            scheme: scheme.to_string(),
            scm: vec![],
            cart: vec![],
            rf: vec![],
            gbt: vec![],
            lasso: vec![],
        }
    }

    /// Adds the next fold's `model`, fitted on a candidate matrix whose row
    /// c is dataset row `rows[c]` (`training_candidates`), under dataset
    /// row numbers.
    fn push(&mut self, model: FoldModel, rows: &[u32]) {
        let row = |f: u32| rows[f as usize];
        let weighted = |list: Vec<WeightedFeature>| -> Vec<WeightedFeature> {
            list.into_iter()
                .map(|w| WeightedFeature {
                    feature: row(w.feature),
                    weight: w.weight,
                })
                .collect()
        };
        match model {
            FoldModel::Scm(mut m) => {
                for rule in &mut m.rules {
                    rule.feature = row(rule.feature);
                }
                self.scm.push(m);
            }
            FoldModel::Cart(mut t) => {
                for split in t.nodes.iter_mut().filter_map(|n| n.split.as_mut()) {
                    split.feature = row(split.feature);
                }
                self.cart.push(t);
            }
            FoldModel::Rf(list) => self.rf.push(weighted(list)),
            FoldModel::Gbt(list) => self.gbt.push(weighted(list)),
            FoldModel::Lasso(mut m) => {
                for (feature, _) in m.iter_mut().flat_map(|m| &mut m.coef) {
                    *feature = row(*feature);
                }
                self.lasso.push(m);
            }
        }
    }

    /// Every feature of the split's scm rules, cart splits and lasso
    /// weights, each once: the learners whose fold model is its whole
    /// feature list. The forest's and boosting's lists are the top of a
    /// ranking, and their features are left out.
    pub fn pattern_features(&self) -> std::collections::BTreeSet<u32> {
        let mut named = std::collections::BTreeSet::new();
        named.extend(self.scm.iter().flat_map(|m| &m.rules).map(|r| r.feature));
        named.extend(
            self.cart
                .iter()
                .flat_map(|t| &t.nodes)
                .filter_map(|n| n.split.as_ref().map(|s| s.feature)),
        );
        named.extend(
            self.lasso
                .iter()
                .flatten()
                .flat_map(|m| &m.coef)
                .map(|c| c.0),
        );
        named
    }
}

/// One outer fold's model as `FoldSplit` keeps it, its features numbered as
/// rows of the matrix it was fitted on.
pub(crate) enum FoldModel {
    Scm(FoldScm),
    Cart(tree::Tree),
    Rf(Vec<WeightedFeature>),
    Gbt(Vec<WeightedFeature>),
    Lasso(Option<lasso::LassoModel>),
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct EvalOutput {
    pub splits: Vec<SplitResult>,
    /// The outer-fold models of each split that ran, in `splits` order.
    /// Empty when read from an output written before they were kept.
    #[serde(default)]
    pub fold_models: Vec<FoldSplit>,
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

pub(crate) struct Scored {
    /// A score for every sample, or None for a null fit, which makes no call.
    pub(crate) scores: Option<Vec<f64>>,
    pub(crate) complexity: f64,
    pub(crate) lasso_fit: Option<LassoFitAudit>,
}

pub(crate) fn fit_scm(
    x: &BitMatrix,
    y: &Bits,
    train: &Bits,
    inner: &[usize],
    cfg: &EvalConfig,
) -> (Scored, FoldScm) {
    let c = scm::select(x, y, train, cfg.scm_max_rules, inner, cfg.inner_folds);
    let m = scm::fit(x, y, train, c.kind, c.p, c.n_rules, false);
    let pred = m.predict(x);
    (
        Scored {
            scores: Some((0..x.cols).map(|i| pred.contains(i) as u8 as f64).collect()),
            complexity: m.rules.len() as f64,
            lasso_fit: None,
        },
        FoldScm {
            kind: m.kind,
            p: m.p,
            rules: m.rules,
        },
    )
}

pub(crate) fn fit_cart(
    x: &BitMatrix,
    y: &Bits,
    train: &Bits,
    inner: &[usize],
    cfg: &EvalConfig,
) -> (Scored, tree::Tree) {
    let (depth, _) = tree::select_depth(
        x,
        y,
        train,
        inner,
        cfg.inner_folds,
        cfg.cart_max_depth,
        cfg.tie_rule,
    );
    let t = tree::fit_tree(
        x,
        y,
        train,
        &TreeParams {
            max_depth: depth,
            min_leaf: 1,
            max_features: None,
            seed: 0,
            tie: cfg.tie_rule,
        },
    );
    let leaves = t.nodes.iter().filter(|n| n.split.is_none()).count();
    (
        Scored {
            scores: Some(t.predict(x, depth)),
            complexity: leaves as f64,
            lasso_fit: None,
        },
        t,
    )
}

/// The forest fits on the same number of chi-square-screened features as
/// boosting and the lasso, screened on its training samples. Its trees name
/// positions in that screen, so importance is returned keyed by rows of `x`.
pub(crate) fn fit_rf(
    x: &BitMatrix,
    y: &Bits,
    train: &Bits,
    cfg: &EvalConfig,
    seed: u64,
) -> (Scored, std::collections::HashMap<u32, f64>) {
    let feats = screen_arm::screen_fit(x, y, train, cfg.screened_features);
    let mut screened = BitMatrix::with_cols(x.cols);
    for &f in &feats {
        screened.push_row(x.row(f as usize));
    }
    let params = ForestParams {
        seed,
        tie: cfg.tie_rule,
        ..cfg.forest
    };
    let forest = tree::fit_forest(&screened, y, train, &params);
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
pub(crate) fn gbt_scores(m: &gbt::Gbt, x: &BitMatrix) -> Option<Vec<f64>> {
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

/// `seed` draws the early-stopping holdout, a partition; `draws` seeds
/// boosting's column samples. They are equal unless
/// `EvalConfig::learner_seed` is set.
#[allow(clippy::too_many_arguments)]
fn fit_gbt(
    x: &BitMatrix,
    y: &Bits,
    train: &Bits,
    groups: Option<&[u32]>,
    unit: &str,
    cfg: &EvalConfig,
    seed: u64,
    draws: u64,
) -> Result<(Scored, gbt::Gbt)> {
    let (fit_m, val) = tuning_holdout(train, y, groups, unit, GBT_HOLDOUT, seed)?;
    let feats = screen_arm::screen_fit(x, y, &fit_m, cfg.screened_features);
    let m = gbt::fit(
        x,
        &Objective::Logistic { y },
        &fit_m,
        &val,
        &feats,
        &GbtParams {
            seed: draws,
            tie: cfg.tie_rule,
            ..cfg.gbt
        },
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
pub(crate) struct LassoContext<'a> {
    pub(crate) dataset: &'a str,
    pub(crate) scheme: &'a str,
    pub(crate) outer_fold: Option<usize>,
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
    let feats = screen_arm::screen_fit(x, y, &fit_m, cfg.screened_features);
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

/// `fit_lasso` for other modules' tests: its scores and selected point.
#[cfg(test)]
pub(crate) fn fit_lasso_for_tests(
    x: &BitMatrix,
    y: &Bits,
    train: &Bits,
    cfg: &EvalConfig,
    seed: u64,
    context: LassoContext<'_>,
) -> (Option<Vec<f64>>, Option<usize>) {
    let (scored, _) = fit_lasso(x, y, train, None, cfg, seed, context).unwrap();
    (
        scored.scores,
        scored.lasso_fit.and_then(|a| a.selected_point),
    )
}

/// `fit_gbt` for other modules' tests: its scores, None for a null fit.
/// `seed` draws both the holdout and the column samples.
#[cfg(test)]
pub(crate) fn fit_gbt_for_tests(
    x: &BitMatrix,
    y: &Bits,
    train: &Bits,
    cfg: &EvalConfig,
    seed: u64,
) -> Option<Vec<f64>> {
    fit_gbt(x, y, train, None, "group", cfg, seed, seed)
        .unwrap()
        .0
        .scores
}

/// Select and refit one already-solved path. Keeping this separate lets the
/// negative tests inject capped statuses without changing solver limits.
pub(crate) fn fit_lasso_from_path(
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

/// The seed of every partition drawn inside outer fold `fold` (its inner
/// folds and tuning holdouts) and, unless `EvalConfig::learner_seed` is set,
/// of the learners' own draws there: `seed + 101 * fold`, wrapping at 2^64
/// (`ablate::plan::fold_seed`). Every seed sum here wraps, because `--seed`
/// takes any 64-bit value.
fn fold_seed(cfg: &EvalConfig, fold: usize) -> u64 {
    crate::ablate::plan::fold_seed(cfg.seed, fold)
}

/// The base seed of the learners' own draws: `EvalConfig::learner_seed`, or
/// `EvalConfig::seed` when it is not set.
pub(crate) fn draw_seed(cfg: &EvalConfig) -> u64 {
    cfg.learner_seed.unwrap_or(cfg.seed)
}

/// The seed of the forest's and boosting's own draws in outer fold `fold`:
/// `fold_seed` on `draw_seed`, so it is `fold_seed` unless
/// `EvalConfig::learner_seed` is set.
fn fold_draw_seed(cfg: &EvalConfig, fold: usize) -> u64 {
    crate::ablate::plan::fold_seed(draw_seed(cfg), fold)
}

/// One learner's fit on outer fold `fold`'s training samples, scoring every
/// sample, with the fitted model as `FoldSplit` keeps it: its features are
/// rows of `x`. Evaluation and any rebuild of a saved score call this alone.
fn fit_fold(
    learner: &str,
    x: &BitMatrix,
    y: &Bits,
    train: &Bits,
    groups: Option<&[u32]>,
    cfg: &EvalConfig,
    context: LassoContext<'_>,
) -> Result<(Scored, FoldModel)> {
    let fold = context.outer_fold.expect("an outer fit needs its fold");
    let seed = fold_seed(cfg, fold);
    let draws = fold_draw_seed(cfg, fold);
    fit_seeded(learner, x, y, train, groups, cfg, seed, draws, context)
}

/// `fit_fold` with its seeds given: `seed` for the fit's inner folds and
/// tuning holdouts, `draws` for the forest's and boosting's own draws. The
/// C9 experiment's inner fits (`inner_pick`) call it with seeds of their own.
#[allow(clippy::too_many_arguments)]
pub(crate) fn fit_seeded(
    learner: &str,
    x: &BitMatrix,
    y: &Bits,
    train: &Bits,
    groups: Option<&[u32]>,
    cfg: &EvalConfig,
    seed: u64,
    draws: u64,
    context: LassoContext<'_>,
) -> Result<(Scored, FoldModel)> {
    let inner = || training_folds(y, train, groups, cfg.inner_folds, seed);
    Ok(match learner {
        "scm" => {
            let (scored, model) = fit_scm(x, y, train, &inner(), cfg);
            (scored, FoldModel::Scm(model))
        }
        "cart" => {
            let (scored, tree) = fit_cart(x, y, train, &inner(), cfg);
            (scored, FoldModel::Cart(tree))
        }
        "rf" => {
            let (scored, importance) = fit_rf(x, y, train, cfg, draws);
            (
                scored,
                FoldModel::Rf(top_weighted(importance, TOP_WEIGHTED)),
            )
        }
        "gbt" => {
            let unit = tuning_unit(context.scheme);
            let (scored, model) = fit_gbt(x, y, train, groups, unit, cfg, seed, draws)?;
            let top = top_weighted(model.importance(), TOP_WEIGHTED);
            (scored, FoldModel::Gbt(top))
        }
        "lasso" => {
            let (scored, model) = fit_lasso(x, y, train, groups, cfg, seed, context)?;
            (scored, FoldModel::Lasso(model))
        }
        _ => unreachable!(),
    })
}

/// Pools one learner's outer-fold fits, in fold order. Sample i takes its
/// score from the fit of fold `test_fold[i]`, the one fit that did not train
/// on it. A null fit predicts its training prevalence, `train_prevalence[f]`
/// (`learn::train_prevalence` of the fold's training samples), for every
/// sample of its fold, and the 0.5 cutoff calls it like any other score, so
/// `metrics`, the ROC points and each fold's balanced accuracy cover every
/// sample. `no_call` counts the samples of null-fit folds by class, and a
/// row with any keeps the metrics over the rest in `called_metrics`.
pub(crate) fn pool_folds(
    learner: &str,
    truth: &[bool],
    test_fold: &[usize],
    fits: &[Scored],
    train_prevalence: &[f64],
    seconds: f64,
) -> Result<LearnerResult> {
    let n = truth.len();
    ensure!(
        test_fold.len() == n && test_fold.iter().all(|&f| f < fits.len()),
        "a sample has no outer test fold"
    );
    ensure!(
        train_prevalence.len() == fits.len()
            && train_prevalence
                .iter()
                .all(|p| p.is_finite() && (0.0..=1.0).contains(p)),
        "{learner} needs one training prevalence in [0, 1] per outer fold"
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
    let null = |i: usize| fits[test_fold[i]].scores.is_none();
    let mut scores = vec![f64::NAN; n];
    for (i, &f) in test_fold.iter().enumerate() {
        scores[i] = match &fits[f].scores {
            Some(s) => {
                ensure!(
                    s[i].is_finite(),
                    "{learner} produced a non-finite test score"
                );
                s[i]
            }
            None => train_prevalence[f],
        };
    }
    let pooled = |keep: &dyn Fn(usize) -> bool| -> (Vec<bool>, Vec<f64>) {
        (0..n)
            .filter(|&i| keep(i))
            .map(|i| (truth[i], scores[i]))
            .unzip()
    };
    let fold_balanced_accuracy = (0..fits.len())
        .map(|f| {
            let (t, s) = pooled(&|i| test_fold[i] == f);
            binary_metrics(&t, &s, 0.5, false).balanced_accuracy
        })
        .collect();
    let in_null_fold = |class: bool| (0..n).filter(|&i| null(i) && truth[i] == class).count();
    let no_call = NoCall {
        resistant: in_null_fold(true),
        susceptible: in_null_fold(false),
    };
    let ranked = !matches!(learner, "scm");
    let (t, s) = pooled(&|_| true);
    let called_metrics = (no_call.resistant + no_call.susceptible > 0).then(|| {
        let (t, s) = pooled(&|i| !null(i));
        (!t.is_empty()).then(|| binary_metrics(&t, &s, 0.5, ranked))
    });
    Ok(LearnerResult {
        learner: learner.to_string(),
        metrics: (!t.is_empty()).then(|| binary_metrics(&t, &s, 0.5, ranked)),
        called_metrics,
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
    // The C4 experiment, off by default (`screen_arm`): its strata are the
    // lineage clusters, each fit reading its own fitting samples' alone.
    let arm = screen_arm::install(groups);
    if let Some(m) = &mic {
        ensure!(m.mic.len() == n, "MIC dimensions differ");
    }
    let truth: Vec<bool> = (0..n).map(|i| y.contains(i)).collect();
    let n_groups = count_groups(groups);
    let mut splits = Vec::new();
    let mut fold_models = Vec::new();

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
        let trains: Vec<Bits> = fold_masks(folds, cfg.folds)
            .into_iter()
            .map(|(train, _)| train)
            .collect();
        let OuterFits {
            fits,
            models,
            seconds,
            prevalence,
            ..
        } = fit_outer_folds(
            dataset,
            scheme,
            x,
            y,
            &trains,
            tuning.as_deref(),
            min_support,
            cfg,
        )?;
        let learners = LEARNERS
            .iter()
            .zip(&fits)
            .zip(seconds)
            .map(|((name, fits), seconds)| {
                pool_folds(name, &truth, folds, fits, &prevalence, seconds)
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
        fold_models.push(models);
    }

    // The C9 experiment, off by default (`inner_pick`), writes only its own
    // file and changes no fit, score or report field.
    if inner_pick::enabled() {
        inner_pick::record_drawn(dataset, x, y, groups, duplicates, min_support, cfg, &splits)?;
    }

    let mic_result = mic.and_then(|m| evaluate_mic(x, y, min_support, &m, cfg));
    let tuning = final_tuning(groups, duplicates);
    let finals = final_models(dataset, x, y, tuning.as_deref(), min_support, cfg)?;
    if arm.is_some() {
        screen_arm::finish(
            dataset,
            "drawn",
            cfg.seed,
            cfg.screened_features,
            min_support,
            n,
            x.rows,
        )?;
    }
    drop(arm);
    Ok(EvalOutput {
        splits,
        fold_models,
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
        check_fold_tuning(scheme, fi, &train, y, tuning, cfg)?;
    }
    Ok(())
}

/// Checks the partitions outer fold `fi`'s fits draw inside its training
/// samples: the inner folds, and both tuning holdouts.
fn check_fold_tuning(
    scheme: &str,
    fi: usize,
    train: &Bits,
    y: &Bits,
    tuning: Option<&[u32]>,
    cfg: &EvalConfig,
) -> Result<()> {
    let seed = fold_seed(cfg, fi);
    let inner = training_folds(y, train, tuning, cfg.inner_folds, seed);
    validate_folds(y, train, &inner, cfg.inner_folds)
        .with_context(|| format!("invalid {scheme} inner evaluation in outer fold {fi}"))?;
    for fraction in [GBT_HOLDOUT, LASSO_HOLDOUT] {
        tuning_holdout(train, y, tuning, tuning_unit(scheme), fraction, seed)
            .with_context(|| format!("invalid {scheme} tuning holdout in outer fold {fi}"))?;
    }
    Ok(())
}

/// The partitions a fit seeded with `seed` draws inside `train`, as
/// `check_fold_tuning` checks an outer fold's: its inner folds and both
/// tuning holdouts. The C9 experiment checks each inner fit with it
/// (`inner_pick`).
pub(crate) fn check_tuning_seeded(
    scheme: &str,
    train: &Bits,
    y: &Bits,
    tuning: Option<&[u32]>,
    cfg: &EvalConfig,
    seed: u64,
) -> Result<()> {
    let inner = training_folds(y, train, tuning, cfg.inner_folds, seed);
    validate_folds(y, train, &inner, cfg.inner_folds)
        .with_context(|| format!("invalid {scheme} inner evaluation"))?;
    for fraction in [GBT_HOLDOUT, LASSO_HOLDOUT] {
        tuning_holdout(train, y, tuning, tuning_unit(scheme), fraction, seed)
            .with_context(|| format!("invalid {scheme} tuning holdout"))?;
    }
    Ok(())
}

/// Every learner's fits on a split's outer folds.
struct OuterFits {
    /// Per learner, in `LEARNERS` order, one fit per outer fold.
    fits: Vec<Vec<Scored>>,
    /// Each fold's model from the fit whose scores `fits` holds for it, its
    /// features numbered as dataset rows (`FoldSplit::push`).
    models: FoldSplit,
    /// Per learner, the seconds its fits took together.
    seconds: Vec<f64>,
    /// Per outer fold, the resistant share of its training samples.
    prevalence: Vec<f64>,
    /// Per outer fold, the candidate features its fits chose from.
    candidates: Vec<usize>,
}

/// Fits every learner on each of `trains`, outer fold `fi` on `trains[fi]`:
/// its candidates come from those samples alone (`training_candidates`) and
/// every draw is seeded from the fold (`fit_fold`: partitions from
/// `fold_seed`, the forest's and boosting's own draws from `fold_draw_seed`).
/// Each fit gives the fold's scores and its model together, so a fold's
/// model is the one that scored it. `evaluate` fits the folds it draws
/// through this and keeps the models as `fold_models`; `evaluate_assigned`
/// fits the folds it is handed and keeps no model, so a supplied partition
/// is fitted exactly as a drawn one. With the null-fold experiment on
/// (`null_refit::enabled`), each boosting or lasso fold that comes out null
/// is refitted after its timer with these same inputs and seeds and recorded
/// (`null_refit::record_fold`); `evaluate_assigned` refuses to run then.
#[allow(clippy::too_many_arguments)]
fn fit_outer_folds(
    dataset: &str,
    scheme: &str,
    x: &BitMatrix,
    y: &Bits,
    trains: &[Bits],
    tuning: Option<&[u32]>,
    min_support: u32,
    cfg: &EvalConfig,
) -> Result<OuterFits> {
    let mut fits: Vec<Vec<Scored>> = LEARNERS.iter().map(|_| Vec::new()).collect();
    let mut models = FoldSplit::new(scheme);
    let mut seconds = vec![0.0; LEARNERS.len()];
    let mut prevalence = Vec::with_capacity(trains.len());
    let mut candidates = Vec::with_capacity(trains.len());
    // The C3 experiment's screen count, off by default: a joined arm's
    // first-block rows in `x`, when `input_arms` made it.
    let joined = crate::input_arms::joined_rows(dataset);
    let mut screens = Vec::new();
    for (fi, train) in trains.iter().enumerate() {
        prevalence.push(train_prevalence(y, train));
        let (cx, rows) = training_candidates(x, train, min_support);
        candidates.push(cx.rows);
        if let Some(first) = joined {
            screens.push(screen_fold(
                fi, &cx, &rows, first, y, train, tuning, scheme, cfg,
            )?);
        }
        // The C4 experiment, off by default (`screen_arm`): counts and fit
        // times for its own file. Without it these do nothing.
        screen_arm::note_fold(scheme, fi, train.count() as usize, cx.rows);
        for (li, learner) in LEARNERS.iter().enumerate() {
            let t0 = Instant::now();
            screen_arm::label(scheme, Some(fi), learner);
            let clock = screen_arm::start_fit();
            let context = LassoContext {
                dataset,
                scheme,
                outer_fold: Some(fi),
            };
            let (fit, model) = fit_fold(learner, &cx, y, train, tuning, cfg, context)?;
            seconds[li] += t0.elapsed().as_secs_f64();
            screen_arm::end_fit(clock, scheme, fi, learner);
            // The null-fold experiment, off by default, writes only its
            // own file and changes no fit, score or fold model.
            if fit.scores.is_none() && null_refit::enabled() {
                null_refit::record_fold(
                    learner,
                    &cx,
                    y,
                    train,
                    tuning,
                    cfg,
                    fold_seed(cfg, fi),
                    fold_draw_seed(cfg, fi),
                    context,
                )?;
            }
            fits[li].push(fit);
            // The fold fitted on its candidates; the report names
            // dataset rows.
            models.push(model, &rows);
        }
    }
    if joined.is_some() {
        crate::input_arms::write_screen(
            dataset,
            scheme,
            cfg.seed,
            cfg.screened_features,
            &screens,
        )?;
    }
    Ok(OuterFits {
        fits,
        models,
        seconds,
        prevalence,
        candidates,
    })
}

/// The C3 experiment's screen count of one outer fold (`input_arms`): the
/// candidates of a joined arm's first block (rows of `x` below `first`)
/// and of K31, and for each learner that screens, the screen it runs on its
/// own fitting genomes: the outer training genomes for the forest
/// (`fit_rf`), and the fitting part of the tuning holdout drawn at the
/// fold's seed for boosting and the lasso (`fit_gbt`, `fit_lasso`). It
/// changes no fit.
#[allow(clippy::too_many_arguments)]
pub(crate) fn screen_fold(
    fold: usize,
    cx: &BitMatrix,
    rows: &[u32],
    first: usize,
    y: &Bits,
    train: &Bits,
    groups: Option<&[u32]>,
    scheme: &str,
    cfg: &EvalConfig,
) -> Result<crate::input_arms::ScreenFold> {
    let in_first = |f: u32| (rows[f as usize] as usize) < first;
    let candidates_first = (0..rows.len() as u32).filter(|&f| in_first(f)).count();
    let seed = fold_seed(cfg, fold);
    let unit = tuning_unit(scheme);
    let mut learners = Vec::with_capacity(3);
    for (learner, fit) in [
        ("rf", train.clone()),
        (
            "gbt",
            tuning_holdout(train, y, groups, unit, GBT_HOLDOUT, seed)?.0,
        ),
        (
            "lasso",
            tuning_holdout(train, y, groups, unit, LASSO_HOLDOUT, seed)?.0,
        ),
    ] {
        let feats = screen(cx, y, &fit, cfg.screened_features);
        let passed_first = feats.iter().filter(|&&f| in_first(f)).count();
        learners.push(crate::input_arms::ScreenLearner {
            learner,
            fitting: fit.count() as usize,
            screened: feats.len(),
            passed_first,
            passed_k31: feats.len() - passed_first,
        });
    }
    Ok(crate::input_arms::ScreenFold {
        fold,
        train: train.count() as usize,
        candidates_first,
        candidates_k31: rows.len() - candidates_first,
        learners,
    })
}

/// What a supplied outer partition makes of one sample (`evaluate_assigned`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Role {
    /// Scored by this outer fold's fits. It trains in every other fold.
    Test(usize),
    /// Trains in every outer fold and is never scored.
    Train,
    /// In no fold: it neither trains nor is scored.
    Unused,
}

/// One outer fold of a supplied partition, as counts.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct FoldCounts {
    pub train: usize,
    pub train_resistant: usize,
    pub test: usize,
    pub test_resistant: usize,
    /// Candidate features the fold's fits chose from (`training_candidates`).
    pub candidates: usize,
}

/// One split scored on a supplied outer partition.
pub struct AssignedSplit {
    /// `test_fold` and each learner's `scores` are over `tested`, in order.
    pub result: SplitResult,
    /// The scored samples, in sample order.
    pub tested: Vec<usize>,
    pub folds: Vec<FoldCounts>,
    /// The support floor every fold's candidates were chosen at
    /// (`supplied_support_floor`).
    pub min_support: u32,
}

/// The per-fit support floor of a cohort of `samples` samples: a candidate
/// feature is present in, and absent from, at least this many of a fold's
/// training samples (`training_candidates`). A report takes it from its
/// whole cohort, a supplied partition from the samples it uses
/// (`supplied_support_floor`). The C4 floor arm, off by default
/// (`screen_arm::floor`), fixes it in place of max(3, n/200).
pub fn support_floor(samples: usize) -> u32 {
    if let Some(floor) = screen_arm::floor() {
        return floor;
    }
    (samples as u32 / 200).max(3)
}

/// The support floor of a supplied partition: `support_floor` of the samples
/// it tests or trains. Unused samples do not count, so listing more of them
/// cannot change a fit.
pub fn supplied_support_floor(roles: &[Role]) -> u32 {
    support_floor(roles.iter().filter(|r| **r != Role::Unused).count())
}

/// The split names `evaluate` writes, which a supplied partition may not take.
pub const DRAWN_SCHEMES: [&str; 2] = ["random", "lineage"];

/// Scores every learner on an outer partition the caller supplies in place
/// of the two `evaluate` draws: `roles[i]` says whether sample i is tested
/// in one of `folds` outer folds, trains in all of them, or is left out.
///
/// Only the outer partition is supplied. Each fold's candidates, inner folds,
/// tuning holdouts and seeds are `evaluate`'s (`fit_outer_folds`), and tuning
/// keeps each k-mer-identical duplicate group whole, as the random split
/// does (`tuning_groups`). So roles that repeat a split `evaluate` drew give
/// that split's scores.
///
/// Unused samples reach no fit. The support floor is that of the tested and
/// training samples alone (`supplied_support_floor`), a duplicate group keeps
/// only its tested and training members before the tuning groups are made,
/// and `n_groups` counts the lineages of those samples. `x` must hold every
/// feature whose support among them lies inside that floor, as a projection
/// of the cohort at that floor does (`matrix::Dataset::project`). Learners
/// read only their training samples, and only tested samples are scored.
///
/// A partition is refused before any fit when a fold tests no sample, when a
/// fold's test or training samples lack a class, when a fold tests a sample
/// and trains on its k-mer-identical copy, when a fold's training samples
/// give no candidate feature, or when an inner fold or tuning holdout lacks
/// a class. `cfg.folds` is not read: `folds` is the count.
/// Final models and the MIC regression are not fitted.
///
/// It refuses to run while the null-fold experiment is on
/// (`null_refit::enabled`): a record names a fold's test samples as every
/// sample outside its training ones, which a partition with unused samples
/// would make wrong, and the records are read against the drawn splits of a
/// dataset report. `grm learn` refuses `--exp-null-refit` beside
/// `--fold-file` for the same reason.
#[allow(clippy::too_many_arguments)]
pub fn evaluate_assigned(
    dataset: &str,
    scheme: &str,
    x: &BitMatrix,
    y: &Bits,
    groups: &[u32],
    duplicates: &[Vec<usize>],
    roles: &[Role],
    folds: usize,
    cfg: &EvalConfig,
) -> Result<AssignedSplit> {
    ensure!(
        !null_refit::enabled(),
        "the null-fold refit records the folds grm learn draws; it does not run on a supplied partition"
    );
    let min_support = supplied_support_floor(roles);
    check_inputs(x, y, groups, duplicates, min_support, cfg)?;
    let n = x.cols;
    ensure!(roles.len() == n, "{} roles for {n} samples", roles.len());
    // The C4 experiment, off by default (`screen_arm`), as in `evaluate`.
    let arm = screen_arm::install(groups);
    let used = |i: usize| roles[i] != Role::Unused;
    let duplicates: Vec<Vec<usize>> = duplicates
        .iter()
        .map(|group| {
            group
                .iter()
                .copied()
                .filter(|&i| used(i))
                .collect::<Vec<_>>()
        })
        .filter(|group| group.len() > 1)
        .collect();
    let duplicates = duplicates.as_slice();
    let used_groups: Vec<u32> = (0..n).filter(|&i| used(i)).map(|i| groups[i]).collect();
    ensure!(
        !scheme.is_empty() && !DRAWN_SCHEMES.contains(&scheme),
        "a supplied partition needs a name of its own, not {scheme:?}"
    );
    ensure!(folds >= 1, "a supplied partition needs an outer fold");
    ensure!(
        roles
            .iter()
            .all(|r| !matches!(r, Role::Test(f) if *f >= folds)),
        "a sample is tested in an outer fold outside 0..{folds}"
    );
    let tuning = tuning_groups(scheme, groups, duplicates);
    let mut trains = Vec::with_capacity(folds);
    let mut counts = Vec::with_capacity(folds);
    for f in 0..folds {
        let test = Bits::from_indices(n, (0..n).filter(|&i| roles[i] == Role::Test(f)));
        let train = Bits::from_indices(
            n,
            (0..n).filter(|&i| match roles[i] {
                Role::Test(g) => g != f,
                Role::Train => true,
                Role::Unused => false,
            }),
        );
        for (side, mask) in [("test", &test), ("training", &train)] {
            let pos = mask.and(y).count();
            ensure!(
                pos > 0 && pos < mask.count(),
                "{scheme} outer fold {f}'s {side} samples lack one or both classes"
            );
        }
        for group in duplicates {
            let tested = group.iter().find(|&&i| test.contains(i));
            let trained = group.iter().find(|&&i| train.contains(i));
            if let (Some(a), Some(b)) = (tested, trained) {
                anyhow::bail!(
                    "{scheme} outer fold {f} tests sample {a} and trains on sample {b}, its k-mer-identical copy"
                );
            }
        }
        // A fit on no candidate would give every sample one score, and scm,
        // cart and the forest would call it as if it had learned it.
        ensure!(
            any_training_candidate(x, &train, min_support),
            "{scheme} outer fold {f}'s {} training samples give no candidate feature: none is present in at least {min_support} of them and absent from at least {min_support}",
            train.count()
        );
        check_fold_tuning(scheme, f, &train, y, tuning.as_deref(), cfg)?;
        counts.push(FoldCounts {
            train: train.count() as usize,
            train_resistant: train.and(y).count() as usize,
            test: test.count() as usize,
            test_resistant: test.and(y).count() as usize,
            candidates: 0,
        });
        trains.push(train);
    }
    // The fold models are dropped: a fold file result keeps its fields.
    let OuterFits {
        fits,
        models: _,
        seconds,
        prevalence,
        candidates,
    } = fit_outer_folds(
        dataset,
        scheme,
        x,
        y,
        &trains,
        tuning.as_deref(),
        min_support,
        cfg,
    )?;
    for (count, candidates) in counts.iter_mut().zip(candidates) {
        count.candidates = candidates;
    }
    let tested: Vec<usize> = (0..n)
        .filter(|&i| matches!(roles[i], Role::Test(_)))
        .collect();
    let truth: Vec<bool> = tested.iter().map(|&i| y.contains(i)).collect();
    let test_fold: Vec<usize> = tested
        .iter()
        .map(|&i| match roles[i] {
            Role::Test(f) => f,
            _ => unreachable!("a tested sample has an outer fold"),
        })
        .collect();
    let learners = LEARNERS
        .iter()
        .zip(fits)
        .zip(seconds)
        .map(|((name, fits), seconds)| {
            let fits: Vec<Scored> = fits
                .into_iter()
                .map(|s| Scored {
                    scores: s.scores.map(|all| tested.iter().map(|&i| all[i]).collect()),
                    ..s
                })
                .collect();
            pool_folds(name, &truth, &test_fold, &fits, &prevalence, seconds)
                .with_context(|| format!("{scheme} split"))
        })
        .collect::<Result<Vec<LearnerResult>>>()?;
    // The C9 experiment, off by default (`inner_pick`), as in `evaluate`.
    if inner_pick::enabled() {
        inner_pick::record_assigned(
            dataset,
            scheme,
            x,
            y,
            groups,
            tuning.as_deref(),
            &trains,
            &tested,
            &test_fold,
            &learners,
            min_support,
            cfg,
        )?;
    }
    if arm.is_some() {
        let used_samples = (0..n).filter(|&i| used(i)).count();
        screen_arm::finish(
            dataset,
            scheme,
            cfg.seed,
            cfg.screened_features,
            min_support,
            used_samples,
            x.rows,
        )?;
    }
    drop(arm);
    Ok(AssignedSplit {
        result: SplitResult {
            scheme: scheme.to_string(),
            folds,
            n_groups: count_groups(&used_groups),
            test_fold,
            learners,
            unavailable_reason: None,
        },
        tested,
        folds: counts,
        min_support,
    })
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

/// The `k` heaviest features of `imp`, heaviest first, and among equal
/// weights the lower feature first.
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
    let (depth, _) = tree::select_depth(
        x,
        y,
        &all,
        &inner,
        cfg.inner_folds,
        cfg.cart_max_depth,
        cfg.tie_rule,
    );
    let cart = tree::fit_tree(
        x,
        y,
        &all,
        &TreeParams {
            max_depth: depth,
            min_leaf: 1,
            max_features: None,
            seed: 0,
            tie: cfg.tie_rule,
        },
    );
    screen_arm::label("final", None, "lasso");
    let (lasso_scored, lasso_model) = final_lasso(
        dataset,
        x,
        y,
        tuning,
        cfg,
        &mut LassoTrace::default(),
        &mut |_, _| {},
    )?;
    let draws = draw_seed(cfg);
    screen_arm::label("final", None, "rf");
    let (_, rf_importance) = fit_rf(x, y, &all, cfg, draws.wrapping_add(17));
    screen_arm::label("final", None, "gbt");
    let (_, gbt_model) = fit_gbt(
        x,
        y,
        &all,
        tuning,
        "group",
        cfg,
        cfg.seed.wrapping_add(13),
        draws.wrapping_add(13),
    )?;
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
        rf_importance: top_weighted(rf_importance, TOP_WEIGHTED),
        gbt_importance: top_weighted(gbt_model.importance(), TOP_WEIGHTED),
    })
}

/// The final models `final_models` fits on every sample of `x`, kept whole
/// and applied to the columns of `scored`, a matrix whose rows are `x`'s
/// features. It gives the final models and, per learner in `LEARNERS`
/// order, a score for every column of `scored`, or None for a null fit
/// (boosting kept no tree, the lasso kept no weight). Each fit is
/// `final_models`' own call with its own inputs and seeds, on the tuning
/// groups `evaluate` gives them (`final_tuning`), so the models equal
/// `final_models`' output: a report keeps only the forest's and boosting's importances, and
/// here the forest and the boosted model themselves score `scored`.
/// `grm learn --exp-cohort` scores a drawn cohort's final models on the
/// dataset's other genomes with it (`crate::cutoff`).
#[allow(clippy::too_many_arguments)]
pub(crate) fn final_scores(
    dataset: &str,
    x: &BitMatrix,
    y: &Bits,
    groups: &[u32],
    duplicates: &[Vec<usize>],
    min_support: u32,
    cfg: &EvalConfig,
    scored: &BitMatrix,
) -> Result<(FinalModels, LearnerScores)> {
    ensure!(
        scored.rows == x.rows,
        "the scored matrix has {} features where the cohort has {}",
        scored.rows,
        x.rows
    );
    check_inputs(x, y, groups, duplicates, min_support, cfg)?;
    let tuning = final_tuning(groups, duplicates);
    let tuning = tuning.as_deref();
    let all = Bits::ones(x.cols);
    let inner = final_plan(x, y, tuning, min_support, cfg)?;
    let choice = scm::select(x, y, &all, cfg.scm_max_rules, &inner, cfg.inner_folds);
    let model = scm::fit(x, y, &all, choice.kind, choice.p, choice.n_rules, true);
    let importance = scm::rule_importance(&model, x, y, &all);
    let scm_scores: Vec<f64> = {
        let called = model.predict(scored);
        (0..scored.cols)
            .map(|i| called.contains(i) as u8 as f64)
            .collect()
    };
    let (depth, _) = tree::select_depth(
        x,
        y,
        &all,
        &inner,
        cfg.inner_folds,
        cfg.cart_max_depth,
        cfg.tie_rule,
    );
    let cart = tree::fit_tree(
        x,
        y,
        &all,
        &TreeParams {
            max_depth: depth,
            min_leaf: 1,
            max_features: None,
            seed: 0,
            tie: cfg.tie_rule,
        },
    );
    let cart_scores = cart.predict(scored, depth);
    let (lasso_scored, lasso_model) = final_lasso(
        dataset,
        x,
        y,
        tuning,
        cfg,
        &mut LassoTrace::default(),
        &mut |_, _| {},
    )?;
    let lasso = lasso_scores(lasso_model.as_ref(), scored);
    // `fit_rf` on every sample, with the forest kept: its trees read the
    // screened rows, so `scored` is screened on the same features.
    let draws = draw_seed(cfg);
    screen_arm::label("final", None, "rf");
    let feats = screen_arm::screen_fit(x, y, &all, cfg.screened_features);
    let (mut screened, mut scored_screened) = (
        BitMatrix::with_cols(x.cols),
        BitMatrix::with_cols(scored.cols),
    );
    for &f in &feats {
        screened.push_row(x.row(f as usize));
        scored_screened.push_row(scored.row(f as usize));
    }
    let forest = tree::fit_forest(
        &screened,
        y,
        &all,
        &ForestParams {
            seed: draws.wrapping_add(17),
            tie: cfg.tie_rule,
            ..cfg.forest
        },
    );
    let rf_importance = forest
        .importance()
        .into_iter()
        .map(|(f, w)| (feats[f as usize], w))
        .collect();
    let rf = forest.predict(&scored_screened);
    let (_, gbt_model) = fit_gbt(
        x,
        y,
        &all,
        tuning,
        "group",
        cfg,
        cfg.seed.wrapping_add(13),
        draws.wrapping_add(13),
    )?;
    let gbt = gbt_scores(&gbt_model, scored);
    Ok((
        FinalModels {
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
            rf_importance: top_weighted(rf_importance, TOP_WEIGHTED),
            gbt_importance: top_weighted(gbt_model.importance(), TOP_WEIGHTED),
        },
        vec![Some(scm_scores), Some(cart_scores), Some(rf), gbt, lasso],
    ))
}

/// Per learner in `LEARNERS` order, a score for every scored column, or
/// None for a null fit (`final_scores`).
pub(crate) type LearnerScores = Vec<Option<Vec<f64>>>;

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
    let inner = training_folds(y, &all, tuning, cfg.inner_folds, cfg.seed.wrapping_add(7));
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
        cfg.seed.wrapping_add(11),
        LassoContext {
            dataset,
            scheme: "final",
            outer_fold: None,
        },
        trace,
        progress,
    )
}

pub(crate) fn evaluate_mic(
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
        let (fit_m, val) = holdout(&train, y, 0.15, cfg.seed.wrapping_add(fi as u64));
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
                seed: draw_seed(cfg).wrapping_add(fi as u64),
                tie: cfg.tie_rule,
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

    #[test]
    fn a_report_records_the_tie_rule_only_when_the_experiment_is_on() {
        let off = serde_json::to_value(EvalConfig::default()).unwrap();
        assert!(off.get("tie_rule").is_none());
        assert!(off["forest"].get("tie").is_none() && off["gbt"].get("tie").is_none());
        let on = serde_json::to_value(EvalConfig {
            tie_rule: TieRule::Exact,
            ..Default::default()
        })
        .unwrap();
        assert_eq!(on["tie_rule"], "exact");
        assert_eq!(
            off.as_object().unwrap().len() + 1,
            on.as_object().unwrap().len()
        );
    }

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
        let (depth, _) = tree::select_depth(
            &x,
            &y,
            &all,
            &inner,
            cfg.inner_folds,
            cfg.cart_max_depth,
            cfg.tie_rule,
        );
        assert_eq!(out.finals.cart_depth, depth);
        let (_, gbt) = fit_gbt(
            &x,
            &y,
            &all,
            tuning.as_deref(),
            "group",
            &cfg,
            cfg.seed + 13,
            cfg.seed + 13,
        )
        .unwrap();
        assert_eq!(
            serde_json::to_value(&out.finals.gbt_importance).unwrap(),
            serde_json::to_value(top_weighted(gbt.importance(), 50)).unwrap()
        );
        // The final lasso's holdout is grouped too: its audit is the grouped
        // fit's, which on this fixture differs from a holdout drawn by class.
        let lasso_audit = |groups: Option<&[u32]>| {
            let (scored, _) = fit_lasso_traced(
                &x,
                &y,
                &all,
                groups,
                &cfg,
                cfg.seed + 11,
                LassoContext {
                    dataset: "fixture__fixturemycin",
                    scheme: "final",
                    outer_fold: None,
                },
                &mut LassoTrace::default(),
                &mut |_, _| {},
            )
            .unwrap();
            serde_json::to_value(&scored.lasso_fit).unwrap()
        };
        let grouped = lasso_audit(tuning.as_deref());
        assert_ne!(grouped, lasso_audit(None));
        assert_eq!(
            serde_json::to_value(&out.finals.lasso_fit).unwrap(),
            grouped
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
                        fit_scm(&cx, &y, &train, &inner, &cfg).0,
                        fit_cart(&cx, &y, &train, &inner, &cfg).0,
                        fit_gbt(&cx, &y, &train, ids, "group", &cfg, seed, seed)
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
                        .unwrap()
                        .0;
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
                // Each fold's training prevalence, counted here from the
                // saved folds and labels: the score of a null fit.
                let prevalence = |f: usize| {
                    let train: Vec<usize> = (0..n).filter(|&i| split.test_fold[i] != f).collect();
                    train.iter().filter(|&&i| truth[i]).count() as f64 / train.len() as f64
                };
                let all: Vec<usize> = (0..n).collect();
                for learner in &split.learners {
                    let ranked = learner.learner != "scm";
                    assert!(learner.scores.iter().all(|s| s.is_finite()));
                    assert_eq!(
                        learner.metrics,
                        Some(metrics_on(&all, &learner.scores, ranked))
                    );
                    let null = |i: usize| {
                        learner.fold_complexity[split.test_fold[i]] == 0.0
                            && matches!(learner.learner.as_str(), "gbt" | "lasso")
                    };
                    let fitted: Vec<usize> = (0..n).filter(|&i| !null(i)).collect();
                    let in_null =
                        |class: bool| (0..n).filter(|&i| null(i) && truth[i] == class).count();
                    assert_eq!(
                        learner.no_call,
                        Some(NoCall {
                            resistant: in_null(true),
                            susceptible: in_null(false)
                        })
                    );
                    assert_eq!(
                        learner.called_metrics,
                        (fitted.len() < n).then(|| {
                            (!fitted.is_empty())
                                .then(|| metrics_on(&fitted, &learner.scores, ranked))
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
                        .unwrap()
                        .0;
                        assert_eq!(learner.fold_complexity[f], refit.complexity);
                        let idx: Vec<usize> = test.iter_ones().collect();
                        match &refit.scores {
                            None => {
                                assert!(idx.iter().all(|&i| learner.scores[i] == prevalence(f)))
                            }
                            Some(s) => assert!(idx.iter().all(|&i| learner.scores[i] == s[i])),
                        }
                        assert_eq!(
                            learner.fold_balanced_accuracy[f],
                            metrics_on(&idx, &learner.scores, false).balanced_accuracy
                        );
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

    fn fitted(scores: &[f64]) -> Scored {
        Scored {
            scores: Some(scores.to_vec()),
            complexity: 3.0,
            lasso_fit: None,
        }
    }

    fn null_fit() -> Scored {
        Scored {
            scores: None,
            complexity: 0.0,
            lasso_fit: None,
        }
    }

    /// Each fold's training prevalence, as `evaluate` counts it.
    fn prevalences(truth: &[bool], test_fold: &[usize], k: usize) -> Vec<f64> {
        let y = Bits::from_indices(truth.len(), (0..truth.len()).filter(|&i| truth[i]));
        fold_masks(test_fold, k)
            .iter()
            .map(|(train, _)| train_prevalence(&y, train))
            .collect()
    }

    #[test]
    fn one_null_fold_of_several_is_scored_at_its_training_prevalence() {
        // Eight samples in three folds; fold 1's fit is null. Its training
        // samples are 0, 1, 2, 5, 6, 7, two of them resistant, so samples 3
        // and 4 score 2/6 = 1/3 and are called susceptible. The test fold's
        // own prevalence (1/2) would call them resistant, and the whole
        // dataset's is 3/8. Each fit also scores samples outside its test
        // fold, which must not be read.
        let truth = [true, true, false, true, false, false, false, false];
        let test_fold = [0, 0, 0, 1, 1, 2, 2, 2];
        let prevalence = prevalences(&truth, &test_fold, 3);
        assert_eq!(prevalence, [1.0 / 5.0, 2.0 / 6.0, 3.0 / 5.0]);
        let fits = [
            fitted(&[0.9, 0.4, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0]),
            null_fit(),
            fitted(&[0.0, 0.0, 0.0, 0.0, 0.0, 0.6, 0.1, 0.3]),
        ];
        let l = pool_folds("gbt", &truth, &test_fold, &fits, &prevalence, 1.0).unwrap();
        let third = 2.0 / 6.0;
        assert_eq!(l.scores, [0.9, 0.4, 0.2, third, third, 0.6, 0.1, 0.3]);
        // Every sample: tp {0}, fn {1, 3}, fp {5}, tn {2, 4, 6, 7}.
        let m = l.metrics.as_ref().unwrap();
        assert_eq!((m.n, m.tp, m.fn_, m.fp, m.tn), (8, 1, 2, 1, 4));
        assert_eq!(m.sensitivity, 1.0 / 3.0);
        assert_eq!(m.specificity, 0.8);
        assert_eq!(m.balanced_accuracy, (1.0 / 3.0 + 0.8) / 2.0);
        // Resistant scores 0.9, 0.4, 1/3 against susceptible 0.2, 1/3, 0.6,
        // 0.1, 0.3: 5 + 4 + 3.5 of 15 pairs, the tie at 1/3 counting half.
        assert!((m.auroc.unwrap() - 12.5 / 15.0).abs() < 1e-12);
        assert_eq!(l.roc, roc_curve(&truth, &l.scores, 60));
        // Fold 0: one of two resistant called, its susceptible sample right.
        // Fold 1: the constant susceptible call, 0 and 1. Fold 2: no
        // resistant sample, so sensitivity 0, and 2 of 3 susceptible right.
        assert_eq!(
            l.fold_balanced_accuracy,
            [0.75, 0.5, (0.0 + 2.0 / 3.0) / 2.0]
        );
        assert_eq!(
            l.no_call,
            Some(NoCall {
                resistant: 1,
                susceptible: 1
            })
        );
        // The figure revision 13 published: samples 3 and 4 left out.
        let c = l.called_metrics.as_ref().unwrap().as_ref().unwrap();
        assert_eq!((c.n, c.tp, c.fn_, c.fp, c.tn), (6, 1, 1, 1, 3));
        assert_eq!(c.balanced_accuracy, (0.5 + 0.75) / 2.0);
        assert!((c.auroc.unwrap() - 7.0 / 8.0).abs() < 1e-12);
        assert_eq!(l.fold_complexity, [3.0, 0.0, 3.0]);
        // Saved and read back, the row keeps both figures.
        let saved = serde_json::to_value(&l).unwrap();
        assert_eq!(saved["called_metrics"]["n"], 6);
        assert_eq!(saved["metrics"]["n"], 8);
        assert!(
            saved["scores"]
                .as_array()
                .unwrap()
                .iter()
                .all(|s| s.is_number())
        );
        let back: LearnerResult = serde_json::from_value(saved).unwrap();
        assert_eq!(back.called_metrics, l.called_metrics);
        assert_eq!(back.metrics, l.metrics);
        assert_eq!(back.scores, l.scores);
    }

    #[test]
    fn a_training_prevalence_of_exactly_one_half_is_a_resistant_call() {
        // Fold 1's training samples are 0 and 1, one resistant: its null fit
        // scores 0.5, which the 0.5 cutoff calls resistant.
        let truth = [true, false, true, false];
        let test_fold = [0, 0, 1, 1];
        let prevalence = prevalences(&truth, &test_fold, 2);
        assert_eq!(prevalence, [0.5, 0.5]);
        let fits = [fitted(&[0.8, 0.3, 0.0, 0.0]), null_fit()];
        let l = pool_folds("gbt", &truth, &test_fold, &fits, &prevalence, 1.0).unwrap();
        assert_eq!(l.scores, [0.8, 0.3, 0.5, 0.5]);
        let m = l.metrics.as_ref().unwrap();
        assert_eq!((m.n, m.tp, m.fn_, m.fp, m.tn), (4, 2, 0, 1, 1));
        assert_eq!(m.balanced_accuracy, 0.75);
        // The constant resistant call: every resistant right, every
        // susceptible wrong.
        assert_eq!(l.fold_balanced_accuracy, [1.0, 0.5]);
        let c = l.called_metrics.as_ref().unwrap().as_ref().unwrap();
        assert_eq!((c.n, c.tp, c.fp, c.tn), (2, 1, 0, 1));
    }

    #[test]
    fn a_row_without_a_null_fit_has_no_called_metrics_and_the_same_numbers() {
        // Every fold fitted: the metrics are those of every pooled score, as
        // before revision 14, and the row has no `called_metrics` field.
        let truth = [true, false, true, false, true, false];
        let test_fold = [0, 0, 1, 1, 2, 2];
        let fits = [
            fitted(&[0.9, 0.2, 0.0, 0.0, 0.0, 0.0]),
            fitted(&[0.0, 0.0, 0.4, 0.7, 0.0, 0.0]),
            fitted(&[0.0, 0.0, 0.0, 0.0, 0.5, 0.6]),
        ];
        let prevalence = prevalences(&truth, &test_fold, 3);
        let l = pool_folds("rf", &truth, &test_fold, &fits, &prevalence, 1.0).unwrap();
        let scores = [0.9, 0.2, 0.4, 0.7, 0.5, 0.6];
        assert_eq!(l.scores, scores);
        assert_eq!(l.metrics, Some(binary_metrics(&truth, &scores, 0.5, true)));
        let m = l.metrics.as_ref().unwrap();
        assert_eq!((m.n, m.tp, m.fn_, m.fp, m.tn), (6, 2, 1, 2, 1));
        assert_eq!(l.called_metrics, None);
        assert_eq!(
            l.no_call,
            Some(NoCall {
                resistant: 0,
                susceptible: 0
            })
        );
        assert_eq!(l.fold_balanced_accuracy, [1.0, 0.0, 0.5]);
        let saved = serde_json::to_value(&l).unwrap();
        assert!(saved.get("called_metrics").is_none());
        // The prevalences are not read when no fit is null.
        let other = pool_folds("rf", &truth, &test_fold, &fits, &[0.0, 1.0, 0.5], 1.0).unwrap();
        assert_eq!(serde_json::to_value(&other).unwrap(), saved);
    }

    #[test]
    fn pooling_refuses_a_missing_or_impossible_training_prevalence() {
        let truth = [true, false];
        let fits = [null_fit(), fitted(&[0.2, 0.7])];
        for bad in [&[0.5][..], &[0.5, f64::NAN], &[1.5, 0.5]] {
            assert!(pool_folds("gbt", &truth, &[0, 1], &fits, bad, 1.0).is_err());
        }
    }

    #[test]
    fn null_fits_are_scored_at_their_folds_training_prevalence_in_evaluate() {
        // With no candidate feature, boosting keeps no tree and the lasso
        // path is empty in every fold: null fits. CART's one leaf and the
        // forest's leaf-only trees score exactly 0.5 too, but they are
        // fitted models, so the 0.5 cutoff calls every sample resistant.
        // Three extra resistant samples (1, 5, 9) make the folds' training
        // prevalences differ from one another, from their test folds' and
        // from the whole dataset's.
        let n = 80;
        let y = Bits::from_indices(n, (0..n).filter(|i| i % 4 == 0 || [1, 5, 9].contains(i)));
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
        let truth: Vec<bool> = (0..n).map(|i| y.contains(i)).collect();
        let none = NoCall {
            resistant: 0,
            susceptible: 0,
        };
        let mut distinct = std::collections::BTreeSet::new();
        for split in &out.splits {
            assert!(split.unavailable_reason.is_none());
            let f = &split.test_fold;
            let share = |keep: &dyn Fn(usize) -> bool| {
                let idx: Vec<usize> = (0..n).filter(|&i| keep(i)).collect();
                idx.iter().filter(|&&i| truth[i]).count() as f64 / idx.len() as f64
            };
            let train: Vec<f64> = (0..cfg.folds).map(|k| share(&|i| f[i] != k)).collect();
            let test: Vec<f64> = (0..cfg.folds).map(|k| share(&|i| f[i] == k)).collect();
            assert!(
                (0..cfg.folds).any(|k| train[k] != test[k]),
                "{}",
                split.scheme
            );
            assert!(train.iter().any(|&p| p != 23.0 / 80.0), "{}", split.scheme);
            distinct.extend(train.iter().map(|p| p.to_bits()));
            for l in &split.learners {
                match l.learner.as_str() {
                    "gbt" | "lasso" => {
                        let want: Vec<f64> = (0..n).map(|i| train[f[i]]).collect();
                        assert_eq!(l.scores, want, "{} {}", split.scheme, l.learner);
                        // Every prevalence is below 0.5: all called susceptible.
                        let m = l.metrics.as_ref().unwrap();
                        assert_eq!((m.n, m.tp, m.fp, m.tn, m.fn_), (n, 0, 0, 57, 23));
                        assert_eq!(m.balanced_accuracy, 0.5);
                        assert_eq!(m.auroc, crate::learn::auroc(&truth, &want));
                        assert_eq!(l.roc, roc_curve(&truth, &want, 60));
                        assert_eq!(
                            l.no_call,
                            Some(NoCall {
                                resistant: 23,
                                susceptible: 57
                            })
                        );
                        assert_eq!(l.called_metrics, Some(None));
                        assert_eq!(l.fold_balanced_accuracy, vec![0.5; cfg.folds]);
                        assert_eq!(l.fold_complexity, vec![0.0; cfg.folds]);
                    }
                    "cart" | "rf" => {
                        assert!(l.scores.iter().all(|&s| s == 0.5));
                        let m = l.metrics.as_ref().unwrap();
                        assert_eq!((m.n, m.tp, m.fp), (n, 23, 57));
                        assert_eq!(l.no_call, Some(none));
                        assert_eq!(l.called_metrics, None);
                    }
                    _ => {
                        assert_eq!(l.no_call, Some(none));
                        assert_eq!(l.called_metrics, None);
                    }
                }
            }
        }
        assert!(distinct.len() > 1);
        // Saved, the scores are numbers and `called_metrics` is null.
        let gbt = &out.splits[0].learners[3];
        let saved = serde_json::to_value(gbt).unwrap();
        assert!(
            saved["scores"]
                .as_array()
                .unwrap()
                .iter()
                .all(|s| s.is_number())
        );
        assert!(saved["called_metrics"].is_null() && saved.get("called_metrics").is_some());
        let back: LearnerResult = serde_json::from_value(saved).unwrap();
        assert_eq!(back.called_metrics, Some(None));
        assert_eq!(back.scores, gbt.scores);
        assert_eq!(back.fold_balanced_accuracy, gbt.fold_balanced_accuracy);
    }

    #[test]
    fn a_learner_with_every_fold_null_predicts_each_folds_prevalence() {
        // E. coli amikacin in revision 5: 1,012 isolates, 30 resistant, and
        // lineage boosting kept no tree in any of its 5 folds. Revision 5
        // scored those isolates 0.5 (all resistant, 0% VME); revisions 7 to
        // 13 left them uncalled with no rates. A null fit predicts its
        // training prevalence: folds 0 and 1 hold 203 isolates, so their
        // training sets are 24 of 809, the others' 24 of 810. All are below
        // 0.5, so every isolate is called susceptible.
        let truth: Vec<bool> = (0..1012).map(|i| i < 30).collect();
        let test_fold: Vec<usize> = (0..1012).map(|i| i % 5).collect();
        let fits: Vec<Scored> = (0..5).map(|_| null_fit()).collect();
        let prevalence = prevalences(&truth, &test_fold, 5);
        let (high, low) = (24.0 / 809.0, 24.0 / 810.0);
        assert_eq!(prevalence, [high, high, low, low, low]);
        let l = pool_folds("gbt", &truth, &test_fold, &fits, &prevalence, 1.0).unwrap();
        let m = l.metrics.as_ref().unwrap();
        assert_eq!((m.n, m.tp, m.fn_, m.fp, m.tn), (1012, 0, 30, 0, 982));
        assert_eq!((m.sensitivity, m.specificity, m.vme), (0.0, 1.0, 1.0));
        assert_eq!(m.balanced_accuracy, 0.5);
        // 12 resistant and 394 susceptible score high, 18 and 588 low:
        // (12 * 588 + (12 * 394 + 18 * 588) / 2) / (30 * 982).
        assert!((m.auroc.unwrap() - 14712.0 / 29460.0).abs() < 1e-12);
        assert!(!l.roc.is_empty());
        assert_eq!(l.fold_balanced_accuracy, vec![0.5; 5]);
        assert_eq!(
            l.no_call,
            Some(NoCall {
                resistant: 30,
                susceptible: 982
            })
        );
        assert_eq!(l.called_metrics, Some(None));
        let saved = serde_json::to_value(&l).unwrap();
        assert_eq!(saved["metrics"]["n"], 1012);
        assert!(saved["called_metrics"].is_null());
        for key in ["scores", "fold_balanced_accuracy"] {
            assert!(saved[key].as_array().unwrap().iter().all(|v| v.is_number()));
        }
        let back: LearnerResult = serde_json::from_value(saved).unwrap();
        assert_eq!(back.metrics, l.metrics);
        assert_eq!(back.called_metrics, Some(None));
    }

    #[test]
    fn rows_written_before_revision_14_read_back_unchanged() {
        // A revision-13 row: fold 1's samples uncalled (null scores, null
        // fold balanced accuracy), metrics over the called samples, and no
        // `called_metrics` field.
        let old = serde_json::json!({
            "learner": "gbt",
            "metrics": {"n": 4, "tp": 2, "fp": 1, "tn": 1, "fn": 0, "accuracy": 0.75,
                        "balanced_accuracy": 0.75, "sensitivity": 1.0, "specificity": 0.5,
                        "precision": 0.6666666666666666, "f1": 0.8, "mcc": 0.5773502691896258,
                        "vme": 0.0, "me": 0.5, "auroc": 0.75, "auprc": 0.8333333333333333},
            "fold_balanced_accuracy": [1.0, null, 0.5],
            "roc": [[0.0, 0.0], [1.0, 1.0]],
            "scores": [0.9, 0.2, null, null, 0.5, 0.6],
            "conformal": [],
            "fold_complexity": [3.0, 0.0, 3.0],
            "seconds": 1.0,
            "no_call": {"resistant": 1, "susceptible": 1}
        });
        let back: LearnerResult = serde_json::from_value(old.clone()).unwrap();
        assert_eq!(back.called_metrics, None);
        assert_eq!(back.metrics.as_ref().unwrap().n, 4);
        assert!(back.scores[2].is_nan() && back.fold_balanced_accuracy[1].is_nan());
        assert_eq!(serde_json::to_value(&back).unwrap(), old);
        // A row with null metrics (every fold null, revisions 7 to 13).
        let mut none = old;
        none["metrics"] = serde_json::Value::Null;
        let back: LearnerResult = serde_json::from_value(none.clone()).unwrap();
        assert_eq!((back.metrics, back.called_metrics), (None, None));
        let again =
            serde_json::to_value(serde_json::from_value::<LearnerResult>(none.clone()).unwrap())
                .unwrap();
        assert_eq!(again, none);
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
                .unwrap()
                .0;
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
        let (scored, model) =
            fit_gbt(&x, &y, &all, None, "group", &small_config(), 42, 42).unwrap();
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
            .0
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
        // Not set, the learner seed is left out, so a release report's
        // `eval_config` keeps its bytes; set, it is recorded.
        assert!(saved.get("learner_seed").is_none());
        let twin = EvalConfig {
            learner_seed: Some(7),
            ..Default::default()
        };
        assert_eq!(serde_json::to_value(twin).unwrap()["learner_seed"], 7);
    }

    /// Synthetic: 120 genomes in 24 lineages of five, a third resistant,
    /// two features that follow the label with noise and four of noise.
    fn draw_fixture() -> (BitMatrix, Bits, Vec<u32>) {
        use rand::{Rng, SeedableRng};
        let n = 120;
        let y = Bits::from_indices(n, (0..n).filter(|i| i % 3 == 0));
        let groups: Vec<u32> = (0..n as u32).map(|i| i / 5).collect();
        let mut rng = rand_chacha::ChaCha8Rng::seed_from_u64(20261114);
        let mut x = BitMatrix::with_cols(n);
        for f in 0..6 {
            let present: Vec<usize> = (0..n)
                .filter(|&i| {
                    if f < 2 {
                        y.contains(i) != rng.random_bool(0.15)
                    } else {
                        rng.random_bool(0.4)
                    }
                })
                .collect();
            x.push_row(Bits::from_indices(n, present).words());
        }
        (x, y, groups)
    }

    fn draw_config() -> EvalConfig {
        EvalConfig {
            folds: 3,
            inner_folds: 2,
            forest: ForestParams {
                n_trees: 8,
                ..Default::default()
            },
            gbt: GbtParams {
                rounds: 12,
                ..Default::default()
            },
            screened_features: 6,
            lasso_lambdas: 6,
            ..Default::default()
        }
    }

    fn learner_scores(out: &EvalOutput, scheme: &str, learner: &str) -> Vec<f64> {
        let split = out.splits.iter().find(|s| s.scheme == scheme).unwrap();
        let row = split.learners.iter().find(|l| l.learner == learner);
        row.unwrap().scores.clone()
    }

    fn test_folds(out: &EvalOutput) -> Vec<Vec<usize>> {
        out.splits.iter().map(|s| s.test_fold.clone()).collect()
    }

    /// Another evaluation seed is another fold draw of the same dataset:
    /// both splits' outer folds are the ones `fold_plan` draws for that
    /// seed, they differ from the release seed's, and nothing about the
    /// dataset (its samples, labels and lineage count) moves with it.
    #[test]
    fn another_seed_redraws_both_splits_outer_folds_on_the_same_dataset() {
        use crate::ablate::plan::{GroupIds, Scheme, fold_plan};
        let (x, y, groups) = draw_fixture();
        let run = |seed| {
            let cfg = EvalConfig {
                seed,
                ..draw_config()
            };
            evaluate("fixture__fixturemycin", &x, &y, &groups, &[], 1, None, &cfg).unwrap()
        };
        let (base, other) = (run(42), run(43));
        assert_eq!(base.splits.len(), 2);
        assert!(base.splits.iter().all(|s| s.unavailable_reason.is_none()));
        assert!(other.splits.iter().all(|s| s.unavailable_reason.is_none()));
        for (a, b) in base.splits.iter().zip(&other.splits) {
            assert_eq!((&a.scheme, a.n_groups), (&b.scheme, b.n_groups));
            assert_eq!(a.test_fold.len(), x.cols);
            assert_ne!(a.test_fold, b.test_fold, "{}", a.scheme);
        }
        let ids = GroupIds::new(groups.clone(), vec![], None).unwrap();
        for (seed, out) in [(42, &base), (43, &other)] {
            for (scheme, split) in [Scheme::Random, Scheme::Lineage]
                .into_iter()
                .zip(&out.splits)
            {
                let plan = fold_plan(&y, &ids, scheme, 3, 2, seed, 0).unwrap();
                assert_eq!(plan.outer, split.test_fold, "seed {seed} {}", split.scheme);
            }
        }
    }

    /// The learner seed moves the forest's and boosting's own draws and
    /// nothing else: every outer fold, and so every score of the SCM, CART
    /// and the lasso, is the one the evaluation seed alone gives. Set to the
    /// evaluation seed it changes no score at all.
    #[test]
    fn the_learner_seed_moves_the_forest_and_boosting_and_no_partition() {
        let (x, y, groups) = draw_fixture();
        let run = |learner_seed| {
            let cfg = EvalConfig {
                learner_seed,
                ..draw_config()
            };
            evaluate("fixture__fixturemycin", &x, &y, &groups, &[], 1, None, &cfg).unwrap()
        };
        let (base, same, twin) = (run(None), run(Some(42)), run(Some(4242)));
        assert_eq!(test_folds(&base), test_folds(&twin));
        for scheme in ["random", "lineage"] {
            for learner in LEARNERS {
                assert_eq!(
                    learner_scores(&base, scheme, learner),
                    learner_scores(&same, scheme, learner),
                    "{scheme} {learner}"
                );
            }
            for learner in ["scm", "cart", "lasso"] {
                assert_eq!(
                    learner_scores(&base, scheme, learner),
                    learner_scores(&twin, scheme, learner),
                    "{scheme} {learner}"
                );
            }
            assert_ne!(
                learner_scores(&base, scheme, "rf"),
                learner_scores(&twin, scheme, "rf"),
                "{scheme}"
            );
        }
        // Boosting draws column samples only when it has columns to drop;
        // its seed still reaches the fit.
        assert_eq!(
            fold_draw_seed(&draw_config(), 2),
            fold_seed(&draw_config(), 2)
        );
        let cfg = EvalConfig {
            learner_seed: Some(4242),
            ..draw_config()
        };
        assert_eq!(fold_seed(&cfg, 2), 42 + 202);
        assert_eq!(fold_draw_seed(&cfg, 2), 4242 + 202);
    }

    fn value<T: Serialize>(v: &T) -> serde_json::Value {
        serde_json::to_value(v).unwrap()
    }

    /// Synthetic: an exact MIC for each of `draw_fixture`'s genomes, two to
    /// four doublings for a resistant genome and -1 to 1 for the rest, so
    /// the MIC regression runs.
    fn draw_mic(y: &Bits) -> Vec<Option<Mic>> {
        (0..y.len())
            .map(|i| {
                let base = if y.contains(i) { 2 } else { -1 };
                Some(Mic {
                    log2: base + (i % 3) as i32,
                    censor: Censor::Exact,
                })
            })
            .collect()
    }

    /// The learner seed reaches each of the forest's and boosting's draws:
    /// boosting's held-out scores on both splits, the final forest's and
    /// boosting's importance lists, and the MIC regression's boosting. SCM's,
    /// CART's and the lasso's final models do not move. Every partition still
    /// comes from `seed`: boosting rebuilt with the evaluation seed's
    /// early-stopping holdout and the learner seed's column samples gives
    /// the fits `evaluate` made, in each outer fold and in the final model,
    /// and a holdout drawn from the learner seed would not.
    /// The MIC regression's boosted models take the tie rule too: on this
    /// synthetic fixture (`wide_fixture(120, 40, 1)`, `draw_mic`) their
    /// out-of-fold predictions differ between the two rules.
    #[test]
    fn the_tie_rule_reaches_the_mic_regression() {
        let (x, y, groups) = wide_fixture(120, 40, 1);
        let mic = draw_mic(&y);
        let predicted = |tie_rule| {
            let cfg = EvalConfig {
                tie_rule,
                screened_features: 40,
                ..draw_config()
            };
            let input = MicInput {
                mic: &mic,
                breakpoint: None,
            };
            let out = evaluate(
                "fixture__fixturemycin",
                &x,
                &y,
                &groups,
                &[],
                1,
                Some(input),
                &cfg,
            )
            .unwrap();
            value(&out.mic.unwrap().predicted)
        };
        assert_ne!(predicted(TieRule::Exact), predicted(TieRule::Float));
    }

    /// `draw_fixture` with `features` features, all but the first two
    /// present at random on half the samples: with many candidates, small
    /// boosting nodes meet features that divide them into the same groups.
    fn wide_fixture(n: usize, features: usize, seed: u64) -> (BitMatrix, Bits, Vec<u32>) {
        use rand::{Rng, SeedableRng};
        let y = Bits::from_indices(n, (0..n).filter(|i| i % 3 == 0));
        let groups: Vec<u32> = (0..n as u32).map(|i| i / 5).collect();
        let mut rng = rand_chacha::ChaCha8Rng::seed_from_u64(seed);
        let mut x = BitMatrix::with_cols(n);
        for f in 0..features {
            let present: Vec<usize> = (0..n)
                .filter(|&i| {
                    if f < 2 {
                        y.contains(i) != rng.random_bool(0.15)
                    } else {
                        rng.random_bool(0.5)
                    }
                })
                .collect();
            x.push_row(Bits::from_indices(n, present).words());
        }
        (x, y, groups)
    }

    #[test]
    fn the_learner_seed_reaches_boosting_the_final_models_and_the_mic_regression() {
        use crate::ablate::plan::fold_seed as seed_of_fold;
        let (x, y, groups) = draw_fixture();
        let mic = draw_mic(&y);
        let run = |learner_seed| {
            let cfg = EvalConfig {
                learner_seed,
                ..draw_config()
            };
            let input = MicInput {
                mic: &mic,
                breakpoint: None,
            };
            evaluate(
                "fixture__fixturemycin",
                &x,
                &y,
                &groups,
                &[],
                1,
                Some(input),
                &cfg,
            )
            .unwrap()
        };
        let (base, twin) = (run(None), run(Some(4242)));
        for scheme in ["random", "lineage"] {
            assert_ne!(
                learner_scores(&base, scheme, "gbt"),
                learner_scores(&twin, scheme, "gbt"),
                "{scheme}"
            );
        }
        let (a, b) = (&base.finals, &twin.finals);
        assert_ne!(value(&a.rf_importance), value(&b.rf_importance));
        assert_ne!(value(&a.gbt_importance), value(&b.gbt_importance));
        let predicted = |out: &EvalOutput| value(&out.mic.as_ref().unwrap().predicted);
        assert_ne!(predicted(&base), predicted(&twin));
        assert_eq!(value(&a.scm), value(&b.scm));
        assert_eq!(
            (value(&a.cart), a.cart_depth),
            (value(&b.cart), b.cart_depth)
        );
        assert_eq!(value(&a.lasso), value(&b.lasso));
        assert_eq!(value(&a.lasso_fit), value(&b.lasso_fit));

        let cfg = EvalConfig {
            learner_seed: Some(4242),
            ..draw_config()
        };
        let mut holdout_matters = false;
        for (scheme, folds) in outer_splits(&y, &groups, &[], &cfg) {
            let tuning = tuning_groups(scheme, &groups, &[]);
            let unit = tuning_unit(scheme);
            for (fi, (train, _)) in fold_masks(&folds, cfg.folds).into_iter().enumerate() {
                let (cx, _) = training_candidates(&x, &train, 1);
                let draws = seed_of_fold(4242, fi);
                let built = |holdout| {
                    let fit = fit_gbt(
                        &cx,
                        &y,
                        &train,
                        tuning.as_deref(),
                        unit,
                        &cfg,
                        holdout,
                        draws,
                    );
                    fit.unwrap().0.scores
                };
                let context = test_context(scheme, fi);
                let fit = fit_fold("gbt", &cx, &y, &train, tuning.as_deref(), &cfg, context);
                let expected = built(seed_of_fold(42, fi));
                assert_eq!(fit.unwrap().0.scores, expected, "{scheme} outer fold {fi}");
                holdout_matters |= built(draws) != expected;
            }
        }
        assert!(holdout_matters);
        let all = Bits::ones(x.cols);
        let importance = |holdout| {
            let fit = fit_gbt(&x, &y, &all, None, "group", &cfg, holdout, 4242 + 13);
            value(&top_weighted(fit.unwrap().1.importance(), 50))
        };
        assert_eq!(value(&b.gbt_importance), importance(42 + 13));
        assert_ne!(importance(42 + 13), importance(4242 + 13));
    }

    /// `--seed` takes any 64-bit value, so every seed sum wraps at 2^64:
    /// the per-fold seeds, the final models' (+7, +11, +13, +17) and the MIC
    /// regression's (+ fold). The release gate builds without overflow
    /// checks, where `+` wraps as well, so this pins the wrapped values; a
    /// debug run of this test is the one that fails on a sum that does not
    /// wrap.
    #[test]
    fn seed_sums_wrap_at_the_largest_seed() {
        let top = EvalConfig {
            seed: u64::MAX,
            ..draw_config()
        };
        assert_eq!(fold_seed(&top, 0), u64::MAX);
        assert_eq!(fold_seed(&top, 1), 100);
        assert_eq!(fold_seed(&top, 4), 403);
        assert_eq!(fold_draw_seed(&top, 1), 100);
        let twin = EvalConfig {
            learner_seed: Some(u64::MAX),
            ..draw_config()
        };
        assert_eq!(fold_seed(&twin, 1), 42 + 101);
        assert_eq!(fold_draw_seed(&twin, 1), 100);
        // Evaluation at the largest seed reaches every sum, the MIC
        // regression's in each of its folds, and its final models are the
        // ones fitted at the wrapped seeds.
        let (x, y, groups) = draw_fixture();
        let mic = draw_mic(&y);
        let input = MicInput {
            mic: &mic,
            breakpoint: None,
        };
        let out = evaluate(
            "fixture__fixturemycin",
            &x,
            &y,
            &groups,
            &[],
            1,
            Some(input),
            &top,
        )
        .unwrap();
        assert!(out.mic.is_some());
        assert!(out.splits.iter().all(|s| s.unavailable_reason.is_none()));
        let all = Bits::ones(x.cols);
        assert_eq!(
            final_plan(&x, &y, None, 1, &top).unwrap(),
            training_folds(&y, &all, None, top.inner_folds, 6)
        );
        let (_, rf) = fit_rf(&x, &y, &all, &top, 16);
        assert_eq!(
            value(&out.finals.rf_importance),
            value(&top_weighted(rf, 50))
        );
        let (_, gbt) = fit_gbt(&x, &y, &all, None, "group", &top, 12, 12).unwrap();
        assert_eq!(
            value(&out.finals.gbt_importance),
            value(&top_weighted(gbt.importance(), 50))
        );
        let (lasso, _) = fit_lasso_traced(
            &x,
            &y,
            &all,
            None,
            &top,
            10,
            LassoContext {
                dataset: "fixture__fixturemycin",
                scheme: "final",
                outer_fold: None,
            },
            &mut LassoTrace::default(),
            &mut |_, _| {},
        )
        .unwrap();
        assert_eq!(value(&out.finals.lasso_fit), value(&lasso.lasso_fit));
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
            called_metrics: None,
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
                    .unwrap()
                    .0;
                    let s = s.scores.expect("every learner fits this matrix");
                    assert!(s.iter().all(|v| v.is_finite()));
                    s[..40].to_vec()
                })
                .collect()
        };
        assert_eq!(fits(&x), fits(&x2));
    }

    /// Synthetic: 120 samples in 30 lineages of four. Every third sample is
    /// resistant, two rows mark resistance imperfectly and three are noise.
    fn assigned_fixture() -> (BitMatrix, Bits, Vec<u32>) {
        let n = 120;
        let rules: [fn(usize) -> bool; 5] = [
            |i| (i % 3 == 0 && i % 10 != 0) || i % 20 == 1,
            |i| i % 3 == 0 && i % 7 != 0,
            |i| i % 5 == 0,
            |i| i % 7 == 1,
            |i| i % 4 < 2,
        ];
        let mut x = BitMatrix::with_cols(n);
        for rule in rules {
            x.push_row(Bits::from_indices(n, (0..n).filter(|&i| rule(i))).words());
        }
        let y = Bits::from_indices(n, (0..n).filter(|i| i % 3 == 0));
        (x, y, (0..n).map(|i| (i / 4) as u32).collect())
    }

    /// `evaluate_assigned` on the synthetic fixture under the name "supplied".
    fn assigned(
        x: &BitMatrix,
        y: &Bits,
        duplicates: &[Vec<usize>],
        roles: &[Role],
        folds: usize,
        cfg: &EvalConfig,
    ) -> Result<AssignedSplit> {
        let groups: Vec<u32> = (0..x.cols).map(|i| (i / 4) as u32).collect();
        evaluate_assigned(
            "fixture__fixturemycin",
            "supplied",
            x,
            y,
            &groups,
            duplicates,
            roles,
            folds,
            cfg,
        )
    }

    #[test]
    fn roles_that_repeat_a_drawn_split_give_that_splits_scores() {
        // The tuning-sensitive synthetic fixture of
        // `scm_cart_boosting_and_lasso_tune_on_the_splits_tuning_groups`:
        // 160 genomes in lineages of four, every other lineage four copies of
        // one genome. Its fits change when tuning splits the copies, so the
        // equality below also holds the supplied path to the random split's
        // tuning groups.
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
        let cfg = EvalConfig {
            forest: ForestParams {
                n_trees: 8,
                ..Default::default()
            },
            ..Default::default()
        };
        let mut tuning_matters = false;
        for duplicates in [&duplicates[..], &[]] {
            let drawn = evaluate(
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
            let random = &drawn.splits[0];
            assert_eq!(random.scheme, "random");
            let roles: Vec<Role> = random.test_fold.iter().map(|&f| Role::Test(f)).collect();
            let supplied = evaluate_assigned(
                "fixture__fixturemycin",
                "supplied",
                &x,
                &y,
                &groups,
                duplicates,
                &roles,
                cfg.folds,
                &cfg,
            )
            .unwrap();
            assert_eq!(supplied.tested, (0..n).collect::<Vec<_>>());
            assert_eq!(supplied.folds.len(), cfg.folds);
            for (f, counts) in supplied.folds.iter().enumerate() {
                let test = random.test_fold.iter().filter(|&&t| t == f).count();
                assert_eq!((counts.test, counts.train), (test, n - test), "fold {f}");
            }
            let renamed = SplitResult {
                scheme: "random".into(),
                ..supplied.result
            };
            assert_eq!(timeless(&renamed), timeless(random));
            // The same outer folds without the copies' tuning groups are
            // other fits, so the supplied path did tune on the groups.
            if !duplicates.is_empty() {
                let untuned = evaluate_assigned(
                    "fixture__fixturemycin",
                    "supplied",
                    &x,
                    &y,
                    &groups,
                    &[],
                    &roles,
                    cfg.folds,
                    &cfg,
                )
                .unwrap();
                tuning_matters = serde_json::to_vec(&untuned.result.learners).unwrap()
                    != serde_json::to_vec(&renamed.learners).unwrap();
            }
        }
        assert!(tuning_matters);
        assert_eq!(
            tuning_groups("supplied", &groups, &duplicates),
            tuning_groups("random", &groups, &duplicates)
        );
        assert_eq!(tuning_unit("supplied"), tuning_unit("random"));
    }

    #[test]
    fn the_features_and_labels_of_unused_samples_and_the_labels_of_tested_ones_cannot_change_a_folds_scores()
     {
        // One fold: samples 0..30 are tested, 30..90 train, 90..120 are in
        // no fold.
        let (x, y, _) = assigned_fixture();
        let n = x.cols;
        let roles: Vec<Role> = (0..n)
            .map(|i| match i {
                0..30 => Role::Test(0),
                30..90 => Role::Train,
                _ => Role::Unused,
            })
            .collect();
        let cfg = small_config();
        let base = assigned(&x, &y, &[], &roles, 1, &cfg).unwrap();
        let train = Bits::from_indices(n, 30..90);
        assert_eq!(
            base.folds,
            [FoldCounts {
                train: 60,
                train_resistant: 20,
                test: 30,
                test_resistant: 10,
                candidates: training_candidates(&x, &train, 3).0.rows,
            }]
        );
        assert!(base.folds[0].candidates > 0);
        assert_eq!(base.tested, (0..30).collect::<Vec<_>>());
        assert_eq!(base.result.test_fold, vec![0; 30]);
        // The lineages of the 90 tested and training samples, 0..90 in
        // lineages of four; the unused ones add none.
        assert_eq!((base.result.folds, base.result.n_groups), (1, 23));
        for learner in &base.result.learners {
            assert_eq!(learner.scores.len(), 30, "{}", learner.learner);
            assert_eq!(learner.metrics.as_ref().unwrap().n, 30);
            assert_eq!(learner.fold_balanced_accuracy.len(), 1);
        }
        // The marker rows fit: every learner has a model and calls both
        // classes, so a changed training set would show in its scores.
        for l in &base.result.learners {
            let m = l.metrics.as_ref().unwrap();
            assert!(
                l.fold_complexity[0] > 0.0 && m.tp > 0 && m.tn > 0,
                "{}",
                l.learner
            );
        }

        // Every feature and the label of each unused sample inverted.
        let mut x2 = BitMatrix::with_cols(n);
        for r in 0..x.rows {
            let row = (0..n).filter(|&i| x.get(r, i) != (i >= 90));
            x2.push_row(Bits::from_indices(n, row).words());
        }
        let y2 = Bits::from_indices(n, (0..n).filter(|&i| y.contains(i) != (i >= 90)));
        let changed = assigned(&x2, &y2, &[], &roles, 1, &cfg).unwrap();
        assert_eq!(timeless(&changed.result), timeless(&base.result));
        assert_eq!(changed.folds, base.folds);
        // The same change does move the scores once those samples train, so
        // the equality above is not a fixture too dull to show a leak.
        let training: Vec<Role> = roles
            .iter()
            .map(|&r| if r == Role::Unused { Role::Train } else { r })
            .collect();
        let leaked = |x: &BitMatrix, y: &Bits| {
            let split = assigned(x, y, &[], &training, 1, &cfg).unwrap();
            assert_eq!(split.folds[0].train, 90);
            timeless(&split.result)
        };
        assert_ne!(leaked(&x, &y), leaked(&x2, &y2));
        assert_ne!(leaked(&x, &y), timeless(&base.result));

        // A test sample's label changes the metrics and never a score.
        let y3 = Bits::from_indices(
            n,
            (0..n).filter(|&i| y.contains(i) != (i < 30 && i % 2 == 0)),
        );
        let relabelled = assigned(&x, &y3, &[], &roles, 1, &cfg).unwrap();
        assert_eq!(relabelled.folds[0].test_resistant, 15);
        let mut metrics_moved = false;
        for (a, b) in relabelled.result.learners.iter().zip(&base.result.learners) {
            assert_eq!(a.scores, b.scores, "{}", a.learner);
            assert_eq!(a.fold_complexity, b.fold_complexity, "{}", a.learner);
            metrics_moved |= a.metrics != b.metrics;
        }
        assert!(metrics_moved);
    }

    #[test]
    fn a_copy_pair_among_unused_samples_cannot_change_a_folds_fits() {
        // One fold: samples 0..30 are tested, 30..90 train, 90..120 are in
        // no fold. Samples 92 and 93 are unused and share a lineage. Naming
        // them a k-mer-identical pair changes no training and no test sample,
        // so no fit may change: a copy group made of unused samples alone
        // must not move the fold onto grouped tuning draws.
        let (x, y, _) = assigned_fixture();
        let n = x.cols;
        let roles: Vec<Role> = (0..n)
            .map(|i| match i {
                0..30 => Role::Test(0),
                30..90 => Role::Train,
                _ => Role::Unused,
            })
            .collect();
        let cfg = EvalConfig {
            forest: ForestParams {
                n_trees: 8,
                ..Default::default()
            },
            ..Default::default()
        };
        let base = assigned(&x, &y, &[], &roles, 1, &cfg).unwrap();
        let copies = assigned(&x, &y, &[vec![92, 93]], &roles, 1, &cfg).unwrap();
        assert_eq!(copies.folds, base.folds);
        let row = |l: &LearnerResult| {
            let mut v = serde_json::to_value(l).unwrap();
            v.as_object_mut().unwrap().remove("seconds");
            v
        };
        let moved: Vec<&str> = base
            .result
            .learners
            .iter()
            .zip(&copies.result.learners)
            .filter(|(a, b)| row(a) != row(b))
            .map(|(a, _)| a.learner.as_str())
            .collect();
        let scores_moved: Vec<&str> = base
            .result
            .learners
            .iter()
            .zip(&copies.result.learners)
            .filter(|(a, b)| a.scores != b.scores)
            .map(|(a, _)| a.learner.as_str())
            .collect();
        assert!(
            moved.is_empty(),
            "an unused copy pair moved the rows of {moved:?} and the scores of {scores_moved:?}"
        );
        // The same pair among training samples does move a fit, so the
        // equality above is not a fixture too dull to show the switch.
        let trained = assigned(&x, &y, &[vec![32, 33]], &roles, 1, &cfg).unwrap();
        assert!(
            base.result
                .learners
                .iter()
                .zip(&trained.result.learners)
                .any(|(a, b)| row(a) != row(b))
        );
    }

    #[test]
    fn each_supplied_fold_chooses_candidates_and_fits_on_its_own_training_samples() {
        // Synthetic: the 120-sample fixture and one more row, present in
        // samples 0..4 only. Fold 0 tests 0..30 and fold 1 tests 30..60, so
        // the row has no support among fold 0's training samples and four
        // among fold 1's.
        let (mut x, y, _) = assigned_fixture();
        let n = x.cols;
        x.push_row(Bits::from_indices(n, 0..4).words());
        let roles: Vec<Role> = (0..n)
            .map(|i| match i {
                0..30 => Role::Test(0),
                30..60 => Role::Test(1),
                _ => Role::Train,
            })
            .collect();
        let cfg = small_config();
        let base = assigned(&x, &y, &[], &roles, 2, &cfg).unwrap();
        let trains = [
            Bits::from_indices(n, 30..n),
            Bits::from_indices(n, (0..30).chain(60..n)),
        ];
        let expected: Vec<usize> = trains
            .iter()
            .map(|train| training_candidates(&x, train, 3).0.rows)
            .collect();
        assert_eq!(expected[1], expected[0] + 1);
        let counted: Vec<usize> = base.folds.iter().map(|f| f.candidates).collect();
        assert_eq!(counted, expected);
        assert_eq!(base.tested, (0..60).collect::<Vec<_>>());

        // Fold 0's test samples train fold 1 and nothing else. With half of
        // them relabelled, fold 0 scores them as before, and fold 1, whose
        // training samples they are, does not.
        let y2 = Bits::from_indices(
            n,
            (0..n).filter(|&i| y.contains(i) != (i < 30 && i % 2 == 0)),
        );
        let relabelled = assigned(&x, &y2, &[], &roles, 2, &cfg).unwrap();
        assert_eq!(
            relabelled.folds[0].train_resistant,
            base.folds[0].train_resistant
        );
        assert_eq!(
            relabelled.folds[1].train_resistant,
            base.folds[1].train_resistant + 5
        );
        let mut fold_1_moved = false;
        for (a, b) in relabelled.result.learners.iter().zip(&base.result.learners) {
            assert_eq!(a.scores[..30], b.scores[..30], "{}", a.learner);
            assert_eq!(a.fold_complexity[0], b.fold_complexity[0], "{}", a.learner);
            fold_1_moved |= a.scores[30..] != b.scores[30..];
        }
        assert!(fold_1_moved);
    }

    #[test]
    fn a_fold_whose_training_samples_give_no_candidate_feature_is_refused() {
        // Synthetic: the 120-sample fixture's labels. Samples 0..30 are
        // tested and the rest train.
        let (fixture, y, _) = assigned_fixture();
        let n = fixture.cols;
        let roles: Vec<Role> = (0..n)
            .map(|i| if i < 30 { Role::Test(0) } else { Role::Train })
            .collect();
        let cfg = small_config();
        // One feature in every sample and one in samples 0..30 and 118..120:
        // all 90 training samples carry the first, and 2 carry the second,
        // under the floor of 3. Neither is a candidate.
        let matrix = |carriers: std::ops::Range<usize>| {
            let mut x = BitMatrix::with_cols(n);
            x.push_row(Bits::ones(n).words());
            x.push_row(Bits::from_indices(n, (0..30).chain(carriers)).words());
            x
        };
        let x = matrix(118..120);
        let train = Bits::from_indices(n, 30..n);
        assert_eq!(training_candidates(&x, &train, 3).0.rows, 0);
        let error = assigned(&x, &y, &[], &roles, 1, &cfg)
            .err()
            .expect("a fold without a candidate feature is refused");
        assert_eq!(
            format!("{error:#}"),
            "supplied outer fold 0's 90 training samples give no candidate feature: none is present in at least 3 of them and absent from at least 3"
        );
        // A third training sample with the second feature makes it a
        // candidate, and the fold is scored.
        let x = matrix(117..120);
        assert_eq!(training_candidates(&x, &train, 3).0.rows, 1);
        let split = assigned(&x, &y, &[], &roles, 1, &cfg).unwrap();
        assert_eq!(split.folds[0].candidates, 1);
        // The fixture's own features give the fold candidates.
        assert!(assigned(&fixture, &y, &[], &roles, 1, &cfg).is_ok());
    }

    #[test]
    fn a_tested_sample_takes_its_own_score_wherever_it_sits_in_the_sample_order() {
        // Synthetic: the 120-sample fixture and one more row that is the
        // label itself, so a one-rule SCM scores every sample by its own
        // label. Samples 80..111 are tested: they are not the first samples
        // of the cohort, their labels are not those of samples 0..31, and
        // they do not read the same in reverse.
        let (mut x, y, _) = assigned_fixture();
        let n = x.cols;
        x.push_row(y.words());
        let roles: Vec<Role> = (0..n)
            .map(|i| match i {
                80..111 => Role::Test(0),
                20..80 => Role::Train,
                _ => Role::Unused,
            })
            .collect();
        let split = assigned(&x, &y, &[], &roles, 1, &small_config()).unwrap();
        assert_eq!(split.tested, (80..111).collect::<Vec<_>>());
        let label = |i: usize| if i.is_multiple_of(3) { 1.0 } else { 0.0 };
        let own: Vec<f64> = (80..111).map(label).collect();
        let first: Vec<f64> = (0..31).map(label).collect();
        let reversed: Vec<f64> = own.iter().rev().copied().collect();
        assert_ne!(own, first);
        assert_ne!(own, reversed);
        let scm = &split.result.learners[0];
        assert_eq!(scm.learner, "scm");
        assert_eq!(scm.scores, own);
        let m = scm.metrics.as_ref().unwrap();
        assert_eq!((m.tp, m.tn, m.fp, m.fn_), (10, 21, 0, 0));
    }

    #[test]
    fn a_supplied_partition_that_cannot_be_scored_fairly_is_refused() {
        let (x, y, _) = assigned_fixture();
        let n = x.cols;
        let cfg = small_config();
        let refusal = |duplicates: &[Vec<usize>], roles: &[Role], folds: usize| {
            let error = assigned(&x, &y, duplicates, roles, folds, &cfg)
                .err()
                .expect("the partition is refused");
            format!("{error:#}")
        };
        let roles = |role: fn(usize) -> Role| (0..n).map(role).collect::<Vec<Role>>();
        let good = roles(|i| if i < 30 { Role::Test(0) } else { Role::Train });
        assert!(assigned(&x, &y, &[], &good, 1, &cfg).is_ok());

        // Its shape.
        assert_eq!(
            refusal(&[], &good, 0),
            "a supplied partition needs an outer fold"
        );
        assert_eq!(refusal(&[], &good[1..], 1), "119 roles for 120 samples");
        let mut outside = good.clone();
        outside[0] = Role::Test(1);
        assert_eq!(
            refusal(&[], &outside, 1),
            "a sample is tested in an outer fold outside 0..1"
        );
        let groups: Vec<u32> = (0..n).map(|i| (i / 4) as u32).collect();
        for name in ["", "random", "lineage"] {
            let error = evaluate_assigned(
                "fixture__fixturemycin",
                name,
                &x,
                &y,
                &groups,
                &[],
                &good,
                1,
                &cfg,
            )
            .err()
            .expect("a drawn split's name is refused");
            assert_eq!(
                format!("{error:#}"),
                format!("a supplied partition needs a name of its own, not {name:?}")
            );
        }

        // A class missing from a fold's test or training samples.
        let susceptible_test = roles(|i| match i {
            i if i < 30 && i % 3 != 0 => Role::Test(0),
            _ => Role::Train,
        });
        assert_eq!(
            refusal(&[], &susceptible_test, 1),
            "supplied outer fold 0's test samples lack one or both classes"
        );
        let resistant_training = roles(|i| match i {
            0..30 => Role::Test(0),
            i if i % 3 == 0 => Role::Train,
            _ => Role::Unused,
        });
        assert_eq!(
            refusal(&[], &resistant_training, 1),
            "supplied outer fold 0's training samples lack one or both classes"
        );
        // A fold the roles test nobody in.
        assert_eq!(
            refusal(&[], &good, 2),
            "supplied outer fold 1's test samples lack one or both classes"
        );
        // One resistant training sample cannot fill two inner folds.
        let one_resistant = roles(|i| match i {
            0..30 => Role::Test(0),
            i if i % 3 != 0 || i == 30 => Role::Train,
            _ => Role::Unused,
        });
        assert!(
            refusal(&[], &one_resistant, 1)
                .starts_with("invalid supplied inner evaluation in outer fold 0"),
            "{}",
            refusal(&[], &one_resistant, 1)
        );

        // Samples 29 and 30 share a lineage and are k-mer-identical copies.
        let copies = [vec![29, 30]];
        let copy_words = "supplied outer fold 0 tests sample 29 and trains on sample 30, its k-mer-identical copy";
        assert_eq!(refusal(&copies, &good, 1), copy_words);
        // The copy in another test fold still trains fold 0.
        let two_folds = roles(|i| match i {
            0..30 => Role::Test(0),
            30..60 => Role::Test(1),
            _ => Role::Train,
        });
        assert_eq!(refusal(&copies, &two_folds, 2), copy_words);
        assert!(assigned(&x, &y, &[], &two_folds, 2, &cfg).is_ok());
        // A copy left out of every fold, or beside its copy, is no leak.
        let mut left_out = good.clone();
        left_out[30] = Role::Unused;
        assert!(assigned(&x, &y, &copies, &left_out, 1, &cfg).is_ok());
        let mut together = good.clone();
        together[30] = Role::Test(0);
        assert!(assigned(&x, &y, &copies, &together, 1, &cfg).is_ok());
        assert!(assigned(&x, &y, &[vec![40, 41]], &good, 1, &cfg).is_ok());
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

    /// Synthetic: 60 samples, half resistant. Rows 0 to 3 are each present
    /// in three samples, the support floor, so a fold that holds one of the
    /// three out drops the row and every later candidate's index falls
    /// behind its dataset row. Rows 4 to 6 follow the label with noise, and
    /// the last nine are noise.
    fn shifted_fixture() -> (BitMatrix, Bits, Vec<u32>) {
        use rand::{Rng, SeedableRng};
        let n = 60;
        let y = Bits::from_indices(n, (0..n).filter(|i| i % 2 == 0));
        let mut rng = rand_chacha::ChaCha8Rng::seed_from_u64(11);
        let mut raw = BitMatrix::with_cols(n);
        for k in 0..4 {
            raw.push_row(Bits::from_indices(n, [3 * k + 1, 3 * k + 2, 3 * k + 3]).words());
        }
        for k in 0..12 {
            let rate = |r: bool| match (k < 3, r) {
                (true, true) => 0.8,
                (true, false) => 0.15,
                _ => 0.4,
            };
            let present = (0..n).filter(|&i| rng.random_bool(rate(y.contains(i))));
            raw.push_row(Bits::from_indices(n, present.collect::<Vec<_>>()).words());
        }
        let (x, kept) = training_candidates(&raw, &Bits::ones(n), 3);
        assert_eq!(kept.len(), raw.rows, "every fixture row is a candidate");
        (x, y, (0..n).map(|i| (i % 15) as u32).collect())
    }

    #[test]
    fn fold_models_give_each_folds_saved_scores_under_dataset_row_numbers() {
        let (x, y, groups) = shifted_fixture();
        let n = x.cols;
        let cfg = EvalConfig {
            scm_max_rules: 2,
            cart_max_depth: 2,
            screened_features: 6,
            ..small_config()
        };
        let out = evaluate("fixture__fixturemycin", &x, &y, &groups, &[], 3, None, &cfg).unwrap();
        assert_eq!(out.fold_models.len(), out.splits.len());
        // Features each learner named under a row number that is not its
        // index among the fold's candidates: what a model left in candidate
        // numbering would get wrong.
        let mut shifted: std::collections::BTreeMap<&str, usize> = Default::default();
        for (split, models) in out.splits.iter().zip(&out.fold_models) {
            assert!(split.unavailable_reason.is_none(), "{split:?}");
            assert_eq!(models.scheme, split.scheme);
            let saved = |name: &str| split.learners.iter().find(|l| l.learner == name).unwrap();
            let tuning = tuning_groups(&split.scheme, &groups, &[]);
            let masks = fold_masks(&split.test_fold, cfg.folds);
            for lists in [models.scm.len(), models.cart.len(), models.rf.len()] {
                assert_eq!(lists, cfg.folds);
            }
            assert_eq!(models.gbt.len(), cfg.folds);
            assert_eq!(models.lasso.len(), cfg.folds);
            for (f, (train, test)) in masks.into_iter().enumerate() {
                let held_out: Vec<usize> = test.iter_ones().collect();
                let prevalence = (0..n)
                    .filter(|&i| train.contains(i) && y.contains(i))
                    .count() as f64
                    / train.count() as f64;
                let (cx, rows) = training_candidates(&x, &train, 3);
                assert!(rows.len() < x.rows, "the fold drops a row");
                let mut note = |name: &'static str, feature: u32| {
                    let c = rows.iter().position(|&r| r == feature);
                    let c = c.unwrap_or_else(|| panic!("{name} names row {feature}, no candidate"));
                    *shifted.entry(name).or_default() += usize::from(c as u32 != feature);
                };
                let refit = |name: &str| {
                    fit_fold(
                        name,
                        &cx,
                        &y,
                        &train,
                        tuning.as_deref(),
                        &cfg,
                        test_context(&split.scheme, f),
                    )
                    .unwrap()
                    .1
                };

                // The rules alone, read on the dataset matrix, give the
                // fold's saved calls.
                let (scm, row) = (&models.scm[f], saved("scm"));
                assert_eq!(scm.rules.len() as f64, row.fold_complexity[f]);
                let calls = scm::ScmModel {
                    kind: scm.kind,
                    p: scm.p,
                    rules: scm.rules.clone(),
                    steps: vec![],
                }
                .predict(&x);
                assert!(
                    held_out
                        .iter()
                        .all(|&i| row.scores[i] == calls.contains(i) as u8 as f64)
                );
                let FoldModel::Scm(again) = refit("scm") else {
                    panic!("another model")
                };
                assert_eq!((scm.kind, scm.p), (again.kind, again.p));
                let renamed: Vec<Rule> = again
                    .rules
                    .iter()
                    .map(|r| Rule {
                        feature: rows[r.feature as usize],
                        presence: r.presence,
                    })
                    .collect();
                assert_eq!(scm.rules, renamed);
                scm.rules.iter().for_each(|r| note("scm", r.feature));

                let (cart, row) = (&models.cart[f], saved("cart"));
                let leaves = cart.nodes.iter().filter(|n| n.split.is_none()).count();
                assert_eq!(leaves as f64, row.fold_complexity[f]);
                let scores = cart.predict(&x, u32::MAX);
                assert!(held_out.iter().all(|&i| row.scores[i] == scores[i]));
                let FoldModel::Cart(again) = refit("cart") else {
                    panic!("another model")
                };
                let split_rows =
                    |t: &tree::Tree, rename: &dyn Fn(u32) -> u32| -> Vec<Option<u32>> {
                        t.nodes
                            .iter()
                            .map(|n| n.split.as_ref().map(|s| rename(s.feature)))
                            .collect()
                    };
                assert_eq!(
                    split_rows(cart, &|f| f),
                    split_rows(&again, &|c| rows[c as usize])
                );
                for feature in split_rows(cart, &|f| f).into_iter().flatten() {
                    note("cart", feature);
                }

                let (model, row) = (&models.lasso[f], saved("lasso"));
                let weights = model.as_ref().map_or(0, |m| m.coef.len());
                assert_eq!(weights as f64, row.fold_complexity[f]);
                match lasso_scores(model.as_ref(), &x) {
                    Some(scores) => assert!(held_out.iter().all(|&i| row.scores[i] == scores[i])),
                    None => assert!(held_out.iter().all(|&i| row.scores[i] == prevalence)),
                }
                let FoldModel::Lasso(again) = refit("lasso") else {
                    panic!("another model")
                };
                assert_eq!(
                    model
                        .as_ref()
                        .map(|m| (m.intercept, m.lambda, m.coef.clone())),
                    again.map(|m| {
                        let coef = m.coef.iter().map(|&(c, w)| (rows[c as usize], w));
                        (m.intercept, m.lambda, coef.collect())
                    })
                );
                for &(feature, _) in model.iter().flat_map(|m| &m.coef) {
                    note("lasso", feature);
                }

                // The forest and boosting keep a ranking: the importance of
                // a fit on the fold's candidates, keyed by dataset row.
                let seed = fold_seed(&cfg, f);
                let draws = fold_draw_seed(&cfg, f);
                let pairs = |list: &[WeightedFeature]| -> Vec<(u32, f64)> {
                    list.iter().map(|w| (w.feature, w.weight)).collect()
                };
                let by_row = |importance: std::collections::HashMap<u32, f64>| {
                    let renamed = importance.into_iter().map(|(c, w)| (rows[c as usize], w));
                    pairs(&top_weighted(renamed.collect(), TOP_WEIGHTED))
                };
                let forest = fit_rf(&cx, &y, &train, &cfg, draws).1;
                assert_eq!(pairs(&models.rf[f]), by_row(forest));
                let unit = tuning_unit(&split.scheme);
                let boosted = fit_gbt(&cx, &y, &train, tuning.as_deref(), unit, &cfg, seed, draws)
                    .unwrap()
                    .1;
                assert_eq!(pairs(&models.gbt[f]), by_row(boosted.importance()));
                assert_eq!(
                    models.gbt[f].is_empty(),
                    saved("gbt").fold_complexity[f] == 0.0
                );
                for (name, list) in [("rf", &models.rf[f]), ("gbt", &models.gbt[f])] {
                    assert!(list.len() <= TOP_WEIGHTED);
                    assert!(in_rank_order(list));
                    list.iter().for_each(|w| note(name, w.feature));
                }
            }
        }
        for name in LEARNERS {
            assert!(
                shifted.get(name).copied().unwrap_or(0) > 0,
                "{name}: {shifted:?}"
            );
        }
    }

    /// Synthetic: 90 samples, a third resistant, and twelve rows of noise
    /// present in about 40% of samples whatever the label, so the lasso has
    /// little to keep.
    fn noise_fixture() -> (BitMatrix, Bits, Vec<u32>) {
        use rand::{Rng, SeedableRng};
        let n = 90;
        let y = Bits::from_indices(n, (0..n).filter(|i| i % 3 == 0));
        let mut rng = rand_chacha::ChaCha8Rng::seed_from_u64(29);
        let mut x = BitMatrix::with_cols(n);
        for _ in 0..12 {
            let present: Vec<usize> = (0..n).filter(|_| rng.random_bool(0.4)).collect();
            x.push_row(Bits::from_indices(n, present).words());
        }
        (x, y, (0..n).map(|i| (i % 30) as u32).collect())
    }

    #[test]
    fn a_weightless_lasso_fold_scored_its_test_samples_at_its_training_prevalence() {
        // docs/report-rule-fields.md: a lasso fold model that is null or has
        // no weights scored its fold's test samples at the fold's training
        // prevalence, read from `test_fold` and the labels alone, and not at
        // the sigmoid of its intercept.
        let (x, y, groups) = noise_fixture();
        let n = x.cols;
        let cfg = EvalConfig {
            screened_features: 12,
            lasso_lambdas: 10,
            ..small_config()
        };
        let out = evaluate("fixture__fixturemycin", &x, &y, &groups, &[], 3, None, &cfg).unwrap();
        let (mut weightless, mut other_call) = (0, 0);
        for (split, models) in out.splits.iter().zip(&out.fold_models) {
            assert_eq!(models.scheme, split.scheme);
            let row = split
                .learners
                .iter()
                .find(|l| l.learner == "lasso")
                .unwrap();
            for (f, model) in models.lasso.iter().enumerate() {
                if model.as_ref().is_some_and(|m| !m.coef.is_empty()) {
                    continue;
                }
                weightless += 1;
                assert_eq!(row.fold_complexity[f], 0.0);
                let train: Vec<usize> = (0..n).filter(|&i| split.test_fold[i] != f).collect();
                let resistant = train.iter().filter(|&&i| y.contains(i)).count();
                let prevalence = resistant as f64 / train.len() as f64;
                assert_ne!(prevalence, 0.5);
                let literal = model.as_ref().map(|m| m.predict_proba(&x));
                for i in (0..n).filter(|&i| split.test_fold[i] == f) {
                    // By the rule: the saved score, exactly.
                    assert_eq!(row.scores[i], prevalence);
                    // Read as written: another score, and here another call.
                    if let Some(literal) = &literal {
                        assert_ne!(row.scores[i], literal[i]);
                        other_call += usize::from((literal[i] >= 0.5) != (prevalence >= 0.5));
                    }
                }
            }
        }
        assert!(weightless > 0, "the fixture has no weightless lasso fold");
        assert!(
            other_call > 0,
            "no weightless fold's intercept gives another call"
        );
    }

    /// Heaviest first, and among equal weights the lower feature first.
    fn in_rank_order(list: &[WeightedFeature]) -> bool {
        list.windows(2).all(|w| {
            let (a, b) = (&w[0], &w[1]);
            a.weight > b.weight || (a.weight == b.weight && a.feature < b.feature)
        })
    }

    #[test]
    fn ranked_lists_put_the_lower_feature_first_among_equal_weights() {
        let importance: std::collections::HashMap<u32, f64> =
            [(7, 0.5), (3, 0.5), (9, 1.0), (1, 0.25), (5, 0.5), (2, 0.25)]
                .into_iter()
                .collect();
        let ranked = |k: usize| -> Vec<(u32, f64)> {
            let list = top_weighted(importance.clone(), k);
            assert!(in_rank_order(&list));
            list.iter().map(|w| (w.feature, w.weight)).collect()
        };
        assert_eq!(
            ranked(6),
            [(9, 1.0), (3, 0.5), (5, 0.5), (7, 0.5), (1, 0.25), (2, 0.25)]
        );
        // A cut inside a tie keeps the lower features.
        assert_eq!(ranked(3), [(9, 1.0), (3, 0.5), (5, 0.5)]);
    }

    #[test]
    fn fold_forest_and_boosting_lists_keep_the_fifty_heaviest_features() {
        // Synthetic: 400 samples, half resistant, and 150 rows that each
        // follow the label with a quarter of the samples flipped, so the
        // forest and boosting spread their weight over many more than 50
        // features.
        use rand::{Rng, SeedableRng};
        let n = 400;
        let y = Bits::from_indices(n, (0..n).filter(|i| i % 2 == 0));
        let mut rng = rand_chacha::ChaCha8Rng::seed_from_u64(41);
        let mut x = BitMatrix::with_cols(n);
        for _ in 0..150 {
            let present: Vec<usize> = (0..n)
                .filter(|&i| y.contains(i) != rng.random_bool(0.25))
                .collect();
            x.push_row(Bits::from_indices(n, present).words());
        }
        let groups: Vec<u32> = (0..n).map(|i| (i / 2) as u32).collect();
        let cfg = EvalConfig {
            forest: ForestParams {
                n_trees: 40,
                max_depth: 4,
                ..Default::default()
            },
            gbt: GbtParams {
                rounds: 60,
                max_depth: 4,
                patience: 60,
                ..Default::default()
            },
            screened_features: 150,
            ..small_config()
        };
        let out = evaluate("fixture__fixturemycin", &x, &y, &groups, &[], 3, None, &cfg).unwrap();
        // The cap the docs give, and the final models' cap.
        assert_eq!(TOP_WEIGHTED, 50);
        assert_eq!(out.finals.rf_importance.len(), 50);
        assert_eq!(out.finals.gbt_importance.len(), 50);
        assert_eq!(out.fold_models.len(), 2);
        for models in &out.fold_models {
            assert_eq!((models.rf.len(), models.gbt.len()), (cfg.folds, cfg.folds));
            for f in 0..cfg.folds {
                let lengths = (models.rf[f].len(), models.gbt[f].len());
                assert_eq!(lengths, (50, 50), "{} fold {f}", models.scheme);
                assert!(in_rank_order(&models.rf[f]) && in_rank_order(&models.gbt[f]));
            }
        }
    }

    #[test]
    fn a_split_that_did_not_run_keeps_no_fold_models() {
        let n = 80;
        let y = Bits::from_indices(n, 0..40);
        let mut x = BitMatrix::with_cols(n);
        x.push_row(y.words());
        // One lineage holds every resistant sample, so the lineage split is
        // unavailable (`unsupported_lineage_plan_preserves_the_random_evaluation`).
        let groups: Vec<u32> = (0..n).map(|i| if i < 40 { 0 } else { i as u32 }).collect();
        let cfg = small_config();
        let out = evaluate("fixture__fixturemycin", &x, &y, &groups, &[], 3, None, &cfg).unwrap();
        assert_eq!(out.splits.len(), 2);
        assert!(out.splits[1].unavailable_reason.is_some());
        let schemes: Vec<&str> = out.fold_models.iter().map(|m| m.scheme.as_str()).collect();
        assert_eq!(schemes, ["random"]);
        // An output written before the fold models were kept reads back
        // with none.
        let mut older = serde_json::to_value(&out).unwrap();
        older.as_object_mut().unwrap().remove("fold_models");
        let older: EvalOutput = serde_json::from_value(older).unwrap();
        assert!(older.fold_models.is_empty());
    }

    #[test]
    fn fold_models_are_written_in_the_final_models_shapes_under_dataset_rows() {
        // Candidate c of this fold is dataset row `rows[c]`.
        let rows = [4, 7, 9, 12];
        let mut split = FoldSplit::new("random");
        let scm = FoldScm {
            kind: Kind::Disjunction,
            p: 0.5,
            rules: vec![
                Rule {
                    feature: 2,
                    presence: true,
                },
                Rule {
                    feature: 0,
                    presence: false,
                },
            ],
        };
        split.push(FoldModel::Scm(scm), &rows);
        let leaf = |value: f64, n_resistant: u32| tree::Node {
            value,
            n: 5,
            n_resistant,
            depth: 1,
            split: None,
        };
        let root = tree::Node {
            value: 0.5,
            n: 10,
            n_resistant: 5,
            depth: 0,
            split: Some(tree::Split {
                feature: 1,
                gain: 0.25,
                absent: 1,
                present: 2,
            }),
        };
        let cart = tree::Tree {
            nodes: vec![root, leaf(0.0, 0), leaf(1.0, 5)],
        };
        split.push(FoldModel::Cart(cart), &rows);
        let weighted = |feature: u32, weight: f64| WeightedFeature { feature, weight };
        split.push(
            FoldModel::Rf(vec![weighted(3, 0.75), weighted(0, 0.25)]),
            &rows,
        );
        split.push(FoldModel::Gbt(vec![]), &rows);
        let lasso = lasso::LassoModel {
            intercept: -1.0,
            coef: vec![(2, 1.5)],
            lambda: 0.125,
        };
        split.push(FoldModel::Lasso(Some(lasso)), &rows);
        split.push(FoldModel::Lasso(None), &rows);
        assert_eq!(
            serde_json::to_value(&split).unwrap(),
            serde_json::json!({
                "scheme": "random",
                "scm": [{
                    "kind": "disjunction",
                    "p": 0.5,
                    "rules": [
                        {"feature": 9, "presence": true},
                        {"feature": 4, "presence": false},
                    ],
                }],
                "cart": [{"nodes": [
                    {
                        "value": 0.5, "n": 10, "n_resistant": 5, "depth": 0,
                        "split": {"feature": 7, "gain": 0.25, "absent": 1, "present": 2},
                    },
                    {"value": 0.0, "n": 5, "n_resistant": 0, "depth": 1, "split": null},
                    {"value": 1.0, "n": 5, "n_resistant": 5, "depth": 1, "split": null},
                ]}],
                "rf": [[{"feature": 12, "weight": 0.75}, {"feature": 4, "weight": 0.25}]],
                "gbt": [[]],
                "lasso": [{"intercept": -1.0, "coef": [[9, 1.5]], "lambda": 0.125}, null],
            })
        );
        // Row 12 is in the forest's ranking alone, so it gets no pattern.
        let patterns: Vec<u32> = split.pattern_features().into_iter().collect();
        assert_eq!(patterns, [4, 7, 9]);
    }

    #[test]
    fn folds_name_the_lowest_row_among_patterns_equal_on_their_training_samples() {
        // Row 1 is the label. Row 0 is the label without sample 0, so the
        // two differ on that one sample.
        let n = 80;
        let y = Bits::from_indices(n, (0..n).filter(|i| i % 2 == 0));
        let mut x = BitMatrix::with_cols(n);
        x.push_row(Bits::from_indices(n, (2..n).filter(|i| i % 2 == 0)).words());
        x.push_row(y.words());
        let groups: Vec<u32> = (0..n).map(|i| (i / 4) as u32).collect();
        let cfg = small_config();
        let out = evaluate("fixture__fixturemycin", &x, &y, &groups, &[], 3, None, &cfg).unwrap();
        let (split, models) = (&out.splits[0], &out.fold_models[0]);
        assert_eq!(
            (split.scheme.as_str(), models.scheme.as_str()),
            ("random", "random")
        );
        let named = |fold: usize| -> Vec<u32> {
            models.scm[fold].rules.iter().map(|r| r.feature).collect()
        };
        // The fold that holds sample 0 out cannot tell the rows apart and
        // names the lower. The other fold trains on sample 0 and names the
        // label itself. The two ids differ for one pattern on the first
        // fold's training samples: equal ids undercount agreement between
        // folds, and the presence patterns settle it.
        let without = split.test_fold[0];
        assert_eq!(named(without), [0]);
        assert_eq!(named(1 - without), [1]);
        let (train, _) = fold_masks(&split.test_fold, cfg.folds).swap_remove(without);
        assert!(train.iter_ones().all(|i| x.get(0, i) == x.get(1, i)));
    }

    /// Synthetic: 120 genomes in 30 lineages of four, every third lineage
    /// resistant, and rows of noise: eight drawn at random, a copy of row 0
    /// with the four genomes of lineage 12 flipped, and a row in three
    /// genomes of three lineages. The matrix is filtered and collapsed at
    /// `min_support` 3, as `evaluate` requires. A fold that tests one of the
    /// three genomes drops the last row from its candidates, and the lineage
    /// fold that tests lineage 12 keeps one copy only, so a fold's
    /// candidates are not the matrix. The noise seed is one under which both
    /// boosting and the lasso have null folds on both splits.
    fn null_fold_fixture() -> (BitMatrix, Bits, Vec<u32>) {
        use rand::{Rng, SeedableRng};
        let n = 120;
        let groups: Vec<u32> = (0..n as u32).map(|i| i / 4).collect();
        let y = Bits::from_indices(n, (0..n).filter(|&i| groups[i].is_multiple_of(3)));
        let mut rng = rand_chacha::ChaCha8Rng::seed_from_u64(20261117);
        let mut x = BitMatrix::with_cols(n);
        for _ in 0..8 {
            let present: Vec<usize> = (0..n).filter(|_| rng.random_bool(0.3)).collect();
            x.push_row(Bits::from_indices(n, present).words());
        }
        let copy = (0..n).filter(|&i| x.get(0, i) != (groups[i] == 12));
        x.push_row(Bits::from_indices(n, copy.collect::<Vec<_>>()).words());
        x.push_row(Bits::from_indices(n, [5, 47, 77]).words());
        assert_eq!(training_candidates(&x, &Bits::ones(n), 3).1.len(), x.rows);
        (x, y, groups)
    }

    fn null_fold_config() -> EvalConfig {
        EvalConfig {
            folds: 3,
            inner_folds: 2,
            scm_max_rules: 1,
            cart_max_depth: 1,
            forest: ForestParams {
                n_trees: 2,
                max_depth: 2,
                ..Default::default()
            },
            gbt: GbtParams {
                rounds: 30,
                max_depth: 2,
                patience: 3,
                ..Default::default()
            },
            screened_features: 6,
            lasso_lambdas: 8,
            ..Default::default()
        }
    }

    /// The null-fold experiment beside a supplied partition (the merge of
    /// `main` after #119). `evaluate_assigned` fits through
    /// `fit_outer_folds`, as `evaluate` does, but a record names every sample
    /// outside a fold's training ones as tested, which an unused sample would
    /// make wrong. So with the experiment on it is refused before any fit
    /// and no record is written; with it off the same call scores the
    /// partition, on which a boosting or lasso fold is null.
    #[test]
    fn a_supplied_partition_is_refused_while_the_null_fold_refit_is_on() {
        let (x, y, groups) = null_fold_fixture();
        let cfg = null_fold_config();
        let [(scheme, random), _] = outer_splits(&y, &groups, &[], &cfg);
        assert_eq!(scheme, "random");
        let roles: Vec<Role> = random.iter().map(|&f| Role::Test(f)).collect();
        let run = || {
            evaluate_assigned(
                "fixture__fixturemycin",
                "supplied",
                &x,
                &y,
                &groups,
                &[],
                &roles,
                cfg.folds,
                &cfg,
            )
        };
        let scored = run().unwrap();
        assert!(
            scored
                .result
                .learners
                .iter()
                .filter(|l| l.learner == "gbt" || l.learner == "lasso")
                .any(|l| l.fold_complexity.contains(&0.0))
        );
        let (refused, lines) = null_refit::capture(run);
        let error = refused.err().expect("refused while the experiment is on");
        assert_eq!(
            format!("{error:#}"),
            "the null-fold refit records the folds grm learn draws; it does not run on a supplied partition"
        );
        assert!(lines.is_empty(), "{lines:?}");
    }

    /// The call site of the null-fold experiment (S2 of the review of 15 Nov
    /// 2026). Each record `evaluate` writes is the arms `null_refit` fits
    /// when called directly with that outer fold's training candidates, its
    /// seed (`fold_seed`), its draws seed (`fold_draw_seed`) and its split's
    /// tuning groups; its counts are the fold's; and there is one record per
    /// null boosting or lasso fold. The test also shows that each of those
    /// inputs moves some record, so a wrong one at the call site fails it.
    #[test]
    fn each_null_fold_record_is_the_arms_of_that_folds_candidates_seeds_and_groups() {
        let (x, y, groups) = null_fold_fixture();
        let min_support = 3;
        let mut moved = [false; 4];
        for learner_seed in [None, Some(4242)] {
            let cfg = EvalConfig {
                learner_seed,
                ..null_fold_config()
            };
            let run = || {
                let dataset = "fixture__fixturemycin";
                evaluate(dataset, &x, &y, &groups, &[], min_support, None, &cfg).unwrap()
            };
            let (out, lines) = null_refit::capture(run);
            let records: Vec<serde_json::Value> = lines
                .iter()
                .map(|l| serde_json::from_str(l).unwrap())
                .collect();
            let mut null_folds = vec![];
            for split in &out.splits {
                for l in split
                    .learners
                    .iter()
                    .filter(|l| l.learner == "gbt" || l.learner == "lasso")
                {
                    for (fold, c) in l.fold_complexity.iter().enumerate() {
                        if *c == 0.0 {
                            null_folds.push((split.scheme.clone(), fold, l.learner.clone()));
                        }
                    }
                }
            }
            let mut keys: Vec<(String, usize, String)> = records
                .iter()
                .map(|r| {
                    let s = |k: &str| r[k].as_str().unwrap().to_string();
                    (
                        s("scheme"),
                        r["fold"].as_u64().unwrap() as usize,
                        s("learner"),
                    )
                })
                .collect();
            null_folds.sort();
            keys.sort();
            assert_eq!(keys, null_folds, "learner seed {learner_seed:?}");
            for (scheme, folds) in outer_splits(&y, &groups, &[], &cfg) {
                let tuning = tuning_groups(scheme, &groups, &[]);
                let unit = tuning_unit(scheme);
                for (fi, (train, _)) in fold_masks(&folds, cfg.folds).into_iter().enumerate() {
                    let (cx, _) = training_candidates(&x, &train, min_support);
                    let test: Vec<usize> = (0..x.cols).filter(|&i| !train.contains(i)).collect();
                    let (seed, draws) = (fold_seed(&cfg, fi), fold_draw_seed(&cfg, fi));
                    let context = test_context(scheme, fi);
                    let lasso = |x: &BitMatrix, groups: Option<&[u32]>, seed| {
                        let arms = null_refit::lasso_arms(
                            x, &y, &train, groups, &cfg, seed, context, &test,
                        );
                        value(&arms.unwrap())
                    };
                    let gbt = |x: &BitMatrix, groups: Option<&[u32]>, seed, draws| {
                        let arms = null_refit::gbt_arms(
                            x, &y, &train, groups, unit, &cfg, seed, draws, &test,
                        );
                        value(&arms.unwrap())
                    };
                    let here = records
                        .iter()
                        .filter(|r| r["scheme"] == scheme && r["fold"] == fi);
                    for r in here {
                        let at = format!(
                            "{scheme} fold {fi} {} learner seed {learner_seed:?}",
                            r["learner"]
                        );
                        assert_eq!(r["dataset"], "fixture__fixturemycin", "{at}");
                        assert_eq!(r["train_n"], train.count(), "{at}");
                        assert_eq!(r["train_resistant"], train.and(&y).count(), "{at}");
                        assert_eq!(r["test"], value(&test), "{at}");
                        let other = fold_seed(&cfg, fi + 1);
                        if r["learner"] == "lasso" {
                            assert!(r["gbt"].is_null(), "{at}");
                            let want = lasso(&cx, tuning.as_deref(), seed);
                            assert_eq!(r["lasso"], want, "{at}");
                            moved[0] |= lasso(&cx, tuning.as_deref(), other) != want;
                            moved[1] |= lasso(&cx, None, seed) != want;
                            moved[2] |= lasso(&x, tuning.as_deref(), seed) != want;
                        } else {
                            assert_eq!(r["learner"], "gbt", "{at}");
                            assert!(r["lasso"].is_null(), "{at}");
                            let want = gbt(&cx, tuning.as_deref(), seed, draws);
                            assert_eq!(r["gbt"], want, "{at}");
                            moved[0] |= gbt(&cx, tuning.as_deref(), other, draws) != want;
                            moved[1] |= gbt(&cx, None, seed, draws) != want;
                            moved[2] |= gbt(&x, tuning.as_deref(), seed, draws) != want;
                            moved[3] |= gbt(&cx, tuning.as_deref(), seed, seed) != want;
                        }
                    }
                }
            }
        }
        // The seed, the tuning groups, the candidates and the draws seed
        // each move some record.
        assert_eq!(moved, [true; 4]);
    }

    /// SYNTHETIC: 60 samples, 30 rows, the first 8 a joined arm's first block.
    fn screen_case(seed: u64) -> (BitMatrix, Bits, Bits) {
        let n = 60;
        let mut s = seed;
        let mut next = move || {
            s ^= s << 13;
            s ^= s >> 7;
            s ^= s << 17;
            s
        };
        let mut x = BitMatrix::with_cols(n);
        for _ in 0..30 {
            let on: Vec<usize> = (0..n).filter(|_| next() % 3 == 0).collect();
            x.push_row(Bits::from_indices(n, on).words());
        }
        let y = Bits::from_indices(n, (0..n).filter(|i| i % 5 < 2));
        let train = Bits::from_indices(n, (0..n).filter(|i| i % 6 != 0));
        (x, y, train)
    }

    #[test]
    fn a_screen_count_reads_each_learners_own_fitting_genomes_only() {
        let (x, y, train) = screen_case(7);
        let cfg = EvalConfig {
            screened_features: 5,
            ..small_config()
        };
        let first = 8;
        let (cx, rows) = training_candidates(&x, &train, 3);
        assert!(cx.rows > 5, "the case screens");
        let got = screen_fold(1, &cx, &rows, first, &y, &train, None, "lineage", &cfg).unwrap();
        let in_first = |f: &u32| (rows[*f as usize] as usize) < first;
        let seed = fold_seed(&cfg, 1);
        let unit = tuning_unit("lineage");
        let fits = [
            train.clone(),
            tuning_holdout(&train, &y, None, unit, GBT_HOLDOUT, seed)
                .unwrap()
                .0,
            tuning_holdout(&train, &y, None, unit, LASSO_HOLDOUT, seed)
                .unwrap()
                .0,
        ];
        assert_eq!((got.fold, got.train), (1, train.count() as usize));
        assert_eq!(
            got.candidates_first,
            (0..cx.rows as u32).filter(|f| in_first(f)).count()
        );
        assert_eq!(got.candidates_first + got.candidates_k31, cx.rows);
        assert_eq!(got.learners.len(), 3);
        for ((l, fit), name) in got.learners.iter().zip(&fits).zip(["rf", "gbt", "lasso"]) {
            let feats = screen(&cx, &y, fit, 5);
            assert_eq!(l.learner, name);
            assert_eq!(l.fitting, fit.count() as usize, "{name}");
            assert_eq!(l.screened, feats.len(), "{name}");
            assert_eq!(
                l.passed_first,
                feats.iter().filter(|f| in_first(f)).count(),
                "{name}"
            );
            assert_eq!(l.passed_first + l.passed_k31, l.screened, "{name}");
        }
        // No outer test genome moves it: their labels and values flipped.
        let test = train.not();
        let flipped_y = y.and(&train).or(&test.and_not(&y));
        let mut flipped_x = BitMatrix::with_cols(x.cols);
        for r in 0..x.rows {
            let row = Bits::from_words(x.cols, x.row(r).to_vec());
            flipped_x.push_row(row.and(&train).or(&test.and_not(&row)).words());
        }
        let (fx, frows) = training_candidates(&flipped_x, &train, 3);
        assert_eq!(frows, rows);
        let again = screen_fold(
            1, &fx, &frows, first, &flipped_y, &train, None, "lineage", &cfg,
        )
        .unwrap();
        assert_eq!(again, got);
    }

    /// A SYNTHETIC cohort for the C4 arms: 96 samples in 12 lineages whose
    /// resistance rates differ, 30 rows, some lineage markers and some
    /// that track resistance inside lineages.
    fn c4_cohort() -> (BitMatrix, Bits, Vec<u32>) {
        let n = 96;
        let mut s = 0x9e37_79b9_7f4a_7c15u64;
        let mut next = move || {
            s ^= s << 13;
            s ^= s >> 7;
            s ^= s << 17;
            (s >> 11) as f64 / (1u64 << 53) as f64
        };
        let groups: Vec<u32> = (0..n).map(|i| (i / 8) as u32).collect();
        let y = Bits::from_indices(
            n,
            (0..n).filter(|&i| next() < 0.2 + 0.6 * (groups[i] % 4) as f64 / 3.0),
        );
        let mut x = BitMatrix::with_cols(n);
        for r in 0..30 {
            let row = Bits::from_indices(
                n,
                (0..n).filter(|&i| {
                    let p = match r % 3 {
                        0 => 0.15 + 0.7 * (groups[i] % 4) as f64 / 3.0,
                        1 => {
                            if y.contains(i) {
                                0.7
                            } else {
                                0.3
                            }
                        }
                        _ => 0.5,
                    };
                    next() < p
                }),
            );
            x.push_row(row.words());
        }
        (x, y, groups)
    }

    fn c4_inputs(dir: &std::path::Path, tag: &str) -> serde_json::Value {
        let path = dir.join(screen_arm::file_name("fixture__fixturemycin", tag, 42));
        serde_json::from_slice(&std::fs::read(&path).unwrap()).unwrap()
    }

    /// With no arm, `screen_fit` is the release's screen and no inputs file
    /// is written; the floor arm at a floor the dataset already has gives
    /// the same bytes; the adjusted screen ranks every forest, boosting and
    /// lasso fit, outer and final, on that fit's own fitting samples, with
    /// strata from the lineages and not from the labels.
    #[test]
    fn the_c4_arms_leave_unflagged_fits_alone_and_screen_every_fit_by_their_rule() {
        use screen_arm::{Rule, Settings};
        let (x, y, groups) = c4_cohort();
        let mut cfg = small_config();
        cfg.screened_features = 6;
        let plain = evaluate("fixture__fixturemycin", &x, &y, &groups, &[], 3, None, &cfg).unwrap();
        let dir = tempfile::tempdir().unwrap();
        let settings = |rule, floor| Settings {
            rule,
            floor,
            out: dir.path().to_path_buf(),
        };
        let floor = screen_arm::with_settings(settings(Rule::ChiSquare, Some(3)), || {
            evaluate("fixture__fixturemycin", &x, &y, &groups, &[], 3, None, &cfg).unwrap()
        });
        for (a, b) in plain.splits.iter().zip(&floor.splits) {
            assert_eq!(timeless(a), timeless(b));
        }
        assert_eq!(
            serde_json::to_vec(&plain.finals).unwrap(),
            serde_json::to_vec(&floor.finals).unwrap()
        );
        let floor_file = c4_inputs(dir.path(), "drawn");
        assert_eq!(
            (
                &floor_file["arm"],
                &floor_file["rule"],
                &floor_file["budget"]
            ),
            (
                &serde_json::json!("floor"),
                &serde_json::json!("chi-square"),
                &serde_json::json!(6)
            )
        );
        std::fs::remove_file(dir.path().join(screen_arm::file_name(
            "fixture__fixturemycin",
            "drawn",
            42,
        )))
        .unwrap();

        screen_arm::with_settings(settings(Rule::Cmh, None), || {
            evaluate("fixture__fixturemycin", &x, &y, &groups, &[], 3, None, &cfg).unwrap()
        });
        let file = c4_inputs(dir.path(), "drawn");
        assert_eq!(
            (&file["arm"], &file["rule"]),
            (&serde_json::json!("cmh"), &serde_json::json!("cmh"))
        );
        assert_eq!(file["dataset_rows"], 30);
        let folds = file["folds"].as_array().unwrap();
        let screens = file["screens"].as_array().unwrap();
        let schemes: Vec<&str> = plain
            .splits
            .iter()
            .filter(|s| s.unavailable_reason.is_none())
            .map(|s| s.scheme.as_str())
            .collect();
        assert_eq!(schemes, ["random", "lineage"]);
        assert_eq!(folds.len(), schemes.len() * cfg.folds);
        // Each outer fit of rf, gbt and lasso is screened once, and the final
        // models three times: lasso, the forest and boosting.
        assert_eq!(screens.len(), folds.len() * 3 + 3);
        assert_eq!(
            file["fit_times"].as_array().unwrap().len(),
            folds.len() * LEARNERS.len()
        );
        for note in screens {
            assert_eq!(note["rule"], "cmh", "{note}");
            assert!(note["kept"].as_u64().unwrap() <= 6);
            if note["split"] == "final" {
                // The final forest fits every sample; boosting and the lasso
                // hold out their tuning samples.
                let fitting = note["fitting"].as_u64().unwrap();
                assert!(
                    if note["learner"] == "rf" {
                        fitting == 96
                    } else {
                        fitting < 96
                    },
                    "{note}"
                );
                continue;
            }
            let fold = folds
                .iter()
                .find(|f| f["split"] == note["split"] && f["fold"] == note["fold"])
                .unwrap();
            assert_eq!(note["candidates"], fold["candidates"]);
            let (fitting, train) = (
                note["fitting"].as_u64().unwrap(),
                fold["train"].as_u64().unwrap(),
            );
            // The forest screens the outer training samples; boosting and the
            // lasso those less their tuning holdout. No test sample is read.
            if note["learner"] == "rf" {
                assert_eq!(fitting, train);
            } else {
                assert!(fitting < train, "{note}");
            }
            // The strata are lineages: twelve of them, each of 8 samples,
            // so a fit of `fitting` samples meets at least fitting / 8.
            let strata = note["strata"].as_u64().unwrap();
            assert!(strata >= fitting.div_ceil(8) && strata <= 12, "{note}");
            assert!(note["informative_strata"].as_u64().unwrap() <= strata);
        }
        for learner in ["rf", "gbt", "lasso"] {
            for split in &schemes {
                for fold in 0..cfg.folds {
                    let n = screens
                        .iter()
                        .filter(|s| {
                            s["learner"] == learner && s["split"] == *split && s["fold"] == fold
                        })
                        .count();
                    assert_eq!(n, 1, "{learner} {split} fold {fold}");
                }
            }
        }
    }

    /// A supplied partition (a fold file, as the study split runs) under the
    /// adjusted arm: each outer fit of rf, gbt and lasso takes the arm's
    /// screen on its own fitting samples, and no tested or unused sample's
    /// lineage reaches a stratum. SYNTHETIC: `c4_cohort`'s 12 lineages of 8;
    /// lineage 11 is unused and lineage g tests in fold g % 2, so fold 0
    /// trains on 5 lineages (40 samples) and fold 1 on 6 (48).
    #[test]
    fn a_supplied_partition_under_the_adjusted_arm_screens_each_fit_on_its_training_lineages() {
        use screen_arm::{Rule, Settings};
        let (x, y, groups) = c4_cohort();
        let mut cfg = small_config();
        cfg.screened_features = 6;
        let roles: Vec<Role> = groups
            .iter()
            .map(|&g| {
                if g == 11 {
                    Role::Unused
                } else {
                    Role::Test(g as usize % cfg.folds)
                }
            })
            .collect();
        let dir = tempfile::tempdir().unwrap();
        let settings = Settings {
            rule: Rule::Cmh,
            floor: None,
            out: dir.path().to_path_buf(),
        };
        screen_arm::with_settings(settings, || {
            evaluate_assigned(
                "fixture__fixturemycin",
                "study-whole",
                &x,
                &y,
                &groups,
                &[],
                &roles,
                cfg.folds,
                &cfg,
            )
            .unwrap()
        });
        let file = c4_inputs(dir.path(), "study-whole");
        assert_eq!(
            (&file["rule"], &file["samples"], &file["min_support"]),
            (
                &serde_json::json!("cmh"),
                &serde_json::json!(88),
                &serde_json::json!(3)
            )
        );
        let folds = file["folds"].as_array().unwrap();
        let train: Vec<u64> = folds.iter().map(|f| f["train"].as_u64().unwrap()).collect();
        assert_eq!(train, [40, 48]);
        let screens = file["screens"].as_array().unwrap();
        assert_eq!(screens.len(), cfg.folds * 3);
        for note in screens {
            assert_eq!(
                (&note["split"], &note["rule"]),
                (&serde_json::json!("study-whole"), &serde_json::json!("cmh")),
                "{note}"
            );
            let fold = note["fold"].as_u64().unwrap() as usize;
            let fitting = note["fitting"].as_u64().unwrap();
            let strata = note["strata"].as_u64().unwrap();
            if note["learner"] == "rf" {
                assert_eq!(fitting, train[fold], "{note}");
            } else {
                assert!(fitting < train[fold], "{note}");
            }
            // Fold 0 trains on 5 lineages, fold 1 on 6: a test or unused
            // sample's lineage in a stratum would push the count past them.
            assert!(
                strata >= fitting.div_ceil(8) && strata <= [5, 6][fold],
                "{note}"
            );
        }
        for learner in ["rf", "gbt", "lasso"] {
            for fold in 0..cfg.folds {
                let n = screens
                    .iter()
                    .filter(|s| s["learner"] == learner && s["fold"] == fold)
                    .count();
                assert_eq!(n, 1, "{learner} fold {fold}");
            }
        }
    }

    /// The floor arm moves the floor of a supplied partition, and the
    /// projection's, to its own; without it the floor is max(3, n/200).
    #[test]
    fn the_floor_arm_fixes_the_support_floor() {
        use screen_arm::{Rule, Settings};
        assert_eq!(support_floor(1600), 8);
        let dir = tempfile::tempdir().unwrap();
        let settings = Settings {
            rule: Rule::ChiSquare,
            floor: Some(3),
            out: dir.path().to_path_buf(),
        };
        let (fixed, supplied) = screen_arm::with_settings(settings, || {
            (
                support_floor(1600),
                supplied_support_floor(&[Role::Train; 1600]),
            )
        });
        assert_eq!((fixed, supplied), (3, 3));
        assert_eq!(support_floor(1600), 8);
        assert_eq!(support_floor(400), 3);
    }
}
