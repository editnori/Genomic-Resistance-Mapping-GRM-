//! No-genome baselines beside the scores a dataset report already holds.
//!
//! `grm baselines` reads each dataset report that `index.json` lists and
//! writes `baselines/<id>.json` beside it, pinned to the report's SHA-256 and
//! evaluation revision. Each file says what two predictors that never read a
//! genome reach on the report's own folds and isolates:
//!
//! - **MIC constants.** For a report with an MIC model, the MIC regression's
//!   folds are rebuilt as `evaluate_mic` draws them: `training_folds` over the
//!   isolates with an MIC, stratified by reported label, with the report's
//!   fold count and seed. These are random folds with no lineage split, and
//!   they are not the random split's folds. Two checks compare the rebuilt
//!   folds with what the fit left in the report, and a contradiction in
//!   either refuses that dataset alone (`Refused`):
//!   - Prediction check (`prediction_check`): on the report's own folds,
//!     groups of isolates with bit-identical held-out predictions lie inside
//!     one fold, and a wrong rebuild splits them. Correct folds can also give
//!     one value to isolates of two folds (a shared start value, or
//!     identical models), so a group that spans folds is evidence, not
//!     proof, of a wrong rebuild. A fold where every isolate gets the same
//!     prediction (a constant fold, as when early stopping keeps no tree and
//!     the model predicts its start value, the upper median of its fit
//!     part's exact MICs) is left out. The check counts the groups of two or
//!     more inside one fold that is not constant (I, `groups`) and the
//!     groups that span folds, at least one of them not constant (X,
//!     `spanning`). X = 0 and I >= 1: the folds are `confirmed`. X > I: the
//!     spanning groups outnumber the others, and the dataset is refused.
//!     Otherwise (0 < X <= I, or no group at all) the file is written
//!     unconfirmed. The rule has no cut read off a release: it asks only
//!     whether the spanning groups outnumber the others.
//!   - Threshold check (`threshold_check`): each rebuilt fold's training
//!     isolates give the
//!     threshold `infer_breakpoint` fits, and that must equal the
//!     threshold the report stores for the fold (`mic.fold_breakpoints`,
//!     every field, agreement bit for bit). It confirms nothing on its own:
//!     where every stored threshold has agreement 1, any stratified draw with
//!     the same fold sizes matches it.
//!
//!   The file also counts the exact-MIC test isolates that have a
//!   k-mer-identical copy (`kset_duplicates`) in their training fold: the
//!   MIC folds are random and can split such a group. Each fold's exact
//!   training MICs give two constants, their lower median and their mode
//!   (`lower_median`, `mode`), and each is scored on the same isolates as the
//!   model: the fold's test isolates with an exact MIC, agreeing when the
//!   constant is within one doubling dilution of the measured value. The
//!   model's own count is recounted from `mic.predicted` under the same rule
//!   (`evaluate_mic` rounds a prediction half away from zero) and must give
//!   the report's `mic.essential_agreement` exactly, or the dataset is
//!   refused.
//! - **Lineage identity.** For each split with folds, inside each of the
//!   report's own outer folds a test genome scores the resistant share of its
//!   lineage among that fold's training genomes, or the fold's training
//!   prevalence (`learn::train_prevalence`) when no training genome shares
//!   its lineage (`lineage_scores`). The report's own scoring conventions
//!   apply, as `eval::pool_folds` scores a learner row: every genome of the
//!   split is scored by the fold that held it out and the scores are pooled;
//!   `learn::binary_metrics` calls a score of 0.5 or more resistant, so a
//!   share of exactly 0.5 is a resistant call, as a null fit's prevalence
//!   is; balanced accuracy is the mean of sensitivity (recall) and
//!   specificity, and a class with no genome counts 0. On a lineage split no
//!   test genome's lineage is in training, so every genome takes its fold's
//!   prevalence: the baseline has nothing of the lineage to read there.
//!
//! The command fits nothing, reads no network and changes no report file,
//! and a second run writes the same bytes. A refused dataset gets no file,
//! and a file it had is removed; every other dataset is written, and the
//! command then fails, naming each refused dataset with its I and X. A
//! report it cannot read, a report of another evaluation revision than this
//! binary's, or an index id that is not a plain file name stops the whole
//! run before anything is written. It also writes `baselines.json` beside
//! `index.json` (`summary`): for the pages that show many datasets, each
//! written dataset's model figure, constant and confirmation, pinned to the
//! SHA-256 of the `index.json` bytes it read and to each report's SHA-256.
//!
//! A baseline says what a predictor reaches without the genome on these
//! cross-validated folds of public data. It does not say why a model does
//! better or worse, and it is not a clinical figure.

use crate::analysis::EVALUATION_REVISION;
use crate::bits::Bits;
use crate::curate::{Breakpoint, Call, Censor, Mic, infer_breakpoint};
use crate::learn::eval::MicResult;
use crate::learn::{binary_metrics, fold_masks, train_prevalence, training_folds};
use anyhow::{Context, Result, bail, ensure};
use serde::Deserialize;
use serde_json::{Value, json};
use std::collections::{BTreeMap, BTreeSet};
use std::path::Path;

/// Version of the baselines file layout.
pub const SCHEMA: u32 = 1;

/// The parts of a dataset report this command reads; serde skips the rest.
#[derive(Deserialize)]
struct Report {
    id: String,
    #[serde(default)]
    evaluation: Option<Evaluation>,
    #[serde(default)]
    eval_config: Option<FoldConfig>,
    samples: Vec<Sample>,
    #[serde(default)]
    splits: Vec<Split>,
    #[serde(default)]
    mic: Option<MicResult>,
    /// Groups of sample indices with the same k-mer set; absent when none.
    #[serde(default)]
    kset_duplicates: Vec<Vec<usize>>,
}

#[derive(Deserialize)]
struct Evaluation {
    #[serde(default)]
    revision: Option<u64>,
}

#[derive(Deserialize)]
struct FoldConfig {
    folds: usize,
    seed: u64,
}

#[derive(Deserialize)]
struct Sample {
    resistant: bool,
    #[serde(default)]
    mic: Option<Mic>,
    lineage: u32,
}

#[derive(Deserialize)]
struct Split {
    scheme: String,
    folds: usize,
    /// Empty when the split is unavailable.
    #[serde(default)]
    test_fold: Vec<usize>,
}

/// The lower median: with m values sorted, the ceil(m/2)-th smallest, so an
/// even count takes the lower of its two middle values. None when empty.
pub fn lower_median(values: &[i32]) -> Option<i32> {
    let mut sorted = values.to_vec();
    sorted.sort_unstable();
    sorted.get(sorted.len().checked_sub(1)? / 2).copied()
}

/// The value with the highest count, a tie going to the lowest, and whether another
/// value was as frequent. None when empty.
pub fn mode(values: &[i32]) -> Option<(i32, bool)> {
    let mut count: BTreeMap<i32, usize> = BTreeMap::new();
    for &v in values {
        *count.entry(v).or_default() += 1;
    }
    let top = count.values().copied().max()?;
    // BTreeMap iterates in ascending order, so the first at the top count is
    // the lowest.
    let (&value, _) = count.iter().find(|&(_, &c)| c == top)?;
    let tied = count.values().filter(|&&c| c == top).count() > 1;
    Some((value, tied))
}

/// Within one doubling dilution, the report's rule for essential agreement.
fn within_one(a: i32, b: i32) -> bool {
    (a - b).abs() <= 1
}

/// One MIC fold: its exact training and test MICs, the two constants its
/// exact training MICs give, and how many of its exact test MICs the model
/// and each constant reach within one dilution.
#[derive(Clone, Debug, PartialEq)]
pub struct FoldConstants {
    pub train_exact: usize,
    pub test_exact: usize,
    /// None when the fold's training isolates have no exact MIC.
    pub median: Option<i32>,
    /// The two middle values differ, so the upper median is another value.
    pub upper_median_differs: bool,
    pub mode: Option<i32>,
    pub mode_tied: bool,
    pub agree_model: usize,
    pub agree_median: usize,
    pub agree_mode: usize,
}

/// Scores the training-fold constants and the model on the exact MICs of
/// each fold's test isolates. `folds` gives each isolate's MIC fold, and
/// `usize::MAX` for an isolate without an MIC, as `training_folds` leaves it.
/// Each constant comes from the exact MICs of the fold's training isolates
/// only, never from its test isolates.
pub fn mic_constants(
    mics: &[Option<Mic>],
    predicted: &[f64],
    folds: &[usize],
    k: usize,
) -> Vec<FoldConstants> {
    let exact = |i: usize| {
        mics[i]
            .filter(|m| m.censor == Censor::Exact)
            .map(|m| m.log2)
    };
    fold_masks(folds, k)
        .into_iter()
        .map(|(train, test)| {
            let train_values: Vec<i32> = train.iter_ones().filter_map(exact).collect();
            let test_isolates: Vec<(usize, i32)> = test
                .iter_ones()
                .filter_map(|i| exact(i).map(|v| (i, v)))
                .collect();
            let median = lower_median(&train_values);
            let mut sorted = train_values.clone();
            sorted.sort_unstable();
            let m = sorted.len();
            let upper_median_differs =
                m >= 2 && m.is_multiple_of(2) && sorted[m / 2 - 1] != sorted[m / 2];
            let (mode, mode_tied) = match mode(&train_values) {
                Some((v, tied)) => (Some(v), tied),
                None => (None, false),
            };
            let reach = |c: Option<i32>| {
                c.map_or(0, |c| {
                    test_isolates
                        .iter()
                        .filter(|&&(_, v)| within_one(c, v))
                        .count()
                })
            };
            FoldConstants {
                train_exact: m,
                test_exact: test_isolates.len(),
                median,
                upper_median_differs,
                mode,
                mode_tied,
                agree_model: test_isolates
                    .iter()
                    .filter(|&&(i, v)| within_one(predicted[i].round() as i32, v))
                    .count(),
                agree_median: reach(median),
                agree_mode: reach(mode),
            }
        })
        .collect()
}

/// What the held-out MIC predictions say about the rebuilt folds.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PredictionCheck {
    /// I: groups of two or more isolates with bit-identical predictions,
    /// each inside one rebuilt fold that is not a constant fold.
    pub groups: usize,
    /// The isolates in those groups.
    pub isolates: usize,
    /// X: groups of two or more isolates with bit-identical predictions that
    /// span rebuilt folds, at least one of them not a constant fold.
    pub spanning: usize,
    /// Rebuilt folds where every isolate gets the same prediction.
    pub constant_folds: Vec<usize>,
}

/// What the prediction check makes of a dataset's rebuilt folds.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Verdict {
    /// X = 0 and I >= 1: the file is written with `confirmed: true`.
    Confirmed,
    /// 0 < X <= I, or no group at all: written with `confirmed: false`.
    Unconfirmed,
    /// X > I: the spanning groups outnumber the others; no file is written.
    Refused,
}

impl PredictionCheck {
    /// The rule asks only whether the spanning groups outnumber the others.
    pub fn verdict(&self) -> Verdict {
        if self.spanning > self.groups {
            Verdict::Refused
        } else if self.spanning == 0 && self.groups > 0 {
            Verdict::Confirmed
        } else {
            Verdict::Unconfirmed
        }
    }
}

/// Checks the rebuilt MIC folds against the held-out predictions. On the
/// report's own folds, groups of isolates with bit-identical predictions lie
/// inside one fold; correct folds can also give one value to isolates of two
/// folds, so a group that spans folds is counted, not taken as proof. A fold
/// where every isolate gets the same prediction (a constant fold, as when
/// early stopping keeps no tree and the model predicts its start value) is
/// left out: a group inside one constant fold or spanning only constant
/// folds counts neither way. `folds` gives each isolate's MIC fold,
/// `usize::MAX` for an isolate without an MIC, which is not read.
pub fn prediction_check(predicted: &[f64], folds: &[usize], k: usize) -> PredictionCheck {
    let mut values: Vec<BTreeSet<u64>> = vec![BTreeSet::new(); k];
    let mut by_value: BTreeMap<u64, Vec<usize>> = BTreeMap::new();
    for (i, &f) in folds.iter().enumerate().filter(|&(_, &f)| f < k) {
        let bits = predicted[i].to_bits();
        values[f].insert(bits);
        by_value.entry(bits).or_default().push(i);
    }
    let constant: Vec<bool> = values.iter().map(|v| v.len() == 1).collect();
    let mut check = PredictionCheck {
        groups: 0,
        isolates: 0,
        spanning: 0,
        constant_folds: (0..k).filter(|&f| constant[f]).collect(),
    };
    for members in by_value.values().filter(|m| m.len() >= 2) {
        let spans: BTreeSet<usize> = members.iter().map(|&i| folds[i]).collect();
        if spans.iter().all(|&f| constant[f]) {
            continue;
        }
        if spans.len() > 1 {
            check.spanning += 1;
            continue;
        }
        check.groups += 1;
        check.isolates += members.len();
    }
    check
}

/// A dataset refused alone: its rebuilt MIC folds contradict the report (the
/// prediction check's spanning groups outnumber the others, or a stored
/// threshold differs), or its saved predictions do not give the report's own
/// figure. No file is written for it, a file it had is removed, and the
/// other datasets are written.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Refused {
    pub dataset: String,
    /// I and X of the prediction check.
    pub inside: usize,
    pub spanning: usize,
    /// What disagrees, with the seed the folds were rebuilt at.
    pub reason: String,
}

impl std::fmt::Display for Refused {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}: refused, no baselines file written (I = {} groups inside one fold, X = {} groups spanning folds): {}",
            self.dataset, self.inside, self.spanning, self.reason
        )
    }
}

impl std::error::Error for Refused {}

/// The exact-MIC isolates with a k-mer-identical copy in their training
/// fold: another member of their `kset_duplicates` group that has an MIC
/// and lies in another rebuilt MIC fold. `folds` is as for
/// `prediction_check`.
pub fn exact_with_training_copy(
    mics: &[Option<Mic>],
    folds: &[usize],
    groups: &[Vec<usize>],
) -> usize {
    let has = |i: usize| folds[i] != usize::MAX;
    groups
        .iter()
        .flat_map(|g| {
            g.iter().filter(move |&&i| {
                mics[i].is_some_and(|m| m.censor == Censor::Exact)
                    && g.iter().any(|&j| j != i && has(j) && folds[j] != folds[i])
            })
        })
        .count()
}

/// The same threshold in every field, agreement bit for bit.
fn same_threshold(a: &Option<Breakpoint>, b: &Option<Breakpoint>) -> bool {
    match (a, b) {
        (None, None) => true,
        (Some(a), Some(b)) => {
            a.r_above == b.r_above && a.n == b.n && a.agreement.to_bits() == b.agreement.to_bits()
        }
        _ => false,
    }
}

/// `agree` of `n` exact MICs, with the share as the report divides it.
fn agreement(agree: usize, n: usize) -> Value {
    json!({ "agree": agree, "n": n, "essential_agreement": agree as f64 / n.max(1) as f64 })
}

/// The MIC part of a baselines file. A `Refused` error names the dataset
/// when the rebuilt folds contradict the held-out predictions or a stored
/// threshold, or the saved predictions do not give the report's own
/// agreement; any other error is a report this command cannot read.
fn mic_part(
    id: &str,
    samples: &[Sample],
    cfg: &FoldConfig,
    mic: &MicResult,
    duplicates: &[Vec<usize>],
) -> Result<Value> {
    let n = samples.len();
    let k = cfg.folds;
    let mics: Vec<Option<Mic>> = samples.iter().map(|s| s.mic).collect();
    let has: Vec<usize> = (0..n).filter(|&i| mics[i].is_some()).collect();
    let n_exact = has
        .iter()
        .filter(|&&i| mics[i].is_some_and(|m| m.censor == Censor::Exact))
        .count();
    ensure!(
        mic.n == has.len() && mic.n_exact == n_exact,
        "{id}: the report's MIC model counts {} isolates with an MIC, {} exact, and its samples hold {}, {n_exact} exact",
        mic.n,
        mic.n_exact,
        has.len()
    );
    ensure!(
        mic.predicted.len() == n && has.iter().all(|&i| mic.predicted[i].is_finite()),
        "{id}: the report's MIC predictions do not give one value per isolate with an MIC"
    );
    ensure!(
        mic.fold_breakpoints.len() == k,
        "{id}: the report stores {} MIC training-fold thresholds for {k} folds",
        mic.fold_breakpoints.len()
    );
    ensure!(
        duplicates.iter().flatten().all(|&i| i < n),
        "{id}: the report's kset_duplicates name a sample beyond its {n}"
    );
    let y = Bits::from_indices(n, (0..n).filter(|&i| samples[i].resistant));
    let with_mic = Bits::from_indices(n, has.iter().copied());
    // The draw `evaluate_mic` makes.
    let folds = training_folds(&y, &with_mic, None, k, cfg.seed);
    let ties = prediction_check(&mic.predicted, &folds, k);
    let refuse = |reason: String| -> anyhow::Error {
        Refused {
            dataset: id.to_string(),
            inside: ties.groups,
            spanning: ties.spanning,
            reason,
        }
        .into()
    };
    if ties.verdict() == Verdict::Refused {
        return Err(refuse(format!(
            "the MIC folds rebuilt at seed {} contradict the report's held-out predictions: more groups of isolates with one bit-identical prediction span rebuilt folds than lie inside one",
            cfg.seed
        )));
    }
    let mut checked = 0;
    for (f, (train, _)) in fold_masks(&folds, k).into_iter().enumerate() {
        let pairs: Vec<(Mic, Call)> = train
            .and(&with_mic)
            .iter_ones()
            .map(|i| {
                (
                    mics[i].unwrap(),
                    if y.contains(i) { Call::R } else { Call::S },
                )
            })
            .collect();
        let rebuilt = infer_breakpoint(&pairs);
        let stored = &mic.fold_breakpoints[f];
        if !same_threshold(&rebuilt, stored) {
            return Err(refuse(format!(
                "the MIC folds rebuilt at seed {} contradict the report's stored training-fold threshold of fold {f} (stored {stored:?}, rebuilt {rebuilt:?})",
                cfg.seed
            )));
        }
        checked += stored.is_some() as usize;
    }
    let per_fold = mic_constants(&mics, &mic.predicted, &folds, k);
    let total = |pick: fn(&FoldConstants) -> usize| per_fold.iter().map(pick).sum::<usize>();
    let model = total(|f| f.agree_model);
    let model_share = model as f64 / n_exact.max(1) as f64;
    if model_share.to_bits() != mic.essential_agreement.to_bits() {
        return Err(refuse(format!(
            "the saved MIC predictions put {model} of {n_exact} exact MICs within one dilution ({model_share}), and the report says {}",
            mic.essential_agreement
        )));
    }
    let defined =
        |pick: fn(&FoldConstants) -> Option<i32>| per_fold.iter().all(|f| pick(f).is_some());
    // A fold whose training isolates hold no exact MIC has no constant, and
    // then the dataset has no figure for it.
    let median = if defined(|f| f.median) {
        agreement(total(|f| f.agree_median), n_exact)
    } else {
        Value::Null
    };
    let mode = if defined(|f| f.mode) {
        agreement(total(|f| f.agree_mode), n_exact)
    } else {
        Value::Null
    };
    Ok(json!({
        "folds": k,
        "seed": cfg.seed,
        "prediction_check": {
            "confirmed": ties.verdict() == Verdict::Confirmed,
            "groups": ties.groups,
            "isolates": ties.isolates,
            "spanning": ties.spanning,
            "constant_folds": ties.constant_folds,
        },
        "threshold_check": {
            "confirmed": checked > 0,
            "checked": checked,
        },
        "exact_with_training_copy": exact_with_training_copy(&mics, &folds, duplicates),
        "n": has.len(),
        "n_exact": n_exact,
        "model": agreement(model, n_exact),
        "median": median,
        "mode": mode,
        "fold_constants": per_fold.iter().map(|f| json!({
            "train_exact": f.train_exact,
            "test_exact": f.test_exact,
            "median": f.median,
            "upper_median_differs": f.upper_median_differs,
            "mode": f.mode,
            "mode_tied": f.mode_tied,
        })).collect::<Vec<_>>(),
    }))
}

/// Lineage-identity scores on a split's own folds, and how many genomes
/// had a lineage no training genome of their fold shares. Inside fold f a
/// test genome scores the resistant share of its lineage among the genomes
/// outside f, or, when none of them shares its lineage, their resistant
/// share (`train_prevalence`). The test genome's own label never counts.
pub fn lineage_scores(
    y: &Bits,
    lineage: &[u32],
    test_fold: &[usize],
    k: usize,
) -> Result<(Vec<f64>, usize)> {
    let n = y.len();
    ensure!(
        lineage.len() == n && test_fold.len() == n && test_fold.iter().all(|&f| f < k),
        "every genome needs a lineage and an outer test fold below {k}"
    );
    let mut scores = vec![f64::NAN; n];
    let mut unseen = 0;
    for (f, (train, test)) in fold_masks(test_fold, k).into_iter().enumerate() {
        ensure!(
            train.count() > 0,
            "outer fold {f} leaves no training genome"
        );
        let mut seen: BTreeMap<u32, (u32, u32)> = BTreeMap::new();
        for i in train.iter_ones() {
            let e = seen.entry(lineage[i]).or_default();
            e.0 += y.contains(i) as u32;
            e.1 += 1;
        }
        let prior = train_prevalence(y, &train);
        for i in test.iter_ones() {
            scores[i] = match seen.get(&lineage[i]) {
                Some(&(r, of)) => r as f64 / of as f64,
                None => {
                    unseen += 1;
                    prior
                }
            };
        }
    }
    Ok((scores, unseen))
}

/// The lineage part of a baselines file: one entry per split with folds.
fn lineage_part(id: &str, samples: &[Sample], splits: &[Split]) -> Result<Value> {
    let n = samples.len();
    let y = Bits::from_indices(n, (0..n).filter(|&i| samples[i].resistant));
    let truth: Vec<bool> = samples.iter().map(|s| s.resistant).collect();
    let lineage: Vec<u32> = samples.iter().map(|s| s.lineage).collect();
    let mut out = Vec::new();
    for split in splits.iter().filter(|s| !s.test_fold.is_empty()) {
        let (scores, unseen) = lineage_scores(&y, &lineage, &split.test_fold, split.folds)
            .with_context(|| format!("{id}: the {} split", split.scheme))?;
        // The learners' cutoff and pooling (`pool_folds`).
        let m = binary_metrics(&truth, &scores, 0.5, false);
        out.push(json!({
            "scheme": split.scheme,
            "folds": split.folds,
            "n": m.n,
            "tp": m.tp,
            "fp": m.fp,
            "tn": m.tn,
            "fn": m.fn_,
            "balanced_accuracy": m.balanced_accuracy,
            "sensitivity": m.sensitivity,
            "specificity": m.specificity,
            "unseen": unseen,
            "unseen_share": unseen as f64 / n.max(1) as f64,
        }));
    }
    Ok(Value::Array(out))
}

/// The baselines file for one dataset report. It refuses a report of
/// another evaluation revision than this binary's.
pub fn for_report(bytes: &[u8]) -> Result<Value> {
    let report: Report = serde_json::from_slice(bytes).context("parse a dataset report")?;
    let id = &report.id;
    let revision = report.evaluation.as_ref().and_then(|e| e.revision);
    if revision != Some(EVALUATION_REVISION as u64) {
        bail!(
            "{id}: the report is of evaluation revision {}, and this grm reads revision {EVALUATION_REVISION}; rebuild the report or use the grm that wrote it",
            revision.map_or_else(|| "none".to_string(), |r| r.to_string())
        );
    }
    let mic = match &report.mic {
        None => Value::Null,
        Some(mic) => {
            let cfg = report
                .eval_config
                .as_ref()
                .with_context(|| format!("{id}: the report records no eval_config"))?;
            mic_part(id, &report.samples, cfg, mic, &report.kset_duplicates)?
        }
    };
    Ok(json!({
        "schema": SCHEMA,
        "dataset": id,
        "report_sha256": crate::research::sha256(bytes),
        "evaluation_revision": EVALUATION_REVISION,
        "mic": mic,
        "lineage": lineage_part(id, &report.samples, &report.splits)?,
    }))
}

/// The summary file for the pages that show many datasets, written beside
/// `index.json` rather than in `baselines/`: no dataset id can name it (an
/// id is a file name inside `baselines/`), and the removal step, which reads
/// `baselines/` alone, never reaches it.
pub const SUMMARY: &str = "baselines.json";

/// Version of the summary file's layout.
pub const SUMMARY_SCHEMA: u32 = 1;

/// The summary of one run: per written dataset with a constant MIC, only
/// what the overview, the datasets table, the methods page and the slide
/// show (the model's figure, the constant's and whether the held-out
/// predictions confirm the folds), with the report's SHA-256, pinned to the
/// SHA-256 of the `index.json` bytes the run read. A new index makes it
/// stale, and the pages then refuse it whole.
pub fn summary(index_bytes: &[u8], files: &[Value]) -> Value {
    let datasets: Vec<Value> = files
        .iter()
        .filter(|f| !f["mic"].is_null() && !f["mic"]["median"].is_null())
        .map(|f| {
            let m = &f["mic"];
            json!({
                "dataset": f["dataset"],
                "report_sha256": f["report_sha256"],
                "model": m["model"]["essential_agreement"],
                "constant": m["median"]["essential_agreement"],
                "confirmed": m["prediction_check"]["confirmed"],
            })
        })
        .collect();
    json!({
        "schema": SUMMARY_SCHEMA,
        "evaluation_revision": EVALUATION_REVISION,
        "index_sha256": crate::research::sha256(index_bytes),
        "datasets": datasets,
    })
}

/// What one run wrote.
#[derive(Debug, PartialEq, Eq)]
pub struct Written {
    pub files: usize,
    /// Files with an MIC part.
    pub mic: usize,
    /// Of those, the ones whose rebuilt folds the held-out predictions
    /// confirm (`prediction_check`).
    pub confirmed: usize,
    /// Exact-MIC test isolates with a k-mer-identical copy in their training
    /// fold, over every file, and the datasets that have one.
    pub copies: usize,
    pub copy_datasets: usize,
    /// The datasets refused alone, in index order.
    pub refused: Vec<Refused>,
}

/// An index id is a file name: not empty, with no path separator and no
/// `..`, so `<id>.json` stays inside `datasets/` and `baselines/`.
fn plain_id(id: &str) -> bool {
    !id.is_empty() && !id.contains(['/', '\\', '\0']) && !id.contains("..")
}

/// Writes `bytes` to `path` through a pending file beside it, as
/// `report::publish` does, and removes the pending file when the write or
/// the rename fails, so a failed run leaves no pending file of its own. The
/// error names the path.
fn publish(path: &Path, bytes: &[u8]) -> Result<()> {
    let pending = path.with_extension(format!("json.{}.tmp", std::process::id()));
    if let Err(e) = std::fs::write(&pending, bytes).and_then(|()| std::fs::rename(&pending, path)) {
        let _ = std::fs::remove_file(&pending);
        return Err(e).with_context(|| format!("write {}", path.display()));
    }
    Ok(())
}

/// Writes `baselines/<id>.json` under `report_dir` for every dataset that
/// its `index.json` lists and the command does not refuse, then `SUMMARY`
/// beside `index.json`, then removes every other regular `*.json` and
/// `*.tmp` file in `baselines/`, as `grm report`'s prune removes its own
/// (a refused dataset's old file among them); a link, a directory or any
/// other entry stays. It refuses a `baselines/` that is a symbolic link or
/// not a directory, and an index id that is not a plain file name, before
/// any report is read or any file written. Every file is computed before
/// any is written, so a report the command cannot read leaves the folder as
/// it was. It holds the report directory's lock, as `grm report` does, so it
/// never reads a build half done.
pub fn write(report_dir: &Path) -> Result<Written> {
    let lock = std::fs::File::open(report_dir)
        .with_context(|| format!("open {}", report_dir.display()))?;
    lock.lock()
        .with_context(|| format!("lock {}", report_dir.display()))?;
    let out = report_dir.join("baselines");
    let exists = match std::fs::symlink_metadata(&out) {
        Ok(m) if m.file_type().is_dir() => true,
        Ok(_) => bail!(
            "{} is a symbolic link or not a directory; nothing was read or written",
            out.display()
        ),
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => false,
        Err(e) => return Err(e).with_context(|| format!("read {}", out.display())),
    };
    let index_path = report_dir.join("index.json");
    let index_bytes = std::fs::read(&index_path)
        .with_context(|| format!("read {}; run `grm report` first", index_path.display()))?;
    let index: Value = serde_json::from_slice(&index_bytes)
        .with_context(|| format!("parse {}", index_path.display()))?;
    let mut ids = Vec::new();
    for row in index["datasets"].as_array().into_iter().flatten() {
        let id = row["id"]
            .as_str()
            .context("index.json lists a dataset with no id")?;
        ensure!(
            plain_id(id),
            "index.json lists the dataset id {id:?}, which holds a path separator or `..`; nothing was written"
        );
        ids.push(id);
    }
    let mut files = Vec::new();
    let mut refused = Vec::new();
    for id in ids {
        let path = report_dir.join("datasets").join(format!("{id}.json"));
        let bytes = std::fs::read(&path).with_context(|| format!("read {}", path.display()))?;
        let file = match for_report(&bytes) {
            Ok(file) => file,
            Err(e) => match e.downcast::<Refused>() {
                Ok(r) => {
                    refused.push(r);
                    continue;
                }
                Err(e) => return Err(e.context(format!("read {}", path.display()))),
            },
        };
        if file["dataset"] != id {
            bail!("{} records the id {}", path.display(), file["dataset"]);
        }
        files.push((format!("{id}.json"), file));
    }
    if !exists {
        std::fs::create_dir(&out).with_context(|| format!("create {}", out.display()))?;
    }
    let mut kept = BTreeSet::new();
    let mut written = Written {
        files: 0,
        mic: 0,
        confirmed: 0,
        copies: 0,
        copy_datasets: 0,
        refused,
    };
    for (name, file) in &files {
        publish(&out.join(name), &serde_json::to_vec(file)?)?;
        written.files += 1;
        if !file["mic"].is_null() {
            written.mic += 1;
            written.confirmed += (file["mic"]["prediction_check"]["confirmed"] == true) as usize;
            let copies = file["mic"]["exact_with_training_copy"]
                .as_u64()
                .unwrap_or(0) as usize;
            written.copies += copies;
            written.copy_datasets += (copies > 0) as usize;
        }
        kept.insert(name.clone());
    }
    let shown: Vec<Value> = files.into_iter().map(|(_, f)| f).collect();
    publish(
        &report_dir.join(SUMMARY),
        &serde_json::to_vec(&summary(&index_bytes, &shown))?,
    )?;
    for entry in std::fs::read_dir(&out).with_context(|| format!("read {}", out.display()))? {
        let entry = entry.with_context(|| format!("read {}", out.display()))?;
        let regular = entry
            .file_type()
            .with_context(|| format!("read {}", entry.path().display()))?
            .is_file();
        let stale = entry
            .file_name()
            .to_str()
            .is_some_and(|n| (n.ends_with(".json") || n.ends_with(".tmp")) && !kept.contains(n));
        if regular && stale {
            std::fs::remove_file(entry.path())
                .with_context(|| format!("remove {}", entry.path().display()))?;
        }
    }
    Ok(written)
}

#[cfg(test)]
mod tests {
    //! Every report here is synthetic: made-up isolates, labels, MICs and
    //! lineages written by the tests, not BV-BRC data.
    use super::*;

    fn exact(log2: i32) -> Option<Mic> {
        Some(Mic {
            log2,
            censor: Censor::Exact,
        })
    }

    #[test]
    fn lower_median_takes_the_lower_middle_value() {
        assert_eq!(lower_median(&[]), None);
        assert_eq!(lower_median(&[3]), Some(3));
        assert_eq!(lower_median(&[5, 1, 3]), Some(3));
        // Even count: the two middle values are 2 and 4; the lower is 2,
        // the ceil(4/2) = 2nd smallest.
        assert_eq!(lower_median(&[4, 1, 2, 9]), Some(2));
        assert_eq!(lower_median(&[0, 0, 1, 1, 1, 1]), Some(1));
        assert_eq!(lower_median(&[-3, -1]), Some(-3));
    }

    #[test]
    fn mode_breaks_a_tie_toward_the_lowest_dilution() {
        assert_eq!(mode(&[]), None);
        assert_eq!(mode(&[2, 2, 5]), Some((2, false)));
        assert_eq!(mode(&[5, 2, 5]), Some((5, false)));
        // 1 and 4 each twice: the lower wins and the tie is recorded.
        assert_eq!(mode(&[4, 1, 4, 1, 7]), Some((1, true)));
        assert_eq!(mode(&[3, -2, 3, -2]), Some((-2, true)));
    }

    #[test]
    fn mic_constants_come_from_the_training_fold_and_score_its_test_isolates() {
        // Synthetic: fold 0 holds the low MICs, fold 1 the high ones, and
        // isolate 8 has no MIC (fold usize::MAX, as `training_folds` leaves
        // it). A constant taken from a fold's own test isolates would land
        // on them; the training fold's lands two or more dilutions away.
        let mics = vec![
            exact(0),
            exact(0),
            exact(1),
            Some(Mic {
                log2: 9,
                censor: Censor::Above,
            }),
            exact(4),
            exact(5),
            exact(5),
            exact(6),
            None,
        ];
        let folds = [0, 0, 0, 0, 1, 1, 1, 1, usize::MAX];
        let predicted = [0.4, 2.5, -0.6, 9.0, 4.0, 3.5, 7.0, 6.0, f64::NAN];
        let got = mic_constants(&mics, &predicted, &folds, 2);
        // Fold 0 trains on fold 1's exact MICs 4, 5, 5, 6 (the censored one
        // is left out): lower median 5, mode 5. Its exact test MICs 0, 0, 1
        // are all more than one dilution from 5.
        assert_eq!(
            got[0],
            FoldConstants {
                train_exact: 4,
                test_exact: 3,
                median: Some(5),
                upper_median_differs: false,
                mode: Some(5),
                mode_tied: false,
                // round(0.4) = 0 agrees with 0; round(2.5) = 3 is three
                // from 0, and round(-0.6) = -1 two from 1.
                agree_model: 1,
                agree_median: 0,
                agree_mode: 0,
            }
        );
        // Fold 1 trains on 0, 0, 1: lower median 0, mode 0; its test MICs
        // 4, 5, 5, 6 are all far from 0.
        assert_eq!(got[1].median, Some(0));
        assert_eq!(got[1].mode, Some(0));
        assert_eq!((got[1].agree_median, got[1].agree_mode), (0, 0));
        // round(4.0) = 4 for 4, round(3.5) = 4 for 5, round(7.0) = 7 for 5
        // (no), round(6.0) = 6 for 6.
        assert_eq!(got[1].agree_model, 3);
    }

    #[test]
    fn mic_constants_use_the_lower_median_and_the_lowest_mode() {
        // Synthetic: fold 1's training MICs are fold 0's: 1, 1, 3, 3. The
        // lower median is 1 and the upper 3; the mode ties 1 and 3 and goes
        // to 1. The test MICs 0, 0 sit within one dilution of 1 and not of 3.
        let mics = vec![exact(1), exact(1), exact(3), exact(3), exact(0), exact(0)];
        let folds = [0, 0, 0, 0, 1, 1];
        let predicted = [0.0; 6];
        let got = &mic_constants(&mics, &predicted, &folds, 2)[1];
        assert_eq!(got.median, Some(1));
        assert!(got.upper_median_differs);
        assert_eq!((got.mode, got.mode_tied), (Some(1), true));
        assert_eq!((got.agree_median, got.agree_mode), (2, 2));
    }

    // ------------------------------------------------- synthetic reports

    const N_MIC: usize = 200;

    /// The 4 isolates whose label disagrees with their MIC in `synthetic`.
    const FLIPS: [usize; 4] = [3, 50, 101, 177];

    /// A synthetic report: `N_MIC` isolates with an MIC (exact, except every
    /// 7th with a low MIC, which is a bound "at most") and 20 without. An
    /// MIC above 2 is resistant, except the 4 `FLIPS`, so each training
    /// fold's threshold has an agreement that depends on which isolates it
    /// holds.
    fn synthetic(seed: u64) -> Value {
        build(seed, &|i, log2| {
            log2.map_or(i % 2 == 0, |v| (v > 2) != FLIPS.contains(&i))
        })
    }

    /// A synthetic report with labels from `resistant(isolate, log2 MIC)`.
    /// The model predicts each MIC exactly, except every 5th isolate, which
    /// it misses by two dilutions, plus a hundredth of a dilution for each
    /// step of its MIC fold at the report's seed: as one boosted model per
    /// fold gives, isolates with the same prediction share a fold. Lineages
    /// run in blocks of 10 isolates; the random split deals isolates round
    /// robin into 5 folds and the lineage split holds out whole blocks.
    fn build(seed: u64, resistant: &dyn Fn(usize, Option<i32>) -> bool) -> Value {
        let n = N_MIC + 20;
        let mut samples = Vec::new();
        for i in 0..n {
            let log2 = (i < N_MIC).then(|| (i % 8) as i32 - 1);
            let mic = log2.map_or(Value::Null, |v| {
                let censor = if i % 7 == 6 && v <= 1 {
                    "AtMost"
                } else {
                    "Exact"
                };
                json!({ "log2": v, "censor": censor })
            });
            samples.push(json!({ "id": format!("g{i}"), "resistant": resistant(i, log2), "mic": mic, "lineage": i / 10 }));
        }
        let mut report = json!({
            "id": "synthetic__fixturemycin",
            "evaluation": { "revision": EVALUATION_REVISION },
            "eval_config": { "folds": 5, "seed": seed },
            "samples": samples,
            "splits": [
                { "scheme": "random", "folds": 5, "n_groups": 22, "test_fold": (0..n).map(|i| i % 5).collect::<Vec<_>>(), "learners": [] },
                { "scheme": "lineage", "folds": 5, "n_groups": 22, "test_fold": (0..n).map(|i| (i / 10) % 5).collect::<Vec<_>>(), "learners": [] },
                { "scheme": "time", "folds": 5, "n_groups": 0, "learners": [], "unavailable_reason": "synthetic" },
            ],
            "mic": {
                "n": N_MIC, "n_exact": 0, "essential_agreement": 0.0, "essential_agreement_all": 0.0,
                "categorical": null, "breakpoint": null, "fold_breakpoints": [],
                "confusion": [], "predicted": [],
            },
        });
        let folds = folds_at(&report, seed);
        report["mic"]["predicted"] = (0..n)
            .map(|i| {
                if i < N_MIC {
                    let v = (i % 8) as f64 - 1.0;
                    let miss = if i % 5 == 0 { 2.0 } else { 0.0 };
                    json!(v + miss + 0.01 * (folds[i] + 1) as f64)
                } else {
                    Value::Null
                }
            })
            .collect();
        report["mic"]["fold_breakpoints"] = json!(stored_thresholds(&report, seed));
        stamp(&mut report);
        report
    }

    /// Stores the model's agreement, as the fit would, from the report's
    /// predictions.
    fn stamp(report: &mut Value) {
        let mics = parse_mics(report);
        let exact: Vec<(usize, i32)> = mics
            .iter()
            .enumerate()
            .filter_map(|(i, m)| m.filter(|m| m.censor == Censor::Exact).map(|m| (i, m.log2)))
            .collect();
        let agree = exact
            .iter()
            .filter(|&&(i, v)| {
                within_one(
                    report["mic"]["predicted"][i].as_f64().unwrap().round() as i32,
                    v,
                )
            })
            .count();
        report["mic"]["n_exact"] = json!(exact.len());
        report["mic"]["essential_agreement"] = json!(agree as f64 / exact.len() as f64);
    }

    fn parse_mics(report: &Value) -> Vec<Option<Mic>> {
        report["samples"]
            .as_array()
            .unwrap()
            .iter()
            .map(|s| serde_json::from_value(s["mic"].clone()).unwrap())
            .collect()
    }

    fn labels(report: &Value) -> Bits {
        let samples = report["samples"].as_array().unwrap();
        Bits::from_indices(
            samples.len(),
            (0..samples.len()).filter(|&i| samples[i]["resistant"] == true),
        )
    }

    /// The MIC folds `evaluate_mic` draws for the report at `seed`.
    fn folds_at(report: &Value, seed: u64) -> Vec<usize> {
        let mics = parse_mics(report);
        let with_mic =
            Bits::from_indices(mics.len(), (0..mics.len()).filter(|&i| mics[i].is_some()));
        training_folds(&labels(report), &with_mic, None, 5, seed)
    }

    fn stored_thresholds(report: &Value, seed: u64) -> Vec<Option<Breakpoint>> {
        let mics = parse_mics(report);
        let y = labels(report);
        let with_mic =
            Bits::from_indices(mics.len(), (0..mics.len()).filter(|&i| mics[i].is_some()));
        fold_masks(&folds_at(report, seed), 5)
            .into_iter()
            .map(|(train, _)| {
                let pairs: Vec<(Mic, Call)> = train
                    .and(&with_mic)
                    .iter_ones()
                    .map(|i| {
                        (
                            mics[i].unwrap(),
                            if y.contains(i) { Call::R } else { Call::S },
                        )
                    })
                    .collect();
                infer_breakpoint(&pairs)
            })
            .collect()
    }

    /// Groups of two or more isolates in folds that `keep` keeps, by fold
    /// and predicted dilution as `build` predicts it, and their isolates.
    fn groups_by_fold(folds: &[usize], keep: impl Fn(usize) -> bool) -> (usize, usize) {
        let mut count: BTreeMap<(usize, usize), usize> = BTreeMap::new();
        for i in (0..N_MIC).filter(|&i| keep(folds[i])) {
            let dilution = i % 8 + if i % 5 == 0 { 2 } else { 0 };
            *count.entry((folds[i], dilution)).or_default() += 1;
        }
        let groups = count.values().filter(|&&c| c >= 2);
        (groups.clone().count(), groups.sum())
    }

    fn predictions(report: &Value) -> Vec<f64> {
        report["mic"]["predicted"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_f64().unwrap_or(f64::NAN))
            .collect()
    }

    fn bytes(v: &Value) -> Vec<u8> {
        serde_json::to_vec(v).unwrap()
    }

    #[test]
    fn baselines_file_pins_the_report_and_confirms_the_rebuilt_folds() {
        let report = synthetic(42);
        let raw = bytes(&report);
        let file = for_report(&raw).unwrap();
        assert_eq!(file["schema"], SCHEMA);
        assert_eq!(file["dataset"], "synthetic__fixturemycin");
        assert_eq!(file["report_sha256"], crate::research::sha256(&raw));
        assert_eq!(file["evaluation_revision"], EVALUATION_REVISION);
        let mic = &file["mic"];
        // Isolates with one prediction share a fold: each predicted
        // dilution in each of the five folds is one group, counted here
        // from the folds and dilutions rather than the predictions' bits.
        let (groups, isolates) = groups_by_fold(&folds_at(&report, 42), |_| true);
        assert!(groups >= 40, "{groups}");
        assert_eq!(
            mic["prediction_check"],
            json!({ "confirmed": true, "groups": groups, "isolates": isolates, "spanning": 0, "constant_folds": [] })
        );
        // Every fold's training isolates support a threshold, so all five
        // stored thresholds are checked.
        assert_eq!(
            mic["threshold_check"],
            json!({ "confirmed": true, "checked": 5 })
        );
        assert_eq!(mic["exact_with_training_copy"], 0);
        assert_eq!(
            (mic["folds"].as_u64(), mic["seed"].as_u64()),
            (Some(5), Some(42))
        );
        assert_eq!(mic["n"], N_MIC);
        let n_exact = mic["n_exact"].as_u64().unwrap() as usize;
        assert_eq!(
            n_exact,
            (0..N_MIC).filter(|i| i % 7 != 6 || i % 8 > 2).count()
        );
        // The model's recount is the report's own figure.
        assert_eq!(
            mic["model"]["essential_agreement"],
            report["mic"]["essential_agreement"]
        );
        assert_eq!(mic["model"]["n"], n_exact);
        // The constants are scored on the same exact MICs, and each count is
        // the sum over the folds of the constant's reach.
        for key in ["median", "mode"] {
            assert_eq!(mic[key]["n"], n_exact);
            let agree = mic[key]["agree"].as_u64().unwrap() as f64;
            assert_eq!(
                mic[key]["essential_agreement"].as_f64().unwrap(),
                agree / n_exact as f64
            );
        }
        // Hand check of the constants: the MICs cycle -1..=6, so every
        // training fold's exact values are spread over eight dilutions with
        // a lower median of 2 or 3, and a constant reaches three of the
        // eight dilutions.
        for f in mic["fold_constants"].as_array().unwrap() {
            let median = f["median"].as_i64().unwrap();
            assert!((2..=3).contains(&median), "{f}");
            assert_eq!(
                f["train_exact"].as_u64().unwrap() + f["test_exact"].as_u64().unwrap(),
                n_exact as u64
            );
        }
        let share = mic["median"]["essential_agreement"].as_f64().unwrap();
        assert!((0.30..0.45).contains(&share), "{share}");
    }

    #[test]
    fn the_folds_are_rebuilt_at_the_report_s_seed_and_a_contradiction_refuses_the_dataset() {
        // The predictions and thresholds come from the folds of seed 42.
        // Read as a report of seed 43, the rebuilt folds contradict them,
        // and the error is the dataset's own refusal, not a stop of the run.
        let mut report = synthetic(42);
        report["eval_config"]["seed"] = json!(43);
        let err = for_report(&bytes(&report)).unwrap_err();
        let refused = err.downcast_ref::<Refused>().expect("a refusal");
        assert_eq!(refused.dataset, "synthetic__fixturemycin");
        assert!(refused.spanning > refused.inside, "{refused}");
        let err = err.to_string();
        assert!(err.contains("synthetic__fixturemycin"), "{err}");
        assert!(err.contains("contradict"), "{err}");
        // A stored threshold one dilution off.
        let mut report = synthetic(42);
        let r = report["mic"]["fold_breakpoints"][2]["r_above"]
            .as_i64()
            .unwrap();
        report["mic"]["fold_breakpoints"][2]["r_above"] = json!(r + 1);
        let err = for_report(&bytes(&report)).unwrap_err();
        let refused = err.downcast_ref::<Refused>().expect("a refusal");
        assert!(refused.reason.contains("fold 2"), "{refused}");
        assert_eq!(refused.spanning, 0);
        // A threshold stored as none where the rebuilt fold fits one.
        let mut report = synthetic(42);
        report["mic"]["fold_breakpoints"][0] = Value::Null;
        assert!(for_report(&bytes(&report)).is_err());
    }

    #[test]
    fn folds_with_no_stored_threshold_are_still_confirmed_by_the_predictions() {
        // Synthetic: every label resistant but 9, so no training fold has the
        // 10 susceptible isolates a threshold needs and every stored
        // threshold is none, which the rebuilt folds reproduce.
        let report = build(42, &|i, _| i >= 9);
        assert!(
            report["mic"]["fold_breakpoints"]
                .as_array()
                .unwrap()
                .iter()
                .all(Value::is_null)
        );
        let file = for_report(&bytes(&report)).unwrap();
        assert_eq!(
            file["mic"]["threshold_check"],
            json!({ "confirmed": false, "checked": 0 })
        );
        assert_eq!(file["mic"]["prediction_check"]["confirmed"], true);
        assert!(file["mic"]["median"]["agree"].is_u64());
    }

    #[test]
    fn a_wrong_seed_fails_the_prediction_check_where_it_passes_the_threshold_check() {
        // Synthetic: no label disagrees with its MIC, so every training
        // fold's threshold separates its isolates without error, and any
        // stratified draw with the same fold sizes fits the same thresholds.
        let report = build(42, &|i, log2| log2.map_or(i % 2 == 0, |v| v > 2));
        let stored = stored_thresholds(&report, 42);
        assert!(stored.iter().all(|b| b.is_some_and(|b| b.agreement == 1.0)));
        assert_eq!(stored_thresholds(&report, 43), stored);
        // The held-out predictions tell the two draws apart.
        let predicted = predictions(&report);
        let right = prediction_check(&predicted, &folds_at(&report, 42), 5);
        let wrong = prediction_check(&predicted, &folds_at(&report, 43), 5);
        assert_eq!((right.verdict(), right.spanning), (Verdict::Confirmed, 0));
        assert_eq!(wrong.verdict(), Verdict::Refused, "{wrong:?}");
        let mut wrong = report.clone();
        wrong["eval_config"]["seed"] = json!(43);
        let err = for_report(&bytes(&wrong)).unwrap_err().to_string();
        assert!(
            err.contains("synthetic__fixturemycin") && err.contains("held-out predictions"),
            "{err}"
        );
        assert_eq!(
            for_report(&bytes(&report)).unwrap()["mic"]["prediction_check"]["confirmed"],
            true
        );
    }

    #[test]
    fn prediction_check_keeps_each_group_in_one_fold_except_across_constant_folds() {
        // Synthetic. Two groups inside folds 0 and 1, and two singletons;
        // isolates without an MIC are not read.
        let folds = [0, 0, 0, 1, 1, 1, usize::MAX, usize::MAX];
        let predicted = [0.5, 0.5, 0.6, 0.7, 0.7, 0.9, f64::NAN, f64::NAN];
        let check = prediction_check(&predicted, &folds, 2);
        assert_eq!(
            check,
            PredictionCheck {
                groups: 2,
                isolates: 4,
                spanning: 0,
                constant_folds: vec![]
            }
        );
        assert_eq!(check.verdict(), Verdict::Confirmed);
        // One value in two folds that each predict more than one value: a
        // spanning group and no group inside a fold, so X > I.
        let check = prediction_check(&[0.5, 0.6, 0.5, 0.7], &[0, 0, 1, 1], 2);
        assert_eq!((check.groups, check.spanning), (0, 1));
        assert_eq!(check.verdict(), Verdict::Refused);
        // Folds 0 and 1 predict one start value each, the same: no tree.
        // Their group is left out, and fold 2's group counts.
        let folds = [0, 0, 1, 1, 2, 2, 2];
        let predicted = [-1.0, -1.0, -1.0, -1.0, 0.25, 0.25, 1.5];
        assert_eq!(
            prediction_check(&predicted, &folds, 3),
            PredictionCheck {
                groups: 1,
                isolates: 2,
                spanning: 0,
                constant_folds: vec![0, 1]
            }
        );
        // The same value in fold 2, which keeps trees: a spanning group.
        let predicted = [-1.0, -1.0, -1.0, -1.0, -1.0, 0.25, 1.5];
        let check = prediction_check(&predicted, &folds, 3);
        assert_eq!((check.groups, check.spanning), (0, 1));
        // A group inside one constant fold is left out too; fold 1's counts.
        let check = prediction_check(&[2.0, 2.0, 3.0, 3.0, 4.0], &[0, 0, 1, 1, 1], 2);
        assert_eq!(
            (check.groups, check.spanning, check.constant_folds),
            (1, 0, vec![0])
        );
        // Only constant folds: no group either way, so the folds are
        // unconfirmed and the file is still written.
        let check = prediction_check(&[2.0, 2.0, 2.0, 2.0], &[0, 0, 1, 1], 2);
        assert_eq!(
            check,
            PredictionCheck {
                groups: 0,
                isolates: 0,
                spanning: 0,
                constant_folds: vec![0, 1]
            }
        );
        assert_eq!(check.verdict(), Verdict::Unconfirmed);
    }

    #[test]
    fn folds_that_predict_one_start_value_are_left_out_and_unconfirmed_without_another_group() {
        // Synthetic: folds 0 and 1 kept no tree and predict -1.0 for every
        // isolate; the other three keep their groups.
        let mut report = synthetic(42);
        let folds = folds_at(&report, 42);
        for (i, &f) in folds.iter().enumerate() {
            if f < 2 {
                report["mic"]["predicted"][i] = json!(-1.0);
            }
        }
        stamp(&mut report);
        let file = for_report(&bytes(&report)).unwrap();
        let check = &file["mic"]["prediction_check"];
        assert_eq!(check["constant_folds"], json!([0, 1]));
        let (groups, _) = groups_by_fold(&folds, |f| f >= 2);
        assert_eq!(
            (check["groups"].as_u64(), check["confirmed"].as_bool()),
            (Some(groups as u64), Some(true))
        );
        // An isolate of fold 2 predicted -1.0 too joins the constant folds'
        // group: one spanning group beside the inside ones, so the file is
        // written unconfirmed.
        let i = folds.iter().position(|&f| f == 2).unwrap();
        report["mic"]["predicted"][i] = json!(-1.0);
        stamp(&mut report);
        let file = for_report(&bytes(&report)).unwrap();
        let check = &file["mic"]["prediction_check"];
        assert_eq!(
            (check["spanning"].as_u64(), check["confirmed"].as_bool()),
            (Some(1), Some(false))
        );
    }

    #[test]
    fn copies_in_training_count_exact_test_isolates_whose_copy_has_an_mic_in_another_fold() {
        let mut report = synthetic(42);
        let folds = folds_at(&report, 42);
        let mics = parse_mics(&report);
        let is_exact = |i: usize| mics[i].is_some_and(|m| m.censor == Censor::Exact);
        let exact_in = |f: usize| -> Vec<usize> {
            (0..N_MIC)
                .filter(|&i| folds[i] == f && is_exact(i))
                .collect()
        };
        let (a, b) = (exact_in(0)[0], exact_in(1)[0]);
        let same = exact_in(2)[..2].to_vec();
        let censored = (0..N_MIC).find(|&i| !is_exact(i) && folds[i] != 3).unwrap();
        let other = exact_in(3)[0];
        // Synthetic groups: a pair split over folds 0 and 1 (both count), a
        // pair inside fold 2 (neither), an exact isolate whose copy has no
        // MIC (none), and an exact isolate whose copy has a censored MIC in
        // another fold (the exact one counts).
        report["kset_duplicates"] =
            json!([[a, b], same, [exact_in(4)[0], N_MIC + 1], [censored, other]]);
        let file = for_report(&bytes(&report)).unwrap();
        assert_eq!(file["mic"]["exact_with_training_copy"], 3);
        report["kset_duplicates"] = json!([[0, N_MIC + 20]]);
        assert!(for_report(&bytes(&report)).is_err());
    }

    #[test]
    fn predictions_that_do_not_give_the_report_s_agreement_are_refused() {
        let mut report = synthetic(42);
        report["mic"]["predicted"][1] = json!(-5.0);
        let err = for_report(&bytes(&report)).unwrap_err().to_string();
        assert!(
            err.contains("synthetic__fixturemycin") && err.contains("report says"),
            "{err}"
        );
    }

    #[test]
    fn a_report_without_an_mic_model_has_no_mic_part() {
        let mut report = synthetic(42);
        report["mic"] = Value::Null;
        let file = for_report(&bytes(&report)).unwrap();
        assert!(file["mic"].is_null());
        assert_eq!(file["lineage"].as_array().unwrap().len(), 2);
    }

    #[test]
    fn a_report_of_another_evaluation_revision_is_refused() {
        for revision in [
            json!(EVALUATION_REVISION - 1),
            json!(EVALUATION_REVISION + 1),
            Value::Null,
        ] {
            let mut report = synthetic(42);
            report["evaluation"]["revision"] = revision;
            let err = for_report(&bytes(&report)).unwrap_err().to_string();
            assert!(err.contains("evaluation revision"), "{err}");
        }
        let mut report = synthetic(42);
        report.as_object_mut().unwrap().remove("evaluation");
        assert!(for_report(&bytes(&report)).is_err());
    }

    #[test]
    fn lineage_scores_read_a_seen_lineage_and_fall_back_on_an_unseen_one() {
        // Synthetic, two folds. Fold 0 tests genomes 0 to 2, fold 1 tests 3
        // to 7.
        //   genome:   0  1  2  3  4  5  6  7
        //   lineage:  7  7  9  7  7  8  8  8
        //   resistant 1  1  0  1  0  0  0  1
        let y = Bits::from_indices(8, [0, 1, 3, 7]);
        let lineage = [7, 7, 9, 7, 7, 8, 8, 8];
        let test_fold = [0, 0, 0, 1, 1, 1, 1, 1];
        let (scores, unseen) = lineage_scores(&y, &lineage, &test_fold, 2).unwrap();
        // Fold 0 trains on 3 to 7. Lineage 7 there: genomes 3 (R) and 4 (S),
        // share 0.5, so genomes 0 and 1 score 0.5. Lineage 9 is unseen, so
        // genome 2 takes the training prevalence, 2 of 5 = 0.4, not 0.
        assert_eq!(scores[..3], [0.5, 0.5, 0.4]);
        // Fold 1 trains on 0 to 2. Lineage 7 there: 0 and 1, both R, share
        // 1.0; the test genome's own label (4 is S) does not count. Lineage 8
        // is unseen, so 5, 6 and 7 take the prevalence 2 of 3.
        assert_eq!(scores[3..5], [1.0, 1.0]);
        assert_eq!(scores[5..], [2.0 / 3.0; 3]);
        assert_eq!(unseen, 4);
        // Called at 0.5 or more: genomes 0, 1, 3, 4, 5, 6, 7 are called
        // resistant and genome 2 susceptible. R: 0, 1, 3, 7 all called, so
        // sensitivity 1; S: 2 right, 4, 5, 6 wrong, so specificity 1/4.
        let truth: Vec<bool> = (0..8).map(|i| y.contains(i)).collect();
        let m = binary_metrics(&truth, &scores, 0.5, false);
        assert_eq!((m.tp, m.fp, m.tn, m.fn_), (4, 3, 1, 0));
        assert_eq!(m.balanced_accuracy, 0.625);
    }

    #[test]
    fn a_lineage_share_never_counts_the_test_genome() {
        // Synthetic: each genome shares its lineage with one genome of the
        // other fold and the other class. Its training share is that
        // genome's label alone, 1.0 or 0.0; counting its own label too would
        // give 0.5 for every genome.
        let y = Bits::from_indices(4, [1, 2]);
        let (scores, unseen) = lineage_scores(&y, &[1, 1, 2, 2], &[0, 1, 0, 1], 2).unwrap();
        assert_eq!(scores, [1.0, 0.0, 0.0, 1.0]);
        assert_eq!(unseen, 0);
    }

    #[test]
    fn lineage_split_has_nothing_to_read() {
        // On the synthetic lineage split every test genome's lineage is out
        // of its training folds, so every genome takes its fold's prevalence.
        let report = synthetic(42);
        let file = for_report(&bytes(&report)).unwrap();
        let splits = file["lineage"].as_array().unwrap();
        assert_eq!(splits.len(), 2, "the unavailable split has no entry");
        let lineage = &splits[1];
        assert_eq!(lineage["scheme"], "lineage");
        assert_eq!(lineage["unseen"], 220);
        assert_eq!(lineage["unseen_share"], 1.0);
        let random = &splits[0];
        assert_eq!(random["scheme"], "random");
        // Lineage blocks of 10 dealt round robin: every lineage is seen.
        assert_eq!(random["unseen"], 0);
        assert_eq!(random["n"], 220);
        let ba = random["balanced_accuracy"].as_f64().unwrap();
        let mean = (random["sensitivity"].as_f64().unwrap()
            + random["specificity"].as_f64().unwrap())
            / 2.0;
        assert_eq!(ba, mean);
    }

    #[test]
    fn write_covers_the_indexed_datasets_writes_the_same_bytes_twice_and_leaves_reports_untouched()
    {
        let temp = tempfile::tempdir().unwrap();
        let dir = temp.path();
        std::fs::create_dir_all(dir.join("datasets")).unwrap();
        std::fs::create_dir_all(dir.join("baselines")).unwrap();
        let raw = bytes(&synthetic(42));
        std::fs::write(dir.join("datasets/synthetic__fixturemycin.json"), &raw).unwrap();
        std::fs::write(
            dir.join("index.json"),
            br#"{"datasets":[{"id":"synthetic__fixturemycin"}]}"#,
        )
        .unwrap();
        std::fs::write(dir.join("baselines/gone__drug.json"), b"{}").unwrap();
        let written = write(dir).unwrap();
        assert_eq!(
            written,
            Written {
                files: 1,
                mic: 1,
                confirmed: 1,
                copies: 0,
                copy_datasets: 0,
                refused: vec![]
            }
        );
        assert_eq!(
            std::fs::read(dir.join("datasets/synthetic__fixturemycin.json")).unwrap(),
            raw
        );
        assert!(!dir.join("baselines/gone__drug.json").exists());
        let first = std::fs::read(dir.join("baselines/synthetic__fixturemycin.json")).unwrap();
        assert_eq!(
            serde_json::from_slice::<Value>(&first).unwrap(),
            for_report(&raw).unwrap()
        );
        write(dir).unwrap();
        assert_eq!(
            std::fs::read(dir.join("baselines/synthetic__fixturemycin.json")).unwrap(),
            first
        );
    }

    #[test]
    fn a_refused_report_leaves_the_folder_as_it_was() {
        let temp = tempfile::tempdir().unwrap();
        let dir = temp.path();
        std::fs::create_dir_all(dir.join("datasets")).unwrap();
        std::fs::create_dir_all(dir.join("baselines")).unwrap();
        let mut first = synthetic(42);
        first["id"] = json!("a__drug");
        std::fs::write(dir.join("datasets/a__drug.json"), bytes(&first)).unwrap();
        let mut old = synthetic(42);
        old["id"] = json!("b__drug");
        old["evaluation"]["revision"] = json!(EVALUATION_REVISION - 1);
        std::fs::write(dir.join("datasets/b__drug.json"), bytes(&old)).unwrap();
        std::fs::write(
            dir.join("index.json"),
            br#"{"datasets":[{"id":"a__drug"},{"id":"b__drug"}]}"#,
        )
        .unwrap();
        std::fs::write(dir.join("baselines/kept__drug.json"), b"{}").unwrap();
        let err = write(dir).unwrap_err().to_string();
        assert!(err.contains("b__drug"), "{err}");
        // The first report passed, and still nothing was written or removed.
        let names: Vec<_> = std::fs::read_dir(dir.join("baselines"))
            .unwrap()
            .map(|e| e.unwrap().file_name().into_string().unwrap())
            .collect();
        assert_eq!(names, ["kept__drug.json"]);
        assert!(!dir.join(SUMMARY).exists());
    }

    // ------------------------------------------- the rule, the files, the folder

    /// A report tree in a temporary folder: `datasets/<id>.json` for each
    /// report and an `index.json` listing them in order.
    fn tree(reports: &[(&str, &Value)]) -> tempfile::TempDir {
        let temp = tempfile::tempdir().unwrap();
        let dir = temp.path();
        std::fs::create_dir_all(dir.join("datasets")).unwrap();
        let mut rows = Vec::new();
        for (id, report) in reports {
            let mut report = (*report).clone();
            report["id"] = json!(id);
            std::fs::write(dir.join(format!("datasets/{id}.json")), bytes(&report)).unwrap();
            rows.push(json!({ "id": id }));
        }
        std::fs::write(dir.join("index.json"), bytes(&json!({ "datasets": rows }))).unwrap();
        temp
    }

    fn names_in(dir: &Path) -> Vec<String> {
        let mut names: Vec<String> = std::fs::read_dir(dir)
            .unwrap()
            .map(|e| e.unwrap().file_name().into_string().unwrap())
            .collect();
        names.sort();
        names
    }

    #[test]
    fn one_spanning_group_beside_three_inside_groups_is_written_unconfirmed() {
        // Synthetic, three folds: three groups inside one fold each (0.25 in
        // fold 0, 0.5 in fold 1, 0.75 in fold 2) and one value, 1.5, shared
        // by fold 0 and fold 1. X = 1 <= I = 3.
        let folds = [0, 0, 0, 1, 1, 1, 2, 2, 2];
        let predicted = [0.25, 0.25, 1.5, 0.5, 0.5, 1.5, 0.75, 0.75, 2.0];
        let check = prediction_check(&predicted, &folds, 3);
        assert_eq!((check.groups, check.spanning), (3, 1));
        assert_eq!(check.verdict(), Verdict::Unconfirmed);
        // Two spanning groups beside one inside: X > I, refused.
        let predicted = [0.25, 0.25, 1.5, 0.5, 2.0, 1.5, 0.75, 2.0, 3.0];
        let check = prediction_check(&predicted, &folds, 3);
        assert_eq!((check.groups, check.spanning), (1, 2));
        assert_eq!(check.verdict(), Verdict::Refused);
        // X = I is not more: written unconfirmed.
        let predicted = [0.25, 0.25, 1.5, 0.5, 0.6, 1.5, 0.75, 0.8, 3.0];
        let check = prediction_check(&predicted, &folds, 3);
        assert_eq!((check.groups, check.spanning), (1, 1));
        assert_eq!(check.verdict(), Verdict::Unconfirmed);
        // On a whole synthetic report: one isolate of fold 1 given the
        // prediction of an isolate of fold 0 makes one spanning group beside
        // the many inside groups, and the file is written with X = 1 and
        // confirmed false.
        let mut report = synthetic(42);
        let f = folds_at(&report, 42);
        let a = (0..N_MIC).find(|&i| f[i] == 0).unwrap();
        let b = (0..N_MIC).find(|&i| f[i] == 1 && i % 8 != a % 8).unwrap();
        report["mic"]["predicted"][b] = report["mic"]["predicted"][a].clone();
        stamp(&mut report);
        let file = for_report(&bytes(&report)).unwrap();
        let check = &file["mic"]["prediction_check"];
        assert_eq!(check["spanning"], 1);
        assert!(check["groups"].as_u64().unwrap() >= 3, "{check}");
        assert_eq!(check["confirmed"], false);
    }

    /// A synthetic dataset run through the repository's own fit
    /// (`learn::eval::evaluate`, default settings, as the review's helper
    /// runs it), written as the report fields this command reads. 200
    /// made-up isolates, not BV-BRC data: 140 susceptible at exact log2 MIC
    /// 0 and 60 resistant at exact log2 MIC 2, lineage `i % 20`, and one
    /// feature present exactly in the resistant ones; with `noise`, five
    /// more features spread over both classes.
    fn fitted(noise: bool) -> Value {
        use crate::bits::BitMatrix;
        use crate::learn::eval::{EvalConfig, MicInput, evaluate};
        let n = 200;
        let resistant = |i: usize| i % 10 < 3;
        let y = Bits::from_indices(n, (0..n).filter(|&i| resistant(i)));
        let mut x = BitMatrix::with_cols(n);
        x.push_row(y.words());
        if noise {
            for (a, b) in [(7, 11), (5, 13), (3, 17), (11, 19), (13, 23)] {
                let row = Bits::from_indices(n, (0..n).filter(|&i| (i * a) % b < b / 3));
                x.push_row(row.words());
            }
        }
        let mics: Vec<Option<Mic>> = (0..n)
            .map(|i| exact(if resistant(i) { 2 } else { 0 }))
            .collect();
        let lineage: Vec<u32> = (0..n).map(|i| (i % 20) as u32).collect();
        let cfg = EvalConfig::default();
        let input = MicInput {
            mic: &mics,
            breakpoint: None,
        };
        let out = evaluate(
            "synthetic__fitmycin",
            &x,
            &y,
            &lineage,
            &[],
            3,
            Some(input),
            &cfg,
        )
        .unwrap();
        let samples: Vec<Value> = (0..n)
            .map(|i| json!({ "id": format!("s{i}"), "resistant": resistant(i), "mic": mics[i], "lineage": lineage[i] }))
            .collect();
        let splits: Vec<Value> = out
            .splits
            .iter()
            .map(|s| json!({ "scheme": s.scheme, "folds": s.folds, "test_fold": s.test_fold }))
            .collect();
        json!({
            "id": "synthetic__fitmycin",
            "evaluation": { "revision": EVALUATION_REVISION },
            "eval_config": cfg,
            "samples": samples,
            "splits": splits,
            "mic": out.mic.expect("an MIC model"),
        })
    }

    #[test]
    fn the_repository_s_own_fit_can_give_one_value_to_two_folds_and_only_that_dataset_is_refused() {
        // With one feature, every fold of the repository's own fit predicts
        // its start value for the susceptible isolates and one shared value
        // for the resistant ones: correct folds, I = 0 and X = 2.
        let one = fitted(false);
        let err = for_report(&bytes(&one)).unwrap_err();
        let refused = err.downcast_ref::<Refused>().expect("a refusal");
        assert_eq!((refused.inside, refused.spanning), (0, 2), "{refused}");
        // With noise features the fit gives each fold its own values, and the
        // same check confirms the same folds.
        let noise = fitted(true);
        let file = for_report(&bytes(&noise)).unwrap();
        let check = &file["mic"]["prediction_check"];
        assert_eq!(
            (check["confirmed"].as_bool(), check["spanning"].as_u64()),
            (Some(true), Some(0))
        );
        assert!(check["groups"].as_u64().unwrap() >= 3, "{check}");
        // In one tree, the refused dataset gets no file and loses the one it
        // had, the other dataset is written, and the run reports the refusal.
        let temp = tree(&[("noise__fitmycin", &noise), ("one__fitmycin", &one)]);
        let dir = temp.path();
        std::fs::create_dir_all(dir.join("baselines")).unwrap();
        std::fs::write(dir.join("baselines/one__fitmycin.json"), b"{}").unwrap();
        let written = write(dir).unwrap();
        assert_eq!((written.files, written.confirmed), (1, 1));
        assert_eq!(written.refused.len(), 1);
        assert_eq!(written.refused[0].dataset, "one__fitmycin");
        assert_eq!(names_in(&dir.join("baselines")), ["noise__fitmycin.json"]);
        let summary: Value =
            serde_json::from_slice(&std::fs::read(dir.join(SUMMARY)).unwrap()).unwrap();
        let ids: Vec<&str> = summary["datasets"]
            .as_array()
            .unwrap()
            .iter()
            .map(|d| d["dataset"].as_str().unwrap())
            .collect();
        assert_eq!(ids, ["noise__fitmycin"]);
    }

    #[test]
    fn constants_and_the_model_score_exact_test_mics_only_and_round_half_away_from_zero() {
        // Synthetic, two folds. Fold 0 holds exact 2, 4, 7 and the censored
        // values "at most 3" and "above 4"; fold 1 holds exact 3, 3, 4 and
        // three values "at most -5".
        let at_most = |log2| {
            Some(Mic {
                log2,
                censor: Censor::AtMost,
            })
        };
        let mics = vec![
            exact(2),
            exact(4),
            exact(7),
            at_most(3),
            Some(Mic {
                log2: 4,
                censor: Censor::Above,
            }),
            exact(3),
            exact(3),
            exact(4),
            at_most(-5),
            at_most(-5),
            at_most(-5),
        ];
        let folds = [0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1];
        // 0.5 and 2.5 round away from zero to 1 and 3, within one of 2 and
        // of 4; rounded half to even they give 0 and 2, two steps away.
        let predicted = [0.5, 2.5, 7.0, 3.0, 4.0, 3.0, 3.0, 4.0, -5.0, -5.0, -5.0];
        let got = mic_constants(&mics, &predicted, &folds, 2);
        // Fold 0 trains on fold 1's exact 3, 3, 4 (the censored -5s left
        // out, which would move the median to -5): median 3, mode 3. Of its
        // exact test MICs 2, 4 and 7, two are within one of 3; the censored
        // "at most 3" and "above 4" would agree too, and are not scored.
        assert_eq!(
            got[0],
            FoldConstants {
                train_exact: 3,
                test_exact: 3,
                median: Some(3),
                upper_median_differs: false,
                mode: Some(3),
                mode_tied: false,
                agree_model: 3,
                agree_median: 2,
                agree_mode: 2,
            }
        );
        // Fold 1 trains on fold 0's exact 2, 4, 7 (the censored 3 and 4 left
        // out, which would make 4 the mode): median 4, mode 2 by the tie to
        // the lowest. Its exact test MICs 3, 3, 4 are all within one of 4,
        // and two of them of 2.
        assert_eq!(
            got[1],
            FoldConstants {
                train_exact: 3,
                test_exact: 3,
                median: Some(4),
                upper_median_differs: false,
                mode: Some(2),
                mode_tied: true,
                agree_model: 3,
                agree_median: 3,
                agree_mode: 2,
            }
        );
    }

    #[test]
    fn a_lineage_share_of_exactly_one_half_is_a_resistant_call() {
        // Synthetic, two outer folds. Genome 0 (fold 0, lineage 1, R) shares
        // its lineage with genomes 2 (R) and 3 (S) of fold 1: share 1 in 2,
        // called resistant. Genome 1 (fold 0, lineage 2, S) is unseen and
        // takes fold 1's prevalence, also 1 in 2. Fold 1's genomes share
        // lineage 1 with genome 0 alone: share 1.
        let sample = |resistant, lineage| Sample {
            resistant,
            mic: None,
            lineage,
        };
        let samples = [
            sample(true, 1),
            sample(false, 2),
            sample(true, 1),
            sample(false, 1),
        ];
        let splits = [Split {
            scheme: "random".into(),
            folds: 2,
            test_fold: vec![0, 0, 1, 1],
        }];
        let part = lineage_part("synthetic__halfmycin", &samples, &splits).unwrap();
        let random = &part[0];
        assert_eq!(
            (
                random["tp"].as_u64(),
                random["fp"].as_u64(),
                random["tn"].as_u64(),
                random["fn"].as_u64()
            ),
            (Some(2), Some(2), Some(0), Some(0))
        );
        assert_eq!(random["sensitivity"], 1.0);
        assert_eq!(random["specificity"], 0.0);
        assert_eq!(random["balanced_accuracy"], 0.5);
        assert_eq!(random["unseen"], 1);
    }

    #[test]
    fn write_waits_while_another_writer_holds_the_report_directory_s_lock() {
        let report = synthetic(42);
        let temp = tree(&[("a__drug", &report)]);
        let dir = temp.path().to_path_buf();
        // The lock `grm report` and `grm intervals` take on the directory.
        let held = std::fs::File::open(&dir).unwrap();
        held.lock().unwrap();
        let (done, finished) = std::sync::mpsc::channel();
        let path = dir.clone();
        let writer = std::thread::spawn(move || {
            let written = write(&path).map(|w| w.files);
            done.send(()).unwrap();
            written
        });
        // While the lock is held, the second writer waits and writes nothing.
        assert!(
            finished
                .recv_timeout(std::time::Duration::from_millis(700))
                .is_err()
        );
        assert!(!dir.join("baselines").exists() && !dir.join(SUMMARY).exists());
        held.unlock().unwrap();
        finished
            .recv_timeout(std::time::Duration::from_secs(60))
            .unwrap();
        assert_eq!(writer.join().unwrap().unwrap(), 1);
        assert_eq!(names_in(&dir.join("baselines")), ["a__drug.json"]);
    }

    #[cfg(unix)]
    #[test]
    fn a_symlinked_baselines_folder_is_refused_before_anything_is_read_or_written() {
        let report = synthetic(42);
        let temp = tree(&[("a__drug", &report)]);
        let dir = temp.path();
        let outside = tempfile::tempdir().unwrap();
        std::fs::write(outside.path().join("notes.txt"), b"kept").unwrap();
        std::fs::write(outside.path().join("data.json"), b"kept").unwrap();
        std::os::unix::fs::symlink(outside.path(), dir.join("baselines")).unwrap();
        let err = write(dir).unwrap_err().to_string();
        assert!(
            err.contains("baselines") && err.contains("symbolic link"),
            "{err}"
        );
        assert_eq!(names_in(outside.path()), ["data.json", "notes.txt"]);
        assert_eq!(
            std::fs::read(outside.path().join("data.json")).unwrap(),
            b"kept"
        );
        assert!(!dir.join(SUMMARY).exists());
        // A plain file under that name is refused the same way.
        std::fs::remove_file(dir.join("baselines")).unwrap();
        std::fs::write(dir.join("baselines"), b"").unwrap();
        assert!(
            write(dir)
                .unwrap_err()
                .to_string()
                .contains("not a directory")
        );
    }

    #[cfg(unix)]
    #[test]
    fn removal_takes_only_regular_json_and_tmp_files() {
        let report = synthetic(42);
        let temp = tree(&[("a__drug", &report)]);
        let dir = temp.path();
        let out = dir.join("baselines");
        std::fs::create_dir_all(out.join("old__drug.json")).unwrap();
        std::fs::write(out.join("gone__drug.json"), b"{}").unwrap();
        std::fs::write(out.join("gone__drug.json.4242.tmp"), b"{").unwrap();
        std::fs::write(out.join("notes.txt"), b"kept").unwrap();
        let outside = tempfile::tempdir().unwrap();
        std::fs::write(outside.path().join("x.json"), b"kept").unwrap();
        std::os::unix::fs::symlink(outside.path().join("x.json"), out.join("link__drug.json"))
            .unwrap();
        write(dir).unwrap();
        // The stale file and the pending write go; the directory, the link
        // and the other file stay, and nothing outside is touched.
        assert_eq!(
            names_in(&out),
            [
                "a__drug.json",
                "link__drug.json",
                "notes.txt",
                "old__drug.json"
            ]
        );
        assert_eq!(
            std::fs::read(outside.path().join("x.json")).unwrap(),
            b"kept"
        );
    }

    #[test]
    fn an_index_id_with_a_path_separator_or_dot_dot_is_refused_before_any_write() {
        for bad in ["../escape", "a/b", "a\\b", "..", ""] {
            let temp = tree(&[]);
            let dir = temp.path();
            std::fs::write(
                dir.join("index.json"),
                bytes(&json!({ "datasets": [{ "id": "a__drug" }, { "id": bad }] })),
            )
            .unwrap();
            let mut report = synthetic(42);
            report["id"] = json!("a__drug");
            std::fs::write(dir.join("datasets/a__drug.json"), bytes(&report)).unwrap();
            let err = write(dir).unwrap_err().to_string();
            assert!(err.contains(&format!("{bad:?}")), "{err}");
            assert_eq!(names_in(dir), ["datasets", "index.json"], "{bad}");
        }
    }

    #[test]
    fn a_failed_write_names_its_path_and_leaves_no_pending_file() {
        let a = synthetic(42);
        let temp = tree(&[("a__drug", &a), ("b__drug", &a)]);
        let dir = temp.path();
        // A non-empty directory where b's file goes: the rename fails.
        std::fs::create_dir_all(dir.join("baselines/b__drug.json/inside")).unwrap();
        let err = format!("{:#}", write(dir).unwrap_err());
        assert!(err.contains("b__drug.json"), "{err}");
        assert!(
            names_in(&dir.join("baselines"))
                .iter()
                .all(|n| !n.ends_with(".tmp")),
            "{:?}",
            names_in(&dir.join("baselines"))
        );
    }

    #[test]
    fn the_summary_pins_the_index_bytes_and_each_report_and_holds_the_median_constant() {
        let mic = synthetic(42);
        let mut plain = synthetic(42);
        plain["mic"] = Value::Null;
        let temp = tree(&[("a__drug", &mic), ("b__drug", &plain)]);
        let dir = temp.path();
        write(dir).unwrap();
        let index = std::fs::read(dir.join("index.json")).unwrap();
        let raw = std::fs::read(dir.join(SUMMARY)).unwrap();
        let summary: Value = serde_json::from_slice(&raw).unwrap();
        let file: Value =
            serde_json::from_slice(&std::fs::read(dir.join("baselines/a__drug.json")).unwrap())
                .unwrap();
        let m = &file["mic"];
        // The fixture's median and mode differ, so a constant read from the
        // mode would show.
        assert_ne!(
            m["median"]["essential_agreement"],
            m["mode"]["essential_agreement"]
        );
        assert_eq!(
            summary,
            json!({
                "schema": SUMMARY_SCHEMA,
                "evaluation_revision": EVALUATION_REVISION,
                "index_sha256": crate::research::sha256(&index),
                "datasets": [{
                    "dataset": "a__drug",
                    "report_sha256": crate::research::sha256(&std::fs::read(dir.join("datasets/a__drug.json")).unwrap()),
                    "model": m["model"]["essential_agreement"],
                    "constant": m["median"]["essential_agreement"],
                    "confirmed": true,
                }],
            })
        );
        // Compact, keys sorted, no final newline; the same bytes again.
        assert_eq!(raw, serde_json::to_vec(&summary).unwrap());
        write(dir).unwrap();
        assert_eq!(std::fs::read(dir.join(SUMMARY)).unwrap(), raw);
        // The removal step never takes it: it lives beside index.json.
        assert!(dir.join(SUMMARY).exists());
    }
}
