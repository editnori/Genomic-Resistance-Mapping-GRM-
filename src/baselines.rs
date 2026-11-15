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
//!   they are not the random split's folds. Each rebuilt fold's training
//!   isolates give the threshold `infer_breakpoint` fits, and that must equal
//!   the threshold the report stores for the fold (`mic.fold_breakpoints`,
//!   every field, agreement bit for bit); one that differs stops the command.
//!   The folds are `confirmed` when at least one fold has a stored threshold
//!   and every fold matches; with no stored threshold to check they are
//!   still rebuilt the same way, and `confirmed` is false. Each fold's exact
//!   training MICs give two constants, their lower median and their mode
//!   (`lower_median`, `mode`), and each is scored on the same isolates as the
//!   model: the fold's test isolates with an exact MIC, agreeing when the
//!   constant is within one doubling dilution of the measured value. The
//!   model's own count is recounted from `mic.predicted` under the same rule
//!   (`evaluate_mic` rounds a prediction half away from zero) and must give
//!   the report's `mic.essential_agreement` exactly.
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
//! and a second run writes the same bytes. It refuses a report of another
//! evaluation revision than this binary's. Every file is computed before any
//! is written, so a refusal leaves the folder as it was.
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

/// The most frequent value, a tie going to the lowest, and whether another
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

/// The MIC part of a baselines file, or an error naming the dataset when the
/// rebuilt folds contradict a stored threshold or the saved predictions do
/// not give the report's own agreement.
fn mic_part(id: &str, samples: &[Sample], cfg: &FoldConfig, mic: &MicResult) -> Result<Value> {
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
    let y = Bits::from_indices(n, (0..n).filter(|&i| samples[i].resistant));
    let with_mic = Bits::from_indices(n, has.iter().copied());
    // The draw `evaluate_mic` makes.
    let folds = training_folds(&y, &with_mic, None, k, cfg.seed);
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
            bail!(
                "{id}: the MIC folds rebuilt at seed {} contradict the report's stored training-fold threshold of fold {f} (stored {stored:?}, rebuilt {rebuilt:?}); no baselines file was written",
                cfg.seed
            );
        }
        checked += stored.is_some() as usize;
    }
    let per_fold = mic_constants(&mics, &mic.predicted, &folds, k);
    let total = |pick: fn(&FoldConstants) -> usize| per_fold.iter().map(pick).sum::<usize>();
    let model = total(|f| f.agree_model);
    let model_share = model as f64 / n_exact.max(1) as f64;
    ensure!(
        model_share.to_bits() == mic.essential_agreement.to_bits(),
        "{id}: the saved MIC predictions put {model} of {n_exact} exact MICs within one dilution ({model_share}), and the report says {}",
        mic.essential_agreement
    );
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
        "confirmed": checked > 0,
        "thresholds_checked": checked,
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
            mic_part(id, &report.samples, cfg, mic)?
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

/// What one run wrote.
#[derive(Debug, PartialEq, Eq)]
pub struct Written {
    pub files: usize,
    /// Files with an MIC part.
    pub mic: usize,
    /// Of those, the ones whose rebuilt folds a stored threshold confirms.
    pub confirmed: usize,
}

/// Writes `baselines/<id>.json` under `report_dir` for every dataset that
/// its `index.json` lists and removes any other file there. Every file is
/// computed first, so a refused report leaves the folder as it was. It holds
/// the report directory's lock, as `grm report` does, so it never reads a
/// build half done.
pub fn write(report_dir: &Path) -> Result<Written> {
    let lock = std::fs::File::open(report_dir)
        .with_context(|| format!("open {}", report_dir.display()))?;
    lock.lock()
        .with_context(|| format!("lock {}", report_dir.display()))?;
    let index: Value = serde_json::from_slice(
        &std::fs::read(report_dir.join("index.json")).context("run `grm report` first")?,
    )?;
    let mut files = Vec::new();
    for row in index["datasets"].as_array().into_iter().flatten() {
        let id = row["id"]
            .as_str()
            .context("index.json lists a dataset with no id")?;
        let path = report_dir.join("datasets").join(format!("{id}.json"));
        let bytes = std::fs::read(&path).with_context(|| format!("read {}", path.display()))?;
        let file = for_report(&bytes)?;
        if file["dataset"] != id {
            bail!("{} records the id {}", path.display(), file["dataset"]);
        }
        files.push((format!("{id}.json"), file));
    }
    let out = report_dir.join("baselines");
    std::fs::create_dir_all(&out)?;
    let mut kept = BTreeSet::new();
    let mut written = Written {
        files: 0,
        mic: 0,
        confirmed: 0,
    };
    for (name, file) in files {
        crate::report::publish(&out.join(&name), &serde_json::to_vec(&file)?)?;
        written.files += 1;
        if !file["mic"].is_null() {
            written.mic += 1;
            written.confirmed += (file["mic"]["confirmed"] == true) as usize;
        }
        kept.insert(name);
    }
    for entry in std::fs::read_dir(&out)? {
        let entry = entry?;
        if !entry.file_name().to_str().is_some_and(|n| kept.contains(n)) {
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

    /// A synthetic report: `N_MIC` isolates with an MIC (exact, except every
    /// 7th with a low MIC, which is a bound "at most") and 20 without. An MIC above 2 is resistant,
    /// except 4 isolates whose label disagrees, so each training fold's
    /// threshold has an agreement that depends on which isolates it holds.
    /// The model predicts each MIC exactly, except every 5th isolate, which
    /// it misses by two dilutions. Lineages run in blocks of 10 isolates;
    /// the random split deals isolates round robin into 5 folds and the
    /// lineage split holds out whole blocks.
    fn synthetic(seed: u64) -> Value {
        let n = N_MIC + 20;
        let mut samples = Vec::new();
        let mut predicted = Vec::new();
        for i in 0..n {
            let (mic, pred, resistant) = if i < N_MIC {
                let log2 = (i % 8) as i32 - 1;
                let censor = if i % 7 == 6 && log2 <= 1 {
                    "AtMost"
                } else {
                    "Exact"
                };
                let flip = [3, 50, 101, 177].contains(&i);
                let pred = if i % 5 == 0 { log2 + 2 } else { log2 };
                (
                    json!({ "log2": log2, "censor": censor }),
                    json!(pred as f64),
                    (log2 > 2) != flip,
                )
            } else {
                (Value::Null, Value::Null, i % 2 == 0)
            };
            samples.push(json!({ "id": format!("g{i}"), "resistant": resistant, "mic": mic, "lineage": i / 10 }));
            predicted.push(pred);
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
                "confusion": [], "predicted": predicted,
            },
        });
        // Store what the fit would have stored: thresholds from the folds
        // drawn at the report's seed, and the model's agreement.
        let thresholds = stored_thresholds(&report, seed);
        let mics: Vec<Option<Mic>> = parse_mics(&report);
        let n_exact = mics
            .iter()
            .flatten()
            .filter(|m| m.censor == Censor::Exact)
            .count();
        let agree = mics
            .iter()
            .enumerate()
            .filter(|(_, m)| m.is_some_and(|m| m.censor == Censor::Exact))
            .filter(|&(i, m)| {
                within_one(
                    report["mic"]["predicted"][i].as_f64().unwrap().round() as i32,
                    m.unwrap().log2,
                )
            })
            .count();
        report["mic"]["n_exact"] = json!(n_exact);
        report["mic"]["essential_agreement"] = json!(agree as f64 / n_exact as f64);
        report["mic"]["fold_breakpoints"] = json!(thresholds);
        report
    }

    fn parse_mics(report: &Value) -> Vec<Option<Mic>> {
        report["samples"]
            .as_array()
            .unwrap()
            .iter()
            .map(|s| serde_json::from_value(s["mic"].clone()).unwrap())
            .collect()
    }

    fn stored_thresholds(report: &Value, seed: u64) -> Vec<Option<Breakpoint>> {
        let mics = parse_mics(report);
        let n = mics.len();
        let samples = report["samples"].as_array().unwrap();
        let y = Bits::from_indices(n, (0..n).filter(|&i| samples[i]["resistant"] == true));
        let with_mic = Bits::from_indices(n, (0..n).filter(|&i| mics[i].is_some()));
        let folds = training_folds(&y, &with_mic, None, 5, seed);
        fold_masks(&folds, 5)
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
        // Every fold's training isolates support a threshold, so all five
        // stored thresholds are checked.
        assert_eq!(mic["confirmed"], true);
        assert_eq!(mic["thresholds_checked"], 5);
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
    fn the_folds_are_rebuilt_at_the_report_s_seed_and_a_contradiction_stops_the_command() {
        // The thresholds were stored from the folds of seed 42. Read as a
        // report of seed 43, the rebuilt folds contradict them.
        let mut report = synthetic(42);
        report["eval_config"]["seed"] = json!(43);
        let err = for_report(&bytes(&report)).unwrap_err().to_string();
        assert!(err.contains("synthetic__fixturemycin"), "{err}");
        assert!(err.contains("contradict"), "{err}");
        // A stored threshold one dilution off.
        let mut report = synthetic(42);
        let r = report["mic"]["fold_breakpoints"][2]["r_above"]
            .as_i64()
            .unwrap();
        report["mic"]["fold_breakpoints"][2]["r_above"] = json!(r + 1);
        assert!(
            for_report(&bytes(&report))
                .unwrap_err()
                .to_string()
                .contains("fold 2")
        );
        // A threshold stored as none where the rebuilt fold fits one.
        let mut report = synthetic(42);
        report["mic"]["fold_breakpoints"][0] = Value::Null;
        assert!(for_report(&bytes(&report)).is_err());
    }

    #[test]
    fn folds_with_no_stored_threshold_are_written_unconfirmed() {
        // Synthetic: every label resistant but 9, so no training fold has the
        // 10 susceptible isolates a threshold needs and every stored
        // threshold is none, which the rebuilt folds reproduce.
        let mut report = synthetic(42);
        for (i, s) in report["samples"]
            .as_array_mut()
            .unwrap()
            .iter_mut()
            .enumerate()
        {
            s["resistant"] = json!(i >= 9);
        }
        report["mic"]["fold_breakpoints"] = json!(stored_thresholds(&report, 42));
        assert!(
            report["mic"]["fold_breakpoints"]
                .as_array()
                .unwrap()
                .iter()
                .all(Value::is_null)
        );
        let file = for_report(&bytes(&report)).unwrap();
        assert_eq!(file["mic"]["confirmed"], false);
        assert_eq!(file["mic"]["thresholds_checked"], 0);
        assert!(file["mic"]["median"]["agree"].is_u64());
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
                confirmed: 1
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
    }
}
