//! Experiment, off by default (plan block C9, docs/served-pick.md): pick the
//! served learner inside each outer training fold, by inner folds grouped
//! by lineage over the five learners, and write what each pick would serve
//! to a file of its own.
//!
//! Nothing here changes a report. `grm learn --exp-pick DIR` turns it on for
//! that process; without the flag `enabled` is false and evaluation gives
//! the same outputs as before. With it, `evaluate` and `evaluate_assigned`
//! fit, score and report exactly as without it, and then write one result
//! per dataset, split and seed, `DIR/<dataset>.<scheme>-s<seed>.pick.json`.
//!
//! In outer fold f, with training genomes T_f, the inner partition is
//! `training_folds(y, T_f, lineage groups, K, fold_seed(seed, f))`: at K = 3
//! the partition the release's SCM and CART tune on in that fold. Inner fit
//! j trains on T_f less inner fold j and is scored on inner fold j. It is
//! fitted as the release fits an outer fold of the same split
//! (`eval::fit_seeded`): candidates from its own training genomes at the
//! dataset's support floor, the split's tuning groups, and every draw seeded
//! with `inner_seed`. No genome outside T_f reaches it. Each learner's calls
//! on the K inner test folds are summed, and the fold's pick is the learner
//! with the highest balanced accuracy of that sum (`report::inner_pick`).
//!
//! A result holds counts, rates, learner names, settings and reasons only:
//! no genome id, no feature and no row per genome.

use super::eval::{
    EvalConfig, LEARNERS, LassoContext, LearnerResult, SplitResult, check_tuning_seeded, fit_seeded,
};
use super::{any_training_candidate, fold_masks, train_prevalence, training_candidates};
use super::{training_folds, validate_folds};
use crate::bits::{BitMatrix, Bits};
use crate::cutoff::{Calls, calls};
use crate::report::{PickTie, inner_pick, release_pick};
use anyhow::{Context, Result, ensure};
use serde::Serialize;
use std::path::{Path, PathBuf};
use std::sync::OnceLock;
use std::time::Instant;

/// The result format. A change to what a result holds changes it.
pub const FORMAT: &str = "grm-c9-pick-1";

/// The experiment's settings, as `enable` set them.
#[derive(Clone, Debug)]
pub struct Settings {
    pub dir: PathBuf,
    /// Inner folds per outer fold; 0 writes the default's pick and the
    /// outer counts and fits no inner fold (the twin's runs).
    pub inner_folds: usize,
    pub tie: PickTie,
}

static SINK: OnceLock<Settings> = OnceLock::new();

/// Tests only: a capture's settings and its results, by file name.
#[cfg(test)]
type Captured = Option<(Settings, Vec<(String, String)>)>;

#[cfg(test)]
thread_local! {
    /// Tests only: while `capture` runs on a thread, that thread's results
    /// go here, by file name, and its settings come from here.
    static CAPTURED: std::cell::RefCell<Captured> = const { std::cell::RefCell::new(None) };
}

/// Turns the experiment on for this process. `dir` is created if missing.
/// A second call fails.
pub fn enable(dir: &Path, inner_folds: usize, tie: PickTie) -> Result<()> {
    ensure!(
        inner_folds == 0 || inner_folds >= 2,
        "--exp-pick-inner takes 0 (no inner fit) or 2 or more folds, not {inner_folds}"
    );
    std::fs::create_dir_all(dir).with_context(|| format!("create {}", dir.display()))?;
    SINK.set(Settings {
        dir: dir.to_path_buf(),
        inner_folds,
        tie,
    })
    .map_err(|_| anyhow::anyhow!("the C9 pick folder is already set"))
}

/// Whether `enable` has run in this process (or, in a unit test, whether
/// `capture` runs on this thread).
pub fn enabled() -> bool {
    SINK.get().is_some() || capturing()
}

#[cfg(not(test))]
fn capturing() -> bool {
    false
}

#[cfg(test)]
fn capturing() -> bool {
    CAPTURED.with(|c| c.borrow().is_some())
}

#[cfg(not(test))]
fn settings() -> Settings {
    SINK.get().expect("enabled").clone()
}

#[cfg(test)]
fn settings() -> Settings {
    CAPTURED
        .with(|c| c.borrow().as_ref().map(|(s, _)| s.clone()))
        .or_else(|| SINK.get().cloned())
        .expect("enabled")
}

/// Tests only: runs `f` with the experiment on for this thread alone, with
/// `inner_folds` and `tie`, and returns `f`'s result with each result file
/// it would have written, as (file name, JSON text).
#[cfg(test)]
pub(crate) fn capture<R>(
    inner_folds: usize,
    tie: PickTie,
    f: impl FnOnce() -> R,
) -> (R, Vec<(String, String)>) {
    let s = Settings {
        dir: PathBuf::new(),
        inner_folds,
        tie,
    };
    CAPTURED.with(|c| *c.borrow_mut() = Some((s, Vec::new())));
    let result = f();
    let files = CAPTURED
        .with(|c| c.borrow_mut().take())
        .map(|(_, v)| v)
        .unwrap_or_default();
    (result, files)
}

fn write(name: &str, text: String) -> Result<()> {
    #[cfg(test)]
    {
        let taken = CAPTURED.with(|c| {
            c.borrow_mut()
                .as_mut()
                .map(|(_, v)| v.push((name.to_string(), text.clone())))
                .is_some()
        });
        if taken {
            return Ok(());
        }
    }
    let dir = settings().dir;
    let path = dir.join(name);
    let tmp = dir.join(format!("{name}.tmp"));
    std::fs::write(&tmp, text).with_context(|| format!("write {}", tmp.display()))?;
    std::fs::rename(&tmp, &path).with_context(|| format!("replace {}", path.display()))
}

/// The seed of inner fit `j` of outer fold `fold` at evaluation seed
/// `seed`: 64-bit FNV-1a of the text `c9-inner/<seed>/<fold>/<j>`. It seeds
/// the fit's own inner folds, tuning holdouts and draws.
pub fn inner_seed(seed: u64, fold: usize, j: usize) -> u64 {
    let mut h = crate::fnv::Fnv::new();
    h.bytes(format!("c9-inner/{seed}/{fold}/{j}").as_bytes());
    h.0
}

/// The default's pick, as the release makes it from the run's splits.
#[derive(Clone, Debug, Serialize)]
pub struct DefaultPick {
    pub learner: String,
    pub random_balanced_accuracy: f64,
    /// Every learner tied with it on the random split, in pipeline order.
    pub tied: Vec<String>,
}

/// One learner's whole split as its report row has it.
#[derive(Clone, Debug, Serialize)]
pub struct SplitLearner {
    pub learner: String,
    pub balanced_accuracy: Option<f64>,
}

/// One learner's calls on one outer test fold, from the release's own
/// outer fit (a null fit's genomes at its training prevalence).
#[derive(Clone, Debug, Serialize)]
pub struct OuterCalls {
    pub learner: String,
    #[serde(flatten)]
    pub calls: Calls,
}

/// One learner over one outer fold's inner folds.
#[derive(Clone, Debug, Serialize)]
pub struct InnerLearner {
    pub learner: String,
    /// Its calls on each inner test fold, in inner fold order; empty when a
    /// fit failed.
    pub folds: Vec<Calls>,
    pub tp: usize,
    pub fp: usize,
    pub tn: usize,
    #[serde(rename = "fn")]
    pub fn_: usize,
    /// Balanced accuracy of the summed calls; None when a fit failed.
    pub balanced_accuracy: Option<f64>,
    /// How many of its inner fits were null (scored at their prevalence).
    pub null_fits: usize,
    pub error: Option<String>,
}

/// The inner pick of one outer fold.
#[derive(Clone, Debug, Serialize)]
pub struct InnerPick {
    /// `picked`, `no pick` or `not run` (inner folds 0).
    pub status: String,
    pub reason: Option<String>,
    /// Per inner fold: its training and test genomes, and resistant ones.
    pub inner_folds: Vec<FoldCounts>,
    pub learners: Vec<InnerLearner>,
    pub tied: Vec<String>,
    pub pick: Option<String>,
}

#[derive(Clone, Copy, Debug, Serialize)]
pub struct FoldCounts {
    pub train: usize,
    pub train_resistant: usize,
    pub test: usize,
    pub test_resistant: usize,
}

/// One outer fold.
#[derive(Clone, Debug, Serialize)]
pub struct FoldRecord {
    pub fold: usize,
    #[serde(flatten)]
    pub counts: FoldCounts,
    pub outer: Vec<OuterCalls>,
    pub inner: InnerPick,
}

/// One result file.
#[derive(Clone, Debug, Serialize)]
pub struct PickResult {
    pub format: &'static str,
    pub evaluation_revision: u32,
    pub dataset: String,
    pub scheme: String,
    pub seed: u64,
    pub inner_folds: usize,
    pub tie: &'static str,
    pub min_support: u32,
    pub samples: usize,
    pub tested: usize,
    pub tested_resistant: usize,
    pub unavailable_reason: Option<String>,
    /// The default's pick; on the lineage split only, from the same run's
    /// random split.
    pub default: Option<DefaultPick>,
    pub learners: Vec<SplitLearner>,
    pub folds: Vec<FoldRecord>,
    pub seconds_inner: f64,
}

pub fn file_name(dataset: &str, scheme: &str, seed: u64) -> String {
    format!("{dataset}.{scheme}-s{seed}.pick.json")
}

/// `evaluate`'s hook: the lineage split of `splits`, with the default's
/// pick from the same splits.
#[allow(clippy::too_many_arguments)]
pub(crate) fn record_drawn(
    dataset: &str,
    x: &BitMatrix,
    y: &Bits,
    groups: &[u32],
    duplicates: &[Vec<usize>],
    min_support: u32,
    cfg: &EvalConfig,
    splits: &[SplitResult],
) -> Result<()> {
    let default =
        release_pick(&serde_json::to_value(splits)?).map(|(learner, ba, tied)| DefaultPick {
            learner,
            random_balanced_accuracy: ba,
            tied,
        });
    let lineage = splits
        .iter()
        .find(|s| s.scheme == "lineage")
        .context("evaluate wrote no lineage split")?;
    let tuning = super::eval::tuning_groups("lineage", groups, duplicates);
    let n = x.cols;
    let all: Vec<usize> = (0..n).collect();
    let masks = if lineage.unavailable_reason.is_none() {
        fold_masks(&lineage.test_fold, lineage.folds)
    } else {
        vec![]
    };
    let trains: Vec<Bits> = masks.into_iter().map(|(train, _)| train).collect();
    let result = build(
        dataset,
        "lineage",
        x,
        y,
        groups,
        tuning.as_deref(),
        &trains,
        &all,
        &lineage.test_fold,
        &lineage.learners,
        lineage.unavailable_reason.clone(),
        default,
        min_support,
        cfg,
    )?;
    write(
        &file_name(dataset, "lineage", cfg.seed),
        serde_json::to_string(&result)?,
    )
}

/// `evaluate_assigned`'s hook: a supplied split. `tested[p]` is the sample
/// that `test_fold[p]` and each learner's `scores[p]` belong to. The
/// default's pick comes from the same species' lineage run.
#[allow(clippy::too_many_arguments)]
pub(crate) fn record_assigned(
    dataset: &str,
    scheme: &str,
    x: &BitMatrix,
    y: &Bits,
    groups: &[u32],
    tuning: Option<&[u32]>,
    trains: &[Bits],
    tested: &[usize],
    test_fold: &[usize],
    learners: &[LearnerResult],
    min_support: u32,
    cfg: &EvalConfig,
) -> Result<()> {
    let result = build(
        dataset,
        scheme,
        x,
        y,
        groups,
        tuning,
        trains,
        tested,
        test_fold,
        learners,
        None,
        None,
        min_support,
        cfg,
    )?;
    write(
        &file_name(dataset, scheme, cfg.seed),
        serde_json::to_string(&result)?,
    )
}

fn counts(train: &Bits, test: &Bits, y: &Bits) -> FoldCounts {
    FoldCounts {
        train: train.count() as usize,
        train_resistant: train.and(y).count() as usize,
        test: test.count() as usize,
        test_resistant: test.and(y).count() as usize,
    }
}

/// The balanced accuracy of summed calls; None when a class is absent.
fn balanced(tp: usize, fp: usize, tn: usize, fn_: usize) -> Option<f64> {
    (tp + fn_ > 0 && tn + fp > 0)
        .then(|| (tp as f64 / (tp + fn_) as f64 + tn as f64 / (tn + fp) as f64) / 2.0)
}

#[allow(clippy::too_many_arguments)]
fn build(
    dataset: &str,
    scheme: &str,
    x: &BitMatrix,
    y: &Bits,
    lineage: &[u32],
    tuning: Option<&[u32]>,
    trains: &[Bits],
    tested: &[usize],
    test_fold: &[usize],
    learners: &[LearnerResult],
    unavailable_reason: Option<String>,
    default: Option<DefaultPick>,
    min_support: u32,
    cfg: &EvalConfig,
) -> Result<PickResult> {
    let started = Instant::now();
    let set = settings();
    let n = x.cols;
    ensure!(
        tested.len() == test_fold.len(),
        "one outer fold per tested sample"
    );
    let truth = |i: usize| y.contains(i);
    let mut result = PickResult {
        format: FORMAT,
        evaluation_revision: crate::analysis::EVALUATION_REVISION,
        dataset: dataset.to_string(),
        scheme: scheme.to_string(),
        seed: cfg.seed,
        inner_folds: set.inner_folds,
        tie: set.tie.name(),
        min_support,
        samples: n,
        tested: if unavailable_reason.is_some() {
            0
        } else {
            tested.len()
        },
        tested_resistant: if unavailable_reason.is_some() {
            0
        } else {
            tested.iter().filter(|&&i| truth(i)).count()
        },
        unavailable_reason,
        default,
        learners: LEARNERS
            .iter()
            .map(|&l| SplitLearner {
                learner: l.to_string(),
                balanced_accuracy: learners
                    .iter()
                    .find(|r| r.learner == l)
                    .and_then(|r| r.metrics.as_ref())
                    .map(|m| m.balanced_accuracy),
            })
            .collect(),
        folds: vec![],
        seconds_inner: 0.0,
    };
    if result.unavailable_reason.is_some() {
        return Ok(result);
    }
    ensure!(
        learners.len() >= LEARNERS.len()
            && LEARNERS.iter().zip(learners).all(|(l, r)| r.learner == *l),
        "a split's learners are not the five in pipeline order"
    );
    for (f, train) in trains.iter().enumerate() {
        // The outer test fold: the tested samples of fold f.
        let positions: Vec<usize> = (0..tested.len()).filter(|&p| test_fold[p] == f).collect();
        let test = Bits::from_indices(n, positions.iter().map(|&p| tested[p]));
        ensure!(
            test.and(train).count() == 0,
            "outer fold {f} trains on a sample it tests"
        );
        let fold_truth: Vec<bool> = positions.iter().map(|&p| truth(tested[p])).collect();
        let outer = LEARNERS
            .iter()
            .zip(learners)
            .map(|(&l, r)| {
                let scores: Vec<f64> = positions.iter().map(|&p| r.scores[p]).collect();
                OuterCalls {
                    learner: l.to_string(),
                    calls: calls(&fold_truth, &scores),
                }
            })
            .collect();
        let inner = inner_fold_pick(
            dataset,
            scheme,
            x,
            y,
            lineage,
            tuning,
            train,
            f,
            min_support,
            cfg,
            &set,
        );
        result.folds.push(FoldRecord {
            fold: f,
            counts: counts(train, &test, y),
            outer,
            inner,
        });
    }
    result.seconds_inner = started.elapsed().as_secs_f64();
    Ok(result)
}

/// The inner folds of outer fold `f` at evaluation seed `seed`: its
/// training samples `train` in `k` folds grouped by `lineage`, drawn at the
/// fold's own seed (`ablate::plan::fold_seed`), the partition the release's
/// SCM and CART tune on there when `k` is 3. A sample outside `train` gets
/// no inner fold (`usize::MAX`).
pub(crate) fn inner_partition(
    y: &Bits,
    train: &Bits,
    lineage: &[u32],
    k: usize,
    seed: u64,
    f: usize,
) -> Vec<usize> {
    training_folds(
        y,
        train,
        Some(lineage),
        k,
        crate::ablate::plan::fold_seed(seed, f),
    )
}

/// The inner pick of outer fold `f`, whose training samples are `train`.
/// It reads no sample outside `train`.
#[allow(clippy::too_many_arguments)]
fn inner_fold_pick(
    dataset: &str,
    scheme: &str,
    x: &BitMatrix,
    y: &Bits,
    lineage: &[u32],
    tuning: Option<&[u32]>,
    train: &Bits,
    f: usize,
    min_support: u32,
    cfg: &EvalConfig,
    set: &Settings,
) -> InnerPick {
    let mut out = InnerPick {
        status: "not run".into(),
        reason: None,
        inner_folds: vec![],
        learners: vec![],
        tied: vec![],
        pick: None,
    };
    let k = set.inner_folds;
    if k == 0 {
        return out;
    }
    out.status = "no pick".into();
    let n = x.cols;
    let inner = inner_partition(y, train, lineage, k, cfg.seed, f);
    if let Err(e) = validate_folds(y, train, &inner, k) {
        out.reason = Some(format!("the inner folds: {e:#}"));
        return out;
    }
    let mut parts = Vec::with_capacity(k);
    for j in 0..k {
        let fit = Bits::from_indices(n, train.iter_ones().filter(|&i| inner[i] != j));
        let score = Bits::from_indices(n, train.iter_ones().filter(|&i| inner[i] == j));
        out.inner_folds.push(counts(&fit, &score, y));
        parts.push((fit, score));
    }
    for (j, (fit, _)) in parts.iter().enumerate() {
        if !any_training_candidate(x, fit, min_support) {
            out.reason = Some(format!("inner fit {j} has no candidate feature"));
            return out;
        }
        if let Err(e) = check_tuning_seeded(scheme, fit, y, tuning, cfg, inner_seed(cfg.seed, f, j))
        {
            out.reason = Some(format!("inner fit {j}: {e:#}"));
            return out;
        }
    }
    let fitted: Vec<_> = parts
        .iter()
        .enumerate()
        .map(|(j, (fit, _))| (j, training_candidates(x, fit, min_support).0))
        .collect();
    for &learner in LEARNERS.iter() {
        let mut row = InnerLearner {
            learner: learner.to_string(),
            folds: vec![],
            tp: 0,
            fp: 0,
            tn: 0,
            fn_: 0,
            balanced_accuracy: None,
            null_fits: 0,
            error: None,
        };
        for ((j, cx), (fit, score)) in fitted.iter().zip(&parts) {
            let seed = inner_seed(cfg.seed, f, *j);
            let context = LassoContext {
                dataset,
                scheme,
                outer_fold: Some(f),
            };
            match fit_seeded(learner, cx, y, fit, tuning, cfg, seed, seed, context) {
                Ok((scored, _model)) => {
                    let members: Vec<usize> = score.iter_ones().collect();
                    let t: Vec<bool> = members.iter().map(|&i| y.contains(i)).collect();
                    let s: Vec<f64> = match &scored.scores {
                        Some(all) => members.iter().map(|&i| all[i]).collect(),
                        None => {
                            row.null_fits += 1;
                            vec![train_prevalence(y, fit); members.len()]
                        }
                    };
                    let c = calls(&t, &s);
                    row.tp += c.tp;
                    row.fp += c.fp;
                    row.tn += c.tn;
                    row.fn_ += c.fn_;
                    row.folds.push(c);
                }
                Err(e) => {
                    row.error = Some(format!("inner fit {j}: {e:#}"));
                    row.folds.clear();
                    break;
                }
            }
        }
        if row.error.is_none() {
            row.balanced_accuracy = balanced(row.tp, row.fp, row.tn, row.fn_);
        } else {
            row.tp = 0;
            row.fp = 0;
            row.tn = 0;
            row.fn_ = 0;
        }
        out.learners.push(row);
    }
    let scores: Vec<Option<f64>> = out.learners.iter().map(|l| l.balanced_accuracy).collect();
    match inner_pick(&scores, set.tie) {
        Some((pick, tied)) => {
            out.status = "picked".into();
            out.tied = tied.iter().map(|&i| LEARNERS[i].to_string()).collect();
            out.pick = Some(LEARNERS[pick].to_string());
        }
        None => out.reason = Some("no learner has an inner score".into()),
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::learn::eval::{Role, evaluate, evaluate_assigned};
    use crate::learn::gbt::GbtParams;
    use crate::learn::tree::ForestParams;
    use serde_json::Value;

    /// SYNTHETIC: 150 genomes in 30 lineages of five, a third resistant,
    /// two features that follow the label with noise and six of noise.
    fn fixture() -> (BitMatrix, Bits, Vec<u32>) {
        use rand::{Rng, SeedableRng};
        let n = 150;
        let y = Bits::from_indices(n, (0..n).filter(|i| i % 3 == 0));
        let groups: Vec<u32> = (0..n as u32).map(|i| i / 5).collect();
        let mut rng = rand_chacha::ChaCha8Rng::seed_from_u64(20261122);
        let mut x = BitMatrix::with_cols(n);
        for f in 0..8 {
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

    fn config() -> EvalConfig {
        EvalConfig {
            folds: 3,
            inner_folds: 2,
            forest: ForestParams {
                n_trees: 6,
                ..Default::default()
            },
            gbt: GbtParams {
                rounds: 8,
                ..Default::default()
            },
            screened_features: 6,
            lasso_lambdas: 5,
            ..Default::default()
        }
    }

    fn run(k: usize, tie: PickTie, cfg: &EvalConfig) -> (Value, Vec<(String, String)>) {
        let (x, y, groups) = fixture();
        let (out, files) = capture(k, tie, || {
            evaluate("synthetic__drug", &x, &y, &groups, &[], 3, None, cfg).unwrap()
        });
        (
            strip_seconds(serde_json::to_value(&out.splits).unwrap()),
            files,
        )
    }

    fn strip_seconds(mut v: Value) -> Value {
        match &mut v {
            Value::Object(m) => {
                m.remove("seconds");
                for x in m.values_mut() {
                    *x = strip_seconds(x.take());
                }
            }
            Value::Array(a) => {
                for x in a.iter_mut() {
                    *x = strip_seconds(x.take());
                }
            }
            _ => {}
        }
        v
    }

    fn one(files: &[(String, String)]) -> Value {
        assert_eq!(files.len(), 1, "one result per dataset, split and seed");
        serde_json::from_str(&files[0].1).unwrap()
    }

    #[test]
    fn the_flag_writes_its_own_file_and_changes_no_split() {
        let (x, y, groups) = fixture();
        let cfg = config();
        let plain = evaluate("synthetic__drug", &x, &y, &groups, &[], 3, None, &cfg).unwrap();
        let (flagged, files) = run(3, PickTie::Later, &cfg);
        assert_eq!(
            strip_seconds(serde_json::to_value(&plain.splits).unwrap()),
            flagged
        );
        assert_eq!(files[0].0, "synthetic__drug.lineage-s42.pick.json");
        let r = one(&files);
        assert_eq!(r["format"], FORMAT);
        assert_eq!(r["inner_folds"], 3);
        assert_eq!(r["tested"], 150);
        // Off, nothing is captured and nothing is written.
        assert!(!enabled());
    }

    #[test]
    fn outer_calls_are_each_learners_report_row_fold_by_fold() {
        let cfg = config();
        let (splits, files) = run(3, PickTie::Later, &cfg);
        let r = one(&files);
        let lineage = splits
            .as_array()
            .unwrap()
            .iter()
            .find(|s| s["scheme"] == "lineage")
            .unwrap();
        for (li, l) in LEARNERS.iter().enumerate() {
            let m = &lineage["learners"][li]["metrics"];
            for key in ["tp", "fp", "tn", "fn"] {
                let sum: u64 = r["folds"]
                    .as_array()
                    .unwrap()
                    .iter()
                    .map(|f| f["outer"][li][key].as_u64().unwrap())
                    .sum();
                assert_eq!(sum, m[key].as_u64().unwrap(), "{l} {key}");
            }
            assert_eq!(
                r["learners"][li]["balanced_accuracy"],
                m["balanced_accuracy"]
            );
        }
        // Each outer fold's counts are its test samples' own.
        let folds = lineage["test_fold"].as_array().unwrap();
        for f in r["folds"].as_array().unwrap() {
            let fi = f["fold"].as_u64().unwrap();
            let tested = folds.iter().filter(|v| v.as_u64() == Some(fi)).count() as u64;
            assert_eq!(f["test"].as_u64().unwrap(), tested);
            assert_eq!(f["outer"][0]["n"].as_u64().unwrap(), tested);
            assert_eq!(f["train"].as_u64().unwrap() + tested, 150);
        }
    }

    /// `x` and `y` with every sample outside `train` changed: each feature
    /// bit flipped and the label flipped.
    fn scrambled(x: &BitMatrix, y: &Bits, train: &Bits) -> (BitMatrix, Bits) {
        let n = x.cols;
        let mut x2 = BitMatrix::with_cols(n);
        for r in 0..x.rows {
            let row = (0..n).filter(|&c| x.get(r, c) != !train.contains(c));
            x2.push_row(Bits::from_indices(n, row).words());
        }
        let y2 = Bits::from_indices(n, (0..n).filter(|&i| y.contains(i) == train.contains(i)));
        (x2, y2)
    }

    #[test]
    fn no_outer_test_genome_reaches_an_inner_fit() {
        let (x, y, groups) = fixture();
        let cfg = config();
        let folds = crate::learn::grouped_folds(&y, &groups, cfg.folds, cfg.seed);
        let set = Settings {
            dir: PathBuf::new(),
            inner_folds: 3,
            tie: PickTie::Later,
        };
        for (f, (train, test)) in fold_masks(&folds, cfg.folds).into_iter().enumerate() {
            let (x2, y2) = scrambled(&x, &y, &train);
            let tuning = Some(groups.as_slice());
            let pick = |x: &BitMatrix, y: &Bits| {
                let p = inner_fold_pick(
                    "d", "lineage", x, y, &groups, tuning, &train, f, 3, &cfg, &set,
                );
                serde_json::to_value(p).unwrap()
            };
            let a = pick(&x, &y);
            assert_eq!(a["status"], "picked", "fold {f}: {a}");
            assert_eq!(
                a,
                pick(&x2, &y2),
                "fold {f}: a test genome changed the inner pick"
            );
            // Every inner fold lies inside the outer training samples.
            let part = inner_partition(&y, &train, &groups, 3, cfg.seed, f);
            assert!(test.iter_ones().all(|i| part[i] == usize::MAX));
            assert!(train.iter_ones().all(|i| part[i] < 3));
            let inner_total: u64 = a["inner_folds"]
                .as_array()
                .unwrap()
                .iter()
                .map(|c| c["test"].as_u64().unwrap())
                .sum();
            assert_eq!(inner_total, train.count() as u64);
        }
    }

    #[test]
    fn inner_folds_keep_each_lineage_whole_and_are_the_release_tuning_partition() {
        let (_, y, groups) = fixture();
        let cfg = config();
        let folds = crate::learn::grouped_folds(&y, &groups, cfg.folds, cfg.seed);
        for (f, (train, _)) in fold_masks(&folds, cfg.folds).into_iter().enumerate() {
            let part = inner_partition(&y, &train, &groups, 3, cfg.seed, f);
            let mut fold_of = std::collections::HashMap::new();
            for i in train.iter_ones() {
                let at = *fold_of.entry(groups[i]).or_insert(part[i]);
                assert_eq!(
                    at, part[i],
                    "lineage {} split across inner folds",
                    groups[i]
                );
            }
            let release = training_folds(
                &y,
                &train,
                Some(&groups),
                3,
                crate::ablate::plan::fold_seed(cfg.seed, f),
            );
            assert_eq!(part, release);
        }
    }

    #[test]
    fn the_seed_keys_every_inner_fit_and_partition() {
        // FNV-1a 64 of the text, computed apart from this code.
        assert_eq!(inner_seed(42, 0, 0), 6374775588177860159);
        assert_eq!(inner_seed(43, 0, 0), 5429478203028429340);
        assert_eq!(inner_seed(42, 1, 0), 6993236585231890516);
        assert_eq!(inner_seed(42, 0, 1), 6374774488666231948);
        let (_, y, groups) = fixture();
        let train = Bits::from_indices(150, 0..120);
        let a = inner_partition(&y, &train, &groups, 3, 42, 0);
        assert_ne!(a, inner_partition(&y, &train, &groups, 3, 43, 0));
        assert_ne!(a, inner_partition(&y, &train, &groups, 3, 42, 1));
        // The written seed is the run's, and another seed writes another file.
        let (_, files) = run(
            3,
            PickTie::Later,
            &EvalConfig {
                seed: 43,
                ..config()
            },
        );
        let r = one(&files);
        assert_eq!(files[0].0, "synthetic__drug.lineage-s43.pick.json");
        assert_eq!(r["seed"], 43);
    }

    #[test]
    fn a_tie_goes_where_the_rule_says_and_is_listed() {
        let tie = |s: &[Option<f64>], t| crate::report::inner_pick(s, t);
        let s = [Some(0.8), Some(0.9), Some(0.9 + 1e-15), Some(0.7), None];
        assert_eq!(tie(&s, PickTie::Later), Some((2, vec![1, 2])));
        assert_eq!(tie(&s, PickTie::Earlier), Some((1, vec![1, 2])));
        let s = [Some(0.6), None, Some(0.6), None, Some(0.6)];
        assert_eq!(tie(&s, PickTie::Later), Some((4, vec![0, 2, 4])));
        assert_eq!(tie(&s, PickTie::Earlier), Some((0, vec![0, 2, 4])));
        assert_eq!(
            tie(&[Some(0.5), Some(0.5000001)], PickTie::Later),
            Some((1, vec![1]))
        );
        assert_eq!(tie(&[None, None], PickTie::Later), None);
        // The file takes the same rule.
        let (_, files) = run(3, PickTie::Earlier, &config());
        let r = one(&files);
        assert_eq!(r["tie"], "earlier");
        for f in r["folds"].as_array().unwrap() {
            let tied = f["inner"]["tied"].as_array().unwrap();
            assert_eq!(f["inner"]["pick"], tied[0]);
        }
    }

    fn row(learner: &str, random: f64, lineage: f64) -> (Value, Value) {
        let m =
            |ba: f64| serde_json::json!({"learner": learner, "metrics": {"balanced_accuracy": ba}});
        (m(random), m(lineage))
    }

    #[test]
    fn the_default_pick_reads_the_random_split_as_the_release_does() {
        let rows = [
            row("scm", 0.70, 0.60),
            row("cart", 0.71, 0.61),
            row("rf", 0.80, 0.62),
            row("gbt", 0.79, 0.90),
            row("lasso", 0.80, 0.63),
        ];
        let splits = serde_json::json!([
            {"scheme": "random", "learners": rows.iter().map(|r| r.0.clone()).collect::<Vec<_>>()},
            {"scheme": "lineage", "unavailable_reason": null, "learners": rows.iter().map(|r| r.1.clone()).collect::<Vec<_>>()},
        ]);
        let (learner, ba, tied) = release_pick(&splits).unwrap();
        assert_eq!(
            learner, "lasso",
            "the random split's best, a tie to the later learner"
        );
        assert_eq!(ba, 0.80);
        assert_eq!(tied, vec!["rf", "lasso"]);
        // And the file's default is the same rule's pick from the same run.
        let (splits, files) = run(3, PickTie::Later, &config());
        let r = one(&files);
        let (learner, _, tied) = release_pick(&splits).unwrap();
        assert_eq!(r["default"]["learner"], learner.as_str());
        assert_eq!(r["default"]["tied"], serde_json::json!(tied));
    }

    #[test]
    fn the_twin_runs_fit_no_inner_fold() {
        let (_, files) = run(0, PickTie::Later, &config());
        let r = one(&files);
        assert_eq!(r["inner_folds"], 0);
        for f in r["folds"].as_array().unwrap() {
            assert_eq!(f["inner"]["status"], "not run");
            assert!(f["inner"]["learners"].as_array().unwrap().is_empty());
            assert_eq!(f["outer"].as_array().unwrap().len(), 5);
        }
        assert!(r["default"]["learner"].is_string());
    }

    #[test]
    fn a_supplied_split_writes_its_pick_with_no_default() {
        let (x, y, groups) = fixture();
        let cfg = config();
        // SYNTHETIC study-like plan: three folds of ten lineages each.
        let roles: Vec<Role> = groups
            .iter()
            .map(|&g| Role::Test((g % 3) as usize))
            .collect();
        let (out, files) = capture(3, PickTie::Later, || {
            evaluate_assigned(
                "synthetic__drug",
                "study-whole",
                &x,
                &y,
                &groups,
                &[],
                &roles,
                3,
                &cfg,
            )
            .unwrap()
        });
        assert_eq!(files[0].0, "synthetic__drug.study-whole-s42.pick.json");
        let r = one(&files);
        assert!(r["default"].is_null());
        assert_eq!(r["tested"], 150);
        for (li, l) in out.result.learners.iter().take(5).enumerate() {
            let m = l.metrics.as_ref().unwrap();
            let tp: u64 = r["folds"]
                .as_array()
                .unwrap()
                .iter()
                .map(|f| f["outer"][li]["tp"].as_u64().unwrap())
                .sum();
            assert_eq!(tp, m.tp as u64);
        }
        for f in r["folds"].as_array().unwrap() {
            assert_eq!(f["inner"]["status"], "picked");
        }
    }
}
