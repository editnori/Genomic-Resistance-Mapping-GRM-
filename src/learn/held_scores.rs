//! Experiment, off by default (plan block C3's run on its datasets of M,
//! `docs/assumption-plan-corrections.md`, entry M run, Scores kept): write
//! each held-out genome's lineage score from a drawn fit to a file of its
//! own, so that a scorer can read the fit again with some genomes left out
//! and refit nothing.
//!
//! Nothing here changes a fit or a report. `grm learn --exp-scores DIR`
//! turns it on for that process, beside `--exp-pick`; without the flag
//! `enabled` is false and nothing is written. With it, `analyze_dataset`
//! fits, scores and reports exactly as without it, and then writes one file
//! per dataset and seed, `DIR/<dataset>.lineage-s<seed>.scores.json`, from
//! the lineage split the fit made: the cohort's genomes in sample order,
//! each with its reported class, its outer test fold and the five learners'
//! out-of-fold scores, the values `cutoff::calls` reads for the pick file's
//! `outer` calls.
//!
//! The file names genomes, so it is written on the run's machine and stays
//! there. A study fit (`--fold-file`) writes no such file: its
//! `FoldFileResult` already holds each tested genome's scores.

use super::eval::{LEARNERS, SplitResult};
use anyhow::{Context, Result, ensure};
use serde::Serialize;
use std::path::{Path, PathBuf};
use std::sync::OnceLock;

/// The file format. A change to what a file holds changes it.
pub const FORMAT: &str = "grm-c3-scores-1";

/// The only split a file holds: the drawn split the pick file's counts are
/// taken on.
pub const SCHEME: &str = "lineage";

static SINK: OnceLock<PathBuf> = OnceLock::new();

#[cfg(test)]
thread_local! {
    /// Tests only: while `capture` runs on a thread, that thread's files go
    /// here, by file name.
    static CAPTURED: std::cell::RefCell<Option<Vec<(String, String)>>> =
        const { std::cell::RefCell::new(None) };
}

/// Turns the experiment on for this process. `dir` is created if missing.
/// A second call fails.
pub fn enable(dir: &Path) -> Result<()> {
    std::fs::create_dir_all(dir).with_context(|| format!("create {}", dir.display()))?;
    SINK.set(dir.to_path_buf())
        .map_err(|_| anyhow::anyhow!("the C3 score folder is already set"))
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

/// Tests only: runs `f` with the experiment on for this thread alone and
/// returns `f`'s result with each file it would have written, as (file
/// name, JSON text).
#[cfg(test)]
pub(crate) fn capture<R>(f: impl FnOnce() -> R) -> (R, Vec<(String, String)>) {
    CAPTURED.with(|c| *c.borrow_mut() = Some(Vec::new()));
    let result = f();
    let files = CAPTURED.with(|c| c.borrow_mut().take()).unwrap_or_default();
    (result, files)
}

fn write(name: &str, text: String) -> Result<()> {
    #[cfg(test)]
    {
        let taken = CAPTURED.with(|c| {
            c.borrow_mut()
                .as_mut()
                .map(|v| v.push((name.to_string(), text.clone())))
                .is_some()
        });
        if taken {
            return Ok(());
        }
    }
    let dir = SINK.get().context("the C3 score folder is not set")?;
    let path = dir.join(name);
    let tmp = dir.join(format!("{name}.tmp"));
    std::fs::write(&tmp, text).with_context(|| format!("write {}", tmp.display()))?;
    std::fs::rename(&tmp, &path).with_context(|| format!("replace {}", path.display()))
}

/// One cohort genome: its id, reported class, outer test fold and the five
/// learners' out-of-fold scores in pipeline order.
#[derive(Clone, Debug, Serialize)]
pub struct HeldGenome {
    pub genome_id: String,
    pub resistant: bool,
    pub fold: usize,
    pub scores: Vec<f64>,
}

/// What `--exp-scores` writes for one dataset and seed.
#[derive(Clone, Debug, Serialize)]
pub struct HeldScores {
    pub format: &'static str,
    pub evaluation_revision: u32,
    pub dataset: String,
    pub scheme: &'static str,
    pub seed: u64,
    /// `k31` for the default and the twin, else the `--exp-inputs` arm.
    pub arm: String,
    pub learners: Vec<String>,
    /// Outer folds; 0 when the split could not be drawn.
    pub folds: usize,
    pub unavailable_reason: Option<String>,
    /// Empty when the split could not be drawn.
    pub genomes: Vec<HeldGenome>,
}

/// The file name for a dataset and seed.
pub fn file_name(dataset: &str, seed: u64) -> String {
    format!("{dataset}.{SCHEME}-s{seed}.scores.json")
}

/// The fit's arm: the `--exp-inputs` arm, or `k31`.
fn arm() -> String {
    crate::input_arms::active()
        .map(|s| s.arm.name().to_string())
        .unwrap_or_else(|| "k31".to_string())
}

/// The file's contents from a drawn fit's splits. `genome_ids` and
/// `resistant` are the cohort's, in the sample order the fit used.
pub fn build(
    dataset: &str,
    seed: u64,
    arm: String,
    genome_ids: &[String],
    resistant: &[bool],
    splits: &[SplitResult],
) -> Result<HeldScores> {
    let lineage = splits
        .iter()
        .find(|s| s.scheme == SCHEME)
        .context("the fit wrote no lineage split")?;
    let mut out = HeldScores {
        format: FORMAT,
        evaluation_revision: crate::analysis::EVALUATION_REVISION,
        dataset: dataset.to_string(),
        scheme: SCHEME,
        seed,
        arm,
        learners: LEARNERS.iter().map(|l| l.to_string()).collect(),
        folds: 0,
        unavailable_reason: lineage.unavailable_reason.clone(),
        genomes: vec![],
    };
    if out.unavailable_reason.is_some() {
        return Ok(out);
    }
    let n = genome_ids.len();
    ensure!(resistant.len() == n, "one class per cohort genome");
    ensure!(
        lineage.test_fold.len() == n,
        "the lineage split tests {} genomes of {n}",
        lineage.test_fold.len()
    );
    ensure!(
        lineage.learners.len() >= LEARNERS.len()
            && LEARNERS
                .iter()
                .zip(&lineage.learners)
                .all(|(l, r)| r.learner == *l && r.scores.len() == n),
        "the lineage split's learners are not the five in pipeline order, one score per genome"
    );
    out.folds = lineage.folds;
    out.genomes = (0..n)
        .map(|i| HeldGenome {
            genome_id: genome_ids[i].clone(),
            resistant: resistant[i],
            fold: lineage.test_fold[i],
            scores: lineage.learners[..LEARNERS.len()]
                .iter()
                .map(|r| r.scores[i])
                .collect(),
        })
        .collect();
    Ok(out)
}

/// Writes the file for a drawn fit, when the experiment is on.
pub(crate) fn record(
    dataset: &str,
    seed: u64,
    genome_ids: &[String],
    resistant: &[bool],
    splits: &[SplitResult],
) -> Result<()> {
    let held = build(dataset, seed, arm(), genome_ids, resistant, splits)?;
    write(&file_name(dataset, seed), serde_json::to_string(&held)?)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::learn::eval::LearnerResult;

    fn learner(name: &str, scores: Vec<f64>) -> LearnerResult {
        LearnerResult {
            learner: name.to_string(),
            metrics: None,
            called_metrics: None,
            fold_balanced_accuracy: vec![],
            roc: vec![],
            scores,
            conformal: vec![],
            fold_complexity: vec![],
            seconds: 0.0,
            no_call: None,
            lasso_fits: None,
        }
    }

    fn split(unavailable: Option<&str>) -> SplitResult {
        SplitResult {
            scheme: SCHEME.to_string(),
            folds: if unavailable.is_some() { 0 } else { 2 },
            n_groups: 2,
            test_fold: if unavailable.is_some() {
                vec![]
            } else {
                vec![0, 1, 0]
            },
            learners: LEARNERS
                .iter()
                .enumerate()
                .map(|(j, l)| learner(l, vec![0.1 * j as f64, 0.5, 0.9]))
                .collect(),
            unavailable_reason: unavailable.map(str::to_string),
        }
    }

    #[test]
    fn held_scores_file_lists_each_genome_with_its_fold_and_five_scores_in_sample_order() {
        let ids: Vec<String> = ["a", "b", "c"].iter().map(|s| s.to_string()).collect();
        let held = build(
            "sp__drug",
            43,
            "k31".into(),
            &ids,
            &[true, false, true],
            &[split(None)],
        )
        .unwrap();
        assert_eq!(
            (held.format, held.scheme, held.seed, held.folds),
            (FORMAT, SCHEME, 43, 2)
        );
        assert_eq!(
            held.learners,
            LEARNERS.iter().map(|l| l.to_string()).collect::<Vec<_>>()
        );
        let rows: Vec<(&str, bool, usize)> = held
            .genomes
            .iter()
            .map(|g| (g.genome_id.as_str(), g.resistant, g.fold))
            .collect();
        assert_eq!(rows, vec![("a", true, 0), ("b", false, 1), ("c", true, 0)]);
        assert_eq!(
            held.genomes[0].scores,
            vec![0.0, 0.1, 0.2, 0.30000000000000004, 0.4]
        );
        assert_eq!(
            file_name("sp__drug", 43),
            "sp__drug.lineage-s43.scores.json"
        );
    }

    #[test]
    fn held_scores_of_an_undrawn_split_give_its_reason_and_no_genome() {
        let ids: Vec<String> = ["a", "b", "c"].iter().map(|s| s.to_string()).collect();
        let held = build(
            "sp__drug",
            42,
            "unitigs".into(),
            &ids,
            &[true, false, true],
            &[split(Some("too few lineages"))],
        )
        .unwrap();
        assert_eq!(held.unavailable_reason.as_deref(), Some("too few lineages"));
        assert_eq!(
            (held.folds, held.genomes.len(), held.arm.as_str()),
            (0, 0, "unitigs")
        );
    }

    #[test]
    fn held_scores_refuse_a_split_whose_genomes_do_not_match_the_cohort() {
        let ids: Vec<String> = ["a", "b"].iter().map(|s| s.to_string()).collect();
        assert!(
            build(
                "sp__drug",
                42,
                "k31".into(),
                &ids,
                &[true, false],
                &[split(None)]
            )
            .is_err()
        );
        let three: Vec<String> = ["a", "b", "c"].iter().map(|s| s.to_string()).collect();
        let mut short = split(None);
        short.learners[4].scores.pop();
        assert!(
            build(
                "sp__drug",
                42,
                "k31".into(),
                &three,
                &[true, false, true],
                &[short]
            )
            .is_err()
        );
    }

    #[test]
    fn held_scores_are_off_unless_enabled_or_captured() {
        assert!(!enabled());
        let (on, files) = capture(enabled);
        assert!(on && files.is_empty());
        assert!(!enabled());
    }
}
