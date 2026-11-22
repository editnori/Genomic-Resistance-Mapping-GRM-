//! Plan block C1 part B (docs/cohort-draw.md): how far an estimate from a
//! small cohort can be trusted. `grm learn --exp-cohort` takes cohorts of r
//! resistant and r susceptible genomes drawn from a dataset
//! (`scripts/c1-cutoff/draw.py` writes them), runs the release's rule on
//! each cohort as if it were a dataset of 2r genomes, and scores the
//! cohort's final models on the dataset's other genomes, the truth the
//! estimate stood for. This module holds what that needs apart from the
//! data root: the cohort file's checks, the scored set, the value of a
//! cohort feature on a genome outside the cohort, the call counts, and the
//! fits of one cohort. `analysis::cohort_draw_species` reads the data root
//! and writes the results.
//!
//! A result holds counts only: no genome id, no row per genome, no k-mer and
//! no feature. Any split here is a retrospective split inside one public
//! database, never external validation.

use crate::bits::{self, BitMatrix, Bits};
use crate::learn::NoCall;
use crate::learn::eval::{self, EvalConfig, LEARNERS, LearnerResult, Role, SplitResult};
use crate::predict::CALL_CUTOFF;
use anyhow::{Result, bail, ensure};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeSet, HashMap, HashSet};
use std::path::{Path, PathBuf};

/// The `format` of a cohort file.
pub const COHORT_FORMAT: &str = "grm-cohort-draw-v1";
/// The `format` of every result.
pub const FORMAT: &str = "grm-cohort-draw-result-v1";
/// Log loss clips a score to [LOSS_CLIP, 1 − LOSS_CLIP], as the time split's
/// record does.
pub const LOSS_CLIP: f64 = 1e-4;
/// The study plan a cohort's study folds come from (docs/study-split.md).
pub const STUDY_SCHEME: &str = "study-whole";
/// Outer folds of the study plan.
pub const STUDY_FOLDS: usize = 5;

/// One drawn cohort, as `draw.py` writes it. It holds genome ids, so it
/// stays where it is made.
#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CohortFile {
    pub format: String,
    pub dataset: String,
    /// Genomes of each class.
    pub r: usize,
    pub draw: usize,
    pub genomes: Vec<String>,
    /// Each genome's reported call, as the draw read it.
    pub resistant: Vec<bool>,
    /// Each genome's lineage number, as the draw read it.
    pub lineage: Vec<u32>,
    /// Each genome's fold in the dataset's `study-whole` plan, for a dataset
    /// in the study set.
    #[serde(default)]
    pub study_fold: Option<Vec<usize>>,
}

/// `<dir>/<dataset>.r<r>.d<draw>.cohort.json`, the name `draw.py` writes.
pub fn cohort_name(dataset: &str, r: usize, draw: usize) -> String {
    format!("{dataset}.r{r}.d{draw}.cohort.json")
}

/// Where the result of a cohort at `seed` goes.
pub fn result_path(out: &Path, dataset: &str, r: usize, draw: usize, seed: u64) -> PathBuf {
    out.join(format!("{dataset}.r{r}.d{draw}.s{seed}.cohort.json"))
}

/// The dataset, r and draw a cohort file's name gives, or None for any
/// other name.
pub fn parse_cohort_name(name: &str) -> Option<(String, usize, usize)> {
    let stem = name.strip_suffix(".cohort.json")?;
    let (rest, draw) = stem.rsplit_once(".d")?;
    let (dataset, r) = rest.rsplit_once(".r")?;
    let ok = |s: &str| !s.is_empty() && s.bytes().all(|b| b.is_ascii_digit());
    (ok(r) && ok(draw) && !dataset.is_empty())
        .then(|| Some((dataset.to_string(), r.parse().ok()?, draw.parse().ok()?)))
        .flatten()
}

/// The dataset rows a cohort file names, ascending, after checking the file
/// against the dataset: its format and dataset, one call, lineage and study
/// fold per genome, r of each class by its own calls, no genome twice, and
/// each genome a row of the dataset with the same call and lineage. A file
/// that fails any check is refused with the reason, so a class short of r is
/// never taken silently.
pub fn cohort_rows(
    file: &CohortFile,
    dataset: &str,
    ids: &[&str],
    resistant: &[bool],
    lineage: &[u32],
) -> Result<Vec<usize>> {
    ensure!(
        file.format == COHORT_FORMAT,
        "the cohort file's format is {:?}, not {COHORT_FORMAT}",
        file.format
    );
    ensure!(
        file.dataset == dataset,
        "the cohort file names {}, not {dataset}",
        file.dataset
    );
    let n = file.genomes.len();
    ensure!(
        file.resistant.len() == n
            && file.lineage.len() == n
            && file.study_fold.as_ref().is_none_or(|f| f.len() == n),
        "the cohort file's lists differ in length"
    );
    ensure!(file.r > 0, "a cohort needs r of at least 1");
    let r_count = file.resistant.iter().filter(|&&r| r).count();
    ensure!(
        r_count == file.r && n - r_count == file.r,
        "the cohort holds {r_count} resistant and {} susceptible genomes, not {} of each",
        n - r_count,
        file.r
    );
    let row_of: HashMap<&str, usize> = ids.iter().enumerate().map(|(i, &g)| (g, i)).collect();
    let mut rows = Vec::with_capacity(n);
    let mut seen = HashSet::new();
    for (k, g) in file.genomes.iter().enumerate() {
        ensure!(seen.insert(g.as_str()), "the cohort lists genome {k} twice");
        let Some(&i) = row_of.get(g.as_str()) else {
            bail!("cohort genome {k} is not a row of the dataset");
        };
        ensure!(
            resistant[i] == file.resistant[k] && lineage[i] == file.lineage[k],
            "cohort genome {k} has another call or lineage in the dataset"
        );
        rows.push(i);
    }
    rows.sort_unstable();
    Ok(rows)
}

/// The dataset rows a cohort's final models are scored on.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ScoredSet {
    /// Rows outside the cohort that are no k-mer-identical copy of a cohort
    /// genome, ascending.
    pub rows: Vec<usize>,
    /// For each of `rows`, whether its lineage holds no cohort genome.
    pub new_lineage: Vec<bool>,
    /// Rows outside the cohort left out as copies of a cohort genome.
    pub copies: Vec<usize>,
}

/// The scored set of `cohort` (dataset rows) among a dataset's `digests`
/// (each row's kset digest) and `lineage`s.
pub fn scored_set(cohort: &[usize], lineage: &[u32], digests: &[&str]) -> ScoredSet {
    let n = lineage.len();
    let inside: HashSet<usize> = cohort.iter().copied().collect();
    let kept: HashSet<&str> = cohort.iter().map(|&i| digests[i]).collect();
    let lineages: HashSet<u32> = cohort.iter().map(|&i| lineage[i]).collect();
    let (mut rows, mut new_lineage, mut copies) = (Vec::new(), Vec::new(), Vec::new());
    for i in 0..n {
        if inside.contains(&i) {
            continue;
        }
        if kept.contains(digests[i]) {
            copies.push(i);
            continue;
        }
        rows.push(i);
        new_lineage.push(!lineages.contains(&lineage[i]));
    }
    ScoredSet {
        rows,
        new_lineage,
        copies,
    }
}

/// The cohort's features on genomes outside it: row f, column j is the bit
/// of feature f's lowest species pattern (`first_pattern[f]`) on species
/// column `columns[j]`. That is whether the genome carries the k-mer the
/// feature's card names (its first k-mer, `annotate`), which is how
/// `grm predict` reads a feature on a genome outside the cohort; every
/// k-mer of one species pattern has the same presence on every genome of
/// the species, so the bit is exact.
pub fn feature_values(species: &BitMatrix, first_pattern: &[u32], columns: &[usize]) -> BitMatrix {
    let mut out = BitMatrix::with_cols(columns.len());
    let mut row = vec![0u64; bits::words_for(columns.len())];
    for &p in first_pattern {
        row.iter_mut().for_each(|w| *w = 0);
        for (j, &c) in columns.iter().enumerate() {
            if species.get(p as usize, c) {
                bits::set(&mut row, j);
            }
        }
        out.push_row(&row);
    }
    out
}

/// Each projected feature's lowest species pattern, the `first_pattern`
/// `feature_values` reads.
pub fn lowest_patterns(ds: &crate::matrix::Dataset) -> Vec<u32> {
    (0..ds.n_features())
        .map(|f| {
            *ds.patterns_of(f)
                .iter()
                .min()
                .expect("a feature covers a pattern")
        })
        .collect()
}

/// Calls at `CALL_CUTOFF` by reported class, and the log loss sum with each
/// score clipped to [LOSS_CLIP, 1 − LOSS_CLIP].
#[derive(Clone, Copy, Debug, Default, PartialEq, Serialize)]
pub struct Calls {
    pub n: usize,
    pub tp: usize,
    pub fp: usize,
    pub tn: usize,
    #[serde(rename = "fn")]
    pub fn_: usize,
    pub loss_sum: f64,
}

/// `Calls` of `scores` against `truth`.
pub fn calls(truth: &[bool], scores: &[f64]) -> Calls {
    assert_eq!(
        truth.len(),
        scores.len(),
        "one score is required per genome"
    );
    let mut c = Calls {
        n: truth.len(),
        ..Calls::default()
    };
    for (&t, &s) in truth.iter().zip(scores) {
        let p = s.clamp(LOSS_CLIP, 1.0 - LOSS_CLIP);
        match (t, s >= CALL_CUTOFF) {
            (true, true) => c.tp += 1,
            (true, false) => c.fn_ += 1,
            (false, true) => c.fp += 1,
            (false, false) => c.tn += 1,
        }
        c.loss_sum -= if t { p.ln() } else { (1.0 - p).ln() };
    }
    c
}

/// Counts by class.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize)]
pub struct Counts {
    pub samples: usize,
    pub resistant: usize,
}

impl Counts {
    pub fn of(resistant: impl IntoIterator<Item = bool>) -> Counts {
        let mut c = Counts::default();
        for r in resistant {
            c.samples += 1;
            c.resistant += r as usize;
        }
        c
    }
}

/// One learner's estimate on one split of the cohort.
#[derive(Clone, Debug, Serialize)]
pub struct EstimateRow {
    pub learner: String,
    /// Over every genome the split scores, from the split's pooled scores.
    pub calls: Calls,
    /// The pooled balanced accuracy the split writes (`metrics`).
    pub balanced_accuracy: f64,
    /// Each outer fold's balanced accuracy, by fold id.
    pub fold_balanced_accuracy: Vec<f64>,
    /// Genomes of null-fit folds, by class (`LearnerResult::no_call`).
    pub no_call: Option<NoCall>,
}

/// One split of the cohort: the release's lineage and random splits, and
/// the study plan's folds for a cohort in the study set.
#[derive(Clone, Debug, Serialize)]
pub struct SplitEstimate {
    pub scheme: String,
    /// Lineages among the genomes the split tests or trains.
    pub n_groups: usize,
    /// Why the release's rule gives no estimate on this split, as it says.
    pub unavailable: Option<String>,
    pub learners: Vec<EstimateRow>,
}

/// One learner's final model on the scored set.
#[derive(Clone, Debug, Serialize)]
pub struct TruthRow {
    pub learner: String,
    /// The final fit kept nothing (boosting no tree, the lasso no weight),
    /// so every scored genome takes the cohort's resistant share.
    pub null_fit: bool,
    pub all: Calls,
    pub new_lineage: Calls,
}

/// The truth, or why there is none.
#[derive(Clone, Debug, Serialize)]
pub struct Truth {
    pub unavailable: Option<String>,
    pub learners: Vec<TruthRow>,
}

/// What one cohort gives, apart from its counts.
pub struct CohortFit {
    pub estimates: Vec<SplitEstimate>,
    pub truth: Truth,
    /// The lineage split's outer fold of each cohort genome, empty when the
    /// split is unavailable. Kept for the tests; no result writes it.
    pub lineage_test_fold: Vec<usize>,
}

/// One cohort's inputs: the cohort's projection at its own support floor
/// and its lineages and copy groups, as `evaluate` reads a dataset's, the
/// study plan's roles when the dataset is in the study set, and the scored
/// set's feature values (`feature_values`), calls and new-lineage flags.
pub struct CohortInputs<'a> {
    pub x: &'a BitMatrix,
    pub y: &'a Bits,
    pub groups: &'a [u32],
    pub duplicates: &'a [Vec<usize>],
    pub min_support: u32,
    pub study_roles: Option<&'a [Role]>,
    pub scored_x: &'a BitMatrix,
    pub scored_y: &'a [bool],
    pub scored_new: &'a [bool],
}

fn estimate_row(l: &LearnerResult, truth: &[bool]) -> EstimateRow {
    EstimateRow {
        learner: l.learner.clone(),
        calls: calls(truth, &l.scores),
        balanced_accuracy: l.metrics.as_ref().map_or(f64::NAN, |m| m.balanced_accuracy),
        fold_balanced_accuracy: l.fold_balanced_accuracy.clone(),
        no_call: l.no_call,
    }
}

fn split_estimate(s: &SplitResult, truth: &[bool], n_groups: usize) -> SplitEstimate {
    SplitEstimate {
        scheme: s.scheme.clone(),
        n_groups,
        unavailable: s.unavailable_reason.clone(),
        learners: s.learners.iter().map(|l| estimate_row(l, truth)).collect(),
    }
}

/// The release's rule on one cohort, and its final models on the scored
/// set. `evaluate` gives the lineage and the random split; when it stops
/// (a random split it cannot draw), both are unavailable with its reason.
/// The study roles go through `evaluate_assigned`, whose refusal makes the
/// study split unavailable. The final models (`eval::final_scores`) score
/// the scored set unless it lacks a class.
pub fn fit_cohort(dataset: &str, c: &CohortInputs, cfg: &EvalConfig) -> CohortFit {
    let truth: Vec<bool> = (0..c.x.cols).map(|i| c.y.contains(i)).collect();
    let n_groups = c.groups.iter().collect::<BTreeSet<_>>().len();
    let mut estimates = Vec::new();
    let mut lineage_test_fold = Vec::new();
    match eval::evaluate(
        dataset,
        c.x,
        c.y,
        c.groups,
        c.duplicates,
        c.min_support,
        None,
        cfg,
    ) {
        Ok(out) => {
            for s in &out.splits {
                if s.scheme == "lineage" {
                    lineage_test_fold = s.test_fold.clone();
                }
                estimates.push(split_estimate(s, &truth, n_groups));
            }
        }
        Err(e) => {
            for scheme in ["random", "lineage"] {
                estimates.push(SplitEstimate {
                    scheme: scheme.into(),
                    n_groups,
                    unavailable: Some(format!("{e:#}")),
                    learners: vec![],
                });
            }
        }
    }
    if let Some(roles) = c.study_roles {
        estimates.push(
            match eval::evaluate_assigned(
                dataset,
                STUDY_SCHEME,
                c.x,
                c.y,
                c.groups,
                c.duplicates,
                roles,
                STUDY_FOLDS,
                cfg,
            ) {
                Ok(a) => {
                    let tested: Vec<bool> = a.tested.iter().map(|&i| truth[i]).collect();
                    let groups: BTreeSet<u32> = (0..roles.len())
                        .filter(|&i| roles[i] != Role::Unused)
                        .map(|i| c.groups[i])
                        .collect();
                    split_estimate(&a.result, &tested, groups.len())
                }
                Err(e) => SplitEstimate {
                    scheme: STUDY_SCHEME.into(),
                    n_groups,
                    unavailable: Some(format!("{e:#}")),
                    learners: vec![],
                },
            },
        );
    }
    let truth = scored_truth(dataset, c, cfg);
    CohortFit {
        estimates,
        truth,
        lineage_test_fold,
    }
}

/// A null final fit's score on every scored genome: the cohort's resistant
/// share, 0.5 for r of each class, never the scored set's.
fn null_fit_score(c: &CohortInputs) -> f64 {
    c.y.count() as f64 / c.x.cols as f64
}

fn scored_truth(dataset: &str, c: &CohortInputs, cfg: &EvalConfig) -> Truth {
    let unavailable = |why: String| Truth {
        unavailable: Some(why),
        learners: vec![],
    };
    let resistant = c.scored_y.iter().filter(|&&r| r).count();
    if resistant == 0 || resistant == c.scored_y.len() {
        return unavailable("the scored genomes lack a class".into());
    }
    let (_, scores) = match eval::final_scores(
        dataset,
        c.x,
        c.y,
        c.groups,
        c.duplicates,
        c.min_support,
        cfg,
        c.scored_x,
    ) {
        Ok(f) => f,
        Err(e) => return unavailable(format!("{e:#}")),
    };
    let share = null_fit_score(c);
    let new: Vec<usize> = (0..c.scored_new.len())
        .filter(|&j| c.scored_new[j])
        .collect();
    let new_y: Vec<bool> = new.iter().map(|&j| c.scored_y[j]).collect();
    let learners = LEARNERS
        .iter()
        .zip(scores)
        .map(|(name, scores)| {
            let null_fit = scores.is_none();
            let scores = scores.unwrap_or_else(|| vec![share; c.scored_y.len()]);
            let new_scores: Vec<f64> = new.iter().map(|&j| scores[j]).collect();
            TruthRow {
                learner: name.to_string(),
                null_fit,
                all: calls(c.scored_y, &scores),
                new_lineage: calls(&new_y, &new_scores),
            }
        })
        .collect();
    Truth {
        unavailable: None,
        learners,
    }
}

/// The cohort's counts.
#[derive(Clone, Debug, Serialize)]
pub struct CohortCounts {
    pub samples: usize,
    pub resistant: usize,
    pub lineages: usize,
    /// k-mer-identical copy groups inside the cohort.
    pub copy_groups: usize,
    pub min_support: u32,
    pub features: usize,
}

/// What `grm learn --exp-cohort` writes for one cohort at one seed. Counts
/// only: no genome id and no row per genome.
#[derive(Clone, Debug, Serialize)]
pub struct CohortResult {
    /// `FORMAT`.
    pub format: String,
    pub evaluation_revision: u32,
    /// The dataset's `input_fingerprint` at the release's seed, so that it
    /// equals the release report's at any seed of this result.
    pub release_fingerprint: String,
    pub eval_config: EvalConfig,
    pub id: String,
    pub species_slug: String,
    pub antibiotic: String,
    pub r: usize,
    pub draw: usize,
    pub seed: u64,
    pub cohort_file_sha256: String,
    /// The dataset's cohort rows (the rows `grm learn` evaluates).
    pub dataset: Counts,
    /// Why the cohort file was refused; nothing is fitted then.
    pub refused: Option<String>,
    pub cohort: Option<CohortCounts>,
    pub scored: Counts,
    pub new_lineage: Counts,
    pub copies_left_out: Counts,
    pub estimates: Vec<SplitEstimate>,
    pub truth: Option<Truth>,
    pub seconds: f64,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::learn::eval::EvalConfig;

    fn file(r: usize, ids: &[&str], res: &[bool], lin: &[u32]) -> CohortFile {
        CohortFile {
            format: COHORT_FORMAT.into(),
            dataset: "species__drug".into(),
            r,
            draw: 0,
            genomes: ids.iter().map(|s| s.to_string()).collect(),
            resistant: res.to_vec(),
            lineage: lin.to_vec(),
            study_fold: None,
        }
    }

    /// SYNTHETIC: six dataset rows g0 to g5.
    const IDS: [&str; 6] = ["g0", "g1", "g2", "g3", "g4", "g5"];
    const RES: [bool; 6] = [true, true, true, false, false, false];
    const LIN: [u32; 6] = [0, 1, 2, 0, 1, 2];

    #[test]
    fn a_cohort_file_must_hold_r_of_each_class_and_rows_of_the_dataset() {
        let ok = file(1, &["g3", "g0"], &[false, true], &[0, 0]);
        assert_eq!(
            cohort_rows(&ok, "species__drug", &IDS, &RES, &LIN).unwrap(),
            vec![0, 3]
        );
        // A class short of r is refused, never taken silently.
        let short = file(2, &["g0", "g1", "g3"], &[true, true, false], &[0, 1, 0]);
        let e = cohort_rows(&short, "species__drug", &IDS, &RES, &LIN).unwrap_err();
        assert!(format!("{e:#}").contains("not 2 of each"), "{e:#}");
        let outside = file(1, &["g0", "g9"], &[true, false], &[0, 0]);
        assert!(cohort_rows(&outside, "species__drug", &IDS, &RES, &LIN).is_err());
        let relabelled = file(1, &["g0", "g1"], &[true, false], &[0, 1]);
        assert!(cohort_rows(&relabelled, "species__drug", &IDS, &RES, &LIN).is_err());
        let twice = file(1, &["g0", "g0"], &[true, false], &[0, 0]);
        assert!(cohort_rows(&twice, "species__drug", &IDS, &RES, &LIN).is_err());
        let other = file(1, &["g0", "g3"], &[true, false], &[0, 0]);
        assert!(cohort_rows(&other, "species__other", &IDS, &RES, &LIN).is_err());
    }

    #[test]
    fn the_scored_set_leaves_out_the_cohort_and_its_copies_and_flags_new_lineages() {
        // SYNTHETIC: rows 0 and 3 are the cohort; row 4 is a copy of row 0
        // (same digest); lineage 2 holds no cohort genome.
        let digests = ["a", "b", "c", "d", "a", "e"];
        let s = scored_set(&[0, 3], &LIN, &digests);
        assert_eq!(s.rows, vec![1, 2, 5]);
        assert_eq!(s.copies, vec![4]);
        assert_eq!(s.new_lineage, vec![true, true, true]);
        for i in [0, 3, 4] {
            assert!(!s.rows.contains(&i), "row {i} must not be scored");
        }
        let s = scored_set(&[0, 1], &LIN, &["a", "b", "c", "d", "e", "f"]);
        assert_eq!(s.rows, vec![2, 3, 4, 5]);
        assert_eq!(s.new_lineage, vec![true, false, false, true]);
    }

    #[test]
    fn a_feature_on_a_genome_outside_the_cohort_reads_its_lowest_pattern() {
        // SYNTHETIC species matrix: 3 patterns over 4 genomes.
        let mut species = BitMatrix::with_cols(4);
        for row in [[1u64 | 4], [2u64 | 8], [1u64 | 2]] {
            species.push_row(&row);
        }
        let x = feature_values(&species, &[0, 2], &[1, 2, 0]);
        assert_eq!((x.rows, x.cols), (2, 3));
        assert_eq!(
            (0..3).map(|j| x.get(0, j)).collect::<Vec<_>>(),
            [false, true, true]
        );
        assert_eq!(
            (0..3).map(|j| x.get(1, j)).collect::<Vec<_>>(),
            [true, false, true]
        );
    }

    #[test]
    fn calls_count_by_class_at_the_cutoff_and_clip_the_log_loss() {
        let c = calls(&[true, true, false, false], &[0.5, 0.0, 0.49, 1.0]);
        assert_eq!((c.n, c.tp, c.fn_, c.tn, c.fp), (4, 1, 1, 1, 1));
        let want = -(0.5f64.ln()) - LOSS_CLIP.ln() - (0.51f64).ln() - LOSS_CLIP.ln();
        assert!((c.loss_sum - want).abs() < 1e-12);
    }

    #[test]
    fn cohort_names_parse_back() {
        let n = cohort_name("escherichia_coli__ampicillin", 25, 3);
        assert_eq!(
            parse_cohort_name(&n),
            Some(("escherichia_coli__ampicillin".into(), 25, 3))
        );
        assert_eq!(parse_cohort_name("x.r1.dz.cohort.json"), None);
        assert_eq!(parse_cohort_name("x.json"), None);
        // A result's name is no cohort file's, so one folder can hold both.
        let result = result_path(Path::new("out"), "escherichia_coli__ampicillin", 25, 3, 42);
        let name = result.file_name().unwrap().to_str().unwrap();
        assert_eq!(name, "escherichia_coli__ampicillin.r25.d3.s42.cohort.json");
        assert_eq!(parse_cohort_name(name), None);
    }

    /// SYNTHETIC cohort: 24 genomes in 12 lineages of two (one R, one S
    /// each), a feature that is the label with three flips, noise features,
    /// and a scored set of 16 genomes, 8 of them in new lineages.
    fn synthetic(seed: u64) -> (BitMatrix, Bits, Vec<u32>, BitMatrix, Vec<bool>, Vec<bool>) {
        let n = 24;
        let mut state = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        let mut next = move || {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (state >> 33) as u32
        };
        let y_vec: Vec<bool> = (0..n).map(|i| i % 2 == 0).collect();
        let groups: Vec<u32> = (0..n as u32).map(|i| i / 2).collect();
        let mut x = BitMatrix::with_cols(n);
        let mut sx = BitMatrix::with_cols(16);
        let sy: Vec<bool> = (0..16).map(|j| j % 2 == 0).collect();
        let mut rows: Vec<(Vec<bool>, Vec<bool>)> = Vec::new();
        let flip = |i: usize| i == 3 || i == 10 || i == 17;
        rows.push((
            (0..n).map(|i| y_vec[i] ^ flip(i)).collect(),
            (0..16).map(|j| sy[j] ^ (j == 5)).collect(),
        ));
        for _ in 0..30 {
            rows.push((
                (0..n).map(|_| next() % 3 == 0).collect(),
                (0..16).map(|_| next() % 3 == 0).collect(),
            ));
        }
        let mut seen = HashSet::new();
        for (a, b) in rows {
            let support = a.iter().filter(|&&v| v).count();
            if !(3..=n - 3).contains(&support) || !seen.insert(a.clone()) {
                continue;
            }
            let mut wa = vec![0u64; bits::words_for(n)];
            a.iter()
                .enumerate()
                .filter(|(_, v)| **v)
                .for_each(|(i, _)| bits::set(&mut wa, i));
            let mut wb = vec![0u64; bits::words_for(16)];
            b.iter()
                .enumerate()
                .filter(|(_, v)| **v)
                .for_each(|(i, _)| bits::set(&mut wb, i));
            x.push_row(&wa);
            sx.push_row(&wb);
        }
        let y = Bits::from_indices(n, (0..n).filter(|&i| y_vec[i]));
        let new: Vec<bool> = (0..16).map(|j| j >= 8).collect();
        (x, y, groups, sx, sy, new)
    }

    fn inputs<'a>(
        s: &'a (BitMatrix, Bits, Vec<u32>, BitMatrix, Vec<bool>, Vec<bool>),
        roles: Option<&'a [Role]>,
    ) -> CohortInputs<'a> {
        CohortInputs {
            x: &s.0,
            y: &s.1,
            groups: &s.2,
            duplicates: &[],
            min_support: 3,
            study_roles: roles,
            scored_x: &s.3,
            scored_y: &s.4,
            scored_new: &s.5,
        }
    }

    #[test]
    fn the_cohort_lineage_split_keeps_each_lineage_in_one_fold() {
        let s = synthetic(7);
        let fit = fit_cohort("species__drug", &inputs(&s, None), &EvalConfig::default());
        let lineage = fit
            .estimates
            .iter()
            .find(|e| e.scheme == "lineage")
            .unwrap();
        assert_eq!(lineage.unavailable, None);
        assert_eq!(lineage.n_groups, 12);
        assert_eq!(fit.lineage_test_fold.len(), 24);
        let mut fold_of: HashMap<u32, usize> = HashMap::new();
        for (i, &f) in fit.lineage_test_fold.iter().enumerate() {
            let at = *fold_of.entry(s.2[i]).or_insert(f);
            assert_eq!(at, f, "lineage {} lies in two folds", s.2[i]);
        }
        // Every learner's calls cover the whole cohort and agree with the
        // split's own metrics.
        for row in &lineage.learners {
            assert_eq!(row.calls.n, 24);
            let ba = (row.calls.tp as f64 / 12.0 + row.calls.tn as f64 / 12.0) / 2.0;
            assert!(
                (ba - row.balanced_accuracy).abs() < 1e-12,
                "{}",
                row.learner
            );
            assert_eq!(row.fold_balanced_accuracy.len(), 5);
        }
    }

    #[test]
    fn the_truth_is_the_final_models_on_the_scored_set() {
        let s = synthetic(11);
        let cfg = EvalConfig::default();
        let fit = fit_cohort("species__drug", &inputs(&s, None), &cfg);
        assert_eq!(fit.truth.unavailable, None);
        // The final models alone, applied to the scored matrix.
        let (finals, direct) =
            eval::final_scores("species__drug", &s.0, &s.1, &s.2, &[], 3, &cfg, &s.3).unwrap();
        let new: Vec<usize> = (0..16).filter(|&j| s.5[j]).collect();
        for (row, scores) in fit.truth.learners.iter().zip(&direct) {
            let scores = scores.clone().unwrap_or_else(|| vec![0.5; 16]);
            assert_eq!(row.all, calls(&s.4, &scores), "{}", row.learner);
            let ny: Vec<bool> = new.iter().map(|&j| s.4[j]).collect();
            let ns: Vec<f64> = new.iter().map(|&j| scores[j]).collect();
            assert_eq!(row.new_lineage, calls(&ny, &ns), "{}", row.learner);
            assert_eq!(row.all.n, 16);
            assert_eq!(row.new_lineage.n, 8);
        }
        // The final models are those `evaluate` reports, so the truth is no
        // fold's model: on the cohort itself the final scm and cart call as
        // the reported final models do.
        let out = eval::evaluate("species__drug", &s.0, &s.1, &s.2, &[], 3, None, &cfg).unwrap();
        assert_eq!(
            serde_json::to_string(&finals).unwrap(),
            serde_json::to_string(&out.finals).unwrap()
        );
        let (_, own) =
            eval::final_scores("species__drug", &s.0, &s.1, &s.2, &[], 3, &cfg, &s.0).unwrap();
        let rules = out.finals.scm.rules.clone();
        let scm = crate::learn::scm::ScmModel {
            kind: out.finals.scm.kind,
            p: out.finals.scm.p,
            rules,
            steps: vec![],
        }
        .predict(&s.0);
        let want: Vec<f64> = (0..24).map(|i| scm.contains(i) as u8 as f64).collect();
        assert_eq!(own[0].as_ref().unwrap(), &want);
        assert_eq!(
            own[1].as_ref().unwrap(),
            &out.finals.cart.predict(&s.0, out.finals.cart_depth)
        );
        // The truth's lasso is the reported final lasso on the scored matrix.
        let lasso = out
            .finals
            .lasso
            .as_ref()
            .filter(|m| !m.coef.is_empty())
            .expect("the synthetic cohort's final lasso keeps a feature");
        assert_eq!(direct[4].as_ref().unwrap(), &lasso.predict_proba(&s.3));
    }

    /// A feature the forest's screen drops cannot move its truth: flipping
    /// every unscreened feature of every scored genome leaves the forest's
    /// scores as they were, so its trees read the screened rows.
    #[test]
    fn the_truths_forest_reads_a_scored_genome_on_the_screened_features() {
        let s = synthetic(11);
        let cfg = EvalConfig {
            screened_features: 4,
            ..EvalConfig::default()
        };
        let feats = crate::learn::screen(&s.0, &s.1, &Bits::ones(24), 4);
        assert_eq!(feats.len(), 4);
        // A screened position whose row is not screened, so a forest that
        // read the scored matrix unscreened would read a flipped row.
        assert!((0..4u32).any(|i| !feats.contains(&i)));
        let mut flipped = BitMatrix::with_cols(16);
        for f in 0..s.3.rows {
            let src = s.3.row(f);
            let mut row = vec![0u64; bits::words_for(16)];
            for j in 0..16 {
                if bits::get(src, j) != !feats.contains(&(f as u32)) {
                    bits::set(&mut row, j);
                }
            }
            flipped.push_row(&row);
        }
        let (_, a) =
            eval::final_scores("species__drug", &s.0, &s.1, &s.2, &[], 3, &cfg, &s.3).unwrap();
        let (_, b) =
            eval::final_scores("species__drug", &s.0, &s.1, &s.2, &[], 3, &cfg, &flipped).unwrap();
        assert_eq!(a[2], b[2]);
    }

    #[test]
    fn a_scored_genome_reads_a_features_lowest_species_pattern() {
        let mut x = BitMatrix::with_cols(2);
        for _ in 0..3 {
            x.push_row(&[0b01]);
        }
        let ds = crate::matrix::Dataset {
            x,
            y: Bits::from_indices(2, [0]),
            columns: vec![0, 1],
            feature_offsets: vec![0, 2, 5, 6],
            feature_patterns: vec![7, 3, 9, 4, 6, 2],
        };
        assert_eq!(lowest_patterns(&ds), vec![3, 4, 2]);
    }

    #[test]
    fn a_null_final_fit_scores_at_the_cohorts_share_not_the_scored_sets() {
        let s = synthetic(11);
        let mut c = inputs(&s, None);
        let sy: Vec<bool> = (0..16).map(|j| j < 5).collect();
        c.scored_y = &sy;
        assert_eq!(null_fit_score(&c), 0.5);
    }

    #[test]
    fn a_scored_set_short_of_a_class_has_no_truth() {
        let mut s = synthetic(5);
        s.4 = vec![true; 16];
        let fit = fit_cohort("species__drug", &inputs(&s, None), &EvalConfig::default());
        assert_eq!(
            fit.truth.unavailable.as_deref(),
            Some("the scored genomes lack a class")
        );
    }

    #[test]
    fn the_study_estimate_scores_the_plans_folds_or_says_why_not() {
        let s = synthetic(3);
        // Each lineage pair (one R, one S) in fold lineage mod 5.
        let roles: Vec<Role> = s.2.iter().map(|&g| Role::Test(g as usize % 5)).collect();
        let fit = fit_cohort(
            "species__drug",
            &inputs(&s, Some(&roles)),
            &EvalConfig::default(),
        );
        let study = fit
            .estimates
            .iter()
            .find(|e| e.scheme == STUDY_SCHEME)
            .unwrap();
        assert_eq!(study.unavailable, None);
        assert_eq!(study.learners.len(), 5);
        assert!(study.learners.iter().all(|l| l.calls.n == 24));
        // A fold with one class only is refused, and the result says why:
        // fold 0's susceptible genomes move to fold 1.
        let bad: Vec<Role> = roles
            .iter()
            .enumerate()
            .map(|(i, &r)| {
                if r == Role::Test(0) && !s.1.contains(i) {
                    Role::Test(1)
                } else {
                    r
                }
            })
            .collect();
        let fit = fit_cohort(
            "species__drug",
            &inputs(&s, Some(&bad)),
            &EvalConfig::default(),
        );
        let study = fit
            .estimates
            .iter()
            .find(|e| e.scheme == STUDY_SCHEME)
            .unwrap();
        assert!(
            study.unavailable.as_deref().unwrap().contains("lack"),
            "{:?}",
            study.unavailable
        );
        assert!(study.learners.is_empty());
    }
}
