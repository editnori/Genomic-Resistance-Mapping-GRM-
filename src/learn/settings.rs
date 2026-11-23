//! Experiment, off by default (plan block C5, docs/settings-arms.md): fit
//! each learner at settings a file names, on a fold file's partition, and
//! write each one's inner scores and calls; and write the inner folds a
//! nested search of those settings runs on.
//!
//! Nothing here changes a report, a fold file result or a fit of a run
//! without its flags. `grm learn --fold-file F --exp-settings S` loads the
//! matrix once and fits each trial of S on F's folds, one learner at a time,
//! with `eval::evaluate_assigned_with`: the fold file's partition checks,
//! candidates, support floor and pooling, with the learner's fit at the
//! trial's setting (`fit_setting`), every draw seeded from S's `seed` as a
//! fold file's fold j is seeded from 42, and tuning that keeps S's `tuning`
//! groups whole. It writes `<fold-out>/<dataset>.<scheme>.t<trial>.json`
//! per trial (`TrialResult`): per learner its balanced accuracy and size per
//! fold and its calls per fold at 0.5 and at each cutoff S names, or the
//! reason its fit failed. A trial result holds no genome id and no score.
//!
//! `grm learn --exp-settings-folds DIR`, beside a drawn run, writes the
//! lineage split's outer folds and, inside each outer fold's training
//! genomes T_f, the inner folds `training_folds(y, T_f, lineage groups, 3,
//! fold_seed(seed, f))`, the partition the release's scm and cart tune on in
//! that fold, so no inner fold splits a lineage. Beside `--fold-file` it
//! writes the file's folds and the same inner folds inside them. These
//! files name genomes, as fold files do, and stay where they are made.
//!
//! A setting of `release` is the release's own fit (`eval::fit_seeded`),
//! with its own inner searches. Any other setting fixes every value of its
//! learner's space in `scripts/hyperopt/protocol.py` and is fitted as given:
//! scm and cart with no inner search, the forest and boosting with their
//! fixed values, boosting's `rounds` as the cap its early stopping runs
//! under, and the lasso's path of `n_lambdas` points with its holdout
//! choosing the penalty. The forest, boosting and the lasso screen their own
//! fitting genomes to `screened_features`, as the release does.

use super::eval::{
    self, AssignedPlan, EvalConfig, FitFn, FoldModel, FoldScm, LEARNERS, LassoContext, Role,
    Scored, SplitResult,
};
use super::scm::{self, Kind};
use super::training_folds;
use super::tree::{self, ForestParams, TreeParams};
use crate::bits::{BitMatrix, Bits};
use crate::cutoff::{Calls, calls_at};
use crate::learn::gbt::GbtParams;
use anyhow::{Context, Result, bail, ensure};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::sync::OnceLock;
use std::time::Instant;

/// The settings file's format.
pub const SETTINGS_FORMAT: &str = "grm-c5-settings-1";
/// A trial result's format. A change to what a result holds changes it.
pub const TRIAL_FORMAT: &str = "grm-c5-trial-1";
/// The folds file's format.
pub const FOLDS_FORMAT: &str = "grm-c5-folds-1";

/// Which groups the fits' own tuning (scm and cart's inner folds, boosting's
/// and the lasso's holdouts) keeps whole.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Deserialize, Serialize)]
#[serde(rename_all = "lowercase")]
pub enum Tuning {
    /// The lineage groups, as the lineage split's fits tune.
    Lineage,
    /// The k-mer-identical copy groups, as a fold file's fits tune.
    Copies,
}

/// One trial of a settings file: a setting per learner it fits, and the
/// cutoff each learner's calls are also counted at.
#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Trial {
    pub id: u32,
    pub learners: BTreeMap<String, serde_json::Value>,
    #[serde(default)]
    pub cutoffs: BTreeMap<String, f64>,
}

/// What `--exp-settings` reads.
#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SettingsFile {
    pub format: String,
    /// The base seed of every fit, in place of the release's 42.
    pub seed: u64,
    pub tuning: Tuning,
    /// Cutoffs every learner's calls are counted at, beside 0.5.
    #[serde(default)]
    pub cutoffs: Vec<f64>,
    pub trials: Vec<Trial>,
}

/// One learner's setting, checked against its space.
#[derive(Clone, Debug, PartialEq)]
pub enum Setting {
    /// The release's own fit and inner searches.
    Release,
    Scm {
        kind: Kind,
        p: f64,
        max_rules: usize,
    },
    Cart {
        max_depth: u32,
        min_leaf: u32,
    },
    Rf {
        n_trees: usize,
        subsample: f64,
        max_depth: u32,
        min_leaf: u32,
    },
    Gbt {
        rounds: usize,
        learning_rate: f64,
        max_depth: u32,
        lambda: f64,
        min_child_weight: f64,
        colsample: f64,
    },
    Lasso {
        n_lambdas: usize,
    },
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct ScmIn {
    kind: Kind,
    p: f64,
    max_rules: usize,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct CartIn {
    max_depth: u32,
    min_leaf: u32,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct RfIn {
    n_trees: usize,
    subsample: f64,
    max_depth: u32,
    min_leaf: u32,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct GbtIn {
    rounds: usize,
    learning_rate: f64,
    max_depth: u32,
    lambda: f64,
    min_child_weight: f64,
    colsample: f64,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct LassoIn {
    n_lambdas: usize,
}

fn within<T: PartialOrd + std::fmt::Display + Copy>(name: &str, v: T, lo: T, hi: T) -> Result<()> {
    ensure!(
        v >= lo && v <= hi,
        "{name} {v} lies outside its space, {lo} to {hi}"
    );
    Ok(())
}

/// `learner`'s setting from `value`: the text `release`, or every value of
/// the learner's space in `scripts/hyperopt/protocol.py`, each inside it.
pub fn parse_setting(learner: &str, value: &serde_json::Value) -> Result<Setting> {
    if value.as_str() == Some("release") {
        return Ok(Setting::Release);
    }
    let v = value.clone();
    let parsed = match learner {
        "scm" => {
            let s: ScmIn = serde_json::from_value(v)?;
            within("scm p", s.p, 0.1, 10.0)?;
            within("scm max_rules", s.max_rules, 1, 10)?;
            Setting::Scm {
                kind: s.kind,
                p: s.p,
                max_rules: s.max_rules,
            }
        }
        "cart" => {
            let s: CartIn = serde_json::from_value(v)?;
            within("cart max_depth", s.max_depth, 1, 30)?;
            within("cart min_leaf", s.min_leaf, 1, 20)?;
            Setting::Cart {
                max_depth: s.max_depth,
                min_leaf: s.min_leaf,
            }
        }
        "rf" => {
            let s: RfIn = serde_json::from_value(v)?;
            within("rf n_trees", s.n_trees, 50, 500)?;
            within("rf subsample", s.subsample, 0.3, 1.0)?;
            within("rf max_depth", s.max_depth, 2, 30)?;
            within("rf min_leaf", s.min_leaf, 1, 10)?;
            Setting::Rf {
                n_trees: s.n_trees,
                subsample: s.subsample,
                max_depth: s.max_depth,
                min_leaf: s.min_leaf,
            }
        }
        "gbt" => {
            let s: GbtIn = serde_json::from_value(v)?;
            within("gbt rounds", s.rounds, 50, 800)?;
            within("gbt learning_rate", s.learning_rate, 0.01, 0.3)?;
            within("gbt max_depth", s.max_depth, 1, 8)?;
            within("gbt lambda", s.lambda, 0.1, 10.0)?;
            within("gbt min_child_weight", s.min_child_weight, 0.1, 10.0)?;
            within("gbt colsample", s.colsample, 0.1, 1.0)?;
            Setting::Gbt {
                rounds: s.rounds,
                learning_rate: s.learning_rate,
                max_depth: s.max_depth,
                lambda: s.lambda,
                min_child_weight: s.min_child_weight,
                colsample: s.colsample,
            }
        }
        "lasso" => {
            let s: LassoIn = serde_json::from_value(v)?;
            within("lasso n_lambdas", s.n_lambdas, 10, 100)?;
            Setting::Lasso {
                n_lambdas: s.n_lambdas,
            }
        }
        other => bail!("{other:?} is not a learner"),
    };
    Ok(parsed)
}

/// A settings file, checked, with each setting parsed.
#[derive(Clone, Debug)]
pub struct Loaded {
    pub file: SettingsFile,
    /// Per trial, its learners in `LEARNERS` order with their settings.
    pub settings: Vec<Vec<(&'static str, Setting)>>,
    pub sha256: String,
}

fn good_cutoff(c: f64) -> bool {
    c.is_finite() && c > 0.0 && c < 1.0
}

/// Reads and checks a settings file: its format, at least one trial, trial
/// ids each once, each learner a learner of the release with a setting
/// inside its space, and each cutoff inside (0, 1) and of a learner the
/// trial fits.
pub fn parse(bytes: &[u8]) -> Result<Loaded> {
    let file: SettingsFile = serde_json::from_slice(bytes).context("settings file")?;
    ensure!(
        file.format == SETTINGS_FORMAT,
        "settings file format {:?}, not {SETTINGS_FORMAT}",
        file.format
    );
    ensure!(!file.trials.is_empty(), "the settings file has no trial");
    ensure!(
        file.cutoffs.iter().all(|&c| good_cutoff(c)),
        "a cutoff lies outside (0, 1)"
    );
    let mut ids = std::collections::BTreeSet::new();
    let mut settings = Vec::with_capacity(file.trials.len());
    for trial in &file.trials {
        ensure!(ids.insert(trial.id), "trial {} is listed twice", trial.id);
        ensure!(
            !trial.learners.is_empty(),
            "trial {} fits no learner",
            trial.id
        );
        for name in trial.learners.keys() {
            ensure!(
                LEARNERS.contains(&name.as_str()),
                "trial {}: {name:?} is not a learner",
                trial.id
            );
        }
        for (name, &c) in &trial.cutoffs {
            ensure!(
                trial.learners.contains_key(name),
                "trial {}: a cutoff for {name:?}, which it does not fit",
                trial.id
            );
            ensure!(
                good_cutoff(c),
                "trial {}: {name}'s cutoff {c} lies outside (0, 1)",
                trial.id
            );
        }
        let mut row = Vec::new();
        for learner in LEARNERS {
            if let Some(value) = trial.learners.get(learner) {
                let s = parse_setting(learner, value)
                    .with_context(|| format!("trial {} {learner}", trial.id))?;
                row.push((learner, s));
            }
        }
        settings.push(row);
    }
    Ok(Loaded {
        file,
        settings,
        sha256: crate::research::sha256(bytes),
    })
}

static TRIALS: OnceLock<Loaded> = OnceLock::new();
static FOLDS: OnceLock<PathBuf> = OnceLock::new();

/// Turns the trials on for this process from the settings file at `path`.
pub fn enable_trials(path: &Path) -> Result<()> {
    let bytes = std::fs::read(path).with_context(|| format!("read {}", path.display()))?;
    let loaded = parse(&bytes).with_context(|| format!("{}", path.display()))?;
    TRIALS
        .set(loaded)
        .map_err(|_| anyhow::anyhow!("the C5 settings file is already set"))
}

/// The settings file `enable_trials` loaded, when on.
pub fn trials() -> Option<&'static Loaded> {
    TRIALS.get()
}

/// Turns the folds writer on for this process; `dir` is created if missing.
pub fn enable_folds(dir: &Path) -> Result<()> {
    std::fs::create_dir_all(dir).with_context(|| format!("create {}", dir.display()))?;
    FOLDS
        .set(dir.to_path_buf())
        .map_err(|_| anyhow::anyhow!("the C5 folds folder is already set"))
}

/// The folder `enable_folds` set, when on.
pub fn folds_dir() -> Option<&'static Path> {
    FOLDS.get().map(PathBuf::as_path)
}

/// One learner's fit at `setting` on one outer fold, seeded as `fit_fold`
/// seeds it: every partition and draw from `fold_seed(cfg.seed, fold)`.
#[allow(clippy::too_many_arguments)]
pub(crate) fn fit_setting(
    setting: &Setting,
    learner: &str,
    x: &BitMatrix,
    y: &Bits,
    train: &Bits,
    groups: Option<&[u32]>,
    cfg: &EvalConfig,
    context: LassoContext<'_>,
) -> Result<(Scored, FoldModel)> {
    let fold = context.outer_fold.expect("an outer fit needs its fold");
    let seed = crate::ablate::plan::fold_seed(cfg.seed, fold);
    let with =
        |c: EvalConfig| eval::fit_seeded(learner, x, y, train, groups, &c, seed, seed, context);
    match (learner, setting) {
        (_, Setting::Release) => {
            eval::fit_seeded(learner, x, y, train, groups, cfg, seed, seed, context)
        }
        ("scm", &Setting::Scm { kind, p, max_rules }) => {
            let m = scm::fit(x, y, train, kind, p, max_rules, false);
            let pred = m.predict(x);
            Ok((
                Scored {
                    scores: Some((0..x.cols).map(|i| pred.contains(i) as u8 as f64).collect()),
                    complexity: m.rules.len() as f64,
                    lasso_fit: None,
                },
                FoldModel::Scm(FoldScm {
                    kind: m.kind,
                    p: m.p,
                    rules: m.rules,
                }),
            ))
        }
        (
            "cart",
            &Setting::Cart {
                max_depth,
                min_leaf,
            },
        ) => {
            let t = tree::fit_tree(
                x,
                y,
                train,
                &TreeParams {
                    max_depth,
                    min_leaf,
                    max_features: None,
                    seed: 0,
                    tie: cfg.tie_rule,
                },
            );
            let leaves = t.nodes.iter().filter(|n| n.split.is_none()).count();
            Ok((
                Scored {
                    scores: Some(t.predict(x, max_depth)),
                    complexity: leaves as f64,
                    lasso_fit: None,
                },
                FoldModel::Cart(t),
            ))
        }
        (
            "rf",
            &Setting::Rf {
                n_trees,
                subsample,
                max_depth,
                min_leaf,
            },
        ) => with(EvalConfig {
            forest: ForestParams {
                n_trees,
                subsample,
                max_depth,
                min_leaf,
                ..cfg.forest
            },
            ..cfg.clone()
        }),
        (
            "gbt",
            &Setting::Gbt {
                rounds,
                learning_rate,
                max_depth,
                lambda,
                min_child_weight,
                colsample,
            },
        ) => with(EvalConfig {
            gbt: GbtParams {
                rounds,
                learning_rate,
                max_depth,
                lambda,
                min_child_weight,
                colsample,
                ..cfg.gbt
            },
            ..cfg.clone()
        }),
        ("lasso", &Setting::Lasso { n_lambdas }) => with(EvalConfig {
            lasso_lambdas: n_lambdas,
            ..cfg.clone()
        }),
        _ => bail!("{learner} cannot take the setting {setting:?}"),
    }
}

/// One learner's calls on each outer fold at one cutoff.
#[derive(Clone, Debug, Serialize)]
pub struct CutoffCalls {
    pub cutoff: f64,
    pub folds: Vec<Calls>,
}

/// One learner of one trial.
#[derive(Clone, Debug, Serialize)]
pub struct TrialLearner {
    pub learner: String,
    /// The setting as the file gave it.
    pub setting: serde_json::Value,
    /// `ok`, or `failed` with `reason`.
    pub status: String,
    pub reason: Option<String>,
    /// The cutoff the trial gave this learner, if any.
    pub cutoff: Option<f64>,
    pub fold_balanced_accuracy: Vec<f64>,
    pub fold_complexity: Vec<f64>,
    /// Test genomes of null-fit folds, by class, scored at their fold's
    /// training prevalence.
    pub no_call_resistant: usize,
    pub no_call_susceptible: usize,
    /// Calls per fold at 0.5 and at every cutoff named, cutoffs ascending.
    pub calls: Vec<CutoffCalls>,
    pub seconds: f64,
}

/// What one trial writes: counts, settings and reasons, no genome id.
#[derive(Clone, Debug, Serialize)]
pub struct TrialResult {
    pub format: &'static str,
    pub evaluation_revision: u32,
    pub dataset: String,
    pub scheme: String,
    pub fold_file_sha256: String,
    pub settings_sha256: String,
    pub seed: u64,
    pub tuning: Tuning,
    pub trial: u32,
    pub min_support: u32,
    pub folds: Vec<eval::FoldCounts>,
    pub learners: Vec<TrialLearner>,
}

/// Each fold's calls of a learner's pooled scores at each of `cutoffs`.
fn fold_calls(
    truth: &[bool],
    test_fold: &[usize],
    scores: &[f64],
    folds: usize,
    cutoffs: &[f64],
) -> Vec<CutoffCalls> {
    cutoffs
        .iter()
        .map(|&cutoff| CutoffCalls {
            cutoff,
            folds: (0..folds)
                .map(|f| {
                    let (t, s): (Vec<bool>, Vec<f64>) = (0..truth.len())
                        .filter(|&p| test_fold[p] == f)
                        .map(|p| (truth[p], scores[p]))
                        .unzip();
                    calls_at(&t, &s, cutoff)
                })
                .collect(),
        })
        .collect()
}

/// The trials of `loaded`, each learner fitted on the supplied partition
/// apart, so one learner's failure is recorded and the others go on. A
/// failure of the partition itself fails every learner with its reason.
#[allow(clippy::too_many_arguments)]
pub fn run_trials(
    dataset: &str,
    scheme: &str,
    x: &BitMatrix,
    y: &Bits,
    groups: &[u32],
    duplicates: &[Vec<usize>],
    roles: &[Role],
    folds: usize,
    base: &EvalConfig,
    loaded: &Loaded,
    fold_file_sha256: &str,
) -> Result<Vec<TrialResult>> {
    let file = &loaded.file;
    let cfg = EvalConfig {
        seed: file.seed,
        ..base.clone()
    };
    let tuning_scheme = match file.tuning {
        Tuning::Lineage => "lineage",
        Tuning::Copies => scheme,
    };
    let mut out = Vec::with_capacity(file.trials.len());
    for (trial, settings) in file.trials.iter().zip(&loaded.settings) {
        let mut learners = Vec::with_capacity(settings.len());
        let mut min_support = eval::supplied_support_floor(roles);
        let mut fold_counts = Vec::new();
        for (learner, setting) in settings {
            let given = trial.cutoffs.get(*learner).copied();
            let mut cutoffs: Vec<f64> = std::iter::once(crate::predict::CALL_CUTOFF)
                .chain(file.cutoffs.iter().copied())
                .chain(given)
                .collect();
            cutoffs.sort_by(f64::total_cmp);
            cutoffs.dedup();
            let fit: FitFn<'_> =
                &|l, x, y, train, g, c, ctx| fit_setting(setting, l, x, y, train, g, c, ctx);
            let names = [*learner];
            let plan = AssignedPlan {
                learners: &names,
                fit,
                tuning_scheme,
            };
            let t0 = Instant::now();
            let result = eval::evaluate_assigned_with(
                dataset,
                scheme,
                x,
                y,
                groups,
                duplicates,
                roles,
                folds,
                &cfg,
                Some(&plan),
            );
            let seconds = t0.elapsed().as_secs_f64();
            let value = trial.learners[*learner].clone();
            learners.push(match result {
                Ok(a) => {
                    min_support = a.min_support;
                    fold_counts = a.folds;
                    let row = &a.result.learners[0];
                    let truth: Vec<bool> = a.tested.iter().map(|&i| y.contains(i)).collect();
                    let (no_r, no_s) = row.no_call.map_or((0, 0), |c| (c.resistant, c.susceptible));
                    TrialLearner {
                        learner: learner.to_string(),
                        setting: value,
                        status: "ok".into(),
                        reason: None,
                        cutoff: given,
                        fold_balanced_accuracy: row.fold_balanced_accuracy.clone(),
                        fold_complexity: row.fold_complexity.clone(),
                        no_call_resistant: no_r,
                        no_call_susceptible: no_s,
                        calls: fold_calls(
                            &truth,
                            &a.result.test_fold,
                            &row.scores,
                            folds,
                            &cutoffs,
                        ),
                        seconds,
                    }
                }
                Err(e) => TrialLearner {
                    learner: learner.to_string(),
                    setting: value,
                    status: "failed".into(),
                    reason: Some(format!("{e:#}")),
                    cutoff: given,
                    fold_balanced_accuracy: vec![],
                    fold_complexity: vec![],
                    no_call_resistant: 0,
                    no_call_susceptible: 0,
                    calls: vec![],
                    seconds,
                },
            });
        }
        out.push(TrialResult {
            format: TRIAL_FORMAT,
            evaluation_revision: crate::analysis::EVALUATION_REVISION,
            dataset: dataset.to_string(),
            scheme: scheme.to_string(),
            fold_file_sha256: fold_file_sha256.to_string(),
            settings_sha256: loaded.sha256.clone(),
            seed: file.seed,
            tuning: file.tuning,
            trial: trial.id,
            min_support,
            folds: fold_counts,
            learners,
        });
    }
    Ok(out)
}

/// `<dir>/<dataset>.<scheme>.t<trial>.json`, where `run_trials`' results go.
pub fn trial_path(dir: &Path, dataset: &str, scheme: &str, trial: u32) -> PathBuf {
    dir.join(format!("{dataset}.{scheme}.t{trial}.json"))
}

/// The outer and inner folds of one split, as `--exp-settings-folds`
/// writes them. It names genomes and stays where it is made.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct FoldsFile {
    pub format: String,
    pub dataset: String,
    pub scheme: String,
    pub seed: u64,
    pub folds: usize,
    pub inner_folds: usize,
    pub unavailable_reason: Option<String>,
    pub genomes: Vec<String>,
    pub resistant: Vec<bool>,
    /// Each genome's lineage number.
    pub lineage: Vec<u32>,
    /// Each genome's outer test fold; null when it trains in every fold or
    /// is unused.
    pub outer: Vec<Option<usize>>,
    /// Whether the split tests or trains the genome.
    pub used: Vec<bool>,
    /// Per outer fold f, each genome's inner fold inside f's training
    /// genomes; null for a genome outside them.
    pub inner: Vec<Vec<Option<usize>>>,
}

/// The folds file of a split: outer folds `outer` over the `used` genomes,
/// and inside each outer fold's training genomes the inner folds grouped by
/// lineage, drawn as the release draws the inner folds of that outer fold
/// at `seed`.
#[allow(clippy::too_many_arguments)]
pub fn folds_file(
    dataset: &str,
    scheme: &str,
    seed: u64,
    folds: usize,
    inner_folds: usize,
    genomes: &[String],
    y: &Bits,
    lineage: &[u32],
    outer: Vec<Option<usize>>,
    used: Vec<bool>,
) -> FoldsFile {
    let n = genomes.len();
    let inner = (0..folds)
        .map(|f| {
            let train = Bits::from_indices(n, (0..n).filter(|&i| used[i] && outer[i] != Some(f)));
            let seed_f = crate::ablate::plan::fold_seed(seed, f);
            training_folds(y, &train, Some(lineage), inner_folds, seed_f)
                .into_iter()
                .map(|j| (j != usize::MAX).then_some(j))
                .collect()
        })
        .collect();
    FoldsFile {
        format: FOLDS_FORMAT.into(),
        dataset: dataset.into(),
        scheme: scheme.into(),
        seed,
        folds,
        inner_folds,
        unavailable_reason: None,
        genomes: genomes.to_vec(),
        resistant: (0..n).map(|i| y.contains(i)).collect(),
        lineage: lineage.to_vec(),
        outer,
        used,
        inner,
    }
}

/// `<dir>/<dataset>.<scheme>-s<seed>.folds.json`.
pub fn folds_path(dir: &Path, dataset: &str, scheme: &str, seed: u64) -> PathBuf {
    dir.join(format!("{dataset}.{scheme}-s{seed}.folds.json"))
}

fn write_json(path: &Path, value: &impl Serialize) -> Result<()> {
    let tmp = path.with_extension("json.tmp");
    std::fs::write(&tmp, serde_json::to_vec(value)?)
        .with_context(|| format!("write {}", tmp.display()))?;
    std::fs::rename(&tmp, path).with_context(|| format!("replace {}", path.display()))
}

/// `analyze_dataset`'s hook: the lineage split of a drawn run.
pub fn write_drawn_folds(
    dir: &Path,
    dataset: &str,
    genomes: &[String],
    y: &Bits,
    lineage: &[u32],
    splits: &[SplitResult],
    cfg: &EvalConfig,
) -> Result<PathBuf> {
    let split = splits
        .iter()
        .find(|s| s.scheme == "lineage")
        .context("evaluate wrote no lineage split")?;
    let n = genomes.len();
    let mut file = if split.unavailable_reason.is_none() {
        ensure!(
            split.test_fold.len() == n,
            "the lineage split scores {} of {n} genomes",
            split.test_fold.len()
        );
        folds_file(
            dataset,
            "lineage",
            cfg.seed,
            split.folds,
            cfg.inner_folds,
            genomes,
            y,
            lineage,
            split.test_fold.iter().map(|&f| Some(f)).collect(),
            vec![true; n],
        )
    } else {
        folds_file(
            dataset,
            "lineage",
            cfg.seed,
            0,
            cfg.inner_folds,
            genomes,
            y,
            lineage,
            vec![None; n],
            vec![false; n],
        )
    };
    file.unavailable_reason = split.unavailable_reason.clone();
    let path = folds_path(dir, dataset, "lineage", cfg.seed);
    write_json(&path, &file)?;
    Ok(path)
}

/// `evaluate_fold_file`'s hook: a fold file's split.
#[allow(clippy::too_many_arguments)]
pub fn write_assigned_folds(
    dir: &Path,
    dataset: &str,
    scheme: &str,
    genomes: &[String],
    y: &Bits,
    lineage: &[u32],
    roles: &[Role],
    folds: usize,
    cfg: &EvalConfig,
) -> Result<PathBuf> {
    let outer = roles
        .iter()
        .map(|r| match r {
            Role::Test(f) => Some(*f),
            _ => None,
        })
        .collect();
    let used = roles.iter().map(|r| *r != Role::Unused).collect();
    let file = folds_file(
        dataset,
        scheme,
        cfg.seed,
        folds,
        cfg.inner_folds,
        genomes,
        y,
        lineage,
        outer,
        used,
    );
    let path = folds_path(dir, dataset, scheme, cfg.seed);
    write_json(&path, &file)?;
    Ok(path)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::learn::eval::{evaluate, evaluate_assigned};

    /// 96 genomes in 24 lineages of 4; 12 features, two of which carry the
    /// label with some noise, the rest drawn from a fixed generator.
    fn fixture() -> (BitMatrix, Bits, Vec<u32>) {
        let n = 96;
        let mut state = 0x2545_f491_4f6c_dd1du64;
        let mut next = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        let y = Bits::from_indices(n, (0..n).filter(|i| (i / 4) % 5 < 2));
        let mut x = BitMatrix::with_cols(n);
        for f in 0..12 {
            let on: Vec<usize> = (0..n)
                .filter(|&i| match f {
                    0 => y.contains(i) != (i % 11 == 0),
                    1 => y.contains(i) != (i % 7 == 3),
                    _ => next() % 3 == 0,
                })
                .collect();
            x.push_row(Bits::from_indices(n, on).words());
        }
        (x, y, (0..n).map(|i| (i / 4) as u32).collect())
    }

    fn settings_json(seed: u64, tuning: &str, trials: &str) -> Vec<u8> {
        format!(
            r#"{{"format":"{SETTINGS_FORMAT}","seed":{seed},"tuning":"{tuning}","cutoffs":[0.3,0.7],"trials":{trials}}}"#
        )
        .into_bytes()
    }

    const RELEASE_ALL: &str = r#"[{"id":0,"learners":{"scm":"release","cart":"release","rf":"release","gbt":"release","lasso":"release"},"cutoffs":{"rf":0.35}}]"#;

    const SAMPLED: &str = r#"[{"id":7,"learners":{
        "scm":{"kind":"disjunction","p":2.5,"max_rules":2},
        "cart":{"max_depth":1,"min_leaf":3},
        "rf":{"n_trees":50,"subsample":0.5,"max_depth":3,"min_leaf":2},
        "gbt":{"rounds":60,"learning_rate":0.2,"max_depth":2,"lambda":1.5,"min_child_weight":0.5,"colsample":0.8},
        "lasso":{"n_lambdas":12}}}]"#;

    fn roles_from(folds: &[usize]) -> Vec<Role> {
        folds.iter().map(|&f| Role::Test(f)).collect()
    }

    fn run(
        x: &BitMatrix,
        y: &Bits,
        groups: &[u32],
        roles: &[Role],
        folds: usize,
        bytes: &[u8],
    ) -> Vec<TrialResult> {
        let loaded = parse(bytes).unwrap();
        run_trials(
            "fixture__fixturemycin",
            "plan",
            x,
            y,
            groups,
            &[],
            roles,
            folds,
            &small(),
            &loaded,
            "0",
        )
        .unwrap()
    }

    /// The release's settings with a screen of 6 of the fixture's 12
    /// features, so a screen that reads a genome outside the fit's own
    /// fitting genomes can change what the forest, boosting and the lasso
    /// see.
    fn small() -> EvalConfig {
        EvalConfig {
            screened_features: 6,
            ..EvalConfig::default()
        }
    }

    /// Without a setting of its own, a trial is the fold file's fit: each
    /// learner's calls per fold at 0.5 equal those of `evaluate_assigned`
    /// on the same partition at the same seed, and the calls at another
    /// cutoff come from the same scores, so their log loss is the same.
    #[test]
    fn the_release_setting_is_the_fold_files_fit_and_a_cutoff_moves_only_the_calls() {
        let (x, y, groups) = fixture();
        let folds: Vec<usize> = (0..96).map(|i| (i / 4) % 3).collect();
        let roles = roles_from(&folds);
        let trial = &run(
            &x,
            &y,
            &groups,
            &roles,
            3,
            &settings_json(42, "copies", RELEASE_ALL),
        )[0];
        let plain = evaluate_assigned(
            "fixture__fixturemycin",
            "plan",
            &x,
            &y,
            &groups,
            &[],
            &roles,
            3,
            &small(),
        )
        .unwrap();
        let truth: Vec<bool> = plain.tested.iter().map(|&i| y.contains(i)).collect();
        for (row, want) in trial.learners.iter().zip(&plain.result.learners) {
            assert_eq!(row.status, "ok", "{}: {:?}", row.learner, row.reason);
            assert_eq!(row.fold_balanced_accuracy, want.fold_balanced_accuracy);
            let at = |c: f64| row.calls.iter().find(|k| k.cutoff == c).unwrap();
            let expect = fold_calls(&truth, &plain.result.test_fold, &want.scores, 3, &[0.5]);
            assert_eq!(at(0.5).folds, expect[0].folds, "{}", row.learner);
            for c in [0.3, 0.7] {
                for (a, b) in at(c).folds.iter().zip(&at(0.5).folds) {
                    assert_eq!(a.loss_sum, b.loss_sum);
                    assert_eq!(a.n, b.n);
                }
            }
        }
        let rf = trial.learners.iter().find(|r| r.learner == "rf").unwrap();
        assert_eq!(rf.cutoff, Some(0.35));
        let cut: Vec<f64> = rf.calls.iter().map(|k| k.cutoff).collect();
        assert_eq!(cut, vec![0.3, 0.35, 0.5, 0.7]);
        let scm = trial.learners.iter().find(|r| r.learner == "scm").unwrap();
        // scm scores 0 or 1: every cutoff gives it the same calls.
        for k in &scm.calls {
            assert_eq!(k.folds, scm.calls[2].folds);
        }
    }

    /// With lineage tuning and the seed `fold_seed(s, f)`, the release
    /// setting fitted on outer fold f of the drawn lineage split (tested)
    /// and every other genome (training) scores fold f's genomes as the
    /// drawn run did.
    #[test]
    fn the_release_setting_on_one_lineage_fold_reproduces_the_drawn_fold() {
        let (x, y, groups) = fixture();
        let cfg = small();
        let out = evaluate("fixture__fixturemycin", &x, &y, &groups, &[], 3, None, &cfg).unwrap();
        let lineage = out.splits.iter().find(|s| s.scheme == "lineage").unwrap();
        assert!(lineage.unavailable_reason.is_none());
        for f in [0, 3] {
            let roles: Vec<Role> = lineage
                .test_fold
                .iter()
                .map(|&g| if g == f { Role::Test(0) } else { Role::Train })
                .collect();
            let seed = crate::ablate::plan::fold_seed(42, f);
            let trial = &run(
                &x,
                &y,
                &groups,
                &roles,
                1,
                &settings_json(seed, "lineage", RELEASE_ALL),
            )[0];
            let tested: Vec<usize> = (0..96).filter(|&i| lineage.test_fold[i] == f).collect();
            let truth: Vec<bool> = tested.iter().map(|&i| y.contains(i)).collect();
            for (row, want) in trial.learners.iter().zip(&lineage.learners) {
                let s: Vec<f64> = tested.iter().map(|&i| want.scores[i]).collect();
                let expect = calls_at(&truth, &s, 0.5);
                let got = &row.calls.iter().find(|k| k.cutoff == 0.5).unwrap().folds[0];
                assert_eq!(got, &expect, "{} fold {f}", row.learner);
                assert_eq!(
                    row.fold_balanced_accuracy[0],
                    want.fold_balanced_accuracy[f]
                );
            }
        }
    }

    /// Genomes listed as unused reach no fit: adding them, with any
    /// features and labels, changes no trial's numbers.
    #[test]
    fn unused_genomes_change_no_trial() {
        let (x, y, groups) = fixture();
        let folds: Vec<usize> = (0..96).map(|i| (i / 4) % 3).collect();
        let base = run(
            &x,
            &y,
            &groups,
            &roles_from(&folds),
            3,
            &settings_json(9, "lineage", SAMPLED),
        );
        // 24 more genomes in 6 more lineages, the opposite label on feature 0.
        let n = 120;
        let mut x2 = BitMatrix::with_cols(n);
        for f in 0..x.rows {
            let on = (0..n).filter(|&i| {
                if i < 96 {
                    x.get(f, i)
                } else {
                    (i + f) % 2 == 0
                }
            });
            x2.push_row(Bits::from_indices(n, on).words());
        }
        let y2 = Bits::from_indices(
            n,
            (0..n).filter(|&i| if i < 96 { y.contains(i) } else { i % 2 == 0 }),
        );
        let groups2: Vec<u32> = (0..n).map(|i| (i / 4) as u32).collect();
        let mut roles2 = roles_from(&folds);
        roles2.extend(std::iter::repeat_n(Role::Unused, 24));
        let more = run(
            &x2,
            &y2,
            &groups2,
            &roles2,
            3,
            &settings_json(9, "lineage", SAMPLED),
        );
        for (a, b) in base[0].learners.iter().zip(&more[0].learners) {
            assert_eq!(a.status, "ok", "{}: {:?}", a.learner, a.reason);
            assert_eq!(
                a.fold_balanced_accuracy, b.fold_balanced_accuracy,
                "{}",
                a.learner
            );
            assert_eq!(a.fold_complexity, b.fold_complexity, "{}", a.learner);
            for (p, q) in a.calls.iter().zip(&b.calls) {
                assert_eq!(p.folds, q.folds, "{}", a.learner);
            }
        }
    }

    /// No test genome's label reaches its own fold's fit: with the labels of
    /// fold 0's test genomes turned round, fold 0's fits call the same number
    /// of its genomes resistant at every cutoff and are as complex, for the
    /// release's setting and a sampled one. A screen that read a test
    /// genome's label could change what the forest, boosting and the lasso
    /// see, and so what they call.
    #[test]
    fn no_test_genomes_label_reaches_its_folds_fit() {
        let (x, y, groups) = fixture();
        let folds: Vec<usize> = (0..96).map(|i| (i / 4) % 3).collect();
        let roles = roles_from(&folds);
        let turned = Bits::from_indices(96, (0..96).filter(|&i| y.contains(i) != (folds[i] == 0)));
        for trials in [RELEASE_ALL, SAMPLED] {
            let a = run(
                &x,
                &y,
                &groups,
                &roles,
                3,
                &settings_json(11, "copies", trials),
            );
            let b = run(
                &x,
                &turned,
                &groups,
                &roles,
                3,
                &settings_json(11, "copies", trials),
            );
            for (p, q) in a[0].learners.iter().zip(&b[0].learners) {
                assert_eq!(p.status, "ok", "{}: {:?}", p.learner, p.reason);
                assert_eq!(q.status, "ok", "{}: {:?}", q.learner, q.reason);
                assert_eq!(p.fold_complexity[0], q.fold_complexity[0], "{}", p.learner);
                for (k, m) in p.calls.iter().zip(&q.calls) {
                    let called = |c: &Calls| c.tp + c.fp;
                    assert_eq!(
                        called(&k.folds[0]),
                        called(&m.folds[0]),
                        "{} at {}",
                        p.learner,
                        k.cutoff
                    );
                }
            }
        }
    }

    /// A sampled setting is fitted as given: scm keeps at most its rules,
    /// cart at most 2^depth leaves, the forest its trees.
    #[test]
    fn a_sampled_setting_is_fitted_as_given() {
        let (x, y, groups) = fixture();
        let folds: Vec<usize> = (0..96).map(|i| (i / 4) % 3).collect();
        let trial = &run(
            &x,
            &y,
            &groups,
            &roles_from(&folds),
            3,
            &settings_json(5, "copies", SAMPLED),
        )[0];
        let get = |l: &str| trial.learners.iter().find(|r| r.learner == l).unwrap();
        assert!(get("scm").fold_complexity.iter().all(|&c| c <= 2.0));
        assert!(get("cart").fold_complexity.iter().all(|&c| c <= 2.0));
        assert!(get("rf").fold_complexity.iter().all(|&c| c == 50.0));
        assert!(get("gbt").fold_complexity.iter().all(|&c| c <= 60.0));
        for row in &trial.learners {
            assert_eq!(row.status, "ok", "{}: {:?}", row.learner, row.reason);
        }
    }

    /// A learner whose fits fail is written as failed with its reason and
    /// the others are still fitted.
    #[test]
    fn a_failed_partition_is_recorded_per_learner() {
        let (x, y, groups) = fixture();
        // Fold 1 tests only susceptible genomes: the partition is refused.
        let folds: Vec<usize> = (0..96)
            .map(|i| if y.contains(i) { i % 2 * 2 } else { i % 3 })
            .collect();
        let trial = &run(
            &x,
            &y,
            &groups,
            &roles_from(&folds),
            3,
            &settings_json(5, "copies", SAMPLED),
        )[0];
        assert_eq!(trial.learners.len(), 5);
        for row in &trial.learners {
            assert_eq!(row.status, "failed");
            assert!(
                row.reason
                    .as_deref()
                    .unwrap()
                    .contains("lack one or both classes")
            );
        }
    }

    /// The inner folds keep every lineage whole, cover exactly the outer
    /// fold's training genomes, and are the partition the release's scm
    /// and cart tune on there.
    #[test]
    fn inner_folds_keep_lineages_whole_inside_the_outer_training_genomes() {
        let (_, y, groups) = fixture();
        let n = 96;
        let genomes: Vec<String> = (0..n).map(|i| format!("g{i}")).collect();
        let outer: Vec<Option<usize>> = (0..n).map(|i| Some((i / 4) % 5)).collect();
        let mut used = vec![true; n];
        used[0] = false;
        let mut outer_u = outer.clone();
        outer_u[0] = None;
        let file = folds_file(
            "d",
            "lineage",
            42,
            5,
            3,
            &genomes,
            &y,
            &groups,
            outer_u.clone(),
            used.clone(),
        );
        for f in 0..5 {
            let train = Bits::from_indices(n, (0..n).filter(|&i| used[i] && outer_u[i] != Some(f)));
            let want = training_folds(
                &y,
                &train,
                Some(&groups),
                3,
                crate::ablate::plan::fold_seed(42, f),
            );
            for (i, (&got, &wanted)) in file.inner[f].iter().zip(&want).enumerate() {
                assert_eq!(got.is_some(), train.contains(i), "genome {i} fold {f}");
                assert_eq!(got.unwrap_or(usize::MAX), wanted);
            }
            let mut of_lineage: BTreeMap<u32, usize> = BTreeMap::new();
            for (&got, &lineage) in file.inner[f].iter().zip(&groups) {
                if let Some(j) = got {
                    assert_eq!(*of_lineage.entry(lineage).or_insert(j), j, "lineage split");
                }
            }
        }
    }

    /// The parser refuses what the driver must never send.
    #[test]
    fn a_settings_file_outside_its_spaces_is_refused() {
        let bad = [
            r#"[{"id":1,"learners":{"scm":{"kind":"conjunction","p":0.01,"max_rules":2}}}]"#,
            r#"[{"id":1,"learners":{"cart":{"max_depth":31,"min_leaf":1}}}]"#,
            r#"[{"id":1,"learners":{"rf":{"n_trees":50,"subsample":0.5,"max_depth":3}}}]"#,
            r#"[{"id":1,"learners":{"lasso":{"n_lambdas":12,"alpha":1}}}]"#,
            r#"[{"id":1,"learners":{"knn":"release"}}]"#,
            r#"[{"id":1,"learners":{"rf":"release"}},{"id":1,"learners":{"rf":"release"}}]"#,
            r#"[{"id":1,"learners":{"rf":"release"},"cutoffs":{"gbt":0.4}}]"#,
            r#"[{"id":1,"learners":{"rf":"release"},"cutoffs":{"rf":1.0}}]"#,
            r#"[]"#,
        ];
        for trials in bad {
            assert!(
                parse(&settings_json(1, "copies", trials)).is_err(),
                "{trials}"
            );
        }
        assert!(parse(br#"{"format":"x","seed":1,"tuning":"copies","trials":[{"id":0,"learners":{"rf":"release"}}]}"#).is_err());
        assert!(parse(&settings_json(1, "copies", SAMPLED)).is_ok());
    }
}
