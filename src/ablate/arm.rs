//! One arm of a full-pipeline input comparison (`grm arm`).
//!
//! An arm is one input configuration of one dataset, already curated,
//! counted and built into a data root: its cap and selection rule chose the
//! genomes, its k counted them, and the arm names its feature blocks and its
//! retained-feature budget (`EvalConfig::screened_features`). This module
//! scores it once per outer fold of each split and repeat and writes the
//! contract §4 prediction rows that J pairs and K reports
//! (`scripts/input-comparison/run.py` drives every arm).
//!
//! Per outer fold, in the order `docs/ablation-runner.md` gives: every block
//! is encoded on the outer training samples (`fit_encoder`, `transform`),
//! filtered there (`block_candidates`), and the blocks are joined
//! (`assemble`) on the same mask; then each learner fits on the assembled
//! matrix with the plan's inner folds and holdouts (`learn::planned`). The
//! learners' own chi-square screen, at the arm's budget on each fit's own
//! fitting samples, is the only step that reads labels.
//!
//! Fold plans follow the pairing rules of DECISIONS V-1:
//! - the reference arm (no `reference`) draws keyed plans from the k = 31
//!   lineages and duplicate groups of the species grouping (`groups`, one
//!   `grm structure` over the union of the comparison's arms' cohorts, at the
//!   arm's own k; without it, its own) and writes them for the others;
//! - `same_plan` and `k31_groups` arms take the reference's plans as they
//!   are, and must have the reference's samples and labels: a k arm is
//!   counted at its own k but split on the k = 31 groups;
//! - a `common_outer_test_set_balanced` arm (another cap or selection rule)
//!   takes a `CohortPlan` from the reference plan: in each outer fold it is
//!   scored on a subsample of the genomes both cohorts hold that the
//!   reference tests there, drawn to the reference fold's class balance, and
//!   trains on its own cohort minus that fold's common test genomes and, on
//!   the lineage split, minus its own genomes of a scored genome's lineage.
//!
//! Every row of every arm carries the reference plan's `fold_plan_sha256`
//! and the reference's lineage cluster, since its `fold` is the reference's;
//! a cohort arm's run records its own plan's hash beside it.
//!
//! Repeats (main's rulings 1 and 2 after VERIFY-H1STATS): the lineage split
//! scores `repeats` keyed plans and the random split one, since only the
//! lineage split decides. A row's `repeat` is its slot, 0 to repeats − 1.
//! The reference draws keyed repeat `slot` for each slot; when a lineage draw
//! is refused, it draws the next keyed repeat past the slots (repeats,
//! repeats + 1, …), at most `MAX_REPLACEMENTS` times per split, keeps each
//! refused draw as `plans/refused-lineage-<repeat>.json` and lists every
//! draw in `repeats.json`. A plan file holds its keyed repeat, so the other
//! arms score the same draws. The twin seeds follow the draw: a learner,
//! balance or tuning seed `s` acts as `plan::repeat_seed(s, repeat)`, so a
//! twin draws anew in every repeat as the default does.

use crate::ablate::blocks::{
    self, BlockSource, ContextSource, CountSource, KmerSource, PanelSource, PredictorBlock,
    PresenceSource, ProbeSource, Provenance,
};
use crate::ablate::encode::{self, EncodeConfig};
use crate::ablate::plan::{self, CohortPlan, FoldPlan, GroupIds, Scheme, SpeciesGroups};
use crate::analysis::{self, ArmInputs};
use crate::bits::Bits;
use crate::catalog::Catalog;
use crate::curate::DatasetSpec;
use crate::dna::KmerSpace;
use crate::kset::KsetReader;
use crate::layout::Layout;
use crate::learn::eval::EvalConfig;
use crate::learn::planned::{self, PlanView};
use anyhow::{Context, Result, bail, ensure};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::time::Instant;

/// How an arm's predictions pair with the reference arm's (DECISIONS V-1).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Pairing {
    /// Same samples and the same fold plan (blocks and budget arms, and the
    /// reference itself).
    SamePlan,
    /// Same samples, counted at another k, split on the k = 31 groups.
    K31Groups,
    /// Another cohort, scored on a class-balanced subsample of the
    /// reference's common outer test set (`plan::CohortPlan`).
    CommonOuterTestSetBalanced,
}

impl Pairing {
    pub fn as_str(self) -> &'static str {
        match self {
            Pairing::SamePlan => "same_plan",
            Pairing::K31Groups => "k31_groups",
            Pairing::CommonOuterTestSetBalanced => "common_outer_test_set_balanced",
        }
    }
}

/// What one arm runs.
#[derive(Clone, Debug)]
pub struct ArmConfig {
    pub arm: String,
    /// Block ids: K, C, H, P (N and X are always unavailable).
    pub blocks: Vec<String>,
    /// `screened_features` is the arm's retained-feature budget.
    pub eval: EvalConfig,
    pub repeats: u32,
    pub pairing: Pairing,
    /// The reference arm's output folder; None for the reference itself.
    pub reference: Option<PathBuf>,
    /// A cohort arm's balanced-subsample seed instead of the reference's
    /// repeat seed (`plan::cohort_plan`).
    pub balance_seed: Option<u64>,
    /// The species grouping's `structure.json` (`load_species_groups`): the
    /// lineages and duplicate groups every arm of a comparison splits on.
    /// A cohort arm needs it; the reference arm without it uses its own.
    pub groups: Option<PathBuf>,
    /// The learners' own seed instead of the plan's (`PlanView::learner_seed`):
    /// the calibration twin of the default.
    pub learner_seed: Option<u64>,
    /// A `same_plan` arm's tuning seed: it scores the reference's outer folds
    /// with inner folds and holdouts drawn again (`plan::retuned_plan`), the
    /// retuned twin that reaches scm, cart and the lasso.
    pub tuning_seed: Option<u64>,
    /// Write `training.jsonl`: a first line with the arm's sample ids, then per
    /// split, repeat and outer fold the genome ids `run_arm` handed the
    /// learners as training samples.
    pub trace_training: bool,
}

/// The splits an arm scores, in `evaluate`'s order.
pub const SCHEMES: [Scheme; 2] = [Scheme::Random, Scheme::Lineage];

/// The replacement draws a reference arm may make for refused lineage
/// repeats, per dataset (main's ruling 2 after VERIFY-H1STATS).
pub const MAX_REPLACEMENTS: u32 = 2;

/// The repeats a split scores: `repeats` on the lineage split, which decides,
/// and one on the random split, which is descriptive.
pub fn split_repeats(cfg: &ArmConfig, scheme: Scheme) -> u32 {
    match scheme {
        Scheme::Lineage => cfg.repeats,
        Scheme::Random => 1,
    }
}

/// A reference arm's plan for one split and repeat, with the samples it is
/// over, their labels and their lineage clusters.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PlanFile {
    pub sample_ids: Vec<String>,
    pub y: Vec<u8>,
    pub cluster: Vec<u32>,
    /// The k-mer duplicate groups as sample indices (`kset_duplicates`).
    pub duplicates: Vec<Vec<usize>>,
    pub plan: FoldPlan,
}

/// One contract §4 row, J's field names.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct PredictionRow {
    pub dataset: String,
    pub arm: String,
    pub learner: String,
    pub split: String,
    pub repeat: u32,
    pub sample_id: String,
    pub y: u8,
    pub fold: usize,
    pub score: f64,
    pub called: bool,
    pub cluster: u32,
    pub fold_plan_sha256: String,
}

/// One (learner, split) of an arm, in the shape of K's run record without
/// its `id`; `run.py` adds the id and the reference arm.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ArmRun {
    pub dataset: String,
    pub arm: String,
    pub blocks: Vec<String>,
    pub learner: String,
    pub split: String,
    pub repeats: u32,
    /// The reference plan of each repeat, the hash every row carries.
    pub fold_plan_sha256: Option<Vec<String>>,
    pub status: String,
    pub reason: Option<String>,
    pub metrics: serde_json::Map<String, serde_json::Value>,
    pub fit_count: usize,
    pub seconds: f64,
    pub predictions: Option<String>,
    pub pairing: Pairing,
    /// A cohort arm's own `CohortPlan` hash per repeat; None otherwise.
    pub arm_plan_sha256: Option<Vec<String>>,
    /// A cohort arm's per-fold drop counts per repeat (`CohortPlan`'s
    /// `balance_dropped` and `lineage_dropped`); None otherwise.
    pub cohort_drops: Option<Vec<CohortDrops>>,
}

/// One repeat's per-outer-fold drop counts of a cohort arm.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct CohortDrops {
    pub repeat: u32,
    /// Common test genomes the class-balanced subsample left out.
    pub balance_dropped: Vec<usize>,
    /// Arm-only genomes of a scored genome's lineage left out of training.
    pub lineage_dropped: Vec<usize>,
    /// Genomes of a scored genome's k-mer duplicate group left out of training.
    pub duplicate_dropped: Vec<usize>,
}

/// What `grm arm` records about the arm's inputs, in `arm.json`.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ArmRecord {
    pub arm: String,
    pub dataset: String,
    pub species: String,
    pub k: usize,
    pub samples: usize,
    pub min_support: u32,
    pub eval_config: serde_json::Value,
    pub selection: Option<crate::curate::Selection>,
    pub pairing: Pairing,
    pub blocks: Vec<blocks::BlockSpec>,
    /// Each requested block that is unavailable, with its reason.
    pub unavailable: Vec<(String, String)>,
}

/// The ksets of the species' matrix genomes, in matrix column order.
fn kset_paths(layout: &Layout, slug: &str, genomes: &[String]) -> Vec<PathBuf> {
    let dir = layout.kset_dir(slug);
    genomes
        .iter()
        .map(|g| dir.join(format!("{g}.kset")))
        .collect()
}

/// A requested block admitted over the arm's samples, or the reason its
/// source has none.
fn block(
    id: &str,
    layout: &Layout,
    spec: &DatasetSpec,
    input: &analysis::SpeciesInputs,
    catalog: Option<&Catalog>,
    samples: &[String],
) -> Result<std::result::Result<PredictorBlock, String>> {
    // `blocks::admit` takes only H's sealed adapters, so each arm below
    // names its concrete source.
    macro_rules! admit {
        ($source:expr, $samples:expr) => {{
            let source = $source;
            match &source.spec().unavailable_reason {
                Some(reason) => Ok(Err(reason.clone())),
                None => blocks::admit(&source, $samples).map(Ok),
            }
        }};
    }
    let genomes = &input.matrix.summary.genomes;
    let k = input.matrix.summary.k;
    let catalog_input = || -> Result<std::result::Result<(&Catalog, Provenance), String>> {
        let Some(cat) = catalog else {
            return Ok(Err("no catalog gene calls were counted".into()));
        };
        Ok(Ok((cat, Provenance::catalog_dir(&layout.catalog_dir())?)))
    };
    match id {
        "K" => admit!(KmerSource::new(&input.matrix)?, samples),
        "C" => {
            let ksets = kset_paths(layout, &spec.species_slug, genomes);
            let source = ProbeSource::from_ksets(
                &spec.species,
                genomes.clone(),
                &ksets,
                &KmerSpace::new(k),
            )?;
            admit!(source, samples)
        }
        "H" => {
            let (cat, provenance) = match catalog_input()? {
                Ok(c) => c,
                Err(why) => return Ok(Err(why)),
            };
            let calls = genomes
                .iter()
                .map(|g| input.hits.get(g).cloned().unwrap_or_default())
                .collect();
            let source = PanelSource::new(
                &spec.antibiotic,
                &provenance,
                &cat.calls_version(),
                genomes.clone(),
                calls,
            )?;
            admit!(source, samples)
        }
        "P" => {
            let (cat, provenance) = match catalog_input()? {
                Ok(c) => c,
                Err(why) => return Ok(Err(why)),
            };
            let ksets = kset_paths(layout, &spec.species_slug, genomes);
            let kmers = ksets.iter().map(|p| KsetReader::open(p)?.read_all());
            admit!(
                PresenceSource::new(cat, &provenance, genomes.clone(), kmers)?,
                samples
            )
        }
        "N" => admit!(CountSource::new(), samples),
        "X" => admit!(ContextSource::new(), samples),
        other => bail!("unknown block {other:?}; blocks are C, H, P, N, X and K"),
    }
}

/// The plan an arm scores one split and repeat on, with the hash its rows
/// carry, its own plan's hash when it differs, and each sample's cluster.
struct ArmPlan {
    fold: Option<FoldPlan>,
    cohort: Option<CohortPlan>,
    rows_sha256: String,
    arm_sha256: Option<String>,
    cluster: Vec<Option<u32>>,
    /// The keyed repeat the reference drew (its `FoldPlan::repeat`).
    repeat: u32,
}

impl ArmPlan {
    fn view(&self) -> PlanView<'_> {
        match (&self.fold, &self.cohort) {
            (Some(p), _) => PlanView::from(p),
            (None, Some(c)) => PlanView::from(c),
            (None, None) => unreachable!("an arm plan has a plan"),
        }
    }
}

/// The partitions `run_arm` scores one plan on: the plan's own (a cohort
/// plan's held-out genomes included), with the arm's learner seed keyed by
/// the plan's repeat.
fn fold_view<'a>(p: &'a ArmPlan, cfg: &ArmConfig) -> PlanView<'a> {
    let mut view = p.view();
    view.learner_seed = cfg.learner_seed.map(|s| plan::repeat_seed(s, p.repeat));
    view
}

#[cfg(test)]
thread_local! {
    /// In test builds, the `train` mask of every `fold_matrix` call on this
    /// thread, in call order, so a test can check what `run_arm` handed each
    /// fold (REVIEW-INPUT18 M1).
    static FOLD_MATRIX_TRAIN: std::cell::RefCell<Vec<Bits>> =
        const { std::cell::RefCell::new(Vec::new()) };
}

/// One outer fold's arm matrix over every sample: each block's encoder is
/// fitted on the fold's `train` samples, the block is transformed, filtered on
/// `train` (`block_candidates`), and the blocks are joined on `train`.
fn fold_matrix(
    raws: &[&blocks::RawBlock],
    samples: &[String],
    train: &Bits,
    ecfg: &EncodeConfig,
) -> Result<crate::bits::BitMatrix> {
    #[cfg(test)]
    FOLD_MATRIX_TRAIN.with(|seen| seen.borrow_mut().push(train.clone()));
    let mut candidates = Vec::with_capacity(raws.len());
    for raw in raws {
        let enc = encode::fit_encoder(raw, samples, train, ecfg)?;
        let block = encode::transform(&enc, raw, samples)?;
        candidates.push(encode::block_candidates(&block, train, ecfg));
    }
    Ok(encode::assemble(&candidates, train).x)
}

fn plan_name(scheme: Scheme, slot: u32) -> String {
    format!("{scheme}-{slot}")
}

/// One split's draws, as `repeats.json` lists them.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct SplitDraws {
    pub split: String,
    /// The keyed repeat each slot scores.
    pub drawn: Vec<u32>,
    /// The reference's refused draws, each with its reason.
    pub refused: Vec<RefusedDraw>,
    /// A retuned twin's refused tuning draws, each with its slot and reason.
    pub refused_tuning: Vec<RefusedTuning>,
    pub max_replacements: u32,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct RefusedTuning {
    pub slot: u32,
    pub draw: u32,
    pub reason: String,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct RefusedDraw {
    pub repeat: u32,
    pub reason: String,
}

/// The plans an arm scores one split on, one per slot, and, when a slot is
/// left without a valid plan after the replacement draws, why. The reference
/// draws keyed repeat `slot` in each slot and, on the lineage split, replaces
/// a refused draw with the next keyed repeat past the slots, at most
/// `MAX_REPLACEMENTS` times; the other arms take the reference's slot plans.
/// A retuned twin draws its tuning partitions again, keyed by draw, when they
/// are refused, at most `MAX_REPLACEMENTS` times per split (`arm_plan`).
fn arm_plans(
    cfg: &ArmConfig,
    a: &ArmInputs,
    species: Option<&SpeciesGroups>,
    scheme: Scheme,
    plans_dir: &Path,
) -> Result<(Vec<ArmPlan>, SplitDraws, Option<String>)> {
    let n = split_repeats(cfg, scheme);
    let replace = cfg.reference.is_none() && scheme == Scheme::Lineage;
    let (mut plans, mut refused, mut spare, mut short) = (Vec::new(), Vec::new(), n, None);
    let mut refused_tuning = Vec::new();
    for slot in 0..n {
        let mut repeat = slot;
        loop {
            let p = arm_plan(
                cfg,
                a,
                species,
                scheme,
                slot,
                repeat,
                plans_dir,
                &mut refused_tuning,
            )?;
            let why = p.fold.as_ref().and_then(|f| f.unavailable_reason.clone());
            match why {
                Some(reason) if replace && (refused.len() as u32) < MAX_REPLACEMENTS => {
                    let name = plan_name(scheme, slot);
                    std::fs::rename(
                        plans_dir.join(format!("{name}.json")),
                        plans_dir.join(format!("refused-{scheme}-{repeat}.json")),
                    )?;
                    refused.push(RefusedDraw { repeat, reason });
                    repeat = spare;
                    spare += 1;
                }
                Some(reason) if replace => {
                    short.get_or_insert(format!(
                        "{scheme} repeat {slot} (keyed repeat {repeat}) is refused: {reason}; with {MAX_REPLACEMENTS} \
                         replacement draws spent, fewer than {n} valid {scheme} repeats remain"
                    ));
                    plans.push(p);
                    break;
                }
                _ => {
                    plans.push(p);
                    break;
                }
            }
        }
    }
    let draws = SplitDraws {
        split: scheme.to_string(),
        drawn: plans.iter().map(|p| p.repeat).collect(),
        refused,
        refused_tuning,
        max_replacements: MAX_REPLACEMENTS,
    };
    Ok((plans, draws, short))
}

fn labels(y: &Bits) -> Vec<u8> {
    (0..y.len()).map(|i| y.contains(i) as u8).collect()
}

fn write_json(path: &Path, value: &impl Serialize) -> Result<()> {
    std::fs::write(path, serde_json::to_vec_pretty(value)?)
        .with_context(|| format!("write {}", path.display()))
}

/// The species grouping a comparison's arms share, from the `structure.json`
/// one `grm structure` wrote over the comparison's species pool (the union of
/// its arms' cohorts), with the k it was counted at.
pub fn load_species_groups(path: &Path) -> Result<(SpeciesGroups, usize)> {
    let s: crate::species::Structure = serde_json::from_slice(
        &std::fs::read(path).with_context(|| format!("read {}", path.display()))?,
    )
    .with_context(|| format!("parse {}", path.display()))?;
    ensure!(
        s.kset_digests.len() == s.genomes.len(),
        "{} records no kset digest per genome; rerun `grm structure`",
        path.display()
    );
    let k = s.k;
    let groups = SpeciesGroups::new(s.genomes, s.lineages.labels, s.kset_digests)
        .with_context(|| format!("species grouping {}", path.display()))?;
    Ok((groups, k))
}

#[allow(clippy::too_many_arguments)]
fn arm_plan(
    cfg: &ArmConfig,
    a: &ArmInputs,
    species: Option<&SpeciesGroups>,
    scheme: Scheme,
    slot: u32,
    repeat: u32,
    plans_dir: &Path,
    refused_tuning: &mut Vec<RefusedTuning>,
) -> Result<ArmPlan> {
    let name = plan_name(scheme, slot);
    let Some(reference) = &cfg.reference else {
        ensure!(
            cfg.pairing == Pairing::SamePlan,
            "the reference arm pairs with itself by same_plan, not {}",
            cfg.pairing.as_str()
        );
        ensure!(
            cfg.tuning_seed.is_none(),
            "a tuning seed retunes the reference's plans, so the reference itself takes none"
        );
        let groups = match species {
            Some(s) => s.groups(&a.sample_ids)?,
            None => GroupIds::new(a.lineage.clone(), a.duplicates.clone(), None)?,
        };
        // Keyed on genome identity (main's ruling on REVIEW-INPUT2 N1); `grm
        // learn` keeps G's index-drawn folds.
        let p = plan::keyed_fold_plan(
            &a.sample_ids,
            &a.y,
            &groups,
            scheme,
            cfg.eval.folds,
            cfg.eval.inner_folds,
            cfg.eval.seed,
            repeat,
        )?;
        write_json(
            &plans_dir.join(format!("{name}.json")),
            &PlanFile {
                sample_ids: a.sample_ids.clone(),
                y: labels(&a.y),
                cluster: groups.lineage().to_vec(),
                duplicates: groups.duplicate().to_vec(),
                plan: p.clone(),
            },
        )?;
        return Ok(ArmPlan {
            rows_sha256: p.sha256.clone(),
            arm_sha256: None,
            cluster: groups.lineage().iter().map(|&c| Some(c)).collect(),
            repeat: p.repeat,
            fold: Some(p),
            cohort: None,
        });
    };
    let path = reference.join("plans").join(format!("{name}.json"));
    let r: PlanFile = serde_json::from_slice(
        &std::fs::read(&path)
            .with_context(|| format!("read the reference plan {}", path.display()))?,
    )
    .with_context(|| format!("parse {}", path.display()))?;
    r.plan.verify()?;
    let drawn = r.plan.repeat;
    ensure!(
        cfg.tuning_seed.is_none() || cfg.pairing == Pairing::SamePlan,
        "a tuning seed retunes the reference's plan under same_plan, not {}",
        cfg.pairing.as_str()
    );
    match cfg.pairing {
        Pairing::SamePlan | Pairing::K31Groups => {
            ensure!(
                r.sample_ids == a.sample_ids && r.y == labels(&a.y),
                "pairing {} needs the reference's samples and labels in its order; this arm has {} samples, the reference {}",
                cfg.pairing.as_str(),
                a.sample_ids.len(),
                r.sample_ids.len()
            );
            let cluster = r.cluster.iter().map(|&c| Some(c)).collect();
            if let (Some(seed), true) = (cfg.tuning_seed, r.plan.is_available()) {
                // The groups the reference keyed its plan on: a species grouping's
                // `groups` carries no source, so its lineages and duplicates rebuild it.
                let groups = GroupIds::new(r.cluster.clone(), r.duplicates.clone(), None)?;
                // A refused draw is kept and replaced by the next keyed draw, at
                // most MAX_REPLACEMENTS times per split; past that its reason stands.
                let mut draw = 0;
                let t = loop {
                    let t = plan::retuned_plan(&r.plan, &r.sample_ids, &a.y, &groups, seed, draw)?;
                    match &t.unavailable_reason {
                        Some(reason) if (refused_tuning.len() as u32) < MAX_REPLACEMENTS => {
                            write_json(
                                &plans_dir
                                    .join(format!("refused-{name}-tuning-{draw}.retuned.json")),
                                &serde_json::json!({ "tuning_seed": seed, "draw": draw, "plan": t }),
                            )?;
                            refused_tuning.push(RefusedTuning {
                                slot,
                                draw,
                                reason: reason.clone(),
                            });
                            draw += 1;
                        }
                        _ => break t,
                    }
                };
                write_json(
                    &plans_dir.join(format!("{name}.retuned.json")),
                    &serde_json::json!({ "tuning_seed": seed, "draw": draw, "plan": t }),
                )?;
                // The run records no plan hash of its own: K keeps arm_plan_sha256 for
                // cohort pairings, and the retuned plan is in its .retuned.json.
                return Ok(ArmPlan {
                    rows_sha256: r.plan.sha256.clone(),
                    arm_sha256: None,
                    cluster,
                    repeat: drawn,
                    fold: Some(t),
                    cohort: None,
                });
            }
            Ok(ArmPlan {
                rows_sha256: r.plan.sha256.clone(),
                arm_sha256: None,
                cluster,
                repeat: drawn,
                fold: Some(r.plan),
                cohort: None,
            })
        }
        Pairing::CommonOuterTestSetBalanced => {
            let at: HashMap<&str, usize> = r
                .sample_ids
                .iter()
                .enumerate()
                .map(|(i, s)| (s.as_str(), i))
                .collect();
            for (i, id) in a.sample_ids.iter().enumerate() {
                if let Some(&j) = at.get(id.as_str()) {
                    ensure!(
                        r.y[j] == a.y.contains(i) as u8,
                        "genome {id} has another label in the reference cohort"
                    );
                }
            }
            if !r.plan.is_available() {
                // The reference cannot score this split, so neither can a
                // cohort scored on its folds; the reason carries over.
                return Ok(ArmPlan {
                    rows_sha256: r.plan.sha256.clone(),
                    arm_sha256: None,
                    cluster: vec![None; a.sample_ids.len()],
                    repeat: drawn,
                    fold: Some(r.plan),
                    cohort: None,
                });
            }
            let species = species.context(
                "a common_outer_test_set_balanced arm needs --groups, the species grouping every arm shares",
            )?;
            let split_on = species.groups(&r.sample_ids)?;
            ensure!(
                split_on.lineage() == r.cluster.as_slice()
                    && split_on.duplicate() == r.duplicates.as_slice(),
                "the reference plan {} was not split on this species grouping",
                path.display()
            );
            let reference_y =
                Bits::from_indices(r.y.len(), (0..r.y.len()).filter(|&j| r.y[j] == 1));
            let c = plan::cohort_plan(
                &r.plan,
                &r.sample_ids,
                &reference_y,
                &a.sample_ids,
                &a.y,
                species,
                cfg.balance_seed.map(|s| plan::repeat_seed(s, drawn)),
            )?;
            write_json(&plans_dir.join(format!("{name}.cohort.json")), &c)?;
            let cluster = a
                .sample_ids
                .iter()
                .map(|id| at.get(id.as_str()).map(|&j| r.cluster[j]))
                .collect();
            Ok(ArmPlan {
                rows_sha256: r.plan.sha256.clone(),
                arm_sha256: Some(c.sha256.clone()),
                cluster,
                repeat: drawn,
                fold: None,
                cohort: Some(c),
            })
        }
    }
}

/// Numeric metrics of a learner result, for K's run record.
fn metrics(
    result: &crate::learn::eval::LearnerResult,
) -> Result<serde_json::Map<String, serde_json::Value>> {
    let mut out = serde_json::Map::new();
    if let serde_json::Value::Object(m) = serde_json::to_value(&result.metrics)? {
        for (k, v) in m {
            if v.is_number() || v.is_null() {
                out.insert(k, v);
            }
        }
    }
    Ok(out)
}

/// Scores one arm of `dataset` in `layout` and writes, under `out`:
/// `predictions.jsonl` (contract §4 rows), `runs.json` (one `ArmRun` per
/// learner and split), `arm.json` (`ArmRecord`) and `plans/` (the reference
/// arm's `PlanFile`s, or a cohort arm's `CohortPlan`s).
pub fn run_arm(
    layout: &Layout,
    slug: &str,
    dataset: &str,
    cfg: &ArmConfig,
    out: &Path,
) -> Result<()> {
    ensure!(cfg.repeats >= 1, "an arm needs at least one repeat");
    ensure!(!cfg.blocks.is_empty(), "an arm needs at least one block");
    let (input, catalog) = analysis::load_inputs(layout, slug)?;
    let spec = input
        .datasets
        .iter()
        .find(|d| d.id == dataset)
        .with_context(|| format!("{slug} has no curated dataset {dataset}"))?;
    let a = analysis::arm_inputs(&input, spec)?
        .with_context(|| format!("{dataset} is short of its class minimum"))?;
    let species = match cfg.groups.as_deref() {
        Some(path) => {
            let (groups, k) = load_species_groups(path)?;
            // Lineages and duplicate groups read k (docs/input-configuration.md),
            // so an arm splits only on a grouping counted at its own k.
            ensure!(
                k == input.matrix.summary.k,
                "the species grouping {} was counted at k = {k}, this arm at k = {}",
                path.display(),
                input.matrix.summary.k
            );
            Some(groups)
        }
        None => None,
    };
    let plans_dir = out.join("plans");
    std::fs::create_dir_all(&plans_dir)?;

    let mut admitted = Vec::new();
    let mut unavailable = Vec::new();
    for id in &cfg.blocks {
        match block(id, layout, spec, &input, catalog.as_ref(), &a.sample_ids)? {
            Ok(b) => admitted.push(b),
            Err(reason) => unavailable.push((id.clone(), reason)),
        }
    }
    write_json(
        &out.join("arm.json"),
        &ArmRecord {
            arm: cfg.arm.clone(),
            dataset: dataset.to_string(),
            species: spec.species.clone(),
            k: input.matrix.summary.k,
            samples: a.sample_ids.len(),
            min_support: a.min_support,
            eval_config: serde_json::to_value(&cfg.eval)?,
            selection: spec.selection.clone(),
            pairing: cfg.pairing,
            blocks: admitted.iter().map(|b| b.spec().clone()).collect(),
            unavailable: unavailable.clone(),
        },
    )?;
    let raws: Vec<&blocks::RawBlock> = admitted.iter().map(|b| b.raw()).collect();
    score_arm(dataset, &a, species.as_ref(), raws, &unavailable, cfg, out)
}

/// `run_arm`'s scoring, after it has loaded the data root and admitted the
/// blocks: every split and repeat's plan (`arm_plan`), then each outer fold's
/// matrix of `raws` on that fold's training samples (`fold_matrix`) and the
/// learners. Writes `predictions.jsonl`, `runs.json` and the plans under
/// `out`. A split is unavailable when any block in `unavailable` is.
fn score_arm(
    dataset: &str,
    a: &ArmInputs,
    species: Option<&SpeciesGroups>,
    raws: Vec<&blocks::RawBlock>,
    unavailable: &[(String, String)],
    cfg: &ArmConfig,
    out: &Path,
) -> Result<()> {
    let plans_dir = out.join("plans");
    std::fs::create_dir_all(&plans_dir)?;
    let ecfg = EncodeConfig {
        min_support: a.min_support,
    };
    let rows_path = out.join("predictions.jsonl");
    let mut rows = std::io::BufWriter::new(std::fs::File::create(&rows_path)?);
    let mut trace = if cfg.trace_training {
        let mut t = std::io::BufWriter::new(std::fs::File::create(out.join("training.jsonl"))?);
        serde_json::to_writer(&mut t, &serde_json::json!({ "sample_ids": a.sample_ids }))?;
        t.write_all(b"\n")?;
        Some(t)
    } else {
        None
    };
    let mut runs: Vec<ArmRun> = Vec::new();
    let mut draws = Vec::new();
    let learners = crate::learn::eval::LEARNERS;
    for scheme in SCHEMES {
        let mut results = Vec::new();
        let mut reason = unavailable
            .first()
            .map(|(id, why)| format!("block {id} is unavailable: {why}"));
        let started = Instant::now();
        let repeats = split_repeats(cfg, scheme);
        let (plans, drawn, short) = arm_plans(cfg, a, species, scheme, &plans_dir)?;
        draws.push(drawn);
        if let Some(why) = short {
            reason.get_or_insert(why);
        }
        for (slot, p) in plans.into_iter().enumerate() {
            let repeat = slot as u32;
            if reason.is_some() {
                results.push((p, None));
                continue;
            }
            let view = fold_view(&p, cfg);
            let split =
                planned::evaluate_plan(dataset, &a.y, &view, &cfg.eval, &mut |f, train| {
                    if let Some(t) = trace.as_mut() {
                        let ids: Vec<&str> = train
                            .iter_ones()
                            .map(|i| a.sample_ids[i].as_str())
                            .collect();
                        let line = serde_json::json!({"split": scheme.to_string(), "repeat": repeat,
                                                  "fold": f, "train": ids});
                        serde_json::to_writer(&mut *t, &line)?;
                        t.write_all(b"\n")?;
                    }
                    fold_matrix(&raws, &a.sample_ids, train, &ecfg)
                })?;
            if let Some(why) = &split.result.unavailable_reason {
                reason.get_or_insert_with(|| why.clone());
            }
            results.push((p, Some(split)));
        }
        let seconds = started.elapsed().as_secs_f64();
        let rows_sha: Vec<String> = results.iter().map(|(p, _)| p.rows_sha256.clone()).collect();
        let arm_sha: Option<Vec<String>> =
            results.iter().map(|(p, _)| p.arm_sha256.clone()).collect();
        let drops: Option<Vec<CohortDrops>> = results
            .iter()
            .enumerate()
            .map(|(r, (p, _))| {
                p.cohort.as_ref().map(|c| CohortDrops {
                    repeat: r as u32,
                    balance_dropped: c.balance_dropped.clone(),
                    lineage_dropped: c.lineage_dropped.clone(),
                    duplicate_dropped: c.duplicate_dropped.clone(),
                })
            })
            .collect();
        for (li, learner) in learners.iter().enumerate() {
            let mut run = ArmRun {
                dataset: dataset.to_string(),
                arm: cfg.arm.clone(),
                blocks: cfg.blocks.clone(),
                learner: learner.to_string(),
                split: scheme.to_string(),
                repeats,
                fold_plan_sha256: Some(rows_sha.clone()),
                status: "ok".into(),
                reason: None,
                metrics: serde_json::Map::new(),
                fit_count: cfg.eval.folds * repeats as usize,
                seconds: seconds / learners.len() as f64,
                predictions: Some("predictions.jsonl".into()),
                pairing: cfg.pairing,
                arm_plan_sha256: arm_sha.clone(),
                cohort_drops: drops.clone(),
            };
            if let Some(why) = &reason {
                run.status = "unavailable".into();
                run.reason = Some(why.clone());
                run.predictions = None;
                run.fit_count = 0;
                runs.push(run);
                continue;
            }
            for (repeat, (p, split)) in results.iter().enumerate() {
                let split = split.as_ref().expect("an available split was scored");
                let result = &split.result.learners[li];
                ensure!(result.learner == *learner, "learner order changed");
                if repeat == 0 {
                    run.metrics = metrics(result)?;
                }
                for (j, &i) in split.tested.iter().enumerate() {
                    let fold = split.result.test_fold[j];
                    let s = result.scores[j];
                    let (score, called) = if s.is_nan() {
                        (split.train_prevalence[fold], false)
                    } else {
                        (s, s >= 0.5)
                    };
                    let row = PredictionRow {
                        dataset: dataset.to_string(),
                        arm: cfg.arm.clone(),
                        learner: learner.to_string(),
                        split: scheme.to_string(),
                        repeat: repeat as u32,
                        sample_id: a.sample_ids[i].clone(),
                        y: a.y.contains(i) as u8,
                        fold,
                        score,
                        called,
                        cluster: p.cluster[i]
                            .context("a scored sample has no reference cluster")?,
                        fold_plan_sha256: p.rows_sha256.clone(),
                    };
                    serde_json::to_writer(&mut rows, &row)?;
                    rows.write_all(b"\n")?;
                }
            }
            runs.push(run);
        }
    }
    rows.flush()?;
    if let Some(t) = trace.as_mut() {
        t.flush()?;
    }
    write_json(&out.join("repeats.json"), &draws)?;
    write_json(&out.join("runs.json"), &runs)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ablate::blocks::{RawBlock, RawColumn};
    use crate::ablate::plan::fold_plan;

    /// A presence block over 20 samples: training samples 0..15, test 15..20.
    fn raw(samples: &[String]) -> RawBlock {
        let n = samples.len();
        let bits = |f: &dyn Fn(usize) -> bool| Bits::from_indices(n, (0..n).filter(|&i| f(i)));
        let present = |f: &dyn Fn(usize) -> bool| RawColumn::Binary {
            present: bits(f),
            assessable: Bits::ones(n),
        };
        // Training counts 0 x6, 2 x3, 4 x6 give cutpoints 2 and 4; adding the
        // test samples' 1s gives 1 and 4 instead, which moves the test bits.
        let count: Vec<Option<u32>> = (0..n)
            .map(|i| {
                Some(if i >= 15 {
                    1
                } else if i < 6 {
                    0
                } else if i < 9 {
                    2
                } else {
                    4
                })
            })
            .collect();
        RawBlock {
            spec: blocks::def("P").unwrap().spec(Provenance {
                version: "fixture".into(),
                source_sha256: "0".repeat(64),
                catalog_date: None,
            }),
            sample_ids: samples.to_vec(),
            names: vec!["P:kept".into(), "P:test_only".into(), "P:count".into()],
            cols: vec![
                present(&|i| i < 8),
                present(&|i| i >= 15),
                RawColumn::count(count),
            ],
        }
    }

    /// `run_arm`'s fold matrix reads only the fold's training samples: the
    /// encoders are fitted on them and the support filter runs on them. Either
    /// step on every sample gives another matrix (the reviewer's M5 and M6).
    #[test]
    fn the_fold_matrix_reads_only_the_training_samples() {
        let samples: Vec<String> = (0..20).map(|i| format!("g{i}")).collect();
        let block = raw(&samples);
        let train = Bits::from_indices(20, 0..15);
        let all = Bits::ones(20);
        let ecfg = EncodeConfig { min_support: 3 };
        let got = fold_matrix(&[&block], &samples, &train, &ecfg).unwrap();
        let build = |fit: &Bits, filter: &Bits| {
            let enc = encode::fit_encoder(&block, &samples, fit, &ecfg).unwrap();
            let x = encode::transform(&enc, &block, &samples).unwrap();
            encode::assemble(&[encode::block_candidates(&x, filter, &ecfg)], &train).x
        };
        let rows =
            |x: &crate::bits::BitMatrix| (0..x.rows).map(|r| x.row(r).to_vec()).collect::<Vec<_>>();
        assert_eq!(rows(&got), rows(&build(&train, &train)));
        assert_ne!(
            rows(&got),
            rows(&build(&train, &all)),
            "filtered on every sample"
        );
        assert_ne!(
            rows(&got),
            rows(&build(&all, &train)),
            "encoder fitted on every sample"
        );
    }

    /// A presence block `id` with rows for `rows` (indices into `samples`):
    /// one always-assessable column per name, present on the samples given.
    fn presence(
        id: &str,
        samples: &[String],
        rows: &[usize],
        cols: &[(&str, Vec<usize>)],
    ) -> RawBlock {
        let n = rows.len();
        RawBlock {
            spec: blocks::def(id).unwrap().spec(Provenance {
                version: "fixture".into(),
                source_sha256: "0".repeat(64),
                catalog_date: None,
            }),
            sample_ids: rows.iter().map(|&i| samples[i].clone()).collect(),
            names: cols
                .iter()
                .map(|(name, _)| format!("{id}:{name}"))
                .collect(),
            cols: cols
                .iter()
                .map(|(_, on)| RawColumn::Binary {
                    present: Bits::from_indices(n, (0..n).filter(|&j| on.contains(&rows[j]))),
                    assessable: Bits::ones(n),
                })
                .collect(),
        }
    }

    /// REVIEW-TL1 I4: a K-only arm can't show a leak in `assemble`'s
    /// cross-block coalescing, because `block_candidates` has already merged
    /// K's own training twins. Here a C column and a K column agree on every
    /// training sample (0..18) and differ on the outer test samples, and C
    /// has no row for test sample 23. `run_arm`'s fold matrix must collapse
    /// them on their training rows alone into one column that takes C's
    /// value where C has a row and K's at 23, and no change to the test
    /// samples (their bits, or which of them C has rows for) may change its
    /// training rows. A collapse keyed on every sample keeps the twins apart
    /// here and merges them once the test bits agree, so it fails both.
    #[test]
    fn a_c_and_k_fold_matrix_collapses_cross_block_twins_on_training_samples_only() {
        let samples: Vec<String> = (0..24).map(|i| format!("g{i}")).collect();
        let train = Bits::from_indices(24, 0..18);
        let ecfg = EncodeConfig { min_support: 3 };
        let all: Vec<usize> = (0..24).collect();
        let test: Vec<usize> = (18..24).collect();
        let on = |training: std::ops::Range<usize>, tested: &[usize]| -> Vec<usize> {
            training.chain(tested.iter().copied()).collect()
        };
        // Test bits of C:twin, C:solo, K:twin and K:other; C's rows.
        let toy = |t: [&[usize]; 4], c_rows: &[usize]| {
            let c = presence(
                "C",
                &samples,
                c_rows,
                &[("twin", on(0..6, t[0])), ("solo", on(6..12, t[1]))],
            );
            let k = presence(
                "K",
                &samples,
                &all,
                &[("twin", on(0..6, t[2])), ("other", on(12..18, t[3]))],
            );
            (c, k)
        };
        let fold =
            |(c, k): &(RawBlock, RawBlock)| fold_matrix(&[c, k], &samples, &train, &ecfg).unwrap();
        let rows =
            |x: &crate::bits::BitMatrix| (0..x.rows).map(|r| x.row(r).to_vec()).collect::<Vec<_>>();
        let training = |x: &crate::bits::BitMatrix| {
            (0..x.rows)
                .map(|r| train.and_slice(x.row(r)).words().to_vec())
                .collect::<Vec<_>>()
        };
        let no_23: Vec<usize> = (0..23).collect();
        let clean = toy([&[18, 19], &[20], &[20, 21, 23], &[22]], &no_23);

        // The toy can show a leak: each block keeps its twin as a candidate,
        // and the twins differ on test samples, so only a collapse on the
        // training rows makes 3 columns of the 4 candidates.
        let candidates: Vec<usize> = [&clean.0, &clean.1]
            .iter()
            .map(|raw| {
                let enc = encode::fit_encoder(raw, &samples, &train, &ecfg).unwrap();
                let x = encode::transform(&enc, raw, &samples).unwrap();
                encode::block_candidates(&x, &train, &ecfg).x.rows
            })
            .collect();
        assert_eq!(candidates, [2, 2]);
        let got = fold(&clean);
        let bits = |on: Vec<usize>| Bits::from_indices(24, on).words().to_vec();
        assert_eq!(
            rows(&got),
            [
                bits(on(0..6, &[18, 19, 23])),
                bits(on(6..12, &[20])),
                bits(on(12..18, &[22])),
            ],
            "C:twin and K:twin collapse to C's value where C has a row and K's at 23"
        );

        for (what, changed) in [
            ("no test bits", toy([&[], &[], &[], &[]], &no_23)),
            ("every test bit", toy([&test, &test, &test, &test], &no_23)),
            (
                "the twins' test bits swapped",
                toy([&[20, 21, 23], &[20], &[18, 19], &[22]], &no_23),
            ),
            (
                "C with a row for every sample",
                toy([&[18, 19], &[20], &[20, 21, 23], &[22]], &all),
            ),
        ] {
            assert_eq!(
                training(&fold(&changed)),
                training(&got),
                "{what} changed the training rows"
            );
        }
    }

    /// REVIEW-INPUT18 M1: the test above calls `fold_matrix` itself, so a leak
    /// where `run_arm` hands it anything but the fold's training samples (every
    /// sample, say) passed every test. This runs `run_arm`'s scoring on a C+K
    /// arm and checks the `train` of every fold matrix it builds: one per
    /// outer fold of each split, in order, each the plan's training samples
    /// for that fold. The toy can show such a leak: K's twin columns are C's
    /// with one sample flipped, so in some fold of each split a matrix on
    /// every sample keeps twins that the fold's own matrix collapses.
    #[test]
    fn run_arm_hands_each_c_and_k_fold_matrix_only_its_training_samples() {
        use rand::{Rng, SeedableRng};
        let n = 60;
        let samples: Vec<String> = (0..n).map(|i| format!("g{i:02}")).collect();
        let every: Vec<usize> = (0..n).collect();
        let mut rng = rand_chacha::ChaCha8Rng::seed_from_u64(18);
        let mut draw = || -> Vec<usize> { (0..n).filter(|_| rng.random_bool(0.4)).collect() };
        let c: Vec<Vec<usize>> = (0..8).map(|_| draw()).collect();
        // Each twin differs from its C column on one sample, 7 apart.
        let mut k: Vec<Vec<usize>> = c
            .iter()
            .enumerate()
            .map(|(i, on)| {
                let s = (7 * i + 3) % n;
                let mut on: Vec<usize> = on.iter().copied().filter(|&x| x != s).collect();
                if !c[i].contains(&s) {
                    on.push(s);
                }
                on
            })
            .collect();
        k.extend((0..4).map(|_| draw()));
        let names: Vec<String> = (0..k.len()).map(|i| format!("f{i}")).collect();
        let block = |id: &str, cols: &[Vec<usize>]| {
            let cols: Vec<(&str, Vec<usize>)> = names
                .iter()
                .zip(cols)
                .map(|(name, on)| (name.as_str(), on.clone()))
                .collect();
            presence(id, &samples, &every, &cols)
        };
        let (cb, kb) = (block("C", &c), block("K", &k));
        let raws = vec![&cb, &kb];
        let a = ArmInputs {
            sample_ids: samples.clone(),
            columns: every.clone(),
            y: Bits::from_indices(n, (0..n).filter(|i| i % 2 == 0)),
            min_support: 3,
            x: crate::bits::BitMatrix::with_cols(n),
            lineage: (0..n as u32).map(|i| i / 3).collect(),
            duplicates: vec![],
        };
        let cfg = ArmConfig {
            arm: "k-c".into(),
            blocks: vec!["K".into(), "C".into()],
            eval: EvalConfig::default(),
            repeats: 1,
            pairing: Pairing::SamePlan,
            reference: None,
            balance_seed: None,
            learner_seed: None,
            tuning_seed: None,
            trace_training: false,
            groups: None,
        };
        let dir = tempfile::tempdir().unwrap();
        FOLD_MATRIX_TRAIN.with(|seen| seen.borrow_mut().clear());
        score_arm("ds", &a, None, raws.clone(), &[], &cfg, dir.path()).unwrap();
        let got = FOLD_MATRIX_TRAIN.with(|seen| seen.take());

        let ecfg = EncodeConfig { min_support: 3 };
        let all = Bits::ones(n);
        let mut want = Vec::new();
        for scheme in SCHEMES {
            let path = dir
                .path()
                .join("plans")
                .join(format!("{}.json", plan_name(scheme, 0)));
            let file: PlanFile = serde_json::from_slice(&std::fs::read(path).unwrap()).unwrap();
            let p = ArmPlan {
                fold: Some(file.plan),
                cohort: None,
                rows_sha256: String::new(),
                arm_sha256: None,
                cluster: vec![],
                repeat: 0,
            };
            let view = fold_view(&p, &cfg);
            assert!(view.unavailable_reason.is_none(), "{scheme}: no plan");
            let mut shows = false;
            for f in 0..view.k {
                let train = view.train(f);
                assert!(
                    train.count() > 0 && (train.count() as usize) < n,
                    "{scheme} fold {f}"
                );
                let rows = |mask: &Bits| fold_matrix(&raws, &samples, mask, &ecfg).unwrap().rows;
                shows |= rows(&train) < rows(&all);
                want.push(train);
            }
            assert!(
                shows,
                "{scheme}: no fold where a matrix on every sample differs"
            );
        }
        assert_eq!(got.len(), want.len(), "fold matrices built");
        for (i, (g, w)) in got.iter().zip(&want).enumerate() {
            assert!(
                g == w,
                "fold matrix {i} got {} samples, not the fold's {} training samples",
                g.count(),
                w.count()
            );
        }
        let runs: Vec<serde_json::Value> =
            serde_json::from_slice(&std::fs::read(dir.path().join("runs.json")).unwrap()).unwrap();
        assert_eq!(runs.len(), 2 * crate::learn::eval::LEARNERS.len());
        assert!(runs.iter().all(|r| r["status"] == "ok"), "{runs:?}");
    }

    /// `run_arm` scores a cohort plan with its held-out genomes out of every
    /// fold's training (the reviewer's M17), and with the arm's learner seed.
    #[test]
    fn a_cohort_plan_keeps_its_held_out_genomes_in_run_arm() {
        let n = 40;
        let lineage: Vec<u32> = (0..n as u32).map(|i| i / 4).collect();
        let y = Bits::from_indices(n, (0..n).filter(|i| i % 2 == 0));
        let ids: Vec<String> = (0..n).map(|i| format!("g{i}")).collect();
        let groups = GroupIds::new(lineage.clone(), vec![], None).unwrap();
        let reference = fold_plan(&y, &groups, Scheme::Random, 5, 3, 42, 0).unwrap();
        // The arm adds a copy of g0 (sample 40), which must never train in g0's fold.
        let mut arm_ids = ids.clone();
        arm_ids.push("copy".into());
        let mut arm_lineage = lineage.clone();
        arm_lineage.push(lineage[0]);
        let arm_y = Bits::from_indices(n + 1, (0..=n).filter(|&i| i == n || y.contains(i)));
        let arm_groups = GroupIds::new(arm_lineage, vec![vec![0, n]], None).unwrap();
        let species = SpeciesGroups::of_cohorts(&[(&ids, &groups), (&arm_ids, &arm_groups)]);
        let cohort =
            plan::cohort_plan(&reference, &ids, &y, &arm_ids, &arm_y, &species, None).unwrap();
        assert!(cohort.held_out.iter().any(|h| !h.is_empty()));
        let p = ArmPlan {
            fold: None,
            cohort: Some(cohort.clone()),
            rows_sha256: reference.sha256.clone(),
            arm_sha256: Some(cohort.sha256.clone()),
            cluster: vec![None; n + 1],
            repeat: 0,
        };
        let cfg = ArmConfig {
            arm: "a".into(),
            blocks: vec!["K".into()],
            eval: EvalConfig::default(),
            repeats: 1,
            pairing: Pairing::CommonOuterTestSetBalanced,
            reference: None,
            balance_seed: None,
            learner_seed: Some(5),
            tuning_seed: None,
            trace_training: false,
            groups: None,
        };
        let view = fold_view(&p, &cfg);
        assert_eq!(view.learner_seed, Some(5));
        for f in 0..cohort.k {
            assert_eq!(view.train(f), cohort.train(f), "outer fold {f}");
            for &i in &cohort.held_out[f] {
                assert!(
                    !view.train(f).contains(i),
                    "outer fold {f} trains held-out sample {i}"
                );
            }
        }
        // Keyed by the plan's repeat: repeat 0 keeps the seed, another draws anew.
        let later = ArmPlan { repeat: 2, ..p };
        assert_eq!(
            fold_view(&later, &cfg).learner_seed,
            Some(plan::repeat_seed(5, 2))
        );
        assert_ne!(plan::repeat_seed(5, 2), 5);
    }

    /// 60 lineages of 4 genomes, every fifth lineage resistant, so some keyed
    /// lineage draws leave an outer fold or a holdout without a resistant
    /// lineage and others don't.
    fn patchy() -> ArmInputs {
        let n = 240;
        ArmInputs {
            sample_ids: (0..n).map(|i| format!("g{i:03}")).collect(),
            columns: (0..n).collect(),
            y: Bits::from_indices(n, (0..n).filter(|i| (i / 4) % 5 == 0)),
            min_support: 1,
            x: crate::bits::BitMatrix::with_cols(n),
            lineage: (0..n as u32).map(|i| i / 4).collect(),
            duplicates: vec![],
        }
    }

    fn reference_cfg(seed: u64, repeats: u32) -> ArmConfig {
        ArmConfig {
            arm: "default".into(),
            blocks: vec!["K".into()],
            eval: EvalConfig {
                seed,
                ..EvalConfig::default()
            },
            repeats,
            pairing: Pairing::SamePlan,
            reference: None,
            balance_seed: None,
            learner_seed: None,
            tuning_seed: None,
            trace_training: false,
            groups: None,
        }
    }

    /// Main's ruling 2 after VERIFY-H1STATS: a refused lineage repeat is
    /// replaced by the next keyed repeat past the slots, at most twice per
    /// split, each refused draw kept and listed; past that the split has
    /// fewer valid repeats than asked and says so. The random split repeats once.
    #[test]
    fn a_refused_lineage_repeat_is_replaced_by_the_next_keyed_draw_at_most_twice() {
        let a = patchy();
        let groups = GroupIds::new(a.lineage.clone(), vec![], None).unwrap();
        let ok = |seed: u64, repeat: u32| {
            plan::keyed_fold_plan(
                &a.sample_ids,
                &a.y,
                &groups,
                Scheme::Lineage,
                5,
                3,
                seed,
                repeat,
            )
            .unwrap()
            .is_available()
        };
        // The first seed whose keyed repeats 0, 1, … are available as `want` says.
        let find = |want: &[bool]| {
            (0..20000u64)
                .find(|&s| (0..want.len()).all(|r| ok(s, r as u32) == want[r]))
                .expect("a seed with this pattern")
        };
        // Slot 1 refused, replaced by keyed repeat 3.
        let seed = find(&[true, false, true, true]);
        let dir = tempfile::tempdir().unwrap();
        let cfg = reference_cfg(seed, 3);
        let (plans, draws, short) = arm_plans(&cfg, &a, None, Scheme::Lineage, dir.path()).unwrap();
        assert_eq!(short, None);
        assert_eq!(draws.drawn, vec![0, 3, 2]);
        assert_eq!(
            draws.refused.iter().map(|d| d.repeat).collect::<Vec<_>>(),
            vec![1]
        );
        assert!(
            plans
                .iter()
                .all(|p| p.fold.as_ref().unwrap().is_available())
        );
        let slot1: PlanFile =
            serde_json::from_slice(&std::fs::read(dir.path().join("lineage-1.json")).unwrap())
                .unwrap();
        assert_eq!(slot1.plan.repeat, 3);
        let kept: PlanFile = serde_json::from_slice(
            &std::fs::read(dir.path().join("refused-lineage-1.json")).unwrap(),
        )
        .unwrap();
        assert!(kept.plan.unavailable_reason.is_some());
        // Slot 1 and both replacements refused: the split is short, named.
        let seed = find(&[true, false, true, false, false]);
        let dir = tempfile::tempdir().unwrap();
        let cfg = reference_cfg(seed, 3);
        let (_, draws, short) = arm_plans(&cfg, &a, None, Scheme::Lineage, dir.path()).unwrap();
        assert_eq!(draws.drawn, vec![0, 4, 2]);
        assert_eq!(draws.refused.len(), MAX_REPLACEMENTS as usize);
        let why = short.expect("short of valid repeats");
        assert!(
            why.contains("2 replacement draws spent, fewer than 3 valid lineage repeats"),
            "{why}"
        );
        // The random split scores one repeat and is never replaced; an arm
        // with a reference takes its slot plans as they are.
        let (plans, draws, _) = arm_plans(&cfg, &a, None, Scheme::Random, dir.path()).unwrap();
        assert_eq!((plans.len(), draws.drawn.clone()), (1, vec![0]));
        let shared = tempfile::tempdir().unwrap();
        std::fs::create_dir(shared.path().join("plans")).unwrap();
        for name in ["lineage-0.json", "lineage-1.json", "lineage-2.json"] {
            std::fs::copy(
                dir.path().join(name),
                shared.path().join("plans").join(name),
            )
            .unwrap();
        }
        let twin = ArmConfig {
            arm: "twin".into(),
            reference: Some(shared.path().to_path_buf()),
            ..reference_cfg(seed, 3)
        };
        let out = tempfile::tempdir().unwrap();
        let (plans, draws, short) =
            arm_plans(&twin, &a, None, Scheme::Lineage, out.path()).unwrap();
        assert_eq!(
            (draws.drawn, draws.refused.len(), short),
            (vec![0, 4, 2], 0, None)
        );
        assert!(!plans[1].fold.as_ref().unwrap().is_available());
    }

    /// Main's ruling after the first calibration failed: a retuned twin whose
    /// tuning draw is refused draws again, keyed by draw, at most twice per
    /// split, keeping each refused draw; past that its split is unavailable.
    #[test]
    fn a_refused_retuned_draw_is_replaced_at_most_twice_per_split() {
        let a = patchy();
        let groups = GroupIds::new(a.lineage.clone(), vec![], None).unwrap();
        let reference_ok = |seed: u64| {
            plan::keyed_fold_plan(&a.sample_ids, &a.y, &groups, Scheme::Lineage, 5, 3, seed, 0)
                .unwrap()
                .is_available()
        };
        let seed = (0..20000u64).find(|&s| reference_ok(s)).unwrap();
        let r = plan::keyed_fold_plan(&a.sample_ids, &a.y, &groups, Scheme::Lineage, 5, 3, seed, 0)
            .unwrap();
        let ok = |tuning: u64, draw: u32| {
            plan::retuned_plan(&r, &a.sample_ids, &a.y, &groups, tuning, draw)
                .unwrap()
                .is_available()
        };
        let reference = tempfile::tempdir().unwrap();
        let plans = reference.path().join("plans");
        std::fs::create_dir(&plans).unwrap();
        arm_plans(&reference_cfg(seed, 1), &a, None, Scheme::Lineage, &plans).unwrap();
        let twin = |tuning: u64| ArmConfig {
            arm: "twin".into(),
            reference: Some(reference.path().to_path_buf()),
            tuning_seed: Some(tuning),
            ..reference_cfg(seed, 1)
        };
        // Draw 0 refused, draw 1 accepted.
        let tuning = (0..20000u64).find(|&t| !ok(t, 0) && ok(t, 1)).unwrap();
        let out = tempfile::tempdir().unwrap();
        let (got, draws, _) =
            arm_plans(&twin(tuning), &a, None, Scheme::Lineage, out.path()).unwrap();
        assert!(got[0].fold.as_ref().unwrap().is_available());
        assert_eq!(draws.refused_tuning.len(), 1);
        assert_eq!(
            (draws.refused_tuning[0].slot, draws.refused_tuning[0].draw),
            (0, 0)
        );
        assert!(
            out.path()
                .join("refused-lineage-0-tuning-0.retuned.json")
                .exists()
        );
        let kept: serde_json::Value = serde_json::from_slice(
            &std::fs::read(out.path().join("lineage-0.retuned.json")).unwrap(),
        )
        .unwrap();
        assert_eq!(kept["draw"], 1);
        // Draws 0, 1 and 2 refused: two replacements spent, the split unavailable.
        let tuning = (0..20000u64)
            .find(|&t| !ok(t, 0) && !ok(t, 1) && !ok(t, 2))
            .unwrap();
        let out = tempfile::tempdir().unwrap();
        let (got, draws, _) =
            arm_plans(&twin(tuning), &a, None, Scheme::Lineage, out.path()).unwrap();
        assert_eq!(draws.refused_tuning.len(), MAX_REPLACEMENTS as usize);
        assert!(!got[0].fold.as_ref().unwrap().is_available());
    }

    /// Main's ruling 1 after VERIFY-H1STATS: a cohort arm's balance seed is
    /// keyed by the reference's draw, so its subsample is drawn anew in each
    /// repeat, as a learner seed is.
    #[test]
    fn a_cohort_arm_s_balance_seed_is_keyed_by_repeat() {
        let a = patchy();
        let groups = GroupIds::new(a.lineage.clone(), vec![], None).unwrap();
        let species = SpeciesGroups::of_cohorts(&[(&a.sample_ids, &groups)]);
        let available = |seed: u64| {
            (0..3).all(|r| {
                plan::keyed_fold_plan(&a.sample_ids, &a.y, &groups, Scheme::Lineage, 5, 3, seed, r)
                    .unwrap()
                    .is_available()
            })
        };
        let seed = (0..20000u64).find(|&s| available(s)).unwrap();
        let reference = tempfile::tempdir().unwrap();
        let plans = reference.path().join("plans");
        std::fs::create_dir(&plans).unwrap();
        arm_plans(
            &reference_cfg(seed, 3),
            &a,
            Some(&species),
            Scheme::Lineage,
            &plans,
        )
        .unwrap();
        let cohort = ArmConfig {
            arm: "cohort-path".into(),
            pairing: Pairing::CommonOuterTestSetBalanced,
            reference: Some(reference.path().to_path_buf()),
            balance_seed: Some(7),
            ..reference_cfg(seed, 3)
        };
        let out = tempfile::tempdir().unwrap();
        let (got, _, _) =
            arm_plans(&cohort, &a, Some(&species), Scheme::Lineage, out.path()).unwrap();
        let seeds: Vec<u64> = got
            .iter()
            .map(|p| p.cohort.as_ref().unwrap().balance_seed)
            .collect();
        assert_eq!(
            seeds,
            (0..3).map(|r| plan::repeat_seed(7, r)).collect::<Vec<_>>()
        );
        assert_eq!(seeds[0], 7);
    }
}
