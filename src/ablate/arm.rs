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
    /// Write `training.jsonl`: a first line with the arm's sample ids, then per
    /// split, repeat and outer fold the genome ids `run_arm` handed the
    /// learners as training samples.
    pub trace_training: bool,
}

/// The splits an arm scores, in `evaluate`'s order.
pub const SCHEMES: [Scheme; 2] = [Scheme::Random, Scheme::Lineage];

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
/// plan's held-out genomes included), with the arm's learner seed.
fn fold_view<'a>(p: &'a ArmPlan, cfg: &ArmConfig) -> PlanView<'a> {
    let mut view = p.view();
    view.learner_seed = cfg.learner_seed;
    view
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
    let mut candidates = Vec::with_capacity(raws.len());
    for raw in raws {
        let enc = encode::fit_encoder(raw, samples, train, ecfg)?;
        let block = encode::transform(&enc, raw, samples)?;
        candidates.push(encode::block_candidates(&block, train, ecfg));
    }
    Ok(encode::assemble(&candidates, train).x)
}

fn plan_name(scheme: Scheme, repeat: u32) -> String {
    format!("{scheme}-{repeat}")
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

fn arm_plan(
    cfg: &ArmConfig,
    a: &ArmInputs,
    species: Option<&SpeciesGroups>,
    scheme: Scheme,
    repeat: u32,
    plans_dir: &Path,
) -> Result<ArmPlan> {
    let name = plan_name(scheme, repeat);
    let Some(reference) = &cfg.reference else {
        ensure!(
            cfg.pairing == Pairing::SamePlan,
            "the reference arm pairs with itself by same_plan, not {}",
            cfg.pairing.as_str()
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
    match cfg.pairing {
        Pairing::SamePlan | Pairing::K31Groups => {
            ensure!(
                r.sample_ids == a.sample_ids && r.y == labels(&a.y),
                "pairing {} needs the reference's samples and labels in its order; this arm has {} samples, the reference {}",
                cfg.pairing.as_str(),
                a.sample_ids.len(),
                r.sample_ids.len()
            );
            Ok(ArmPlan {
                rows_sha256: r.plan.sha256.clone(),
                arm_sha256: None,
                cluster: r.cluster.iter().map(|&c| Some(c)).collect(),
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
                cfg.balance_seed,
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
    let learners = crate::learn::eval::LEARNERS;
    for scheme in SCHEMES {
        let mut results = Vec::new();
        let mut reason = unavailable
            .first()
            .map(|(id, why)| format!("block {id} is unavailable: {why}"));
        let started = Instant::now();
        for repeat in 0..cfg.repeats {
            let p = arm_plan(cfg, &a, species.as_ref(), scheme, repeat, &plans_dir)?;
            if reason.is_some() {
                results.push((p, None));
                continue;
            }
            let view = fold_view(&p, cfg);
            let raws: Vec<&blocks::RawBlock> = admitted.iter().map(|b| b.raw()).collect();
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
                repeats: cfg.repeats,
                fold_plan_sha256: Some(rows_sha.clone()),
                status: "ok".into(),
                reason: None,
                metrics: serde_json::Map::new(),
                fit_count: cfg.eval.folds * cfg.repeats as usize,
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
    }
}
