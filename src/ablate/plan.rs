//! Fold plans for feature-block ablations, and the fits they cost.
//!
//! A fold plan fixes every partition one (dataset, split, repeat) is scored
//! on: the outer folds, each outer fold's inner folds, and its two tuning
//! holdouts. Every ablation arm of that dataset, split and repeat uses the
//! same plan, so arms differ only in their feature blocks. Repeat 0 with the
//! pipeline's seed has the outer folds `learn::eval::evaluate` scores at
//! c90ebaf; later repeats draw from a seed derived from the base seed
//! (`repeat_seed`).
//!
//! Inner folds and tuning holdouts follow the duplicate-safe policy of
//! `tuning_groups`: a lineage split keeps every lineage whole, and a random
//! split keeps every k-mer-identical duplicate group whole.
//!
//! `FoldPlan::sha256` identifies a plan across the Rust and Python tasks;
//! `canonical_json` documents its exact bytes. `plan_compute` counts the fits
//! a protocol would run and refuses one that exceeds a tier's budget.

use crate::bits::Bits;
use crate::learn::{
    fold_masks, grouped_folds, isolate_groups, scm, stratified_folds, training_folds,
    tuning_holdout, validate_folds,
};
use anyhow::{Context, Result, bail, ensure};
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::fmt::Write as _;

/// The fraction of outer-training samples boosting holds out to stop early:
/// `learn::eval`'s own share, so plans and evaluations can't drift apart.
pub const GBT_HOLDOUT: f64 = crate::learn::eval::GBT_HOLDOUT;
/// The fraction of outer-training samples the lasso holds out to choose its
/// penalty: `learn::eval`'s own share.
pub const LASSO_HOLDOUT: f64 = crate::learn::eval::LASSO_HOLDOUT;
/// The `format` of every plan this module writes. It is part of the hashed
/// bytes, so a change of layout changes every hash.
pub const PLAN_FORMAT: &str = "grm-fold-plan-v1";
/// Outer fold `f` seeds its inner draws with `repeat_seed + 101 f`, as
/// `learn::eval::fold_seed` does with the pipeline's seed.
const FOLD_SEED_STRIDE: u64 = 101;

/// How outer folds are drawn.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Scheme {
    /// Stratified k-fold, each duplicate group moved whole into one fold.
    Random,
    /// Whole lineages held out together.
    Lineage,
}

impl Scheme {
    pub fn as_str(self) -> &'static str {
        match self {
            Scheme::Random => "random",
            Scheme::Lineage => "lineage",
        }
    }
}

impl std::fmt::Display for Scheme {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}

/// How samples are related: lineage clusters, k-mer-identical duplicate
/// groups in `kset_duplicates` format (two or more ascending sample indices
/// per group, no sample in two groups), and optionally a source such as the
/// submitting study.
///
/// These ids decide splits, inner folds and inference units only. They are
/// never predictors: the type has no conversion to a feature matrix or to
/// any block type, and its fields are private, readable only by reference.
/// The two split schemes read `lineage` and `duplicate`; `source` is carried
/// for inference units and later schemes.
///
/// ```
/// use grm::ablate::plan::GroupIds;
/// let groups = GroupIds::new(vec![0, 0, 1], vec![vec![0, 1]], None).unwrap();
/// assert_eq!(groups.lineage(), &[0, 0, 1]);
/// ```
///
/// ```compile_fail
/// use grm::ablate::plan::GroupIds;
/// use grm::bits::BitMatrix;
/// let groups = GroupIds::new(vec![0, 0, 1], vec![vec![0, 1]], None).unwrap();
/// let _: BitMatrix = groups.into();
/// ```
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GroupIds {
    lineage: Vec<u32>,
    duplicate: Vec<Vec<usize>>,
    source: Option<Vec<u32>>,
}

impl GroupIds {
    /// Checks the ids as `learn::eval::evaluate` checks its inputs: one
    /// lineage per sample, and each duplicate group two or more ascending
    /// in-range indices, disjoint from the others, inside one lineage. A
    /// source, when given, has one id per sample.
    pub fn new(
        lineage: Vec<u32>,
        duplicate: Vec<Vec<usize>>,
        source: Option<Vec<u32>>,
    ) -> Result<Self> {
        let n = lineage.len();
        ensure!(
            u32::try_from(n).is_ok(),
            "{n} samples is more than u32 sample ids can name"
        );
        let mut grouped = Bits::zeros(n);
        for group in &duplicate {
            ensure!(
                group.len() >= 2
                    && group.windows(2).all(|w| w[0] < w[1])
                    && group[group.len() - 1] < n,
                "a duplicate group needs two or more ascending sample indices below {n}: {group:?}"
            );
            for &i in group {
                ensure!(
                    !grouped.contains(i),
                    "sample {i} is in two duplicate groups"
                );
                grouped.insert(i);
                ensure!(
                    lineage[i] == lineage[group[0]],
                    "duplicate samples {} and {i} are in different lineages",
                    group[0]
                );
            }
        }
        if let Some(source) = &source {
            ensure!(
                source.len() == n,
                "{} source ids for {n} samples",
                source.len()
            );
        }
        Ok(GroupIds {
            lineage,
            duplicate,
            source,
        })
    }

    pub fn len(&self) -> usize {
        self.lineage.len()
    }

    pub fn is_empty(&self) -> bool {
        self.lineage.is_empty()
    }

    pub fn lineage(&self) -> &[u32] {
        &self.lineage
    }

    pub fn duplicate(&self) -> &[Vec<usize>] {
        &self.duplicate
    }

    pub fn source(&self) -> Option<&[u32]> {
        self.source.as_deref()
    }
}

/// The groups an outer fold's inner folds and tuning holdouts keep whole, or
/// None for the plain stratified draws:
/// - lineage: the lineages, as at c90ebaf;
/// - random without duplicates: None, the c90ebaf stratified path;
/// - random with duplicates: composite ids, where every member of a
///   duplicate group takes its first member's index and every other sample
///   keeps its own index.
///
/// This is `learn::eval::tuning_groups`, the ids `evaluate` tunes on, so
/// plans and evaluations can't diverge. The first member is the lowest when
/// the group is ascending, as `GroupIds::new` and `evaluate` require, but
/// this function takes raw slices and doesn't check: an unsorted group takes
/// its first member's index, not its lowest.
pub fn tuning_groups(
    scheme: Scheme,
    lineage: &[u32],
    duplicates: &[Vec<usize>],
) -> Option<Vec<u32>> {
    crate::learn::eval::tuning_groups(scheme.as_str(), lineage, duplicates)
        .map(std::borrow::Cow::into_owned)
}

/// The seed repeat `repeat` draws its outer folds from. Repeat 0 uses `seed`
/// itself, so it reproduces the pipeline. Repeat r > 0 uses the first eight
/// bytes, read big-endian, of the SHA-256 of the ASCII text
/// `grm-fold-plan-repeat:{seed}:{r}` with both numbers in decimal; in Python,
/// `int.from_bytes(hashlib.sha256(f"grm-fold-plan-repeat:{seed}:{r}".encode()).digest()[:8], "big")`.
/// A hash, unlike `seed + r`, cannot hand one repeat another repeat's
/// per-fold seeds (`fold_seed`).
pub fn repeat_seed(seed: u64, repeat: u32) -> u64 {
    if repeat == 0 {
        return seed;
    }
    let tag = format!("grm-fold-plan-repeat:{seed}:{repeat}");
    let digest = ring::digest::digest(&ring::digest::SHA256, tag.as_bytes());
    let mut head = [0u8; 8];
    head.copy_from_slice(&digest.as_ref()[..8]);
    u64::from_be_bytes(head)
}

/// The seed of outer fold `fold`'s inner folds and both tuning holdouts:
/// `repeat_seed + 101 * fold`, wrapping at 2^64. With the pipeline's seed at
/// repeat 0 this is `learn::eval::fold_seed`.
pub fn fold_seed(repeat_seed: u64, fold: usize) -> u64 {
    repeat_seed.wrapping_add((fold as u64).wrapping_mul(FOLD_SEED_STRIDE))
}

/// One outer fold's tuning holdouts as their validation samples, ascending.
/// The fitting samples are the rest of that outer fold's training samples.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TuningHoldouts {
    /// `GBT_HOLDOUT` (0.15): boosting's early-stopping samples.
    pub gbt: Vec<usize>,
    /// `LASSO_HOLDOUT` (0.2): the lasso's penalty-selection samples.
    pub lasso: Vec<usize>,
}

/// Every partition one (dataset, split, repeat) is scored on. Ablation arms
/// share it. Its serde JSON is the plan file; `sha256` covers every other
/// key of it (`canonical_json`). A plan file with a key this type does not
/// name is refused, so Rust and Python hash the same keys.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FoldPlan {
    /// `PLAN_FORMAT`.
    pub format: String,
    pub scheme: Scheme,
    pub repeat: u32,
    /// Outer folds.
    pub k: usize,
    /// Inner folds per outer fold.
    pub inner_k: usize,
    /// The base seed, `EvalConfig::seed` for the pipeline's own plan.
    pub seed: u64,
    /// `repeat_seed(seed, repeat)`, which drew `outer`.
    pub repeat_seed: u64,
    /// Each sample's outer test fold.
    pub outer: Vec<usize>,
    /// Per outer fold, each sample's inner fold, None for that outer fold's
    /// test samples. Empty when the plan is unavailable.
    pub inner: Vec<Vec<Option<usize>>>,
    /// Per outer fold, its tuning holdouts. Empty when the plan is
    /// unavailable.
    pub holdouts: Vec<TuningHoldouts>,
    /// Why the nested plan cannot be scored, or None when it can.
    pub unavailable_reason: Option<String>,
    /// Lowercase hex SHA-256 of `canonical_json`.
    pub sha256: String,
}

impl FoldPlan {
    /// True when the plan can be scored: it has no unavailable reason and
    /// the shape `check_structure` requires.
    pub fn is_available(&self) -> bool {
        self.unavailable_reason.is_none() && self.check_structure().is_ok()
    }

    /// Fails unless the plan has the shape `fold_plan` writes:
    /// - `format` is `PLAN_FORMAT`, both fold counts are at least two, and
    ///   `repeat_seed` is `repeat_seed(seed, repeat)`;
    /// - every outer fold is below `k`;
    /// - an unavailable plan has a non-empty reason and no inner folds or
    ///   holdouts;
    /// - an available plan has `k` inner plans and `k` holdout pairs. In
    ///   outer fold `f`, a sample has an inner fold below `inner_k` exactly
    ///   when it trains in `f`, and both holdouts are ascending training
    ///   samples of `f`.
    pub fn check_structure(&self) -> Result<()> {
        ensure!(
            self.format == PLAN_FORMAT || self.format == KEYED_PLAN_FORMAT,
            "fold plan format {:?} is not {PLAN_FORMAT} or {KEYED_PLAN_FORMAT}",
            self.format
        );
        ensure!(
            self.k >= 2 && self.inner_k >= 2,
            "fold plan needs at least two outer and two inner folds"
        );
        ensure!(
            self.repeat_seed == repeat_seed(self.seed, self.repeat),
            "fold plan repeat_seed does not follow from its seed and repeat"
        );
        let n = self.outer.len();
        ensure!(
            self.outer.iter().all(|&f| f < self.k),
            "fold plan has an outer fold outside 0..{}",
            self.k
        );
        if let Some(reason) = &self.unavailable_reason {
            ensure!(
                !reason.is_empty(),
                "unavailable fold plan has an empty reason"
            );
            ensure!(
                self.inner.is_empty() && self.holdouts.is_empty(),
                "unavailable fold plan has inner folds or holdouts"
            );
            return Ok(());
        }
        ensure!(
            self.inner.len() == self.k && self.holdouts.len() == self.k,
            "available fold plan needs {} inner plans and holdout pairs, not {} and {}",
            self.k,
            self.inner.len(),
            self.holdouts.len()
        );
        for (f, (inner, held)) in self.inner.iter().zip(&self.holdouts).enumerate() {
            ensure!(
                inner.len() == n,
                "outer fold {f}'s inner plan covers {} of {n} samples",
                inner.len()
            );
            for (i, fold) in inner.iter().enumerate() {
                let trains = self.outer[i] != f;
                ensure!(
                    match fold {
                        Some(j) => trains && *j < self.inner_k,
                        None => !trains,
                    },
                    "outer fold {f}'s inner fold for sample {i} does not match its outer fold"
                );
            }
            for (name, val) in [("gbt", &held.gbt), ("lasso", &held.lasso)] {
                ensure!(
                    val.windows(2).all(|w| w[0] < w[1])
                        && val.iter().all(|&i| i < n && self.outer[i] != f),
                    "outer fold {f}'s {name} holdout is not ascending training samples"
                );
            }
        }
        Ok(())
    }

    /// The bytes `sha256` covers: this plan's serde JSON without its
    /// `sha256` key, in the canonical form `canonical_json` defines. In
    /// Python, with `plan` the parsed plan file:
    /// `json.dumps({k: v for k, v in plan.items() if k != "sha256"}, sort_keys=True, separators=(",", ":"), ensure_ascii=True)`.
    pub fn canonical_json(&self) -> Result<String> {
        let mut value = serde_json::to_value(self)?;
        if let Value::Object(map) = &mut value {
            map.remove("sha256");
        }
        canonical_json(&value)
    }

    /// The hash `sha256` should hold for the rest of the plan.
    pub fn compute_sha256(&self) -> Result<String> {
        Ok(sha256_hex(self.canonical_json()?.as_bytes()))
    }

    /// Fails when the plan's shape is wrong (`check_structure`) or `sha256`
    /// does not match the rest of the plan.
    pub fn verify(&self) -> Result<()> {
        self.check_structure()?;
        ensure!(
            self.compute_sha256()? == self.sha256,
            "fold plan sha256 does not match its contents"
        );
        Ok(())
    }
}

/// The fold plan of one split scheme and repeat.
///
/// Outer folds: random is `isolate_groups(stratified_folds(y, k, s),
/// duplicates, k)` and lineage is `grouped_folds(y, lineage, k, s)`, with
/// `s = repeat_seed(seed, repeat)`. At repeat 0 with `EvalConfig::default()`
/// (k 5, inner_k 3, seed 42) these are the outer folds `evaluate` scores.
///
/// Inner folds and tuning holdouts of outer fold `f` are drawn from its
/// training samples alone, seeded `fold_seed(s, f)`, by `training_folds` and
/// `tuning_holdout` with `tuning_groups(scheme, ..)`. Outer-test labels
/// cannot change them.
///
/// A plan that cannot be scored is still returned, with its outer folds,
/// empty inner folds and holdouts, and `unavailable_reason`: a lineage split
/// with fewer than 2k lineages, any outer, inner or holdout partition with a
/// side lacking a class, in the order and wording of
/// `learn::eval::check_split_plan`. Errors are kept for inputs that are
/// wrong rather than unlucky: mismatched lengths or fewer than two folds.
pub fn fold_plan(
    y: &Bits,
    groups: &GroupIds,
    scheme: Scheme,
    k: usize,
    inner_k: usize,
    seed: u64,
    repeat: u32,
) -> Result<FoldPlan> {
    ensure!(
        y.len() == groups.len(),
        "fold plan dimensions differ: {} labels, {} samples with group ids",
        y.len(),
        groups.len()
    );
    ensure!(
        k >= 2 && inner_k >= 2,
        "cross-validation needs at least two folds"
    );
    let drawn = repeat_seed(seed, repeat);
    let outer = match scheme {
        Scheme::Random => isolate_groups(&stratified_folds(y, k, drawn), &groups.duplicate, k),
        Scheme::Lineage => grouped_folds(y, &groups.lineage, k, drawn),
    };
    let (inner, holdouts, unavailable_reason) =
        match nested_plan(y, groups, scheme, &outer, k, inner_k, drawn) {
            Ok((inner, holdouts)) => (inner, holdouts, None),
            Err(error) => (vec![], vec![], Some(format!("{error:#}"))),
        };
    let mut plan = FoldPlan {
        format: PLAN_FORMAT.to_string(),
        scheme,
        repeat,
        k,
        inner_k,
        seed,
        repeat_seed: drawn,
        outer,
        inner,
        holdouts,
        unavailable_reason,
        sha256: String::new(),
    };
    plan.sha256 = plan.compute_sha256()?;
    Ok(plan)
}

/// The `format` of every cohort plan.
pub const COHORT_PLAN_FORMAT: &str = "grm-cohort-plan-v4";

/// The plan of an arm whose cohort differs from the reference arm's, such
/// as another cap or selection rule (DECISIONS V-1 and the rulings after
/// lane INPUT: the `common_outer_test_set_balanced` pairing).
///
/// Per outer fold f, the arm is scored on a class-balanced subsample of the
/// genomes in both cohorts that the reference tests in f: every class keeps
/// as many genomes as it can while the subsample's resistant-to-susceptible
/// ratio matches the reference's whole test fold f (`reference_test_classes`),
/// so a cohort that the selection rule tilted towards one class is scored on
/// the population the reference represents. The subsample is drawn by a
/// ChaCha8 shuffle of each class, in sample order, seeded `balance_seed(
/// balance_seed, f)` (`balance_seed` is the reference's repeat seed) and
/// ranked by `keyed_hash(genome id, that seed, "balance")`, so the draw is keyed on genome
/// identity.
///
/// In fold f the arm trains on its own cohort minus `held_out[f]` and minus
/// the scored genomes: `held_out[f]` holds the common genomes the reference
/// tests in f that the subsample left out (`balance_dropped[f]` of them),
/// and, on the lineage split, every genome only the arm has whose lineage
/// (the arm's own k = 31 lineages) is a lineage of fold f's scored genomes
/// (`lineage_dropped[f]`), and, on every split, every other genome of a
/// scored genome's k-mer duplicate group (`duplicate_dropped[f]`), so no exact
/// copy of a scored genome trains in its fold. A common genome the subsample leaves out is never
/// scored (`outer` None), and a genome only the arm has trains in every other
/// fold. Inner folds and both tuning holdouts are keyed (`keyed_tuning`),
/// seeded `fold_seed(repeat_seed, f)`, on unit keys taken through the
/// reference (`SpeciesGroups::unit_keys`). The arm's lineages and duplicate
/// groups are the species grouping's, not its own. So a genome both cohorts
/// hold has the reference arm's inner folds and holdouts in every outer fold
/// whose training samples it is in, whichever other genomes the arm holds.
/// Which genomes it trains beside, and so its scores, still change.
///
/// `sha256` covers every other key (`canonical_json`, as for `FoldPlan`),
/// and `reference_sha256` names the reference plan, so the arm's prediction
/// rows carry the reference plan's hash (their `fold` is the reference's)
/// while its run records this plan's.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CohortPlan {
    /// `COHORT_PLAN_FORMAT`.
    pub format: String,
    /// The reference `FoldPlan`'s `sha256`.
    pub reference_sha256: String,
    pub scheme: Scheme,
    pub repeat: u32,
    pub k: usize,
    pub inner_k: usize,
    pub seed: u64,
    pub repeat_seed: u64,
    /// The seed of the class-balanced subsample: the reference's repeat seed
    /// unless `cohort_plan` was given another.
    pub balance_seed: u64,
    /// Per outer fold, the reference's whole test fold as [resistant,
    /// susceptible] counts: the balance the subsample matches.
    pub reference_test_classes: Vec<[usize; 2]>,
    /// Each sample's outer test fold: the fold that scores it, or None for a
    /// genome that is never scored.
    pub outer: Vec<Option<usize>>,
    /// Per outer fold, the samples that neither train nor are scored there,
    /// ascending.
    pub held_out: Vec<Vec<usize>>,
    /// Per outer fold, the common test genomes the balanced subsample left out.
    pub balance_dropped: Vec<usize>,
    /// Per outer fold, the arm-only genomes dropped from training for sharing
    /// a lineage with a scored genome (lineage split only; 0 on random).
    pub lineage_dropped: Vec<usize>,
    /// Per outer fold, the genomes dropped from training for being in the
    /// k-mer duplicate group of a scored genome (every split).
    pub duplicate_dropped: Vec<usize>,
    /// As `FoldPlan::inner`: per outer fold, each training sample's inner
    /// fold.
    pub inner: Vec<Vec<Option<usize>>>,
    pub holdouts: Vec<TuningHoldouts>,
    pub unavailable_reason: Option<String>,
    pub sha256: String,
}

/// The seed of outer fold `fold`'s class-balanced subsample: the first eight
/// bytes, big-endian, of the SHA-256 of `grm-cohort-balance:{seed}:{fold}`.
pub fn balance_seed(seed: u64, fold: usize) -> u64 {
    let tag = format!("grm-cohort-balance:{seed}:{fold}");
    let digest = ring::digest::digest(&ring::digest::SHA256, tag.as_bytes());
    let mut head = [0u8; 8];
    head.copy_from_slice(&digest.as_ref()[..8]);
    u64::from_be_bytes(head)
}

impl CohortPlan {
    /// Outer fold `f`'s training samples: every sample neither scored in `f`
    /// nor held out of it.
    pub fn train(&self, f: usize) -> Bits {
        let n = self.outer.len();
        let mut train = Bits::from_indices(n, (0..n).filter(|&i| self.outer[i] != Some(f)));
        for &i in self.held_out.get(f).map(Vec::as_slice).unwrap_or(&[]) {
            train.remove(i);
        }
        train
    }

    /// Fails unless `sha256` matches the rest of the plan and its shape is
    /// the one `cohort_plan` writes.
    pub fn verify(&self) -> Result<()> {
        ensure!(
            self.format == COHORT_PLAN_FORMAT,
            "cohort plan format {:?} is not {COHORT_PLAN_FORMAT}",
            self.format
        );
        ensure!(
            self.repeat_seed == repeat_seed(self.seed, self.repeat),
            "cohort plan repeat_seed does not follow from its seed and repeat"
        );
        ensure!(
            self.outer.iter().flatten().all(|&f| f < self.k),
            "cohort plan has an outer fold outside 0..{}",
            self.k
        );
        let n = self.outer.len();
        ensure!(
            self.held_out.len() == self.k
                && self.balance_dropped.len() == self.k
                && self.lineage_dropped.len() == self.k
                && self.duplicate_dropped.len() == self.k
                && self.reference_test_classes.len() == self.k,
            "cohort plan needs one held-out list, drop count and reference balance per outer fold"
        );
        for (f, held) in self.held_out.iter().enumerate() {
            ensure!(
                held.windows(2).all(|w| w[0] < w[1])
                    && held.iter().all(|&i| i < n && self.outer[i] != Some(f))
                    && held.len()
                        == self.balance_dropped[f]
                            + self.lineage_dropped[f]
                            + self.duplicate_dropped[f],
                "outer fold {f}'s held-out samples are not ascending unscored samples matching its drop counts"
            );
        }
        if self.unavailable_reason.is_none() {
            ensure!(
                self.inner.len() == self.k && self.holdouts.len() == self.k,
                "available cohort plan needs {} inner plans and holdout pairs",
                self.k
            );
            for (f, inner) in self.inner.iter().enumerate() {
                let train = self.train(f);
                ensure!(
                    inner.len() == n
                        && inner
                            .iter()
                            .enumerate()
                            .all(|(i, fold)| fold.is_some() == train.contains(i)),
                    "outer fold {f}'s inner plan does not match its training samples"
                );
            }
        }
        let mut value = serde_json::to_value(self)?;
        if let Value::Object(map) = &mut value {
            map.remove("sha256");
        }
        ensure!(
            sha256_hex(canonical_json(&value)?.as_bytes()) == self.sha256,
            "cohort plan sha256 does not match its contents"
        );
        Ok(())
    }
}

/// How many resistant and susceptible genomes a class-balanced subsample of
/// `have` = [resistant, susceptible] keeps to match `want`'s ratio: as many as
/// possible, the scarcer class (relative to `want`) kept whole.
fn balanced_counts(have: [usize; 2], want: [usize; 2]) -> [usize; 2] {
    let ([r, s], [wr, ws]) = (have, want);
    if wr == 0 || ws == 0 {
        return [if wr == 0 { 0 } else { r }, if ws == 0 { 0 } else { s }];
    }
    if (r as u128) * (ws as u128) <= (s as u128) * (wr as u128) {
        // Resistant genomes are the scarcer: keep all, match the susceptible.
        let keep_s = ((r as f64) * (ws as f64) / (wr as f64)).round() as usize;
        [r, keep_s.min(s)]
    } else {
        let keep_r = ((s as f64) * (wr as f64) / (ws as f64)).round() as usize;
        [keep_r.min(r), s]
    }
}

/// The cohort plan of an arm with samples `ids` and labels `y`, against the
/// available `reference` plan over `reference_ids` with labels
/// `reference_y`, both grouped by the species grouping `species` (which
/// must hold every genome of both); `subsample_seed` replaces the
/// reference's repeat seed as the subsample's `balance_seed`. A plan that
/// cannot be scored, because some outer fold's scored genomes or its
/// training samples lack a class, or an inner or holdout partition does, is
/// returned with its reason, as `fold_plan` does.
pub fn cohort_plan(
    reference: &FoldPlan,
    reference_ids: &[String],
    reference_y: &Bits,
    ids: &[String],
    y: &Bits,
    species: &SpeciesGroups,
    subsample_seed: Option<u64>,
) -> Result<CohortPlan> {
    reference.verify()?;
    ensure!(
        reference.is_available(),
        "the reference plan cannot be scored: {}",
        reference.unavailable_reason.as_deref().unwrap_or("")
    );
    ensure!(
        reference_ids.len() == reference.outer.len() && reference_y.len() == reference.outer.len(),
        "{} reference ids and {} labels for a plan over {} samples",
        reference_ids.len(),
        reference_y.len(),
        reference.outer.len()
    );
    let n = ids.len();
    ensure!(
        n == y.len(),
        "cohort plan dimensions differ: {n} ids, {} labels",
        y.len()
    );
    let groups = &species.groups(ids)?;
    let units = species.unit_keys(reference_ids, ids, reference.scheme)?;
    let fold_of: HashMap<&str, usize> = reference_ids
        .iter()
        .map(String::as_str)
        .zip(reference.outer.iter().copied())
        .collect();
    ensure!(
        fold_of.len() == reference_ids.len(),
        "a reference id appears twice"
    );
    let common: Vec<Option<usize>> = ids
        .iter()
        .map(|id| fold_of.get(id.as_str()).copied())
        .collect();
    let (k, inner_k, drawn) = (reference.k, reference.inner_k, reference.repeat_seed);
    let subsample = subsample_seed.unwrap_or(drawn);
    let reference_test_classes: Vec<[usize; 2]> = (0..k)
        .map(|f| {
            let test = (0..reference.outer.len()).filter(|&j| reference.outer[j] == f);
            let r = test.clone().filter(|&j| reference_y.contains(j)).count();
            [r, test.count() - r]
        })
        .collect();

    let mut outer = vec![None; n];
    let mut held_out = vec![Vec::new(); k];
    let mut balance_dropped = vec![0; k];
    for f in 0..k {
        let class = |resistant: bool| -> Vec<usize> {
            (0..n)
                .filter(|&i| common[i] == Some(f) && y.contains(i) == resistant)
                .collect()
        };
        let (mut r, mut s) = (class(true), class(false));
        let [keep_r, keep_s] = balanced_counts([r.len(), s.len()], reference_test_classes[f]);
        // Keyed on genome id, so a genome's place in the draw doesn't depend on
        // which other genomes the arm holds.
        let rank = balance_seed(subsample, f);
        r.sort_by_key(|&i| (keyed_hash(&ids[i], rank, "balance"), i));
        s.sort_by_key(|&i| (keyed_hash(&ids[i], rank, "balance"), i));
        for &i in r[..keep_r].iter().chain(&s[..keep_s]) {
            outer[i] = Some(f);
        }
        held_out[f].extend(r[keep_r..].iter().chain(&s[keep_s..]));
        balance_dropped[f] = held_out[f].len();
    }
    let mut lineage_dropped = vec![0; k];
    if reference.scheme == Scheme::Lineage {
        for f in 0..k {
            let scored: BTreeSet<u32> = (0..n)
                .filter(|&i| outer[i] == Some(f))
                .map(|i| groups.lineage[i])
                .collect();
            let leaking: Vec<usize> = (0..n)
                .filter(|&i| common[i].is_none() && scored.contains(&groups.lineage[i]))
                .collect();
            lineage_dropped[f] = leaking.len();
            held_out[f].extend(leaking);
        }
    }
    // An exact copy of a scored genome never trains in its fold, on any split
    // (G's fold plans keep a duplicate group in one fold; a cohort arm's own
    // genomes need this).
    let mut duplicate_dropped = vec![0; k];
    for f in 0..k {
        let held: BTreeSet<usize> = held_out[f].iter().copied().collect();
        let copies: BTreeSet<usize> = groups
            .duplicate
            .iter()
            .filter(|group| group.iter().any(|&i| outer[i] == Some(f)))
            .flat_map(|group| group.iter().copied())
            .filter(|&i| outer[i] != Some(f) && !held.contains(&i))
            .collect();
        duplicate_dropped[f] = copies.len();
        held_out[f].extend(copies);
    }
    for held in &mut held_out {
        held.sort_unstable();
    }

    let mut plan = CohortPlan {
        format: COHORT_PLAN_FORMAT.to_string(),
        reference_sha256: reference.sha256.clone(),
        scheme: reference.scheme,
        repeat: reference.repeat,
        k,
        inner_k,
        seed: reference.seed,
        repeat_seed: drawn,
        balance_seed: subsample,
        reference_test_classes,
        outer,
        held_out,
        balance_dropped,
        lineage_dropped,
        duplicate_dropped,
        inner: vec![],
        holdouts: vec![],
        unavailable_reason: None,
        sha256: String::new(),
    };
    let nested = || -> Result<NestedPlan> {
        let (mut inner, mut holdouts) = (Vec::with_capacity(k), Vec::with_capacity(k));
        for f in 0..k {
            let test = Bits::from_indices(n, (0..n).filter(|&i| plan.outer[i] == Some(f)));
            let train = plan.train(f);
            for (side, mask) in [("scored test", &test), ("training", &train)] {
                let pos = mask.and(y).count();
                ensure!(
                    pos > 0 && pos < mask.count(),
                    "outer fold {f}'s {side} samples lack one or both classes"
                );
            }
            let (folds, held) = keyed_tuning(
                &units,
                y,
                &train,
                inner_k,
                fold_seed(drawn, f),
                plan.scheme,
                f,
            )?;
            inner.push(folds);
            holdouts.push(held);
        }
        Ok((inner, holdouts))
    };
    match nested() {
        Ok((inner, holdouts)) => {
            plan.inner = inner;
            plan.holdouts = holdouts;
        }
        Err(error) => plan.unavailable_reason = Some(format!("{error:#}")),
    }
    let mut value = serde_json::to_value(&plan)?;
    if let Value::Object(map) = &mut value {
        map.remove("sha256");
    }
    plan.sha256 = sha256_hex(canonical_json(&value)?.as_bytes());
    Ok(plan)
}

/// The `format` of a keyed fold plan (`keyed_fold_plan`).
pub const KEYED_PLAN_FORMAT: &str = "grm-fold-plan-keyed-v1";

/// A genome-keyed draw: the first eight bytes, big-endian, of the SHA-256 of
/// the ASCII text `grm-keyed:{purpose}:{seed}:{key}`. It depends on the key,
/// the seed and the purpose alone, never on which other genomes are drawn.
pub fn keyed_hash(key: &str, seed: u64, purpose: &str) -> u64 {
    let tag = format!("grm-keyed:{purpose}:{seed}:{key}");
    let digest = ring::digest::digest(&ring::digest::SHA256, tag.as_bytes());
    let mut head = [0u8; 8];
    head.copy_from_slice(&digest.as_ref()[..8]);
    u64::from_be_bytes(head)
}

/// Each sample's draw unit under `scheme`, by key: on the lineage split its
/// lineage, keyed by the lowest genome id of it in `ids`; on the random split
/// its k-mer duplicate group, keyed by the group's lowest genome id in `ids`,
/// or the genome's own id when it has no duplicate. A unit moves whole. The
/// keys depend on `ids`: a cohort that gains or loses a unit's lowest id
/// re-keys the unit, which is why a cohort arm keys its units through the
/// reference (`SpeciesGroups::unit_keys`). For the reference itself, grouped
/// by `SpeciesGroups::groups`, the two give the same keys.
pub fn unit_keys(ids: &[String], groups: &GroupIds, scheme: Scheme) -> Vec<String> {
    let mut key: Vec<String> = ids.to_vec();
    match scheme {
        Scheme::Lineage => {
            let mut lowest: HashMap<u32, &str> = HashMap::new();
            for (i, id) in ids.iter().enumerate() {
                let e = lowest.entry(groups.lineage[i]).or_insert(id.as_str());
                if id.as_str() < *e {
                    *e = id.as_str();
                }
            }
            for (i, k) in key.iter_mut().enumerate() {
                *k = lowest[&groups.lineage[i]].to_string();
            }
        }
        Scheme::Random => {
            for group in &groups.duplicate {
                let low = group
                    .iter()
                    .map(|&i| ids[i].as_str())
                    .min()
                    .unwrap_or_default()
                    .to_string();
                for &i in group {
                    key[i] = low.clone();
                }
            }
        }
    }
    key
}

/// The k = 31 grouping every arm of one input comparison shares (main's
/// rulings on REVIEW-INPUT3 N2 and REVIEW-INPUT4 M1): each genome of the
/// species pool, the union of the comparison's arms' cohorts, with its
/// lineage and its kset digest, as one `grm structure` over that pool gave
/// them, before any `grm arm`. An arm's lineages and duplicate
/// groups are this grouping restricted to its genomes (`groups`), so two
/// arms never disagree on whether two genomes share a lineage or are
/// copies; a genome outside the pool is refused.
#[derive(Clone, Debug)]
pub struct SpeciesGroups {
    ids: Vec<String>,
    lineage: Vec<u32>,
    digest: Vec<String>,
    at: HashMap<String, usize>,
}

impl SpeciesGroups {
    /// Fails unless `ids` are distinct, every genome has a lineage and a
    /// digest, and genomes with one digest share a lineage.
    pub fn new(ids: Vec<String>, lineage: Vec<u32>, digest: Vec<String>) -> Result<Self> {
        ensure!(
            ids.len() == lineage.len() && ids.len() == digest.len(),
            "the species grouping has {} genomes, {} lineages and {} kset digests",
            ids.len(),
            lineage.len(),
            digest.len()
        );
        let mut at = HashMap::with_capacity(ids.len());
        let mut lineage_of: HashMap<&str, u32> = HashMap::new();
        for (i, id) in ids.iter().enumerate() {
            ensure!(
                at.insert(id.clone(), i).is_none(),
                "genome {id} appears twice in the species grouping"
            );
            ensure!(!digest[i].is_empty(), "genome {id} has no kset digest");
            let l = *lineage_of.entry(digest[i].as_str()).or_insert(lineage[i]);
            ensure!(
                l == lineage[i],
                "genome {id} has another genome's kset digest but not its lineage"
            );
        }
        Ok(SpeciesGroups {
            ids,
            lineage,
            digest,
            at,
        })
    }

    fn index(&self, ids: &[String]) -> Result<Vec<usize>> {
        ids.iter()
            .map(|id| {
                self.at
                    .get(id)
                    .copied()
                    .with_context(|| format!("genome {id} is not in the species grouping"))
            })
            .collect()
    }

    /// `ids`' group ids: each genome's lineage in the grouping, and as
    /// duplicate groups the genomes of `ids` that share a kset digest.
    pub fn groups(&self, ids: &[String]) -> Result<GroupIds> {
        let at = self.index(ids)?;
        let mut by_digest: BTreeMap<&str, Vec<usize>> = BTreeMap::new();
        for (i, &p) in at.iter().enumerate() {
            by_digest
                .entry(self.digest[p].as_str())
                .or_default()
                .push(i);
        }
        let mut duplicate: Vec<Vec<usize>> =
            by_digest.into_values().filter(|g| g.len() >= 2).collect();
        duplicate.sort_unstable();
        GroupIds::new(
            at.iter().map(|&p| self.lineage[p]).collect(),
            duplicate,
            None,
        )
    }

    /// Each of `ids`' draw units under `scheme` (its lineage, or its kset
    /// digest's genomes), keyed by the unit's lowest genome id in the
    /// reference cohort `reference_ids`, or, for a unit the reference doesn't
    /// hold, by its lowest genome id in the pool. A key so depends on the
    /// unit, the reference cohort and the pool alone, never on which other
    /// genomes `ids` holds. The pool is the union of the run's arms, so an
    /// arm-only unit's key can change when an arm is added to the run.
    pub fn unit_keys(
        &self,
        reference_ids: &[String],
        ids: &[String],
        scheme: Scheme,
    ) -> Result<Vec<String>> {
        let unit = |p: usize| match scheme {
            Scheme::Lineage => format!("lineage:{}", self.lineage[p]),
            Scheme::Random => format!("digest:{}", self.digest[p]),
        };
        let lowest = |members: &mut dyn Iterator<Item = usize>| {
            let mut low: HashMap<String, &str> = HashMap::new();
            for p in members {
                let e = low.entry(unit(p)).or_insert(self.ids[p].as_str());
                if self.ids[p].as_str() < *e {
                    *e = self.ids[p].as_str();
                }
            }
            low
        };
        let in_reference = lowest(&mut self.index(reference_ids)?.into_iter());
        let in_pool = lowest(&mut (0..self.ids.len()));
        Ok(self
            .index(ids)?
            .into_iter()
            .map(|p| {
                let u = unit(p);
                in_reference.get(&u).unwrap_or(&in_pool[&u]).to_string()
            })
            .collect())
    }

    /// A grouping of the union of `cohorts` (genome ids with their group
    /// ids), for tests: one lineage label space, and genomes in one
    /// duplicate group of any cohort share a digest.
    #[cfg(test)]
    pub(crate) fn of_cohorts(cohorts: &[(&[String], &GroupIds)]) -> SpeciesGroups {
        let mut root: BTreeMap<String, String> = BTreeMap::new();
        let mut lineage: BTreeMap<String, u32> = BTreeMap::new();
        fn find(root: &BTreeMap<String, String>, id: &str) -> String {
            let mut id = id.to_string();
            while root[&id] != id {
                id = root[&id].clone();
            }
            id
        }
        for (ids, groups) in cohorts {
            for (i, id) in ids.iter().enumerate() {
                root.entry(id.clone()).or_insert_with(|| id.clone());
                lineage.insert(id.clone(), groups.lineage()[i]);
            }
            for group in groups.duplicate() {
                for &i in &group[1..] {
                    let (a, b) = (find(&root, &ids[group[0]]), find(&root, &ids[i]));
                    if a != b {
                        let (low, high) = if a < b { (a, b) } else { (b, a) };
                        root.insert(high, low);
                    }
                }
            }
        }
        let ids: Vec<String> = root.keys().cloned().collect();
        let digest = ids.iter().map(|id| find(&root, id)).collect();
        let lineage = ids.iter().map(|id| lineage[id]).collect();
        SpeciesGroups::new(ids, lineage, digest).unwrap()
    }
}

/// One outer fold's keyed inner folds and tuning holdouts over its `train`
/// samples: a unit's inner fold is `keyed_hash(unit, seed, "inner") %
/// inner_k`; it is in boosting's early-stopping holdout when
/// `keyed_hash(unit, seed, "gbt") % 7 == 0` (round(1 / 0.15)) and in the
/// lasso's when `keyed_hash(unit, seed, "lasso") % 5 == 0` (round(1 / 0.2)).
/// Every inner fold's training and validation samples, and both sides of
/// both holdouts, must hold both classes, or the outer fold is refused.
fn keyed_tuning(
    units: &[String],
    y: &Bits,
    train: &Bits,
    inner_k: usize,
    seed: u64,
    scheme: Scheme,
    fi: usize,
) -> Result<(Vec<Option<usize>>, TuningHoldouts)> {
    let n = y.len();
    let folds: Vec<usize> = (0..n)
        .map(|i| {
            if train.contains(i) {
                (keyed_hash(&units[i], seed, "inner") % inner_k as u64) as usize
            } else {
                usize::MAX
            }
        })
        .collect();
    validate_folds(y, train, &folds, inner_k)
        .with_context(|| format!("invalid {scheme} inner evaluation in outer fold {fi}"))?;
    let holdout = |purpose: &str, frac: f64| -> Result<Vec<usize>> {
        let m = ((1.0 / frac).round() as u64).max(2);
        let val = Bits::from_indices(
            n,
            train
                .iter_ones()
                .filter(|&i| keyed_hash(&units[i], seed, purpose).is_multiple_of(m)),
        );
        let fit = train.and_not(&val);
        for mask in [&fit, &val] {
            let pos = mask.and(y).count();
            ensure!(
                pos > 0 && pos < mask.count(),
                "invalid {scheme} tuning holdout in outer fold {fi}: the {purpose} holdout lacks one or both classes on one side"
            );
        }
        Ok(val.iter_ones().collect())
    };
    let held = TuningHoldouts {
        gbt: holdout("gbt", GBT_HOLDOUT)?,
        lasso: holdout("lasso", LASSO_HOLDOUT)?,
    };
    let inner = (0..n)
        .map(|i| train.contains(i).then_some(folds[i]))
        .collect();
    Ok((inner, held))
}

/// The input comparison's fold plan: G's `FoldPlan` shape, with every
/// partition keyed on genome identity instead of drawn over sample
/// indices (main's ruling on REVIEW-INPUT2 N1). A unit (`unit_keys`) goes to
/// outer fold `keyed_hash(unit, repeat_seed, "outer") % k`, and each outer
/// fold's inner folds and holdouts follow `keyed_tuning` seeded
/// `fold_seed(repeat_seed, f)`. A genome's folds depend on its unit's key
/// and on no genome's position. With `groups` from `SpeciesGroups::groups`,
/// a unit's key is its lowest genome id in `ids`, so another cohort keeps a
/// shared genome's folds as long as it holds that genome's unit's lowest id
/// (a cohort arm keys through the reference instead: `cohort_plan`). Folds
/// are balanced only approximately:
/// the plan is refused (`unavailable_reason`) unless the lineage split has
/// at least 2k lineages, every outer fold's test and training samples hold
/// both classes, and every inner fold and holdout passes `keyed_tuning`.
#[allow(clippy::too_many_arguments)]
pub fn keyed_fold_plan(
    ids: &[String],
    y: &Bits,
    groups: &GroupIds,
    scheme: Scheme,
    k: usize,
    inner_k: usize,
    seed: u64,
    repeat: u32,
) -> Result<FoldPlan> {
    ensure!(
        ids.len() == y.len() && ids.len() == groups.len(),
        "keyed fold plan dimensions differ: {} ids, {} labels, {} samples with group ids",
        ids.len(),
        y.len(),
        groups.len()
    );
    ensure!(
        k >= 2 && inner_k >= 2,
        "cross-validation needs at least two folds"
    );
    let drawn = repeat_seed(seed, repeat);
    let units = unit_keys(ids, groups, scheme);
    let outer: Vec<usize> = units
        .iter()
        .map(|u| (keyed_hash(u, drawn, "outer") % k as u64) as usize)
        .collect();
    let nested = || -> Result<NestedPlan> {
        if scheme == Scheme::Lineage {
            let n_groups = groups.lineage.iter().collect::<BTreeSet<_>>().len();
            ensure!(
                n_groups >= k * 2,
                "fewer than twice as many lineages as requested folds"
            );
        }
        validate_folds(y, &Bits::ones(y.len()), &outer, k)
            .with_context(|| format!("invalid {scheme} outer evaluation"))?;
        let (mut inner, mut holdouts) = (Vec::with_capacity(k), Vec::with_capacity(k));
        for (fi, (train, _)) in fold_masks(&outer, k).into_iter().enumerate() {
            let (folds, held) =
                keyed_tuning(&units, y, &train, inner_k, fold_seed(drawn, fi), scheme, fi)?;
            inner.push(folds);
            holdouts.push(held);
        }
        Ok((inner, holdouts))
    };
    let (inner, holdouts, unavailable_reason) = match nested() {
        Ok((inner, holdouts)) => (inner, holdouts, None),
        Err(error) => (vec![], vec![], Some(format!("{error:#}"))),
    };
    let mut plan = FoldPlan {
        format: KEYED_PLAN_FORMAT.to_string(),
        scheme,
        repeat,
        k,
        inner_k,
        seed,
        repeat_seed: drawn,
        outer,
        inner,
        holdouts,
        unavailable_reason,
        sha256: String::new(),
    };
    plan.sha256 = plan.compute_sha256()?;
    Ok(plan)
}

type NestedPlan = (Vec<Vec<Option<usize>>>, Vec<TuningHoldouts>);

/// The inner folds and tuning holdouts of every outer fold, after the checks
/// `learn::eval::check_split_plan` makes before any fit.
fn nested_plan(
    y: &Bits,
    groups: &GroupIds,
    scheme: Scheme,
    outer: &[usize],
    k: usize,
    inner_k: usize,
    drawn: u64,
) -> Result<NestedPlan> {
    if scheme == Scheme::Lineage {
        let n_groups = groups.lineage.iter().collect::<BTreeSet<_>>().len();
        ensure!(
            n_groups >= k * 2,
            "fewer than twice as many lineages as requested folds"
        );
    }
    validate_folds(y, &Bits::ones(y.len()), outer, k)
        .with_context(|| format!("invalid {scheme} outer evaluation"))?;
    let tuning = tuning_groups(scheme, &groups.lineage, &groups.duplicate);
    let mut inner = Vec::with_capacity(k);
    let mut holdouts = Vec::with_capacity(k);
    for (fi, (train, _)) in fold_masks(outer, k).into_iter().enumerate() {
        let (folds, held) = fold_tuning_plan(
            y,
            &train,
            tuning.as_deref(),
            inner_k,
            fold_seed(drawn, fi),
            scheme,
            fi,
        )?;
        inner.push(folds);
        holdouts.push(held);
    }
    Ok((inner, holdouts))
}

/// One outer fold's inner folds and tuning holdouts. Each is drawn from the
/// samples in `train` and their labels only. A failure reads as
/// `check_split_plan` words it for outer fold `fi` of `scheme`.
fn fold_tuning_plan(
    y: &Bits,
    train: &Bits,
    tuning: Option<&[u32]>,
    inner_k: usize,
    seed: u64,
    scheme: Scheme,
    fi: usize,
) -> Result<(Vec<Option<usize>>, TuningHoldouts)> {
    let folds = training_folds(y, train, tuning, inner_k, seed);
    validate_folds(y, train, &folds, inner_k)
        .with_context(|| format!("invalid {scheme} inner evaluation in outer fold {fi}"))?;
    let holdout = |frac| {
        let unit = crate::learn::eval::tuning_unit(scheme.as_str());
        tuning_holdout(train, y, tuning, unit, frac, seed)
            .with_context(|| format!("invalid {scheme} tuning holdout in outer fold {fi}"))
    };
    let (_, gbt) = holdout(GBT_HOLDOUT)?;
    let (_, lasso) = holdout(LASSO_HOLDOUT)?;
    let inner = folds
        .iter()
        .enumerate()
        .map(|(i, &f)| train.contains(i).then_some(f))
        .collect();
    Ok((
        inner,
        TuningHoldouts {
            gbt: gbt.iter_ones().collect(),
            lasso: lasso.iter_ones().collect(),
        },
    ))
}

/// The canonical byte form of the shared fold hash, identical to Python's
/// `json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)`:
/// - object keys sorted by code point, and no whitespace between tokens;
/// - integers in plain decimal, of any size up to u64; floats are refused,
///   so float formatting never enters the hash;
/// - strings escape `"` as `\"` and `\` as `\\`, write `\b \f \n \r \t` for
///   those controls, and write every other character outside U+0020..=U+007E
///   as `\uXXXX` in lowercase hex (a UTF-16 surrogate pair above U+FFFF), so
///   the output is ASCII.
///
/// The fold hash is the lowercase hex SHA-256 of these ASCII bytes.
pub fn canonical_json(value: &Value) -> Result<String> {
    let mut out = String::new();
    write_canonical(value, &mut out)?;
    Ok(out)
}

fn write_canonical(value: &Value, out: &mut String) -> Result<()> {
    match value {
        Value::Null => out.push_str("null"),
        Value::Bool(b) => out.push_str(if *b { "true" } else { "false" }),
        Value::Number(n) => {
            if let Some(u) = n.as_u64() {
                write!(out, "{u}")?;
            } else if let Some(i) = n.as_i64() {
                write!(out, "{i}")?;
            } else {
                bail!("canonical JSON holds integers only, not {n}");
            }
        }
        Value::String(s) => write_json_string(s, out)?,
        Value::Array(items) => {
            out.push('[');
            for (j, item) in items.iter().enumerate() {
                if j > 0 {
                    out.push(',');
                }
                write_canonical(item, out)?;
            }
            out.push(']');
        }
        Value::Object(map) => {
            let mut keys: Vec<&String> = map.keys().collect();
            keys.sort();
            out.push('{');
            for (j, key) in keys.into_iter().enumerate() {
                if j > 0 {
                    out.push(',');
                }
                write_json_string(key, out)?;
                out.push(':');
                write_canonical(&map[key], out)?;
            }
            out.push('}');
        }
    }
    Ok(())
}

fn write_json_string(s: &str, out: &mut String) -> Result<()> {
    out.push('"');
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\u{8}' => out.push_str("\\b"),
            '\u{c}' => out.push_str("\\f"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            ' '..='~' => out.push(c),
            _ => {
                let mut units = [0u16; 2];
                for unit in c.encode_utf16(&mut units) {
                    write!(out, "\\u{unit:04x}")?;
                }
            }
        }
    }
    out.push('"');
    Ok(())
}

fn sha256_hex(bytes: &[u8]) -> String {
    ring::digest::digest(&ring::digest::SHA256, bytes)
        .as_ref()
        .iter()
        .map(|b| format!("{b:02x}"))
        .collect()
}

// ---------------------------------------------------------------- compute

/// One learner's tuning cost in fit units.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct LearnerCost {
    pub learner: String,
    /// Model fits one outer fold spends choosing hyperparameters, before its
    /// final fit on the outer fold's training samples.
    pub inner_candidates: u64,
}

impl LearnerCost {
    /// The pipeline's tuning fits per outer fold at c90ebaf, with `inner_k`
    /// inner folds. A fit unit is one model fit, whatever its size:
    /// - scm: 2 rule kinds x 7 penalties (`scm::P_GRID`) x inner_k, one
    ///   longest-rule model per inner fold whose prefixes are all scored;
    /// - cart: inner_k, one deepest tree per inner fold whose depths are all
    ///   scored;
    /// - rf: 0, no tuning;
    /// - gbt: 0, its one early-stopped fit on the 0.15 holdout is the final
    ///   model;
    /// - lasso: 1, the penalty path on the 0.2 holdout's fitting samples,
    ///   before the final refit path on all training samples.
    ///
    /// Fails for fewer than two inner folds, as the pipeline does, and when
    /// the count overflows u64.
    pub fn pipeline(learner: &str, inner_k: usize) -> Result<Self> {
        ensure!(inner_k >= 2, "cross-validation needs at least two folds");
        let inner_k = u64::try_from(inner_k).context("inner folds overflow u64")?;
        let scm_fits = 2 * scm::P_GRID.len() as u64;
        let inner_candidates = match learner {
            "scm" => scm_fits.checked_mul(inner_k),
            "cart" => Some(inner_k),
            "rf" | "gbt" => Some(0),
            "lasso" => Some(1),
            _ => bail!("no pipeline tuning cost for learner {learner:?}"),
        }
        .context("fit units overflow u64")?;
        Ok(LearnerCost {
            learner: learner.to_string(),
            inner_candidates,
        })
    }
}

/// One tier of an ablation protocol: what it runs and the most it may spend.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Tier {
    pub name: String,
    /// The most fit units this tier may spend; `plan_compute` refuses more.
    pub budget: u64,
    /// Ablation arms, by name; Task H defines their blocks.
    pub arms: Vec<String>,
    pub learners: Vec<LearnerCost>,
    pub schemes: Vec<Scheme>,
    /// Repeats 0..repeats are run.
    pub repeats: u32,
    /// Outer folds per plan.
    pub k: usize,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Protocol {
    pub tiers: Vec<Tier>,
}

/// One dataset and the tiers it is run in.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct InventoryDataset {
    pub name: String,
    pub tiers: Vec<String>,
    /// Schemes whose fold plan is unavailable for this dataset; they are
    /// not run, so they cost nothing.
    #[serde(default)]
    pub unavailable: Vec<Scheme>,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Inventory {
    pub datasets: Vec<InventoryDataset>,
}

/// The fit units of one (tier, arm, learner, split, repeat, dataset).
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct ComputeRow {
    pub tier: String,
    pub arm: String,
    pub learner: String,
    pub scheme: Scheme,
    pub repeat: u32,
    pub dataset: String,
    pub units: u64,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct TierCompute {
    pub tier: String,
    pub units: u64,
    pub budget: u64,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct ComputePlan {
    /// Rows in protocol order: tier, arm, learner, scheme, repeat, then
    /// dataset in inventory order.
    pub rows: Vec<ComputeRow>,
    /// One total per tier, in protocol order.
    pub tiers: Vec<TierCompute>,
    pub total_units: u64,
}

/// Counts the fit units a protocol would run on an inventory, and refuses
/// when any tier's units exceed its budget. One row's units are
///
/// `k x (inner_candidates + 1)`
///
/// that is outer folds x (tuning fits + the final fit) for that tier's `k`
/// and that learner's `inner_candidates`. A dataset runs in the tiers it
/// names, for every scheme it does not list as unavailable.
///
/// Ablation evidence comes from outer-fold predictions alone, so ablation
/// runs fit no full-data final model and no MIC regression. The formula has
/// no term for them because no such fit is run.
pub fn plan_compute(protocol: &Protocol, inventory: &Inventory) -> Result<ComputePlan> {
    unique_names("tier", protocol.tiers.iter().map(|t| t.name.as_str()))?;
    unique_names(
        "dataset",
        inventory.datasets.iter().map(|d| d.name.as_str()),
    )?;
    let known: BTreeSet<&str> = protocol.tiers.iter().map(|t| t.name.as_str()).collect();
    for dataset in &inventory.datasets {
        for tier in &dataset.tiers {
            ensure!(
                known.contains(tier.as_str()),
                "dataset {} names unknown tier {tier}",
                dataset.name
            );
        }
    }
    let mut rows = Vec::new();
    let mut tiers = Vec::new();
    for tier in &protocol.tiers {
        ensure!(tier.k >= 2, "tier {} needs at least two folds", tier.name);
        ensure!(tier.repeats >= 1, "tier {} runs no repeat", tier.name);
        unique_names("arm", tier.arms.iter().map(String::as_str))?;
        unique_names("learner", tier.learners.iter().map(|l| l.learner.as_str()))?;
        ensure!(
            !tier.schemes.is_empty()
                && tier.schemes.iter().collect::<BTreeSet<_>>().len() == tier.schemes.len(),
            "tier {} needs one or more distinct schemes",
            tier.name
        );
        let datasets: Vec<&InventoryDataset> = inventory
            .datasets
            .iter()
            .filter(|d| d.tiers.contains(&tier.name))
            .collect();
        let mut units = 0u64;
        for arm in &tier.arms {
            for learner in &tier.learners {
                let per_plan = learner
                    .inner_candidates
                    .checked_add(1)
                    .and_then(|fits| fits.checked_mul(tier.k as u64))
                    .context("fit units overflow u64")?;
                for &scheme in &tier.schemes {
                    for repeat in 0..tier.repeats {
                        for dataset in &datasets {
                            if dataset.unavailable.contains(&scheme) {
                                continue;
                            }
                            units = units
                                .checked_add(per_plan)
                                .context("fit units overflow u64")?;
                            rows.push(ComputeRow {
                                tier: tier.name.clone(),
                                arm: arm.clone(),
                                learner: learner.learner.clone(),
                                scheme,
                                repeat,
                                dataset: dataset.name.clone(),
                                units: per_plan,
                            });
                        }
                    }
                }
            }
        }
        tiers.push(TierCompute {
            tier: tier.name.clone(),
            units,
            budget: tier.budget,
        });
    }
    let over: Vec<String> = tiers
        .iter()
        .filter(|t| t.units > t.budget)
        .map(|t| {
            format!(
                "tier {} needs {} fit units, over its budget of {}",
                t.tier, t.units, t.budget
            )
        })
        .collect();
    ensure!(over.is_empty(), "compute plan refused: {}", over.join("; "));
    let total_units = tiers
        .iter()
        .try_fold(0u64, |sum, t| sum.checked_add(t.units))
        .context("fit units overflow u64")?;
    Ok(ComputePlan {
        rows,
        tiers,
        total_units,
    })
}

fn unique_names<'a>(what: &str, names: impl Iterator<Item = &'a str>) -> Result<()> {
    let mut seen = BTreeSet::new();
    let mut any = false;
    for name in names {
        any = true;
        ensure!(!name.is_empty(), "a {what} has an empty name");
        ensure!(seen.insert(name), "{what} {name} is named twice");
    }
    ensure!(any, "no {what} is given");
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::Scheme::{Lineage, Random};
    use super::*;

    /// The re-key case of REVIEW-INPUT3 N2: 640 genomes `g0000`.. in 80
    /// lineages of 8 (a third resistant in each), `g0001` and `g0002` one
    /// duplicate group, `g0010` and `a-dup` another; the pool adds `a00`..`a09`
    /// (ids below every `g`) to lineages 10 to 19.
    fn rekey_pool() -> (Vec<String>, Bits, SpeciesGroups) {
        let n = 640;
        let mut ids: Vec<String> = (0..n).map(|i| format!("g{i:04}")).collect();
        let mut lineage: Vec<u32> = (0..n as u32).map(|i| i / 8).collect();
        let mut digest: Vec<String> = ids.clone();
        digest[2] = digest[1].clone();
        for l in 0..10u32 {
            ids.push(format!("a{l:02}"));
            lineage.push(10 + l);
            digest.push(format!("a{l:02}"));
        }
        ids.push("a-dup".into());
        lineage.push(lineage[10]);
        digest.push(digest[10].clone());
        let reference: Vec<String> = ids[..n].to_vec();
        let y = Bits::from_indices(n, (0..n).filter(|i| i % 3 == 0));
        (
            reference,
            y,
            SpeciesGroups::new(ids, lineage, digest).unwrap(),
        )
    }

    /// A cohort arm that drops the lowest-id genome of lineages 0 to 4 and of
    /// the duplicate group {g0001, g0002}, and adds genomes below the lowest
    /// id of lineages 10 to 19 and of g0010's group, keeps every shared
    /// genome's inner fold and holdout memberships, on both splits and three
    /// seeds (main's ruling on REVIEW-INPUT3 N2). Keyed on the arm's own
    /// lowest ids (`unit_keys`), the same arm re-keys those units.
    #[test]
    fn a_cohort_arm_keeps_shared_folds_when_it_drops_or_adds_a_unit_minimum() {
        let (ref_ids, ref_y, species) = rekey_pool();
        let n = ref_ids.len();
        let dropped: Vec<usize> = (0..5).map(|l| l * 8).chain([1]).collect();
        let kept: Vec<usize> = (0..n).filter(|i| !dropped.contains(i)).collect();
        let added: Vec<String> = (0..10)
            .map(|l| format!("a{l:02}"))
            .chain(["a-dup".to_string()])
            .collect();
        let mut arm_ids: Vec<String> = kept.iter().map(|&i| ref_ids[i].clone()).collect();
        arm_ids.extend(added.iter().cloned());
        let order = {
            let mut o: Vec<usize> = (0..arm_ids.len()).collect();
            o.sort_by(|&a, &b| arm_ids[a].cmp(&arm_ids[b]));
            o
        };
        let arm_ids: Vec<String> = order.iter().map(|&a| arm_ids[a].clone()).collect();
        let resistant = |id: &str| match ref_ids.iter().position(|r| r == id) {
            Some(i) => ref_y.contains(i),
            None => id.ends_with('0') || id.ends_with('5'),
        };
        let arm_y = Bits::from_indices(
            arm_ids.len(),
            (0..arm_ids.len()).filter(|&a| resistant(&arm_ids[a])),
        );
        let ref_groups = species.groups(&ref_ids).unwrap();
        let arm_groups = species.groups(&arm_ids).unwrap();
        let at: HashMap<&str, usize> = ref_ids
            .iter()
            .enumerate()
            .map(|(i, s)| (s.as_str(), i))
            .collect();
        for scheme in [Lineage, Random] {
            // The case does re-key under the arm's own lowest ids.
            let own = unit_keys(&arm_ids, &arm_groups, scheme);
            let through = species.unit_keys(&ref_ids, &arm_ids, scheme).unwrap();
            let ref_keys = unit_keys(&ref_ids, &ref_groups, scheme);
            let rekeyed = (0..arm_ids.len())
                .filter(|&a| at.contains_key(arm_ids[a].as_str()))
                .filter(|&a| own[a] != ref_keys[at[arm_ids[a].as_str()]])
                .count();
            assert!(rekeyed >= 2, "{scheme}: {rekeyed} shared genomes re-keyed");
            for (a, id) in arm_ids.iter().enumerate() {
                if let Some(&i) = at.get(id.as_str()) {
                    assert_eq!(through[a], ref_keys[i], "{scheme} {id}");
                }
            }
            for seed in [1, 2, 3] {
                let r =
                    keyed_fold_plan(&ref_ids, &ref_y, &ref_groups, scheme, 5, 3, seed, 0).unwrap();
                assert!(
                    r.is_available(),
                    "{scheme} seed {seed}: {:?}",
                    r.unavailable_reason
                );
                let c =
                    cohort_plan(&r, &ref_ids, &ref_y, &arm_ids, &arm_y, &species, None).unwrap();
                assert!(
                    c.unavailable_reason.is_none(),
                    "{scheme} seed {seed}: {:?}",
                    c.unavailable_reason
                );
                let mut compared = 0;
                for f in 0..5 {
                    let train = c.train(f);
                    for (a, id) in arm_ids.iter().enumerate() {
                        let Some(&i) = at.get(id.as_str()) else {
                            continue;
                        };
                        if let Some(o) = c.outer[a] {
                            assert_eq!(o, r.outer[i], "{scheme} seed {seed} {id}");
                        }
                        if r.outer[i] == f {
                            continue;
                        }
                        assert!(train.contains(a), "{scheme} seed {seed} {id} trains in {f}");
                        assert_eq!(
                            c.inner[f][a], r.inner[f][i],
                            "{scheme} seed {seed} inner {id} in {f}"
                        );
                        for (ca, ra) in [
                            (&c.holdouts[f].gbt, &r.holdouts[f].gbt),
                            (&c.holdouts[f].lasso, &r.holdouts[f].lasso),
                        ] {
                            assert_eq!(
                                ca.binary_search(&a).is_ok(),
                                ra.binary_search(&i).is_ok(),
                                "{scheme} seed {seed} holdout {id} in {f}"
                            );
                        }
                        compared += 1;
                    }
                }
                assert_eq!(compared, (n - dropped.len()) * 4, "{scheme} seed {seed}");
            }
        }
    }

    /// A genome missing from the species grouping is refused, and so is a
    /// reference genome: an arm can't be grouped or keyed off the pool.
    #[test]
    fn the_species_grouping_refuses_a_genome_outside_the_pool() {
        let (ref_ids, _, species) = rekey_pool();
        let mut ids = ref_ids[..5].to_vec();
        ids.push("z-outside".into());
        let why = format!("{:#}", species.groups(&ids).unwrap_err());
        assert!(
            why.contains("z-outside is not in the species grouping"),
            "{why}"
        );
        assert!(species.unit_keys(&ids, &ref_ids, Lineage).is_err());
        let dup = SpeciesGroups::new(
            vec!["a".into(), "b".into()],
            vec![0, 1],
            vec!["d".into(), "d".into()],
        );
        assert!(dup.is_err(), "copies in two lineages");
    }

    /// The keyed plan's refusals, one per branch (REVIEW-INPUT3 L2).
    #[test]
    fn a_keyed_lineage_plan_with_fewer_than_twice_k_lineages_is_refused() {
        let n = 180;
        let ids: Vec<String> = (0..n).map(|i| format!("g{i:03}")).collect();
        let y = Bits::from_indices(n, (0..n).filter(|i| i % 2 == 0));
        let groups = GroupIds::new((0..n as u32).map(|i| i % 9).collect(), vec![], None).unwrap();
        let p = keyed_fold_plan(&ids, &y, &groups, Lineage, 5, 3, 42, 0).unwrap();
        let why = p
            .unavailable_reason
            .clone()
            .expect("9 lineages for 5 folds");
        assert!(
            why.contains("fewer than twice as many lineages as requested folds"),
            "{why}"
        );
        p.verify().unwrap();
    }

    #[test]
    fn keyed_tuning_refuses_an_inner_fold_or_holdout_side_without_a_class() {
        let n = 300;
        let units: Vec<String> = (0..n).map(|i| format!("g{i:03}")).collect();
        let train = Bits::ones(n);
        let seed = 7;
        let hash = |i: usize, purpose: &str| keyed_hash(&units[i], seed, purpose);
        let tune = |y: &Bits| keyed_tuning(&units, y, &train, 3, seed, Random, 0);
        // Every genome of inner fold 0 susceptible.
        let y = Bits::from_indices(
            n,
            (0..n).filter(|&i| hash(i, "inner") % 3 != 0 && i % 2 == 0),
        );
        let why = format!("{:#}", tune(&y).unwrap_err());
        assert!(
            why.contains("invalid random inner evaluation in outer fold 0"),
            "{why}"
        );
        // Every genome of the boosting holdout susceptible, the inner folds fine.
        let y = Bits::from_indices(
            n,
            (0..n).filter(|&i| !hash(i, "gbt").is_multiple_of(7) && i % 2 == 0),
        );
        let why = format!("{:#}", tune(&y).unwrap_err());
        assert!(
            why.contains("the gbt holdout lacks one or both classes"),
            "{why}"
        );
        // The same for the lasso's, with boosting's holdout holding both.
        let y = Bits::from_indices(
            n,
            (0..n).filter(|&i| !hash(i, "lasso").is_multiple_of(5) && i % 2 == 0),
        );
        let why = format!("{:#}", tune(&y).unwrap_err());
        assert!(
            why.contains("the lasso holdout lacks one or both classes"),
            "{why}"
        );
        let y = Bits::from_indices(n, (0..n).filter(|i| i % 2 == 0));
        assert!(tune(&y).is_ok());
    }

    /// The cohort plan's class-balanced subsample is keyed on genome id
    /// (REVIEW-INPUT3 L3): an arm with 30 more genomes of its own, all sorting
    /// before the shared ones, scores the same shared genomes in every fold.
    #[test]
    fn the_balanced_subsample_keeps_its_shared_genomes_when_the_arm_grows() {
        let n = 400;
        let ref_ids: Vec<String> = (0..n).map(|i| format!("g{i:03}")).collect();
        let ref_y = Bits::from_indices(n, (0..n).filter(|i| i % 4 == 0));
        let extra: Vec<String> = (0..30).map(|j| format!("a{j:02}")).collect();
        let mut pool_ids = extra.clone();
        pool_ids.extend(ref_ids.iter().cloned());
        let species = SpeciesGroups::new(
            pool_ids.clone(),
            (0..pool_ids.len() as u32).collect(),
            pool_ids.clone(),
        )
        .unwrap();
        let r = keyed_fold_plan(
            &ref_ids,
            &ref_y,
            &species.groups(&ref_ids).unwrap(),
            Random,
            5,
            3,
            42,
            0,
        )
        .unwrap();
        // Every resistant genome, every third susceptible: the resistant ones
        // are cut to the reference's ratio.
        let small: Vec<String> = (0..n)
            .filter(|&i| ref_y.contains(i) || i % 3 == 0)
            .map(|i| ref_ids[i].clone())
            .collect();
        let mut big = extra.clone();
        big.extend(small.iter().cloned());
        let y_of = |ids: &[String]| {
            Bits::from_indices(
                ids.len(),
                (0..ids.len()).filter(|&a| {
                    ids[a].starts_with('a')
                        || ref_y.contains(ref_ids.binary_search(&ids[a]).unwrap())
                }),
            )
        };
        let a = cohort_plan(&r, &ref_ids, &ref_y, &small, &y_of(&small), &species, None).unwrap();
        let b = cohort_plan(&r, &ref_ids, &ref_y, &big, &y_of(&big), &species, None).unwrap();
        assert!(
            a.balance_dropped.iter().all(|&d| d > 0),
            "{:?}",
            a.balance_dropped
        );
        assert_eq!(a.balance_dropped, b.balance_dropped);
        for (j, id) in small.iter().enumerate() {
            assert_eq!(a.outer[j], b.outer[j + extra.len()], "{id}");
        }
    }
    use crate::bits::BitMatrix;
    use crate::learn::eval::{EvalConfig, evaluate};

    /// 160 samples in 40 lineages of four. Every lineage holds both classes,
    /// and every third lineage from the second holds a duplicate pair, some
    /// of them of mixed class.
    fn fixture() -> (Bits, GroupIds) {
        let n = 160;
        let y = Bits::from_indices(
            n,
            (0..n).filter(|&i| {
                let (j, m) = (i / 4, i % 4);
                m == 0 || (m == 1 && j % 2 == 0) || (m == 2 && j % 3 == 0)
            }),
        );
        let lineage = (0..n as u32).map(|i| i / 4).collect();
        let duplicate = (0..n / 4)
            .filter(|j| j % 3 == 1)
            .map(|j| vec![4 * j + 1, 4 * j + 2])
            .collect();
        (y, GroupIds::new(lineage, duplicate, None).unwrap())
    }

    fn flip(y: &mut Bits, i: usize) {
        if y.contains(i) {
            y.remove(i);
        } else {
            y.insert(i);
        }
    }

    /// The first group that some partition of `plan` splits, if any.
    fn crossing(groups: &[Vec<usize>], plan: &FoldPlan) -> Option<String> {
        let split = |group: &[usize], at: &dyn Fn(usize) -> Option<usize>| {
            group.iter().any(|&i| at(i) != at(group[0]))
        };
        for group in groups.iter().map(Vec::as_slice) {
            if split(group, &|i| Some(plan.outer[i])) {
                return Some(format!("{group:?} crosses outer folds"));
            }
            for (fi, inner) in plan.inner.iter().enumerate() {
                if split(group, &|i| inner[i]) {
                    return Some(format!("{group:?} crosses outer fold {fi}'s inner folds"));
                }
                let held = &plan.holdouts[fi];
                for (name, val) in [("gbt", &held.gbt), ("lasso", &held.lasso)] {
                    if split(group, &|i| Some(val.binary_search(&i).is_ok() as usize)) {
                        return Some(format!(
                            "{group:?} crosses outer fold {fi}'s {name} holdout"
                        ));
                    }
                }
            }
        }
        None
    }

    #[test]
    fn repeat_zero_reproduces_the_pipeline_outer_folds() {
        let cfg = EvalConfig::default();
        let (y, g) = fixture();
        let plan =
            |scheme| fold_plan(&y, &g, scheme, cfg.folds, cfg.inner_folds, cfg.seed, 0).unwrap();
        let random = plan(Random);
        let stratified = stratified_folds(&y, cfg.folds, cfg.seed);
        assert_eq!(
            random.outer,
            isolate_groups(&stratified, g.duplicate(), cfg.folds)
        );
        // The fixture's duplicates move samples, so the isolation is tested.
        assert_ne!(random.outer, stratified);
        let lineage = plan(Lineage);
        assert_eq!(
            lineage.outer,
            grouped_folds(&y, g.lineage(), cfg.folds, cfg.seed)
        );
        for p in [&random, &lineage] {
            assert_eq!(p.unavailable_reason, None);
            assert_eq!((p.k, p.inner_k, p.seed, p.repeat_seed), (5, 3, 42, 42));
            assert_eq!((p.inner.len(), p.holdouts.len()), (5, 5));
        }
    }

    #[test]
    fn repeat_zero_inner_plans_follow_the_pipeline_where_its_policy_is_unchanged() {
        let cfg = EvalConfig::default();
        let (y, g) = fixture();
        let no_dups = GroupIds::new(g.lineage().to_vec(), vec![], None).unwrap();
        let cases: [(Scheme, &GroupIds, Option<&[u32]>); 2] =
            [(Lineage, &g, Some(g.lineage())), (Random, &no_dups, None)];
        for (scheme, groups, tuning) in cases {
            let plan =
                fold_plan(&y, groups, scheme, cfg.folds, cfg.inner_folds, cfg.seed, 0).unwrap();
            assert!(plan.is_available());
            for (fi, (train, _)) in fold_masks(&plan.outer, cfg.folds).into_iter().enumerate() {
                // learn::eval::fold_seed
                let seed = cfg.seed + fi as u64 * 101;
                let inner = training_folds(&y, &train, tuning, cfg.inner_folds, seed);
                let expected: Vec<Option<usize>> = (0..y.len())
                    .map(|i| train.contains(i).then_some(inner[i]))
                    .collect();
                assert_eq!(plan.inner[fi], expected, "{scheme} outer fold {fi}");
                let held = &plan.holdouts[fi];
                for (frac, got) in [(0.15, &held.gbt), (0.2, &held.lasso)] {
                    let unit = crate::learn::eval::tuning_unit(scheme.as_str());
                    let (_, val) = tuning_holdout(&train, &y, tuning, unit, frac, seed).unwrap();
                    assert_eq!(*got, val.iter_ones().collect::<Vec<_>>());
                }
            }
        }
    }

    #[test]
    fn tuning_groups_give_composite_ids_only_to_random_splits_with_duplicates() {
        let lineage = [7, 7, 7, 8, 8, 8];
        assert_eq!(tuning_groups(Random, &lineage, &[]), None);
        assert_eq!(
            tuning_groups(Random, &lineage, &[vec![1, 2], vec![3, 5]]),
            Some(vec![0, 1, 1, 3, 4, 3])
        );
        assert_eq!(
            tuning_groups(Lineage, &lineage, &[vec![1, 2]]),
            Some(lineage.to_vec())
        );
    }

    #[test]
    fn no_duplicate_or_lineage_group_crosses_any_partition() {
        let (y, g) = fixture();
        let lineages: Vec<Vec<usize>> = (0..40).map(|j| (4 * j..4 * j + 4).collect()).collect();
        for repeat in 0..3 {
            let random = fold_plan(&y, &g, Random, 5, 3, 42, repeat).unwrap();
            assert!(random.is_available(), "{:?}", random.unavailable_reason);
            assert_eq!(crossing(g.duplicate(), &random), None);
            let lineage = fold_plan(&y, &g, Lineage, 5, 3, 42, repeat).unwrap();
            assert!(lineage.is_available(), "{:?}", lineage.unavailable_reason);
            assert_eq!(crossing(g.duplicate(), &lineage), None);
            assert_eq!(crossing(&lineages, &lineage), None);
        }
        // Control: c90ebaf's stratified inner folds split some duplicate
        // pair, and `crossing` sees it.
        let mut old = fold_plan(&y, &g, Random, 5, 3, 42, 0).unwrap();
        for (fi, (train, _)) in fold_masks(&old.outer, 5).into_iter().enumerate() {
            let inner = training_folds(&y, &train, None, 3, fold_seed(42, fi));
            old.inner[fi] = (0..y.len())
                .map(|i| train.contains(i).then_some(inner[i]))
                .collect();
        }
        assert!(crossing(g.duplicate(), &old).is_some());
    }

    #[test]
    fn the_hash_is_stable_and_changes_with_any_input() {
        let (y, g) = fixture();
        let plan = |y: &Bits, g: &GroupIds, scheme, k, inner_k, seed, repeat| {
            let p = fold_plan(y, g, scheme, k, inner_k, seed, repeat).unwrap();
            p.verify().unwrap();
            p
        };
        let base = plan(&y, &g, Random, 5, 3, 42, 0);
        assert_eq!(base, plan(&y, &g, Random, 5, 3, 42, 0));
        assert_eq!(base.sha256.len(), 64);
        let mut flipped = y.clone();
        flip(&mut flipped, 0);
        let fewer_dups =
            GroupIds::new(g.lineage().to_vec(), g.duplicate()[1..].to_vec(), None).unwrap();
        let merged = GroupIds::new(
            g.lineage().iter().map(|l| l / 2).collect(),
            g.duplicate().to_vec(),
            None,
        )
        .unwrap();
        let variants = [
            plan(&y, &g, Lineage, 5, 3, 42, 0),
            plan(&y, &g, Random, 4, 3, 42, 0),
            plan(&y, &g, Random, 5, 2, 42, 0),
            plan(&y, &g, Random, 5, 3, 43, 0),
            plan(&y, &g, Random, 5, 3, 42, 1),
            plan(&flipped, &g, Random, 5, 3, 42, 0),
            plan(&y, &fewer_dups, Random, 5, 3, 42, 0),
            plan(&y, &merged, Lineage, 5, 3, 42, 0),
        ];
        let mut seen = BTreeSet::from([base.sha256.clone()]);
        for (j, v) in variants.iter().enumerate() {
            assert!(seen.insert(v.sha256.clone()), "variant {j} repeats a hash");
        }
        let mut tampered = base.clone();
        tampered.outer[0] = (tampered.outer[0] + 1) % 5;
        assert!(tampered.verify().is_err());
    }

    #[test]
    fn the_hash_bytes_match_python_json_dumps() {
        let plan = FoldPlan {
            format: PLAN_FORMAT.to_string(),
            scheme: Lineage,
            repeat: 1,
            k: 2,
            inner_k: 2,
            seed: 42,
            repeat_seed: u64::MAX,
            outer: vec![0, 1, 0, 1],
            inner: vec![
                vec![None, Some(0), None, Some(1)],
                vec![Some(1), None, Some(0), None],
            ],
            holdouts: vec![
                TuningHoldouts {
                    gbt: vec![1],
                    lasso: vec![3],
                },
                TuningHoldouts {
                    gbt: vec![0],
                    lasso: vec![2],
                },
            ],
            unavailable_reason: Some("\u{e9} \"q\"\n\u{1f600}\u{7f}".to_string()),
            sha256: "ignored".to_string(),
        };
        // Python: json.dumps(obj, sort_keys=True, separators=(",", ":"),
        // ensure_ascii=True) and hashlib.sha256 of its ASCII bytes.
        assert_eq!(
            plan.canonical_json().unwrap(),
            r#"{"format":"grm-fold-plan-v1","holdouts":[{"gbt":[1],"lasso":[3]},{"gbt":[0],"lasso":[2]}],"inner":[[null,0,null,1],[1,null,0,null]],"inner_k":2,"k":2,"outer":[0,1,0,1],"repeat":1,"repeat_seed":18446744073709551615,"scheme":"lineage","seed":42,"unavailable_reason":"\u00e9 \"q\"\n\ud83d\ude00\u007f"}"#
        );
        assert_eq!(
            plan.compute_sha256().unwrap(),
            "0ee6b051c60f64ac3366b392b1f845f18138d27590f980e32eec2e543d8a4cbf"
        );
        assert!(canonical_json(&serde_json::json!({"x": 0.5})).is_err());
    }

    #[test]
    fn too_few_lineages_leave_the_lineage_split_unavailable() {
        let n = 40;
        let y = Bits::from_indices(n, (0..n).filter(|i| i % 3 == 0));
        let with = |lineages: u32| {
            let g =
                GroupIds::new((0..n as u32).map(|i| i % lineages).collect(), vec![], None).unwrap();
            fold_plan(&y, &g, Lineage, 5, 3, 42, 0).unwrap()
        };
        let nine = with(9);
        let reason = nine.unavailable_reason.as_deref().unwrap();
        assert!(
            reason.contains("fewer than twice as many lineages as requested folds"),
            "{reason}"
        );
        assert!(nine.inner.is_empty() && nine.holdouts.is_empty());
        nine.verify().unwrap();
        let ten = with(10);
        assert!(
            !ten.unavailable_reason
                .as_deref()
                .is_some_and(|r| r.contains("fewer than twice")),
            "{:?}",
            ten.unavailable_reason
        );
        let g = GroupIds::new((0..n as u32).map(|i| i % 9).collect(), vec![], None).unwrap();
        assert!(
            fold_plan(&y, &g, Random, 5, 3, 42, 0)
                .unwrap()
                .is_available()
        );
    }

    #[test]
    fn a_fold_missing_a_class_leaves_the_split_unavailable() {
        // Three resistant samples cannot reach five folds.
        let n = 40;
        let y = Bits::from_indices(n, [0, 1, 2]);
        let g = GroupIds::new((0..n as u32).map(|i| i / 2).collect(), vec![], None).unwrap();
        for scheme in [Random, Lineage] {
            let plan = fold_plan(&y, &g, scheme, 5, 3, 42, 0).unwrap();
            let reason = plan.unavailable_reason.as_deref().unwrap();
            assert!(
                reason.contains(&format!("invalid {scheme} outer evaluation"))
                    && reason.contains("lacks one or both classes"),
                "{reason}"
            );
        }
    }

    #[test]
    fn repeats_differ_and_are_deterministic() {
        assert_eq!(repeat_seed(42, 0), 42);
        // Python: int.from_bytes(hashlib.sha256(
        //     b"grm-fold-plan-repeat:42:1").digest()[:8], "big")
        assert_eq!(repeat_seed(42, 1), 10622993013839509822);
        assert_eq!(repeat_seed(42, 2), 6315510376195533293);
        assert_eq!(fold_seed(42, 3), 345);
        assert_eq!(fold_seed(u64::MAX, 1), 100);
        let (y, g) = fixture();
        for scheme in [Random, Lineage] {
            let plans: Vec<FoldPlan> = (0..3)
                .map(|r| fold_plan(&y, &g, scheme, 5, 3, 42, r).unwrap())
                .collect();
            for (r, p) in plans.iter().enumerate() {
                assert!(p.is_available(), "{scheme} {r}: {:?}", p.unavailable_reason);
                assert_eq!(*p, fold_plan(&y, &g, scheme, 5, 3, 42, r as u32).unwrap());
            }
            for (a, b) in [(0, 1), (0, 2), (1, 2)] {
                assert_ne!(plans[a].outer, plans[b].outer, "{scheme} {a} {b}");
                assert_ne!(plans[a].inner, plans[b].inner, "{scheme} {a} {b}");
            }
        }
    }

    #[test]
    fn outer_test_labels_cannot_change_the_inner_plans() {
        let (y, g) = fixture();
        let no_dups = GroupIds::new(g.lineage().to_vec(), vec![], None).unwrap();
        for (scheme, groups) in [(Random, &g), (Random, &no_dups), (Lineage, &g)] {
            let plan = fold_plan(&y, groups, scheme, 5, 3, 42, 0).unwrap();
            let tuning = tuning_groups(scheme, groups.lineage(), groups.duplicate());
            for (fi, (train, test)) in fold_masks(&plan.outer, 5).into_iter().enumerate() {
                let seed = fold_seed(plan.repeat_seed, fi);
                let mut flipped = y.clone();
                for i in test.iter_ones() {
                    flip(&mut flipped, i);
                }
                for relabeled in [flipped, y.or(&test), y.and_not(&test)] {
                    let (inner, held) = fold_tuning_plan(
                        &relabeled,
                        &train,
                        tuning.as_deref(),
                        3,
                        seed,
                        scheme,
                        fi,
                    )
                    .unwrap();
                    assert_eq!(inner, plan.inner[fi], "{scheme} outer fold {fi}");
                    assert_eq!(held, plan.holdouts[fi], "{scheme} outer fold {fi}");
                }
            }
            // Control: a training label does reach the inner plan.
            let (train, _) = fold_masks(&plan.outer, 5).swap_remove(0);
            let mut relabeled = y.clone();
            for i in train.iter_ones().take(8) {
                flip(&mut relabeled, i);
            }
            let changed = fold_tuning_plan(
                &relabeled,
                &train,
                tuning.as_deref(),
                3,
                fold_seed(42, 0),
                scheme,
                0,
            );
            if let Ok((inner, held)) = changed {
                assert!(inner != plan.inner[0] || held != plan.holdouts[0]);
            }
        }
    }

    #[test]
    fn group_ids_reject_malformed_groups() {
        let lineage = vec![0, 0, 0, 1];
        // A singleton, a descending pair, an index out of range, a sample in
        // two groups, and a pair across lineages.
        for duplicate in [
            vec![vec![1]],
            vec![vec![2, 1]],
            vec![vec![0, 4]],
            vec![vec![0, 1], vec![1, 2]],
            vec![vec![2, 3]],
        ] {
            assert!(GroupIds::new(lineage.clone(), duplicate, None).is_err());
        }
        assert!(GroupIds::new(lineage.clone(), vec![], Some(vec![0, 0, 0])).is_err());
        let ok = GroupIds::new(lineage.clone(), vec![vec![0, 2]], Some(vec![5, 5, 6, 6])).unwrap();
        assert_eq!(ok.source(), Some(&[5, 5, 6, 6][..]));
        let (y, _) = fixture();
        assert!(fold_plan(&y, &ok, Random, 5, 3, 42, 0).is_err());
    }

    #[test]
    fn pipeline_learner_costs_follow_the_tuning_loops() {
        let cost = |l| LearnerCost::pipeline(l, 3).unwrap().inner_candidates;
        assert_eq!(
            [
                cost("scm"),
                cost("cart"),
                cost("rf"),
                cost("gbt"),
                cost("lasso")
            ],
            [42, 3, 0, 0, 1]
        );
        assert!(LearnerCost::pipeline("svm", 3).is_err());
        for (learner, inner_k) in [("scm", 0), ("cart", 1), ("lasso", 1)] {
            assert!(
                LearnerCost::pipeline(learner, inner_k).is_err(),
                "{learner} {inner_k}"
            );
        }
        let error = LearnerCost::pipeline("scm", usize::MAX).unwrap_err();
        assert!(
            format!("{error:#}").contains("fit units overflow u64"),
            "{error:#}"
        );
        assert_eq!(
            LearnerCost::pipeline("cart", usize::MAX)
                .unwrap()
                .inner_candidates,
            usize::MAX as u64
        );
    }

    fn protocol(budget: u64) -> Protocol {
        Protocol {
            tiers: vec![Tier {
                name: "smoke".to_string(),
                budget,
                arms: vec!["all".to_string(), "no_amr".to_string()],
                learners: vec![
                    LearnerCost::pipeline("scm", 3).unwrap(),
                    LearnerCost::pipeline("rf", 3).unwrap(),
                ],
                schemes: vec![Random, Lineage],
                repeats: 2,
                k: 5,
            }],
        }
    }

    fn inventory() -> Inventory {
        let dataset = |name: &str, unavailable| InventoryDataset {
            name: name.to_string(),
            tiers: vec!["smoke".to_string()],
            unavailable,
        };
        Inventory {
            datasets: vec![dataset("d1", vec![]), dataset("d2", vec![Lineage])],
        }
    }

    #[test]
    fn plan_compute_matches_a_hand_count() {
        // Per (arm, learner, repeat): random runs d1 and d2, lineage only d1,
        // so 3 plans. scm: 5 x (42 + 1) = 215 a plan; rf: 5 x (0 + 1) = 5.
        // 2 arms x 2 repeats x 3 plans x (215 + 5) = 2,640.
        let plan = plan_compute(&protocol(2640), &inventory()).unwrap();
        assert_eq!(plan.rows.len(), 2 * 2 * 2 * 3);
        assert_eq!(plan.total_units, 2640);
        assert_eq!(
            plan.tiers,
            vec![TierCompute {
                tier: "smoke".to_string(),
                units: 2640,
                budget: 2640
            }]
        );
        let scm: Vec<&ComputeRow> = plan.rows.iter().filter(|r| r.learner == "scm").collect();
        assert_eq!(scm.len(), 12);
        assert!(scm.iter().all(|r| r.units == 215));
        assert!(
            !plan
                .rows
                .iter()
                .any(|r| r.dataset == "d2" && r.scheme == Lineage)
        );
        assert_eq!(
            plan.rows[0],
            ComputeRow {
                tier: "smoke".to_string(),
                arm: "all".to_string(),
                learner: "scm".to_string(),
                scheme: Random,
                repeat: 0,
                dataset: "d1".to_string(),
                units: 215,
            }
        );
    }

    #[test]
    fn plan_compute_refuses_over_budget_and_malformed_protocols() {
        let error = plan_compute(&protocol(2639), &inventory()).unwrap_err();
        assert!(
            format!("{error:#}")
                .contains("tier smoke needs 2640 fit units, over its budget of 2639"),
            "{error:#}"
        );
        let mut unknown = inventory();
        unknown.datasets[0].tiers.push("full".to_string());
        assert!(plan_compute(&protocol(10_000), &unknown).is_err());
        let mut twice = protocol(10_000);
        twice.tiers[0].arms.push("all".to_string());
        assert!(plan_compute(&twice, &inventory()).is_err());
    }

    /// The reason `learn::eval::evaluate` gives for `scheme` failing at the
    /// pipeline's seed. Both come from `check_split_plan`: a random failure
    /// is `evaluate`'s error, and a lineage failure is its split's
    /// unavailable reason.
    ///
    /// The toy matrix meets `evaluate`'s final-model precondition: its rows
    /// differ from each other and from `y`, and each has support of at least
    /// `MIN_SUPPORT` on both sides over all samples.
    fn pipeline_reason(y: &Bits, g: &GroupIds, scheme: Scheme, cfg: &EvalConfig) -> String {
        const MIN_SUPPORT: u32 = 3;
        let n = y.len();
        let patterns: [&dyn Fn(usize) -> bool; 4] = [
            &|i| i % 5 < 2,
            &|i| i % 7 < 3,
            &|i| (i / 3) % 2 == 0,
            &|i| i % 11 < 4,
        ];
        let mut rows: Vec<Bits> = vec![];
        for pattern in patterns {
            let row = Bits::from_indices(n, (0..n).filter(|&i| pattern(i)));
            let support = row.count() as usize;
            assert!(
                row != *y
                    && !rows.contains(&row)
                    && support >= MIN_SUPPORT as usize
                    && support + MIN_SUPPORT as usize <= n,
                "toy feature row is a duplicate or outside the support range"
            );
            rows.push(row);
        }
        let mut x = BitMatrix::with_cols(n);
        for row in &rows {
            x.push_row(row.words());
        }
        assert_eq!(
            crate::learn::training_candidates(&x, &Bits::ones(n), MIN_SUPPORT)
                .1
                .len(),
            x.rows
        );
        match evaluate(
            "plan_reason",
            &x,
            y,
            g.lineage(),
            g.duplicate(),
            MIN_SUPPORT,
            None,
            cfg,
        ) {
            Err(error) => {
                assert_eq!(
                    scheme, Random,
                    "evaluate failed before it could report the lineage split: {error:#}"
                );
                format!("{error:#}")
            }
            Ok(out) => out
                .splits
                .iter()
                .find(|s| s.scheme == scheme.as_str())
                .and_then(|s| s.unavailable_reason.clone())
                .unwrap_or_else(|| panic!("the pipeline's {scheme} split is available")),
        }
    }

    #[test]
    fn unavailable_reasons_are_the_pipelines_words() {
        // Small models: the lineage cases fit the random split first.
        let mut cfg = EvalConfig {
            scm_max_rules: 1,
            cart_max_depth: 1,
            screened_features: 2,
            lasso_lambdas: 3,
            ..Default::default()
        };
        cfg.forest.n_trees = 2;
        cfg.forest.max_depth = 2;
        cfg.gbt.rounds = 3;
        cfg.gbt.max_depth = 1;
        let labels =
            |n: usize, r: &dyn Fn(usize) -> bool| Bits::from_indices(n, (0..n).filter(|&i| r(i)));
        let lineages = |n: usize, l: &dyn Fn(usize) -> u32| {
            GroupIds::new((0..n).map(l).collect(), vec![], None).unwrap()
        };
        let cases = [
            // Nine lineages for five folds.
            (
                labels(40, &|i| i % 3 == 0),
                lineages(40, &|i| (i % 9) as u32),
                Lineage,
                3,
                "fewer than twice as many lineages as requested folds",
            ),
            // Three resistant samples for five folds.
            (
                labels(40, &|i| i < 3),
                lineages(40, &|i| (i / 2) as u32),
                Random,
                3,
                "invalid random outer evaluation: ",
            ),
            // Eight resistant training samples for nine inner folds.
            (
                labels(40, &|i| i % 4 == 0),
                lineages(40, &|i| i as u32),
                Random,
                9,
                "invalid random inner evaluation in outer fold 0: ",
            ),
            // Two all-resistant lineages for five folds.
            (
                labels(120, &|i| i < 20),
                lineages(120, &|i| (i / 10) as u32),
                Lineage,
                3,
                "invalid lineage outer evaluation: ",
            ),
            // At most ten training lineages for twelve inner folds.
            (
                labels(120, &|i| i % 2 == 0),
                lineages(120, &|i| (i / 10) as u32),
                Lineage,
                12,
                "invalid lineage inner evaluation in outer fold 0: ",
            ),
        ];
        for (y, g, scheme, inner_k, prefix) in cases {
            cfg.inner_folds = inner_k;
            let plan = fold_plan(&y, &g, scheme, cfg.folds, inner_k, cfg.seed, 0).unwrap();
            let reason = plan.unavailable_reason.clone().unwrap();
            assert!(reason.starts_with(prefix), "{reason}");
            assert_eq!(reason, pipeline_reason(&y, &g, scheme, &cfg));
            assert!(plan.inner.is_empty() && plan.holdouts.is_empty());
            assert!(!plan.is_available());
            plan.verify().unwrap();
        }
    }

    #[test]
    fn plan_files_with_unknown_keys_are_refused() {
        let (y, g) = fixture();
        let plan = fold_plan(&y, &g, Random, 5, 3, 42, 0).unwrap();
        let value = serde_json::to_value(&plan).unwrap();
        let back: FoldPlan = serde_json::from_value(value.clone()).unwrap();
        back.verify().unwrap();
        assert_eq!(back, plan);
        let mut extra = value.clone();
        extra["note"] = serde_json::json!("x");
        assert!(serde_json::from_value::<FoldPlan>(extra).is_err());
        let mut nested = value;
        nested["holdouts"][0]["note"] = serde_json::json!("x");
        assert!(serde_json::from_value::<FoldPlan>(nested).is_err());
    }

    #[test]
    fn verify_and_is_available_check_the_plan_structure() {
        let (y, g) = fixture();
        let available = fold_plan(&y, &g, Random, 5, 3, 42, 0).unwrap();
        let n = 40;
        let few = GroupIds::new((0..n as u32).map(|i| i % 9).collect(), vec![], None).unwrap();
        let unavailable = fold_plan(
            &Bits::from_indices(n, (0..n).filter(|i| i % 3 == 0)),
            &few,
            Lineage,
            5,
            3,
            42,
            0,
        )
        .unwrap();
        for plan in [&available, &unavailable] {
            plan.verify().unwrap();
        }
        assert!(available.is_available() && !unavailable.is_available());
        // Every edit is re-hashed, so only the structure check can catch it.
        let broken = |plan: &FoldPlan, what: &str, edit: &dyn Fn(&mut FoldPlan)| {
            let mut plan = plan.clone();
            edit(&mut plan);
            plan.sha256 = plan.compute_sha256().unwrap();
            assert!(plan.verify().is_err(), "{what}");
            assert!(!plan.is_available(), "{what}");
        };
        let test0 = available.outer.iter().position(|&f| f == 0).unwrap();
        let train0 = available.outer.iter().position(|&f| f != 0).unwrap();
        broken(&unavailable, "reason nulled", &|p| {
            p.unavailable_reason = None
        });
        broken(&unavailable, "reason emptied", &|p| {
            p.unavailable_reason = Some(String::new())
        });
        broken(&available, "reason added", &|p| {
            p.unavailable_reason = Some("x".to_string())
        });
        broken(&available, "an inner plan dropped", &|p| {
            p.inner.pop();
        });
        broken(&available, "a holdout pair dropped", &|p| {
            p.holdouts.pop();
        });
        broken(&available, "an inner plan shortened", &|p| {
            p.inner[0].pop();
        });
        broken(&available, "a test sample given an inner fold", &|p| {
            p.inner[0][test0] = Some(0)
        });
        broken(&available, "a training sample's inner fold dropped", &|p| {
            p.inner[0][train0] = None
        });
        broken(&available, "an inner fold past inner_k", &|p| {
            p.inner[0][train0] = Some(3)
        });
        broken(&available, "a test sample in a holdout", &|p| {
            p.holdouts[0].lasso.push(test0);
            p.holdouts[0].lasso.sort();
        });
        broken(&available, "a holdout out of order", &|p| {
            p.holdouts[0].gbt.reverse()
        });
        broken(&available, "an outer fold past k", &|p| p.outer[0] = 5);
        broken(&available, "repeat_seed changed", &|p| p.repeat_seed += 1);
        broken(&available, "format changed", &|p| {
            p.format = "grm-fold-plan-v0".to_string()
        });
    }

    fn cart_tier(name: &str, budget: u64) -> Tier {
        Tier {
            name: name.to_string(),
            budget,
            arms: vec!["all".to_string()],
            learners: vec![LearnerCost::pipeline("cart", 3).unwrap()],
            schemes: vec![Random],
            repeats: 1,
            k: 5,
        }
    }

    fn one_dataset(tiers: &[&str]) -> Inventory {
        Inventory {
            datasets: vec![InventoryDataset {
                name: "d".to_string(),
                tiers: tiers.iter().map(|t| t.to_string()).collect(),
                unavailable: vec![],
            }],
        }
    }

    #[test]
    fn plan_compute_names_every_tier_over_budget_and_refuses_bad_counts() {
        // Each cart plan costs 5 x (3 + 1) = 20 units.
        let three = Protocol {
            tiers: vec![cart_tier("a", 19), cart_tier("b", 20), cart_tier("c", 1)],
        };
        let error = format!(
            "{:#}",
            plan_compute(&three, &one_dataset(&["a", "b", "c"])).unwrap_err()
        );
        assert!(
            error.contains("tier a needs 20 fit units, over its budget of 19")
                && error.contains("tier c needs 20 fit units, over its budget of 1")
                && !error.contains("tier b"),
            "{error}"
        );
        let refused = |expected: &str, edit: &dyn Fn(&mut Tier)| {
            let mut protocol = Protocol {
                tiers: vec![cart_tier("a", u64::MAX)],
            };
            assert!(plan_compute(&protocol, &one_dataset(&["a"])).is_ok());
            edit(&mut protocol.tiers[0]);
            let error = plan_compute(&protocol, &one_dataset(&["a"])).unwrap_err();
            assert!(
                format!("{error:#}").contains(expected),
                "{expected}: {error:#}"
            );
        };
        refused("needs at least two folds", &|t| t.k = 1);
        refused("runs no repeat", &|t| t.repeats = 0);
        refused("a arm has an empty name", &|t| t.arms = vec![String::new()]);
        refused("no arm is given", &|t| t.arms.clear());
        refused("a learner has an empty name", &|t| {
            t.learners[0].learner.clear()
        });
        refused("fit units overflow u64", &|t| {
            t.learners[0].inner_candidates = u64::MAX
        });
        // Each plan costs u64::MAX, so a second repeat overflows the tier.
        refused("fit units overflow u64", &|t| {
            t.learners[0].inner_candidates = u64::MAX / 5 - 1;
            t.repeats = 2;
        });
        // Two tiers of u64::MAX units each overflow the total.
        let mut full = cart_tier("a", u64::MAX);
        full.learners[0].inner_candidates = u64::MAX / 5 - 1;
        let both = Protocol {
            tiers: vec![
                full.clone(),
                Tier {
                    name: "b".to_string(),
                    ..full
                },
            ],
        };
        let error = plan_compute(&both, &one_dataset(&["a", "b"])).unwrap_err();
        assert!(
            format!("{error:#}").contains("fit units overflow u64"),
            "{error:#}"
        );
        let unnamed = Inventory {
            datasets: vec![InventoryDataset {
                name: String::new(),
                ..one_dataset(&["a"]).datasets.remove(0)
            }],
        };
        let error = plan_compute(
            &Protocol {
                tiers: vec![cart_tier("a", 100)],
            },
            &unnamed,
        )
        .unwrap_err();
        assert!(
            format!("{error:#}").contains("a dataset has an empty name"),
            "{error:#}"
        );
    }

    /// The tuning ids and holdout shares G used before integration, kept as
    /// the reference that Task A's functions must reproduce.
    fn pre_integration_tuning_groups(
        scheme: Scheme,
        lineage: &[u32],
        duplicates: &[Vec<usize>],
    ) -> Option<Vec<u32>> {
        match scheme {
            Scheme::Lineage => Some(lineage.to_vec()),
            Scheme::Random if duplicates.is_empty() => None,
            Scheme::Random => {
                let mut ids: Vec<u32> = (0..lineage.len() as u32).collect();
                for group in duplicates {
                    if let Some(&lowest) = group.iter().min() {
                        for &i in group {
                            ids[i] = lowest as u32;
                        }
                    }
                }
                Some(ids)
            }
        }
    }

    #[test]
    fn tuning_groups_and_holdout_shares_are_task_a_s() {
        assert_eq!((GBT_HOLDOUT, LASSO_HOLDOUT), (0.15, 0.2));
        // xorshift64, so the cases need no extra crate.
        let mut state = 0x9e37_79b9_7f4a_7c15_u64;
        let mut next = |m: usize| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state % m as u64) as usize
        };
        for case in 0..500 {
            let n = 2 + next(40);
            let lineage: Vec<u32> = (0..n).map(|_| next(5) as u32).collect();
            let mut duplicate = Vec::new();
            for l in 0..5 {
                let members: Vec<usize> = (0..n).filter(|&i| lineage[i] == l).collect();
                let mut rest = members.as_slice();
                while rest.len() >= 2 {
                    let take = 1 + next(rest.len());
                    if take >= 2 && next(2) == 0 {
                        duplicate.push(rest[..take].to_vec());
                    }
                    rest = &rest[take..];
                }
            }
            let groups = GroupIds::new(lineage, duplicate, None).unwrap();
            for scheme in [Random, Lineage] {
                assert_eq!(
                    tuning_groups(scheme, groups.lineage(), groups.duplicate()),
                    pre_integration_tuning_groups(scheme, groups.lineage(), groups.duplicate()),
                    "case {case}, {scheme}"
                );
            }
        }
    }

    #[test]
    fn a_tuning_holdout_failure_reads_as_check_split_plan_words_it() {
        // Task A's `a_grouped_holdout_that_fails_after_its_inner_folds_pass_
        // stops_the_plan`, through the planner: outer fold 0 tests two
        // resistant and two susceptible genomes, and its inner folds pass,
        // but six groups of 100 fail boosting's grouped holdout and four fail
        // the lasso's. The plan gives `evaluate`'s reason, naming the unit.
        // Both holdouts fail in the same words, so each case also draws both
        // on their own: only the named one fails, so the plan's error is its.
        for (big, failing) in [(6, "gbt"), (4, "lasso")] {
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
            let mut ids: Vec<u32> = (0..n as u32).collect();
            for g in &duplicates {
                for &i in g {
                    ids[i] = g[0] as u32;
                }
            }
            let outer: Vec<usize> = (0..n).map(|i| usize::from(i >= 4)).collect();
            for (scheme, listed, reason) in [
                (Random, duplicates, "no group tuning partition"),
                (Lineage, vec![], "no lineage tuning partition"),
            ] {
                let groups = GroupIds::new(ids.clone(), listed, None).unwrap();
                let error = nested_plan(&y, &groups, scheme, &outer, 2, 2, 42).unwrap_err();
                assert_eq!(
                    format!("{error:#}"),
                    format!(
                        "invalid {scheme} tuning holdout in outer fold 0: \
                         {reason} contains both classes"
                    ),
                    "{big} groups of 100 fail the {failing} holdout"
                );
                let train = Bits::from_indices(n, 4..n);
                let tuning = tuning_groups(scheme, groups.lineage(), groups.duplicate());
                let unit = crate::learn::eval::tuning_unit(scheme.as_str());
                let seed = fold_seed(42, 0);
                for (name, frac) in [("gbt", GBT_HOLDOUT), ("lasso", LASSO_HOLDOUT)] {
                    let held = tuning_holdout(&train, &y, tuning.as_deref(), unit, frac, seed);
                    match held {
                        Err(error) => assert_eq!(
                            (name, format!("{error:#}")),
                            (failing, format!("{reason} contains both classes")),
                            "{scheme}: {big} groups of 100 fail the {name} holdout"
                        ),
                        Ok(_) => assert_ne!(
                            name, failing,
                            "{scheme}: {big} groups of 100 pass the {name} holdout"
                        ),
                    }
                }
            }
        }
    }

    #[test]
    fn check_structure_refuses_fewer_than_two_outer_or_inner_folds() {
        // An unavailable plan has no inner folds or holdouts, so only the
        // fold-count floor refuses a re-hashed k = 1 or inner_k = 1.
        let (y, g) = fixture();
        let mut base = fold_plan(&y, &g, Random, 5, 3, 42, 0).unwrap();
        base.unavailable_reason = Some("a reason".to_string());
        base.inner.clear();
        base.holdouts.clear();
        base.sha256 = base.compute_sha256().unwrap();
        base.check_structure().unwrap();
        for (k, inner_k) in [(1, 3), (5, 1)] {
            let mut plan = base.clone();
            plan.k = k;
            plan.inner_k = inner_k;
            if k == 1 {
                plan.outer = vec![0; plan.outer.len()];
            }
            plan.sha256 = plan.compute_sha256().unwrap();
            let error = plan.check_structure().unwrap_err();
            assert_eq!(
                format!("{error:#}"),
                "fold plan needs at least two outer and two inner folds",
                "k {k}, inner_k {inner_k}"
            );
        }
    }
}
