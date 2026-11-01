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
use std::collections::BTreeSet;
use std::fmt::Write as _;

/// The fraction of outer-training samples boosting holds out to stop early.
pub const GBT_HOLDOUT: f64 = 0.15;
/// The fraction of outer-training samples the lasso holds out to choose its
/// penalty.
pub const LASSO_HOLDOUT: f64 = 0.2;
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
///   duplicate group takes its lowest member's index and every other sample
///   keeps its own index.
///
/// To be unified with Task A's `learn::eval::tuning_groups` at integration:
/// both must give the same ids, or plans and evaluations diverge.
pub fn tuning_groups(
    scheme: Scheme,
    lineage: &[u32],
    duplicates: &[Vec<usize>],
) -> Option<Vec<u32>> {
    match scheme {
        Scheme::Lineage => Some(lineage.to_vec()),
        Scheme::Random if duplicates.is_empty() => None,
        Scheme::Random => {
            let id = |i: usize| u32::try_from(i).expect("sample indices fit in u32");
            let mut ids: Vec<u32> = (0..lineage.len()).map(id).collect();
            for group in duplicates {
                if let Some(&lowest) = group.iter().min() {
                    for &i in group {
                        ids[i] = id(lowest);
                    }
                }
            }
            Some(ids)
        }
    }
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
            self.format == PLAN_FORMAT,
            "fold plan format {:?} is not {PLAN_FORMAT}",
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
        tuning_holdout(train, y, tuning, frac, seed)
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
                    let (_, val) = tuning_holdout(&train, &y, tuning, frac, seed).unwrap();
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
}
