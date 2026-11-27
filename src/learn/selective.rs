//! Nested selective calling and abstention (`grm learn --exp-selective <DIR>`,
//! docs/selective-calling.md).
//!
//! Off by default: without `--exp-selective`, `enabled()` is false and
//! `grm learn` writes the exact same reports, fold-file results, and fits as
//! before. With it, each evaluated split (`random` and `lineage` on drawn runs,
//! or the supplied split on `--fold-file` runs) calibrates four auditable
//! selection/abstention policies plus class-conditional Mondrian conformal
//! prediction sets strictly inside each outer training fold `T_f`, with zero
//! reads from outer test fold `E_f`.
//!
//! Two evaluation layers are reported separately in `grm-selective-v1` JSON:
//! - `calibrated_operating_point`: applies thresholds `(t_lo, t_hi, tau_d)` and
//!   joint risk-model parameters fit strictly inside `T_f` directly to `E_f`,
//!   where achieved coverage may shift under lineage or study shift.
//! - `matched_coverage`: retains exactly `floor(c * n_f)` of the same `n_f`
//!   outer-fold rows on every outer fold `f` for every policy using only that
//!   policy's outer confidence scores and a frozen deterministic, label-blind
//!   tie-breaking rule (`tie_rank`).
//!
//! Pairwise Mash distances are computed strictly from cached MinHash sketches
//! (`species/<slug>/sketches.bin` when present, else
//! `species-check/sketches/<genome_id>.json`) via `cluster::distance_matrix`,
//! never from `structure.json`'s species-wide `nearest`:
//! - Calibration distance for `i in T_f`:
//!   `d_inner(i) = min_{j in T_f, inner_fold[j] != inner_fold[i]} D[i, j]`.
//! - Outer test distance for `u in E_f`:
//!   `d_outer(u) = min_{j in T_f} D[u, j]`.
//!
//! When sketches are unavailable, score-only policies run and distance-based
//! policies are marked unavailable with a zero-ID reason.
//!
//! Emitted JSON holds aggregate counts, thresholds, rates, intervals, and
//! metadata only: zero genome IDs and zero per-genome predictions.

use super::eval::{EvalConfig, LEARNERS, LearnerResult, Role, SplitResult, tuning_groups};
use super::fold_masks;
use super::inner_pick::{InnerOof, collect_inner_oof};
use crate::bits::{BitMatrix, Bits};
use crate::cluster;
use crate::fnv::Fnv;
use crate::intervals::{self, ALPHA as INTERVAL_ALPHA};
use crate::layout::Layout;
use crate::species_check::{self, GenomeSketch};
use anyhow::{Context, Result, ensure};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::sync::OnceLock;
use std::time::Instant;

/// Result JSON format identifier.
pub const FORMAT: &str = "grm-selective-v1";
/// Pre-committed primary target coverage for cross-method contrasts.
pub const PRIMARY_COVERAGE: f64 = 0.85;
/// Pre-committed coverage ladder evaluated on every dataset, split, learner, and policy.
pub const COVERAGE_LADDER: [f64; 8] = [1.00, 0.95, 0.90, 0.85, 0.80, 0.70, 0.60, 0.50];
/// Pre-committed error levels for class-conditional Mondrian conformal calibration.
pub const CONFORMAL_ALPHAS: [f64; 6] = [0.01, 0.02, 0.05, 0.10, 0.15, 0.20];
/// Minimum class-conditional denominator below which rates are flagged sparse.
pub const SPARSE_DENOMINATOR_MIN: usize = 25;
/// Two-sided significance level for exact Clopper-Pearson confidence intervals (95%).
pub const CI_ALPHA: f64 = INTERVAL_ALPHA;
/// L2 penalty for the strictly convex monotone joint margin+distance logistic model.
pub const JOINT_L2_LAMBDA: f64 = 0.05;
/// Deterministic, label-blind tie-breaking rule identifier.
pub const TIE_RULE: &str = "fnv1a64(selective-tie/<seed>/<scheme>/<fold>/<sample_pos>)";
/// The four auditable selection/abstention policies in fixed output order.
pub const POLICIES: [&str; 4] = [
    "outcome_blind",
    "score_margin",
    "mash_distance",
    "margin_plus_distance",
];
/// Binary header for optional `species/<slug>/sketches.bin`.
const SKETCHES_BIN_MAGIC: &[u8; 8] = b"GRMSK001";

static SINK: OnceLock<PathBuf> = OnceLock::new();

#[cfg(test)]
type Captured = Option<Vec<(String, String)>>;

#[cfg(test)]
thread_local! {
    static CAPTURED: std::cell::RefCell<Captured> = const { std::cell::RefCell::new(None) };
}

/// Turns the selective-calling experiment on for this process. `dir` is created
/// if missing. A second call fails.
pub fn enable(dir: &Path) -> Result<()> {
    std::fs::create_dir_all(dir).with_context(|| format!("create {}", dir.display()))?;
    SINK.set(dir.to_path_buf())
        .map_err(|_| anyhow::anyhow!("the selective-calling output folder is already set"))
}

/// Whether `enable` has run in this process (or, in a unit test, `capture` is
/// active on this thread).
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

#[cfg(test)]
pub(crate) fn capture<R>(f: impl FnOnce() -> R) -> (R, Vec<(String, String)>) {
    CAPTURED.with(|c| *c.borrow_mut() = Some(Vec::new()));
    let result = f();
    let files = CAPTURED.with(|c| c.borrow_mut().take()).unwrap_or_default();
    (result, files)
}

fn write_result(name: &str, text: String) -> Result<()> {
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
    let dir = SINK.get().expect("selective enabled");
    let path = dir.join(name);
    let tmp = dir.join(format!("{name}.tmp"));
    std::fs::write(&tmp, text).with_context(|| format!("write {}", tmp.display()))?;
    std::fs::rename(&tmp, &path).with_context(|| format!("replace {}", path.display()))
}

/// Output filename for one dataset, split scheme, and evaluation seed.
pub fn file_name(dataset: &str, scheme: &str, seed: u64) -> String {
    format!("{dataset}.{scheme}-s{seed}.selective.json")
}

/// Pre-committed constants embedded in every `grm-selective-v1` file.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct SelectiveConstants {
    pub primary_coverage: f64,
    pub coverage_ladder: Vec<f64>,
    pub conformal_alphas: Vec<f64>,
    pub sparse_denominator_min: usize,
    pub ci_alpha: f64,
    pub joint_l2_lambda: f64,
    pub tie_rule: String,
}

impl SelectiveConstants {
    pub fn frozen() -> Self {
        Self {
            primary_coverage: PRIMARY_COVERAGE,
            coverage_ladder: COVERAGE_LADDER.to_vec(),
            conformal_alphas: CONFORMAL_ALPHAS.to_vec(),
            sparse_denominator_min: SPARSE_DENOMINATOR_MIN,
            ci_alpha: CI_ALPHA,
            joint_l2_lambda: JOINT_L2_LAMBDA,
            tie_rule: TIE_RULE.to_string(),
        }
    }
}

// ------------------------------------------------ sketches & Mash distances

/// Symmetric row-major `n x n` Mash distance matrix over cohort sample indices.
/// Inactive samples (e.g. `Role::Unused` in a fold file) store `f32::NAN` so
/// any accidental read outside active folds is caught immediately.
#[derive(Clone, Debug)]
pub struct DistanceMatrix {
    pub n: usize,
    pub k: usize,
    pub sketch_size: usize,
    pub values: Vec<f32>,
}

impl DistanceMatrix {
    /// Builds an `n x n` symmetric Mash distance matrix from per-sample bottom-s sketches.
    pub fn from_sketches(sketches: &[Vec<u64>], sketch_size: usize, k: usize) -> Self {
        let n = sketches.len();
        let values = cluster::distance_matrix(sketches, sketch_size, k);
        Self {
            n,
            k,
            sketch_size,
            values,
        }
    }

    /// Builds a distance matrix directly from a row-major `n x n` slice (for tests).
    pub fn from_dense(n: usize, k: usize, sketch_size: usize, values: Vec<f32>) -> Result<Self> {
        ensure!(
            values.len() == n * n,
            "distance matrix expects {} entries, got {}",
            n * n,
            values.len()
        );
        Ok(Self {
            n,
            k,
            sketch_size,
            values,
        })
    }

    #[inline]
    pub fn get(&self, i: usize, j: usize) -> f32 {
        self.values[i * self.n + j]
    }

    #[inline]
    pub fn set(&mut self, i: usize, j: usize, v: f32) {
        self.values[i * self.n + j] = v;
    }
}

/// Availability and provenance of the cohort Mash distance matrix.
/// Never records any genome ID.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct DistanceStatus {
    pub available: bool,
    pub source: Option<String>,
    pub k: Option<usize>,
    pub sketch_size: Option<usize>,
    pub unavailable_reason: Option<String>,
}

/// Loaded cohort distance matrix together with its zero-ID status metadata.
#[derive(Clone, Debug)]
pub struct LoadedDistances {
    pub matrix: Option<DistanceMatrix>,
    pub status: DistanceStatus,
}

impl LoadedDistances {
    pub fn unavailable(reason: impl Into<String>) -> Self {
        Self {
            matrix: None,
            status: DistanceStatus {
                available: false,
                source: None,
                k: None,
                sketch_size: None,
                unavailable_reason: Some(reason.into()),
            },
        }
    }

    pub fn from_matrix(matrix: DistanceMatrix, source: impl Into<String>) -> Self {
        let status = DistanceStatus {
            available: true,
            source: Some(source.into()),
            k: Some(matrix.k),
            sketch_size: Some(matrix.sketch_size),
            unavailable_reason: None,
        };
        Self {
            matrix: Some(matrix),
            status,
        }
    }
}

/// Writes an optional `species/<slug>/sketches.bin` sidecar.
pub fn write_sketches_bin(
    path: &Path,
    k: usize,
    sketch_size: usize,
    sketches: &[(String, Vec<u64>)],
) -> Result<()> {
    ensure!(
        k > 0 && sketch_size > 0,
        "sketches.bin requires positive k and sketch_size"
    );
    let mut buf = Vec::new();
    buf.extend_from_slice(SKETCHES_BIN_MAGIC);
    buf.extend_from_slice(&(u32::try_from(k)?).to_le_bytes());
    buf.extend_from_slice(&(u32::try_from(sketch_size)?).to_le_bytes());
    buf.extend_from_slice(&(u32::try_from(sketches.len())?).to_le_bytes());
    for (id, hashes) in sketches {
        let id_bytes = id.as_bytes();
        buf.extend_from_slice(&(u32::try_from(id_bytes.len())?).to_le_bytes());
        buf.extend_from_slice(id_bytes);
        buf.extend_from_slice(&(u32::try_from(hashes.len())?).to_le_bytes());
        for &h in hashes {
            buf.extend_from_slice(&h.to_le_bytes());
        }
    }
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent).with_context(|| format!("create {}", parent.display()))?;
    }
    let tmp = path.with_extension("bin.tmp");
    std::fs::write(&tmp, &buf).with_context(|| format!("write {}", tmp.display()))?;
    std::fs::rename(&tmp, path).with_context(|| format!("replace {}", path.display()))
}

/// Reads an optional `species/<slug>/sketches.bin` sidecar.
pub fn read_sketches_bin(path: &Path) -> Result<(usize, usize, HashMap<String, Vec<u64>>)> {
    let bytes = std::fs::read(path).with_context(|| format!("read {}", path.display()))?;
    ensure!(
        bytes.len() >= 20 && &bytes[..8] == SKETCHES_BIN_MAGIC,
        "sketches sidecar has invalid magic header"
    );
    let mut pos = 8usize;
    let read_u32 = |pos: &mut usize| -> Result<u32> {
        ensure!(*pos + 4 <= bytes.len(), "truncated sketches sidecar");
        let v = u32::from_le_bytes(bytes[*pos..*pos + 4].try_into().unwrap());
        *pos += 4;
        Ok(v)
    };
    let k = read_u32(&mut pos)? as usize;
    let sketch_size = read_u32(&mut pos)? as usize;
    let count = read_u32(&mut pos)? as usize;
    ensure!(
        k > 0 && sketch_size > 0,
        "sketches sidecar has invalid k or sketch_size"
    );
    let mut map = HashMap::with_capacity(count);
    for _ in 0..count {
        let id_len = read_u32(&mut pos)? as usize;
        ensure!(pos + id_len <= bytes.len(), "truncated sketches sidecar id");
        let id = std::str::from_utf8(&bytes[pos..pos + id_len])
            .context("non-UTF-8 genome id in sketches sidecar")?
            .to_string();
        pos += id_len;
        let h_len = read_u32(&mut pos)? as usize;
        let byte_len = h_len
            .checked_mul(8)
            .context("overflow in sketches sidecar length")?;
        ensure!(
            pos + byte_len <= bytes.len(),
            "truncated sketches sidecar hashes"
        );
        let mut hashes = Vec::with_capacity(h_len);
        for &chunk in bytes[pos..pos + byte_len].as_chunks::<8>().0 {
            hashes.push(u64::from_le_bytes(chunk));
        }
        pos += byte_len;
        ensure!(
            hashes.windows(2).all(|w| w[0] < w[1]),
            "sketches sidecar holds unsorted or duplicate hashes"
        );
        map.insert(id, hashes);
    }
    ensure!(pos == bytes.len(), "trailing bytes in sketches sidecar");
    Ok((k, sketch_size, map))
}

fn read_cached_species_check_sketch(layout: &Layout, genome_id: &str) -> Option<Vec<u64>> {
    let path = species_check::sketch_path(layout, genome_id);
    let bytes = std::fs::read(path).ok()?;
    let sketch: GenomeSketch = serde_json::from_slice(&bytes).ok()?;
    if sketch.k == species_check::K
        && sketch.sketch_size == species_check::SKETCH
        && sketch.genome_id == genome_id
        && sketch.sketch.windows(2).all(|w| w[0] < w[1])
    {
        Some(sketch.sketch)
    } else {
        None
    }
}

/// Loads cohort sketches from `species/<slug>/sketches.bin` (if present) or
/// cached `species-check/sketches/<genome_id>.json` (read-only; never sketches
/// FASTAs), and computes the symmetric pairwise Mash distance matrix via
/// `cluster::distance_matrix` over all active cohort genomes (`active`).
/// Missing sketches are handled gracefully by returning `available = false`
/// with a zero-ID reason.
pub fn load_cohort_distances(
    layout: &Layout,
    slug: &str,
    sample_ids: &[String],
    active: &Bits,
) -> Result<LoadedDistances> {
    let n = sample_ids.len();
    ensure!(
        active.len() == n,
        "active mask length {} differs from cohort size {n}",
        active.len()
    );
    let active_indices: Vec<usize> = active.iter_ones().collect();
    if active_indices.is_empty() {
        return Ok(LoadedDistances::unavailable("no active cohort genomes"));
    }

    let sidecar = layout.species_dir(slug).join("sketches.bin");
    let mut k = species_check::K;
    let mut sketch_size = species_check::SKETCH;
    let mut bin_map: Option<HashMap<String, Vec<u64>>> = None;
    let mut used_bin = false;
    let mut used_cache = false;

    if sidecar.exists() {
        match read_sketches_bin(&sidecar) {
            Ok((bk, bs, map)) => {
                k = bk;
                sketch_size = bs;
                bin_map = Some(map);
            }
            Err(_) => {
                return Ok(LoadedDistances::unavailable(
                    "species sketches.bin sidecar is unreadable or malformed",
                ));
            }
        }
    }

    let can_fallback_cache = k == species_check::K && sketch_size == species_check::SKETCH;
    let mut active_sketches = Vec::with_capacity(active_indices.len());
    let mut missing = 0usize;

    for &idx in &active_indices {
        let id = &sample_ids[idx];
        if let Some(h) = bin_map.as_ref().and_then(|m| m.get(id)) {
            active_sketches.push(h.clone());
            used_bin = true;
        } else if can_fallback_cache && let Some(h) = read_cached_species_check_sketch(layout, id) {
            active_sketches.push(h);
            used_cache = true;
        } else {
            missing += 1;
        }
    }

    if missing > 0 {
        return Ok(LoadedDistances::unavailable(format!(
            "{missing} of {} active cohort genomes lack a cached MinHash sketch",
            active_indices.len()
        )));
    }

    let compact = cluster::distance_matrix(&active_sketches, sketch_size, k);
    let m = active_indices.len();
    let mut full = vec![f32::NAN; n * n];
    for (ai, &i) in active_indices.iter().enumerate() {
        for (aj, &j) in active_indices.iter().enumerate() {
            full[i * n + j] = compact[ai * m + aj];
        }
    }
    let source = match (used_bin, used_cache) {
        (true, false) => "species_sketches_bin",
        (false, true) => "species_check_cache",
        (true, true) => "species_sketches_bin+species_check_cache",
        (false, false) => unreachable!(),
    };
    Ok(LoadedDistances::from_matrix(
        DistanceMatrix {
            n,
            k,
            sketch_size,
            values: full,
        },
        source,
    ))
}

/// Fold-safe calibration distances inside outer training fold `train`:
/// `d_inner(i) = min_{j in train, inner_fold[j] != inner_fold[i]} D[i, j]`.
/// Strictly never reads any entry involving samples outside `train`.
pub fn inner_nearest_distances(
    dist: &DistanceMatrix,
    train: &Bits,
    inner_fold: &[usize],
) -> Result<Vec<f64>> {
    ensure!(
        dist.n == train.len() && dist.n == inner_fold.len(),
        "distance matrix and fold dimensions differ"
    );
    let train_idx: Vec<usize> = train.iter_ones().collect();
    let mut out = vec![f64::NAN; dist.n];
    for &i in &train_idx {
        let fi = inner_fold[i];
        ensure!(fi != usize::MAX, "training sample lacks an inner fold");
        let mut best = f32::INFINITY;
        for &j in &train_idx {
            if inner_fold[j] != fi {
                let d = dist.get(i, j);
                ensure!(d.is_finite(), "non-finite distance inside training fold");
                if d < best {
                    best = d;
                }
            }
        }
        ensure!(
            best.is_finite(),
            "no out-of-inner-fold training reference for calibration distance"
        );
        out[i] = best as f64;
    }
    Ok(out)
}

/// Fold-safe outer test distances for `u in test` strictly against outer
/// training references `train`:
/// `d_outer(u) = min_{j in train} D[u, j]`.
/// Strictly never reads `D[u, v]` for `v in test` or unused samples.
pub fn outer_nearest_distances(
    dist: &DistanceMatrix,
    train: &Bits,
    test: &Bits,
) -> Result<Vec<f64>> {
    ensure!(
        dist.n == train.len() && dist.n == test.len(),
        "distance matrix and mask dimensions differ"
    );
    ensure!(
        train.and(test).count() == 0,
        "outer test fold overlaps outer training fold"
    );
    let train_idx: Vec<usize> = train.iter_ones().collect();
    ensure!(
        !train_idx.is_empty(),
        "outer training fold is empty for distance calculation"
    );
    let mut out = vec![f64::NAN; dist.n];
    for u in test.iter_ones() {
        let mut best = f32::INFINITY;
        for &j in &train_idx {
            let d = dist.get(u, j);
            ensure!(d.is_finite(), "non-finite distance between test and train");
            if d < best {
                best = d;
            }
        }
        out[u] = best as f64;
    }
    Ok(out)
}

// ------------------------------------------------ tie-breaking & retention

/// Frozen deterministic, label-blind tie-breaking rank derived solely from
/// `(seed, scheme, fold, sample_pos)`. Smaller values win ties first.
pub fn tie_rank(seed: u64, scheme: &str, fold: usize, sample_pos: usize) -> u64 {
    let mut h = Fnv::new();
    h.bytes(format!("selective-tie/{seed}/{scheme}/{fold}/{sample_pos}").as_bytes());
    h.0
}

/// Deterministic label-blind uniform score in `[0.0, 1.0)` derived from
/// `tie_rank`, higher values retained first (`1.0 - u`, matching smaller
/// `tie_rank` retained first).
pub fn blind_confidence(seed: u64, scheme: &str, fold: usize, sample_pos: usize) -> f64 {
    const SCALE: f64 = 1.0 / ((1u64 << 53) as f64);
    let top53 = tie_rank(seed, scheme, fold, sample_pos) >> 11;
    1.0 - (top53 as f64) * SCALE
}

/// Exact integer `floor(c * n)` avoiding floating-point underflow on ladder
/// percentages (e.g. `0.95 * 140 = 133`, not `132`).
pub fn retention_count(coverage: f64, n: usize) -> usize {
    if coverage >= 1.0 {
        return n;
    }
    if coverage <= 0.0 || n == 0 {
        return 0;
    }
    let pct = (coverage * 100.0).round();
    if (coverage * 100.0 - pct).abs() < 1e-12 && (0.0..=100.0).contains(&pct) {
        ((pct as u128 * n as u128) / 100) as usize
    } else {
        (coverage * n as f64).floor().clamp(0.0, n as f64) as usize
    }
}

/// Selects the exact `floor(c * candidates.len())` kept entries in
/// `matched_coverage` mode using `(confidence DESC, tie_rank ASC, sample_pos ASC)`.
/// Returns a boolean mask aligned with `candidates`.
pub fn select_matched_coverage(
    candidates: &[(usize, f64)], // (sample_pos, confidence)
    target_coverage: f64,
    seed: u64,
    scheme: &str,
    fold: usize,
) -> Vec<bool> {
    let n = candidates.len();
    let keep = retention_count(target_coverage, n);
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| {
        let (pos_a, conf_a) = candidates[a];
        let (pos_b, conf_b) = candidates[b];
        conf_b
            .total_cmp(&conf_a)
            .then_with(|| {
                tie_rank(seed, scheme, fold, pos_a).cmp(&tie_rank(seed, scheme, fold, pos_b))
            })
            .then_with(|| pos_a.cmp(&pos_b))
    });
    let mut kept = vec![false; n];
    for &idx in order.iter().take(keep) {
        kept[idx] = true;
    }
    kept
}

// ------------------------------------------ calibration & joint risk model

/// Threshold calibrated on descending confidence scores over `T_f` to retain
/// `retention_count(c, values.len())` calibration samples.
fn calibrate_confidence_threshold(values: &[f64], target_coverage: f64) -> f64 {
    let m = values.len();
    let k_keep = retention_count(target_coverage, m);
    if m == 0 || k_keep >= m {
        return f64::NEG_INFINITY;
    }
    if k_keep == 0 {
        return f64::INFINITY;
    }
    let mut sorted = values.to_vec();
    sorted.sort_by(|a, b| b.total_cmp(a));
    if sorted[k_keep - 1] > sorted[k_keep] {
        0.5 * (sorted[k_keep - 1] + sorted[k_keep])
    } else {
        sorted[k_keep - 1]
    }
}

/// Threshold `tau_d` calibrated on ascending Mash distances over `T_f` so
/// `d <= tau_d` retains `retention_count(c, values.len())` calibration samples.
fn calibrate_distance_threshold(distances: &[f64], target_coverage: f64) -> f64 {
    let m = distances.len();
    let k_keep = retention_count(target_coverage, m);
    if target_coverage >= 1.0 || m == 0 || k_keep >= m {
        return 1.0;
    }
    if k_keep == 0 {
        return -1.0;
    }
    let mut sorted = distances.to_vec();
    sorted.sort_by(f64::total_cmp);
    if sorted[k_keep - 1] < sorted[k_keep] {
        0.5 * (sorted[k_keep - 1] + sorted[k_keep])
    } else {
        sorted[k_keep - 1]
    }
}

/// Finite-sample Mondrian conformal quantile at significance `alpha` on sorted
/// ascending nonconformity scores. `None` represents `+infinity` (always
/// include the class in the prediction set).
pub fn conformal_quantile(sorted_nc: &[f64], alpha: f64) -> Option<f64> {
    let n = sorted_nc.len();
    if n == 0 {
        return None;
    }
    let rank = (((n + 1) as f64) * (1.0 - alpha)).ceil() as usize;
    if rank == 0 {
        Some(-1.0)
    } else if rank > n {
        None
    } else {
        Some(sorted_nc[rank - 1])
    }
}

/// Parameters of the strictly convex monotone L2-regularized logistic
/// error-risk model fit strictly on outer training fold `T_f`:
/// `P(error | m, d) = sigmoid(intercept - weight_margin * z_m + weight_distance * z_d)`
/// with `weight_margin >= 0` and `weight_distance >= 0`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct JointModelParams {
    pub intercept: f64,
    pub weight_margin: f64,
    pub weight_distance: f64,
    pub mean_margin: f64,
    pub sd_margin: f64,
    pub mean_distance: f64,
    pub sd_distance: f64,
    pub l2_lambda: f64,
    pub fallback_equal_weights: bool,
}

impl JointModelParams {
    /// Standardized confidence score `-eta(m, d)` (higher = lower predicted
    /// error risk = more confident to call).
    pub fn confidence(&self, margin: f64, distance: f64) -> f64 {
        let z_m = if self.sd_margin > 1e-12 {
            (margin - self.mean_margin) / self.sd_margin
        } else {
            0.0
        };
        let z_d = if self.sd_distance > 1e-12 {
            (distance - self.mean_distance) / self.sd_distance
        } else {
            0.0
        };
        self.weight_margin * z_m - self.weight_distance * z_d - self.intercept
    }

    /// Predicted error probability `sigmoid(-confidence(m, d))`.
    pub fn error_prob(&self, margin: f64, distance: f64) -> f64 {
        sigmoid(-self.confidence(margin, distance))
    }

    fn z_margin(&self, margin: f64) -> f64 {
        if self.sd_margin > 1e-12 {
            (margin - self.mean_margin) / self.sd_margin
        } else {
            0.0
        }
    }

    fn z_distance(&self, distance: f64) -> f64 {
        if self.sd_distance > 1e-12 {
            (distance - self.mean_distance) / self.sd_distance
        } else {
            0.0
        }
    }
}

#[inline]
fn sigmoid(x: f64) -> f64 {
    if x >= 0.0 {
        1.0 / (1.0 + (-x).exp())
    } else {
        let e = x.exp();
        e / (1.0 + e)
    }
}

#[inline]
fn log1p_exp(x: f64) -> f64 {
    if x > 0.0 {
        x + (-x).exp().ln_1p()
    } else {
        x.exp().ln_1p()
    }
}

/// Fits the monotone L2-regularized logistic error-risk model strictly on
/// `(errors, margins, distances)` of outer training fold `T_f`.
///
/// Objective:
/// `(1/m) sum_i [ log(1 + exp(eta_i)) - e_i * eta_i ] + (lambda / 2) (w_m^2 + w_d^2)`
/// where `eta_i = b0 - w_m * z_{m,i} + w_d * z_{d,i}`, subject to `w_m >= 0, w_d >= 0`.
pub fn fit_joint_risk_model(
    errors: &[bool],
    margins: &[f64],
    distances: &[f64],
    lambda: f64,
) -> JointModelParams {
    let m = errors.len();
    assert_eq!(m, margins.len());
    assert_eq!(m, distances.len());
    let mf = m.max(1) as f64;

    let mean_margin = margins.iter().sum::<f64>() / mf;
    let var_margin = margins
        .iter()
        .map(|&v| (v - mean_margin).powi(2))
        .sum::<f64>()
        / mf;
    let sd_margin = var_margin.sqrt();

    let mean_distance = distances.iter().sum::<f64>() / mf;
    let var_distance = distances
        .iter()
        .map(|&v| (v - mean_distance).powi(2))
        .sum::<f64>()
        / mf;
    let sd_distance = var_distance.sqrt();

    let has_m = sd_margin > 1e-12;
    let has_d = sd_distance > 1e-12;
    let z_m: Vec<f64> = margins
        .iter()
        .map(|&v| {
            if has_m {
                (v - mean_margin) / sd_margin
            } else {
                0.0
            }
        })
        .collect();
    let z_d: Vec<f64> = distances
        .iter()
        .map(|&v| {
            if has_d {
                (v - mean_distance) / sd_distance
            } else {
                0.0
            }
        })
        .collect();
    let err: Vec<f64> = errors.iter().map(|&e| if e { 1.0 } else { 0.0 }).collect();
    let err_rate = (err.iter().sum::<f64>() / mf).clamp(1e-6, 1.0 - 1e-6);

    let mut b0 = (err_rate / (1.0 - err_rate)).ln();
    let mut w_m = 0.0f64;
    let mut w_d = 0.0f64;

    let objective = |b0: f64, w_m: f64, w_d: f64| -> f64 {
        let mut loss = 0.0;
        for i in 0..m {
            let eta = b0 - w_m * z_m[i] + w_d * z_d[i];
            loss += log1p_exp(eta) - err[i] * eta;
        }
        loss / mf + 0.5 * lambda * (w_m * w_m + w_d * w_d)
    };

    if m > 0 && (has_m || has_d) {
        let mut cur_obj = objective(b0, w_m, w_d);
        for _ in 0..200 {
            let (prev_b0, prev_wm, prev_wd) = (b0, w_m, w_d);

            // 1. Intercept b0 (unconstrained)
            let (mut g0, mut h0) = (0.0, 0.0);
            for i in 0..m {
                let p = sigmoid(b0 - w_m * z_m[i] + w_d * z_d[i]);
                g0 += p - err[i];
                h0 += p * (1.0 - p);
            }
            let step0 = (g0 / mf) / (h0 / mf + 1e-9);
            let mut alpha = 1.0;
            while alpha >= 1e-6 {
                let cand = b0 - alpha * step0;
                let obj = objective(cand, w_m, w_d);
                if obj <= cur_obj + 1e-15 {
                    b0 = cand;
                    cur_obj = obj;
                    break;
                }
                alpha *= 0.5;
            }

            // 2. Margin weight w_m >= 0 (eta decreases with margin)
            if has_m {
                let (mut gm, mut hm) = (0.0, 0.0);
                for i in 0..m {
                    let p = sigmoid(b0 - w_m * z_m[i] + w_d * z_d[i]);
                    gm += -(p - err[i]) * z_m[i];
                    hm += p * (1.0 - p) * z_m[i] * z_m[i];
                }
                let grad = gm / mf + lambda * w_m;
                let hess = hm / mf + lambda;
                let raw_target = (w_m - grad / hess).max(0.0);
                let delta = raw_target - w_m;
                let mut alpha = 1.0;
                while alpha >= 1e-6 {
                    let cand = (w_m + alpha * delta).max(0.0);
                    let obj = objective(b0, cand, w_d);
                    if obj <= cur_obj + 1e-15 {
                        w_m = cand;
                        cur_obj = obj;
                        break;
                    }
                    alpha *= 0.5;
                }
            }

            // 3. Distance weight w_d >= 0 (eta increases with distance)
            if has_d {
                let (mut gd, mut hd) = (0.0, 0.0);
                for i in 0..m {
                    let p = sigmoid(b0 - w_m * z_m[i] + w_d * z_d[i]);
                    gd += (p - err[i]) * z_d[i];
                    hd += p * (1.0 - p) * z_d[i] * z_d[i];
                }
                let grad = gd / mf + lambda * w_d;
                let hess = hd / mf + lambda;
                let raw_target = (w_d - grad / hess).max(0.0);
                let delta = raw_target - w_d;
                let mut alpha = 1.0;
                while alpha >= 1e-6 {
                    let cand = (w_d + alpha * delta).max(0.0);
                    let obj = objective(b0, w_m, cand);
                    if obj <= cur_obj + 1e-15 {
                        w_d = cand;
                        cur_obj = obj;
                        break;
                    }
                    alpha *= 0.5;
                }
            }

            let max_change = (b0 - prev_b0)
                .abs()
                .max((w_m - prev_wm).abs())
                .max((w_d - prev_wd).abs());
            if max_change < 1e-12 {
                break;
            }
        }
    }

    // If neither feature received positive weight (e.g. zero inner errors or
    // purely degenerate inner error labels), fall back deterministically to
    // equal positive standardized weights on whichever features vary on T_f.
    let mut fallback_equal_weights = false;
    if w_m <= 1e-12 && w_d <= 1e-12 && (has_m || has_d) {
        w_m = if has_m { 1.0 } else { 0.0 };
        w_d = if has_d { 1.0 } else { 0.0 };
        fallback_equal_weights = true;
    }

    JointModelParams {
        intercept: b0,
        weight_margin: w_m,
        weight_distance: w_d,
        mean_margin,
        sd_margin,
        mean_distance,
        sd_distance,
        l2_lambda: lambda,
        fallback_equal_weights,
    }
}

// ------------------------------------------------ accounting & JSON schema

/// Proportion `x / n` with exact two-sided 95% Clopper-Pearson confidence
/// interval (`crate::intervals::clopper_pearson`) and sparse-denominator flag.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct RateCi {
    pub x: usize,
    pub n: usize,
    pub rate: Option<f64>,
    pub ci_lo: Option<f64>,
    pub ci_hi: Option<f64>,
    pub sparse_denominator: bool,
}

impl RateCi {
    pub fn compute(x: usize, n: usize) -> Self {
        let sparse_denominator = n < SPARSE_DENOMINATOR_MIN;
        if n == 0 {
            return Self {
                x,
                n,
                rate: None,
                ci_lo: None,
                ci_hi: None,
                sparse_denominator,
            };
        }
        let (lo, hi) = intervals::clopper_pearson(x as u64, n as u64, CI_ALPHA)
            .expect("x <= n and n > 0 for Clopper-Pearson interval");
        Self {
            x,
            n,
            rate: Some(x as f64 / n as f64),
            ci_lo: Some(lo),
            ci_hi: Some(hi),
            sparse_denominator,
        }
    }
}

/// Breakdown of why samples were abstained. Always sums to
/// `abstain_r + abstain_s`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct AbstainReasons {
    pub null_fit: usize,
    pub ood_mash_distance: usize,
    pub uncertain_score: usize,
    pub conformal_both: usize,
    pub conformal_empty: usize,
    pub random_thinning: usize,
}

impl AbstainReasons {
    pub fn total(&self) -> usize {
        self.null_fit
            + self.ood_mash_distance
            + self.uncertain_score
            + self.conformal_both
            + self.conformal_empty
            + self.random_thinning
    }

    pub fn add_assign(&mut self, other: &Self) {
        self.null_fit += other.null_fit;
        self.ood_mash_distance += other.ood_mash_distance;
        self.uncertain_score += other.uncertain_score;
        self.conformal_both += other.conformal_both;
        self.conformal_empty += other.conformal_empty;
        self.random_thinning += other.random_thinning;
    }
}

/// Per-100 clinical accounting among all eligible resistant isolates.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Per100R {
    pub caught: Option<f64>,
    pub missed: Option<f64>,
    pub deferred: Option<f64>,
}

/// Per-100 clinical accounting among all eligible susceptible isolates.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Per100S {
    pub cleared: Option<f64>,
    pub overcalled: Option<f64>,
    pub deferred: Option<f64>,
}

/// Complete, zero-ID selective calling accounting over a split's tested rows.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct SelectiveAccounting {
    pub n: usize,
    pub n_r: usize,
    pub n_s: usize,
    pub called: usize,
    pub called_r: usize,
    pub called_s: usize,
    pub abstain_r: usize,
    pub abstain_s: usize,
    pub tp: usize,
    pub fp: usize,
    pub tn: usize,
    #[serde(rename = "fn")]
    pub fn_: usize,
    pub abstain_reasons: AbstainReasons,
    pub coverage: Option<f64>,
    pub coverage_r: Option<f64>,
    pub coverage_s: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub set_true_label_coverage: Option<RateCi>,
    pub vme_called: RateCi,
    pub me_called: RateCi,
    pub npv: RateCi,
    pub for_s: RateCi,
    pub ppv: RateCi,
    pub balanced_accuracy_called: Option<f64>,
    pub vme_all: RateCi,
    pub me_all: RateCi,
    pub per_100_r: Per100R,
    pub per_100_s: Per100S,
    /// True when either called-class denominator (`called_r` or `called_s`) is
    /// below `SPARSE_DENOMINATOR_MIN`.
    pub sparse_denominator: bool,
    pub sparse_denominator_r: bool,
    pub sparse_denominator_s: bool,
}

/// Per-fold aggregate counts (zero genome IDs, zero per-genome predictions).
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct FoldSelectiveCounts {
    pub fold: usize,
    pub n: usize,
    pub n_r: usize,
    pub n_s: usize,
    pub called: usize,
    pub called_r: usize,
    pub called_s: usize,
    pub abstain_r: usize,
    pub abstain_s: usize,
    pub tp: usize,
    pub fp: usize,
    pub tn: usize,
    #[serde(rename = "fn")]
    pub fn_: usize,
    pub abstain_reasons: AbstainReasons,
}

/// Raw accumulator for building `SelectiveAccounting` and `FoldSelectiveCounts`.
#[derive(Clone, Debug, Default)]
pub struct CallTally {
    pub n_r: usize,
    pub n_s: usize,
    pub tp: usize,
    pub fp: usize,
    pub tn: usize,
    pub fn_: usize,
    pub abstain_r: usize,
    pub abstain_s: usize,
    pub abstain_reasons: AbstainReasons,
    pub conformal_covered: usize,
    pub conformal_total: usize,
}

impl CallTally {
    pub fn record_call(&mut self, truth_r: bool, pred_r: bool) {
        if truth_r {
            self.n_r += 1;
            if pred_r {
                self.tp += 1;
            } else {
                self.fn_ += 1;
            }
        } else {
            self.n_s += 1;
            if pred_r {
                self.fp += 1;
            } else {
                self.tn += 1;
            }
        }
    }

    pub fn record_abstain(&mut self, truth_r: bool, reason: AbstainReason) {
        if truth_r {
            self.n_r += 1;
            self.abstain_r += 1;
        } else {
            self.n_s += 1;
            self.abstain_s += 1;
        }
        match reason {
            AbstainReason::NullFit => self.abstain_reasons.null_fit += 1,
            AbstainReason::OodMashDistance => self.abstain_reasons.ood_mash_distance += 1,
            AbstainReason::UncertainScore => self.abstain_reasons.uncertain_score += 1,
            AbstainReason::ConformalBoth => self.abstain_reasons.conformal_both += 1,
            AbstainReason::ConformalEmpty => self.abstain_reasons.conformal_empty += 1,
            AbstainReason::RandomThinning => self.abstain_reasons.random_thinning += 1,
        }
    }

    pub fn record_conformal(&mut self, truth_r: bool, has_r: bool, has_s: bool) {
        self.conformal_total += 1;
        let covered = if truth_r { has_r } else { has_s };
        if covered {
            self.conformal_covered += 1;
        }
        match (has_r, has_s) {
            (true, false) => self.record_call(truth_r, true),
            (false, true) => self.record_call(truth_r, false),
            (true, true) => self.record_abstain(truth_r, AbstainReason::ConformalBoth),
            (false, false) => self.record_abstain(truth_r, AbstainReason::ConformalEmpty),
        }
    }

    pub fn add_assign(&mut self, other: &Self) {
        self.n_r += other.n_r;
        self.n_s += other.n_s;
        self.tp += other.tp;
        self.fp += other.fp;
        self.tn += other.tn;
        self.fn_ += other.fn_;
        self.abstain_r += other.abstain_r;
        self.abstain_s += other.abstain_s;
        self.abstain_reasons.add_assign(&other.abstain_reasons);
        self.conformal_covered += other.conformal_covered;
        self.conformal_total += other.conformal_total;
    }

    pub fn fold_counts(&self, fold: usize) -> FoldSelectiveCounts {
        let called_r = self.tp + self.fn_;
        let called_s = self.tn + self.fp;
        FoldSelectiveCounts {
            fold,
            n: self.n_r + self.n_s,
            n_r: self.n_r,
            n_s: self.n_s,
            called: called_r + called_s,
            called_r,
            called_s,
            abstain_r: self.abstain_r,
            abstain_s: self.abstain_s,
            tp: self.tp,
            fp: self.fp,
            tn: self.tn,
            fn_: self.fn_,
            abstain_reasons: self.abstain_reasons,
        }
    }

    pub fn finish(&self, include_conformal: bool) -> SelectiveAccounting {
        let n = self.n_r + self.n_s;
        let called_r = self.tp + self.fn_;
        let called_s = self.tn + self.fp;
        let called = called_r + called_s;
        debug_assert_eq!(called_r + self.abstain_r, self.n_r);
        debug_assert_eq!(called_s + self.abstain_s, self.n_s);
        debug_assert_eq!(
            self.abstain_reasons.total(),
            self.abstain_r + self.abstain_s
        );

        let ratio = |num: usize, den: usize| (den > 0).then(|| num as f64 / den as f64);
        let pct = |num: usize, den: usize| (den > 0).then(|| 100.0 * num as f64 / den as f64);

        let balanced_accuracy_called = (called_r > 0 && called_s > 0)
            .then(|| 0.5 * (self.tp as f64 / called_r as f64 + self.tn as f64 / called_s as f64));
        let set_true_label_coverage = include_conformal
            .then(|| RateCi::compute(self.conformal_covered, self.conformal_total));

        let called_pred_s = self.tn + self.fn_;
        let called_pred_r = self.tp + self.fp;
        let sparse_r = called_r < SPARSE_DENOMINATOR_MIN;
        let sparse_s = called_s < SPARSE_DENOMINATOR_MIN;

        SelectiveAccounting {
            n,
            n_r: self.n_r,
            n_s: self.n_s,
            called,
            called_r,
            called_s,
            abstain_r: self.abstain_r,
            abstain_s: self.abstain_s,
            tp: self.tp,
            fp: self.fp,
            tn: self.tn,
            fn_: self.fn_,
            abstain_reasons: self.abstain_reasons,
            coverage: ratio(called, n),
            coverage_r: ratio(called_r, self.n_r),
            coverage_s: ratio(called_s, self.n_s),
            set_true_label_coverage,
            vme_called: RateCi::compute(self.fn_, called_r),
            me_called: RateCi::compute(self.fp, called_s),
            npv: RateCi::compute(self.tn, called_pred_s),
            for_s: RateCi::compute(self.fn_, called_pred_s),
            ppv: RateCi::compute(self.tp, called_pred_r),
            balanced_accuracy_called,
            vme_all: RateCi::compute(self.fn_, self.n_r),
            me_all: RateCi::compute(self.fp, self.n_s),
            per_100_r: Per100R {
                caught: pct(self.tp, self.n_r),
                missed: pct(self.fn_, self.n_r),
                deferred: pct(self.abstain_r, self.n_r),
            },
            per_100_s: Per100S {
                cleared: pct(self.tn, self.n_s),
                overcalled: pct(self.fp, self.n_s),
                deferred: pct(self.abstain_s, self.n_s),
            },
            sparse_denominator: sparse_r || sparse_s,
            sparse_denominator_r: sparse_r,
            sparse_denominator_s: sparse_s,
        }
    }
}

/// Reason a single sample was deferred/abstained.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AbstainReason {
    NullFit,
    OodMashDistance,
    UncertainScore,
    ConformalBoth,
    ConformalEmpty,
    RandomThinning,
}

/// Thresholds calibrated strictly on outer training fold `T_f` for one target
/// coverage rung.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct CoverageThresholds {
    pub target_coverage: f64,
    pub tau_blind: f64,
    pub t_lo: f64,
    pub t_hi: f64,
    pub tau_margin: f64,
    pub tau_distance: Option<f64>,
    pub tau_joint: Option<f64>,
}

/// Class-conditional Mondrian conformal thresholds calibrated strictly on
/// `s_inner[T_f]` at one significance level `alpha`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ConformalThresholds {
    pub alpha: f64,
    pub q_r: Option<f64>,
    pub q_s: Option<f64>,
    pub t_lo: Option<f64>,
    pub t_hi: Option<f64>,
}

/// Calibration artifacts fit strictly inside outer training fold `T_f` for one
/// learner. Mutating outer test labels or outer test distances must leave every
/// field bit-for-bit identical.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct FoldCalibration {
    pub fold: usize,
    pub train_n: usize,
    pub train_r: usize,
    pub train_s: usize,
    pub inner_null_fits: usize,
    pub inner_unavailable_reason: Option<String>,
    pub thresholds: Vec<CoverageThresholds>,
    pub conformal: Vec<ConformalThresholds>,
    pub joint_model: Option<JointModelParams>,
}

/// Evaluation of one policy at one target coverage rung.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct CoveragePoint {
    pub target_coverage: f64,
    #[serde(flatten)]
    pub accounting: SelectiveAccounting,
    pub folds: Vec<FoldSelectiveCounts>,
}

/// Evaluation of class-conditional Mondrian conformal sets at one error level.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ConformalOperatingPoint {
    pub alpha: f64,
    #[serde(flatten)]
    pub accounting: SelectiveAccounting,
    pub folds: Vec<FoldSelectiveCounts>,
}

/// Dual-layer evaluation of one selective calling policy for one learner.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct PolicySelective {
    pub policy: String,
    pub available: bool,
    pub unavailable_reason: Option<String>,
    /// Out-of-sample evaluation on `E_f` applying thresholds fit strictly on `T_f`.
    pub calibrated_operating_point: Vec<CoveragePoint>,
    /// Exact label-blind matched-coverage evaluation retaining `floor(c * n_f)`
    /// rows on every outer fold `f`.
    pub matched_coverage: Vec<CoveragePoint>,
}

/// Complete selective calling evaluation for one learner on one split.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct LearnerSelective {
    pub learner: String,
    pub outer_null_folds: usize,
    pub calibrations: Vec<FoldCalibration>,
    pub policies: Vec<PolicySelective>,
    pub mondrian_conformal_available: bool,
    pub mondrian_conformal_unavailable_reason: Option<String>,
    pub mondrian_conformal: Vec<ConformalOperatingPoint>,
}

/// Top-level `grm-selective-v1` output file (zero genome IDs, zero per-genome predictions).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct SelectiveResult {
    pub format: String,
    pub evaluation_revision: u32,
    pub dataset: String,
    pub species_slug: String,
    pub antibiotic: String,
    pub scheme: String,
    pub seed: u64,
    pub folds: usize,
    pub inner_folds: usize,
    pub min_support: u32,
    pub samples: usize,
    pub tested: usize,
    pub tested_resistant: usize,
    pub tested_susceptible: usize,
    pub fold_file_sha256: Option<String>,
    pub unavailable_reason: Option<String>,
    pub constants: SelectiveConstants,
    pub distances: DistanceStatus,
    pub learners: Vec<LearnerSelective>,
    pub seconds: f64,
}

// ------------------------------------------------ core fold calibration & eval

fn is_null_outer_fold(learner: &str, result: &LearnerResult, fold: usize) -> bool {
    matches!(learner, "gbt" | "lasso")
        && result.fold_complexity.get(fold).is_some_and(|&c| c == 0.0)
}

/// Calibrates every threshold and the joint margin+distance error-risk model
/// strictly on outer training fold `T_f`, reading ZERO data from outer test
/// fold `E_f`.
#[allow(clippy::too_many_arguments)]
fn calibrate_fold(
    fold: usize,
    scheme: &str,
    seed: u64,
    y: &Bits,
    train: &Bits,
    inner_scores: Option<&[f64]>,
    inner_null_fits: usize,
    inner_error: Option<String>,
    d_inner: Option<&[f64]>,
) -> FoldCalibration {
    let train_idx: Vec<usize> = train.iter_ones().collect();
    let train_n = train_idx.len();
    let train_r = train_idx.iter().filter(|&&i| y.contains(i)).count();
    let train_s = train_n - train_r;

    let blind_vals: Vec<f64> = train_idx
        .iter()
        .map(|&i| blind_confidence(seed, scheme, fold, i))
        .collect();

    let dist_vals: Option<Vec<f64>> = d_inner.map(|d| train_idx.iter().map(|&i| d[i]).collect());

    let margins: Option<Vec<f64>> =
        inner_scores.map(|scores| train_idx.iter().map(|&i| (scores[i] - 0.5).abs()).collect());

    // Joint margin + distance monotone risk model fit strictly on T_f.
    let joint_model = match (inner_scores, &margins, &dist_vals) {
        (Some(scores), Some(m_vals), Some(d_vals)) => {
            let errors: Vec<bool> = train_idx
                .iter()
                .map(|&i| (scores[i] >= 0.5) != y.contains(i))
                .collect();
            Some(fit_joint_risk_model(
                &errors,
                m_vals,
                d_vals,
                JOINT_L2_LAMBDA,
            ))
        }
        _ => None,
    };

    let joint_confs: Option<Vec<f64>> = match (&joint_model, &margins, &dist_vals) {
        (Some(jm), Some(m_vals), Some(d_vals)) => Some(
            m_vals
                .iter()
                .zip(d_vals)
                .map(|(&m, &d)| jm.confidence(m, d))
                .collect(),
        ),
        _ => None,
    };

    let thresholds: Vec<CoverageThresholds> = COVERAGE_LADDER
        .iter()
        .map(|&c| {
            let tau_blind = if c >= 1.0 {
                0.0
            } else {
                let t = calibrate_confidence_threshold(&blind_vals, c);
                if t.is_finite() { t } else { 1.0 }
            };

            let tau_margin = match &margins {
                Some(m_vals) => {
                    if c >= 1.0 {
                        0.0
                    } else {
                        let t = calibrate_confidence_threshold(m_vals, c);
                        if t.is_finite() {
                            t.clamp(0.0, 0.5)
                        } else {
                            0.5
                        }
                    }
                }
                None => 0.0,
            };
            let t_lo = (0.5 - tau_margin).max(0.0);
            let t_hi = (0.5 + tau_margin).min(1.0);

            let tau_distance = dist_vals
                .as_ref()
                .map(|d_vals| calibrate_distance_threshold(d_vals, c));

            let tau_joint = joint_confs.as_ref().map(|j_vals| {
                if c >= 1.0 {
                    j_vals.iter().copied().fold(f64::INFINITY, f64::min) - 1.0
                } else {
                    let t = calibrate_confidence_threshold(j_vals, c);
                    if t.is_finite() {
                        t
                    } else {
                        j_vals.iter().copied().fold(f64::NEG_INFINITY, f64::max) + 1.0
                    }
                }
            });

            CoverageThresholds {
                target_coverage: c,
                tau_blind,
                t_lo,
                t_hi,
                tau_margin,
                tau_distance,
                tau_joint,
            }
        })
        .collect();

    // Class-conditional Mondrian conformal calibration strictly on s_inner[T_f].
    let conformal = if let Some(scores) = inner_scores {
        let mut nc_r: Vec<f64> = train_idx
            .iter()
            .filter(|&&i| y.contains(i))
            .map(|&i| 1.0 - scores[i])
            .collect();
        let mut nc_s: Vec<f64> = train_idx
            .iter()
            .filter(|&&i| !y.contains(i))
            .map(|&i| scores[i])
            .collect();
        nc_r.sort_by(f64::total_cmp);
        nc_s.sort_by(f64::total_cmp);
        CONFORMAL_ALPHAS
            .iter()
            .map(|&alpha| {
                let q_r = conformal_quantile(&nc_r, alpha);
                let q_s = conformal_quantile(&nc_s, alpha);
                ConformalThresholds {
                    alpha,
                    q_r,
                    q_s,
                    t_lo: q_r.map(|q| 1.0 - q),
                    t_hi: q_s,
                }
            })
            .collect()
    } else {
        vec![]
    };

    FoldCalibration {
        fold,
        train_n,
        train_r,
        train_s,
        inner_null_fits,
        inner_unavailable_reason: inner_error,
        thresholds,
        conformal,
        joint_model,
    }
}

/// Evaluates all four policies + Mondrian conformal calibration on one split,
/// given precomputed outer splits, outer learner results, and optional cohort
/// Mash distance matrix.
#[allow(clippy::too_many_arguments)]
pub fn evaluate_split(
    dataset: &str,
    species_slug: &str,
    antibiotic: &str,
    scheme: &str,
    x: &BitMatrix,
    y: &Bits,
    groups: &[u32],
    duplicates: &[Vec<usize>],
    trains: &[Bits],
    tested: &[usize],
    test_fold: &[usize],
    learners: &[LearnerResult],
    unavailable_reason: Option<String>,
    fold_file_sha256: Option<String>,
    min_support: u32,
    cfg: &EvalConfig,
    distances: &LoadedDistances,
) -> Result<SelectiveResult> {
    let started = Instant::now();
    let n = x.cols;
    let mut base = SelectiveResult {
        format: FORMAT.to_string(),
        evaluation_revision: crate::analysis::EVALUATION_REVISION,
        dataset: dataset.to_string(),
        species_slug: species_slug.to_string(),
        antibiotic: antibiotic.to_string(),
        scheme: scheme.to_string(),
        seed: cfg.seed,
        folds: trains.len(),
        inner_folds: cfg.inner_folds,
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
            tested.iter().filter(|&&i| y.contains(i)).count()
        },
        tested_susceptible: if unavailable_reason.is_some() {
            0
        } else {
            tested.iter().filter(|&&i| !y.contains(i)).count()
        },
        fold_file_sha256,
        unavailable_reason,
        constants: SelectiveConstants::frozen(),
        distances: distances.status.clone(),
        learners: vec![],
        seconds: 0.0,
    };

    if base.unavailable_reason.is_some() {
        base.seconds = started.elapsed().as_secs_f64();
        return Ok(base);
    }

    ensure!(
        tested.len() == test_fold.len(),
        "selective evaluation requires one outer fold per tested sample"
    );
    ensure!(
        learners.len() >= LEARNERS.len()
            && LEARNERS.iter().zip(learners).all(|(l, r)| r.learner == *l),
        "split learners are not the five release learners in pipeline order"
    );

    let tuning = tuning_groups(scheme, groups, duplicates);
    let n_folds = trains.len();

    // Per outer fold: collect test positions, validate disjointness from train,
    // collect inner OOF scores strictly inside T_f, and compute fold-safe Mash
    // distances d_inner (inside T_f only) and d_outer (E_f vs T_f only).
    let mut fold_test_pos: Vec<Vec<usize>> = Vec::with_capacity(n_folds);
    let mut fold_oofs: Vec<InnerOof> = Vec::with_capacity(n_folds);
    let mut fold_d_inner: Vec<Option<Vec<f64>>> = Vec::with_capacity(n_folds);
    let mut fold_d_outer: Vec<Option<Vec<f64>>> = Vec::with_capacity(n_folds);
    let mut dist_fold_error: Option<String> = None;

    for (f, train) in trains.iter().enumerate() {
        let positions: Vec<usize> = (0..tested.len()).filter(|&p| test_fold[p] == f).collect();
        let test = Bits::from_indices(n, positions.iter().map(|&p| tested[p]));
        ensure!(
            test.and(train).count() == 0,
            "outer fold {f} trains on a sample it tests"
        );
        let oof = collect_inner_oof(
            dataset,
            scheme,
            x,
            y,
            tuning.as_deref(),
            tuning.as_deref(),
            train,
            f,
            cfg.inner_folds,
            min_support,
            cfg,
        );
        if let Some(dm) = distances.matrix.as_ref() {
            if oof.reason.is_none() {
                match (
                    inner_nearest_distances(dm, train, &oof.inner_fold),
                    outer_nearest_distances(dm, train, &test),
                ) {
                    (Ok(di), Ok(do_)) => {
                        fold_d_inner.push(Some(di));
                        fold_d_outer.push(Some(do_));
                    }
                    (Err(e), _) | (_, Err(e)) => {
                        dist_fold_error.get_or_insert_with(|| format!("fold {f}: {e:#}"));
                        fold_d_inner.push(None);
                        fold_d_outer.push(None);
                    }
                }
            } else {
                // Even if inner fits failed on candidates/tuning, compute d_inner
                // if inner_fold itself is valid for distance-only mash_distance.
                let valid_inner = train
                    .iter_ones()
                    .all(|i| oof.inner_fold[i] < cfg.inner_folds)
                    && (0..cfg.inner_folds)
                        .filter(|&j| train.iter_ones().any(|i| oof.inner_fold[i] == j))
                        .count()
                        >= 2;
                if valid_inner {
                    match (
                        inner_nearest_distances(dm, train, &oof.inner_fold),
                        outer_nearest_distances(dm, train, &test),
                    ) {
                        (Ok(di), Ok(do_)) => {
                            fold_d_inner.push(Some(di));
                            fold_d_outer.push(Some(do_));
                        }
                        (Err(e), _) | (_, Err(e)) => {
                            dist_fold_error.get_or_insert_with(|| format!("fold {f}: {e:#}"));
                            fold_d_inner.push(None);
                            fold_d_outer.push(None);
                        }
                    }
                } else {
                    dist_fold_error.get_or_insert_with(|| {
                        oof.reason
                            .clone()
                            .unwrap_or_else(|| format!("fold {f}: inner partition invalid"))
                    });
                    fold_d_inner.push(None);
                    fold_d_outer.push(None);
                }
            }
        } else {
            fold_d_inner.push(None);
            fold_d_outer.push(None);
        }
        fold_test_pos.push(positions);
        fold_oofs.push(oof);
    }

    let mut learner_results = Vec::with_capacity(LEARNERS.len());
    for (li, &learner_name) in LEARNERS.iter().enumerate() {
        let outer_res = &learners[li];
        let mut outer_null_folds = 0usize;
        let mut calibrations = Vec::with_capacity(n_folds);
        let mut inner_score_unavailable: Option<String> = None;

        for f in 0..n_folds {
            if is_null_outer_fold(learner_name, outer_res, f) {
                outer_null_folds += 1;
            }
            let oof = &fold_oofs[f];
            let (scores_opt, null_fits, err_opt) = if let Some(ref reason) = oof.reason {
                (None, 0, Some(reason.clone()))
            } else {
                let row = &oof.learners[li];
                if let Some(ref e) = row.error {
                    (None, row.null_fits, Some(e.clone()))
                } else {
                    (Some(row.scores.as_slice()), row.null_fits, None)
                }
            };
            if let Some(ref e) = err_opt {
                inner_score_unavailable.get_or_insert_with(|| format!("fold {f}: {e}"));
            }
            calibrations.push(calibrate_fold(
                f,
                scheme,
                cfg.seed,
                y,
                &trains[f],
                scores_opt,
                null_fits,
                err_opt,
                fold_d_inner[f].as_deref(),
            ));
        }

        let distance_unavailable: Option<String> = if !distances.status.available {
            distances.status.unavailable_reason.clone()
        } else {
            dist_fold_error.clone()
        };

        let mut policies = Vec::with_capacity(POLICIES.len());
        for &policy in &POLICIES {
            let unavail_reason = match policy {
                "outcome_blind" => None,
                "score_margin" => inner_score_unavailable.clone(),
                "mash_distance" => distance_unavailable.clone(),
                "margin_plus_distance" => inner_score_unavailable
                    .clone()
                    .or_else(|| distance_unavailable.clone()),
                _ => unreachable!(),
            };
            if let Some(reason) = unavail_reason {
                policies.push(PolicySelective {
                    policy: policy.to_string(),
                    available: false,
                    unavailable_reason: Some(reason),
                    calibrated_operating_point: vec![],
                    matched_coverage: vec![],
                });
                continue;
            }

            let mut cal_points = Vec::with_capacity(COVERAGE_LADDER.len());
            let mut matched_points = Vec::with_capacity(COVERAGE_LADDER.len());

            for (ci, &target_c) in COVERAGE_LADDER.iter().enumerate() {
                let mut cal_total = CallTally::default();
                let mut cal_folds = Vec::with_capacity(n_folds);
                let mut match_total = CallTally::default();
                let mut match_folds = Vec::with_capacity(n_folds);

                for f in 0..n_folds {
                    let positions = &fold_test_pos[f];
                    let outer_null = is_null_outer_fold(learner_name, outer_res, f);
                    let cal = &calibrations[f];
                    let thr = &cal.thresholds[ci];

                    let mut cal_fold = CallTally::default();
                    let mut match_fold = CallTally::default();

                    // Compute label-blind outer confidence per test sample in E_f.
                    let mut candidates: Vec<(usize, f64)> = Vec::with_capacity(positions.len());
                    for &p in positions {
                        let u = tested[p];
                        let s = outer_res.scores[p];
                        let conf = match policy {
                            "outcome_blind" => blind_confidence(cfg.seed, scheme, f, u),
                            "score_margin" => (s - 0.5).abs(),
                            "mash_distance" => {
                                let d = fold_d_outer[f].as_ref().unwrap()[u];
                                -d
                            }
                            "margin_plus_distance" => {
                                let d = fold_d_outer[f].as_ref().unwrap()[u];
                                let jm = cal.joint_model.as_ref().unwrap();
                                jm.confidence((s - 0.5).abs(), d)
                            }
                            _ => unreachable!(),
                        };
                        candidates.push((u, conf));
                    }

                    // Layer 1: calibrated_operating_point (thresholds strictly from T_f).
                    for (local_i, &p) in positions.iter().enumerate() {
                        let u = tested[p];
                        let truth_r = y.contains(u);
                        let s = outer_res.scores[p];
                        let pred_r = s >= 0.5;
                        let (_, conf) = candidates[local_i];

                        match policy {
                            "outcome_blind" => {
                                if target_c >= 1.0 || conf >= thr.tau_blind {
                                    cal_fold.record_call(truth_r, pred_r);
                                } else {
                                    cal_fold.record_abstain(truth_r, AbstainReason::RandomThinning);
                                }
                            }
                            "score_margin" => {
                                if outer_null {
                                    cal_fold.record_abstain(truth_r, AbstainReason::NullFit);
                                } else if target_c >= 1.0 || (s - 0.5).abs() >= thr.tau_margin {
                                    cal_fold.record_call(truth_r, pred_r);
                                } else {
                                    cal_fold.record_abstain(truth_r, AbstainReason::UncertainScore);
                                }
                            }
                            "mash_distance" => {
                                let d = fold_d_outer[f].as_ref().unwrap()[u];
                                let tau_d = thr.tau_distance.unwrap();
                                if target_c >= 1.0 || d <= tau_d {
                                    cal_fold.record_call(truth_r, pred_r);
                                } else {
                                    cal_fold
                                        .record_abstain(truth_r, AbstainReason::OodMashDistance);
                                }
                            }
                            "margin_plus_distance" => {
                                if outer_null {
                                    cal_fold.record_abstain(truth_r, AbstainReason::NullFit);
                                } else {
                                    let tau_j = thr.tau_joint.unwrap();
                                    if target_c >= 1.0 || conf >= tau_j {
                                        cal_fold.record_call(truth_r, pred_r);
                                    } else {
                                        let reason = classify_joint_abstain(
                                            (s - 0.5).abs(),
                                            fold_d_outer[f].as_ref().unwrap()[u],
                                            thr.tau_margin,
                                            thr.tau_distance.unwrap(),
                                            cal.joint_model.as_ref().unwrap(),
                                        );
                                        cal_fold.record_abstain(truth_r, reason);
                                    }
                                }
                            }
                            _ => unreachable!(),
                        }
                    }

                    // Layer 2: matched_coverage (exact floor(c * n_f) per outer fold).
                    let kept = select_matched_coverage(&candidates, target_c, cfg.seed, scheme, f);
                    for (local_i, &p) in positions.iter().enumerate() {
                        let u = tested[p];
                        let truth_r = y.contains(u);
                        let s = outer_res.scores[p];
                        let pred_r = s >= 0.5;
                        if kept[local_i] {
                            match_fold.record_call(truth_r, pred_r);
                        } else {
                            let reason = match policy {
                                "outcome_blind" => AbstainReason::RandomThinning,
                                "score_margin" => {
                                    if outer_null {
                                        AbstainReason::NullFit
                                    } else {
                                        AbstainReason::UncertainScore
                                    }
                                }
                                "mash_distance" => AbstainReason::OodMashDistance,
                                "margin_plus_distance" => {
                                    if outer_null {
                                        AbstainReason::NullFit
                                    } else {
                                        classify_joint_abstain(
                                            (s - 0.5).abs(),
                                            fold_d_outer[f].as_ref().unwrap()[u],
                                            thr.tau_margin,
                                            thr.tau_distance.unwrap(),
                                            cal.joint_model.as_ref().unwrap(),
                                        )
                                    }
                                }
                                _ => unreachable!(),
                            };
                            match_fold.record_abstain(truth_r, reason);
                        }
                    }

                    cal_total.add_assign(&cal_fold);
                    cal_folds.push(cal_fold.fold_counts(f));
                    match_total.add_assign(&match_fold);
                    match_folds.push(match_fold.fold_counts(f));
                }

                cal_points.push(CoveragePoint {
                    target_coverage: target_c,
                    accounting: cal_total.finish(false),
                    folds: cal_folds,
                });
                matched_points.push(CoveragePoint {
                    target_coverage: target_c,
                    accounting: match_total.finish(false),
                    folds: match_folds,
                });
            }

            policies.push(PolicySelective {
                policy: policy.to_string(),
                available: true,
                unavailable_reason: None,
                calibrated_operating_point: cal_points,
                matched_coverage: matched_points,
            });
        }

        // Class-conditional Mondrian conformal sets calibrated strictly on s_inner[T_f].
        let (mc_avail, mc_unavail, mc_points) = if let Some(reason) = inner_score_unavailable {
            (false, Some(reason), vec![])
        } else {
            let mut pts = Vec::with_capacity(CONFORMAL_ALPHAS.len());
            for (ai, &alpha) in CONFORMAL_ALPHAS.iter().enumerate() {
                let mut total = CallTally::default();
                let mut folds_out = Vec::with_capacity(n_folds);
                for f in 0..n_folds {
                    let positions = &fold_test_pos[f];
                    let outer_null = is_null_outer_fold(learner_name, outer_res, f);
                    let ct = &calibrations[f].conformal[ai];
                    let mut fold_tally = CallTally::default();
                    for &p in positions {
                        let u = tested[p];
                        let truth_r = y.contains(u);
                        if outer_null {
                            fold_tally.record_abstain(truth_r, AbstainReason::NullFit);
                        } else {
                            let s = outer_res.scores[p];
                            let has_r = ct.q_r.is_none_or(|q| (1.0 - s) <= q);
                            let has_s = ct.q_s.is_none_or(|q| s <= q);
                            fold_tally.record_conformal(truth_r, has_r, has_s);
                        }
                    }
                    total.add_assign(&fold_tally);
                    folds_out.push(fold_tally.fold_counts(f));
                }
                pts.push(ConformalOperatingPoint {
                    alpha,
                    accounting: total.finish(true),
                    folds: folds_out,
                });
            }
            (true, None, pts)
        };

        learner_results.push(LearnerSelective {
            learner: learner_name.to_string(),
            outer_null_folds,
            calibrations,
            policies,
            mondrian_conformal_available: mc_avail,
            mondrian_conformal_unavailable_reason: mc_unavail,
            mondrian_conformal: mc_points,
        });
    }

    base.learners = learner_results;
    base.seconds = started.elapsed().as_secs_f64();
    Ok(base)
}

/// Attributes an abstention by `margin_plus_distance` to genomic novelty
/// (`OodMashDistance`) vs score uncertainty (`UncertainScore`) using only the
/// `T_f`-calibrated thresholds and joint model weights.
fn classify_joint_abstain(
    margin: f64,
    distance: f64,
    tau_margin: f64,
    tau_distance: f64,
    jm: &JointModelParams,
) -> AbstainReason {
    let ood = distance > tau_distance;
    let low_margin = margin < tau_margin;
    if ood && !low_margin {
        AbstainReason::OodMashDistance
    } else if low_margin && !ood {
        AbstainReason::UncertainScore
    } else {
        // Both or neither individually breached: attribute to whichever
        // standardized risk term contributes more upward risk.
        let dist_risk = jm.weight_distance * jm.z_distance(distance);
        let margin_risk = -jm.weight_margin * jm.z_margin(margin);
        if dist_risk > margin_risk {
            AbstainReason::OodMashDistance
        } else {
            AbstainReason::UncertainScore
        }
    }
}

// ------------------------------------------------ hooks for analysis.rs

/// Hook for drawn splits (`random` and `lineage`): writes one
/// `grm-selective-v1` file per split (`<dataset>.<scheme>-s<seed>.selective.json`).
#[allow(clippy::too_many_arguments)]
pub fn record_drawn(
    layout: &Layout,
    dataset: &str,
    species_slug: &str,
    antibiotic: &str,
    sample_ids: &[String],
    x: &BitMatrix,
    y: &Bits,
    groups: &[u32],
    duplicates: &[Vec<usize>],
    min_support: u32,
    cfg: &EvalConfig,
    splits: &[SplitResult],
) -> Result<()> {
    let n = x.cols;
    let active = Bits::ones(n);
    let distances = load_cohort_distances(layout, species_slug, sample_ids, &active)?;
    let all: Vec<usize> = (0..n).collect();
    for split in splits {
        let trains: Vec<Bits> = if split.unavailable_reason.is_none() {
            fold_masks(&split.test_fold, split.folds)
                .into_iter()
                .map(|(train, _)| train)
                .collect()
        } else {
            vec![]
        };
        let result = evaluate_split(
            dataset,
            species_slug,
            antibiotic,
            &split.scheme,
            x,
            y,
            groups,
            duplicates,
            &trains,
            &all,
            &split.test_fold,
            &split.learners,
            split.unavailable_reason.clone(),
            None,
            min_support,
            cfg,
            &distances,
        )?;
        write_result(
            &file_name(dataset, &split.scheme, cfg.seed),
            serde_json::to_string(&result)?,
        )?;
    }
    Ok(())
}

/// Hook for `--fold-file` supplied splits (`study-whole`, `time`, etc.): writes
/// `<dataset>.<scheme>-s<seed>.selective.json`.
#[allow(clippy::too_many_arguments)]
pub fn record_assigned(
    layout: &Layout,
    dataset: &str,
    species_slug: &str,
    antibiotic: &str,
    scheme: &str,
    fold_file_sha256: &str,
    sample_ids: &[String],
    x: &BitMatrix,
    y: &Bits,
    groups: &[u32],
    duplicates: &[Vec<usize>],
    roles: &[Role],
    folds: usize,
    tested: &[usize],
    test_fold: &[usize],
    learners: &[LearnerResult],
    min_support: u32,
    cfg: &EvalConfig,
) -> Result<()> {
    let n = x.cols;
    let used = |i: usize| roles[i] != Role::Unused;
    let active = Bits::from_indices(n, (0..n).filter(|&i| used(i)));
    let distances = load_cohort_distances(layout, species_slug, sample_ids, &active)?;

    let used_duplicates: Vec<Vec<usize>> = duplicates
        .iter()
        .map(|group| {
            group
                .iter()
                .copied()
                .filter(|&i| used(i))
                .collect::<Vec<_>>()
        })
        .filter(|group| group.len() > 1)
        .collect();

    let mut trains = Vec::with_capacity(folds);
    for f in 0..folds {
        let train = Bits::from_indices(
            n,
            (0..n).filter(|&i| match roles[i] {
                Role::Test(g) => g != f,
                Role::Train => true,
                Role::Unused => false,
            }),
        );
        trains.push(train);
    }

    let result = evaluate_split(
        dataset,
        species_slug,
        antibiotic,
        scheme,
        x,
        y,
        groups,
        &used_duplicates,
        &trains,
        tested,
        test_fold,
        learners,
        None,
        Some(fold_file_sha256.to_string()),
        min_support,
        cfg,
        &distances,
    )?;
    write_result(
        &file_name(dataset, scheme, cfg.seed),
        serde_json::to_string(&result)?,
    )
}

/// Checks whether a split's selective output file already exists in the
/// enabled sink directory (if any).
pub fn outputs_exist_drawn(dataset: &str, seed: u64) -> bool {
    let Some(dir) = SINK.get() else {
        return false;
    };
    ["random", "lineage"]
        .iter()
        .all(|scheme| dir.join(file_name(dataset, scheme, seed)).exists())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::species_check::GenomeSketch;

    #[test]
    fn retention_count_is_exact_floor_without_float_underflow() {
        // 0.95 * 140 = 133.0 mathematically, though 0.95_f64 * 140.0 = 132.99999999999997.
        assert_eq!(retention_count(0.95, 140), 133);
        assert_eq!(retention_count(0.85, 20), 17);
        assert_eq!(retention_count(0.85, 25), 21);
        assert_eq!(retention_count(0.70, 30), 21);
        assert_eq!(retention_count(1.00, 37), 37);
        assert_eq!(retention_count(0.50, 37), 18);
        assert_eq!(retention_count(0.0, 37), 0);
        assert_eq!(retention_count(0.85, 0), 0);
    }

    #[test]
    fn sketches_bin_round_trips_and_loads_with_cache_fallback() {
        let dir = tempfile::tempdir().unwrap();
        let layout = Layout::new(dir.path());
        let ids: Vec<String> = (0..4).map(|i| format!("9999.{i}")).collect();
        let sketches = [
            vec![1u64, 2, 3, 4],
            vec![1u64, 2, 3, 5],
            vec![1u64, 2, 6, 7],
            vec![10u64, 20, 30, 40],
        ];
        // Write 3 sketches to sketches.bin and 1 to species-check cache.
        let bin_path = layout.species_dir("sp").join("sketches.bin");
        let bin_entries: Vec<(String, Vec<u64>)> = ids[..3]
            .iter()
            .cloned()
            .zip(sketches[..3].iter().cloned())
            .collect();
        write_sketches_bin(
            &bin_path,
            species_check::K,
            species_check::SKETCH,
            &bin_entries,
        )
        .unwrap();
        let (bk, bs, map) = read_sketches_bin(&bin_path).unwrap();
        assert_eq!(
            (bk, bs, map.len()),
            (species_check::K, species_check::SKETCH, 3)
        );

        // Without the 4th sketch in cache, distances report unavailable without leaking any ID.
        let all = Bits::ones(4);
        let missing = load_cohort_distances(&layout, "sp", &ids, &all).unwrap();
        assert!(!missing.status.available);
        let reason = missing.status.unavailable_reason.unwrap();
        assert!(reason.contains("1 of 4"), "{reason}");
        for id in &ids {
            assert!(!reason.contains(id), "genome id leaked in reason: {reason}");
        }

        // Adding the 4th sketch to species-check cache completes the cohort cleanly.
        let cache_path = species_check::sketch_path(&layout, &ids[3]);
        std::fs::create_dir_all(cache_path.parent().unwrap()).unwrap();
        let gs = GenomeSketch {
            genome_id: ids[3].clone(),
            k: species_check::K,
            sketch_size: species_check::SKETCH,
            length: 1000,
            gc: 0.5,
            sketch: sketches[3].clone(),
        };
        std::fs::write(&cache_path, serde_json::to_vec(&gs).unwrap()).unwrap();

        let loaded = load_cohort_distances(&layout, "sp", &ids, &all).unwrap();
        assert!(loaded.status.available);
        assert_eq!(
            loaded.status.source.as_deref(),
            Some("species_sketches_bin+species_check_cache")
        );
        let dm = loaded.matrix.unwrap();
        assert_eq!(dm.get(0, 0), 0.0);
        assert_eq!(dm.get(0, 1), dm.get(1, 0));
        assert!(dm.get(0, 1) < dm.get(0, 3));
    }

    #[test]
    fn fold_safe_distances_never_read_outer_test_or_same_inner_fold() {
        // 6 samples: 0..4 in train (inner folds [0, 0, 1, 1]), 4..6 in test.
        let mut dm = DistanceMatrix::from_dense(
            6,
            21,
            1000,
            vec![
                // 0     1     2     3     4     5
                0.00, 0.01, 0.20, 0.30, 0.001, 0.001, // 0 (inner 0)
                0.01, 0.00, 0.25, 0.15, 0.001, 0.001, // 1 (inner 0)
                0.20, 0.25, 0.00, 0.02, 0.001, 0.001, // 2 (inner 1)
                0.30, 0.15, 0.02, 0.00, 0.001, 0.001, // 3 (inner 1)
                0.001, 0.001, 0.001, 0.001, 0.00, 0.0001, // 4 (test)
                0.001, 0.001, 0.001, 0.001, 0.0001, 0.00, // 5 (test)
            ],
        )
        .unwrap();
        let train = Bits::from_indices(6, 0..4);
        let test = Bits::from_indices(6, 4..6);
        let inner_fold = [0, 0, 1, 1, usize::MAX, usize::MAX];

        let d_in = inner_nearest_distances(&dm, &train, &inner_fold).unwrap();
        // Sample 0 (inner 0): nearest in inner 1 is sample 2 (0.20), ignoring sample 1 (0.01) and test (0.001).
        assert!((d_in[0] - 0.20).abs() < 1e-6);
        // Sample 1 (inner 0): nearest in inner 1 is sample 3 (0.15).
        assert!((d_in[1] - 0.15).abs() < 1e-6);
        // Sample 2 (inner 1): nearest in inner 0 is sample 0 (0.20).
        assert!((d_in[2] - 0.20).abs() < 1e-6);
        // Sample 3 (inner 1): nearest in inner 0 is sample 1 (0.15).
        assert!((d_in[3] - 0.15).abs() < 1e-6);
        assert!(d_in[4].is_nan() && d_in[5].is_nan());

        // Overwrite test rows/cols with NaN; inner_nearest_distances must still succeed bit-for-bit!
        for u in 4..6 {
            for j in 0..6 {
                dm.set(u, j, f32::NAN);
                dm.set(j, u, f32::NAN);
            }
        }
        let d_in2 = inner_nearest_distances(&dm, &train, &inner_fold).unwrap();
        for i in 0..4 {
            assert_eq!(d_in[i].to_bits(), d_in2[i].to_bits());
        }

        // For d_outer: test-to-test entry (4,5) = NaN must never be read; only test-to-train.
        dm.set(4, 0, 0.08);
        dm.set(4, 1, 0.12);
        dm.set(4, 2, 0.05);
        dm.set(4, 3, 0.09);
        dm.set(5, 0, 0.40);
        dm.set(5, 1, 0.35);
        dm.set(5, 2, 0.30);
        dm.set(5, 3, 0.25);
        let d_out = outer_nearest_distances(&dm, &train, &test).unwrap();
        assert!((d_out[4] - 0.05).abs() < 1e-6);
        assert!((d_out[5] - 0.25).abs() < 1e-6);
        assert!(d_out[0..4].iter().all(|v| v.is_nan()));
    }

    #[test]
    fn joint_risk_model_is_monotone_and_handles_constant_margin_scm() {
        // Synthetic inner fold: errors occur when margin is small or distance is large.
        let margins = [0.48, 0.45, 0.42, 0.40, 0.05, 0.08, 0.10, 0.45, 0.40, 0.04];
        let dists = [0.01, 0.02, 0.01, 0.03, 0.02, 0.03, 0.09, 0.15, 0.14, 0.16];
        let errors = [
            false, false, false, false, true, true, true, true, true, true,
        ];
        let jm = fit_joint_risk_model(&errors, &margins, &dists, JOINT_L2_LAMBDA);
        assert!(jm.weight_margin > 0.0, "{jm:?}");
        assert!(jm.weight_distance > 0.0, "{jm:?}");
        assert!(!jm.fallback_equal_weights);
        // Higher margin -> higher confidence; higher distance -> lower confidence.
        assert!(jm.confidence(0.45, 0.02) > jm.confidence(0.05, 0.02));
        assert!(jm.confidence(0.45, 0.02) > jm.confidence(0.45, 0.15));

        // Constant margin (binary SCM): sd_margin == 0 => weight_margin == 0, distance still works!
        let scm_margins = [0.5; 10];
        let jm_scm = fit_joint_risk_model(&errors, &scm_margins, &dists, JOINT_L2_LAMBDA);
        assert_eq!(jm_scm.weight_margin, 0.0);
        assert!(jm_scm.weight_distance > 0.0, "{jm_scm:?}");
        assert!(jm_scm.confidence(0.5, 0.01) > jm_scm.confidence(0.5, 0.15));

        // Zero inner errors: triggers deterministic equal-weight fallback without NaN.
        let no_errors = [false; 10];
        let jm_zero = fit_joint_risk_model(&no_errors, &margins, &dists, JOINT_L2_LAMBDA);
        assert!(jm_zero.fallback_equal_weights);
        assert_eq!(jm_zero.weight_margin, 1.0);
        assert_eq!(jm_zero.weight_distance, 1.0);
        assert!(jm_zero.confidence(0.45, 0.01).is_finite());
    }

    #[test]
    fn hand_computed_accounting_and_conformal_both_vs_empty() {
        // Test Case A: conformal_both ({R, S}) vs conformal_empty ({}) separated cleanly.
        // q_r = 0.25 (has_r iff 1-s <= 0.25 <=> s >= 0.75)
        // q_s = 0.30 (has_s iff s <= 0.30)
        // So s in (0.30, 0.75) -> empty {} (conformal_empty, NOT covering true label!).
        let mut t_empty = CallTally::default();
        let samples = [
            (true, 0.90),  // {R} -> TP, covered
            (true, 0.80),  // {R} -> TP, covered
            (true, 0.20),  // {S} -> FN, NOT covered
            (true, 0.50),  // {}  -> conformal_empty (R), NOT covered
            (false, 0.10), // {S} -> TN, covered
            (false, 0.25), // {S} -> TN, covered
            (false, 0.85), // {R} -> FP, NOT covered
            (false, 0.55), // {}  -> conformal_empty (S), NOT covered
        ];
        for (truth_r, s) in samples {
            t_empty.record_conformal(truth_r, (1.0 - s) <= 0.25, s <= 0.30);
        }
        let acc = t_empty.finish(true);
        assert_eq!((acc.n, acc.n_r, acc.n_s), (8, 4, 4));
        assert_eq!((acc.called, acc.called_r, acc.called_s), (6, 3, 3));
        assert_eq!((acc.abstain_r, acc.abstain_s), (1, 1));
        assert_eq!((acc.tp, acc.fp, acc.tn, acc.fn_), (2, 1, 2, 1));
        assert_eq!(acc.abstain_reasons.conformal_empty, 2);
        assert_eq!(acc.abstain_reasons.conformal_both, 0);
        // True label covered: 4 of 8 (2 TP + 2 TN; empty sets do NOT cover true label).
        let cov = acc.set_true_label_coverage.as_ref().unwrap();
        assert_eq!((cov.x, cov.n), (4, 8));
        assert!((cov.rate.unwrap() - 0.5).abs() < 1e-12);
        assert!((acc.coverage.unwrap() - 0.75).abs() < 1e-12);
        assert!((acc.coverage_r.unwrap() - 0.75).abs() < 1e-12);
        assert!((acc.coverage_s.unwrap() - 0.75).abs() < 1e-12);
        assert!((acc.vme_called.rate.unwrap() - 1.0 / 3.0).abs() < 1e-12);
        assert!((acc.vme_all.rate.unwrap() - 1.0 / 4.0).abs() < 1e-12);
        assert!((acc.me_called.rate.unwrap() - 1.0 / 3.0).abs() < 1e-12);
        assert!((acc.me_all.rate.unwrap() - 1.0 / 4.0).abs() < 1e-12);
        assert!((acc.npv.rate.unwrap() - 2.0 / 3.0).abs() < 1e-12);
        assert!((acc.for_s.rate.unwrap() - 1.0 / 3.0).abs() < 1e-12);
        assert!((acc.ppv.rate.unwrap() - 2.0 / 3.0).abs() < 1e-12);
        assert!(acc.sparse_denominator);

        // Test Case B: q_r = 0.70 (s >= 0.30), q_s = 0.65 (s <= 0.65)
        // So s in [0.30, 0.65] -> {R, S} (conformal_both, COVERS true label!).
        let mut t_both = CallTally::default();
        for (truth_r, s) in samples {
            t_both.record_conformal(truth_r, (1.0 - s) <= 0.70, s <= 0.65);
        }
        let acc_b = t_both.finish(true);
        assert_eq!(acc_b.abstain_reasons.conformal_both, 2);
        assert_eq!(acc_b.abstain_reasons.conformal_empty, 0);
        // True label covered: 2 TP + 2 TN + 2 {R, S} = 6 of 8!
        let cov_b = acc_b.set_true_label_coverage.as_ref().unwrap();
        assert_eq!((cov_b.x, cov_b.n), (6, 8));
        assert!((cov_b.rate.unwrap() - 0.75).abs() < 1e-12);
    }

    #[test]
    fn capture_records_drawn_and_assigned_selective_json_files() {
        let dir = tempfile::tempdir().unwrap();
        let layout = Layout::new(dir.path());
        let n = 24;
        let y = Bits::from_indices(n, (0..n).filter(|i| i % 2 == 0));
        let groups: Vec<u32> = (0..n as u32).map(|i| i / 2).collect();
        let mut x = BitMatrix::with_cols(n);
        x.push_row(Bits::from_indices(n, (0..n).filter(|i| i % 2 == 0)).words());
        let sample_ids: Vec<String> = (0..n).map(|i| format!("562.{}", 70000 + i)).collect();
        let cfg = EvalConfig {
            folds: 2,
            inner_folds: 2,
            ..Default::default()
        };
        let unavail_split = SplitResult {
            scheme: "lineage".into(),
            folds: 2,
            n_groups: 12,
            test_fold: vec![],
            learners: vec![],
            unavailable_reason: Some("synthetic unavailable".into()),
        };
        let (res, files) = capture(|| {
            assert!(enabled());
            record_drawn(
                &layout,
                "sp__drug",
                "sp",
                "drug",
                &sample_ids,
                &x,
                &y,
                &groups,
                &[],
                3,
                &cfg,
                &[unavail_split],
            )
        });
        res.unwrap();
        assert!(!enabled());
        assert_eq!(files.len(), 1);
        assert_eq!(files[0].0, "sp__drug.lineage-s42.selective.json");
        let parsed: SelectiveResult = serde_json::from_str(&files[0].1).unwrap();
        assert_eq!(parsed.format, FORMAT);
        assert_eq!(
            parsed.unavailable_reason.as_deref(),
            Some("synthetic unavailable")
        );
        assert!(parsed.learners.is_empty());
        for id in &sample_ids {
            assert!(!files[0].1.contains(id));
        }
    }
}
