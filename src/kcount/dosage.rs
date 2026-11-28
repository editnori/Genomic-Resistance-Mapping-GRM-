//! Modality-verified core-panel depth normalization (`D_core`), target/plasmid
//! dosage extraction, and training-fold core-panel construction.
//!
//! What this module checks in code:
//! - Core-panel schema version (`1`), `k` range (`5..=64`), non-empty unique
//!   canonical k-mers, and mandatory canonical `panel_sha256` computed over
//!   schema version, `k`, `input_type`, and sorted canonical k-mers.
//! - Stored `.ktab` header integrity (SHA-256) and stored `InputModality`
//!   (`reads` vs `assembly`) directly against `panel.input_type`, rejecting
//!   `InputModality::Unknown` and any cross-modality mismatch even when a CLI
//!   `--input-type` flag agrees with the panel. If `--input-type` is supplied,
//!   both the flag and the table header modality must match `panel.input_type`.
//! - Optional `--expected-fold-id` equality against `panel.fold_id`.
//! - Exact `k` agreement across table, core panel, and optional target panel.
//! - Core-coverage cutoffs (`d_core >= min_core_depth`, default `10.0`, and
//!   `core_breadth >= min_core_breadth`, default `0.80`), emitting
//!   `dosage_status: "insufficient_core_coverage"` and JSON `null` (`None`)
//!   for every normalized dosage value when coverage is insufficient.
//! - Optional QC vector (`DosageQcVector`) where every field is `Option<...>`,
//!   serializing as JSON `null` when `--qc-json` is omitted or a field is
//!   absent.
//!
//! What the caller must supply:
//! - When constructing a training-fold core panel with
//!   [`build_training_fold_core_panel`] (`grm kcount build-core-panel`), the
//!   caller is responsible for excluding held-out test genomes from the
//!   training input list.

use crate::kcount::format::{
    AnyTableReader, DEFAULT_SAMPLE_NAME, InputModality, KtabReader, hex_bytes,
};
use crate::kcount::key::{KmerKey, KmerSpace64, KmerSpace128, MAX_K, MAX_K_U64, MIN_K};
use crate::kset::{BUCKET_BITS, BUCKETS};
use crate::matrix::{REVIEWED_COLUMNS, ensure_reviewed};
use anyhow::{Context, Result, bail, ensure};
use ring::digest::{Context as Sha256Context, SHA256};
use rustc_hash::FxHashMap;
use serde::{Deserialize, Serialize};
use std::cmp::Reverse;
use std::collections::{BTreeMap, BinaryHeap};
use std::fs;
use std::path::{Path, PathBuf};

const CORE_PANEL_SCHEMA_VERSION: u32 = 1;
const CORE_PANEL_DIGEST_DOMAIN: &[u8] = b"GRM-CORE-PANEL-V1\n";

/// Input sequencing modality recorded on a core panel and matched against the
/// `.ktab` header's stored [`InputModality`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum InputType {
    Reads,
    Assembly,
}

impl InputType {
    #[inline]
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Reads => "reads",
            Self::Assembly => "assembly",
        }
    }

    #[inline]
    pub fn to_modality(self) -> InputModality {
        match self {
            Self::Reads => InputModality::Reads,
            Self::Assembly => InputModality::Assembly,
        }
    }

    pub fn from_modality(modality: InputModality) -> Result<Self> {
        match modality {
            InputModality::Reads => Ok(Self::Reads),
            InputModality::Assembly => Ok(Self::Assembly),
            InputModality::Unknown => {
                bail!("table header records unknown input modality; expected reads or assembly")
            }
        }
    }
}

/// Computes the canonical SHA-256 digest of a core panel over its domain tag,
/// `schema_version`, `k`, `input_type`, and sorted canonical k-mers.
pub fn compute_panel_sha256(
    schema_version: u32,
    k: usize,
    input_type: InputType,
    sorted_canonical_kmers: &[String],
) -> String {
    let mut ctx = Sha256Context::new(&SHA256);
    ctx.update(CORE_PANEL_DIGEST_DOMAIN);
    ctx.update(format!("schema_version={schema_version}\n").as_bytes());
    ctx.update(format!("k={k}\n").as_bytes());
    ctx.update(format!("input_type={}\n", input_type.as_str()).as_bytes());
    for s in sorted_canonical_kmers {
        ctx.update(s.as_bytes());
        ctx.update(b"\n");
    }
    hex_bytes(ctx.finish().as_ref())
}

/// Single-copy core k-mer panel with mandatory canonical `panel_sha256`.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct CorePanel {
    pub schema_version: u32,
    pub k: usize,
    pub input_type: InputType,
    pub fold_id: Option<String>,
    pub panel_sha256: String,
    pub kmers: Vec<String>,
}

impl CorePanel {
    /// Builds a validated `CorePanel` with canonicalized, sorted, deduplicated
    /// k-mers and a mandatory canonical SHA-256 digest covering schema
    /// version, `k`, `input_type`, and sorted k-mers.
    pub fn new(
        k: usize,
        input_type: InputType,
        fold_id: Option<String>,
        kmers: Vec<String>,
    ) -> Result<Self> {
        ensure!((MIN_K..=MAX_K).contains(&k), "unsupported k {k}");
        ensure!(
            !kmers.is_empty(),
            "core panel must contain at least 1 k-mer"
        );
        let canonical_kmers = canonicalize_panel_strings(k, &kmers)?;
        let panel_sha256 =
            compute_panel_sha256(CORE_PANEL_SCHEMA_VERSION, k, input_type, &canonical_kmers);
        Ok(Self {
            schema_version: CORE_PANEL_SCHEMA_VERSION,
            k,
            input_type,
            fold_id,
            panel_sha256,
            kmers: canonical_kmers,
        })
    }

    /// Loads and validates a core panel JSON file, rejecting missing, empty,
    /// or mismatched `panel_sha256`.
    pub fn load(path: &Path) -> Result<Self> {
        let bytes =
            fs::read(path).with_context(|| format!("read core panel {}", path.display()))?;
        let panel: Self = serde_json::from_slice(&bytes)
            .with_context(|| format!("parse core panel {}", path.display()))?;
        panel.validate()?;
        Ok(panel)
    }

    /// Validates schema version, `k`, non-empty unique canonical k-mers, and
    /// mandatory `panel_sha256`.
    pub fn validate(&self) -> Result<()> {
        ensure!(
            self.schema_version == CORE_PANEL_SCHEMA_VERSION,
            "unsupported core panel schema_version {}",
            self.schema_version
        );
        ensure!(
            (MIN_K..=MAX_K).contains(&self.k),
            "unsupported core panel k {}",
            self.k
        );
        ensure!(
            !self.kmers.is_empty(),
            "core panel must contain at least 1 k-mer"
        );
        ensure!(
            !self.panel_sha256.trim().is_empty(),
            "core panel requires a non-empty panel_sha256"
        );
        let canonical = canonicalize_panel_strings(self.k, &self.kmers)?;
        ensure!(
            canonical.len() == self.kmers.len(),
            "core panel contains duplicate canonical k-mers"
        );
        let actual_sha =
            compute_panel_sha256(self.schema_version, self.k, self.input_type, &canonical);
        ensure!(
            actual_sha == self.panel_sha256,
            "core panel sha256 mismatch (expected {}, got {actual_sha})",
            self.panel_sha256
        );
        Ok(())
    }
}

/// Named sets of target (plasmid or AMR) k-mers for normalized dosage scoring.
#[derive(Clone, Debug, Default, Serialize, Deserialize, PartialEq, Eq)]
pub struct TargetPanel {
    pub k: usize,
    pub targets: BTreeMap<String, Vec<String>>,
}

impl TargetPanel {
    pub fn load(path: &Path) -> Result<Self> {
        let bytes =
            fs::read(path).with_context(|| format!("read target panel {}", path.display()))?;
        let panel: Self = serde_json::from_slice(&bytes)
            .with_context(|| format!("parse target panel {}", path.display()))?;
        ensure!(
            (MIN_K..=MAX_K).contains(&panel.k),
            "unsupported target panel k {}",
            panel.k
        );
        Ok(panel)
    }
}

/// Dosage callability status; when `InsufficientCoreCoverage`, every
/// normalized dosage field serializes as JSON `null` (`None`).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum DosageStatus {
    Ok,
    InsufficientCoreCoverage,
}

impl DosageStatus {
    #[inline]
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Ok => "ok",
            Self::InsufficientCoreCoverage => "insufficient_core_coverage",
        }
    }
}

/// Configuration thresholds and verification options for
/// [`extract_sample_dosage`] and [`evaluate_table_core_coverage`].
#[derive(Clone, Debug)]
pub struct DosageConfig {
    pub sample_name: String,
    /// Optional explicit input modality check in addition to the `.ktab`
    /// header's stored modality. When set, both this value and the `.ktab`
    /// header's stored modality must match `panel.input_type`.
    pub expected_input_type: Option<InputType>,
    /// Optional fold identifier verified against `panel.fold_id`.
    pub expected_fold_id: Option<String>,
    pub min_count: u32,
    pub min_core_depth: f64,
    pub min_core_breadth: f64,
}

impl Default for DosageConfig {
    fn default() -> Self {
        Self {
            sample_name: DEFAULT_SAMPLE_NAME.to_string(),
            expected_input_type: None,
            expected_fold_id: None,
            min_count: 1,
            min_core_depth: 10.0,
            min_core_breadth: 0.80,
        }
    }
}

/// Optional read/contig QC metrics parsed independently from `--qc-json`.
///
/// Every field is `Option<...>` so omitted `--qc-json` or absent fields
/// serialize as JSON `null`, never `0` or `0.0`.
#[derive(Clone, Debug, Default, Serialize, Deserialize, PartialEq)]
pub struct DosageQcVector {
    #[serde(default)]
    pub total_reads: Option<u64>,
    #[serde(default)]
    pub total_bases: Option<u64>,
    #[serde(default)]
    pub valid_bases: Option<u64>,
    #[serde(default)]
    pub masked_bases: Option<u64>,
    #[serde(default)]
    pub gc_fraction: Option<f64>,
    #[serde(default)]
    pub usable_read_len_median: Option<f64>,
    #[serde(default)]
    pub usable_read_len_mean: Option<f64>,
    #[serde(default)]
    pub usable_read_len_p10: Option<f64>,
    #[serde(default)]
    pub usable_read_len_p90: Option<f64>,
    #[serde(default)]
    pub n50: Option<u64>,
    #[serde(default)]
    pub singleton_fraction: Option<f64>,
}

impl DosageQcVector {
    fn validate_finite(&self) -> Result<()> {
        for (name, val) in [
            ("gc_fraction", self.gc_fraction),
            ("usable_read_len_median", self.usable_read_len_median),
            ("usable_read_len_mean", self.usable_read_len_mean),
            ("usable_read_len_p10", self.usable_read_len_p10),
            ("usable_read_len_p90", self.usable_read_len_p90),
            ("singleton_fraction", self.singleton_fraction),
        ] {
            if let Some(v) = val {
                ensure!(v.is_finite(), "QC field {name} must be finite");
            }
        }
        Ok(())
    }
}

/// Parses optional QC metrics from either a `CountSummary`-shaped JSON file
/// (`{"qc": { ... }}`) or a standalone `ReadQcMetrics`-shaped JSON object.
/// Present fields decode to `Some(...)`; absent fields remain `None` (JSON
/// `null`).
pub fn load_qc_metrics(path: &Path) -> Result<DosageQcVector> {
    let bytes = fs::read(path).with_context(|| format!("read qc json {}", path.display()))?;
    let val: serde_json::Value = serde_json::from_slice(&bytes)
        .with_context(|| format!("parse qc json {}", path.display()))?;
    let obj = val
        .as_object()
        .with_context(|| format!("qc json {} must be a JSON object", path.display()))?;
    let target = if let Some(inner) = obj.get("qc") {
        ensure!(
            inner.is_object(),
            "qc json {} has non-object \"qc\" field",
            path.display()
        );
        inner.clone()
    } else {
        val
    };
    let qc: DosageQcVector = serde_json::from_value(target)
        .with_context(|| format!("decode qc fields from {}", path.display()))?;
    qc.validate_finite()?;
    Ok(qc)
}

/// Per-k-mer raw count and normalized dosage entry.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct KmerDosageEntry {
    pub kmer: String,
    pub raw_count: u32,
    pub dosage: Option<f64>,
}

/// Per-target (e.g. plasmid replicon or AMR locus) dosage summary.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct TargetDosageSummary {
    pub target: String,
    pub n_kmers: usize,
    pub detected_kmers: usize,
    pub breadth: f64,
    pub median_raw_count: f64,
    pub mean_raw_count: f64,
    /// Median raw count across all target k-mers divided by `d_core`; `None`
    /// (JSON `null`) whenever `dosage_status == InsufficientCoreCoverage`.
    pub normalized_dosage: Option<f64>,
    /// Median raw count across detected target k-mers divided by `d_core`;
    /// `None` (JSON `null`) whenever `dosage_status == InsufficientCoreCoverage`
    /// or no target k-mers are detected.
    pub detected_normalized_dosage: Option<f64>,
    pub kmers: Vec<KmerDosageEntry>,
}

/// Complete per-sample dosage and optional QC report. Contains no input paths.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct SampleDosageReport {
    pub sample: String,
    pub k: usize,
    pub input_type: InputType,
    pub table_modality: InputModality,
    pub fold_id: Option<String>,
    pub core_panel_sha256: String,
    pub dosage_status: DosageStatus,
    pub d_core: f64,
    pub d_core_all_median: f64,
    pub core_breadth: f64,
    pub core_panel_size: usize,
    pub core_detected: usize,
    pub min_core_depth: f64,
    pub min_core_breadth: f64,
    pub qc: DosageQcVector,
    pub targets: Vec<TargetDosageSummary>,
}

fn validate_dosage_config(cfg: &DosageConfig) -> Result<()> {
    ensure!(cfg.min_count >= 1, "min_count must be >= 1");
    ensure!(
        cfg.min_core_depth.is_finite() && cfg.min_core_depth > 0.0,
        "min_core_depth must be finite and > 0.0"
    );
    ensure!(
        cfg.min_core_breadth.is_finite() && (0.0..=1.0).contains(&cfg.min_core_breadth),
        "min_core_breadth must be finite and within [0.0, 1.0]"
    );
    Ok(())
}

fn verify_modality_and_fold(
    table_modality: Option<InputModality>,
    core_panel: &CorePanel,
    cfg: &DosageConfig,
) -> Result<InputModality> {
    validate_dosage_config(cfg)?;
    if let Some(expected_fold) = &cfg.expected_fold_id {
        ensure!(
            core_panel.fold_id.as_deref() == Some(expected_fold.as_str()),
            "core panel fold_id mismatch: panel has {:?}, expected {:?}",
            core_panel.fold_id,
            expected_fold
        );
    }
    if let Some(flag_type) = cfg.expected_input_type {
        ensure!(
            flag_type == core_panel.input_type,
            "input_type flag mismatch: flag is {}, core panel requires {}",
            flag_type.as_str(),
            core_panel.input_type.as_str()
        );
    }
    let stored = table_modality.context(
        "input table does not record an InputModality header; dosage requires a verified .ktab table with reads or assembly modality",
    )?;
    let stored_type = InputType::from_modality(stored).with_context(|| {
        format!(
            "table header modality must be {} to match core panel",
            core_panel.input_type.as_str()
        )
    })?;
    ensure!(
        stored_type == core_panel.input_type,
        "table modality mismatch: table header records {}, core panel requires {}",
        stored.as_str(),
        core_panel.input_type.as_str()
    );
    Ok(stored)
}

/// Extracts modality-verified `D_core` and optional target/plasmid normalized
/// dosages from a verified counted `.ktab` table.
pub fn extract_sample_dosage(
    table_path: &Path,
    core_panel: &CorePanel,
    targets: Option<&TargetPanel>,
    read_qc: Option<&DosageQcVector>,
    cfg: &DosageConfig,
) -> Result<SampleDosageReport> {
    core_panel.validate()?;
    if let Some(q) = read_qc {
        q.validate_finite()?;
    }
    let reader = KtabReader::open_verified(table_path)?;
    ensure!(
        reader.k == core_panel.k,
        "k mismatch: table has k={}, core panel has k={}",
        reader.k,
        core_panel.k
    );
    let table_modality = verify_modality_and_fold(Some(reader.modality), core_panel, cfg)?;
    if let Some(tp) = targets {
        ensure!(
            tp.k == reader.k,
            "k mismatch: table has k={}, target panel has k={}",
            reader.k,
            tp.k
        );
    }

    let PanelLookup {
        core_counts,
        target_results,
    } = if reader.k <= MAX_K_U64 {
        lookup_panel_counts::<u64>(&reader, core_panel, targets)?
    } else {
        lookup_panel_counts::<u128>(&reader, core_panel, targets)?
    };

    let detected_core: Vec<u32> = core_counts
        .iter()
        .copied()
        .filter(|&c| c >= cfg.min_count)
        .collect();
    let core_detected = detected_core.len();
    let core_panel_size = core_counts.len();
    let core_breadth = core_detected as f64 / core_panel_size as f64;
    let d_core = median_u32(&detected_core);
    let d_core_all_median = median_u32(&core_counts);

    let dosage_status =
        if d_core >= cfg.min_core_depth && core_breadth >= cfg.min_core_breadth && d_core > 0.0 {
            DosageStatus::Ok
        } else {
            DosageStatus::InsufficientCoreCoverage
        };

    let mut target_summaries = Vec::new();
    for (tname, kmer_counts) in target_results {
        let n_kmers = kmer_counts.len();
        let raw_vals: Vec<u32> = kmer_counts.iter().map(|(_, c)| *c).collect();
        let det_vals: Vec<u32> = raw_vals
            .iter()
            .copied()
            .filter(|&c| c >= cfg.min_count)
            .collect();
        let detected_kmers = det_vals.len();
        let breadth = if n_kmers == 0 {
            0.0
        } else {
            detected_kmers as f64 / n_kmers as f64
        };
        let median_raw_count = median_u32(&raw_vals);
        let mean_raw_count = if n_kmers == 0 {
            0.0
        } else {
            raw_vals.iter().map(|&c| c as f64).sum::<f64>() / n_kmers as f64
        };
        let normalized_dosage = if dosage_status == DosageStatus::Ok {
            Some(median_raw_count / d_core)
        } else {
            None
        };
        let detected_normalized_dosage =
            if dosage_status == DosageStatus::Ok && !det_vals.is_empty() {
                Some(median_u32(&det_vals) / d_core)
            } else {
                None
            };
        let kmers = kmer_counts
            .into_iter()
            .map(|(kmer, raw_count)| KmerDosageEntry {
                kmer,
                raw_count,
                dosage: if dosage_status == DosageStatus::Ok {
                    Some(raw_count as f64 / d_core)
                } else {
                    None
                },
            })
            .collect();
        target_summaries.push(TargetDosageSummary {
            target: tname,
            n_kmers,
            detected_kmers,
            breadth,
            median_raw_count,
            mean_raw_count,
            normalized_dosage,
            detected_normalized_dosage,
            kmers,
        });
    }

    let qc = read_qc.cloned().unwrap_or_default();

    Ok(SampleDosageReport {
        sample: cfg.sample_name.clone(),
        k: reader.k,
        input_type: core_panel.input_type,
        table_modality,
        fold_id: core_panel.fold_id.clone(),
        core_panel_sha256: core_panel.panel_sha256.clone(),
        dosage_status,
        d_core,
        d_core_all_median,
        core_breadth,
        core_panel_size,
        core_detected,
        min_core_depth: cfg.min_core_depth,
        min_core_breadth: cfg.min_core_breadth,
        qc,
        targets: target_summaries,
    })
}

/// Evaluates modality-verified `D_core`, `core_breadth`, and [`DosageStatus`]
/// for a counted `.ktab` reader (used by `grm kcount matrix-export --mode
/// dosage`).
pub fn evaluate_table_core_coverage(
    reader: &AnyTableReader,
    core_panel: &CorePanel,
    cfg: &DosageConfig,
) -> Result<(DosageStatus, f64, f64)> {
    core_panel.validate()?;
    let tab = match reader {
        AnyTableReader::Tab(t) => t,
        _ => bail!(
            "dosage mode requires counted .ktab inputs (got {})",
            reader.format().magic_label()
        ),
    };
    ensure!(
        tab.k == core_panel.k,
        "k mismatch: table has k={}, core panel has k={}",
        tab.k,
        core_panel.k
    );
    verify_modality_and_fold(Some(tab.modality), core_panel, cfg)?;

    let counts = if tab.k <= MAX_K_U64 {
        lookup_core_counts::<u64>(tab, core_panel)?
    } else {
        lookup_core_counts::<u128>(tab, core_panel)?
    };
    let detected: Vec<u32> = counts
        .iter()
        .copied()
        .filter(|&c| c >= cfg.min_count)
        .collect();
    let breadth = detected.len() as f64 / counts.len() as f64;
    let d_core = median_u32(&detected);
    let status = if d_core >= cfg.min_core_depth && breadth >= cfg.min_core_breadth && d_core > 0.0
    {
        DosageStatus::Ok
    } else {
        DosageStatus::InsufficientCoreCoverage
    };
    Ok((status, d_core, breadth))
}

struct PanelLookup {
    core_counts: Vec<u32>,
    target_results: Vec<(String, Vec<(String, u32)>)>,
}

fn lookup_panel_counts<K: KmerKey>(
    reader: &KtabReader,
    core_panel: &CorePanel,
    targets: Option<&TargetPanel>,
) -> Result<PanelLookup> {
    let k = reader.k;
    let space = K::new_space(k);

    let mut bucket_queries: BTreeMap<usize, Vec<K>> = BTreeMap::new();
    let mut core_keys = Vec::with_capacity(core_panel.kmers.len());
    for s in &core_panel.kmers {
        let canon = K::canonical(&space, s)
            .with_context(|| format!("invalid k={k} core panel k-mer: {s}"))?;
        let mixed = K::mix(&space, canon);
        core_keys.push(mixed);
        bucket_queries
            .entry(mixed.bucket_of(k, BUCKET_BITS))
            .or_default()
            .push(mixed);
    }

    let mut target_keys: Vec<(String, Vec<(String, K)>)> = Vec::new();
    if let Some(tp) = targets {
        for (tname, kmers) in &tp.targets {
            let mut list = Vec::with_capacity(kmers.len());
            for s in kmers {
                let canon = K::canonical(&space, s)
                    .with_context(|| format!("invalid k={k} target k-mer in {tname}: {s}"))?;
                let mixed = K::mix(&space, canon);
                list.push((s.clone(), mixed));
                bucket_queries
                    .entry(mixed.bucket_of(k, BUCKET_BITS))
                    .or_default()
                    .push(mixed);
            }
            target_keys.push((tname.clone(), list));
        }
    }

    let mut found = FxHashMap::<u128, u32>::default();
    let mut kbuf = Vec::new();
    let mut cbuf = Vec::new();
    let mut keys = Vec::new();
    let mut counts = Vec::new();
    for (b, qkeys) in bucket_queries {
        keys.clear();
        counts.clear();
        reader.read_bucket(b, &mut kbuf, &mut cbuf, &mut keys, &mut counts)?;
        for q in qkeys {
            if let Ok(idx) = keys.binary_search(&q) {
                found.insert(q.to_u128(), counts[idx]);
            }
        }
    }

    let core_counts: Vec<u32> = core_keys
        .into_iter()
        .map(|key| *found.get(&key.to_u128()).unwrap_or(&0))
        .collect();
    let mut target_out = Vec::with_capacity(target_keys.len());
    for (tname, list) in target_keys {
        let entries = list
            .into_iter()
            .map(|(seq, key)| (seq, *found.get(&key.to_u128()).unwrap_or(&0)))
            .collect();
        target_out.push((tname, entries));
    }
    Ok(PanelLookup {
        core_counts,
        target_results: target_out,
    })
}

fn lookup_core_counts<K: KmerKey>(reader: &KtabReader, core_panel: &CorePanel) -> Result<Vec<u32>> {
    let k = reader.k;
    let space = K::new_space(k);
    let mut bucket_queries: BTreeMap<usize, Vec<K>> = BTreeMap::new();
    let mut core_keys = Vec::with_capacity(core_panel.kmers.len());
    for s in &core_panel.kmers {
        let canon = K::canonical(&space, s)
            .with_context(|| format!("invalid k={k} core panel k-mer: {s}"))?;
        let mixed = K::mix(&space, canon);
        core_keys.push(mixed);
        bucket_queries
            .entry(mixed.bucket_of(k, BUCKET_BITS))
            .or_default()
            .push(mixed);
    }

    let mut found = FxHashMap::<u128, u32>::default();
    let mut kbuf = Vec::new();
    let mut cbuf = Vec::new();
    let mut keys = Vec::new();
    let mut counts = Vec::new();
    for (b, qkeys) in bucket_queries {
        keys.clear();
        counts.clear();
        reader.read_bucket(b, &mut kbuf, &mut cbuf, &mut keys, &mut counts)?;
        for q in qkeys {
            if let Ok(idx) = keys.binary_search(&q) {
                found.insert(q.to_u128(), counts[idx]);
            }
        }
    }
    Ok(core_keys
        .into_iter()
        .map(|key| *found.get(&key.to_u128()).unwrap_or(&0))
        .collect())
}

fn canonicalize_panel_strings(k: usize, kmers: &[String]) -> Result<Vec<String>> {
    let mut out = Vec::with_capacity(kmers.len());
    if k <= MAX_K_U64 {
        let sp = KmerSpace64::new(k);
        for s in kmers {
            let c = sp
                .canonical(s)
                .with_context(|| format!("invalid k={k} DNA k-mer in core panel: {s}"))?;
            out.push(sp.decode(c));
        }
    } else {
        let sp = KmerSpace128::new(k);
        for s in kmers {
            let c = sp
                .canonical(s)
                .with_context(|| format!("invalid k={k} DNA k-mer in core panel: {s}"))?;
            out.push(sp.decode(c));
        }
    }
    out.sort();
    out.dedup();
    Ok(out)
}

pub(crate) fn median_u32(vals: &[u32]) -> f64 {
    if vals.is_empty() {
        return 0.0;
    }
    let mut sorted = vals.to_vec();
    sorted.sort_unstable();
    let n = sorted.len();
    if n % 2 == 1 {
        sorted[n / 2] as f64
    } else {
        (sorted[n / 2 - 1] as f64 + sorted[n / 2] as f64) / 2.0
    }
}

/// Computes the exact median k-mer count of one table from a sorted count
/// histogram without retaining all k-mer counts in memory.
fn table_median_count<K: KmerKey>(reader: &AnyTableReader) -> Result<f64> {
    let total = reader.n_kmers();
    if total == 0 {
        return Ok(0.0);
    }
    match reader {
        AnyTableReader::V1(_) | AnyTableReader::Gr(_) => return Ok(1.0),
        AnyTableReader::Tab(t) if t.presence_only => return Ok(1.0),
        AnyTableReader::Tab(_) => {}
    }
    let mut hist: BTreeMap<u32, u64> = BTreeMap::new();
    let mut kbuf = Vec::new();
    let mut cbuf = Vec::new();
    let mut keys: Vec<K> = Vec::new();
    let mut counts: Vec<u32> = Vec::new();
    for b in 0..BUCKETS {
        keys.clear();
        counts.clear();
        reader.read_bucket(b, &mut kbuf, &mut cbuf, &mut keys, &mut counts)?;
        for &c in &counts {
            *hist.entry(c).or_insert(0) += 1;
        }
    }
    Ok(median_from_histogram(&hist, total))
}

fn median_from_histogram(hist: &BTreeMap<u32, u64>, total: u64) -> f64 {
    if total == 0 {
        return 0.0;
    }
    let lo_rank = (total - 1) / 2;
    let hi_rank = total / 2;
    let mut seen = 0u64;
    let mut lo_val = None;
    let mut hi_val = None;
    for (&val, &freq) in hist {
        let next = seen + freq;
        if lo_val.is_none() && next > lo_rank {
            lo_val = Some(val as f64);
        }
        if hi_val.is_none() && next > hi_rank {
            hi_val = Some(val as f64);
            break;
        }
        seen = next;
    }
    match (lo_val, hi_val) {
        (Some(a), Some(b)) => (a + b) / 2.0,
        _ => 0.0,
    }
}

/// Configuration for [`build_training_fold_core_panel`].
#[derive(Clone, Debug)]
pub struct CorePanelBuildConfig {
    pub k: usize,
    pub input_type: InputType,
    pub min_prevalence: f64,
    pub max_copy_factor: f64,
    pub max_kmers: usize,
    pub fold_id: Option<String>,
}

/// Builds a training-fold single-copy core k-mer panel from 1 to
/// `REVIEWED_COLUMNS` verified training tables.
///
/// A k-mer is retained when:
/// 1. Its prevalence across the supplied training tables is at least
///    `cfg.min_prevalence`.
/// 2. In every training table carrying that k-mer, its count is at most
///    `sample_median_depth * cfg.max_copy_factor` (excluding multicopy/mobile
///    or plasmid k-mers).
///
/// The caller is responsible for passing only training-fold tables and
/// excluding held-out test genomes.
pub fn build_training_fold_core_panel(
    tables: &[PathBuf],
    cfg: &CorePanelBuildConfig,
) -> Result<CorePanel> {
    ensure!(!tables.is_empty(), "at least one training table required");
    ensure_reviewed(tables.len(), REVIEWED_COLUMNS)?;
    ensure!((MIN_K..=MAX_K).contains(&cfg.k), "unsupported k {}", cfg.k);
    ensure!(
        cfg.min_prevalence.is_finite() && cfg.min_prevalence > 0.0 && cfg.min_prevalence <= 1.0,
        "min_prevalence must be in (0.0, 1.0]"
    );
    ensure!(
        cfg.max_copy_factor.is_finite() && cfg.max_copy_factor >= 1.0,
        "max_copy_factor must be finite and >= 1.0"
    );
    ensure!(cfg.max_kmers >= 1, "max_kmers must be >= 1");

    let want_modality = cfg.input_type.to_modality();
    let mut readers = Vec::with_capacity(tables.len());
    for (idx, p) in tables.iter().enumerate() {
        let r = AnyTableReader::open(p)?;
        r.verify_integrity()?;
        ensure!(
            r.k() == cfg.k,
            "k mismatch: training table #{idx} has k={}, expected {}",
            r.k(),
            cfg.k
        );
        let stored = r.modality().with_context(|| {
            format!(
                "training table #{idx} has no InputModality header; expected {}",
                cfg.input_type.as_str()
            )
        })?;
        ensure!(
            stored == want_modality,
            "table modality mismatch: training table #{idx} records {}, expected {}",
            stored.as_str(),
            want_modality.as_str()
        );
        readers.push(r);
    }

    if cfg.k <= MAX_K_U64 {
        build_panel_generic::<u64>(&readers, cfg)
    } else {
        build_panel_generic::<u128>(&readers, cfg)
    }
}

fn build_panel_generic<K: KmerKey>(
    readers: &[AnyTableReader],
    cfg: &CorePanelBuildConfig,
) -> Result<CorePanel> {
    let n_tables = readers.len();
    let space = K::new_space(cfg.k);

    let mut max_allowed = Vec::with_capacity(n_tables);
    for r in readers {
        let med = table_median_count::<K>(r)?;
        max_allowed.push(med * cfg.max_copy_factor);
    }

    let mut selected = Vec::new();
    let mut kbuf = Vec::new();
    let mut cbuf = Vec::new();
    let mut per_keys: Vec<Vec<K>> = (0..n_tables).map(|_| Vec::new()).collect();
    let mut per_counts: Vec<Vec<u32>> = (0..n_tables).map(|_| Vec::new()).collect();
    let mut cursors = vec![0usize; n_tables];
    let mut heap: BinaryHeap<Reverse<(K, usize)>> = BinaryHeap::with_capacity(n_tables);

    for b in 0..BUCKETS {
        if selected.len() >= cfg.max_kmers {
            break;
        }
        heap.clear();
        for i in 0..n_tables {
            per_keys[i].clear();
            per_counts[i].clear();
            readers[i].read_bucket(
                b,
                &mut kbuf,
                &mut cbuf,
                &mut per_keys[i],
                &mut per_counts[i],
            )?;
            cursors[i] = 0;
            if let Some(&first_key) = per_keys[i].first() {
                heap.push(Reverse((first_key, i)));
            }
        }

        while let Some(Reverse((key, mut idx))) = heap.pop() {
            let mut support = 0usize;
            let mut exceeded_copy_cap = false;

            loop {
                let pos = cursors[idx];
                let c = per_counts[idx][pos];
                cursors[idx] = pos + 1;
                if let Some(&next_key) = per_keys[idx].get(pos + 1) {
                    heap.push(Reverse((next_key, idx)));
                }
                support += 1;
                if (c as f64) > max_allowed[idx] {
                    exceeded_copy_cap = true;
                }

                match heap.peek() {
                    Some(&Reverse((next_k, _))) if next_k == key => {
                        idx = heap.pop().unwrap().0.1;
                    }
                    _ => break,
                }
            }

            let prevalence = support as f64 / n_tables as f64;
            if prevalence >= cfg.min_prevalence && !exceeded_copy_cap {
                selected.push(K::decode(&space, K::unmix(&space, key)));
                if selected.len() >= cfg.max_kmers {
                    break;
                }
            }
        }
    }

    CorePanel::new(cfg.k, cfg.input_type, cfg.fold_id.clone(), selected)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn median_u32_and_histogram_are_exact_for_even_odd_and_empty() {
        assert_eq!(median_u32(&[]), 0.0);
        assert_eq!(median_u32(&[42]), 42.0);
        assert_eq!(median_u32(&[40, 10, 20]), 20.0);
        assert_eq!(median_u32(&[40, 10, 30, 20]), 25.0);
        assert_eq!(median_u32(&[1, 2, 3, 4]), 2.5);

        let mut hist = BTreeMap::new();
        hist.insert(10, 1);
        hist.insert(20, 1);
        hist.insert(30, 1);
        hist.insert(40, 1);
        assert_eq!(median_from_histogram(&hist, 4), 25.0);

        hist.insert(50, 1);
        assert_eq!(median_from_histogram(&hist, 5), 30.0);
    }

    #[test]
    fn core_panel_requires_valid_mandatory_sha256_over_metadata_and_kmers() {
        let dir = tempfile::tempdir().unwrap();
        let panel = CorePanel::new(
            21,
            InputType::Reads,
            Some("fold-0".to_string()),
            vec![
                "ACGTACGTACGTACGTACGTA".to_string(),
                "TGCATGCATGCATGCATGCAT".to_string(),
            ],
        )
        .unwrap();
        assert_eq!(panel.panel_sha256.len(), 64);
        panel.validate().unwrap();

        // Missing panel_sha256 in JSON fails deserialization.
        let mut val = serde_json::to_value(&panel).unwrap();
        val.as_object_mut().unwrap().remove("panel_sha256");
        let missing_path = dir.path().join("missing_sha.json");
        fs::write(&missing_path, serde_json::to_vec(&val).unwrap()).unwrap();
        assert!(CorePanel::load(&missing_path).is_err());

        // Empty panel_sha256 fails validate().
        let mut empty_panel = panel.clone();
        empty_panel.panel_sha256 = String::new();
        assert!(empty_panel.validate().is_err());

        // Tampered panel_sha256 fails validate().
        let mut bad_sha = panel.clone();
        bad_sha.panel_sha256 = "0".repeat(64);
        let err = bad_sha.validate().unwrap_err();
        assert!(format!("{err:#}").contains("sha256 mismatch"));

        // Editing input_type without recomputing panel_sha256 fails validate().
        let mut swapped_mod = panel;
        swapped_mod.input_type = InputType::Assembly;
        let err = swapped_mod.validate().unwrap_err();
        assert!(format!("{err:#}").contains("sha256 mismatch"));
    }

    #[test]
    fn qc_vector_serializes_all_null_when_omitted_and_parses_partial_or_wrapped_json() {
        let dir = tempfile::tempdir().unwrap();

        // Default DosageQcVector must serialize EVERY field as JSON null, zero 0 / 0.0 values.
        let empty = DosageQcVector::default();
        let v = serde_json::to_value(&empty).unwrap();
        let obj = v.as_object().unwrap();
        assert_eq!(obj.len(), 11);
        for (field, fv) in obj {
            assert!(
                fv.is_null(),
                "omitted QC field {field} serialized as {fv}, expected null"
            );
        }

        // Wrapped CountSummary-shaped JSON.
        let wrapped = dir.path().join("wrapped.json");
        fs::write(
            &wrapped,
            r#"{
                "sample": "ignored",
                "k": 31,
                "qc": {
                    "total_reads": 500,
                    "gc_fraction": 0.52,
                    "usable_read_len_median": 150.0,
                    "singleton_fraction": 0.12,
                    "is_fastq": true
                }
            }"#,
        )
        .unwrap();
        let parsed = load_qc_metrics(&wrapped).unwrap();
        assert_eq!(parsed.total_reads, Some(500));
        assert_eq!(parsed.gc_fraction, Some(0.52));
        assert_eq!(parsed.usable_read_len_median, Some(150.0));
        assert_eq!(parsed.singleton_fraction, Some(0.12));
        assert_eq!(parsed.total_bases, None);
        assert_eq!(parsed.n50, None);
        let pv = serde_json::to_value(&parsed).unwrap();
        assert!(pv["total_bases"].is_null());
        assert!(pv["n50"].is_null());

        // Flat ReadQcMetrics-shaped JSON.
        let flat = dir.path().join("flat.json");
        fs::write(&flat, r#"{"total_reads": 42, "n50": 12000}"#).unwrap();
        let parsed_flat = load_qc_metrics(&flat).unwrap();
        assert_eq!(parsed_flat.total_reads, Some(42));
        assert_eq!(parsed_flat.n50, Some(12_000));
        assert_eq!(parsed_flat.gc_fraction, None);
    }
}
