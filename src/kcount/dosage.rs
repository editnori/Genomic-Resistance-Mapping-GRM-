//! Fold-safe `D_core` depth normalization and plasmid/target dosage extractor.
//!
//! Enforces statistical anti-leakage guardrails:
//! - The single-copy core k-mer panel must match the exact `k` and input type
//!   (`reads` vs `assembly`) of the counted sample table.
//! - If `d_core < min_core_depth` (default `10.0`) or
//!   `core_breadth < min_core_breadth` (default `0.80`), the sample is marked
//!   `insufficient_core_coverage` and all normalized dosage outputs are `None`
//!   (JSON `null`), never silently zero-filled.
//! - Emits a per-sample QC vector alongside dosage summaries so downstream
//!   models can run QC-only control checks.

use crate::kcount::format::{AnyTableReader, KtabReader};
use crate::kcount::key::{KmerKey, KmerSpace64, KmerSpace128, MAX_K, MAX_K_U64, MIN_K};
use crate::kcount::superkmer::{CountSummary, ReadQcMetrics, hex_bytes};
use crate::kset::{BUCKET_BITS, BUCKETS};
use anyhow::{Context, Result, ensure};
use ring::digest::{Context as Sha256Context, SHA256};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::fs;
use std::path::Path;

/// Input sequencing modality enforced to match between core panel and sample.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum InputType {
    Reads,
    Assembly,
}

/// Pinned or training-fold single-copy core k-mer panel.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct CorePanel {
    pub schema_version: u32,
    pub k: usize,
    pub input_type: InputType,
    pub fold_id: Option<String>,
    pub panel_sha256: Option<String>,
    pub kmers: Vec<String>,
}

impl CorePanel {
    /// Builds a validated `CorePanel` with a deterministic SHA-256 over its
    /// sorted canonical k-mer strings.
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
        let sha = panel_digest(&canonical_kmers);
        Ok(Self {
            schema_version: 1,
            k,
            input_type,
            fold_id,
            panel_sha256: Some(sha),
            kmers: canonical_kmers,
        })
    }

    pub fn load(path: &Path) -> Result<Self> {
        let bytes =
            fs::read(path).with_context(|| format!("read core panel {}", path.display()))?;
        let panel: Self = serde_json::from_slice(&bytes)
            .with_context(|| format!("parse core panel {}", path.display()))?;
        panel.validate()?;
        Ok(panel)
    }

    pub fn validate(&self) -> Result<()> {
        ensure!(
            self.schema_version == 1,
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
            "core panel must contain at least one k-mer"
        );
        let canonical = canonicalize_panel_strings(self.k, &self.kmers)?;
        if let Some(expected_sha) = &self.panel_sha256 {
            let actual_sha = panel_digest(&canonical);
            ensure!(
                actual_sha == *expected_sha,
                "core panel sha256 mismatch (expected {expected_sha}, got {actual_sha})"
            );
        }
        Ok(())
    }
}

/// Named sets of target (plasmid or AMR) k-mers for normalized dosage scoring.
#[derive(Clone, Debug, Default, Serialize, Deserialize, PartialEq)]
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

/// Dosage callability status; low-depth samples NEVER silently zero-fill.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum DosageStatus {
    Ok,
    InsufficientCoreCoverage,
}

/// Configuration thresholds for [`extract_sample_dosage`].
#[derive(Clone, Copy, Debug)]
pub struct DosageConfig {
    pub input_type: InputType,
    pub min_count: u32,
    pub min_core_depth: f64,
    pub min_core_breadth: f64,
}

impl Default for DosageConfig {
    fn default() -> Self {
        Self {
            input_type: InputType::Reads,
            min_count: 1,
            min_core_depth: 10.0,
            min_core_breadth: 0.80,
        }
    }
}

/// Per-sample QC vector emitted alongside dosage calls for QC-only controls.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct SampleQcVector {
    pub sample: String,
    pub usable_read_len_median: f64,
    pub gc_fraction: f64,
    pub singleton_fraction: f64,
    pub d_core: f64,
    pub core_breadth: f64,
    pub total_reads: u64,
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
    /// whenever `dosage_status == InsufficientCoreCoverage`.
    pub normalized_dosage: Option<f64>,
    /// Median raw count across detected target k-mers divided by `d_core`;
    /// `None` whenever `dosage_status == InsufficientCoreCoverage` or no
    /// target k-mers are detected.
    pub detected_normalized_dosage: Option<f64>,
    pub kmers: Vec<KmerDosageEntry>,
}

/// Complete per-sample dosage and QC report.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct SampleDosageReport {
    pub sample: String,
    pub k: usize,
    pub input_type: InputType,
    pub fold_id: Option<String>,
    pub core_panel_sha256: Option<String>,
    pub dosage_status: DosageStatus,
    pub d_core: f64,
    pub d_core_all_median: f64,
    pub core_breadth: f64,
    pub core_panel_size: usize,
    pub core_detected: usize,
    pub min_core_depth: f64,
    pub min_core_breadth: f64,
    pub qc: SampleQcVector,
    pub targets: Vec<TargetDosageSummary>,
}

/// Extracts fold-safe `D_core` and optional target/plasmid normalized dosages
/// from a counted `.ktab` table.
pub fn extract_sample_dosage(
    table_path: &Path,
    core_panel: &CorePanel,
    targets: Option<&TargetPanel>,
    read_qc: Option<&ReadQcMetrics>,
    cfg: &DosageConfig,
) -> Result<SampleDosageReport> {
    core_panel.validate()?;
    let reader = KtabReader::open(table_path)?;
    ensure!(
        reader.k == core_panel.k,
        "k mismatch: table has k={}, core panel has k={}",
        reader.k,
        core_panel.k
    );
    ensure!(
        cfg.input_type == core_panel.input_type,
        "input_type mismatch: sample is {:?}, core panel is {:?}",
        cfg.input_type,
        core_panel.input_type
    );
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

    let qc = SampleQcVector {
        sample: reader.name.clone(),
        usable_read_len_median: read_qc.map_or(0.0, |q| q.usable_read_len_median),
        gc_fraction: read_qc.map_or(0.0, |q| q.gc_fraction),
        singleton_fraction: read_qc.map_or(0.0, |q| q.singleton_fraction),
        d_core,
        core_breadth,
        total_reads: read_qc.map_or(0, |q| q.total_reads),
    };

    Ok(SampleDosageReport {
        sample: reader.name,
        k: reader.k,
        input_type: cfg.input_type,
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

/// Evaluates `D_core`, `core_breadth`, and [`DosageStatus`] for any table
/// reader (used by `op matrix-export --mode dosage`).
pub fn evaluate_table_core_coverage(
    reader: &AnyTableReader,
    core_panel: &CorePanel,
    cfg: &DosageConfig,
) -> Result<(DosageStatus, f64, f64)> {
    core_panel.validate()?;
    ensure!(
        reader.k() == core_panel.k,
        "k mismatch: table '{}' has k={}, core panel has k={}",
        reader.name(),
        reader.k(),
        core_panel.k
    );
    ensure!(
        cfg.input_type == core_panel.input_type,
        "input_type mismatch: requested {:?}, core panel is {:?}",
        cfg.input_type,
        core_panel.input_type
    );
    let counts = if reader.k() <= MAX_K_U64 {
        lookup_core_any::<u64>(reader, core_panel)?
    } else {
        lookup_core_any::<u128>(reader, core_panel)?
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

/// Loads [`ReadQcMetrics`] from a metadata JSON file (either a raw
/// [`ReadQcMetrics`] or a full [`CountSummary`]).
pub fn load_qc_metrics(path: &Path) -> Result<ReadQcMetrics> {
    let bytes = fs::read(path).with_context(|| format!("read qc json {}", path.display()))?;
    if let Ok(summary) = serde_json::from_slice::<CountSummary>(&bytes) {
        return Ok(summary.qc);
    }
    serde_json::from_slice::<ReadQcMetrics>(&bytes)
        .with_context(|| format!("parse qc json {}", path.display()))
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

    // Map every query mixed key to its bucket so we read each touched bucket once.
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

    let mut found = rustc_hash::FxHashMap::<u128, u32>::default();
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

fn lookup_core_any<K: KmerKey>(
    reader: &AnyTableReader,
    core_panel: &CorePanel,
) -> Result<Vec<u32>> {
    let k = reader.k();
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

    let mut found = rustc_hash::FxHashMap::<u128, u32>::default();
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

fn panel_digest(sorted_kmers: &[String]) -> String {
    let mut ctx = Sha256Context::new(&SHA256);
    for s in sorted_kmers {
        ctx.update(s.as_bytes());
        ctx.update(b"\n");
    }
    hex_bytes(ctx.finish().as_ref())
}

fn median_u32(vals: &[u32]) -> f64 {
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

/// Builds a training-fold single-copy core k-mer panel from a set of
/// training-fold tables (present in at least `min_prevalence` fraction of
/// training samples with per-sample count <= `max_copy_factor * sample_median`).
pub fn build_training_fold_core_panel(
    tables: &[&Path],
    k: usize,
    input_type: InputType,
    min_prevalence: f64,
    max_kmers: usize,
    fold_id: Option<String>,
) -> Result<CorePanel> {
    ensure!(!tables.is_empty(), "at least one training table required");
    ensure!(
        (0.0..=1.0).contains(&min_prevalence) && min_prevalence > 0.0,
        "min_prevalence must be in (0, 1]"
    );
    let readers: Vec<AnyTableReader> = tables
        .iter()
        .map(|p| AnyTableReader::open(p))
        .collect::<Result<Vec<_>>>()?;
    for r in &readers {
        ensure!(
            r.k() == k,
            "training table '{}' has k={}, expected {k}",
            r.name(),
            r.k()
        );
    }
    let min_support = ((readers.len() as f64) * min_prevalence).ceil() as usize;

    let mut selected = Vec::new();
    let mut kbuf = Vec::new();
    let mut cbuf = Vec::new();
    let sp64 = (k <= MAX_K_U64).then(|| KmerSpace64::new(k));
    let sp128 = (k > MAX_K_U64).then(|| KmerSpace128::new(k));

    for b in 0..BUCKETS {
        if selected.len() >= max_kmers {
            break;
        }
        let mut counts_across = rustc_hash::FxHashMap::<u128, usize>::default();
        for r in &readers {
            let mut keys = Vec::new();
            let mut cnts = Vec::new();
            r.read_bucket::<u128>(b, &mut kbuf, &mut cbuf, &mut keys, &mut cnts)?;
            for key in keys {
                *counts_across.entry(key).or_insert(0) += 1;
            }
        }
        let mut bucket_hits: Vec<u128> = counts_across
            .into_iter()
            .filter(|&(_, sup)| sup >= min_support)
            .map(|(key, _)| key)
            .collect();
        bucket_hits.sort_unstable();
        for key in bucket_hits {
            if selected.len() >= max_kmers {
                break;
            }
            let seq = if let Some(sp) = sp64 {
                sp.decode(sp.unmix(key as u64))
            } else {
                let sp = sp128.unwrap();
                sp.decode(sp.unmix(key))
            };
            selected.push(seq);
        }
    }
    CorePanel::new(k, input_type, fold_id, selected)
}
