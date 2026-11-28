//! Streaming bucket-by-bucket set algebra (`union`, `intersect`, `diff`,
//! `symdiff`) and ML-ready `matrix-export` (`pa`, `counts`, `dosage`) over
//! any mix of `.kset` v1 (`GRMKSET1`), counted `.ktab` (`GRMKTB1\0`), and
//! Golomb-Rice sorted-key (`GRMKGR1\0`) tables.

use crate::fnv::Fnv;
use crate::kcount::dosage::{CorePanel, DosageConfig, DosageStatus, evaluate_table_core_coverage};
use crate::kcount::format::{
    AnyTableReader, encode_bucket_counts, encode_bucket_keys, write_gr_kset_from_buckets,
    write_ktab_from_buckets,
};
use crate::kcount::key::{KmerKey, MAX_K_U64};
use crate::kcount::superkmer::{hex_bytes, write_json_atomic};
use crate::kset::BUCKETS;
use anyhow::{Context, Result, ensure};
use ring::digest::{Context as Sha256Context, SHA256};
use serde::{Deserialize, Serialize};
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

/// Set algebra operation applied across 2 or more k-mer tables.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SetOp {
    /// Keys present in at least one input; counts are summed (saturating `u32`).
    Union,
    /// Keys present in every input; counts take the minimum across inputs.
    Intersect,
    /// Keys present in input 0 and absent from all subsequent inputs.
    Diff,
    /// Keys present in an odd number of inputs (exact symmetric difference for 2 inputs).
    SymDiff,
}

/// Summary returned by [`run_set_op`].
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct SetOpSummary {
    pub op: SetOp,
    pub k: usize,
    pub n_inputs: usize,
    pub n_kmers_out: u64,
    pub max_count_out: u32,
    pub presence_only: bool,
    pub fnv64: String,
    pub keys_sha256: String,
    pub output_path: String,
}

/// Feature value encoding mode for [`export_matrix`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum MatrixMode {
    Pa,
    Counts,
    Dosage,
}

/// Output file layout for [`export_matrix`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum MatrixFormat {
    Sparse,
    Dense,
}

/// Configuration for [`export_matrix`].
#[derive(Clone, Debug)]
pub struct MatrixExportConfig {
    pub mode: MatrixMode,
    pub format: MatrixFormat,
    pub min_support: usize,
    pub max_support: Option<usize>,
    pub core_panel: Option<CorePanel>,
    pub dosage: DosageConfig,
}

/// Per-sample missingness metadata recorded in `missing_mask.tsv` and `manifest.json`.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct SampleMaskEntry {
    pub sample_idx: usize,
    pub sample: String,
    pub missing: bool,
    pub dosage_status: Option<DosageStatus>,
    pub d_core: Option<f64>,
    pub core_breadth: Option<f64>,
}

/// Deterministic manifest written to `<out_dir>/manifest.json`.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct MatrixManifest {
    pub schema_version: u32,
    pub k: usize,
    pub mode: MatrixMode,
    pub format: MatrixFormat,
    pub n_samples: usize,
    pub n_features: u64,
    pub n_missing_samples: usize,
    pub samples: Vec<SampleMaskEntry>,
    pub sha256_features: String,
    pub sha256_matrix: String,
    pub sha256_missing_mask: String,
}

/// Runs a streaming bucket-by-bucket merge-join set operation (`union`,
/// `intersect`, `diff`, or `symdiff`) across 2 or more tables.
pub fn run_set_op(
    op: SetOp,
    inputs: &[PathBuf],
    out_path: &Path,
    name: &str,
    presence_only_out: bool,
) -> Result<SetOpSummary> {
    ensure!(
        inputs.len() >= 2,
        "set operation requires at least 2 inputs"
    );
    let readers: Vec<AnyTableReader> = inputs
        .iter()
        .map(|p| AnyTableReader::open(p))
        .collect::<Result<Vec<_>>>()?;
    let k = readers[0].k();
    for r in &readers[1..] {
        ensure!(
            r.k() == k,
            "k mismatch across inputs: '{}' has k={}, expected k={k}",
            r.name(),
            r.k()
        );
    }
    if k <= MAX_K_U64 {
        run_set_op_generic::<u64>(op, &readers, k, out_path, name, presence_only_out)
    } else {
        run_set_op_generic::<u128>(op, &readers, k, out_path, name, presence_only_out)
    }
}

fn run_set_op_generic<K: KmerKey>(
    op: SetOp,
    readers: &[AnyTableReader],
    k: usize,
    out_path: &Path,
    name: &str,
    presence_only_out: bool,
) -> Result<SetOpSummary> {
    let n_inputs = readers.len();
    let mut bucket_counts = [0u32; BUCKETS];
    let mut key_payloads = vec![Vec::new(); BUCKETS];
    let mut count_payloads = vec![Vec::new(); BUCKETS];

    let mut fnv = Fnv::new();
    let mut keys_sha = Sha256Context::new(&SHA256);
    let mut n_kmers_out = 0u64;
    let mut max_count_out = 0u32;
    let mut all_ones = true;

    let mut kbuf = Vec::new();
    let mut cbuf = Vec::new();
    let mut per_input_keys: Vec<Vec<K>> = (0..n_inputs).map(|_| Vec::new()).collect();
    let mut per_input_counts: Vec<Vec<u32>> = (0..n_inputs).map(|_| Vec::new()).collect();
    let mut cursors = vec![0usize; n_inputs];
    let mut out_keys: Vec<K> = Vec::new();
    let mut out_counts: Vec<u32> = Vec::new();

    for b in 0..BUCKETS {
        for i in 0..n_inputs {
            per_input_keys[i].clear();
            per_input_counts[i].clear();
            readers[i].read_bucket(
                b,
                &mut kbuf,
                &mut cbuf,
                &mut per_input_keys[i],
                &mut per_input_counts[i],
            )?;
            cursors[i] = 0;
        }
        out_keys.clear();
        out_counts.clear();

        // Streaming k-way merge-join over sorted keys in bucket b.
        loop {
            let mut min_key: Option<K> = None;
            for i in 0..n_inputs {
                if let Some(&key) = per_input_keys[i].get(cursors[i]) {
                    min_key = Some(match min_key {
                        None => key,
                        Some(cur) => cur.min(key),
                    });
                }
            }
            let Some(key) = min_key else {
                break;
            };

            let mut present_count = 0usize;
            let mut in_first = false;
            let mut first_count = 0u32;
            let mut sum_count = 0u32;
            let mut min_count = u32::MAX;
            for i in 0..n_inputs {
                if per_input_keys[i].get(cursors[i]) == Some(&key) {
                    let c = per_input_counts[i][cursors[i]];
                    cursors[i] += 1;
                    present_count += 1;
                    if i == 0 {
                        in_first = true;
                        first_count = c;
                    }
                    sum_count = sum_count.saturating_add(c);
                    min_count = min_count.min(c);
                }
            }

            let emitted: Option<u32> = match op {
                SetOp::Union => Some(sum_count),
                SetOp::Intersect if present_count == n_inputs => Some(min_count),
                SetOp::Diff if in_first && present_count == 1 => Some(first_count),
                SetOp::SymDiff if present_count % 2 == 1 => Some(sum_count),
                _ => None,
            };

            if let Some(c) = emitted {
                let c_final = if presence_only_out { 1 } else { c };
                out_keys.push(key);
                out_counts.push(c_final);
                key.with_le_bytes(|b| {
                    fnv.bytes(b);
                    keys_sha.update(b);
                });
                n_kmers_out += 1;
                if c_final > max_count_out {
                    max_count_out = c_final;
                }
                if c_final != 1 {
                    all_ones = false;
                }
            }
        }

        bucket_counts[b] = out_keys.len() as u32;
        encode_bucket_keys(k, b, &out_keys, &mut key_payloads[b]);
        if !presence_only_out && !out_counts.is_empty() {
            encode_bucket_counts(&out_counts, &mut count_payloads[b]);
        }
    }

    if let Some(parent) = out_path.parent()
        && !parent.as_os_str().is_empty()
    {
        fs::create_dir_all(parent)?;
    }
    if presence_only_out {
        write_gr_kset_from_buckets(out_path, name, k, &bucket_counts, &key_payloads)?;
    } else {
        let counts_arg = if all_ones {
            vec![Vec::new(); BUCKETS]
        } else {
            count_payloads
        };
        write_ktab_from_buckets(
            out_path,
            name,
            k,
            1,
            max_count_out,
            all_ones,
            &bucket_counts,
            &key_payloads,
            &counts_arg,
        )?;
    }

    Ok(SetOpSummary {
        op,
        k,
        n_inputs,
        n_kmers_out,
        max_count_out,
        presence_only: presence_only_out || all_ones,
        fnv64: fnv.hex(),
        keys_sha256: hex_bytes(keys_sha.finish().as_ref()),
        output_path: out_path.display().to_string(),
    })
}

/// Streams a deterministic ML-ready feature matrix (`pa`, `counts`, or
/// `dosage`) plus explicit sample missingness mask under `out_dir`.
pub fn export_matrix(
    inputs: &[PathBuf],
    out_dir: &Path,
    cfg: &MatrixExportConfig,
) -> Result<MatrixManifest> {
    ensure!(
        !inputs.is_empty(),
        "matrix-export requires at least 1 input"
    );
    let readers: Vec<AnyTableReader> = inputs
        .iter()
        .map(|p| AnyTableReader::open(p))
        .collect::<Result<Vec<_>>>()?;
    let k = readers[0].k();
    for r in &readers[1..] {
        ensure!(
            r.k() == k,
            "k mismatch across inputs: '{}' has k={}, expected k={k}",
            r.name(),
            r.k()
        );
    }
    if cfg.mode == MatrixMode::Dosage {
        ensure!(
            cfg.core_panel.is_some(),
            "--core-panel is required for --mode dosage"
        );
    }
    if k <= MAX_K_U64 {
        export_matrix_generic::<u64>(&readers, k, out_dir, cfg)
    } else {
        export_matrix_generic::<u128>(&readers, k, out_dir, cfg)
    }
}

fn export_matrix_generic<K: KmerKey>(
    readers: &[AnyTableReader],
    k: usize,
    out_dir: &Path,
    cfg: &MatrixExportConfig,
) -> Result<MatrixManifest> {
    fs::create_dir_all(out_dir)
        .with_context(|| format!("create matrix export dir {}", out_dir.display()))?;
    let n_samples = readers.len();
    let space = K::new_space(k);

    let mut samples = Vec::with_capacity(n_samples);
    for (idx, r) in readers.iter().enumerate() {
        if cfg.mode == MatrixMode::Dosage {
            let panel = cfg.core_panel.as_ref().unwrap();
            let (status, d_core, breadth) = evaluate_table_core_coverage(r, panel, &cfg.dosage)?;
            samples.push(SampleMaskEntry {
                sample_idx: idx,
                sample: r.name().to_string(),
                missing: status == DosageStatus::InsufficientCoreCoverage,
                dosage_status: Some(status),
                d_core: Some(d_core),
                core_breadth: Some(breadth),
            });
        } else {
            samples.push(SampleMaskEntry {
                sample_idx: idx,
                sample: r.name().to_string(),
                missing: false,
                dosage_status: None,
                d_core: None,
                core_breadth: None,
            });
        }
    }
    let n_missing_samples = samples.iter().filter(|s| s.missing).count();

    // Write explicit missingness mask first.
    let mask_path = out_dir.join("missing_mask.tsv");
    let mut mask_sha = Sha256Context::new(&SHA256);
    {
        let mut mw = BufWriter::new(File::create(&mask_path)?);
        let header = "sample_idx\tsample\tmissing\tdosage_status\td_core\tcore_breadth\n";
        mw.write_all(header.as_bytes())?;
        mask_sha.update(header.as_bytes());
        for s in &samples {
            let status_str = match s.dosage_status {
                Some(DosageStatus::Ok) => "ok",
                Some(DosageStatus::InsufficientCoreCoverage) => "insufficient_core_coverage",
                None => "na",
            };
            let d_str = s
                .d_core
                .map_or_else(|| "NA".to_string(), |v| format!("{v:.6}"));
            let b_str = s
                .core_breadth
                .map_or_else(|| "NA".to_string(), |v| format!("{v:.6}"));
            let line = format!(
                "{}\t{}\t{}\t{}\t{}\t{}\n",
                s.sample_idx,
                s.sample,
                if s.missing { 1 } else { 0 },
                status_str,
                d_str,
                b_str
            );
            mw.write_all(line.as_bytes())?;
            mask_sha.update(line.as_bytes());
        }
        mw.flush()?;
    }

    let feat_path = out_dir.join("features.tsv");
    let matrix_file_name = match cfg.format {
        MatrixFormat::Sparse => "matrix.coo.tsv",
        MatrixFormat::Dense => "matrix.tsv",
    };
    let mat_path = out_dir.join(matrix_file_name);

    let mut feat_sha = Sha256Context::new(&SHA256);
    let mut mat_sha = Sha256Context::new(&SHA256);
    let mut feat_w = BufWriter::new(File::create(&feat_path)?);
    let mut mat_w = BufWriter::new(File::create(&mat_path)?);

    let feat_hdr = "feature_idx\tmixed_key_hex\tkmer_seq\tsupport\n";
    feat_w.write_all(feat_hdr.as_bytes())?;
    feat_sha.update(feat_hdr.as_bytes());

    match cfg.format {
        MatrixFormat::Sparse => {
            let hdr = "feature_idx\tsample_idx\tvalue\tis_missing\n";
            mat_w.write_all(hdr.as_bytes())?;
            mat_sha.update(hdr.as_bytes());
        }
        MatrixFormat::Dense => {
            let mut hdr = String::from("feature_idx");
            for s in &samples {
                hdr.push('\t');
                hdr.push_str(&s.sample);
            }
            hdr.push('\n');
            mat_w.write_all(hdr.as_bytes())?;
            mat_sha.update(hdr.as_bytes());
        }
    }

    let mut kbuf = Vec::new();
    let mut cbuf = Vec::new();
    let mut per_keys: Vec<Vec<K>> = (0..n_samples).map(|_| Vec::new()).collect();
    let mut per_counts: Vec<Vec<u32>> = (0..n_samples).map(|_| Vec::new()).collect();
    let mut cursors = vec![0usize; n_samples];
    let mut row_counts = vec![0u32; n_samples];
    let mut feature_idx = 0u64;

    for b in 0..BUCKETS {
        for i in 0..n_samples {
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
        }
        loop {
            let mut min_key: Option<K> = None;
            for i in 0..n_samples {
                if let Some(&key) = per_keys[i].get(cursors[i]) {
                    min_key = Some(match min_key {
                        None => key,
                        Some(cur) => cur.min(key),
                    });
                }
            }
            let Some(key) = min_key else {
                break;
            };
            let mut support = 0usize;
            for i in 0..n_samples {
                if per_keys[i].get(cursors[i]) == Some(&key) {
                    row_counts[i] = per_counts[i][cursors[i]];
                    cursors[i] += 1;
                    support += 1;
                } else {
                    row_counts[i] = 0;
                }
            }
            if support < cfg.min_support || cfg.max_support.is_some_and(|mx| support > mx) {
                continue;
            }
            let seq = K::decode(&space, K::unmix(&space, key));
            let feat_line = format!("{feature_idx}\t{}\t{seq}\t{support}\n", key.to_hex());
            feat_w.write_all(feat_line.as_bytes())?;
            feat_sha.update(feat_line.as_bytes());

            match cfg.format {
                MatrixFormat::Sparse => {
                    for (i, &c) in row_counts.iter().enumerate() {
                        if samples[i].missing {
                            // Explicitly record missingness rather than omitting as silent zero.
                            let line = format!("{feature_idx}\t{i}\tNA\t1\n");
                            mat_w.write_all(line.as_bytes())?;
                            mat_sha.update(line.as_bytes());
                            continue;
                        }
                        if c == 0 {
                            continue;
                        }
                        let line = match cfg.mode {
                            MatrixMode::Pa => format!("{feature_idx}\t{i}\t1\t0\n"),
                            MatrixMode::Counts => format!("{feature_idx}\t{i}\t{c}\t0\n"),
                            MatrixMode::Dosage => {
                                let d = samples[i].d_core.unwrap();
                                format!("{feature_idx}\t{i}\t{:.6}\t0\n", c as f64 / d)
                            }
                        };
                        mat_w.write_all(line.as_bytes())?;
                        mat_sha.update(line.as_bytes());
                    }
                }
                MatrixFormat::Dense => {
                    let mut line = format!("{feature_idx}");
                    for (i, &c) in row_counts.iter().enumerate() {
                        line.push('\t');
                        if samples[i].missing {
                            line.push_str("NA");
                        } else {
                            match cfg.mode {
                                MatrixMode::Pa => line.push(if c > 0 { '1' } else { '0' }),
                                MatrixMode::Counts => line.push_str(&c.to_string()),
                                MatrixMode::Dosage => {
                                    let d = samples[i].d_core.unwrap();
                                    line.push_str(&format!("{:.6}", c as f64 / d));
                                }
                            }
                        }
                    }
                    line.push('\n');
                    mat_w.write_all(line.as_bytes())?;
                    mat_sha.update(line.as_bytes());
                }
            }
            feature_idx += 1;
        }
    }
    feat_w.flush()?;
    mat_w.flush()?;

    let manifest = MatrixManifest {
        schema_version: 1,
        k,
        mode: cfg.mode,
        format: cfg.format,
        n_samples,
        n_features: feature_idx,
        n_missing_samples,
        samples,
        sha256_features: hex_bytes(feat_sha.finish().as_ref()),
        sha256_matrix: hex_bytes(mat_sha.finish().as_ref()),
        sha256_missing_mask: hex_bytes(mask_sha.finish().as_ref()),
    };
    write_json_atomic(&out_dir.join("manifest.json"), &manifest)?;
    Ok(manifest)
}
