//! Streaming bucket-by-bucket set algebra (`union`, `intersect`, `diff`,
//! `symdiff`) and bounded `matrix-export` (`pa`, `counts`, `sparse`, `dense`)
//! across `.kset` v1 (`GRMKSET1`), `.grk` (`GRMKGR1\0`), and `.ktab`
//! (`GRMKTB1\0`) tables.
//!
//! All entry points enforce `crate::matrix::REVIEWED_COLUMNS` (4,096 inputs)
//! before opening any input file, verify input SHA-256 digests when present,
//! enforce output path extensions, refuse silent `u32` count clamping unless
//! `--allow-saturation` is explicitly enabled, and default sample labels to
//! `"sample"` (for set ops) or `"s000001"`.. (for `matrix-export`) so input
//! filenames and paths are never copied into outputs.

use crate::confine;
use crate::fnv::Fnv;
use crate::kcount::format::{
    AnyTableReader, DEFAULT_SAMPLE_NAME, InputModality, KtabWriteMeta, MAX_COUNT, TableFormat,
    check_output_extension, encode_bucket_counts, encode_bucket_keys, hex_bytes,
    write_gr_kset_from_buckets, write_ktab_from_buckets,
};
use crate::kcount::key::{KmerKey, MAX_K_U64};
use crate::kcount::write_json_atomic;
use crate::kset::BUCKETS;
use crate::matrix::{REVIEWED_COLUMNS, ensure_reviewed};
use anyhow::{Context, Result, bail, ensure};
use ring::digest::{Context as Sha256Context, SHA256};
use rustc_hash::FxHashSet;
use serde::{Deserialize, Serialize};
use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

/// Set algebra operation applied across 2 to `REVIEWED_COLUMNS` k-mer tables.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SetOp {
    /// Keys present in at least one input; counts are summed with checked
    /// `u64` arithmetic and capped at `2^32 - 1` only under
    /// `--allow-saturation`.
    Union,
    /// Keys present in every input; counts take the minimum across inputs.
    Intersect,
    /// Keys present in input 0 and absent from all subsequent inputs.
    Diff,
    /// Keys present in an odd number of inputs (exact symmetric difference for
    /// 2 inputs; counts summed across carrying inputs).
    SymDiff,
}

/// Configuration for [`run_set_op`].
#[derive(Clone, Debug)]
pub struct SetOpConfig {
    pub sample_name: String,
    pub presence_only_out: bool,
    pub allow_saturation: bool,
}

impl Default for SetOpConfig {
    fn default() -> Self {
        Self {
            sample_name: DEFAULT_SAMPLE_NAME.to_string(),
            presence_only_out: false,
            allow_saturation: false,
        }
    }
}

/// Summary returned by [`run_set_op`]. Contains no input or output paths.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct SetOpSummary {
    pub op: SetOp,
    pub format: String,
    pub k: usize,
    pub sample: String,
    pub modality: InputModality,
    pub n_inputs: usize,
    pub n_kmers_out: u64,
    pub max_count_out: u32,
    pub saturated_counts: u64,
    pub counts_saturated: bool,
    pub presence_only: bool,
    pub fnv64: String,
    pub keys_sha256: String,
}

/// Feature value encoding mode for [`export_matrix`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum MatrixMode {
    Pa,
    Counts,
}

/// Output matrix file layout for [`export_matrix`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum MatrixLayout {
    Sparse,
    Dense,
}

/// Configuration for [`export_matrix`].
#[derive(Clone, Debug)]
pub struct MatrixExportConfig {
    pub mode: MatrixMode,
    pub layout: MatrixLayout,
    pub min_support: usize,
    pub max_support: Option<usize>,
    pub min_prevalence: Option<f64>,
    pub max_prevalence: Option<f64>,
    pub sample_names: Option<Vec<String>>,
}

impl Default for MatrixExportConfig {
    fn default() -> Self {
        Self {
            mode: MatrixMode::Pa,
            layout: MatrixLayout::Sparse,
            min_support: 1,
            max_support: None,
            min_prevalence: None,
            max_prevalence: None,
            sample_names: None,
        }
    }
}

/// Per-sample missingness metadata recorded in `missing_mask.tsv` and
/// `manifest.json`.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct SampleMaskEntry {
    pub sample_idx: usize,
    pub sample: String,
    pub missing_core_coverage: bool,
    pub dosage_status: Option<String>,
    pub d_core: Option<f64>,
    pub core_breadth: Option<f64>,
}

/// Deterministic manifest written to `<out_dir>/manifest.json`. Contains no
/// input file stems or paths.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct MatrixManifest {
    pub schema_version: u32,
    pub k: usize,
    pub mode: MatrixMode,
    pub layout: MatrixLayout,
    pub n_samples: usize,
    pub n_features: u64,
    pub n_missing_samples: usize,
    pub samples: Vec<SampleMaskEntry>,
    pub sha256_features: String,
    pub sha256_matrix: String,
    pub sha256_missing_mask: String,
}

struct HashedPartWriter {
    path: PathBuf,
    part: PathBuf,
    writer: Option<BufWriter<File>>,
    sha: Sha256Context,
}

impl HashedPartWriter {
    fn create(path: PathBuf) -> Result<Self> {
        confine::regular_or_absent(&path)?;
        let part = confine::suffixed(&path, ".part");
        let file = confine::create_part(&part)?;
        Ok(Self {
            path,
            part,
            writer: Some(BufWriter::new(file)),
            sha: Sha256Context::new(&SHA256),
        })
    }

    fn write_all(&mut self, bytes: &[u8]) -> Result<()> {
        self.writer
            .as_mut()
            .expect("part writer still open")
            .write_all(bytes)?;
        self.sha.update(bytes);
        Ok(())
    }

    fn finish(mut self) -> Result<String> {
        if let Some(mut w) = self.writer.take() {
            w.flush()?;
        }
        confine::replace_with(&self.part, &self.path)?;
        Ok(hex_bytes(self.sha.clone().finish().as_ref()))
    }
}

impl Drop for HashedPartWriter {
    fn drop(&mut self) {
        self.writer.take();
        if confine::regular_or_absent(&self.part).unwrap_or(false) {
            let _ = fs::remove_file(&self.part);
        }
    }
}

/// Propagates `InputModality::Reads` or `InputModality::Assembly` only when
/// every input table agrees on that modality; otherwise returns `Unknown`.
fn resolve_output_modality(readers: &[AnyTableReader]) -> InputModality {
    let mut shared: Option<InputModality> = None;
    for r in readers {
        let m = r.modality().unwrap_or(InputModality::Unknown);
        match m {
            InputModality::Unknown => return InputModality::Unknown,
            InputModality::Assembly | InputModality::Reads => match shared {
                None => shared = Some(m),
                Some(prev) if prev == m => {}
                Some(_) => return InputModality::Unknown,
            },
        }
    }
    shared.unwrap_or(InputModality::Unknown)
}

fn resolve_sample_names(n_samples: usize, supplied: Option<&[String]>) -> Result<Vec<String>> {
    match supplied {
        None => Ok((0..n_samples).map(|i| format!("s{:06}", i + 1)).collect()),
        Some(names) => {
            ensure!(
                names.len() == n_samples,
                "--sample-names provided {} names for {n_samples} inputs",
                names.len()
            );
            let mut seen = FxHashSet::default();
            for name in names {
                ensure!(!name.is_empty(), "sample name must not be empty");
                ensure!(
                    !name.contains(['\t', '\n', '\r']),
                    "sample name {name:?} must not contain tabs or newlines"
                );
                ensure!(seen.insert(name.as_str()), "duplicate sample name {name:?}");
            }
            Ok(names.to_vec())
        }
    }
}

fn open_and_verify_tables(inputs: &[PathBuf]) -> Result<(Vec<AnyTableReader>, usize)> {
    let mut readers = Vec::with_capacity(inputs.len());
    for p in inputs {
        let r = AnyTableReader::open(p)?;
        r.verify_integrity()?;
        readers.push(r);
    }
    let k = readers[0].k();
    for (idx, r) in readers.iter().enumerate().skip(1) {
        ensure!(
            r.k() == k,
            "k mismatch across inputs: input #{idx} has k={}, expected k={k}",
            r.k()
        );
    }
    Ok((readers, k))
}

/// Runs a streaming bucket-by-bucket merge-join set operation (`union`,
/// `intersect`, `diff`, or `symdiff`) across 2 to `REVIEWED_COLUMNS` tables.
pub fn run_set_op(
    op: SetOp,
    inputs: &[PathBuf],
    out_path: &Path,
    cfg: &SetOpConfig,
) -> Result<SetOpSummary> {
    ensure!(
        inputs.len() >= 2,
        "set operation requires at least 2 inputs"
    );
    ensure_reviewed(inputs.len(), REVIEWED_COLUMNS)?;
    let out_fmt = if cfg.presence_only_out {
        TableFormat::GrKset
    } else {
        TableFormat::Ktab
    };
    check_output_extension(out_path, out_fmt)?;

    let (readers, k) = open_and_verify_tables(inputs)?;
    let modality = resolve_output_modality(&readers);
    if k <= MAX_K_U64 {
        run_set_op_generic::<u64>(op, &readers, k, modality, out_path, cfg)
    } else {
        run_set_op_generic::<u128>(op, &readers, k, modality, out_path, cfg)
    }
}

fn run_set_op_generic<K: KmerKey>(
    op: SetOp,
    readers: &[AnyTableReader],
    k: usize,
    modality: InputModality,
    out_path: &Path,
    cfg: &SetOpConfig,
) -> Result<SetOpSummary> {
    let n_inputs = readers.len();
    let mut bucket_counts = [0u32; BUCKETS];
    let mut key_payloads = vec![Vec::new(); BUCKETS];
    let mut count_payloads = vec![Vec::new(); BUCKETS];

    let mut fnv = Fnv::new();
    let mut keys_sha = Sha256Context::new(&SHA256);
    let mut n_kmers_out = 0u64;
    let mut max_count_out = 0u32;
    let mut saturated_counts = 0u64;
    let mut all_ones = true;

    let mut kbuf = Vec::new();
    let mut cbuf = Vec::new();
    let mut per_input_keys: Vec<Vec<K>> = (0..n_inputs).map(|_| Vec::new()).collect();
    let mut per_input_counts: Vec<Vec<u32>> = (0..n_inputs).map(|_| Vec::new()).collect();
    let mut cursors = vec![0usize; n_inputs];
    let mut heap: BinaryHeap<Reverse<(K, usize)>> = BinaryHeap::with_capacity(n_inputs);
    let mut out_keys: Vec<K> = Vec::new();
    let mut out_counts: Vec<u32> = Vec::new();

    for b in 0..BUCKETS {
        heap.clear();
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
            if let Some(&first_key) = per_input_keys[i].first() {
                heap.push(Reverse((first_key, i)));
            }
        }
        out_keys.clear();
        out_counts.clear();

        // Streaming k-way merge-join over sorted keys in bucket b.
        while let Some(Reverse((key, mut idx))) = heap.pop() {
            let mut present_count = 0usize;
            let mut in_first = false;
            let mut first_count = 0u32;
            let mut sum_count_u64 = 0u64;
            let mut min_count = u32::MAX;

            loop {
                let pos = cursors[idx];
                let c = per_input_counts[idx][pos];
                cursors[idx] = pos + 1;
                if let Some(&next_key) = per_input_keys[idx].get(pos + 1) {
                    heap.push(Reverse((next_key, idx)));
                }
                present_count += 1;
                if idx == 0 {
                    in_first = true;
                    first_count = c;
                }
                sum_count_u64 = sum_count_u64
                    .checked_add(c as u64)
                    .context("k-mer count sum overflowed u64")?;
                min_count = min_count.min(c);

                match heap.peek() {
                    Some(&Reverse((next_k, _))) if next_k == key => {
                        idx = heap.pop().unwrap().0.1;
                    }
                    _ => break,
                }
            }

            let emitted_u64: Option<u64> = match op {
                SetOp::Union => Some(sum_count_u64),
                SetOp::Intersect if present_count == n_inputs => Some(min_count as u64),
                SetOp::Diff if in_first && present_count == 1 => Some(first_count as u64),
                SetOp::SymDiff if present_count % 2 == 1 => Some(sum_count_u64),
                _ => None,
            };

            if let Some(raw_c) = emitted_u64 {
                let c_final = if cfg.presence_only_out {
                    1
                } else if raw_c > MAX_COUNT as u64 {
                    if !cfg.allow_saturation {
                        bail!(
                            "k-mer count sum {raw_c} exceeds u32 cap {MAX_COUNT} (2^32 - 1); pass --allow-saturation to clamp counts at 2^32 - 1"
                        );
                    }
                    saturated_counts = saturated_counts
                        .checked_add(1)
                        .context("saturated count overflow")?;
                    MAX_COUNT
                } else {
                    raw_c as u32
                };
                out_keys.push(key);
                out_counts.push(c_final);
                key.with_le_bytes(|bytes| {
                    fnv.bytes(bytes);
                    keys_sha.update(bytes);
                });
                n_kmers_out = n_kmers_out
                    .checked_add(1)
                    .context("output k-mer count overflow")?;
                if c_final > max_count_out {
                    max_count_out = c_final;
                }
                if c_final != 1 {
                    all_ones = false;
                }
            }
        }

        bucket_counts[b] = u32::try_from(out_keys.len()).context("bucket count exceeds u32")?;
        encode_bucket_keys(k, b, &out_keys, &mut key_payloads[b]);
        if !cfg.presence_only_out && !out_counts.is_empty() {
            encode_bucket_counts(&out_counts, &mut count_payloads[b]);
        }
    }

    if let Some(parent) = out_path.parent()
        && !parent.as_os_str().is_empty()
    {
        fs::create_dir_all(parent)?;
    }
    let out_fmt = if cfg.presence_only_out {
        write_gr_kset_from_buckets(
            out_path,
            &cfg.sample_name,
            k,
            modality,
            &bucket_counts,
            &key_payloads,
        )?;
        TableFormat::GrKset
    } else {
        let presence_only = all_ones;
        let counts_arg = if presence_only {
            vec![Vec::new(); BUCKETS]
        } else {
            count_payloads
        };
        let meta = KtabWriteMeta {
            name: &cfg.sample_name,
            k,
            modality,
            min_count: 1,
            max_count_observed: max_count_out,
            saturated_counts,
            presence_only,
        };
        write_ktab_from_buckets(out_path, &meta, &bucket_counts, &key_payloads, &counts_arg)?;
        TableFormat::Ktab
    };

    Ok(SetOpSummary {
        op,
        format: out_fmt.magic_label().to_string(),
        k,
        sample: cfg.sample_name.clone(),
        modality,
        n_inputs,
        n_kmers_out,
        max_count_out,
        saturated_counts,
        counts_saturated: saturated_counts > 0,
        presence_only: cfg.presence_only_out || all_ones,
        fnv64: fnv.hex(),
        keys_sha256: hex_bytes(keys_sha.finish().as_ref()),
    })
}

/// Streams a bounded feature matrix (`pa` or `counts`, `sparse` or `dense`)
/// plus `features.tsv`, `missing_mask.tsv`, and `manifest.json` under
/// `out_dir`.
pub fn export_matrix(
    inputs: &[PathBuf],
    out_dir: &Path,
    cfg: &MatrixExportConfig,
) -> Result<MatrixManifest> {
    ensure!(
        !inputs.is_empty(),
        "matrix-export requires at least 1 input"
    );
    ensure_reviewed(inputs.len(), REVIEWED_COLUMNS)?;
    ensure!(cfg.min_support >= 1, "min_support must be >= 1");
    if let Some(mx) = cfg.max_support {
        ensure!(
            mx >= cfg.min_support,
            "max_support ({mx}) must be >= min_support ({})",
            cfg.min_support
        );
    }
    if let Some(mp) = cfg.min_prevalence {
        ensure!(
            mp.is_finite() && (0.0..=1.0).contains(&mp),
            "min_prevalence must be within [0.0, 1.0]"
        );
    }
    if let Some(mxp) = cfg.max_prevalence {
        ensure!(
            mxp.is_finite() && (0.0..=1.0).contains(&mxp),
            "max_prevalence must be within [0.0, 1.0]"
        );
    }
    if let (Some(mp), Some(mxp)) = (cfg.min_prevalence, cfg.max_prevalence) {
        ensure!(
            mxp >= mp,
            "max_prevalence ({mxp}) must be >= min_prevalence ({mp})"
        );
    }
    let sample_names = resolve_sample_names(inputs.len(), cfg.sample_names.as_deref())?;

    let (readers, k) = open_and_verify_tables(inputs)?;
    if k <= MAX_K_U64 {
        export_matrix_generic::<u64>(&readers, &sample_names, k, out_dir, cfg)
    } else {
        export_matrix_generic::<u128>(&readers, &sample_names, k, out_dir, cfg)
    }
}

fn export_matrix_generic<K: KmerKey>(
    readers: &[AnyTableReader],
    sample_names: &[String],
    k: usize,
    out_dir: &Path,
    cfg: &MatrixExportConfig,
) -> Result<MatrixManifest> {
    fs::create_dir_all(out_dir)
        .with_context(|| format!("create matrix export dir {}", out_dir.display()))?;
    let n_samples = readers.len();
    let space = K::new_space(k);

    let samples: Vec<SampleMaskEntry> = sample_names
        .iter()
        .enumerate()
        .map(|(idx, name)| SampleMaskEntry {
            sample_idx: idx,
            sample: name.clone(),
            missing_core_coverage: false,
            dosage_status: None,
            d_core: None,
            core_breadth: None,
        })
        .collect();
    let n_missing_samples = samples.iter().filter(|s| s.missing_core_coverage).count();

    let sha256_missing_mask = {
        let mut mw = HashedPartWriter::create(out_dir.join("missing_mask.tsv"))?;
        mw.write_all(
            b"sample_idx\tsample\tmissing_core_coverage\tdosage_status\td_core\tcore_breadth\n",
        )?;
        for s in &samples {
            let status_str = s.dosage_status.as_deref().unwrap_or("na");
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
                if s.missing_core_coverage { 1 } else { 0 },
                status_str,
                d_str,
                b_str
            );
            mw.write_all(line.as_bytes())?;
        }
        mw.finish()?
    };

    let mut feat_w = HashedPartWriter::create(out_dir.join("features.tsv"))?;
    let mut mat_w = HashedPartWriter::create(out_dir.join("matrix.tsv"))?;

    feat_w.write_all(b"feature_idx\tmixed_key_hex\tkmer_seq\tsupport\n")?;
    match cfg.layout {
        MatrixLayout::Sparse => {
            mat_w.write_all(b"feature_idx\tsample_idx\tvalue\n")?;
        }
        MatrixLayout::Dense => {
            let mut hdr = String::from("feature_idx");
            for s in &samples {
                hdr.push('\t');
                hdr.push_str(&s.sample);
            }
            hdr.push('\n');
            mat_w.write_all(hdr.as_bytes())?;
        }
    }

    let mut kbuf = Vec::new();
    let mut cbuf = Vec::new();
    let mut per_keys: Vec<Vec<K>> = (0..n_samples).map(|_| Vec::new()).collect();
    let mut per_counts: Vec<Vec<u32>> = (0..n_samples).map(|_| Vec::new()).collect();
    let mut cursors = vec![0usize; n_samples];
    let mut heap: BinaryHeap<Reverse<(K, usize)>> = BinaryHeap::with_capacity(n_samples);
    let mut present_row: Vec<(usize, u32)> = Vec::with_capacity(n_samples);
    let mut dense_row = vec![0u32; n_samples];
    let mut feature_idx = 0u64;

    for b in 0..BUCKETS {
        heap.clear();
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
            if let Some(&first_key) = per_keys[i].first() {
                heap.push(Reverse((first_key, i)));
            }
        }

        while let Some(Reverse((key, mut idx))) = heap.pop() {
            present_row.clear();
            loop {
                let pos = cursors[idx];
                let c = per_counts[idx][pos];
                cursors[idx] = pos + 1;
                if let Some(&next_key) = per_keys[idx].get(pos + 1) {
                    heap.push(Reverse((next_key, idx)));
                }
                present_row.push((idx, c));

                match heap.peek() {
                    Some(&Reverse((next_k, _))) if next_k == key => {
                        idx = heap.pop().unwrap().0.1;
                    }
                    _ => break,
                }
            }

            let support = present_row.len();
            if support < cfg.min_support || cfg.max_support.is_some_and(|mx| support > mx) {
                continue;
            }
            let prevalence = support as f64 / n_samples as f64;
            if cfg.min_prevalence.is_some_and(|mp| prevalence < mp)
                || cfg.max_prevalence.is_some_and(|mxp| prevalence > mxp)
            {
                continue;
            }

            let seq = K::decode(&space, K::unmix(&space, key));
            let feat_line = format!("{feature_idx}\t{}\t{seq}\t{support}\n", key.to_hex());
            feat_w.write_all(feat_line.as_bytes())?;

            match cfg.layout {
                MatrixLayout::Sparse => {
                    for &(i, c) in &present_row {
                        if samples[i].missing_core_coverage || c == 0 {
                            continue;
                        }
                        let line = match cfg.mode {
                            MatrixMode::Pa => format!("{feature_idx}\t{i}\t1\n"),
                            MatrixMode::Counts => format!("{feature_idx}\t{i}\t{c}\n"),
                        };
                        mat_w.write_all(line.as_bytes())?;
                    }
                }
                MatrixLayout::Dense => {
                    for &(i, c) in &present_row {
                        dense_row[i] = c;
                    }
                    let mut line = format!("{feature_idx}");
                    for (i, &c) in dense_row.iter().enumerate() {
                        line.push('\t');
                        if samples[i].missing_core_coverage {
                            line.push_str("NA");
                        } else {
                            match cfg.mode {
                                MatrixMode::Pa => line.push(if c > 0 { '1' } else { '0' }),
                                MatrixMode::Counts => line.push_str(&c.to_string()),
                            }
                        }
                    }
                    line.push('\n');
                    mat_w.write_all(line.as_bytes())?;
                    for &(i, _) in &present_row {
                        dense_row[i] = 0;
                    }
                }
            }
            feature_idx = feature_idx
                .checked_add(1)
                .context("feature index overflow")?;
        }
    }

    let sha256_features = feat_w.finish()?;
    let sha256_matrix = mat_w.finish()?;

    let manifest = MatrixManifest {
        schema_version: 1,
        k,
        mode: cfg.mode,
        layout: cfg.layout,
        n_samples,
        n_features: feature_idx,
        n_missing_samples,
        samples,
        sha256_features,
        sha256_matrix,
        sha256_missing_mask,
    };
    write_json_atomic(&out_dir.join("manifest.json"), &manifest)?;
    Ok(manifest)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kcount::format::{GrKsetReader, KtabReader, write_gr_kset, write_ktab};
    use crate::kcount::key::{KmerSpace64, KmerSpace128};

    fn make_ktab_u64(path: &Path, k: usize, modality: InputModality, items: &[(&str, u32)]) {
        let sp = KmerSpace64::new(k);
        let mut entries: Vec<(u64, u32)> = items
            .iter()
            .map(|&(s, c)| (sp.mix(sp.canonical(s).unwrap()), c))
            .collect();
        entries.sort_unstable_by_key(|&(key, _)| key);
        let sat = entries.iter().filter(|&&(_, c)| c == MAX_COUNT).count() as u64;
        write_ktab(path, DEFAULT_SAMPLE_NAME, k, modality, 1, sat, &entries).unwrap();
    }

    #[test]
    fn union_checked_summation_and_opt_in_saturation() {
        let dir = tempfile::tempdir().unwrap();
        let k = 21;
        let k1 = "ACGTACGTACGTACGTACGTA";
        let k2 = "TGCATGCATGCATGCATGCAT";
        let a = dir.path().join("a.ktab");
        let b = dir.path().join("b.ktab");
        make_ktab_u64(&a, k, InputModality::Reads, &[(k1, 5), (k2, 3_000_000_000)]);
        make_ktab_u64(&b, k, InputModality::Reads, &[(k1, 7), (k2, 2_000_000_000)]);

        // Without --allow-saturation, 3e9 + 2e9 > u32::MAX MUST return an error.
        let out_fail = dir.path().join("u_fail.ktab");
        let err = run_set_op(
            SetOp::Union,
            &[a.clone(), b.clone()],
            &out_fail,
            &SetOpConfig::default(),
        )
        .unwrap_err();
        assert!(
            format!("{err:#}").contains("allow-saturation"),
            "unexpected error: {err:#}"
        );
        assert!(!out_fail.exists());

        // With --allow-saturation, exact sum 5 + 7 = 12 and 3e9 + 2e9 clamps to u32::MAX
        // with saturated_counts recorded in both .ktab header and JSON summary.
        let out_sat = dir.path().join("u_sat.ktab");
        let summary = run_set_op(
            SetOp::Union,
            &[a, b],
            &out_sat,
            &SetOpConfig {
                allow_saturation: true,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(summary.n_kmers_out, 2);
        assert_eq!(summary.max_count_out, u32::MAX);
        assert_eq!(summary.saturated_counts, 1);
        assert!(summary.counts_saturated);
        assert_eq!(summary.modality, InputModality::Reads);

        let reader = KtabReader::open_verified(&out_sat).unwrap();
        assert_eq!(reader.saturated_counts, 1);
        assert_eq!(reader.max_count_observed, u32::MAX);
        assert_eq!(reader.modality, InputModality::Reads);
        let sp = KmerSpace64::new(k);
        let got: std::collections::HashMap<u64, u32> =
            reader.read_all::<u64>().unwrap().into_iter().collect();
        assert_eq!(got[&sp.mix(sp.canonical(k1).unwrap())], 12);
        assert_eq!(got[&sp.mix(sp.canonical(k2).unwrap())], u32::MAX);
    }

    #[test]
    fn reviewed_columns_cap_rejects_before_opening_files() {
        let dir = tempfile::tempdir().unwrap();
        let bogus_inputs: Vec<PathBuf> = (0..=REVIEWED_COLUMNS)
            .map(|i| dir.path().join(format!("does_not_exist_{i}.ktab")))
            .collect();

        for op in [SetOp::Union, SetOp::Intersect, SetOp::Diff, SetOp::SymDiff] {
            let err = run_set_op(
                op,
                &bogus_inputs,
                &dir.path().join("out.ktab"),
                &SetOpConfig::default(),
            )
            .unwrap_err();
            let msg = format!("{err:#}");
            assert!(
                msg.contains("4096"),
                "op {op:?} did not enforce REVIEWED_COLUMNS before opening files: {msg}"
            );
        }

        let err = export_matrix(
            &bogus_inputs,
            &dir.path().join("mat"),
            &MatrixExportConfig::default(),
        )
        .unwrap_err();
        assert!(
            format!("{err:#}").contains("4096"),
            "matrix-export did not enforce REVIEWED_COLUMNS before opening files: {err:#}"
        );
    }

    #[test]
    fn modality_propagation_and_extension_validation() {
        let dir = tempfile::tempdir().unwrap();
        let k = 45;
        let sp = KmerSpace128::new(k);
        let s1 = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTA";
        let s2 = "TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCAT";
        let m1 = sp.mix(sp.canonical(s1).unwrap());
        let m2 = sp.mix(sp.canonical(s2).unwrap());
        let mut keys = vec![m1, m2];
        keys.sort_unstable();

        let asm_a = dir.path().join("asm_a.grk");
        let asm_b = dir.path().join("asm_b.grk");
        let reads_c = dir.path().join("reads_c.ktab");
        write_gr_kset(&asm_a, "sample", k, InputModality::Assembly, &keys).unwrap();
        write_gr_kset(&asm_b, "sample", k, InputModality::Assembly, &keys[..1]).unwrap();
        write_ktab(
            &reads_c,
            "sample",
            k,
            InputModality::Reads,
            1,
            0,
            &[(keys[0], 4)],
        )
        .unwrap();

        // Matching Assembly inputs -> Assembly output.
        let out_gr = dir.path().join("same_mod.grk");
        let sum_same = run_set_op(
            SetOp::Union,
            &[asm_a.clone(), asm_b],
            &out_gr,
            &SetOpConfig {
                presence_only_out: true,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(sum_same.modality, InputModality::Assembly);
        assert_eq!(
            GrKsetReader::open_verified(&out_gr).unwrap().modality,
            InputModality::Assembly
        );

        // Mixed Assembly + Reads inputs -> Unknown output.
        let out_mixed = dir.path().join("mixed_mod.ktab");
        let sum_mix = run_set_op(
            SetOp::Union,
            &[asm_a.clone(), reads_c.clone()],
            &out_mixed,
            &SetOpConfig::default(),
        )
        .unwrap();
        assert_eq!(sum_mix.modality, InputModality::Unknown);
        assert_eq!(
            KtabReader::open_verified(&out_mixed).unwrap().modality,
            InputModality::Unknown
        );

        // Extension mismatches refused cleanly.
        assert!(
            run_set_op(
                SetOp::Union,
                &[asm_a.clone(), reads_c.clone()],
                &dir.path().join("wrong.grk"),
                &SetOpConfig::default(),
            )
            .is_err()
        );
        assert!(
            run_set_op(
                SetOp::Union,
                &[asm_a, reads_c],
                &dir.path().join("wrong.ktab"),
                &SetOpConfig {
                    presence_only_out: true,
                    ..Default::default()
                },
            )
            .is_err()
        );
    }
}
