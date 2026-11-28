//! Pure-Rust k-mer key primitives, verified binary formats, streaming set
//! algebra, bounded matrix export, and modality-verified core-panel dosage
//! extraction (`grm kcount`).
//!
//! Submodules:
//! - [`key`]: dual-width canonical k-mer spaces (`u64` for `k = 5..=32`,
//!   `u128` for `k = 33..=64`) with branchless O(1) SWAR reverse complement
//!   and invertible `2k`-bit mixers.
//! - [`mod@format`]: verified Golomb-Rice sorted-key table (`.grk`,
//!   `GRMKGR1\0`), counted table (`.ktab`, `GRMKTB1\0`), per-sequence count
//!   profile (`.prof.bin`, `GRMPRF1\0`), and `.kset` v1 (`GRMKSET1`) unified
//!   reader.
//! - [`ops`]: streaming bucket-by-bucket set algebra (`union`, `intersect`,
//!   `diff`, `symdiff`) and bounded `matrix-export` (`pa`, `counts`, `dosage`,
//!   `sparse`, `dense`).
//! - [`dosage`]: modality-verified core-panel depth normalization (`D_core`),
//!   target/plasmid dosage extraction, and training-fold core-panel
//!   construction.

pub mod dosage;
pub mod format;
pub mod key;
pub mod ops;

use crate::confine;
use crate::fnv::Fnv;
use crate::kcount::dosage::{
    CorePanel, CorePanelBuildConfig, DosageConfig, InputType, TargetPanel,
    build_training_fold_core_panel, extract_sample_dosage, load_qc_metrics,
};
use crate::kcount::format::{
    AnyTableReader, DEFAULT_SAMPLE_NAME, InputModality, MAGIC_PROF1, TableFormat, convert_to_gr,
    convert_to_ktab, convert_to_v1, detect_format, hex_bytes, read_profile_header,
    read_profiles_verified,
};
use crate::kcount::key::{KmerKey, MAX_K, MAX_K_U64, MIN_K};
use crate::kcount::ops::{
    MatrixExportConfig, MatrixLayout, MatrixMode, SetOp, SetOpConfig, export_matrix, run_set_op,
};
use crate::kset::BUCKETS;
use anyhow::{Context, Result, bail, ensure};
use clap::builder::RangedU64ValueParser;
use clap::{Args, Subcommand, ValueEnum};
use ring::digest::{Context as Sha256Context, SHA256};
use serde::Serialize;
use std::fs::{self, File};
use std::io::Read;
use std::path::{Path, PathBuf};

/// CLI arguments for `grm kcount`.
#[derive(Args, Debug)]
pub struct KcountArgs {
    #[command(subcommand)]
    pub command: KcountCommand,
}

#[derive(Subcommand, Debug)]
pub enum KcountCommand {
    /// Inspect metadata, bucket statistics, sketch prefixes, or verified
    /// SHA-256 and key digests for `.kset` v1, `.grk`, `.ktab`, or `.prof.bin`
    /// files.
    Inspect {
        /// File to inspect.
        #[arg(long)]
        input: PathBuf,
        /// Optional bottom-`s` MinHash sketch prefix length to include.
        #[arg(long)]
        sketch: Option<usize>,
        /// Recompute stored header SHA-256 digests and decode every bucket or
        /// sequence record.
        #[arg(long)]
        verify: bool,
        /// Optional JSON output path (prints to stdout when omitted).
        #[arg(long)]
        out: Option<PathBuf>,
    },
    /// Convert between `.kset` v1 (`GRMKSET1`), `.grk` (`GRMKGR1`), and
    /// `.ktab` (`GRMKTB1`) tables, refusing output extension mismatches.
    ConvertV1 {
        /// Source table (`.kset` v1, `.grk`, or `.ktab`).
        #[arg(long)]
        input: PathBuf,
        /// Output table path; extension must match the format written
        /// (`.grk`, `.ktab`, or `.kset`).
        #[arg(long)]
        out: PathBuf,
        /// Target conversion direction (`to-gr`, `to-ktab`, or `to-v1`;
        /// inferred from source format and `--out` extension when omitted).
        #[arg(long, value_enum)]
        direction: Option<ConvertDirectionArg>,
        /// Sample name written into the output header (defaults to `"sample"`
        /// so input file names are never copied into output headers).
        #[arg(long)]
        name: Option<String>,
        /// Input modality recorded in `.grk`/`.ktab` headers (`unknown`,
        /// `assembly`, or `reads`; preserves source modality when omitted for
        /// `.grk`/`.ktab` inputs, defaults to `unknown` for `.kset` v1).
        #[arg(long, value_enum)]
        modality: Option<InputModalityArg>,
    },
    /// Streaming bucket-by-bucket set algebra (`union`, `intersect`, `diff`,
    /// `symdiff`) across 2 to 4,096 `.kset` v1, `.grk`, or `.ktab` tables.
    Op {
        #[command(subcommand)]
        action: OpCommand,
    },
    /// Stream a bounded sparse or dense feature matrix (`pa`, `counts`, or
    /// `dosage`) across 1 to 4,096 tables, defaulting column names to
    /// `"s000001"`..
    MatrixExport {
        /// One to 4,096 input tables (`.kset` v1, `.grk`, or `.ktab`; `--mode
        /// dosage` requires verified counted `.ktab` tables matching the core
        /// panel modality).
        #[arg(long, required = true)]
        input: Vec<PathBuf>,
        /// Output directory receiving `features.tsv`, `matrix.tsv`,
        /// `missing_mask.tsv`, and `manifest.json`.
        #[arg(long)]
        out_dir: PathBuf,
        /// Feature value encoding (`pa`, `counts`, or `dosage`).
        #[arg(long, value_enum, default_value_t = MatrixModeArg::Pa)]
        mode: MatrixModeArg,
        /// Matrix output layout (`sparse` or `dense`).
        #[arg(long, alias = "format", value_enum, default_value_t = MatrixLayoutArg::Sparse)]
        layout: MatrixLayoutArg,
        /// Minimum number of samples carrying a k-mer (`>= 1`).
        #[arg(long, default_value_t = 1)]
        min_support: usize,
        /// Optional maximum number of samples carrying a k-mer.
        #[arg(long)]
        max_support: Option<usize>,
        /// Optional minimum sample prevalence fraction in `[0.0, 1.0]`.
        #[arg(long)]
        min_prevalence: Option<f64>,
        /// Optional maximum sample prevalence fraction in `[0.0, 1.0]`.
        #[arg(long)]
        max_prevalence: Option<f64>,
        /// Optional explicit sample column names (one per input, comma- or
        /// flag-separated; defaults to `"s000001"`, `"s000002"`, ... so input
        /// file names are never copied into outputs).
        #[arg(long, value_delimiter = ',')]
        sample_names: Vec<String>,
        /// Core k-mer panel JSON with canonical `panel_sha256` (required for
        /// `--mode dosage`).
        #[arg(long)]
        core_panel: Option<PathBuf>,
        /// Optional input modality assertion (`reads` or `assembly`) checked
        /// alongside each table's stored `.ktab` header modality.
        #[arg(long, value_enum)]
        input_type: Option<InputTypeArg>,
        /// Optional fold identifier verified against `panel.fold_id`.
        #[arg(long)]
        expected_fold_id: Option<String>,
        /// Minimum count required for a core k-mer to count as detected in
        /// `--mode dosage`.
        #[arg(long, default_value_t = 1)]
        min_count: u32,
        /// Minimum median core depth `D_core` required in `--mode dosage`.
        #[arg(long, default_value_t = 10.0)]
        min_core_depth: f64,
        /// Minimum fraction of core k-mers detected in `--mode dosage`.
        #[arg(long, default_value_t = 0.80)]
        min_core_breadth: f64,
    },
    /// Build a single-copy core k-mer panel JSON with canonical `panel_sha256`
    /// from 1 to 4,096 training-fold tables (caller must exclude held-out test
    /// genomes).
    BuildCorePanel {
        /// One to 4,096 training-fold `.ktab` or `.grk` tables.
        #[arg(long, required = true)]
        input: Vec<PathBuf>,
        /// K-mer length (`5..=64`), validated against every input table.
        #[arg(
            long,
            value_parser = RangedU64ValueParser::<usize>::new().range(MIN_K as u64..=MAX_K as u64)
        )]
        k: usize,
        /// Input modality (`reads` or `assembly`), validated against every
        /// input table header.
        #[arg(long, value_enum)]
        input_type: InputTypeArg,
        /// Minimum fraction of training tables carrying a k-mer (`(0, 1]`).
        #[arg(long, default_value_t = 0.95)]
        min_prevalence: f64,
        /// Maximum ratio of a k-mer's count to its training table's median
        /// depth (`>= 1.0`); k-mers exceeding this in any training table are
        /// excluded.
        #[arg(long, default_value_t = 1.5)]
        max_copy_factor: f64,
        /// Maximum number of core k-mers to retain in the panel.
        #[arg(long, default_value_t = 10_000)]
        max_kmers: usize,
        /// Optional fold identifier recorded in the panel metadata.
        #[arg(long)]
        fold_id: Option<String>,
        /// Optional output core panel JSON path (prints to stdout when
        /// omitted).
        #[arg(long)]
        out: Option<PathBuf>,
    },
    /// Modality-verified `D_core` depth normalization and target/plasmid
    /// dosage extraction from a verified counted `.ktab` table.
    Dosage {
        /// Counted `.ktab` (`GRMKTB1\0`) file for the sample.
        #[arg(long)]
        table: PathBuf,
        /// Single-copy core k-mer panel JSON with mandatory `panel_sha256`.
        #[arg(long)]
        core_panel: PathBuf,
        /// Sample label written into the JSON report (defaults to `"sample"`
        /// so input file paths are never copied into outputs).
        #[arg(long)]
        name: Option<String>,
        /// Optional input modality assertion (`reads` or `assembly`); when
        /// set, both this flag and the `.ktab` header's stored modality must
        /// match `panel.input_type`.
        #[arg(long, value_enum)]
        input_type: Option<InputTypeArg>,
        /// Optional fold identifier verified against `panel.fold_id`.
        #[arg(long)]
        expected_fold_id: Option<String>,
        /// Optional target/plasmid k-mer panel JSON.
        #[arg(long)]
        targets: Option<PathBuf>,
        /// Optional metadata or QC JSON (matching `CountSummary` or
        /// `ReadQcMetrics` fields). When omitted, every QC field serializes as
        /// JSON `null`.
        #[arg(long)]
        qc_json: Option<PathBuf>,
        /// Minimum count required for a k-mer to count as detected.
        #[arg(long, default_value_t = 1)]
        min_count: u32,
        /// Minimum median core depth `D_core` required for valid dosage calls.
        #[arg(long, default_value_t = 10.0)]
        min_core_depth: f64,
        /// Minimum fraction of core k-mers detected for valid dosage calls.
        #[arg(long, default_value_t = 0.80)]
        min_core_breadth: f64,
        /// Optional JSON output path (prints to stdout when omitted).
        #[arg(long)]
        out: Option<PathBuf>,
    },
}

#[derive(Subcommand, Debug)]
pub enum OpCommand {
    /// Union of 2 to 4,096 tables (summing counts with checked `u64`
    /// arithmetic; capped at `2^32 - 1` only under `--allow-saturation`).
    Union {
        #[arg(long, required = true)]
        input: Vec<PathBuf>,
        /// Output table path (`.ktab`, or `.grk` with `--presence-only`).
        #[arg(long)]
        out: PathBuf,
        /// Sample name written into the output header (defaults to `"sample"`).
        #[arg(long)]
        name: Option<String>,
        /// Write presence-only `.grk` output instead of counted `.ktab`.
        #[arg(long)]
        presence_only: bool,
        /// Clamp summed counts exceeding `2^32 - 1` to `2^32 - 1` and record
        /// saturation in the `.ktab` header and JSON summary instead of
        /// failing.
        #[arg(long)]
        allow_saturation: bool,
    },
    /// Intersection of 2 to 4,096 tables (minimum count across inputs).
    Intersect {
        #[arg(long, required = true)]
        input: Vec<PathBuf>,
        /// Output table path (`.ktab`, or `.grk` with `--presence-only`).
        #[arg(long)]
        out: PathBuf,
        /// Sample name written into the output header (defaults to `"sample"`).
        #[arg(long)]
        name: Option<String>,
        /// Write presence-only `.grk` output instead of counted `.ktab`.
        #[arg(long)]
        presence_only: bool,
    },
    /// Set difference: keys in input 0 absent from all subsequent inputs.
    Diff {
        #[arg(long, required = true)]
        input: Vec<PathBuf>,
        /// Output table path (`.ktab`, or `.grk` with `--presence-only`).
        #[arg(long)]
        out: PathBuf,
        /// Sample name written into the output header (defaults to `"sample"`).
        #[arg(long)]
        name: Option<String>,
        /// Write presence-only `.grk` output instead of counted `.ktab`.
        #[arg(long)]
        presence_only: bool,
    },
    /// Symmetric difference: keys present in an odd number of inputs (counts
    /// summed with checked `u64` arithmetic; capped at `2^32 - 1` only under
    /// `--allow-saturation`).
    Symdiff {
        #[arg(long, required = true)]
        input: Vec<PathBuf>,
        /// Output table path (`.ktab`, or `.grk` with `--presence-only`).
        #[arg(long)]
        out: PathBuf,
        /// Sample name written into the output header (defaults to `"sample"`).
        #[arg(long)]
        name: Option<String>,
        /// Write presence-only `.grk` output instead of counted `.ktab`.
        #[arg(long)]
        presence_only: bool,
        /// Clamp summed counts exceeding `2^32 - 1` to `2^32 - 1` and record
        /// saturation in the `.ktab` header and JSON summary instead of
        /// failing.
        #[arg(long)]
        allow_saturation: bool,
    },
}

#[derive(Clone, Copy, Debug, ValueEnum)]
pub enum MatrixModeArg {
    Pa,
    Counts,
    Dosage,
}

impl From<MatrixModeArg> for MatrixMode {
    fn from(v: MatrixModeArg) -> Self {
        match v {
            MatrixModeArg::Pa => MatrixMode::Pa,
            MatrixModeArg::Counts => MatrixMode::Counts,
            MatrixModeArg::Dosage => MatrixMode::Dosage,
        }
    }
}

#[derive(Clone, Copy, Debug, ValueEnum)]
pub enum InputTypeArg {
    Reads,
    Assembly,
}

impl From<InputTypeArg> for InputType {
    fn from(v: InputTypeArg) -> Self {
        match v {
            InputTypeArg::Reads => InputType::Reads,
            InputTypeArg::Assembly => InputType::Assembly,
        }
    }
}

#[derive(Clone, Copy, Debug, ValueEnum)]
pub enum MatrixLayoutArg {
    Sparse,
    Dense,
}

impl From<MatrixLayoutArg> for MatrixLayout {
    fn from(v: MatrixLayoutArg) -> Self {
        match v {
            MatrixLayoutArg::Sparse => MatrixLayout::Sparse,
            MatrixLayoutArg::Dense => MatrixLayout::Dense,
        }
    }
}

#[derive(Clone, Copy, Debug, ValueEnum)]
pub enum ConvertDirectionArg {
    ToGr,
    ToKtab,
    ToV1,
}

#[derive(Clone, Copy, Debug, ValueEnum)]
pub enum InputModalityArg {
    Unknown,
    Assembly,
    Reads,
}

impl From<InputModalityArg> for InputModality {
    fn from(v: InputModalityArg) -> Self {
        match v {
            InputModalityArg::Unknown => InputModality::Unknown,
            InputModalityArg::Assembly => InputModality::Assembly,
            InputModalityArg::Reads => InputModality::Reads,
        }
    }
}

#[derive(Serialize)]
struct InspectReport {
    format: String,
    file_bytes: u64,
    k: usize,
    sample: Option<String>,
    modality: Option<InputModality>,
    n_kmers: Option<u64>,
    bytes_per_kmer: Option<f64>,
    total_bases: Option<u64>,
    n_sequences: Option<usize>,
    presence_only: Option<bool>,
    max_count_observed: Option<u32>,
    saturated_counts: Option<u64>,
    verified: bool,
    header_sha256: Option<String>,
    fnv64: Option<String>,
    keys_sha256: Option<String>,
    sketch_hex: Option<Vec<String>>,
}

#[derive(Serialize)]
struct ConvertReport {
    format: &'static str,
    k: usize,
    sample: String,
    modality: Option<InputModality>,
    n_kmers: u64,
    output_bytes: u64,
    counts_dropped: bool,
}

/// Writes `val` as pretty-printed JSON atomically via `<path>.part` and rename.
pub fn write_json_atomic<T: Serialize>(path: &Path, val: &T) -> Result<()> {
    let mut bytes = serde_json::to_vec_pretty(val)?;
    bytes.push(b'\n');
    confine::write_atomic(path, &bytes, false)
}

/// Entry point wired from `Command::Kcount` in `src/main.rs`.
pub fn run_cli(args: KcountArgs) -> Result<()> {
    match args.command {
        KcountCommand::Inspect {
            input,
            sketch,
            verify,
            out,
        } => {
            let report = inspect_file(&input, sketch, verify)?;
            emit_json(&report, out.as_deref())
        }
        KcountCommand::ConvertV1 {
            input,
            out,
            direction,
            name,
            modality,
        } => {
            let fmt = detect_format(&input)?;
            let dir = direction.unwrap_or_else(|| infer_direction(fmt, &out));
            let sample_name = name.as_deref().unwrap_or(DEFAULT_SAMPLE_NAME);
            let modality_override = modality.map(InputModality::from);
            let outcome = match dir {
                ConvertDirectionArg::ToGr => {
                    convert_to_gr(&input, &out, sample_name, modality_override)?
                }
                ConvertDirectionArg::ToKtab => {
                    convert_to_ktab(&input, &out, sample_name, modality_override)?
                }
                ConvertDirectionArg::ToV1 => convert_to_v1(&input, &out, sample_name)?,
            };
            let output_bytes = fs::metadata(&out)?.len();
            let report = ConvertReport {
                format: outcome.format.magic_label(),
                k: outcome.k,
                sample: outcome.sample,
                modality: outcome.modality,
                n_kmers: outcome.n_kmers,
                output_bytes,
                counts_dropped: outcome.counts_dropped,
            };
            println!("{}", serde_json::to_string_pretty(&report)?);
            Ok(())
        }
        KcountCommand::Op { action } => {
            let (op, input, out, name, presence_only, allow_saturation) = match action {
                OpCommand::Union {
                    input,
                    out,
                    name,
                    presence_only,
                    allow_saturation,
                } => (
                    SetOp::Union,
                    input,
                    out,
                    name,
                    presence_only,
                    allow_saturation,
                ),
                OpCommand::Intersect {
                    input,
                    out,
                    name,
                    presence_only,
                } => (SetOp::Intersect, input, out, name, presence_only, false),
                OpCommand::Diff {
                    input,
                    out,
                    name,
                    presence_only,
                } => (SetOp::Diff, input, out, name, presence_only, false),
                OpCommand::Symdiff {
                    input,
                    out,
                    name,
                    presence_only,
                    allow_saturation,
                } => (
                    SetOp::SymDiff,
                    input,
                    out,
                    name,
                    presence_only,
                    allow_saturation,
                ),
            };
            let cfg = SetOpConfig {
                sample_name: name.unwrap_or_else(|| DEFAULT_SAMPLE_NAME.to_string()),
                presence_only_out: presence_only,
                allow_saturation,
            };
            let summary = run_set_op(op, &input, &out, &cfg)?;
            println!("{}", serde_json::to_string_pretty(&summary)?);
            Ok(())
        }
        KcountCommand::MatrixExport {
            input,
            out_dir,
            mode,
            layout,
            min_support,
            max_support,
            min_prevalence,
            max_prevalence,
            sample_names,
            core_panel,
            input_type,
            expected_fold_id,
            min_count,
            min_core_depth,
            min_core_breadth,
        } => {
            let names = if sample_names.is_empty() {
                None
            } else {
                Some(sample_names)
            };
            let panel = match &core_panel {
                Some(p) => Some(CorePanel::load(p)?),
                None => None,
            };
            let cfg = MatrixExportConfig {
                mode: mode.into(),
                layout: layout.into(),
                min_support,
                max_support,
                min_prevalence,
                max_prevalence,
                sample_names: names,
                core_panel: panel,
                dosage: DosageConfig {
                    sample_name: DEFAULT_SAMPLE_NAME.to_string(),
                    expected_input_type: input_type.map(InputType::from),
                    expected_fold_id,
                    min_count,
                    min_core_depth,
                    min_core_breadth,
                },
            };
            let manifest = export_matrix(&input, &out_dir, &cfg)?;
            println!("{}", serde_json::to_string_pretty(&manifest)?);
            Ok(())
        }
        KcountCommand::BuildCorePanel {
            input,
            k,
            input_type,
            min_prevalence,
            max_copy_factor,
            max_kmers,
            fold_id,
            out,
        } => {
            let cfg = CorePanelBuildConfig {
                k,
                input_type: input_type.into(),
                min_prevalence,
                max_copy_factor,
                max_kmers,
                fold_id,
            };
            let panel = build_training_fold_core_panel(&input, &cfg)?;
            emit_json(&panel, out.as_deref())
        }
        KcountCommand::Dosage {
            table,
            core_panel,
            name,
            input_type,
            expected_fold_id,
            targets,
            qc_json,
            min_count,
            min_core_depth,
            min_core_breadth,
            out,
        } => {
            let panel = CorePanel::load(&core_panel)?;
            let tp = match &targets {
                Some(p) => Some(TargetPanel::load(p)?),
                None => None,
            };
            let qc = match &qc_json {
                Some(p) => Some(load_qc_metrics(p)?),
                None => None,
            };
            let cfg = DosageConfig {
                sample_name: name.unwrap_or_else(|| DEFAULT_SAMPLE_NAME.to_string()),
                expected_input_type: input_type.map(InputType::from),
                expected_fold_id,
                min_count,
                min_core_depth,
                min_core_breadth,
            };
            let report = extract_sample_dosage(&table, &panel, tp.as_ref(), qc.as_ref(), &cfg)?;
            emit_json(&report, out.as_deref())
        }
    }
}

fn infer_direction(src_format: TableFormat, out: &Path) -> ConvertDirectionArg {
    match src_format {
        TableFormat::KsetV1 => match out.extension().and_then(|s| s.to_str()) {
            Some("ktab") => ConvertDirectionArg::ToKtab,
            _ => ConvertDirectionArg::ToGr,
        },
        TableFormat::GrKset | TableFormat::Ktab => ConvertDirectionArg::ToV1,
    }
}

fn emit_json<T: Serialize>(val: &T, out: Option<&Path>) -> Result<()> {
    if let Some(p) = out {
        write_json_atomic(p, val)
    } else {
        println!("{}", serde_json::to_string_pretty(val)?);
        Ok(())
    }
}

fn inspect_file(path: &Path, sketch: Option<usize>, verify: bool) -> Result<InspectReport> {
    let file_bytes = fs::metadata(path)
        .with_context(|| format!("stat {}", path.display()))?
        .len();

    let mut magic = [0u8; 8];
    let n_read = File::open(path)
        .with_context(|| format!("open {}", path.display()))?
        .read(&mut magic)?;
    if n_read == 8 && &magic == MAGIC_PROF1 {
        if verify {
            let (hdr, profs) = read_profiles_verified(path)?;
            let total_bases = profs
                .iter()
                .try_fold(0u64, |acc, p| acc.checked_add(p.seq_len))
                .context("profile total_bases overflow")?;
            let max_obs = profs
                .iter()
                .flat_map(|p| p.counts.iter().copied())
                .max()
                .unwrap_or(0);
            return Ok(InspectReport {
                format: "profile (GRMPRF1)".to_string(),
                file_bytes,
                k: hdr.k,
                sample: None,
                modality: None,
                n_kmers: None,
                bytes_per_kmer: None,
                total_bases: Some(total_bases),
                n_sequences: Some(profs.len()),
                presence_only: None,
                max_count_observed: Some(max_obs),
                saturated_counts: None,
                verified: true,
                header_sha256: Some(hex_bytes(&hdr.digest_sha256)),
                fnv64: None,
                keys_sha256: None,
                sketch_hex: None,
            });
        }
        let hdr = read_profile_header(path)?;
        let n_seq = usize::try_from(hdr.n_sequences).context("profile n_sequences overflow")?;
        return Ok(InspectReport {
            format: "profile (GRMPRF1)".to_string(),
            file_bytes,
            k: hdr.k,
            sample: None,
            modality: None,
            n_kmers: None,
            bytes_per_kmer: None,
            total_bases: Some(hdr.total_bases),
            n_sequences: Some(n_seq),
            presence_only: None,
            max_count_observed: Some(hdr.max_count_observed),
            saturated_counts: None,
            verified: false,
            header_sha256: Some(hex_bytes(&hdr.digest_sha256)),
            fnv64: None,
            keys_sha256: None,
            sketch_hex: None,
        });
    }

    let reader = AnyTableReader::open(path)?;
    let fmt_str = match reader.format() {
        TableFormat::KsetV1 => "kset-v1 (GRMKSET1)",
        TableFormat::GrKset => "gr-kset (GRMKGR1)",
        TableFormat::Ktab => "ktab (GRMKTB1)",
    };
    let n_kmers = reader.n_kmers();
    let (presence_only, max_obs, saturated_counts) = match &reader {
        AnyTableReader::V1(_) | AnyTableReader::Gr(_) => {
            (Some(true), Some(if n_kmers == 0 { 0 } else { 1 }), None)
        }
        AnyTableReader::Tab(r) => (
            Some(r.presence_only),
            Some(r.max_count_observed),
            Some(r.saturated_counts),
        ),
    };
    let bpk = if n_kmers == 0 {
        0.0
    } else {
        file_bytes as f64 / n_kmers as f64
    };
    let sketch_hex = match sketch {
        Some(s) => {
            if reader.k() <= MAX_K_U64 {
                Some(
                    reader
                        .sketch::<u64>(s)?
                        .into_iter()
                        .map(|v| v.to_hex())
                        .collect(),
                )
            } else {
                Some(
                    reader
                        .sketch::<u128>(s)?
                        .into_iter()
                        .map(|v| v.to_hex())
                        .collect(),
                )
            }
        }
        None => None,
    };
    let (verified, fnv64, keys_sha256) = if verify {
        let digest_matched = reader.verify_integrity()?;
        let (f, s) = if reader.k() <= MAX_K_U64 {
            verify_table_digests::<u64>(&reader)?
        } else {
            verify_table_digests::<u128>(&reader)?
        };
        (digest_matched, Some(f), Some(s))
    } else {
        (false, None, None)
    };
    Ok(InspectReport {
        format: fmt_str.to_string(),
        file_bytes,
        k: reader.k(),
        sample: Some(reader.name().to_string()),
        modality: reader.modality(),
        n_kmers: Some(n_kmers),
        bytes_per_kmer: Some(bpk),
        total_bases: None,
        n_sequences: None,
        presence_only,
        max_count_observed: max_obs,
        saturated_counts,
        verified,
        header_sha256: reader.digest_sha256().map(|d| hex_bytes(&d)),
        fnv64,
        keys_sha256,
        sketch_hex,
    })
}

fn verify_table_digests<K: KmerKey>(reader: &AnyTableReader) -> Result<(String, String)> {
    let mut fnv = Fnv::new();
    let mut sha = Sha256Context::new(&SHA256);
    let mut kbuf = Vec::new();
    let mut cbuf = Vec::new();
    let mut keys = Vec::new();
    let mut counts = Vec::new();
    let mut observed_max = 0u32;
    for b in 0..BUCKETS {
        keys.clear();
        counts.clear();
        reader.read_bucket::<K>(b, &mut kbuf, &mut cbuf, &mut keys, &mut counts)?;
        if keys.len() != counts.len() {
            bail!("bucket {b} key/count length mismatch");
        }
        if let Some(&m) = counts.iter().max() {
            observed_max = observed_max.max(m);
        }
        for &key in &keys {
            key.with_le_bytes(|b| {
                fnv.bytes(b);
                sha.update(b);
            });
        }
    }
    if let AnyTableReader::Tab(r) = reader {
        ensure!(
            observed_max == r.max_count_observed,
            "ktab header max_count_observed ({}) != decoded max ({observed_max})",
            r.max_count_observed
        );
    }
    Ok((fnv.hex(), hex_bytes(sha.finish().as_ref())))
}
