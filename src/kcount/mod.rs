//! Next-generation pure-Rust k-mer counting engine (`grm kcount`).
//!
//! Submodules:
//! - [`key`]: dual-width canonical k-mer spaces (`u64` for `k = 5..=32`,
//!   `u128` for `k = 33..=64`) with branchless O(1) SWAR reverse complement
//!   and invertible `2k`-bit mixers.
//! - [`superkmer`]: minimizer-bucketed 2-bit super-k-mer streaming counter
//!   with bounded-memory spilling, intra-genome Rayon parallelism, exact
//!   uncapped `u32` counts, single-pass histogram/digest/sketch emission, and
//!   optional FastK-style per-sequence count profiles.
//! - [`mod@format`]: compact counted `.ktab` (`GRMKTB1\0`), Golomb-Rice sorted-key
//!   baseline (`GRMKGR1\0`), profile (`GRMPRF1\0`), and `.kset` v1 (`GRMKSET1`)
//!   unified reader.
//! - [`spss`]: experimental 2-bit Spectrum-Preserving String Set greedy
//!   simplitig archive (`spss-archive`, `GRMSPS1\0`) with lossless round-trip
//!   to `.kset` v1 (`GRMKSET1`).
//! - [`dosage`]: fold-safe `D_core` normalization, plasmid dosage extractor,
//!   explicit missingness gate, and per-sample QC vector.
//! - [`preflight`]: `ntHash`-style single-pass HyperLogLog + adaptive
//!   low-abundance preflight estimator.
//! - [`ops`]: streaming bucket-by-bucket set algebra (`union`, `intersect`,
//!   `diff`, `symdiff`) and ML-ready `matrix-export`.

pub mod dosage;
pub mod format;
pub mod key;
pub mod ops;
pub mod preflight;
pub mod spss;
pub mod superkmer;

use crate::fnv::Fnv;
use crate::kcount::dosage::{
    CorePanel, DosageConfig, InputType, TargetPanel, extract_sample_dosage, load_qc_metrics,
};
use crate::kcount::format::{
    AnyTableReader, MAGIC_PROF1, TableFormat, convert_to_v1, convert_v1_to_gr, read_profiles,
};
use crate::kcount::key::{KmerKey, MAX_K, MAX_K_U64, MIN_K};
use crate::kcount::ops::{
    MatrixExportConfig, MatrixFormat, MatrixMode, SetOp, export_matrix, run_set_op,
};
use crate::kcount::preflight::preflight_inputs;
use crate::kcount::spss::{
    archive_auto, expand_to_kset_v1, expand_to_ktab, is_spss_file, read_spss_header, verify_spss,
};
use crate::kcount::superkmer::{CountConfig, count_inputs, hex_bytes, write_json_atomic};
use crate::kset::BUCKETS;
use anyhow::{Context, Result, bail};
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
    /// Count canonical k-mers (`k = 5..=64`) from plain or `.gz` FASTA/FASTQ
    /// using minimizer-bucketed 2-bit super-k-mers and uncapped `u32` counts.
    Count {
        /// One or more FASTA/FASTQ files (plain or `.gz`).
        #[arg(long, required = true)]
        input: Vec<PathBuf>,
        /// K-mer length (`5..=64`).
        #[arg(
            long,
            default_value_t = 31,
            value_parser = RangedU64ValueParser::<usize>::new().range(MIN_K as u64..=MAX_K as u64)
        )]
        k: usize,
        /// Sample name stored in output headers (defaults to first input stem).
        #[arg(long)]
        name: Option<String>,
        /// Primary output table path (`.ktab`, or `.grk` with `--presence-only`).
        #[arg(long)]
        out: PathBuf,
        /// Optional Golomb-Rice presence-only sidecar output (`.grk`).
        #[arg(long)]
        gr_kset_out: Option<PathBuf>,
        /// Optional abundance histogram JSON output (`.hist.json`).
        #[arg(long)]
        hist_out: Option<PathBuf>,
        /// Optional full summary, digests, sketch, and QC JSON output.
        #[arg(long)]
        meta_out: Option<PathBuf>,
        /// Optional FastK-style per-sequence k-mer count profile (`.prof.bin`).
        #[arg(long)]
        profile_out: Option<PathBuf>,
        /// Minimum exact k-mer count to retain (`1` for assemblies, `>= 3` for raw reads).
        #[arg(long, default_value_t = 1)]
        min_count: u32,
        /// Optional upper count ceiling.
        #[arg(long)]
        max_count: Option<u32>,
        /// Phred base-quality floor for FASTQ inputs; bases below break windows like `N`.
        #[arg(long, default_value_t = 0)]
        min_qual: u8,
        /// Bottom-`s` MinHash sketch size.
        #[arg(long, default_value_t = 1000)]
        sketch_size: usize,
        /// Memory budget in MiB before spilling super-k-mer partitions to disk.
        #[arg(long, default_value_t = 512)]
        max_memory_mb: usize,
        /// Power-of-2 partition count (`16..=1024`).
        #[arg(long, default_value_t = 256)]
        partitions: usize,
        /// Optional minimizer length override (`m <= k`).
        #[arg(long)]
        minimizer_len: Option<usize>,
        /// Write presence-only Golomb-Rice table instead of counted `.ktab`.
        #[arg(long)]
        presence_only: bool,
        /// Optional temporary directory for spilled super-k-mer partitions.
        #[arg(long)]
        tmp_dir: Option<PathBuf>,
    },
    /// Pack a `.kset` v1, `.ktab`, `.grk`, or FASTA/FASTQ into a 2-bit
    /// Spectrum-Preserving String Set (`spss-archive`, `GRMSPS1\0`).
    Archive {
        /// Input `.kset` v1, `.ktab`, `.grk`, or FASTA/FASTQ file.
        #[arg(long)]
        input: PathBuf,
        /// Output `spss-archive` path (`.spss`).
        #[arg(long)]
        out: PathBuf,
        /// K-mer length when archiving directly from FASTA/FASTQ (`5..=64`).
        #[arg(
            long,
            value_parser = RangedU64ValueParser::<usize>::new().range(MIN_K as u64..=MAX_K as u64)
        )]
        k: Option<usize>,
        /// Sample name when archiving directly from FASTA/FASTQ.
        #[arg(long)]
        name: Option<String>,
    },
    /// Expand an `spss-archive` (`GRMSPS1\0`) losslessly back to `.kset` v1
    /// (`GRMKSET1`, `k <= 32`) or presence-only `.ktab` (`GRMKTB1\0`),
    /// verifying its `fnv64` digest.
    Expand {
        /// Input `spss-archive` path (`.spss`).
        #[arg(long)]
        input: PathBuf,
        /// Output `.kset` or `.ktab` path.
        #[arg(long)]
        out: PathBuf,
        /// Target format (`kset-v1` default for `k <= 32`, or `ktab` for any `k`).
        #[arg(long, value_enum, default_value_t = ExpandTargetArg::KsetV1)]
        format: ExpandTargetArg,
    },
    /// Single-pass streaming preflight estimating `F0`, `f1`, `f2+`, and disk
    /// sizes without sorting.
    Preflight {
        /// One or more FASTA/FASTQ files (plain or `.gz`).
        #[arg(long, required = true)]
        input: Vec<PathBuf>,
        /// K-mer length (`5..=64`).
        #[arg(
            long,
            default_value_t = 31,
            value_parser = RangedU64ValueParser::<usize>::new().range(MIN_K as u64..=MAX_K as u64)
        )]
        k: usize,
        /// Phred base-quality floor for FASTQ inputs.
        #[arg(long, default_value_t = 0)]
        min_qual: u8,
        /// Optional JSON output path (prints to stdout when omitted).
        #[arg(long)]
        out: Option<PathBuf>,
    },
    /// Fold-safe `D_core` depth normalization and plasmid/target dosage
    /// extraction with explicit missingness on low core coverage.
    Dosage {
        /// Counted `.ktab` (`GRMKTB1\0`) file for the sample.
        #[arg(long)]
        table: PathBuf,
        /// Single-copy core k-mer panel JSON (pinned or training-fold derived).
        #[arg(long)]
        core_panel: PathBuf,
        /// Input modality (`reads` or `assembly`), validated against the core panel.
        #[arg(long, value_enum, default_value_t = InputTypeArg::Reads)]
        input_type: InputTypeArg,
        /// Optional target/plasmid k-mer panel JSON.
        #[arg(long)]
        targets: Option<PathBuf>,
        /// Optional metadata or QC JSON emitted by `grm kcount count`.
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
    /// Streaming bucket-by-bucket set algebra and matrix export across 2+
    /// `.kset` v1, `.ktab`, or `.grk` tables.
    Op {
        #[command(subcommand)]
        action: OpCommand,
    },
    /// Inspect metadata, bucket statistics, sketch prefixes, or full verified
    /// digests for `.kset` v1, `.ktab`, `.grk`, `.spss`, or `.prof.bin` files.
    Inspect {
        /// File to inspect.
        #[arg(long)]
        input: PathBuf,
        /// Optional bottom-`s` MinHash sketch prefix length to include.
        #[arg(long)]
        sketch: Option<usize>,
        /// Decode every bucket or simplitig and verify exact digests.
        #[arg(long)]
        verify: bool,
        /// Optional JSON output path (prints to stdout when omitted).
        #[arg(long)]
        out: Option<PathBuf>,
    },
    /// Lossless conversion between `.kset` v1 (`GRMKSET1`) and `.grk` / `.ktab`
    /// for `k <= 32`.
    ConvertV1 {
        /// Source file (`.kset` v1, `.grk`, or `.ktab`).
        #[arg(long)]
        input: PathBuf,
        /// Output file path.
        #[arg(long)]
        out: PathBuf,
        /// Target direction (`to-gr` or `to-v1`; auto-detects opposite when omitted).
        #[arg(long, value_enum)]
        direction: Option<ConvertDirectionArg>,
    },
}

#[derive(Subcommand, Debug)]
pub enum OpCommand {
    /// Union of 2+ tables (summing counts).
    Union {
        #[arg(long, required = true)]
        input: Vec<PathBuf>,
        #[arg(long)]
        out: PathBuf,
        #[arg(long, default_value = "union")]
        name: String,
        #[arg(long)]
        presence_only: bool,
    },
    /// Intersection of 2+ tables (taking minimum count across inputs).
    Intersect {
        #[arg(long, required = true)]
        input: Vec<PathBuf>,
        #[arg(long)]
        out: PathBuf,
        #[arg(long, default_value = "intersect")]
        name: String,
        #[arg(long)]
        presence_only: bool,
    },
    /// Difference: keys in first input absent from all subsequent inputs.
    Diff {
        #[arg(long, required = true)]
        input: Vec<PathBuf>,
        #[arg(long)]
        out: PathBuf,
        #[arg(long, default_value = "diff")]
        name: String,
        #[arg(long)]
        presence_only: bool,
    },
    /// Symmetric difference: keys present in an odd number of inputs.
    Symdiff {
        #[arg(long, required = true)]
        input: Vec<PathBuf>,
        #[arg(long)]
        out: PathBuf,
        #[arg(long, default_value = "symdiff")]
        name: String,
        #[arg(long)]
        presence_only: bool,
    },
    /// Stream a deterministic sparse or dense ML feature matrix (`pa`, `counts`,
    /// or `dosage`) with explicit missingness mask.
    MatrixExport {
        #[arg(long, required = true)]
        input: Vec<PathBuf>,
        #[arg(long)]
        out_dir: PathBuf,
        #[arg(long, value_enum, default_value_t = MatrixModeArg::Pa)]
        mode: MatrixModeArg,
        #[arg(long, value_enum, default_value_t = MatrixFormatArg::Sparse)]
        format: MatrixFormatArg,
        #[arg(long, default_value_t = 1)]
        min_support: usize,
        #[arg(long)]
        max_support: Option<usize>,
        #[arg(long)]
        core_panel: Option<PathBuf>,
        #[arg(long, value_enum, default_value_t = InputTypeArg::Reads)]
        input_type: InputTypeArg,
        #[arg(long, default_value_t = 1)]
        min_count: u32,
        #[arg(long, default_value_t = 10.0)]
        min_core_depth: f64,
        #[arg(long, default_value_t = 0.80)]
        min_core_breadth: f64,
    },
}

#[derive(Clone, Copy, Debug, ValueEnum)]
pub enum ExpandTargetArg {
    KsetV1,
    Ktab,
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
pub enum MatrixFormatArg {
    Sparse,
    Dense,
}

impl From<MatrixFormatArg> for MatrixFormat {
    fn from(v: MatrixFormatArg) -> Self {
        match v {
            MatrixFormatArg::Sparse => MatrixFormat::Sparse,
            MatrixFormatArg::Dense => MatrixFormat::Dense,
        }
    }
}

#[derive(Clone, Copy, Debug, ValueEnum)]
pub enum ConvertDirectionArg {
    ToGr,
    ToV1,
}

#[derive(Serialize)]
struct InspectReport {
    path: String,
    format: String,
    file_bytes: u64,
    k: usize,
    sample: Option<String>,
    n_kmers: Option<u64>,
    bytes_per_kmer: Option<f64>,
    n_simplitigs: Option<u64>,
    total_bases: Option<u64>,
    n_sequences: Option<usize>,
    presence_only: Option<bool>,
    max_count_observed: Option<u32>,
    verified: bool,
    fnv64: Option<String>,
    keys_sha256: Option<String>,
    sketch_hex: Option<Vec<String>>,
}

/// Entry point wired from `Command::Kcount` in `src/main.rs`.
pub fn run_cli(args: KcountArgs) -> Result<()> {
    match args.command {
        KcountCommand::Count {
            input,
            k,
            name,
            out,
            gr_kset_out,
            hist_out,
            meta_out,
            profile_out,
            min_count,
            max_count,
            min_qual,
            sketch_size,
            max_memory_mb,
            partitions,
            minimizer_len,
            presence_only,
            tmp_dir,
        } => {
            let sample_name = name.unwrap_or_else(|| {
                input[0]
                    .file_stem()
                    .and_then(|s| s.to_str())
                    .unwrap_or("sample")
                    .to_string()
            });
            let cfg = CountConfig {
                k,
                sample_name,
                min_count,
                max_count,
                min_qual,
                sketch_size,
                max_memory_mb,
                partitions,
                minimizer_len,
                presence_only_out: presence_only,
                out_path: out,
                gr_kset_out,
                hist_out,
                meta_out,
                profile_out,
                tmp_dir,
            };
            let summary = count_inputs(&input, &cfg)?;
            println!("{}", serde_json::to_string_pretty(&summary)?);
            Ok(())
        }
        KcountCommand::Archive {
            input,
            out,
            k,
            name,
        } => {
            let summary = archive_auto(&input, &out, k, name.as_deref())?;
            println!("{}", serde_json::to_string_pretty(&summary)?);
            Ok(())
        }
        KcountCommand::Expand { input, out, format } => {
            let summary = match format {
                ExpandTargetArg::KsetV1 => expand_to_kset_v1(&input, &out)?,
                ExpandTargetArg::Ktab => expand_to_ktab(&input, &out)?,
            };
            println!("{}", serde_json::to_string_pretty(&summary)?);
            Ok(())
        }
        KcountCommand::Preflight {
            input,
            k,
            min_qual,
            out,
        } => {
            let report = preflight_inputs(&input, k, min_qual)?;
            emit_json(&report, out.as_deref())
        }
        KcountCommand::Dosage {
            table,
            core_panel,
            input_type,
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
                input_type: input_type.into(),
                min_count,
                min_core_depth,
                min_core_breadth,
            };
            let report = extract_sample_dosage(&table, &panel, tp.as_ref(), qc.as_ref(), &cfg)?;
            emit_json(&report, out.as_deref())
        }
        KcountCommand::Op { action } => match action {
            OpCommand::Union {
                input,
                out,
                name,
                presence_only,
            } => {
                let s = run_set_op(SetOp::Union, &input, &out, &name, presence_only)?;
                println!("{}", serde_json::to_string_pretty(&s)?);
                Ok(())
            }
            OpCommand::Intersect {
                input,
                out,
                name,
                presence_only,
            } => {
                let s = run_set_op(SetOp::Intersect, &input, &out, &name, presence_only)?;
                println!("{}", serde_json::to_string_pretty(&s)?);
                Ok(())
            }
            OpCommand::Diff {
                input,
                out,
                name,
                presence_only,
            } => {
                let s = run_set_op(SetOp::Diff, &input, &out, &name, presence_only)?;
                println!("{}", serde_json::to_string_pretty(&s)?);
                Ok(())
            }
            OpCommand::Symdiff {
                input,
                out,
                name,
                presence_only,
            } => {
                let s = run_set_op(SetOp::SymDiff, &input, &out, &name, presence_only)?;
                println!("{}", serde_json::to_string_pretty(&s)?);
                Ok(())
            }
            OpCommand::MatrixExport {
                input,
                out_dir,
                mode,
                format,
                min_support,
                max_support,
                core_panel,
                input_type,
                min_count,
                min_core_depth,
                min_core_breadth,
            } => {
                let panel = match &core_panel {
                    Some(p) => Some(CorePanel::load(p)?),
                    None => None,
                };
                let cfg = MatrixExportConfig {
                    mode: mode.into(),
                    format: format.into(),
                    min_support,
                    max_support,
                    core_panel: panel,
                    dosage: DosageConfig {
                        input_type: input_type.into(),
                        min_count,
                        min_core_depth,
                        min_core_breadth,
                    },
                };
                let manifest = export_matrix(&input, &out_dir, &cfg)?;
                println!("{}", serde_json::to_string_pretty(&manifest)?);
                Ok(())
            }
        },
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
        } => {
            let fmt = crate::kcount::format::detect_format(&input)?;
            let dir = direction.unwrap_or(match fmt {
                TableFormat::KsetV1 => ConvertDirectionArg::ToGr,
                TableFormat::GrKset | TableFormat::Ktab => ConvertDirectionArg::ToV1,
            });
            let n = match dir {
                ConvertDirectionArg::ToGr => convert_v1_to_gr(&input, &out)?,
                ConvertDirectionArg::ToV1 => convert_to_v1(&input, &out)?,
            };
            let out_bytes = fs::metadata(&out)?.len();
            println!(
                "{}",
                serde_json::to_string_pretty(&serde_json::json!({
                    "input": input.display().to_string(),
                    "output": out.display().to_string(),
                    "n_kmers": n,
                    "output_bytes": out_bytes,
                }))?
            );
            Ok(())
        }
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
    if is_spss_file(path)? {
        let h = if verify {
            verify_spss(path)?
        } else {
            read_spss_header(path)?
        };
        let bpk = if h.n_kmers == 0 {
            0.0
        } else {
            file_bytes as f64 / h.n_kmers as f64
        };
        return Ok(InspectReport {
            path: path.display().to_string(),
            format: "spss-archive (GRMSPS1)".to_string(),
            file_bytes,
            k: h.k,
            sample: Some(h.name),
            n_kmers: Some(h.n_kmers),
            bytes_per_kmer: Some(bpk),
            n_simplitigs: Some(h.n_simplitigs),
            total_bases: Some(h.total_bases),
            n_sequences: None,
            presence_only: Some(true),
            max_count_observed: None,
            verified: verify,
            fnv64: Some(format!("{:016x}", h.fnv64_digest)),
            keys_sha256: None,
            sketch_hex: None,
        });
    }

    let mut magic = [0u8; 8];
    let n_read = File::open(path)?.read(&mut magic)?;
    if n_read == 8 && &magic == MAGIC_PROF1 {
        let (k, profs) = read_profiles(path)?;
        return Ok(InspectReport {
            path: path.display().to_string(),
            format: "profile (GRMPRF1)".to_string(),
            file_bytes,
            k,
            sample: None,
            n_kmers: None,
            bytes_per_kmer: None,
            n_simplitigs: None,
            total_bases: Some(profs.iter().map(|p| p.seq_len).sum()),
            n_sequences: Some(profs.len()),
            presence_only: None,
            max_count_observed: profs.iter().flat_map(|p| p.counts.iter().copied()).max(),
            verified: true,
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
    let (presence_only, max_obs) = match &reader {
        AnyTableReader::V1(_) | AnyTableReader::Gr(_) => (Some(true), Some(1)),
        AnyTableReader::Tab(r) => (Some(r.presence_only), Some(r.max_count_observed)),
    };
    let n_kmers = reader.n_kmers();
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
    let (fnv64, keys_sha256) = if verify {
        let (f, s) = if reader.k() <= MAX_K_U64 {
            verify_table_digests::<u64>(&reader)?
        } else {
            verify_table_digests::<u128>(&reader)?
        };
        (Some(f), Some(s))
    } else {
        (None, None)
    };
    Ok(InspectReport {
        path: path.display().to_string(),
        format: fmt_str.to_string(),
        file_bytes,
        k: reader.k(),
        sample: Some(reader.name().to_string()),
        n_kmers: Some(n_kmers),
        bytes_per_kmer: Some(bpk),
        n_simplitigs: None,
        total_bases: None,
        n_sequences: None,
        presence_only,
        max_count_observed: max_obs,
        verified: verify,
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
    for b in 0..BUCKETS {
        keys.clear();
        counts.clear();
        reader.read_bucket::<K>(b, &mut kbuf, &mut cbuf, &mut keys, &mut counts)?;
        if keys.len() != counts.len() {
            bail!("bucket {b} key/count length mismatch");
        }
        for &key in &keys {
            key.with_le_bytes(|b| {
                fnv.bytes(b);
                sha.update(b);
            });
        }
    }
    Ok((fnv.hex(), hex_bytes(sha.finish().as_ref())))
}
