//! Pure-Rust k-mer counting engine, key primitives, and verified binary
//! formats (`grm kcount`).
//!
//! Submodules:
//! - [`key`]: dual-width canonical k-mer spaces (`u64` for `k = 5..=32`,
//!   `u128` for `k = 33..=64`) with branchless O(1) SWAR reverse complement
//!   and invertible `2k`-bit mixers.
//! - [`mod@format`]: verified Golomb-Rice sorted-key table (`.grk`,
//!   `GRMKGR1\0`), counted table (`.ktab`, `GRMKTB1\0`), per-sequence count
//!   profile (`.prof.bin`, `GRMPRF1\0`), and `.kset` v1 (`GRMKSET1`) unified
//!   reader.
//! - [`superkmer`]: bounded-memory parallel minimizer-bucketed super-k-mer
//!   counter (`grm kcount count`) with `u32` counts capped at `2^32 - 1`,
//!   exact histograms and QC metrics, and indexed parallel `.prof.bin` writing.

pub mod format;
pub mod key;
pub mod superkmer;

use crate::confine;
use crate::fnv::Fnv;
use crate::kcount::format::{
    AnyTableReader, DEFAULT_SAMPLE_NAME, InputModality, MAGIC_PROF1, TableFormat, convert_to_gr,
    convert_to_ktab, convert_to_v1, detect_format, hex_bytes, read_profile_header,
    read_profiles_verified,
};
use crate::kcount::key::{KmerKey, MAX_K, MAX_K_U64, MIN_K};
use crate::kcount::superkmer::{CountConfig, InputTypeChoice, count_inputs};
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
    /// Count canonical k-mers (`k = 5..=64`) from plain or `.gz` FASTA/FASTQ
    /// using minimizer-bucketed 2-bit super-k-mers and `u32` counts capped at
    /// `2^32 - 1`.
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
        /// Sample name stored in output headers (defaults to `"sample"` so
        /// input file names are never copied into output headers).
        #[arg(long)]
        name: Option<String>,
        /// Input modality (`auto`, `reads`, or `assembly`). `auto` detects
        /// FASTQ as `reads` and FASTA as `assembly`; an explicit value errors
        /// when the input records disagree with it.
        #[arg(long, value_enum, default_value_t = InputTypeArg::Auto)]
        input_type: InputTypeArg,
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
        /// Optional per-sequence k-mer count profile (`.prof.bin`).
        #[arg(long)]
        profile_out: Option<PathBuf>,
        /// Retain original FASTA/FASTQ record IDs in `.prof.bin` instead of
        /// anonymized sequential labels (`r000001`, `r000002`, ...).
        #[arg(long)]
        keep_record_ids: bool,
        /// Minimum exact k-mer count to retain (`1` for assemblies, `>= 3` for raw reads).
        #[arg(long, default_value_t = 1)]
        min_count: u32,
        /// Optional upper count ceiling (`<= 2^32 - 1`).
        #[arg(long)]
        max_count: Option<u32>,
        /// Allow counts exceeding `2^32 - 1` (`4,294,967,295`) to clamp at
        /// `u32::MAX` and record `saturated_counts` instead of failing.
        #[arg(long)]
        allow_saturation: bool,
        /// Phred base-quality floor for FASTQ inputs; bases below break windows like `N`.
        #[arg(long, default_value_t = 0)]
        min_qual: u8,
        /// Bottom-`s` MinHash sketch size.
        #[arg(long, default_value_t = 1000)]
        sketch_size: usize,
        /// Target working-memory budget in MiB for bounded partitioning, merge,
        /// and spilling.
        #[arg(long, default_value_t = 512)]
        max_memory_mb: usize,
        /// Power-of-2 minimizer partition count (`16..=1024`).
        #[arg(long, default_value_t = 256)]
        partitions: usize,
        /// Optional minimizer length override (`1 <= m <= min(k, 31)`).
        #[arg(long)]
        minimizer_len: Option<usize>,
        /// Write presence-only Golomb-Rice `.grk` table instead of counted `.ktab`.
        #[arg(long)]
        presence_only: bool,
        /// Optional temporary directory for spilled partitions and shards.
        #[arg(long)]
        tmp_dir: Option<PathBuf>,
    },
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
}

#[derive(Clone, Copy, Debug, ValueEnum)]
pub enum InputTypeArg {
    Auto,
    Reads,
    Assembly,
}

impl From<InputTypeArg> for InputTypeChoice {
    fn from(v: InputTypeArg) -> Self {
        match v {
            InputTypeArg::Auto => InputTypeChoice::Auto,
            InputTypeArg::Reads => InputTypeChoice::Reads,
            InputTypeArg::Assembly => InputTypeChoice::Assembly,
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
        KcountCommand::Count {
            input,
            k,
            name,
            input_type,
            out,
            gr_kset_out,
            hist_out,
            meta_out,
            profile_out,
            keep_record_ids,
            min_count,
            max_count,
            allow_saturation,
            min_qual,
            sketch_size,
            max_memory_mb,
            partitions,
            minimizer_len,
            presence_only,
            tmp_dir,
        } => {
            let sample_name = name.unwrap_or_else(|| DEFAULT_SAMPLE_NAME.to_string());
            let cfg = CountConfig {
                k,
                sample_name,
                input_type: input_type.into(),
                min_count,
                max_count,
                allow_saturation,
                min_qual,
                sketch_size,
                max_memory_mb,
                partitions,
                minimizer_len,
                presence_only_out: presence_only,
                keep_record_ids,
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
