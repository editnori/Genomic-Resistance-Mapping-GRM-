# Checked 2-bit SPSS simplitig archive and streaming preflight (`grm kcount`)

This page documents the checked 2-bit Spectrum-Preserving String Set archive format (`src/kcount/spss.rs`), the single-pass rolling-hash HyperLogLog preflight estimator (`src/kcount/preflight.rs`), and the `grm kcount archive`, `grm kcount expand`, and `grm kcount preflight` CLI subcommands. Existing `grm count`, `grm matrix`, `.kset` v1 (`GRMKSET1`), and Revision 14 learners are unchanged.

## 1. Checked 2-bit SPSS simplitig archive (`src/kcount/spss.rs`)

At bacterial density near 5 million distinct 31-mers, 5.15 bytes is the floor for a code that uses only the sorted values. A presence-only k-mer set still carries `(k - 1)`-base overlaps between adjacent k-mers, so greedy ProphAsm-style simplitigs packed at 2 bits per base reduce archival storage well below the sorted-value floor while keeping `.kset` v1, `.grk`, and `.ktab` untouched for random access and multi-genome merges.

### Binary layout (`GRMSPS1\0`, 80-byte fixed header)

- `[0..8]`: magic `GRMSPS1\0`.
- `[8]`: `k` (`5..=64`).
- `[9]`: 1-byte `InputModality` (`0 = unknown`, `1 = assembly`, `2 = reads`), preserved across `.grk` and `.ktab` archiving and expansion.
- `[10..12]`: reserved bytes (must be zero).
- `[12..16]`: `name_len` (`u32`, little-endian).
- `[16..24]`: `n_kmers` (`u64`, little-endian).
- `[24..32]`: `n_simplitigs` (`u64`, little-endian).
- `[32..40]`: `total_bases` (`u64`, little-endian), required to equal `n_kmers + (k - 1) * n_simplitigs`.
- `[40..48]`: `fnv64_digest` (`u64`, little-endian) over the sorted mixed keys.
- `[48..80]`: 32-byte SHA-256 digest covering `[0..48]`, the UTF-8 sample name, the 8-byte `lengths_byte_len`, the ULEB128 simplitig excess-length stream (`slen - k` per simplitig), and the LSB-first 2-bit packed base stream (`ceil(total_bases / 4)` bytes, unused high padding bits required to be zero).

### Integrity, bounds checks, and count-loss guard

- **Checked header and payload bounds:** Both `read_spss_header` (`grm kcount inspect`) and `read_spss_body` (`grm kcount expand` and `inspect --verify`) validate every length and count (`name_len`, `n_kmers`, `n_simplitigs`, `total_bases`, `lengths_byte_len`, and packed byte length) with checked arithmetic against the remaining unread bytes of the file and allocate with fallible `try_reserve` before touching heap memory.
- **Full verification:** Every `expand` and `inspect --verify` recomputes the 32-byte SHA-256 digest, checks zero trailing padding bits, verifies strict uniqueness of decoded keys across all 1,024 buckets, and checks both `n_kmers` and `fnv64_digest`. `inspect` without `--verify` reports `verified: false`; `inspect --verify` sets `verified: true` only after matching the stored SHA-256 and `fnv64` digests.
- **Counted `.ktab` guard:** Because `.spss` archives presence-only k-mer sets, `grm kcount archive` refuses any counted `.ktab` (non-presence-only or any count above 1) unless `--drop-counts` is passed explicitly. Expanding `.spss` to `.ktab` writes a presence-only `.ktab` with every count equal to `1`. Expanding `.spss` to `.kset` v1 (`k <= 32`) or `.grk` (`k = 5..=64`) reproduces the exact sorted k-mer set, byte for byte when sample name and modality match.
- **Default sample name:** When `--name` is omitted on `grm kcount archive`, the archive header records `"sample"`. Input file paths and filename stems are never copied into `.spss` headers or JSON output.

### Locally measured synthetic archive sizes (`k = 31`)

Measured locally on synthetic inputs built in Rust integration tests (Ubuntu 24.04):

- Deterministic 500-contig synthetic assembly (`tests/kcount_archive.rs`, 439,226 distinct 31-mers with repeats and N runs): `.kset` v1 takes 2,935,132 bytes (6.6825 B/k-mer), `.grk` takes 2,471,268 bytes (5.6264 B/k-mer), and `.spss` takes 115,076 bytes (550 simplitigs, 455,726 total bases, 0.2620 B/k-mer).
- Synthetic 1%-error FASTQ read sets (`tests/kcount_preflight.rs`, 25 kb genome at ~20x coverage, 150 bp reads, 1% substitution errors): `.spss` archives above 0.40 B/k-mer because sequencing errors create short isolated singleton simplitigs, compared to under 0.27 B/k-mer on clean multi-contig assemblies.

## 2. Single-pass rolling-hash preflight estimator (`src/kcount/preflight.rs`)

`grm kcount preflight` streams plain or `.gz` FASTA/FASTQ inputs in a single pass without sorting to estimate distinct k-mers (`F0`), singletons (`f1`), solid k-mers (`f2+`), peak coverage, GC fraction, and expected `.kset` v1, `.grk`, and `.spss` file sizes.

- **Rolling canonical hash and seed provenance:** `NtHash64` computes a strand-invariant O(1) rolling 64-bit canonical hash for any `k = 5..=64`, breaking windows on non-ACGT bytes or FASTQ bases strictly below `--min-qual` (Phred+33; bases at exactly `min_qual` are retained). The four odd 64-bit base seeds (`A=0, C=1, G=2, T=3`) are derived deterministically from `SHA-256("GRM-KCOUNT-NTHASH-SEED-V1:0")` through `":3"` (first 8 digest bytes interpreted little-endian with the low bit set to 1) and verified in unit and integration tests using `ring::digest::SHA256`. The 64-bit output finalizer uses the two public-domain SplitMix64 multipliers (`0xbf58_476d_1ce4_e5b9` and `0x94d0_49bb_1331_11eb`, Steele, Lea, and Flood, OOPSLA 2014).
- **Sketch table memory (distinct from process RSS):** `HyperLogLog14` uses 16,384 1-byte registers (`p = 14`, 16 KiB) with LinearCounting in the small-cardinality regime (`raw <= 2.5 * M`). `AdaptiveFreqSketch` caps its hash-prefix subsampled `FxHashMap<u64, u32>` at 65,536 entries (~2.1 MiB table capacity), returning exact `F0`, `f1`, and `f2+` when `sample_exact == true` (`F0 <= 65,536`) and scaling by `2^subsample_shift` above that threshold. Total process RSS is larger than these sketch tables because it also includes binary memory, allocator overhead, the 1 MiB `BufReader` input buffer, and reusable per-record line buffers.
- **Modality-calibrated `.spss` size estimate:** Rather than applying one fixed constant per k-mer, preflight selects an `"assembly"` or `"reads"` profile (auto-detected from FASTQ presence, mean k-mer multiplicity, and mean record length, or set explicitly with `--modality`). For assemblies, it models `total_records` simplitigs at `(F0 + (k - 1) * simplitigs) / 4` bytes plus length and header overhead; for reads, it scales per-k-mer bytes with the estimated singleton fraction to account for sequencing-error tips.
- **Default sample label:** When `--name` is omitted on `grm kcount preflight`, the JSON report sets `sample` to `"sample"` and omits input file paths and filename stems.

## 3. CLI commands

```text
grm kcount archive   --input <TABLE|FASTA/FASTQ[.gz]> --out <PATH.spss> [--k <5..=64>] [--name <NAME>] [--modality unknown|assembly|reads] [--drop-counts]
grm kcount expand    --input <PATH.spss> --out <PATH.kset|PATH.grk|PATH.ktab> [--format kset-v1|grk|ktab]
grm kcount preflight --input <FASTA/FASTQ[.gz]>... [--k <5..=64>] [--min-qual <Q>] [--name <NAME>] [--modality unknown|assembly|reads] [--out <PATH>]
grm kcount inspect   --input <FILE> [--sketch <S>] [--verify] [--out <PATH>]
```

- `grm kcount expand` enforces that `--out`'s file extension matches the target format (`.kset` for `kset-v1`, `.grk` for `grk`, `.ktab` for `ktab`) and infers `--format` from `--out` when omitted.
