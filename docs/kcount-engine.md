# Next-Generation K-mer Counter Engine (`grm kcount`)

`grm kcount` (`src/kcount/`) is an additive, pure-Rust k-mer counting, archiving, preflight estimation, fold-safe dosage normalization, and streaming set-algebra engine. Existing `grm count`, `grm matrix`, `.kset` v1 (`GRMKSET1`), and Revision 14 learners remain untouched.

## Architecture

### 1. Dual-Width Canonical Keys (`src/kcount/key.rs`)

- **Fast `u64` path (`k = 5..=32`):** `KmerSpace64` uses the exact 2-bit encoding (`A=0, C=1, G=2, T=3`) and SplitMix64-derived bijective `2k`-bit mixer of `crate::dna::KmerSpace`, preserving 100% bit-for-bit parity with `.kset` v1 digests and bottom-`s` sketches.
- **Wide `u128` path (`k = 33..=64`):** `KmerSpace128` packs `2k` bits into a `u128` with a 3-round odd-multiply and right-xorshift invertible bijection across `[0, 2^(2k))`, supporting `k = 61` mobile-element context resolution without truncation.
- **Branchless O(1) SWAR reverse complement:** Both `u64` and `u128` reverse complement via 2-bit pair swaps, nibble swaps, `swap_bytes`, bitwise complement, and logical right shift.
- **Rolling canonical scanner:** Resets cleanly on any non-ACGT byte (`N`, IUPAC codes, whitespace) and optionally on FASTQ bases below `--min-qual` (Phred+33).

### 2. Minimizer-Bucketed Super-k-mer Counter (`src/kcount/superkmer.rs`)

- **Streaming input:** Zero per-read heap allocations on plain and `.gz` FASTA and FASTQ inputs after buffer warmup, computing read/contig QC metrics (`total_reads`, `usable_read_len_median`, `gc_fraction`, `singleton_fraction`, `n50`) in the same pass.
- **Super-k-mer partitioning:** Contiguous valid ACGT runs compute canonical `m`-mer minimizers via an O(1) amortized monotonic deque and pack consecutive k-mers sharing a minimizer as 2-bit super-k-mers (length up to 255 bases) into `P` power-of-2 partition buffers (`P = 256` default). Because a k-mer and its reverse complement share the exact same set of canonical `m`-mers, every occurrence of a canonical k-mer routes to one unique partition.
- **Bounded memory spilling:** When buffered super-k-mers exceed half of `--max-memory-mb`, partition buffers spill to temporary files in an RAII directory and are removed immediately after processing.
- **Intra-genome Rayon parallelism:** Partitions unpack, radix-sort, and aggregate in parallel into exact uncapped `u32` counts (no `u8` 255 saturation and no `u16` 32,767 cap), emitting in one run:
  - Sorted counted k-mer table (`.ktab`, `GRMKTB1\0`)
  - Exact abundance histogram (`.hist.json`: `f1`, `f2+`, `F0`, total occurrences, peak coverage, full bins)
  - Deterministic digests (`fnv64` matching `kset::digest` for `k <= 32`, `keys_sha256`, `counts_sha256`)
  - Bottom-`s` MinHash sketch (`s = 1000` default) from the sorted bucket prefix
  - Optional FastK-style per-sequence k-mer count profile (`.prof.bin`, `GRMPRF1\0`, run-length + zigzag-delta encoded)

### 3. Storage Formats & Experimental SPSS Archive (`src/kcount/format.rs`, `src/kcount/spss.rs`)

Sorted uniformly mixed 62-bit integers (~5M per bacterial genome) have an information-theoretic lower bound of `log2(2^62 / N) + log2(e)` (~41.2 bits, or ~5.14 B/k-mer). Two complementary binary representations address random access and archival compactness separately:

- **Merge / random-access forms:**
  - `.kset` v1 (`GRMKSET1`, untouched in `src/kset.rs`): 1,024-bucket LEB128 gap encoding (~5.97 B/k-mer at `k = 31`), fast for multi-genome merges.
  - `.grk` (`GRMKGR1\0`): 1,024-bucket per-bucket bit-optimal Golomb-Rice gap baseline approaching the sorted-integer bound (~5.2 B/k-mer at `k = 31` for bacterial density).
  - `.ktab` (`GRMKTB1\0`): 1,024-bucket Golomb-Rice keys paired with exact uncapped `u32` counts (0 extra bytes when all counts equal 1; per-bucket adaptive tiered 1B/3B/5B or `u16` + `0xFFFF` escape to `u32`).
- **Experimental archive form (`spss-archive`, `GRMSPS1\0`):**
  - `grm kcount archive` and `grm kcount expand` exploit `(k - 1)`-base overlaps via deterministic greedy ProphAsm-style simplitigs packed into a contiguous 2-bit base stream plus LEB128 excess lengths and a header `fnv64` digest.
  - Typical bacterial k-mer sets archive at **~0.25 to 0.35 B/k-mer** (~18x to 23x smaller than `.kset` v1) and expand losslessly back to byte-for-byte identical `.kset` v1 (`GRMKSET1`) verified against `fnv64`.

### 4. Fold-Safe `D_core` Normalization & Plasmid Dosage (`src/kcount/dosage.rs`)

`grm kcount dosage` enforces statistical anti-leakage guardrails:

- Core panel (`--core-panel`) must match the sample table's exact `k` and input modality (`reads` vs `assembly`).
- Computes `d_core` (median count across detected core k-mers) and `core_breadth` (fraction detected at `>= min_count`).
- **Explicit missingness:** When `d_core < min_core_depth` (default `10.0`) or `core_breadth < min_core_breadth` (default `0.80`), emits `dosage_status: "insufficient_core_coverage"` and `null` (`None`) for every normalized dosage value rather than zero-filling.
- Emits a per-sample QC vector (`usable_read_len_median`, `gc_fraction`, `singleton_fraction`, `d_core`, `core_breadth`, `total_reads`) alongside target summaries so downstream pipelines can run QC-only controls.

### 5. Single-Pass Preflight & Streaming Set Algebra (`src/kcount/preflight.rs`, `src/kcount/ops.rs`)

- **`grm kcount preflight`:** `ntHash`-style O(1) rolling canonical 64-bit hash feeding a HyperLogLog (`p = 14`, 16 KB) and an adaptive hash-prefix subsampled exact frequency map (~2 MB bounded RAM) estimating `F0`, `f1`, `f2+`, peak coverage, and expected `.kset` v1 / `.grk` / `.spss` disk footprints without sorting.
- **`grm kcount op`:** Streaming bucket-by-bucket merge-join across 2+ `.kset` v1, `.ktab`, or `.grk` tables:
  - `union`, `intersect`, `diff`, `symdiff`
  - `matrix-export --mode pa|counts|dosage --format sparse|dense`, writing deterministic `manifest.json`, `features.tsv`, `missing_mask.tsv`, and `matrix.coo.tsv` or `matrix.tsv` with explicit `NA` / `is_missing = 1` entries for low-core-coverage samples.

## CLI Subcommands

```text
grm kcount count      --input <FASTA/FASTQ[.gz]>... --k <5..=64> --out <PATH> [OPTIONS]
grm kcount archive    --input <KSET|KTAB|FASTA> --out <PATH.spss> [--k <5..=64>]
grm kcount expand     --input <PATH.spss> --out <PATH.kset|PATH.ktab> [--format kset-v1|ktab]
grm kcount preflight  --input <FASTA/FASTQ[.gz]>... --k <5..=64> [--min-qual <Q>]
grm kcount dosage     --table <PATH.ktab> --core-panel <PANEL.json> [--targets <TARGETS.json>]
grm kcount op         <union|intersect|diff|symdiff|matrix-export> ...
grm kcount inspect    --input <FILE> [--sketch <S>] [--verify]
grm kcount convert-v1 --input <FILE> --out <FILE> [--direction to-gr|to-v1]
```
