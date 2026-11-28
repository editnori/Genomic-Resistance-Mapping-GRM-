# Bounded-memory parallel super-k-mer counter (`grm kcount count`)

This page documents the bounded-memory parallel super-k-mer counting engine (`src/kcount/superkmer.rs`) and the `grm kcount count` CLI command. Existing `grm count`, `grm matrix`, `.kset` v1 (`GRMKSET1`), and Revision 14 learners are unchanged.

## 1. Three-stage bounded-memory pipeline

`grm kcount count` supports `k = 5..=64` (dispatching to the `u64` path for `k <= 32` and the `u128` path for `k = 33..=64`) over plain or `.gz` (including multi-member gzip) FASTA and FASTQ inputs:

1. **Stage 1 (Parallel record-aligned ingest and minimizer partitioning):** Reads FASTA/FASTQ inputs in record-aligned chunks sized to `--max-memory-mb` (splitting single FASTA contigs that exceed the chunk size with `k - 1` base overlap while tracking full contig lengths). Worker threads reuse amortized scratch buffers to extract canonical `m`-mer minimizer super-k-mers (packed 2-bit, up to 255 bases per super-k-mer) using the public-domain SplitMix64 finalizer multipliers (Steele, Lea, and Flood, OOPSLA 2014) for strand-invariant minimizer ordering. Read/contig QC metrics (median, p10, p90, mean usable length, GC fraction, and N50) are accumulated in exact length-frequency histograms rather than per-read vectors. Whenever buffered super-k-mer bytes reach the ingest budget, partitions spill to `grm-kcount-<pid>-<nonce>/p0000.skm`.. under `--tmp-dir`.
2. **Stage 2 (Bounded-concurrency partition counting):** When the estimated in-memory footprint exceeds half of `--max-memory-mb` (or Stage 1 spilled), all partitions stream from disk and workers are capped in concurrency and per-chunk key expansion to stay within `--max-memory-mb`. Skewed partitions that exceed a worker's key cap are sub-chunked, sorted, collapsed to `(key, u64)` runs, and k-way merged with checked `u64` addition. Retained `(mixed_key, count)` pairs are routed into 256 super-bucket shard files (`s000.shd`..`s255.shd`, each covering 4 of the 1,024 output buckets) and flushed whenever shard buffers reach budget. Each partition spill file is unlinked as soon as it is read.
3. **Stage 3 (Single-super-bucket streaming finalization):** Loads one super-bucket (`1/256` of retained k-mers) at a time, sorts its 4 output buckets, computes FNV-1a64, `keys_sha256`, `counts_sha256`, and bottom-`s` MinHash sketch, streams Golomb-Rice key and count payloads directly into verified `.ktab` (`GRMKTB1\0`) and/or `.grk` (`GRMKGR1\0`) `.part` files, writes the full header SHA-256 digest in a 64 KiB streaming pass, and atomically renames into place.

## 2. Optional indexed parallel `--profile-out`

When `--profile-out` is requested, after Stages 1 to 3 finish and release their working buffers, the counter loads a compact sorted-key lookup index (`keys: Vec<K>`, optional `counts: Vec<u32>` omitted for presence-only tables, and an 18-bit sub-bucket directory for `k >= 9`, or 10-bit for `k < 9`). Inputs are streamed in record-aligned chunks, positional k-mer counts are looked up and run-length + zigzag-delta ULEB128 encoded on the fly in parallel without allocating per-read `Vec<u32>` buffers, and encoded records are written in exact input order to a verified `.prof.bin` (`GRMPRF1\0`) file. Note that `--profile-out` holds the compact index of retained k-mers in memory (about 12 bytes per retained `u64` k-mer for counted tables, or 8 bytes when presence-only), in addition to the chunk buffer.

## 3. Contracts and guardrails

- **Count ceiling (`u32` capped at `2^32 - 1`):** Exact `u64` occurrence counts are checked against `4,294,967,295` (`u32::MAX`). Exceeding `2^32 - 1` fails with an error by default; passing `--allow-saturation` clamps overflowing counts to `u32::MAX` and records the exact number of clamped retained k-mers in `saturated_counts` (both in the `.ktab` header and JSON summary) and all observed saturated k-mers in `histogram.saturated_total`.
- **Input modality (`--input-type auto|reads|assembly`):** Default `auto` detects FASTQ inputs as `reads`, FASTA inputs as `assembly`, and empty inputs as `unknown`, storing the resulting `InputModality` byte in every written `.ktab` and `.grk` header. Passing `--input-type reads` on FASTA, `--input-type assembly` on FASTQ, or mixing FASTA and FASTQ records in one run fails with an error.
- **No input filename or record-ID leaks by default:** When `--name` is omitted, headers and JSON reports record `"sample"`. In `.prof.bin`, sequence records use anonymized sequential labels (`r000001`, `r000002`, ...) unless `--keep-record-ids` is explicitly passed.
- **Minimizer length and output extensions:** `--minimizer-len` requires `1 <= m <= min(k, 31)` (help text and validation agree). Output extensions must match the format written (`.ktab` for counted tables, `.grk` for `--presence-only` and `--gr-kset-out`).
- **Spill cleanup on signals and startup:** On Unix, an async-signal-safe `SIGINT`/`SIGTERM` handler sets an atomic cancellation flag checked at chunk, spill, partition, and super-bucket boundaries so early exit unwinds RAII spill directories and `.part` files cleanly before exiting non-zero. At startup, stale `grm-kcount-<old_pid>-<nonce>` directories in `--tmp-dir` whose owning `<old_pid>` is dead (`ESRCH`) are removed without touching any live process's spill directory.
- **Determinism and `.kset` v1 parity:** Output `.ktab`, `.grk`, `.prof.bin`, `.hist.json`, and summary JSON files are byte-for-byte SHA-256 identical across thread counts (1 vs 8 threads) and across spill vs in-memory runs. For `k = 5..=32`, keys, `n_kmers`, `fnv64`, bottom-1000 sketches, and `convert-v1` `.kset` bytes match `grm count` bit for bit.

## 4. CLI usage

```text
grm kcount count --input <FASTA/FASTQ[.gz]> [--input <FILE> ...] \
  --k <5..=64> --out <OUT.ktab|OUT.grk> \
  [--name <SAMPLE>] [--input-type auto|reads|assembly] \
  [--gr-kset-out <OUT.grk>] [--hist-out <OUT.hist.json>] \
  [--meta-out <OUT.meta.json>] [--profile-out <OUT.prof.bin>] [--keep-record-ids] \
  [--min-count <N>] [--max-count <N>] [--allow-saturation] [--min-qual <Q>] \
  [--sketch-size <S>] [--max-memory-mb <MB>] [--partitions <16..=1024>] \
  [--minimizer-len <1..=min(k,31)>] [--presence-only] [--tmp-dir <DIR>]
```

## 5. Comparison with `grm count` and synthetic benchmarks

`grm count` processes multiple small bacterial assemblies concurrently in memory inside a single process and writes presence-only `.kset` v1 files for `k <= 32`; spawning one `grm kcount count` process per small assembly carries higher per-process and partitioning overhead. `grm kcount count` targets large raw read sets, `k = 5..=64`, bounded working memory, exact `u32` counts, abundance histograms, read QC metrics, and verified `.ktab`/`.grk`/`.prof.bin` outputs.

Measured on a synthetic 50 Mbp noisy FASTQ (`333,334` reads of 150 bp from a 2 Mbp synthetic sequence with about 1% random substitution noise, `12,359,806` distinct 31-mers at `--min-count 1`, `2,038,478` retained 31-mers at `--min-count 3`, release build on Ubuntu 24.04):

| Run | `--min-count` | `--max-memory-mb` | Spilled | Wall time | Peak RSS |
| --- | ---: | ---: | --- | ---: | ---: |
| Table only (`.ktab`) | 1 | 64 | yes | 10.88 s | 60,204 KiB (58.8 MiB) |
| Table only (`.ktab`) | 1 | 256 | yes | 9.00 s | 135,540 KiB (132.4 MiB) |
| Table only (`.ktab`) | 1 | 1,024 | yes | 4.20 s | 136,148 KiB (133.0 MiB) |
| Table only (`.ktab`) | 3 | 512 | yes | 2.94 s | 136,784 KiB (133.6 MiB) |
| Table + `--profile-out` | 3 | 512 | yes | 4.87 s | 135,880 KiB (132.7 MiB) |

All three `--min-count 1` `.ktab` outputs (`64`, `256`, and `1,024` MiB) matched byte for byte by SHA-256, left zero temporary files in `--tmp-dir`, and `grm kcount inspect --verify` validated the `35,372,238`-byte `.prof.bin` (`333,334` sequences, `50,000,100` bases) in `1.14 s`.
