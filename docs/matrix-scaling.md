# Matrix merge beyond the packed-key ceiling and the disk budget

`grm matrix` (`src/matrix.rs`) merges one kset per genome into the species
pattern matrix. This note covers two limits on large MTB cohorts and the
software built for them:

1. The packed merge key, which capped cohorts at 4,096 genomes (k = 31) or
   1,024 (k = 32). A wide-key fallback lifts it.
2. Disk and open files: the stock merge needs every genome's kset on disk
   and open at once. A prototype out-of-core merge over column-block
   partials needs only one block of ksets at a time.

Status, 2026-10-27: software and synthetic tests only. No biological input
was merged, no 4k or 8k MTB matrix was built, and nothing in the pipeline
(`species.rs`, `main.rs`) calls the out-of-core path yet. There is no frozen
cap-4000 cohort; the MTB follow-up waits for a newly frozen one and for
fresh independent review of both the matrix change and the disk strategy.

The wide key removes the packed **algorithmic** ceiling, not the downstream
safety gate. Every public build/merge entry point, and `FeatureMatrix::load`,
refuses more than `REVIEWED_COLUMNS` (4,096) columns. `Dataset::project`, the
learners, annotation, and reports have not been reviewed past that size.
Only private functions in this module's tests can merge larger synthetic
cohorts; there is no environment variable or production configuration that
lifts the cap.

## 1. The packed key and its ceiling

A kset splits the 2k-bit value range into 2^10 buckets (`BUCKET_BITS`). The
merge reads the same bucket from every genome and sorts it by k-mer, which
groups each k-mer's genomes together. The original merge packs each entry
into one u64,

    (k-mer bits below the bucket prefix) << ceil(log2 n) | genome index

and radix sorts it. That needs `2k - 10 + ceil(log2 n) <= 64`, so at most
2^(74 - 2k) genomes:

| k | genomes the packed key holds |
|---|---|
| 32 | 1,024 |
| 31 | 4,096 |
| 30 | 16,384 |
| 27 | 1,048,576 |

Past that, the old code refused to merge (the guard) rather than truncate.
The wide key lifts only this key-size ceiling. The separate 4,096-column
review gate still applies to every public construction path and to loading
a saved matrix.

### Wide-key fallback

`merge_key(k, n)` picks the key for each merge:

- **Packed**, the unchanged fast path, whenever the inequality holds.
- **Wide** otherwise. Each entry is a `(u64 k-mer bits, u32 genome index)`
  tuple, sorted with `sort_unstable`. A bucket's k-mers are visited in the
  same ascending order and grouped by the same code, so patterns, their
  numbering, support, k-mer lists, and summary counts are what the packed
  key would give if it had the bits.

What did not change: k, genome (column) order, the presence and support
rules (`min_support <= support <= n - min_support`), pattern deduplication
and numbering, the file format, and `Dataset::project`.

Guards, all hard errors and never truncation:

- `build`, `write_partial`, and `build_from_partials` refuse any block or
  cohort ending past 4,096 columns before opening their inputs.
  `FeatureMatrix::load` applies the same gate before a saved matrix can reach
  downstream stages. Tests use the private `build_with`,
  `write_partial_with` (with a test `Limits`), and
  `build_from_partials_with`; production code has no uncapped public option.
- Mixed k and duplicate genome names are rejected before the merge. Partial
  columns carry names and per-kset FNV digests, and may not duplicate a name.
- A zero gap after a kset's first value in a bucket is a repeated k-mer and
  is rejected. The shared grouping code also rejects a repeated
  (k-mer, genome) entry. Full partial rows are derived from their actual
  presence bits, never from an entry count alone.
- The packed key is never used past its ceiling. Forcing it there (only
  possible from tests) is an error, and the packed bucket merge rechecks
  its bit budget.
- n must fit the u32 genome index and support counters (n <= 2^32 - 1).
- Pattern ids stay below 2^32 - 1, globally and per bucket.
- When opening inputs runs out of file descriptors (EMFILE/ENFILE), the
  error names the `ulimit -n` needed and points to the column-block merge.
  The partial hint counts open blocks, not the original number of genomes.

When the wide key is used, `build` prints one stderr line saying so;
supported cohorts print nothing new.

### Cost of the in-memory merge (`build`)

Per bucket b, let T_b be the number of (k-mer, genome) entries. Mixed values
are spread evenly, so T_b is about T / 1024.

| | packed | wide |
|---|---|---|
| sort | LSD radix, at most 8 byte passes, O(T_b) | pdqsort, O(T_b log T_b) |
| bucket working memory | 8 B key + 8 B radix scratch = 16 B/entry | 16 B/entry tuple, in place |

Unchanged around that: up to `threads` buckets sort at once (`species.rs`
passes `threads_chunk` = 2x the rayon threads). Every kset stays open for
the whole merge, which costs one fd and about 12.3 KB of bucket tables per
genome. The process also needs descriptors for its executable, standard
streams, and other files, so the soft limit must exceed the genome count by
sufficient headroom; exactly `n` is not enough. In particular, the
pipeline's `ulimit -n 8192` does **not** make an 8,193-kset stock build work
by default (and the public 4,096-column gate refuses it first). A private
benchmark past that fd limit must deliberately raise the limit; a future
production path would need the same explicit resource change or column
blocks. A pattern row is `ceil(n / 64)` words (1,032 B at 8,193 genomes),
held twice during the merge. Each kept k-mer costs 16 B in the pair list
plus 8 B in `kmers`.

## 2. Disk: the column-block (out-of-core) merge

### The constraint

The MTB scale-up estimates about 96 GB of temporary ksets for a stock
4,000-genome merge, i.e. about 24 MB per genome. About 47 GB is free and a
30 GiB (32.2 GB) reserve must remain, so about 14.8 GB is usable. The stock
merge needs every kset to exist when the merge starts, so it cannot fit.
It also needs about 4,000 open files at once.

### Assessment

A kset spends about 6 bytes (one LEB128 gap) per (genome, k-mer). An MTB
cohort is clonal: most k-mers occur in nearly every genome, so the same
information as a presence bitset costs 1 bit per (genome, k-mer) plus one
k-mer value per block. Storing blocks of genomes that way, then deleting
their ksets, bounds temporary disk to one block of ksets plus the much
smaller blocks already folded. Merging the blocks k-mer by k-mer needs one
open file per block, and no genome-index key at all, so the packed ceiling
does not apply either.

### Prototype (in `src/matrix.rs`)

- `write_partial(ksets, columns, first_col, out, threads_chunk, reserve_bytes)`
  reads one block of genome ksets (columns
  `first_col..first_col + len`, where `first_col` is a multiple of 64).
  `columns` is the expected genome name and `kset::digest` for every column,
  as recorded in `counted.json`; every kset is read and checked against both
  before a rewrite begins. The partial holds each distinct k-mer once, in
  ascending order per bucket, as a LEB128 `gap << 1 | full`. Unless `full`
  is set (all of the block's actual presence bits are set), the row's
  `ceil(len / 64)` words follow. Nothing is filtered, because support is only
  known once every block is merged.
- The version-2 partial header stores `n_kmers` and `full_rows`, followed by
  each column's name and kset digest. The adjacent
  `<out>.receipt.json` repeats those fields and records the partial's byte
  length and SHA-256. The partial and receipt are written through temporary
  files, `sync_all`ed, renamed, and followed by a parent-directory sync.
  Success therefore means both are durable before a caller deletes ksets.
  Temporary suffixes append to the whole file name, so `x.a` and `x.b`
  cannot collide.
- A rewrite retires the old receipt and partial before writing the new one.
  An error after creation removes the temporary file; a missing receipt, a
  leftover `.part`, or a partial without its matching receipt is not a
  resumable block. A failed rewrite consequently leaves no stale block that
  a later merge can accept.
- `build_from_partials(partials, columns, min_support, threads_chunk)`
  requires the expected genome/kset-digest list. It opens one fd per block,
  verifies every receipt's SHA and byte count against its file, checks the
  receipt against the header, and then checks every column against the
  expected name and digest. Thus a catalog change, re-assembly/recount, file
  replacement, edited receipt, or interrupted rewrite is refused rather than
  silently reused. Blocks must share k and tile columns 0..n in order,
  without gaps, overlaps, or duplicates. For each bucket it k-way merges the
  blocks by k-mer, joins their rows at word offsets, and uses the same
  pattern collector and global assembly as `build`.
- The partial reader validates magic/version, reserved zero bytes, k, bucket
  width, block alignment and cap, bounded name and width allocations,
  per-bucket counts and offsets, strictly ascending in-bucket k-mers, clear
  tail bits, non-full row support, and the header's full-row total during the
  merge. The supported-cohort result is byte-identical to `build`, including
  genome order. As with `build`, the pipeline caller must fill the matrix
  summary's `kset_digests` from its counted manifest.

The disk-sparse sequence a caller would run, which the small equivalence
test and benchmark simulate, is:

    derive expected names and kset digests from the counted manifest
    for each block of W genomes, in column order:
        ensure_free_space(block dir, estimated new bytes, reserve)
        count the block's genomes to ksets
        write_partial(block ksets, expected columns, first_col, reserve)
        only after success, delete the block's ksets
    build_from_partials(all partials, all expected columns)
    save_with_reserve(matrix dir, reserve)

Resume must reconstruct the same expected list and rewrite every missing or
stale block. Receipts make that decision explicit; they do not replace the
caller's counted-manifest bookkeeping.

### Resource model

- Temporary disk: about W x (kset size) for the block being folded, plus
  the partials so far, plus the final matrix. A block's partial is at most
  D_block x (ceil(W/64) x 8 + ~7) bytes for D_block distinct k-mers in the
  block, and less where rows are full. Over the cohort that is about
  n x D x (1/8 + 7/W) bytes when no rows are full.
- Open files: W ksets while writing a partial; n/W partials for the final
  merge.
- Memory: writing a partial is `build` over W genomes (16 B per entry of a
  bucket, a few buckets at once). The final merge decodes one bucket of
  every block: about D_b x (8 + W/8) bytes per block, then the same global
  assembly as `build`.

`ensure_free_space(dir, need, reserve)` provides the caller's pre-count
check on 64-bit Linux (`statvfs` available bytes). `write_partial` checks
its header/table and each encoded chunk against the reserve while writing;
`FeatureMatrix::save_with_reserve` checks the complete matrix save, with
files from the save it replaces counted as reusable space. The merge itself
only reads partials; its new output is the subsequent matrix save. These
checks are best-effort preflights, not space reservations: another writer,
metadata allocation, or delayed allocation can still cause ENOSPC. Write
errors propagate and partial temporaries are cleaned up, but a real full
filesystem and power-loss recovery have not been fault-injected.

For the stated 4,000-genome budget (about 24 MB per kset, from the 96 GB
estimate), and assuming D_block of 4.4M to 8.8M distinct 31-mers per block
(not measured):

| W | block ksets | all partials, no full rows, D_block 8.8M | peak temp. (excl. matrix) |
|---|---|---|---|
| 1,024 | ~24.6 GB | ~4.6 GB | over budget |
| 256 | ~6.1 GB | ~5.4 GB | ~11.5 GB |
| 128 | ~3.1 GB | ~6.3 GB | ~9.4 GB |

With D_block 4.4M or a high full-row share, the partials shrink
accordingly. Only W of about 256 or below fits the ~14.8 GB headroom, and
only if the final matrix and other outputs fit the rest. That has not been
estimated (see below).

## Evidence

### Tests (`cargo test --locked --offline`, all in `src/matrix.rs`)

Wide key:

- `merge_key_is_packed_exactly_when_the_key_fits_64_bits`: named
  boundaries (k = 31: 4,096 packed, 4,097 / 8,193 / 2^20 wide; k = 32:
  1,024 packed, 1,025 / 8,193 wide; k = 30 at 16,384 / 16,385), every k
  from 5 to 32 against an independently computed inequality, and rejection
  of n = 2^32.
- `supported_cohorts_keep_their_digests_under_both_keys`: SHA-256 digests
  of five saved synthetic matrices (k = 5 to 32, up to 4,096 genomes) were
  recorded with the code before this change and are pinned. The packed
  path and the forced wide key reproduce them byte for byte.
- `k31_merges_at_the_packed_ceiling_of_4096`, `..._past_the_ceiling_at_4097`,
  `k31_merges_8193_genomes`, `k32_merges_at_the_packed_ceiling_of_1024`,
  `..._past_the_ceiling_at_1025`, `k32_merges_4097_genomes`: file-backed
  synthetic cohorts against a plain BTreeMap oracle (patterns, numbering,
  support, k-mer lists, distinct/rare/core counts, genome names and order).
  Cases past 4,096 use only private test entry points. They include k-mers
  at both ends of the value space and at bucket edges, and k-mers private to
  the last genome, whose index needs the top bit. Each also checks that the
  forced wide key and partial path give the same saved bytes, and that
  forcing the packed key past its algorithmic ceiling errors.
- `merges_refuse_more_columns_than_reviewed_unless_a_test_lifts_the_cap`
  proves public stock, partial-write, and partial-merge paths refuse 4,097
  columns before opening inputs, and that a synthetic oversized matrix made
  through private test code is refused by `FeatureMatrix::load`.
- `rejects_no_ksets_and_mixed_k_at_any_size`; `small_and_empty_cohorts_...`
  (k = 5, 7, 31, 32; 1, 2, 3, 64, 65 genomes; full, one empty, all empty;
  min_support 0, 1, 2; packed, wide, and 64-genome partials all agree);
  `duplicate_patterns_merge_across_buckets_under_both_keys`;
  `digest_and_genome_order_do_not_depend_on_scheduling` (`threads_chunk`
  1, 7, 1024; reversed input order); `bucket_merge_past_both_ceilings_...`
  (70,000 genomes in memory, private test code only).

Out-of-core integrity and scaling:

- `partial_receipts_describe_their_blocks`: names, per-column kset digests,
  k-mers, full rows, header totals, bytes and SHA-256 match independent
  counts and the saved receipts. Partial bytes are identical for
  `threads_chunk` 1, 7, 1024. A block of identical genomes stores no rows
  and is more than 32x smaller than its ksets.
- `partials_must_tile_the_columns_in_order_with_one_k`: swapped, missing,
  duplicated, leading-gap and reordered-tail partials are each rejected;
  so are incomplete column lists, an empty list, a block start off a
  64-column boundary, wrong expected-column counts, and mixed k.
- `partial_files_reject_corruption`: bad magic/version, reserved bytes, k,
  bucket width, first column and cap, zero/huge width, k-mer and full-row
  totals, overlong names, bucket counts one more/fewer, trailing padding,
  and six truncation points. It also exercises file/receipt mismatches and
  forged receipts, including a false full-row total.
- `receipts_bind_each_partial_to_its_file_and_columns`: missing, malformed,
  edited, swapped and forged receipts; same-genome but changed-kset file
  replacement; leftover `.part` inputs; successful rewrite; and distinct
  temporary/receipt names for `x.a` and `x.b`.
- `stale_partials_are_refused_after_catalog_or_assembly_changes`: an
  assembly/recount and a catalog replacement invalidate exactly the affected
  columns; failed rewrites cannot resurrect old content; valid rewrites
  again match stock `build`; duplicate expected names are refused.
- `a_kset_that_repeats_a_kmer_is_refused_by_every_path`: the reviewed case
  where g0 repeats v and g1 lacks it is rejected by packed, wide, and
  partial paths and by the shared grouper. The clean two-genome version
  keeps g1's bit clear and agrees with the oracle and stock build.
- `dense_buckets_and_uneven_blocks_match_build_under_any_schedule`: crowded
  buckets with real non-full rows, uneven and short final blocks, write and
  merge chunk sizes 1, 7, and 1024, all compared with packed/wide `build`.
- `free_space_is_measured_and_kept_for_partials_and_saves`: statvfs against
  `df`, reserve refusal before output, exact replacement-size accounting,
  and no summary after a refused save. The stale-partial test injects a
  mid-write free-space failure and checks temporary/receipt cleanup and
  refusal on resume.
- `the_disk_sparse_sequence_makes_the_matrix_build_makes`: the same block
  helper used by `bench_block_case` deletes each synthetic block's ksets,
  leaves only partials/receipts, and matches stock `build`.
- `partials_merge_under_an_open_file_limit_that_stops_build`: reruns the
  test binary under `ulimit -Sn 96`. There, `build` on 200 genomes fails
  with the `ulimit -n` hint, while 64-genome partials complete. The child's
  matrix digest equals `build`'s digest from the parent, which has no limit.

Tests that hold thousands of ksets open take turns behind a mutex and check
`/proc/self/limits` first. They need a soft open-file limit of at least
8,449; this machine has 10,240. Debug `cargo test --lib matrix::` takes
about 3 minutes on this shared machine, mostly writing and opening small
synthetic ksets.

### Resource benchmark

    cargo test --release --lib matrix::tests::bench_merge_resources -- --ignored --nocapture

Each case runs in a fresh child process under an explicit `ulimit -Sn`
(the cohort plus 64 for `build`, the block width plus 64 for the
out-of-core merge) and `ulimit -Sv` of 4 GiB, so each peak RSS (VmHWM) is
that case's own. These are earlier synthetic measurements, not a biological
run and not a default-resource claim. The ignored benchmark deliberately
uses private test-only uncapped functions and raises fd limits; it was not
rerun after the integrity fixes. It is not a pipeline command or evidence
that production may exceed 4,096 columns. Original benchmark environment:
12-thread WSL2 machine, release build, 2026-10-27. Wall times vary between
runs on this shared machine: a later run measured 0.4 to 1.6 s for the
`files` rows below, with RSS unchanged.

`files`: tiny cohorts merged by `build` (1,024 pooled values, each genome
keeping about half). At this size, per-(genome, bucket) overhead dominates,
not sorting.

| k | n | key | ulimit -n | entries | merge s | peak RSS MiB |
|---|---|---|---|---|---|---|
| 31 | 1,024 | packed | 1,088 | 524,766 | 0.315 | 18.1 |
| 31 | 1,024 | wide | 1,088 | 524,766 | 0.300 | 18.2 |
| 31 | 4,096 | packed | 4,160 | 2,097,546 | 0.887 | 56.4 |
| 31 | 4,096 | wide | 4,160 | 2,097,546 | 0.912 | 56.9 |
| 31 | 4,097 | wide | 4,161 | 2,098,062 | 1.031 | 56.8 |
| 31 | 8,193 | wide | 8,257 | 4,196,107 | 2.590 | 108.5 |
| 32 | 1,024 | packed | 1,088 | 524,164 | 0.298 | 18.4 |
| 32 | 1,025 | wide | 1,089 | 524,674 | 0.267 | 18.2 |
| 32 | 4,097 | wide | 4,161 | 2,098,723 | 0.880 | 57.3 |

`memory`: one bucket in memory at an assumed MTB-like density scaled down
8x (540 values per genome in the bucket, 90% kept), single thread, to
isolate sort-and-group cost per entry.

| k | n | key | entries | ns/entry | peak RSS MiB (setup) |
|---|---|---|---|---|---|
| 31 | 1,024 | packed | 552,470 | 83 | 21.4 (12.5) |
| 31 | 1,024 | wide | 552,470 | 125 | 21.2 (12.4) |
| 31 | 4,096 | packed | 2,211,992 | 78 | 71.1 (36.6) |
| 31 | 4,096 | wide | 2,211,992 | 128 | 71.1 (36.5) |
| 31 | 8,193 | wide | 4,424,289 | 117 | 137.6 (68.8) |
| 32 | 1,024 | packed | 552,616 | 89 | 21.2 (12.5) |
| 32 | 1,024 | wide | 552,616 | 152 | 21.4 (12.6) |
| 32 | 8,193 | wide | 4,424,370 | 135 | 137.5 (68.6) |

The wide key costs about 1.5 to 1.7x the packed key per entry, and both
add the same ~16 B/entry over setup (35 MiB for 2.2M entries, 69 MiB for
4.4M), as modelled.

`blocks`: the disk-sparse sequence above, run end to end. The cohort is
clonal-looking: 4,096 pooled values, nine in ten core (each genome misses
one with the given dropout), the rest in about half the genomes. "All
ksets" is what the stock merge would need on disk at once. "Peak temp." is
the most that existed at once here: one block's ksets plus the finished
partials.

| k | n | W | dropout | ulimit -n | all ksets MiB | peak temp. MiB | partials MiB | full rows | write s | merge s | peak RSS MiB |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 31 | 4,097 | 512 | 0.2% | 576 | 162.5 | 22.0 | 1.77 | 14,413 | 1.81 | 0.03 | 14.4 |
| 31 | 8,193 | 1,024 | 0.2% | 1,088 | 325.0 | 44.6 | 3.99 | 7,685 | 3.95 | 0.09 | 21.4 |
| 31 | 8,193 | 1,024 | 0 | 1,088 | 325.4 | 41.5 | 0.85 | 33,382 | 5.19 | 0.04 | 21.0 |
| 31 | 8,193 | 256 | 0.2% | 320 | 325.0 | 13.4 | 3.31 | 74,533 | 4.35 | 0.07 | 17.1 |
| 32 | 4,097 | 1,024 | 0.2% | 1,088 | 167.4 | 43.8 | 2.02 | 5,763 | 2.06 | 0.04 | 20.8 |

On this synthetic data, the partials are 80x to 380x smaller than the
ksets they replace, and peak temporary disk follows the block width. Its
byte equivalence with `build` is established by the tests above, not by
this benchmark.

### Packed path against the code before this change

To check that the shared refactor did not slow the packed path, the code
before this change and the current code were built in scratch copies
outside the repository, with the same benchmark appended.

End to end (`build` on tiny file-backed cohorts: k = 31 at 1,024 / 4,096,
and 2,048 genomes with a 16,384-value pool; k = 32 at 1,024; 7 repeats per
run, old and new interleaved for 3 to 4 rounds at 4 and 12 rayon threads),
both copies produced identical outputs. The timings were dominated by noise
on this shared machine. The new/old ratio of median times ranged from 0.85
to 1.33, and its sign flipped between runs (one run had new slower at 4
threads, the next at 12). These 0.2 to 3 s runs cannot resolve a
difference of a few percent.

In memory, single-threaded, the pre-change bucket merge was copied verbatim
(reading from memory instead of ksets) and timed against the new
`merge_entries` path on identical buckets. There were 15 alternating
rounds per case, two runs, and outputs were asserted equal:

| case | entries | new/old median, run 1 | run 2 |
|---|---|---|---|
| k = 31, 4,096 genomes | 2,211,609 | 0.955 | 1.117 |
| k = 31, 2,048 genomes, 4,000-value pool | 4,096,325 | 1.013 | 0.975 |
| k = 32, 1,024 genomes | 553,115 | 1.021 | 1.034 |
| k = 31, 64 genomes, 20,000-value pool | 640,314 | 1.027 | 0.964 |

That is no consistent slowdown: 6 of 8 ratios are within +/-4%, and the
outlier's run 2 also had the new path's fastest single round (53.5 vs 59.1
ns/entry). A dedicated benchmark on a quiet machine would be needed to rule
out a difference of a few percent.

## What this enables

- Within the reviewed 4,096-column gate, `build` keeps the packed path and
  its bytes wherever that key fits, and can use the wide algorithm when a
  packed ceiling below that gate is crossed (for example k = 32 past 1,024
  genomes). The five pinned supported-cohort digests remain unchanged.
- Private synthetic tests show that the wide and partial algorithms can be
  correct past their old key and fd limits. That is not production
  permission to exceed 4,096 columns.
- `write_partial` + `build_from_partials` now provide content-bound,
  durable, resumable building blocks for a future disk-sparse integration.
  They are not wired into the pipeline, and their synthetic resource
  measurements do not prove a real cap-4000 MTB merge fits.

## What remains to be proven or built

Fresh independent review of this integrity/scaling change comes before any
biological run. After that, at minimum:

- **Pipeline integration (not built).** Nothing calls the partial path.
  The count -> partial -> delete loop belongs in `species.rs` / `main.rs`
  and needs its own change and review. It must persist and reconstruct the
  expected names and kset digests, use the free-space/reserve APIs with a
  defensible estimate, delete ksets only after durable partial+receipt
  success, and rewrite missing/stale blocks on resume. The current tests
  exercise library calls, not pipeline wiring or manifest bookkeeping.
- **Other stages read ksets.** The count cache treats an existing kset as a
  cache hit. `structure` takes MinHash sketches from ksets. Probe panels
  (`analysis.rs` -> `probes::build_probes`) read ksets and fall back when
  they are unreadable. Deleting ksets block by block means storing
  sketches first, accepting recounts, and deciding how probes are built.
  None of this is done.
- **Downstream review before any cap lift.** `Dataset::project`, learners,
  annotation, and reports were not reviewed or tested past 4,096 columns.
  Lifting the gate requires an explicit reviewed code change and downstream
  evidence, not a test helper, benchmark flag, or larger `ulimit`.
- **Real-data pilot before 4k/8k.** On a newly frozen, small biological
  cohort, first compare wired partial output with stock `build`, then test
  catalog/assembly invalidation and interrupted resume. A cap-4000 pilot
  follows only after those checks and resource measurement. No biological
  pilot has run, and no real 4k or 8k resource claim is established.
- **Real sizes.** Per-genome kset size (24 MB is the scale-up estimate),
  distinct k-mers per block, full-row share, and partial size on real MTB
  ksets are unmeasured. The disk table above rests on assumed D_block
  values; the old synthetic block benchmark was not rerun after these
  integrity changes.
- **Final matrix size.** `patterns.bin` is n_patterns x ceil(n/64) x 8
  bytes and `kmers.bin` is 8 B per kept k-mer. How many patterns and kept
  k-mers real 4k or 8k MTB data gives is unknown, so the matrix's own disk
  and RAM, which must also fit the headroom, are not estimated.
- **Runtime and memory at real scale.** Benchmarked buckets reach 4.4M
  entries; a real 8k `build` bucket could be much larger. Sort cost, memory
  bandwidth across concurrent buckets, kset I/O, and partial merge runtime
  are extrapolated, not measured. The earlier 120 to 140 ns/entry synthetic
  wide-key figures and 16 B x T_b per concurrent bucket are not real-data
  guarantees. An 8,193-input stock build additionally needs more than 8,192
  available fds plus headroom, deliberate limit changes, and downstream
  review; it does not work under the pipeline defaults.
- **Failure evidence.** Tests inject free-space exhaustion and corrupt or
  stale files, but do not fill a real filesystem, kill a process mid-write,
  or verify power-loss recovery on the deployment filesystem. fsync and
  directory-sync behavior therefore remains a deployment validation item.
