# K-mer keys and verified binary formats (`grm kcount`)

This page documents the dual-width canonical k-mer key primitives (`src/kcount/key.rs`), the verified binary table and profile formats (`src/kcount/format.rs`), and the `grm kcount convert-v1` and `grm kcount inspect` CLI commands. Existing `grm count`, `grm matrix`, `.kset` v1 (`GRMKSET1`), and Revision 14 learners are unchanged.

## 1. Dual-width canonical keys (`src/kcount/key.rs`)

- **Fast `u64` path (`k = 5..=32`):** `KmerSpace64` delegates 2-bit encoding (`A=0, C=1, G=2, T=3`) and mixing to `crate::dna::KmerSpace`, matching `.kset` v1 keys, FNV-1a64 digests, and bottom-`s` MinHash sketches bit for bit.
- **Wide `u128` path (`k = 33..=64`):** `KmerSpace128` packs `2k` bits into a `u128` and applies a 3-round odd-multiply and right-xorshift bijection across `[0, 2^(2k))`. For `k <= 32`, `KmerSpace128` delegates to `KmerSpace64`.
- **Branchless SWAR reverse complement:** `reverse_complement_u64` and `reverse_complement_u128` reverse 2-bit bases in O(1) using pair swaps, nibble swaps, `swap_bytes`, bitwise complement, and right shift.
- **Mixer constant provenance:** Round 0 of `KmerSpace128` joins the two public-domain SplitMix64 finalizer multipliers (`0xbf58_476d_1ce4_e5b9` and `0x94d0_49bb_1331_11eb`, Steele, Lea, and Flood, OOPSLA 2014). Rounds 1 and 2 use four odd 64-bit halves derived deterministically from `SHA-256("GRM-KCOUNT-MIX128-V1:0")` through `":3"` (first 8 digest bytes interpreted little-endian with the low bit set to 1), asserted in unit tests with `ring::digest::SHA256` alongside fixed golden `mix`/`unmix` vectors for `k = 33, 45, 61, 64` and 511 test values per `k` across `k = 33..=64` (11 boundary values plus 500 pseudorandom draws per `k`).

## 2. Verified binary formats (`src/kcount/format.rs`)

At bacterial density near 5 million distinct 31-mers, 5.15 bytes per k-mer is the lower bound for a code that uses only the sorted values. The formats below provide random-access 1,024-bucket storage and per-sequence count profiles:

- **`.grk` (`GRMKGR1\0`, 72-byte fixed header):** Presence-only sorted mixed keys for `k = 5..=64`, encoded per bucket with Golomb-Rice gap coding (`best_r` chosen around `round(log2(mean_gap * ln 2))`). Stores `k`, `bucket_bits` (10), `key_width` (8 or 16), 1-byte `InputModality` (`0 = unknown`, `1 = assembly`, `2 = reads`), `n_kmers`, `name_len`, and a 32-byte SHA-256 digest covering all header metadata outside the digest slot, the sample name, bucket tables, and key payloads.
- **`.ktab` (`GRMKTB1\0`, 80-byte fixed header):** Sorted mixed keys paired with `u32` counts capped at `4,294,967,295` (`2^32 - 1`). Stores `k`, `bucket_bits`, `key_width`, flags (`0x01` presence-only, `0x02` saturated counts present), 1-byte `InputModality`, `min_count`, `max_count_observed`, `n_kmers`, `saturated_counts` (number of k-mers clamped at `2^32 - 1`), `name_len`, and a 32-byte SHA-256 digest covering header metadata, sample name, bucket tables, key payloads, and count payloads. Per-bucket count encoding selects between a 1-byte all-ones mode, tiered 1B/3B/5B encoding, or 2-byte `u16` with `0xFFFF` escape to `u32`. When all counts equal 1, `.ktab` sets the presence-only flag and writes zero count-payload bytes; its fixed header and 1,025-entry `count_offsets` table add 8,208 bytes of overhead relative to `.grk`.
- **`.prof.bin` (`GRMPRF1\0`, 72-byte fixed header):** Per-sequence positional `u32` k-mer counts encoded as run-length and zigzag-delta ULEB128 streams, storing `k`, `max_count_observed`, `n_sequences`, `total_bases`, and a 32-byte SHA-256 digest covering header metadata and all sequence records.
- **`.kset` v1 (`GRMKSET1`, unchanged in `src/kset.rs`):** Readable through `AnyTableReader` alongside `.grk` and `.ktab`.

### Integrity and bounds checks

- Every length and count field read from `.grk`, `.ktab`, and `.prof.bin` is validated with checked arithmetic against the remaining unread bytes of the file before any `Vec` or `String` allocation, and dynamic allocations use fallible `try_reserve`.
- Every full read (`GrKsetReader::read_all`, `KtabReader::read_all`, `read_profiles`, `convert-v1`, and `inspect --verify`) recomputes the 32-byte SHA-256 digest and fails if any byte in the header, metadata, bucket tables, or body differs.

## 3. CLI commands

```text
grm kcount convert-v1 --input <FILE> --out <FILE> [--direction to-gr|to-ktab|to-v1] [--name <NAME>] [--modality unknown|assembly|reads]
grm kcount inspect    --input <FILE> [--sketch <S>] [--verify] [--out <PATH>]
```

- **Default sample name:** When `--name` is omitted on `convert-v1`, the output header records `"sample"`. Input file paths and filename stems are never copied into output table headers or JSON reports by default.
- **Output extension check:** `convert-v1` refuses to write when `--out`'s extension disagrees with the format being written (`.grk` for `GRMKGR1`, `.ktab` for `GRMKTB1`, `.kset` for `GRMKSET1`). Converting a counted `.ktab` to `.kset` v1 or `.grk` retains the sorted key set only and sets `counts_dropped: true` in the JSON output when any count exceeded 1.
- **Strict `verified` flag on `inspect`:** Without `--verify`, `inspect` reads header metadata only and reports `verified: false` for all four formats (`.kset`, `.grk`, `.ktab`, `.prof.bin`). With `--verify`, `inspect` recomputes the stored header SHA-256 digest and decodes every bucket or sequence record, setting `verified: true` only when a stored digest was compared and matched. Because `.kset` v1 carries no stored digest in its header, `inspect --verify` decodes all buckets and reports `fnv64` and `keys_sha256` with `verified: false`.
