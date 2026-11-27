# K-mer counting and multi-genome matrix technical review

This review audits GRM's Rust k-mer counting, storage, multi-genome presence-matrix merge, unitig compaction, and sketch-clustering pipeline (`src/dna.rs`, `src/fnv.rs`, `src/fasta.rs`, `src/kset.rs`, `src/matrix.rs`, `src/input_build.rs`, `src/cluster.rs`, `src/species.rs`, `src/learn/mod.rs`) against 13 external k-mer and colored de Bruijn graph systems, cites measured repository baselines, specifies a two-tier optimization roadmap, and defines the reproducible benchmark harness under `scripts/kmer-bench/`.

This document and the benchmark harness change no production Rust code (`src/`), no `Cargo.toml`/`Cargo.lock` pin, and no Revision 14 output. Tier-1 items below are designed to preserve bit-exact Revision 14 matrices and reports; Tier-2 items introduce new representations or inputs that require separate evaluation before any release change.

## 1. GRM's current k-mer pipeline

GRM builds its species feature matrix in four stages with zero external binary dependencies in the default path:

1. **Canonical 2-bit encoding and invertible mixing (`src/dna.rs`).** `KmerSpace::new` (`src/dna.rs:48-76`) supports `MIN_K = 5` to `MAX_K = 32` (`src/dna.rs:17-18`), packing bases (`A=0, C=1, G=2, T=3` via `CODE` at `src/dna.rs:20-34`) into the low `2k` bits of a `u64`. `for_each_canonical` (`src/dna.rs:103-121`) maintains rolling forward and reverse-complement words, resets the window on any non-ACGT byte (so no k-mer spans `N` or IUPAC codes), and emits `fwd.min(rev)`. Every stored value passes through `KmerSpace::mix` (`src/dna.rs:83-90`), a bijective SplitMix64 finalizer scaled to `2k` bits with exact inverse `unmix` (`src/dna.rs:92-99`). Mixing solves two problems at once: raw canonical values concentrate in the lower half of the `2k`-bit range (>70% below midpoint, tested at `src/dna.rs:267-294`), which would skew prefix bucketing, and sorted mixed values directly form a bottom-`s` MinHash sketch (`KsetReader::sketch`, `src/kset.rs:395-406`) without a separate hashing pass.
2. **Sequence reading and assembly QC (`src/fasta.rs`).** `fasta::open` (`src/fasta.rs:11-24`) checks the first two bytes for gzip magic (`0x1f 0x8b`) and wraps `MultiGzDecoder` or a plain `File` in a 1 MiB `BufReader`. `for_each_record` (`src/fasta.rs:28-72`) streams multi-line FASTA and 4-line FASTQ records. `StatsBuilder` (`src/fasta.rs:96-143`) tracks contig lengths, GC, N50, and ambiguous base counts.
3. **Per-genome counting, radix sort, and `.kset` v1 storage (`src/kset.rs`, `src/fnv.rs`).** `kset::count_into` (`src/kset.rs:41-57`) reuses caller-owned `values` and `scratch` buffers per Rayon worker (`src/species.rs:371-391`), pushes mixed k-mers into `values`, sorts with an 8-bit LSD radix sort (`radix_sort`, `src/kset.rs:61-93`, skipping constant byte positions), and collapses duplicates at `min_abundance` (`dedup_with_min_count`, `src/kset.rs:97-117`). `kset::encode` (`src/kset.rs:202-221`) partitions values into `BUCKETS = 1024` (`BUCKET_BITS = 10`, `src/kset.rs:22-23`) by top 10 bits and writes unsigned LEB128 gaps from each bucket floor (`GRMKSET1`, ~5.9 to 6.0 B/k-mer on 5 Mbp bacterial genomes). `kset::digest` (`src/kset.rs:265-271`) hashes the sorted mixed `u64` stream with deterministic 64-bit FNV-1a (`src/fnv.rs:6-32`).
4. **Multi-genome presence matrix merge (`src/matrix.rs`).**
   - Stock `matrix::build` (`src/matrix.rs:219-226`) merges one bucket at a time across Rayon chunks (`threads_chunk = 2 * threads`, `src/species.rs:598`). When `2k - 10 + ceil(log2 n) <= 64` (up to 4,096 genomes at `k = 31`, 1,024 at `k = 32`), `MergeKey::Packed` (`src/matrix.rs:182-206`, `453-472`) packs `(low_bits << gbits) | genome_index` into one `u64` and LSD-radix-sorts it in `O(T_b)` time. Beyond that ceiling, `MergeKey::Wide` (`src/matrix.rs:473-485`) falls back to `Vec<(u64, u32)>` sorted with `sort_unstable` (`O(T_b log T_b)`).
   - `group_sorted` (`src/matrix.rs:494-528`) builds each k-mer's presence bitset row (`words = ceil(n / 64)` `u64` words), drops k-mers outside `[min_support, n - min_support]`, and collapses identical presence rows into unique *patterns* inside `BucketMerge::add` (`src/matrix.rs:153-171`) and globally in `assemble` (`src/matrix.rs:300-397`). Public entry points enforce `REVIEWED_COLUMNS = 4096` (`src/matrix.rs:46-56`) because downstream stages (`Dataset::project`, learners, annotation, reports) have not been reviewed past 4,096 columns.
   - Out-of-core column-block partials (`GRMKPRT2`, `src/matrix.rs:530-1666`) fold `W`-genome blocks (default `W = 256` under `--disk-sparse`, `src/main.rs:777`), eliding all-present ("full") rows, binding every block with SHA-256 receipts and per-column FNV digests, and enforcing a 30 GiB free-space reserve (`src/species.rs:29`).
5. **Native unitig compaction (`src/input_build.rs`) and lineage clustering (`src/cluster.rs`).** For plan block C3's unitig arm, `input_build::union_values` (`src/input_build.rs:390-409`) takes the bucket-wise union of all canonical 31-mers across the species' `.kset` files (single-genome k-mers included), indexes them in a prefix-indexed sorted `Table` (`src/input_build.rs:336-383`), computes unambiguous left/right de Bruijn links in parallel (`src/input_build.rs:446-502`), and runs union-find compaction (`src/input_build.rs:503-533`). `cluster::distance_matrix` (`src/cluster.rs:43-60`) computes all pairwise Mash distances from bottom-1,000 sketches for UPGMA lineage clustering.

## 2. Measured repository baselines

The following measurements are recorded in committed repository files and serve as reference baselines for the benchmark suite.

### 2.1 Per-species counting, matrix, and C3 input-build wall times

From `reference/input-arms-m/first-run-times.tsv` (16-CPU, 32 GB machines, whole-second step logs) and `docs/input-arms-result.md` (cohort sizes):

| Species | Genomes | `count` (s) | `matrix` (s) | `structure` (s) | `inputs_counts` (s) | `inputs_unitigs` (s) | `inputs_aa` (s) | `prodigal` (s) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| *Escherichia coli* | 2,000 | 49 | 59 | 2 | 26 | 105 | 30 | 429 |
| *Klebsiella pneumoniae* | 2,000 | 92 | 110 | 2 | 55 | 143 | 42 | 899 |
| *Mycobacterium tuberculosis* | 2,000 | 97 | 83 | 1 | 83 | 93 | 36 | 747 |
| *Staphylococcus aureus* | 2,000 | 101 | 225 | 2 | 82 | 108 | 36 | 514 |
| *Enterococcus faecium* | 2,000 | 66 | 81 | 2 | 59 | 74 | 25 | 290 |
| *Neisseria gonorrhoeae* | 2,000 | 51 | 43 | 1 | 42 | 26 | 11 | 214 |
| *Streptococcus pneumoniae* | 2,000 | 49 | 43 | 1 | 33 | 34 | 12 | 205 |
| *Acinetobacter baumannii* | 1,702 | 134 | 194 | 3 | 98 | 174 | 48 | 549 |
| *Pseudomonas aeruginosa* | 1,421 | 127 | 140 | 1 | 89 | 147 | 49 | 666 |

Key observation: on 2,000 *E. coli* assemblies (10.232 Gbp total sequence), stock `grm count` (49 s) plus `grm matrix` (59 s) completes in **108 s** end to end (~208 Mbp/s counting throughput across 16 cores), roughly 4x faster than running Prodigal gene calling alone (429 s) and far below the downstream learner fit times (~719 to 805 s per draw for `fit-k31`).

### 2.2 Native unitig compaction vs. `bcalm2 2.2.3`

From `docs/input-arms.md:32` and `docs/assumption-plan-corrections.md:444-452`:

- ***Klebsiella aerogenes*** (93 genomes; 31,083,227 distinct 31-mers; 1,370,388 unitigs): GRM and `bcalm2 2.2.3` (`-kmer-size 31 -abundance-min 1`) produced identical unitig sets (1,370,388 in both, 0 in one only). GRM's unitig build plus comparison took 14.81 s versus 55.28 s for `bcalm2` alone (2,404,604 KB peak RSS).
- ***Escherichia coli*** (2,000 genomes, 10.232 Gbp; 104,518,454 distinct 31-mers; 5,493,343 unitigs): exact set equality (5,493,343 in both, 0 in one only). GRM standalone `exp-input-build unitigs` ran in **214.22 s** at **2,930,764 KB** (2.80 GiB) peak RSS (232.75 s / 4,268,980 KB when also verifying against `bcalm2`). `bcalm2 2.2.3` alone took **689.13 s** at **9,445,936 KB** (9.01 GiB) peak RSS, before GRM even parsed `bcalm2`'s output FASTA. GRM's native Rust compaction was **3.22x faster** and used **3.22x lower peak RAM** than `bcalm2 2.2.3`.

### 2.3 Disk-sparse column-block merge and wide-key scaling

From `docs/matrix-scaling.md:626-639` and `docs/matrix-scaling.md:728-735`:

- **4,000-genome *M. tuberculosis* (`--disk-sparse --block-width 128`)**: four matrices plus checks and duplicate screening took 57 min 10 s; peak disk growth during the full matrix build was **6,226,591,744 B (5.80 GiB)** versus ~96 GB estimated to hold all 4,000 `.kset` files simultaneously, and the kept full matrix was **367,757,922 B (350.7 MiB)**.
- **Packed vs. Wide sort cost (single-thread bucket microbenchmark, `docs/matrix-scaling.md:626-639`)**: at 4,096 genomes (`k = 31`, 2,211,992 entries), `MergeKey::Packed` took **78 ns/entry** vs. **128 ns/entry** for `MergeKey::Wide` (1.64x slower, 1.5x to 1.7x across cases), directly attributable to `sort_unstable` (`O(T_b log T_b)`) replacing LSD radix sort (`O(T_b)`).

### 2.4 C3 empirical finding on assembly-only k-mer counts

From `docs/input-arms-result.md:17,38-42,102-106` and `docs/input-arms-m-count.md:15-18,56-60,86-90`:

- Adding assembly-only `[count >= 2]` 31-mer features (`k31-counts`) carried **0 of 5 learners** on both the first run (232 datasets) and the combined count with M (267 datasets first reading, 244 second reading). Median lineage balanced accuracy differences were 0.000 (`scm`, `cart`, `rf`, `lasso`) to 0.001 (`gbt`), with Holm-adjusted sign-test p-values 0.270 to 1.000 on the first run (0.204 to 1.000 on all datasets).
- **Technical explanation**: short-read assemblers (SPAdes, SKESA, Velvet) collapse high-copy plasmids and tandem repeats into single-copy contigs in FASTA output. Counting `>= 2` inside an assembly therefore captures paralogous insertion sequences and duplicated rRNA operons rather than plasmid copy-number dosage. Recovering true gene/plasmid copy-number amplification (for example beta-lactamase or efflux copy dosage) requires counting **raw sequencing reads (FASTQ)** and normalizing per-k-mer depth against single-copy chromosomal median coverage.

## 3. Comparison with external k-mer and colored de Bruijn graph tools

The table below compares GRM against 13 external systems across single-sample k-mer counting, streaming estimation, multi-sample presence matrices, compacted de Bruijn graphs, and colored k-mer indexes.

| Tool | Primary Role | Core Algorithm / Data Structure | `k` Support | Abundance Counts | Disk / RAM Profile | Multi-Genome Matrix / Colors | Relevance to GRM |
| --- | --- | --- | --- | --- | --- | --- | --- |
| **GRM** (`src/`) | Assembly k-mer sets, collapsed presence `BitMatrix`, unitigs | 2-bit canonical + invertible SplitMix64 mixer, 8-bit LSD radix sort, 1,024-bucket delta-LEB128 `.kset`, packed/wide bucket merge + pattern collapse | `5..=32` (`u64`) | Floor filter only (`min_abundance`); binary presence stored | ~5.9 B/k-mer `.kset`; 5.80 GiB peak temp on 4k MTB (`--disk-sparse`); 2.80 GiB RAM on 2k *E. coli* unitigs | Exact pattern-collapsed `BitMatrix` (capped at `REVIEWED_COLUMNS = 4096`) | Baseline: fast for assemblies <= 4,096 genomes, zero external deps, free bottom-`s` MinHash sketches |
| **FastK** | High-speed single-sample counter & profiles | Minimizer super-k-mer binning + multi-threaded radix sort + prefix-compressed `.ktab` / positional `.prof` | Arbitrary (optimized `<= 64`) | `u8` (capped at 255 in `.ktab`; exact histogram `.hist`) | Very low peak RAM (bounded `-M`), compact prefix-stripped `.ktab` (~4-5 B/k-mer) | Single-sample (external set-ops via `Logex` / `Fastmerge`) | Top reference for Tier-2 raw FASTQ counting (T2) and prefix-stripped `.kset` v2 layout |
| **KMC3** | Disk-based exact counter for large read sets | Two-stage signature/minimizer super-k-mer disk bins + parallel radix sort | `1..=256` | Exact up to `u32` (configurable `-cs`, e.g. 65535) | Strict RAM cap (`-m`), heavy temp disk on raw reads (~1x input FASTQ), compact `.kmc_pre/.kmc_suf` | Single-sample (`kmc_tools` pairwise set ops) | Gold-standard exact abundance reference for T1/T2 and arbitrary `k > 32` |
| **Jellyfish 2** | In-memory exact k-mer counter | Lock-free CAS reprobing hash table with key quotienting; spills to disk merge if undersized | Arbitrary (`u64` fast path `<= 32`) | Exact arbitrary-width via overflow bits | High RAM on noisy raw FASTQ unless pre-sized or Bloom-filtered; fast on clean assemblies | Single-sample (`jellyfish merge`) | Classic exact baseline for T1/T2 parity and lock-free hashing trade-offs |
| **DSK** | Memory-bounded exact disk counter | Multi-pass hash partitioning into disk slices within fixed RAM/disk budget | Arbitrary (compile-time multi-word) | Exact (`u32`) | Minimal RAM (can run in <1 GB), higher wall time from multiple passes | Single-sample | Reference for extreme memory-constrained counting; slower than KMC3/FastK on NVMe nodes |
| **Gerbil** | High-`k` & optional GPU counter | Two-stage minimizer super-k-mer disk partitioning + hash-table counting (CPU or CUDA) | Large `k` (fast `k >= 32`) | Exact (`u32`) | Low RAM, competitive disk footprint on huge read sets | Single-sample | Strong reference if Tier-2 scales raw FASTQ counting to `k > 32` or GPU nodes |
| **ntCard** | Streaming cardinality & histogram estimator | ntHash rolling hash + sampling sketch estimating `F0` and `f1..f64` | Arbitrary | Approximate frequency spectrum only (no per-k-mer output) | O(1) memory (~MBs), single pass at I/O speed | Single/multi-stream cardinality only | Useful for fast preflight sizing of hash tables/disk reserves; cannot replace exact counting |
| **kmtricks** | Multi-sample raw-FASTQ k-mer matrix builder | Minimizer super-k-mer partitioning + Bloom matrices + cross-sample k-mer rescue | Typically `<= 64` | Binary presence or quantized/exact abundance matrices | Tunable RAM/disk via partition count; avoids full assembly per sample | Native multi-sample presence/abundance matrix (T3, T5) | Direct external competitor for T3 and Tier-2 raw-FASTQ matrix (T5) |
| **GGCAT** | Parallel compacted & colored dBG builder (Rust) | Bucketed super-k-mers, multi-threaded unitig compaction, run-length/bitvector color compression | Arbitrary (optimized `<= 64`) | Abundance threshold filtering (`-s`) | Low RAM and fast NVMe-backed buckets; native Rust implementation | Native colored compacted unitigs (T3, T4, T5) scaling to 100k+ bacterial genomes | Strongest architectural reference for T4 and Tier-2 colored unitigs beyond 4,096 genomes |
| **Cuttlefish 2** | Ultra-low-RAM dBG compaction | Deterministic finite automaton (DFA) state transitions over BBHash MPHF of k-mers | Typically `<= 64` | Reference mode (`--ref`) or read mode via KMC3 | Very low RAM per distinct k-mer (~3-4 B/k-mer in DFA state table) | Uncolored unitigs (T4) | Strong memory efficiency reference for T4 unitig compaction |
| **Fulgor** | Compressed colored unitig index | SSHash SPSS k-mer dictionary + hybrid color-set encoding (list, dense bitvector, runs, Elias-Fano) | `<= 32` (SSHash fast path) | Binary color presence (extensions for weights) | ~1-2 bytes per colored k-mer including index + color classes | Directly mirrors GRM's pattern-collapsed `BitMatrix` at 10k-1M genomes | Best algorithmic blueprint for replacing GRM's uncompressed `patterns.bin` + `kmers.bin` at >4,096 genomes |
| **Themisto** | Bacterial colored k-mer index & pseudoaligner | Spectral Burrows-Wheeler Transform (SBWT) + hybrid sparse/dense color-set compression | `<= 32` | Binary color presence | Compact SBWT (~5 bits/k-mer) + compressed color matrix | Native bacterial pan-genome color sets (100k+ genomes) | Proven at 100k+ bacterial isolate scale; strong benchmark for colored presence queries |
| **MetaGraph** | Petabase-scale annotated dBG | Succinct BOSS de Bruijn graph + row-diff / multi-BRWT compressed annotation matrices | Arbitrary | Supports binary presence, counts, and coordinates | Highest compression ratios via row-diff / multi-BRWT; heavier index build cost | Supports arbitrary-scale color + abundance matrices | Reference architecture for massive cross-species repositories with count annotations |
| **SSHash** | Minimal-space exact k-mer dictionary | Sparse and skew hashing over minimizer-indexed unitig spectra (SPSS) | `<= 32` | Dictionary mapping k-mer -> unitig/position (plus optional weights) | ~6-8 bits/k-mer (vs. GRM's 8 B/k-mer in `kmers.bin` and `Table`) | Foundation for Fulgor color indexing | Directly applicable to shrinking GRM's 8 B/k-mer `kmers.bin` and `input_build::Table` by ~8x |

### Architectural synthesis

- **Where GRM's pipeline is already well matched to its workload:** For curated bacterial assembly cohorts up to 4,096 genomes (2,000 in Revision 14), GRM's in-memory LSD radix sort + invertible SplitMix64 mixing + 1,024-bucket merge achieves ~108 s end-to-end (`count` + `matrix`) on 2,000 *E. coli* genomes without external binaries, delivers exact bottom-1,000 MinHash sketches for free, outputs an ML-ready row-major `BitMatrix` directly for AVX2 popcount learners, and compacts unitigs 3.2x faster at 3.2x lower RAM than `bcalm2 2.2.3`.
- **Where GRM lags external systems:**
  1. **Raw FASTQ reads (T2, T5):** GRM pushes every k-mer occurrence (including 80-90% singleton sequencing errors at 30x-50x coverage) into a flat `Vec<u64>` before sorting (`src/kset.rs:50-54`), causing `O(coverage * genome_len)` RAM spikes instead of using super-k-mer minimizer binning + Bloom singleton rejection (as in FastK, KMC3, kmtricks, GGCAT).
  2. **Pattern-table allocation overhead during merge & projection:** GRM allocates a heap `Box<[u64]>` on every lookup in `BucketMerge`, `assemble`, `Dataset::project`, and `learn::training_candidates`.
  3. **Scale past `REVIEWED_COLUMNS = 4096`:** Uncompressed `patterns.bin` (`n_patterns * ceil(n/64) * 8` B) and `kmers.bin` (8 B/k-mer) scale linearly with `n * n_patterns` and do not exploit unitig spectrum compression (SSHash/GGCAT) or color-class compression (Fulgor/Themisto/MetaGraph).

## 4. Tier-1 optimization roadmap (Revision-14 bit-exact Rust changes)

Every Tier-1 optimization below preserves byte-for-byte Revision 14 `.kset`, `matrix/`, `structure.json`, and `learn` outputs (verifiable against `cargo test --release --locked`, `scripts/learner-crosscheck/check_fits.sh`, and pinned matrix SHA-256 digests in `src/matrix.rs:2342-2373`).

| Priority | Target (`file:lines`) | Current Bottleneck | Proposed Change | Expected Impact |
| --- | --- | --- | --- | --- |
| **1.1** | `src/matrix.rs:131,165`, `309,333`, `1958,1975`, `1986,1992`; `src/learn/mod.rs:458,468` | `FxHashMap<Box<[u64]>, u32>` / `FxHashSet<Box<[u64]>>` allocates a heap `Box<[u64]>` on *every* lookup (`entry(row.into())`, `seen.insert(masked.as_slice().into())`), even when the pattern is already present (>95% hit rate), and keeps a duplicate copy of every pattern row in RAM. | Replace with an index-only hash table (`hashbrown::HashTable<u32>` or flat open-addressed `Vec<u32>`) hashing `&[u64]` slices directly against the flat arena (`out.rows`, `patterns.data`, `sub.data`). Zero allocations on cache hits; single copy of row data. | Cuts pattern-dedup heap allocations by >95%, reduces `matrix` and `Dataset::project` peak RAM by ~30-45%, and speeds up `assemble`, `Dataset::project`, and `training_candidates`. |
| **1.2** | `src/matrix.rs:473-485` | `MergeKey::Wide` sorts `Vec<(u64, u32)>` using comparison `entries.sort_unstable()` (`src/matrix.rs:481`), costing `O(T_b log T_b)` (1.5x to 1.7x slower per entry than `MergeKey::Packed`, `docs/matrix-scaling.md:637`). | Replace with stable 8-bit LSD radix sort over the `u64` k-mer key alone. Because `fill(g)` inserts genomes in ascending order `0..plan.n` (`src/matrix.rs:476-480`), a stable LSD radix sort on the `2k - 10` key bits (7 byte passes at `k = 31`, 8 at `k = 32`) preserves ascending genome order per k-mer automatically (keeping repeated-genome checks in `group_sorted` exact). | Eliminates the 1.5x-1.7x penalty of `MergeKey::Wide`, giving linear `O(T_b)` merge time for `k = 32` past 1,024 genomes and `k = 31` past 4,096 genomes. |
| **1.3** | `src/fasta.rs:11-24`, `50-62`, `105-117`; `src/kset.rs:48-56` | (a) `fasta::open` opens every file twice (`src/fasta.rs:12,15`) to check 2 gzip magic bytes; (b) FASTQ loop allocates `let mut skip = Vec::new()` per read (`src/fasta.rs:57`); (c) `count_into` scans `seq` twice: once in `stats.add(seq)` (`src/kset.rs:51`) and once in `for_each_canonical` (`src/kset.rs:52`). | (a) Single `File::open` with `Read::chain` over the 2 peeked bytes (or `seek(Start(0))`); (b) hoist `skip` buffer outside the read loop; (c) fuse `StatsBuilder` base counting into the `CODE[b]` lookup loop during canonical k-mer extraction. | Halves per-file `open()` syscalls, eliminates millions of per-read heap allocations on FASTQ, and halves L1/DRAM scan traffic during `grm count`. |
| **1.4** | `src/dna.rs:133-141` | `KmerSpace::reverse_complement` runs a `self.k`-iteration scalar loop (31 iterations per call). Called 8 times per distinct 31-mer during `input_build::compact` (`src/input_build.rs:415`, ~836M calls on *E. coli*). | Replace with branchless O(1) SWAR bit-pair reversal: `let c = (!kmer).swap_bytes();` + 4-bit nibble swap (`0x0F0F...`) + 2-bit pair swap (`0x3333...`) + `c >> (64 - self.bits)`. Exact for all `k` in `5..=32` (verified in `scripts/kmer-bench/parity.py`). | ~5x-8x faster `reverse_complement`, meaningfully speeding up `input_build::compact` and `read_unitigs` with ~6 ALU instructions instead of 31 loop iterations. |
| **1.5** | `src/input_build.rs:423` | `successors` computes `let mask = (1u64 << space.bits()) - 1;`. When `k = 32`, `space.bits() == 64`, causing shift overflow (panic in debug mode, `mask = 0` in release mode on x86-64). | Expose `KmerSpace::mask(&self) -> u64` from `src/dna.rs:41` (which already guards `bits == 64 => u64::MAX` at `src/dna.rs:54-58`) and use `space.mask()` in `input_build.rs:423`. | Fixes latent `k = 32` correctness bug in unitig compaction with zero runtime cost. |
| **1.6** | `src/cluster.rs:43-60` | `distance_matrix` computes `mash_distance(jaccard(&sketches[i], &sketches[j], s), k)` for all `(i, j)` pairs including both upper and lower triangles (`n^2 - n` evaluations). | Compute only `i < j` in parallel (using triangular indexing or chunked upper rows) and mirror `dist[j * n + i] = dist[i * n + j]`. | Exact 2x speedup in pairwise Mash distance computation for `grm structure` (important at 4k-10k genomes where `n^2` grows quadratically). |
| **1.7** | `src/kset.rs:202-221`, `351-383` | `.kset` v1 (`GRMKSET1`) uses byte-aligned LEB128 gaps across 1,024 buckets (~5.9-6.0 B/k-mer at 5M 31-mers because mean gap `~2^39.7` requires 6 bytes = 48 bits). | Introduce opt-in `.kset` v2 (`GRMKSET2`) with per-bucket Golomb-Rice (or prefix-stripped Elias-Fano) gap coding at `b = floor(log2(2^(2k - BUCKET_BITS) / max(1, count)))`, keeping `GRMKSET1` readable for backward compatibility and preserving `KsetReader::sketch` and `kset::digest`. | Reduces `.kset` disk footprint and `matrix` read I/O by ~13-15% (~5.1-5.2 B/k-mer) while retaining random bucket access. |

## 5. Tier-2 roadmap (scientific and beyond-4,096 scale capabilities)

Tier-2 items extend GRM's scientific scope or scale beyond Revision 14's contracts and therefore require explicit evaluation plans and downstream review before production adoption.

### 5.1 Raw FASTQ streaming and plasmid copy-number normalization

To test whether real sequencing dosage overcomes the C3 assembly-count null result (0/5 learners carried for assembly `[count >= 2]`, Section 2.4):

1. **Streaming FASTQ counting with singleton rejection:** Partition incoming reads into `B` minimizer super-k-mer buckets (FastK/KMC3 design) or filter singletons via a compact Bloom/quotient filter before sorting, preventing 30x-50x sequencing error clouds (typically 80-90% of distinct raw k-mers) from dominating RAM. Retain `u16` saturating counts (`1..=65535`) for k-mers meeting `min_abundance >= 2` (or adaptive valley cutoff above the error peak).
2. **Chromosomal depth normalization:** Estimate each isolate's haploid chromosomal depth `c_chrom` from the robust median/mode of single-copy core k-mer counts (those present across core species patterns). Normalize every k-mer's abundance to relative copy number `r_i = c_i / c_chrom`, and emit calibrated binary dosage indicators (for example `r_i >= 2.5` for low-copy plasmid/tandem duplication and `r_i >= 8.0` for high-copy plasmid amplification) alongside presence features so learners can distinguish high-dosage resistance alleles without confounding by library sequencing depth.

### 5.2 Colored compacted unitigs beyond `REVIEWED_COLUMNS = 4096`

Scaling past `REVIEWED_COLUMNS = 4096` (`src/matrix.rs:46`) to 10,000-100,000+ genomes hits two bottlenecks: uncompressed `patterns.bin` (`n_patterns * ceil(n/64) * 8` bytes) and `kmers.bin` (8 bytes per kept k-mer).

- **Graph + color representation:** Adopt a colored compacted de Bruijn graph layout modeled on **GGCAT** (construction) and **Fulgor** / **Themisto** (indexing): store distinct k-mers as an **SSHash** unitig spectrum (~6-8 bits/k-mer, ~8x smaller than `kmers.bin`) and compress pattern bitvectors using hybrid color-class encodings (sorted sparse lists for rare patterns, run-length/Elias-Fano or bitvectors for dense patterns).
- **Downstream gate:** Any lift of `REVIEWED_COLUMNS` still requires an explicit end-to-end review of `Dataset::project`, learners, annotation, and report schemas (`docs/matrix-scaling.md:724-727`).

### 5.3 Support for `k > 32`

Extending `KmerSpace` (`src/dna.rs:18`) from `u64` (`k <= 32`) to a `u128` specialization (`33 <= k <= 64`) enables benchmarking longer sequence contexts (flanking mobile elements, promoter mutations) while keeping the `u64` path untouched for `k <= 32`.

## 6. Reproducible benchmark suite (`scripts/kmer-bench/`)

`scripts/kmer-bench/` provides a zero-dependency Python 3.12 benchmark harness covering five workloads, exact reference parity verification, `/usr/bin/time -v` telemetry parsing, disk-footprint sampling, binary SHA-256 provenance capture, and zero-genome-ID hygiene checks (`grm-kmer-bench-v1`).

### 6.1 Benchmark workloads

| Workload | Input | Task | Tools Covered | Parity Check |
| --- | --- | --- | --- | --- |
| **T1** | Multi-contig FASTA assemblies | Per-genome canonical k-mer counting (`k = 31`, `min_abundance = 1`) | `grm`, `fastk`, `kmc3`, `jellyfish`, `ntcard` | Exact canonical set + FNV-1a64 (`kset::digest`) + SHA-256 + abundance spectrum (`ntcard` marked cardinality-only) |
| **T2** | Paired-end noisy FASTQ reads (with multi-copy plasmids) | Raw read counting and abundance recovery (`k = 31`, `min_abundance >= 2`) | `grm`, `fastk`, `kmc3`, `jellyfish`, `ntcard` | Exact canonical set + abundance spectrum (with optional saturation cap, e.g. 255 for FastK `u8`) + plasmid depth ratio |
| **T3** | `N` FASTA assemblies | Multi-genome presence matrix (`min_support = 2`, stock / `--column-blocks` / `--disk-sparse` variants) | `grm`, `kmtricks`, `ggcat` | Exact per-k-mer presence bitset and collapsed pattern count |
| **T4** | `N` FASTA assemblies | Compacted de Bruijn graph unitigs over species union (`min_abundance = 1`) | `grm`, `bcalm2`, `ggcat`, `cuttlefish2` | Exact set-of-sets unitig equality (`src/input_build.rs:575-626`) |
| **T5** | `N` paired-end FASTQ read sets | Direct-from-FASTQ multi-genome presence matrix (`min_abundance >= 2`, `min_support = 2`) | `grm`, `kmtricks`, `ggcat` | Exact per-k-mer presence bitset over `min_abundance >= 2` k-mers |

### 6.2 Metrics and output hygiene (`grm-kmer-bench-v1`)

Every benchmark run writes `bench.json` (`format: "grm-kmer-bench-v1"`) and `summary.tsv` containing:

- **Compute metrics (`/usr/bin/time -v`):** wall-clock seconds, user seconds, system seconds, peak resident set size (`max_rss_kb`), exit status, and optional CPU percentage / page faults / filesystem I/O counters.
- **Disk metrics (`DiskSampler`):** `peak_temp_bytes`, `final_artifact_bytes`, `peak_total_bytes`, `kmers_counted`, and `bytes_per_kmer`.
- **Provenance:** binary basename, binary SHA-256 digest, and version line (with non-neutral home paths redacted).
- **Hygiene:** aggregate cohort counts only (zero genome IDs, zero Windows/private home paths, zero em dashes), checked against `scripts/c7-outside/scan.py` token rules before output is written.

### 6.3 Exact CLI commands (local WSL2 / Linux and Lambda runs)

1. **Generate a deterministic synthetic bacterial cohort (FASTA + noisy paired-end FASTQ + ground truth):**

```bash
python3 -I -B scripts/kmer-bench/cli.py generate \
  --out work/kmer-bench/cohort-small \
  --seed 42 --n-genomes 16 --k 31 \
  --core-len 50000 --n-plasmids 2 --plasmid-copies 4,8 \
  --coverage 20 --read-length 150 --insert-size 300 --error-rate 0.005
```

2. **Preview command plan or run deterministic CI mock suite across T1-T5:**

```bash
python3 -I -B scripts/kmer-bench/cli.py plan \
  --cohort-dir work/kmer-bench/cohort-small \
  --out-dir work/kmer-bench/plan-out \
  --workloads all --tools all --include-grm-variants

python3 -I -B scripts/kmer-bench/cli.py run \
  --cohort-dir work/kmer-bench/cohort-small \
  --out-dir work/kmer-bench/mock-out \
  --workloads all --tools all --include-grm-variants --mock
```

3. **Run real benchmarks on a Linux/Lambda host (gracefully skipping optional uninstalled binaries):**

```bash
python3 -I -B scripts/kmer-bench/cli.py run \
  --cohort-dir work/kmer-bench/cohort-large \
  --out-dir work/kmer-bench/real-out \
  --workloads T1,T2,T3,T4,T5 \
  --tools grm,fastk,kmc3,jellyfish,ntcard,kmtricks,ggcat,cuttlefish2,bcalm2 \
  --threads 16 --mem-gb 32 --k 31 --include-grm-variants
```

4. **Verify exact parity on any candidate `.kset`, text k-mer dump, or unitigs FASTA:**

```bash
python3 -I -B scripts/kmer-bench/cli.py verify-parity \
  --mode set --k 31 --min-abundance 1 \
  --reference work/kmer-bench/cohort-small/assemblies/syn_g0000.fna \
  --candidate-kset work/kmer-bench/real-out/work/T1_grm_default/data/species/synthetic_sp/ksets/syn_g0000.kset

python3 -I -B scripts/kmer-bench/cli.py verify-parity \
  --mode abundance --k 31 --min-abundance 2 --max-cap 255 \
  --reference work/kmer-bench/cohort-small/reads/syn_g0000_R1.fastq \
              work/kmer-bench/cohort-small/reads/syn_g0000_R2.fastq \
  --candidate-dump work/kmer-bench/real-out/fastk_tabex.tsv
```

5. **Aggregate one or more `bench.json` reports into combined TSV and Markdown tables:**

```bash
python3 -I -B scripts/kmer-bench/cli.py summarize \
  --inputs work/kmer-bench/real-out/bench.json \
  --tsv-out work/kmer-bench/combined.tsv \
  --md-out work/kmer-bench/combined.md
```
