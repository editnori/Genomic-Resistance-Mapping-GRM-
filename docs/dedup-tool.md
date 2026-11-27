# Cross-database deduplication tool (`grm dedup`)

`grm dedup` builds compact versioned reference indices (`grm-dedup-index-v1`) and screens external or newly acquired genome batches against one or more reference indices, while performing within-batch connected-component deduplication and longitudinal patient grouping.

## Epistemic scope and limits

Every screen summary written by `grm dedup screen` carries two explicit status fields:

- `claim_scope`: `"zero_detected_overlap_under_stated_checks"`
- `sealed_set_status`: `"unqueried_unknown"`

Deduplication cannot prove absolute absence of overlap across all public and private databases. In particular:

- **Unopened sealed records are never queried by default.** Screening checks only the reference indices explicitly passed via `--reference`.
- **Accession gaps make metadata overlap counts lower bounds.** Reference and query isolates without assembly, BioSample, or SRA accessions cannot be linked by accession alone (`ACCESSION_GAP` is counted for both references and queries).
- **Missing sequences prevent clean certification.** Queries without FASTA or sketch fingerprints cannot be checked at the sequence level and are placed in `UNRESOLVED_OVERLAP` rather than `certified_clean`.
- **Study overlap is not isolate identity.** Sharing a BioProject (`SHARED_BIOPROJECT`) or publication (`SHARED_PMID`) is tagged separately and never conflated with isolate duplication.
- **Close relatedness is not same-specimen identity.** Genomic near-duplicates (`<= near_dup_threshold`) are separated from close relatives (`(near_dup_threshold, close_relative_threshold]`), and thresholds should be validated per species on permitted development records rather than fixed blindly across species.

## Streaming FASTA fingerprinting

`grm::kset::fingerprint_fasta(path, k_digest, k_sketch, sketch_size)` reads plain or gzip-compressed FASTA assemblies directly in memory without writing `.kset` files, computing:

| Field | Definition |
| --- | --- |
| `fnv64_k31` | 16-char lowercase hex FNV-1a64 over sorted mixed canonical `k_digest`-mers (default `k_digest = 31`), matching `kset::digest` and `counted.json` `kset_digest`. |
| `sha256_k31` | 64-char lowercase hex SHA-256 over sorted mixed canonical `k_digest`-mers encoded as little-endian `u64` bytes, eliminating 64-bit FNV collision risk. |
| `assembly_sha256` | 64-char lowercase hex SHA-256 over strand-canonical uppercase `ACGTN` contig sequences sorted lexicographically and length-delimited (`GRM-ASM-V1`). Invariant to FASTA headers, line wrapping, base case, contig order, and contig strand, while staying distinct from canonical k-mer-set identity (for example, when a contig is duplicated or split). |
| `n_kmers`, `contigs`, `bp` | Distinct canonical `k_digest`-mer count, contig count, and assembly base-pair length. |
| `sketch_k21` | Bottom-`sketch_size` MinHash sketch at `k_sketch` (default `k_sketch = 21`, `sketch_size = 1000`), matching `species-check/sketches/*.json`. |
| `sketch_k31` | Optional bottom-`sketch_size` MinHash sketch at `k_digest`. |

## Building a reference index (`grm dedup index`)

`grm dedup index` builds a deterministic versioned index directory containing `records.jsonl` and `manifest.json` (`schema: "grm-dedup-index-v1"`), pinned by `records_sha256` and `index_sha256`. It reads metadata files strictly read-only and never modifies `bvbrc/genomes.jsonl`.

```sh
grm dedup index \
  --metadata data/bvbrc/genomes.jsonl \
  --provenance data/bvbrc/genome_provenance.jsonl \
  --counted data/species \
  --sketches-dir data/species-check/sketches \
  --source-label bvbrc-rev14 \
  --out dedup-ref/bvbrc-rev14
```

Optional inputs:

- `--metadata <PATH>` (repeatable): genome metadata JSONL files (for example `bvbrc/genomes.jsonl`).
- `--provenance <PATH>`: sidecar provenance JSONL (`genome_provenance.jsonl`) carrying `biosample_accession`, `assembly_accession`, `sra_accession`, `bioproject_accession`, `strain`, `pmid`/`publication`, and `subject_id` without altering `genomes.jsonl`.
- `--genome-ids <PATH>`: optional allowlist of genome ids to index.
- `--fasta-dir <DIR>`: directory of FASTA assemblies fingerprinted in memory via `kset::fingerprint_fasta`.
- `--counted <PATH>` (repeatable): `counted.json` file(s) or a `species/` directory holding `*/counted.json`. When both `--fasta-dir` and `--counted` are provided, `fnv64_k31` is checked against `counted.json` `kset_digest`.
- `--sketches-dir <DIR>`: existing `species-check/sketches/` directory of `<id>.json` sketches (`k = 21`, `s = 1000`).

On every `grm dedup screen` run, `load_index` verifies `manifest.json` schema, `records_bytes`, `records_sha256`, and `index_sha256` before screening.

## Screening a query batch (`grm dedup screen`)

```sh
grm dedup screen \
  --reference dedup-ref/bvbrc-rev14 \
  --query external/queries.jsonl \
  --query-fasta-dir external/fasta \
  --synonym-map external/synonyms.tsv \
  --thresholds external/thresholds.json \
  --out external/dedup-screen
```

### Identity and overlap layers

Each query record in `decisions.jsonl` reports every matching layer independently:

1. **Accession and alias normalization (`D1`, `D2`, `STRAIN`):**
   - Assembly keys (`D1`) strip and uppercase, and map `GCA_`/`GCF_` accessions to their numeric core (`ASM:<digits>`, ignoring prefix and `.version`).
   - BioSample and SRA/ENA/DDBJ keys (`D2`) split multi-valued `sra_accession` fields on `,`, `;`, `|`, and whitespace, and resolve alias equivalence classes loaded from `--synonym-map` (TSV or JSONL linking `SAMEA`<->`ERS`, `SAMN`<->`SRS`, `SAMD`<->`DRS`, run-to-sample links `ERR`/`SRR`/`DRR` <-> BioSample, and `PRJEB`<->`ERP` / `PRJNA`<->`SRP` study aliases).
   - Strain keys (`STRAIN`) uppercase with collapsed whitespace and match only within the same normalized two-word species; strain-only matches are counted separately in `layer_counts.STRAIN_ref_only`.
2. **Normalized assembly content (`EXACT_ASSEMBLY`):** exact `assembly_sha256` match against a reference record.
3. **Canonical k-mer-set identity (`EXACT_KSET`):** matches on `sha256_k31` when both query and reference carry it (rejecting any 64-bit FNV collision where `fnv64_k31` agrees but `sha256_k31` differs), or on `fnv64_k31` with `n_kmers` and optional `sketch_k31` consistency when one side came from `counted.json`.
4. **Genomic Mash distance (`NEAR_DUP` vs `CLOSE_RELATIVE`):**
   - `SketchIndex` builds an inverted posting index over the lowest `index_hashes` (default 64) values of each reference bottom-`s` sketch, filtering candidates before exact `cluster::mash_distance` verification rather than running dense pairwise comparisons.
   - `mash_distance <= near_dup_threshold` is recorded as `NEAR_DUP` (`D3_mash_near_dup`).
   - `near_dup_threshold < mash_distance <= close_relative_threshold` is recorded as `CLOSE_RELATIVE`.
5. **Within-batch connected components (`D4_within_batch_dup`):**
   - Surviving queries are linked via Union-Find across normalized/alias-resolved assembly accessions, BioSample/SRA accessions, exact `assembly_sha256`, exact k-mer digest, and within-batch Mash distance `<= min(near_dup_a, near_dup_b)`.
   - Each component elects the deterministic smallest `genome_id` under R0 `(taxon, suffix)` / lexicographical order as representative and excludes the remaining members (`components.json`).
6. **Longitudinal patient grouping (`subject_id`):**
   - Serial isolates sharing a normalized `subject_id` are grouped in `subject_groups.json` and share `split_group_id: "subject:<id>"` so downstream splits keep serial patient samples together without collapsing distinct serial isolates as `D4`.
   - A query sharing `subject_id` with a reference genome is flagged as `UNRESOLVED_OVERLAP` (`shared_subject_with_reference`).
7. **Unresolved overlap gating (`UNRESOLVED_OVERLAP`):**
   - Surviving representatives are excluded from `certified_clean` as `UNRESOLVED_OVERLAP` when they lack a sequence fingerprint (`missing_sequence_fingerprint`), share a reference `subject_id`, fall in the `CLOSE_RELATIVE` band with missing query accessions (`missing_accessions_with_close_relative`) or against an unaccessioned reference (`close_relative_to_unaccessioned_reference`), or have missing accessions against a metadata-only reference index.
8. **Study-overlap tagging (`SHARED_BIOPROJECT`, `SHARED_PMID`):**
   - Tagged independently on each decision (after `--synonym-map` study alias resolution) without excluding distinct isolates or conflating study overlap with isolate identity.

### Species-specific thresholds (`--thresholds`)

`--thresholds <json>` accepts a `grm-dedup-thresholds-v1` document mapping normalized species slugs (or two-word species names) to `(near_dup_threshold, close_relative_threshold)` pairs validated on permitted development records:

```json
{
  "schema": "grm-dedup-thresholds-v1",
  "default": {
    "near_dup_threshold": 0.0005,
    "close_relative_threshold": 0.001,
    "basis": "fallback when species override is absent"
  },
  "species": {
    "escherichia_coli": {
      "near_dup_threshold": 0.0005,
      "close_relative_threshold": 0.0015,
      "basis": "validated on permitted development records"
    },
    "mycobacterium_tuberculosis": {
      "near_dup_threshold": 0.0001,
      "close_relative_threshold": 0.0005,
      "basis": "clonal species; tighter threshold"
    }
  }
}
```

### Funnel codes and outputs

`grm dedup screen` writes four deterministic files under `--out`: `summary.json`, `decisions.jsonl`, `components.json`, and `subject_groups.json`.

In `summary.json`, the sequential exclusion stages partition `total_queries` exactly:

```text
total_queries
  = C4_ref_metadata
  + EXACT_ASSEMBLY
  + EXACT_KSET
  + D3_mash_near_dup
  + D4_within_batch_dup
  + UNRESOLVED_OVERLAP
  + certified_clean
```

| Funnel / tag code | Meaning |
| --- | --- |
| `total_queries` | Distinct query `genome_id`s screened. |
| `C4_ref_metadata` | Excluded at stage 1 by reference metadata (`D1`, `D2`, or same-species `STRAIN`). |
| `D1_ref_assembly` | Queries matching a reference assembly key (`GCA_`/`GCF_` core or assembly alias). |
| `D2_ref_biosample_sra` | Queries matching a reference BioSample, SRA run, or INSDC/ENA alias. |
| `EXACT_ASSEMBLY` | Remaining queries excluded at stage 2 by exact `assembly_sha256`. |
| `EXACT_KSET` | Remaining queries excluded at stage 3 by exact canonical k-mer-set digest (`sha256_k31` or verified `fnv64_k31`). |
| `D3_mash_near_dup` | Remaining queries excluded at stage 4 by genomic Mash distance `<= near_dup_threshold`. |
| `CLOSE_RELATIVE` | Non-exclusive count of queries with a reference hit in `(near_dup_threshold, close_relative_threshold]`. |
| `D4_within_batch_dup` | Non-representative members of within-batch connected components excluded at stage 5. |
| `UNRESOLVED_OVERLAP` | Remaining queries excluded at stage 6 because sequence checks were absent or close-relative/patient overlap could not be resolved. |
| `SHARED_BIOPROJECT` | Non-exclusive count of queries sharing a BioProject or study alias with a reference genome. |
| `SHARED_PMID` | Non-exclusive count of queries sharing a publication PMID with a reference genome. |
| `ACCESSION_GAP` | Non-exclusive count of queries lacking assembly, BioSample, and SRA accessions. |
| `certified_clean` | Queries surviving all exclusion stages (`claim_scope: "zero_detected_overlap_under_stated_checks"`). |

Marginal (non-exclusive) hit totals for every layer are reported separately in `summary.json` under `layer_counts`, and tags within `certified_clean` alone are reported under `certified_clean_tags`.

## External portal metadata normalizers (`scripts/dedup/portals.py`)

`scripts/dedup/portals.py` uses only the Python standard library to convert public AMR portal exports (TSV, CSV, or JSONL, plain or `.gz`, including multi-row per-antibiotic AST tables) into canonical `grm dedup` query JSONL and INSDC/ENA/NCBI synonym TSV files:

| `--portal` | Source portal | Extracted alias links |
| --- | --- | --- |
| `ebi_amr` (`cabbage`) | EMBL-EBI AMR Portal / CABBAGE | `SAMEA`<->`ERS`, `ERR`/`SRR`<->BioSample, `PRJEB`<->`ERP` |
| `cryptic` | CRyPTIC *Mycobacterium tuberculosis* compendium | `BIOSAMPLE`<->`ENA_SAMPLE`, `ENA_RUN`<->sample, `BIOPROJECT`<->`ENA_STUDY`, `SITEID:SUBJID` patient grouping |
| `ncbi_pd` (`narms`) | NARMS / NCBI Pathogen Detection | `SAMN`<->`SRS`, `SRR`<->`SAMN`, `GCA_`<->`GCF_`, `PRJNA`<->`SRP` |
| `enterobase` | EnteroBase strain/assembly tables | `Sample ID`<->`Secondary Sample Accession`, `Accession`<->sample, `Project ID`<->`Study Accession` |
| `pubmlst` | PubMLST isolate tables | `biosample_accession`<->`ena_sample`, `ena_run_accession`<->sample, `GCA_`<->`GCF_`, `bioproject_accession`<->`ena_study` |

Example conversion and synonym merging:

```sh
python3 -I scripts/dedup/portals.py \
  --portal ebi_amr \
  --input ebi_amr_export.tsv.gz \
  --out external/ebi_queries.jsonl \
  --synonym-out external/ebi_synonyms.tsv

python3 -I scripts/dedup/portals.py \
  --merge-synonyms external/ebi_synonyms.tsv external/ncbi_synonyms.tsv \
  --synonym-out external/all_synonyms.tsv
```

Unit tests for the portal normalizers live in `tests/dedup_portals`:

```sh
python3 -I -B -W error -m unittest discover -s tests/dedup_portals -t . -v
```
