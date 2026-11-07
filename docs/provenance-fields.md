# Source, year and measurement fields from fetch to report

GRM fetches the source, year and measurement fields of every BV-BRC AST record and genome, and most of them stop at `grm curate`. Only the consolidated MIC and the genome's collection year reach a dataset report. The method and standard reach `corpus.json` only as corpus-wide counts of the raw records. Work-program row A06 records this as "source/year/measurement provenance is not retained and used throughout" (`docs/work-program-20261027.md:133`). Its first step adds side files that carry them, unchanged and as provenance only, to `report/provenance/` ("The side files", below); the trace up to there describes the pinned files, which that step leaves as they were.

The line numbers below are those of `b90be1ca`, read from the record types, the code that reads them and the synthetic test fixtures. This trace proposes no new field and no use of any field as a predictor, and H's proxy refusal stays in force (see the last section).

## The path

| Step | Reads | Writes |
|---|---|---|
| `grm fetch amr`, `grm fetch genomes` | BV-BRC `genome_amr` (laboratory-method records only) and `genome` | `bvbrc/genome_amr.jsonl`, `bvbrc/genomes.jsonl` |
| `grm curate` | both raw files | `curated/phenotypes.jsonl`, `datasets.json`, `funnel.json` |
| `grm learn` | `phenotypes.jsonl`, `datasets.json`, `genomes.jsonl` | `species/<slug>/results/<id>.json`, one report per dataset |
| `grm report` | the curated files and both raw files | `report/index.json`, `corpus.json`, `species/<slug>.json`, and `datasets/<id>.json` copied as `grm learn` wrote it |

Fetch is `src/bvbrc.rs:229-251` and `:255-275`, and the paths are `src/layout.rs:15-21`, `:31-33` and `:57-59`. Curate reads at `src/main.rs:524-525` and writes at `:578-590`. Learn reads at `src/analysis.rs:331-341` and writes at `:941-952`. Report reads at `src/report.rs:524-535`, copies at `:570-579` and writes `corpus.json` at `:779`.

The raw files hold exactly the fields of `AmrRecord` (`src/bvbrc.rs:22-59`) and `GenomeRecord` (`:65-120`). Fetch selects those fields (`AMR_FIELDS`, `:61-63`; `GENOME_FIELDS`, `:122-125`) and writes each record back through its type (`:237-240`, `:264-267`). A number becomes text and a list is joined with `;` (`:402-417`). The fixture `parses_mixed_type_fields` (`:433-443`) stores a `testing_standard_year` of 2021 as `"2021"` and a `pmid` of `[26686880,1]` as `"26686880;1"`.

## The fields

| Kind | Field | Record, `src/bvbrc.rs` | Last kept in | Where it stops |
|---|---|---|---|---|
| Source | `source` | AMR, `:54` | raw file | no code reads it |
| Source | `pmid` | AMR, `:56` | raw file | no code reads it |
| Source | `isolation_source` | genome, `:103` | raw file | no code reads it |
| Year | `testing_standard_year` | AMR, `:52` | raw file | no code reads it |
| Year | `date_inserted` | AMR, `:58` | raw file | no code reads it |
| Year | `collection_year` | genome, `:99` | report | reaches the report, per genome |
| Measurement | `measurement_value`, `measurement` | AMR, `:40`, `:36` | report, as one log2 step per pair | `grm curate` rounds it and keeps one per pair |
| Measurement | `measurement_sign` | AMR, `:38` | report, as the MIC's censor | `grm curate` turns it into the censor |
| Measurement | `measurement_unit` | AMR, `:42` | raw file | `grm curate` checks it, then drops it |
| Measurement | `laboratory_typing_platform` | AMR, `:46` | raw file | `grm curate` checks it for E-test, then drops it |
| Measurement | `vendor` | AMR, `:48` | raw file | no code reads it |
| Measurement | `laboratory_typing_method` | AMR, `:44` | `phenotypes.jsonl`, as a set per pair | `grm learn` and `grm report` do not read the set |
| Measurement | `testing_standard` | AMR, `:50` | `phenotypes.jsonl`, as an uppercased set per pair | `grm learn` and `grm report` do not read the set |

The genome record's other metadata (country, host, MLST, sequencing platform and the assembly fields) is outside this trace.

### Stopped at curation

`curate` reads each AMR record's `genome_id`, `antibiotic` and `resistant_phenotype`, its measurement fields through `parse_mic`, and its `laboratory_typing_method` and `testing_standard` (`src/curate.rs:421-455`). From the genome record it reads the species, name, length, quality, contig count and CheckM values (`:310-339`). `Phenotype`, the one row per genome and drug that `grm curate` writes, has no field for anything else (`:201-220`).

So `source`, `pmid`, `testing_standard_year`, `date_inserted`, `vendor` and `isolation_source` never leave the raw files. No other code in `src/` reads them, and the fixture above is the only other use of `testing_standard_year` and `pmid`.

### Measurement fields, reduced to one MIC

`parse_mic` (`src/curate.rs:58-119`) turns one record into an MIC or nothing. The unit must be `mg/L`, or absent under a method that is not disk diffusion, and under disk diffusion only an E-test platform passes (`:59-79`). The value is `measurement_value`, else `measurement`. The second component of a combination such as `2/38` is dropped, and the value is rounded to the nearest log2 step (`:80-96`). `measurement_sign`, else a sign written in the value, becomes the censor: exact, at most or above (`:84-89`, `:97-118`). The unit and the platform are not kept.

The fixture `parses_mics_and_rejects_zone_diameters` (`:968-1006`) shows each step: `<=` 0.25 mg/L becomes step -2 at most, `>2/38` becomes step 1 above, 0.06 becomes step -4 exact (0.0625 mg/L), a zone diameter gives nothing, and an E-test strip under disk diffusion gives an MIC.

Each genome and drug pair then keeps one MIC: the median exact step, or else the tightest bound when every censored value points the same way (`consolidate_mics`, `:577-621`). It also keeps `mic_spread`, the range of the exact steps (`:470`, `:585`), and `n_records`, the number of records (`:449`, `:483`). The record-level values stop here.

The pair also keeps every `laboratory_typing_method` and every uppercased `testing_standard` its records had, as sets (`:452-454`, `:481-482`), and a set does not record which call or MIC each one came with.

### Written to `phenotypes.jsonl`, then unread

`grm curate` writes `methods`, `standards`, `mic_spread` and `n_records` into `phenotypes.jsonl` (`src/main.rs:578-582`). `grm learn` and `grm report` read that file (`src/analysis.rs:331-335`, `src/report.rs:529`) but use none of the four. No code in `src/` outside `curate.rs` names them except a test fixture (`src/analysis.rs:1965-1970`), and no Python under `scripts/` or `tests/` names them.

### What reaches the report

- **MIC.** Each sample's consolidated MIC, as a log2 step and censor, is in the dataset report's `samples` (`src/analysis.rs:681`), which `grm report` copies unchanged (`src/report.rs:570-579`). The MIC regression reads the same value (`src/analysis.rs:516`).
- **Collection year.** `collection_year`, as BV-BRC stored it, is in each dataset report's `samples[].year` (`src/analysis.rs:684`), in each species file's `genomes[].year` (`src/report.rs:417`, `:486`), and in `corpus.json`'s `years`, counted over genomes with a label (`src/report.rs:731-741`, `:772`). It is a genome field and is not linked to any AST record.
- **Method and standard, corpus-wide.** `corpus.json`'s `methods` and `standards` count `laboratory_typing_method` and `testing_standard` over every raw AMR record (`src/report.rs:535`, `:748-757`), including the records `grm curate` set aside. They are not per dataset or per label, and the standards are counted as submitted, not uppercased.

`grm learn` reads the genome metadata only to fill each `Sample` row (`src/analysis.rs:676-686`), and `input.genomes` has no other use in `src/analysis.rs`. Its input fingerprint leaves the metadata out because it changes "only descriptive report fields" (`src/analysis.rs:1074-1075`).

## The side files (register row A06, step 1)

Three side files now carry these fields to the report directory. No pinned file changes: `genomes.jsonl`, `GENOME_FIELDS`, `Phenotype` and every curated file `grm curate` wrote before are byte for byte as they were, and no dataset report changes. No side file enters an input fingerprint, a feature block or a learner.

| Step | Writes | Pinned by |
|---|---|---|
| `grm fetch genome-dates` | `bvbrc/genome_dates.jsonl`: one row per genome of `species/*/genomes.txt` that BV-BRC returns, in genome id order | `bvbrc/genome_dates.manifest.json`: each list's SHA-256, the requested, returned and missing genomes, and the file's size and SHA-256 |
| `grm curate` | `curated/provenance.jsonl`: one row per genome and drug, beside each `phenotypes.jsonl` row and in the same order | the SHA-256 that `grm provenance` records |
| `grm provenance` | `report/provenance/<id>.json` for each dataset `index.json` lists: the report's samples in its order, each with its genome's dates and its pair's row | the report's SHA-256 and evaluation revision, and both inputs' SHA-256; `grm provenance --check` fails when any differs or a file was edited |

The code is `src/provenance.rs`, `curate::PairProvenance` and `bvbrc::GenomeDate`.

### What each field reaches now

| Field | From | Side file field | Reaches | Proxy check by name |
|---|---|---|---|---|
| `date_inserted` (genome) | BV-BRC `genome` | `genome_dates.jsonl` `date_inserted` | `report/provenance`, per sample | refused (insertion date) |
| `completion_date` (genome) | BV-BRC `genome` | `genome_dates.jsonl` `completion_date` | `report/provenance`, per sample | refused (date) |
| `date_inserted` (AST record) | `genome_amr.jsonl` | `date_inserted_first`, `date_inserted_last`: the earliest and latest over the pair's records | `report/provenance`, per sample | refused (insertion date) |
| `testing_standard_year` | `genome_amr.jsonl` | `testing_standard_years`, a set | `report/provenance`, per sample | refused (year) |
| `pmid` | `genome_amr.jsonl` | `pmids`, a set, a list split at `;` | `report/provenance`, per sample | refused (record field) |
| `source` | `genome_amr.jsonl` | `sources`, a set | `report/provenance`, per sample | refused (record field) |
| `laboratory_typing_method` | `genome_amr.jsonl` | `methods`, as in `phenotypes.jsonl` | `report/provenance`, per sample | refused (record field) |
| `testing_standard` | `genome_amr.jsonl` | `standards`, as in `phenotypes.jsonl` | `report/provenance`, per sample | refused (record field) |
| `n_records`, `mic_spread` | `grm curate` | as in `phenotypes.jsonl` | `report/provenance`, per sample | refused (record field) |
| `vendor`, `isolation_source`, `measurement_unit`, `laboratory_typing_platform` | raw files | none | raw file only, as before | |

A pair's row covers the records its phenotype consolidates: the records `grm curate` sets aside (an unresolved drug name, a genome that fails QC, a record with neither a call nor an MIC) are in neither. The record dates are ordered by time, the whole seconds and then the fraction, since BV-BRC writes some without milliseconds and `13.5Z` sorts before `13Z` as text. They stay as BV-BRC wrote them.

A sample whose genome has no row in the dates file has `genome_dates` null, and one with no curated pair has `pair` null; each file counts both. A genome BV-BRC returns without a date has a row with that date null, and the manifest counts those too.

### The proxy check

`ablate::blocks::admit` refuses a block with any side-file field as a column, block id or version name. The dates and years match the date rules. `pmid`, `pmids`, `source`, `sources`, `methods`, `standards`, `n_records` and `mic_spread` are an AST record's study and lab-method fields (`RECORD_FIELDS`, kind "study or lab-method field of an AST record"). They are refused only as a whole name, or as the part after a block prefix (`C:pmid`), ignoring case, so a genotype name that merely contains one, such as `P:resource`, still passes.

- `admit_refuses_a_block_of_the_provenance_side_file_fields` names each side-file field with the rule that refuses it, so the test fails when a side file gains a field until the field has a rule.
- `admit_refuses_record_fields_by_exact_name_only` checks each record field as a column and as a version, and genotype-shaped names that carry one inside a longer name.
- `the_real_catalog_s_gene_groups_pass_the_deny_list` (ignored; it needs the catalog in `GRM_CATALOG_DIR`) runs the deny-list over every gene group of a real AMRFinderPlus catalog and every panel name. No name in catalog 2026-08-07.1 contains "source".

Neither side file reaches a block in any case: no adapter reads them, and `admit` accepts only this module's adapters.

### Kept out of copies and hashes

What reads, ships or hashes data, and where the side files stay out:
- `grm learn`, `grm report` and every input fingerprint read named curated files only, never `curated/provenance.jsonl` or the genome dates.
- The release and publish scripts (`scripts/deployment/`, `ops/cloudflare/`) never read `curated/`, but they publish exactly what a tree's `SOURCE.sha256` lists. `report/provenance/` is not part of a publication, so a release `SOURCE.sha256` must leave it out, even when `grm provenance` ran in the same root (listing all of `report/`, as for the intervals, would take it in).
- The box kit's `stage-manifest.py` (outside this repo) stages every `curated/*` line of a master manifest, so it sends `curated/provenance.jsonl` to a box until main installs the version that skips it (`stage-manifest.py.next`, sha256 `68cf7b4c…`, from lane PROV). A push-data manifest built by hand must leave out `curated/provenance.jsonl` and `bvbrc/genome_dates*`: the runner's allowlist accepts any `curated/*.jsonl` line.
- `scripts/input-comparison/run.py` removes `curated/provenance.jsonl` from each data root right after it curates (about 310 MB per root on the real root, and no arm reads it). The species pool never copies it, and the disk estimate never counts it. RAW's input manifest leaves out `bvbrc/genome_dates.jsonl` and its manifest, so fetching the dates into RAW doesn't refuse a WORK made before.

## H's proxy refusal stays in force

This trace proposes no field and no predictor. H, the ablation task that wrote the feature-block registry (`docs/ablation-runner.md:3`, `src/ablate/plan.rs:654`; not block H, the hierarchical panel calls), refuses source proxies in `src/ablate/blocks.rs`, and that refusal is unchanged:

> Source proxies. Lineage and cluster ids, BioProject, study, submitter or centre, country, collection and insertion dates, sequencing platform and assembly QC say where a genome came from, not what it carries, and a model can learn resistance prevalence from them. `proxy_kind` names them, and `admit` refuses any block whose id, version or column names match one. They can only ever be reported as a named diagnostic.

That is `src/ablate/blocks.rs:20-25`. `check_spec` applies it to a block's id and version (`:376`) and `check_block` to its column names (`:903`). `admit` (`:843-860`) is the only way to make a predictor block, and its sealed bound accepts only that module's adapters (`:30-32`).

Block X, genome context, stays unavailable because "the context fields in genome metadata are source proxies and are refused" (`:205`, `:253`). The deny-list's test names `collection_year`, `Date Inserted` and `sequencing_platform` among the refused names (`:1750`, `:1794-1802`). The module doc states the limits of a name check and names the seal on `admit` as the real guard (`:44-48`).
