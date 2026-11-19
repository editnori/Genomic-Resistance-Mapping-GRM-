# Frozen rules: candidate population for a retrospective test on newly deposited studies from the same database

These rules were written and hashed before any per-dataset class count was computed from the new labels. The hashed text has SHA-256 `7882b6f9dba4b81b4965709d15e67286bd8aee3d6890505e2cfda7d2ed525aeb` and was hashed at 2026-11-11T19:01:54Z. This copy differs from it in the title, in this paragraph and in the Inputs paragraph, where the locations of the local files are replaced by what the files hold. The rules R0 to R9 are unchanged. File names such as `new-genome-ids.txt` are files of the local cache of BV-BRC API responses the tables were built from; the cache is not in this repository. The tables the rules produced are beside this file: `candidates.tsv.gz`, `pairs.tsv.gz` and `datasets.tsv`, pinned by `../tables.sha256` on their uncompressed bytes.

Written and hashed (RULES.sha256) before any per-dataset class count was computed from the new labels.
Scope: genomes whose BV-BRC lab records were inserted 5 to 6 Nov 2026 (20,108 genome ids,
`bvbrc-check/new-genome-ids.txt`). This is a retrospective test on newly deposited studies from
the same database. It is not a lockbox: the batch's metadata and label totals were already read by
an earlier review. The insertion date is database ingestion time; it is neither collection time nor
the time a label became available. Collection fields are carried as cached and are not used as a filter.

Inputs (read only, a local cache of BV-BRC API responses and of the published site data):
new-batch metadata `bvbrc-check/new-genome-meta.jsonl.gz`; new-batch labelled records
`bvbrc-check/new-amr-labelled.jsonl.gz` (rows with a phenotype word only; MIC-only rows were never
fetched); held-FASTA metadata `bvbrc-check/master-genome-meta.jsonl.gz` (35,151 ids);
older-record ids `bvbrc-check/old-genome-ids.txt` and metadata `old-genome-meta.jsonl.gz`;
published drug-name table `live-check/corpus.json`; published dataset list `live-check/index.json`;
cohort genomes `live-check/species/*.json` and the 274 served reports
(`datasets/*.json` of the published report tree; sample ids and R/S per dataset).
Source rules copied from `src/curate.rs`, `src/analysis.rs` at commit aaa46ab.

## R0. Identifiers

- Genome id order: split `taxon.suffix` at the first dot and compare as the integer pair
  (taxon, suffix). "Smallest id" means smallest under this order.
- Accession normalisation: strip, upper case. `sra_accession` is split on `,` `;` and spaces;
  each run is its own key. `assembly_accession` keys on the digits of `GCA_`/`GCF_` (version
  dropped, GCA and GCF share a key). `strain` keys on upper case with runs of whitespace collapsed,
  and is matched only between genomes of the same species (R3).
- Cohort = the union of genome ids in `live-check/species/*.json` and of every served report's
  `samples[].id`. Held FASTA = the 35,151 ids of `master-fasta-ids.txt` / master metadata.
  Cohort metadata (accessions, studies) is read from the held-FASTA metadata.

## R1. Not new to the database

A new-batch id that already had a lab record inserted before 1 Jul 2025 (`old-genome-ids.txt`) is
dropped ("older lab record"). Reads no label.

## R2. Species scope (the cohort's own definition)

Species = the first two words of `species`, else of `genome_name`, when there are two words and the
first starts with an upper-case letter (copy of `species_of` in `src/curate.rs`). Slug as
`curate::slug`. No species -> dropped ("no species"). A species whose slug is not one of the
published species slugs in `index.json` -> dropped ("species outside the published 23").
Salmonella: the cohort defines Salmonella as the two-word species `Salmonella enterica`, which
includes typhoidal serovars (the cohort's own Salmonella holds Typhi and Paratyphi genomes). The
same definition is applied: typhoidal genomes stay in scope. Each Salmonella candidate is tagged
`typhoidal` (genome_name or species matches `\b(Typhi|Paratyphi)\b`), `non-typhoidal` (another
serovar named) or `serovar unknown`, and counts are reported by tag. The tag never drops a genome.

## R3. Quality rule, as far as the cache allows

Copy of `qc_failure` in `src/curate.rs` with the default thresholds, first failing reason wins, in
this order: `genome_quality == "Poor"`; contigs > 1000; CheckM completeness < 90; CheckM
contamination > 10; no genome length; |length / species median - 1| > 0.25 (skipped when the
species has no median). Missing quality, contig or CheckM values pass, as in the source.

Cached inputs: `genome_quality` (filled on 5,102 of 20,108 new ids), `contigs`,
`checkm_completeness`, `checkm_contamination`, `genome_length`.
Unavailable inputs: (a) the curation reference for the species median length (the medians the
curation run computed over all genome records it fetched; not in the cache); (b) the sequence
species check (`grm species-check`, `SPECIES_CHECK_FAILED`), which needs the FASTA.
Substitute for (a), fixed now: the upper median (`v[n/2]` of sorted lengths, as the source) over
the 35,151 held-FASTA genomes of the same species. Sensitivity, reported and not applied: the same
median over the cohort genomes only. (b) is not applied; every count below is therefore provisional
on the species check.

## R4. Cross-training duplicates (dropped)

A new id that shares any of these with any cohort or held-FASTA genome (a different genome id) is
dropped: BioSample, any SRA run, assembly accession (R0 key), or strain record within the same species.
Cached fields: `biosample_accession`, `sra_accession`, `assembly_accession`, `strain` on both sides.
Not cached and not checked: any NCBI-side link (BioSample to SRA to assembly) beyond these fields,
and isolates resubmitted under a new BioSample. The strain match is applied as written, with no
exception for short or generic strain names; strain-only drops are counted apart so their effect
is visible. Repeats of an older-record (pre-Jul-2025) genome outside the cohort and held set are
flagged (`repeats_older_pool_sample`) and not dropped: those genomes trained no published model.

## R5. Within-batch duplicates (dropped, label-blind)

Among new ids still in after R1 to R4, ids are joined into groups when they share a BioSample, an
SRA run or an assembly key (connected components over the three fields). In each group of two or
more, the smallest genome id (R0) is kept and the rest are dropped ("within-batch duplicate of
<kept id>"). The rule reads no label. Dropped ids are listed in `candidates.tsv`.

## R6. Label rule (copy of the source; reported calls only)

- Drug name: `antibiotic`, stripped and lower-cased, is looked up in the published table
  `corpus.json` `drug_names[].raw` (also stripped and lower-cased). kind `drug` -> its `canonical`;
  kind `class` or `excluded` -> not a drug, record dropped (as the source); kind `unresolved` or no
  entry -> unmapped, record dropped and the spelling listed with its record count. No spelling is
  mapped by hand, before or after.
- Call: `parse_call` of `resistant_phenotype`: resistant / nonsusceptible / non-susceptible -> R;
  susceptible -> S; intermediate / susceptible-dose dependent / reduced susceptibility -> I; any
  other word -> no call.
- Per genome and canonical drug: R and S both present -> conflict, no label. Else R if any R, else S
  if any S, else I if any I (intermediate only: no label). Only unconflicted R or S is a label
  (`Phenotype::reported_label`). MIC-only records are not in the cache and would never be a label
  (`analysis.rs` `cohort`: reported calls only; MIC-derived labels excluded).
- Pair to dataset: dataset id `<species slug>__<slug(canonical)>`; only the 274 published dataset ids count.

## R7. Same-study handling (flagged, reported apart, never pooled)

Study keys of a new genome: its `bioproject_accession`, every PubMed id in its `publication`
field and every `pmid` on its label records. Cohort study keys: the BioProjects and PubMed ids
(`publication`) of the cohort genomes in the held-FASTA metadata. A new genome with any key in the
cohort keys is `same-study` and goes to a separate population with its own counts. Others are
`other-study`, the primary population. A second flag, `study_in_held_set_only` (key in held-FASTA
genomes outside the cohort), is recorded and reported, and does not move the genome. Study identity
is approximate (one project can hold several papers and the reverse).

## R8. Dataset eligibility

A published dataset is eligible in a population when, after R1 to R7, it has at least 25 genomes
with label R and at least 25 with label S in that population. Populations: primary (`other-study`)
and same-study, never summed. A genome may sit in several datasets (one per drug); datasets are
not independent and genome counts are always distinct ids, never sums over datasets.

## R9. Order and funnel

Genome stages in order: R1, R2 (no species, then out of scope), R3, R4, R5, then "no label on any
published dataset", then the R7 split, then "on no eligible dataset" per population. Every stage
subtracts from the one before. What is frozen: these rules, then `candidates.tsv`, `pairs.tsv`,
`datasets.tsv` with sha256 in `SHA256SUMS`. Nothing here selects a learner, threshold or metric.
