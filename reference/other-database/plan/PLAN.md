# Lane outside: plan for a test on data from another public database

Written 15 Nov 2026 before any susceptibility result of a release species was read. Its sha256 is recorded in `REPORT.md` and in `results/PLAN.sha256` before the first label is fetched. Any decision made after a label is read is labelled exploratory (rule 6).

## 0. What was read before this plan

- From NCBI Pathogen Detection, label-blind: AST row fields id, biosample_acc, bioproject_acc, target_acc, scientific_name, antibiotic, standard, platform, vendor, reagent, target_creation_date; isolate fields target_acc, biosample_acc, bioproject_acc, asm_acc, Run, scientific_name, asm_level, asm_stats_n_contig, asm_stats_length_bp, target_creation_date, sra_release_date. No phenotype, MIC, zone diameter, measurement sign, host or other person field was asked for or returned (the fetch script asserts it). Counts are in `results/survey-by-species.tsv` and `results/drug-name-map.tsv`.
- One interface test at 12:00:24Z returned three AST rows with phenotypes for *Enterobacter soli*, a species the release does not cover.
- From the release: `index.json`, the field names of the reports, and, while reading the field names, one printed line of a training-cohort cross-validation metric (Salmonella amoxicillin/clavulanate, SCM, random split). No outside label was read.
- A first isolate fetch lost records (the backend returns some isolate records twice and the record limit cut answers short); `scripts/ncbi_iso_fill.py` fetched the missing ones at 12:23Z to 12:28Z, after which every AST target of every species has its isolate record. The survey was rerun at 12:28:58Z.
- Two BV-BRC API tests (12:14:38Z, 12:16:03Z) returned accession fields for three public genomes. One NCBI Datasets JSON test (12:27:15Z) for two public assemblies carried a nested BioSample block; that file was deleted unread and only the TSV form, which returns the columns asked for, is used.

## 1. Source rules (S1 to S8)

A source passes when it meets all eight. Each surveyed source is scored against them in `REPORT.md`.

- S1. Public interface, no account and no key, reachable from this machine.
- S2. Its terms could be read on its own pages and allow this use.
- S3. Rule 1 of the program: the data are held at NCBI or BV-BRC.
- S4. Laboratory susceptibility results (a reported call from a laboratory test), not predictions from the genome.
- S5. Not taken from BV-BRC: the isolates used are outside BV-BRC's labelled set (step C4 below).
- S6. A genome assembly for each isolate can be fetched as FASTA from a public interface.
- S7. Susceptibility and genome fields can be fetched without any field tied to a person's record.
- S8. Enough data: for at least one release species, at least 100 candidate genomes in the label-blind survey (steps C1 to C3, C4 (a) and (b), D2, and the assembly part of D1: column `of_which_on_a_release_drug`), and at least one release drug with at least 50 of them (`drugs_with_at_least_50_candidates` at least 1).

## 2. Pilot species and size

- P1. The pilot species is the release species with the most candidate genome x release-drug pairs (`candidate_genome_drug_pairs` in `results/survey-by-species.tsv`, the run after the isolate gap fill); a tie goes to more drugs with at least 50 candidates.
- P2. If more than 4,000 candidates remain after step D4, the first 4,000 by sha256 of `"20261115:" + BioSample` (hex, ascending) are kept. The draw reads no label.
- By P1, on the survey of 12:28:58Z, the pilot species is *Escherichia coli*: 10,378 candidates on a release drug, 136,274 candidate genome x drug pairs, 34 release drugs with at least 50 candidates, 36 release datasets. Next is *Salmonella enterica* with 65,610 pairs.

## 3. Candidate funnel, in order, each step counted

- C1. NCBI AST isolates of the species: scientific name equal to the species or starting with it and a space, joined by BioSample.
- C2. At least one AST row naming a drug of a release dataset for the species. Names map by, in order: exact equality with the release's canonical name; the release's own raw-to-canonical name list (`curated/drug_names.json` of the data root, case-insensitive); the fixed table in `scripts/survey_overlap.py` (`MANUAL`). No other name maps.
- C3. A GenBank assembly (GCA_) in the isolate record. If an isolate lists several, the one with the fewest contigs, then the smallest accession.
- C4. Outside BV-BRC's labelled set: (a) the BioSample and (b) the assembly key (digits of GCA_/GCF_, version dropped) match no genome of the BV-BRC labelled set in the program's data root (`/home/user/grm-h1raw/bvbrc/`, genome_amr genome ids joined to the genome table). BV-BRC is not asked live: its pages answered HTTP 403 with a browser check at 12:14:38Z and 12:26:38Z, so its terms could not be read from here (two interface tests of its API at 12:14:38Z and 12:16:03Z were made before that was known, and no further BV-BRC request is made). Results BV-BRC took in after the data root was copied are therefore not excluded; this is a stated limit.
- D1. Duplicate by accession: the candidate shares an assembly key or any SRA run with a training genome. Training genomes are the union of `samples[].id` over the species' release reports. Their BioSample and assembly accession come from the data root's BV-BRC genome table; their SRA runs, further assembly keys and BioProjects come from NCBI (Pathogen Detection isolate records by BioSample, and the NCBI Datasets v2 genome table by assembly accession), by `scripts/training_meta.py`. Training genomes NCBI returns nothing for are counted.
- D2. Duplicate by sample: the candidate's BioSample is a training genome's BioSample.
- D4. Within the batch, candidates that share an assembly key or an SRA run are joined; the smallest BioSample accession is kept.
- P2 cap (section 2).
- F1. FASTA fetched from NCBI's public FTP site over HTTPS (`ftp.ncbi.nlm.nih.gov/genomes/all/`), the path built from the assembly accession and the assembly name given by NCBI Datasets v2; at most 3 downloads at a time, files kept in scratch. A failed download is tried once more after 30 s and then counted as failed.
- Q1. The QC of the prediction contract (`QcConfig::default`, `src/curate.rs:192`, and `qc_failure`, `src/curate.rs:410`): at most 1,000 contigs (counted from the FASTA); CheckM completeness at least 90 and contamination at most 10 (NCBI Datasets v2 values; a missing value passes, as in curate); length within 25% of the species median, `|len/median - 1| <= 0.25`, the median computed as curate does over the data root's BV-BRC genome table.
- Q2. Species check: `grm species-check` on a scratch data root holding the release's recorded check (`species-check/species_check.json`, unchanged references), the recorded references' FASTA files, the training genomes and the candidates. A candidate judged `fail` is out; `pass` and `borderline` stay, as in curate.
- D3. Duplicate by sequence: Mash distance at most 0.0005 to any training genome, from `grm`'s own k=21, s=1000 sketches and its Jaccard and distance formula (`src/cluster.rs:12`, `:35`).
- S. Scored: `grm predict` with each species report, contract grm-predict-1, on the candidates that have an AST row naming that report's drug.
- L. Labelled: section 4. Fetched only after every step above is written.

The overlap left after the funnel is reported in numbers: the nearest-training-genome distance distribution of the scored genomes (bands (0.0005, 0.001], (0.001, 0.005], above 0.005), the scored genomes whose BioProject also holds a training genome, and the overlap with the outside-study candidates (a retrospective test on newly deposited studies from the same database, kept separate from this test).

## 4. Labels

- The AST row field `phenotype` is mapped as `parse_call` (`src/curate.rs:32`): resistant, nonsusceptible, non-susceptible to R; susceptible to S; intermediate, susceptible-dose dependent, reduced susceptibility to I; anything else to no call.
- One label per genome and drug, as curate consolidates: R and S both present is a conflict and no label; otherwise R if any R, else S if any S; I only gives no label. Only R and S labels are scored.
- No label is inferred from an MIC or zone. The `mic`, `disk_diffusion` and `measurement_sign` fields are not fetched. The `standard` field and the method fields are counted per dataset; NCBI records no standard year in these fields, so years are reported as not recorded.
- Label fields fetched: id, biosample_acc, target_acc, antibiotic, phenotype, standard, platform. No host or person field.

## 5. Models

- The pilot species' reports of release `rev14-274-20261112`, each pinned by sha256 in `results/pilot-reports.sha256` before scoring.
- SCM and CART always; the lasso where `grm predict` finds it complete; where the lasso is the null fit, `grm predict` makes no call and it is reported as not scored. RF and GBT are published as importances only and are not scored; every table says so.
- Replay check before any outside genome is scored: for the first 10 training genomes in report order with a FASTA file, `grm predict` on every pilot report must find each model feature exactly as the report's `bits` give for that sample. Any mismatch stops the pilot.

## 6. Measures per dataset and learner

A dataset is estimable when at least 25 R and 25 S labelled genomes are scored; otherwise only counts are given.

- Recall (TP / (TP + FN)) with a 95% Clopper-Pearson interval, beside log loss.
- Log loss: mean of -[y ln p + (1 - y) ln(1 - p)], p clipped to [1e-4, 1 - 1e-4], with p the SCM call (1 or 0), the CART leaf value, or the lasso probability; 95% bootstrap percentile interval from 2,000 resamples of genomes, seed 1.
- Also specificity, balanced accuracy, very major error rate FN / R, major error rate FP / S, with the cutoff 0.5 of the contract.
- Reference: the same recall and log loss from the report's lineage split, from its held-out `scores` and `samples[].resistant`; the random split is shown beside it.
- A row with a null score (no call) is left out of the measures and counted, on both sides.

Readings:

- Recall is below when the pilot's upper Clopper-Pearson bound is under the lineage-split recall, above when its lower bound is over it, otherwise consistent.
- Log loss is worse when the pilot's lower bootstrap bound is over the lineage-split log loss, better when its upper bound is under it, otherwise consistent.
- A learner holds on a dataset when recall is not below and log loss is not worse.

Pass rule per learner: the learner's lineage-split figures carry over to this source for the pilot species when it holds on more than half of its estimable datasets and at least 3 are estimable. With fewer than 3 estimable datasets the learner is not assessable.

## 7. Which word the result may carry

"External validation" is allowed only when all four hold:

- W1. The source is independent of BV-BRC, whose data trained the models: for every release species with at least 100 isolates in the source, BV-BRC's labelled set holds at most 5% of them by BioSample or assembly.
- W2. No candidate is a duplicate of a training genome by accession, sample or sequence (D1 to D3 applied).
- W3. No scored genome comes from a BioProject that also holds a training genome.
- W4. At least 3 datasets are estimable.

The label-blind survey already shows W1 fails for NCBI Pathogen Detection: BV-BRC's labelled set holds 88.8% of NCBI's *M. tuberculosis* AST isolates, 80.2% of its *S. pneumoniae* and 64.8% of its *A. baumannii* AST isolates by BioSample or assembly (`results/survey-by-species.tsv`), so the two databases hold results for the same isolates and the source is not independent of the training database. The 5% bound was set after these counts were seen, but before any label; no bound at or above 64.8% would be a reading of "independent". Whatever the outcome, this pilot is "a test on data from another public database", with the overlap left stated in numbers, and never external validation. It is also not the outside-study test, which is a retrospective test on newly deposited studies from the same database.

## 8. Pre-specified descriptions (not part of the pass rule)

- The measures for scored genomes more than 0.005 from every training genome, and for those whose BioProject holds no training genome, each where 25 R and 25 S remain.
- Counts by testing standard and method per dataset.
- Anything else is exploratory.

## 9. Limits fixed in advance

- At most 5 threads while the other lanes build or fit, at most 12 GB, at most 3 concurrent downloads, one request at a time to each NCBI interface, no BV-BRC request.
- On HTTP 429, 401, or 503 `auth_unavailable`, work stops and the error is written to `REPORT.md`.
