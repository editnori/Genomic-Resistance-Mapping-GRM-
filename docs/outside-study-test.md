# Outside-study test

This plan is for a retrospective test on newly deposited studies from the same database. The genomes are public BV-BRC genomes whose lab records the database added on 5 to 6 Nov 2026, from studies (BioProject or PubMed id) the cohort does not hold. For 260 of the 7,258 genomes, the genome record was inserted before July 2025 (120 of the 803 *S. aureus* records, 17 of those in 2014): there the sequence record is older and the lab record is new. It is not a lockbox, because the batch's metadata and label totals were read before these rules were frozen. It is not external validation either, because the genomes come from the same database, under the same labelling practice, as the cohort.

The kit has four parts:

- this plan;
- `scripts/outside-study/run.py`, the run script, which drives `grm` stage by stage;
- `scripts/outside-study/score.py`, the scoring script;
- `tests/outside_study/`, the tests.

The frozen population is in `scripts/outside-study/frozen/`: `candidates.tsv.gz`, `pairs.tsv.gz` and `datasets.tsv`, with the rules that produced them in `RULES.md`.

## The frozen gate

`scripts/outside-study/FROZEN.sha256` pins the SHA-256 of six files: this plan, `score.py`, `run.py`, `reports.sha256`, `tables.sha256` and `frozen/RULES.md`. Every stage of `run.py` and every run of `score.py` first passes `score.frozen_gate`. It refuses unless all of these hold:

- git is the first of `/usr/bin/git`, `/bin/git` and `/usr/local/bin/git` that exists. It is never looked up on `PATH`. Every call passes `-C <kit root>` and `--no-replace-objects`, with an environment that holds only `PATH=/usr/bin:/bin`, `HOME`, `LC_ALL=C`, `GIT_CONFIG_NOSYSTEM=1`, `GIT_CONFIG_GLOBAL=/dev/null` and `GIT_OPTIONAL_LOCKS=0`. No inherited `GIT_*` variable reaches it.
- The kit root is the folder two levels above the running script. No option points it or a pin file elsewhere. `git rev-parse --show-toplevel` must be that folder.
- The running `score.py` is the kit's own file, and `run.py` checks the same for itself and for the `score.py` it loaded.
- `FROZEN.sha256` as committed at HEAD has the same bytes as the file on disk. It pins exactly the six files.
- Each frozen file, read from the commit at HEAD (`git cat-file blob HEAD:<path>`), has its pinned SHA-256, and so does the file on disk. No frozen path, or folder above it, is a symbolic link. The gate reads the committed content itself, so index flags (`--assume-unchanged`, `--skip-worktree`), a staged-only edit and another index file cannot hide a change.
- HEAD is the commit given as `--expect-commit` (40 lower-case hex characters). A real run requires that option; a rehearsal may leave it out.

The gate returns the commit, the expected commit, the SHA-256 of `FROZEN.sha256` and the pinned hashes. Every receipt and `results.json` record them.

The gate stops an accident or a casual edit. It cannot stop someone who rewrites the gate and commits that, and no check inside the kit can. The check that holds is outside the kit: the commit is public before the first genome is fetched, the launcher passes that commit as `--expect-commit`, the outputs name it, and `score.py` at that commit can be run again on the saved work directory to compare the outputs byte for byte.

## Inputs

| Input | Where | Bound by |
|---|---|---|
| The frozen tables | `scripts/outside-study/frozen/`: `candidates.tsv.gz`, `pairs.tsv.gz`, `datasets.tsv` | `scripts/outside-study/tables.sha256`, on the uncompressed bytes, checked by `score.load_tables` each time either script reads them |
| The served revision-12 dataset reports | the `datasets/` folder of the published report tree, one `<dataset>.json` each | `scripts/outside-study/reports.sha256`, 274 lines written by `sha256sum *.json` over that folder and checked by `score.load_report` |
| The master data root | read only; `run.py` never writes it (see `setup`) | the receipts, which list what `setup` read |
| The `grm` binary | built from the kit's commit with `cargo build --release --locked` | its SHA-256, recorded in every receipt and compared between stages |

Both pin files are frozen files. `--tables` names a folder of stand-in tables and is accepted only with `--rehearsal`. A real run reads the kit's own tables.

## The items

### 1. Population

The population is the rows of the frozen `candidates.tsv` with `species_slug` equal to the species, `population` equal to `primary (other-study)`, and `on_eligible_dataset` equal to `yes` (`score.population`). The datasets are the rows of `datasets.tsv` for that species with `eligible_primary` equal to `yes` (`score.eligible`). In a real run `run.py` derives both lists itself, with these functions, for `--species` and `--population`. `setup` writes them to `ids.txt` and `datasets.txt` in the work directory, and its receipt records them.

For *Staphylococcus aureus* that is 803 genomes on 7 datasets: ciprofloxacin, erythromycin, gentamicin, levofloxacin, oxacillin, sulfamethoxazole/trimethoprim and tetracycline. The 803 genomes are in 3 BioProjects: PRJNA1050907 has 598 genomes, PRJEB24782 has 172 and PRJNA819273 has 33. Across all species, the frozen tables list 7,258 primary genomes on an eligible dataset, in 7 species, with 82 eligible datasets.

The same-study set uses `population` equal to `same-study`, `eligible_same_study` and the `in_same_study_set` flag (`score.POPULATIONS`). It is scored in its own run, with `--population same-study`, and reported on its own. The two are never pooled. The frozen tables list 3,155 same-study genomes on 49 eligible datasets. *S. aureus* has none.

### 2. Labels

A label is a frozen `pairs.tsv` row whose `call` is `R` or `S`, on an eligible dataset, for a genome in the population (`score.labels`). Pairs with `call` equal to `I`, `conflict` or `none` carry no label. Such a genome is still predicted, but it is not counted on that dataset. Nothing is relabelled from MIC.

`score.labels` refuses in four cases:

- an R or S pair that is not flagged `in_frozen_test` (or `in_same_study_set` for the same-study run);
- a flagged pair that is not R or S;
- a flagged pair whose genome is outside the population;
- a pair listed twice.

Any of these would mean the two populations, or labelled and unlabelled pairs, had been mixed.

### 3. Exclusions

Exclusions are label-blind and are applied in this order (`score.EXCLUSION_ORDER`, `score.exclusions`):

1. **No sequence.** BV-BRC answered and the answer held no sequence. `grm fetch fasta` writes each failure to `fasta_failures.tsv` as the genome id and a message, and the message ends in "empty FASTA response" exactly when the database answered with no sequence (`fetch_fasta` in `src/bvbrc.rs`). Only that answer is this exclusion. Any other failure (a timeout, a connection error, an HTTP error status), or a genome with neither a FASTA file nor a recorded failure, stops the fetch stage without a receipt. The stage can then be run again in the same work directory, and each attempt is logged.
2. **Species check failed.** `grm species-check` against the master root's recorded references in `species-check/species_check.json` gave `fail`, under the default thresholds (fail distance 0.13, nearer-species margin 0.05). A `borderline` genome stays in and is flagged in `genomes.tsv` and in the counts.
3. **Not counted.** `grm count` wrote no k-mer set. This rule is label-blind, and without it a genome would have no prediction. The reason is the line in `counted.json`'s `failed` list.

A genome takes the first reason that applies. Every excluded genome is listed in `genomes.tsv` with its reason and the tool's message.

Without `--reset`, `grm species-check` keeps the recorded references and judges the cohort and the new genomes together. It recomputes each species' GC median and spread over all of them, so a new genome can move a cohort genome's GC flag between `pass` and `borderline`. A fail never depends on GC, and a borderline genome is kept either way, so such a move changes neither which new genome is scored nor any cohort report. The species-check stage compares every recorded reference and verdict with the new run (`run.species_check`). A recorded genome whose verdict moved between `pass` and `borderline` while its distance, its nearest species and that distance, and its `inside_species` all stayed as recorded is listed in `species-check/cohort-gc-moves.tsv`, counted in the receipt, and the run goes on. Any other change stops the run with `species-check/cohort-changes.tsv`: a changed reference, a recorded genome the new check leaves out, a verdict that moves to or from `fail`, or a flag that moves with a changed distance. This rule was fixed before the first new genome was fetched.

`grm count` records a digest of each genome's k-mer set (`kset_digest` in `counted.json`). `duplicates.tsv` lists every scored genome whose digest equals a cohort genome's or another scored genome's, with the other id and where it is. It is a list only, and nobody is excluded by it: the frozen duplicate rules (R4 and R5 of `frozen/RULES.md`) work on accessions. `results.json` and `reading.txt` give the two counts.

### 4. Models

Each dataset's served revision-12 report is pinned by SHA-256 (`reports.sha256`). `grm predict` reads it under contract `grm-predict-1`, with calls at the 0.5 cutoff (docs/prediction-contract.md).

The scored learners are SCM and CART always, and the lasso when its published form is a complete model: every feature it reads is listed with a readable k-mer (`score.complete_learners`). For *S. aureus*, the lasso is complete for gentamicin, levofloxacin and oxacillin. It is incomplete for ciprofloxacin, erythromycin, sulfamethoxazole/trimethoprim and tetracycline, so it is not scored there and `grm predict` says why on stderr. Forests and boosted trees are not scored. No model is refitted or tuned, and no learner is chosen by its result.

The scoring checks every row (`score.check_row`). It refuses a row if any of these differs:

- the contract;
- the evaluation revision;
- the report SHA-256;
- the dataset;
- k or the input fingerprint;
- the learner;
- the model's features, or their k-mers;
- the call or score the frozen model gives on the row's own feature values.

It also refuses:

- a genome outside the frozen list;
- a row for an excluded genome;
- two rows for one genome and learner;
- a scored genome without a row;
- rows for a dataset that is not eligible in the population;
- a `predict/survivors.txt` that differs from the survivors the scoring derives from the availability, species-check and count files.

The `predict` stage first runs `grm predict` on the cohort's own genomes for every report and judges all of those rows. Every feature value must equal the report's published bit, and every call must equal the final model's call on those bits (`score.replay_problems`). Otherwise the stage refuses, and no `predict/*.jsonl` exists. Only after every replay passes does it submit a job for the new genomes.

### 5. Metrics

Metrics are computed per dataset and learner (`score.confusion`, `score.estimates`), over the scored genomes with a label on that dataset:

- TP, FN, FP and TN;
- sensitivity TP/(TP+FN), specificity TN/(TN+FP), balanced accuracy (the mean of the two) and precision TP/(TP+FP);
- VME, which is FN/(TP+FN) (FN over R), and ME, which is FP/(FP+TN) (FP over S);
- exact Clopper-Pearson 95% intervals on sensitivity and specificity (`score.clopper_pearson`), the same algorithm as `src/intervals.rs`, tested against its 30 scipy values to 1e-12.

Sensitivity and specificity are over the genomes a learner called. A genome with no call is counted in its own column (`no_call_resistant`, `no_call_susceptible`) and is in neither denominator.

A dataset with fewer than 25 scored R genomes, or fewer than 25 scored S genomes, after exclusions is reported as counts only: no estimate, no interval and no reading. For *S. aureus*, the frozen tables count these labelled genomes before exclusions:

| Dataset | R genomes | S genomes |
|---|---|---|
| ciprofloxacin | 254 | 476 |
| erythromycin | 209 | 574 |
| gentamicin | 92 | 624 |
| levofloxacin | 60 | 346 |
| oxacillin | 659 | 56 |
| sulfamethoxazole/trimethoprim | 94 | 661 |
| tetracycline | 104 | 627 |

These are `primary_R` and `primary_S` in `datasets.tsv`.

### 6. Beside each estimate

Each row of `metrics.tsv` carries the report's own figures from both cross-validation splits (`score.reference_row`): the lineage-split and random-split sensitivity and specificity, with their n, TP, FN, FP and TN. The known-gene rule is set beside the report's `catalog` row.

The two kinds of figure measure different things:

- The report's figures for the learners come from models refitted in each fold on part of the cohort and scored on held-out cohort genomes.
- The new figures come from one final model, fitted once on the whole cohort, applied unchanged to genomes from other studies.

The known-gene rule is different again. Its rule comes from the drug's name and class alone (`catalog::drug_rule`, called in `src/analysis.rs`), with no use of the cohort's labels, so its figure in the report is a fixed rule applied to the whole cohort: nothing was fitted, and there is no cross-validation figure for it. The `reference_kind` column of `metrics.tsv` says which kind each row has ("cross-validation, lineage split" or "fixed rule, whole cohort"), and `reference_at_0_or_1` flags a row whose reference sensitivity or specificity is exactly 0 or 1. `reading.txt` says this in plain words.

### 7. Fixed reading

The reading compares the report's lineage-split estimate with the new 95% interval (`score.reading`):

- "below" when the lineage-split estimate is above the new interval's upper bound;
- "above" when it is under the lower bound;
- "consistent" otherwise.

It is read for sensitivity and for specificity. There is no multiplicity correction and no hypothesis test. `reading_counts.tsv` counts the readings per learner. The datasets share genomes and studies, so they are not independent of each other, and `reading.txt` says so. It also says that precision depends on the share of resistant genomes, which differs between the cohort and the new genomes.

A report whose lineage split is unavailable has no lineage-split estimate. Its reading is "not read: no lineage-split estimate", and the random split is not substituted. For *S. aureus* it applies to levofloxacin, whose lineage split is unavailable ("invalid lineage outer evaluation: fold 4 validation lacks one or both classes"). So both of levofloxacin's readings are "not read" for every learner. The random-split figure is still shown beside the estimate.

### 8. By study

`groups.tsv` gives the counts per learner and dataset for three groupings:

- by BioProject (`candidates.tsv` column `bioproject`);
- by country (`isolation_country`);
- by whether the genome's MLST type occurs in the cohort.

A genome with an empty BioProject is grouped under "no BioProject recorded". That group is never left out as a study and is not counted as one. An empty country is grouped under "no country recorded". The frozen primary population has neither: 0 of the 7,258 genomes has an empty `bioproject`, and 0 an empty `isolation_country`.

The MLST grouping compares the whole value of the `candidates.tsv` column `mlst` as written, scheme name and ST, with the set of `samples[].mlst` values over every pinned report of the species (`score.cohort_sts`). For *S. aureus* that is 18 reports, 2,000 cohort genomes and 169 STs. A value that is empty, or whose ST field (after the last dot) is `-` or empty, is "no ST" (`score.no_st`). For the primary population:

| Species | ST in the cohort | ST not in the cohort | no ST |
|---|---|---|---|
| E. coli | 419 | 143 | 39 |
| K. pneumoniae | 1,442 | 497 | 75 |
| M. tuberculosis | 1,093 | 35 | 24 |
| Salmonella enterica | 1,800 | 287 | 20 |
| S. flexneri | 77 | 1 | 0 |
| S. aureus | 726 | 53 | 24 |
| S. pneumoniae | 396 | 83 | 24 |

So 779 of the 803 *S. aureus* genomes have an ST. There is no table of scheme aliases. BV-BRC writes a few genomes under an older scheme name, so by ST number alone 3 more *K. pneumoniae* genomes and 15 more *Salmonella enterica* genomes would count as in the cohort.

`leave_one_study_out.tsv` drops each BioProject in turn and recomputes the confusion counts on the rest. It gives sensitivity, specificity and balanced accuracy only when the rest still has at least 25 R and 25 S genomes. This is the same threshold as item 5, and it stops a range being drawn from two or three genomes. `metrics.tsv` gives the range of those estimates (min and max) and how many there are. There is no study-level interval. *S. aureus* has 3 BioProjects, so each dataset has at most 3 leave-one-out estimates.

For each dataset, the scoring's `datasets.tsv` gives the BioProject with the most scored R genomes, with its count and share, the same for S, the number of studies, and how many leave-one-study-out folds keep at least 25 R and 25 S genomes. `reading.txt` prints those figures on one line per dataset. Where no fold is reportable, it says that no study-to-study range can be shown and that the estimate is in effect one study's. `datasets.tsv` also says, for each of those two studies, whether the fold that leaves it out is reportable. Where some fold is reportable but that one is not, `reading.txt` says that no estimate can be shown without that study and that every leave-one-study-out range for the dataset includes it. This is description only; nothing is decided by it.

### 9. Known-gene comparator

The comparator uses the gene calls `grm count` already writes (`species/<slug>/counted.json`, `genomes[].genes`, catalog version `<version>|group-calls-v2`). It does not use AMRFinderPlus. The rule is the report's own `catalog.rule`: every rule group must have a hit whose class (`C:`) or subclass (`S:`) token matches, split on "/" (`score.token_match` and `score.explain`, ported from `catalog::token_match` and `catalog::explain`).

Before any new genome is called, the `known-gene` stage counts every cohort genome of the datasets' reports and replays the published catalog for each report (`score.known_gene_replay`). It checks:

- each sample's `groups` (as names, through `catalog_groups`);
- each sample's `catalog_resistant`;
- the whole `catalog.groups` table, including `counts_for_rule`, with its order;
- the catalog's TP, FN, FP and TN.

Any mismatch writes `known-gene/replay-problems.txt` and refuses. The comparator is then exact or absent. A report with no catalog rule gets "no rule published" and no comparator row.

The stage also refuses if the scratch root's catalog version differs from a report's `catalog.version` (2026-08-07.1 for every *S. aureus* report), or if any gene call came from another catalog.

The dry run had no AMRFinderPlus catalog, so the stage wrote "not available" there, which only a rehearsal accepts. The tests run the stage with a stand-in `grm` that writes gene calls, through its pass and its refusals. On the real run, the replay against the real cohort is the check that the port is exact.

### 10. Feature carriage

For each model's features, `features.tsv` gives the share of new R genomes and of new S genomes that carry the published k-mer (the row's `found`). Beside them it gives the cohort's `n_resistant` and `n_susceptible` from the report, over its `counts.resistant` and `counts.susceptible`. These are counts of carriage; they say nothing about cause.

### 11. Beside the known-gene rule

`paired.tsv` sets each learner beside the known-gene rule, genome by genome, on each dataset where the rule is scored. Over the scored labelled genomes that have both a learner call and a known-gene call, it counts among the R genomes those both call resistant, those only the learner calls resistant, those only the known-gene rule calls resistant, and those neither does, and the same four among the S genomes. These are counts only, with no test. They show, on genomes from other studies, how many resistant genomes the k-mer model calls that the known-gene rule misses, and how many susceptible genomes it calls resistant to do so.

### 12. Not run, and said so

These are not run: MinHash novelty, the `--matrix` partial-feature counts, calibration, any hypothesis test, forests and boosted trees. `reading.txt` ends by saying so.

### 13. Other species

Other species run unchanged once the machinery checks pass on them: the species-check verdicts, the final-model replay on the cohort, and the known-gene replay. Use `--species <slug>`; the run script derives that species' ids and eligible datasets from the frozen tables. Whether they run does not depend on the *S. aureus* scores.

## Commands

At the kit's published commit `$SHA`. `$MASTER` is the master data root, `$REPORTS` the served `datasets/` folder and `$WORK` a new directory outside the master root.

```sh
cd "$KIT"
git checkout "$SHA"
cargo build --release --locked
python3 -I scripts/outside-study/run.py all --grm target/release/grm --master "$MASTER" \
  --work "$WORK" --species staphylococcus_aureus --population primary --reports "$REPORTS" \
  --expect-commit "$SHA" --jobs 8
python3 -I scripts/outside-study/score.py --work "$WORK" --reports "$REPORTS" \
  --species staphylococcus_aureus --population primary --expect-commit "$SHA" --out "$WORK/score"
```

Each stage can also be run alone, in order: `setup`, `fetch`, `species-check`, `count`, `known-gene`, `predict`.

## Receipts and stages

Each stage writes `$WORK/receipts/<stage>.json`. It holds the stage's arguments, the `grm` path and SHA-256, what the gate returned with the hashes of the two scripts, the SHA-256 of each input and output, the counts it prints, and its seconds. Files inside the work directory are named `work/<path>`, so a saved work directory can be scored anywhere. Large folders are pinned as a listing of `name<TAB>size`.

A stage refuses when its own receipt already exists: a stage never runs twice in one work directory, and there is no option to rerun one. To start again, start a new work directory. `setup` also refuses when `$WORK/root` exists and is not empty. A stage refuses when an earlier receipt is missing, or when any file it lists, the `grm` binary, the kit commit, the scripts or the arguments have changed. Every start and every end of a stage (done, refused or failed, with the message) is appended to `receipts/attempts.jsonl`, and `results.json` gives the number of attempts per stage.

`score.py` reads the six receipts. It refuses unless all six exist; each one names the commit at HEAD; the script hashes are the pinned ones and the `grm` SHA-256 is the same in all six; the species, population and rehearsal flag match its own arguments; and every work file it reads has the SHA-256 the stage that wrote it recorded.

| Stage | Writes |
|---|---|
| setup | `ids.txt` and `datasets.txt`, the lists it derived; `root/`, the scratch data root. `root/` holds per-file links to the master's FASTA and sketches, copies of its catalog, `curated/species.json` and `species_check.json`, and the lists `outside__<species>` (frozen ids) and `replay__<species>` (the cohort genomes of the datasets' reports, 1,988 for the 7 *S. aureus* datasets). |
| fetch | `fetch/availability.tsv` |
| species-check | `species-check/verdicts.tsv`, `species-check/cohort-gc-moves.tsv` (and `cohort-changes.tsv` on refusal) |
| count | `count/counted.tsv` (with each genome's `kset_digest`), `count/cohort-digests.tsv`, plus k-mer sets and gene calls under `root/species/` |
| known-gene | `known-gene/datasets.tsv`, `known-gene/calls.tsv` (or `replay-problems.txt` on refusal) |
| predict | `predict-replay/<dataset>.jsonl` (or `predict-replay/replay-problems.txt` on refusal), then `predict/survivors.txt` and `predict/<dataset>.jsonl` |

## Outputs

The scoring script writes these files into `--out`:

- `genomes.tsv`
- `datasets.tsv`, with the per-study figures of item 8
- `metrics.tsv`
- `reading_counts.tsv`
- `groups.tsv`
- `leave_one_study_out.tsv`
- `features.tsv`
- `paired.tsv`
- `duplicates.tsv`
- `reading.txt`
- `results.json`, with every number, the input hashes, what the gate returned, the attempts per stage, the `grm` SHA-256 and whether the table pins were checked

In a rehearsal every file name starts with `rehearsal-`.

Three runs of `score.py` on the same work directory, each in its own `python3 -I` process with its own random hash seed, give byte-identical outputs. A run on the same inputs with their rows shuffled gives the same outputs, except for the input hashes in `results.json`.

## Rehearsal

`--rehearsal` on both scripts runs the kit on cohort genomes only, to show that it works. A rehearsal adds `--tables DIR` and, for `run.py`, `--ids FILE --datasets FILE`; each of the three is refused without `--rehearsal`. Both scripts refuse any non-cohort genome in a rehearsal, and refuse any cohort genome outside one. `scripts/outside-study/rehearsal.py` writes stand-in tables for such a list. Their labels are the reports' own labels, and they include an I pair and a conflict pair. In a rehearsal, the scoring script also checks each genome's rows against the report's published bits and final-model calls, and it does not check the stand-in tables against `tables.sha256` (`results.json` records that). A rehearsal measures nothing, because those genomes trained the models, and `reading.txt` says so on its first line.

## Tests

```sh
python3 -I -m unittest discover -s tests/outside_study -t . -v
python3 -I tests/outside_study/mutants.py
```

The suite runs in CI (docs/ci.md) and holds 195 tests. It needs `git`, and it passes in a copy of the tree with no `.git`, because the tests make their own repositories. Its fixtures are hand-built: a made-up species with three reports and 126 frozen genomes whose counts, intervals, readings, groupings, paired counts and leave-one-out values are worked out by hand in `tests/outside_study/fixtures.py`. The learners give different calls on at least one R and one S genome, and two datasets label two genomes the other way round, so a learner or dataset swap changes a number. The suite covers:

- every refusal of both scripts;
- the label rules and the exclusion order;
- the Clopper-Pearson values;
- the frozen gate, on temporary committed copies of the kit, against each way around it the review tried (index flags, a staged-only edit, `GIT_DIR`, `GIT_WORK_TREE`, `GIT_INDEX_FILE`, a stand-in `git` on `PATH`, an edited pin file, a modified copy of `score.py` run from elsewhere, and a commit other than `--expect-commit`) and a few more (replace objects, a work tree elsewhere, a copy inside another repository, symbolic links, a missing or malformed `--expect-commit`);
- every stage of `run.py` and its receipts, with a stand-in `grm` (`tests/outside_study/fake_grm.py`) that fetches, checks species, counts with gene calls and digests, and predicts, including the refusals of each stage;
- the receipts `score.py` reads;
- byte-identical outputs from three separate processes, each with its own random hash seed, and from shuffled input rows.

`mutants.py` breaks the scripts in 73 ways, one at a time, in a temporary copy, and fails unless the suite catches every one.

## Limits

- The genomes come from the same database as the cohort. Their labels follow the same practices, and the metadata and label totals of the batch were read before this plan was frozen.
- Each new estimate has one model per dataset and no refit. Its interval covers only the sampling of these genomes, not study-to-study variation.
- For *S. aureus*, one study dominates the resistant genomes. Before exclusions, PRJNA1050907 holds 74.5% to 100% of the R genomes on every dataset: 252 of 254 on ciprofloxacin, 195 of 209 on erythromycin, 78 of 92 on gentamicin, 60 of 60 on levofloxacin, 491 of 659 on oxacillin, 94 of 94 on sulfamethoxazole/trimethoprim and 96 of 104 on tetracycline. Leaving it out leaves a reportable fold only on oxacillin (168 R and 34 S genomes). On the other six datasets no estimate can be shown without that study, and each estimate is in effect that one study's. Levofloxacin has no other study. On the other five, one or two folds stay reportable (they leave out a smaller study), so `metrics.tsv` shows a range, but every fold in it includes PRJNA1050907. On oxacillin the susceptible genomes are few: PRJNA819273 holds 33 of the 56, and leaving it out leaves 23.
- `setup` needs a FASTA file in the master root for every cohort genome of the datasets' reports; the replays depend on it.
