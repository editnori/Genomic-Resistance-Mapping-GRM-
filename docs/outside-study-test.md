# Outside-study test

This plan is for a retrospective test on newly deposited studies from the same database. The genomes are public BV-BRC genomes that were deposited after the cohort was built and that come from studies the cohort does not hold. It is not a lockbox, because the batch's metadata and label totals were read before these rules were frozen. It is not external validation either, because the genomes come from the same database, under the same labelling practice, as the cohort.

The kit has four parts:

- this plan;
- `scripts/outside-study/run.py`, the box-side run script;
- `scripts/outside-study/score.py`, the scoring script;
- `tests/outside_study/`, the tests.

`scripts/outside-study/FROZEN.sha256` pins the SHA-256 of the plan, `score.py` and `run.py`. The `known-gene` and `predict` stages and the scoring script refuse to run unless all three files match those pins and git holds them unchanged at HEAD (`score.frozen_gate`).

## Inputs

| Input | Where | Bound by |
|---|---|---|
| W7's frozen tables | `candidates.tsv`, `pairs.tsv` and `datasets.tsv` from W7's report folder | `scripts/outside-study/w7.sha256`, checked by `score.load_inputs` |
| The served revision-12 dataset reports | the `datasets/` folder of the published report tree, one `<dataset>.json` each | `scripts/outside-study/reports.sha256`, 274 lines written by `sha256sum *.json` over that folder and checked by `score.load_report` |
| The master data root | `/home/user/grm-h1raw` on the box, read only | `run.py`, which never writes it (see `setup`) |
| The `grm` binary | built from the kit's commit with `cargo build --release --locked` | its SHA-256, recorded in every receipt |

No option moves either pin file. `--pins` and `--kit-root` do not exist; only the tests set other pin files, through `main(..., overrides=...)`.

## The twelve items

### 1. Population

The population is the rows of W7's `candidates.tsv` with `species_slug` equal to the species, `population` equal to `primary (other-study)`, and `on_eligible_dataset` equal to `yes` (`score.population`). The datasets are the rows of `datasets.tsv` for that species with `eligible_primary` equal to `yes` (`score.eligible`).

For *Staphylococcus aureus* that is 803 genomes on 7 datasets: ciprofloxacin, erythromycin, gentamicin, levofloxacin, oxacillin, sulfamethoxazole/trimethoprim and tetracycline. The 803 genomes are in 3 BioProjects: PRJNA1050907 has 598 genomes, PRJEB24782 has 172 and PRJNA819273 has 33. Across all species, W7 lists 7,258 primary genomes on an eligible dataset, in 7 species, with 82 eligible datasets.

The same-study set uses `population` equal to `same-study`, `eligible_same_study` and the `in_same_study_set` flag (`score.POPULATIONS`). It is scored in its own run, with `--population same-study`, and reported on its own. The two are never pooled. W7 lists 3,155 same-study genomes on 49 eligible datasets. *S. aureus* has none.

### 2. Labels

A label is a W7 `pairs.tsv` row whose `call` is `R` or `S`, on an eligible dataset, for a genome in the population (`score.labels`). Pairs with `call` equal to `I`, `conflict` or `none` carry no label. Such a genome is still predicted, but it is not counted on that dataset. Nothing is relabelled from MIC.

`score.labels` refuses in four cases:

- an R or S pair that is not flagged `in_frozen_test` (or `in_same_study_set` for the same-study run);
- a flagged pair that is not R or S;
- a flagged pair whose genome is outside the population;
- a pair listed twice.

Any of these would mean the two populations, or labelled and unlabelled pairs, had been mixed.

### 3. Exclusions

Exclusions are label-blind and are applied in this order (`score.EXCLUSION_ORDER`, `score.exclusions`):

1. **No sequence.** `grm fetch fasta` wrote no FASTA file. The reason is the line in `fasta_failures.tsv`.
2. **Species check failed.** `grm species-check` against the master root's recorded references in `species-check/species_check.json` gave `fail`, under the default thresholds (fail distance 0.13, nearer-species margin 0.05). A `borderline` genome stays in and is flagged in `genomes.tsv` and in the counts.
3. **Not counted.** `grm count` wrote no k-mer set. This third rule is not in the brief. It is label-blind, and without it a genome would have no prediction. The reason is the line in `counted.json`'s `failed` list.

A genome takes the first reason that applies. Every excluded genome is listed in `genomes.tsv` with its reason and the tool's message.

Without `--reset`, `grm species-check` keeps the recorded references and judges the cohort and the new genomes together. It recomputes each species' GC median and spread over all of them, so a new genome can move a cohort genome's GC flag between `pass` and `borderline`. A fail never depends on GC, and a borderline genome is kept either way, so such a move changes neither which new genome is scored nor any cohort report. The species-check stage compares every recorded reference and verdict with the new run (`run.species_check`). A recorded genome whose verdict moved between `pass` and `borderline` while its distance, its nearest species and that distance, and its `inside_species` all stayed as recorded is listed in `species-check/cohort-gc-moves.tsv`, counted in the receipt, and the run goes on. Any other change stops the run with `species-check/cohort-changes.tsv`: a changed reference, a recorded genome the new check leaves out, a verdict that moves to or from `fail`, or a flag that moves with a changed distance. This rule was fixed before the first new genome was fetched.

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
- rows for a dataset that is not eligible in the population.

Before any new genome is predicted, the `predict` stage runs `grm predict` on the cohort's own genomes for each report. Every feature value must equal the report's published bit, and every call must equal the final model's call on those bits (`score.replay_problems`). Otherwise the stage refuses.

### 5. Metrics

Metrics are computed per dataset and learner (`score.confusion`, `score.estimates`), over the scored genomes with a label on that dataset:

- TP, FN, FP and TN;
- sensitivity TP/(TP+FN), specificity TN/(TN+FP), balanced accuracy (the mean of the two) and precision TP/(TP+FP);
- VME, which is FN/(TP+FN) (FN over R), and ME, which is FP/(FP+TN) (FP over S);
- exact Clopper-Pearson 95% intervals on sensitivity and specificity (`score.clopper_pearson`), the same algorithm as `src/intervals.rs`, tested against its 30 scipy values to 1e-12.

A dataset with fewer than 25 scored R genomes, or fewer than 25 scored S genomes, after exclusions is reported as counts only: no estimate, no interval and no reading. For *S. aureus*, W7 counts these labelled genomes before exclusions:

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

- The report's figures come from models refitted in each fold on part of the cohort and scored on held-out cohort genomes.
- The new figures come from one final model, fitted once on the whole cohort, applied unchanged to genomes from other studies.

`reading.txt` says this in plain words.

### 7. Fixed reading

The reading compares the report's lineage-split estimate with the new 95% interval (`score.reading`):

- "below" when the lineage-split estimate is above the new interval's upper bound;
- "above" when it is under the lower bound;
- "consistent" otherwise.

It is read for sensitivity and for specificity. There is no multiplicity correction and no hypothesis test. `reading_counts.tsv` counts the readings per learner.

A report whose lineage split is unavailable has no lineage-split estimate. Its reading is "not read: no lineage-split estimate", and the random split is not substituted. This differs from the brief, which reads every dataset. For *S. aureus* it applies to levofloxacin, whose lineage split is unavailable ("invalid lineage outer evaluation: fold 4 validation lacks one or both classes"). So both of levofloxacin's readings are "not read" for every learner. The random-split figure is still shown beside the estimate.

### 8. By study

`groups.tsv` gives the counts per learner and dataset for three groupings:

- by BioProject (`candidates.tsv` column `bioproject`);
- by country (`isolation_country`);
- by whether the genome's MLST type occurs in the cohort.

The MLST grouping compares the `candidates.tsv` column `mlst` with the set of `samples[].mlst` values over every pinned report of the species (`score.cohort_sts`). For *S. aureus* that is 18 reports, 2,000 cohort genomes and 169 STs. A genome with an empty `mlst` is "no ST". Of the 803 *S. aureus* genomes, 793 have an ST.

`leave_one_study_out.tsv` drops each BioProject in turn and recomputes the confusion counts on the rest. It gives sensitivity, specificity and balanced accuracy only when the rest still has at least 25 R and 25 S genomes. This 25-genome threshold is not in the brief; it is the same threshold as item 5, and it stops a range being drawn from two or three genomes. `metrics.tsv` gives the range of those estimates (min and max) and how many there are. There is no study-level interval. *S. aureus* has 3 BioProjects, so each dataset has at most 3 leave-one-out estimates.

### 9. Known-gene comparator

The comparator uses the gene calls `grm count` already writes (`species/<slug>/counted.json`, `genomes[].genes`, catalog version `<version>|group-calls-v2`). It does not use AMRFinderPlus. The rule is the report's own `catalog.rule`: every rule group must have a hit whose class (`C:`) or subclass (`S:`) token matches, split on "/" (`score.token_match` and `score.explain`, ported from `catalog::token_match` and `catalog::explain`).

Before any new genome is called, the `known-gene` stage counts every cohort genome of the datasets' reports and replays the published catalog for each report (`score.known_gene_replay`). It checks:

- each sample's `groups` (as names, through `catalog_groups`);
- each sample's `catalog_resistant`;
- the whole `catalog.groups` table, including `counts_for_rule`, with its order;
- the catalog's TP, FN, FP and TN.

Any mismatch writes `known-gene/replay-problems.txt` and refuses. The comparator is then exact or absent. A report with no catalog rule gets "no rule published" and no comparator row.

The stage also refuses if the scratch root's catalog version differs from a report's `catalog.version` (2026-08-07.1 for every *S. aureus* report), or if any gene call came from another catalog.

The dry run on this host has no AMRFinderPlus catalog, and fetching one was out of scope. The stage therefore wrote "not available" there, which only a rehearsal accepts. The pass and fail paths are covered by unit tests on hand-built calls (`tests/outside_study/test_models.py`). On the box, the replay against the real cohort is the check that it is exact.

### 10. Feature carriage

For each model's features, `features.tsv` gives the share of new R genomes and of new S genomes that carry the published k-mer (the row's `found`). Beside them it gives the cohort's `n_resistant` and `n_susceptible` from the report, over its `counts.resistant` and `counts.susceptible`.

### 11. Not run, and said so

These are not run: MinHash novelty, the `--matrix` partial-feature counts, calibration, any hypothesis test, forests and boosted trees. `reading.txt` ends by saying so.

### 12. Other species

Other species run unchanged once the machinery checks pass on them: the species-check verdicts, the final-model replay on the cohort, and the known-gene replay. Use `--species <slug>`, that species' ids and that species' eligible datasets. Whether they run does not depend on the *S. aureus* scores.

## Commands

On a Boat box, at the kit's commit, the master root at `/home/user/grm-h1raw`. `$KIT` is the checkout, `$W7` the folder with W7's three tables, `$REPORTS` the served `datasets/` folder and `$WORK` a new directory outside the master root.

```sh
cd "$KIT"
cargo build --release --locked
sha256sum -c scripts/outside-study/FROZEN.sha256
# The frozen ids and the eligible datasets, from W7's tables.
python3 -I - "$W7" > "$WORK.ids" <<'EOF'
import sys
rows = open(sys.argv[1] + "/candidates.tsv").read().split("\n")
h = rows[0].split("\t")
for line in rows[1:]:
    if line:
        r = dict(zip(h, line.split("\t")))
        if r["species_slug"] == "staphylococcus_aureus" and r["population"] == "primary (other-study)" and r["on_eligible_dataset"] == "yes":
            print(r["genome_id"])
EOF
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)c[$i]=i;next} $c["species_slug"]=="staphylococcus_aureus" && $c["eligible_primary"]=="yes"{print $c["dataset_id"]}' "$W7/datasets.tsv" > "$WORK.datasets"
python3 -I scripts/outside-study/run.py all --grm target/release/grm --master /home/user/grm-h1raw \
  --work "$WORK" --species staphylococcus_aureus --ids "$WORK.ids" --reports "$REPORTS" \
  --datasets "$WORK.datasets" --jobs 8
python3 -I scripts/outside-study/score.py --work "$WORK" --w7 "$W7" --reports "$REPORTS" \
  --species staphylococcus_aureus --population primary --out "$WORK/score"
```

Each stage can also be run alone in order (`setup`, `fetch`, `species-check`, `count`, `known-gene`, `predict`). Each one checks the receipts before it.

## Outputs

`$WORK/receipts/<stage>.json` holds the stage's arguments, the `grm` path and SHA-256, the kit commit and script hashes, the SHA-256 of each input and output, the counts it prints, and its seconds. Large folders are pinned as a listing of `name<TAB>size`. A stage refuses when an earlier receipt is missing, or when any file it lists, the `grm` binary, the scripts or the arguments have changed.

| Stage | Writes |
|---|---|
| setup | `root/`, the scratch data root. It holds per-file links to the master's FASTA and sketches, and copies of its catalog, `curated/species.json` and `species_check.json`. It also holds the lists `outside__<species>` (frozen ids) and `replay__<species>` (the cohort genomes of the datasets' reports, 1,988 for the 7 *S. aureus* datasets). |
| fetch | `fetch/availability.tsv` |
| species-check | `species-check/verdicts.tsv`, `species-check/cohort-gc-moves.tsv` (and `cohort-changes.tsv` on refusal) |
| count | `count/counted.tsv`, plus k-mer sets and gene calls under `root/species/` |
| known-gene | `known-gene/datasets.tsv`, `known-gene/calls.tsv` (or `replay-problems.txt` on refusal) |
| predict | `predict/<dataset>.jsonl`, `predict/survivors.txt`, `predict-replay/<dataset>.jsonl` |

The scoring script writes these files into `--out`:

- `genomes.tsv`
- `datasets.tsv`
- `metrics.tsv`
- `reading_counts.tsv`
- `groups.tsv`
- `leave_one_study_out.tsv`
- `features.tsv`
- `reading.txt`
- `results.json`, with every number and the input hashes

The same inputs give byte-identical outputs.

## Rehearsal

`--rehearsal` on both scripts runs the kit on cohort genomes only, to show that it works. Both scripts refuse any non-cohort genome in a rehearsal, and refuse any cohort genome outside one. `scripts/outside-study/rehearsal.py` writes stand-in W7 tables for such a list. Their labels are the reports' own labels, and they include an I pair and a conflict pair. In a rehearsal, the scoring script also checks each genome's rows against the report's published bits and final-model calls. A rehearsal measures nothing, because those genomes trained the models, and `reading.txt` says so on its first line.

## Tests

```sh
python3 -I -m unittest discover -s tests/outside_study -t . -v
python3 -I tests/outside_study/mutants.py
```

The suite runs in CI (docs/ci.md). Its fixtures are hand-built: a made-up species with three reports and 126 frozen genomes whose counts, intervals, readings, groupings and leave-one-out values are worked out by hand in `tests/outside_study/fixtures.py`. The suite covers:

- every refusal of both scripts;
- the label rules and the exclusion order;
- the Clopper-Pearson values;
- the frozen gate, on temporary git repositories;
- the receipts, with a stand-in `grm` for fetch, species check and count;
- byte-identical reruns at three hash seeds.

`mutants.py` breaks the scripts in 20 ways, one at a time, in a temporary copy, and fails unless the suite catches every one.

## Limits

- The genomes come from the same database as the cohort. Their labels follow the same practices, and the metadata and label totals of the batch were read before this plan was frozen.
- Each new estimate has one model per dataset and no refit. Its interval covers only the sampling of these genomes, not study-to-study variation; the leave-one-BioProject-out range shows that variation without an interval.
- With three BioProjects for *S. aureus*, one study holds 598 of the 803 genomes.
- `setup` needs a FASTA file in the master root for every cohort genome of the datasets' reports; the replays depend on it.
