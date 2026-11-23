# Cap arms: more than 2,000 genomes a species

Plan block C2 of `docs/assumption-plan.md` asks whether fitting more than 2,000 genomes a species helps, scored on the same held-out genomes. This page gives the arms, how held-out genomes are kept out of every wider arm's training, the scripts under `scripts/c2-cap/` that list, fit and score one species, and their tests. Nothing here changes `src/` or the release; with no new option, the release's output bytes stay the same. The release keeps its cap of 2,000 genomes a species. A split by study or lineage is a retrospective split inside one database, never external validation.

## The arms

| Arm | Cohort per species | Pairing |
|---|---|---|
| `cap-2000`, the default | the release's 2,000 genomes (`species/<slug>.json` of `rev14-274-20261112`) | the reference |
| `cap-4000` | the top 4,000 by the release's selection rule on the rebuilt revision 14 table | `common_outer_test_set_balanced` against `cap-2000` |
| `cap-all` | the top 4,096 by the same rule, or every genome the rule ranks when fewer | the same |
| `twin` | `cap-2000` again on slots of its own | scored as an arm of no change |

The rule is `curate::select_genomes` (`src/curate.rs`): a genome with an unconflicted reported R or S call on a drug of the species' datasets, not inferred from a MIC, scores the sum of 1 / (size of that call's class); genomes sort by score, highest first, then by id. The ranking does not depend on the cap, so the cohorts nest. `cohorts.py` ports the rule, cuts the ranking at each cap, and refuses a cap past 4,096 (the reviewed matrix gate, `REVIEWED_COLUMNS`), a genome the rebuilt table does not rank, and a wider cohort that drops a genome of the release's 2,000. No genome comes from a new BV-BRC query.

| Species | Ranked by the rule | Added at 4,000 | Added at 4,096 |
|---|---:|---:|---:|
| *Escherichia coli* | 8,704 | 2,000 | 2,096 |
| *Enterococcus faecium* | 2,907 | 907 | 907 |
| *Klebsiella pneumoniae* | 6,159 | 2,000 | 2,096 |
| *Mycobacterium tuberculosis* | 21,295 | 2,000 | 2,096 |
| *Neisseria gonorrhoeae* | 5,396 | 2,000 | 2,096 |
| *Salmonella enterica* | 8,170 | 2,000 | 2,096 |
| *Staphylococcus aureus* | 5,045 | 2,000 | 2,096 |
| *Streptococcus pneumoniae* | 15,530 | 2,000 | 2,096 |

The release's 2,000 lie in the top 4,000 in all 8 species, so each wider cohort is the 2,000 plus the next genomes of the ranking. The default is the release's list, not the rebuilt table's top 2,000, which differs from it by 4, 68, 2 and 1 genomes in 4 species. *Enterococcus faecium* ranks 2,907 of its 3,138 labelled genomes; the other 231 hold no call the rule counts and can enter no dataset's cohort, so its `cap-4000` and `cap-all` are one cohort, fitted once as `cap-all`; `draws.tsv` names `cap-all` in the `fit` column of its `cap-4000` rows.

## Datasets

The CAP datasets: every dataset of L in the 8 species whose release cohort is 2,000, 166 in all (`scripts/c2-cap/lists.tsv`, `species.tsv`). 23 are in M; the other 143 are the first run, 39 of them in S. Of the 143, 115 depend on none of the four label rules and 28 on the most-resistant call against the other three. Training labels come from the rebuilt revision 14 table, with each dataset's breakpoints set to its release report's (`scripts/c6-mic/release_specs.py`). `lists_check.py` refuses a run's dataset that is not first run, is another species', is listed twice, or is marked M by the label count's `marks.tsv`.

## The held-out genomes

**Lineage split.** Every arm runs through `grm arm` with K blocks and `--groups` set to the `structure.json` of one `grm structure` over the `cap-all` root, the union of every arm's genomes, so all arms split on one species grouping. `cap-2000` is the reference: keyed lineage plans, 5 outer folds. A wider arm reads the reference's plan for each of its slots and builds a cohort plan (`plan::cohort_plan`). In outer fold f:

- the arm is scored on the genomes both cohorts hold that the reference tests in f, subsampled by class to the reference fold's ratio. Here the 2,000 lie whole inside each wider cohort with the same labels, so the subsample keeps the reference's whole test fold;
- the arm trains on its cohort less the scored genomes and less every added genome whose lineage is a scored genome's lineage, and every added genome in a scored genome's k-mer duplicate group;
- an added genome is never scored.

`check_folds.py lineage` reads each fit's plans, cohort plans and `--trace-training` file on the run box and stops the species on any breach: an outer test fold that differs between arms, an added genome in an outer test fold, a scored genome in training, an added genome of a held-out lineage in training, an added copy of a scored genome in training, a default genome that trains otherwise than in the default, a cohort plan whose lineage count differs from the grouping's, or a balance drop. It writes counts per fold only.

**Study split.** `grm arm` has no study split, so each arm's study fit is `grm learn --fold-file` at seed 42 on the first-run datasets in S. The default takes the reviewed `study-whole` plan as it is. `study.py` writes a wider arm's fold file: the default's test folds, `train` and `unused` unchanged, and each added genome of the dataset's cohort in `train` when its BioProject is known and no tested genome has it, and no tested genome has its k-mer digest; else in `unused`, counted by reason. An added genome is never tested, and trains in every fold or in none. When no added genome trains, the dataset is not fitted for that arm and stays out of the study count. `check_folds.py study` holds each file to the default's test folds and every added genome in `train` to a known BioProject and digest no tested genome has.

BioProjects of added genomes come from `master-genome-meta.jsonl.gz`, which holds none of the genomes added in *Escherichia coli*, *Enterococcus faecium*, *Klebsiella pneumoniae* and *Streptococcus pneumoniae*. So a wider arm adds no study genome on 33 of the 39 datasets in S unless those records are fetched by id, which this block leaves open.

## Draws, the twin and the count

A `grm arm` fit takes evaluation seed 42; its lineage draws are its keyed slots. `cap-2000` runs 6 slots: 0 to 2 are the draws, 3 to 5 the twin. Each wider arm runs slots 0 to 2 against the reference's. `expect.py` lists every row a species must score.

`score.py` runs on the run box and writes under `score/`, never `logs/`. Before any score it checks, per row: the dataset is first run; the plan and the prediction rows exist, five learners, one row per scored genome per fold (a missing draw stops the scorer and is never read as zero); a twin row is the default's own fit at a slot no draw has; a wider arm's scored genomes in each fold are the default's; a study row's `--exp-pick` file is the dataset's `study-whole` fit at seed 42 with no inner fold, and a wider arm's test folds hold the default's tested and resistant counts, fold by fold, with at least the default's training genomes. Per dataset, arm, split, slot and learner it sums the outer calls into the measures of `scripts/c9-pick/score.py`. The count and conditions come from its own instance of `scripts/c3-inputs/score.py` with the family `cap-4000`, `cap-all` and `twin` against `cap-2000`: Holm over the 15 lineage rows of the primary measure, 3 arms by 5 learners. Study rows are one fit each, read beside the lineage count.

## The fetch and the matrices

`grm fetch fasta` fetches every `cap-all` genome. `fetch_check.py` stops the species when a release genome fails or its length, contig count or N50 differs from the release's. It counts, without stopping, release genomes whose bytes differ from the master list, as `scripts/c3-inputs/run_species.sh` does, since a FASTA whose header text changed differs in bytes, not in sequence. It records the sha256 of added genomes the master list does not hold, and takes a failed added genome out of every wider cohort. Only `cap-all` counts ksets; `cap-2000` and `cap-4000` take hard links to them through `scripts/input-comparison/run.py`'s sharing (`roots.py`). Each root runs `grm structure` and `grm matrix` with no network.

## Running

```
bash scripts/c2-cap/run_species.sh KIT SPECIES DATASETS OUT [JOBS] [DRAWS] [TWIN] [ARMS] [SPLITS] [SCORE]
```

The full run is `16 3 3 cap-4000,cap-all lineage,study yes`, one species to one 16-CPU, 32 GB machine. Each step writes a time record (wall, CPU, peak memory), a log scrubbed of genome ids and a flag, so a restart skips finished steps. The end scan (`scrub.py --scan`) counts genome ids of the release and every cohort list in `results/`, `score/` and `logs/` and fails on any. `SCORE no` moves each fit's prediction and score files to `work/unread/`, which nothing reads. The cohort lists, fold files, fits and scores stay on the run machine.

On *Enterococcus faecium* (vancomycin, ampicillin and gentamicin, one lineage slot, `cap-2000` and `cap-all`, scores unread), the fetch took 439.58 s for 2,907 genomes at 8 workers with no failure, the `cap-all` count 129.53 s, its matrix 105.81 s at 3.1 GB peak, and each fit 38 to 55 s at 2.4 to 3.6 GB peak. Every fold check passed, with no balance drop and no added copy of a scored genome. Over the 5 outer folds, `cap-all`'s added genomes trained in 3,644 of 4,505 fold places on vancomycin (901 labelled added genomes) and were kept out by lineage in 861; on ampicillin, 753 of 940 (188) and 187. None of the 907 added genomes carries a gentamicin call, so on gentamicin `cap-all`'s fit repeats the default's.

## Tests and planted faults

`tests/c2_cap` holds SYNTHETIC tests of the cohorts, roots, lists, expected rows, both fold checks, the study files, the scorer and count, the fetch check, the scrub, the time records and the runner's text; the study split is proved on SYNTHETIC fold plans only. `scripts/c2-cap/planted_faults.py` plants each fault below in a copy of the last commit, runs the one test named for it, and puts the text back, on `scripts/c9-pick/planted_faults.py`'s runner. Each test passes with no fault planted and fails with its fault.

| Fault | Caught by |
|---|---|
| 1. An added genome of a held-out lineage in training | `Lineage.test_an_added_genome_of_a_held_out_lineage_in_training_breaches` |
| 2. An added genome of a held-out BioProject in training, written by `study.py` (2a) or passed by the fold check (2b) | `Study.test_an_added_genome_of_a_held_out_bioproject_stays_out_of_training`; `Study.test_the_study_check_breaches_on_a_held_out_bioproject_in_training` |
| 3. An outer test fold that differs between arms, passed by the fold check (3a) or the scorer (3b) | `Lineage.test_an_outer_test_fold_that_differs_between_arms_breaches`; `Checks.test_an_arm_scored_on_other_genomes_than_the_defaults_outer_test_folds_fails` |
| 4. An added genome in an outer test fold | `Lineage.test_an_added_genome_in_an_outer_test_fold_breaches` |
| 5. A genome drawn that the rebuilt revision 14 table does not hold | `Gate.test_a_genome_the_rebuilt_table_does_not_hold_is_refused` |
| 6. A cohort past 4,096 | `Gate.test_a_cohort_past_4096_is_refused` |
| 7. A wider cohort that drops a genome of the 2,000 | `Gate.test_a_wider_cohort_that_drops_a_default_genome_is_refused` |
| 8. The twin scored as an arm: at a draw's slot or from a wider arm's fit | `Checks.test_a_twin_at_a_draw_slot_or_from_a_wider_fit_stops_the_scorer` |
| 9. A missing draw read as a zero | `Checks.test_a_missing_draw_stops_the_scorer_and_is_never_a_zero` |
| 10. A dataset of M counted in the first run, by the list check (10a) or the scorer (10b) | `ListsCheck.test_a_dataset_of_m_is_not_counted_in_the_first_run`; `Checks.test_a_dataset_of_m_or_outside_the_first_run_is_refused` |

## What it does not do

It reads no sealed, held-out, clinical or spent test set, and no later record: the rebuilt table is built from the records before the cut. It builds nothing past 4,096 genomes and finds, names or ranks no marker. No row per genome, genome id, k-mer, sequence, model or feature weight leaves the run machine or enters the repository; a cohort list is a row per genome and stays there. `grm arm` saves no final model, so a cap arm that passes would need a final fit on its cohort per dataset before block T1 can score it.
