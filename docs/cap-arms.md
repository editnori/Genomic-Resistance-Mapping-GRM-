# Cap arms: more than 2,000 genomes a species

Plan block C2 of `docs/assumption-plan.md` asks whether fitting more than 2,000 genomes a species helps, scored on the same held-out genomes. This page gives the arms, how held-out genomes are kept out of every wider arm's training, the scripts under `scripts/c2-cap/` that list, fit and score one species, and their tests. Nothing here changes `src/` or the release; with no new option, the release's output bytes stay the same. The release keeps its cap of 2,000 genomes a species. A split by study or lineage is a retrospective split inside one database, never external validation. Every arm, the default included, is fitted here by `grm arm`, and on the study split by `grm learn --fold-file`, so the default's figures here are not the release report's figures and are never set beside them as if they were.

## The arms

| Arm | Cohort per species | Pairing |
|---|---|---|
| `cap-2000`, the default | the release's 2,000 genomes (`species/<slug>.json` of `rev14-274-20261112`) | the reference |
| `cap-4000`, withdrawn | the top 4,000 by the release's selection rule on the rebuilt revision 14 table | none: not fitted |
| `cap-all` | the top 4,096 by the same rule, or every genome the rule ranks when fewer | `common_outer_test_set_balanced` against `cap-2000` |
| `twin` | `cap-2000` again on slots of its own | scored as an arm of no change |

`cap-4000` was withdrawn at 09:22Z on 23 Nov 2026, before any fit: it differs from `cap-all` by 96 genomes in 7 species and is the same cohort in the eighth, so one result would be paid for twice. The runner builds no `cap-4000` root unless `ARMS` names it, and the arm's code and tests stay. The count keeps the Holm family at 15: `score.py --count --withdrawn cap-4000` holds `cap-4000`'s five lineage pairs at p = 1, and its tables name it withdrawn.

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

The release's 2,000 lie in the top 4,000 in all 8 species, so each wider cohort is the 2,000 plus the next genomes of the ranking. The default is the release's list, not the rebuilt table's top 2,000, which differs from it by 4, 68, 2 and 1 genomes in 4 species. *Enterococcus faecium* ranks 2,907 of its 3,138 labelled genomes; the other 231 hold no call the rule counts and can enter no dataset's cohort, so its `cap-4000` and `cap-all` are one cohort, fitted and counted under `cap-all` only. When `ARMS` names `cap-4000`, its `cap-4000` rows in `expect.tsv` are `same-cohort` and never scored; the run's `ARMS`, `cap-all`, lists no `cap-4000` row.

## Added genomes with a label

`cohorts.py` writes `added-labelled.tsv`: per dataset of the species and per wider arm, the species and antibiotic in full, the cohort's labelled genomes, its added genomes, and the added genomes with a label by class (resistant, susceptible). A wider arm with no added genome labelled on a dataset is `none-added` there: it is not fitted, and the dataset is out of that arm's count, its denominator included, on both splits. No threshold is set, so a dataset with one added genome labelled is fitted.

Under `cap-all`, 27 of the 143 first-run datasets are none-added, and the 143 hold 105,906 added labelled genomes:

| Species | None-added of first-run datasets | The none-added datasets |
|---|---|---|
| *Enterococcus faecium* | 8 of 11 | clindamycin, dalfopristin/quinupristin, doxycycline, erythromycin, gentamicin, kanamycin, streptomycin, teicoplanin |
| *Escherichia coli* | 8 of 29 | cefpodoxime, cefpodoxime/clavulanate, ceftiofur, colistin, erythromycin, nitrofurantoin, streptomycin, sulfamethoxazole |
| *Klebsiella pneumoniae* | 4 of 32 | cefotetan, ceftazidime/avibactam, ofloxacin, streptomycin |
| *Mycobacterium tuberculosis* | 1 of 20 | nicotinamide |
| *Neisseria gonorrhoeae* | 1 of 6 | erythromycin |
| *Salmonella enterica* | 3 of 24 | ampicillin/sulbactam, cefazolin, pefloxacin |
| *Staphylococcus aureus* | 0 of 13 | none |
| *Streptococcus pneumoniae* | 2 of 8 | azithromycin, cefuroxime |

Another 19 of the 116 fitted datasets add 1 to 19 labelled genomes each, 105 in all. No threshold drops them, so each enters `cap-all`'s count and denominator with the weight of any other dataset, though at most 19 labelled genomes join its cohort:

| Species | Antibiotic | Added labelled | Resistant | Susceptible |
|---|---|---|---|---|
| *Enterococcus faecium* | tetracycline | 6 | 6 | 0 |
| *Escherichia coli* | ampicillin/sulbactam | 1 | 1 | 0 |
| *Escherichia coli* | cefazolin | 2 | 2 | 0 |
| *Escherichia coli* | ceftriaxone | 3 | 2 | 1 |
| *Escherichia coli* | chloramphenicol | 7 | 5 | 2 |
| *Escherichia coli* | doxycycline | 2 | 2 | 0 |
| *Escherichia coli* | levofloxacin | 2 | 1 | 1 |
| *Klebsiella pneumoniae* | cefpodoxime | 6 | 1 | 5 |
| *Klebsiella pneumoniae* | cefpodoxime/ceftriaxone | 3 | 2 | 1 |
| *Klebsiella pneumoniae* | cefpodoxime/clavulanate | 5 | 0 | 5 |
| *Klebsiella pneumoniae* | chloramphenicol | 11 | 8 | 3 |
| *Klebsiella pneumoniae* | doripenem | 2 | 1 | 1 |
| *Klebsiella pneumoniae* | doxycycline | 5 | 0 | 5 |
| *Klebsiella pneumoniae* | fosfomycin | 2 | 0 | 2 |
| *Klebsiella pneumoniae* | ticarcillin/clavulanate | 5 | 0 | 5 |
| *Klebsiella pneumoniae* | trimethoprim | 15 | 10 | 5 |
| *Salmonella enterica* | sulfamethoxazole | 1 | 0 | 1 |
| *Staphylococcus aureus* | daptomycin | 14 | 0 | 14 |
| *Streptococcus pneumoniae* | meropenem | 13 | 0 | 13 |

## Datasets

The CAP datasets: every dataset of L in the 8 species whose release cohort is 2,000, 166 in all (`scripts/c2-cap/lists.tsv`, `species.tsv`). 23 are in M; the other 143 are the first run, 39 of them in S. Of the 143, 115 depend on none of the four label rules and 28 on the most-resistant call against the other three. Training labels come from the rebuilt revision 14 table, with each dataset's breakpoints set to its release report's (`scripts/c6-mic/release_specs.py`). `lists_check.py` refuses a run's dataset that is not first run, is another species', is listed twice, or is marked M by the label count's `marks.tsv`.

## The held-out genomes

**Lineage split.** Every arm runs through `grm arm` with K blocks and `--groups` set to the `structure.json` of one `grm structure` over the `cap-all` root, the union of every arm's genomes, so all arms split on one species grouping. `cap-2000` is the reference: keyed lineage plans, 5 outer folds. A wider arm reads the reference's plan for each of its slots and builds a cohort plan (`plan::cohort_plan`). In outer fold f:

- the arm is scored on the genomes both cohorts hold that the reference tests in f, subsampled by class to the reference fold's ratio. Here the 2,000 lie whole inside each wider cohort with the same labels, so the subsample keeps the reference's whole test fold;
- the arm trains on its cohort less the scored genomes and less every added genome whose lineage is a scored genome's lineage, and every added genome in a scored genome's k-mer duplicate group;
- an added genome is never scored.

`check_folds.py lineage` reads each fit's plans, cohort plans and `--trace-training` file on the run box and stops the species on any breach: an outer test fold that differs between arms, an added genome in an outer test fold, a scored genome in training, an added genome of a held-out lineage in training, an added copy of a scored genome in training, a default genome that trains otherwise than in the default, a cohort plan whose lineage count differs from the grouping's, or a balance drop. It writes counts per fold only.

**Study split.** `grm arm` has no study split, so each arm's study fit is `grm learn --fold-file` at seed 42 on the first-run datasets in S. The default takes the reviewed `study-whole` plan as it is. `study.py` writes a wider arm's fold file: the default's test folds, `train` and `unused` unchanged, and each added genome of the dataset's cohort in `train` when its BioProject is known and no tested genome has it, and no tested genome has its k-mer digest; else in `unused`, counted by reason. An added genome is never tested, and trains in every fold or in none. When no added genome trains, or none is labelled, the dataset is none-added for that arm: not fitted, and out of the study count. `check_folds.py study` holds each file to the default's test folds and every added genome in `train` to a known BioProject and digest no tested genome has.

**BioProjects.** `master-genome-meta.jsonl.gz` holds no BioProject for the genomes added in *Escherichia coli*, *Enterococcus faecium*, *Klebsiella pneumoniae* and *Streptococcus pneumoniae*, and the first and third have 33 of the 39 datasets in S. So `bioprojects.py` asks BV-BRC's genome table for each `cap-all` added genome's record, by genome id, two fields only:

```
POST https://www.bv-brc.org/api/genome/
in(genome_id,(<ids>))&select(genome_id,bioproject_accession)&limit(<n>)
```

It asks in batches of 500 and reads no other table, no other field and no lab record. A query of another table or field is refused before it is sent, and an answer with a field not asked, a genome not asked or a genome twice stops the fetch, so no new genome comes in. The answers stay on the run box beside the cohort lists. `results/bioprojects.json` holds the endpoint, the fetch's UTC time and counts only. `study.py` and `check_folds.py study` read the master rows and the answers together; a genome the two name differently counts as unknown, and an added genome whose record has no BioProject stays `unused`, counted. It runs only for a species with a dataset in S.

Fetched on 23 Nov 2026 from 09:45:37Z to 09:45:53Z, 5 batches a species, with 2,096 added genomes asked and 2,096 answered in each:

| Species | With a BioProject | Without | Same as the master rows | Other than the master rows |
|---|---:|---:|---:|---:|
| *Escherichia coli* | 985 | 1,111 | 0 | 0 |
| *Klebsiella pneumoniae* | 2,006 | 90 | 0 | 0 |
| *Mycobacterium tuberculosis* | 1,771 | 325 | 1,704 | 0 |
| *Salmonella enterica* | 1,748 | 348 | 1,689 | 0 |
| *Staphylococcus aureus* | 1,137 | 959 | 1,137 | 0 |

Before the k-mer digest check, which needs the run's ksets, `cap-all`'s added genomes could train on 27 of the 39 datasets in S. Summed over the 39, 25,853 added genomes are labelled: 21,039 share a tested genome's BioProject, 2,836 have none known, and 1,978 could train. The most are *Klebsiella pneumoniae* ceftazidime (430 of 1,988 added labelled), ciprofloxacin (418 of 1,834) and sulfamethoxazole/trimethoprim (402 of 1,735). *Klebsiella pneumoniae* cefepime (0 of 1,163) and ceftriaxone (0 of 1,315) add none, as do *Escherichia coli* nitrofurantoin, streptomycin and sulfamethoxazole and *Mycobacterium tuberculosis* nicotinamide, which have no added genome labelled.

## Draws, the twin and the count

A `grm arm` fit takes evaluation seed 42; its lineage draws are its keyed slots. `cap-2000` runs 6 slots: 0 to 2 are the draws, 3 to 5 the twin. Each wider arm runs slots 0 to 2 against the reference's. The plan names `grm learn` seeds 42 to 44 for the draws and 45 to 47 for the twin; `grm arm` takes no seed of its own, and its keyed slots are its draws, so the same roles here are slots 0 to 2 and 3 to 5 at seed 42 (`docs/assumption-plan-corrections.md`). `expect.py` lists every row a species must score, with its role: `draw`, `twin` or `study`; `none-added` (no added genome labelled, or on the study split none that trains); or `same-cohort` (*Enterococcus faecium*'s `cap-4000`). Rows of the last two are not scored.

**Refused draws.** When the default's lineage draw at a slot is refused, `grm arm` draws the next keyed repeat past the slots (6, then 7), at most twice a dataset (`src/ablate/arm.rs`, `arm_plans`). The slot's plan and each wider arm's cohort plan then hold the keyed repeat, while the trace lines and prediction rows hold the slot, so `check_folds.py` and `score.py` read both by slot. When the replacements are refused too, every lineage run of the dataset is unavailable and none is produced, as H-1's ruling has it (`docs/input-comparison-runner.md`, A refused repeat). `check_folds.py` then prints the dataset's lineage split as unavailable and breaches only on a trace line of the default or a wider arm's cohort plan at a refused slot; the run goes on. `score.py` scores each of the dataset's lineage rows 0, with why `lineage split unavailable`, on every arm and the twin, so the dataset leaves each arm's lineage count. Its study rows are scored as usual.

`score.py` runs on the run box and writes under `score/`, never `logs/`. Before any score it checks, per row: the dataset is first run; the plan and the prediction rows exist, five learners, one row per scored genome per fold (a missing draw stops the scorer and is never read as zero); a twin row is the default's own fit at a slot no draw has; a wider arm's scored genomes in each fold are the default's; a study row's `--exp-pick` file is the dataset's `study-whole` fit at seed 42 with no inner fold, and a wider arm's test folds hold the default's tested and resistant counts, fold by fold, with at least the default's training genomes. Per dataset, arm, split, slot and learner it sums the outer calls into the measures of `scripts/c9-pick/score.py`. The count and conditions come from its own instance of `scripts/c3-inputs/score.py` with the family `cap-4000`, `cap-all` and `twin` against `cap-2000`: Holm over the 15 lineage rows of the primary measure, 3 arms by 5 learners. The count is taken once, over every species' `score/datasets.tsv`, with `cap-4000` withdrawn:

```
python3 -I scripts/c2-cap/score.py --count count --withdrawn cap-4000 <each species' score/datasets.tsv>
```

A withdrawn arm may hold no row; its five lineage pairs stay in the family at p = 1, `count.tsv`'s `withdrawn` column is 1 on its rows, and `conditions.tsv` names it withdrawn and reads none of its conditions. A none-added dataset has no row in `datasets.tsv`, so it is out of the arm's denominator. Study rows are one fit each, read beside the lineage count.

## The fetch and the matrices

`grm fetch fasta` fetches every `cap-all` genome. `fetch_check.py` stops the species when a release genome fails or its length, contig count or N50 differs from the release's. It counts, without stopping, release genomes whose bytes differ from the master list, as `scripts/c3-inputs/run_species.sh` does, since a FASTA whose header text changed differs in bytes, not in sequence. It records the sha256 of added genomes the master list does not hold, and takes a failed added genome out of every wider cohort. Only `cap-all` counts ksets; `cap-2000`, and `cap-4000` when `ARMS` names it, take hard links to them through `scripts/input-comparison/run.py`'s sharing (`roots.py`). Each root runs `grm structure` and `grm matrix` with no network.

## Running

```
bash scripts/c2-cap/run_species.sh KIT SPECIES DATASETS OUT [JOBS] [DRAWS] [TWIN] [ARMS] [SPLITS] [SCORE] [AT_ONCE]
```

The full run is `16 3 3 cap-all lineage,study yes N`, one species to one 16-CPU, 32 GB machine, with N, at most 4, from the species' scaled peak fit memory with 4 GB kept free. Each step writes a time record (wall, CPU, peak memory), a log scrubbed of genome ids and a flag, so a restart skips finished steps. The end scan (`scrub.py --scan`) counts genome ids of the release and every cohort list in `results/`, `score/` and `logs/` and fails on any. `SCORE no` moves each fit's prediction and score files to `work/unread/`, which nothing reads. The cohort lists, fold files, fits and scores stay on the run machine.

On *Enterococcus faecium* (vancomycin, ampicillin and gentamicin, one lineage slot, `cap-2000` and `cap-all`, scores unread), the fetch took 439.58 s for 2,907 genomes at 8 workers with no failure, the `cap-all` count 129.53 s, its matrix 105.81 s at 3.1 GiB peak, and each fit 38 to 55 s at 2.4 to 3.5 GiB peak. Every fold check passed, with no balance drop and no added copy of a scored genome. Over the 5 outer folds, `cap-all`'s added genomes trained in 3,644 of 4,505 fold places on vancomycin (901 labelled added genomes) and were kept out by lineage in 861; on ampicillin, 753 of 940 (188) and 187. None of the 907 added genomes carries a gentamicin call, so on gentamicin `cap-all`'s fit repeats the default's; such a dataset is now none-added and not fitted.

**Fits at once.** Each dataset's lineage chain (the default's fit, each wider arm's, the fold check) is one job. The chains run largest first, by the `cap-all` cohort's labelled genomes on the dataset; the first runs alone, and after it N chains run at once: the smaller of `AT_ONCE` and (the box's memory less 4 GiB) / (the largest fit peak so far), at least 1, read again before each chain starts and written to `results/at-once.txt`. A fit killed for memory fails the run once the running chains end, and the next start runs at half the N in force then. Study fits run one at a time.

No limit in the fit path depends on wall time. Every clock read in `grm arm`'s path (`src/ablate/arm.rs`, `src/learn/planned.rs`, `src/learn/eval.rs`, `src/learn/lasso.rs`, `src/learn/inner_pick.rs`, `src/analysis.rs`, `src/main.rs`) only writes a `seconds` field; the plan's `budget` (`src/ablate/plan.rs`) counts fit units, not time; the only timeouts outside tests are the fetch's (`src/bvbrc.rs`). On the pilot's data, the six fits of *Enterococcus faecium* vancomycin, ampicillin and gentamicin (`cap-2000` and `cap-all`, slot 0) each ran alone, then in two batches of four at once, 8 fits in all, 16 threads each. Each batch fit's files equal the same fit's run alone, 56 of 56: JSON parsed, with its `seconds` values set aside (80 in the alone files), and every other file byte for byte. Four at once ran 3.58 times the pace of one at a time over the 8 fits, and 3.26 times in the slower batch.

## Tests and planted faults

`tests/c2_cap` holds SYNTHETIC tests of the cohorts, roots, lists, expected rows, both fold checks, the study files, the scorer and count, the fetch check, the scrub, the time records and the runner's text; the study split is proved on SYNTHETIC fold plans only. `scripts/c2-cap/planted_faults.py` plants each fault below in a copy of the last commit, runs the one test named for it, and puts the text back, on `scripts/c9-pick/planted_faults.py`'s runner. Each test passes with no fault planted and fails with its fault, and `test_faults.py` holds each fault's text to one occurrence in its file, since the runner stops at a text it does not find once.

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
| 11. The BioProject fetch asks another field of the genome table, `collection_year` (11a), or another table, `genome_amr` (11b) | `Fetch.test_the_fetch_asks_the_genome_table_for_genome_id_and_bioproject_accession_only` |
| 12. A wider arm with no added genome labelled fitted and counted on the lineage split | `Expect.test_an_arm_with_no_added_genome_labelled_is_none_added_on_both_splits` |
| 13. A withdrawn arm's conditions read as a fitted arm's | `Count.test_a_withdrawn_arm_keeps_its_five_pairs_at_p_one_and_the_tables_name_it_withdrawn` |
| 14. *Enterococcus faecium*'s one wider cohort scored under both wider arms | `Checks.test_a_none_added_or_same_cohort_lineage_row_is_not_scored` |

## Limits

**The study split decides nothing.** Its fits, one per dataset and arm, are kept beside the lineage count, and no condition reads them (the lead's ruling of 23 Nov 2026).

**Few added genomes may train on the study split.** Before the k-mer digest check, which can only lower these counts, 1,978 added genomes may train on 27 of the 39 datasets in S: 1,921 in *Klebsiella pneumoniae*, 30 in *Escherichia coli*, 23 in *Staphylococcus aureus* and 4 in *Salmonella enterica*. On the other 12, no added genome trains, and the arm is none-added there. Of the 25,853 added genomes labelled on the 39, 21,039 share a tested genome's BioProject and 2,836 have none known. Per dataset, with the default's tested genomes and those of them with no BioProject known:

| Species | Antibiotic | Added labelled | Share a tested genome's BioProject | BioProject unknown | May train | Tested | Tested, no BioProject |
|---|---|---:|---:|---:|---:|---:|---:|
| *Escherichia coli* | amikacin | 284 | 96 | 186 | 2 | 1,125 | 36 |
| *Escherichia coli* | ampicillin/sulbactam | 1 | 1 | 0 | 0 | 283 | 24 |
| *Escherichia coli* | aztreonam | 349 | 114 | 233 | 2 | 699 | 38 |
| *Escherichia coli* | cefazolin | 2 | 2 | 0 | 0 | 414 | 23 |
| *Escherichia coli* | cefepime | 100 | 96 | 2 | 2 | 1,124 | 77 |
| *Escherichia coli* | cefoxitin | 94 | 91 | 0 | 3 | 1,066 | 186 |
| *Escherichia coli* | ceftriaxone | 3 | 2 | 1 | 0 | 923 | 53 |
| *Escherichia coli* | ertapenem | 108 | 98 | 7 | 3 | 820 | 55 |
| *Escherichia coli* | levofloxacin | 2 | 2 | 0 | 0 | 403 | 34 |
| *Escherichia coli* | nitrofurantoin | 0 | 0 | 0 | 0 | 482 | 24 |
| *Escherichia coli* | piperacillin/tazobactam | 1,739 | 882 | 855 | 2 | 979 | 145 |
| *Escherichia coli* | streptomycin | 0 | 0 | 0 | 0 | 386 | 21 |
| *Escherichia coli* | sulfamethoxazole | 0 | 0 | 0 | 0 | 341 | 37 |
| *Escherichia coli* | tetracycline | 65 | 28 | 22 | 15 | 782 | 94 |
| *Escherichia coli* | tobramycin | 490 | 95 | 394 | 1 | 1,070 | 199 |
| *Klebsiella pneumoniae* | ampicillin/sulbactam | 1,244 | 1,171 | 0 | 73 | 744 | 0 |
| *Klebsiella pneumoniae* | aztreonam | 1,529 | 1,470 | 32 | 27 | 1,050 | 0 |
| *Klebsiella pneumoniae* | cefazolin | 1,348 | 1,276 | 0 | 72 | 1,141 | 24 |
| *Klebsiella pneumoniae* | cefepime | 1,163 | 1,163 | 0 | 0 | 1,233 | 21 |
| *Klebsiella pneumoniae* | cefoxitin | 1,371 | 1,287 | 0 | 84 | 1,181 | 0 |
| *Klebsiella pneumoniae* | ceftazidime | 1,988 | 1,526 | 32 | 430 | 1,310 | 23 |
| *Klebsiella pneumoniae* | ceftriaxone | 1,315 | 1,283 | 32 | 0 | 1,208 | 24 |
| *Klebsiella pneumoniae* | ciprofloxacin | 1,834 | 1,379 | 37 | 418 | 1,706 | 117 |
| *Klebsiella pneumoniae* | colistin | 50 | 2 | 32 | 16 | 336 | 0 |
| *Klebsiella pneumoniae* | doripenem | 2 | 2 | 0 | 0 | 243 | 0 |
| *Klebsiella pneumoniae* | ertapenem | 533 | 400 | 31 | 102 | 952 | 0 |
| *Klebsiella pneumoniae* | imipenem | 1,221 | 1,131 | 2 | 88 | 1,158 | 40 |
| *Klebsiella pneumoniae* | levofloxacin | 1,239 | 1,150 | 0 | 89 | 871 | 0 |
| *Klebsiella pneumoniae* | nitrofurantoin | 576 | 540 | 0 | 36 | 562 | 0 |
| *Klebsiella pneumoniae* | polymyxin b | 55 | 2 | 53 | 0 | 100 | 11 |
| *Klebsiella pneumoniae* | sulfamethoxazole/trimethoprim | 1,735 | 1,294 | 39 | 402 | 1,546 | 118 |
| *Klebsiella pneumoniae* | tigecycline | 212 | 148 | 0 | 64 | 924 | 0 |
| *Klebsiella pneumoniae* | tobramycin | 1,036 | 1,016 | 0 | 20 | 1,230 | 0 |
| *Mycobacterium tuberculosis* | nicotinamide | 0 | 0 | 0 | 0 | 166 | 0 |
| *Salmonella enterica* | sulfisoxazole | 1,439 | 1,240 | 195 | 4 | 544 | 78 |
| *Staphylococcus aureus* | clindamycin | 386 | 335 | 44 | 7 | 1,543 | 124 |
| *Staphylococcus aureus* | erythromycin | 1,382 | 857 | 519 | 6 | 1,957 | 384 |
| *Staphylococcus aureus* | penicillin | 811 | 764 | 44 | 3 | 1,508 | 99 |
| *Staphylococcus aureus* | sulfamethoxazole/trimethoprim | 147 | 96 | 44 | 7 | 1,389 | 29 |
| All 39 | | **25,853** | **21,039** | **2,836** | **1,978** | **35,499** | **2,138** |

**Tested genomes with no BioProject.** 2,138 of the 35,499 tested genomes, summed over the 39 datasets, have no BioProject known. An added genome of the same project, recorded only on the added genome, is not held out by the BioProject rule; only the k-mer digest check reaches it, and only when it is an exact copy.

**A copy is an equal digest on the study split.** The study split holds out an added genome whose k-mer digest equals a tested genome's. A near copy from another BioProject may train. On the lineage split, a near copy shares the scored genome's lineage and is held out with it.

**A refused lineage split is not produced.** A dataset whose lineage draws `grm arm` refuses past its two replacements leaves every arm's lineage count (Draws, the twin and the count, above).

## What it does not do

It reads no sealed, held-out, clinical or spent test set, and no later record: the rebuilt table is built from the records before the cut. Beside the FASTA fetch, it asks BV-BRC only the genome table, by the ids of genomes already ranked, for `genome_id` and `bioproject_accession`. It builds nothing past 4,096 genomes and finds, names or ranks no marker. No row per genome, genome id, k-mer, sequence, model or feature weight leaves the run machine or enters the repository; a cohort list is a row per genome and stays there. `grm arm` saves no final model, so a cap arm that passes would need a final fit on its cohort per dataset before block T1 can score it.
