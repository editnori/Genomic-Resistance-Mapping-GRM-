# Input arms: other inputs beside 31-mer presence

An experiment, off by default, for plan block C3 (`docs/assumption-plan.md`). The release fits its five learners on the presence of 31-mers in each genome (K31). C3 asks whether other inputs separate the classes better: unitigs, K31 plus k-mer counts, K31 plus amino-acid k-mers of called genes, the AMRFinderPlus program's calls alone, and K31 plus the program's calls. All five arms are built. The two program-call arms came first; unitigs, K31 plus counts and K31 plus amino-acid k-mers follow the lead's rulings of 00:47Z on 23 Nov 2026 (below). By the lead's ruling 17, one run with one binary fits the default, the twin and all five arms.

`grm learn --exp-inputs ARM` fits one arm in place of K31. Without the flag nothing new is called and `grm learn` is unchanged, byte for byte. `EVALUATION_REVISION` stays 14 and the report schema stays 2. An arm writes no report and removes none.

Any split here is a retrospective split inside one public database. A split by lineage or by study is not external validation. GRM's known-gene rule is GRM's own reader of the AMRFinderPlus catalog; "the program" on this page is the AMRFinderPlus program itself, a separate thing.

## The arms

| Arm | Columns | State |
|---|---|---|
| K31, the default | the release's, `grm learn` with no `--exp-inputs` | runs today |
| `program` | one per element symbol of the program's Type AMR rows | built |
| `k31-program` | the program's columns, then the dataset's K31 rows less those whose cohort pattern a program column has | built |
| `k31-counts` | [count ≥ 2] of each 31-mer in a genome, then the K31 rows less those whose cohort pattern a count row has | built |
| `unitigs` | one per distinct unitig pattern of the graph over every 31-mer of the species' genomes, in place of K31 | built |
| `k31-aa` | amino-acid 10-mers of the genes Prodigal calls, then the K31 rows less those whose cohort pattern an amino-acid row has | built |

The program is AMRFinderPlus 4.2.7 with database 2026-08-07.1 and `--plus`, run from the AMRFinderPlus kit as `scripts/refit/run_species.sh` runs it. Each table must equal the sha256 the release recorded for it before any fit reads one. A genome has a column when its table has a row of Type AMR with that element symbol, either scope and any method. The program arms run on A, 237 datasets in 21 species: database 2026-08-07.1 has no organism option for *Mycobacterium tuberculosis* or *Staphylococcus haemolyticus*.

The three matrix arms follow the lead's rulings of 00:47Z on 23 Nov 2026 on questions 1 to 6. The unitig graph is built over every 31-mer of the species' genomes, those held by one genome included, as the plan's words have it ("a compacted graph of the cohort's 31-mers"), and not over the matrix's 31-mers alone. GRM's own compaction may build it only if, on the pilot species, its unitigs equal bcalm2's or Bifrost's on the same genomes as a set, with that program's version recorded and counts printed only; otherwise that program builds it. A genome holds a unitig when its kset holds every 31-mer of it, and the inputs file also counts the genome and unitig pairs held in part. There is one graph per species, built with no label. Genes are called by Prodigal 2.6.3 in single mode at table 11, with amino-acid k of 10. The counts come from GRM's own counter with one bin, a count of 2 or more, checked on SYNTHETIC sequence against a recount written apart.

## The three matrix arms

`grm exp-input-build counts|aa|unitigs --species SP --out DIR` builds one matrix per species, once, over every genome of the species matrix and in its column order. It reads no label, phenotype, fold or seed. Each genome's values go into a kset and the release's `matrix::build` merges them at support 2, so identical patterns collapse as they do for the species matrix. `DIR/matrix` is in the release's matrix format; `DIR/genomes.tsv` has a row per genome and stays on the box; `DIR/build.json`, format `grm-c3-build-1`, holds counts only.

- `counts`: the canonical 31-mers a genome's assembly holds 2 or more times, by the release's own counter (`kset::count_file`) with its floor at 2. A SYNTHETIC test holds it to a recount from strings written apart.
- `aa`: Prodigal 2.6.3 calls genes on each assembly in single mode at table 11 (`-p single -g 11`). GRM reads each protein, finds its contig by the sequence number in its `ID=<sequence>_<gene>` (Prodigal numbers an assembly's records from 1 in file order, so two records with the same first word, or one with an empty header, are still told apart), reads that contig as Prodigal reads a sequence line (the bytes `A` to `z` only, so a space, a `-` or a carriage return is dropped), holds the protein to its own table 11 translation at Prodigal's coordinates, and takes every amino-acid 10-mer. By the lead's ruling 20, a protein unequal to that translation stops the build. A stop or an ambiguous residue ends a window. GRM runs Prodigal as a separate program and links none of its code. The run script installs it from the kit, records its version and sha256, and refuses at the kit step when it is absent or differs.
- `unitigs`: the compacted de Bruijn graph over every 31-mer of the species' genomes, those held by one genome included, built by GRM (`--graph grm`) or read from an outside build's unitigs FASTA (`--graph outside --outside-unitigs FILE`). A genome holds a unitig when its kset holds every 31-mer of it; `build.json` also counts the genome and unitig pairs held in part. Given both, `--graph grm --outside-unitigs FILE` builds GRM's graph and writes in `equality` the 31-mers and unitigs in each build, in both and in one only, a unitig being its set of 31-mers.

GRM builds the graph. On the pilot species, *Klebsiella aerogenes* (93 genomes), GRM's graph and that of bcalm2 2.2.3 (`-kmer-size 31 -abundance-min 1`, its release binaries pinned by sha256 in the run script) hold the same 31,083,227 31-mers and 1,370,388 unitigs each: 1,370,388 in both and 0 in one only. On *Escherichia coli* (2,000 genomes, 10.232 Gbp), the second species, they hold the same 104,518,454 31-mers and 5,493,343 unitigs each: 5,493,343 in both and 0 in one only. So, by the lead's ruling on question 1, GRM's compaction builds it, and `GRAPH` defaults to `grm`. On *Escherichia coli* GRM's build, with the comparison, took 232.75 s at a peak of 4,268,980 KB. bcalm2 alone took 689.13 s at 9,445,936 KB, and GRM still had to read its unitigs after that. `GRAPH` check runs both on any species: GRM's graph is used, the equality counts go into the build's record, and a unitig or 31-mer in one build only stops the species, which then runs again in a fresh `OUT` with `GRAPH` bcalm2.

`grm learn --exp-inputs ARM --exp-input-matrix DIR` projects the build onto each cohort as `Dataset::project` projects K31: at the dataset's floor, with every label set to S, identical cohort patterns collapsed to the first. `k31-counts` and `k31-aa` put the added rows first and then the K31 rows, less a K31 row whose cohort pattern an added row has, which the inputs file counts (questions 9 and 10). `unitigs` uses the added rows in K31's place.

## Inside the folds

An arm hands `evaluate` a new matrix and nothing else. `analysis::project_cohort`, which both the drawn splits and `--fold-file` use, passes the dataset's K31 projection and its genome ids to `input_arms::apply`. The release's per-fit code then runs over the arm's matrix as it runs over K31: the support floor of max(3, n/200) and the collapse of identical patterns on each fit's training genomes, and the chi-square screen to 4,000 for forest, boosted trees and lasso on each fit's own fitting genomes.

What an arm builds once per dataset reads no label. Its one cohort-wide filter is the one `Dataset::project` applies to K31: a program column or an added row is kept when its support over the cohort lies in [floor, n − floor], and identical cohort patterns collapse to the first (in byte order of the symbol for a program column). Every training subset implies that range, so the filter can neither add nor drop a training candidate, whatever the outer test genomes hold. A cohort genome with no table stops the run on that dataset; no row is ever filled.

The folds come from the labels, the lineages and the seed, never from a column, so every arm at seed s is scored on the default's folds at seed s. The scorer checks that.

## The files

`grm learn --exp-inputs ARM --exp-program-tables DIR --exp-pick OUT --exp-pick-inner 0 --force` reads `DIR/<genome id>.tsv` and writes, per dataset, split and seed, C9's result file at no inner fold (`docs/served-pick.md`, The result file): each outer fold's genomes by class and each learner's calls, with no genome id. Beside it, per dataset, it writes `OUT/<dataset>.<tag>.inputs.json`, format `grm-c3-inputs-1`, where the tag is `drawn` for the drawn splits and the fold file's scheme for a study run. That file holds counts only: the tables read, the program names seen, out of range, collapsed and kept, the K31 rows, those left out as a program column's pattern, and the arm's rows; for a matrix arm, its build's kind and species patterns, the added rows, the K31 rows left out as an added row's pattern, and for unitigs the cohort's genome and unitig pairs held whole and in part. It holds no column name, genome id, k-mer, unitig, protein k-mer, sequence or position.

Every joined arm, `k31-program` among them, also writes each split's screen count, `OUT/<dataset>.<split>-s<seed>.screen.json`, format `grm-c3-screen-1`. Per outer fold it gives the fold's training genomes, the candidates from each block after the fold's floor and collapse, and for the forest, boosted trees and lasso in that order the genomes it fits on (boosted trees and lasso hold some out to tune), the number the chi-square screen keeps (4,000 or every candidate) and how many of each block pass. One candidate pass and one screen run over both blocks (question 9); the count only reads what that screen kept.

`--exp-program-column` takes `symbol` (the default), `subclass` or `class`, and `--exp-program-rows` takes `amr` (the default) or `all`, so the lead's other answers to questions 7 and 8 need no new code. The flag needs `--force`, so a kept report can never skip an arm, and it refuses beside `--learner-seed`, the `--exp-tie-*` options, `--exp-null-refit`, `--exp-mic` and `--exp-cohort`.

## The scorer

`scripts/c3-inputs/score.py` checks every result `expect.py` lists and stops before scoring on any failure:

- C9's checks at no inner fold: format, revision, dataset, scheme, seed, the five learners in order, each fold's calls equal to its test genomes, and the calls summed to the split's balanced accuracy;
- each arm's folds equal to the default's at the same dataset, split and seed, by class;
- the default's seed-42 lineage figures and served learner equal to the release report's;
- the twin at seeds of its own, and neither the default nor the twin with an inputs file;
- each inputs file of its arm, dataset and tag, with a table for every cohort genome and its rows adding up;
- each matrix arm's inputs file naming its build and holding no program field, its added rows at most the build's species patterns, its rows its blocks' each once, and only the unitig arm counting the pairs held whole and in part;
- each joined arm's screen count, of its run's format, revision, dataset, scheme, seed, arm and first block: one row per outer fold with the fold's training genomes, each learner's fitting genomes inside them, forest, boosted trees and lasso in that order, and each block's passed count inside its candidates;
- an arm's study result refused for no candidate as a refusal record of its own run and that reason, beside a scored default.

It scores each learner on its outer calls summed over the folds: balanced accuracy, recall, specificity, the very major and the major error rate, each with its count, denominator and Wilson 95% interval, and log loss. A dataset's value is the mean of its draws.

Only `grm learn --fold-file` refuses a run because an outer fold's training genomes give no candidate feature, and it refuses that dataset's whole study run; in a drawn split such a fold gets the release's fits on no feature and is scored as any other. By the lead's ruling 16, an arm refused this way where the default is not is scored as the constant model on the default's folds: every test genome of a fold gets the majority class of the fold's training genomes, R on a tie, and a score equal to their R share, clipped as log loss clips. Its balanced accuracy is 0.5 when every fold's training majority is the same class. `draws.tsv` marks those rows in `constant`, and `datasets.tsv` counts them in `constant_draws`.

`score.py --count` takes several species' tables once every arm of the family has results or is withdrawn, and refuses before that. By the lead's addition to question 12, the family stays 30 pairs whatever happens: an arm withdrawn before it runs is named with `--withdrawn ARM`, keeps its five pairs at p = 1 (`withdrawn` in `count.tsv`), and the entry "C3 design" names the arm and the time. The count refuses a withdrawn arm that has results, and the twin or the default as withdrawn. Per arm and learner it counts the datasets better, worse and the same by 0.02 against K31 on the arm's own datasets. Every difference is the arm's value less the default's, at 12 decimals; lower is better for log loss and both error rates. Beside the count stand the median difference with its interval from 4,000 species resamples and the one-sided sign test. Holm runs over 30 lineage pairs: 6 arm rows (unitigs, counts, amino-acid, program, K31 plus program and the twin) by 5 learners. An input arm touches all five learners, so `conditions.tsv` reads the plan's conditions by the lead's ruling 15: a learner carries the change when it meets 1 and 2 both, and at least 3 of 5 must carry it; condition 3 follows the plan's words (the twin meets neither 1 nor 2), so it holds when the twin meets 1 in fewer than 3 of 5 learners and 2 in fewer than 3 of 5, each read apart; 4 holds for all five learners; 5 holds for each learner that carries. Beside it, and deciding nothing, the file writes the stricter reading: 1 for at least 3 of 5, 2 and 5 for all five, 3 the twin meeting neither, 4 for all five. Datasets with a constant-model result enter the count, and `count.tsv` gives their number per row in `constant`; `count-no-constant.tsv` is the count with them left out, with Holm over its own 30 rows, and decides nothing.

## Datasets

`scripts/c3-inputs/lists.py` writes `scripts/c3-inputs/lists.tsv` from `scripts/c9-pick/l.tsv` and the AMRFinderPlus comparison's record, with each species and antibiotic in full.

| Set | Datasets | Species | Less M, the first run | First run in S |
|---|---:|---:|---:|---:|
| L | 268 | 23 | 233 | 55 |
| A | 237 | 21 | 204 | 54 |

On the datasets not in M, three of the four label rules give the same labels. Of L's 233, 32 hold an S and R pair and depend on the most resistant call against the other three; of A's 204, 21 do. The rest depend on none of the four. On M, the entry Label rule of `docs/assumption-plan-corrections.md` now sets `r-or-s-wins`, the release's rule, and each claim on M takes two readings; the run on M has not happened.

## Running

`scripts/c3-inputs/run_species.sh KIT SPECIES DATASETS OUT [JOBS] [SEEDS] [TWIN] [ARMS] [TABLES] [GRAPH]` runs one species on a box of its own, with no prompt. It checks its kit by sha256 (the release, the data root's list, the curated table, the study plans, the label count's table, the AMRFinderPlus kit and the recorded tables' sha256 list). It refuses a dataset in M, and refuses at the kit step a box without numpy, `cargo`, `objcopy`, `unshare`, `flock`, `/usr/bin/time`, `curl`, `tar`, `xz` or `zstd`. It builds `grm` and the species' data root as `scripts/c9-pick/run_species.sh` does, installs the program and proves it on the kit's fixture, runs it on every genome of the species with a recorded table, and holds each table to its recorded sha256, again at every start before any fit, since a restart skips the program's step. Given `TABLES`, a folder of the program's tables made before as `TABLES/<species>/<genome id>.tsv`, it installs and runs no program: it copies the cohort's tables in and holds each to its recorded sha256, and a table missing or unequal stops the run with the counts.

By the lead's ruling 17, one run with one binary covers the default, the twin and all five arms, since the default is refitted in the same run with the same binary. `ARMS` defaults to all five. Before the fits it builds the three matrices with no label: `inputs_counts`; `prodigal`, then `inputs_aa`, which removes the proteins; and, for `GRAPH` bcalm2 or check, `bcalm`, then `inputs_unitigs`, which removes bcalm2's files. The kit step refuses a box when Prodigal 2.6.3 (`prodigal.linux`) is absent from the kit or differs from its pinned sha256, and the same for bcalm2 when `GRAPH` needs it. The tools step records each version and sha256 in `binary.txt`. Each build's files get a sha256 list in `logs/`, held again at every start, and its counts-only `build.json` is copied to `results/builds/`. A step or fit the kernel kills for memory writes `OUT/FAILED` naming it, and `logs/run.log` says so (signal 9, which GNU time records for a step as the shell's exit status 137); a restart with a smaller `JOBS`, which holds fewer genomes and buckets at once (32 GB is the largest box), redoes that step and skips the rest.

It fits the default at each seed and checks its seed-42 reports against the release's with `scripts/refit/refit_check.py` before any arm runs. It then fits each arm at each seed, the twin, and the study plans, scores into `score/` (the scorer's output too, as `score/score.log`, so nothing in `logs/` or `logs/run.log` holds a score), writes `logs/timing.tsv` of each step's wall time and peak memory, and writes `OUT/DONE` with the sha256 of its list of output files. Each log is copied with every genome id written as `[genome id]`, and the end fails if a file under `results/`, `check/`, `score/` or `logs/` names one. Run again, it skips finished steps and finished results. When `grm learn` refuses an arm's study fit for no candidate and the default's study result exists, `scripts/c3-inputs/refused.py` writes `results/picks/<arm>/<dataset>.study-whole-s42.refused.json` in place of the result, with the reason and the fold's counts, and the run goes on. Any other failure, that refusal of the default among them, writes `OUT/FAILED` and exits 1.

The first pilot, one draw of *Klebsiella aerogenes* with ertapenem, imipenem and meropenem (93 genomes) on a 16-CPU, 32 GB box, took 8 min 44 s end to end with the two program arms: the build 42.26 s, the program's environment 38.73 s, the program 372.69 s, and the fits 17.98 s for K31, 2.74 s for the program alone and 11.79 s for K31 plus the program. The peak was 2,087,888 KB, in the count, and the out folder ended at 6.77 GB.

The pilot of all five arms ran the same draw at seed 42, with no twin, `GRAPH` check and `TABLES` (the program's 93 tables of an earlier pilot, each equal to its recorded sha256), alone on a 16-CPU, 32 GB box from 04:13:58Z to 04:19:16Z on 23 Nov 2026, 5 min 18 s. The builds took 4.52 s for the counts; 34.86 s for Prodigal and 1.57 s for the amino-acid build; and 55.28 s for bcalm2 and 14.81 s for the unitig build with the comparison. The fits took 32.51 s for K31, 3.82 s for the program alone, 28.72 s for K31 plus the program, 14.00 s for K31 plus counts, 11.37 s for K31 plus amino-acid k-mers and 16.75 s for unitigs, each under 1,000,000 KB at peak. The peak was 2,404,604 KB, in bcalm2; without it, 2,123,940 KB, in the count. The out folder ended at 5.17 GB. The two graphs were equal again, 1,370,388 unitigs in each and 0 in one only. K31's seed-42 reports equal the release's in all 3 datasets, 513 of 513 checks pass, and the files copied back hold 0 genome ids. This run used 1.10 to 2.05 times the CPU seconds of an earlier run of the same steps on the same box, at an earlier commit (each step over 5 s), so a box's speed varies by that much; the cost in the C3 arms entry of `docs/assumption-plan-corrections.md` takes this slower run.

## Tests and planted faults

`input_arms` holds 14 SYNTHETIC tests, `input_build` 8, `learn::eval` one on the screen count and `src/main.rs` 2 on the flags. `tests/c3_inputs` holds 347 tests: 152 SYNTHETIC ones, 71 of the scorer, the refusal record, the lists, the expected results and the run script's rules and 81 of the result page's table script, the run's checks and their fault runner (`docs/input-arms-result.md`, Method); and 195 of the run on the datasets of M (below). `scripts/c3-inputs/planted_faults.py` plants each fault below in a copy of the tree and runs the test meant to catch it, which must fail.

| Fault | Caught by |
|---|---|
| 1. A cohort filter that reads the outer test genomes' calls | `no_outer_test_genome_changes_a_fits_candidates`: emptying or filling the outer test genomes' calls changes no fit's candidates |
| 2. A joined arm that drops K31 (2a), holds the program block twice (2b) or keeps a K31 copy of a program column (2d), or a scorer that passes one (2c) | `the_joined_arm_holds_each_block_once_the_program_first`; `a_k31_row_equal_to_a_program_column_enters_once_as_the_programs`; `Checks.test_a_joined_arm_that_drops_or_doubles_a_block_fails` |
| 3. A program column built from the label | `an_arm_matrix_does_not_change_with_the_labels` |
| 4. An arm scored on other folds (4a) or paired with the default at another seed (4b) | `Checks.test_an_arm_on_other_folds_fails`; `Checks.test_an_arm_is_held_to_the_default_at_its_own_seed` |
| 5. A twin at a draw's seed, in the scorer (5a) or the list (5b), or fitted with an arm (5c) | `Checks.test_a_twin_at_a_draw_seed_fails`; `Expect.test_a_twin_at_a_draw_seed_is_refused`; `Checks.test_a_twin_with_an_inputs_file_fails` |
| 6. A default not held to the release | `Checks.test_a_default_that_differs_from_the_release_fails` |
| 7. A difference taken as the default's value less the arm's | `Count.test_every_difference_is_the_arms_value_less_the_defaults` |
| 8. Holm over one arm's five learners | `Count.test_holm_runs_over_the_thirty_lineage_pairs_of_the_primary_measure` |
| 9. A count taken on part of the family | `Count.test_the_count_refuses_part_of_the_family` |
| 10. An arm run without `--force` | `learn_takes_the_c3_input_arm_only_when_asked` |
| 11. A matrix arm's added block built from the label (11a), its floor taken on the cohort's counts with the outer test genomes in (11b), or boosting's screen count taken on a holdout drawn from every sample (11c) | `a_matrix_arm_does_not_change_with_the_labels`; `no_outer_test_genome_changes_a_matrix_arms_candidates`; `a_screen_count_reads_each_learners_own_fitting_genomes_only` |
| 12. A joined matrix arm that drops K31 (12a), holds its added block twice (12b) or keeps a K31 copy of an added row (12c), or a scorer that passes one (12d) | `a_matrix_arm_holds_each_block_once_the_added_first`; `MatrixArms.test_a_matrix_arm_that_drops_or_doubles_a_block_fails` |
| 13. A matrix arm scored on other folds | `MatrixArms.test_a_matrix_arm_on_other_folds_fails` |
| 14. A count column that is presence, every 31-mer seen once kept | `the_count_arm_keeps_the_31mers_seen_twice_as_a_recount_from_strings_does` |
| 15. A unitig held on any one of its 31-mers | `a_unitig_is_held_only_when_every_kmer_of_it_is` |
| 16. A protein read in another frame (16a) or with table 4, TGA as W (16b), or an amino-acid 10-mer across a stop or an ambiguous residue (16c) | `a_protein_must_be_the_table_11_translation_at_its_coordinates`; `no_amino_acid_10mer_spans_a_stop_or_an_ambiguous_residue` |
| 17. A graph that drops the 31-mers of one genome only | `the_graph_holds_single_genome_kmers_and_breaks_at_their_branch` |
| 18. A screen count taken on all samples, the outer test genomes too | `a_screen_count_reads_each_learners_own_fitting_genomes_only` |

## The run on the datasets of M

The first run left out the 35 datasets of M, in 9 species. A dataset is in M when, for its antibiotic, some genome's records pair an intermediate call with S or with R. The entry M run of `docs/assumption-plan-corrections.md` says how they run once, read two ways. This part says what block C3's code for that run does.

**Not run.** The run has not happened: no dataset of M has been fitted, and no figure below comes from one. These are retrospective splits inside one database, never external validation, and "no difference shown" is never "equal".

For block C3, `run_m.sh` and `m_reading.py` admit the datasets of M; `run_species.sh` and `page_tables.py` stay as the first run left them and still refuse one, because blocks C4 and C7 read them and Same results needs the first run's scorer unchanged.

`scripts/c3-inputs/run_m.sh KIT SPECIES OUT [JOBS] [TABLES]` runs a species' datasets of M at the first run's seeds (42, 43 and 44), twin (45, 46 and 47), arms (all five, the program arms on the 33 in A) and graph (GRM's own), and the study split on the 14 in S. Its data steps, builds, program and fit functions are `run_species.sh`'s own, read from it; a test holds the 27 functions it reads and 10 of its 11 pinned files to that script.

- **Binary.** It builds two binaries from `git archive` with rustc 1.98.1, `--locked` and cargo's home at `/home/user/.cargo`, the first run's. The first is the first run's commit, `8bfd4ca`. The second is the head of `run-m/input-arms-grm`, and every fit uses it.
- **The change.** The second is `8bfd4ca` with one change, `grm learn --exp-scores DIR`, off by default: each lineage fit writes its held-out genomes with their fold and the five learners' scores to `DIR/<dataset>.lineage-s<seed>.scores.json` (format `grm-c3-scores-1`), and no fit changes.
- **Diff.** On the box, before any build, the run script refuses the second commit if its diff from `8bfd4ca` touches any file but the six it names. A test makes the same check where the clone holds both commits.
- **Same build.** The first binary's code sha256, the binary less its build id, must be the species' row of `reference/input-arms/page-binaries.tsv` (`m_checks.py build`). Cargo's home is built into the binary, so each build's record names it, and so does a difference.
- **Check species.** The species with the fewest first-run datasets among those that hold a dataset of M and a first-run dataset, the first by name on a tie: *Neisseria gonorrhoeae*, with 6.
- **Same results.** The check species' box fits its 6 first-run datasets again with the second binary. `score.py`'s `datasets.tsv` must give back its rows of the committed `reference/input-arms/page-datasets.tsv` at 12 decimals, that file equal to its sha256 in `SHA256SUMS` (`m_checks.py results`).
- **Checked.** Its `check/m-checked.txt` names the second binary's code sha256, and every other species' box refuses to fit until its own binary matches that record.
- **Same counts.** For every fit, each fold's calls and log loss sum must come back from the genomes' scores by GRM's cutoff (R at 0.5 and above): a lineage fit's from its score file, the default's and the twin's also from their reports, which the run keeps as each fit writes them, and a study fit's from its fold result (`m_checks.py counts`). It runs on the check species' first-run fits before any dataset of M, and on every fit of M before scoring.

On a difference in any of the three the run names what differs, writes `OUT/FAILED` and stops. The lead, who rules on the plan's entries in `docs/assumption-plan-corrections.md`, then rules before another build, machine or genome set is used.

**List.** After the matrix and before any fit, `m_list.py` builds the list. The pinned table `curated.tar` holds each pair's consolidated call and no record's call. So `m_list.py` reads the calls from the kit's `genome_amr.jsonl`, the 1,285,111 lab records that table was rebuilt from on 21 Nov 2026, which the run pins by its own sha256. `m_list.py` prints and flushes that sha256 before any check, so a stop still shows it.

**Records.** Every genome and antibiotic pair of the records must consolidate to the table's call under r-or-s-wins (S with an intermediate call is S, an intermediate call with R is R, R with S is a conflict), and the records' name counts must equal `curated/drug_names.json`'s, or it stops.

**Listed genomes.** A dataset's list is the genomes of its cohort (the matrix's genomes with an unconflicted R or S call not taken from an MIC, as `lists.tsv` counts them) whose records for its antibiotic hold an intermediate call. The list stays in `work/m/list/` on the box.

**Floor.** `results/m/floor.tsv` gives each dataset's cohort, its listed genomes and the genomes left for the second reading, each by class, and whether the dataset passes the floor (the column `enters`): 25 genomes of each class at least, the release's class minimum. Small lists, below, says what else the second reading needs.

**Floor form.** `m_reading.py floor` writes `m-floor.tsv` and `m-left-out.tsv` from the species' `floor.tsv` files alone, before any score is read. `m-floor.tsv` gives each dataset's second reading as in it, left out under the floor, or left out for a short list (Small lists, below). `m_reading.py count` stops unless both files are, byte for byte, the ones the floor form writes from the same species' files.

**Left-out table.** `m-left-out.tsv` names each dataset left out of the second reading, for either reason, in full, with the reason and its counts by class. It goes on the page as `reference/input-arms-m/m-left-out.tsv` before any score is read.

**Second reading.** `second.py` reads every result again with no refit. Each dataset's listed genomes leave every outer test fold of every arm, the default and the twin, at every seed and on both splits. A fold that holds a listed genome gets each learner's calls taken again from the scores; a fold that holds none keeps its counts; the training counts stay as fitted. An arm's study fit refused for no candidate is the constant model on the default's study folds less the list (ruling 16), and a split that could not be drawn stays unscored. Before it reads, every result's counts must come back from its scores, every listed genome must be in each lineage fit's cohort, and the genomes each lineage fit keeps must give `floor.tsv`'s counts. A dataset under the floor is not scored in this reading. It writes `score/second/` in `score.py`'s columns.

**Counts.** `m_reading.py count` takes three counts in each reading, each its own call of `score.py`'s count, so condition 2's Holm family of 30 is the block's own for each count and each reading:

- the first run as read: the rows of the committed record, `reference/input-arms/page-datasets.tsv`, on balanced accuracy and recall, the only measures it holds; the same in both readings;
- M alone: the datasets of M on the six measures; it decides nothing (the entry Pick rule);
- all datasets: the record's rows and those of M, on balanced accuracy and recall.

Before it counts, `m_reading.py check` must give back from the record every row of the first run's published count on those two measures, 110 rows, and its 55 condition rows (`page-count.tsv` and `page-conditions.tsv`, each equal to its sha256 in `SHA256SUMS`); it stops on any difference. No count leaves out the constant-model datasets, since the counts are read once with no other subset. A dataset left out of the second reading, for either reason, is not in that reading's counts.

**Small lists.** A dataset enters the second reading only when it passes the floor and each class lists no genome or at least 5 (`SMALLEST_LISTED`, read from `floor.tsv`'s listed counts by class). Every other dataset that passes the floor is left out of it for a short list and read once, in the first.

**Smallest cell.** Five is a convention for a smallest cell. It does not prove that no genome's call can be worked out: a class of 5 genomes all called alike still shows each one's call.

**Dropped rows.** `second.py` still writes the rows of a dataset left out for a short list. `m_reading.py count` stops without them, then drops them before any count, table or line reads them.

**Cause.** In a dataset that lists 1 genome, the first reading less the second gives back that genome's class, calls and scores. A table that counts a dataset and hides its figures can be worked back: the re-check of this code, on made-up species directories, named 44 hidden figures from the medians and interval ends shown, and 42 of the 44 were true. Leaving the dataset out ends that, since no table then holds any figure of its second reading.

**Tables.** `m_reading.py count` writes five tables, each of counts per dataset and no genome:

- `m-count.tsv` and `m-conditions.tsv`: each count and reading;
- `m-both.tsv`: per arm, conditions 1 to 5 of ruling 15 on the count over all datasets in each reading, and whether both readings meet each, as a claim on all datasets needs;
- `m-size.tsv`: per reading, arm (the twin too) and learner, the datasets entered, better, worse and the same, the counts the Pick rule's size rule reads; it decides no pick;
- `m-datasets.tsv`: per reading, dataset, split, pair and learner, with each dataset's listed counts by class and its second reading's state; the second reading gives only the datasets in it.

**Page tables.** `scripts/c3-inputs/m_page.py` keeps the record and prints the page's tables. `floor` copies `m-floor.tsv` and `m-left-out.tsv` into `reference/input-arms-m` before any score is read, and `record` the other 5 tables after the count. `print` gives each table in `page_tables.py`'s number forms; a dataset left out reads its reason in place of its second reading's figures. `onpage` checks the left-out table's rows on a page. Each stops on a file missing or changed, a count over the datasets read, or a genome id's shape.

**Result page.** `m_page.py page` prints what `docs/input-arms-result.md` takes after the count, in four tables: per arm, conditions 1 to 5 in both readings, then the datasets read by species, the count over all datasets and the count on M alone. Every other table stays in the record's files, among them each count on every measure and one row per dataset, arm and learner. `test_m_record.py` fails on any page that holds one of those tables' header rows.

**Pack.** Only `OUT/pack.tar` and `OUT/DONE` leave the box. `pack.py` packs every file under OUT but `work/` and `flags/`: the result files (counts per fold), the build records, `floor.tsv`, the checks, both readings' tables and the logs. The pack holds both readings, as the count needs them. So it goes unopened to the one machine that takes the count over all 9 species, and to no other; it is no page file.

**Pack check.** `pack.py check` refuses the pack, naming the path and why with no genome id, when it would carry a link, a score file (by its name or its format), a study fit's fold result, a report or a list of samples, a list file of genomes, a compressed file, a file that is not UTF-8 text or holds a NUL byte, or a genome id of the species in a path or a file.

**Ids.** An id is read after "_" or "." and inside a longer run of points and digits, and each log is copied with every such id written as `[genome id]`. The tar holds the checked files alone, in name order, with no user, group or time stamp. Only the count tables and the page tables leave the count machine.

**Cost.** `cost.py measured` estimates the run, before it starts, from the first run's own step times in `reference/input-arms-m/first-run-times.tsv`: 592 steps of the 9 species, in whole seconds, each row naming its box's `logs/timing.tsv` by sha256. C3's budget is 45 xlarge-hours.

- With the program's tables reused on the 8 species with a dataset of M in set A, whose first-run boxes hold them: 12.32 xlarge-hours, 197.1 core-hours, $2.46. *Mycobacterium tuberculosis* has no dataset of M in A and needs no table.
- With the program run again on those 8: 17.02 xlarge-hours, 272.3 core-hours, $3.40.
- `cost.py model`, the entry C3 arms' cost model, gives 44.76 xlarge-hours central and 102.77 xlarge-hours upper. It stands only as an upper bound: its arm figures take the slower of two pilots, and each matrix arm the larger of its samples' and genomes' share.

**Scaling.** Each lineage fit is scaled by its arm's samples of M over its first-run samples, and each study fit by the same counts in set S. *Enterococcus faecium* has datasets of M in S and no first-run study fit, so its study fits are taken as its seed-42 lineage fit per sample. The grm build counts twice, and the check species' first-run fits run again.

**Left out.** The estimate names as left out `m_list`, Same counts, `score.py`, `second.py`, the pack, the `git archive` copies, and box start-up and idle time.

**Old figure.** The `cost.py` the review read gave 12.468 xlarge-hours ($2.49) from the same 9 timing tables. Five changes take 0.148 off it, to 12.320:

- *Neisseria gonorrhoeae*'s program counts as tables copied, 6 s, not run again, 981 s: less 0.271;
- each study fit is scaled by its samples in set S, not all samples of M: plus 0.136;
- a species with no dataset of M in S runs no study fit: less 0.018;
- *Enterococcus faecium*'s study fits are counted: plus 0.007;
- whole seconds: less 0.002.

**Tests.** `tests/c3_inputs/test_m_*.py` holds 195 tests of this run, on SYNTHETIC fits, score files, reports, lists, species directories, tables, logs and timing copies, and on the committed `lists.tsv`, `page-binaries.tsv`, `page-datasets.tsv`, `page-count.tsv`, `page-conditions.tsv`, `first-run-times.tsv`, `reference/input-arms-m`, C3 arms entry, this page and `docs/input-arms-result.md`. `scripts/c3-inputs/m_faults.py` plants 170 faults, at least one for each new rule, in a copy of the tree and runs the test meant to catch each.

## What it does not do

It reads no sealed, held-out, clinical or spent test set: neither the outside-study genomes nor the other-database genomes. It finds, names and ranks no resistance marker; a program column is named only as the program's table names it, and no column name, k-mer, unitig, protein k-mer or gene call leaves a fit or a build. It changes no report and no release rule. A winning program-call arm would enter only with the program's version pinned, and the maintainer would rule how GRM may name its calls.
