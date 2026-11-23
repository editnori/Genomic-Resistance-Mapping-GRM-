# Outside-learner result

The release fits GRM's own learners. This page asks whether learners from outside GRM, fitted on the same folds, call the same held-out genomes better: XGBoost 3.0.5 and LightGBM 4.6.0 on the 4,000 rows `grm export` screens, each read against GRM's boosted trees (`gbt`), and a scikit-learn 1.7.2 elastic net on unitigs, read against GRM's lasso (`lasso`), beside a twin. It reads block C7 of the assumption plan (`docs/assumption-plan.md`), with the method of `docs/outside-learners.md`.

The design was fixed at 07:37Z on 23 Nov 2026, the lead's rulings on its questions at 09:05Z and on the elastic net's stop at 11:26Z, before any C7 score was read (`docs/assumption-plan-corrections.md`, entry "C7 design"). The table script, its gates and this page's form were committed before any score of the run was opened.

<!-- part 2: the clock of the reading, as served-pick-result.md gives it. -->

Every split here is a retrospective split inside one public database. A split by lineage or by study is never external validation. The margin of 0.02 is an adapted research screen. It is no applicable FDA validation criterion.

"No difference shown" is the rule's reading and no statement that two learners are equal. Recall stands beside log loss wherever log loss is given.

## Reading

<!-- part 2: the result per arm, conditions 1 to 5 with their figures, from page_tables.py print [conditions] and [count]. -->

**Arms.** All three arms ran: `xgboost`, `lightgbm` and `elasticnet`.

<!-- part 2: if the elastic net is withdrawn, the sentence above alone changes, to: "The elastic net was withdrawn at <time> before any score was read, for <reason>; its pair stays in the family of five at p = 1 and holds no result." No table changes shape: its rows read withdrawn. -->

**Pairs.** Each arm touches one learner. `xgboost` and `lightgbm` are read against `gbt`, and `elasticnet` against `lasso`, each against the default's own learner at the arm's seeds, 42, 43 and 44 (`k31-at-arm-seeds`). The twin, the default at seeds 45, 46 and 47, is read against the default on each learner, `twin-gbt` and `twin-lasso`.

**Conditions.** An arm carries the change when its pair meets condition 1 (the lineage median difference's interval lies wholly above 0) and condition 2 (better beats worse by the one-sided sign test at p < 0.05 after Holm over the 5 lineage pairs). Condition 3 reads both twin pairs for the whole block: the twin meets 1 in neither pair and 2 in neither, each read apart. Condition 4: the arm's recall interval does not lie wholly below -0.02. Condition 5: its study median has its lineage median's sign. Condition 6 is the maintainer's.

**Settings.** Each booster runs GRM's settings where GRM has them, on the rows GRM's screen kept, so the pair reads the library and leaves the screen as it is. The elastic net stops at tol 1e-4 or 2,000 epochs, and every elastic net result carries that limit as `docs/outside-learners.md` words it: saga's tol of 1e-4 stops on the change in the weights and bounds no distance to the optimum, and a difference between the elastic net and `lasso` reads saga's stop as well as the learner.

The release keeps its five learners and their settings whatever the count says, until the maintainer rules. Nothing changes in the release or on the site.

## Data

**Datasets.** The first run covers the 233 datasets of L less M in 23 species, each species and antibiotic named in full in the section Datasets. The 35 datasets of M wait for the label rule. A dataset's genomes are its genome count in `scripts/c7-outside/lists.tsv`, which counts each genome once.

**Draws.** The default and the twin come from plan block C3's run, the release's rule refitted with the same binary. Each arm runs three lineage draws a dataset, at fold seeds 42, 43 and 44, on the folds `grm export --seed` writes for that seed, checked against the default's folds by class and by genome before any fit. A dataset's value is the mean of its draws, and a dataset with a draw not scored is not scored.

**Study split.** The 55 first-run datasets in S get one study run per arm at seed 42 on the reviewed `study-whole` fold file, with no twin.

**Passes.** A species whose run started with the elastic net left out restarted with it on the same run folder, keeping its finished steps and booster fits. Its first copy, moved aside when it restarted, is named and set aside by its time and never read. If the elastic net is withdrawn after some second passes ran, every species is read from the copy whose `DONE` says the elastic net was left out: the species' own folder where no second pass ran, its first copy where one did. Each second-pass copy is named with its species, counted and set aside, never opened, and `page-withdrawn.tsv` gives how many species ran a second pass.

<!-- part 2: species, results and datasets not scored, with their counts, from checks.txt and page_tables.py print [checks]. -->

## Count

The primary measure is each pair's lineage balanced accuracy on the outer test folds. A dataset is better, worse or the same by more than 0.02, compared at 12 decimals. Every difference is the arm's value less the default's; for the very major error rate, the major error rate and log loss, lower is better, so better and worse turn round. The interval is a 95% interval of the median over 4,000 resamples of species. The sign test is one-sided, and Holm runs over the 5 lineage pairs on balanced accuracy, withdrawn or not.

<!-- part 2: the family's 5 rows, from page_tables.py print [count]. -->

## Per arm

<!-- part 2: per arm, split and measure, the median with better and worse, recall beside log loss, from page_tables.py print [measures]. -->

**Folds.** Per arm and split, counts of folds over the folds of its reports: for each booster the folds that kept one tree; for the elastic net the refits that reached 2,000 epochs, the folds whose fit kept no weight, and the folds scored at their training share.

<!-- part 2: page_tables.py print [folds]. -->

**Libraries.** The versions each species' run read back, as one row when every species agrees, or a row per species.

<!-- part 2: page_tables.py print [libraries] and [binaries]. -->

**Time.** Per arm, the fits timed, with the sums and medians of their wall and CPU seconds.

<!-- part 2: page_tables.py print [times]. -->

## Beside it

These counts stand beside the count. They are in no Holm family and decide nothing.

**Own twin pair.** Each arm against its own learner's twin pair alone, condition 3 read on that learner.

**Without constant folds.** The count with every dataset that holds a draw scored at its training share left out (`count-no-constant.tsv`).

**Without unconverged refits.** The count with every dataset whose elastic net draws hold a refit that reached 2,000 epochs left out (`count-no-unconverged.tsv`).

<!-- part 2: page_tables.py print [conditions] beside rows, [no-constant] and [no-unconverged]. -->

## Seen after

These points are seen in the tables after the reading. They change no rule and no reading above, and decide nothing.

**Species.** Per species, named in full, the rows are dataset and pair counts, with their denominator. They decide nothing.

<!-- part 2: page_tables.py print [species], and what else is seen after the reading. -->

## Limits

**One database.** Every split is a retrospective split inside one public database. A split by lineage or by study is not external validation.

**Elastic net stop.** The elastic net stops at tol 1e-4 or 2,000 epochs. saga's tol of 1e-4 stops on the change in the weights and bounds no distance to the optimum, so a difference between the elastic net and `lasso` reads saga's stop as well as the learner.

**Screen.** The boosters fit the rows GRM's screen kept, 4,000 per outer fold, so a pair reads the library on GRM's screen and says nothing of another screen.

**Study split.** One study run per dataset at seed 42, and no study twin.

**Labels.** The 35 datasets of M wait for the label rule.

**Markers.** No feature, unitig or k-mer is named or ranked. No elastic net weight and no booster's feature importance leaves the box. The page holds aggregates over datasets, species, splits, arms and folds only.

## Method

**Checks of the run.** Before any score was read, `scripts/c7-outside/run_checks.py check` read each species' copy and printed counts only: the species, each `SHA256SUMS`, each `DONE` with its elastic net pass, `expect.py` against each species' results, the genome-id count and the count of unitigs and matrix column names over every file that is read or published. `score/score.log` stays out of the record.

<!-- part 2: the checks' table, from checks.txt. -->

**Reading.** `scripts/c7-outside/run_checks.py gather` copied each species' `score/` tables, as its run wrote them, and its run files into the reading, and `scripts/c7-outside/score.py --count` counted them once.

<!-- part 2: the commands, as reading.txt lists them with the sha256 of each file they wrote. -->

**Tables.** `scripts/c7-outside/page_tables.py write` checked the reading against 17 gates and wrote `reference/outside-learners/`. The four count files are `score.py`'s own tables: the script counted them again, copied them and summed them.

| Gate | What it holds |
|---|---|
| 1. reading | the reading holds `reading.txt`, `checks.txt` with every check passed, `withdrawn.tsv`, `count/` with its four files, `score/` and `run/`, and nothing else; each withdrawn arm is named with a reason and a time |
| 2. species | `score/` and `run/` each hold every species of the first run once; a folder named `<species slug>.<UTC time>` is a first copy set aside and never read as a species |
| 3. M | no dataset of M, or outside the first run, stands in any table |
| 4. pairs | `xgboost` and `lightgbm` stand on `gbt` only and `elasticnet` on `lasso` only, each read against `k31-at-arm-seeds` |
| 5. twin | the twin's rows are lineage draws at seeds 45, 46 and 47 only, and the conditions never read the twin as an arm |
| 6. study | a study row is a study run at seed 42 on the study split, and is never counted on the lineage split |
| 7. family | each count holds the five lineage pairs, withdrawn or not, each with its Holm p |
| 8. withdrawn | a withdrawn arm holds no row, no report and no fold, and its pair stays at p = 1 |
| 9. passes | a species whose `DONE` gives its first pass alone holds no elastic net result |
| 10. rows | each species holds exactly its first-run datasets' rows |
| 11. missing | a draw not scored holds no figure and is never read as zero |
| 12. draws | `datasets.tsv` is the mean of `draws.tsv`'s draws |
| 13. reports | one report per library arm draw and no other, each of its listed keys only, its fold counts the sums over its folds |
| 14. checks | every `score.py` check and every `plans.py` check passed |
| 15. binary | one code sha256 and one version of each library in every species, or each named with its species |
| 16. beside | `count-no-constant.tsv`, `count-no-unconverged.tsv` and each arm's own twin pair are written as deciding nothing |
| 17. count | the four count files are what `score.py --count` writes from the `datasets.tsv` files |

The fold counts that read 0 by construction, a booster's unconverged refits and null fits and the elastic net's one-tree folds, stay in `page-folds.tsv` and off this page.

**Tests.** `tests/c7_outside/test_page_tables.py` and `tests/c7_outside/test_run_checks.py` hold SYNTHETIC tests of the table script and the run's checks, on three made-up species. `scripts/c7-outside/page_faults.py` plants each fault of its table in a copy of the tree and runs the test meant to catch it.

<!-- part 2: the planted faults' table, caught of planted, from page_faults.py's faults.tsv. -->

**Files.** The reading's tables are in `reference/outside-learners/`. `scripts/c7-outside/page_tables.py print` prints this page's tables from those files, rounded to 3 decimals.

## Datasets

The first run's 233 datasets, per species: the datasets, those in S, and each antibiotic with the dataset's genomes. Per-dataset figures are in `reference/outside-learners/page-datasets.tsv`, with species and antibiotic in full.

| Species | Datasets | In S | Antibiotics, with genomes |
|---|---:|---:|---|
| *Acinetobacter baumannii* | 7 | 2 | cefepime (437), cefotaxime (808), doripenem (213), doxycycline (167), minocycline (139), moxifloxacin (91), sulfamethoxazole/trimethoprim (1,114) |
| *Campylobacter coli* | 8 | 0 | azithromycin (389), ciprofloxacin (391), clindamycin (389), erythromycin (391), gentamicin (391), nalidixic acid (375), telithromycin (374), tetracycline (391) |
| *Campylobacter jejuni* | 4 | 0 | ciprofloxacin (508), clindamycin (508), nalidixic acid (506), tetracycline (508) |
| *Clostridioides difficile* | 4 | 0 | azithromycin (203), ceftriaxone (136), clarithromycin (203), moxifloxacin (204) |
| *Corynebacterium diphtheriae* | 3 | 0 | sulfamethoxazole/trimethoprim (236), tetracycline (243), trimethoprim (229) |
| *Enterobacter cloacae* | 9 | 8 | cefepime (103), ciprofloxacin (110), gentamicin (104), imipenem (86), levofloxacin (104), meropenem (99), sulfamethoxazole/trimethoprim (80), tobramycin (99), trimethoprim (124) |
| *Enterococcus faecium* | 11 | 0 | ampicillin (1,680), clindamycin (901), dalfopristin/quinupristin (1,118), doxycycline (345), erythromycin (1,151), gentamicin (1,345), kanamycin (944), streptomycin (1,239), teicoplanin (1,192), tetracycline (958), vancomycin (1,741) |
| *Escherichia coli* | 29 | 15 | amikacin (1,125), ampicillin/sulbactam (283), aztreonam (699), cefazolin (414), cefepime (1,124), cefotaxime (1,149), cefoxitin (1,066), cefpodoxime (175), cefpodoxime/clavulanate (136), ceftazidime (1,466), ceftiofur (203), ceftriaxone (923), cefuroxime (548), cephalothin (472), chloramphenicol (709), colistin (239), doxycycline (227), ertapenem (820), erythromycin (205), imipenem (614), levofloxacin (403), nalidixic acid (406), nitrofurantoin (482), piperacillin/tazobactam (979), streptomycin (386), sulfamethoxazole (341), tetracycline (782), tobramycin (1,070), trimethoprim (861) |
| *Klebsiella aerogenes* | 3 | 0 | ertapenem (83), imipenem (68), meropenem (78) |
| *Klebsiella oxytoca* | 4 | 0 | cefotaxime (98), cefoxitin (99), cefpodoxime (74), ertapenem (112) |
| *Klebsiella pneumoniae* | 32 | 18 | ampicillin/sulbactam (744), aztreonam (1,050), cefazolin (1,141), cefepime (1,233), cefotetan (61), cefoxitin (1,181), cefpodoxime (221), cefpodoxime/ceftriaxone (56), cefpodoxime/clavulanate (207), ceftazidime (1,310), ceftazidime/avibactam (259), ceftriaxone (1,208), cefuroxime (565), chloramphenicol (409), ciprofloxacin (1,706), colistin (336), doripenem (243), doxycycline (141), ertapenem (952), fosfomycin (294), imipenem (1,158), levofloxacin (871), nitrofurantoin (562), ofloxacin (71), polymyxin b (100), streptomycin (80), sulfamethoxazole/trimethoprim (1,546), tetracycline (741), ticarcillin/clavulanate (250), tigecycline (924), tobramycin (1,230), trimethoprim (258) |
| *Klebsiella quasipneumoniae* | 3 | 0 | cefotaxime (101), ceftazidime (102), sulfamethoxazole/trimethoprim (101) |
| *Mycobacterium tuberculosis* | 20 | 1 | amikacin (1,813), amoxicillin (57), capreomycin (670), ciprofloxacin (101), cycloserine (401), delamanid (1,248), ethambutol (1,947), ethionamide (1,774), kanamycin (1,800), levofloxacin (1,248), linezolid (1,248), moxifloxacin (1,524), nicotinamide (166), ofloxacin (643), para-aminosalicylic acid (1,650), prothionamide (170), pyrazinamide (492), rifabutin (1,317), rifampin (1,967), streptomycin (698) |
| *Neisseria gonorrhoeae* | 6 | 0 | cefixime (600), ceftriaxone (1,926), ciprofloxacin (1,746), erythromycin (179), penicillin (1,354), tetracycline (1,371) |
| *Pseudomonas aeruginosa* | 12 | 6 | amikacin (708), aztreonam (235), cefepime (245), ceftazidime (833), ceftazidime/avibactam (309), ceftolozane/tazobactam (259), gentamicin (346), imipenem (447), levofloxacin (702), meropenem (883), piperacillin/tazobactam (210), tobramycin (645) |
| *Salmonella enterica* | 24 | 1 | amoxicillin/clavulanate (799), ampicillin (1,997), ampicillin/sulbactam (183), azithromycin (1,355), cefazolin (185), cefotaxime (941), cefoxitin (1,080), ceftazidime (783), ceftiofur (689), ceftriaxone (1,155), chloramphenicol (1,989), ciprofloxacin (1,503), gentamicin (1,675), kanamycin (584), nalidixic acid (1,599), pefloxacin (348), spectinomycin (396), streptomycin (1,269), sulfamethoxazole (275), sulfamethoxazole/trimethoprim (1,145), sulfathiazole (225), sulfisoxazole (544), tetracycline (1,869), trimethoprim (634) |
| *Shigella flexneri* | 4 | 0 | ampicillin (358), azithromycin (361), ciprofloxacin (361), trimethoprim (362) |
| *Shigella sonnei* | 20 | 0 | amikacin (580), ampicillin (1,674), azithromycin (1,307), cefoperazone (257), cefotaxime (593), cefoxitin (593), cefpirome (333), ceftazidime (582), ceftriaxone (1,163), ciprofloxacin (1,397), colistin (268), ertapenem (338), gentamicin (1,288), nalidixic acid (854), streptomycin (595), sulfamethoxazole/trimethoprim (617), temocillin (333), tetracycline (878), tobramycin (541), trimethoprim (1,104) |
| *Staphylococcus aureus* | 13 | 4 | cefoxitin (1,407), chloramphenicol (1,038), ciprofloxacin (1,624), clindamycin (1,543), daptomycin (904), erythromycin (1,957), fusidic acid (1,101), gentamicin (1,913), mupirocin (736), penicillin (1,508), rifampin (931), sulfamethoxazole/trimethoprim (1,389), tetracycline (1,647) |
| *Staphylococcus haemolyticus* | 3 | 0 | ciprofloxacin (124), fusidic acid (114), tetracycline (104) |
| *Streptococcus agalactiae* | 3 | 0 | clindamycin (141), erythromycin (141), tetracycline (141) |
| *Streptococcus pneumoniae* | 8 | 0 | azithromycin (1,297), cefuroxime (1,732), chloramphenicol (1,859), clindamycin (1,894), erythromycin (1,991), meropenem (1,808), sulfamethoxazole/trimethoprim (1,777), tetracycline (1,971) |
| *Vibrio cholerae* | 3 | 0 | nalidixic acid (91), sulfamethoxazole (74), trimethoprim (74) |
| All | 233 | 55 | |
