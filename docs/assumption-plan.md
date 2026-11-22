# Assumption plan

Written on 21 Nov 2026. Part 1, the inventory, lists each assumption a GRM result depends on. Part 2 fixes how each one gets tested, before any fit. Nothing was fitted or compared to write it. The commit that adds this page fixes its time.

Amended once, on 21 Nov 2026 at 23:48Z, before any fit: block T1, what the largest arm leaves out, Holm's family, the twin, and the label count. No arm, set, measure or cut point changes. This page is the public copy of the working plan, sha256 `2852db4bef1be2743a7f4aa7986a82707de310f284aa69c31425281556180ad9`. Roles replace names, quoted typing is tidied, and the study-split row is corrected to its source.

**All data.** The maintainer asked to "use all data". The largest arm stops at 4,096 genomes a species, since the matrix stages past that are not reviewed. 7 species hold more, in 158 of the 274 datasets. At 4,096, 41,627 of 80,967 genomes with a scorable label are left out, and 94 of the 274 datasets lose at least one. The release's cap of 2,000 leaves out 57,206.

It reads `main` at `77341aa` (merge of pull request #141) and release `rev14-274-20261112`, evaluation revision 14. The release holds 274 dataset reports in 23 species, and 268 of those reports have a lineage figure.

Every result in it is a retrospective analysis of public BV-BRC genomes. A split by study inside one database is never external validation. A threshold here is an adapted research screen. It is never an applicable FDA validation criterion, because no clinical claim rests on it.

Every arm reads only the release's public genomes and the AMRFinderPlus program's tables made from them. No sealed, held-out or clinical data enters an arm. No new marker is found, named or ranked, and every output is an aggregate. A gene is named only as a public catalogue names it.

The maintainer's questions, 21 Nov 2026: "so we only did 31-mers? … fix all our bad assumptions, like we probably borrowed things from Kover that we probably shouldn't have? the 25 S to 25 R thing was one? … should we improve all of this to more SOTA techniques? and utilize all data sets? and build hyperparameters … there were interpretability elements I wanted to solve".

The ruling this plan follows: each idea is a hypothesis with a measured answer. The release keeps its rule until a comparison shows a difference and the maintainer says to change it. Runs past 2,000 genomes are allowed inside a comparison, and the release pipeline keeps its cap.

## Inventory

The inventory holds 22 assumptions, A1 to A22. "Defended" means a page gives a reason or a record has measured it. 9 of the 22 assumptions have neither and are marked **none**: A1, A2, A6, A8, A12 on real genomes, A15, A16, counted in pull request #144 and not yet compared, A18 on lineage and study splits, and A20. Two more are defended only in part. No lineage-adjusted screen has been measured (A9), and no page gives a reason for picking the served learner on the random split (A19). "Here" means a choice made in this repository with no outside source.

| # | Assumption | Set at | Stated on | Came from | Measured | Defended |
|---|---|---|---|---|---|---|
| A1 | 25 resistant and 25 susceptible genomes per dataset | `src/curate.rs:196`, `src/main.rs:40` | README, `docs/input-configuration.md` | Here. The legacy app at `7711459` used 50 R and 50 S records (`legacy/src/app.py:3494-3501`). No commit links the two | nothing | **none** |
| A2 | The 25 is one number at two inclusion gates, with no separate reporting rule | `src/curate.rs:666`, then `src/analysis.rs:1311-1312` after download | README | here | comparison records reuse 25 to score; the time split used 20 and 25 | **none** |
| A3 | At most 2,000 genomes per species, not per dataset | `src/main.rs:36` | `docs/input-comparison-runner.md` | here | H-1: cap-1000 worse on 25 of 33 datasets; cap-4000 better on 31 and worse on 18 of 49 (sign p 0.0854) | measured, 4 species |
| A4 | Capped genomes ranked by summed 1/(class size) over their reported calls, ties by id | `src/curate.rs:753` | code comment | here | H-1 grr-v1: 15 better, 18 worse of 33 | measured, 4 species |
| A5 | k = 31 | `src/main.rs:695` | `docs/kover-comparison.md` | Kover's default and the legacy app's | H-1: k21 20 better, 29 worse of 49 (p 0.2529); k25 19 and 30 (p 0.1524) | measured |
| A6 | Canonical k-mers, presence only, no counts | `src/dna.rs`, `src/ablate/blocks.rs:207` | `docs/input-configuration.md` | Kover | nothing | **none** |
| A7 | Identical presence patterns become one feature | `src/matrix.rs` header | the same | here | nothing | stated, exact for presence learners |
| A8 | Floors: a k-mer in under 2 genomes is dropped; a fit's candidates need support of max(3, n/200) | `src/main.rs:698`, `src/learn/eval.rs:1513` | `docs/input-configuration.md:69` | Kover for the first, here for the second | nothing | **none** |
| A9 | Forest, boosted trees and lasso see the top 4,000 features by a 2x2 chi-square on their fitting samples, unadjusted for lineage | `src/learn/mod.rs:487`, `src/main.rs:701` | `docs/input-configuration.md` | here | H-1: budget-1000 worse on 37 of 49; budget-16000 no direction. A median dataset has 479,983 features | budget measured, adjustment **none** |
| A10 | SCM: at most 6 rules, 7 p values, both kinds, chosen by inner 3 folds | `src/learn/eval.rs:105`, `src/learn/scm.rs:218` | `docs/kover-comparison.md` | Kover, except 6 rules where Kover's default is 10 | Kover minus GRM scm −0.052 (−0.094 to −0.018) on 12 datasets | measured |
| A11 | CART: depth 1 to 8 by inner folds, leaf 1 | `src/learn/eval.rs:106` | `docs/hyperopt-selection.md` | here | synthetic cohorts only | synthetic only |
| A12 | Forest: 200 trees, subsample 0.632, depth 30, leaf 1, about 64 of 4,000 features a split. Only the screen is chosen | `src/learn/tree.rs:413-417`, `:431` | `docs/hyperopt-selection.md` | here | synthetic only; 20 of H-1's 23 firm results are rf | **none** on real genomes |
| A13 | Boosted trees: fixed settings, rounds by early stopping on a 15% holdout | `src/learn/gbt.rs:74`, `src/learn/eval.rs:643` | `docs/hyperopt-selection.md` | here | synthetic only | synthetic only |
| A14 | Lasso: penalty from a 30-point path on a 20% holdout, at most 400 nonzero | `src/learn/eval.rs:110`, `:645`, `src/learn/lasso.rs:272` | `docs/input-configuration.md` | here | synthetic only | synthetic only |
| A15 | A score of 0.5 or more is a resistant call | `src/learn/eval.rs:37-45` | `docs/prediction-contract.md` | here | nothing | **none** |
| A16 | Reported R and S only. S with I is S, I with R is R, I alone is set aside, Nonsusceptible is R | `src/curate.rs:34-36`, `:568`, `:227` | the methods page (`web/src/views/MethodsView.tsx:202-204`), silent on S with I and I with R | here | Pull request #144: 1,320 of 1,064,823 pairs mix I with S or R, and they give 853 of the reports' 217,514 labels | **none**: counted, not compared |
| A17 | A label read from an MIC is never scored. Its breakpoint is fitted on the whole cohort and used at 95% agreement. No breakpoint table | `src/curate.rs:151-174` | `src/curate.rs` header | here | 295,219 pairs labelled from an MIC, none scored | stated |
| A18 | MIC models scored on random stratified folds only | `src/learn/eval.rs:1979-1992` | code only | here | `grm baselines` on the release: EA-exact median 0.848 against 0.722 for the training median, random folds | **none** on lineage or study |
| A19 | Served learner: best random-split balanced accuracy, a tie to the later learner | `src/report.rs:143` | `docs/fold-draws.md` | here | pick stable on 24 of 78; optimism 0.0071; best of five on lineage in 100 of 235 | measured, reason **none** |
| A20 | Lineages: Mash on a 1,000-value sketch, average linkage, the coarsest cut with 20 or more clusters and the largest at 25% or less | `src/species.rs:28`, `:502`, `src/cluster.rs:235-239` | code only | here | nothing | **none** |
| A21 | QC: at most 1,000 contigs, completeness 90 or more, contamination 10 or less, length within 25% | `src/curate.rs:189-199` | `docs/prediction-contract.md:46` | here | nothing | stated, no reason |
| A22 | One fold draw: 5 outer, 3 inner, seed 42 | `src/learn/eval.rs:103-104` | `docs/fold-draws.md` | here | lineage SD over draws 0.0083 to 0.0140 per learner | measured |

## Kover items

| Item | Kover 2.0.7 default | GRM | Source |
|---|---|---|---|
| k | 31 | 31 | `legacy/bin/kover/kover --kmer-size` |
| Canonical k-mers, presence only | yes | yes | `src/dna.rs` |
| Singletons left out | yes | yes, floor 2 | `src/main.rs:698` |
| SCM with utility covered − p·lost | yes | yes | `src/learn/scm.rs` header |
| Both kinds, conjunction and disjunction | yes | yes | `src/learn/scm.rs` |
| p values | 8, with 999999 | the same 7 without 999999 | `src/learn/scm.rs:218`; the legacy app listed 10 |
| Choice by cross-validation | yes | yes, inner 3 folds | `docs/kover-comparison.md` |
| Most rules | 10 | 6 | `src/learn/eval.rs:105` |

Not from Kover: the 25/25 line, the cap, the pattern grouping, the screen and the other four learners. `docs/kover-comparison.md` lists what else differs from Kover: the features offered, the rule tie, the inner folds, the selection criterion and the fewest-rules rule.

## Measured records

| Record | Datasets | Finding | Spent |
|---|---|---|---|
| H-1 input comparison | 72 in 4 species, 11 arms, 2 twins | No arm beats the default across the four species. P1, cohort size, failed | no |
| Study split | 70 in 9 species | Study folds cost specificity and log loss for 5 of 5 learners, and recall for 2 | no |
| Time split | 24 in 4 species | Training on the past alone reads worse on log loss and balanced accuracy for 5 of 5 | no |
| Other public database | 13,377 genomes, 8 species | SCM meets the rule in 6 to 8 of 8 species by reading, CART 5 to 7, lasso 3 to 4 of 4 | **yes, read once** |
| Outside-study test | 7,258 genomes, 81 readable datasets | SCM holds on 32 of 81, CART on 23, lasso on 7 of 25 | **yes, read once** |
| Kover | 12, one per species | Kover minus GRM scm −0.052 (−0.094 to −0.018) | no |
| Dense network | 24 | −0.059 (−0.132 to −0.044) against the best of five in hindsight | no |
| TBProfiler | 14 graded M. tuberculosis reports | Higher than scm on 8, cart 9, rf 13, gbt 10, lasso 10 | no |
| AMRFinderPlus | 235, primary-new-235 | Ahead of GRM's own reader of the AMRFinderPlus catalog by 0.10 or more on 48, 20 of them quinolones. Each learner's median gap to the program is within 0.007 | no |
| Fold draws | 78 | The pick is stable on 24 of 78, with optimism 0.0071 | no |

The two spent test sets are separate from each other and from the release's cohorts: 13,377 genomes from another public database and 7,258 newly deposited BV-BRC genomes from studies the cohorts do not hold.

## Rules

These rules hold in every block below. A block says "Rules'" where it uses them unchanged.

| Rule | Plan |
|---|---|
| Default arm | The release's rule at `77341aa`, refitted in the same run with the same binary |
| Twin | The default at seeds 45, 46 and 47. Its count is the count with no change. H-1's tuning twin moved 102 of 230 contrasts past the margin |
| Draws | Each arm at fold seeds 42, 43 and 44. A dataset's value is the mean of the three |
| Splits | Lineage on every dataset with a lineage figure. Study on those passing `docs/study-split.md`'s entry rule. The random split is shown and decides nothing |
| Inside the folds | Every screen, setting, cutoff and learner choice uses the outer fold's training samples only |
| Primary measure | Lineage balanced accuracy, per learner |
| Beside it | Recall, specificity, the very major error rate (resistant called susceptible), the major error rate (susceptible called resistant), each with a Wilson 95% interval, and log loss |
| Interval | Paired by dataset. The median difference gets a 95% interval from 4,000 species resamples, read at the 101st and 3,900th values, as the AMRFinderPlus record does |
| Same | Within 0.02 either way in the primary measure, log loss and EA-exact, and within 0.05 in overlap |
| Margin basis | One draw's SD is 0.0083 to 0.0140, so a difference of two three-draw means has an SD near 0.007 to 0.011 |
| Spent test sets | No comparison reads the outside-study or other-database genomes. A read after a ruling would be a second look, and the maintainer decides it |
| Labels | Training labels come from the rebuilt revision 14 table alone. M marks a dataset holding a pair with S and I or with I and R: 36 of the 274 |
| Running | Long fits run unattended on rented machines, from a script that writes its result and then stops its machine |

The release changes only when all six hold:

1. On the lineage split, the median difference's interval lies wholly above 0 for each learner the change touches, or for at least 3 of 5 when it touches all.
2. Better beats worse in the count, by a one-sided sign test at p < 0.05 after Holm over the block's arm and learner pairs, the twin's included. Blocks are separate families.
3. The twin meets neither 1 nor 2. It is the default fitted again with only the fold seeds changed, so its count is what moves with no change. If it meets either, seeds alone look like a gain, and the block decides nothing.
4. No touched learner's recall interval lies wholly below −0.02.
5. The study-split median has the same sign for the touched learners.
6. The maintainer rules the change.

**Labels.** Fits start on the 238 datasets that are not M, where three of the four label rules give the same labels. The 36 wait for the maintainer's rule. Each result names the rules it does not depend on. Under the most resistant call, 33 of the 238 hold a pair with R and S and are refitted where labels change.

The release does not change for a gain on the random split only, on log loss only, on one species, on a set chosen by its result, or on single contrasts the count does not carry.

| Set | Rule | Datasets | Species |
|---|---|---:|---:|
| L | reports with a lineage figure | 268 | 23 |
| S | `docs/study-split.md` entry rule: 5 or more BioProjects, 0.8 of genomes in one, none over 0.5 | 70, 69 in L | 9 |
| CAP | L in the 8 species whose cohort reached 2,000 | 166, 50 in S | 8 |
| BIG | L with 200 or more genomes in the smaller class | 108, 36 in S | 11 |
| A | L with program tables and an organism option | 237 | 21 |
| MIC | reports with an MIC model and a lineage figure | 132 of 134, 43 in S | 12 |
| MISS | A where the program beats a learner by 0.10 or more on lineage | 75 | 16 |

Each set splits into M and not M, and each block gives its M count. The sets overlap: S, CAP, BIG, A, MIC and MISS each lie inside L, apart from the one S dataset with no lineage figure.

## Diagnostic

D0 runs before any comparison. It adds no feature and no fit. It asks why a learner misses where the program does not, and its answer orders C3, C4 and C5.

The AMRFinderPlus record's 48 datasets are where the program beats GRM's own reader of the AMRFinderPlus catalog. The learners are a separate comparison. On those 48, they trail the program by 0.10 or more on only 6 (scm), 4 (cart), 15 (rf), 5 (gbt) and 3 (lasso). So D0 runs on each learner's own miss set. The 48 get a second question: which catalogue element GRM's known-gene rule does not read there, such as a point mutation, a subclass or a family.

On each lineage fold where a learner missed, D0 puts the dataset in one of four causes:

1. **Absent.** No feature matches the program's call vector on the fold's training genomes at 0.95 agreement or more.
2. **Removed.** A feature matches, but its support is under the floor (A8), or its chi-square rank on the fitting samples is past 4,000 (A9, for forest, boosted trees and lasso).
3. **Lineage.** A matching feature is offered, but the learner chose a feature that agrees with it at 0.95 or more on the training genomes and lies in at most two lineage clusters.
4. **Unused.** The learner had the input and did not choose it.

| D0 | Diagnostic |
|---|---|
| Question | Where a learner misses and the program does not, why? Asked: "fix all our bad assumptions" |
| Arms | None. It reads the release's fold models, matrices and screens |
| Datasets | MISS, 75 in 16 species, each learner on its own set: scm 26, cart 34, rf 59, gbt 48, lasso 37 |
| Splits | The release's lineage split only |
| Measure | Datasets in each cause per learner, and each cause's share of the learner's gap. A dataset goes to the cause covering most of its missed folds, and a tie counts in both |
| Release | Changes nothing. Its sets are chosen by their result, so they never test an arm |
| Ordering | The cause with the most datasets over the five learners sets the order if it holds 40% or more. Absent puts C3 first, Removed C4's floor arm, Lineage C4's adjusted arms, Unused C5. Otherwise C4, C3, C5 |
| Spent sets | Reads neither |
| Cost | Matrices for 16 species, at most 3.44 xlarge-hours, and the matching, under 0.5. About 4 xlarge-hours, 64 core-hours, $0.80 |
| Code | Needed: a reader matching the program's per-genome calls to matrix features and ranking them in each fold's screen. Exists: the program's tables, `fold_models`, `grm export` |
| Order | Now. It fits nothing, so it does not wait for the label count |

## Comparisons

Each block has the same rows. Costs come from the units in section Cost.

**C1 cutoff**

| C1 | Cutoff |
|---|---|
| Question | What does a lower line let in, and how far can an estimate from a small cohort be trusted? Asked: "the 25 S to 25 R thing was one?" |
| Arms | Part A counts what 5, 10, 15 and 20 per class admit against 25, and 50, the legacy line. Part B draws cohorts of r resistant and r susceptible genomes, r = 10, 15, 20, 25 (default) and 50, five draws each |
| Datasets | Part A: the whole curated corpus. Part B: BIG, 108 datasets in 11 species, marked M as the label count finds |
| Splits | The small cohort's own lineage split gives the estimate. The truth is that model scored on the dataset's other genomes, all of them and those in lineages the cohort lacks. Study split on the 36 BIG datasets in S |
| Measure | The estimate's error, its absolute difference from the truth in balanced accuracy, per learner. Beside it, how often the interval covers the truth, and Rules' measures |
| Count | For each line against 25: datasets where the error is smaller by more than 0.02 (better), larger by more (worse), or within it (same) |
| Release | A lower line is a candidate when the median of its error minus 25's has a species interval wholly below 0.02 and its coverage is within 0.05 of 25's. The maintainer then decides whether to fetch and fit what part A counts |
| Spent sets | Reads neither |
| Cost | Part A under 0.1 xlarge-hours. Part B: 2,700 evaluations of two fits at about 0.0027 core-hours each, the study run's smallest fit (2.39 s at 4 threads), plus matrices, at most 3.44 xlarge-hours. About 4.5 xlarge-hours, 72 core-hours, $0.90 |
| Code | Needed: a subsample arm for `grm arm` drawing r per class within lineage groups, and a scorer for the left-out genomes. Exists: `grm curate --min-class`, final-model scoring from the outside-study path |
| Order | Part A now, with D0. Part B after the label count |

**C2 cap**

| C2 | Cap |
|---|---|
| Question | Does using more of each species' genomes help, scored on the same held-out genomes? Asked: "utilize all data sets?" |
| Arms | 2,000 (default), 4,000, and all labelled genomes up to 4,096. Past 4,096 waits for the review in `docs/matrix-scaling.md`. The ranking (A4) is not rerun, since grr-v1 showed no direction |
| Datasets | CAP, 166 in 8 species. 7 of the 8 species have more than 4,096 labelled genomes, M. tuberculosis the most at 32,058. E. faecium, at 3,138, fits whole |
| Splits | Lineage and study, 50 CAP datasets in S. Outer test folds come from the 2,000 cohort and stay the same in every arm (`common_outer_test_set_balanced`). Added genomes join training only, never from a held-out lineage |
| Measure | Rules' primary measure and those beside it |
| Count | Rules' count, per learner, for 4,000 and for all against 2,000 |
| Release | Rules' change rule on CAP. H-1 found cap-4000 better on 14 of 18 Salmonella datasets and 17 to 14 on the other three species, so one species alone changes nothing |
| Spent sets | Reads neither |
| Cost | Default 166 × 3 × 0.022 = 11.0 xlarge-hours, twin 11.0, 4,000 about 21.9, all to 4,096 about 22.5. 66 xlarge-hours at H-1's clean rate and 99 at its full rate, 1,060 to 1,590 core-hours, $13 to $20 |
| Code | Needed: the review past 4,096 for the full arm. Exists: the cap arms and common test-set pairing in `scripts/input-comparison/run.py`, and `--disk-sparse` |
| Order | After the label count. Does not wait on D0 |

**C3 inputs**

| C3 | Inputs |
|---|---|
| Question | Do other encodings separate the classes better than 31-mer presence? Asked: "so we only did 31-mers?" and "more SOTA techniques" |
| Arms | K31 (default). Unitigs from a compacted graph of the cohort's 31-mers. K31 plus k-mer counts. K31 plus amino-acid k-mers of called genes. The program's calls alone. K31 plus the program's calls. k 21 and 25 are not rerun |
| Datasets | L, 268. The two program-call arms run on A, 237 |
| Splits | Lineage on all, study on the 69 in S |
| Measure | Rules', per learner |
| Count | Rules', each arm against K31 |
| Release | Rules' change rule. A winning program-call arm enters only with the program's version pinned, and the maintainer rules how GRM may name its calls |
| Spent sets | Reads neither |
| Cost | Default and twin 9.8 xlarge-hours. Program-call arms 4.9 each. Each new encoding 3.44 for matrices plus 4.9. About 45 xlarge-hours, 720 core-hours, $9. A one-species pilot fixes the build times first, since no record holds them |
| Code | Needed: unitig, count and amino-acid sources for `grm arm` (N is unavailable today, `src/ablate/blocks.rs:207`), and a source reading the program tables. Exists: the block framework, K, C, H and P, the program's tables |
| Order | Set by D0. Its unitig source comes before C7's elastic net |

**C4 screen**

| C4 | Screen and floor |
|---|---|
| Question | Does the screen or the floor remove the inputs that separate the classes, or keep lineage markers? Asked: "fix all our bad assumptions" |
| Arms | Chi-square top 4,000 (default). Lineage-adjusted top 4,000, a Cochran-Mantel-Haenszel test with the training genomes' lineage clusters as strata. Lineage-adjusted top 16,000. Candidate floor fixed at 3 in place of max(3, n/200) |
| Datasets | L, 268 |
| Splits | Lineage, and study on the 69 in S. The strata come from sketches, which carry no label |
| Measure | Rules'. Screen arms: forest, boosted trees and lasso, with SCM and CART refitted to check they do not move. Floor arm: all five |
| Count | Rules', each arm against the default |
| Release | Rules' change rule |
| Spent sets | Reads neither |
| Cost | Default and twin 9.8 xlarge-hours, three arms at 4.9 each. About 25 xlarge-hours, 400 core-hours, $5 |
| Code | Needed: an adjusted screen in `src/learn/mod.rs` behind a flag, and a floor flag. Exists: `--screened-features`, and an earlier lineage-adjusted scan, a guide only since it failed its own negative control (28 of 200 permutations, cut 15) |
| Order | Set by D0 |

**C5 settings**

| C5 | Settings |
|---|---|
| Question | Does a search over each learner's settings, made inside the training folds, beat the release's rule? Asked: "build hyperparameters" |
| Arms | The release's rule (default), which already chooses some settings inside each fit. A nested search of 30 random settings per learner from the spaces frozen in `scripts/hyperopt/protocol.py`, by inner 3-fold balanced accuracy. The call cutoff chosen on inner scores in place of 0.5 |
| Datasets | The 69 datasets in both S and L. The search costs about 90 times the release's rule, so it runs where both splits exist |
| Splits | Lineage and study, with inner folds drawn as the release draws them |
| Measure | Rules', per learner |
| Count | Rules', each arm against the default |
| Release | Rules' change rule. On synthetic cohorts the default-anchored rule lost 0.0126 (0.0088 to 0.0163) for SCM on the lineage split, so SCM's result is read beside that |
| Spent sets | Reads neither |
| Cost | 69 × 2 splits × 3 draws × 90 × 0.019 core-hours, about 708 core-hours or 44 xlarge-hours. Default, twin and cutoff add about 1.5. About 46 xlarge-hours, 735 core-hours, $9 |
| Code | Needed: a nested-search driver over `grm learn --fold-file` that writes its choices, and a cutoff option. Exists: the frozen spaces, the fold files, `--fold-file` |
| Order | Set by D0, else after C4 and C3 |

**C6 labels and MIC**

| C6 | Classes and MIC |
|---|---|
| Question | Does keeping Intermediate as a class, or modelling the MIC, hold on lineage and study splits? Asked: "fix all our bad assumptions" |
| Arms | A: two classes, I set aside (default); three ordered classes as two binary fits with the five learners; I as R. B: the released censored MIC model, now on lineage and study folds; an ordered model over dilution steps; the training median (`docs/baselines.md`) |
| Datasets | Part A: L datasets with 25 or more genomes whose only call is I. Their count needs the label count's table. Part B: MIC, 132 |
| Splits | Lineage and study, 43 MIC datasets in S. The breakpoint for category agreement is fitted on training folds only |
| Measure | A: balanced accuracy of the R against S call on R and S genomes only, so all arms score the same genomes, with the share of I genomes called R beside it. B: EA-exact, with category agreement and both error rates beside it |
| Count | Rules', each arm against its default |
| Release | A: Rules' change rule, after the maintainer rules the label pattern. B: if the model's lineage EA-exact beats the training median by under 0.02 on most datasets, the site shows the lineage reading beside the random one. The maintainer rules that wording |
| Spent sets | Reads neither |
| Cost | A: at most 2 arms × 2 fits × 4.9, plus 9.8, so 29 xlarge-hours over all of L. B: at most 3 runs × 3 draws × 1.63, so 15. At most 44 xlarge-hours, 700 core-hours, $9 |
| Code | Needed: a three-class label mode, the MIC model on lineage and study fold plans, an ordered MIC model. Exists: `evaluate_mic`, and the constants and scorer of `grm baselines` |
| Order | B first among fits, since the gap is in a figure the release already shows. A after the label count |

**C7 outside learners**

| C7 | Outside learners |
|---|---|
| Question | Do public learners do better on GRM's own folds and screens? Asked: "more SOTA techniques" |
| Arms | GRM's five (default). XGBoost and LightGBM, versions pinned, on the export's 4,000 screened rows, with rounds by early stopping on a 15% inner holdout as GRM's boosted trees do. An elastic net on unitigs, with mix and penalty chosen by inner 3 folds |
| Datasets | L, 268 |
| Splits | Lineage, and study on the 69 in S, from `grm export`'s fold plans |
| Measure | Rules'. Each library is paired with GRM's boosted trees, and the elastic net with GRM's lasso |
| Count | Rules' |
| Release | Rules' change rule, and the maintainer's ruling on adding a dependency |
| Spent sets | Reads neither |
| Cost | Default and twin 9.8 xlarge-hours. Libraries: 337 split-datasets × 5 folds × 3 draws × 2 × 0.0026 core-hours, 26 core-hours. Exports at most 3.44. About 15 xlarge-hours, 240 core-hours, $3, plus the elastic net after a pilot |
| Code | Needed: a runner that reads exports and writes reports in GRM's form with pinned libraries. Exists: `grm export`, `--fold-file`, the dense-network runner |
| Order | Libraries any time after the label count. The elastic net after C3's unitig source |

**C8 stability**

| C8 | Stability |
|---|---|
| Question | Do the learners choose the same inputs across folds and draws, and how many lie in catalogue-named genes? Asked: "there were interpretability elements I wanted to solve" |
| Arms | The release's rule from its saved `fold_models`, with no fit. Then each arm above, from its own reports |
| Datasets | All 274 for the release, then each arm's own set |
| Splits | The five lineage and five random folds, and the three draws for arms |
| Measure | Per dataset and learner, the median pairwise overlap (Jaccard) of chosen feature sets across lineage folds. Beside it, the count of distinct features chosen and the count whose notes carry a catalogue family. Only catalogue family names appear |
| Count | Each arm against the release: datasets with overlap higher, lower or within 0.05 |
| Release | Nothing alone. An arm that meets Rules' change rule while its overlap falls on most datasets goes to the maintainer with that said |
| Spent sets | Reads neither |
| Cost | Under 0.1 xlarge-hours for the release, and inside each arm's cost after |
| Code | Needed: a reader over `fold_models` and the notes' `catalog_families`. Exists: both, in every report |
| Order | Now, with D0 |

**C9 pick**

| C9 | Served learner |
|---|---|
| Question | Is the random-split pick, with ties to the later learner, the right learner to serve? Asked: "fix all our bad assumptions" |
| Arms | The release's pick (default). A pick made inside each outer training fold by inner lineage-grouped folds over the five learners |
| Datasets | L, 268 |
| Splits | Lineage, and study on the 69 in S |
| Measure | Lineage balanced accuracy of the picked learner on the outer test folds, with Rules' measures beside it |
| Count | Rules' |
| Release | Rules' change rule, and the maintainer's ruling, since the site names the served learner |
| Spent sets | Reads neither |
| Cost | Three inner fits of a draw's lineage half, 3 × 0.82 = 2.4 xlarge-hours a draw, 7.3 for three. Default and twin 9.8. About 17 xlarge-hours, 270 core-hours, $3.40 |
| Code | Needed: an inner pick behind a flag, beside `best_learner` in `src/report.rs`. Exists: `best_learner` and the fold draws' optimism reader |
| Order | After the label count. Does not wait on D0 |

## Later-records test

| T1 | Later-records test |
|---|---|
| Question | Do the release's rule and each passing arm hold on labels deposited after the release's download? |
| Set | The 252,148 BV-BRC lab records inserted at or after 2026-11-05T19:10:00. File `later-records.jsonl`, sha256 `2efc5328ce4267a7963b93578899545758035d4b6bf004034ac41596d3f6a7dc`, sealed on 21 Nov 2026 at 23:18:12Z |
| Holds | 155,448 later R or S labels under today's rule, 101,593 of them in 166 of the 274 datasets. None is on a genome in a release report. Counts: `docs/assumption-plan-later-records.md` |
| Closed | Until the scoring, no fit, screen, settings search or pick reads the file, a table built from it, or a genome whose labels all come from it |
| Scored | Once, after every block is ruled: the release's rule and each arm that met conditions 1 to 5, on the same genomes and measure |
| Datasets | Those of the 274 with 25 or more later labels in each class: 102 under each label rule. The label rule picked applies here too |
| Measure | Balanced accuracy per dataset and learner, with Rules' other measures. An arm is read as its paired difference from the release's rule |
| Reported | Whatever it says, with the counts scored and not scored. A test on later records of the same public database, never external validation |
| Release | Changes nothing alone. A second read is the maintainer's to decide |
| Cost | 12,617 genomes, 54.85 Gbp, fetched at the end and not yet in the cost table |
| Code | Needed: scoring of saved final models. The release publishes no forest or boosted-tree model, so the default arm's refit saves them, named as a refit |

## Parked

| Item | Reason | Need |
|---|---|---|
| Pretrained protein or DNA sequence models as inputs | none is on hand | a GPU machine, model weights, and a licence that allows this use |
| A20, the lineage rule | another rule makes another test, so no paired count exists | a separate plan, with the study split as the check beside it |
| A21, QC thresholds | changing them changes the cohorts, and nothing asks for it | a separate plan |
| A4, the ranking | grr-v1 showed no direction on H-1 | nothing now |
| k 21 and 25 | H-1 showed no direction | nothing now |
| All genomes past 4,096 | the matrix stages past 4,096 are not reviewed | the review in `docs/matrix-scaling.md` |

## Order

1. Fits start on the datasets that are not M. The M datasets wait for the maintainer's label rule.
2. Now, with no fit: D0, C8 on the release, and C1 part A.
3. The code each comparison needs, each piece reviewed before it runs.
4. C6 part B, first among fits.
5. The input group in D0's order, by default C4, C3 and C5, then C7.
6. On their own machines, not waiting on D0: C1 part B, C2, C6 part A and C9.
7. The maintainer rules each result.
8. T1, once, after every block is ruled.

## Cost

One xlarge machine is 16 CPUs at $0.20 an hour. Every cost is an estimate from these units unless it says measured. Comparisons on the same set and split that run together could share one default and one twin. Each block counts its own, and the block costs sum to the total, so the total is an upper bound.

| Unit | Value | Source |
|---|---|---|
| One learning draw, 274 datasets, two splits, five learners | 1.63 xlarge-hours, 26.1 core-hours, $0.33 | `docs/fold-draws.md`: 4.892 hours for three draws |
| One fold plan on one study dataset, five learners | 0.019 core-hours | the study-split draws: 20.00 core-hours for 1,050 fits |
| One H-1 dataset-arm with its own data root, capped species | 0.022 to 0.033 xlarge-hours | H-1: 22.148 clean to 33.189 in all, over 72 datasets × 14 runs |
| Matrices for the whole release, no learning | at most 3.44 xlarge-hours | the release build's 21,232 s less its 8,831.6 s of learning |
| AMRFinderPlus program, one genome | 0.0065 core-hours | the AMRFinderPlus run: 153.88 core-hours for 23,751 genomes |
| One small outside fit on 4,000 screened rows | 0.0026 core-hours | the dense-network run: 1.84 CPU-hours for 720 fits |

| Block | xlarge-hours | Core-hours | Dollars |
|---|---:|---:|---:|
| D0 | 4 | 64 | 0.80 |
| C1 | 4.5 | 72 | 0.90 |
| C2 | 66 to 99 | 1,060 to 1,590 | 13 to 20 |
| C3 | 45 | 720 | 9 |
| C4 | 25 | 400 | 5 |
| C5 | 46 | 735 | 9 |
| C6 | 44 at most | 700 | 9 |
| C7 | 15 | 240 | 3 |
| C8 | 0.1 | 2 | 0.02 |
| C9 | 17 | 270 | 3.40 |
| Total | 267 to 300 | 4,260 to 4,790 | 53 to 60 |

## Literature

Literature chose arms and measures here. No result rests on it. Every item was read as an abstract only, through PubMed's E-utilities, and none was read in full.

| Item | Used for |
|---|---|
| Lees et al. 2018, Bioinformatics, pyseer (PMID 30535304), abstract only | the lineage-adjusted screen in C4 |
| Do et al. 2024, Comput Struct Biotechnol J, unitig pan-genome learning (PMID 38707536), abstract only | the unitig arm in C3 |
| ValizadehAslani et al. 2020, Biology (Basel), amino-acid k-mers (PMID 33126516), abstract only | the amino-acid arm in C3 |
| Nguyen et al. 2019, J Clin Microbiol, MIC prediction for Salmonella (PMID 30333126), abstract only | the MIC arms in C6, the boosting libraries in C7 |
| Hu et al. 2024, Brief Bioinform, 78-dataset benchmark (PMID 38706320), abstract only | the program-call arms in C3, lineage over random splits |
| Lüftinger et al. 2021, Front Cell Infect Microbiol, limited data (PMID 33659219), abstract only | small-cohort trust in C1, genome-distance folds |
| Hicks et al. 2019, PLoS Comput Biol, parameters for SCM and forests (PMID 31479500), abstract only | the settings arms in C5 |
