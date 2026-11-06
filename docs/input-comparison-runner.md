# Full-pipeline input comparison: the runner

Layth specified a full-pipeline comparison over the input data: genome and cohort count, which genomes are included, k-mer length, feature families and the retained-feature budget. Each input configuration goes through the whole pipeline and is compared with today's default. Layth delegated the open choices to main on 5 Nov 2026; main's decisions are DECISIONS V-1 (the values and pairing rules below), the rulings after lane INPUT (the balanced common test set, the lineage drop, MLST groups, block screens, and the seed and margin as decided values) and H-1 (real data comes later, on a 4-species panel). Nothing here is Layth's own acceptance.

**Built and proved on synthetic data only.** `scripts/input-comparison/run.py` runs every arm through `grm` and hands the predictions to J and K. Its end-to-end test (`tests/input_comparison/test_run.py`, `EndToEndTests`) runs all 11 arms plus the default on a synthetic data root with a planted effect. No real cohort has been run; that is H-1's, after revision 11.

[input-configuration.md](input-configuration.md) traces each input as the code has it now.

## The decided arms (V-1)

One factor at a time against the default; 11 arms plus the default. `scripts/input-comparison/spec.json` holds them, and `run.py` refuses an arm that changes no input or more than one.

| Factor | Default | Arms | Pairing rule |
|---|---|---|---|
| Genome and cohort count | cap 2,000, class minimum 25 | `cap-1000`, `cap-4000` (class minimum stays 25) | `common_outer_test_set_balanced` |
| Which genomes | `legacy` (`curate::select_genomes`) | `grr-v1`: `group-round-robin-v1` (`curate::plan_group_sampling`), seed 1, each genome grouped by its BV-BRC MLST sequence type | `common_outer_test_set_balanced` |
| k-mer length | 31 | `k21`, `k25` | `k31_groups` |
| Feature families | K | `k-c`, `k-h`, `k-p`, `k-c-h-p` (blocks K+C, K+H, K+P, K+C+H+P); N and X stay unavailable, and no source is built for them | `same_plan` |
| Retained-feature budget | `screened_features` 4,000 | `budget-1000`, `budget-16000` | `same_plan` |

Decisions stay K's `improves`, `worse` or `inconclusive`; there is no equivalence claim (contract §15). Main's rulings after lane INPUT (5 Nov 2026, under the same delegation) made two spec values decided: the group draw's seed is 1, and the log-loss margin is 0.02. The margin decides nothing: it feeds only K's descriptive within-margin note.

## Pairing rules (V-1 and the rulings after lane INPUT)

- **`same_plan`**: the same samples and the same fold plan, J's rule before V-1. Blocks and budget arms score on the default arm's data root and its plans.
- **`k31_groups`**: a k arm counts the default cohort again at its own k, but its fold plans are the default's, drawn on the k = 31 lineages and duplicate groups. So its samples, folds, plan hashes and clusters equal the default's.
- **`common_outer_test_set_balanced`**: a cap or selection arm curates its own cohort, and `ablate::plan::cohort_plan` (format `grm-cohort-plan-v2`) gives its plan.
  - **Scored set.** In each outer fold f, the arm is scored on a class-balanced subsample of the genomes both cohorts hold that the default tests in f. Each class keeps as many genomes as it can while the subsample's resistant-to-susceptible ratio matches the default's whole test fold f (`reference_test_classes`). The scarcer class relative to that ratio is kept whole. The draw is a ChaCha8 shuffle of each class, seeded from the default plan's repeat seed (`balance_seed`, recorded in the plan) through `plan::balance_seed`.
  - **Training.** The arm trains in fold f on its own cohort minus the scored genomes and minus `held_out[f]`. `held_out[f]` holds the common genomes of the default's fold f that the subsample left out (`balance_dropped[f]`). On the lineage split it also holds every genome only the arm has whose lineage (the arm's own k = 31 lineages) is a lineage of a scored genome (`lineage_dropped[f]`), so no test lineage reaches the arm's training.
  - **The rest.** A genome only the arm has trains in every other fold and is never scored. Inner folds and tuning holdouts are drawn on each fold's training samples with the arm's own lineages and duplicate groups. J restricts the default arm to the arm's scored genomes.
  - **Recorded.** `grm arm` writes the plan to `plans/<split>-<repeat>.cohort.json` and the per-fold drop counts to each run (`cohort_drops`).
- **`common_outer_test_set`**, the unbalanced rule of V-1, stays in J and K only so records written before the rulings still validate; the runner no longer uses it.

Every row of every arm carries the default plan's `fold_plan_sha256` and the default cohort's lineage cluster, since its `fold` is the default's. A cohort arm's own `CohortPlan` hash is recorded on its run (`pairing.arm_plan_sha256`). J names the rule in every contrast record (`pairing`), and K checks it against the runs and refuses any pairing outside the rules ([ablation-report.md](ablation-report.md#pairing-rules)).

## What runs

`python3 -B scripts/input-comparison/run.py SPEC --raw RAW --work WORK --grm GRM --source-commit SHA [--jobs N]`

1. **The species pool.** `run.py` curates every arm's root first. The pool is the union of their curated cohorts (genomes.txt), counted at the default's k and structured once in its own root (`roots/pool-<hash>`), before any `grm arm`. Its `structure.json` (lineages and kset digests) is the one grouping the default and the cohort arms split on (`grm arm --groups`); the other arms take the default's plans. `run.py` refuses, before counting, a genome of the union without a FASTA (naming it) and a pool whose kset estimate exceeds the disk share (below); after `grm structure`, it refuses a pool that doesn't hold exactly the union, or that was counted at another k. The grouping is recorded in `WORK/species-groups.json`, `summary.json` and the manifest's `species_pool`, and K's report header names it.
2. **Data roots.** One per distinct (cohort, k): the default root, one per cap or selection arm, and one per k arm. Each links `RAW/bvbrc`, `RAW/fasta` and `RAW/catalog`, then runs `grm curate --cap --min-class --selection [--selection-seed --stratum]`, `grm count --k`, `grm structure` and `grm matrix`. Blocks and budget arms reuse the default root.
3. **Arms.** `grm arm` (`src/ablate/arm.rs`, `run_arm`) scores one arm of one dataset.
   - It admits the arm's blocks through H's adapters: `KmerSource`, `ProbeSource::from_ksets`, `PanelSource::new` over the counted gene calls, and `PresenceSource::new` over the ksets.
   - For each split (random, lineage), repeat and outer fold, it runs `fit_encoder` and `transform` on the fold's training samples, `block_candidates` there, `assemble` on the same mask, and every learner through `learn::planned::evaluate_plan` with the plan's inner folds and holdouts. Each arm is scored once per outer fold.
   - The default arm runs first and writes the fold plans the others read (`plans/<split>-<repeat>.json`: sample ids, labels, the species grouping's clusters and duplicate groups, and the keyed `FoldPlan`).
4. **Rows.** `grm arm` writes contract §4 rows (`predictions.jsonl`): `fold` from the plan, `called` = score ≥ 0.5, and a no-call row scoring its outer fold's training prevalence with `called` false. It also writes `runs.json` (K's run record without its id, plus a cohort arm's drop counts) and `arm.json` (the arm's k, sample count, support floor, eval config, selection record, block specs and unavailable blocks).
5. **J.** For each dataset, arm, learner and split, `stats.pair_arms(..., pairing=rule)` and `stats.compare(family="C", endpoints=["log_loss"])`, then one `stats.finalize` (BH over all planned family C contrasts). A pairing J refuses, or an arm with no result, is not produced, with its reason.
6. **K.** `run.py` builds the run manifest (each paired run with `pairing`: rule, reference arm `default`, and its own plan hashes) and calls `report.build_report` with the prediction root, which writes `WORK/report/ablation-report.json` and `.md`. `WORK/summary.json` lists each arm's decision per learner and split.

## Groups for the group draw

`grm curate --selection group-round-robin-v1 --selection-seed 1` groups each eligible genome by the `mlst` field of its BV-BRC genome record (`bvbrc/genomes.jsonl`) as `mlst:<sequence type>`. A genome without one is its own singleton group, `singleton:<genome id>`, so no genome is refused. The `Selection` record names the source (`group_source` `bvbrc-genome-mlst`, `curate::MLST_GROUP_SOURCE`) and counts the species' singleton genomes (`singletons`). The `--groups` TSV of the first build is gone.

## Block screens (ruling 4)

H's K block declares a screen cap equal to the default budget, and `docs/ablation-runner.md` lists `screen_block` once per fit. `grm arm` applies no separate block screen. The learners' own chi-square screen, at the arm's `screened_features` on each fit's own fitting samples, is the one label-dependent screen, over the assembled arm matrix. Main ruled that this satisfies contract §16, since labels choose columns only on the rows of the fit that uses them. It also keeps a budget arm a one-factor change; a separate block screen would change two settings at once.

## The synthetic end-to-end test

`scripts/input-comparison/synthetic_root.py` writes a stdlib-generated data root with nothing that reads the network:
- 3,000 Salmonella enterica genomes in 40 clades, some of them exact duplicates;
- one meropenem dataset, a catalog and empty BV-BRC feature caches;
- an MLST sequence type per clade, missing on about 10% of genomes;
- the planted effect: resistance is carriage of any of 300 unrelated blaKPC alleles. Each allele is too rare for K's per-fit support floor, while H's carbapenemase panels and P's group presence call the family;
- the null effect C sees: every genome carries the gyrA sequence of H's `salmonella_qrdr` assay, and 30% carry a QRDR substitution drawn independently of the phenotype;
- 2% of labels flipped.

`GRM_INPUT_E2E_BINARY=target/release/grm python3 -B -m unittest discover -s tests/input_comparison -t .` generates the root, runs `run.py` on `spec.json` and requires:
- all 120 planned contrasts produced (11 arms and the sentinel);
- the lasso's decision `improves` on the lineage split for `k-h`, `k-p` and `k-c-h-p`, and every random-split contrast descriptive (family D, an estimate but no p);
- the sentinel present in the manifest, and a learner flagged unreliable exactly when the sentinel's lineage-split contrast for it isn't `inconclusive`;
- no `improves` anywhere for `k-c`;
- every contrast naming its arm's rule;
- the report reloading against its prediction files.

Without the variable only the spec tests run.

On the rerun after the rulings (Boat xlarge, 5 Nov 2026, seed 11) the test passed 7/7. All 110 contrasts were produced; `K+H`, `K+P` and `K+C+H+P` were `improves` for every learner on both splits, and `K+C` had no `improves`.

## The null study for the cohort arms

`scripts/input-comparison/null_study.py` runs the whole decided comparison once per cohort seed and pools the decisions. Each seed's data root has the planted H and P effect and no k-mer signal. Over 20 seeds (101 to 120, Boat xlarge, 5 Nov 2026), each arm had 200 family C contrasts (5 learners × 2 splits × 20). Shares of `improves`, with Clopper-Pearson 95% intervals:

| Arm | Pairing | `improves` | `worse` |
|---|---|---|---|
| cap-1000 | balanced | 13/200 = 6.5% [3.5, 10.9] | 133 |
| cap-4000 | balanced | 45/200 = 22.5% [16.9, 28.9] | 40 |
| grr-v1 | balanced | 35/200 = 17.5% [12.5, 23.5] | 55 |
| cohort arms pooled | balanced | 93/600 = 15.5% [12.7, 18.6] | 228 |
| k21, k25, budget-1000, budget-16000 pooled | same samples and plans | 141/800 = 17.6% | |
| K+C | same_plan | 2/200 = 1.0% | |
| K+H, K+P, K+C+H+P | same_plan | 200/200 each | 0 |

Main's ruling on this study (5 Nov 2026) withdrew the first criterion: with no k-mer signal, two different models still differ in expected log loss, so "`improves` at no more than 5%" can't hold for an arm that changes the model. The table is now context (ruling (a)): the cohort arms' no-signal `improves` rate, 15.5%, against the same-sample arms' 17.6%. The decision procedure's calibration is the next section.
- **Same-sample arms do it too.** Arms with the default's exact samples and fold plans (k and budget) score `improves` at 17.6%, more often than the cohort arms. With no k-mer signal, two different models fitted to noise really do differ in log loss (how far they overfit, how well they calibrate), and the paired test finds those differences. "No k-mer signal" doesn't make log loss equal, so no pairing rule can bring a model-changing arm to the nominal rate here. Only K+C, which adds at most two columns and leaves the model nearly unchanged, is near it.
- **The balance artefact is gone.** In the first build, `cap-1000` was `improves` for gbt, lasso and rf from the 65%-resistant test set. Under the balanced rule those learners are `worse` (rf 40 of 40), as expected for an arm trained on half the genomes at another prevalence and scored at the default's balance.
- **Contrasts within a seed are correlated** (one cohort, one test set), so the per-contrast intervals are optimistic; every seed had at least one `improves` for cap-4000 and grr-v1.

Receipts: `logs/null-study-20seeds.json` and each seed's `summary.json` and `contrasts.json` in the lane folder.

## Calibration of the decision procedure (ruling (b))

`scripts/input-comparison/calibration.json` defines two arms that are exchangeable with the default, so every `improves` or `worse` on them is a false positive:
- **(b1) `cohort-path`**: the default's inputs sent through the cohort-arm path. Its cohort is the default's, under `common_outer_test_set_balanced`, with its own subsample seed (`grm arm --balance-seed 7`).
- **(b2) `seed-twin`**: a `same_plan` twin of the default whose forest and boosting draw from another seed (`grm arm --learner-seed 1001`, `PlanView::learner_seed`). Every partition is unchanged.

`null_study.py run --spec calibration.json` ran both on seeds 101 to 120, fixed in advance, with each seed's own BH family of 20 contrasts. Rates are cluster-robust with the seed as the cluster: the mean of per-seed shares, with a 95% Student t interval over the 20 independent seeds.

| Arm | `improves` | `worse` | Contrasts |
|---|---|---|---|
| cohort-path (b1) | 0.0% [0.0, 0.0] | 0.0% [0.0, 0.0] | 0/200 and 0/200 |
| seed-twin (b2) | 2.0% [0.1, 3.9] | 3.0% [0.8, 5.2] | 4/200 and 6/200 |

Both are at or below 5%, so the procedure passes as ruled. Two caveats belong with that:
- **(b1) is exact by construction.** When the cohort equals the default's, the balanced subsample keeps every common genome whatever its seed, and the inner draws are the default's. So the arm reproduces the default bit for bit (Rust test `the_calibration_arms_are_exchangeable_with_the_default`), and its contrasts are exactly zero. It shows the cohort path adds nothing of its own. It doesn't test the paired test's calibration.
- **(b2)'s per-arm rate is diluted.** scm, cart and the lasso draw nothing random, so 6 of the twin's 10 contrasts per seed are exactly zero. On the four that can differ:
  - rf: 0 of 40;
  - gbt on the lineage split: 1 of 20 `improves`;
  - **gbt on the random split: 3 of 20 `improves` and 6 of 20 `worse`**, 45% against 5%.

  So for boosting on the random split, J's test is anti-conservative. The likely cause is fold dependence. On the random split every lineage cluster spans all five outer folds, and the twin-minus-default difference has a per-fold part (each fold's boosted model stops at its own tree count) that every cluster in that fold shares. The cluster jackknife and the sign flip treat clusters as independent. The lineage split nests clusters in folds, and rf averages its seed away over many trees, which fits both being near nominal.

  The fix would change J's inference unit for the random split: resample outer folds, or fold-by-cluster blocks. With five folds, the exact sign flip can't reach p < 0.05 (smallest p 2/32). That is a method change, and main hasn't decided it, so J is unchanged here.

Receipts: `logs/calibration-20seeds.json` and each seed's `summary.json` and `contrasts.json` (`logs/cal1/cal`, `logs/cal2/cal`).

## Which split decides, and the sentinel (main's ruling after the calibration study)

The calibration above (seed twin over 20 seeds: 2.0% `improves` and 3.0% `worse` overall; gbt on the random split non-inconclusive in 9 of 20; gbt on the lineage split 1 of 20; rf 0 of 40) led main to three rules, which the runner and K now apply.

1. **The random split is descriptive only.** `run.py` plans every random-split contrast in family D, the contract's family with no p and no direction. J still computes the estimate, its interval and the AUROC delta. K prints these contrasts in their own block, "Random split: descriptive only", with the reason: with five outer folds a fold-aware test can't reach p < 0.05, and the cluster-level test is anti-conservative on this split. Using family D, rather than a new `descriptive` decision value, keeps K's closed decision set, J's `decide` and the BH families unchanged. Family D's fixed `inconclusive` is never printed as a test result: the block says it decides nothing, and summary.json shows these rows as `descriptive`.
2. **Decisions come from the lineage split only.** It is the only family C split, so BH adjusts over the lineage contrasts alone.
3. **A sentinel in every comparison.**
   - `run.py` adds `sentinel-seed-twin` to every spec run, once per species. It is the default's learner-seed twin (learner seed 1001, `same_plan`), exchangeable with the default, and the manifest names it in `sentinel_arms`.
   - K reads its lineage-split log-loss contrast for each learner. If that contrast is anything but `inconclusive`, or was not produced, every contrast of that learner on that species is marked unreliable. The report's `unreliable` section lists the species, learner, sentinel contrast, reason and flagged contrast ids, and each flagged line ends "Unreliable (rule 10)."
   - The records are kept as J wrote them. summary.json shows them as `unreliable` (J's value is in `raw_decisions`) and leaves them out of `any_improves` and `any_worse`.
   - K's reading rules 9 and 10 state this and the calibration numbers in every report, H-1's included. K's tests force a sentinel decision on the fixture and check the flags.

## After REVIEW-INPUT (main's rulings, 5 Nov 2026)

- **Canonical genome order (M1).** The review found that genome order alone moved decisions. An arm with the default's genome set in another order came out non-inconclusive on 10 of 50 lineage contrasts, scm, cart and the lasso included. So `run.py` rewrites every data root's genomes.txt sorted by genome id right after `grm curate`, before count, structure, matrix, fold plans and evaluation. The default arm is included.
  - Two arms with the same genome set therefore give bit-identical predictions.
  - The drawn order's sha256 is kept in the root's `curated.json`, and a group draw's order stays in `selection.json`.
  - `grm learn`'s production path is unchanged (no revision change). The default-arm identity check is now: the default arm equals `grm learn` run on the same genomes in canonical order. The e2e test `test_the_default_arm_is_grm_learn_on_the_canonical_order` checks it bit for bit, and `test_a_shuffled_genome_list_gives_byte_identical_arm_output` shows that a shuffled genomes.txt gives byte-identical predictions and plans.
- **Identical cohorts.** A cohort arm (cap or selection) whose genome set equals the default's is not run. Its runs are `identical_cohort`, with the reason; K lists them under "Identical cohorts, not run" (reading rule 11), and every planned contrast of the arm is not produced with that reason. `--run-identical-cohorts` runs such an arm anyway, for the order-twin proof only.
- **H-1 panel floor.** `spec.json`'s `panel_floor` requires more than 2,000 eligible genomes per species (`min_eligible_genomes` 2001), so that cap-4000 and grr-v1 change the genome set. `run.py` refuses a species below it unless `--below-panel-floor` is given (synthetic studies only). It counts eligible genomes as curation does: an unconflicted reported R or S call on one of the species' dataset drugs.
- **Duplicates on every split (L1).** `cohort_plan` (format `grm-cohort-plan-v3`) holds out of fold f's training every genome in the k-mer duplicate group of a genome scored in f, on both splits, and records the count (`duplicate_dropped`, in each run's `cohort_drops`). Test: `a_cohort_arm_never_trains_on_a_copy_of_a_scored_genome`.
- **Where block C exists (L3).** H's probe assays are defined for M. tuberculosis, N. gonorrhoeae, E. coli and Salmonella enterica only (`src/probes.rs`), and only two of them can call today:
  - **Callable:** N. gonorrhoeae (`gc_gyra`) and Salmonella (`salmonella_qrdr`).
  - **Not callable:** the MTB assays have no verified wild-type sequence, and E. coli's `ecoli_qrdr` is marked unsupported (synonymous changes read as mutations).
  - **Effect on H-1:** the panel prefers species with block C. On any other species, K+C and K+C+H+P are reported as not produced, with C's reason, because an arm with an unavailable block makes every one of its runs unavailable.
- **The sentinel's reach (I1).** The learner-seed twin probes only the learners that draw randomness, rf and gbt. scm, cart and the lasso draw nothing, so their twin contrast is zero by construction. Canonical order removes genome-order effects for them. K's reading rule 10 says so.
- **MLST groups (I3).** The group draw takes a genome record's BV-BRC `mlst` field whole, after trimming, as its group (`mlst:<value>`). A record that joins more than one scheme's sequence types (for E. coli, Achtman and Pasteur) therefore forms a group per combination, which is finer than one sequence type. A genome without the field is its own singleton.
- **Reading rules (L2).** K's rules are numbered 1 to 11 in order, and a test renders the Markdown and checks that every "rule N" points at the rule that says it.
- **Mutants (L4).** New tests kill the reviewer's surviving mutants:
  - M5 and M6: `the_fold_matrix_reads_only_the_training_samples`, on `ablate::arm::fold_matrix`;
  - M17: `a_cohort_plan_keeps_its_held_out_genomes_in_run_arm`, on `fold_view`;
  - M14: `test_summary_keeps_unreliable_and_descriptive_rows_out_of_its_decisions`.

## After REVIEW-INPUT2 (main's rulings, 5 Nov 2026)

- **Folds keyed on genome identity (N1).** Canonical order made two equal genome sets index-identical, but one extra genome shifted every later index, so a near-identical arm's inner folds and holdouts were re-drawn. Five extra genomes changed 84% of prediction rows. Inside the comparison every fold assignment is now keyed on genome identity:
  - the outer folds of the default arm's plan (`plan::keyed_fold_plan`, format `grm-fold-plan-keyed-v1`), and every arm's inner folds and tuning holdouts (`keyed_tuning`);
  - the cohort plan's class-balanced subsample (`cohort_plan`, format `grm-cohort-plan-v4`).
- **How a unit is keyed.** A draw unit is a lineage on the lineage split, keyed by its lowest genome id; on the random split it is a k-mer duplicate group, keyed by its lowest genome id, or a single genome by its id. Each unit's assignment is `sha256("grm-keyed:{purpose}:{seed}:{key}")`'s first eight bytes, taken modulo:
  - **outer fold:** `% k`, with the repeat seed;
  - **inner fold:** `% inner_k`;
  - **boosting holdout:** `% 7 == 0`;
  - **lasso holdout:** `% 5 == 0`.

  The last three use outer fold f's seed (`repeat_seed + 101 f`). Tests: `keyed_folds_do_not_depend_on_the_other_genomes` (Rust) and `test_the_default_arm_folds_match_an_independent_id_keyed_reference` (Python, from the plan file's ids and groups).

  REVIEW-INPUT3 found that this keying was not membership-free: a cohort arm keyed a unit by its own lowest id, so dropping or adding that one genome re-keyed the whole unit, and each arm's root clustered its own lineages. The rule as it now stands is in the next section.
- **Balance and refusal.** Stratification gives way to approximate balance. A keyed plan is refused, with its reason, unless:
  - the lineage split has at least 2k lineages;
  - every outer fold's test and training samples hold both classes;
  - every inner fold's training and validation samples hold both classes;
  - both sides of both holdouts hold both classes.

  Tests: `a_keyed_plan_is_refused_when_a_fold_lacks_a_class`, and since REVIEW-INPUT3 one per branch (next section).
- **The comparison's protocol differs from production.** `grm learn` and `grm pipeline` keep G's index-drawn, stratified folds over genomes.txt order, unchanged (no revision change). The comparison differs in two ways: its fold draw is keyed on genome id, and its genome order is canonical. So the comparison's default arm is not `grm learn`'s evaluation. It is checked against the id-keyed reference above, and its scores are not production scores.
- **What remains index-based:** the forest's and boosting's own draws inside a fit, such as bootstrap samples and feature subsets. These are what the learner-seed sentinel probes.
- **The near-identical twin, before and after (N1's proof).** 2,005 genomes; cap-4000 adds 5 genomes and grr-v1 swaps 5; seeds 201–205, shipped flags.
  - **Before (a2ebe1b9):** shared genomes kept their outer folds but only 57% (cap-4000) and 38% (grr-v1) of their inner folds and 84% and 77% of their holdout memberships. 84% and 83% of prediction rows changed, and 13 of 50 lineage contrasts were non-inconclusive.
  - **After (keyed):** shared genomes keep 100% of their outer folds, inner folds and holdout memberships. Prediction rows still change (79% and 78%), because 5 more training genomes change what every learner fits. 7 of 50 lineage contrasts are non-inconclusive (14.0%; per-seed 0.2, 0, 0, 0.2, 0.3; 95% t interval over seeds −2.7 to 30.7, clamped at 0).
  - Receipts: `logs/neartwin` in the lane folder.
- **Base rate (N1).** The near-identical twin's no-signal decision rate was printed beside every cohort arm's contrasts. REVIEW-INPUT3 I-c replaced it with a rate per arm kind (next section).
- **M17 at `run_arm` (L4).** `grm arm --trace-training` writes `training.jsonl`: the arm's sample ids, then each outer fold's training ids exactly as `run_arm` handed them to the learners. The e2e test `test_run_arm_never_trains_a_cohort_arm_on_its_held_out_or_scored_genomes` checks them against the cohort plans `run_arm` wrote.
- **Rule 10 (L-a)** now says that canonical order removes genome-order effects only between arms with equal genome sets, and that keyed folds keep a shared genome's folds. Neither makes a cohort arm's decision free of the change in its training samples.
- **The 20.7% to 30.0% change (L-b), corrected.** The first build's 31 of 150 and the review fixes' 45 of 150 lineage `improves` (seeds 101–110) were attributed above to identical-set arms no longer running and to BH covering lineage contrasts only. REVIEW-INPUT2 shows that neither explains the change:
  - no seed had an identical-set arm;
  - re-adjusting the old p values by BH over lineage contrasts alone gives exactly the old counts;
  - the old family was 110 contrasts per seed, not 120.

  Every one of the 150 contrasts changed its estimate, because canonical order re-drew the default's partitions too. So the two figures are two draws of one noisy base rate (per-seed SD about 12–13 points), not an effect of either rule.

## After REVIEW-INPUT3 (main's rulings, 5 Nov 2026)

- **One grouping, keys through the default (N2).** A unit's fold key is a function of the unit's identity and the default cohort only:
  - **Grouping.** Every arm's lineages and k-mer duplicate groups are the species pool's (`plan::SpeciesGroups`, read from the pool's `structure.json` by `grm arm --groups`), restricted to the arm's genomes. Before this, each arm's root ran its own `grm structure`, whose average-linkage lineages and threshold depend on the cohort. So two arms now never disagree on whether two genomes share a lineage or are copies. A cohort arm needs `--groups`, and it refuses a reference plan that wasn't split on the same grouping.
  - **Keys.** A lineage or duplicate group is keyed by its lowest-id member in the default cohort. A unit with no member in the default is keyed by its lowest-id member in the species pool, fixed once per species before any arm. For the default arm this is its own lowest id, as before (`unit_keys`); a cohort arm keys through the reference (`SpeciesGroups::unit_keys`).
  - **What holds.** A genome both the default and a cohort arm hold is scored only in its default outer fold (or not at all, when the balanced subsample leaves it out), and in every outer fold where it trains it has the default's inner fold and gbt and lasso holdout memberships, whichever genomes the arm drops or adds. Its scores still change, because its training neighbours do.
  - **What doesn't.** Two default arms on different roots (two separate comparisons) each key on their own default cohort, and two pools give two groupings. The rule is per comparison.
  - **Tests.** `a_cohort_arm_keeps_shared_folds_when_it_drops_or_adds_a_unit_minimum` (Rust: drop five lineage minima and a duplicate group's lower id, add genomes below ten lineages' lowest ids and below a group's, three seeds, both splits; it asserts that keying on the arm's own ids would re-key those units), `the_species_grouping_refuses_a_genome_outside_the_pool`, and in the e2e `test_every_cohort_arm_keeps_the_default_folds_of_its_shared_genomes` and `test_every_arm_splits_on_one_species_grouping`.
  - **Cost.** Superseded by the union pool (next section).
- **Refusal tests (L2, L3).** One test per branch: `a_keyed_lineage_plan_with_fewer_than_twice_k_lineages_is_refused` (K3), `keyed_tuning_refuses_an_inner_fold_or_holdout_side_without_a_class` (inner fold, gbt and lasso holdout sides; K4), and `the_balanced_subsample_keeps_its_shared_genomes_when_the_arm_grows` (K5: 30 arm-only genomes sorting before the shared ones leave the scored set unchanged).
- **Refused plans (I-a).** A keyed holdout or fold lacking a class still refuses the plan. `summary.json` lists every refused plan with its reason (`refused_plans`).
- **The re-keying twin (N2's proof and I-c's rates).** `tools/boat/rekeytwin.py` in the lane folder, seeds 301–310 fixed in advance, Boat xlarge. Per seed: a synthetic pool of 2,020 genomes; a default cohort of the pool minus 5 lineage minima and 10 other genomes drawn across the id range; three twins of the default, each built, counted and scored through `grm arm` against it with `--groups` the pool:
  - **cap-down:** drops 5 genomes, 2 of them their lineage's lowest id in the default;
  - **cap-up:** adds 5 of the left-out genomes, 2 of them below their lineage's lowest default id;
  - **selection-swap:** both.

  | Kind | Shared genomes the twin's own lowest ids would re-key (lineage, 10 seeds) | Folds kept (outer, inner, holdouts; both splits) | Lineage contrasts non-inconclusive, per seed | Rate (95% t interval over 10 seeds) |
  |---|---|---|---|---|
  | cap-down | 964 | 100% (20,000 / 80,000 / 160,000 per split) | 0, .2, 0, .2, .4, 0, 0, 0, .6, 0 | 14.0% [0, 29.2] (lower bound clamped from −1.2) |
  | cap-up | 1,364 | 100% (20,000 / 80,200 / 160,400) | 0, .2, 0, .4, 0, .2, .2, 0, .4, .2 | 16.0% [4.7, 27.3] |
  | selection-swap | 2,434 | 100% (20,000 / 80,000 / 160,000) | .2, 0, .6, .6, 0, 0, .2, 0, .2, 0 | 18.0% [0.9, 35.1] |

  No default or twin plan was refused (0 of 80). Each rate was measured on 5-genome changes. A real cap-1000 arm changes about 1,000 genomes; how the rate of a larger cohort change relates to these is unmeasured, in either direction. Each rate's source string says so, and K prints it. Receipts: `logs/twin8` in the lane folder.
- **Base rates per kind (I-c).** See the twin study above; `spec.json` `base_rates` holds one record per kind (`cap-down`, `cap-up`, `selection-swap`). `run.py` gives each cohort arm its kind (a cap below the default's, above it, or another selection rule) and K prints that kind's rate beside the arm's contrasts, or "not measured" when the spec has none.

## After REVIEW-INPUT4 (main's rulings, 5 Nov 2026)

Main revised the "full eligible pool" ruling: that pool can't be built on real data. The master root has FASTAs only for selected genomes, and Salmonella's 14,521 eligible genomes need about 390 GiB of ksets against 260 GiB free on an xlarge.

- **The pool is the union of the arms' cohorts (M1).** Every arm's root is curated first. The pool is the union of their curated genomes.txt, counted and structured once, before any `grm arm`. A unit with no member in the default is keyed by its lowest-id member in that union.
  - **Refusals.** A genome of the union without a FASTA refuses the run, naming it; it is never dropped. The pool's count and structure read only FASTAs, not annotations. After `grm structure`, the pool must hold exactly the union. Before this, `grm count` skipped a genome without a FASTA and the run went on with a smaller pool.
  - **Disk.** Superseded by the whole-comparison check (next section).
  - **One union per comparison run.** The union, and with it the lineages and the keys of units outside the default, depends on which arms the run holds. Adding an arm means rerunning the whole comparison. K's report header says so ("Species pool: the union of the cohorts of this run's N arms …").
- **The default's lineage split (I3).** Since REVIEW-INPUT3, the default arm's lineages come from clustering the pool, not its own cohort. `choose_lineages` depends on the genomes and their tree, so on real data the default's lineage split can differ from the one its own `grm structure` would give. On the synthetic seeds the two partitions agreed.
- **k (I5).** `grm arm` refuses `--groups` counted at another k than the arm's, and `run.py` passes `--groups` only to the default and the cohort arms. The k arms take the default's k = 31 plans (`k31_groups`) and read no grouping. Test: `test_an_arm_refuses_a_grouping_counted_at_another_k`.
- **Refused-plan listing (L2).** `test_refused_plans_lists_an_unavailable_plan_with_its_reason` gives `refused_plans` a plan folder with a refused fold plan and a refused cohort plan.
- **Tests (M1).** `test_the_pool_refuses_a_genome_of_the_union_without_a_fasta`, `test_the_pool_refuses_a_structure_that_is_not_the_union`, `test_the_pool_s_kset_estimate_refuses_more_than_its_disk_share`, and in the e2e `test_every_arm_splits_on_one_species_grouping` (pool = union of every root's genomes; `--groups` on exactly the default and cohort arms; the manifest's `species_pool` and K's header line).

## After REVIEW-INPUT5 (main's rulings, 5 Nov 2026)

- **The disk check sizes the whole comparison (M1(a)).** Every data root is curated first, so the genome count of each is known before any count. The estimate is:
  - the pool's ksets: the union's genome count times the per-genome size;
  - plus, for each distinct data root the run counts (`counted_roots`: the default's, and each cap, selection and k arm's; an identical cohort isn't built; blocks, budget and calibration arms reuse the default's), its genome count times the per-genome size, plus `pool_budget.root_overhead_share` for its matrix and sidecars (5% in round 10, 10% since REVIEW-INPUT6).
  - The per-genome size is the mean of the species' own counted ksets in RAW (`species/<slug>/ksets`), or the spec's stated default (`pool_budget.default_kset_bytes`, 27.5 MiB, a Salmonella k31 kset). A k21 or k25 root is sized as a k31 one, which is an upper bound: on the synthetic e2e their ksets are 0.63 and 0.80 times a k31 kset.
  - It refuses above 80% (`pool_budget.max_disk_share`) of the free disk, before any count.
- **Where it runs.**
  - **Before any box:** `run.py SPEC --raw RAW --work WORK --grm GRM --source-commit SHA --estimate-only --free-disk-bytes N` curates every arm, prints the estimate (pool, each root, total) against N, the box's free disk, and exits 1 when it doesn't fit. With `--free-disk-bytes`, this machine's disk is not read.
  - **On the box:** the same check runs before the pool is counted. It takes the smaller of `--free-disk-bytes` (when given) and WORK's measured free space (I2), and writes `WORK/disk-estimate.json`.
- **Ksets are counted once per root (M1(b), as measured).** `grm count` reuses a kset only when the same root's own kset folder already holds it for the same genome and k with a matching sidecar (`species::count_one`). Every comparison root has its own kset folder, so a genome in the pool and in three arm roots is counted four times, and the estimate counts it four times. No counting behaviour changed. Test: `test_every_root_counts_its_own_ksets` (e2e: the pool's and the default's kset of one genome are two files with the same bytes).
- **What it means for the spec's species.** For Salmonella at 27.5 MiB per kset, the spec's roots are 13,000 genome-ksets: default 2,000, cap-1000, cap-4000 4,000, grr-v1 2,000, k21 2,000 and k25 2,000. They come to 366.6 GiB with the margin. With the pool's 4,000–6,000 genomes (107.4–161.1 GiB), the total is 474.0–527.7 GiB, against 207.4 GiB (80% of an xlarge's 259.3 GiB). So the comparison refuses on one xlarge. Splitting a comparison across boxes is not this round; H-1's sizing plans it. Test: `test_a_salmonella_sized_comparison_refuses_on_one_xlarge`.
- **The e2e's kset size.** At the spec's 27.5 MiB, the e2e's own comparison (a 3,000-genome pool and 12,000 genome-ksets in its six roots) needs about 419 GiB, so it refuses on an xlarge before any count. On the first round-10 gates run it did exactly that. The e2e runs a copy of the spec stating 1 MiB per kset (`spec-e2e.json` in its output folder), an upper bound for its 23 kb genomes.
- **Measured against the estimate.** `test_the_disk_estimate_covers_every_root_the_run_counted` takes the e2e's own per-genome size (the mean of its pool's ksets). It requires the estimate to name exactly the roots the run counted, to be at least the disk the run used, and to be within 1.5 times it. It writes `WORK/disk-measured.json`. The lane report gives the figures.
- **Tests (L1, L2).** `test_estimate_only_sizes_the_box_it_is_given_and_exits_1_when_it_does_not_fit` and `test_the_pool_is_the_arms_union_not_the_eligible_set` run `run.main` with a stub `grm` that only curates (50 eligible genomes, the arms selecting 30).
- **Reused WORK (I5).** `refused_plans` reads only this run's arms. A folder kept for an arm no longer in the spec is named on stderr and in `summary.json` (`stale_arm_folders`), and not read.
- **Base-rate sources (I6).** The rates were measured with a full 2,020-genome pool. REVIEW-INPUT5 reran them on the union pool and got identical rates (30 of 30 seeds), and the source strings now say so.

## After REVIEW-INPUT6 (main's rulings, 6 Nov 2026)

- **Every run checks disk before counting anything (M1).** Before this, the on-box check ran only while the pool's `structure.json` was missing. So a rerun after a failure, or an arm added to a reused WORK, counted new roots unchecked: the reviewer's k27 arm at 10^15 B per kset counted 2,000 ksets. Now `check_disk` runs on every run, after curation and before the pool or any root is counted. It estimates only what is still to be written:
  - **Pool:** its genomes without a finished kset (the kset and its sidecar both written); none once its `structure.json` exists.
  - **Each data root:** its genomes without a finished kset; plus, unless its `built.json` exists, its matrix and sidecars at `root_overhead_share` (10%) of all its genomes' ksets.
  - **Arm outputs** (`WORK/arms`, rewritten on every run): `arm_output_share` (5%) of every root's ksets.
  - **Fresh WORK:** everything is pending, so the estimate is the full one.

  It refuses above 80% of the free disk measured at that moment: the smaller of WORK's volume and `--free-disk-bytes` when given. It writes `WORK/disk-estimate.json` with each root's `pending` and `built`. Test: `test_a_reused_work_checks_disk_before_counting_a_new_root`. That e2e test reproduces the reviewer's case on a link copy of the finished e2e WORK: a k27 arm at 10^15 B per kset refuses, with the k27 root 2,000 pending, every other root and the pool 0 pending, and no k27 kset written.
- **Margins (L1).** The matrix margin is 10% of a root's ksets (synthetic roots measured 4.0–5.9%; the MTB 4k matrix was about 6% of its ksets). The arm outputs are 5% of the roots' ksets (synthetic: 3%). Both, like the 27.5 MiB kset, are stated figures until H-1 measures real per-species ksets and matrices, and `pool_budget.reason` says so. The 80% share is headroom above them, not the margin. Code defaults apply only to a spec without `pool_budget`; both shipped specs state their own.
- **The box's free disk (L2).** `test_the_box_takes_the_smaller_of_the_given_and_the_measured_free_disk` gives `--free-disk-bytes` 10^15 while WORK's volume reads 1,000 B free; the check refuses on the 1,000.
- **Stale state (L3).**
  - `run.py` deletes an arm's `out/plans` before it runs that arm, or records it as an identical cohort. Test: `test_an_arm_s_plans_from_an_earlier_run_are_cleared_before_it_runs`.
  - WORK records its inputs in `WORK/inputs.json`: the sha256 of RAW's input manifest (every file under `bvbrc/`, `fasta/` and `catalog/` with its sha256 and size), the file count, and the grm binary's sha256. A reused WORK whose RAW or grm changed is refused, naming which; so is a WORK with data roots but no record. Hashing RAW reads every input file once per run. Test: `test_a_reused_work_refuses_another_raw_or_another_grm`.
- **Measured against the estimate, whole WORK.** `test_the_disk_estimate_covers_every_root_the_run_counted` now measures the whole WORK: roots, pool, arm outputs, logs, `j/`, `k/` and the report. It compares that with pool + roots × 1.10 + roots × 0.05 at the e2e's own mean kset size, and requires measured ≤ estimate ≤ 1.5 × measured. The lane report gives the figures.

## After REVIEW-INPUT7 (main's rulings, 6 Nov 2026)

- **Roots curate as revision 12's S0 did (M1).** Each data root, the pool included, links RAW's `species-check/` beside `bvbrc/`, `fasta/` and `catalog/`, so `grm curate` drops the genomes the species check fails. `grm curate` also gets `--sigdrug-vocabulary`: `run.py --sigdrug-vocabulary DIR` names GRM's own vocabulary copy, the one S0 used (a copy of `grm-r9-curate/vocabulary`; S0 `REPORT.md`, section 2).
  - **Before this,** the roots curated without either. On the revision-12 master root they would have kept genomes revision 12 drops as another species: 69 K. pneumoniae, 4 E. coli, 2 MTB and 1 Salmonella.
  - **The spec requires both** (`curate_inputs`). `run.py` refuses a RAW without `species-check/species_check.json`, and a run without `--sigdrug-vocabulary`.
  - **The revision-12 gate.** Superseded in REVIEW-INPUT8 (next section): the gate compares the default's drawn genome order and its `datasets.json` without `selection`.
  - **Synthetic RAW.** `synthetic_root.py` now writes `species-check/species_check.json`, which fails the lowest-id resistant genome (one the default would select). It also writes a SigDrug-format `vocabulary/` whose brand fold turns `synthmerrem` (the name every 150th genome's record uses) into meropenem. `truth.json` names both. GRM's own resolver can't reach that name (its typo match stops at edit distance 2), so only the vocabulary labels those genomes. A typo such as `meropenam` would not do: the first box run showed GRM resolves it without the vocabulary.
  - **Tests.** `test_the_default_curates_as_revision_12_did` (e2e) checks three things: the default drops the failed genome; every misspelled genome is labelled meropenem; and the default's `datasets.json` and genome set equal a plain `grm curate --cap 2000 --min-class 25 --sigdrug-vocabulary` over a root laid out as S0's. Unit tests: `test_the_revision12_gate_refuses_another_default_cohort` and `test_curate_inputs_are_required_when_the_spec_asks`.
- **Check placement (L1).** `test_main_checks_disk_before_it_counts_anything` runs `run.main` past `--estimate-only`, with a stub grm that records and fails every count and a free disk that refuses. The refusal comes first and nothing is counted. That test kills a check moved after the pool, moved after the default root, or removed.
- **Names and folders (L2).**
  - Arm names and dataset ids must match `^[a-z0-9][a-z0-9._-]*$`, and a dataset listed twice is refused.
  - `arm()` refuses an output folder that doesn't resolve inside `WORK/arms`, and a `plans` folder that is a symlink. It clears `plans` with a plain `shutil.rmtree`, so a failure stops the run.
  - **The four block arms are renamed** `k-c`, `k-h`, `k-p` and `k-c-h-p` (they were `K+C`, `K+H`, `K+P`, `K+C+H+P`), since `+` and capitals fall outside that pattern. Their blocks are unchanged.
  - Tests: `test_arm_and_dataset_names_are_plain_and_unique` (`../../escape` among them), `test_an_arm_folder_outside_work_or_a_linked_plans_folder_is_refused` and `test_an_identical_cohort_s_old_plans_are_cleared_too`.
- **Inputs (L3).**
  - **What `inputs.json` hashes:** every file under RAW's `bvbrc/`, `fasta/`, `catalog/` and `species-check/`, and each arm's stratum file, read in 1 MiB blocks; the SigDrug vocabulary folder; and the grm binary.
  - **How RAW is walked:** without following links. A symlink inside the hashed parts is refused (since REVIEW-INPUT8, a symlink anywhere in RAW).
  - **`--grm`** is resolved with the same PATH lookup `subprocess` uses, and that file is hashed and run.
  - **Tests:** `test_the_input_record_covers_non_fasta_inputs_and_refuses_links` changes an annotation, the species check and the vocabulary in turn, and adds a symlink; `test_a_bare_grm_name_is_hashed_as_path_finds_it`. The full hash still reads the hashed parts once per run, 61.2 GB in 51,545 files on the master root (REVIEW-INPUT8 I1).
- **Disk terms (L4).**
  - The estimate now sizes the pool's copy of the default's `curated/` (its measured size in disk blocks, until the pool has it) and WORK-level files: logs, `inputs.json`, `disk-estimate.json`, `species-groups.json`, `j/`, `k/` and the report, at a stated 64 MiB (`pool_budget.work_files_bytes`).
  - The e2e now measures disk in allocated blocks, as `df` and the check count it (`run.disk_bytes`). The lane report gives the figures.

## After REVIEW-INPUT8 (main's rulings, 6 Nov 2026)

- **The revision-12 gate compares two values (M1).** Diff 01 makes `grm curate` record a `selection` on every dataset, so this tree's `datasets.json` (`00aaacf2…` on the real root) could never equal revision 12's `5519c61d…`, and H-1 would always have refused. Before anything is counted, and in `--estimate-only` too, `check_revision12` now requires both of:
  - **the default's genome order** as curate drew it (`curated.json` `drawn_order_sha256`), equal to revision 12's `species/salmonella_enterica/genomes.txt` (`5864bd5d…`). This value pins the cohort.
  - **its `datasets.json` with every `selection` record removed** (`datasets_view`: two-space indent, no trailing newline, as grm writes it), equal to `5519c61d…`. This value pins the labels.

  `--estimate-only` prints the gate's values with its estimate.
  - **Pinned through the real binary.** The e2e runs with the gate on, set to the synthetic default's own values (`SYNTHETIC_GATE` in the test file: drawn order `2299eda7…`, `datasets.json` view `fa73a622…`). A change in curate's output or its record then fails the e2e before H-1. Tests: `test_the_gate_held_the_pinned_synthetic_values`, `test_the_gate_s_datasets_view_drops_only_selection_records` and the stub gate test, which covers both keys and `--estimate-only`.
  - **On the real root:** a laptop-only, read-only `--estimate-only` passed the gate (lane report, round 13).
- **The binary K records (L2).** K's manifest hashes the binary the run resolved and ran (`Runner.grm`), the same one `inputs.json` records. The e2e passes a bare `--grm` name found on PATH, with a decoy named `grm` in run.py's working directory. Test: `test_the_manifest_hashes_the_binary_the_run_ran`.
- **The vocabulary variable (L3).** Every grm child runs without `GRM_SIGDRUG_VOCABULARY`, so a vocabulary reaches curate only through `--sigdrug-vocabulary`, which `inputs.json` hashes. A spec must state `curate_inputs`; `calibration.json` now says species check yes, vocabulary no. Tests: `test_grm_children_never_see_the_vocabulary_variable` and `test_a_spec_must_state_its_curate_inputs`.
- **Every root curates alike (L4).** `test_every_root_curates_with_the_species_check_and_the_vocabulary` (e2e) checks every curated root, cohort arms included: the `species-check` link, `--sigdrug-vocabulary` on its curate line, and no genome the check fails.
- **Links and strata (L5).** RAW is walked whole, without following links, and any symlink anywhere in it is refused. A stratum must lie inside RAW; it is hashed by content. Tests: `test_a_symlink_anywhere_in_raw_is_refused` (a linked part, a linked file, a link outside the hashed parts) and `test_a_stratum_file_is_hashed_and_must_lie_inside_raw`.
- **Names (I2).** Matched whole (`fullmatch`), at most 64 characters, so a trailing newline or a 300-character name is refused.

## After REVIEW-INPUT9 (main's rulings, 6 Nov 2026)

- **No real run goes ungated (L1).** A spec without `revision12_gate` refuses, in `--estimate-only` too, unless it says `"synthetic": true`. K's manifest records the flag (`synthetic`), and the report header says "Synthetic run: yes" or "no".
  - The e2e's spec copy and `calibration.json` say `"synthetic": true`.
  - `null_study.py` runs a synthetic copy of its spec: no gate, `"synthetic": true`, and a 1 MiB kset size. Without it, the real-root gate and the 27.5 MiB disk default would have refused every seed.
- **One genome order per species.** `revision12_gate.genome_order` maps each H-1 panel species to the sha256 of the master root's `species/<slug>/genomes.txt`. The gate compares every species the map names: the run's own species through its drawn order (`curated.json`), and the others through the default root's untouched `species/<slug>/genomes.txt`. The run's own species must be named. The `datasets.json` view still covers all 277 datasets.
- **The panel**, applying DECISIONS H-1 to revision 12's run: more than 2,000 eligible genomes, then the four with the most datasets whose pipeline stage took under one xlarge-hour (`grm-r12-run-20261105/budget.tsv`).
  - K. pneumoniae (39 datasets, R9 1.698 h) and E. coli (37, R8 1.351 h) are over the hour.
  - The panel is Salmonella (25 datasets, R10 0.856 h), MTB (22, R7 0.584 h), S. aureus (19, R3 0.959 h) and A. baumannii (18, R6 0.782 h; 2,482 eligible).

  | Species | `genome_order` (master `species/<slug>/genomes.txt`) |
  |---|---|
  | salmonella_enterica | `5864bd5d9a7cb6662127386e7218d608e43febffc45a5a0b15e4c24fbf6e4a8c` |
  | mycobacterium_tuberculosis | `b104d8bd7ffaa2d0dc3032fe30ca0c0b32670a7099a02bcd49860e5c77f59a64` (S0's MTB 2k binding) |
  | staphylococcus_aureus | `3269d7b9d182e65e77d47709d928659183380545dcbf5b9d5d42ee4ed8042c56` |
  | acinetobacter_baumannii | `f0fc6b442c5b3d139263dea41071b194f10ee38ac3d7b713d4021e79b6f379e2` (replaced; see below) |

  A species the map doesn't name refuses. Adding one is one line: its master `genomes.txt` sha256.
- **Tests (L2).**
  - `test_the_gate_s_datasets_view_drops_only_selection_records` uses a three-entry `datasets.json`.
  - `test_a_symlink_anywhere_in_raw_is_refused` links a folder outside the hashed parts.
  - `test_the_revision12_gate_refuses_another_default_cohort` uses wrong values that differ from the right ones only in the last hex digit. It also covers a species the map doesn't name, a named species whose file is missing, and a second species that differs.
  - `test_a_real_run_without_the_gate_refuses` covers `synthetic` false, absent and `"true"` (a string).

## After INPUT round 14: the panel (main's ruling, 6 Nov 2026)

- **A. baumannii is out.** Its default selects 1,702 of its 2,482 eligible genomes, so its cap-4000 and grr-v1 arms would select the default's genomes.
- **The floor now reads:** the default selection is capped at 2,000, and more than 2,000 genomes are eligible.
- **The fourth species** is the next by dataset count that meets the floor and whose revision-12 pipeline stage took under one xlarge-hour, block C species first. That is N. gonorrhoeae.
- **Candidates** (master root `curated/species.json`; stage times from `grm-r12-run-20261105/budget.tsv`, stage members from its NOTES and `logs/drive-R1.out`):

  | Species | Datasets | Eligible | Selected | Stage, xlarge-h | Floor | Under 1 h | Block C |
  |---|---:|---:|---:|---|---|---|---|
  | K. pneumoniae | 39 | 8,257 | 2,000 | R9 1.698 | yes | no | |
  | E. coli | 37 | 12,411 | 2,000 | R8 1.351 | yes | no | |
  | **Salmonella** | 25 | 14,521 | 2,000 | R10 0.856 | yes | yes | yes |
  | **MTB** | 22 | 32,058 | 2,000 | R7 0.584 | yes | yes | |
  | Shigella sonnei | 20 | 1,994 | 1,984 | R5 0.596 | no | yes | |
  | **S. aureus** | 19 | 5,374 | 2,000 | R3 0.959 | yes | yes | |
  | A. baumannii | 18 | 2,482 | 1,702 | R6 0.782 | no | yes | |
  | E. faecium | 14 | 3,138 | 2,000 | R2 0.618 | yes | yes | |
  | S. pneumoniae | 12 | 17,837 | 2,000 | R4 1.538 (with P. aeruginosa) | yes | not shown | |
  | **N. gonorrhoeae** | 7 | 10,223 | 2,000 | R1 0.646 (with 12 others) | yes | yes | yes |

  - **Every other species** fails the floor: 1,462 or fewer eligible, or a selection under 2,000.
  - **E. faecium** qualifies too and has more datasets. N. gonorrhoeae comes first because it is a block C species (block C calls only for N. gonorrhoeae and Salmonella).
  - **N. gonorrhoeae's own time** is at most its stage's total, 0.646 h for 13 species.
  - **S. pneumoniae** can't be shown to finish under the hour, because its stage's 1.538 h includes P. aeruginosa.
- **The gate:** `genome_order` names `neisseria_gonorrhoeae` = `03985e841ff3f9a931c98d9881ee99b0ae1891470e2b8235d4dd471b42279c3b` (master `species/neisseria_gonorrhoeae/genomes.txt`, 2,000 lines) in place of A. baumannii's.

## Log loss and discrimination

A log-loss `improves` can come from better calibration or less overfitting, not only from better discrimination. Every log-loss contrast record therefore carries J's AUROC delta (treatment minus baseline) and its 95% cluster jackknife-t interval (`auroc_delta`, `auroc_ci`). K prints them beside the decision as descriptive only; they decide nothing. K's reading rule 8 says the same in every report, including H-1's.

## Open questions

Main ruled on all four questions from the first build: the class balance (the balanced rule above), lineage leakage (the lineage drop), groups (MLST) and block screens (kept, ruling 4). After the null study, main restated the proof as the calibration above. One question is open: J's fold dependence for boosting on the random split ([Calibration](#calibration-of-the-decision-procedure-ruling-b)). The synthetic phenotype carries no k-mer signal, so the arms that change only the K model (k, budget) compare models fitted to noise. Their occasional `improves` or `worse` is a real log-loss difference between two such models, not a planted truth.

## Cost

An arm costs a curate plus count, structure and matrix for each new data root (cap, selection and k arms; blocks and budget arms reuse the default root), and one `grm arm` per arm and dataset: two splits × five outer folds × five learners, with no final models, MIC regression or annotation. On the synthetic tree (2,000 genomes of 23 kb, about 30,000 K patterns, one dataset, Boat xlarge, three arms at a time) a root took 1 to 5 s and an arm 7 to 50 s; the whole comparison took 176 s. Real genomes are about 200 times longer and real matrices hold far more patterns, so these numbers don't scale to a real species; a real arm costs about what `grm learn`'s cross-validation costs for that dataset, plus a count, structure and matrix of the species for each cap, selection or k arm. H-1 sizes the real panel.
