# Full-pipeline input comparison: the runner

Layth specified a full-pipeline comparison over the input data: genome and cohort count, which genomes are included, k-mer length, feature families and the retained-feature budget. Each input configuration goes through the whole pipeline and is compared with today's default. Layth delegated the open choices to main on 5 Nov 2026; main's decisions are DECISIONS V-1 (the values and pairing rules below), the rulings after lane INPUT (the balanced common test set, the lineage drop, MLST groups, block screens, and the seed and margin as decided values) and H-1 (real data comes later, on a 4-species panel). Nothing here is Layth's own acceptance.

**Built and proved on synthetic data only.** `scripts/input-comparison/run.py` runs every arm through `grm` and hands the predictions to J and K. Its end-to-end test (`tests/input_comparison/test_run.py`, `EndToEndTests`) runs all 11 arms plus the default on a synthetic data root with a planted effect. No real cohort has been run; that is H-1's, after revision 11.

[input-configuration.md](input-configuration.md) traces each input as the code has it now.

## The decision method H-1 uses: a test per learner (pre-registered, round 3, 7 Nov 2026)

Main's ruling "STATSFIX round 2 calibration, and the per-learner test". This section was committed before any round-3 calibration run, and it replaces the test of every section below.

**The unit.** Unit u is one (repeat r, outer fold f) of the lineage split's 3 keyed repeats and 5 outer folds, so J = 15 units, slots and replacement draws as in the sections below.
- m_u is the mean, over the genomes paired in u, of log loss(default) − log loss(arm), with probabilities clipped to [1e-4, 1 − 1e-4].
- m̄ = (1/J) Σ_u m_u and s² = Σ_u (m_u − m̄)² / (J − 1).

**The test by learner** (`stats.LEARNER_TESTS`):

| Learner | Test | SE | `ci_method` |
|---|---|---|---|
| scm, cart, rf | paired t over the unit means | √(s² / J) | `repeat_fold_t` |
| gbt, lasso | Nadeau and Bengio's corrected resampled t | √((1/J + ρ) · s²) | `repeat_fold_corrected_t` |

- **ρ** = (1/J) Σ_u n_test(u) / n_train(u). Here n_test(u) is the genomes paired in u, and n_train(u) is the default's genomes in repeat r outside fold f.
- **Statistic and intervals.** Both tests use t = m̄ / SE with df = J − 1 = 14. p is two-sided from Student's t(14). The 95% interval is m̄ ± t(0.975, 14) · SE and the 90% interval is m̄ ± t(0.95, 14) · SE.
- **The AUROC delta** takes the same test over the units' AUROC deltas.
- **The random split** stays descriptive, with the same test at J = 5 and df = 4.
- **Decisions are unchanged:** `improves` or `worse` needs the 95% interval off zero and p_adj < 0.05. BH runs over the run's non-sentinel lineage contrasts, and each sentinel is judged on its own p.
- **Why this split.** It rests only on synthetic null calibration. The paired t held size for scm, cart and rf with full power (round 1). For gbt and the lasso its rows failed, while the corrected t held their size (round 2).

**Fresh-seed confirmation.** It uses H-1's exact decision code above.
- **Seeds.** Calibration seeds 9401 to 9500, then reserve seeds 9501 upward, at most 20, while a judged row has fewer than 100 contrasts. The ruling named seeds 9301 to 9400 as never used, but round 2's power curve used 9301 to 9320, so the confirmation starts at 9401 to keep its condition that no seed was used before.
- **Size.** Every judged (twin, learner) row needs at least 100 contrasts, and a rate (improves + worse) of at most 5% with a Clopper-Pearson upper bound of at most 11%. Each direction is reported.
- **Power.** The 1× planted arms (k-h, k-p and k-c-h-p) must be `improves` on at least 95% of seeds for every learner, as the run decided them.
- **Power curve.** k-h and k-p on seeds 9601 to 9620, at visible shares 0.125, 0.25, 0.5 and 1. For each learner, the minimum detectable effect is the smallest share detected on at least 80% of seeds, with the mean log-loss gain at that share.
- **A learner that fails either bar** is reported, and H-1 treats it as exploratory: its decisions are reported but left out of the answer. Nothing is tuned further.

### Round 3's results (7 Nov 2026)

**Round 3 passes, so H-1 may decide with all five learners.** Every judged sentinel row meets the size bar on 100 new seeds, and every learner calls the three planted arms `improves` on all 100. No reserve seed was needed, and nothing was tuned after the first seed ran.

**Setup.** The runs used commit f7fc489, whose run code is the pre-registration's (db7ab2b). Nine Boat xlarge boxes with 16 vCPUs each built grm `cc19f66a…` from the same `git archive` at the same path. Calibration seeds 9401 to 9500 ran 20 per box on `bx_da8n6ht5` (EPYC 4585PX), `bx_j85dctp5`, `bx_ckuk23uz`, `bx_z9xu2u2u` and `bx_t2jnk9j8` (9950X). Power seeds 9601 to 9620 ran 5 per box on `bx_vrkfcpn8` (EPYC 4585PX), `bx_7t68uwtt` (EPYC 9455), `bx_37cmxyae` (9700X) and `bx_8ne8d6w4` (9950X).

**False positives.** Rows are the sentinels on their own p, with Clopper-Pearson 95% intervals:

| Twin | Learner | Test | Hits / contrasts | Rate [95% CI] | Improves | Worse | Result |
|---|---|---|---|---|---|---|---|
| learner-seed | rf | paired t | 2/100 | 2.0% [0.2, 7.0] | 1 | 1 | pass |
| learner-seed | gbt | corrected t | 0/100 | 0.0% [0.0, 3.6] | 0 | 0 | pass |
| retuned | scm | paired t | 2/100 | 2.0% [0.2, 7.0] | 1 | 1 | pass |
| retuned | cart | paired t | 4/100 | 4.0% [1.1, 9.9] | 4 | 0 | pass |
| retuned | gbt | corrected t | 0/100 | 0.0% [0.0, 3.6] | 0 | 0 | pass |
| retuned | lasso | corrected t | 0/100 | 0.0% [0.0, 3.6] | 0 | 0 | pass |

Retuned cart comes closest to the bar: 4 of 100, all `improves`, with an upper bound of 9.9% against 11%. The other rows are zero by construction and not judged.

**Planted arms (pass).** k-h, k-p and k-c-h-p are `improves` on 100 of 100 seeds for every learner, against the 95% bar. In round 2, cart missed 33 of 100 under the corrected t. Under the paired t it misses none.

**Power curve.** The minimum detectable effect is the smallest visible share detected on at least 80% of seeds. k-p gives the same table:

| Learner | 0.125 | 0.25 | 0.5 | 1 | Minimum detectable (share, mean log-loss gain) |
|---|---|---|---|---|---|
| scm | 20/20 | 20/20 | 20/20 | 20/20 | 0.125, 1.86 |
| cart | 16/20 | 20/20 | 20/20 | 20/20 | 0.125, 1.04 |
| rf | 0/20 | 4/20 | 17/20 | 20/20 | 0.5, 0.012 |
| gbt | 0/20 | 0/20 | 10/20 | 20/20 | 1, 0.47 |
| lasso | 1/20 | 12/20 | 19/20 | 20/20 | 0.5, 0.16 |

gbt still needs every carrier to show the family. At a quarter, its mean gain of 0.018 was detected on 0 of 20 seeds. So a gbt gain near stage 2's 0.02 margin would usually read as no measured difference.

**H-1's learners.** All five pass both pre-registered bars, so none is exploratory. The power curve is not a bar.

**Receipts.** `grm-collect/statsfix-r3` on the Boat controller holds `reports/r1/` (the calibration and power reports) and `boxes/<box>/` (every seed's files), with `SHA256SUMS`.

### After REVIEW-STATSFIX-R3 (7 Nov 2026)

REVIEW-STATSFIX-R3 approved round 3 with three Lows, and none of them changes its results.

- **L1.** `calibrate()` found the zero rows in the data. So a judged row with no contrast on any seed passed, and so did a row that should vary but came out all zero. It now names the nine zero rows the calibration's rule lists (`ZERO_ROWS`) and judges every other row:
  - A zero row fails on any nonzero contrast.
  - A judged row fails when every contrast is exactly zero.
  - An undefined contrast counts toward no row's 100, and each row reports its `undefined` count.
- **L2.** A planted record must carry its learner's test, as a twin's record already did.
- **L3.** The calibration report's `test` field gives `stats.LEARNER_TESTS`. Both specs, `stats.py` and this doc's AUROC section no longer name one test, an every-seed bar or the cluster jackknife for it.

Rerun over round 3's 540 receipt files, the fixed `null_study.py report` gives byte-identical `calibration.txt`, `power.json` and `power.txt`. `calibration.json` differs only in its `test` field and the 15 rows' new `undefined` counts, all 0. The calibration still passes, with no failed row and no missed planted arm.

The four tests for these fixes fail on round 3's `null_study.py`, and seven mutants of the new checks are each caught by the test they name.

## The decision method (main's rulings after VERIFY-H1STATS, 7 Nov 2026)

REVIEW-ICODE-PY's H1 and VERIFY-H1STATS showed that the lineage-cluster test was anti-conservative: one pair of fitted models scores every lineage of an outer fold, so a fold-level shift moves all of that fold's clusters together, and gbt's learner-seed twin rejected in 26 to 46% of seeds. The method below replaces it. H-1 runs it only after its calibration passes, and [the first calibration](#the-calibrations-results-7-nov-2026) failed for gbt and the lasso. The null-study and calibration figures in the sections after this one came from the old test and are **superseded**.

1. **The test.** `grm arm` scores the spec's `repeats` keyed lineage plans (`spec.json` says 3) and one random plan. J pairs each arm with the default per repeat (`stats.pair_units`) and gives each (repeat, outer fold) unit the mean log-loss delta of its paired test genomes. The decision is a paired t over the 15 unit means, df 14 (`stats.compare_units`, `ci_method` `repeat_fold_t`, unit `repeat_fold`): `improves` or `worse` needs the 95% t interval off zero and p_adj < 0.05, as before. The AUROC delta is the mean of the units' AUROC deltas with its 95% t interval. The random split runs once and stays descriptive (family D: estimate and interval over its 5 units, no p).
2. **Seeds keyed by repeat.** Each repeat's plan is keyed (`plan::repeat_seed(seed, repeat)`). A twin's seed `s` (learner, balance or tuning) acts as `repeat_seed(s, repeat)` in each repeat, so a twin draws anew in every repeat as the default does; at repeat 0 it is `s`, as before.
3. **Cohort arms pair per repeat.** A cohort arm's balanced subsample can differ by repeat, so J pairs it with the default's rows of the same repeat. The refusal of repeats > 1 with a cohort arm (PR #69's review) is gone. Under `common_outer_test_set_balanced`, J also checks the subsample is the balanced draw of the default's whole test fold: every default fold of every repeat is scored, with both classes, in counts `grm arm`'s `balanced_counts` gives (REVIEW-ICODE-PY L1).
4. **A refused repeat.** When the default's lineage draw for a repeat is refused, `grm arm` draws the next keyed repeat past the slots (3, then 4), at most 2 replacements per dataset. Each refused draw is kept as `plans/refused-lineage-<repeat>.json` (summary.json's `refused_plans` lists it) and `repeats.json` lists every slot's keyed repeat. With fewer than 3 valid repeats after that, every lineage run of the arm is unavailable with that reason, and its contrasts are not produced. **A dataset with no decision-bearing lineage contrast can't be tested this way, and the run still finishes:** summary.json gives `decisions_made` per dataset and names it in `no_decision_datasets`, the report lists it under "Datasets the per-dataset test can't answer" with its planned, not produced and undefined counts, and `run.py` prints a `NO LINEAGE ANSWER` line and exits 0. This is Layth's ruling on H-1's pilot (2026-11-09), where Salmonella ciprofloxacin's lineage draws were refused. It replaced the exit 3 of REVIEW-ICODE-PY M3, which stopped the pilot after its report was written.
5. **Sentinels.** Every run adds two arms exchangeable with the default. `sentinel-seed-twin` (learner seed 1001) judges rf and gbt, the learners its seed reaches. `sentinel-tuning-twin` (tuning seed 1002) judges scm, cart and the lasso: it keeps the default's outer folds and draws every outer fold's inner folds and both tuning holdouts again (`plan::retuned_plan`). Those draws are the only randomness those learners have. The default's own draw is the same procedure at its own seed (Rust test `a_retuned_plan_keeps_the_outer_folds_and_redraws_tuning`), so the twin is exchangeable with it. VERIFY-H1STATS listed a bootstrap or subsample of each training fold as other options; both train on other data than the default, so they aren't exchangeable with it and would flag real differences. The manifest records `sentinel_learners`. A sentinel's family C contrast is judged on its own p, outside BH (`adjust` none, `p_adj` = p; the manifest's `sentinel_adjust`), so its sensitivity no longer depends on the other arms (VERIFY's M1). When it isn't `inconclusive` for a learner it judges, or isn't produced, K marks that learner's decisions on the species unreliable.
6. **BH's family.** One family over every non-sentinel lineage contrast of the run's species (REVIEW-ICODE-PY L4). For H-1's panel the family is every lineage contrast of every species run: `run.pending_contrasts` gives each run's records before adjustment, and `stats.finalize_sets` adjusts several runs' sets as one family, for `run.py combine`. K recomputes BH over the planned non-sentinel family C ids.
7. **Inputs.** Cap and selection are two inputs, so an arm changing both is refused (REVIEW-ICODE-PY L2). A RAW holding `truth.json` (a synthetic root's marker) refuses to run unless the spec says `"synthetic": true`, so a synthetic run can't be recorded as real (L3).
8. **Stage 2's entry, D3 as amended (ruling 5).** `scripts/input-comparison/stage2.py`: a stage-1 setting enters on a dataset when at least one reliable learner is `improves` on the lineage unit t with its log-loss interval's lower end above the 0.02 margin and its AUROC-delta interval above 0, and no reliable learner is `worse`. Blocks that enter combine by union; of two entering values of one other input, the larger mean log-loss gain wins, then the value nearer the default, then the smaller budget; k joins only after the k add-on. A spec that says `"stage": 2` may hold an arm changing the blocks and one other input, paired by that input's rule (never two inputs besides the blocks, so never k with a cohort change). Stage 2 itself runs later and is exploratory.

### The calibration H-1 must pass (ruling 6)

Pre-registered here before the run. `null_study.py run --spec calibration.json` with H-1's exact decision code (the unit t at 3 repeats, the sentinels on their own p, BH over the run's other lineage contrasts) on synthetic seeds 9101 to 9200. `calibration.json`'s arms: `cohort-path` (balance seed 7) and the planted `k-h`, `k-p` and `k-c-h-p`; the runner adds both sentinels. `null_study.py report --calibration` then judges:
- **Rows.** One row per twin (learner-seed twin, retuned twin, cohort path) and learner, counting that twin's lineage contrast once per seed, so a row's contrasts are independent. `improves` and `worse` are counted apart, each with a Clopper-Pearson 95% interval. A row whose every contrast is exactly zero is zero by construction (the twin doesn't reach that learner: the learner-seed twin for scm, cart and the lasso, the retuned twin for rf, the cohort path for all) and is reported, not judged.
- **Bar.** Every judged row needs at least 100 contrasts, improves + worse at most 5%, and the interval's upper bound at most 11%. Every learner needs a judged row.
- **Power.** Every planted arm is `improves` for every learner on every seed, as the run decided it.
- **Refused plans.** A seed whose run made no lineage decision (its summary's `no_decision_datasets`) is listed and counts toward no row. If any judged row ends with fewer than 100 contrasts, reserve seeds 9201 upward run in order, at most 20, until every row has 100; every seed run is reported.

### The calibration's results (7 Nov 2026)

**It fails, so H-1 can't run on this method as it stands.** Boat xlarge `bx_wybets94` (9950X), grm `b9c532d8…` built at `2cfa4278`, seeds 9101 to 9200 and then reserve seeds 9201 to 9207. The reserve was needed because the retuned twin went unproduced on 7 seeds, so its rows reached 100 contrasts only at seed 9207. All 107 seeds made lineage decisions. Rows are the sentinels on their own p, improves + worse, with Clopper-Pearson 95% intervals:

| Twin | Learner | Hits / contrasts | Rate [95% CI] | improves, worse | Result |
|---|---|---|---|---|---|
| learner-seed | rf | 3/107 | 2.8% [0.6, 8.0] | 0, 3 | pass |
| learner-seed | gbt | 9/107 | **8.4% [3.9, 15.4]** | 6, 3 | **fail** |
| learner-seed | scm, cart, lasso | 0/107 | zero by construction | – | not judged |
| retuned | scm | 2/100 | 2.0% [0.2, 7.0] | 1, 1 | pass |
| retuned | cart | 0/100 | 0.0% [0.0, 3.6] | 0, 0 | pass |
| retuned | gbt | 2/100 | 2.0% [0.2, 7.0] | 1, 1 | pass |
| retuned | lasso | 9/100 | **9.0% [4.2, 16.4]** | 7, 2 | **fail** |
| retuned | rf | 0/100 | zero by construction | – | not judged |
| cohort path | all five | 0/107 | zero by construction | – | not judged |

- **Power.** k-h, k-p and k-c-h-p were `improves` for every learner on all 107 seeds (BH over each run's lineage contrasts).
- **Not produced.** The retuned twin's redrawn tuning partitions were refused on 7 of 107 seeds: on 6 the gbt holdout lacked a class (outer fold 0, 1 or 3, two seeds each), and on 1 an inner fold did. On real data such a seed would flag scm, cart and the lasso unreliable on the species.
- **The two failing rows** are both learners whose fits move with a holdout draw: gbt's trees on the learner seed, the lasso's penalty on its holdout. Their hits mostly have p between 0.02 and 0.05. That pattern points to a test somewhat anti-conservative for those learners rather than to a few extreme seeds; VERIFY-H1STATS measured gbt at 6.2% [2.5, 12.5] at R = 3 and attributed the creep to units that share samples across repeats.
- **Receipts.** `logs/cal-box/grm-logs/cal/calibration.json` and `calibration.txt`, every seed's `summary.json` and `contrasts.json` in `cal-seeds.tgz`, and a laptop recomputation from those files that equals the box's report (`logs/calibration-local.json`). All in the STATSFIX lane folder.

### Round 2, pre-registered before any rerun (main's ruling "STATSFIX calibration fails", 7 Nov 2026)

This section fixes the test, the retuned twin's draws, the power curve and the rerun before any of them run. Where it differs from the sections above, it replaces them.

**The test: the corrected resampled t (Nadeau and Bengio), for every learner.**
- **Units.** Unit u is one (repeat r, outer fold f) of the lineage split's 3 keyed repeats and 5 outer folds, so J = 15 units, with the slots of the replacement rule above.
- **Unit value.** m_u is the mean, over the genomes paired in u, of log loss(default) − log loss(arm), with probabilities clipped to [1e-4, 1 − 1e-4].
- **Size ratio.** n_test(u) is the number of genomes paired in u. n_train(u) is the number of the default's genomes in repeat r outside fold f, its training samples in that fold. ρ = (1/J) Σ_u n_test(u) / n_train(u), which is about 0.25 for the default's own folds.
- **Statistic.** m̄ = (1/J) Σ_u m_u and s² = Σ_u (m_u − m̄)² / (J − 1). SE = √((1/J + ρ) · s²), t = m̄ / SE, df = J − 1 = 14.
- **p and intervals.** p is two-sided from Student's t with 14 df. The 95% interval is m̄ ± t(0.975, 14) · SE and the 90% interval is m̄ ± t(0.95, 14) · SE.
- **Decisions** are unchanged: `improves` or `worse` needs the 95% interval off zero and p_adj < 0.05, BH over the run's non-sentinel lineage contrasts, and a sentinel on its own p.
- **AUROC delta.** Its interval uses the same SE formula over the units' AUROC deltas.
- **Random split.** It stays descriptive and uses the same formula with J = 5 and df = 4.
- **Degenerate units.** Unit means equal to rounding give p = 1 when they are all zero; otherwise the contrast is undefined.
- **Record.** `ci_method` is `repeat_fold_corrected_t` and the unit is `repeat_fold`.

**The retuned twin's replacement draws.** When a slot's retuned tuning partitions are refused, `grm arm` draws them again from `repeat_seed(repeat_seed(tuning seed, keyed repeat), d)` for d = 1, then 2; d = 0 is today's draw. That allows at most 2 replacements per split, each refused draw kept as `plans/refused-<split>-<slot>-tuning-<d>.retuned.json` and listed in `repeats.json`. Past that, the twin's split is unavailable with the reason.

**The power curve.** The 2× size in the ruling's example can't be built: the planted effect already explains every resistant genome, since resistance is carriage of a blaKPC allele that H and P call. So the sizes are shares of the resistant genomes whose blaKPC allele is visible:
- **Visible shares.** v = 0.125, 0.25, 0.5 and 1 (1 is today's planted effect). `synthetic_root.py --visible v` keeps each carrier resistant but writes its allele contig only with probability v, drawn from its own random stream. At v = 1 the root is byte-identical to today's.
- **Arms and seeds.** For each v, `power.json` runs `k-h` and `k-p` against the default on seeds 9301 to 9320, with 3 repeats and the corrected t. The runner's sentinels run too.
- **Detection.** A learner detects arm a at size v on a seed when a's lineage contrast is `improves` in that run.
- **Minimum detectable effect.** For each learner and arm, the smallest v detected on at least 80% of the seeds, with its observed share and Clopper-Pearson interval. If no tested size reaches 80%, it is reported as above 1×. This is the minimum detectable effect that H-1's "no measured difference" is stated against, as a visible share and as the mean log-loss gain the arms showed at that size.

**The rerun.** The calibration reruns as pre-registered above:
- seeds 9101 to 9200, then reserve seeds 9201 upward while a judged row has fewer than 100 contrasts, at most 20;
- the same bar: at least 100 contrasts, rate at most 5% and upper bound at most 11% per judged row, each direction reported;
- the planted arms `improves` for every learner on every seed.

If any row still fails, the lane stops and reports, and nothing is tuned further.

### Round 2's results (7 Nov 2026)

**The rerun fails on power, so H-1 still can't run.** The false-positive rate is now controlled, but cart misses the planted arms on 33 of 100 seeds. Per the ruling, nothing was tuned after this.

**Setup.** Boat xlarge `bx_y9r6ck25` (9950X), grm `13ef526a…` built at `c7ee4225`, seeds 9101 to 9200. No reserve seed was needed: the retuned twin was produced on all 100.

**False positives.** Rows are the sentinels on their own p, improves + worse, with Clopper-Pearson 95% intervals:

| Twin | Learner | Hits / contrasts | Rate [95% CI] | Result |
|---|---|---|---|---|
| learner-seed | rf | 0/100 | 0.0% [0.0, 3.6] | pass |
| learner-seed | gbt | 0/100 | 0.0% [0.0, 3.6] | pass |
| retuned | scm | 0/100 | 0.0% [0.0, 3.6] | pass |
| retuned | cart | 0/100 | 0.0% [0.0, 3.6] | pass |
| retuned | gbt | 0/100 | 0.0% [0.0, 3.6] | pass |
| retuned | lasso | 0/100 | 0.0% [0.0, 3.6] | pass |
| zero by construction | (as in round 1) | 0/100 | – | not judged |

Every judged row is 0 of 100 in both directions. The corrected t holds its size with room to spare and may be conservative.

**Planted arms (fails).** k-h, k-p and k-c-h-p are `improves` on all 100 seeds for scm, rf, gbt and the lasso, but for cart on only 67 of 100, with the same 33 seeds missed by all three arms. On a missed seed, cart's contrast typically has its interval off zero and p just under 0.05, but p_adj over 0.05 after BH (seed 9104: estimate 1.39, 95% interval [0.012, 2.76], p 0.048, p_adj 0.064). So the miss comes from the corrected t's wider SE (×2.18 at ρ = 0.25) combined with BH, not from missing signal. On this root the three planted arms give identical contrasts: H and P each add one column calling the same blaKPC family, so they are one test, not three.

**Power curve** (`bx_txpqa8m9`, seeds 9301 to 9320 at each size, the corrected t). The minimum detectable effect is the smallest visible share detected on at least 80% of seeds:

| Learner | 0.125 | 0.25 | 0.5 | 1 | Minimum detectable (share; mean log-loss gain there) |
|---|---|---|---|---|---|
| scm | 4/20 | 9/20 | 19/20 | 20/20 | 0.5; 2.95 |
| cart | 0/20 | 1/20 | 3/20 | 16/20 | 1 (exactly 80%); 1.67 |
| rf | 0/20 | 0/20 | 5/20 | 20/20 | 1; 0.075 |
| gbt | 0/20 | 0/20 | 11/20 | 20/20 | 1; 0.47 |
| lasso | 0/20 | 11/20 | 20/20 | 20/20 | 0.5; 0.16 |

k-p gives the same table. So under the corrected t, rf, gbt and cart detect the planted family reliably only when every carrier shows it; at half visibility they mostly don't. For rf, a mean log-loss gain of about 0.012 (half visibility) was detected on only 5 of 20 seeds, and only the full effect, about 0.075, reached 80%. Both are on the scale of stage 2's 0.02 margin, so a gain near that margin would usually read as no measured difference for rf.

**Receipts.** In the STATSFIX lane folder:
- `logs/r2/cal-box/grm-logs/cal/` (`calibration.json`, `calibration.txt`, `cal-seeds.tgz`);
- `logs/r2/gates-box/grm-logs/gates/` (`power.json`, `power.txt`, `power-seeds.tgz`).

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
2. **Data roots.** One per distinct (cohort, k): the default root, one per cap or selection arm, and one per k arm. Each links `RAW/bvbrc`, `RAW/fasta` and `RAW/catalog`, then runs `grm curate --cap --min-class --selection [--selection-seed --stratum]`, `grm count --k`, `grm structure` and `grm matrix`. Blocks and budget arms reuse the default root. A root counted at the pool's k takes the pool's ksets as hard links before its count ([H-1 support](#h-1-support-h1code-7-nov-2026)).
3. **Arms.** `grm arm` (`src/ablate/arm.rs`, `run_arm`) scores one arm of one dataset.
   - It admits the arm's blocks through H's adapters: `KmerSource`, `ProbeSource::from_ksets`, `PanelSource::new` over the counted gene calls, and `PresenceSource::new` over the ksets.
   - For each split (random, lineage), repeat and outer fold, it runs `fit_encoder` and `transform` on the fold's training samples, `block_candidates` there, `assemble` on the same mask, and every learner through `learn::planned::evaluate_plan` with the plan's inner folds and holdouts. Each arm is scored once per outer fold.
   - The default arm runs first and writes the fold plans the others read (`plans/<split>-<repeat>.json`: sample ids, labels, the species grouping's clusters and duplicate groups, and the keyed `FoldPlan`).
4. **Rows.** `grm arm` writes contract §4 rows (`predictions.jsonl`): `fold` from the plan, `called` = score ≥ 0.5, and a no-call row scoring its outer fold's training prevalence with `called` false. It also writes `runs.json` (K's run record without its id, plus a cohort arm's drop counts) and `arm.json` (the arm's k, sample count, support floor, eval config, selection record, block specs and unavailable blocks).
5. **J.** For each dataset, arm, learner and split, `stats.pair_units(..., pairing=rule)` and `stats.compare_units` (family C on the lineage split, D on the random split), then one `stats.finalize`: BH over every planned non-sentinel family C contrast, the sentinels on their own p ([the decision method](#the-decision-method-mains-rulings-after-verify-h1stats-7-nov-2026)). A pairing J refuses, or an arm with no result, is not produced, with its reason.
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
- all 130 planned contrasts produced (11 arms and the two sentinels);
- the lasso's decision `improves` on the lineage split for `k-h`, `k-p` and `k-c-h-p`, and every random-split contrast descriptive (family D, an estimate but no p);
- 3 lineage repeats and 1 random in every run, every lineage contrast on 15 (repeat, fold) units, the retuned twin keeping the default's outer folds and drawing other inner folds in each repeat, and the sentinels outside BH;
- both sentinels present in the manifest, and a learner flagged unreliable exactly when the lineage-split contrast of the sentinel that judges it isn't `inconclusive`;
- no `improves` anywhere for `k-c`;
- every contrast naming its arm's rule;
- the report reloading against its prediction files.

Without the variable only the spec tests run.

On the rerun after the rulings (Boat xlarge, 5 Nov 2026, seed 11) the test passed 7/7. All 110 contrasts were produced; `K+H`, `K+P` and `K+C+H+P` were `improves` for every learner on both splits, and `K+C` had no `improves`.

## The null study for the cohort arms

**Superseded (7 Nov 2026):** every figure in this section and the next two came from the lineage-cluster test, which [the decision method](#the-decision-method-mains-rulings-after-verify-h1stats-7-nov-2026) replaced. They aren't H-1's calibration and no longer describe what the runner decides; they are kept as the record of what was measured. The cohort kinds' base rates in `spec.json` were measured the same way.

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

Main's ruling on this study (5 Nov 2026) withdrew the first criterion: with no k-mer signal, two different models still differ in expected log loss, so "`improves` at no more than 5%" can't hold for an arm that changes the model. The table is now context (ruling (a)): the cohort arms' no-signal `improves` rate, 15.5%, against the same-sample arms' 17.6%. The decision procedure's calibration is the next section. **Superseded:** these rates are over both splits. On the lineage split alone, the only split that decides now, the order reverses (21.0% against 13.0%, a point estimate with overlapping intervals); see [After PR #69's Codex review](#after-pr-69s-codex-review-mains-rulings-6-nov-2026).
- **Same-sample arms do it too.** Arms with the default's exact samples and fold plans (k and budget) score `improves` at 17.6%, more often than the cohort arms over both splits (not on the lineage split; see the pointer above). With no k-mer signal, two different models fitted to noise really do differ in log loss (how far they overfit, how well they calibrate), and the paired test finds those differences. "No k-mer signal" doesn't make log loss equal, so no pairing rule can bring a model-changing arm to the nominal rate here. Only K+C, which adds at most two columns and leaves the model nearly unchanged, is near it.
- **The balance artefact is gone.** In the first build, `cap-1000` was `improves` for gbt, lasso and rf from the 65%-resistant test set. Under the balanced rule those learners are `worse` (rf 40 of 40), as expected for an arm trained on half the genomes at another prevalence and scored at the default's balance.
- **Contrasts within a seed are correlated** (one cohort, one test set), so the per-contrast intervals are optimistic; every seed had at least one `improves` for cap-4000 and grr-v1.

Receipts: `logs/null-study-20seeds.json` and each seed's `summary.json` and `contrasts.json` in the lane folder.

## Calibration of the decision procedure (ruling (b))

**Superseded (7 Nov 2026)** by [the calibration H-1 must pass](#the-calibration-h-1-must-pass-ruling-6): this one judged BH-adjusted decisions averaged over all five learners, which hid gbt's rate (REVIEW-ICODE-PY H1). `calibration.json` now holds the new calibration's arms.

The calibration of 5 Nov defined two arms that are exchangeable with the default, so every `improves` or `worse` on them is a false positive:
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

**Superseded in part (7 Nov 2026):** the random split stays descriptive and the lineage split decides, but the test, the sentinels and BH's family are now [the decision method](#the-decision-method-mains-rulings-after-verify-h1stats-7-nov-2026)'s.

The calibration above (seed twin over 20 seeds: 2.0% `improves` and 3.0% `worse` over both splits, 1.0% and 0.0% over the lineage split alone; gbt on the random split non-inconclusive in 9 of 20; gbt on the lineage split 1 of 20; rf 0 of 40) led main to three rules, which the runner and K now apply.

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
- **Short cohorts (9 Nov 2026).** A cohort arm's own cohort can leave a dataset with fewer genomes in a class than the class minimum, 25. `grm arm` then refuses it ("short of its class minimum"). On real data, H-1's NG-1 stopped there: grr-v1's 2,000 genomes, drawn round robin over MLST groups, hold too few cefixime-resistant ones. `run.py` now records such an arm's runs on that dataset as `unavailable`, with the reason `short_of_class_minimum`, and prints a `SHORT COHORT` line. K lists them under "Unavailable runs", every planned contrast of the arm on that dataset is not produced with that reason, and the run goes on. This applies Layth's ruling on H-1's pilot ("Report it, keep going") to a cohort arm. The same refusal of the default, or of an arm scored on the default's cohort, still stops the run. So does any other `grm arm` failure. Tests: `test_a_cohort_arm_short_of_a_class_minimum_is_unavailable_and_the_run_goes_on`, `test_the_default_s_class_minimum_refusal_still_stops_the_run`, `test_only_grm_s_own_refusal_of_a_cohort_arm_counts_as_short` and `test_a_short_cohort_arm_s_contrasts_are_not_produced_and_the_rest_decided`.
- **A block some datasets refuse (9 Nov 2026).** H's panel source refuses H on any drug without panel targets ("no panel targets <drug>"). On H-1's SAL-1 that left k-h and k-c-h-p without H on 11 of 13 datasets. K keeps one block list for the whole run, and the part's first dataset was one that refused H.

  `run.py` made that refusal the list's H entry. K's schema then refused every ok run with H on the 2 datasets with a panel, and the run stopped at K's report with exit 1.

  Now a block that any dataset admits keeps an admitted spec, and each refusing dataset's runs stay `unavailable` with the block's reason. H's spec is per drug, so the list holds the first admitting dataset's, and each dataset's own spec is in its arm.json. A block no dataset admits is unavailable with each distinct reason, in dataset order. K's schema is unchanged. Test: `test_a_block_one_dataset_refuses_stays_available_for_another_dataset_s_runs`, which fails on 63c5058 with SAL-1's error.
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
- **Ksets are counted once per root (M1(b), as measured).** `grm count` reuses a kset only when the same root's own kset folder already holds it for the same genome and k with a matching sidecar (`species::count_one`). Every comparison root has its own kset folder, so a genome in the pool and in three arm roots is counted four times, and the estimate counts it four times. No counting behaviour changed. Test: `test_every_root_counts_its_own_ksets` (e2e: the pool's and the default's kset of one genome are two files with the same bytes). **Superseded for k31 roots** by kset sharing ([H-1 support](#h-1-support-h1code-7-nov-2026)); `--no-share-ksets` keeps this behaviour for the identity proof.
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

## After PR #69's Codex review (main's rulings, 6 Nov 2026)

- **The null study counts decision-bearing lineage contrasts only.** `null_study.py report` read every row of `summary.json` `decisions`. Since main's ruling after the calibration study, the random split's rows read `descriptive`, so they sat in every denominator.
  - **Reproduced** (`logs/round16/repro-before.json` in the lane folder): a planted arm improving on all its lineage contrasts scored 0.5, so the study could never pass, and a true lineage false-positive rate of 0.2 read 0.1.
  - **Now:** rates and acceptance use the lineage split's contrasts. A contrast not produced is counted apart.
- **Calibration sees what the sentinel masks (P1).** K marks a learner's rows `unreliable` when the sentinel fires, and the study counted those as neither `improves` nor `worse`.
  - **Reproduced:** a sentinel firing on every seed passed calibration with 0 false positives.
  - **Now:** the study reads J's own decisions (`raw_decisions`; a summary written before the mask has only `decisions`, which are then J's), and reports the sentinel's firings apart (`sentinel_firings`). The same case now fails, with 40 firings in 20 seeds.
- **`repeats` above 1 with a cohort arm is refused.** A cohort arm's balanced subsample can differ by repeat, and `stats.pair_arms` refuses that ("covers different samples in different repeats"); reproduced. Both decided specs use `repeats` 1, so `load_spec` refuses such a spec up front.
- **Tests:** `test_null_study_rates_count_decision_bearing_lineage_contrasts_only`, `test_null_study_calibration_sees_what_the_sentinel_masks` and `test_a_spec_with_repeats_and_a_cohort_arm_is_refused`.
- **The rates, re-measured** from the stored per-seed summaries. The old column is `logs/round16/remeasure.json` in the lane folder. The new column is lane INPUT17's recomputation, `tools/rebh17.py` → `logs/rebh17.json` (lineage-only Benjamini-Hochberg, below); `remeasure.json`'s lineage figures keep the stored decisions.
  - **Why they change:** the 20-seed null study and the calibration ran before random-split contrasts became descriptive and before the sentinel existed. Their old rates were over both splits, as both decided then, and no summary held a `descriptive` or `unreliable` row. The new rates are over the lineage split alone, the only split that decides now, with each lineage contrast decided as the code decides it now: family C holds the lineage split only, so Benjamini-Hochberg adjusts each seed's 55 lineage contrasts, not the 110 of both splits the stored decisions were adjusted over (REVIEW-NULL16 L3). Seeds 101–120; Clopper-Pearson 95% intervals.

  | Arm | Old, both splits | New, lineage only | `worse` old → new |
  |---|---|---|---|
  | cap-1000 | 13/200 = 6.5% | 6/100 = 6.0% [2.2, 12.6] | 133 → 65 |
  | cap-4000 | 45/200 = 22.5% | 33/100 = 33.0% [23.9, 43.1] | 40 → 10 |
  | grr-v1 | 35/200 = 17.5% | 24/100 = 24.0% [16.0, 33.6] | 55 → 12 |
  | cohort arms pooled | 93/600 = 15.5% | 63/300 = 21.0% [16.5, 26.1] | |
  | k21, k25, budget-1000, budget-16000 pooled | 141/800 = 17.6% | 52/400 = 13.0% [9.9, 16.7] | |
  | K+C (now `k-c`) | 2/200 = 1.0% | 2/100 = 2.0% | 1 → 0 |
  | K+H, K+P, K+C+H+P | 200/200 each | 100/100 each | 0 |
  | calibration, seed twin (b2) | 4/200 = 2.0% improves, 6/200 = 3.0% worse | 1/100 = 1.0% improves, 0/100 worse | |
  | calibration, cohort path (b1) | 0/200 | 0/100 | |

  - **Lineage-only Benjamini-Hochberg moves 4 of 1,100 decisions** (`logs/rebh17.json` in lane INPUT17's folder; J's own `stats.bh` reproduces all 2,200 stored `p_adj` over both splits first). k21 loses two `improves` (20 → 18), k25 one (19 → 18), and grr-v1 gains one (23 → 24). With the stored decisions the lineage figures were grr-v1 23/100, cohort pooled 62/300 = 20.7% and same-sample pooled 55/400 = 13.8%. No `worse` moves, and no calibration decision moves.
  - **One sentence above no longer holds:** on the lineage split, the cohort arms' no-signal `improves` rate (21.0%) is above the same-sample arms' (13.0%), not below it. That is a point estimate, not a clear difference: the two Clopper-Pearson intervals overlap ([16.5, 26.1] and [9.9, 16.7]), and both are optimistic because contrasts within a seed are correlated; cap-1000 alone (6.0%) is below the same-sample pool. The point stands that model-changing arms differ in log loss with no k-mer signal.
  - **Calibration still passes,** with both arms at or below 5% on the lineage split. The stored seeds predate the sentinel (no sentinel arm ran and no summary has an `unreliable` list), so they say nothing about how often it fires.
  - **The per-kind base rates in `spec.json`** (cap-down 14.0%, cap-up 16.0%, selection-swap 18.0%) didn't come from the null study. The re-keying twin counted them from J's lineage-split records, and they recompute unchanged.
  - **K's rule 9** now gives the seed twin's lineage figures (1.0% `improves`, 0.0% `worse`) beside the both-split ones.

## Log loss and discrimination

A log-loss `improves` can come from better calibration or less overfitting, not only from better discrimination. Every log-loss contrast record therefore carries J's AUROC delta (treatment minus baseline) and its 95% interval by the record's own test (`auroc_delta`, `auroc_ci`): the learner's unit test in the input comparison, the cluster jackknife-t elsewhere. K prints them beside the decision as descriptive only; they decide nothing. K's reading rule 8 says the same in every report, including H-1's.

## Open questions

Main ruled on all four questions from the first build: the class balance (the balanced rule above), lineage leakage (the lineage drop), groups (MLST) and block screens (kept, ruling 4). After the null study, main restated the proof as the calibration above. One question is open: J's fold dependence for boosting on the random split ([Calibration](#calibration-of-the-decision-procedure-ruling-b)). The synthetic phenotype carries no k-mer signal, so the arms that change only the K model (k, budget) compare models fitted to noise. Their occasional `improves` or `worse` is a real log-loss difference between two such models, not a planted truth.

## After REVIEW-INTEG7 (main's rulings, 6 Nov 2026)

- **"synthetic": true needs a synthetic RAW (L2).** `check_synthetic` runs right after the curate-input check, before RAW is hashed or any root curates, in `--estimate-only` too. It refuses unless RAW's `truth.json` has synthetic_root.py's format (`grm-input-comparison-synthetic`) and lists exactly the genome ids of RAW's `bvbrc/genomes.jsonl`.
  - **Why this marker.** The format string alone could arrive in a copied file. Tying it to RAW's own genome list means the marker describes this data root, and a real root's BV-BRC metadata can't match a generator's list by accident. A dedicated flag file would be weaker still: it carries nothing about the data.
  - A synthetic spec that also carries a gate still runs the gate.
- **The gate decision on stderr.** At the start of every run, before anything is counted, run.py prints `revision-12 gate: PASSED: …` (what it compared, and "the spec also says synthetic; the gate ran" when it does) or `revision-12 gate: NOT REQUIRED: …` for a synthetic spec with no gate.
- **K records it.** The manifest's optional `revision12_gate` is `passed` or `not_required`, and the schema refuses `not_required` unless `synthetic` is true. The report header says which: "Synthetic run: yes. The spec also carries a revision-12 gate, and it passed.", "Synthetic run: yes, on a synthetic data root. The spec carries no revision-12 gate, so none ran; a real run needs it." or "Synthetic run: no. The revision-12 gate passed." A manifest without the field (written before this ruling) keeps the old wording and says it doesn't record whether a gate ran.
- **Tests.**
  - `test_a_synthetic_spec_refuses_a_raw_that_is_not_a_synthetic_root`: no `truth.json`, another format, a marker from another synthetic root, and one genome short, in both modes, with no root curated.
  - `test_the_gate_decision_is_printed_on_stderr`: both lines.
  - `test_the_synthetic_flag_is_named_in_the_header` (K): the header for each value, and the schema refusals.
  - The e2e checks that its manifest records `passed`.
  - The other stub tests use RAWs that aren't synthetic roots, so they stand the check down and test what comes after it.
- **Added after REVIEW-INTEG7B:** the refusal test also asserts that WORK holds no `inputs.json`, which pins the check before RAW is hashed (L1), and `test_the_manifest_records_whether_the_gate_ran` runs `run.main` past the manifest with the counting, arm and contrast steps stubbed: `not_required` for an ungated synthetic run, `passed` with the gate (L2).

## After REVIEW-NULL16 (main's rulings, 6 Nov 2026)

- **The null study needs 20 seeds with a lineage contrast per arm (L1).** `pass` counted summary files, so a calibration arm with lineage contrasts in 1 of 20 seeds (95 not produced) passed on "improves 0/5". Now every arm the study judges needs a lineage contrast in at least 20 seeds. The calibration judges the arms `calibration.json` names (`cohort-path`, `seed-twin`) and any other arm its summaries hold (the runner's sentinel); the null study judges the cohort and planted arms. `report` prints each arm's not-produced count and its seeds with a lineage contrast, and the total not produced.
- **Old arm names (I3).** `report` reads the stored seeds' `K+C`, `K+H`, `K+P` and `K+C+H+P` as `k-c`, `k-h`, `k-p` and `k-c-h-p`, so it runs on them without a patch.
- **The repeats refusal names the cohort-paired arms (I4).**
- **The lineage-only figures** in the section above now use Benjamini-Hochberg over the lineage contrasts alone (L3), and say what the stored seeds can't show about the sentinel (L4).
- **Tests:** `test_null_study_passes_only_with_a_lineage_contrast_in_20_seeds_per_arm` (the reviewer's 1-of-20 case, the 20/19 boundary, and a null study with one cohort arm in 19 seeds), `test_null_study_reads_summaries_written_before_the_mask` and `test_null_study_counts_a_contrast_not_produced_apart_and_refuses_an_arm_without_any` (L2: the reviewer's surviving mutants B2 and C1–C3), `test_null_study_reads_the_arm_names_from_before_the_rename`, and the arm names in `test_a_spec_with_repeats_and_a_cohort_arm_is_refused`.
- **After REVIEW-INPUT17 (main's rulings, 6 Nov 2026).**
  - **M1:** the first version of the gate judged only the arms it found, so an empty or mistyped folder printed `seeds 0; pass True` and a calibration arm absent from every seed passed. Now no seeds, a named arm in no seed, or any arm under 20 seeds fails the calibration, and the report names the arm. A folder that doesn't exist is refused, and the null study names a missing arm instead of raising `KeyError`.
  - **L1:** tests with two calibration arms, one of them thin, and with a planted arm thinned.
  - **L2:** seeds are counted, not files. Summaries are keyed by seed folder: a folder given twice, or the same seed in two folders, counts once when the files are byte-identical and is refused, naming both, when they differ.
  - **Tests:** `test_the_calibration_judges_every_arm_its_spec_names`, `test_every_judged_arm_needs_20_seeds_not_any` and `test_a_seed_counts_once_however_many_folders_hold_it`.
- **After REVIEW-INPUT17B (main's rulings, 6 Nov 2026).**
  - **L1:** `test_the_calibration_judges_the_runner_s_sentinel_too` runs `calibration.json`'s arms with the sentinel. A sentinel in 19 of 20 seeds fails the calibration, and the report names it.
  - **L2:** a seed is the one its `cohort.json` records (`one_seed` writes it beside `summary.json`), not its folder's name. So `seed-0101` and `seed-101.bak` beside `seed-101`, or a copy of seeds 101–110 renamed 111–120, count as the 10 seeds they are.
    - Copies of one seed count once when their summaries are byte-identical. When they differ, they are refused, naming both.
    - A summary whose `cohort.json` is missing or unreadable, has no `seed`, or has a seed that isn't an integer is refused.
    - Each of the 40 stored seeds (`cal1`, `cal2`, `nullA`, `nullB`) records the seed its folder names, so their reports are unchanged.
  - **Tests:** `test_a_seed_is_the_one_its_cohort_json_records` covers the reviewer's three cases, a 20-seed control, a differing copy and each unreadable seed. The test helpers now write `cohort.json` as `one_seed` does.
- **After REVIEW-INPUT18 (main's rulings, 7 Nov 2026).**
  - **L2:** `one_seed` wrote `summary.json` before `cohort.json`, and a rerun skips a seed folder that has a summary, so a seed killed between the two writes was never rerun and `report` refused it. It now writes `cohort.json` and `contrasts.json` first, then `summary.json` whole: a copy to `summary.json.part`, then a rename. `test_a_seed_killed_between_its_writes_reruns` kills `one_seed` after each of its writes in turn. Each time, a folder with a summary has its `cohort.json`, and a rerun leaves one seed that `report` reads.
  - **F1:** byte-identical summaries recorded under two seeds are refused, naming both files. Every run records its own `seconds`, so only a copy whose `cohort.json` was edited matches. The 40 stored summaries are byte-distinct, so their reports are unchanged.
  - **L3:** tests for a sentinel that improves, or worsens, on every lineage contrast; a thin first arm beside two full ones; a differing copy of the same length; and seeds 100 apart. The test helpers' summaries now carry `seconds`, as real ones do.
- **A multi-block fold matrix (REVIEW-TL1 I4).** A K-only arm can't show a leak in `assemble`'s cross-block coalescing, because `block_candidates` has already merged K's own training twins. `ablate::arm::tests::a_c_and_k_fold_matrix_collapses_cross_block_twins_on_training_samples_only` runs `run_arm`'s fold matrix on a C block and a K block.
  - Their twin columns agree on the training samples and differ on the outer test samples, and C has no row for one test sample.
  - It checks that the twins collapse into one column, with C's value where C has a row and K's where C has none.
  - It also checks that no change to the test samples changes the matrix's training rows: none of their bits set, all of them set, the twins' test bits swapped, or a C row for every sample.
  - A collapse keyed on every sample fails it, and so does assembling on every sample, ignoring availability, or letting K overwrite C.
  - **After REVIEW-INPUT18 M1:** that test calls `fold_matrix` itself, so a leak where `run_arm` hands it anything but the fold's training samples passed every test. `run_arm` now loads the data root and admits the blocks, then calls `score_arm`, which holds its scoring unchanged. In test builds, `fold_matrix` records each `train` it receives.
    - `run_arm_hands_each_c_and_k_fold_matrix_only_its_training_samples` runs `score_arm` on a C+K arm of 60 samples.
    - It checks that `fold_matrix` is called once per outer fold of each split, each time with exactly that fold's training samples from the plan, and that every run is scored.
    - It also checks that the toy can show the leak: in some fold of each split, a matrix on every sample keeps twins that the fold's own matrix collapses.
    - `run_arm` passing every sample, or the test samples, to `fold_matrix` fails it.

## H-1 support (H1CODE, 7 Nov 2026)

H1PLAN (`grm-h1plan-20261106`) found that only N. gonorrhoeae fits one xlarge and that Salmonella's 25 datasets don't fit one TTL. Main's rulings after H1PLAN (D1–D6) gave lane H1CODE gaps G1, G2(a), G6 and G10. This section describes them as built.

### Kset sharing (G1)

Every k31 root's genomes are pool genomes counted at the pool's k, so their ksets are the pool's, byte for byte (`test_every_root_counts_its_own_ksets` measured it).
- **Linking.** After the pool is counted, `Runner.root` hard-links the pool's finished kset (`<id>.kset`) of each of a k31 root's genomes into the root's kset folder and copies its sidecar (`<id>.json`, a few kilobytes), then runs `grm count`. `species::count_one` reuses a kset whose k matches and whose sidecar has the same catalog version and a digest, so the root counts nothing. A genome the root already holds a kset or sidecar for keeps its own. A link that fails stops the run.
- **Why the sidecar is a copy (REVIEW-H1CODE L1).** `grm count` writes a new kset through a `.part` file and rename, which replaces a link and never writes through it. But it rewrites a reused sidecar in place (`species.rs`, when its catalog version differs or it has no kset digest), and through a link that write would reach the pool and every other root.
- **Checked after the count.** After the count, every shared kset must still be the pool's inode and every sidecar copy must hold the pool's bytes, or the run stops: a recount or a sidecar rewritten for another catalog means the root counted otherwise than the pool. The root's `shared.json` records `genomes`, `shared` (ksets that are the pool's file), `own` and `linked_now`.
- **k21 and k25 roots** count their own, as before.
- **Disk.** `counted_roots` gives each root `shared`, the pending genomes whose kset it will link, and leaves those out of `pending`. They are every union genome while the pool is still to be counted (the pool's pending ksets count them), else the pool's finished ksets. The matrix margin (10%) and arm outputs (5%) are still taken on all of a root's genomes. `disk_bytes` counts a file with several hard links in the folder once.
- **`--no-share-ksets`** counts every root's own ksets, as before G1. It exists for the identity proof only.
- **Tests.** `test_the_k31_roots_share_the_pool_s_ksets_and_count_none_anew` (e2e: the 4 k31 roots' every kset is the pool's inode and every sidecar a byte-equal copy after the count, `own` 0; k21 and k25 hold their own) and `test_sharing_changes_no_output` (e2e: the whole comparison again with `--no-share-ksets`; predictions, plans and training traces byte-identical, J's records byte-identical, K's manifest and report equal after timings, and the shared WORK's roots smaller on disk). The measured-against-estimate test now takes the k31 roots' ksets as shared: estimate = pool + pending ksets (k21, k25) + 10% of all roots' ksets + 5% of them + the pool's `curated/` + 64 MiB, still measured ≤ estimate ≤ 1.5 × measured. Unit tests in `test_h1.py` (`KsetSharingTests`) cover the links, a kept own kset, a refused recount, a sidecar rewritten in place (refused, and the pool's files unchanged), `--no-share-ksets` and the disk check with and without the pool counted.
- **What it does for H-1** (H1PLAN's figures with G1, no k arms): Salmonella 183.5 GB, MTB 182.7, S. aureus 94.6, N. gonorrhoeae 77.7 and the pilot 120.6, each under its box's limit of 214–219 GB.

### `run.py combine` (G2(a))

`python3 -B run.py combine SPEC --work-from WORK [--work-from WORK ...] --master-manifest FILE --out DIR`

Splitting a species by arms would change its union and so its grouping, except for k arms, whose cohort is the default's. Splitting by datasets keeps the union, the grouping and each dataset's default arm. So a part holds some of SPEC's datasets and every arm, or every arm but the k arms, which a k add-on part runs (G2(b), below). `combine` then puts the parts back together.
- **What a run now keeps for it:** `WORK/spec.json` (the spec as canonical JSON, so its sha256 is K's `protocol_sha256`) and `WORK/raw-manifest.txt` (the RAW file list `inputs.json` hashes).
- **The pinned master manifest (REVIEW-H1CODE L6).** SPEC's `data_root.manifest_sha256` names it; `--master-manifest` must hash to it, and a spec without it refuses. `load_spec` refuses a `data_root` whose `manifest_sha256` isn't a sha256.
- **Refusals,** before any arm output is read:
  - a part's spec differs from SPEC in anything but `datasets` and `arms` (naming the keys), or K's manifest wasn't made from it;
  - a part's arms aren't SPEC's arms in SPEC's order, as SPEC states them, with only k arms left out; a k add-on part (it has `reference.json`) holds anything but k arms as SPEC states them;
  - a part lists a dataset SPEC doesn't (or with another partition), a dataset is in two main parts, or a dataset of SPEC is in none;
  - an arm on a dataset is in two parts, or an arm of SPEC on a dataset is in none;
  - a k add-on part takes its datasets from more than one main part, or its `reference.json` names another run than the main part that holds them;
  - the parts differ in the species pool (`species-groups.json` without its structure path, so two boxes' WORK paths may differ, and without its arm list, since a main part that leaves out the k arms lists fewer arms over the same union), the grm binary or the vocabulary (`inputs.json`), the source commit, `synthetic` or the gate (K's manifest), or eligible genomes and the panel floor;
  - K's `binary_sha256` isn't the grm the part's `inputs.json` hashed;
  - a part's `raw-manifest.txt` isn't the manifest `inputs.json` hashed, or any of its RAW files is missing from the master manifest (sha256sum lines) with that sha256;
  - `--out` is not new or empty, or a part lacks an arm's output on one of its datasets.
- **Each arm output is checked before it is copied (REVIEW-H1CODE L2)** against its part's `resume.json`: no folder on its path from the part's `arms/` (`arms/`, `arms/<arm>`, `arms/<arm>/<dataset>`) is a symlink (REVIEW-H1CODE-B L-a); the record names that arm and dataset, was made from the part's `inputs.json` and spec, records both run flags (a record without them refuses, naming the file: L-b), lists `runs.json`, `predictions.jsonl` and `arm.json`, and the folder holds exactly the files it lists, none a symlink, each with its recorded sha256. Every arm of every part must have run with the same `--trace-training` and `--run-identical-cohorts`.
- **Output.** It copies each part's `arms/<arm>/<dataset>/`, then runs J and K as a run does (`report_run`): one `stats.finalize`, so one BH family per species over every part's lineage contrasts; K's manifest and report with SPEC's `protocol_sha256`; and `summary.json`. `combine.json` and `summary.json`'s `combine` record the master manifest and its sha256, and each part's WORK, spec sha256, datasets, RAW manifest sha256 and species grouping sha256. K's manifest is a closed schema, so the part hashes are not in it.
- **Timings.** J's records are byte-identical to the unsplit run's. K's manifest and report differ only in timings: each run's `seconds`, `compute.actual_seconds`, the report's sha256 of the manifest (`inputs_sha256.manifest`), and the Markdown's seconds and input-hash lines. `tests/input_comparison/test_run.py` `timeless` is that comparison.
- **Tests.** `CombineEndToEndTests` (e2e): `synthetic_root.py --second-drug` adds a ciprofloxacin dataset (QRDR carriers resistant, flips from their own stream; the default root is byte for byte the same without it, `test_the_second_drug_adds_its_records_and_changes_nothing_else`). A comparison over both datasets (1,200 genomes, cap-600, k-c and k21 and the sentinels; H has no ciprofloxacin panel, so k-h would not be produced there) runs whole and as one-dataset parts, one of them with a k add-on (G2(b), below). Combining the parts equals the whole run after timings, its predictions and summary decisions included. Combining the whole run alone gives its own J and K (H1-PLAN §7 acceptance 11). `CombineRefusalTests` checks every refusal on fabricated parts, the master-manifest pin and the arm-output checks (a truncated predictions file, an extra file, a missing record, a link, a record of another arm, other inputs or another spec, and mixed flags) among them.

### The k add-on (`--reference-from`, G2(b))

`python3 -B run.py SPEC --raw RAW --work WORK --grm GRM --source-commit SHA --reference-from MAIN [...]`

Layth lifted register §13 G2 on 8 Nov 2026, so H-1 runs k21 and k25. A k arm's cohort is the default's, so it changes neither the union nor the grouping, and it is scored on the default arm's plans (`k31_groups`). So its runs can go on boxes of their own: a main part runs every arm but the k arms, and a k add-on runs only k arms against the main part's finished WORK, MAIN.
- **Checks, before anything is counted.** SPEC must hold k arms only, and each k arm's cohort must be the default's. MAIN must be a finished run, read as `combine` reads a part, and not an add-on itself. Then:
  - its spec may differ from SPEC only in `arms` and `datasets`, and K's manifest must have been made from it;
  - it ran each of SPEC's datasets with the same partition and none of SPEC's k arms;
  - it ran the same grm binary and vocabulary, source commit, `synthetic` flag, gate outcome, eligible genomes and panel floor;
  - every RAW file of this run is in its `raw-manifest.txt` with that size and sha256.
- **Copies.** For each dataset, the default arm's and both sentinels' outputs. Each is checked against MAIN's `resume.json` as `combine` checks an arm output, must have run with this run's flags, and is copied whole, `resume.json` too, then compared file for file. A linked destination refuses. The add-on writes `WORK/reference.json` (MAIN, the sha256s of its spec and `inputs.json`, the grouping and each copy's `resume.json` sha256), takes MAIN's `species-groups.json` as its own, and records MAIN's grouping sha256 in each k arm's `resume.json`. `summary.json`'s `reference_from` is `reference.json`.
- **Runs.** No species pool, no default arm and no sentinel: only the k arms' roots, counted at their own k, and the k arms. The disk check counts only those roots (`pool_bytes` 0). A rerun copies the default and sentinels again and resumes a finished k arm. With `--estimate-only`, `--reference-from` plans the k roots alone and MAIN is not read.
- **In `combine`.** An add-on part's k arms must have been scored on the main part's default plans (their `resume.json` `reference_plans`), and its default and sentinel copies must equal the main part's byte for byte. `combine` takes the k arms from the add-on and every other arm from the main part. `combine.json` records each part's arms and, for an add-on, the main part it ran against.
- **Tests.**
  - `ReferenceTests` (stub grm): an add-on takes the default and sentinels and runs only k21, with one `grm arm` call, count, structure and matrix at k 21 only, and no pool root. A rerun resumes k21 on a fresh copy, and `--estimate-only` plans the k21 root alone.
  - `ReferenceTests`' refusals: 20 references that differ each refuse before anything is counted. They are a cohort arm or no arm in SPEC, another margin, another partition, k21 already run, another grm, vocabulary, binary, commit, `synthetic` flag, gate, eligible count or floor, a K manifest not made from the spec, a RAW file MAIN didn't hash, other run flags, a reference that is itself an add-on, a damaged default, a missing sentinel and an unfinished run. A linked copy destination refuses too, with nothing written through it.
  - `CombineRefusalTests`: two main parts without k21 and k25 and three add-ons (k21 and k25 on one, one each on the other) reach the report in either order, and missing add-ons leave named gaps.
  - `CombineRefusalTests`' refusals: an arm in two parts, an add-on with a cohort arm or a restated, repeated or empty k arm list, an add-on over two main parts, a `reference.json` naming another run or no record, another pool, a changed, missing or linked copy, other reference plans and a damaged k arm. A main part that leaves out, reorders, restates or adds a non-k arm refuses too.
  - `CombineEndToEndTests` (e2e) now runs cap-600, k-c and k21: part1 runs every arm on meropenem, part2 every arm but k21 on ciprofloxacin, and `part2-k21`, a k add-on against part2, runs k21 there. Combining the three equals the whole run after timings. The add-on counted no pool, its copies equal part2's, and its k21 predictions are the whole run's byte for byte.

### Resume (G6)

A reused WORK skips an arm on a dataset when its `resume.json` still matches. That file is written after the arm finished and holds:
- the sha256 of `inputs.json`, of the spec (canonical) and of the species grouping;
- the reference plans (the default arm's `plans/*.json` for that dataset, which is all `grm arm --reference` reads; none for the default itself);
- `--trace-training` and `--run-identical-cohorts`;
- every output file with its sha256.

The folder must hold exactly the files listed, `runs.json`, `predictions.jsonl` and `arm.json` among them, none a symlink, each with its sha256. If the dataset folder, its arm folder or `arms/` itself is a symlink, the run refuses (REVIEW-H1CODE-B L-a): skipping would take another folder's output, and rerunning would write through the link. Anything else reruns the arm, after `resume.json` is removed and `plans/` cleared. Roots (`built.json`), the pool (`structure.json`) and curation were already reused. J and K always rerun. `summary.json` lists the resumed outputs (`resumed_arms`), and stderr names them. A WORK whose RAW, grm or vocabulary changed is still refused as a whole.
- **Tests.** `ResumeTests` (stub grm): a finished arm resumes; a changed spec, `--trace-training`, grouping, `inputs.json` or reference plans reruns it; `--run-identical-cohorts` toggled, a damaged, missing, extra or linked output, or a missing or unreadable `resume.json`, reruns it; a rerun whose `grm arm` fails leaves no `resume.json`; a linked dataset or arm folder refuses and nothing is written through it. `test_a_reused_work_resumes_every_finished_arm_and_reruns_a_damaged_one` (e2e): the finished WORK, copied, with one arm's predictions truncated, reruns that arm only (12 of 13 resumed) and gives the same J and K after timings.

### The FASTA check in `--estimate-only`

`--estimate-only` now also reports the union's genomes without a FASTA in RAW (`union_without_fasta`, the first 20 in `union_without_fasta_examples`) and says on stderr that the run would refuse at the pool. This is the same check as the pool's refusal. The exit code is still the disk's: H-1's gate before any box needs both rc 0 and `union_without_fasta` 0 (H1-PLAN §5).

### The revision-12 gate tests substance (main's ruling after H1DATA)

H-1's species check of the newly fetched genomes fails 16 Salmonella Typhi. Revision 12's Salmonella default keeps the same 2,000 genomes, but its legacy score reweights. So 456 genomes change position (drawn order `c60015cf…`, not `5864bd5d…`), and the `datasets.json` view moves (`2a64c19c…`, not `5519c61d…`), because 9 datasets' eligible-pool counts change. The order gate of REVIEW-INPUT8 would refuse every H-1 spec on that. Main ruled that the gate tests substance instead.
- **What it compares.** `run.substance(root, slugs)`, the rule of H1DATA's `tools/substantive_check.py`, for each species `genome_set` names (the run's own must be one):
  - `genome_set[<slug>]`: the sha256 of the default's `species/<slug>/genomes.txt` sorted, one id per line;
  - `datasets[<slug>][<id>]`: for each dataset revision 12 has, the sha256 of the canonical JSON of {its `datasets.json` entry without `n_resistant`, `n_susceptible`, `n_intermediate` and `selection`; the `curated/phenotypes.jsonl` row of every selected genome for its antibiotic}. The whole row fixes the learned label (`Phenotype::reported_label`) and what curate recorded with it.

  The default must give every pinned value, hold every pinned dataset and hold no other dataset of those species (main's ruling after REVIEW-H1DATA: an extra dataset would be learned and reported though revision 12 has none). A species without pinned datasets refuses. A duplicate phenotype row for one genome and drug refuses.
- **What it records and doesn't gate:** the drawn order (the run's own species from `curated.json`, the others from their untouched `genomes.txt`) and the `datasets.json` view. Each is compared with the spec's `recorded` values, and stderr says whether they are as recorded. The gate's record goes to `WORK/revision12-gate.json`, `summary.json`'s `revision12_gate` and `--estimate-only`'s output. K's manifest keeps `revision12_gate: passed`, since K's manifest schema is closed.
- **The pins.** `spec.json`'s gate holds revision 12's values for the four panel species (73 datasets: 25, 22, 19 and 7), computed from its master root with `run.substance` (lane H1CODE's `tools/pin_gate.py`). It also holds the sha256 of each source file (`reference`) and revision 12's own order and view as `recorded`. The H-1 specs carry the same pins, and record `c60015cf…` and `2a64c19c…` as "revision 12's cohort with the H-1 species check".
- **Checked on real data.** A default curated from the master root with its own species check, and one with H1DATA's staged `species_check.json`, both pass, by `run.substance` against the pins and by H1DATA's checker. The second records `c60015cf` and `2a64c19c`.
- **Tests.**
  - `SubstanceGateTests` (stub curate): revision 12 itself passes as recorded; the same set in another order passes and is recorded as not as the spec records it; one genome swapped refuses; one selected label changed refuses; an extra dataset refuses; a duplicate phenotype row refuses; changed eligible-pool counts and an unselected genome's row decide nothing.
  - `MultiSpeciesGateTests` (REVIEW-H1CODE M1): the stub curates a second panel species beside Salmonella, and the gate pins both; in a Salmonella run, the other species' changed genome set refuses, and so does a changed selected label in its dataset.
  - `test_the_revision12_gate_refuses_another_default_cohort`: a one-digit change of either pin, a missing pinned dataset, an unnamed own species, a missing `genomes.txt` and a species without pinned datasets, in both modes.
  - The e2e pins the synthetic default's substance (`SYNTHETIC_GATE`, from the real binary), and checks the gate record in WORK and `summary.json`.

### The H-1 specs (G10)

`scripts/input-comparison/h1/` holds the specs H-1 runs. Each is `spec.json` with only these changed:
- `decisions`: the H-1 rulings;
- `species`;
- `datasets`: revision 12's reported datasets, all confirmatory;
- `arms`: all of `spec.json`'s in the species specs (Layth lifted register §13 G2 on 8 Nov 2026), every arm but k21 and k25 in a main part, one k arm in a k add-on part;
- `pool_budget.default_kset_bytes`, with its source in `pool_budget.reason` (D2): revision 12's stage kset total ÷ 2,000, rounded up;
- `revision12_gate.recorded`: H-1's drawn order and view (above), with the H-1 root's species check `c2e7c389…` named in its note;
- `data_root` (main's ruling after REVIEW-H1DATA): the H-1 data root `/home/layth/review/grm-h1-data-root-20261107/` (hard links to the master root's files plus H1DATA's 9,394 FASTAs and species check, read only; the master root stays untouched) and its manifest, H1DATA's `MASTER-MANIFEST-h1.sha256`, pinned at `7dd00779…` for `combine`. The run takes it, or a box's copy pushed from that manifest, as `--raw`; it writes only under `--work`.

| File | Species | Datasets | Arms | kset bytes per genome (stage) |
|---|---|---:|---|---:|
| `pilot.json` | Salmonella | 1: `salmonella_enterica__ciprofloxacin` (D1) | grr-v1, k21, k-c-h-p, budget-16000 | 29,147,870 (R10) |
| `h1-salmonella_enterica.json` | Salmonella | 25 | 11 | 29,147,870 (R10) |
| `h1-salmonella_enterica-part1.json`, `-part2.json` | Salmonella | 13 and 12 (SAL-1, SAL-2: alternate ranks of revision 12's per-dataset time) | 9 | 29,147,870 |
| `h1-salmonella_enterica-part1-k21.json`, `-part1-k25.json`, `-part2-k21.json`, `-part2-k25.json` | Salmonella | part1's 13 or part2's 12 | 1 | 29,147,870 |
| `h1-mycobacterium_tuberculosis.json` | MTB | 22 | 11 | 25,881,860 (R7) |
| `h1-mycobacterium_tuberculosis-part1.json`, `-part1-k21.json`, `-part1-k25.json` | MTB | 22 | 9, 1 and 1 | 25,881,860 |
| `h1-staphylococcus_aureus.json` | S. aureus | 18 | 11 | 16,812,689 (R3) |
| `h1-staphylococcus_aureus-part1.json`, `-part1-k21.json`, `-part1-k25.json` | S. aureus | 18 | 9, 1 and 1 | 16,812,689 |
| `h1-neisseria_gonorrhoeae.json` | N. gonorrhoeae | 7 | 11 | 12,750,902 (R1) |
| `h1-neisseria_gonorrhoeae-part1.json`, `-part1-k21.json`, `-part1-k25.json` | N. gonorrhoeae | 7 | 9, 1 and 1 | 12,750,902 |

The 11 stage-1 arms are cap-1000, cap-4000, grr-v1, k21, k25, k-c, k-h, k-p, k-c-h-p, budget-1000 and budget-16000, and `run.py` adds the sentinels.

A species runs as its main parts, each with the nine arms that aren't k arms, and a k21 and a k25 add-on part per main part, each run with `--reference-from` that main part's WORK. `combine` of a species' full spec takes them all. Every part equals its species spec but for `datasets` and `arms`. The pilot runs k21 in its own run, with no add-on part, so it times a k21 root on real data.

The gate's pins, panel floor, base rates, curate inputs, learners, splits, repeats, margin and seed are `spec.json`'s; only the gate's `recorded` values are H-1's (above). `H1SpecTests` checks all of this, and feeds each species' parts to `combine`'s checks as finished parts. The lane report has each spec's `--estimate-only` on the master root.

## Cost

An arm costs a curate plus count, structure and matrix for each new data root (cap, selection and k arms; blocks and budget arms reuse the default root), and one `grm arm` per arm and dataset: two splits × five outer folds × five learners, with no final models, MIC regression or annotation. On the synthetic tree (2,000 genomes of 23 kb, about 30,000 K patterns, one dataset, Boat xlarge, three arms at a time) a root took 1 to 5 s and an arm 7 to 50 s; the whole comparison took 176 s. Real genomes are about 200 times longer and real matrices hold far more patterns, so these numbers don't scale to a real species; a real arm costs about what `grm learn`'s cross-validation costs for that dataset, plus a count, structure and matrix of the species for each cap, selection or k arm. H-1 sizes the real panel.
