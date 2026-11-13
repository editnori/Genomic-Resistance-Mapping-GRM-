# H-1 reading plan

I wrote this on 11 Nov 2026 at 06:10 UTC and tightened the wording at 06:15 UTC. The five main parts were still running. No part had written a report, and I had opened none of the new run's checkpoint files. Private working note. The ledger holds its sha256 (row L25).

The report at `31b7edb` prints as it is: BH inside each species run, rule 10 as written, every decision listed. I change no file of the run and no number in its report. This plan adds a second reading beside it and fixes that reading before the results exist.

## Seen so far

I had seen two things. The first is a Salmonella report over 25 datasets, 9 arms and the two twins, which I rebuilt from the arm outputs the earlier stages saved on 9 and 10 Nov (`rehearsal/pulled/keep/sal-c`). The genomes are real and public. A rebuilt report is diagnostic. The second is the N. gonorrhoeae counts and medians from the stage that stopped on 9 Nov, 7 datasets and the same arms without grr-v1 (`rehearsal/pulled/keep/ng-1`).

I had seen nothing for S. aureus or M. tuberculosis, since their rebuilt reports failed before they wrote a contrast. I had seen nothing for k21 or k25 on any species. So the 40 S. aureus and M. tuberculosis datasets are unseen for every arm, and k21 and k25 are unseen everywhere.

## Rule 10

Rule 10 marks every decision of a learner on a species unreliable when one twin contrast of that learner rejects or is missing. The real data shows two faults in it.

The first is the missing twin. The report lists a twin contrast as "not produced" on any dataset whose default has no lineage result. Salmonella has 7 such datasets of 25 and N. gonorrhoeae 3 of 7. So the rule flags every learner on both species, 383 contrasts each on Salmonella, before it looks at a twin.

The second is the threshold. The rule judges each twin contrast at 0.05 on its own, and an exactly calibrated test still rejects at least once in 18 datasets 60 times in 100.

On Salmonella the twins rejected 2 of 90 judged contrasts (2.2%, exact 95% interval 0.3% to 7.8%). A calibrated test gives 5%. The rule still flags all 1,915 contrasts.

## Added reading

**Calibration, per species and learner.** I count k, the lineage contrasts of the twin that judges the learner which reject (p < 0.05 and the 95% interval off zero), out of n produced. I flag the learner on the species when P(k or more | n, 0.05) < 0.05. A twin contrast that is missing because the default has no lineage result stays out of the count, since that dataset has no decision to protect. A twin contrast missing for any other reason counts as a rejection. I give the same count per learner across the four species.

**Family.** I read each decision twice: under BH inside its species (the code) and under BH across the four species together (the doc). A firm decision holds under both and sits on a learner the calibration count leaves unflagged.

**Size of the change.** Per arm I print the median and quartiles of the lineage log-loss change, the contrasts whose 90% interval lies inside the 0.02 margin, the contrasts that are exactly zero, and the median AUROC change. The twins go in the same table, and I read each arm against its twin.

**By dataset.** Per arm I count the datasets where the mean change over the five learners is above zero and below zero, with an exact sign test.

## Predictions

Salmonella showed one pattern. More genomes helped: cap-4000 was better on 14 of 18 datasets and cap-1000 worse on 9 of 10. Every feature-side arm sat with the twins.

P1 is the one prediction with a direction. On S. aureus and M. tuberculosis together, cap-4000 is better on more than half the datasets and cap-1000 is worse on more than half. I test each with a one-sided exact sign test at 0.05, ties dropped. A pass on both supports the pattern. A fail on either means it does not carry past Salmonella, and I report that.

- **P2, feature side.** For k-c, k-p, budget-1000 and budget-16000 I predict no direction. Two-sided sign test, and the share of contrasts inside the margin.
- **P3, k.** For k21 and k25 I predict no direction. Two-sided, all four species.
- **P4, selection.** For grr-v1 I predict no direction. Two-sided.

I report N. gonorrhoeae and Salmonella beside these and leave them out of the P1 test.

## Limits

- Every interval is conditional on these cohorts. Nothing here is external validation.
- A lower log loss can come from calibration alone, so the AUROC change prints beside it.
- The sign test treats datasets as separate, and datasets of one species share genomes. Its p guides and does not decide.
- A cohort arm scores on the test genomes both cohorts hold, drawn to the default's class balance. Its change is the effect of the training cohort on those genomes.
- The calibration count has little power on a small species. With 4 datasets it takes 2 rejections to flag a learner, and with 18 it takes 4.
- A retune alone moves the tree learners by more than the 0.02 margin (Salmonella, retuned twin, scm: quartiles -0.080 and +0.081). Few of their contrasts can sit inside the margin, whatever the arm does.

Script: `analysis/h1_contrasts.py`.
