# First-public-day split

This page compares two ways of choosing training genomes. The past arm trains only on genomes public in BV-BRC by each dataset's cut day. The random arm trains on a random draw of the same size and class balance. Both are scored on the same later genomes of studies absent from training, in 24 datasets of 4 species at evaluation revision 14. This is a split by first-public day inside one database, BV-BRC. It is not a prospective test and not external validation.

Five learners fit each arm: scm, cart, rf, gbt and lasso. As the test genomes are later genomes of new studies, this is a retrospective test on newly deposited studies from the same database. Every figure on this page comes from the rerun after the review, which the section Review describes.

## Reading

In 24 datasets of 4 species, scored on later genomes of studies absent from training, training only on genomes public by each dataset's cut day read worse than a random draw of the same size on log loss for 5 of 5 learners and on balanced accuracy for 5 of 5, and on recall for 2 of 5 (rf and lasso; not shown for scm, cart and gbt); no learner read better on any measure or contrast.

Median differences over the 24 datasets, past minus random, in recall and log loss were -0.0306 and +0.2107 for scm, -0.0116 and +0.1727 for cart, -0.0173 and +0.0905 for rf, -0.0381 and +0.0921 for gbt, and -0.0311 and +0.0677 for lasso.

scm scores only 0 or 1, so its log loss is a count of wrong calls at the clip: 738 of 3,390 tested pairs in the past arm against 567 of 3,390 in the random arm.

The plan's rule, fixed before any fit, reads cost shown, with a cost for 5 of 5 learners; the interval resamples only the 4 species, so it is coarse.

## Plan and order

The plan was hashed with its analysis script (sha256 `d415dbd4`) before any fold file, fit or result of the experiment existed, by its hash list's own note. It fixed the design, the rule that admits a dataset, recall and log loss per learner on three contrasts, the species interval and the rule for a cost.

The plan named specificity and balanced accuracy as supporting figures and gave them no rule. The table layout script was hashed after the first result and, by its hash list's own note, before any score was read.

The rerun's plan (sha256 `5b886192`) was written after the results of 14 Nov and the review's counts had been read. It fixed what runs and what is compared, and adds no rule. The 14 Nov rule applies to the rerun as written.

| Event | Time (UTC) | Source |
|---|---|---|
| The plan of 14 Nov and its analysis script hashed | 2026-11-14 15:16:56Z | the plan's hash list |
| First result of the run of 14 Nov | 2026-11-14 15:43:54Z | the run log of 14 Nov |
| The table layout script hashed | 2026-11-14 15:47:19Z | the layout script's hash list |
| Last result of the run of 14 Nov | 2026-11-14 16:11:39Z | the run log of 14 Nov |
| The rerun's plan hashed | 2026-11-15 14:23:05Z | the rerun plan's hash list |
| First fit of the rerun started | 2026-11-15 14:23:53Z | the rerun's run log |
| Last result of the rerun | 2026-11-15 14:37:00Z | the rerun's run log |

**Exploratory.** Every figure the plan did not fix is exploratory and enters no rule: the share of log loss from wrong calls at the clip and balanced accuracy beside recall and log loss in the main table, which the rerun's plan added, the comparison of the study and time contrasts, the share of past genomes in each arm's training genomes as a reading of the study contrast, and the interval over every ordered species draw.

## Method

A dataset is one species and one antibiotic, and its cohort is the genomes `grm learn` evaluates for it. A genome's day is the day its record entered BV-BRC, so it is first public inside this database. A genome's study is its BV-BRC BioProject.

- **Cut day.** The earliest day by which at least 60% of the cohort was public. Past genomes have a day on or before it, later genomes a day after it.
- **Test studies.** A study is new when none of its genomes of the species is public by the cut day. A keyed hash of the accession picks half of the new studies with a later genome; the choice depends on no dataset and no label.
- **Test set.** A keyed hash picks two in three of the later genomes of the test studies. The rest are the reserve. All three arms are scored on exactly this test set.
- **Arms.** past trains on the past genomes. study-out draws from the past genomes and later genomes of other studies, with no genome of a test study. random draws from the study-out pool and the reserve, so it holds genomes of the test genomes' own studies and period.
- **Size.** Each arm trains on as many resistant and as many susceptible genomes as the past arm.
- **Left out.** A later genome with no recorded study is in no pool and no test set. A copy of a test genome, with identical inputs, is in no pool.
- **Entry.** A dataset enters when its test set holds at least 20 genomes of each class and each arm trains on at least 25 of each class.
- **Datasets.** 24 datasets in 4 species enter: 6 of *Escherichia coli*, 8 of *Klebsiella pneumoniae*, 2 of *Pseudomonas aeruginosa* and 8 of *Salmonella enterica*.
- **Fits.** Each arm is one fold file for `grm learn --fold-file` ([fold-file.md](fold-file.md)), with the default settings of `grm learn` (seed 42, 4,000 screened features, 3 inner folds). Each fit selects candidate features and tunes its settings inside its own training genomes.

Three contrasts compare the arms, each a paired difference on one dataset's test genomes. total is past minus random, the plan's question. study is study-out minus random: the test genomes' own studies absent from training against present. time is past minus study-out: later genomes of other studies absent against present. The plan gave the study and time contrasts no rule.

A score of 0.5 or more calls resistant. Log loss clips each score to [0.0001, 1 - 0.0001], so one wrong call at the clip adds 9.21 to a sum over the test genomes. scm scores only 0 or 1, so its log loss is a count of wrong calls at the clip.

Each contrast is read as its median over the 24 datasets, with a 95% interval: 10,000 draws of the 4 species with replacement (seed 20261114), each the median over the datasets of the drawn species, from the 251st to the 9,750th sorted draw.

4 species give only 35 distinct draws (256 in order), so the interval excludes zero only when nearly every mix of species agrees in sign. Over all 256 ordered draws the same two quantiles give the same reading in 60 of 60 readings.

A median reads worse when its interval lies wholly on the worse side of zero, better when it lies wholly on the other side, and not shown otherwise. Worse is lower recall, balanced accuracy or specificity, or a higher log loss.

For one learner the total contrast shows a cost when recall or log loss reads worse and neither reads better. The experiment shows a cost when at least 3 of 5 learners show one and none reads better, on at least 12 datasets in 3 species.

## Total contrast

Past minus random. Each cell gives the median over 24 datasets, the interval and the reading.

| Learner | Recall | Log loss | Balanced accuracy | Specificity | Rule |
|---|---|---|---|---|---|
| scm | -0.0306 [-0.0551, +0.0000] not shown | +0.2107 [+0.1226, +0.4016] worse | -0.0258 [-0.0622, -0.0042] worse | +0.0033 [+0.0000, +0.0385] not shown | cost |
| cart | -0.0116 [-0.0377, +0.0000] not shown | +0.1727 [+0.0689, +0.2972] worse | -0.0527 [-0.0704, -0.0151] worse | +0.0000 [-0.0197, +0.0085] not shown | cost |
| rf | -0.0173 [-0.1062, -0.0111] worse | +0.0905 [+0.0770, +0.1161] worse | -0.0608 [-0.0932, -0.0355] worse | +0.0109 [-0.0185, +0.0505] not shown | cost |
| gbt | -0.0381 [-0.1496, +0.0233] not shown | +0.0921 [+0.0447, +0.1702] worse | -0.0314 [-0.0748, -0.0037] worse | -0.0043 [-0.0419, +0.0000] not shown | cost |
| lasso | -0.0311 [-0.0648, -0.0111] worse | +0.0677 [+0.0345, +0.1494] worse | -0.0172 [-0.0440, -0.0047] worse | +0.0000 [+0.0000, +0.0000] not shown | cost |

## Study and time contrasts

Study is study-out minus random, and time is past minus study-out.

| Learner | Contrast | Recall | Log loss | Balanced accuracy | Specificity |
|---|---|---|---|---|---|
| scm | study | +0.0000 [-0.0278, +0.0000] not shown | +0.1185 [+0.0487, +0.3684] worse | -0.0156 [-0.0195, -0.0125] worse | +0.0000 [+0.0000, +0.0385] not shown |
| scm | time | +0.0000 [-0.0161, +0.0000] not shown | +0.0000 [+0.0000, +0.0501] not shown | +0.0000 [-0.0053, +0.0017] not shown | +0.0000 [+0.0000, +0.0000] not shown |
| cart | study | -0.0093 [-0.0417, +0.0000] not shown | +0.0726 [+0.0208, +0.1344] worse | -0.0098 [-0.0251, +0.0122] not shown | +0.0000 [-0.0192, +0.0035] not shown |
| cart | time | +0.0000 [-0.0538, +0.0000] not shown | +0.0506 [-0.0207, +0.2026] not shown | -0.0210 [-0.0642, +0.0000] not shown | +0.0000 [-0.1957, +0.0000] not shown |
| rf | study | -0.0155 [-0.0424, -0.0073] worse | +0.1035 [+0.0286, +0.1294] worse | -0.0096 [-0.0321, +0.0065] not shown | +0.0097 [+0.0000, +0.0283] not shown |
| rf | time | +0.0000 [-0.0233, +0.0000] not shown | -0.0030 [-0.0208, +0.0518] not shown | -0.0171 [-0.0276, -0.0110] worse | -0.0035 [-0.0185, -0.0035] worse |
| gbt | study | -0.0253 [-0.0331, +0.0000] not shown | +0.0880 [+0.0110, +0.1182] worse | -0.0144 [-0.0275, +0.0116] not shown | -0.0144 [-0.0354, +0.0000] not shown |
| gbt | time | +0.0000 [-0.0315, +0.0139] not shown | +0.0033 [-0.0121, +0.0515] not shown | -0.0145 [-0.0377, +0.0133] not shown | +0.0035 [-0.0741, +0.0346] not shown |
| lasso | study | -0.0392 [-0.0660, +0.0000] not shown | +0.0666 [+0.0150, +0.0797] worse | -0.0127 [-0.0478, +0.0019] not shown | +0.0000 [-0.0323, +0.0194] not shown |
| lasso | time | +0.0000 [-0.0446, +0.0278] not shown | +0.0163 [+0.0119, +0.0699] worse | -0.0071 [-0.0348, -0.0031] worse | +0.0000 [-0.0221, +0.0000] not shown |

Exploratory, in the review's narrower form: on log loss the study contrast's median is larger than the time contrast's for 5 of 5 learners. The time contrast reads worse on log loss for lasso alone, on balanced accuracy for rf and lasso and on specificity for rf, and not shown everywhere else; on 4 species, not shown is weak evidence either way.

The study contrast takes the test genomes' own studies out of training, and it also raises the share of past genomes in training from 76.19% in the random arm to 81.47% in the study-out arm. So it moves study membership and the mix of periods together. Later genomes also differ from past ones in place, laboratory and method, and the design separates none of these.

| Arm | Training genomes | Past | Later, other studies | Later, test studies | Share past |
|---|---|---|---|---|---|
| past | 20,968 | 20,968 | 0 | 0 | 100.00% |
| random | 20,968 | 15,976 | 3,747 | 1,245 | 76.19% |
| study-out | 20,968 | 17,082 | 3,886 | 0 | 81.47% |

## Per arm

Test genomes of the 24 datasets added up per learner and arm, as the table layout script adds them; a genome tested in several datasets of a species counts once per dataset. These sums are a description, and the plan's summary is the median over datasets above. The last column is exploratory: the share of the pooled log loss from wrong calls at the clip. scm's share is 100.0% at one decimal, as a right call at the clip adds about 0.0001 to the sum.

| Learner | Arm | Resistant called resistant | Recall | Susceptible called susceptible | Specificity | Wrong calls | Log loss from wrong calls at the clip |
|---|---|---|---|---|---|---|---|
| scm | past | 1,251 of 1,729 | 0.7235 | 1,401 of 1,661 | 0.8435 | 738 | 100.0% |
| scm | random | 1,419 of 1,729 | 0.8207 | 1,404 of 1,661 | 0.8453 | 567 | 100.0% |
| scm | study-out | 1,300 of 1,729 | 0.7519 | 1,377 of 1,661 | 0.8290 | 713 | 100.0% |
| cart | past | 1,238 of 1,729 | 0.7160 | 1,270 of 1,661 | 0.7646 | 882 | 52.4% |
| cart | random | 1,404 of 1,729 | 0.8120 | 1,381 of 1,661 | 0.8314 | 605 | 36.3% |
| cart | study-out | 1,321 of 1,729 | 0.7640 | 1,348 of 1,661 | 0.8116 | 721 | 44.9% |
| rf | past | 1,121 of 1,729 | 0.6484 | 1,452 of 1,661 | 0.8742 | 817 | 2.5% |
| rf | random | 1,344 of 1,729 | 0.7773 | 1,481 of 1,661 | 0.8916 | 565 | 7.2% |
| rf | study-out | 1,173 of 1,729 | 0.6784 | 1,503 of 1,661 | 0.9049 | 714 | 4.7% |
| gbt | past | 1,224 of 1,729 | 0.7079 | 1,444 of 1,661 | 0.8694 | 722 | 1.0% |
| gbt | random | 1,373 of 1,729 | 0.7941 | 1,471 of 1,661 | 0.8856 | 546 | 0.0% |
| gbt | study-out | 1,259 of 1,729 | 0.7282 | 1,467 of 1,661 | 0.8832 | 664 | 0.0% |
| lasso | past | 1,285 of 1,729 | 0.7432 | 1,447 of 1,661 | 0.8712 | 658 | 8.8% |
| lasso | random | 1,388 of 1,729 | 0.8028 | 1,449 of 1,661 | 0.8724 | 553 | 0.0% |
| lasso | study-out | 1,302 of 1,729 | 0.7530 | 1,446 of 1,661 | 0.8706 | 642 | 7.6% |

## Test sets

The test sets hold 3,390 genome and dataset pairs, 1,729 of them resistant, from 828 distinct genomes. A genome tested in several datasets of its species is one pair per dataset and one distinct genome. A test set holds 54 to 230 genomes (median 138.5), and each arm trains on 478 to 1,207 genomes (median 815.5), the same number in all three arms of a dataset.

The tested genomes come from 22 distinct studies. The per-species counts below add up to 26, as 4 studies have tested genomes in two species. A test set holds genomes of 1 to 6 studies (1 to 7 drawn), and 3 datasets hold one. In 6 datasets a drawn study has no tested genome, since all its later genomes went to the reserve.

| Species | Datasets | Tested pairs | Resistant pairs | Distinct tested genomes | Test studies | Genomes of the largest test study |
|---|---|---|---|---|---|---|
| *Escherichia coli* | 6 | 645 | 253 | 257 | 6 | 84 |
| *Klebsiella pneumoniae* | 8 | 1,148 | 736 | 269 | 10 | 107 |
| *Pseudomonas aeruginosa* | 2 | 136 | 63 | 72 | 6 | 26 |
| *Salmonella enterica* | 8 | 1,461 | 677 | 230 | 4 | 155 |

## Review

An independent review of the run of 14 Nov found that genomes a fold file lists as unused reached the fits by two paths. The support floor, the fewest training genomes a candidate feature must be present in and absent from, came from the whole cohort, and a copy group of unused genomes alone moved tuning onto grouped draws. With `grm learn --fold-file` fixed to take both from the tested and training genomes alone, the 72 fold files ran again unchanged.

47 of the 360 fits changed a score (3,224 of 50,850 scores), in 28 of the 72 result files and 20 of the 24 datasets. 6 of the 60 readings changed, none of them on the total contrast, and the rule's answer stayed cost shown. cart's and gbt's study recall, gbt's study specificity and cart's time log loss now read not shown; lasso's time log loss and balanced accuracy now read worse.

The review also found that the test sets hold genomes of 22 distinct studies where the first report said 26, which counted each species' studies apart, and 1 to 6 studies per dataset where it said 1 to 7. The tables on this page are the rerun's alone.

| Learner | Fits changed | Scores changed | Readings changed |
|---|---|---|---|
| scm | 10 of 72 | 136 of 10,170 | 0 of 12 |
| cart | 21 of 72 | 1,412 of 10,170 | 2 of 12 |
| rf | 0 of 72 | 0 of 10,170 | 0 of 12 |
| gbt | 8 of 72 | 838 of 10,170 | 2 of 12 |
| lasso | 8 of 72 | 838 of 10,170 | 2 of 12 |

## Limits

- One database, BV-BRC, with one curation of its labels and studies. Training and test genomes take their labels from the same kind of record.
- The day is the day a genome entered BV-BRC, which can differ from the day its sample became public in the sequence archives. A split on the archive day could order genomes differently.
- Later genomes differ from past ones in study, place, laboratory and method all at once. The design separates study membership from the rest and nothing finer, and the study contrast also changes the mix of periods.
- Past genomes with no recorded study train in every arm, and later ones are left out.
- One draw per arm, one cut per dataset and 4 species: read the direction as shown and the size as coarse.
- The interval resamples species only. It holds nothing of the draw, the seed or other species and datasets.
- Recall and specificity rest on one cutoff: a score of 0.5 or more calls resistant.
- scm's log loss is a count of wrong calls at the clip.
- The page gives no figure for a single dataset, and it names no feature and no study.

## Files

- `reference/time-split/` holds the record: 18 files that `scripts/time-split/record.py` writes from the input set. Its README lists each file and its columns, and `SHA256SUMS` gives each file's sha256.
- The record holds rows per dataset, species, arm and learner, and single counts. It holds no genome id, label or score per genome, no fold file and no list of test-set studies.
- `checks.tsv` gives 0 as the count of the input set's genome ids found in the record.
- The input set holds the result files, whose scores are per genome, so the repository does not hold it. `--check` alone derives every derived file from the row files.
- `tests/time_split` rebuilds the section Reading and every table of this page from the record, and finds every other number of the page in the record.

```sh
python3 -I scripts/time-split/record.py --check reference/time-split
python3 -I scripts/time-split/record.py --inputs <dir> --check reference/time-split
python3 -I -m unittest discover -s tests/time_split -t .
```
