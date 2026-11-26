# Staggered split by time

Not run yet. This page sets out the frozen plan of 25 Nov 2026 and the five rulings of 26 Nov that read it: R33, R35, R38, R40 and R41. No fold file, fit or score of this split exists.

## Question

A learner is fitted on the genomes public by one day, then used for months or years. Two questions, per learner:

1. Does a fit made at an earlier day score worse on a later window of genomes than a fit made just before that window, on the same test genomes?
2. Is that cost still there when both fits train on the same number of genomes?

This is a retrospective split inside one database, BV-BRC. It is never a prospective test and never external validation. It ranks no feature and names none.

The split of 14 Nov used one cut day per dataset. This one cuts each dataset at several days and compares fits made at each.

## Clocks

Each dataset is cut on two clocks, read apart:

- Public day: the day a genome entered BV-BRC. Every genome has one.
- Collection year: the four-digit year BV-BRC records for the sample's collection. A genome with no such year is in no training set and no test set on this clock.

Each clock gets its own design, tables and verdicts. No figure pools the two.

## Design

Example: shares of 40, 60 and 80% on a dataset of 500 dated genomes put the cut points at the days of genomes 200, 300 and 400 in date order.

- Cut point at share s of n dated genomes: the day (or year) of the genome at place ceil(s x n) in date order. Two shares that give one day leave the windows between them empty.
- Window j runs from after cut point j to cut point j + 1; the last window runs to the end.
- Test set of window j, from window 2 on: every genome of the window whose study has no genome of the species available by cut point j. A genome with no recorded study is tested nowhere.
- An undated genome is in no training set and no test set, and makes no study known.
- all(k): every genome dated by cut point k. recent(k): per class, the latest genomes of all(k), as many as all(1) holds. So recent(1) is all(1).
- Copies: a genome whose feature set is identical to a test genome's of window j trains nowhere on window j. A test genome of an earlier window can train a later cut point.
- A tie in date inside recent(k) goes to the smaller of the plan's keyed hashes, which name the clock, the dataset and the genome.
- One fold file per dataset, clock, cut point, arm and window, so both fits of a pair score the same test genomes. recent(1) and all(1) share one file.
- Each fit runs at the evaluator's default settings: seed 42, 4,000 screened features, 3 inner folds, five learners.

## Contrasts

| Contrast | Older fit | Fit at the window's own cut point | Scored on |
|---|---|---|---|
| Stale all | all(k) | all(j) | window j |
| Stale same size | recent(k) | recent(j) | window j |

The horizon is j - k. The difference is the older fit minus the fit made at the window's own cut point. A pair needs the test minimum of each class in its test set and at least 25 genomes of each class in both training sets. A dataset with several pairs at one horizon gives their mean. The main reading is stale all, horizon 1, on the public-day clock.

## Share sets

Before any fit, a count on metadata alone tries 4 share sets (40, 60, 80; 50, 70; 40, 60; 50, 75) and 3 test minimums per class (20, 15, 10): 12 options per clock.

Each clock gets one option. Among options with 12 datasets in 3 species at horizon 1 on stale all, the plan takes the most cut points, then the highest test minimum, then the most datasets, then the first in the order above. If no option reaches 12 datasets in 3 species, the plan takes the one with the most datasets, and that clock is descriptive only.

## Rule

Recall and log loss per learner, as on 14 Nov: cutoff 0.5, scores clipped to 0.0001 and 0.9999. Specificity and balanced accuracy stand beside recall and log loss and decide nothing.

Per learner, contrast, horizon and measure: the median over datasets, with a 95% interval from 10,000 draws of species with replacement, seed 20261125. An interval wholly on the worse side of zero reads worse, wholly on the other side better, and otherwise not shown.

A learner shows a cost when recall or log loss reads worse and neither reads better. A contrast then gets one verdict:

| Verdict | Condition |
|---|---|
| Cost | At least 3 of 5 learners show a cost and no learner reads better, on at least 12 datasets in 3 species |
| Mixed | Any learner reads better on recall or log loss |
| Not shown | Neither row above applies |
| Descriptive only | Under 12 datasets or 3 species |

Not shown means the interval includes zero. It never means the two fits score alike.

## Predictions

The plan wrote three before any fold file:

1. Stale same size, horizon 1, recall: not shown for at least 4 of 5 learners.
2. Stale all, horizon 1: a cost through log loss for at least 3 of 5.
3. Collection-year clock: not shown on every contrast.

The analysis marks each met, not met or not read, by a rule fixed in its code before any result. It reads the first two on the public-day clock. It reads the third on both contrasts at every horizon the collection-year clock has. The third is met only when horizon 1 has a verdict on both contrasts and every verdict at every horizon is not shown. A later horizon that is descriptive only, or has no dataset, gives no verdict.

Descriptive only is never met. A prediction is met only when both readings are made and both meet it. It is not met when a reading that was made fails it, and otherwise it is not read.

## Second reading

5 of the 24 datasets of the 14 Nov split list genomes with an intermediate call. They stay in the design with the release's labels. Each contrast is then read twice: first on every test genome, then with the listed genomes taken out of a pair's test genomes for both fits, with no refit.

The plan gives one verdict from one reading. The second reading, and the rule that a verdict stands only in both, are additions of ruling R33, made before any fold file.

- A pair enters the second reading only when its test set lists no genome, or at least 5 of each class, and each class keeps the test minimum once they are out.
- A dataset with no such pair at a horizon is left out of the second reading there and named with its reason. A pair left out gives no figure.
- A verdict stands only when both readings give it.
- The lists stay on the computer that makes the run.
- With lists given, the analysis stops when the list directory has no list file, when a list file names no dataset of the design, or when a list file names no test genome of its dataset. Without lists, no second reading is made.

Example, at a test minimum of 15: a test set of 20 resistant and 20 susceptible genomes that lists 5 of each is read again on 15 and 15. One that lists 3 resistant genomes and no susceptible one leaves the second reading.

## Checks

- Before any fit, the design at the single share 60%, with the 14 Nov rules, gives the cut day of each of the 14 Nov split's 24 datasets. The check stops unless its table holds those 24, each once.
- No path of the sealed later records is opened.
- Before any reading, the four checks of 14 Nov: each result names its fold file's sha256, both fits of a pair scored the same genomes with the same labels, training counts match the design, and recall worked out again matches the evaluator's.
- 101 planted faults, one a rule, each with the test that fails with it.

## Limits

- One database and one curation. Later genomes differ in study, place, laboratory and method at once.
- BV-BRC adds genomes in batches, so windows are uneven.
- The same-size arm trains on fewer genomes than a user would have, over a narrower period, so it can hold fewer studies.
- Species are few, so the interval is coarse.
- The record will hold rows per dataset, species, clock, contrast, horizon, arm and learner. It holds no genome id and names no study.

Record: frozen plan `FREEZE-STAGGER-20261125.md`, sha256 ceac8b5bf286fb5e94f05dad69d6ade622ad6234ee0a6fd0ac2469dfc182e39c; rulings R33, R35, R38, R40 and R41 of 26 Nov 2026; `scripts/time-split/stagger_design.py`, `scripts/time-split/stagger_analyze.py` and `scripts/time-split/stagger_faults.py`; the split of 14 Nov, [time-split.md](time-split.md).
