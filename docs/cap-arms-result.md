# Cap arms result

The release fits each species on 2,000 genomes picked by its selection rule. This page asks whether a model fitted on the top 4,096 genomes a species by the same rule, `cap-all`, does better than the release's 2,000, `cap-2000`, on the default's held-out genomes, beside a twin of the default. It reads block C2 of the assumption plan (`docs/assumption-plan.md`), with the method of `docs/cap-arms.md`, on 8 species: *Enterococcus faecium*, *Escherichia coli*, *Klebsiella pneumoniae*, *Mycobacterium tuberculosis*, *Neisseria gonorrhoeae*, *Salmonella enterica*, *Staphylococcus aureus* and *Streptococcus pneumoniae*.

The reading's tables and its gates were committed before any C2 score was read (`docs/assumption-plan-corrections.md`, entries "C2 design", "C2 rulings" and "C2 review"). The run's copies were read once, on 23 Nov 2026, after `scripts/c2-cap/run_checks.py` printed 8 of 8 checks passed (`reference/cap-arms/checks.txt`), and the table script wrote `reference/cap-arms/` after its 14 of 14 gates passed (`reference/cap-arms/gates.txt`).

Every split here is a retrospective split inside one public database. A split by lineage or by study is never external validation. The margin of 0.02 is an adapted research screen. It is no applicable FDA validation criterion.

An added genome only trains and is never scored, so every arm is read on the default's held-out genomes. `cap-4000`, the top 4,000, was withdrawn before any fit and holds no result. Every arm, the default included, is fitted here by `grm arm`, so the default's figures here are not the release report's figures and are never set beside them as if they were.

## Reading

**Result.** `cap-all` does not meet conditions 1 and 2: 0 of the 5 learners carry the change, so the rule reads no difference shown for `cap-all` against `cap-2000`. It meets condition 3, the twin's, and condition 4, recall, for all 5 learners. Condition 5 is shown and decides nothing. Condition 6 is the maintainer's. `cap-4000` was withdrawn before any fit, and none of its conditions is read (`reference/cap-arms/page-reading.tsv`, `reference/cap-arms/page-conditions.tsv`).

"No difference shown" is the rule's reading. It is not a statement that two caps are equal.

Conditions 1 and 2 are read under the lead's ruling 15: a learner carries the change when it meets conditions 1 and 2 both, and at least 3 of the 5 learners must carry it. Condition 3 holds when the twin meets 1 in fewer than 3 of 5 learners and 2 in fewer than 3 of 5, each read apart. Condition 4 reads recall. `conditions.tsv` also writes condition 5, which reads the study split's median; the study split decides nothing in C2, so condition 5 is shown and not read. Condition 6 is the maintainer's. The stricter reading `conditions.tsv` writes beside decides nothing.

`cap-all`'s deciding conditions, as `reference/cap-arms/page-reading.tsv` holds them. Learners counts the learners that meet the condition.

| Arm | Condition | Reads | Met | Learners |
|---|---|---|---|---|
| `cap-all` | 1 and 2 | a learner carries when its lineage median difference's interval lies wholly above 0 and better beats worse by the one-sided sign test at p < 0.05 after Holm over 15; at least 3 of 5 carry | not met | 0 |
| `cap-all` | 3 | the twin meets 1 in fewer than 3 of 5 learners and 2 in fewer than 3 of 5, each read apart (learners: the larger of the twin's two counts) | met | 0 |
| `cap-all` | 4 | no learner's recall interval lies wholly below -0.02, for all five | met | 5 |
| `cap-all` | all deciding | conditions 1 and 2, 3, 4 all met | not met | none |

Every condition per arm and reading, as `reference/cap-arms/page-conditions.tsv` holds them.

| Arm | Reading | Condition | Reads | Met | Learners | Decides |
|---|---|---|---|---|---|---|
| `cap-4000` | ruling 15 | 1 and 2 | a learner carries when its lineage median difference's interval lies wholly above 0 and better beats worse by the one-sided sign test at p < 0.05 after Holm over 15; at least 3 of 5 carry | not read | none | withdrawn before any fit: not read |
| `cap-4000` | ruling 15 | 3 | the twin meets 1 in fewer than 3 of 5 learners and 2 in fewer than 3 of 5, each read apart (learners: the larger of the twin's two counts) | not read | none | withdrawn before any fit: not read |
| `cap-4000` | ruling 15 | 4 | no learner's recall interval lies wholly below -0.02, for all five | not read | none | withdrawn before any fit: not read |
| `cap-4000` | ruling 15 | 5 | each learner that carries has a study median of its lineage median's sign | not read | none | withdrawn before any fit: not read |
| `cap-4000` | ruling 15 | 6 | the maintainer rules the change | the maintainer's | none | withdrawn before any fit: not read |
| `cap-4000` | stricter | 1 | the lineage median difference's interval lies wholly above 0, for at least 3 of 5 learners | not read | none | withdrawn before any fit: not read |
| `cap-4000` | stricter | 2 | better beats worse by the one-sided sign test at p < 0.05 after Holm over 15, for all five | not read | none | withdrawn before any fit: not read |
| `cap-4000` | stricter | 3 | the twin meets neither 1 nor 2 as this reading reads them | not read | none | withdrawn before any fit: not read |
| `cap-4000` | stricter | 4 | no learner's recall interval lies wholly below -0.02, for all five | not read | none | withdrawn before any fit: not read |
| `cap-4000` | stricter | 5 | each of the five learners has a study median of its lineage median's sign | not read | none | withdrawn before any fit: not read |
| `cap-4000` | stricter | 6 | the maintainer rules the change | the maintainer's | none | withdrawn before any fit: not read |
| `cap-all` | ruling 15 | 1 and 2 | a learner carries when its lineage median difference's interval lies wholly above 0 and better beats worse by the one-sided sign test at p < 0.05 after Holm over 15; at least 3 of 5 carry | not met | 0 | decides |
| `cap-all` | ruling 15 | 3 | the twin meets 1 in fewer than 3 of 5 learners and 2 in fewer than 3 of 5, each read apart (learners: the larger of the twin's two counts) | met | 0 | decides |
| `cap-all` | ruling 15 | 4 | no learner's recall interval lies wholly below -0.02, for all five | met | 5 | decides |
| `cap-all` | ruling 15 | 5 | each learner that carries has a study median of its lineage median's sign | not met | 0 | decides nothing: it reads the study split |
| `cap-all` | ruling 15 | 6 | the maintainer rules the change | the maintainer's | none | the maintainer's |
| `cap-all` | stricter | 1 | the lineage median difference's interval lies wholly above 0, for at least 3 of 5 learners | not met | 0 | decides nothing: the stricter reading |
| `cap-all` | stricter | 2 | better beats worse by the one-sided sign test at p < 0.05 after Holm over 15, for all five | not met | 0 | decides nothing: the stricter reading |
| `cap-all` | stricter | 3 | the twin meets neither 1 nor 2 as this reading reads them | met | none | decides nothing: the stricter reading |
| `cap-all` | stricter | 4 | no learner's recall interval lies wholly below -0.02, for all five | met | 5 | decides nothing: the stricter reading |
| `cap-all` | stricter | 5 | each of the five learners has a study median of its lineage median's sign | not met | 2 | decides nothing: the stricter reading |
| `cap-all` | stricter | 6 | the maintainer rules the change | the maintainer's | none | decides nothing: the stricter reading |

The release keeps its cap of 2,000 genomes a species and its pipeline's rule whatever the count says, until the maintainer rules. Nothing changes in the release or on the site.

## Data

**Datasets.** The first run: the 143 CAP datasets of the 8 species that are not in M, 39 of them in S (`scripts/c2-cap/lists.tsv`). The datasets of M wait for the label rule. A dataset's genomes are its genomes in the release's lists.

**Arms.** `cap-2000` is the release's 2,000 genomes a species. `cap-all` is the top 4,096 by the release's selection rule on the rebuilt revision 14 table, or every genome the rule ranks when fewer. The twin is `cap-2000` again on slots of its own, scored as an arm of no change. `cap-4000` was withdrawn before any fit.

`cap-all` holds 4,096 genomes in 7 species and 2,907 in *Enterococcus faecium*, where the rule ranks 2,907. It adds 2,096 genomes in each of the 7 and 907 in *Enterococcus faecium*, and no added genome failed the fetch (`reference/cap-arms/page-species.tsv`).

**Ties at the cut.** Genomes of equal score sort by id, so where the 4,096 cut falls inside a tie of the selection rule, which of the tied genomes are added rests on their ids and not on their calls.

**None-added.** A dataset with no added genome labelled under `cap-all` is none-added: not fitted, and out of `cap-all`'s count, its denominator included, on both splits. The section Datasets names each. 27 of the 143 datasets are none-added. On the study split, 8 more datasets of S have added genomes labelled but none that trains, and are out of `cap-all`'s study rows (`reference/cap-arms/page-not-fitted.tsv`).

**Draws.** The default's and `cap-all`'s draws are `grm arm` slots 0 to 2 at evaluation seed 42, and the twin's are the default's slots 3 to 5. A dataset's value is the mean of its draws. A missing draw stops the scorer and is never read as zero.

**Unavailable lineage split.** A dataset whose lineage draws `grm arm` refuses past both replacements has no lineage test for any arm. It leaves the lineage count of every arm and the twin alike, and the section Datasets names it with its reason. Its study rows are scored as before. 58 of the 143 datasets have no lineage split, with 243 default slots refused; 18 of them are also none-added (`reference/cap-arms/checks.txt`, `reference/cap-arms/page-not-fitted.tsv`).

**Study split.** The datasets in S get one study fit a dataset and arm, `grm learn --fold-file` at seed 42 on the default's `study-whole` test folds. An added genome trains there in every fold or in none, and only when its BioProject is known and no tested genome has it or its k-mer digest. `cap-all` has a study fit on 27 of the 39 datasets in S, and its study fits train on 1,978 added genomes, counted once a dataset (`reference/cap-arms/page-study.tsv`).

Per species, as `reference/cap-arms/page-species.tsv` holds it. The pairs are a dataset and a learner on the lineage split, out of those scored, and decide nothing.

| Species | Datasets | In S | `cap-2000` genomes | `cap-all` genomes | Added | Short of 4,096 | None-added | `cap-all` better / worse / same | Twin better / worse / same |
|---|---|---|---|---|---|---|---|---|---|
| *Enterococcus faecium* | 11 | 0 | 2,000 | 2,907 | 907 | the rule ranks fewer genomes | 8 | 0 / 0 / 10 of 10 | 14 / 4 / 22 of 40 |
| *Escherichia coli* | 29 | 15 | 2,000 | 4,096 | 2,096 | none | 8 | 3 / 5 / 42 of 50 | 3 / 9 / 43 of 55 |
| *Klebsiella pneumoniae* | 32 | 18 | 2,000 | 4,096 | 2,096 | none | 4 | 10 / 26 / 59 of 95 | 17 / 14 / 64 of 95 |
| *Mycobacterium tuberculosis* | 20 | 1 | 2,000 | 4,096 | 2,096 | none | 1 | 17 / 15 / 58 of 90 | 3 / 10 / 77 of 90 |
| *Neisseria gonorrhoeae* | 6 | 0 | 2,000 | 4,096 | 2,096 | none | 1 | 4 / 4 / 7 of 15 | 4 / 4 / 7 of 15 |
| *Salmonella enterica* | 24 | 1 | 2,000 | 4,096 | 2,096 | none | 3 | 35 / 6 / 49 of 90 | 13 / 22 / 55 of 90 |
| *Staphylococcus aureus* | 13 | 4 | 2,000 | 4,096 | 2,096 | none | 0 | 10 / 0 / 5 of 15 | 2 / 1 / 12 of 15 |
| *Streptococcus pneumoniae* | 8 | 0 | 2,000 | 4,096 | 2,096 | none | 2 | 4 / 2 / 9 of 15 | 3 / 3 / 19 of 25 |

## Count

The primary measure is the pooled lineage balanced accuracy on the default's outer test folds, per learner. A dataset is better, worse or the same by more than 0.02, compared at 12 decimals. The family is `cap-4000`, `cap-all` and the twin against `cap-2000`, Holm over 15 lineage pairs, 3 arms by 5 learners. `cap-4000`'s five pairs stay in the family at p = 1, named withdrawn.

The datasets that add fewer than 20 labelled genomes under `cap-all` stay in the count: the count was fixed with no threshold before the cohorts were listed (the lead's ruling of 12:19Z on 23 Nov 2026). 19 datasets add fewer than 20 (`reference/cap-arms/page-few-added.tsv`).

`cap-all` enters 76 of its 116 fitted datasets on the lineage split, and the twin 85 of 143. No learner's `cap-all` pair meets condition 1 or condition 2. The smallest Holm p is 0.479, for `scm`, and `cart`'s interval starts at 0.000, not above it. The family, as `reference/cap-arms/page-count.tsv` holds its lineage balanced accuracy rows:

| Pair | Learner | Condition 1 | Condition 2 | Better | Worse | Same | Median | Interval | Sign p | Holm p |
|---|---|---|---|---|---|---|---|---|---|---|
| cap-4000 less cap-2000 | `scm` | withdrawn | withdrawn | none | none | none | none | none | 1.000 | 1.000 |
| cap-4000 less cap-2000 | `cart` | withdrawn | withdrawn | none | none | none | none | none | 1.000 | 1.000 |
| cap-4000 less cap-2000 | `rf` | withdrawn | withdrawn | none | none | none | none | none | 1.000 | 1.000 |
| cap-4000 less cap-2000 | `gbt` | withdrawn | withdrawn | none | none | none | none | none | 1.000 | 1.000 |
| cap-4000 less cap-2000 | `lasso` | withdrawn | withdrawn | none | none | none | none | none | 1.000 | 1.000 |
| cap-all less cap-2000 | `scm` | not met | not met | 17 of 76 | 7 of 76 | 52 of 76 | 0.002 | -0.003 to 0.005 | 0.032 | 0.479 |
| cap-all less cap-2000 | `cart` | not met | not met | 17 of 76 | 13 of 76 | 46 of 76 | 0.001 | 0.000 to 0.007 | 0.292 | 1.000 |
| cap-all less cap-2000 | `rf` | not met | not met | 20 of 76 | 10 of 76 | 46 of 76 | 0.002 | -0.002 to 0.017 | 0.049 | 0.691 |
| cap-all less cap-2000 | `gbt` | not met | not met | 13 of 76 | 14 of 76 | 49 of 76 | -0.002 | -0.005 to 0.008 | 0.649 | 1.000 |
| cap-all less cap-2000 | `lasso` | not met | not met | 16 of 76 | 14 of 76 | 46 of 76 | -0.002 | -0.005 to 0.004 | 0.428 | 1.000 |
| twin less cap-2000 | `scm` | not met | not met | 12 of 85 | 13 of 85 | 60 of 85 | 0.000 | -0.001 to 0.001 | 0.655 | 1.000 |
| twin less cap-2000 | `cart` | not met | not met | 13 of 85 | 15 of 85 | 57 of 85 | 0.000 | -0.002 to 0.001 | 0.714 | 1.000 |
| twin less cap-2000 | `rf` | not met | not met | 8 of 85 | 11 of 85 | 66 of 85 | -0.001 | -0.003 to 0.001 | 0.820 | 1.000 |
| twin less cap-2000 | `gbt` | not met | not met | 13 of 85 | 13 of 85 | 59 of 85 | -0.001 | -0.003 to 0.001 | 0.577 | 1.000 |
| twin less cap-2000 | `lasso` | not met | not met | 13 of 85 | 15 of 85 | 57 of 85 | 0.000 | -0.001 to 0.001 | 0.714 | 1.000 |

Every row of the count, each measure on both splits, as `reference/cap-arms/page-count.tsv` holds it. The study rows decide nothing.

| Split | Pair | Learner | Measure | Better is | Entered | Better | Worse | Same | Median | Interval | Species | Reading |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| lineage | cap-all less cap-2000 | `scm` | balanced accuracy | higher | 76 of 116 | 17 | 7 | 52 | 0.002 | -0.003 to 0.005 | 8 | the count |
| lineage | cap-all less cap-2000 | `scm` | recall | higher | 76 of 116 | 26 | 18 | 32 | 0.005 | -0.001 to 0.009 | 8 | the count |
| lineage | cap-all less cap-2000 | `scm` | specificity | higher | 76 of 116 | 25 | 15 | 36 | 0.000 | -0.003 to 0.007 | 8 | the count |
| lineage | cap-all less cap-2000 | `scm` | very major error rate | lower | 76 of 116 | 26 | 18 | 32 | -0.005 | -0.009 to 0.001 | 8 | the count |
| lineage | cap-all less cap-2000 | `scm` | major error rate | lower | 76 of 116 | 25 | 15 | 36 | 0.000 | -0.007 to 0.003 | 8 | the count |
| lineage | cap-all less cap-2000 | `scm` | log loss | lower | 76 of 116 | 43 | 25 | 8 | -0.026 | -0.054 to 0.027 | 8 | the count |
| lineage | cap-all less cap-2000 | `cart` | balanced accuracy | higher | 76 of 116 | 17 | 13 | 46 | 0.001 | 0.000 to 0.007 | 8 | the count |
| lineage | cap-all less cap-2000 | `cart` | recall | higher | 76 of 116 | 30 | 21 | 25 | 0.002 | -0.003 to 0.017 | 8 | the count |
| lineage | cap-all less cap-2000 | `cart` | specificity | higher | 76 of 116 | 15 | 14 | 47 | 0.000 | -0.001 to 0.002 | 8 | the count |
| lineage | cap-all less cap-2000 | `cart` | very major error rate | lower | 76 of 116 | 30 | 21 | 25 | -0.002 | -0.017 to 0.003 | 8 | the count |
| lineage | cap-all less cap-2000 | `cart` | major error rate | lower | 76 of 116 | 15 | 14 | 47 | 0.000 | -0.002 to 0.001 | 8 | the count |
| lineage | cap-all less cap-2000 | `cart` | log loss | lower | 76 of 116 | 31 | 21 | 24 | -0.010 | -0.018 to 0.001 | 8 | the count |
| lineage | cap-all less cap-2000 | `rf` | balanced accuracy | higher | 76 of 116 | 20 | 10 | 46 | 0.002 | -0.002 to 0.017 | 8 | the count |
| lineage | cap-all less cap-2000 | `rf` | recall | higher | 76 of 116 | 16 | 23 | 37 | -0.002 | -0.014 to 0.003 | 8 | the count |
| lineage | cap-all less cap-2000 | `rf` | specificity | higher | 76 of 116 | 31 | 10 | 35 | 0.009 | -0.001 to 0.034 | 8 | the count |
| lineage | cap-all less cap-2000 | `rf` | very major error rate | lower | 76 of 116 | 16 | 23 | 37 | 0.002 | -0.003 to 0.014 | 8 | the count |
| lineage | cap-all less cap-2000 | `rf` | major error rate | lower | 76 of 116 | 31 | 10 | 35 | -0.009 | -0.033 to 0.000 | 8 | the count |
| lineage | cap-all less cap-2000 | `rf` | log loss | lower | 76 of 116 | 18 | 22 | 36 | 0.001 | -0.009 to 0.014 | 8 | the count |
| lineage | cap-all less cap-2000 | `gbt` | balanced accuracy | higher | 76 of 116 | 13 | 14 | 49 | -0.002 | -0.005 to 0.008 | 8 | the count |
| lineage | cap-all less cap-2000 | `gbt` | recall | higher | 76 of 116 | 14 | 30 | 32 | -0.002 | -0.015 to 0.001 | 8 | the count |
| lineage | cap-all less cap-2000 | `gbt` | specificity | higher | 76 of 116 | 22 | 7 | 47 | 0.005 | 0.001 to 0.011 | 8 | the count |
| lineage | cap-all less cap-2000 | `gbt` | very major error rate | lower | 76 of 116 | 14 | 30 | 32 | 0.002 | -0.001 to 0.015 | 8 | the count |
| lineage | cap-all less cap-2000 | `gbt` | major error rate | lower | 76 of 116 | 22 | 7 | 47 | -0.005 | -0.011 to -0.001 | 8 | the count |
| lineage | cap-all less cap-2000 | `gbt` | log loss | lower | 76 of 116 | 19 | 16 | 41 | 0.001 | -0.013 to 0.012 | 8 | the count |
| lineage | cap-all less cap-2000 | `lasso` | balanced accuracy | higher | 76 of 116 | 16 | 14 | 46 | -0.002 | -0.005 to 0.004 | 8 | the count |
| lineage | cap-all less cap-2000 | `lasso` | recall | higher | 76 of 116 | 17 | 20 | 39 | -0.001 | -0.011 to 0.005 | 8 | the count |
| lineage | cap-all less cap-2000 | `lasso` | specificity | higher | 76 of 116 | 26 | 14 | 36 | 0.003 | 0.000 to 0.008 | 8 | the count |
| lineage | cap-all less cap-2000 | `lasso` | very major error rate | lower | 76 of 116 | 17 | 20 | 39 | 0.001 | -0.005 to 0.011 | 8 | the count |
| lineage | cap-all less cap-2000 | `lasso` | major error rate | lower | 76 of 116 | 26 | 14 | 36 | -0.003 | -0.008 to 0.000 | 8 | the count |
| lineage | cap-all less cap-2000 | `lasso` | log loss | lower | 76 of 116 | 22 | 18 | 36 | 0.002 | -0.015 to 0.011 | 8 | the count |
| lineage | twin less cap-2000 | `scm` | balanced accuracy | higher | 85 of 143 | 12 | 13 | 60 | 0.000 | -0.001 to 0.001 | 8 | the count |
| lineage | twin less cap-2000 | `scm` | recall | higher | 85 of 143 | 21 | 26 | 38 | -0.005 | -0.013 to 0.000 | 8 | the count |
| lineage | twin less cap-2000 | `scm` | specificity | higher | 85 of 143 | 18 | 17 | 50 | 0.000 | -0.003 to 0.004 | 8 | the count |
| lineage | twin less cap-2000 | `scm` | very major error rate | lower | 85 of 143 | 21 | 26 | 38 | 0.005 | 0.000 to 0.014 | 8 | the count |
| lineage | twin less cap-2000 | `scm` | major error rate | lower | 85 of 143 | 18 | 17 | 50 | 0.000 | -0.004 to 0.003 | 8 | the count |
| lineage | twin less cap-2000 | `scm` | log loss | lower | 85 of 143 | 30 | 39 | 16 | 0.005 | -0.007 to 0.034 | 8 | the count |
| lineage | twin less cap-2000 | `cart` | balanced accuracy | higher | 85 of 143 | 13 | 15 | 57 | 0.000 | -0.002 to 0.001 | 8 | the count |
| lineage | twin less cap-2000 | `cart` | recall | higher | 85 of 143 | 17 | 19 | 49 | 0.000 | -0.002 to 0.001 | 8 | the count |
| lineage | twin less cap-2000 | `cart` | specificity | higher | 85 of 143 | 13 | 16 | 56 | 0.000 | -0.001 to 0.001 | 8 | the count |
| lineage | twin less cap-2000 | `cart` | very major error rate | lower | 85 of 143 | 17 | 19 | 49 | 0.000 | 0.000 to 0.002 | 8 | the count |
| lineage | twin less cap-2000 | `cart` | major error rate | lower | 85 of 143 | 13 | 16 | 56 | 0.000 | -0.001 to 0.001 | 8 | the count |
| lineage | twin less cap-2000 | `cart` | log loss | lower | 85 of 143 | 24 | 34 | 27 | 0.004 | 0.000 to 0.006 | 8 | the count |
| lineage | twin less cap-2000 | `rf` | balanced accuracy | higher | 85 of 143 | 8 | 11 | 66 | -0.001 | -0.003 to 0.001 | 8 | the count |
| lineage | twin less cap-2000 | `rf` | recall | higher | 85 of 143 | 16 | 12 | 57 | 0.001 | -0.003 to 0.002 | 8 | the count |
| lineage | twin less cap-2000 | `rf` | specificity | higher | 85 of 143 | 9 | 14 | 62 | 0.000 | -0.005 to 0.004 | 8 | the count |
| lineage | twin less cap-2000 | `rf` | very major error rate | lower | 85 of 143 | 16 | 12 | 57 | -0.001 | -0.002 to 0.003 | 8 | the count |
| lineage | twin less cap-2000 | `rf` | major error rate | lower | 85 of 143 | 9 | 14 | 62 | 0.000 | -0.004 to 0.005 | 8 | the count |
| lineage | twin less cap-2000 | `rf` | log loss | lower | 85 of 143 | 5 | 11 | 69 | 0.003 | 0.000 to 0.007 | 8 | the count |
| lineage | twin less cap-2000 | `gbt` | balanced accuracy | higher | 85 of 143 | 13 | 13 | 59 | -0.001 | -0.003 to 0.001 | 8 | the count |
| lineage | twin less cap-2000 | `gbt` | recall | higher | 85 of 143 | 25 | 20 | 40 | 0.001 | -0.001 to 0.003 | 8 | the count |
| lineage | twin less cap-2000 | `gbt` | specificity | higher | 85 of 143 | 13 | 24 | 48 | -0.002 | -0.010 to -0.001 | 8 | the count |
| lineage | twin less cap-2000 | `gbt` | very major error rate | lower | 85 of 143 | 25 | 20 | 40 | -0.001 | -0.003 to 0.001 | 8 | the count |
| lineage | twin less cap-2000 | `gbt` | major error rate | lower | 85 of 143 | 13 | 24 | 48 | 0.002 | 0.001 to 0.010 | 8 | the count |
| lineage | twin less cap-2000 | `gbt` | log loss | lower | 85 of 143 | 14 | 21 | 50 | 0.000 | -0.003 to 0.008 | 8 | the count |
| lineage | twin less cap-2000 | `lasso` | balanced accuracy | higher | 85 of 143 | 13 | 15 | 57 | 0.000 | -0.001 to 0.001 | 8 | the count |
| lineage | twin less cap-2000 | `lasso` | recall | higher | 85 of 143 | 26 | 16 | 43 | 0.001 | -0.004 to 0.012 | 8 | the count |
| lineage | twin less cap-2000 | `lasso` | specificity | higher | 85 of 143 | 11 | 25 | 49 | -0.001 | -0.007 to 0.002 | 8 | the count |
| lineage | twin less cap-2000 | `lasso` | very major error rate | lower | 85 of 143 | 26 | 16 | 43 | -0.001 | -0.012 to 0.004 | 8 | the count |
| lineage | twin less cap-2000 | `lasso` | major error rate | lower | 85 of 143 | 11 | 25 | 49 | 0.001 | -0.002 to 0.007 | 8 | the count |
| lineage | twin less cap-2000 | `lasso` | log loss | lower | 85 of 143 | 23 | 26 | 36 | 0.000 | -0.007 to 0.012 | 8 | the count |
| study | cap-all less cap-2000 | `scm` | balanced accuracy | higher | 27 of 27 | 8 | 3 | 16 | 0.000 | -0.004 to 0.008 | 4 | the study split: decides nothing |
| study | cap-all less cap-2000 | `scm` | recall | higher | 27 of 27 | 5 | 6 | 16 | -0.004 | -0.009 to 0.000 | 4 | the study split: decides nothing |
| study | cap-all less cap-2000 | `scm` | specificity | higher | 27 of 27 | 8 | 3 | 16 | 0.003 | -0.003 to 0.010 | 4 | the study split: decides nothing |
| study | cap-all less cap-2000 | `scm` | very major error rate | lower | 27 of 27 | 5 | 6 | 16 | 0.004 | 0.000 to 0.009 | 4 | the study split: decides nothing |
| study | cap-all less cap-2000 | `scm` | major error rate | lower | 27 of 27 | 8 | 3 | 16 | -0.003 | -0.010 to 0.003 | 4 | the study split: decides nothing |
| study | cap-all less cap-2000 | `scm` | log loss | lower | 27 of 27 | 10 | 11 | 6 | 0.000 | -0.032 to 0.071 | 4 | the study split: decides nothing |
| study | cap-all less cap-2000 | `cart` | balanced accuracy | higher | 27 of 27 | 3 | 8 | 16 | -0.009 | -0.013 to 0.000 | 4 | the study split: decides nothing |
| study | cap-all less cap-2000 | `cart` | recall | higher | 27 of 27 | 5 | 9 | 13 | 0.000 | -0.007 to 0.000 | 4 | the study split: decides nothing |
| study | cap-all less cap-2000 | `cart` | specificity | higher | 27 of 27 | 3 | 10 | 14 | 0.000 | -0.008 to 0.001 | 4 | the study split: decides nothing |
| study | cap-all less cap-2000 | `cart` | very major error rate | lower | 27 of 27 | 5 | 9 | 13 | 0.000 | 0.000 to 0.007 | 4 | the study split: decides nothing |
| study | cap-all less cap-2000 | `cart` | major error rate | lower | 27 of 27 | 3 | 10 | 14 | 0.000 | -0.001 to 0.008 | 4 | the study split: decides nothing |
| study | cap-all less cap-2000 | `cart` | log loss | lower | 27 of 27 | 10 | 10 | 7 | 0.001 | -0.018 to 0.005 | 4 | the study split: decides nothing |
| study | cap-all less cap-2000 | `rf` | balanced accuracy | higher | 27 of 27 | 2 | 4 | 21 | 0.000 | -0.005 to 0.008 | 4 | the study split: decides nothing |
| study | cap-all less cap-2000 | `rf` | recall | higher | 27 of 27 | 2 | 2 | 23 | -0.002 | -0.005 to 0.000 | 4 | the study split: decides nothing |
| study | cap-all less cap-2000 | `rf` | specificity | higher | 27 of 27 | 5 | 4 | 18 | 0.001 | 0.000 to 0.014 | 4 | the study split: decides nothing |
| study | cap-all less cap-2000 | `rf` | very major error rate | lower | 27 of 27 | 2 | 2 | 23 | 0.002 | 0.000 to 0.005 | 4 | the study split: decides nothing |
| study | cap-all less cap-2000 | `rf` | major error rate | lower | 27 of 27 | 5 | 4 | 18 | -0.001 | -0.014 to 0.000 | 4 | the study split: decides nothing |
| study | cap-all less cap-2000 | `rf` | log loss | lower | 27 of 27 | 3 | 6 | 18 | -0.008 | -0.010 to 0.005 | 4 | the study split: decides nothing |
| study | cap-all less cap-2000 | `gbt` | balanced accuracy | higher | 27 of 27 | 2 | 4 | 21 | -0.001 | -0.007 to 0.001 | 4 | the study split: decides nothing |
| study | cap-all less cap-2000 | `gbt` | recall | higher | 27 of 27 | 9 | 6 | 12 | -0.003 | -0.004 to 0.014 | 4 | the study split: decides nothing |
| study | cap-all less cap-2000 | `gbt` | specificity | higher | 27 of 27 | 7 | 5 | 15 | 0.000 | -0.012 to 0.008 | 4 | the study split: decides nothing |
| study | cap-all less cap-2000 | `gbt` | very major error rate | lower | 27 of 27 | 9 | 6 | 12 | 0.003 | -0.014 to 0.004 | 4 | the study split: decides nothing |
| study | cap-all less cap-2000 | `gbt` | major error rate | lower | 27 of 27 | 7 | 5 | 15 | 0.000 | -0.008 to 0.013 | 4 | the study split: decides nothing |
| study | cap-all less cap-2000 | `gbt` | log loss | lower | 27 of 27 | 5 | 12 | 10 | 0.006 | -0.013 to 0.031 | 4 | the study split: decides nothing |
| study | cap-all less cap-2000 | `lasso` | balanced accuracy | higher | 27 of 27 | 4 | 3 | 20 | 0.002 | -0.006 to 0.003 | 4 | the study split: decides nothing |
| study | cap-all less cap-2000 | `lasso` | recall | higher | 27 of 27 | 4 | 7 | 16 | -0.004 | -0.009 to -0.003 | 4 | the study split: decides nothing |
| study | cap-all less cap-2000 | `lasso` | specificity | higher | 27 of 27 | 8 | 6 | 13 | 0.000 | -0.007 to 0.006 | 4 | the study split: decides nothing |
| study | cap-all less cap-2000 | `lasso` | very major error rate | lower | 27 of 27 | 4 | 7 | 16 | 0.004 | 0.003 to 0.009 | 4 | the study split: decides nothing |
| study | cap-all less cap-2000 | `lasso` | major error rate | lower | 27 of 27 | 8 | 6 | 13 | 0.000 | -0.006 to 0.007 | 4 | the study split: decides nothing |
| study | cap-all less cap-2000 | `lasso` | log loss | lower | 27 of 27 | 11 | 7 | 9 | -0.005 | -0.025 to 0.008 | 4 | the study split: decides nothing |

## Per arm

Recall stands beside log loss wherever log loss is given.

Per learner on the lineage split, balanced accuracy and recall side by side, with log loss, from `reference/cap-arms/page-count.tsv`. For `cap-all`, the recall medians run from -0.002 (`rf`, `gbt`) to 0.005 (`scm`), and the lowest end of a recall interval is -0.015, for `gbt`, so no recall interval lies wholly below -0.02.

| Pair | Learner | Entered | Balanced accuracy: better / worse / same | Balanced accuracy: median (interval) | Recall: better / worse / same | Recall: median (interval) | Log loss: median (interval) |
|---|---|---|---|---|---|---|---|
| cap-all less cap-2000 | `scm` | 76 of 116 | 17 / 7 / 52 | 0.002 (-0.003 to 0.005) | 26 / 18 / 32 | 0.005 (-0.001 to 0.009) | -0.026 (-0.054 to 0.027) |
| cap-all less cap-2000 | `cart` | 76 of 116 | 17 / 13 / 46 | 0.001 (0.000 to 0.007) | 30 / 21 / 25 | 0.002 (-0.003 to 0.017) | -0.010 (-0.018 to 0.001) |
| cap-all less cap-2000 | `rf` | 76 of 116 | 20 / 10 / 46 | 0.002 (-0.002 to 0.017) | 16 / 23 / 37 | -0.002 (-0.014 to 0.003) | 0.001 (-0.009 to 0.014) |
| cap-all less cap-2000 | `gbt` | 76 of 116 | 13 / 14 / 49 | -0.002 (-0.005 to 0.008) | 14 / 30 / 32 | -0.002 (-0.015 to 0.001) | 0.001 (-0.013 to 0.012) |
| cap-all less cap-2000 | `lasso` | 76 of 116 | 16 / 14 / 46 | -0.002 (-0.005 to 0.004) | 17 / 20 / 39 | -0.001 (-0.011 to 0.005) | 0.002 (-0.015 to 0.011) |
| twin less cap-2000 | `scm` | 85 of 143 | 12 / 13 / 60 | 0.000 (-0.001 to 0.001) | 21 / 26 / 38 | -0.005 (-0.013 to 0.000) | 0.005 (-0.007 to 0.034) |
| twin less cap-2000 | `cart` | 85 of 143 | 13 / 15 / 57 | 0.000 (-0.002 to 0.001) | 17 / 19 / 49 | 0.000 (-0.002 to 0.001) | 0.004 (0.000 to 0.006) |
| twin less cap-2000 | `rf` | 85 of 143 | 8 / 11 / 66 | -0.001 (-0.003 to 0.001) | 16 / 12 / 57 | 0.001 (-0.003 to 0.002) | 0.003 (0.000 to 0.007) |
| twin less cap-2000 | `gbt` | 85 of 143 | 13 / 13 / 59 | -0.001 (-0.003 to 0.001) | 25 / 20 / 40 | 0.001 (-0.001 to 0.003) | 0.000 (-0.003 to 0.008) |
| twin less cap-2000 | `lasso` | 85 of 143 | 13 / 15 / 57 | 0.000 (-0.001 to 0.001) | 26 / 16 / 43 | 0.001 (-0.004 to 0.012) | 0.000 (-0.007 to 0.012) |

## Beside it

These counts decide nothing.

**Without the few added.** The count over `cap-all`'s datasets less those that add fewer than 20 labelled genomes, made with the scorer's own counting code and Holm over its own rows. It enters 71 of 97 datasets on the lineage split, and no learner meets condition 1 or condition 2 there; the smallest Holm p is 0.260, for `scm` (`reference/cap-arms/page-count-without-few.tsv`).

| Split | Pair | Learner | Measure | Better is | Entered | Better | Worse | Same | Median | Interval | Species | Reading |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| lineage | cap-all less cap-2000 | `scm` | balanced accuracy | higher | 71 of 97 | 17 | 6 | 48 | 0.003 | -0.001 to 0.007 | 8 | beside: decides nothing |
| lineage | cap-all less cap-2000 | `scm` | recall | higher | 71 of 97 | 26 | 16 | 29 | 0.007 | -0.001 to 0.015 | 8 | beside: decides nothing |
| lineage | cap-all less cap-2000 | `cart` | balanced accuracy | higher | 71 of 97 | 16 | 13 | 42 | 0.001 | 0.000 to 0.006 | 8 | beside: decides nothing |
| lineage | cap-all less cap-2000 | `cart` | recall | higher | 71 of 97 | 29 | 20 | 22 | 0.002 | -0.005 to 0.021 | 8 | beside: decides nothing |
| lineage | cap-all less cap-2000 | `rf` | balanced accuracy | higher | 71 of 97 | 20 | 10 | 41 | 0.002 | -0.002 to 0.019 | 8 | beside: decides nothing |
| lineage | cap-all less cap-2000 | `rf` | recall | higher | 71 of 97 | 16 | 23 | 32 | -0.002 | -0.019 to 0.004 | 8 | beside: decides nothing |
| lineage | cap-all less cap-2000 | `gbt` | balanced accuracy | higher | 71 of 97 | 13 | 12 | 46 | -0.001 | -0.005 to 0.008 | 8 | beside: decides nothing |
| lineage | cap-all less cap-2000 | `gbt` | recall | higher | 71 of 97 | 14 | 27 | 30 | -0.001 | -0.012 to 0.003 | 8 | beside: decides nothing |
| lineage | cap-all less cap-2000 | `lasso` | balanced accuracy | higher | 71 of 97 | 16 | 12 | 43 | -0.002 | -0.004 to 0.004 | 8 | beside: decides nothing |
| lineage | cap-all less cap-2000 | `lasso` | recall | higher | 71 of 97 | 17 | 19 | 35 | -0.002 | -0.013 to 0.006 | 8 | beside: decides nothing |
| lineage | twin less cap-2000 | `scm` | balanced accuracy | higher | 85 of 143 | 12 | 13 | 60 | 0.000 | -0.001 to 0.001 | 8 | beside: decides nothing |
| lineage | twin less cap-2000 | `scm` | recall | higher | 85 of 143 | 21 | 26 | 38 | -0.005 | -0.013 to 0.000 | 8 | beside: decides nothing |
| lineage | twin less cap-2000 | `cart` | balanced accuracy | higher | 85 of 143 | 13 | 15 | 57 | 0.000 | -0.002 to 0.001 | 8 | beside: decides nothing |
| lineage | twin less cap-2000 | `cart` | recall | higher | 85 of 143 | 17 | 19 | 49 | 0.000 | -0.002 to 0.001 | 8 | beside: decides nothing |
| lineage | twin less cap-2000 | `rf` | balanced accuracy | higher | 85 of 143 | 8 | 11 | 66 | -0.001 | -0.003 to 0.001 | 8 | beside: decides nothing |
| lineage | twin less cap-2000 | `rf` | recall | higher | 85 of 143 | 16 | 12 | 57 | 0.001 | -0.003 to 0.002 | 8 | beside: decides nothing |
| lineage | twin less cap-2000 | `gbt` | balanced accuracy | higher | 85 of 143 | 13 | 13 | 59 | -0.001 | -0.003 to 0.001 | 8 | beside: decides nothing |
| lineage | twin less cap-2000 | `gbt` | recall | higher | 85 of 143 | 25 | 20 | 40 | 0.001 | -0.001 to 0.003 | 8 | beside: decides nothing |
| lineage | twin less cap-2000 | `lasso` | balanced accuracy | higher | 85 of 143 | 13 | 15 | 57 | 0.000 | -0.001 to 0.001 | 8 | beside: decides nothing |
| lineage | twin less cap-2000 | `lasso` | recall | higher | 85 of 143 | 26 | 16 | 43 | 0.001 | -0.004 to 0.012 | 8 | beside: decides nothing |
| study | cap-all less cap-2000 | `scm` | balanced accuracy | higher | 27 of 27 | 8 | 3 | 16 | 0.000 | -0.004 to 0.008 | 4 | beside: decides nothing |
| study | cap-all less cap-2000 | `scm` | recall | higher | 27 of 27 | 5 | 6 | 16 | -0.004 | -0.009 to 0.000 | 4 | beside: decides nothing |
| study | cap-all less cap-2000 | `cart` | balanced accuracy | higher | 27 of 27 | 3 | 8 | 16 | -0.009 | -0.013 to 0.000 | 4 | beside: decides nothing |
| study | cap-all less cap-2000 | `cart` | recall | higher | 27 of 27 | 5 | 9 | 13 | 0.000 | -0.007 to 0.000 | 4 | beside: decides nothing |
| study | cap-all less cap-2000 | `rf` | balanced accuracy | higher | 27 of 27 | 2 | 4 | 21 | 0.000 | -0.005 to 0.008 | 4 | beside: decides nothing |
| study | cap-all less cap-2000 | `rf` | recall | higher | 27 of 27 | 2 | 2 | 23 | -0.002 | -0.005 to 0.000 | 4 | beside: decides nothing |
| study | cap-all less cap-2000 | `gbt` | balanced accuracy | higher | 27 of 27 | 2 | 4 | 21 | -0.001 | -0.007 to 0.001 | 4 | beside: decides nothing |
| study | cap-all less cap-2000 | `gbt` | recall | higher | 27 of 27 | 9 | 6 | 12 | -0.003 | -0.004 to 0.014 | 4 | beside: decides nothing |
| study | cap-all less cap-2000 | `lasso` | balanced accuracy | higher | 27 of 27 | 4 | 3 | 20 | 0.002 | -0.006 to 0.003 | 4 | beside: decides nothing |
| study | cap-all less cap-2000 | `lasso` | recall | higher | 27 of 27 | 4 | 7 | 16 | -0.004 | -0.009 to -0.003 | 4 | beside: decides nothing |

The 19 datasets that add fewer than 20 labelled genomes, with the added genomes by class (`reference/cap-arms/page-few-added.tsv`):

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

**Study split.** The study split decides nothing. Its fits, one a dataset and arm, are counted beside the lineage count. `cap-all` enters 27 of 27 study datasets in 4 species; its rows stand in the count table above, marked "the study split: decides nothing".

## Seen after

These points were seen in the tables after the reading. They change no rule and no reading above, and decide nothing.

**The largest drop.** `cap-all`'s largest drop in lineage balanced accuracy is *Mycobacterium tuberculosis*, amoxicillin, `rf`: 0.655 for `cap-2000` and 0.508 for `cap-all`, on 57 genomes, where all 130 added labelled genomes are susceptible; its recall goes from 0.560 to 0.027. The twin, a refit of the same 2,000 genomes at other slots, has a drop as large: *Salmonella enterica*, sulfisoxazole, `lasso`, 0.876 for `cap-2000` and 0.725 for the twin. The largest gain is *Klebsiella pneumoniae*, nitrofurantoin, `scm`: 0.548 and 0.680 (`reference/cap-arms/page-datasets.tsv`).

**Recall moves more often.** Of `cap-all`'s 380 lineage rows of a dataset and a learner, recall is better by more than 0.02 in 103 and worse in 112, and balanced accuracy is better in 83 and worse in 58 (`reference/cap-arms/page-datasets.tsv`).

**Per species.** In *Salmonella enterica*, `cap-all` is better in 35 and worse in 6 of 90 lineage pairs, and the twin better in 13 and worse in 22 of 90. In *Klebsiella pneumoniae*, `cap-all` is better in 10 and worse in 26 of 95, and the twin better in 17 and worse in 14 of 95. In *Staphylococcus aureus*, `cap-all` is better in 10 and worse in 0 of 15 (`reference/cap-arms/page-species.tsv`).

**Cost.** `cap-all`'s lineage fits took a median of 193.05 s of wall time against 87.91 s for `cap-2000`'s, with a median peak memory of 5,722,228 KB against 2,715,660 KB and a largest of 17,952,344 KB (`reference/cap-arms/page-times.tsv`).

## Limits

**One database.** Every split is a retrospective split inside one public database. A split by lineage or by study is not external validation.

**No difference shown.** "No difference shown" is the rule's reading. It is not a statement that two caps are equal.

**Added genomes.** An added genome only trains and is never scored, so every arm is read on the default's held-out genomes.

**Withdrawn.** `cap-4000` was withdrawn before any fit. It holds no result, and its five lineage pairs stand in the family at p = 1.

**Ties at the cut.** The 4,096 cut can fall inside a tie of the selection rule, and the tied genomes added are then picked by id. The genomes that share the rule's score at the 4,096th place number 974 of the 5,396 the rule ranks in *Neisseria gonorrhoeae*, 946 of 5,045 in *Staphylococcus aureus*, 538 of 8,170 in *Salmonella enterica*, 233 of 8,704 in *Escherichia coli*, 121 of 6,159 in *Klebsiella pneumoniae*, 70 of 21,295 in *Mycobacterium tuberculosis* and 12 of 15,530 in *Streptococcus pneumoniae*, and *Enterococcus faecium*'s rule ranks 2,907 genomes, so it has no cut (`reference/cap-arms-checks/ties.tsv`).

**Unavailable lineage splits.** 58 of the 143 datasets have no lineage split on any arm, so the lineage count reads 76 of `cap-all`'s 116 fitted datasets and 85 of 143 for the twin.

**The study split decides nothing.** No condition that decides reads it.

**Labels.** The datasets of M wait for the label rule.

**No final model.** No final model is built here. `grm arm` saves no final model, so a cap arm that passed would need a final fit on its cohort per dataset before block T1 could score it.

**The release.** The release keeps its cap of 2,000 genomes a species and its pipeline's rule whatever the count says, until the maintainer rules.

**Markers.** No marker is found, named or ranked. The page holds aggregates over datasets, species, splits, arms and learners only.

## Method

**Checks of the run.** Before any score was read, `scripts/c2-cap/run_checks.py` read each species' copy and printed counts only: the species, every `SHA256SUMS`, `DONE` and `binary.txt`, one binary, `expect.py`'s rows, the cohorts against the selection rule's counts, genome ids and BioProject accessions, and the refused lineage draws. It reads the tables of each species' `score.tar` in place, without unpacking it, and `score/score.log` is never read or hashed by it. A copy that lacks `lists.txt` takes it rebuilt from `lists.tsv` and the driver's arguments, and the rebuilt text must match its `SHA256SUMS` line. All 8 checks passed: 257 of 257 `SHA256SUMS` lines hold, 8 of 8 `expect.tsv` files are `expect.py`'s with 1,365 rows, 8 of 8 cohorts equal the selection rule's counts, and 0 genome ids and 0 BioProject accessions stand in 1,306 files (`reference/cap-arms/checks.txt`).

**Genome ids of the added genomes.** Before the reading, the genome ids in the files each species' run wrote were counted where that species ran, against the release's ids and the ids the wider cohort adds, 15,579 beyond the release over the 8 species (`reference/cap-arms-checks/added-genome-ids.tsv`, named per species), and 0 were found in 8 of 8 species, so `scripts/c2-cap/run_checks.py` ran with the release directory and no id list, since the cohort lists stay where each species ran.

**Tables.** `scripts/c2-cap/take_tables.py` takes `score/checks.tsv`, `score/draws.tsv` and `score/datasets.tsv` out of each species' `score.tar` by name, and no other member: it walks the archive's headers and seeks past every other member's bytes, so `score/score.log` is never read or written to disk. `scripts/c2-cap/page_tables.py write` then takes `run_checks.py`'s printout, refuses to write unless it holds a pass line for each of `run_checks.py`'s checks and no line starting FAIL (gate 1), checks the reading against its other gates and writes `reference/cap-arms/`. `scripts/c2-cap/page_tables.py print` prints this page's tables from those files, rounded to 3 decimals. The reading ran each step once: `run_checks.py`, `take_tables.py` on the 8 species' `score.tar`, `score.py --count` with `cap-4000` withdrawn over the 8 species' `datasets.tsv` (2,340 dataset rows), `page_tables.py write` and `page_tables.py print`.

**Ties.** `scripts/c2-cap/ties.py` counts, per species, the genomes that share the selection rule's score at the 4,096th place, scoring each genome as `cohorts.py` does and holding its order to `cohorts.py`'s ranking. It writes one row per species, counts only. It ran on the ranking's inputs as the cohorts were listed, and its `ranked` column equals `docs/cap-arms.md`'s selection counts in 8 of 8 species (`reference/cap-arms-checks/ties.tsv`).

**Tests.** `tests/c2_cap/test_page_tables.py`, `tests/c2_cap/test_run_checks.py`, `tests/c2_cap/test_review_checks.py`, `tests/c2_cap/test_review_faults.py`, `tests/c2_cap/test_take_tables.py` and `tests/c2_cap/test_ties.py` hold SYNTHETIC tests of the four scripts. `scripts/c2-cap/page_faults.py` plants each fault of the table script in a copy of the tree and runs the test meant to catch it.

A column that reads 0 by construction stays off the printed page: C2's `score.py` writes `constant` 0 on every draw, so that column is not printed, and its count writes no count without constant learners.

The scorer's checks, passed of made (`reference/cap-arms/page-checks.tsv`):

| Check | Passed | Species |
|---|---|---|
| cart has rows | 738 of 738 | 8 |
| cart: one row per scored genome in each of 5 folds | 738 of 738 | 8 |
| gbt has rows | 738 of 738 | 8 |
| gbt: one row per scored genome in each of 5 folds | 738 of 738 | 8 |
| lasso has rows | 738 of 738 | 8 |
| lasso: one row per scored genome in each of 5 folds | 738 of 738 | 8 |
| rf has rows | 738 of 738 | 8 |
| rf: one row per scored genome in each of 5 folds | 738 of 738 | 8 |
| scm has rows | 738 of 738 | 8 |
| scm: one row per scored genome in each of 5 folds | 738 of 738 | 8 |
| the lineage split is unavailable (a default plan refused past grm arm's replacements); not produced | 468 of 468 | 8 |
| the pick file is the dataset's study-whole fit at seed 42 with no inner fold | 66 of 66 | 5 |
| the plan and the prediction rows exist | 738 of 738 | 8 |
| the scored genomes of each outer fold are the default's | 228 of 228 | 8 |
| the test folds are the default's, with at least its training genomes | 27 of 27 | 4 |
| the twin is the default's fit at a slot no draw has | 255 of 255 | 8 |

The binary, one commit and one code sha256 for the 8 species (`reference/cap-arms/page-binaries.tsv`):

| Commit | Code sha256 | Species | Binary sha256s |
|---|---|---|---|
| `ee04e21` | `cc878fced544` | Enterococcus faecium, Escherichia coli, Klebsiella pneumoniae, Mycobacterium tuberculosis, Neisseria gonorrhoeae, Salmonella enterica, Staphylococcus aureus, Streptococcus pneumoniae | `8cc865103d84`: Enterococcus faecium, Escherichia coli, Klebsiella pneumoniae, Mycobacterium tuberculosis, Neisseria gonorrhoeae, Salmonella enterica, Staphylococcus aureus, Streptococcus pneumoniae |

The fits' times, 325 fits in all (`reference/cap-arms/page-times.tsv`):

| Arm | Split | Fits | Wall s, sum | Wall s, median | CPU s, sum | CPU s, median | Peak KB, median | Peak KB, largest |
|---|---|---|---|---|---|---|---|---|
| `cap-2000` | lineage | 143 | 23163.47 | 87.91 | 65434.69 | 210.9 | 2,715,660 | 8,761,564 |
| `cap-2000` | study | 39 | 940.69 | 21.22 | 5071.21 | 129.97 | 2,729,208 | 3,762,248 |
| `cap-all` | lineage | 116 | 28508.18 | 193.05 | 93048.08 | 462.04 | 5,722,228 | 17,952,344 |
| `cap-all` | study | 27 | 1048.19 | 34.69 | 5681.09 | 181.37 | 5,427,848 | 8,123,680 |

## Datasets

Per species, the datasets of the first run, those in S, the cohort sizes, the added genomes, and per dataset the default's and `cap-all`'s figures, with each dataset not fitted named with its reason.

Each dataset an arm does not score, with its reason (`reference/cap-arms/page-not-fitted.tsv`):

| Species | Antibiotic | Arm | Split | Reason |
|---|---|---|---|---|
| *Enterococcus faecium* | clindamycin | cap-all | lineage | none-added: no added genome labelled |
| *Enterococcus faecium* | dalfopristin/quinupristin | cap-all | lineage | none-added: no added genome labelled |
| *Enterococcus faecium* | doxycycline | cap-all | lineage | none-added: no added genome labelled |
| *Enterococcus faecium* | erythromycin | cap-all | lineage | none-added: no added genome labelled |
| *Enterococcus faecium* | gentamicin | cap-all | lineage | none-added: no added genome labelled |
| *Enterococcus faecium* | kanamycin | cap-all | lineage | none-added: no added genome labelled |
| *Enterococcus faecium* | streptomycin | cap-all | lineage | none-added: no added genome labelled |
| *Enterococcus faecium* | teicoplanin | cap-all | lineage | none-added: no added genome labelled |
| *Enterococcus faecium* | doxycycline | every arm and the twin | lineage | lineage split unavailable: the default's lineage draws refused past grm arm's two replacements |
| *Enterococcus faecium* | teicoplanin | every arm and the twin | lineage | lineage split unavailable: the default's lineage draws refused past grm arm's two replacements |
| *Enterococcus faecium* | vancomycin | every arm and the twin | lineage | lineage split unavailable: the default's lineage draws refused past grm arm's two replacements |
| *Escherichia coli* | cefpodoxime | cap-all | lineage | none-added: no added genome labelled |
| *Escherichia coli* | cefpodoxime/clavulanate | cap-all | lineage | none-added: no added genome labelled |
| *Escherichia coli* | ceftiofur | cap-all | lineage | none-added: no added genome labelled |
| *Escherichia coli* | colistin | cap-all | lineage | none-added: no added genome labelled |
| *Escherichia coli* | erythromycin | cap-all | lineage | none-added: no added genome labelled |
| *Escherichia coli* | nitrofurantoin | cap-all | lineage and study | none-added: no added genome labelled |
| *Escherichia coli* | streptomycin | cap-all | lineage and study | none-added: no added genome labelled |
| *Escherichia coli* | sulfamethoxazole | cap-all | lineage and study | none-added: no added genome labelled |
| *Escherichia coli* | ampicillin/sulbactam | cap-all | study | none-added: no added genome trains on the study split |
| *Escherichia coli* | cefazolin | cap-all | study | none-added: no added genome trains on the study split |
| *Escherichia coli* | ceftriaxone | cap-all | study | none-added: no added genome trains on the study split |
| *Escherichia coli* | levofloxacin | cap-all | study | none-added: no added genome trains on the study split |
| *Escherichia coli* | amikacin | every arm and the twin | lineage | lineage split unavailable: the default's lineage draws refused past grm arm's two replacements |
| *Escherichia coli* | ampicillin/sulbactam | every arm and the twin | lineage | lineage split unavailable: the default's lineage draws refused past grm arm's two replacements |
| *Escherichia coli* | aztreonam | every arm and the twin | lineage | lineage split unavailable: the default's lineage draws refused past grm arm's two replacements |
| *Escherichia coli* | cefazolin | every arm and the twin | lineage | lineage split unavailable: the default's lineage draws refused past grm arm's two replacements |
| *Escherichia coli* | cefepime | every arm and the twin | lineage | lineage split unavailable: the default's lineage draws refused past grm arm's two replacements |
| *Escherichia coli* | cefpodoxime | every arm and the twin | lineage | lineage split unavailable: the default's lineage draws refused past grm arm's two replacements |
| *Escherichia coli* | cefpodoxime/clavulanate | every arm and the twin | lineage | lineage split unavailable: the default's lineage draws refused past grm arm's two replacements |
| *Escherichia coli* | ceftiofur | every arm and the twin | lineage | lineage split unavailable: the default's lineage draws refused past grm arm's two replacements |
| *Escherichia coli* | cefuroxime | every arm and the twin | lineage | lineage split unavailable: the default's lineage draws refused past grm arm's two replacements |
| *Escherichia coli* | chloramphenicol | every arm and the twin | lineage | lineage split unavailable: the default's lineage draws refused past grm arm's two replacements |
| *Escherichia coli* | colistin | every arm and the twin | lineage | lineage split unavailable: the default's lineage draws refused past grm arm's two replacements |
| *Escherichia coli* | doxycycline | every arm and the twin | lineage | lineage split unavailable: the default's lineage draws refused past grm arm's two replacements |
| *Escherichia coli* | erythromycin | every arm and the twin | lineage | lineage split unavailable: the default's lineage draws refused past grm arm's two replacements |
| *Escherichia coli* | imipenem | every arm and the twin | lineage | lineage split unavailable: the default's lineage draws refused past grm arm's two replacements |
| *Escherichia coli* | nalidixic acid | every arm and the twin | lineage | lineage split unavailable: the default's lineage draws refused past grm arm's two replacements |
| *Escherichia coli* | nitrofurantoin | every arm and the twin | lineage | lineage split unavailable: the default's lineage draws refused past grm arm's two replacements |
| *Escherichia coli* | streptomycin | every arm and the twin | lineage | lineage split unavailable: the default's lineage draws refused past grm arm's two replacements |
| *Escherichia coli* | tobramycin | every arm and the twin | lineage | lineage split unavailable: the default's lineage draws refused past grm arm's two replacements |
| *Klebsiella pneumoniae* | cefotetan | cap-all | lineage | none-added: no added genome labelled |
| *Klebsiella pneumoniae* | ceftazidime/avibactam | cap-all | lineage | none-added: no added genome labelled |
| *Klebsiella pneumoniae* | ofloxacin | cap-all | lineage | none-added: no added genome labelled |
| *Klebsiella pneumoniae* | streptomycin | cap-all | lineage | none-added: no added genome labelled |
| *Klebsiella pneumoniae* | cefepime | cap-all | study | none-added: no added genome trains on the study split |
| *Klebsiella pneumoniae* | ceftriaxone | cap-all | study | none-added: no added genome trains on the study split |
| *Klebsiella pneumoniae* | doripenem | cap-all | study | none-added: no added genome trains on the study split |
| *Klebsiella pneumoniae* | polymyxin b | cap-all | study | none-added: no added genome trains on the study split |
| *Klebsiella pneumoniae* | cefotetan | every arm and the twin | lineage | lineage split unavailable: the default's lineage draws refused past grm arm's two replacements |
| *Klebsiella pneumoniae* | cefpodoxime/ceftriaxone | every arm and the twin | lineage | lineage split unavailable: the default's lineage draws refused past grm arm's two replacements |
| *Klebsiella pneumoniae* | cefpodoxime/clavulanate | every arm and the twin | lineage | lineage split unavailable: the default's lineage draws refused past grm arm's two replacements |
| *Klebsiella pneumoniae* | ceftazidime/avibactam | every arm and the twin | lineage | lineage split unavailable: the default's lineage draws refused past grm arm's two replacements |
| *Klebsiella pneumoniae* | chloramphenicol | every arm and the twin | lineage | lineage split unavailable: the default's lineage draws refused past grm arm's two replacements |
| *Klebsiella pneumoniae* | colistin | every arm and the twin | lineage | lineage split unavailable: the default's lineage draws refused past grm arm's two replacements |
| *Klebsiella pneumoniae* | doripenem | every arm and the twin | lineage | lineage split unavailable: the default's lineage draws refused past grm arm's two replacements |
| *Klebsiella pneumoniae* | doxycycline | every arm and the twin | lineage | lineage split unavailable: the default's lineage draws refused past grm arm's two replacements |
| *Klebsiella pneumoniae* | fosfomycin | every arm and the twin | lineage | lineage split unavailable: the default's lineage draws refused past grm arm's two replacements |
| *Klebsiella pneumoniae* | ofloxacin | every arm and the twin | lineage | lineage split unavailable: the default's lineage draws refused past grm arm's two replacements |
| *Klebsiella pneumoniae* | polymyxin b | every arm and the twin | lineage | lineage split unavailable: the default's lineage draws refused past grm arm's two replacements |
| *Klebsiella pneumoniae* | streptomycin | every arm and the twin | lineage | lineage split unavailable: the default's lineage draws refused past grm arm's two replacements |
| *Klebsiella pneumoniae* | ticarcillin/clavulanate | every arm and the twin | lineage | lineage split unavailable: the default's lineage draws refused past grm arm's two replacements |
| *Mycobacterium tuberculosis* | nicotinamide | cap-all | lineage and study | none-added: no added genome labelled |
| *Mycobacterium tuberculosis* | ciprofloxacin | every arm and the twin | lineage | lineage split unavailable: the default's lineage draws refused past grm arm's two replacements |
| *Mycobacterium tuberculosis* | nicotinamide | every arm and the twin | lineage | lineage split unavailable: the default's lineage draws refused past grm arm's two replacements |
| *Neisseria gonorrhoeae* | erythromycin | cap-all | lineage | none-added: no added genome labelled |
| *Neisseria gonorrhoeae* | cefixime | every arm and the twin | lineage | lineage split unavailable: the default's lineage draws refused past grm arm's two replacements |
| *Neisseria gonorrhoeae* | ceftriaxone | every arm and the twin | lineage | lineage split unavailable: the default's lineage draws refused past grm arm's two replacements |
| *Neisseria gonorrhoeae* | erythromycin | every arm and the twin | lineage | lineage split unavailable: the default's lineage draws refused past grm arm's two replacements |
| *Salmonella enterica* | ampicillin/sulbactam | cap-all | lineage | none-added: no added genome labelled |
| *Salmonella enterica* | cefazolin | cap-all | lineage | none-added: no added genome labelled |
| *Salmonella enterica* | pefloxacin | cap-all | lineage | none-added: no added genome labelled |
| *Salmonella enterica* | ampicillin/sulbactam | every arm and the twin | lineage | lineage split unavailable: the default's lineage draws refused past grm arm's two replacements |
| *Salmonella enterica* | cefazolin | every arm and the twin | lineage | lineage split unavailable: the default's lineage draws refused past grm arm's two replacements |
| *Salmonella enterica* | ceftazidime | every arm and the twin | lineage | lineage split unavailable: the default's lineage draws refused past grm arm's two replacements |
| *Salmonella enterica* | pefloxacin | every arm and the twin | lineage | lineage split unavailable: the default's lineage draws refused past grm arm's two replacements |
| *Salmonella enterica* | sulfamethoxazole | every arm and the twin | lineage | lineage split unavailable: the default's lineage draws refused past grm arm's two replacements |
| *Salmonella enterica* | sulfathiazole | every arm and the twin | lineage | lineage split unavailable: the default's lineage draws refused past grm arm's two replacements |
| *Staphylococcus aureus* | cefoxitin | every arm and the twin | lineage | lineage split unavailable: the default's lineage draws refused past grm arm's two replacements |
| *Staphylococcus aureus* | chloramphenicol | every arm and the twin | lineage | lineage split unavailable: the default's lineage draws refused past grm arm's two replacements |
| *Staphylococcus aureus* | ciprofloxacin | every arm and the twin | lineage | lineage split unavailable: the default's lineage draws refused past grm arm's two replacements |
| *Staphylococcus aureus* | daptomycin | every arm and the twin | lineage | lineage split unavailable: the default's lineage draws refused past grm arm's two replacements |
| *Staphylococcus aureus* | fusidic acid | every arm and the twin | lineage | lineage split unavailable: the default's lineage draws refused past grm arm's two replacements |
| *Staphylococcus aureus* | gentamicin | every arm and the twin | lineage | lineage split unavailable: the default's lineage draws refused past grm arm's two replacements |
| *Staphylococcus aureus* | mupirocin | every arm and the twin | lineage | lineage split unavailable: the default's lineage draws refused past grm arm's two replacements |
| *Staphylococcus aureus* | rifampin | every arm and the twin | lineage | lineage split unavailable: the default's lineage draws refused past grm arm's two replacements |
| *Staphylococcus aureus* | sulfamethoxazole/trimethoprim | every arm and the twin | lineage | lineage split unavailable: the default's lineage draws refused past grm arm's two replacements |
| *Staphylococcus aureus* | tetracycline | every arm and the twin | lineage | lineage split unavailable: the default's lineage draws refused past grm arm's two replacements |
| *Streptococcus pneumoniae* | azithromycin | cap-all | lineage | none-added: no added genome labelled |
| *Streptococcus pneumoniae* | cefuroxime | cap-all | lineage | none-added: no added genome labelled |
| *Streptococcus pneumoniae* | chloramphenicol | every arm and the twin | lineage | lineage split unavailable: the default's lineage draws refused past grm arm's two replacements |
| *Streptococcus pneumoniae* | clindamycin | every arm and the twin | lineage | lineage split unavailable: the default's lineage draws refused past grm arm's two replacements |
| *Streptococcus pneumoniae* | meropenem | every arm and the twin | lineage | lineage split unavailable: the default's lineage draws refused past grm arm's two replacements |

Per dataset, split and learner, `cap-all`'s balanced accuracy less `cap-2000`'s (`reference/cap-arms/page-datasets.tsv`, which also holds each arm's balanced accuracy and recall, the twin's figures and the draws scored):

| Species | Antibiotic | Split | `cap-all` role | Genomes | Added labelled | `scm` difference | `cart` difference | `rf` difference | `gbt` difference | `lasso` difference |
|---|---|---|---|---|---|---|---|---|---|---|
| *Enterococcus faecium* | ampicillin | lineage | draw | 1,680 | 188 (188 R, 0 S) | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| *Enterococcus faecium* | clindamycin | lineage | none-added | 901 | 0 (0 R, 0 S) | none | none | none | none | none |
| *Enterococcus faecium* | dalfopristin/quinupristin | lineage | none-added | 1,118 | 0 (0 R, 0 S) | none | none | none | none | none |
| *Enterococcus faecium* | doxycycline | lineage | lineage split unavailable | 345 | 0 (0 R, 0 S) | none | none | none | none | none |
| *Enterococcus faecium* | erythromycin | lineage | none-added | 1,151 | 0 (0 R, 0 S) | none | none | none | none | none |
| *Enterococcus faecium* | gentamicin | lineage | none-added | 1,345 | 0 (0 R, 0 S) | none | none | none | none | none |
| *Enterococcus faecium* | kanamycin | lineage | none-added | 944 | 0 (0 R, 0 S) | none | none | none | none | none |
| *Enterococcus faecium* | streptomycin | lineage | none-added | 1,239 | 0 (0 R, 0 S) | none | none | none | none | none |
| *Enterococcus faecium* | teicoplanin | lineage | lineage split unavailable | 1,192 | 0 (0 R, 0 S) | none | none | none | none | none |
| *Enterococcus faecium* | tetracycline | lineage | draw | 958 | 6 (6 R, 0 S) | 0.001 | 0.007 | 0.001 | -0.002 | -0.005 |
| *Enterococcus faecium* | vancomycin | lineage | lineage split unavailable | 1,741 | 901 (373 R, 528 S) | none | none | none | none | none |
| *Escherichia coli* | amikacin | lineage | lineage split unavailable | 1,125 | 284 (0 R, 284 S) | none | none | none | none | none |
| *Escherichia coli* | amikacin | study | study | 1,125 | 284 (0 R, 284 S) | 0.035 | 0.045 | -0.033 | 0.002 | -0.020 |
| *Escherichia coli* | ampicillin/sulbactam | lineage | lineage split unavailable | 283 | 1 (1 R, 0 S) | none | none | none | none | none |
| *Escherichia coli* | ampicillin/sulbactam | study | none-added | 283 | 1 (1 R, 0 S) | none | none | none | none | none |
| *Escherichia coli* | aztreonam | lineage | lineage split unavailable | 699 | 349 (3 R, 346 S) | none | none | none | none | none |
| *Escherichia coli* | aztreonam | study | study | 699 | 349 (3 R, 346 S) | -0.044 | -0.022 | -0.012 | -0.029 | -0.015 |
| *Escherichia coli* | cefazolin | lineage | lineage split unavailable | 414 | 2 (2 R, 0 S) | none | none | none | none | none |
| *Escherichia coli* | cefazolin | study | none-added | 414 | 2 (2 R, 0 S) | none | none | none | none | none |
| *Escherichia coli* | cefepime | lineage | lineage split unavailable | 1,124 | 100 (1 R, 99 S) | none | none | none | none | none |
| *Escherichia coli* | cefepime | study | study | 1,124 | 100 (1 R, 99 S) | -0.004 | -0.001 | 0.000 | 0.006 | 0.004 |
| *Escherichia coli* | cefotaxime | lineage | draw | 1,149 | 2,012 (456 R, 1556 S) | -0.011 | 0.004 | -0.007 | 0.006 | 0.008 |
| *Escherichia coli* | cefoxitin | lineage | draw | 1,066 | 94 (1 R, 93 S) | -0.004 | 0.002 | 0.001 | 0.005 | 0.010 |
| *Escherichia coli* | cefoxitin | study | study | 1,066 | 94 (1 R, 93 S) | -0.013 | -0.038 | -0.005 | 0.033 | 0.005 |
| *Escherichia coli* | cefpodoxime | lineage | lineage split unavailable | 175 | 0 (0 R, 0 S) | none | none | none | none | none |
| *Escherichia coli* | cefpodoxime/clavulanate | lineage | lineage split unavailable | 136 | 0 (0 R, 0 S) | none | none | none | none | none |
| *Escherichia coli* | ceftazidime | lineage | draw | 1,466 | 1,976 (349 R, 1627 S) | -0.003 | 0.000 | 0.021 | 0.006 | 0.004 |
| *Escherichia coli* | ceftiofur | lineage | lineage split unavailable | 203 | 0 (0 R, 0 S) | none | none | none | none | none |
| *Escherichia coli* | ceftriaxone | lineage | draw | 923 | 3 (2 R, 1 S) | -0.003 | 0.000 | -0.006 | -0.007 | 0.000 |
| *Escherichia coli* | ceftriaxone | study | none-added | 923 | 3 (2 R, 1 S) | none | none | none | none | none |
| *Escherichia coli* | cefuroxime | lineage | lineage split unavailable | 548 | 1,801 (405 R, 1396 S) | none | none | none | none | none |
| *Escherichia coli* | cephalothin | lineage | draw | 472 | 216 (81 R, 135 S) | 0.013 | 0.002 | -0.013 | -0.012 | -0.039 |
| *Escherichia coli* | chloramphenicol | lineage | lineage split unavailable | 709 | 7 (5 R, 2 S) | none | none | none | none | none |
| *Escherichia coli* | colistin | lineage | lineage split unavailable | 239 | 0 (0 R, 0 S) | none | none | none | none | none |
| *Escherichia coli* | doxycycline | lineage | lineage split unavailable | 227 | 2 (2 R, 0 S) | none | none | none | none | none |
| *Escherichia coli* | ertapenem | lineage | draw | 820 | 108 (1 R, 107 S) | -0.011 | -0.052 | -0.002 | -0.011 | -0.005 |
| *Escherichia coli* | ertapenem | study | study | 820 | 108 (1 R, 107 S) | -0.007 | -0.113 | -0.002 | -0.011 | 0.041 |
| *Escherichia coli* | erythromycin | lineage | lineage split unavailable | 205 | 0 (0 R, 0 S) | none | none | none | none | none |
| *Escherichia coli* | imipenem | lineage | lineage split unavailable | 614 | 839 (0 R, 839 S) | none | none | none | none | none |
| *Escherichia coli* | levofloxacin | lineage | draw | 403 | 2 (1 R, 1 S) | -0.079 | 0.000 | 0.005 | -0.049 | 0.014 |
| *Escherichia coli* | levofloxacin | study | none-added | 403 | 2 (1 R, 1 S) | none | none | none | none | none |
| *Escherichia coli* | nalidixic acid | lineage | lineage split unavailable | 406 | 479 (85 R, 394 S) | none | none | none | none | none |
| *Escherichia coli* | nitrofurantoin | lineage | lineage split unavailable | 482 | 0 (0 R, 0 S) | none | none | none | none | none |
| *Escherichia coli* | nitrofurantoin | study | none-added | 482 | 0 (0 R, 0 S) | none | none | none | none | none |
| *Escherichia coli* | piperacillin/tazobactam | lineage | draw | 979 | 1,739 (179 R, 1560 S) | -0.019 | -0.019 | -0.057 | -0.004 | -0.018 |
| *Escherichia coli* | piperacillin/tazobactam | study | study | 979 | 1,739 (179 R, 1560 S) | 0.030 | 0.005 | 0.003 | 0.013 | 0.002 |
| *Escherichia coli* | streptomycin | lineage | lineage split unavailable | 386 | 0 (0 R, 0 S) | none | none | none | none | none |
| *Escherichia coli* | streptomycin | study | none-added | 386 | 0 (0 R, 0 S) | none | none | none | none | none |
| *Escherichia coli* | sulfamethoxazole | lineage | none-added | 341 | 0 (0 R, 0 S) | none | none | none | none | none |
| *Escherichia coli* | sulfamethoxazole | study | none-added | 341 | 0 (0 R, 0 S) | none | none | none | none | none |
| *Escherichia coli* | tetracycline | lineage | draw | 782 | 65 (38 R, 27 S) | -0.007 | -0.006 | -0.002 | -0.019 | -0.006 |
| *Escherichia coli* | tetracycline | study | study | 782 | 65 (38 R, 27 S) | -0.007 | -0.015 | 0.001 | -0.024 | 0.010 |
| *Escherichia coli* | tobramycin | lineage | lineage split unavailable | 1,070 | 490 (22 R, 468 S) | none | none | none | none | none |
| *Escherichia coli* | tobramycin | study | study | 1,070 | 490 (22 R, 468 S) | 0.004 | 0.002 | 0.000 | -0.006 | 0.003 |
| *Escherichia coli* | trimethoprim | lineage | draw | 861 | 308 (75 R, 233 S) | 0.039 | 0.023 | 0.016 | -0.005 | 0.015 |
| *Klebsiella pneumoniae* | ampicillin/sulbactam | lineage | draw | 744 | 1,244 (1210 R, 34 S) | 0.020 | 0.015 | -0.011 | -0.006 | -0.046 |
| *Klebsiella pneumoniae* | ampicillin/sulbactam | study | study | 744 | 1,244 (1210 R, 34 S) | 0.052 | -0.097 | 0.034 | 0.018 | -0.006 |
| *Klebsiella pneumoniae* | aztreonam | lineage | draw | 1,050 | 1,529 (1462 R, 67 S) | 0.034 | 0.012 | -0.003 | 0.014 | 0.012 |
| *Klebsiella pneumoniae* | aztreonam | study | study | 1,050 | 1,529 (1462 R, 67 S) | -0.021 | 0.009 | 0.011 | -0.002 | -0.026 |
| *Klebsiella pneumoniae* | cefazolin | lineage | draw | 1,141 | 1,348 (1308 R, 40 S) | 0.006 | 0.020 | -0.013 | -0.003 | -0.035 |
| *Klebsiella pneumoniae* | cefazolin | study | study | 1,141 | 1,348 (1308 R, 40 S) | -0.010 | -0.020 | -0.039 | -0.012 | -0.044 |
| *Klebsiella pneumoniae* | cefepime | lineage | draw | 1,233 | 1,163 (924 R, 239 S) | 0.008 | -0.030 | -0.042 | -0.056 | -0.066 |
| *Klebsiella pneumoniae* | cefepime | study | none-added | 1,233 | 1,163 (924 R, 239 S) | none | none | none | none | none |
| *Klebsiella pneumoniae* | cefotetan | lineage | lineage split unavailable | 61 | 0 (0 R, 0 S) | none | none | none | none | none |
| *Klebsiella pneumoniae* | cefoxitin | lineage | draw | 1,181 | 1,371 (911 R, 460 S) | 0.060 | 0.021 | -0.029 | -0.114 | -0.111 |
| *Klebsiella pneumoniae* | cefoxitin | study | study | 1,181 | 1,371 (911 R, 460 S) | 0.001 | 0.006 | -0.004 | -0.025 | -0.022 |
| *Klebsiella pneumoniae* | cefpodoxime | lineage | draw | 221 | 6 (1 R, 5 S) | -0.003 | 0.090 | -0.002 | 0.009 | -0.064 |
| *Klebsiella pneumoniae* | cefpodoxime/ceftriaxone | lineage | lineage split unavailable | 56 | 3 (2 R, 1 S) | none | none | none | none | none |
| *Klebsiella pneumoniae* | cefpodoxime/clavulanate | lineage | lineage split unavailable | 207 | 5 (0 R, 5 S) | none | none | none | none | none |
| *Klebsiella pneumoniae* | ceftazidime | lineage | draw | 1,310 | 1,988 (1597 R, 391 S) | 0.001 | -0.011 | 0.013 | 0.014 | -0.002 |
| *Klebsiella pneumoniae* | ceftazidime | study | study | 1,310 | 1,988 (1597 R, 391 S) | 0.022 | -0.012 | 0.008 | 0.015 | 0.066 |
| *Klebsiella pneumoniae* | ceftazidime/avibactam | lineage | lineage split unavailable | 259 | 0 (0 R, 0 S) | none | none | none | none | none |
| *Klebsiella pneumoniae* | ceftriaxone | lineage | draw | 1,208 | 1,315 (1273 R, 42 S) | 0.000 | -0.009 | -0.010 | -0.003 | -0.021 |
| *Klebsiella pneumoniae* | ceftriaxone | study | none-added | 1,208 | 1,315 (1273 R, 42 S) | none | none | none | none | none |
| *Klebsiella pneumoniae* | cefuroxime | lineage | draw | 565 | 1,387 (1387 R, 0 S) | -0.051 | -0.007 | -0.025 | -0.036 | -0.035 |
| *Klebsiella pneumoniae* | chloramphenicol | lineage | lineage split unavailable | 409 | 11 (8 R, 3 S) | none | none | none | none | none |
| *Klebsiella pneumoniae* | ciprofloxacin | lineage | draw | 1,706 | 1,834 (1420 R, 414 S) | -0.002 | -0.034 | -0.044 | -0.005 | 0.023 |
| *Klebsiella pneumoniae* | ciprofloxacin | study | study | 1,706 | 1,834 (1420 R, 414 S) | 0.028 | -0.021 | 0.009 | 0.008 | 0.002 |
| *Klebsiella pneumoniae* | colistin | lineage | lineage split unavailable | 336 | 50 (10 R, 40 S) | none | none | none | none | none |
| *Klebsiella pneumoniae* | colistin | study | study | 336 | 50 (10 R, 40 S) | -0.018 | 0.174 | -0.020 | 0.029 | 0.025 |
| *Klebsiella pneumoniae* | doripenem | lineage | lineage split unavailable | 243 | 2 (1 R, 1 S) | none | none | none | none | none |
| *Klebsiella pneumoniae* | doripenem | study | none-added | 243 | 2 (1 R, 1 S) | none | none | none | none | none |
| *Klebsiella pneumoniae* | doxycycline | lineage | lineage split unavailable | 141 | 5 (0 R, 5 S) | none | none | none | none | none |
| *Klebsiella pneumoniae* | ertapenem | lineage | draw | 952 | 533 (309 R, 224 S) | 0.012 | 0.006 | 0.008 | -0.012 | 0.036 |
| *Klebsiella pneumoniae* | ertapenem | study | study | 952 | 533 (309 R, 224 S) | 0.008 | -0.028 | -0.112 | 0.002 | -0.017 |
| *Klebsiella pneumoniae* | fosfomycin | lineage | lineage split unavailable | 294 | 2 (0 R, 2 S) | none | none | none | none | none |
| *Klebsiella pneumoniae* | imipenem | lineage | draw | 1,158 | 1,221 (485 R, 736 S) | 0.003 | -0.075 | 0.007 | -0.005 | 0.007 |
| *Klebsiella pneumoniae* | imipenem | study | study | 1,158 | 1,221 (485 R, 736 S) | -0.007 | -0.020 | -0.011 | 0.006 | 0.020 |
| *Klebsiella pneumoniae* | levofloxacin | lineage | draw | 871 | 1,239 (1168 R, 71 S) | -0.013 | 0.001 | 0.003 | -0.069 | -0.016 |
| *Klebsiella pneumoniae* | levofloxacin | study | study | 871 | 1,239 (1168 R, 71 S) | 0.008 | 0.000 | -0.010 | -0.001 | -0.004 |
| *Klebsiella pneumoniae* | nitrofurantoin | lineage | draw | 562 | 576 (558 R, 18 S) | 0.133 | 0.005 | -0.037 | -0.062 | 0.079 |
| *Klebsiella pneumoniae* | nitrofurantoin | study | study | 562 | 576 (558 R, 18 S) | 0.053 | -0.009 | 0.016 | -0.005 | 0.011 |
| *Klebsiella pneumoniae* | ofloxacin | lineage | lineage split unavailable | 71 | 0 (0 R, 0 S) | none | none | none | none | none |
| *Klebsiella pneumoniae* | polymyxin b | lineage | lineage split unavailable | 100 | 55 (55 R, 0 S) | none | none | none | none | none |
| *Klebsiella pneumoniae* | polymyxin b | study | none-added | 100 | 55 (55 R, 0 S) | none | none | none | none | none |
| *Klebsiella pneumoniae* | streptomycin | lineage | lineage split unavailable | 80 | 0 (0 R, 0 S) | none | none | none | none | none |
| *Klebsiella pneumoniae* | sulfamethoxazole/trimethoprim | lineage | draw | 1,546 | 1,735 (1196 R, 539 S) | -0.005 | -0.006 | 0.002 | -0.011 | -0.009 |
| *Klebsiella pneumoniae* | sulfamethoxazole/trimethoprim | study | study | 1,546 | 1,735 (1196 R, 539 S) | 0.025 | -0.017 | 0.008 | -0.013 | -0.017 |
| *Klebsiella pneumoniae* | tetracycline | lineage | draw | 741 | 1,047 (559 R, 488 S) | -0.001 | -0.043 | -0.019 | -0.056 | -0.027 |
| *Klebsiella pneumoniae* | ticarcillin/clavulanate | lineage | lineage split unavailable | 250 | 5 (0 R, 5 S) | none | none | none | none | none |
| *Klebsiella pneumoniae* | tigecycline | lineage | draw | 924 | 212 (12 R, 200 S) | 0.026 | 0.041 | 0.009 | 0.007 | 0.010 |
| *Klebsiella pneumoniae* | tigecycline | study | study | 924 | 212 (12 R, 200 S) | 0.024 | -0.009 | 0.029 | -0.023 | 0.009 |
| *Klebsiella pneumoniae* | tobramycin | lineage | draw | 1,230 | 1,036 (699 R, 337 S) | 0.009 | 0.019 | -0.001 | 0.001 | 0.005 |
| *Klebsiella pneumoniae* | tobramycin | study | study | 1,230 | 1,036 (699 R, 337 S) | -0.021 | -0.025 | 0.019 | 0.004 | 0.008 |
| *Klebsiella pneumoniae* | trimethoprim | lineage | draw | 258 | 15 (10 R, 5 S) | -0.012 | -0.012 | 0.011 | -0.032 | -0.021 |
| *Mycobacterium tuberculosis* | amikacin | lineage | draw | 1,813 | 1,959 (135 R, 1824 S) | -0.012 | 0.001 | 0.011 | -0.012 | -0.013 |
| *Mycobacterium tuberculosis* | amoxicillin | lineage | draw | 57 | 130 (0 R, 130 S) | 0.107 | 0.081 | -0.147 | 0.036 | 0.108 |
| *Mycobacterium tuberculosis* | capreomycin | lineage | draw | 670 | 1,157 (131 R, 1026 S) | -0.017 | -0.003 | 0.015 | -0.022 | -0.023 |
| *Mycobacterium tuberculosis* | ciprofloxacin | lineage | lineage split unavailable | 101 | 125 (0 R, 125 S) | none | none | none | none | none |
| *Mycobacterium tuberculosis* | cycloserine | lineage | draw | 401 | 383 (0 R, 383 S) | 0.033 | 0.003 | 0.034 | 0.037 | -0.036 |
| *Mycobacterium tuberculosis* | delamanid | lineage | draw | 1,248 | 774 (117 R, 657 S) | -0.032 | -0.031 | -0.033 | -0.014 | -0.026 |
| *Mycobacterium tuberculosis* | ethambutol | lineage | draw | 1,947 | 2,034 (1191 R, 843 S) | 0.011 | -0.019 | 0.023 | 0.017 | 0.002 |
| *Mycobacterium tuberculosis* | ethionamide | lineage | draw | 1,774 | 1,022 (377 R, 645 S) | 0.018 | 0.024 | 0.011 | -0.006 | -0.006 |
| *Mycobacterium tuberculosis* | kanamycin | lineage | draw | 1,800 | 1,135 (53 R, 1082 S) | -0.001 | -0.001 | 0.009 | 0.005 | -0.003 |
| *Mycobacterium tuberculosis* | levofloxacin | lineage | draw | 1,248 | 775 (774 R, 1 S) | 0.000 | -0.012 | -0.095 | -0.033 | -0.013 |
| *Mycobacterium tuberculosis* | linezolid | lineage | draw | 1,248 | 774 (1 R, 773 S) | -0.018 | 0.018 | -0.008 | -0.023 | -0.003 |
| *Mycobacterium tuberculosis* | moxifloxacin | lineage | draw | 1,524 | 1,542 (464 R, 1078 S) | -0.006 | -0.001 | -0.004 | 0.001 | -0.002 |
| *Mycobacterium tuberculosis* | nicotinamide | lineage | lineage split unavailable | 166 | 0 (0 R, 0 S) | none | none | none | none | none |
| *Mycobacterium tuberculosis* | nicotinamide | study | none-added | 166 | 0 (0 R, 0 S) | none | none | none | none | none |
| *Mycobacterium tuberculosis* | ofloxacin | lineage | draw | 643 | 1,105 (256 R, 849 S) | -0.004 | 0.000 | 0.011 | 0.002 | -0.001 |
| *Mycobacterium tuberculosis* | para-aminosalicylic acid | lineage | draw | 1,650 | 1,157 (6 R, 1151 S) | 0.002 | -0.001 | -0.017 | -0.002 | -0.012 |
| *Mycobacterium tuberculosis* | prothionamide | lineage | draw | 170 | 821 (323 R, 498 S) | -0.028 | -0.028 | -0.018 | 0.046 | 0.014 |
| *Mycobacterium tuberculosis* | pyrazinamide | lineage | draw | 492 | 1,055 (428 R, 627 S) | 0.093 | 0.099 | 0.130 | 0.073 | 0.112 |
| *Mycobacterium tuberculosis* | rifabutin | lineage | draw | 1,317 | 866 (361 R, 505 S) | 0.003 | 0.005 | 0.002 | 0.002 | 0.005 |
| *Mycobacterium tuberculosis* | rifampin | lineage | draw | 1,967 | 2,063 (1492 R, 571 S) | -0.002 | -0.003 | -0.013 | -0.026 | -0.005 |
| *Mycobacterium tuberculosis* | streptomycin | lineage | draw | 698 | 1,255 (989 R, 266 S) | 0.019 | 0.031 | 0.030 | -0.022 | -0.004 |
| *Neisseria gonorrhoeae* | cefixime | lineage | lineage split unavailable | 600 | 795 (0 R, 795 S) | none | none | none | none | none |
| *Neisseria gonorrhoeae* | ceftriaxone | lineage | lineage split unavailable | 1,926 | 1,990 (0 R, 1990 S) | none | none | none | none | none |
| *Neisseria gonorrhoeae* | ciprofloxacin | lineage | draw | 1,746 | 1,857 (381 R, 1476 S) | 0.004 | 0.003 | 0.056 | 0.005 | 0.004 |
| *Neisseria gonorrhoeae* | erythromycin | lineage | lineage split unavailable | 179 | 0 (0 R, 0 S) | none | none | none | none | none |
| *Neisseria gonorrhoeae* | penicillin | lineage | draw | 1,354 | 1,126 (0 R, 1126 S) | 0.058 | -0.056 | 0.023 | 0.091 | -0.001 |
| *Neisseria gonorrhoeae* | tetracycline | lineage | draw | 1,371 | 1,191 (0 R, 1191 S) | -0.069 | -0.023 | -0.079 | -0.018 | 0.002 |
| *Salmonella enterica* | amoxicillin/clavulanate | lineage | draw | 799 | 1,468 (1 R, 1467 S) | 0.000 | 0.001 | 0.000 | 0.004 | -0.006 |
| *Salmonella enterica* | ampicillin | lineage | draw | 1,997 | 2,091 (824 R, 1267 S) | 0.003 | -0.002 | 0.106 | 0.096 | 0.065 |
| *Salmonella enterica* | ampicillin/sulbactam | lineage | lineage split unavailable | 183 | 0 (0 R, 0 S) | none | none | none | none | none |
| *Salmonella enterica* | azithromycin | lineage | draw | 1,355 | 1,738 (0 R, 1738 S) | -0.002 | 0.014 | -0.008 | -0.004 | 0.003 |
| *Salmonella enterica* | cefazolin | lineage | lineage split unavailable | 185 | 0 (0 R, 0 S) | none | none | none | none | none |
| *Salmonella enterica* | cefotaxime | lineage | draw | 941 | 308 (85 R, 223 S) | -0.036 | 0.073 | -0.011 | 0.015 | 0.068 |
| *Salmonella enterica* | cefoxitin | lineage | draw | 1,080 | 1,521 (4 R, 1517 S) | 0.015 | 0.016 | -0.006 | 0.010 | 0.000 |
| *Salmonella enterica* | ceftazidime | lineage | lineage split unavailable | 783 | 149 (1 R, 148 S) | none | none | none | none | none |
| *Salmonella enterica* | ceftiofur | lineage | draw | 689 | 1,481 (3 R, 1478 S) | 0.003 | -0.022 | 0.020 | 0.008 | -0.012 |
| *Salmonella enterica* | ceftriaxone | lineage | draw | 1,155 | 1,593 (10 R, 1583 S) | 0.031 | 0.049 | 0.050 | 0.011 | 0.021 |
| *Salmonella enterica* | chloramphenicol | lineage | draw | 1,989 | 2,085 (402 R, 1683 S) | 0.008 | 0.024 | 0.021 | 0.036 | 0.058 |
| *Salmonella enterica* | ciprofloxacin | lineage | draw | 1,503 | 1,988 (130 R, 1858 S) | 0.054 | 0.008 | 0.025 | 0.016 | 0.034 |
| *Salmonella enterica* | gentamicin | lineage | draw | 1,675 | 1,974 (122 R, 1852 S) | 0.037 | -0.010 | 0.069 | 0.097 | 0.066 |
| *Salmonella enterica* | kanamycin | lineage | draw | 584 | 805 (10 R, 795 S) | 0.005 | -0.007 | 0.014 | -0.017 | -0.008 |
| *Salmonella enterica* | nalidixic acid | lineage | draw | 1,599 | 1,880 (206 R, 1674 S) | -0.018 | -0.020 | 0.019 | 0.019 | 0.030 |
| *Salmonella enterica* | pefloxacin | lineage | lineage split unavailable | 348 | 0 (0 R, 0 S) | none | none | none | none | none |
| *Salmonella enterica* | spectinomycin | lineage | draw | 396 | 386 (191 R, 195 S) | 0.050 | -0.023 | 0.034 | -0.015 | 0.041 |
| *Salmonella enterica* | streptomycin | lineage | draw | 1,269 | 1,956 (820 R, 1136 S) | 0.007 | 0.071 | 0.007 | 0.024 | -0.012 |
| *Salmonella enterica* | sulfamethoxazole | lineage | lineage split unavailable | 275 | 1 (0 R, 1 S) | none | none | none | none | none |
| *Salmonella enterica* | sulfamethoxazole/trimethoprim | lineage | draw | 1,145 | 1,582 (141 R, 1441 S) | 0.002 | -0.017 | -0.012 | -0.011 | 0.008 |
| *Salmonella enterica* | sulfathiazole | lineage | lineage split unavailable | 225 | 158 (125 R, 33 S) | none | none | none | none | none |
| *Salmonella enterica* | sulfisoxazole | lineage | draw | 544 | 1,439 (298 R, 1141 S) | 0.065 | 0.055 | 0.071 | 0.051 | -0.006 |
| *Salmonella enterica* | sulfisoxazole | study | study | 544 | 1,439 (298 R, 1141 S) | 0.000 | -0.003 | -0.006 | -0.007 | -0.014 |
| *Salmonella enterica* | tetracycline | lineage | draw | 1,869 | 2,035 (1016 R, 1019 S) | 0.007 | 0.010 | 0.117 | 0.041 | 0.082 |
| *Salmonella enterica* | trimethoprim | lineage | draw | 634 | 443 (62 R, 381 S) | -0.004 | 0.016 | 0.071 | -0.079 | -0.057 |
| *Staphylococcus aureus* | cefoxitin | lineage | lineage split unavailable | 1,407 | 612 (519 R, 93 S) | none | none | none | none | none |
| *Staphylococcus aureus* | chloramphenicol | lineage | lineage split unavailable | 1,038 | 455 (0 R, 455 S) | none | none | none | none | none |
| *Staphylococcus aureus* | ciprofloxacin | lineage | lineage split unavailable | 1,624 | 1,312 (583 R, 729 S) | none | none | none | none | none |
| *Staphylococcus aureus* | clindamycin | lineage | draw | 1,543 | 386 (187 R, 199 S) | 0.073 | 0.099 | 0.059 | 0.020 | 0.084 |
| *Staphylococcus aureus* | clindamycin | study | study | 1,543 | 386 (187 R, 199 S) | 0.000 | -0.015 | -0.001 | -0.008 | -0.018 |
| *Staphylococcus aureus* | daptomycin | lineage | lineage split unavailable | 904 | 14 (0 R, 14 S) | none | none | none | none | none |
| *Staphylococcus aureus* | erythromycin | lineage | draw | 1,957 | 1,382 (505 R, 877 S) | 0.034 | 0.073 | 0.098 | 0.029 | 0.067 |
| *Staphylococcus aureus* | erythromycin | study | study | 1,957 | 1,382 (505 R, 877 S) | 0.000 | 0.061 | 0.000 | 0.007 | -0.002 |
| *Staphylococcus aureus* | fusidic acid | lineage | lineage split unavailable | 1,101 | 1,246 (60 R, 1186 S) | none | none | none | none | none |
| *Staphylococcus aureus* | gentamicin | lineage | lineage split unavailable | 1,913 | 1,385 (8 R, 1377 S) | none | none | none | none | none |
| *Staphylococcus aureus* | mupirocin | lineage | lineage split unavailable | 736 | 865 (0 R, 865 S) | none | none | none | none | none |
| *Staphylococcus aureus* | penicillin | lineage | draw | 1,508 | 811 (811 R, 0 S) | 0.000 | -0.001 | -0.002 | -0.001 | -0.004 |
| *Staphylococcus aureus* | penicillin | study | study | 1,508 | 811 (811 R, 0 S) | 0.000 | 0.008 | 0.005 | -0.007 | 0.006 |
| *Staphylococcus aureus* | rifampin | lineage | lineage split unavailable | 931 | 1,339 (0 R, 1339 S) | none | none | none | none | none |
| *Staphylococcus aureus* | sulfamethoxazole/trimethoprim | lineage | lineage split unavailable | 1,389 | 147 (6 R, 141 S) | none | none | none | none | none |
| *Staphylococcus aureus* | sulfamethoxazole/trimethoprim | study | study | 1,389 | 147 (6 R, 141 S) | 0.000 | 0.000 | 0.000 | 0.019 | 0.010 |
| *Staphylococcus aureus* | tetracycline | lineage | lineage split unavailable | 1,647 | 1,339 (92 R, 1247 S) | none | none | none | none | none |
| *Streptococcus pneumoniae* | azithromycin | lineage | none-added | 1,297 | 0 (0 R, 0 S) | none | none | none | none | none |
| *Streptococcus pneumoniae* | cefuroxime | lineage | none-added | 1,732 | 0 (0 R, 0 S) | none | none | none | none | none |
| *Streptococcus pneumoniae* | chloramphenicol | lineage | lineage split unavailable | 1,859 | 1,767 (510 R, 1257 S) | none | none | none | none | none |
| *Streptococcus pneumoniae* | clindamycin | lineage | lineage split unavailable | 1,894 | 148 (0 R, 148 S) | none | none | none | none | none |
| *Streptococcus pneumoniae* | erythromycin | lineage | draw | 1,991 | 2,001 (1358 R, 643 S) | 0.011 | 0.022 | 0.009 | 0.007 | -0.005 |
| *Streptococcus pneumoniae* | meropenem | lineage | lineage split unavailable | 1,808 | 13 (0 R, 13 S) | none | none | none | none | none |
| *Streptococcus pneumoniae* | sulfamethoxazole/trimethoprim | lineage | draw | 1,777 | 1,610 (1146 R, 464 S) | 0.041 | 0.036 | 0.023 | -0.004 | -0.005 |
| *Streptococcus pneumoniae* | tetracycline | lineage | draw | 1,971 | 1,993 (1418 R, 575 S) | -0.040 | -0.053 | 0.000 | -0.002 | -0.002 |
