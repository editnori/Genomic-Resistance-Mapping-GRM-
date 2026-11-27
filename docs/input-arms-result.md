# Input-arm result

The release fits its five learners on the presence of 31-mers in each genome (K31). This page asks whether five other inputs separate the classes better: unitigs, K31 plus k-mer counts, K31 plus amino-acid k-mers of called genes, the AMRFinderPlus program's calls alone, and K31 plus the program's calls. It reads block C3 of the assumption plan (`docs/assumption-plan.md`), with the method of `docs/input-arms.md`.

The reading's rules were fixed before any C3 score was read (`docs/assumption-plan-corrections.md`, entries "C3 design" and "C3 arms"), and the scripts that read it at 14:56:32Z on 23 Nov 2026. The count was read once, from 14:59:17Z to 15:00:54Z on 23 Nov 2026, with no other option, subset or measure.

Every split here is a retrospective split inside one public database. A split by lineage or by study is never external validation. The margin of 0.02 is an adapted research screen. It is no applicable FDA validation criterion.

GRM's known-gene rule is GRM's own reader of the AMRFinderPlus catalog. "The program" on this page is the AMRFinderPlus program itself, a separate thing.

## Reading

**Result.** By ruling 15, one arm of the five meets conditions 1 to 5: the AMRFinderPlus program's calls alone, on 203 of A's 204 datasets in 21 species. Four of its five learners carry the change, `cart`, `rf`, `gbt` and `lasso`, and `scm` does not. The twin meets condition 1 in 0 of 5 learners and condition 2 in 0 of 5. No learner's recall interval lies wholly below -0.02. Each of the four that carry has a study median of its lineage median's sign. Condition 6 is the maintainer's, so nothing changes until the maintainer rules.

**K31 plus amino-acid k-mers.** 2 of 5 learners carry, `rf` and `gbt`, fewer than 3, so conditions 1 and 2 are not met. Better beats worse after Holm for all five learners, and the median difference's interval lies wholly above 0 for `rf` and `gbt` only.

**The other three arms.** Unitigs, K31 plus k-mer counts and K31 plus the AMRFinderPlus program's calls: 0 of 5 learners carry, so conditions 1 and 2 are not met for any of them.

**Ruling 15.** An input arm touches all five learners, so the plan's conditions are read per learner. A learner carries the change when it meets condition 1 and condition 2 both, and at least 3 of the 5 learners must carry it.

Condition 3 follows the plan's words, the twin meets neither 1 nor 2: it holds when the twin meets 1 in fewer than 3 of 5 learners and 2 in fewer than 3 of 5, each read apart. Condition 4 holds for all five learners. Condition 5 holds for each learner that carries. Condition 6 is the maintainer's.

| Condition | Reads by ruling 15, which decides |
|---|---|
| 1 and 2 | a learner carries when its lineage median difference's interval lies wholly above 0 and better beats worse by the one-sided sign test at p < 0.05 after Holm over 30; at least 3 of 5 carry |
| 3 | the twin meets 1 in fewer than 3 of 5 learners and 2 in fewer than 3 of 5, each read apart |
| 4 | no learner's recall interval lies wholly below -0.02, for all five |
| 5 | each learner that carries has a study median of its lineage median's sign |
| 6 | the maintainer rules the change |

| Arm | Condition | Reads | Read | Learners |
|---|---|---|---|---|
| unitigs | 1 and 2 | a learner carries when its lineage median difference's interval lies wholly above 0 and better beats worse by the one-sided sign test at p < 0.05 after Holm over 30; at least 3 of 5 carry | not met | 0 |
| unitigs | 3 | the twin meets 1 in fewer than 3 of 5 learners and 2 in fewer than 3 of 5, each read apart (learners: the larger of the twin's two counts) | met | 0 |
| unitigs | 4 | no learner's recall interval lies wholly below -0.02, for all five | met | 5 |
| unitigs | 5 | each learner that carries has a study median of its lineage median's sign | not met | 0 |
| unitigs | 6 | the maintainer rules the change | the maintainer's |  |
| K31 plus k-mer counts | 1 and 2 | a learner carries when its lineage median difference's interval lies wholly above 0 and better beats worse by the one-sided sign test at p < 0.05 after Holm over 30; at least 3 of 5 carry | not met | 0 |
| K31 plus k-mer counts | 3 | the twin meets 1 in fewer than 3 of 5 learners and 2 in fewer than 3 of 5, each read apart (learners: the larger of the twin's two counts) | met | 0 |
| K31 plus k-mer counts | 4 | no learner's recall interval lies wholly below -0.02, for all five | met | 5 |
| K31 plus k-mer counts | 5 | each learner that carries has a study median of its lineage median's sign | not met | 0 |
| K31 plus k-mer counts | 6 | the maintainer rules the change | the maintainer's |  |
| K31 plus amino-acid k-mers | 1 and 2 | a learner carries when its lineage median difference's interval lies wholly above 0 and better beats worse by the one-sided sign test at p < 0.05 after Holm over 30; at least 3 of 5 carry | not met | 2 |
| K31 plus amino-acid k-mers | 3 | the twin meets 1 in fewer than 3 of 5 learners and 2 in fewer than 3 of 5, each read apart (learners: the larger of the twin's two counts) | met | 0 |
| K31 plus amino-acid k-mers | 4 | no learner's recall interval lies wholly below -0.02, for all five | met | 5 |
| K31 plus amino-acid k-mers | 5 | each learner that carries has a study median of its lineage median's sign | met | 2 |
| K31 plus amino-acid k-mers | 6 | the maintainer rules the change | the maintainer's |  |
| the AMRFinderPlus program's calls alone | 1 and 2 | a learner carries when its lineage median difference's interval lies wholly above 0 and better beats worse by the one-sided sign test at p < 0.05 after Holm over 30; at least 3 of 5 carry | met | 4 |
| the AMRFinderPlus program's calls alone | 3 | the twin meets 1 in fewer than 3 of 5 learners and 2 in fewer than 3 of 5, each read apart (learners: the larger of the twin's two counts) | met | 0 |
| the AMRFinderPlus program's calls alone | 4 | no learner's recall interval lies wholly below -0.02, for all five | met | 5 |
| the AMRFinderPlus program's calls alone | 5 | each learner that carries has a study median of its lineage median's sign | met | 4 |
| the AMRFinderPlus program's calls alone | 6 | the maintainer rules the change | the maintainer's |  |
| K31 plus the AMRFinderPlus program's calls | 1 and 2 | a learner carries when its lineage median difference's interval lies wholly above 0 and better beats worse by the one-sided sign test at p < 0.05 after Holm over 30; at least 3 of 5 carry | not met | 0 |
| K31 plus the AMRFinderPlus program's calls | 3 | the twin meets 1 in fewer than 3 of 5 learners and 2 in fewer than 3 of 5, each read apart (learners: the larger of the twin's two counts) | met | 0 |
| K31 plus the AMRFinderPlus program's calls | 4 | no learner's recall interval lies wholly below -0.02, for all five | met | 5 |
| K31 plus the AMRFinderPlus program's calls | 5 | each learner that carries has a study median of its lineage median's sign | not met | 0 |
| K31 plus the AMRFinderPlus program's calls | 6 | the maintainer rules the change | the maintainer's |  |

**The stricter reading.** Beside ruling 15 stands a stricter reading: condition 1 for at least 3 of 5 learners, 2 and 5 for all five, 3 the twin meeting neither, 4 for all five. It decides nothing, and it stands under Beside it.

"No difference shown" is the rule's reading and no statement that two inputs are equal.

The release keeps its input whatever the count says, until the maintainer rules. Nothing changes in the release or on the site.

## Data

**Datasets.** The first run is L's datasets less those in M, which wait for the maintainer's label rule: 233 datasets in 23 species, 55 of them in S with a study split (`docs/input-arms.md`, Datasets). The two program arms run on A's first run, 204 datasets in 21 species, 54 of them in S: database 2026-08-07.1 of the program has no organism option for *Mycobacterium tuberculosis* or *Staphylococcus haemolyticus*. The section Datasets names every species.

**Entered.** On the lineage split, 232 of L's 233 datasets in 23 species enter the count of the default, the twin and the three matrix arms, and 203 of A's 204 in 21 species enter the count of the two program arms. The one left out is named under Not scored. On the study split, 55 of 55 enter on L and 54 of 54 on A. No species was set aside. *Vibrio cholerae*'s run left a `FAILED` file at its build at 07:29:21Z, older than its `DONE` at 07:49:38Z; the check of the run named that file and set it aside, and the species entered.

**Arms.** Each arm is fitted in place of K31 with the release's five learners and scored against K31 refitted in the same run, on the same folds, with the same binary (ruling 17). Each arm is counted on its own datasets only: the three matrix arms, the default and the twin on L's first run, the two program arms on A's.

**Draws.** The default and each arm come from the same runs at fold seeds 42, 43 and 44 on the lineage split. A dataset's value is the mean of its draws. The twin is the default at seeds 45, 46 and 47. It is never read as an arm.

**Study split.** Each dataset in S gets one study run per arm and for the default, on its reviewed study plan, at the release's seeds, with no twin. A study row is never counted on the lineage split.

**Refused study runs.** Only a study run can be refused because an outer fold's training genomes give no candidate feature. By ruling 16, an arm's study run refused this way, where the default's is not, is scored as the constant model on the default's folds: every test genome of a fold gets the majority class of the fold's training genomes, R on a tie. That dataset is kept, enters the count, and is marked in `constant`.

The count with those datasets left out stands under Beside it and decides nothing.

No study run was refused for no candidate: 0 for each arm (`page-arms.tsv`), and the check of the run found 0 refusal records. No dataset was scored as the constant model.

**Not scored.** A dataset with a draw not scored is not scored. A missing draw is never read as zero.

*Klebsiella pneumoniae* with polymyxin b leaves the lineage count, so 232 of 233 enter it on L and 203 of 204 on A. Its lineage split did not run at fold seed 43, for the default and for each of the five arms, nor at the twin's seed 46, and a dataset with a draw not scored is not scored. Each of those 7 results records: `invalid lineage outer evaluation: fold 4 validation lacks one or both classes`. Its study split ran and is counted.

**Runs.** Each species ran on a machine of its own, one run with one binary for the default, the twin and all five arms, and wrote its scorer's tables there. The count is the scorer's count over those tables, copied as the runs wrote them.

## Count

The primary measure is each learner's pooled lineage balanced accuracy on the outer test folds. A dataset is better, worse or the same by more than 0.02, compared at 12 decimals. Every difference is the arm's value less the default's.

The median difference comes with a 95% interval from 4,000 species resamples. The one-sided sign test of better against worse runs with Holm over 30 lineage pairs on the primary measure: 6 arm rows (unitigs, K31 plus counts, K31 plus amino-acid k-mers, the program alone, K31 plus the program and the twin) by 5 learners. An arm withdrawn before it ran keeps its five pairs at p = 1.

| Arm | Learner | Better | Worse | Same | Median | Interval | Sign p | Holm p | Constant | Condition 1 | Condition 2 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| unitigs | `scm` | 19 of 232 | 29 of 232 | 184 of 232 | 0.000 | -0.001 to 0.000 | 0.944 | 1.000 | 0 | not met | not met |
| unitigs | `cart` | 24 of 232 | 19 of 232 | 189 of 232 | 0.000 | 0.000 to 0.000 | 0.271 | 1.000 | 0 | not met | not met |
| unitigs | `rf` | 10 of 232 | 12 of 232 | 210 of 232 | 0.000 | -0.001 to 0.000 | 0.738 | 1.000 | 0 | not met | not met |
| unitigs | `gbt` | 34 of 232 | 15 of 232 | 183 of 232 | 0.000 | 0.000 to 0.001 | 0.005 | 0.099 | 0 | not met | not met |
| unitigs | `lasso` | 12 of 232 | 11 of 232 | 209 of 232 | 0.000 | -0.001 to 0.001 | 0.500 | 1.000 | 0 | not met | not met |
| K31 plus k-mer counts | `scm` | 19 of 232 | 18 of 232 | 195 of 232 | 0.000 | 0.000 to 0.000 | 0.500 | 1.000 | 0 | not met | not met |
| K31 plus k-mer counts | `cart` | 27 of 232 | 23 of 232 | 182 of 232 | 0.000 | 0.000 to 0.000 | 0.336 | 1.000 | 0 | not met | not met |
| K31 plus k-mer counts | `rf` | 18 of 232 | 12 of 232 | 202 of 232 | 0.000 | 0.000 to 0.000 | 0.181 | 1.000 | 0 | not met | not met |
| K31 plus k-mer counts | `gbt` | 35 of 232 | 18 of 232 | 179 of 232 | 0.001 | 0.000 to 0.002 | 0.014 | 0.270 | 0 | met | not met |
| K31 plus k-mer counts | `lasso` | 19 of 232 | 11 of 232 | 202 of 232 | 0.000 | 0.000 to 0.000 | 0.100 | 1.000 | 0 | not met | not met |
| K31 plus amino-acid k-mers | `scm` | 48 of 232 | 23 of 232 | 161 of 232 | 0.001 | 0.000 to 0.003 | 0.002 | 0.045 | 0 | not met | met |
| K31 plus amino-acid k-mers | `cart` | 54 of 232 | 25 of 232 | 153 of 232 | 0.000 | 0.000 to 0.003 | < 0.001 | 0.017 | 0 | not met | met |
| K31 plus amino-acid k-mers | `rf` | 39 of 232 | 8 of 232 | 185 of 232 | 0.002 | 0.001 to 0.005 | < 0.001 | < 0.001 | 0 | met | met |
| K31 plus amino-acid k-mers | `gbt` | 51 of 232 | 16 of 232 | 165 of 232 | 0.003 | 0.001 to 0.005 | < 0.001 | < 0.001 | 0 | met | met |
| K31 plus amino-acid k-mers | `lasso` | 46 of 232 | 16 of 232 | 170 of 232 | 0.001 | 0.000 to 0.004 | < 0.001 | 0.002 | 0 | not met | met |
| the AMRFinderPlus program's calls alone | `scm` | 70 of 203 | 53 of 203 | 80 of 203 | 0.005 | 0.000 to 0.013 | 0.074 | 1.000 | 0 | not met | not met |
| the AMRFinderPlus program's calls alone | `cart` | 85 of 203 | 42 of 203 | 76 of 203 | 0.009 | 0.005 to 0.018 | < 0.001 | 0.002 | 0 | met | met |
| the AMRFinderPlus program's calls alone | `rf` | 121 of 203 | 28 of 203 | 54 of 203 | 0.036 | 0.025 to 0.053 | < 0.001 | < 0.001 | 0 | met | met |
| the AMRFinderPlus program's calls alone | `gbt` | 105 of 203 | 26 of 203 | 72 of 203 | 0.021 | 0.013 to 0.030 | < 0.001 | < 0.001 | 0 | met | met |
| the AMRFinderPlus program's calls alone | `lasso` | 102 of 203 | 32 of 203 | 69 of 203 | 0.020 | 0.011 to 0.027 | < 0.001 | < 0.001 | 0 | met | met |
| K31 plus the AMRFinderPlus program's calls | `scm` | 6 of 203 | 9 of 203 | 188 of 203 | 0.000 | 0.000 to 0.000 | 0.849 | 1.000 | 0 | not met | not met |
| K31 plus the AMRFinderPlus program's calls | `cart` | 5 of 203 | 5 of 203 | 193 of 203 | 0.000 | 0.000 to 0.000 | 0.623 | 1.000 | 0 | not met | not met |
| K31 plus the AMRFinderPlus program's calls | `rf` | 4 of 203 | 1 of 203 | 198 of 203 | 0.000 | 0.000 to 0.000 | 0.188 | 1.000 | 0 | not met | not met |
| K31 plus the AMRFinderPlus program's calls | `gbt` | 10 of 203 | 10 of 203 | 183 of 203 | 0.000 | 0.000 to 0.000 | 0.588 | 1.000 | 0 | not met | not met |
| K31 plus the AMRFinderPlus program's calls | `lasso` | 1 of 203 | 1 of 203 | 201 of 203 | 0.000 | 0.000 to 0.000 | 0.750 | 1.000 | 0 | not met | not met |
| the twin, the default at seeds 45, 46 and 47 | `scm` | 29 of 232 | 27 of 232 | 176 of 232 | 0.000 | 0.000 to 0.001 | 0.447 | 1.000 | 0 | not met | not met |
| the twin, the default at seeds 45, 46 and 47 | `cart` | 27 of 232 | 30 of 232 | 175 of 232 | 0.000 | 0.000 to 0.000 | 0.702 | 1.000 | 0 | not met | not met |
| the twin, the default at seeds 45, 46 and 47 | `rf` | 18 of 232 | 22 of 232 | 192 of 232 | -0.001 | -0.003 to 0.000 | 0.785 | 1.000 | 0 | not met | not met |
| the twin, the default at seeds 45, 46 and 47 | `gbt` | 40 of 232 | 33 of 232 | 159 of 232 | 0.000 | 0.000 to 0.001 | 0.241 | 1.000 | 0 | not met | not met |
| the twin, the default at seeds 45, 46 and 47 | `lasso` | 28 of 232 | 19 of 232 | 185 of 232 | 0.000 | -0.001 to 0.000 | 0.121 | 1.000 | 0 | not met | not met |

## Per arm

**Unitigs.** One column per distinct unitig pattern of the graph over every 31-mer of the species' genomes, in place of K31. It runs on L.

**K31 plus k-mer counts.** [count ≥ 2] of each 31-mer in a genome, so a genome has the row when the 31-mer occurs in it at least twice, then the K31 rows less those whose cohort pattern a count row has. It runs on L.

**K31 plus amino-acid k-mers.** Amino-acid 10-mers of the genes Prodigal 2.6.3 calls, then the K31 rows less those whose cohort pattern an amino-acid row has. It runs on L.

**The AMRFinderPlus program's calls alone.** One column per element symbol of the program's Type AMR rows, as the program's table names it. It runs on A. The program is AMRFinderPlus 4.2.7 with database 2026-08-07.1 and `--plus`.

**K31 plus the AMRFinderPlus program's calls.** The program's columns, then the dataset's K31 rows less those whose cohort pattern a program column has. It runs on A.

A program-call arm would enter only with the program's version pinned, and the maintainer would rule how GRM may name its calls.

**Reading per arm.** Learners that carry the change by ruling 15: unitigs 0 of 5; K31 plus k-mer counts 0 of 5; K31 plus amino-acid k-mers 2 of 5, `rf` and `gbt`; the AMRFinderPlus program's calls alone 4 of 5, `cart`, `rf`, `gbt` and `lasso`; K31 plus the AMRFinderPlus program's calls 0 of 5. Every arm ran, and none was withdrawn.

| Arm | Set | State | Datasets | Species | Study datasets | Study runs refused |
|---|---|---|---|---|---|---|
| K31, the default | L | ran | 233 | 23 | 55 | 0 |
| unitigs | L | ran | 233 | 23 | 55 | 0 |
| K31 plus k-mer counts | L | ran | 233 | 23 | 55 | 0 |
| K31 plus amino-acid k-mers | L | ran | 233 | 23 | 55 | 0 |
| the AMRFinderPlus program's calls alone | A | ran | 204 | 21 | 54 | 0 |
| K31 plus the AMRFinderPlus program's calls | A | ran | 204 | 21 | 54 | 0 |
| the twin, the default at seeds 45, 46 and 47 | L | ran | 233 | 23 | 0 | 0 |

**Builds.** Each species built the three matrix arms' inputs once, and each record holds counts only.

| Species | Build | Genomes | Values summed | Patterns | Proteins | Unitigs | Pairs held in part |
|---|---|---|---|---|---|---|---|
| *Acinetobacter baumannii* | aa | 1,702 | 1,901,796,518 | 440,154 | 6,552,874 |  |  |
| *Acinetobacter baumannii* | counts | 1,702 | 42,982,498 | 78,652 |  |  |  |
| *Acinetobacter baumannii* | unitigs | 1,702 | 1,084,530,272 | 1,569,889 |  | 2,899,863 | 64,818 |
| *Campylobacter coli* | aa | 391 | 204,299,294 | 86,977 | 714,742 |  |  |
| *Campylobacter coli* | counts | 391 | 2,092,278 | 11,742 |  |  |  |
| *Campylobacter coli* | unitigs | 391 | 46,390,988 | 222,241 |  | 368,389 | 13,310 |
| *Campylobacter jejuni* | aa | 508 | 261,575,942 | 90,611 | 912,294 |  |  |
| *Campylobacter jejuni* | counts | 508 | 2,664,633 | 10,441 |  |  |  |
| *Campylobacter jejuni* | unitigs | 508 | 78,715,651 | 264,997 |  | 456,474 | 10,160 |
| *Clostridioides difficile* | aa | 204 | 225,585,493 | 53,770 | 784,319 |  |  |
| *Clostridioides difficile* | counts | 204 | 5,341,957 | 18,475 |  |  |  |
| *Clostridioides difficile* | unitigs | 204 | 55,508,819 | 200,221 |  | 703,397 | 16,367 |
| *Corynebacterium diphtheriae* | aa | 243 | 169,730,671 | 156,198 | 559,282 |  |  |
| *Corynebacterium diphtheriae* | counts | 243 | 989,780 | 7,242 |  |  |  |
| *Corynebacterium diphtheriae* | unitigs | 243 | 66,071,381 | 556,020 |  | 856,628 | 11,248 |
| *Enterobacter cloacae* | aa | 270 | 383,185,052 | 338,250 | 1,269,052 |  |  |
| *Enterobacter cloacae* | counts | 270 | 3,604,887 | 16,244 |  |  |  |
| *Enterobacter cloacae* | unitigs | 270 | 163,239,569 | 1,222,404 |  | 3,927,570 | 15,840 |
| *Enterococcus faecium* | aa | 2,000 | 1,581,423,541 | 315,529 | 5,827,108 |  |  |
| *Enterococcus faecium* | counts | 2,000 | 34,168,097 | 61,636 |  |  |  |
| *Enterococcus faecium* | unitigs | 2,000 | 712,173,077 | 814,056 |  | 1,493,776 | 100,988 |
| *Escherichia coli* | aa | 2,000 | 2,866,173,068 | 1,008,020 | 9,795,417 |  |  |
| *Escherichia coli* | counts | 2,000 | 53,447,222 | 189,522 |  |  |  |
| *Escherichia coli* | unitigs | 2,000 | 2,089,259,406 | 3,351,087 |  | 5,493,343 | 150,437 |
| *Klebsiella aerogenes* | aa | 93 | 138,711,652 | 90,641 | 455,919 |  |  |
| *Klebsiella aerogenes* | counts | 93 | 1,643,521 | 7,596 |  |  |  |
| *Klebsiella aerogenes* | unitigs | 93 | 43,935,938 | 577,629 |  | 1,370,388 | 803 |
| *Klebsiella oxytoca* | aa | 112 | 192,731,465 | 84,452 | 639,505 |  |  |
| *Klebsiella oxytoca* | counts | 112 | 1,648,648 | 7,550 |  |  |  |
| *Klebsiella oxytoca* | unitigs | 112 | 45,279,878 | 311,720 |  | 1,749,258 | 10,252 |
| *Klebsiella pneumoniae* | aa | 2,000 | 3,129,460,178 | 765,004 | 10,520,926 |  |  |
| *Klebsiella pneumoniae* | counts | 2,000 | 35,994,549 | 88,228 |  |  |  |
| *Klebsiella pneumoniae* | unitigs | 2,000 | 2,128,201,382 | 2,584,874 |  | 4,426,252 | 96,938 |
| *Klebsiella quasipneumoniae* | aa | 109 | 169,510,077 | 75,921 | 566,993 |  |  |
| *Klebsiella quasipneumoniae* | counts | 109 | 1,319,703 | 6,694 |  |  |  |
| *Klebsiella quasipneumoniae* | unitigs | 109 | 43,939,445 | 400,476 |  | 1,158,295 | 10,351 |
| *Mycobacterium tuberculosis* | aa | 2,000 | 2,505,160,281 | 159,567 | 8,167,681 |  |  |
| *Mycobacterium tuberculosis* | counts | 2,000 | 33,066,659 | 21,081 |  |  |  |
| *Mycobacterium tuberculosis* | unitigs | 2,000 | 560,601,791 | 275,719 |  | 558,497 | 170,812 |
| *Neisseria gonorrhoeae* | aa | 2,000 | 1,134,678,586 | 103,055 | 4,255,825 |  |  |
| *Neisseria gonorrhoeae* | counts | 2,000 | 30,387,471 | 26,882 |  |  |  |
| *Neisseria gonorrhoeae* | unitigs | 2,000 | 303,566,186 | 202,190 |  | 353,872 | 193,847 |
| *Pseudomonas aeruginosa* | aa | 1,421 | 2,735,703,912 | 752,488 | 8,885,356 |  |  |
| *Pseudomonas aeruginosa* | counts | 1,421 | 33,226,726 | 82,507 |  |  |  |
| *Pseudomonas aeruginosa* | unitigs | 1,421 | 1,772,395,920 | 2,401,921 |  | 4,164,312 | 167,894 |
| *Salmonella enterica* | aa | 2,000 | 2,749,730,989 | 277,069 | 9,347,836 |  |  |
| *Salmonella enterica* | counts | 2,000 | 49,798,931 | 51,224 |  |  |  |
| *Salmonella enterica* | unitigs | 2,000 | 1,258,426,533 | 1,098,067 |  | 1,941,439 | 114,177 |
| *Shigella flexneri* | aa | 362 | 442,703,705 | 17,922 | 1,631,095 |  |  |
| *Shigella flexneri* | counts | 362 | 3,172,906 | 3,280 |  |  |  |
| *Shigella flexneri* | unitigs | 362 | 39,193,976 | 29,177 |  | 212,961 | 196,765 |
| *Shigella sonnei* | aa | 1,984 | 2,507,684,273 | 153,843 | 9,078,208 |  |  |
| *Shigella sonnei* | counts | 1,984 | 23,754,653 | 27,241 |  |  |  |
| *Shigella sonnei* | unitigs | 1,984 | 661,887,852 | 473,032 |  | 1,070,586 | 733,188 |
| *Staphylococcus aureus* | aa | 2,000 | 1,503,719,561 | 260,754 | 5,218,222 |  |  |
| *Staphylococcus aureus* | counts | 2,000 | 19,680,858 | 32,087 |  |  |  |
| *Staphylococcus aureus* | unitigs | 2,000 | 885,620,251 | 765,269 |  | 1,291,757 | 41,649 |
| *Staphylococcus haemolyticus* | aa | 126 | 87,185,931 | 21,642 | 316,877 |  |  |
| *Staphylococcus haemolyticus* | counts | 126 | 1,061,879 | 4,204 |  |  |  |
| *Staphylococcus haemolyticus* | unitigs | 126 | 7,685,210 | 44,428 |  | 170,422 | 2,779 |
| *Streptococcus agalactiae* | aa | 141 | 87,778,323 | 5,106 | 304,148 |  |  |
| *Streptococcus agalactiae* | counts | 141 | 474,623 | 1,115 |  |  |  |
| *Streptococcus agalactiae* | unitigs | 141 | 6,131,463 | 8,648 |  | 85,039 | 1,019 |
| *Streptococcus pneumoniae* | aa | 2,000 | 1,126,211,172 | 202,021 | 4,126,155 |  |  |
| *Streptococcus pneumoniae* | counts | 2,000 | 31,402,846 | 35,617 |  |  |  |
| *Streptococcus pneumoniae* | unitigs | 2,000 | 566,708,472 | 524,146 |  | 814,228 | 97,135 |
| *Vibrio cholerae* | aa | 91 | 102,692,050 | 36,212 | 327,075 |  |  |
| *Vibrio cholerae* | counts | 91 | 1,227,630 | 6,190 |  |  |  |
| *Vibrio cholerae* | unitigs | 91 | 26,676,162 | 235,254 |  | 756,851 | 3,963 |

## Beside it

These tables are scored beside the reading. They decide nothing.

**The stricter reading.** No arm meets it: the AMRFinderPlus program's calls alone meet condition 2 in 4 of 5 learners, not all five, and K31 plus amino-acid k-mers meet condition 1 in 2 of 5, fewer than 3.

| Arm | Condition | Reads | Read | Learners |
|---|---|---|---|---|
| unitigs | 1 | the lineage median difference's interval lies wholly above 0, for at least 3 of 5 learners | not met | 0 |
| unitigs | 2 | better beats worse by the one-sided sign test at p < 0.05 after Holm over 30, for all five | not met | 0 |
| unitigs | 3 | the twin meets neither 1 nor 2 as this reading reads them | met |  |
| unitigs | 4 | no learner's recall interval lies wholly below -0.02, for all five | met | 5 |
| unitigs | 5 | each of the five learners has a study median of its lineage median's sign | not met | 0 |
| unitigs | 6 | the maintainer rules the change | the maintainer's |  |
| K31 plus k-mer counts | 1 | the lineage median difference's interval lies wholly above 0, for at least 3 of 5 learners | not met | 1 |
| K31 plus k-mer counts | 2 | better beats worse by the one-sided sign test at p < 0.05 after Holm over 30, for all five | not met | 0 |
| K31 plus k-mer counts | 3 | the twin meets neither 1 nor 2 as this reading reads them | met |  |
| K31 plus k-mer counts | 4 | no learner's recall interval lies wholly below -0.02, for all five | met | 5 |
| K31 plus k-mer counts | 5 | each of the five learners has a study median of its lineage median's sign | not met | 1 |
| K31 plus k-mer counts | 6 | the maintainer rules the change | the maintainer's |  |
| K31 plus amino-acid k-mers | 1 | the lineage median difference's interval lies wholly above 0, for at least 3 of 5 learners | not met | 2 |
| K31 plus amino-acid k-mers | 2 | better beats worse by the one-sided sign test at p < 0.05 after Holm over 30, for all five | met | 5 |
| K31 plus amino-acid k-mers | 3 | the twin meets neither 1 nor 2 as this reading reads them | met |  |
| K31 plus amino-acid k-mers | 4 | no learner's recall interval lies wholly below -0.02, for all five | met | 5 |
| K31 plus amino-acid k-mers | 5 | each of the five learners has a study median of its lineage median's sign | met | 5 |
| K31 plus amino-acid k-mers | 6 | the maintainer rules the change | the maintainer's |  |
| the AMRFinderPlus program's calls alone | 1 | the lineage median difference's interval lies wholly above 0, for at least 3 of 5 learners | met | 4 |
| the AMRFinderPlus program's calls alone | 2 | better beats worse by the one-sided sign test at p < 0.05 after Holm over 30, for all five | not met | 4 |
| the AMRFinderPlus program's calls alone | 3 | the twin meets neither 1 nor 2 as this reading reads them | met |  |
| the AMRFinderPlus program's calls alone | 4 | no learner's recall interval lies wholly below -0.02, for all five | met | 5 |
| the AMRFinderPlus program's calls alone | 5 | each of the five learners has a study median of its lineage median's sign | met | 5 |
| the AMRFinderPlus program's calls alone | 6 | the maintainer rules the change | the maintainer's |  |
| K31 plus the AMRFinderPlus program's calls | 1 | the lineage median difference's interval lies wholly above 0, for at least 3 of 5 learners | not met | 0 |
| K31 plus the AMRFinderPlus program's calls | 2 | better beats worse by the one-sided sign test at p < 0.05 after Holm over 30, for all five | not met | 0 |
| K31 plus the AMRFinderPlus program's calls | 3 | the twin meets neither 1 nor 2 as this reading reads them | met |  |
| K31 plus the AMRFinderPlus program's calls | 4 | no learner's recall interval lies wholly below -0.02, for all five | met | 5 |
| K31 plus the AMRFinderPlus program's calls | 5 | each of the five learners has a study median of its lineage median's sign | not met | 0 |
| K31 plus the AMRFinderPlus program's calls | 6 | the maintainer rules the change | the maintainer's |  |

**Without the constant model.** The count with the constant-model datasets left out, Holm over its own 30 rows. No dataset was scored as the constant model, so `count-no-constant.tsv` is the same file as `count.tsv`, byte for byte, and its table is the count's above.

**Other measures.** Recall, specificity, the very major and the major error rate and log loss on the lineage split, and the primary measure on the study split. Lower is better for log loss and both error rates. Recall stands beside log loss wherever log loss is given.

| Split | Arm | Learner | Measure | Better is | Entered | Better | Worse | Same | Median | Interval |
|---|---|---|---|---|---|---|---|---|---|---|
| lineage | unitigs | `scm` | recall | higher | 232 of 233 | 33 | 39 | 160 | 0.000 | 0.000 to 0.000 |
| lineage | unitigs | `scm` | specificity | higher | 232 of 233 | 35 | 36 | 161 | 0.000 | 0.000 to 0.000 |
| lineage | unitigs | `scm` | very major error rate | lower | 232 of 233 | 33 | 39 | 160 | 0.000 | 0.000 to 0.000 |
| lineage | unitigs | `scm` | major error rate | lower | 232 of 233 | 35 | 36 | 161 | 0.000 | 0.000 to 0.000 |
| lineage | unitigs | `scm` | log loss | lower | 232 of 233 | 75 | 77 | 80 | 0.000 | 0.000 to 0.000 |
| lineage | unitigs | `cart` | recall | higher | 232 of 233 | 38 | 36 | 158 | 0.000 | 0.000 to 0.000 |
| lineage | unitigs | `cart` | specificity | higher | 232 of 233 | 30 | 24 | 178 | 0.000 | 0.000 to 0.000 |
| lineage | unitigs | `cart` | very major error rate | lower | 232 of 233 | 38 | 36 | 158 | 0.000 | 0.000 to 0.000 |
| lineage | unitigs | `cart` | major error rate | lower | 232 of 233 | 30 | 24 | 178 | 0.000 | 0.000 to 0.000 |
| lineage | unitigs | `cart` | log loss | lower | 232 of 233 | 81 | 59 | 92 | 0.000 | -0.003 to 0.000 |
| lineage | unitigs | `rf` | recall | higher | 232 of 233 | 17 | 14 | 201 | 0.000 | 0.000 to 0.000 |
| lineage | unitigs | `rf` | specificity | higher | 232 of 233 | 8 | 16 | 208 | 0.000 | 0.000 to 0.000 |
| lineage | unitigs | `rf` | very major error rate | lower | 232 of 233 | 17 | 14 | 201 | 0.000 | 0.000 to 0.000 |
| lineage | unitigs | `rf` | major error rate | lower | 232 of 233 | 8 | 16 | 208 | 0.000 | 0.000 to 0.000 |
| lineage | unitigs | `rf` | log loss | lower | 232 of 233 | 3 | 7 | 222 | 0.000 | 0.000 to 0.001 |
| lineage | unitigs | `gbt` | recall | higher | 232 of 233 | 45 | 30 | 157 | 0.000 | -0.001 to 0.000 |
| lineage | unitigs | `gbt` | specificity | higher | 232 of 233 | 30 | 26 | 176 | 0.000 | 0.000 to 0.001 |
| lineage | unitigs | `gbt` | very major error rate | lower | 232 of 233 | 45 | 30 | 157 | 0.000 | 0.000 to 0.001 |
| lineage | unitigs | `gbt` | major error rate | lower | 232 of 233 | 30 | 26 | 176 | 0.000 | -0.001 to 0.000 |
| lineage | unitigs | `gbt` | log loss | lower | 232 of 233 | 26 | 19 | 187 | -0.001 | -0.002 to 0.001 |
| lineage | unitigs | `lasso` | recall | higher | 232 of 233 | 25 | 17 | 190 | 0.000 | 0.000 to 0.000 |
| lineage | unitigs | `lasso` | specificity | higher | 232 of 233 | 23 | 20 | 189 | 0.000 | 0.000 to 0.000 |
| lineage | unitigs | `lasso` | very major error rate | lower | 232 of 233 | 25 | 17 | 190 | 0.000 | 0.000 to 0.000 |
| lineage | unitigs | `lasso` | major error rate | lower | 232 of 233 | 23 | 20 | 189 | 0.000 | 0.000 to 0.000 |
| lineage | unitigs | `lasso` | log loss | lower | 232 of 233 | 26 | 28 | 178 | 0.000 | 0.000 to 0.002 |
| lineage | K31 plus k-mer counts | `scm` | recall | higher | 232 of 233 | 34 | 34 | 164 | 0.000 | 0.000 to 0.000 |
| lineage | K31 plus k-mer counts | `scm` | specificity | higher | 232 of 233 | 29 | 31 | 172 | 0.000 | 0.000 to 0.000 |
| lineage | K31 plus k-mer counts | `scm` | very major error rate | lower | 232 of 233 | 34 | 34 | 164 | 0.000 | 0.000 to 0.000 |
| lineage | K31 plus k-mer counts | `scm` | major error rate | lower | 232 of 233 | 29 | 31 | 172 | 0.000 | 0.000 to 0.000 |
| lineage | K31 plus k-mer counts | `scm` | log loss | lower | 232 of 233 | 70 | 68 | 94 | 0.000 | 0.000 to 0.000 |
| lineage | K31 plus k-mer counts | `cart` | recall | higher | 232 of 233 | 35 | 33 | 164 | 0.000 | 0.000 to 0.000 |
| lineage | K31 plus k-mer counts | `cart` | specificity | higher | 232 of 233 | 23 | 34 | 175 | 0.000 | 0.000 to 0.000 |
| lineage | K31 plus k-mer counts | `cart` | very major error rate | lower | 232 of 233 | 35 | 33 | 164 | 0.000 | 0.000 to 0.000 |
| lineage | K31 plus k-mer counts | `cart` | major error rate | lower | 232 of 233 | 23 | 34 | 175 | 0.000 | 0.000 to 0.000 |
| lineage | K31 plus k-mer counts | `cart` | log loss | lower | 232 of 233 | 67 | 61 | 104 | 0.000 | 0.000 to 0.000 |
| lineage | K31 plus k-mer counts | `rf` | recall | higher | 232 of 233 | 34 | 13 | 185 | 0.000 | 0.000 to 0.000 |
| lineage | K31 plus k-mer counts | `rf` | specificity | higher | 232 of 233 | 14 | 18 | 200 | 0.000 | 0.000 to 0.000 |
| lineage | K31 plus k-mer counts | `rf` | very major error rate | lower | 232 of 233 | 34 | 13 | 185 | 0.000 | 0.000 to 0.000 |
| lineage | K31 plus k-mer counts | `rf` | major error rate | lower | 232 of 233 | 14 | 18 | 200 | 0.000 | 0.000 to 0.000 |
| lineage | K31 plus k-mer counts | `rf` | log loss | lower | 232 of 233 | 12 | 3 | 217 | 0.000 | -0.001 to 0.001 |
| lineage | K31 plus k-mer counts | `gbt` | recall | higher | 232 of 233 | 43 | 30 | 159 | 0.000 | 0.000 to 0.000 |
| lineage | K31 plus k-mer counts | `gbt` | specificity | higher | 232 of 233 | 36 | 23 | 173 | 0.000 | 0.000 to 0.001 |
| lineage | K31 plus k-mer counts | `gbt` | very major error rate | lower | 232 of 233 | 43 | 30 | 159 | 0.000 | -0.001 to 0.000 |
| lineage | K31 plus k-mer counts | `gbt` | major error rate | lower | 232 of 233 | 36 | 23 | 173 | 0.000 | -0.001 to 0.000 |
| lineage | K31 plus k-mer counts | `gbt` | log loss | lower | 232 of 233 | 33 | 18 | 181 | -0.001 | -0.002 to 0.002 |
| lineage | K31 plus k-mer counts | `lasso` | recall | higher | 232 of 233 | 31 | 22 | 179 | 0.000 | 0.000 to 0.000 |
| lineage | K31 plus k-mer counts | `lasso` | specificity | higher | 232 of 233 | 23 | 19 | 190 | 0.000 | 0.000 to 0.000 |
| lineage | K31 plus k-mer counts | `lasso` | very major error rate | lower | 232 of 233 | 31 | 22 | 179 | 0.000 | 0.000 to 0.000 |
| lineage | K31 plus k-mer counts | `lasso` | major error rate | lower | 232 of 233 | 23 | 19 | 190 | 0.000 | 0.000 to 0.000 |
| lineage | K31 plus k-mer counts | `lasso` | log loss | lower | 232 of 233 | 33 | 24 | 175 | 0.000 | 0.000 to 0.001 |
| lineage | K31 plus amino-acid k-mers | `scm` | recall | higher | 232 of 233 | 55 | 45 | 132 | 0.000 | 0.000 to 0.005 |
| lineage | K31 plus amino-acid k-mers | `scm` | specificity | higher | 232 of 233 | 58 | 38 | 136 | 0.000 | 0.000 to 0.001 |
| lineage | K31 plus amino-acid k-mers | `scm` | very major error rate | lower | 232 of 233 | 55 | 45 | 132 | 0.000 | -0.005 to 0.000 |
| lineage | K31 plus amino-acid k-mers | `scm` | major error rate | lower | 232 of 233 | 58 | 38 | 136 | 0.000 | -0.001 to 0.000 |
| lineage | K31 plus amino-acid k-mers | `scm` | log loss | lower | 232 of 233 | 101 | 78 | 53 | 0.000 | -0.016 to 0.000 |
| lineage | K31 plus amino-acid k-mers | `cart` | recall | higher | 232 of 233 | 66 | 36 | 130 | 0.000 | 0.000 to 0.001 |
| lineage | K31 plus amino-acid k-mers | `cart` | specificity | higher | 232 of 233 | 52 | 29 | 151 | 0.000 | 0.000 to 0.001 |
| lineage | K31 plus amino-acid k-mers | `cart` | very major error rate | lower | 232 of 233 | 66 | 36 | 130 | 0.000 | -0.001 to 0.000 |
| lineage | K31 plus amino-acid k-mers | `cart` | major error rate | lower | 232 of 233 | 52 | 29 | 151 | 0.000 | -0.001 to 0.000 |
| lineage | K31 plus amino-acid k-mers | `cart` | log loss | lower | 232 of 233 | 94 | 57 | 81 | -0.003 | -0.012 to 0.000 |
| lineage | K31 plus amino-acid k-mers | `rf` | recall | higher | 232 of 233 | 48 | 19 | 165 | 0.001 | 0.000 to 0.004 |
| lineage | K31 plus amino-acid k-mers | `rf` | specificity | higher | 232 of 233 | 28 | 15 | 189 | 0.000 | 0.000 to 0.000 |
| lineage | K31 plus amino-acid k-mers | `rf` | very major error rate | lower | 232 of 233 | 48 | 19 | 165 | -0.001 | -0.004 to 0.000 |
| lineage | K31 plus amino-acid k-mers | `rf` | major error rate | lower | 232 of 233 | 28 | 15 | 189 | 0.000 | -0.001 to 0.000 |
| lineage | K31 plus amino-acid k-mers | `rf` | log loss | lower | 232 of 233 | 44 | 6 | 182 | -0.005 | -0.007 to -0.003 |
| lineage | K31 plus amino-acid k-mers | `gbt` | recall | higher | 232 of 233 | 55 | 36 | 141 | 0.001 | 0.000 to 0.004 |
| lineage | K31 plus amino-acid k-mers | `gbt` | specificity | higher | 232 of 233 | 52 | 27 | 153 | 0.002 | 0.000 to 0.003 |
| lineage | K31 plus amino-acid k-mers | `gbt` | very major error rate | lower | 232 of 233 | 55 | 36 | 141 | -0.001 | -0.004 to 0.000 |
| lineage | K31 plus amino-acid k-mers | `gbt` | major error rate | lower | 232 of 233 | 52 | 27 | 153 | -0.002 | -0.004 to 0.000 |
| lineage | K31 plus amino-acid k-mers | `gbt` | log loss | lower | 232 of 233 | 59 | 32 | 141 | -0.001 | -0.006 to 0.001 |
| lineage | K31 plus amino-acid k-mers | `lasso` | recall | higher | 232 of 233 | 53 | 33 | 146 | 0.000 | 0.000 to 0.006 |
| lineage | K31 plus amino-acid k-mers | `lasso` | specificity | higher | 232 of 233 | 47 | 33 | 152 | 0.000 | 0.000 to 0.002 |
| lineage | K31 plus amino-acid k-mers | `lasso` | very major error rate | lower | 232 of 233 | 53 | 33 | 146 | 0.000 | -0.005 to 0.000 |
| lineage | K31 plus amino-acid k-mers | `lasso` | major error rate | lower | 232 of 233 | 47 | 33 | 152 | 0.000 | -0.002 to 0.000 |
| lineage | K31 plus amino-acid k-mers | `lasso` | log loss | lower | 232 of 233 | 69 | 37 | 126 | -0.003 | -0.009 to 0.000 |
| lineage | the AMRFinderPlus program's calls alone | `scm` | recall | higher | 203 of 204 | 66 | 62 | 75 | 0.000 | -0.001 to 0.007 |
| lineage | the AMRFinderPlus program's calls alone | `scm` | specificity | higher | 203 of 204 | 77 | 47 | 79 | 0.004 | 0.000 to 0.014 |
| lineage | the AMRFinderPlus program's calls alone | `scm` | very major error rate | lower | 203 of 204 | 66 | 62 | 75 | 0.000 | -0.007 to 0.001 |
| lineage | the AMRFinderPlus program's calls alone | `scm` | major error rate | lower | 203 of 204 | 77 | 47 | 79 | -0.004 | -0.014 to 0.000 |
| lineage | the AMRFinderPlus program's calls alone | `scm` | log loss | lower | 203 of 204 | 96 | 80 | 27 | -0.016 | -0.060 to 0.000 |
| lineage | the AMRFinderPlus program's calls alone | `cart` | recall | higher | 203 of 204 | 79 | 63 | 61 | 0.002 | 0.000 to 0.014 |
| lineage | the AMRFinderPlus program's calls alone | `cart` | specificity | higher | 203 of 204 | 85 | 30 | 88 | 0.008 | 0.001 to 0.022 |
| lineage | the AMRFinderPlus program's calls alone | `cart` | very major error rate | lower | 203 of 204 | 79 | 63 | 61 | -0.002 | -0.014 to 0.000 |
| lineage | the AMRFinderPlus program's calls alone | `cart` | major error rate | lower | 203 of 204 | 85 | 30 | 88 | -0.008 | -0.022 to -0.001 |
| lineage | the AMRFinderPlus program's calls alone | `cart` | log loss | lower | 203 of 204 | 117 | 62 | 24 | -0.054 | -0.092 to -0.024 |
| lineage | the AMRFinderPlus program's calls alone | `rf` | recall | higher | 203 of 204 | 104 | 38 | 61 | 0.023 | 0.011 to 0.043 |
| lineage | the AMRFinderPlus program's calls alone | `rf` | specificity | higher | 203 of 204 | 78 | 48 | 77 | 0.009 | 0.000 to 0.017 |
| lineage | the AMRFinderPlus program's calls alone | `rf` | very major error rate | lower | 203 of 204 | 104 | 38 | 61 | -0.023 | -0.048 to -0.011 |
| lineage | the AMRFinderPlus program's calls alone | `rf` | major error rate | lower | 203 of 204 | 78 | 48 | 77 | -0.009 | -0.018 to 0.000 |
| lineage | the AMRFinderPlus program's calls alone | `rf` | log loss | lower | 203 of 204 | 119 | 35 | 49 | -0.039 | -0.063 to -0.017 |
| lineage | the AMRFinderPlus program's calls alone | `gbt` | recall | higher | 203 of 204 | 88 | 49 | 66 | 0.010 | 0.002 to 0.026 |
| lineage | the AMRFinderPlus program's calls alone | `gbt` | specificity | higher | 203 of 204 | 88 | 42 | 73 | 0.010 | 0.004 to 0.021 |
| lineage | the AMRFinderPlus program's calls alone | `gbt` | very major error rate | lower | 203 of 204 | 88 | 49 | 66 | -0.010 | -0.026 to -0.002 |
| lineage | the AMRFinderPlus program's calls alone | `gbt` | major error rate | lower | 203 of 204 | 88 | 42 | 73 | -0.010 | -0.021 to -0.005 |
| lineage | the AMRFinderPlus program's calls alone | `gbt` | log loss | lower | 203 of 204 | 95 | 54 | 54 | -0.017 | -0.031 to -0.002 |
| lineage | the AMRFinderPlus program's calls alone | `lasso` | recall | higher | 203 of 204 | 82 | 50 | 71 | 0.012 | 0.004 to 0.018 |
| lineage | the AMRFinderPlus program's calls alone | `lasso` | specificity | higher | 203 of 204 | 87 | 35 | 81 | 0.014 | 0.004 to 0.024 |
| lineage | the AMRFinderPlus program's calls alone | `lasso` | very major error rate | lower | 203 of 204 | 82 | 50 | 71 | -0.012 | -0.018 to -0.004 |
| lineage | the AMRFinderPlus program's calls alone | `lasso` | major error rate | lower | 203 of 204 | 87 | 35 | 81 | -0.014 | -0.024 to -0.004 |
| lineage | the AMRFinderPlus program's calls alone | `lasso` | log loss | lower | 203 of 204 | 109 | 48 | 46 | -0.024 | -0.041 to -0.012 |
| lineage | K31 plus the AMRFinderPlus program's calls | `scm` | recall | higher | 203 of 204 | 19 | 10 | 174 | 0.000 | 0.000 to 0.000 |
| lineage | K31 plus the AMRFinderPlus program's calls | `scm` | specificity | higher | 203 of 204 | 6 | 10 | 187 | 0.000 | 0.000 to 0.000 |
| lineage | K31 plus the AMRFinderPlus program's calls | `scm` | very major error rate | lower | 203 of 204 | 19 | 10 | 174 | 0.000 | 0.000 to 0.000 |
| lineage | K31 plus the AMRFinderPlus program's calls | `scm` | major error rate | lower | 203 of 204 | 6 | 10 | 187 | 0.000 | 0.000 to 0.000 |
| lineage | K31 plus the AMRFinderPlus program's calls | `scm` | log loss | lower | 203 of 204 | 41 | 27 | 135 | 0.000 | 0.000 to 0.000 |
| lineage | K31 plus the AMRFinderPlus program's calls | `cart` | recall | higher | 203 of 204 | 9 | 7 | 187 | 0.000 | 0.000 to 0.000 |
| lineage | K31 plus the AMRFinderPlus program's calls | `cart` | specificity | higher | 203 of 204 | 2 | 8 | 193 | 0.000 | 0.000 to 0.000 |
| lineage | K31 plus the AMRFinderPlus program's calls | `cart` | very major error rate | lower | 203 of 204 | 9 | 7 | 187 | 0.000 | 0.000 to 0.000 |
| lineage | K31 plus the AMRFinderPlus program's calls | `cart` | major error rate | lower | 203 of 204 | 2 | 8 | 193 | 0.000 | 0.000 to 0.000 |
| lineage | K31 plus the AMRFinderPlus program's calls | `cart` | log loss | lower | 203 of 204 | 33 | 33 | 137 | 0.000 | 0.000 to 0.000 |
| lineage | K31 plus the AMRFinderPlus program's calls | `rf` | recall | higher | 203 of 204 | 10 | 6 | 187 | 0.000 | 0.000 to 0.000 |
| lineage | K31 plus the AMRFinderPlus program's calls | `rf` | specificity | higher | 203 of 204 | 3 | 5 | 195 | 0.000 | 0.000 to 0.000 |
| lineage | K31 plus the AMRFinderPlus program's calls | `rf` | very major error rate | lower | 203 of 204 | 10 | 6 | 187 | 0.000 | 0.000 to 0.000 |
| lineage | K31 plus the AMRFinderPlus program's calls | `rf` | major error rate | lower | 203 of 204 | 3 | 5 | 195 | 0.000 | 0.000 to 0.000 |
| lineage | K31 plus the AMRFinderPlus program's calls | `rf` | log loss | lower | 203 of 204 | 5 | 0 | 198 | 0.000 | -0.001 to 0.000 |
| lineage | K31 plus the AMRFinderPlus program's calls | `gbt` | recall | higher | 203 of 204 | 19 | 18 | 166 | 0.000 | 0.000 to 0.000 |
| lineage | K31 plus the AMRFinderPlus program's calls | `gbt` | specificity | higher | 203 of 204 | 12 | 14 | 177 | 0.000 | 0.000 to 0.000 |
| lineage | K31 plus the AMRFinderPlus program's calls | `gbt` | very major error rate | lower | 203 of 204 | 19 | 18 | 166 | 0.000 | 0.000 to 0.000 |
| lineage | K31 plus the AMRFinderPlus program's calls | `gbt` | major error rate | lower | 203 of 204 | 12 | 14 | 177 | 0.000 | 0.000 to 0.000 |
| lineage | K31 plus the AMRFinderPlus program's calls | `gbt` | log loss | lower | 203 of 204 | 8 | 12 | 183 | 0.000 | 0.000 to 0.000 |
| lineage | K31 plus the AMRFinderPlus program's calls | `lasso` | recall | higher | 203 of 204 | 3 | 3 | 197 | 0.000 | 0.000 to 0.000 |
| lineage | K31 plus the AMRFinderPlus program's calls | `lasso` | specificity | higher | 203 of 204 | 2 | 0 | 201 | 0.000 | 0.000 to 0.000 |
| lineage | K31 plus the AMRFinderPlus program's calls | `lasso` | very major error rate | lower | 203 of 204 | 3 | 3 | 197 | 0.000 | 0.000 to 0.000 |
| lineage | K31 plus the AMRFinderPlus program's calls | `lasso` | major error rate | lower | 203 of 204 | 2 | 0 | 201 | 0.000 | 0.000 to 0.000 |
| lineage | K31 plus the AMRFinderPlus program's calls | `lasso` | log loss | lower | 203 of 204 | 5 | 3 | 195 | 0.000 | 0.000 to 0.000 |
| lineage | the twin, the default at seeds 45, 46 and 47 | `scm` | recall | higher | 232 of 233 | 44 | 43 | 145 | 0.000 | 0.000 to 0.000 |
| lineage | the twin, the default at seeds 45, 46 and 47 | `scm` | specificity | higher | 232 of 233 | 45 | 36 | 151 | 0.000 | 0.000 to 0.001 |
| lineage | the twin, the default at seeds 45, 46 and 47 | `scm` | very major error rate | lower | 232 of 233 | 44 | 43 | 145 | 0.000 | 0.000 to 0.000 |
| lineage | the twin, the default at seeds 45, 46 and 47 | `scm` | major error rate | lower | 232 of 233 | 45 | 36 | 151 | 0.000 | -0.001 to 0.000 |
| lineage | the twin, the default at seeds 45, 46 and 47 | `scm` | log loss | lower | 232 of 233 | 85 | 79 | 68 | 0.000 | -0.003 to 0.002 |
| lineage | the twin, the default at seeds 45, 46 and 47 | `cart` | recall | higher | 232 of 233 | 40 | 44 | 148 | 0.000 | 0.000 to 0.000 |
| lineage | the twin, the default at seeds 45, 46 and 47 | `cart` | specificity | higher | 232 of 233 | 37 | 38 | 157 | 0.000 | 0.000 to 0.000 |
| lineage | the twin, the default at seeds 45, 46 and 47 | `cart` | very major error rate | lower | 232 of 233 | 40 | 44 | 148 | 0.000 | 0.000 to 0.000 |
| lineage | the twin, the default at seeds 45, 46 and 47 | `cart` | major error rate | lower | 232 of 233 | 37 | 38 | 157 | 0.000 | 0.000 to 0.000 |
| lineage | the twin, the default at seeds 45, 46 and 47 | `cart` | log loss | lower | 232 of 233 | 81 | 56 | 95 | 0.000 | -0.003 to 0.000 |
| lineage | the twin, the default at seeds 45, 46 and 47 | `rf` | recall | higher | 232 of 233 | 25 | 39 | 168 | 0.000 | -0.001 to 0.000 |
| lineage | the twin, the default at seeds 45, 46 and 47 | `rf` | specificity | higher | 232 of 233 | 21 | 28 | 183 | 0.000 | -0.001 to 0.000 |
| lineage | the twin, the default at seeds 45, 46 and 47 | `rf` | very major error rate | lower | 232 of 233 | 25 | 39 | 168 | 0.000 | 0.000 to 0.001 |
| lineage | the twin, the default at seeds 45, 46 and 47 | `rf` | major error rate | lower | 232 of 233 | 21 | 28 | 183 | 0.000 | 0.000 to 0.001 |
| lineage | the twin, the default at seeds 45, 46 and 47 | `rf` | log loss | lower | 232 of 233 | 10 | 18 | 204 | 0.001 | 0.000 to 0.002 |
| lineage | the twin, the default at seeds 45, 46 and 47 | `gbt` | recall | higher | 232 of 233 | 53 | 47 | 132 | 0.000 | 0.000 to 0.000 |
| lineage | the twin, the default at seeds 45, 46 and 47 | `gbt` | specificity | higher | 232 of 233 | 37 | 38 | 157 | 0.000 | 0.000 to 0.001 |
| lineage | the twin, the default at seeds 45, 46 and 47 | `gbt` | very major error rate | lower | 232 of 233 | 53 | 47 | 132 | 0.000 | 0.000 to 0.000 |
| lineage | the twin, the default at seeds 45, 46 and 47 | `gbt` | major error rate | lower | 232 of 233 | 37 | 38 | 157 | 0.000 | -0.001 to 0.000 |
| lineage | the twin, the default at seeds 45, 46 and 47 | `gbt` | log loss | lower | 232 of 233 | 36 | 48 | 148 | 0.001 | -0.001 to 0.003 |
| lineage | the twin, the default at seeds 45, 46 and 47 | `lasso` | recall | higher | 232 of 233 | 41 | 32 | 159 | 0.000 | 0.000 to 0.001 |
| lineage | the twin, the default at seeds 45, 46 and 47 | `lasso` | specificity | higher | 232 of 233 | 33 | 31 | 168 | 0.000 | -0.003 to 0.000 |
| lineage | the twin, the default at seeds 45, 46 and 47 | `lasso` | very major error rate | lower | 232 of 233 | 41 | 32 | 159 | 0.000 | -0.001 to 0.000 |
| lineage | the twin, the default at seeds 45, 46 and 47 | `lasso` | major error rate | lower | 232 of 233 | 33 | 31 | 168 | 0.000 | 0.000 to 0.003 |
| lineage | the twin, the default at seeds 45, 46 and 47 | `lasso` | log loss | lower | 232 of 233 | 52 | 54 | 126 | 0.000 | -0.002 to 0.006 |
| study | unitigs | `scm` | balanced accuracy | higher | 55 of 55 | 8 | 6 | 41 | 0.000 | 0.000 to 0.001 |
| study | unitigs | `cart` | balanced accuracy | higher | 55 of 55 | 7 | 6 | 42 | 0.000 | 0.000 to 0.002 |
| study | unitigs | `rf` | balanced accuracy | higher | 55 of 55 | 2 | 2 | 51 | -0.001 | -0.001 to 0.000 |
| study | unitigs | `gbt` | balanced accuracy | higher | 55 of 55 | 4 | 6 | 45 | 0.000 | -0.004 to 0.001 |
| study | unitigs | `lasso` | balanced accuracy | higher | 55 of 55 | 2 | 1 | 52 | 0.000 | -0.001 to 0.002 |
| study | K31 plus k-mer counts | `scm` | balanced accuracy | higher | 55 of 55 | 9 | 7 | 39 | 0.000 | -0.001 to 0.004 |
| study | K31 plus k-mer counts | `cart` | balanced accuracy | higher | 55 of 55 | 11 | 11 | 33 | 0.000 | 0.000 to 0.003 |
| study | K31 plus k-mer counts | `rf` | balanced accuracy | higher | 55 of 55 | 5 | 6 | 44 | 0.000 | -0.003 to 0.004 |
| study | K31 plus k-mer counts | `gbt` | balanced accuracy | higher | 55 of 55 | 8 | 10 | 37 | 0.001 | 0.000 to 0.009 |
| study | K31 plus k-mer counts | `lasso` | balanced accuracy | higher | 55 of 55 | 8 | 9 | 38 | 0.000 | -0.001 to 0.001 |
| study | K31 plus amino-acid k-mers | `scm` | balanced accuracy | higher | 55 of 55 | 17 | 11 | 27 | 0.003 | 0.000 to 0.019 |
| study | K31 plus amino-acid k-mers | `cart` | balanced accuracy | higher | 55 of 55 | 17 | 11 | 27 | 0.005 | 0.002 to 0.012 |
| study | K31 plus amino-acid k-mers | `rf` | balanced accuracy | higher | 55 of 55 | 6 | 3 | 46 | 0.002 | 0.000 to 0.008 |
| study | K31 plus amino-acid k-mers | `gbt` | balanced accuracy | higher | 55 of 55 | 13 | 7 | 35 | 0.004 | 0.002 to 0.012 |
| study | K31 plus amino-acid k-mers | `lasso` | balanced accuracy | higher | 55 of 55 | 11 | 8 | 36 | 0.002 | 0.000 to 0.017 |
| study | the AMRFinderPlus program's calls alone | `scm` | balanced accuracy | higher | 54 of 54 | 23 | 10 | 21 | 0.015 | -0.001 to 0.043 |
| study | the AMRFinderPlus program's calls alone | `cart` | balanced accuracy | higher | 54 of 54 | 21 | 11 | 22 | 0.007 | 0.002 to 0.027 |
| study | the AMRFinderPlus program's calls alone | `rf` | balanced accuracy | higher | 54 of 54 | 36 | 7 | 11 | 0.040 | 0.027 to 0.061 |
| study | the AMRFinderPlus program's calls alone | `gbt` | balanced accuracy | higher | 54 of 54 | 25 | 8 | 21 | 0.016 | 0.007 to 0.040 |
| study | the AMRFinderPlus program's calls alone | `lasso` | balanced accuracy | higher | 54 of 54 | 32 | 5 | 17 | 0.030 | 0.009 to 0.051 |
| study | K31 plus the AMRFinderPlus program's calls | `scm` | balanced accuracy | higher | 54 of 54 | 3 | 2 | 49 | 0.000 | 0.000 to 0.000 |
| study | K31 plus the AMRFinderPlus program's calls | `cart` | balanced accuracy | higher | 54 of 54 | 1 | 1 | 52 | 0.000 | 0.000 to 0.000 |
| study | K31 plus the AMRFinderPlus program's calls | `rf` | balanced accuracy | higher | 54 of 54 | 3 | 2 | 49 | 0.000 | -0.002 to 0.000 |
| study | K31 plus the AMRFinderPlus program's calls | `gbt` | balanced accuracy | higher | 54 of 54 | 0 | 3 | 51 | 0.000 | -0.002 to 0.001 |
| study | K31 plus the AMRFinderPlus program's calls | `lasso` | balanced accuracy | higher | 54 of 54 | 1 | 0 | 53 | 0.000 | 0.000 to 0.000 |

## Seen after

These points are seen in the tables after the reading. They change no rule and no reading above, and decide nothing.

**Joined to K31.** The program's calls joined to K31 score as K31 alone on most datasets: K31 plus the AMRFinderPlus program's calls is the same as the default on 183 to 201 of 203 lineage datasets per learner, where the program's calls alone are the same on 54 to 80 of 203. These tables do not say why.

**Seed alone.** The twin moves 70 of 1,160 lineage dataset and learner pairs by more than 0.05 in balanced accuracy, from -0.141 (*Klebsiella pneumoniae* with ticarcillin/clavulanate, `gbt`) to 0.169 (*Mycobacterium tuberculosis* with nicotinamide, `gbt`). The AMRFinderPlus program's calls alone move 407 of 1,015 by more than 0.05, from -0.365 (*Klebsiella pneumoniae* with ceftazidime/avibactam, `cart`, 0.895 to 0.530) to 0.363 (*Corynebacterium diphtheriae* with sulfamethoxazole/trimethoprim, `rf`, 0.500 to 0.863). K31 plus amino-acid k-mers move 142 of 1,160, from -0.362 (*Klebsiella pneumoniae* with ceftazidime/avibactam, `scm`, 0.930 to 0.568) to 0.349 (*Salmonella enterica* with sulfamethoxazole, `gbt`, 0.603 to 0.952).

**Species.** Summed over the five learners as dataset and learner pairs, the program's calls alone are better on more pairs than worse in 20 of A's 21 species. In *Klebsiella aerogenes* they are worse on 12 of 15 pairs and better on 1.

**Before Holm.** Unitigs with `gbt` give a sign p of 0.005 for 34 better against 15 worse of 232, and 0.099 after Holm.

## Limits

**One database.** Every split is a retrospective split inside one public database. A split by lineage or by study is not external validation.

**Labels.** The 35 datasets of M wait for the maintainer's label rule. The count over them is taken when that rule lets them run, with these rules, this scorer and these arguments unchanged. After the floor `m_page.py left` writes the left-out page, `docs/input-arms-m-left-out.md`, and after the count `m_page.py page` writes the count page, `docs/input-arms-m-count.md`; this page holds no table of that run.

**Study split.** One study run per dataset at the release's seeds, and no study twin.

**Contig order.** The figures of K31 plus amino-acid k-mers are for the contig order BV-BRC served on 23 Nov 2026, since Prodigal's gene calls change with it. On a second fetch on 27 Nov, 98 of the 2,000 *Neisseria gonorrhoeae* files gave their contigs in another order. That arm's difference from K31 moved on 13 of the species' 30 lineage rows, on 2 of them by more than 0.02 (`docs/input-arms.md`, Stop of 27 Nov 2026). The count on this page was not taken again.

**Constant model.** A study run refused for no candidate is scored as the constant model (ruling 16). It enters the count and is marked. The count without it decides nothing.

**The program.** A program-call arm would enter only with the program's version pinned. A program column is named only as the AMRFinderPlus program's table names it.

**The catalogue.** The AMRFinderPlus program reads a catalogue curated from published resistance work. A lineage or a study held out of a fit is still inside that catalogue. So the program's calls bring knowledge from outside the training folds, and 31-mer presence brings none. These tables do not measure how much of the difference comes from that.

**Binaries.** One binary: every species' `binary.txt` gives the commit `8bfd4ca` and the code sha256 `d3de72fcc10f6ef9503eb43f8cb5fd1449bca137d5f0a1594caf6fda0c614363`, 23 of 23 species, with no restart.

**Unread files.** A file under `score/` that neither `SHA256SUMS` lists nor the run's rule covers stays unchecked and unread. `score/score.log` in each `score.tar` was never opened, hashed or copied out.

**Markers.** No marker is found, named or ranked. The page holds aggregates over datasets, species, splits, arms and learners only. No file holds a genome, a k-mer, a unitig, a sequence, a position or a column name of a matrix.

**The release.** The release keeps its input whatever the count says, until the maintainer rules.

## Method

**Checks of the run.** Before any score is read, `scripts/c3-inputs/run_checks.py` reads each species' copy and prints counts only. Each check can fail, and the count is taken only when all six pass.

| Check | Holds |
|---|---|
| species | the copies are exactly the species of `scripts/c3-inputs/lists.tsv`'s first run |
| `SHA256SUMS` | every line of each species' `SHA256SUMS` holds, the files under `score/` read in `score.tar` when it is not unpacked, `DONE` names its sha256, and the files under `results/` and `check/` are its lines |
| `DONE` and `binary.txt` | each `DONE` names its species, binary, code, result count and `SHA256SUMS` as `binary.txt` does; no `FAILED` but one named, set aside and older than its species' `DONE` |
| `expect.py` | each `expect.tsv` is `expect.py`'s from the first run, `lists.txt` names the design's seeds, twin seeds and arms, and the result files are exactly its rows |
| refusals | each refusal record is an arm's study run refused for no candidate, with five constant-model rows in `draws.tsv`, one per learner, and no others |
| genome ids | no genome id of the release's species lists in any file of the copies that is read or published |

All six passed. `checks.txt` in `reference/input-arms/` holds the printout.

| Check | Result |
|---|---|
| species | 23 copies, the 23 species of `lists.tsv`'s first run |
| `SHA256SUMS` | holds in 23 of 23 species, 10,953 of 10,953 lines; `results/` and `check/` are the listed files in 23 of 23 |
| `DONE` and `binary.txt` | 23 of 23 `DONE` lines name their species, `SHA256SUMS`, result count, binary and code as `binary.txt` does; 23 species on commit `8bfd4ca` and code `d3de72fcc10f`, no restart; one `FAILED` file, *Vibrio cholerae*'s at its build at 07:29:21Z, older than its `DONE` at 07:49:38Z, named and set aside |
| `expect.py` | 23 of 23 `expect.tsv` files are `expect.py`'s from the first run; 5,047 of 5,047 results found and none more; 233 datasets, 55 in S, 23 species |
| refusals | 0 refusal records and 0 constant-model rows |
| genome ids | 0 in 13,362 files, against the 23,757 genome ids of the release's 23 species lists |

**Reading.** The answers read are the run's own: each species' `score/` tables, copied as its run wrote them. The count runs once over them:

```
python3 -I scripts/c3-inputs/run_checks.py <back dir> scripts/c3-inputs/lists.tsv <release dir> --set-aside vibrio_cholerae > checks.txt
# gather: per species, score/checks.tsv, score/draws.tsv and score/datasets.tsv out of score.tar into read/score/<species>/,
#         binary.txt, DONE and results/builds/<kind>.build.json into read/run/<species>/
python3 -I scripts/c3-inputs/score.py --count read/count read/score/<species>/datasets.tsv   # 23 files, once
python3 -I scripts/c3-inputs/page_tables.py write read scripts/c3-inputs/lists.tsv checks.txt reference/input-arms
python3 -I scripts/c3-inputs/page_tables.py print reference/input-arms
```

`reading.txt` lists every command with the sha256 of each file it wrote.

The gather step read each species' `checks.tsv`, `draws.tsv` and `datasets.tsv` out of its `score.tar` in place, 69 tables, each matched to its line in the species' `SHA256SUMS`, and copied `binary.txt`, `DONE` and the three build records of each species as they are, 115 files. It never opened `score/score.log`. The count printed `count: 9505 dataset rows from 23 files`, and the table script `gates: 13 of 13 pass; wrote 16 files`. Before the tables were copied out, a count of the release's 23,757 genome ids over them, the checks printout and the steps' logs found 0 in 34 files.

**Tables.** `scripts/c3-inputs/page_tables.py write` checks the reading against thirteen gates and stops at the first that fails, naming it. All thirteen passed. Gate 1 refuses to write unless the run's checks printout says that every check passed. The count files are the scorer's own tables: the program counts them again from the species' `datasets.tsv` files (gate 13), copies them as they are, and sums them with every other file in `SHA256SUMS`.

| Gate | Passes when |
|---|---|
| 1. reading | the run's checks printout ends `CHECKS n of n`, n being `run_checks.py`'s own count of checks, with a PASS line for each and no line starting FAIL; the reading holds `reading.txt`, `count/` with `count.tsv`, `count-no-constant.tsv` and `conditions.tsv`, `score/` and `run/`, and nothing else; an arm `count.tsv` marks withdrawn is one of the five |
| 2. species | `score/` and `run/` each hold exactly the species of the first run, each once, with its files, and each species' rows name it alone |
| 3. twin | the twin's rows are lineage draws at seeds 45, 46 and 47 only, no arm's row is at a twin seed, and `conditions.tsv` reads the five arms and never the twin |
| 4. study | a study row is a study run at seed 42 on the study split, and every lineage row is a draw at 42, 43 or 44 or the twin's |
| 5. M | no dataset of M, or outside the first run, stands in any table; species and antibiotic are `lists.tsv`'s |
| 6. own datasets | each arm's rows are on its own set's datasets only: the default, the twin and the three matrix arms on L's first run, the two program arms on A's, the study split on those in S, and a withdrawn arm on none |
| 7. draws | `datasets.tsv` is the mean of `draws.tsv`'s draws, and a dataset with a draw not scored is not scored, never read as zero |
| 8. constant | a constant-model draw is an arm's study run refused for no candidate, marked in `draws.tsv` and counted in `datasets.tsv`, and its dataset is kept |
| 9. checks | each species' `checks.tsv` holds checks and every one passed |
| 10. builds | each species' builds are `grm-c3-build-1` records of its species, one per kind of the matrix arms run, holding counts only |
| 11. binary | each species' `binary.txt` and `DONE` name the same binary and code; one code sha256 in every species, or the page names each binary with its species |
| 12. conditions | ruling 15's rows decide and the stricter reading's decide nothing, and conditions 1 and 2 as read here off `count.tsv` equal `conditions.tsv`'s |
| 13. count | `count.tsv`, `count-no-constant.tsv` and `conditions.tsv` are what `score.py --count` writes from the `datasets.tsv` files, with the arms `count.tsv` marks withdrawn |

**Tests.** `tests/c3_inputs/test_page_tables.py` holds 53 SYNTHETIC tests of the table script and `tests/c3_inputs/test_run_checks.py` 24 of the run's checks, on three made-up species, and `tests/c3_inputs/test_page_faults.py` 4 of the fault runner. `scripts/c3-inputs/page_faults.py` plants each fault below in a copy of the tree and runs the test meant to catch it, which must fail; `tests/c3_inputs/test_page_faults.py` checks that each fault's text occurs once in its file and that its test exists. All 29 were caught.

| Fault | Caught by |
|---|---|
| 1. A folder outside the reading passed (1), or a reading written though the run's checks did not all pass (1b) | `Gates.test_a_folder_outside_the_reading_is_refused`; `Gates.test_a_checks_printout_with_a_failed_check_is_refused` |
| 2. A species read twice | `Gates.test_a_species_twice_is_refused` |
| 3. The twin read as an arm in `conditions.tsv` | `Gates.test_conditions_that_read_the_twin_as_an_arm_are_refused` |
| 4. A study row counted on the lineage split, by the gate (4a) or in the species table (4b), or conditions 1 and 2 read on the study rows (4c) | `Gates.test_a_study_run_on_the_lineage_split_is_refused`; `Tables.test_a_study_row_stays_on_the_study_split`; `Tables.test_the_count_rows_carry_conditions_one_and_two_on_the_lineage_primary_rows` |
| 5. A dataset of M counted | `Gates.test_a_dataset_of_m_is_refused` |
| 6. An arm's rows passed on datasets not its own (6a), or a program arm counted on L's datasets in the arms table (6b) | `Gates.test_a_program_arm_outside_a_is_refused`; `Tables.test_each_arm_stands_on_its_own_datasets` |
| 7. A value passed for a dataset with a draw not scored (7a), or a missing draw read as zero (7b) | `Gates.test_a_value_for_a_dataset_with_a_draw_not_scored_is_refused`; `Tables.test_a_draw_not_scored_leaves_its_dataset_unscored_not_zero` |
| 8. A constant-model draw passed without its reason (8a) or off an arm's study run (8b), or a study run refused for no candidate left unmarked (8c) | `Gates.test_a_constant_model_draw_without_its_reason_is_refused`; `Gates.test_a_constant_model_draw_off_an_arms_study_run_is_refused`; `Tables.test_a_constant_model_row_is_kept_and_marked` |
| 9. A failed check of the scorer passed | `Gates.test_a_failed_check_is_refused` |
| 10. A build record holding more than its counts passed | `Gates.test_a_build_record_beyond_its_counts_is_refused` |
| 11. A `DONE` line naming another code passed (11a), a page that leaves a species' binary out passed (11b), or two binaries written as one (11c) | `Gates.test_a_done_line_naming_another_code_is_refused`; `Gates.test_a_page_that_leaves_a_species_binary_out_is_refused`; `Tables.test_two_binaries_are_each_named_with_their_species` |
| 12. The stricter reading written as deciding (12a), or `conditions.tsv` passed without reading it again off `count.tsv` (12b) | `Gates.test_a_page_that_writes_the_stricter_reading_as_deciding_is_refused`; `Gates.test_conditions_that_differ_from_the_count_are_refused` |
| 13. A count passed that is not the scorer's | `Gates.test_a_count_that_is_not_the_scorers_is_refused` |
| R1. A file changed since its sum passed | `RunChecks.test_a_file_changed_since_its_sum_fails` |
| R2. A `FAILED` file not set aside passed | `RunChecks.test_a_failed_file_not_set_aside_fails_and_one_set_aside_is_named` |
| R3. A result file more than `expect.py`'s rows passed | `RunChecks.test_a_result_file_more_than_its_rows_fails` |
| R4. A refused default passed | `RunChecks.test_a_refused_default_fails` |
| R5. A refusal for another reason passed | `RunChecks.test_a_refusal_for_another_reason_fails` |
| R6. `score.tar` left out of the genome-id count | `RunChecks.test_a_genome_id_in_score_tar_is_counted` |

**Files.** The reading's tables are in `reference/input-arms/`. `scripts/c3-inputs/page_tables.py print` prints this page's tables from those files, rounded to 3 decimals.

## Datasets

Per species, named in full: the datasets of the first run, those in A and in S, those of M waiting, and better, worse and same on the primary measure by arm on the lineage split. A species' median decides nothing.

| Species | Datasets | In A | In S | M waiting | unitigs: better / worse / same | K31 plus k-mer counts: better / worse / same | K31 plus amino-acid k-mers: better / worse / same | the AMRFinderPlus program's calls alone: better / worse / same | K31 plus the AMRFinderPlus program's calls: better / worse / same | the twin, the default at seeds 45, 46 and 47: better / worse / same |
|---|---|---|---|---|---|---|---|---|---|---|
| *Acinetobacter baumannii* | 7 | 7 | 2 | 11 | 4 / 3 / 28 | 7 / 2 / 26 | 10 / 3 / 22 | 32 / 1 / 2 | 1 / 0 / 34 | 12 / 6 / 17 |
| *Campylobacter coli* | 8 | 7 | 0 | 0 | 1 / 2 / 37 | 2 / 1 / 37 | 2 / 2 / 36 | 11 / 0 / 24 | 1 / 0 / 34 | 2 / 4 / 34 |
| *Campylobacter jejuni* | 4 | 4 | 0 | 0 | 1 / 1 / 18 | 1 / 1 / 18 | 1 / 1 / 18 | 6 / 1 / 13 | 1 / 0 / 19 | 3 / 1 / 16 |
| *Clostridioides difficile* | 4 | 4 | 0 | 0 | 0 / 0 / 20 | 0 / 0 / 20 | 0 / 0 / 20 | 11 / 0 / 9 | 0 / 0 / 20 | 0 / 2 / 18 |
| *Corynebacterium diphtheriae* | 3 | 3 | 0 | 0 | 2 / 4 / 9 | 1 / 1 / 13 | 7 / 0 / 8 | 10 / 1 / 4 | 0 / 0 / 15 | 3 / 3 / 9 |
| *Enterobacter cloacae* | 9 | 9 | 8 | 0 | 12 / 8 / 25 | 13 / 5 / 27 | 24 / 8 / 13 | 32 / 12 / 1 | 6 / 1 / 38 | 11 / 5 / 29 |
| *Enterococcus faecium* | 11 | 11 | 0 | 2 | 3 / 0 / 52 | 6 / 5 / 44 | 4 / 1 / 50 | 18 / 10 / 27 | 1 / 1 / 53 | 11 / 2 / 42 |
| *Escherichia coli* | 29 | 29 | 15 | 7 | 6 / 12 / 127 | 15 / 11 / 119 | 25 / 7 / 113 | 69 / 27 / 49 | 0 / 7 / 138 | 9 / 10 / 126 |
| *Klebsiella aerogenes* | 3 | 3 | 0 | 0 | 3 / 4 / 8 | 2 / 4 / 9 | 7 / 2 / 6 | 1 / 12 / 2 | 0 / 1 / 14 | 4 / 5 / 6 |
| *Klebsiella oxytoca* | 4 | 4 | 0 | 0 | 7 / 5 / 8 | 4 / 5 / 11 | 7 / 7 / 6 | 17 / 2 / 1 | 1 / 1 / 18 | 8 / 7 / 5 |
| *Klebsiella pneumoniae* | 32 | 31 | 18 | 6 | 21 / 12 / 122 | 23 / 11 / 121 | 35 / 19 / 101 | 78 / 36 / 36 | 3 / 4 / 143 | 18 / 26 / 111 |
| *Klebsiella quasipneumoniae* | 3 | 3 | 0 | 0 | 0 / 0 / 15 | 2 / 0 / 13 | 5 / 0 / 10 | 7 / 1 / 7 | 0 / 0 / 15 | 6 / 0 / 9 |
| *Mycobacterium tuberculosis* | 20 | 0 | 1 | 2 | 2 / 7 / 91 | 1 / 8 / 91 | 3 / 4 / 93 | 0 / 0 / 0 | 0 / 0 / 0 | 10 / 11 / 79 |
| *Neisseria gonorrhoeae* | 6 | 6 | 0 | 1 | 4 / 3 / 23 | 4 / 3 / 23 | 1 / 1 / 28 | 12 / 10 / 8 | 1 / 1 / 28 | 4 / 4 / 22 |
| *Pseudomonas aeruginosa* | 12 | 12 | 6 | 1 | 11 / 3 / 46 | 11 / 6 / 43 | 23 / 5 / 32 | 32 / 1 / 27 | 0 / 3 / 57 | 8 / 8 / 44 |
| *Salmonella enterica* | 24 | 23 | 1 | 0 | 8 / 11 / 101 | 12 / 9 / 99 | 35 / 10 / 75 | 54 / 18 / 43 | 2 / 1 / 112 | 14 / 18 / 88 |
| *Shigella flexneri* | 4 | 4 | 0 | 0 | 2 / 1 / 17 | 1 / 1 / 18 | 6 / 1 / 13 | 9 / 1 / 10 | 1 / 1 / 18 | 2 / 0 / 18 |
| *Shigella sonnei* | 20 | 18 | 0 | 0 | 2 / 2 / 96 | 1 / 2 / 97 | 10 / 2 / 88 | 30 / 20 / 40 | 0 / 1 / 89 | 3 / 4 / 93 |
| *Staphylococcus aureus* | 13 | 12 | 4 | 1 | 7 / 5 / 53 | 6 / 3 / 56 | 20 / 5 / 40 | 23 / 11 / 26 | 3 / 2 / 55 | 3 / 6 / 56 |
| *Staphylococcus haemolyticus* | 3 | 0 | 0 | 0 | 1 / 3 / 11 | 5 / 2 / 8 | 2 / 3 / 10 | 0 / 0 / 0 | 0 / 0 / 0 | 4 / 3 / 8 |
| *Streptococcus agalactiae* | 3 | 3 | 0 | 0 | 0 / 0 / 15 | 0 / 0 / 15 | 2 / 0 / 13 | 10 / 3 / 2 | 0 / 0 / 15 | 0 / 1 / 14 |
| *Streptococcus pneumoniae* | 8 | 8 | 0 | 4 | 0 / 0 / 40 | 0 / 2 / 38 | 4 / 6 / 30 | 15 / 11 / 14 | 0 / 0 / 40 | 5 / 3 / 32 |
| *Vibrio cholerae* | 3 | 3 | 0 | 0 | 2 / 0 / 13 | 1 / 0 / 14 | 5 / 1 / 9 | 6 / 3 / 6 | 5 / 2 / 8 | 2 / 2 / 11 |
