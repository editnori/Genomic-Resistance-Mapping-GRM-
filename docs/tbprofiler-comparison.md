# TBProfiler comparison

This page counts the calls of TBProfiler 6.6.3, a public catalogue-based caller, against the phenotypes recorded in GRM's 22 M. tuberculosis reports of evaluation revision 14, beside the five learners' calls on the same genomes. This is a retrospective comparison on assemblies of public BV-BRC genomes. It is not a prospective test and not external validation.

This is a diagnostic comparison of a published catalogue-based caller: a call is a marker read against a phenotype and establishes no mechanism. TBProfiler is built for sequencing reads and was run here on assemblies, which cannot show a mixed population, and it gave no call in 819 of 21,737 genome-drug rows where an assembly did not cover a target region in full; these figures say nothing about the tool on reads.

## Reading

GRM ran TBProfiler 6.6.3, a public catalogue-based caller, with its shipped database on assemblies of the 2,000 genomes in its 22 revision 14 M. tuberculosis reports, and counted its calls against the recorded phenotypes on the 16 reports whose drug it names, beside each of the five learners' calls on the same genomes held out of that learner's fit.

On the 16 M. tuberculosis reports whose drug it names, TBProfiler made a call on 20,918 of 21,737 genome-drug rows and none on the other 819, with median sensitivity 0.779 (0.000 to 0.961) beside median specificity 0.848 (0.684 to 1.000); on the 14 within these 16 whose drug its database grades (none for the other 2), with a call on 18,869 of 19,686 rows, they were 0.802 (0.009 to 0.961) and 0.829 (0.684 to 0.997).

On the lineage split and the same genomes, TBProfiler's balanced accuracy was higher than scm's on 8 of the 14 M. tuberculosis reports whose drug its database grades, cart's on 9, rf's on 13, gbt's on 10 and lasso's on 10, and lower on the rest; its catalogue may include some of these genomes, which would favour it.

TBProfiler's figures are retrospective counts from assemblies, which cannot show a mixed population and left 817 of 19,686 rows of the 14 graded reports without a call, so they say nothing about the tool on sequencing reads; its catalogue may include some of these genomes, which would favour it, and each learner saw only genomes held out of its fit.

## Plan and order

The plan (sha256 `a12a03b3`) was hashed at 16:42:57Z on 14 Nov and its addendum (sha256 `099af373`) at 16:51:41Z, both before any recorded phenotype was read against a call, by the hash list's own note. The scoring run checked both hashes at 17:28:38Z, before the evaluator read any phenotype.

The plan fixed `who_assoc` as the primary way of reading, the rule that scores a report, sensitivity beside specificity per report, and the medians over the 16 scored reports whose drug the tool names, with the same for each learner on the same genomes, each learner minus the tool, and the lineage interval of each.

The shipped database grades no entry for cycloserine or para-aminosalicylic acid (`grades.tsv`), so the tool cannot call either resistant; this was written into the plan's addendum, hashed at 16:51:41Z, before the evaluator first read a recorded phenotype after 17:28:42Z. The addendum added the 14 scored reports whose drug the database grades as a second set beside the plan's 16, and changed no figure the plan defines.

| Event | Time (UTC) | Source |
|---|---|---|
| The plan hashed | 2026-11-14 16:42:57Z | the plan's hash list |
| The plan's addendum hashed | 2026-11-14 16:51:41Z | the plan's hash list |
| Both hashes checked by the scoring run | 2026-11-14 17:28:38Z | the log of the scoring run |
| Last step before the evaluator read a phenotype | 2026-11-14 17:28:42Z | the log of the scoring run |
| Scoring done | 2026-11-14 17:29:39Z | the log of the scoring run |

## Method

TBProfiler read each assembly with the database it ships, `tbdb` commit `3ed1cb99`, whose graded entries carry the source label "WHO catalogue v2". The tool names 18 drugs. A report is primary when the tool has a drug of the same name, with rifampin read as rifampicin. The four ways of reading its output:

- **who_assoc**, the primary way: resistant when the tool reports an entry graded Assoc w R or Assoc w R - Interim for the drug.
- **who_assoc_only**: resistant on an entry graded Assoc w R alone.
- **tool_default**: resistant on any entry the tool reports for the drug, the database's own ungraded and Uncertain significance entries included, as its printed report does.
- **who_assoc_small**: as who_assoc, leaving out an insertion or deletion of 50 bases or more or the loss of a whole coding region, which a contig end of an assembly can mimic.

Each genome and report gets one call: resistant, susceptible or no call, the last where the way of reading finds no resistant entry and an assembly did not cover a target region in full; a resistant call stands whatever the coverage. A no call is never read as susceptible. Each learner and the catalog row are scored on the same genomes, those where the tool has a call.

A learner calls resistant at a report score of 0.5 or more, on the lineage split and the random split. Each learner was scored only on genomes held out of its fit, and production selects its settings inside each fit.

A report is scored when at least 25 resistant and 25 susceptible genomes have a call. Sensitivity is the share of resistant genomes called resistant, specificity the share of susceptible genomes called susceptible, and balanced accuracy their mean.

Each of the 2,000 resamples draws 881 group numbers with Python's `random.Random(20261114).randrange(881)`, one call per draw, with the groups numbered in the text order of their smallest genome id. The interval's ends are the 51st and the 1,950th of the 2,000 sorted values.

An interval needs at least 40 resamples with a value; every interval here had a value in all 2,000.

The groups are the reports' own lineage ids, joined across the 22 reports: two genomes share a group when any report puts them in one lineage. A resample takes every genome of each group drawn. A figure has no value in a resample where a class has no genome, and a median skips a report without one. Brackets on this page give this 95% lineage interval.

Every figure on this page is rounded once from the full value in the record.

The shipped database's entries of the two counted grades, per tool drug:

| Tool drug | Assoc w R | Assoc w R - Interim | Counted entries |
|---|---|---|---|
| amikacin | 2 | 2 | 4 |
| bedaquiline | 10 | 86 | 96 |
| capreomycin | 10 | 64 | 74 |
| clofazimine | 6 | 61 | 67 |
| cycloserine | 0 | 0 | 0 |
| delamanid | 0 | 54 | 54 |
| ethambutol | 13 | 0 | 13 |
| ethionamide | 10 | 281 | 291 |
| isoniazid | 12 | 135 | 147 |
| kanamycin | 6 | 2 | 8 |
| levofloxacin | 12 | 6 | 18 |
| linezolid | 1 | 7 | 8 |
| moxifloxacin | 10 | 8 | 18 |
| para-aminosalicylic acid | 0 | 0 | 0 |
| pretomanid | 0 | 30 | 30 |
| pyrazinamide | 144 | 202 | 346 |
| rifampicin | 26 | 113 | 139 |
| streptomycin | 19 | 144 | 163 |

## Primary reports

TBProfiler under `who_assoc` on the 16 primary reports, with the lineage interval in brackets. Graded says whether the database grades an entry for the report's drug.

| Report | Graded | Rows with a call | Resistant called resistant | Sensitivity | Susceptible called susceptible | Specificity | Balanced accuracy |
|---|---|---|---|---|---|---|---|
| amikacin | yes | 1,683 of 1,813 | 432 of 542 | 0.797 [0.746, 0.845] | 1,099 of 1,141 | 0.963 [0.947, 0.979] | 0.880 [0.854, 0.906] |
| capreomycin | yes | 652 of 670 | 347 of 408 | 0.850 [0.792, 0.892] | 205 of 244 | 0.840 [0.774, 0.890] | 0.845 [0.805, 0.879] |
| clofazimine | yes | 1,249 of 1,303 | 5 of 188 | 0.027 [0.007, 0.052] | 1,048 of 1,061 | 0.988 [0.980, 0.994] | 0.507 [0.496, 0.520] |
| cycloserine | no | 400 of 401 | 0 of 119 | 0.000 [0.000, 0.000] | 281 of 281 | 1.000 [1.000, 1.000] | 0.500 [0.500, 0.500] |
| delamanid | yes | 1,217 of 1,248 | 4 of 441 | 0.009 [0.002, 0.020] | 774 of 776 | 0.997 [0.993, 1.000] | 0.503 [0.499, 0.509] |
| ethambutol | yes | 1,871 of 1,947 | 991 of 1,301 | 0.762 [0.707, 0.801] | 390 of 570 | 0.684 [0.624, 0.743] | 0.723 [0.696, 0.748] |
| ethionamide | yes | 1,767 of 1,774 | 656 of 914 | 0.718 [0.646, 0.766] | 696 of 853 | 0.816 [0.737, 0.875] | 0.767 [0.719, 0.804] |
| isoniazid | yes | 1,939 of 1,954 | 1,311 of 1,381 | 0.949 [0.928, 0.966] | 478 of 558 | 0.857 [0.800, 0.907] | 0.903 [0.878, 0.925] |
| kanamycin | yes | 1,677 of 1,800 | 489 of 606 | 0.807 [0.755, 0.850] | 1,027 of 1,071 | 0.959 [0.942, 0.976] | 0.883 [0.857, 0.907] |
| levofloxacin | yes | 1,171 of 1,248 | 442 of 722 | 0.612 [0.491, 0.706] | 354 of 449 | 0.788 [0.717, 0.868] | 0.700 [0.662, 0.725] |
| linezolid | yes | 1,183 of 1,248 | 24 of 193 | 0.124 [0.057, 0.202] | 982 of 990 | 0.992 [0.986, 0.997] | 0.558 [0.525, 0.596] |
| moxifloxacin | yes | 1,447 of 1,524 | 597 of 722 | 0.827 [0.768, 0.866] | 593 of 725 | 0.818 [0.772, 0.862] | 0.822 [0.797, 0.842] |
| para-aminosalicylic acid | no | 1,649 of 1,650 | 0 of 429 | 0.000 [0.000, 0.000] | 1,220 of 1,220 | 1.000 [1.000, 1.000] | 0.500 [0.500, 0.500] |
| pyrazinamide | yes | 486 of 492 | 321 of 377 | 0.851 [0.798, 0.898] | 83 of 109 | 0.761 [0.660, 0.874] | 0.806 [0.747, 0.868] |
| rifampin | yes | 1,830 of 1,967 | 1,293 of 1,346 | 0.961 [0.940, 0.974] | 395 of 484 | 0.816 [0.756, 0.866] | 0.888 [0.862, 0.910] |
| streptomycin | yes | 697 of 698 | 495 of 538 | 0.920 [0.871, 0.951] | 117 of 159 | 0.736 [0.598, 0.871] | 0.828 [0.761, 0.893] |

## Proxy reports

The plan read 4 reports from another drug with a shared name stem, and the call is resistant when any drug read is. They enter no set and no median.

| Report | Read from | Rows with a call | Resistant called resistant | Sensitivity | Susceptible called susceptible | Specificity | Balanced accuracy |
|---|---|---|---|---|---|---|---|
| ciprofloxacin | levofloxacin and moxifloxacin | 101 of 101 | 65 of 73 | 0.890 [0.785, 0.971] | 22 of 28 | 0.786 [0.529, 0.947] | 0.838 [0.707, 0.927] |
| ofloxacin | levofloxacin and moxifloxacin | 643 of 643 | 321 of 355 | 0.904 [0.867, 0.935] | 266 of 288 | 0.924 [0.884, 0.956] | 0.914 [0.887, 0.935] |
| prothionamide | ethionamide | 168 of 170 | 66 of 91 | 0.725 [0.578, 0.829] | 52 of 77 | 0.675 [0.495, 0.850] | 0.700 [0.616, 0.792] |
| rifabutin | rifampicin | 1,182 of 1,317 | 645 of 663 | 0.973 [0.948, 0.987] | 379 of 519 | 0.730 [0.665, 0.789] | 0.852 [0.822, 0.875] |

The tool names no drug for 2 reports, amoxicillin (57 genomes) and nicotinamide (166 genomes). They are counted here and enter no table of calls and no median.

## Medians

Each cell gives the median over the set's reports, the range in parentheses and the lineage interval in brackets. The plan fixed the set of 16, and the addendum added the 14. Each median is taken on its own, so a median sensitivity and a median specificity need not come from one report; the tables above give both per report.

The catalogue behind TBProfiler was built from public isolates with recorded phenotypes, and some of these 2,000 genomes may be among them, which would favour the tool; the overlap was not checked, and each learner was scored only on genomes held out of its fit.

| Caller | Set | Sensitivity | Specificity | Balanced accuracy |
|---|---|---|---|---|
| TBProfiler | the plan's 16 | 0.779 (0.000 to 0.961) [0.734, 0.806] | 0.848 (0.684 to 1.000) [0.821, 0.892] | 0.787 (0.500 to 0.903) [0.744, 0.809] |
| TBProfiler | the addendum's 14 | 0.802 (0.009 to 0.961) [0.762, 0.830] | 0.829 (0.684 to 0.997) [0.804, 0.873] | 0.814 (0.503 to 0.903) [0.777, 0.833] |
| scm, lineage split | the plan's 16 | 0.800 (0.487 to 0.949) [0.724, 0.829] | 0.785 (0.223 to 0.952) [0.738, 0.839] | 0.697 (0.521 to 0.898) [0.677, 0.725] |
| scm, lineage split | the addendum's 14 | 0.813 (0.576 to 0.949) [0.755, 0.840] | 0.785 (0.223 to 0.952) [0.737, 0.837] | 0.755 (0.521 to 0.898) [0.714, 0.776] |
| cart, lineage split | the plan's 16 | 0.731 (0.479 to 0.941) [0.683, 0.796] | 0.764 (0.438 to 0.956) [0.718, 0.815] | 0.716 (0.542 to 0.900) [0.696, 0.736] |
| cart, lineage split | the addendum's 14 | 0.731 (0.479 to 0.941) [0.677, 0.796] | 0.793 (0.513 to 0.956) [0.749, 0.847] | 0.750 (0.542 to 0.900) [0.706, 0.771] |
| rf, lineage split | the plan's 16 | 0.767 (0.000 to 0.960) [0.716, 0.793] | 0.802 (0.266 to 0.998) [0.755, 0.838] | 0.701 (0.499 to 0.868) [0.672, 0.718] |
| rf, lineage split | the addendum's 14 | 0.777 (0.000 to 0.960) [0.740, 0.822] | 0.788 (0.266 to 0.998) [0.741, 0.830] | 0.727 (0.499 to 0.868) [0.697, 0.756] |
| gbt, lineage split | the plan's 16 | 0.724 (0.372 to 0.951) [0.682, 0.768] | 0.747 (0.526 to 0.951) [0.719, 0.810] | 0.730 (0.506 to 0.901) [0.702, 0.749] |
| gbt, lineage split | the addendum's 14 | 0.765 (0.372 to 0.951) [0.721, 0.801] | 0.747 (0.526 to 0.951) [0.720, 0.796] | 0.743 (0.506 to 0.901) [0.713, 0.773] |
| lasso, lineage split | the plan's 16 | 0.713 (0.027 to 0.917) [0.616, 0.773] | 0.815 (0.477 to 0.975) [0.764, 0.850] | 0.723 (0.501 to 0.899) [0.702, 0.749] |
| lasso, lineage split | the addendum's 14 | 0.770 (0.027 to 0.917) [0.683, 0.807] | 0.821 (0.477 to 0.975) [0.776, 0.850] | 0.766 (0.501 to 0.899) [0.726, 0.783] |
| scm, random split | the plan's 16 | 0.794 (0.507 to 0.950) [0.735, 0.826] | 0.808 (0.310 to 0.948) [0.778, 0.852] | 0.713 (0.542 to 0.901) [0.691, 0.737] |
| scm, random split | the addendum's 14 | 0.814 (0.507 to 0.950) [0.760, 0.839] | 0.830 (0.310 to 0.948) [0.787, 0.866] | 0.771 (0.542 to 0.901) [0.747, 0.790] |
| cart, random split | the plan's 16 | 0.769 (0.624 to 0.945) [0.745, 0.808] | 0.742 (0.382 to 0.954) [0.705, 0.807] | 0.724 (0.531 to 0.900) [0.705, 0.743] |
| cart, random split | the addendum's 14 | 0.805 (0.624 to 0.945) [0.764, 0.832] | 0.786 (0.382 to 0.954) [0.739, 0.842] | 0.772 (0.531 to 0.900) [0.750, 0.790] |
| rf, random split | the plan's 16 | 0.774 (0.021 to 0.952) [0.740, 0.811] | 0.801 (0.303 to 0.994) [0.769, 0.839] | 0.726 (0.508 to 0.879) [0.697, 0.747] |
| rf, random split | the addendum's 14 | 0.797 (0.021 to 0.952) [0.764, 0.836] | 0.792 (0.303 to 0.994) [0.758, 0.828] | 0.758 (0.508 to 0.879) [0.726, 0.784] |
| gbt, random split | the plan's 16 | 0.790 (0.298 to 0.945) [0.737, 0.820] | 0.800 (0.563 to 0.949) [0.759, 0.838] | 0.753 (0.536 to 0.900) [0.730, 0.771] |
| gbt, random split | the addendum's 14 | 0.813 (0.298 to 0.945) [0.768, 0.837] | 0.807 (0.563 to 0.949) [0.765, 0.840] | 0.793 (0.536 to 0.900) [0.766, 0.804] |
| lasso, random split | the plan's 16 | 0.800 (0.165 to 0.912) [0.736, 0.825] | 0.802 (0.394 to 0.946) [0.775, 0.847] | 0.753 (0.518 to 0.902) [0.726, 0.772] |
| lasso, random split | the addendum's 14 | 0.827 (0.165 to 0.912) [0.771, 0.841] | 0.805 (0.394 to 0.946) [0.787, 0.851] | 0.785 (0.518 to 0.902) [0.756, 0.797] |

## Learners beside the tool

On the lineage split and the same genomes, TBProfiler's balanced accuracy was ahead of or within 0.02 of the best of the five learners on 12 of 16 endpoints and 12 of the 14 graded ones; the 0.02 margin is an adapted research screen the plan proposed and no one has adopted, and it is never an applicable FDA validation criterion.

The catalogue behind TBProfiler was built from public isolates with recorded phenotypes, and some of these 2,000 genomes may be among them, which would favour the tool; the overlap was not checked, and each learner was scored only on genomes held out of its fit.

Per learner and split: on how many reports the tool's balanced accuracy is higher or lower, on how many the tool is ahead or the learner at most 0.02 ahead, the median of the learner minus the tool with its lineage interval, and the genome-drug rows only one of the two gets right, summed over the set.

| Learner | Set | Tool higher | Tool lower | Tool ahead or within 0.02 | Learner minus tool, median | Right only for the learner | Right only for the tool |
|---|---|---|---|---|---|---|---|
| scm, lineage split | the plan's 16 | 8 | 8 | 13 | +0.001 [-0.009, +0.011] | 1,323 | 2,751 |
| scm, lineage split | the addendum's 14 | 8 | 6 | 13 | -0.006 [-0.012, +0.006] | 1,049 | 2,507 |
| cart, lineage split | the plan's 16 | 9 | 7 | 12 | -0.004 [-0.010, +0.007] | 1,389 | 2,772 |
| cart, lineage split | the addendum's 14 | 9 | 5 | 12 | -0.007 [-0.013, +0.002] | 969 | 2,202 |
| rf, lineage split | the plan's 16 | 13 | 3 | 13 | -0.042 [-0.046, -0.022] | 1,134 | 1,483 |
| rf, lineage split | the addendum's 14 | 13 | 1 | 13 | -0.044 [-0.052, -0.027] | 931 | 1,323 |
| gbt, lineage split | the plan's 16 | 10 | 6 | 13 | -0.003 [-0.013, +0.002] | 1,330 | 2,198 |
| gbt, lineage split | the addendum's 14 | 10 | 4 | 13 | -0.005 [-0.018, -0.001] | 994 | 1,881 |
| lasso, lineage split | the plan's 16 | 10 | 6 | 13 | -0.005 [-0.013, +0.005] | 1,143 | 1,848 |
| lasso, lineage split | the addendum's 14 | 10 | 4 | 13 | -0.006 [-0.017, +0.001] | 812 | 1,553 |
| scm, random split | the plan's 16 | 8 | 8 | 11 | +0.003 [-0.004, +0.014] | 1,364 | 2,746 |
| scm, random split | the addendum's 14 | 8 | 6 | 11 | -0.001 [-0.007, +0.009] | 1,045 | 2,428 |
| cart, random split | the plan's 16 | 7 | 9 | 10 | +0.008 [-0.002, +0.013] | 1,500 | 2,853 |
| cart, random split | the addendum's 14 | 7 | 7 | 10 | +0.003 [-0.005, +0.009] | 1,105 | 2,406 |
| rf, random split | the plan's 16 | 12 | 4 | 13 | -0.019 [-0.028, -0.009] | 1,302 | 1,398 |
| rf, random split | the addendum's 14 | 12 | 2 | 13 | -0.024 [-0.032, -0.012] | 1,038 | 1,249 |
| gbt, random split | the plan's 16 | 5 | 11 | 10 | +0.005 [-0.002, +0.017] | 1,402 | 1,973 |
| gbt, random split | the addendum's 14 | 5 | 9 | 10 | +0.003 [-0.005, +0.012] | 1,036 | 1,643 |
| lasso, random split | the plan's 16 | 9 | 7 | 12 | -0.001 [-0.008, +0.012] | 1,252 | 1,666 |
| lasso, random split | the addendum's 14 | 9 | 5 | 12 | -0.001 [-0.011, +0.007] | 873 | 1,345 |

The best lineage learner on each report of the plan's set is picked after the fact, which favours the learners. Gap is its balanced accuracy minus the tool's.

| Report | Graded | TBProfiler | Best lineage learner | Its balanced accuracy | Gap | Tool ahead or within 0.02 |
|---|---|---|---|---|---|---|
| pyrazinamide | yes | 0.806 | gbt | 0.640 | -0.1667 | yes |
| ethionamide | yes | 0.767 | gbt | 0.732 | -0.0344 | yes |
| streptomycin | yes | 0.828 | lasso | 0.813 | -0.0146 | yes |
| levofloxacin | yes | 0.700 | lasso | 0.695 | -0.0056 | yes |
| ethambutol | yes | 0.723 | lasso | 0.718 | -0.0048 | yes |
| moxifloxacin | yes | 0.822 | lasso | 0.819 | -0.0034 | yes |
| isoniazid | yes | 0.903 | gbt | 0.901 | -0.0020 | yes |
| capreomycin | yes | 0.845 | gbt | 0.851 | +0.0058 | yes |
| rifampin | yes | 0.888 | scm | 0.898 | +0.0096 | yes |
| amikacin | yes | 0.880 | scm | 0.890 | +0.0100 | yes |
| linezolid | yes | 0.558 | lasso | 0.573 | +0.0148 | yes |
| kanamycin | yes | 0.883 | lasso | 0.899 | +0.0161 | yes |
| clofazimine | yes | 0.507 | cart | 0.558 | +0.0508 | no |
| delamanid | yes | 0.503 | gbt | 0.572 | +0.0692 | no |
| cycloserine | no | 0.500 | lasso | 0.613 | +0.1133 | no |
| para-aminosalicylic acid | no | 0.500 | cart | 0.729 | +0.2286 | no |

GRM's own reader of the AMRFinderPlus catalog has a row on 6 of the 16 reports. On them, the tool's balanced accuracy minus the catalog row's has median 0.325. Of the 7,858 genome-drug rows with a call there, the two calls agree on 4,931, only the tool is right on 2,386 and only the catalog row on 541, which sum to 7,858.

## Ways of reading

The catalogue behind TBProfiler was built from public isolates with recorded phenotypes, and some of these 2,000 genomes may be among them, which would favour the tool; the overlap was not checked, and each learner was scored only on genomes held out of its fit.

The tool's medians and ranges under each way of reading over each set, and the screen's count. `who_assoc` is the primary way and the only one the plan's screen reads; the other three show how the calls move with the grades counted or with large changes left out; their counts in the last column are exploratory, and the plan named no screen for them. The screen and its label stand in the section Learners beside the tool.

| Way of reading | Set | TBProfiler sensitivity | TBProfiler specificity | Ahead of or within 0.02 of the best lineage learner |
|---|---|---|---|---|
| who_assoc | the plan's 16 | 0.779 (0.000 to 0.961) | 0.848 (0.684 to 1.000) | 12 of 16 |
| who_assoc | the addendum's 14 | 0.802 (0.009 to 0.961) | 0.829 (0.684 to 0.997) | 12 of 14 |
| who_assoc_only | the plan's 16 | 0.763 (0.000 to 0.954) | 0.857 (0.684 to 1.000) | 10 of 16 |
| who_assoc_only | the addendum's 14 | 0.777 (0.000 to 0.954) | 0.847 (0.684 to 1.000) | 10 of 14 |
| tool_default | the plan's 16 | 0.813 (0.009 to 0.961) | 0.847 (0.640 to 0.997) | 13 of 16 |
| tool_default | the addendum's 14 | 0.830 (0.009 to 0.961) | 0.828 (0.640 to 0.997) | 12 of 14 |
| who_assoc_small | the plan's 16 | 0.779 (0.000 to 0.961) | 0.848 (0.684 to 1.000) | 12 of 16 |
| who_assoc_small | the addendum's 14 | 0.798 (0.005 to 0.961) | 0.829 (0.684 to 1.000) | 12 of 14 |

## Disagreements by grade

Under `who_assoc`, per primary report: the resistant calls on a susceptible genome by the grade they rest on, the susceptible calls on a resistant genome by what the tool reported for the drug, and the genomes with no call by label. An entry of another grade is one the tool reports that `who_assoc` does not count.

| Report | Resistant call on Assoc w R, label susceptible | Resistant call on Interim alone, label susceptible | Susceptible call with an entry of another grade, label resistant | Susceptible call with nothing reported, label resistant | No call, label resistant | No call, label susceptible |
|---|---|---|---|---|---|---|
| amikacin | 30 | 12 | 0 | 110 | 6 | 124 |
| capreomycin | 39 | 0 | 0 | 61 | 14 | 4 |
| clofazimine | 12 | 1 | 0 | 183 | 5 | 49 |
| cycloserine | 0 | 0 | 42 | 77 | 0 | 1 |
| delamanid | 0 | 2 | 0 | 437 | 7 | 24 |
| ethambutol | 180 | 0 | 85 | 225 | 37 | 39 |
| ethionamide | 124 | 33 | 11 | 247 | 1 | 6 |
| isoniazid | 79 | 1 | 10 | 60 | 7 | 8 |
| kanamycin | 44 | 0 | 15 | 102 | 11 | 112 |
| levofloxacin | 92 | 3 | 3 | 277 | 31 | 46 |
| linezolid | 6 | 2 | 0 | 169 | 9 | 56 |
| moxifloxacin | 126 | 6 | 7 | 118 | 11 | 66 |
| para-aminosalicylic acid | 0 | 0 | 153 | 276 | 0 | 1 |
| pyrazinamide | 21 | 5 | 20 | 36 | 4 | 2 |
| rifampin | 88 | 1 | 0 | 53 | 19 | 118 |
| streptomycin | 41 | 1 | 13 | 30 | 1 | 0 |

## Exploratory

**Exploratory.** The figures in this section were read after the outcomes and are in no plan.

Counting each no call as wrong or as right gives bounds on every genome. Reading a no call as susceptible, as the tool's printed report does, gives the lower sensitivity bound and the upper specificity bound.

| Set | Sensitivity, no call counted wrong | Sensitivity, no call counted right | Specificity, no call counted wrong | Specificity, no call counted right |
|---|---|---|---|---|
| the plan's 16 | 0.764 (0.000 to 0.947) | 0.784 (0.000 to 0.961) | 0.836 (0.640 to 0.999) | 0.855 (0.704 to 1.000) |
| the addendum's 14 | 0.790 (0.009 to 0.947) | 0.805 (0.025 to 0.961) | 0.818 (0.640 to 0.968) | 0.847 (0.704 to 0.998) |

**Exploratory.** On the kanamycin report, the tool's rule left its note in 4 result files and set aside 1 entry. As run, the tool called 489 of 606 resistant genomes resistant and 1,027 of 1,071 susceptible genomes susceptible (sensitivity 0.8069, specificity 0.9589). With the entry counted it calls 1,026 of 1,071 susceptible genomes susceptible (sensitivity 0.8069, specificity 0.9580), and the best lineage learner minus the tool moves from +0.0161 to +0.0165.

The catalogue behind TBProfiler was built from public isolates with recorded phenotypes, and some of these 2,000 genomes may be among them, which would favour the tool; the overlap was not checked, and each learner was scored only on genomes held out of its fit.

These counts rest on the tool's raw output alone and cannot be rebuilt from the call file.

## Review

An independent review recounted every row of the first report's table from the call file and the reports and found no blocker. It narrowed the first report in these ways, and this page takes each in the review's form.

- The plan's words alone do not fix the lineage interval. The review read the evaluator's two choices, which the section Method states, and this record's script gives the review's ends.
- Some cells of the first report were rounded twice. Every figure here is rounded once.
- The word the first report used for the screen's outcome overstates a screen no one has adopted. The section Learners beside the tool gives the review's wording, and the screen's count has no interval here.
- The first report's account of how the two drugs without a graded entry were found rests on an unlogged run. This page gives the fact from the database's grade table.
- The catalogue's possible overlap with these genomes goes beside every figure that sets the tool beside a learner.

## Limits

- The tool's raw output and the assemblies it read no longer exist; the per-genome call file, kept outside this repository, is the lowest layer these counts can be rebuilt from.
- The tool read assemblies. It is built for sequencing reads, and an assembly cannot show a mixed population.
- The overlap of the catalogue's source isolates with these genomes was not counted.
- The best of the five learners is picked per report after the fact, which favours the learners.
- A no call leaves a genome out for every caller. The section Exploratory gives the tool's bounds on every genome.
- One cutoff: a learner's score of 0.5 or more calls resistant.
- The record leaves out these figures the plan named: the exact per-report interval, which treats related genomes as independent, the McNemar and sign tests, the quartiles of the differences, the per-fold measures, and sensitivity and specificity over all genome-drug rows taken together.
- One database release, `tbdb` commit `3ed1cb99`, and one tool version.
- The page gives no figure for a single genome, and it names no genome.

## Files

- `reference/tbprofiler-comparison/` holds the record: 20 files that `scripts/tbprofiler-comparison/record.py` writes from the input set. Its README lists each file and its columns, and `SHA256SUMS` gives each file's sha256.
- The record holds rows per report and caller, per way of reading, per tool drug, per grade and per learner, and single counts. It holds no genome id, label, call or score per genome.
- `checks.tsv` gives 0 as the count of the input set's genome ids found in the record.
- The input set holds the call file and the reports' scores, which are per genome, so the repository does not hold it. `--check` alone derives every derived file from the row files.
- `tests/tbprofiler_comparison` rebuilds the section Reading, three paragraphs and every table from the record, holds every other paragraph and list item word for word, and checks that each other number stands in the record.
