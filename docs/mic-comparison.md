# MIC comparison

This page asks how well three ways of predicting a genome's MIC hold when related genomes are kept apart. It reads block C6 part B of the assumption plan (`docs/assumption-plan.md`). It is an experiment beside the release, and it changes nothing in the release or on the site.

The three MIC models were fitted again on lineage folds and on study folds. The genomes are public BV-BRC records of release `rev14-274-20261112`, evaluation revision 14. The main measure is EA-exact: exact MICs whose predicted step is within one dilution, over the exact MICs tested.

The reading's rules were fixed at 14:43Z on 22 Nov 2026. All 12 species had run by then, and no score had been read.

## Reading

**Release rule.** On the lineage split, the released model's EA-exact minus the training median's is under 0.02 on 53 of the 100 datasets. That is more than half, so the condition of part B's release rule holds. The maintainer rules what the site says.

**Lineage pairs.** The released model minus the training median is better on 47 of 100 datasets, the same on 19 and worse on 34. The ordered model minus the released one is better on 41 of the 100 datasets, the same on 31 and worse on 28. No lineage interval lies wholly above 0, and no p value is under 0.05 after Holm.

**Twin.** The released model at other seeds is better on 22 of the 100 datasets, the same on 56 and worse on 22. Its median difference is 0.000.

**Study split.** By study, the released model is worse than the training median on 17 of 28 datasets. The median difference is -0.055 (-0.133 to -0.027). The study split decides nothing here.

## Data

The 12 species are *Acinetobacter baumannii*, *Campylobacter coli*, *Campylobacter jejuni*, *Enterococcus faecium*, *Escherichia coli*, *Klebsiella pneumoniae*, *Mycobacterium tuberculosis*, *Neisseria gonorrhoeae*, *Pseudomonas aeruginosa*, *Salmonella enterica*, *Staphylococcus aureus* and *Streptococcus pneumoniae*.

**Lineage split.** 100 datasets: the plan's 132 MIC datasets less the 32 marked M. They test 66,116 MICs, one per genome and antibiotic. Each dataset ran at seeds 42 to 47, so the split holds 600 results.

**Study split.** 28 datasets in 6 species, with 13,224 MICs and 168 results. That is the MIC list of the study split less *Klebsiella pneumoniae* with ampicillin/sulbactam, whose seed-47 study plan failed its check (`docs/assumption-plan-corrections.md`, MIC list).

The 768 results are the 600 and the 168 together, and all passed the scorer's checks. Each dataset and split holds 15 score rows, 1,920 in all.

| Arm | Model | Seeds |
|---|---|---|
| released | the release's censored MIC regression, on the plan's folds | 42, 43, 44 |
| ordered | one boosted fit per dilution threshold | 42, 43, 44 |
| median | each fold's training median, the lower median of its exact MICs | 42, 43, 44 |
| twin | the released model again | 45, 46, 47 |

A dataset's value for an arm is the mean of its three draws. On the lineage split, the twin changes only the seed, for the folds and for each fit. On the study split, the fits keep the release's seed and the study plan changes from seed to seed. The seed-47 plan is coarser than the other five.

Each pair is the first arm minus the second: P1 released minus median, P2 ordered minus released, P3 twin minus released. A difference above 0.02 counts as better, one below -0.02 as worse, and the rest as the same.

## Pairs

| Split | Pair | Compared | Better | Same | Worse | Median | Interval | Sign p | Holm p |
|---|---|---:|---:|---:|---:|---:|---|---:|---:|
| lineage | P1 released minus median | 100 of 100 | 47 | 19 | 34 | 0.011 | -0.010 to 0.063 | 0.091 | 0.222 |
| lineage | P2 ordered minus released | 100 of 100 | 41 | 31 | 28 | 0.006 | -0.011 to 0.024 | 0.074 | 0.222 |
| lineage | P3 twin minus released | 100 of 100 | 22 | 56 | 22 | 0.000 | -0.004 to 0.004 | 0.560 | 0.560 |
| study | P1 released minus median | 28 of 28 | 10 | 1 | 17 | -0.055 | -0.133 to -0.027 | 0.939 | |
| study | P2 ordered minus released | 28 of 28 | 9 | 3 | 16 | -0.025 | -0.046 to 0.024 | 0.946 | |
| study | P3 twin minus released | 28 of 28 | 10 | 12 | 6 | 0.007 | 0.002 to 0.051 | 0.227 | |

**Median.** The median difference over the datasets compared. On the lineage split it is above 0 for P1 and P2 and 0 for P3. On the study split it is below 0 for P1 and P2 and above 0 for P3.

**Interval.** A 95% interval from 4,000 resamples of species, read at the 101st and the 3,900th sorted medians. It uses the method and generator of the AMRFinderPlus record (`docs/amrfinderplus-comparison.md`). The lineage split has 12 species and the study split 6.

**Sign test.** One-sided and exact: is better more often than worse? The same are left out. Holm runs over the three lineage pairs. The study rows are shown beside them and decide nothing.

**Release rule.** The plan's rule for part B: if the model's lineage EA-exact beats the training median by under 0.02 on most datasets, the site shows the lineage reading beside the random one.

On the lineage split, 53 of the 100 datasets compared have a P1 difference under 0.02: the 34 worse and the 19 the same. 53 is more than half of 100, so the condition holds. This page changes nothing in the release or on the site. The maintainer rules that wording.

## Beside it

These read the R and S labels through each fold's training breakpoint. A genome is called resistant when its predicted step is above its fold's cut-off.

**Category agreement.** Genomes whose call equals the reported call, over the MIC genomes tested.

**Very major error rate.** Resistant genomes called susceptible, over the resistant MIC genomes tested.

**Major error rate.** Susceptible genomes called resistant, over the susceptible MIC genomes tested.

Each cell gives the median over datasets of the three-draw mean. Then, in brackets, the count pooled over the datasets and the three draws, over its pooled denominator, so a genome counts once a draw.

Lineage split, defined on 91 of 100 datasets for every arm:

| Arm | Agreement | Very major | Major |
|---|---|---|---|
| released | 0.867 (163,679 of 189,510) | 0.147 (13,598 of 75,798) | 0.087 (12,233 of 113,712) |
| ordered | 0.886 (166,942 of 189,510) | 0.158 (13,727 of 75,798) | 0.047 (8,841 of 113,712) |
| median | 0.641 (107,470 of 189,510) | 1.000 (46,107 of 75,798) | 0.000 (35,933 of 113,712) |
| twin | 0.867 (164,341 of 189,510) | 0.154 (13,488 of 75,798) | 0.084 (11,681 of 113,712) |

Study split, defined on 24 of 28 datasets, and 23 of 28 for the twin:

| Arm | Agreement | Very major | Major |
|---|---|---|---|
| released | 0.744 (29,183 of 36,342) | 0.124 (2,817 of 22,212) | 0.272 (4,342 of 14,130) |
| ordered | 0.772 (29,554 of 36,342) | 0.098 (2,436 of 22,212) | 0.298 (4,352 of 14,130) |
| median | 0.536 (19,132 of 36,342) | 0.989 (12,866 of 22,212) | 0.088 (4,344 of 14,130) |
| twin | 0.755 (28,661 of 35,709) | 0.098 (2,694 of 22,089) | 0.367 (4,354 of 13,620) |

**Left out of a pair.** None. On both splits every draw of every arm was available, left no fold undefined and tested at least one exact MIC.

**Left out beside it.** The three figures are not defined where a scored fold's training genomes gave no breakpoint. On the lineage split that holds for 9 of the 100 datasets. They are *Acinetobacter baumannii* with cefepime, *Campylobacter coli* with gentamicin, *Campylobacter jejuni* with clindamycin, *Escherichia coli* with nitrofurantoin, *Klebsiella pneumoniae* with tigecycline, *Neisseria gonorrhoeae* with ciprofloxacin, *Salmonella enterica* with ceftazidime, and *Staphylococcus aureus* with clindamycin and with daptomycin.

On the study split it holds for 4 of the 28 datasets: *Acinetobacter baumannii* with cefepime, *Escherichia coli* with nitrofurantoin, *Klebsiella pneumoniae* with tigecycline and *Staphylococcus aureus* with clindamycin. The twin also lacks them for *Staphylococcus aureus* with erythromycin.

## Limits

**Labels.** No dataset here is marked M, so each result is the same under label rules (a), (b) and (d). EA-exact reads MIC values only. Category agreement and both error rates read the R and S labels and each fold's training breakpoint.

**Rule (c).** Under the most resistant call, a dataset holding a pair with an R and an S record is checked again. This page does not name those datasets.

**Datasets marked M.** The other 32 MIC datasets wait for the maintainer's label rule.

**Random folds.** This reading has no random split. It does not compare its figures with the release's random-fold MIC figures.

**Study twin.** On the study split, the twin changes the study plan and keeps the seed. Its count is not the count of seeds alone, and its interval lies above 0.

**Wilson intervals.** The reading's rules give none for the figures beside EA-exact.

**Causes.** It does not say why a model reads better or worse on a dataset. It finds, names or ranks no resistance marker.

**Other databases.** Every figure here comes from a lineage split and a retrospective split by study inside one database, BV-BRC.

## Method

`scripts/c6-mic/compare.py` reads each species' `scores.tsv` and `checks.tsv` from `score.py`, with the run's `expect.tsv`. It stops on a missing row, a row held twice or a check that is not ok. Values are exact fractions, so a difference of exactly 0.02 counts as the same. The interval's generator is seeded with the string `c6read-20261122|<split>|<pair>|ea_exact`, fixed at 14:48Z on 22 Nov 2026, before any score was read. The outputs hold one row per dataset and split and one row per pair and split.

```
python3 -I scripts/c6-mic/compare.py --out OUT --results 768 --lineage 100 --study 28 --species 12 SCORE_DIR...
```

## Datasets

Each row is one dataset, an aggregate over its genomes. MIC genomes are the genomes with an MIC that the plan tests. Exact MICs are EA-exact's denominator at each draw. The four arm columns give EA-exact as the mean of three draws. P1 to P3 give the class of each pair.

**Lineage split.** 100 datasets.

| Species | Antibiotic | MIC genomes | Exact MICs | Released | Ordered | Median | Twin | P1 | P2 | P3 |
|---|---|---:|---:|---:|---:|---:|---:|---|---|---|
| *Acinetobacter baumannii* | cefepime | 162 | 53 | 0.899 | 0.711 | 0.195 | 0.937 | better | worse | better |
| *Acinetobacter baumannii* | cefotaxime | 697 | 72 | 0.491 | 0.250 | 0.722 | 0.546 | worse | worse | better |
| *Acinetobacter baumannii* | doripenem | 213 | 107 | 0.919 | 0.963 | 0.981 | 0.822 | worse | better | worse |
| *Acinetobacter baumannii* | sulfamethoxazole/trimethoprim | 783 | 229 | 0.852 | 0.847 | 0.860 | 0.898 | same | same | better |
| *Campylobacter coli* | azithromycin | 308 | 266 | 0.921 | 0.969 | 0.974 | 0.916 | worse | better | same |
| *Campylobacter coli* | ciprofloxacin | 309 | 295 | 0.949 | 0.974 | 0.695 | 0.938 | better | better | same |
| *Campylobacter coli* | clindamycin | 308 | 304 | 0.861 | 0.894 | 0.715 | 0.870 | better | better | same |
| *Campylobacter coli* | erythromycin | 309 | 271 | 0.909 | 0.937 | 0.900 | 0.919 | same | better | same |
| *Campylobacter coli* | gentamicin | 309 | 297 | 0.967 | 0.982 | 0.983 | 0.972 | same | same | same |
| *Campylobacter coli* | nalidixic acid | 294 | 165 | 0.984 | 0.986 | 0.933 | 0.986 | better | same | same |
| *Campylobacter coli* | telithromycin | 294 | 278 | 0.920 | 0.875 | 0.845 | 0.922 | better | worse | same |
| *Campylobacter coli* | tetracycline | 309 | 134 | 0.789 | 0.784 | 0.470 | 0.821 | better | same | better |
| *Campylobacter jejuni* | ciprofloxacin | 477 | 470 | 0.918 | 0.950 | 0.787 | 0.904 | better | better | same |
| *Campylobacter jejuni* | clindamycin | 477 | 449 | 0.915 | 0.940 | 0.958 | 0.933 | worse | better | same |
| *Campylobacter jejuni* | nalidixic acid | 475 | 84 | 0.964 | 0.956 | 0.905 | 0.972 | better | same | same |
| *Campylobacter jejuni* | tetracycline | 477 | 258 | 0.863 | 0.907 | 0.570 | 0.827 | better | better | worse |
| *Enterococcus faecium* | ampicillin | 1219 | 353 | 0.844 | 0.873 | 0.502 | 0.811 | better | better | worse |
| *Enterococcus faecium* | dalfopristin/quinupristin | 804 | 765 | 0.925 | 0.931 | 0.906 | 0.923 | same | same | same |
| *Enterococcus faecium* | teicoplanin | 902 | 117 | 0.556 | 0.533 | 0.479 | 0.570 | better | worse | same |
| *Enterococcus faecium* | vancomycin | 917 | 86 | 0.938 | 0.930 | 0.988 | 0.922 | worse | same | same |
| *Escherichia coli* | amikacin | 472 | 87 | 0.759 | 0.617 | 0.828 | 0.785 | worse | worse | better |
| *Escherichia coli* | aztreonam | 184 | 54 | 0.630 | 0.654 | 0.833 | 0.728 | worse | better | better |
| *Escherichia coli* | cefepime | 435 | 86 | 0.543 | 0.612 | 0.651 | 0.581 | worse | better | better |
| *Escherichia coli* | cefotaxime | 175 | 59 | 0.576 | 0.571 | 0.305 | 0.661 | better | same | better |
| *Escherichia coli* | ceftazidime | 443 | 128 | 0.557 | 0.589 | 0.521 | 0.552 | better | better | same |
| *Escherichia coli* | nitrofurantoin | 287 | 63 | 0.704 | 0.688 | 0.952 | 0.661 | worse | same | worse |
| *Escherichia coli* | piperacillin/tazobactam | 265 | 62 | 0.586 | 0.871 | 0.935 | 0.817 | worse | better | better |
| *Escherichia coli* | tetracycline | 310 | 67 | 0.473 | 0.736 | 0.448 | 0.443 | better | better | worse |
| *Klebsiella pneumoniae* | ampicillin/sulbactam | 555 | 85 | 0.839 | 0.694 | 0.965 | 0.835 | worse | worse | same |
| *Klebsiella pneumoniae* | aztreonam | 802 | 108 | 0.660 | 0.546 | 0.667 | 0.651 | same | worse | same |
| *Klebsiella pneumoniae* | cefazolin | 742 | 54 | 0.574 | 0.451 | 0.759 | 0.568 | worse | worse | same |
| *Klebsiella pneumoniae* | cefepime | 789 | 83 | 0.627 | 0.558 | 0.639 | 0.618 | same | worse | same |
| *Klebsiella pneumoniae* | cefoxitin | 766 | 75 | 0.844 | 0.831 | 0.893 | 0.880 | worse | same | better |
| *Klebsiella pneumoniae* | ceftazidime | 844 | 178 | 0.625 | 0.599 | 0.646 | 0.603 | worse | worse | worse |
| *Klebsiella pneumoniae* | ciprofloxacin | 1003 | 78 | 0.739 | 0.774 | 0.833 | 0.744 | worse | better | same |
| *Klebsiella pneumoniae* | colistin | 275 | 95 | 0.309 | 0.302 | 0.256 | 0.379 | better | same | better |
| *Klebsiella pneumoniae* | ertapenem | 535 | 62 | 0.344 | 0.172 | 0.118 | 0.263 | better | worse | worse |
| *Klebsiella pneumoniae* | fosfomycin | 220 | 145 | 0.667 | 0.807 | 0.841 | 0.687 | worse | better | better |
| *Klebsiella pneumoniae* | imipenem | 903 | 174 | 0.718 | 0.707 | 0.621 | 0.755 | better | same | better |
| *Klebsiella pneumoniae* | nitrofurantoin | 356 | 54 | 0.599 | 0.284 | 0.685 | 0.574 | worse | worse | worse |
| *Klebsiella pneumoniae* | sulfamethoxazole/trimethoprim | 981 | 63 | 0.778 | 0.640 | 0.815 | 0.698 | worse | worse | worse |
| *Klebsiella pneumoniae* | tetracycline | 704 | 175 | 0.869 | 0.827 | 0.829 | 0.865 | better | worse | same |
| *Klebsiella pneumoniae* | tigecycline | 449 | 288 | 0.876 | 0.903 | 0.934 | 0.840 | worse | better | worse |
| *Klebsiella pneumoniae* | tobramycin | 754 | 140 | 0.783 | 0.800 | 0.643 | 0.781 | better | same | same |
| *Mycobacterium tuberculosis* | amikacin | 1317 | 646 | 0.913 | 0.918 | 0.918 | 0.911 | same | same | same |
| *Mycobacterium tuberculosis* | delamanid | 1248 | 1248 | 0.846 | 0.833 | 0.862 | 0.846 | same | same | same |
| *Mycobacterium tuberculosis* | ethambutol | 1298 | 1087 | 0.777 | 0.787 | 0.547 | 0.783 | better | same | same |
| *Mycobacterium tuberculosis* | ethionamide | 1317 | 1142 | 0.790 | 0.805 | 0.722 | 0.796 | better | same | same |
| *Mycobacterium tuberculosis* | kanamycin | 1317 | 728 | 0.914 | 0.900 | 0.933 | 0.908 | same | same | same |
| *Mycobacterium tuberculosis* | levofloxacin | 1248 | 1084 | 0.771 | 0.795 | 0.671 | 0.774 | better | better | same |
| *Mycobacterium tuberculosis* | linezolid | 1248 | 1225 | 0.686 | 0.695 | 0.676 | 0.690 | same | same | same |
| *Mycobacterium tuberculosis* | moxifloxacin | 1296 | 1017 | 0.780 | 0.770 | 0.562 | 0.775 | better | same | same |
| *Mycobacterium tuberculosis* | para-aminosalicylic acid | 1317 | 834 | 0.586 | 0.568 | 0.585 | 0.577 | same | same | same |
| *Mycobacterium tuberculosis* | rifabutin | 1317 | 583 | 0.728 | 0.675 | 0.611 | 0.739 | better | worse | same |
| *Mycobacterium tuberculosis* | rifampin | 1317 | 527 | 0.719 | 0.703 | 0.787 | 0.713 | worse | same | same |
| *Neisseria gonorrhoeae* | cefixime | 561 | 446 | 0.800 | 0.783 | 0.411 | 0.835 | better | same | better |
| *Neisseria gonorrhoeae* | ceftriaxone | 816 | 734 | 0.811 | 0.812 | 0.442 | 0.769 | better | same | worse |
| *Neisseria gonorrhoeae* | ciprofloxacin | 602 | 532 | 0.454 | 0.547 | 0.148 | 0.452 | better | better | same |
| *Neisseria gonorrhoeae* | erythromycin | 179 | 174 | 0.808 | 0.837 | 0.630 | 0.762 | better | better | worse |
| *Neisseria gonorrhoeae* | penicillin | 278 | 253 | 0.651 | 0.734 | 0.696 | 0.715 | worse | better | better |
| *Neisseria gonorrhoeae* | tetracycline | 359 | 356 | 0.900 | 0.869 | 0.539 | 0.889 | better | worse | same |
| *Pseudomonas aeruginosa* | amikacin | 544 | 401 | 0.638 | 0.673 | 0.623 | 0.652 | same | better | same |
| *Pseudomonas aeruginosa* | aztreonam | 192 | 113 | 0.705 | 0.702 | 0.726 | 0.717 | worse | same | same |
| *Pseudomonas aeruginosa* | cefepime | 203 | 136 | 0.684 | 0.645 | 0.654 | 0.684 | better | worse | same |
| *Pseudomonas aeruginosa* | ceftazidime | 672 | 545 | 0.415 | 0.435 | 0.317 | 0.420 | better | better | same |
| *Pseudomonas aeruginosa* | ceftazidime/avibactam | 282 | 274 | 0.457 | 0.513 | 0.383 | 0.437 | better | better | worse |
| *Pseudomonas aeruginosa* | ceftolozane/tazobactam | 234 | 234 | 0.363 | 0.511 | 0.329 | 0.370 | better | better | same |
| *Pseudomonas aeruginosa* | imipenem | 223 | 127 | 0.475 | 0.528 | 0.391 | 0.470 | better | better | same |
| *Pseudomonas aeruginosa* | levofloxacin | 539 | 322 | 0.598 | 0.618 | 0.410 | 0.623 | better | same | better |
| *Pseudomonas aeruginosa* | meropenem | 838 | 685 | 0.378 | 0.414 | 0.194 | 0.373 | better | better | same |
| *Pseudomonas aeruginosa* | piperacillin/tazobactam | 185 | 103 | 0.573 | 0.605 | 0.650 | 0.576 | worse | better | same |
| *Pseudomonas aeruginosa* | tobramycin | 505 | 347 | 0.476 | 0.522 | 0.429 | 0.500 | better | better | better |
| *Salmonella enterica* | amoxicillin/clavulanate | 399 | 124 | 0.815 | 0.890 | 0.435 | 0.812 | better | better | same |
| *Salmonella enterica* | ampicillin/sulbactam | 183 | 139 | 0.758 | 0.825 | 0.849 | 0.719 | worse | better | worse |
| *Salmonella enterica* | azithromycin | 865 | 659 | 0.906 | 0.923 | 0.912 | 0.931 | same | same | better |
| *Salmonella enterica* | cefazolin | 185 | 139 | 0.604 | 0.552 | 0.453 | 0.583 | better | worse | worse |
| *Salmonella enterica* | cefotaxime | 268 | 76 | 0.829 | 0.092 | 0.961 | 0.930 | worse | worse | better |
| *Salmonella enterica* | cefoxitin | 600 | 320 | 0.893 | 0.826 | 0.444 | 0.900 | better | worse | same |
| *Salmonella enterica* | ceftazidime | 267 | 53 | 0.692 | 0.063 | 0.849 | 0.698 | worse | worse | same |
| *Salmonella enterica* | ceftiofur | 402 | 196 | 0.944 | 0.963 | 0.367 | 0.944 | better | same | same |
| *Salmonella enterica* | ceftriaxone | 493 | 295 | 0.886 | 0.871 | 0.881 | 0.875 | same | same | same |
| *Salmonella enterica* | chloramphenicol | 720 | 486 | 0.938 | 0.973 | 0.986 | 0.924 | worse | better | same |
| *Salmonella enterica* | ciprofloxacin | 721 | 269 | 0.720 | 0.695 | 0.392 | 0.679 | better | worse | worse |
| *Salmonella enterica* | gentamicin | 629 | 241 | 0.746 | 0.809 | 0.755 | 0.721 | same | better | worse |
| *Salmonella enterica* | nalidixic acid | 641 | 485 | 0.814 | 0.880 | 0.951 | 0.838 | worse | better | better |
| *Salmonella enterica* | streptomycin | 427 | 137 | 0.893 | 0.815 | 0.445 | 0.861 | better | worse | worse |
| *Salmonella enterica* | sulfamethoxazole/trimethoprim | 728 | 58 | 0.776 | 0.552 | 0.828 | 0.764 | worse | worse | same |
| *Salmonella enterica* | sulfisoxazole | 402 | 67 | 0.706 | 0.622 | 0.970 | 0.751 | worse | worse | better |
| *Salmonella enterica* | tetracycline | 658 | 97 | 0.412 | 0.509 | 0.773 | 0.368 | worse | better | worse |
| *Staphylococcus aureus* | chloramphenicol | 166 | 115 | 1.000 | 1.000 | 1.000 | 1.000 | same | same | same |
| *Staphylococcus aureus* | clindamycin | 212 | 165 | 0.846 | 0.976 | 0.988 | 0.863 | worse | better | same |
| *Staphylococcus aureus* | daptomycin | 181 | 172 | 0.891 | 0.859 | 0.890 | 0.895 | same | worse | same |
| *Staphylococcus aureus* | erythromycin | 211 | 155 | 0.968 | 0.991 | 0.994 | 0.927 | worse | better | worse |
| *Streptococcus pneumoniae* | cefuroxime | 1621 | 66 | 0.934 | 0.859 | 0.561 | 0.929 | better | worse | same |
| *Streptococcus pneumoniae* | chloramphenicol | 1674 | 1469 | 0.995 | 0.983 | 0.982 | 0.991 | same | same | same |
| *Streptococcus pneumoniae* | clindamycin | 1770 | 152 | 0.943 | 0.967 | 0.954 | 0.919 | same | better | worse |
| *Streptococcus pneumoniae* | erythromycin | 1798 | 186 | 0.677 | 0.747 | 0.306 | 0.647 | better | better | worse |
| *Streptococcus pneumoniae* | meropenem | 1714 | 50 | 0.380 | 0.107 | 0.820 | 0.393 | worse | worse | same |
| *Streptococcus pneumoniae* | sulfamethoxazole/trimethoprim | 1632 | 227 | 0.824 | 0.900 | 0.410 | 0.918 | better | better | better |
| *Streptococcus pneumoniae* | tetracycline | 1805 | 132 | 0.902 | 0.929 | 0.419 | 0.917 | better | better | same |

**Study split.** 28 datasets.

| Species | Antibiotic | MIC genomes | Exact MICs | Released | Ordered | Median | Twin | P1 | P2 | P3 |
|---|---|---:|---:|---:|---:|---:|---:|---|---|---|
| *Acinetobacter baumannii* | cefepime | 162 | 53 | 0.497 | 0.509 | 0.170 | 0.453 | better | same | worse |
| *Escherichia coli* | amikacin | 472 | 87 | 0.609 | 0.521 | 0.766 | 0.667 | worse | worse | better |
| *Escherichia coli* | aztreonam | 184 | 54 | 0.346 | 0.321 | 0.833 | 0.457 | worse | worse | better |
| *Escherichia coli* | cefepime | 435 | 86 | 0.283 | 0.240 | 0.287 | 0.291 | same | worse | same |
| *Escherichia coli* | nitrofurantoin | 287 | 63 | 0.746 | 0.683 | 0.952 | 0.831 | worse | worse | better |
| *Escherichia coli* | piperacillin/tazobactam | 265 | 62 | 0.629 | 0.565 | 0.935 | 0.419 | worse | worse | worse |
| *Escherichia coli* | tetracycline | 310 | 67 | 0.398 | 0.662 | 0.199 | 0.517 | better | better | better |
| *Klebsiella pneumoniae* | aztreonam | 802 | 108 | 0.710 | 0.549 | 0.667 | 0.602 | better | worse | worse |
| *Klebsiella pneumoniae* | cefazolin | 742 | 54 | 0.333 | 0.142 | 0.296 | 0.340 | better | worse | same |
| *Klebsiella pneumoniae* | cefepime | 789 | 83 | 0.510 | 0.406 | 0.639 | 0.514 | worse | worse | same |
| *Klebsiella pneumoniae* | cefoxitin | 766 | 75 | 0.773 | 0.729 | 0.893 | 0.711 | worse | worse | worse |
| *Klebsiella pneumoniae* | ceftazidime | 844 | 178 | 0.519 | 0.511 | 0.646 | 0.552 | worse | same | better |
| *Klebsiella pneumoniae* | ciprofloxacin | 1003 | 78 | 0.662 | 0.748 | 0.833 | 0.675 | worse | better | same |
| *Klebsiella pneumoniae* | colistin | 275 | 95 | 0.267 | 0.344 | 0.182 | 0.274 | better | better | same |
| *Klebsiella pneumoniae* | ertapenem | 535 | 62 | 0.290 | 0.151 | 0.102 | 0.349 | better | worse | better |
| *Klebsiella pneumoniae* | imipenem | 903 | 174 | 0.452 | 0.718 | 0.621 | 0.517 | worse | better | better |
| *Klebsiella pneumoniae* | nitrofurantoin | 356 | 54 | 0.593 | 0.321 | 0.389 | 0.444 | better | worse | worse |
| *Klebsiella pneumoniae* | sulfamethoxazole/trimethoprim | 981 | 63 | 0.444 | 0.492 | 0.741 | 0.439 | worse | better | same |
| *Klebsiella pneumoniae* | tigecycline | 449 | 288 | 0.863 | 0.888 | 0.934 | 0.866 | worse | better | same |
| *Klebsiella pneumoniae* | tobramycin | 754 | 140 | 0.786 | 0.764 | 0.136 | 0.752 | better | worse | worse |
| *Pseudomonas aeruginosa* | aztreonam | 192 | 113 | 0.693 | 0.619 | 0.726 | 0.676 | worse | worse | same |
| *Pseudomonas aeruginosa* | cefepime | 203 | 136 | 0.627 | 0.578 | 0.654 | 0.642 | worse | worse | same |
| *Pseudomonas aeruginosa* | ceftazidime/avibactam | 282 | 274 | 0.602 | 0.656 | 0.380 | 0.605 | better | better | same |
| *Pseudomonas aeruginosa* | imipenem | 223 | 127 | 0.383 | 0.383 | 0.304 | 0.407 | better | same | better |
| *Pseudomonas aeruginosa* | piperacillin/tazobactam | 185 | 103 | 0.401 | 0.346 | 0.650 | 0.408 | worse | worse | same |
| *Salmonella enterica* | sulfisoxazole | 402 | 67 | 0.861 | 0.836 | 0.970 | 0.846 | worse | worse | same |
| *Staphylococcus aureus* | clindamycin | 212 | 165 | 0.937 | 0.988 | 0.988 | 0.988 | worse | better | better |
| *Staphylococcus aureus* | erythromycin | 211 | 155 | 0.935 | 0.994 | 0.996 | 0.974 | worse | better | better |
