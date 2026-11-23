# Mechanism split result

This page reads blocks M1 and M2 of the mechanism split once, by the rules `docs/mechanism-split.md` fixed before any run. M1 asks whether each learner calls resistant genomes whose catalogue calls are point mutations as often as those whose calls are acquired genes. M2 asks whether the mix of catalogue calls among deposited genomes moves between collection years. Neither block fits anything: both read the AMRFinderPlus program's tables beside the release's reports.

The rules of the reading are the method page's sections "M1, recall by mechanism", "M2, the mechanism mix over collection years" and "What a table can and cannot say", with the rulings on questions 1 to 14 in `docs/assumption-plan-corrections.md` (entries "Mechanism design" and "Mechanism split, first run"). No measure, subset or cut was added. The reading ran once on 23 Nov 2026, with no other option, over the runs of the 21 species at main `11c2996`.

Three callers stand apart on this page. The AMRFinderPlus program is the program itself, AMRFinderPlus 4.2.7 with database 2026-08-07.1. GRM's known-gene rule is GRM's own reader of the AMRFinderPlus catalogue, the catalog row of each report. The five learners are `scm`, `cart`, `rf`, `gbt` and `lasso`.

Every split here is a retrospective split inside one public database. A split by lineage or by study is never external validation. A share is no prevalence, and a trend describes deposits. The margin of 0.02 is an adapted research screen. Nothing here is a change comparison: no default, no released file and no sentence on the site changes.

## Reading

**M1, in plain words.** Recall is the share of resistant genomes a caller calls resistant. A resistant genome is gene only when the AMRFinderPlus program's rows that meet the antibiotic's rule are acquired genes alone, and point only when they are point mutations alone. On the lineage split of the primary set, 19 of the 204 datasets, in 5 species, hold 10 or more resistant genomes in both. On those 19, each learner's recall is higher on gene only, by more than 0.02, on 13 or 14 datasets, and higher on point only on 5 or 6, with a median difference of 0.206 to 0.492. With 5 species there is no species interval, since the method gives none under 6. The sign test after Holm over the five learners gives 0.318 for each, so each learner reads "no difference shown".

| Row | Datasets | Species | Median gene less point | Species interval | Gene ahead | Point ahead | Same | Sign p | Holm p | Reading |
|---|---|---|---|---|---|---|---|---|---|---|
| set covering machine (`scm`) | 19 | 5 | 0.492 | none | 14 | 5 | 0 | 0.064 | 0.318 | no difference shown |
| decision tree (`cart`) | 19 | 5 | 0.341 | none | 13 | 6 | 0 | 0.167 | 0.318 | no difference shown |
| random forest (`rf`) | 19 | 5 | 0.206 | none | 13 | 5 | 1 | 0.096 | 0.318 | no difference shown |
| boosted trees (`gbt`) | 19 | 5 | 0.338 | none | 14 | 5 | 0 | 0.064 | 0.318 | no difference shown |
| lasso (`lasso`) | 19 | 5 | 0.396 | none | 14 | 5 | 0 | 0.064 | 0.318 | no difference shown |
| GRM's known-gene rule | 19 | 5 | 0.984 | none | 17 | 0 | 2 | < 0.001 | beside | decides nothing |
| the AMRFinderPlus program | 19 | 5 | 0.000 | none | 0 | 0 | 19 | 1.000 | beside | decides nothing |

**GRM's known-gene rule** stands beside the learners and decides nothing. Its recall is higher on gene only on 17 of the 19 datasets and the same on 2, with a median difference of 0.984. Its recall on gene only is 1.000 on 17 of the 19 and 0.984 or more on all 19. On point only it is 0.000 on 10 of the 19.

**The AMRFinderPlus program** stands beside and decides nothing. Its recall is 1.000 on both strata in each of the 19 datasets, so its difference is 0.000 on 19 of 19. The strata are made from the program's own rows, so this reads how the strata are made, not how well the program calls.

**Seen after the reading, deciding nothing.** On every learner, the 5 quinolone datasets of the 19 are point ahead, and the gene-ahead datasets are all among the 14 others (section Seen after).

**M2, in plain words.** Per species and drug class, as the program's Class column names it, M2 compares genomes collected in 2010 to 2014 with those collected in 2015 to 2019: the share with a point mutation row of the class, and among genomes with an acquired gene row, the mix of gene families. A BioProject with 5 or more genomes in a period gets one vote, and a cell is read only with 3 or more voting BioProjects in each period. 110 of the 285 species and class cells pass the vote, all in 6 of the 21 species: *Acinetobacter baumannii*, *Enterococcus faecium*, *Escherichia coli*, *Klebsiella pneumoniae*, *Salmonella enterica* and *Staphylococcus aureus*. No cell moves. On the point share 0 of the 110 move, and the smallest Holm p is 0.330. On the family mix 0 of the 76 passing cells with a mix in both periods move, every Holm p 1.000. On 2000 to 2009 against 2010 to 2019, 89 of the 285 pass and none moves.

**Guards.** No guard moves a cell, on either pair of periods or either statistic (section Guards). The guards decide nothing.

Nothing in the release or on the site changes.

## Data

**Datasets.** A: the 237 datasets of the 21 species with the program's tables (`scripts/c3-inputs/lists.tsv`), 204 in the primary set and 33 in M, which wait for the maintainer's label rule. All 237 entered. The rule equals `drug_rules.tsv`'s on 237 of 237, the hard check passed on 237 of 237 with 2,844 of 2,844 learner and catalog rows equal to the reports' metrics, the program's call equals the comparison adapter's on 237 of 237, and the program's counts equal the AMRFinderPlus comparison's record on 237 of 237.

**M1.** On the lineage split of the primary set, 19 datasets in 5 species hold 10 or more resistant genomes in both gene only and point only: *Escherichia coli* with cefuroxime, cephalothin and ertapenem; *Klebsiella pneumoniae* with ampicillin/sulbactam, cefazolin, cefoxitin, ceftazidime, ceftriaxone, cefuroxime, ciprofloxacin, doripenem, ertapenem, imipenem and levofloxacin; *Pseudomonas aeruginosa* with amikacin and tobramycin; *Salmonella enterica* with nalidixic acid; *Shigella sonnei* with ciprofloxacin and nalidixic acid. In M, 5 datasets in 2 species do: *Escherichia coli* with amoxicillin/clavulanate, ampicillin and meropenem, and *Klebsiella pneumoniae* with meropenem and piperacillin/tazobactam.

**M2.** 21,625 genomes with a program table appear in the release's reports of the 21 species, and 6 genomes of *Staphylococcus aureus* in a report have none. Of the 21,625, 4,528 were collected in 2010 to 2014 and 6,563 in 2015 to 2019, 2,626 in 2000 to 2009, 359 before 2000 and 258 from 2020 on, and 7,291 have no year. A cell's genomes in a period are those of its voting BioProjects. There are 285 species and class cells, and 204 dataset cells, one per dataset of the primary set, for the secondary reading.

**Lineage clusters.** In each of the 21 species the reports agree on the lineage clusters, with 0 genome pairs differing, so the lineage guard is given in 21 of 21.

**Runs.** Each species ran `scripts/mechanism/run_species.sh` at main `11c2996` from the program's tables made before and the database 2026-08-07.1, with nothing installed. The reading is `combine.py` over the 21 runs' tables as they wrote them.

## M1

M1 reads each learner's recall on the resistant genomes of each stratum. A genome's stratum comes from the program's counted rows that meet the dataset's own rule, never from its label: point only (a `POINT` or `POINT_DISRUPT` row and no `AMR` row), gene only (an `AMR` row and no point row), both, and neither. A `POINT_DISRUPT` row is read with point (question 1).

The reading is, per learner on the lineage split of the primary set, recall on gene only less recall on point only, per dataset where both strata hold 10 or more resistant genomes. The median is given with a species interval from 4,000 resamples of species, read at the 101st and 3,900th values. A dataset is gene ahead or point ahead by more than 0.02. The two-sided sign test of gene ahead against point ahead runs with Holm over the five learners. A learner reads "gene only higher" or "point only higher" at a Holm p under 0.05, and "no difference shown" otherwise.

GRM's known-gene rule and the AMRFinderPlus program stand beside the learners and decide nothing. The random split stands beside, labelled "training and test share lineages", and the datasets in M stand apart. M1 gives recall alone, as the method fixes it (question 11).

**Beside: M and the random split.** On the 5 datasets in M, on the lineage split, each learner is gene ahead on 5 of 5, with a median of 0.269 to 0.630 and a sign p of 0.062. M stands apart and enters no Holm family. On the random split, where training and test share lineages, each learner is gene ahead on 13 or 14 of the 19 primary datasets, with a median of 0.276 to 0.452 and a sign p of 0.031 to 0.167. It stands beside, with no Holm, and decides nothing.

| Split | Set | Row | Datasets | Species | Median gene less point | Species interval | Gene ahead | Point ahead | Same | Sign p |
|---|---|---|---|---|---|---|---|---|---|---|
| lineage | M | set covering machine (`scm`) | 5 | 2 | 0.620 | none | 5 | 0 | 0 | 0.062 |
| lineage | M | decision tree (`cart`) | 5 | 2 | 0.630 | none | 5 | 0 | 0 | 0.062 |
| lineage | M | random forest (`rf`) | 5 | 2 | 0.269 | none | 5 | 0 | 0 | 0.062 |
| lineage | M | boosted trees (`gbt`) | 5 | 2 | 0.583 | none | 5 | 0 | 0 | 0.062 |
| lineage | M | lasso (`lasso`) | 5 | 2 | 0.570 | none | 5 | 0 | 0 | 0.062 |
| lineage | M | GRM's known-gene rule | 5 | 2 | 0.640 | none | 4 | 0 | 1 | 0.125 |
| lineage | M | the AMRFinderPlus program | 5 | 2 | 0.000 | none | 0 | 0 | 5 | 1.000 |
| random (training and test share lineages) | primary | set covering machine (`scm`) | 19 | 5 | 0.452 | none | 14 | 4 | 1 | 0.031 |
| random (training and test share lineages) | primary | decision tree (`cart`) | 19 | 5 | 0.314 | none | 14 | 4 | 1 | 0.031 |
| random (training and test share lineages) | primary | random forest (`rf`) | 19 | 5 | 0.288 | none | 14 | 5 | 0 | 0.064 |
| random (training and test share lineages) | primary | boosted trees (`gbt`) | 19 | 5 | 0.280 | none | 13 | 5 | 1 | 0.096 |
| random (training and test share lineages) | primary | lasso (`lasso`) | 19 | 5 | 0.276 | none | 13 | 6 | 0 | 0.167 |
| random (training and test share lineages) | primary | GRM's known-gene rule | 19 | 5 | 0.984 | none | 17 | 0 | 2 | < 0.001 |
| random (training and test share lineages) | primary | the AMRFinderPlus program | 19 | 5 | 0.000 | none | 0 | 0 | 19 | 1.000 |
| random (training and test share lineages) | M | set covering machine (`scm`) | 5 | 2 | 0.726 | none | 5 | 0 | 0 | 0.062 |
| random (training and test share lineages) | M | decision tree (`cart`) | 5 | 2 | 0.319 | none | 4 | 1 | 0 | 0.375 |
| random (training and test share lineages) | M | random forest (`rf`) | 5 | 2 | 0.407 | none | 5 | 0 | 0 | 0.062 |
| random (training and test share lineages) | M | boosted trees (`gbt`) | 5 | 2 | 0.500 | none | 5 | 0 | 0 | 0.062 |
| random (training and test share lineages) | M | lasso (`lasso`) | 5 | 2 | 0.745 | none | 5 | 0 | 0 | 0.062 |
| random (training and test share lineages) | M | GRM's known-gene rule | 5 | 2 | 0.640 | none | 4 | 0 | 1 | 0.125 |
| random (training and test share lineages) | M | the AMRFinderPlus program | 5 | 2 | 0.000 | none | 0 | 0 | 5 | 1.000 |

**Neither.** Over the 92 primary datasets with 10 or more resistant genomes in neither stratum, resistant genomes with no program row that meets the rule, the learners' median recall on the lineage split is 0.167 to 0.253. The median recall of GRM's known-gene rule is 0.000, and that of the AMRFinderPlus program 0.000.

| Split | Set | Row | Datasets | Median recall on neither |
|---|---|---|---|---|
| lineage | primary | set covering machine (`scm`) | 92 | 0.193 |
| lineage | primary | decision tree (`cart`) | 92 | 0.230 |
| lineage | primary | random forest (`rf`) | 92 | 0.167 |
| lineage | primary | boosted trees (`gbt`) | 92 | 0.253 |
| lineage | primary | lasso (`lasso`) | 92 | 0.218 |
| lineage | primary | GRM's known-gene rule | 92 | 0.000 |
| lineage | primary | the AMRFinderPlus program | 92 | 0.000 |
| lineage | M | set covering machine (`scm`) | 19 | 0.069 |
| lineage | M | decision tree (`cart`) | 19 | 0.189 |
| lineage | M | random forest (`rf`) | 19 | 0.150 |
| lineage | M | boosted trees (`gbt`) | 19 | 0.125 |
| lineage | M | lasso (`lasso`) | 19 | 0.097 |
| lineage | M | GRM's known-gene rule | 19 | 0.000 |
| lineage | M | the AMRFinderPlus program | 19 | 0.000 |

**Strata.** Genome counts per stratum on the lineage split, summed over the entered datasets, a genome once per dataset. Of the 62,878 resistant genomes of the primary set, 39,814 are gene only, 7,963 point only, 9,066 both and 6,035 neither. Point only holds 10 or more resistant genomes in 47 of the 204 datasets and gene only in 156. The page writes no recall pooled over datasets.

| Set | Stratum | Datasets | With a resistant genome | With 10 or more resistant | Resistant genomes | Susceptible genomes |
|---|---|---|---|---|---|---|
| primary | point only | 204 | 78 | 47 | 7,963 | 6,791 |
| primary | gene only | 204 | 171 | 156 | 39,814 | 11,041 |
| primary | both | 204 | 85 | 68 | 9,066 | 1,479 |
| primary | neither | 204 | 150 | 92 | 6,035 | 64,480 |
| primary | partial only:point only | 204 | 0 | 0 | 0 | 0 |
| primary | partial only:gene only | 204 | 59 | 8 | 329 | 243 |
| primary | partial only:both | 204 | 0 | 0 | 0 | 0 |
| primary | with point disrupt:point only | 204 | 38 | 11 | 637 | 657 |
| primary | with point disrupt:gene only | 204 | 0 | 0 | 0 | 0 |
| primary | with point disrupt:both | 204 | 51 | 43 | 4,967 | 598 |
| primary | point disrupt alone:point only | 204 | 32 | 5 | 177 | 433 |
| primary | point disrupt alone:both | 204 | 45 | 33 | 2,877 | 267 |
| primary | all genomes | 204 | 204 | 204 | 62,878 | 83,791 |
| M | point only | 33 | 15 | 12 | 3,513 | 2,248 |
| M | gene only | 33 | 24 | 21 | 11,383 | 3,208 |
| M | both | 33 | 18 | 14 | 2,707 | 205 |
| M | neither | 33 | 28 | 19 | 1,002 | 17,476 |
| M | partial only:point only | 33 | 0 | 0 | 0 | 0 |
| M | partial only:gene only | 33 | 9 | 3 | 115 | 20 |
| M | partial only:both | 33 | 0 | 0 | 0 | 0 |
| M | with point disrupt:point only | 33 | 7 | 3 | 108 | 220 |
| M | with point disrupt:gene only | 33 | 0 | 0 | 0 | 0 |
| M | with point disrupt:both | 33 | 7 | 5 | 1,519 | 69 |
| M | point disrupt alone:point only | 33 | 7 | 2 | 39 | 171 |
| M | point disrupt alone:both | 33 | 6 | 5 | 996 | 65 |
| M | all genomes | 33 | 33 | 33 | 18,605 | 23,137 |

## M2

A cell is a species with a drug class as the program's Class column names it. Per period, the cell gives the share of genomes with a point row of the class, the share with a gene row (Subtype `AMR`), and among genomes with a gene row the mix of catalogue families over their Subtype `AMR` rows alone. The primary periods are 2010 to 2014 against 2015 to 2019, and the secondary periods 2000 to 2009 against 2010 to 2019. A year is four digits in `samples[].year`, and nothing is imputed.

One vote per BioProject: a BioProject votes in a period with 5 or more genomes in it, and a cell passes with 3 or more voting BioProjects in each period. A period's share is the mean over its voting BioProjects of each one's share. A permutation of period labels over BioProjects, 1,000 times, gives the p of the point share's difference and of the family mix's total variation distance. Holm runs over the passing cells, each statistic apart, and a cell moves at a Holm p under 0.05. A cell that does not pass writes no share, mix or p. A family is named only where, in each period, 5 or more genomes of the voting BioProjects carry it or none do; the rest share one row, `other families`.

Only the BioProject vote on the primary periods over species and class cells decides. Everything else in this section and the next stands beside.

**The deciding cells.** On the primary periods 110 of the 285 cells pass, in 6 species: *Acinetobacter baumannii* 15 of 15, *Enterococcus faecium* 17 of 17, *Escherichia coli* 21 of 21, *Klebsiella pneumoniae* 19 of 19, *Salmonella enterica* 17 of 17 and *Staphylococcus aureus* 21 of 21. No cell of the other 15 species passes. On the point share 0 of the 110 move. 76 of the 110 have a family mix in both periods, and 0 of the 76 move. The table gives the 110, with the voting BioProjects and their genomes per period, earlier then later. The family mix of each passing cell, by the database's family names over 5 or more genomes, is `page-m2-families.tsv`.

| Species | Class | BioProjects voting | Genomes | Point share | Holm p | Gene share | Mix distance | Holm p | Moves |
|---|---|---|---|---|---|---|---|---|---|
| *Acinetobacter baumannii* | AMINOGLYCOSIDE | 7 / 5 | 152 / 169 | 0.000 to 0.000 | 1.000 | 0.990 to 1.000 | 0.261 | 1.000 | no |
| *Acinetobacter baumannii* | BETA-LACTAM | 7 / 5 | 152 / 169 | 0.018 to 0.379 | 0.330 | 1.000 to 1.000 | 0.046 | 1.000 | no |
| *Acinetobacter baumannii* | BLEOMYCIN | 7 / 5 | 152 / 169 | 0.000 to 0.000 | 1.000 | 0.000 to 0.008 | none | none | no |
| *Acinetobacter baumannii* | COLISTIN | 7 / 5 | 152 / 169 | 0.029 to 0.055 | 1.000 | 0.000 to 0.000 | none | none | no |
| *Acinetobacter baumannii* | EFFLUX | 7 / 5 | 152 / 169 | 0.012 to 0.025 | 1.000 | 0.000 to 0.000 | none | none | no |
| *Acinetobacter baumannii* | FOSFOMYCIN | 7 / 5 | 152 / 169 | 0.000 to 0.000 | 1.000 | 0.000 to 0.000 | none | none | no |
| *Acinetobacter baumannii* | MACROLIDE | 7 / 5 | 152 / 169 | 0.000 to 0.000 | 1.000 | 0.536 to 0.489 | 0.000 | 1.000 | no |
| *Acinetobacter baumannii* | PHENICOL | 7 / 5 | 152 / 169 | 0.000 to 0.000 | 1.000 | 0.296 to 0.091 | 0.264 | 1.000 | no |
| *Acinetobacter baumannii* | QUINOLONE | 7 / 5 | 152 / 169 | 0.856 to 0.985 | 1.000 | 0.000 to 0.000 | none | none | no |
| *Acinetobacter baumannii* | RIFAMYCIN | 7 / 5 | 152 / 169 | 0.000 to 0.005 | 1.000 | 0.043 to 0.010 | 0.000 | 1.000 | no |
| *Acinetobacter baumannii* | STREPTOGRAMIN | 7 / 5 | 152 / 169 | 0.000 to 0.000 | 1.000 | 0.536 to 0.489 | 0.000 | 1.000 | no |
| *Acinetobacter baumannii* | STREPTOTHRICIN | 7 / 5 | 152 / 169 | 0.000 to 0.000 | 1.000 | 0.000 to 0.020 | none | none | no |
| *Acinetobacter baumannii* | SULFONAMIDE | 7 / 5 | 152 / 169 | 0.000 to 0.000 | 1.000 | 0.847 to 0.756 | 0.334 | 1.000 | no |
| *Acinetobacter baumannii* | TETRACYCLINE | 7 / 5 | 152 / 169 | 0.000 to 0.000 | 1.000 | 0.504 to 0.683 | 0.260 | 1.000 | no |
| *Acinetobacter baumannii* | TRIMETHOPRIM | 7 / 5 | 152 / 169 | 0.000 to 0.000 | 1.000 | 0.214 to 0.065 | 0.917 | 1.000 | no |
| *Enterococcus faecium* | AMINOGLYCOSIDE | 4 / 9 | 237 / 545 | 0.000 to 0.000 | 1.000 | 0.995 to 1.000 | 0.076 | 1.000 | no |
| *Enterococcus faecium* | AVILAMYCIN | 4 / 9 | 237 / 545 | 0.000 to 0.000 | 1.000 | 0.000 to 0.000 | none | none | no |
| *Enterococcus faecium* | BETA-LACTAM | 4 / 9 | 237 / 545 | 0.965 to 0.936 | 1.000 | 0.000 to 0.000 | none | none | no |
| *Enterococcus faecium* | FOSFOMYCIN | 4 / 9 | 237 / 545 | 0.000 to 0.000 | 1.000 | 0.000 to 0.000 | none | none | no |
| *Enterococcus faecium* | GLYCOPEPTIDE | 4 / 9 | 237 / 545 | 0.000 to 0.000 | 1.000 | 0.710 to 0.776 | 0.060 | 1.000 | no |
| *Enterococcus faecium* | LINCOSAMIDE | 4 / 9 | 237 / 545 | 0.000 to 0.000 | 1.000 | 0.750 to 0.827 | 0.199 | 1.000 | no |
| *Enterococcus faecium* | LIPOPEPTIDE | 4 / 9 | 237 / 545 | 0.769 to 0.636 | 1.000 | 0.000 to 0.000 | none | none | no |
| *Enterococcus faecium* | MACROLIDE | 4 / 9 | 237 / 545 | 0.000 to 0.000 | 1.000 | 1.000 to 1.000 | 0.091 | 1.000 | no |
| *Enterococcus faecium* | OXAZOLIDINONE | 4 / 9 | 237 / 545 | 0.065 to 0.075 | 1.000 | 0.008 to 0.116 | 0.626 | 1.000 | no |
| *Enterococcus faecium* | PHENICOL | 4 / 9 | 237 / 545 | 0.000 to 0.000 | 1.000 | 0.015 to 0.196 | 0.563 | 1.000 | no |
| *Enterococcus faecium* | PLEUROMUTILIN | 4 / 9 | 237 / 545 | 0.287 to 0.273 | 1.000 | 0.005 to 0.105 | 0.667 | 1.000 | no |
| *Enterococcus faecium* | QUINOLONE | 4 / 9 | 237 / 545 | 0.912 to 0.887 | 1.000 | 0.000 to 0.000 | none | none | no |
| *Enterococcus faecium* | RIFAMYCIN | 4 / 9 | 237 / 545 | 0.078 to 0.044 | 1.000 | 0.000 to 0.000 | none | none | no |
| *Enterococcus faecium* | STREPTOGRAMIN | 4 / 9 | 237 / 545 | 0.000 to 0.000 | 1.000 | 1.000 to 1.000 | 0.103 | 1.000 | no |
| *Enterococcus faecium* | STREPTOTHRICIN | 4 / 9 | 237 / 545 | 0.000 to 0.000 | 1.000 | 0.406 to 0.284 | 0.000 | 1.000 | no |
| *Enterococcus faecium* | TETRACYCLINE | 4 / 9 | 237 / 545 | 0.125 to 0.027 | 1.000 | 0.764 to 0.848 | 0.100 | 1.000 | no |
| *Enterococcus faecium* | TRIMETHOPRIM | 4 / 9 | 237 / 545 | 0.000 to 0.000 | 1.000 | 0.681 to 0.712 | 0.058 | 1.000 | no |
| *Escherichia coli* | AMINOGLYCOSIDE | 8 / 6 | 351 / 375 | 0.000 to 0.000 | 1.000 | 0.635 to 0.697 | 0.178 | 1.000 | no |
| *Escherichia coli* | BETA-LACTAM | 8 / 6 | 351 / 375 | 0.086 to 0.139 | 1.000 | 0.710 to 0.769 | 0.152 | 1.000 | no |
| *Escherichia coli* | BLEOMYCIN | 8 / 6 | 351 / 375 | 0.000 to 0.000 | 1.000 | 0.000 to 0.006 | none | none | no |
| *Escherichia coli* | COLISTIN | 8 / 6 | 351 / 375 | 0.020 to 0.003 | 1.000 | 0.102 to 0.024 | 0.000 | 1.000 | no |
| *Escherichia coli* | FOSFOMYCIN | 8 / 6 | 351 / 375 | 0.435 to 0.370 | 1.000 | 0.008 to 0.000 | none | none | no |
| *Escherichia coli* | FOSMIDOMYCIN | 8 / 6 | 351 / 375 | 0.062 to 0.127 | 1.000 | 0.000 to 0.000 | none | none | no |
| *Escherichia coli* | GLYCOPEPTIDE | 8 / 6 | 351 / 375 | 0.000 to 0.000 | 1.000 | 0.005 to 0.000 | none | none | no |
| *Escherichia coli* | ISONIAZID | 8 / 6 | 351 / 375 | 0.000 to 0.000 | 1.000 | 0.000 to 0.000 | none | none | no |
| *Escherichia coli* | LINCOSAMIDE | 8 / 6 | 351 / 375 | 0.000 to 0.000 | 1.000 | 0.002 to 0.035 | 0.750 | 1.000 | no |
| *Escherichia coli* | MACROLIDE | 8 / 6 | 351 / 375 | 0.000 to 0.000 | 1.000 | 0.342 to 0.323 | 0.231 | 1.000 | no |
| *Escherichia coli* | MULTIDRUG | 8 / 6 | 351 / 375 | 0.096 to 0.125 | 1.000 | 0.000 to 0.000 | none | none | no |
| *Escherichia coli* | NITROFURAN | 8 / 6 | 351 / 375 | 0.063 to 0.190 | 1.000 | 0.000 to 0.000 | none | none | no |
| *Escherichia coli* | PHENICOL | 8 / 6 | 351 / 375 | 0.000 to 0.000 | 1.000 | 0.413 to 0.394 | 0.101 | 1.000 | no |
| *Escherichia coli* | QUINOLONE | 8 / 6 | 351 / 375 | 0.535 to 0.736 | 1.000 | 0.237 to 0.196 | 0.325 | 1.000 | no |
| *Escherichia coli* | RIFAMYCIN | 8 / 6 | 351 / 375 | 0.000 to 0.000 | 1.000 | 0.005 to 0.003 | 0.667 | 1.000 | no |
| *Escherichia coli* | STREPTOGRAMIN | 8 / 6 | 351 / 375 | 0.000 to 0.000 | 1.000 | 0.001 to 0.008 | 0.000 | 1.000 | no |
| *Escherichia coli* | STREPTOTHRICIN | 8 / 6 | 351 / 375 | 0.000 to 0.000 | 1.000 | 0.042 to 0.060 | 0.250 | 1.000 | no |
| *Escherichia coli* | SULFONAMIDE | 8 / 6 | 351 / 375 | 0.000 to 0.000 | 1.000 | 0.622 to 0.580 | 0.162 | 1.000 | no |
| *Escherichia coli* | TETRACYCLINE | 8 / 6 | 351 / 375 | 0.000 to 0.000 | 1.000 | 0.534 to 0.575 | 0.162 | 1.000 | no |
| *Escherichia coli* | TRICLOSAN | 8 / 6 | 351 / 375 | 0.000 to 0.007 | 1.000 | 0.000 to 0.000 | none | none | no |
| *Escherichia coli* | TRIMETHOPRIM | 8 / 6 | 351 / 375 | 0.000 to 0.000 | 1.000 | 0.513 to 0.575 | 0.330 | 1.000 | no |
| *Klebsiella pneumoniae* | AMINOGLYCOSIDE | 8 / 5 | 615 / 156 | 0.000 to 0.000 | 1.000 | 0.889 to 0.910 | 0.151 | 1.000 | no |
| *Klebsiella pneumoniae* | BETA-LACTAM | 8 / 5 | 615 / 156 | 0.552 to 0.368 | 1.000 | 0.998 to 0.993 | 0.175 | 1.000 | no |
| *Klebsiella pneumoniae* | BLEOMYCIN | 8 / 5 | 615 / 156 | 0.000 to 0.000 | 1.000 | 0.132 to 0.105 | 0.667 | 1.000 | no |
| *Klebsiella pneumoniae* | COLISTIN | 8 / 5 | 615 / 156 | 0.028 to 0.021 | 1.000 | 0.000 to 0.007 | none | none | no |
| *Klebsiella pneumoniae* | FOSFOMYCIN | 8 / 5 | 615 / 156 | 0.000 to 0.000 | 1.000 | 1.000 to 1.000 | 0.017 | 1.000 | no |
| *Klebsiella pneumoniae* | LINCOSAMIDE | 8 / 5 | 615 / 156 | 0.000 to 0.000 | 1.000 | 0.114 to 0.129 | 0.392 | 1.000 | no |
| *Klebsiella pneumoniae* | MACROLIDE | 8 / 5 | 615 / 156 | 0.000 to 0.000 | 1.000 | 0.571 to 0.373 | 0.267 | 1.000 | no |
| *Klebsiella pneumoniae* | MULTIDRUG | 8 / 5 | 615 / 156 | 0.043 to 0.010 | 1.000 | 0.000 to 0.000 | none | none | no |
| *Klebsiella pneumoniae* | NITROFURAN | 8 / 5 | 615 / 156 | 0.108 to 0.141 | 1.000 | 0.889 to 0.837 | 0.088 | 1.000 | no |
| *Klebsiella pneumoniae* | OXAZOLIDINONE | 8 / 5 | 615 / 156 | 0.000 to 0.000 | 1.000 | 0.000 to 0.014 | none | none | no |
| *Klebsiella pneumoniae* | PHENICOL | 8 / 5 | 615 / 156 | 0.013 to 0.013 | 1.000 | 0.961 to 0.991 | 0.127 | 1.000 | no |
| *Klebsiella pneumoniae* | PLEUROMUTILIN | 8 / 5 | 615 / 156 | 0.000 to 0.000 | 1.000 | 0.000 to 0.014 | none | none | no |
| *Klebsiella pneumoniae* | QUINOLONE | 8 / 5 | 615 / 156 | 0.648 to 0.411 | 1.000 | 0.973 to 1.000 | 0.095 | 1.000 | no |
| *Klebsiella pneumoniae* | RIFAMYCIN | 8 / 5 | 615 / 156 | 0.000 to 0.000 | 1.000 | 0.187 to 0.393 | 0.231 | 1.000 | no |
| *Klebsiella pneumoniae* | STREPTOGRAMIN | 8 / 5 | 615 / 156 | 0.000 to 0.000 | 1.000 | 0.186 to 0.115 | 0.243 | 1.000 | no |
| *Klebsiella pneumoniae* | STREPTOTHRICIN | 8 / 5 | 615 / 156 | 0.000 to 0.000 | 1.000 | 0.030 to 0.006 | 0.000 | 1.000 | no |
| *Klebsiella pneumoniae* | SULFONAMIDE | 8 / 5 | 615 / 156 | 0.000 to 0.000 | 1.000 | 0.833 to 0.833 | 0.032 | 1.000 | no |
| *Klebsiella pneumoniae* | TETRACYCLINE | 8 / 5 | 615 / 156 | 0.102 to 0.102 | 1.000 | 0.968 to 0.943 | 0.127 | 1.000 | no |
| *Klebsiella pneumoniae* | TRIMETHOPRIM | 8 / 5 | 615 / 156 | 0.000 to 0.000 | 1.000 | 0.814 to 0.808 | 0.255 | 1.000 | no |
| *Salmonella enterica* | AMINOGLYCOSIDE | 11 / 6 | 550 / 453 | 0.000 to 0.000 | 1.000 | 0.736 to 0.743 | 0.230 | 1.000 | no |
| *Salmonella enterica* | BETA-LACTAM | 11 / 6 | 550 / 453 | 0.000 to 0.000 | 1.000 | 0.734 to 0.628 | 0.384 | 1.000 | no |
| *Salmonella enterica* | BLEOMYCIN | 11 / 6 | 550 / 453 | 0.000 to 0.000 | 1.000 | 0.091 to 0.142 | 0.232 | 1.000 | no |
| *Salmonella enterica* | COLISTIN | 11 / 6 | 550 / 453 | 0.000 to 0.000 | 1.000 | 0.009 to 0.013 | 0.500 | 1.000 | no |
| *Salmonella enterica* | FOSFOMYCIN | 11 / 6 | 550 / 453 | 0.000 to 0.000 | 1.000 | 0.121 to 0.150 | 0.395 | 1.000 | no |
| *Salmonella enterica* | LINCOSAMIDE | 11 / 6 | 550 / 453 | 0.000 to 0.000 | 1.000 | 0.012 to 0.022 | 0.782 | 1.000 | no |
| *Salmonella enterica* | MACROLIDE | 11 / 6 | 550 / 453 | 0.003 to 0.001 | 1.000 | 0.142 to 0.167 | 0.292 | 1.000 | no |
| *Salmonella enterica* | MULTIDRUG | 11 / 6 | 550 / 453 | 0.023 to 0.072 | 1.000 | 0.000 to 0.000 | none | none | no |
| *Salmonella enterica* | NITROFURAN | 11 / 6 | 550 / 453 | 0.000 to 0.000 | 1.000 | 0.038 to 0.099 | 0.011 | 1.000 | no |
| *Salmonella enterica* | PHENICOL | 11 / 6 | 550 / 453 | 0.000 to 0.000 | 1.000 | 0.464 to 0.473 | 0.304 | 1.000 | no |
| *Salmonella enterica* | QUINOLONE | 11 / 6 | 550 / 453 | 0.320 to 0.603 | 1.000 | 0.120 to 0.253 | 0.311 | 1.000 | no |
| *Salmonella enterica* | RIFAMYCIN | 11 / 6 | 550 / 453 | 0.000 to 0.000 | 1.000 | 0.098 to 0.170 | 0.231 | 1.000 | no |
| *Salmonella enterica* | STREPTOGRAMIN | 11 / 6 | 550 / 453 | 0.000 to 0.000 | 1.000 | 0.003 to 0.030 | 0.583 | 1.000 | no |
| *Salmonella enterica* | SULFONAMIDE | 11 / 6 | 550 / 453 | 0.000 to 0.000 | 1.000 | 0.725 to 0.722 | 0.038 | 1.000 | no |
| *Salmonella enterica* | TETRACYCLINE | 11 / 6 | 550 / 453 | 0.000 to 0.000 | 1.000 | 0.648 to 0.539 | 0.147 | 1.000 | no |
| *Salmonella enterica* | TRICLOSAN | 11 / 6 | 550 / 453 | 0.114 to 0.211 | 1.000 | 0.000 to 0.000 | none | none | no |
| *Salmonella enterica* | TRIMETHOPRIM | 11 / 6 | 550 / 453 | 0.000 to 0.000 | 1.000 | 0.242 to 0.479 | 0.400 | 1.000 | no |
| *Staphylococcus aureus* | AMINOCOUMARIN | 4 / 6 | 365 / 815 | 0.161 to 0.037 | 1.000 | 0.000 to 0.000 | none | none | no |
| *Staphylococcus aureus* | AMINOGLYCOSIDE | 4 / 6 | 365 / 815 | 0.000 to 0.000 | 1.000 | 0.714 to 0.292 | 0.151 | 1.000 | no |
| *Staphylococcus aureus* | BETA-LACTAM | 4 / 6 | 365 / 815 | 0.000 to 0.039 | 1.000 | 0.968 to 0.772 | 0.235 | 1.000 | no |
| *Staphylococcus aureus* | BLEOMYCIN | 4 / 6 | 365 / 815 | 0.000 to 0.000 | 1.000 | 0.035 to 0.005 | 0.000 | 1.000 | no |
| *Staphylococcus aureus* | FOSFOMYCIN | 4 / 6 | 365 / 815 | 0.479 to 0.315 | 1.000 | 0.717 to 0.688 | 0.032 | 1.000 | no |
| *Staphylococcus aureus* | FUSIDIC ACID | 4 / 6 | 365 / 815 | 0.000 to 0.000 | 1.000 | 0.013 to 0.013 | 0.000 | 1.000 | no |
| *Staphylococcus aureus* | FUSIDIC_ACID | 4 / 6 | 365 / 815 | 0.010 to 0.001 | 1.000 | 0.000 to 0.000 | none | none | no |
| *Staphylococcus aureus* | LINCOSAMIDE | 4 / 6 | 365 / 815 | 0.000 to 0.000 | 1.000 | 0.640 to 0.356 | 0.382 | 1.000 | no |
| *Staphylococcus aureus* | LIPOPEPTIDE | 4 / 6 | 365 / 815 | 0.056 to 0.062 | 1.000 | 0.000 to 0.000 | none | none | no |
| *Staphylococcus aureus* | MACROLIDE | 4 / 6 | 365 / 815 | 0.000 to 0.000 | 1.000 | 0.638 to 0.330 | 0.397 | 1.000 | no |
| *Staphylococcus aureus* | MUPIROCIN | 4 / 6 | 365 / 815 | 0.056 to 0.000 | 1.000 | 0.085 to 0.001 | 0.000 | 1.000 | no |
| *Staphylococcus aureus* | OXAZOLIDINONE | 4 / 6 | 365 / 815 | 0.000 to 0.000 | 1.000 | 0.000 to 0.000 | none | none | no |
| *Staphylococcus aureus* | PHENICOL | 4 / 6 | 365 / 815 | 0.000 to 0.000 | 1.000 | 0.301 to 0.241 | 0.064 | 1.000 | no |
| *Staphylococcus aureus* | PLEUROMUTILIN | 4 / 6 | 365 / 815 | 0.000 to 0.000 | 1.000 | 0.000 to 0.000 | none | none | no |
| *Staphylococcus aureus* | QUINOLONE | 4 / 6 | 365 / 815 | 0.695 to 0.320 | 1.000 | 0.000 to 0.000 | none | none | no |
| *Staphylococcus aureus* | RIFAMYCIN | 4 / 6 | 365 / 815 | 0.029 to 0.060 | 1.000 | 0.000 to 0.000 | none | none | no |
| *Staphylococcus aureus* | STREPTOGRAMIN | 4 / 6 | 365 / 815 | 0.000 to 0.000 | 1.000 | 0.386 to 0.325 | 0.246 | 1.000 | no |
| *Staphylococcus aureus* | STREPTOTHRICIN | 4 / 6 | 365 / 815 | 0.000 to 0.000 | 1.000 | 0.262 to 0.015 | 0.000 | 1.000 | no |
| *Staphylococcus aureus* | SULFONAMIDE | 4 / 6 | 365 / 815 | 0.415 to 0.046 | 1.000 | 0.000 to 0.000 | none | none | no |
| *Staphylococcus aureus* | TETRACYCLINE | 4 / 6 | 365 / 815 | 0.090 to 0.172 | 1.000 | 0.645 to 0.399 | 0.283 | 1.000 | no |
| *Staphylococcus aureus* | TRIMETHOPRIM | 4 / 6 | 365 / 815 | 0.085 to 0.047 | 1.000 | 0.387 to 0.275 | 0.162 | 1.000 | no |

**Secondary periods.** On 2000 to 2009 against 2010 to 2019, 89 of the 285 cells pass, in 5 species: *Acinetobacter baumannii*, *Escherichia coli*, *Pseudomonas aeruginosa*, *Salmonella enterica* and *Staphylococcus aureus*. 59 of the 89 have a family mix in both periods, and none moves on either statistic.

## Guards

The guards stand beside the BioProject vote and decide nothing: the share per genome with no vote, the largest BioProject left out, flagged BioProjects left out (0.95 or more resistant labels for the class, or one family in 0.95 or more of their gene genomes, read over every year), one vote per lineage cluster where the species' reports agree on the clusters, the clusters voting in both periods, countries present in both periods only, and rows of no partial method only. Resistant genomes per dataset are a secondary reading. The time-split link is exploratory.

**Species and class cells.** Under every guard, on both pairs of periods and both statistics, 0 cells move. The guards change how many cells pass. On the primary periods, among the guards with a vote, that runs from 34 cells with a mix, with flagged BioProjects left out, to 145 cells on the point share, with one vote per lineage cluster. The share per genome has no vote and no p, so by its own rule it passes and moves no cell. Its shares are written only where the cell's BioProject vote passes.

| Periods | Guard | Statistic | Cells | Pass | Move | Decides |
|---|---|---|---|---|---|---|
| 2010 to 2014 against 2015 to 2019 | one vote per BioProject (decides) | point | 285 | 110 | 0 | yes |
| 2010 to 2014 against 2015 to 2019 | one vote per BioProject (decides) | tvd | 285 | 76 | 0 | yes |
| 2010 to 2014 against 2015 to 2019 | countries in both periods only | point | 285 | 110 | 0 | no |
| 2010 to 2014 against 2015 to 2019 | countries in both periods only | tvd | 285 | 73 | 0 | no |
| 2010 to 2014 against 2015 to 2019 | flagged BioProjects left out | point | 285 | 76 | 0 | no |
| 2010 to 2014 against 2015 to 2019 | flagged BioProjects left out | tvd | 285 | 34 | 0 | no |
| 2010 to 2014 against 2015 to 2019 | rows of no partial method only | point | 285 | 110 | 0 | no |
| 2010 to 2014 against 2015 to 2019 | rows of no partial method only | tvd | 285 | 75 | 0 | no |
| 2010 to 2014 against 2015 to 2019 | largest BioProject left out | point | 285 | 110 | 0 | no |
| 2010 to 2014 against 2015 to 2019 | largest BioProject left out | tvd | 285 | 74 | 0 | no |
| 2010 to 2014 against 2015 to 2019 | one vote per lineage cluster | point | 285 | 145 | 0 | no |
| 2010 to 2014 against 2015 to 2019 | one vote per lineage cluster | tvd | 285 | 86 | 0 | no |
| 2010 to 2014 against 2015 to 2019 | lineage clusters voting in both periods | point | 285 | 104 | 0 | no |
| 2010 to 2014 against 2015 to 2019 | lineage clusters voting in both periods | tvd | 285 | 68 | 0 | no |
| 2010 to 2014 against 2015 to 2019 | share per genome, no vote | point | 285 | 0 | 0 | no |
| 2010 to 2014 against 2015 to 2019 | share per genome, no vote | tvd | 285 | 0 | 0 | no |
| 2000 to 2009 against 2010 to 2019 | one vote per BioProject (decides) | point | 285 | 89 | 0 | no |
| 2000 to 2009 against 2010 to 2019 | one vote per BioProject (decides) | tvd | 285 | 59 | 0 | no |
| 2000 to 2009 against 2010 to 2019 | countries in both periods only | point | 285 | 89 | 0 | no |
| 2000 to 2009 against 2010 to 2019 | countries in both periods only | tvd | 285 | 55 | 0 | no |
| 2000 to 2009 against 2010 to 2019 | flagged BioProjects left out | point | 285 | 50 | 0 | no |
| 2000 to 2009 against 2010 to 2019 | flagged BioProjects left out | tvd | 285 | 19 | 0 | no |
| 2000 to 2009 against 2010 to 2019 | rows of no partial method only | point | 285 | 89 | 0 | no |
| 2000 to 2009 against 2010 to 2019 | rows of no partial method only | tvd | 285 | 58 | 0 | no |
| 2000 to 2009 against 2010 to 2019 | largest BioProject left out | point | 285 | 59 | 0 | no |
| 2000 to 2009 against 2010 to 2019 | largest BioProject left out | tvd | 285 | 38 | 0 | no |
| 2000 to 2009 against 2010 to 2019 | one vote per lineage cluster | point | 285 | 141 | 0 | no |
| 2000 to 2009 against 2010 to 2019 | one vote per lineage cluster | tvd | 285 | 84 | 0 | no |
| 2000 to 2009 against 2010 to 2019 | lineage clusters voting in both periods | point | 285 | 126 | 0 | no |
| 2000 to 2009 against 2010 to 2019 | lineage clusters voting in both periods | tvd | 285 | 68 | 0 | no |
| 2000 to 2009 against 2010 to 2019 | share per genome, no vote | point | 285 | 0 | 0 | no |
| 2000 to 2009 against 2010 to 2019 | share per genome, no vote | tvd | 285 | 0 | 0 | no |

**Dataset cells.** Resistant genomes per dataset are the secondary reading, over the 204 datasets of the primary set. On the primary periods 30 dataset cells pass on the point share and 29 on the mix, on the secondary periods 11 and 11, and none moves.

| Periods | Guard | Statistic | Dataset cells | Pass | Move |
|---|---|---|---|---|---|
| 2010 to 2014 against 2015 to 2019 | one vote per BioProject (decides) | point | 204 | 30 | 0 |
| 2010 to 2014 against 2015 to 2019 | one vote per BioProject (decides) | tvd | 204 | 29 | 0 |
| 2010 to 2014 against 2015 to 2019 | share per genome, no vote | point | 204 | 0 | 0 |
| 2010 to 2014 against 2015 to 2019 | share per genome, no vote | tvd | 204 | 0 | 0 |
| 2000 to 2009 against 2010 to 2019 | one vote per BioProject (decides) | point | 204 | 11 | 0 |
| 2000 to 2009 against 2010 to 2019 | one vote per BioProject (decides) | tvd | 204 | 11 | 0 |
| 2000 to 2009 against 2010 to 2019 | share per genome, no vote | point | 204 | 0 | 0 |
| 2000 to 2009 against 2010 to 2019 | share per genome, no vote | tvd | 204 | 0 | 0 |

**Time-split link, exploratory.** Over the 24 datasets of `docs/time-split.md`, Spearman's rank correlation of the share of resistant later test genomes with a catalogue call no resistant training genome carries (a point row is one call, "point"), against each learner's past-less-random difference. On recall it runs from -0.257 (`scm`) to 0.268 (`rf`), and on log loss it runs from -0.186 (`gbt`) to 0.371 (`cart`). It is a retrospective split inside one database and decides nothing.

| Learner | Datasets | Recall difference, Spearman | Log loss difference, Spearman |
|---|---|---|---|
| set covering machine (`scm`) | 24 | -0.257 | 0.235 |
| decision tree (`cart`) | 24 | 0.184 | 0.371 |
| random forest (`rf`) | 24 | 0.268 | -0.047 |
| boosted trees (`gbt`) | 24 | 0.151 | -0.186 |
| lasso (`lasso`) | 24 | -0.060 | -0.092 |

## Seen after

These points were seen in the tables after the reading. The rules did not ask for them. They change no reading above and decide nothing. The entry "Mechanism split, result" in `docs/assumption-plan-corrections.md` records the rule in two parts, the two names of one class and the scan of the record.

**Quinolones.** Of the 19 datasets of the M1 reading, 5 have the drug class quinolone: *Klebsiella pneumoniae* with ciprofloxacin and levofloxacin, *Salmonella enterica* with nalidixic acid, and *Shigella sonnei* with ciprofloxacin and nalidixic acid. On each of the 5, every learner is point ahead, and GRM's known-gene rule has a recall of 1.000 on gene only and 0.000 on point only. On the 14 others, of the classes cephalosporin (1st-2nd gen), cephalosporin (3rd gen), carbapenem, penicillin + inhibitor and aminoglycoside, `scm`, `gbt` and `lasso` are gene ahead on 14 of 14. `cart` is gene ahead on 13 and point ahead on 1, *Klebsiella pneumoniae* with doripenem (0.577 against 0.600). `rf` is gene ahead on 13 and the same on that one (1.000 against 1.000). The direction of the learners' difference follows the drug class here, and a median over the 19 mixes the two.

**Sulfamethoxazole/trimethoprim.** Its rule has two parts, `C:TRIMETHOPRIM & C:SULFONAMIDE`, the only rule of the 237 with two. A genome's stratum reads its rows that meet either part, and the program's call needs both. So on 9 of the 10 sulfamethoxazole/trimethoprim datasets the program does not call every resistant genome of a stratum. On the primary set it calls 38,673 of the 39,814 gene-only resistant genomes and 7,959 of the 7,963 point-only ones, and the 1,141 and 4 it does not call are all on sulfamethoxazole/trimethoprim. On every other dataset it calls every resistant genome of every stratum. No sulfamethoxazole/trimethoprim dataset is among the 19 of the reading or the 5 in M, so no M1 reading changes. The counts of the section Strata include them.

| Species | Set | Gene only, called by the program | Point only, called by the program |
|---|---|---|---|
| *Acinetobacter baumannii* | primary | 65 of 905 | 0 of 0 |
| *Corynebacterium diphtheriae* | primary | 10 of 27 | 0 of 0 |
| *Enterobacter cloacae* | primary | 39 of 41 | 0 of 0 |
| *Escherichia coli* | M | 545 of 564 | 0 of 0 |
| *Klebsiella pneumoniae* | primary | 935 of 1,023 | 0 of 0 |
| *Klebsiella quasipneumoniae* | primary | 53 of 55 | 0 of 0 |
| *Salmonella enterica* | primary | 228 of 272 | 0 of 0 |
| *Shigella sonnei* | primary | 492 of 522 | 0 of 0 |
| *Staphylococcus aureus* | primary | 0 of 118 | 2 of 6 |
| *Streptococcus pneumoniae* | primary | 0 of 0 | 0 of 0 |

**Fusidic acid.** In *Staphylococcus aureus* the program's Class column names fusidic acid two ways, `FUSIDIC ACID` and `FUSIDIC_ACID`, so the cell rule makes two cells. Among the voting genomes, on both pairs of periods, `FUSIDIC ACID` has a point share of 0.000 in each period and `FUSIDIC_ACID` a gene share of 0.000. No other two of the 32 class names differ only by a space or an underscore. Neither cell moves.

**Before Holm.** On the primary periods, 2 of the 110 point p values and 4 of the 76 mix p values are under 0.05 before Holm. The smallest point p is 0.003, *Acinetobacter baumannii* with `BETA-LACTAM` (a point share of 0.018 against 0.379, 7 and 5 voting BioProjects), and its Holm p is 0.330. A share describes deposits.

## Limits

**What a table can say.** A table says which catalogue calls the resistant genomes carry. It does not say what causes resistance in a genome, and it names no new marker. A gene family is named only as the AMRFinderPlus database names it, over 5 or more genomes.

**Deposits.** The genomes are BV-BRC deposits chosen for having lab results, capped at 2,000 a species. A share is no prevalence, and a trend over collection years describes deposits, not resistance in any population.

**One database.** A split by lineage or by study is a retrospective split inside one database, never external validation.

**Few datasets.** The M1 reading stands on 19 datasets in 5 species, 11 of them in *Klebsiella pneumoniae*. Datasets of one species can share genomes, and the sign test counts them apart. Under 6 species the method gives no species interval.

**A rule in two parts.** On sulfamethoxazole/trimethoprim a stratum counts a row that meets one part of the rule (section Seen after). The maintainer has not ruled on it.

**Labels.** The datasets in M wait for the maintainer's label rule. They are marked and never join the primary reading.

**M4.** The lineage curve for Mycobacterium tuberculosis has not run. TBProfiler 6.6.3 must run again first on the 2,000 genomes.

**Markers.** No marker is found, named or ranked. The page holds aggregates over datasets, species, cells and learners only: no genome id, no BioProject accession, no sequence and no position.

## Method

**Checks of the run.** Before any table was read, a check program read each species' copy and printed counts only. All six checks passed. `checks.txt` in `reference/mechanism-split/` holds its printout.

| Check | Result |
|---|---|
| species | 21 copies, the 21 species of A |
| `SHA256SUMS` | holds in 21 of 21 species, 210 of 210 lines, every file listed and none more |
| `DONE` | 21 of 21 name their species, commit `11c299674c68`, 10 results and the sha256 of their `SHA256SUMS`; no `FAILED` file |
| first run | 19 first-run copies with `DONE`, at commit `253d60e12759`; 19 of 19 `SHA256SUMS` equal the second run's byte for byte; *Streptococcus agalactiae* and *Streptococcus pneumoniae* not in the first run |
| scan | 0 files with a hit in 609 files (189 tables) over 21 species, against each species' own genome ids and against the 23,757 genome ids of the release's 23 species files |
| hard check | 237 of 237 datasets pass, 2,844 of 2,844 rows equal, 237 entered |
| program | 21 species read the program's tables made before and 0 ran it again, and none installed it; 21,625 tables, each equal to its recorded sha256; 21 of 21 database logs give 43 files equal to the kit's list |

**Reading.** The reading ran once, over each species' `DONE` and `results/` as its run wrote them:

```
python3 -I scripts/mechanism/combine.py --run runs/<species slug> (21, in sorted order) --contrasts reference/time-split/contrasts.tsv --out combined
python3 -I scripts/mechanism/page_tables.py write <the reading> scripts/c3-inputs/lists.tsv reference/time-split/contrasts.tsv 11c299674c681581bbdaea0717058ed47034cd68 reference/mechanism-split
python3 -I scripts/mechanism/scan.py --ids <every genome id of the release> reference/mechanism-split
python3 -I scripts/mechanism/page_tables.py print reference/mechanism-split
```

`combine.py` printed `combine: 21 species runs; m1 rows 28, m2 rows 10752, time-split link rows 10`. `reading.txt` lists the commands with the sha256 of each table `combine.py` wrote, and `SHA256SUMS` gives every other file's.

**Tables.** `scripts/mechanism/page_tables.py write` checked the reading against six gates and wrote `reference/mechanism-split/`. All six passed:

| Gate | What passed |
|---|---|
| 1. layout | the reading holds reading.txt, checks.txt, combine.py's five tables under combined/ and runs/ |
| 2. species | runs/ holds exactly the species of A in lists.tsv, none dropped and none more |
| 3. runs | each run's DONE names its species, the reading's commit, 10 results and the sha256 of its SHA256SUMS; SHA256SUMS lists the ten tables and each file equals its line; no run is a first-run copy |
| 4. datasets | each run's datasets.tsv holds exactly its species' datasets of A, species and antibiotic as lists.tsv names them, set M exactly where lists.tsv has in_M 1; strata.tsv and pooled.tsv hold only entered datasets, each with the set datasets.tsv gives it |
| 5. floors | every run table passes scan.py's table checks: a recall only over 10 or more resistant genomes, an M2 share or mix only in a cell whose BioProject vote passes, a family named only over 5 or more genomes |
| 6. combined | combined/ is byte for byte what combine.py writes from runs/, species in sorted order, and the time split's contrasts |

**Scan of the record.** `scan.py` over the record's 20 files found 0 genome ids, 0 BioProject accessions, 0 sequences and 0 positions, against the 23,757 genome ids of the release. Its table rules know a run's tables by name, so it flags each of the record's 15 tables as a table it does not know and returns 1. The floors hold on the run tables those 15 are made from (gate 5).

**Tests.** `tests/mechanism/test_page_tables.py` holds 18 SYNTHETIC tests of the table script, on three made-up species, and `tests/mechanism/test_page.py` holds 19 tests of this page against `reference/mechanism-split/`. `scripts/mechanism/page_faults.py` plants each fault below in the table script and runs the test meant to catch it. All 7 were caught.

| Fault | Caught by |
|---|---|
| P1. A dataset of M counted in the primary reading | `Gates.test_a_dataset_of_m_marked_primary_is_refused` |
| P2. A species dropped | `Gates.test_a_species_missing_from_the_reading_is_refused` |
| P3. A study row read as a lineage row | `Tables.test_a_lineage_guard_row_keeps_its_own_figures` |
| P4. A cell under its floor written as a share | `Gates.test_a_share_in_a_cell_under_its_floor_is_refused` |
| P5. A guard read as the deciding figure | `Tables.test_the_deciding_count_reads_the_bioproject_vote_alone` |
| P6. A first-run table mixed in | `Gates.test_a_first_run_copy_is_refused` |
| P7. A species named by its slug, not in full | `Tables.test_species_and_datasets_are_named_in_full` |

**Files.** The reading's tables are in `reference/mechanism-split/`. `scripts/mechanism/page_tables.py print` prints this page's tables from those files, rounded to 3 decimals, and `tests/mechanism/test_page.py` checks that each printed row stands here.

## Datasets

Per species: the datasets of A, those in M, those entered, the learner and catalog rows equal to the reports, the M1 datasets of the primary set with both strata at 10 or more, the M2 genomes with a table, the classes, the cells passing on each pair of periods, and the cells moving.

| Species | Datasets | In M | Entered | Rows equal | M1 primary, both 10 or more | M2 genomes | Classes | Pass, primary | Pass, secondary | Point moves | Mix moves |
|---|---|---|---|---|---|---|---|---|---|---|---|
| *Acinetobacter baumannii* | 18 | 11 | 18 | 216 of 216 | 0 | 1,702 | 15 | 15 | 15 | 0 | 0 |
| *Campylobacter coli* | 7 | 0 | 7 | 84 of 84 | 0 | 391 | 8 | 0 | 0 | 0 | 0 |
| *Campylobacter jejuni* | 4 | 0 | 4 | 48 of 48 | 0 | 508 | 9 | 0 | 0 | 0 | 0 |
| *Clostridioides difficile* | 4 | 0 | 4 | 48 of 48 | 0 | 204 | 15 | 0 | 0 | 0 | 0 |
| *Corynebacterium diphtheriae* | 3 | 0 | 3 | 36 of 36 | 0 | 243 | 11 | 0 | 0 | 0 | 0 |
| *Enterobacter cloacae* | 9 | 0 | 9 | 108 of 108 | 0 | 270 | 15 | 0 | 0 | 0 | 0 |
| *Enterococcus faecium* | 13 | 2 | 13 | 156 of 156 | 0 | 2,000 | 17 | 17 | 0 | 0 | 0 |
| *Escherichia coli* | 36 | 7 | 36 | 432 of 432 | 3 | 2,000 | 21 | 21 | 21 | 0 | 0 |
| *Klebsiella aerogenes* | 3 | 0 | 3 | 36 of 36 | 0 | 93 | 13 | 0 | 0 | 0 | 0 |
| *Klebsiella oxytoca* | 4 | 0 | 4 | 48 of 48 | 0 | 112 | 11 | 0 | 0 | 0 | 0 |
| *Klebsiella pneumoniae* | 37 | 6 | 37 | 444 of 444 | 11 | 2,000 | 19 | 19 | 0 | 0 | 0 |
| *Klebsiella quasipneumoniae* | 3 | 0 | 3 | 36 of 36 | 0 | 109 | 14 | 0 | 0 | 0 | 0 |
| *Neisseria gonorrhoeae* | 7 | 1 | 7 | 84 of 84 | 0 | 2,000 | 9 | 0 | 0 | 0 | 0 |
| *Pseudomonas aeruginosa* | 13 | 1 | 13 | 156 of 156 | 2 | 1,421 | 15 | 0 | 15 | 0 | 0 |
| *Salmonella enterica* | 23 | 0 | 23 | 276 of 276 | 1 | 2,000 | 17 | 17 | 17 | 0 | 0 |
| *Shigella flexneri* | 4 | 0 | 4 | 48 of 48 | 0 | 362 | 13 | 0 | 0 | 0 | 0 |
| *Shigella sonnei* | 18 | 0 | 18 | 216 of 216 | 2 | 1,984 | 19 | 0 | 0 | 0 | 0 |
| *Staphylococcus aureus* | 13 | 1 | 13 | 156 of 156 | 0 | 1,994 | 21 | 21 | 21 | 0 | 0 |
| *Streptococcus agalactiae* | 3 | 0 | 3 | 36 of 36 | 0 | 141 | 5 | 0 | 0 | 0 | 0 |
| *Streptococcus pneumoniae* | 12 | 4 | 12 | 144 of 144 | 0 | 2,000 | 9 | 0 | 0 | 0 | 0 |
| *Vibrio cholerae* | 3 | 0 | 3 | 36 of 36 | 0 | 91 | 9 | 0 | 0 | 0 | 0 |

Per dataset of A, species and antibiotic in full: the set, whether it entered, the resistant genomes per stratum on the lineage split, and each learner's recall on gene only less recall on point only where both hold 10 or more resistant genomes, `none` otherwise. A dataset's figure decides nothing alone.

| Species | Antibiotic | Set | Entered | Gene only | Point only | Both | Neither | `scm` | `cart` | `rf` | `gbt` | `lasso` |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| *Acinetobacter baumannii* | amikacin | M | yes | 635 | 0 | 0 | 150 | none | none | none | none | none |
| *Acinetobacter baumannii* | ampicillin/sulbactam | M | yes | 496 | 0 | 59 | 0 | none | none | none | none | none |
| *Acinetobacter baumannii* | cefepime | primary | yes | 351 | 0 | 61 | 0 | none | none | none | none | none |
| *Acinetobacter baumannii* | cefotaxime | primary | yes | 767 | 0 | 11 | 0 | none | none | none | none | none |
| *Acinetobacter baumannii* | ceftazidime | M | yes | 1017 | 0 | 70 | 0 | none | none | none | none | none |
| *Acinetobacter baumannii* | ceftriaxone | M | yes | 949 | 0 | 62 | 0 | none | none | none | none | none |
| *Acinetobacter baumannii* | ciprofloxacin | M | yes | 0 | 1309 | 1 | 9 | none | none | none | none | none |
| *Acinetobacter baumannii* | doripenem | primary | yes | 61 | 0 | 0 | 0 | none | none | none | none | none |
| *Acinetobacter baumannii* | doxycycline | primary | yes | 90 | 0 | 0 | 17 | none | none | none | none | none |
| *Acinetobacter baumannii* | gentamicin | M | yes | 1165 | 0 | 0 | 39 | none | none | none | none | none |
| *Acinetobacter baumannii* | imipenem | M | yes | 705 | 0 | 88 | 1 | none | none | none | none | none |
| *Acinetobacter baumannii* | levofloxacin | M | yes | 0 | 759 | 0 | 1 | none | none | none | none | none |
| *Acinetobacter baumannii* | meropenem | M | yes | 393 | 0 | 62 | 0 | none | none | none | none | none |
| *Acinetobacter baumannii* | minocycline | primary | yes | 20 | 0 | 0 | 8 | none | none | none | none | none |
| *Acinetobacter baumannii* | moxifloxacin | primary | yes | 0 | 60 | 0 | 0 | none | none | none | none | none |
| *Acinetobacter baumannii* | sulfamethoxazole/trimethoprim | primary | yes | 905 | 0 | 0 | 57 | none | none | none | none | none |
| *Acinetobacter baumannii* | tetracycline | M | yes | 667 | 1 | 0 | 91 | none | none | none | none | none |
| *Acinetobacter baumannii* | tobramycin | M | yes | 475 | 0 | 0 | 206 | none | none | none | none | none |
| *Campylobacter coli* | azithromycin | primary | yes | 1 | 77 | 0 | 0 | none | none | none | none | none |
| *Campylobacter coli* | ciprofloxacin | primary | yes | 0 | 132 | 0 | 0 | none | none | none | none | none |
| *Campylobacter coli* | clindamycin | primary | yes | 11 | 0 | 0 | 83 | none | none | none | none | none |
| *Campylobacter coli* | erythromycin | primary | yes | 1 | 77 | 0 | 0 | none | none | none | none | none |
| *Campylobacter coli* | gentamicin | primary | yes | 64 | 0 | 0 | 1 | none | none | none | none | none |
| *Campylobacter coli* | nalidixic acid | primary | yes | 0 | 116 | 0 | 3 | none | none | none | none | none |
| *Campylobacter coli* | tetracycline | primary | yes | 301 | 0 | 0 | 3 | none | none | none | none | none |
| *Campylobacter jejuni* | ciprofloxacin | primary | yes | 0 | 93 | 16 | 4 | none | none | none | none | none |
| *Campylobacter jejuni* | clindamycin | primary | yes | 1 | 0 | 0 | 24 | none | none | none | none | none |
| *Campylobacter jejuni* | nalidixic acid | primary | yes | 0 | 92 | 15 | 3 | none | none | none | none | none |
| *Campylobacter jejuni* | tetracycline | primary | yes | 355 | 0 | 0 | 2 | none | none | none | none | none |
| *Clostridioides difficile* | azithromycin | primary | yes | 22 | 0 | 0 | 69 | none | none | none | none | none |
| *Clostridioides difficile* | ceftriaxone | primary | yes | 3 | 0 | 87 | 0 | none | none | none | none | none |
| *Clostridioides difficile* | clarithromycin | primary | yes | 22 | 0 | 0 | 70 | none | none | none | none | none |
| *Clostridioides difficile* | moxifloxacin | primary | yes | 0 | 86 | 0 | 0 | none | none | none | none | none |
| *Corynebacterium diphtheriae* | sulfamethoxazole/trimethoprim | primary | yes | 27 | 0 | 0 | 1 | none | none | none | none | none |
| *Corynebacterium diphtheriae* | tetracycline | primary | yes | 40 | 0 | 0 | 6 | none | none | none | none | none |
| *Corynebacterium diphtheriae* | trimethoprim | primary | yes | 10 | 0 | 0 | 20 | none | none | none | none | none |
| *Enterobacter cloacae* | cefepime | primary | yes | 70 | 0 | 0 | 0 | none | none | none | none | none |
| *Enterobacter cloacae* | ciprofloxacin | primary | yes | 11 | 0 | 56 | 0 | none | none | none | none | none |
| *Enterobacter cloacae* | gentamicin | primary | yes | 26 | 0 | 0 | 1 | none | none | none | none | none |
| *Enterobacter cloacae* | imipenem | primary | yes | 53 | 0 | 0 | 5 | none | none | none | none | none |
| *Enterobacter cloacae* | levofloxacin | primary | yes | 7 | 0 | 51 | 0 | none | none | none | none | none |
| *Enterobacter cloacae* | meropenem | primary | yes | 36 | 0 | 0 | 4 | none | none | none | none | none |
| *Enterobacter cloacae* | sulfamethoxazole/trimethoprim | primary | yes | 41 | 0 | 0 | 0 | none | none | none | none | none |
| *Enterobacter cloacae* | tobramycin | primary | yes | 25 | 0 | 0 | 20 | none | none | none | none | none |
| *Enterobacter cloacae* | trimethoprim | primary | yes | 31 | 0 | 0 | 7 | none | none | none | none | none |
| *Enterococcus faecium* | ampicillin | primary | yes | 0 | 1288 | 1 | 2 | none | none | none | none | none |
| *Enterococcus faecium* | clindamycin | primary | yes | 543 | 0 | 0 | 76 | none | none | none | none | none |
| *Enterococcus faecium* | dalfopristin/quinupristin | primary | yes | 112 | 0 | 0 | 0 | none | none | none | none | none |
| *Enterococcus faecium* | daptomycin | M | yes | 0 | 47 | 0 | 8 | none | none | none | none | none |
| *Enterococcus faecium* | doxycycline | primary | yes | 133 | 0 | 0 | 3 | none | none | none | none | none |
| *Enterococcus faecium* | erythromycin | primary | yes | 1061 | 0 | 0 | 2 | none | none | none | none | none |
| *Enterococcus faecium* | gentamicin | primary | yes | 396 | 0 | 0 | 0 | none | none | none | none | none |
| *Enterococcus faecium* | kanamycin | primary | yes | 612 | 0 | 0 | 0 | none | none | none | none | none |
| *Enterococcus faecium* | linezolid | M | yes | 38 | 5 | 0 | 27 | none | none | none | none | none |
| *Enterococcus faecium* | streptomycin | primary | yes | 484 | 0 | 0 | 4 | none | none | none | none | none |
| *Enterococcus faecium* | teicoplanin | primary | yes | 0 | 0 | 0 | 516 | none | none | none | none | none |
| *Enterococcus faecium* | tetracycline | primary | yes | 663 | 0 | 1 | 53 | none | none | none | none | none |
| *Enterococcus faecium* | vancomycin | primary | yes | 896 | 0 | 0 | 0 | none | none | none | none | none |
| *Escherichia coli* | amikacin | primary | yes | 23 | 0 | 0 | 7 | none | none | none | none | none |
| *Escherichia coli* | amoxicillin/clavulanate | M | yes | 389 | 20 | 75 | 64 | 0.656 | 0.630 | 0.553 | 0.583 | 0.570 |
| *Escherichia coli* | ampicillin | M | yes | 1175 | 16 | 376 | 29 | 0.909 | 0.866 | 0.732 | 0.835 | 0.876 |
| *Escherichia coli* | ampicillin/sulbactam | primary | yes | 106 | 1 | 24 | 11 | none | none | none | none | none |
| *Escherichia coli* | aztreonam | primary | yes | 190 | 3 | 18 | 14 | none | none | none | none | none |
| *Escherichia coli* | cefazolin | primary | yes | 259 | 3 | 46 | 14 | none | none | none | none | none |
| *Escherichia coli* | cefepime | primary | yes | 168 | 0 | 40 | 3 | none | none | none | none | none |
| *Escherichia coli* | cefotaxime | primary | yes | 521 | 2 | 84 | 8 | none | none | none | none | none |
| *Escherichia coli* | cefoxitin | primary | yes | 318 | 8 | 69 | 4 | none | none | none | none | none |
| *Escherichia coli* | cefpodoxime | primary | yes | 79 | 0 | 24 | 0 | none | none | none | none | none |
| *Escherichia coli* | cefpodoxime/clavulanate | primary | yes | 91 | 0 | 10 | 0 | none | none | none | none | none |
| *Escherichia coli* | ceftazidime | primary | yes | 511 | 4 | 78 | 8 | none | none | none | none | none |
| *Escherichia coli* | ceftiofur | primary | yes | 104 | 1 | 0 | 0 | none | none | none | none | none |
| *Escherichia coli* | ceftriaxone | primary | yes | 371 | 3 | 55 | 9 | none | none | none | none | none |
| *Escherichia coli* | cefuroxime | primary | yes | 211 | 12 | 26 | 23 | 0.863 | 0.858 | 0.957 | 0.882 | 0.867 |
| *Escherichia coli* | cephalothin | primary | yes | 226 | 13 | 23 | 16 | 0.889 | 0.881 | 0.833 | 0.679 | 0.675 |
| *Escherichia coli* | chloramphenicol | primary | yes | 208 | 3 | 24 | 15 | none | none | none | none | none |
| *Escherichia coli* | ciprofloxacin | M | yes | 6 | 494 | 316 | 10 | none | none | none | none | none |
| *Escherichia coli* | colistin | primary | yes | 43 | 0 | 2 | 5 | none | none | none | none | none |
| *Escherichia coli* | doxycycline | primary | yes | 138 | 0 | 0 | 9 | none | none | none | none | none |
| *Escherichia coli* | ertapenem | primary | yes | 156 | 11 | 21 | 17 | 0.949 | 0.638 | 0.821 | 0.722 | 0.917 |
| *Escherichia coli* | erythromycin | primary | yes | 65 | 0 | 0 | 6 | none | none | none | none | none |
| *Escherichia coli* | gentamicin | M | yes | 431 | 0 | 0 | 24 | none | none | none | none | none |
| *Escherichia coli* | imipenem | primary | yes | 74 | 5 | 27 | 7 | none | none | none | none | none |
| *Escherichia coli* | levofloxacin | primary | yes | 1 | 183 | 108 | 2 | none | none | none | none | none |
| *Escherichia coli* | meropenem | M | yes | 54 | 12 | 8 | 8 | 0.620 | 0.130 | 0.056 | 0.620 | 0.852 |
| *Escherichia coli* | nalidixic acid | primary | yes | 2 | 97 | 62 | 3 | none | none | none | none | none |
| *Escherichia coli* | nitrofurantoin | primary | yes | 1 | 23 | 0 | 22 | none | none | none | none | none |
| *Escherichia coli* | norfloxacin | M | yes | 3 | 34 | 3 | 1 | none | none | none | none | none |
| *Escherichia coli* | piperacillin/tazobactam | primary | yes | 209 | 2 | 62 | 19 | none | none | none | none | none |
| *Escherichia coli* | streptomycin | primary | yes | 240 | 0 | 0 | 14 | none | none | none | none | none |
| *Escherichia coli* | sulfamethoxazole | primary | yes | 153 | 0 | 0 | 8 | none | none | none | none | none |
| *Escherichia coli* | sulfamethoxazole/trimethoprim | M | yes | 564 | 0 | 0 | 16 | none | none | none | none | none |
| *Escherichia coli* | tetracycline | primary | yes | 441 | 0 | 0 | 23 | none | none | none | none | none |
| *Escherichia coli* | tobramycin | primary | yes | 158 | 0 | 0 | 185 | none | none | none | none | none |
| *Escherichia coli* | trimethoprim | primary | yes | 410 | 0 | 0 | 21 | none | none | none | none | none |
| *Klebsiella aerogenes* | ertapenem | primary | yes | 10 | 4 | 0 | 34 | none | none | none | none | none |
| *Klebsiella aerogenes* | imipenem | primary | yes | 9 | 2 | 0 | 15 | none | none | none | none | none |
| *Klebsiella aerogenes* | meropenem | primary | yes | 4 | 2 | 0 | 20 | none | none | none | none | none |
| *Klebsiella oxytoca* | cefotaxime | primary | yes | 56 | 0 | 1 | 0 | none | none | none | none | none |
| *Klebsiella oxytoca* | cefoxitin | primary | yes | 47 | 0 | 0 | 0 | none | none | none | none | none |
| *Klebsiella oxytoca* | cefpodoxime | primary | yes | 46 | 0 | 1 | 0 | none | none | none | none | none |
| *Klebsiella oxytoca* | ertapenem | primary | yes | 76 | 0 | 0 | 2 | none | none | none | none | none |
| *Klebsiella pneumoniae* | amikacin | M | yes | 272 | 0 | 0 | 46 | none | none | none | none | none |
| *Klebsiella pneumoniae* | amoxicillin/clavulanate | M | yes | 288 | 4 | 298 | 1 | none | none | none | none | none |
| *Klebsiella pneumoniae* | ampicillin/sulbactam | primary | yes | 173 | 28 | 310 | 25 | 0.527 | 0.296 | 0.072 | 0.228 | 0.473 |
| *Klebsiella pneumoniae* | aztreonam | primary | yes | 498 | 7 | 18 | 170 | none | none | none | none | none |
| *Klebsiella pneumoniae* | cefazolin | primary | yes | 295 | 26 | 414 | 28 | 0.725 | 0.543 | 0.610 | 0.571 | 0.803 |
| *Klebsiella pneumoniae* | cefepime | primary | yes | 215 | 8 | 394 | 2 | none | none | none | none | none |
| *Klebsiella pneumoniae* | cefotetan | primary | yes | 3 | 1 | 21 | 0 | none | none | none | none | none |
| *Klebsiella pneumoniae* | cefoxitin | primary | yes | 173 | 20 | 340 | 17 | 0.457 | 0.341 | 0.320 | 0.276 | 0.153 |
| *Klebsiella pneumoniae* | cefpodoxime | primary | yes | 118 | 0 | 62 | 0 | none | none | none | none | none |
| *Klebsiella pneumoniae* | cefpodoxime/clavulanate | primary | yes | 121 | 0 | 60 | 1 | none | none | none | none | none |
| *Klebsiella pneumoniae* | ceftazidime | primary | yes | 302 | 18 | 510 | 30 | 0.492 | 0.429 | 0.325 | 0.407 | 0.421 |
| *Klebsiella pneumoniae* | ceftazidime/avibactam | primary | yes | 5 | 3 | 71 | 0 | none | none | none | none | none |
| *Klebsiella pneumoniae* | ceftriaxone | primary | yes | 359 | 12 | 453 | 19 | 0.630 | 0.611 | 0.311 | 0.369 | 0.516 |
| *Klebsiella pneumoniae* | cefuroxime | primary | yes | 214 | 14 | 215 | 25 | 0.868 | 0.349 | 0.286 | 0.338 | 0.396 |
| *Klebsiella pneumoniae* | chloramphenicol | primary | yes | 276 | 0 | 26 | 0 | none | none | none | none | none |
| *Klebsiella pneumoniae* | ciprofloxacin | primary | yes | 186 | 18 | 892 | 0 | -0.656 | -0.774 | -0.306 | -0.565 | -0.559 |
| *Klebsiella pneumoniae* | colistin | primary | yes | 0 | 24 | 0 | 88 | none | none | none | none | none |
| *Klebsiella pneumoniae* | doripenem | primary | yes | 26 | 10 | 178 | 3 | 0.131 | -0.023 | 0.000 | 0.169 | 0.100 |
| *Klebsiella pneumoniae* | doxycycline | primary | yes | 36 | 0 | 10 | 0 | none | none | none | none | none |
| *Klebsiella pneumoniae* | ertapenem | primary | yes | 213 | 22 | 429 | 7 | 0.850 | 0.580 | 0.071 | 0.413 | 0.070 |
| *Klebsiella pneumoniae* | fosfomycin | primary | yes | 83 | 0 | 0 | 0 | none | none | none | none | none |
| *Klebsiella pneumoniae* | gentamicin | M | yes | 656 | 0 | 0 | 33 | none | none | none | none | none |
| *Klebsiella pneumoniae* | imipenem | primary | yes | 113 | 11 | 476 | 15 | 0.865 | 0.865 | 0.652 | 0.752 | 0.588 |
| *Klebsiella pneumoniae* | levofloxacin | primary | yes | 11 | 16 | 522 | 0 | -0.750 | -1.000 | -0.750 | -1.000 | -0.812 |
| *Klebsiella pneumoniae* | meropenem | M | yes | 150 | 23 | 549 | 14 | 0.556 | 0.760 | 0.269 | 0.319 | 0.416 |
| *Klebsiella pneumoniae* | nitrofurantoin | primary | yes | 242 | 3 | 110 | 20 | none | none | none | none | none |
| *Klebsiella pneumoniae* | norfloxacin | M | yes | 9 | 0 | 20 | 0 | none | none | none | none | none |
| *Klebsiella pneumoniae* | ofloxacin | primary | yes | 12 | 0 | 32 | 0 | none | none | none | none | none |
| *Klebsiella pneumoniae* | piperacillin/tazobactam | M | yes | 319 | 25 | 512 | 12 | 0.600 | 0.050 | 0.141 | 0.375 | 0.434 |
| *Klebsiella pneumoniae* | polymyxin b | primary | yes | 0 | 6 | 0 | 42 | none | none | none | none | none |
| *Klebsiella pneumoniae* | streptomycin | primary | yes | 34 | 0 | 0 | 20 | none | none | none | none | none |
| *Klebsiella pneumoniae* | sulfamethoxazole/trimethoprim | primary | yes | 1023 | 0 | 0 | 19 | none | none | none | none | none |
| *Klebsiella pneumoniae* | tetracycline | primary | yes | 342 | 0 | 40 | 5 | none | none | none | none | none |
| *Klebsiella pneumoniae* | ticarcillin/clavulanate | primary | yes | 48 | 1 | 57 | 1 | none | none | none | none | none |
| *Klebsiella pneumoniae* | tigecycline | primary | yes | 211 | 1 | 28 | 17 | none | none | none | none | none |
| *Klebsiella pneumoniae* | tobramycin | primary | yes | 559 | 0 | 0 | 71 | none | none | none | none | none |
| *Klebsiella pneumoniae* | trimethoprim | primary | yes | 106 | 0 | 0 | 6 | none | none | none | none | none |
| *Klebsiella quasipneumoniae* | cefotaxime | primary | yes | 62 | 0 | 0 | 0 | none | none | none | none | none |
| *Klebsiella quasipneumoniae* | ceftazidime | primary | yes | 59 | 0 | 0 | 0 | none | none | none | none | none |
| *Klebsiella quasipneumoniae* | sulfamethoxazole/trimethoprim | primary | yes | 55 | 0 | 0 | 2 | none | none | none | none | none |
| *Neisseria gonorrhoeae* | azithromycin | M | yes | 0 | 404 | 2 | 2 | none | none | none | none | none |
| *Neisseria gonorrhoeae* | cefixime | primary | yes | 0 | 103 | 0 | 0 | none | none | none | none | none |
| *Neisseria gonorrhoeae* | ceftriaxone | primary | yes | 0 | 39 | 0 | 1 | none | none | none | none | none |
| *Neisseria gonorrhoeae* | ciprofloxacin | primary | yes | 0 | 1083 | 0 | 3 | none | none | none | none | none |
| *Neisseria gonorrhoeae* | erythromycin | primary | yes | 0 | 97 | 1 | 0 | none | none | none | none | none |
| *Neisseria gonorrhoeae* | penicillin | primary | yes | 0 | 481 | 406 | 1 | none | none | none | none | none |
| *Neisseria gonorrhoeae* | tetracycline | primary | yes | 0 | 237 | 122 | 0 | none | none | none | none | none |
| *Pseudomonas aeruginosa* | amikacin | primary | yes | 63 | 22 | 21 | 21 | 0.332 | 0.241 | 0.206 | 0.512 | 0.630 |
| *Pseudomonas aeruginosa* | aztreonam | primary | yes | 91 | 0 | 23 | 0 | none | none | none | none | none |
| *Pseudomonas aeruginosa* | cefepime | primary | yes | 0 | 0 | 122 | 0 | none | none | none | none | none |
| *Pseudomonas aeruginosa* | ceftazidime | primary | yes | 2 | 0 | 370 | 0 | none | none | none | none | none |
| *Pseudomonas aeruginosa* | ceftazidime/avibactam | primary | yes | 0 | 54 | 85 | 0 | none | none | none | none | none |
| *Pseudomonas aeruginosa* | ceftolozane/tazobactam | primary | yes | 0 | 44 | 86 | 2 | none | none | none | none | none |
| *Pseudomonas aeruginosa* | ciprofloxacin | M | yes | 0 | 360 | 11 | 20 | none | none | none | none | none |
| *Pseudomonas aeruginosa* | gentamicin | primary | yes | 64 | 4 | 16 | 17 | none | none | none | none | none |
| *Pseudomonas aeruginosa* | imipenem | primary | yes | 0 | 158 | 77 | 0 | none | none | none | none | none |
| *Pseudomonas aeruginosa* | levofloxacin | primary | yes | 0 | 309 | 3 | 18 | none | none | none | none | none |
| *Pseudomonas aeruginosa* | meropenem | primary | yes | 0 | 371 | 116 | 1 | none | none | none | none | none |
| *Pseudomonas aeruginosa* | piperacillin/tazobactam | primary | yes | 0 | 0 | 87 | 0 | none | none | none | none | none |
| *Pseudomonas aeruginosa* | tobramycin | primary | yes | 69 | 34 | 39 | 52 | 0.192 | 0.177 | 0.148 | 0.134 | 0.148 |
| *Salmonella enterica* | amoxicillin/clavulanate | primary | yes | 441 | 0 | 0 | 30 | none | none | none | none | none |
| *Salmonella enterica* | ampicillin | primary | yes | 1338 | 0 | 0 | 68 | none | none | none | none | none |
| *Salmonella enterica* | ampicillin/sulbactam | primary | yes | 17 | 0 | 0 | 109 | none | none | none | none | none |
| *Salmonella enterica* | azithromycin | primary | yes | 157 | 1 | 0 | 26 | none | none | none | none | none |
| *Salmonella enterica* | cefazolin | primary | yes | 16 | 0 | 0 | 103 | none | none | none | none | none |
| *Salmonella enterica* | cefotaxime | primary | yes | 279 | 0 | 0 | 40 | none | none | none | none | none |
| *Salmonella enterica* | cefoxitin | primary | yes | 385 | 0 | 0 | 10 | none | none | none | none | none |
| *Salmonella enterica* | ceftazidime | primary | yes | 21 | 0 | 0 | 9 | none | none | none | none | none |
| *Salmonella enterica* | ceftiofur | primary | yes | 456 | 0 | 0 | 5 | none | none | none | none | none |
| *Salmonella enterica* | ceftriaxone | primary | yes | 560 | 0 | 0 | 8 | none | none | none | none | none |
| *Salmonella enterica* | chloramphenicol | primary | yes | 654 | 0 | 2 | 57 | none | none | none | none | none |
| *Salmonella enterica* | ciprofloxacin | primary | yes | 9 | 211 | 106 | 27 | none | none | none | none | none |
| *Salmonella enterica* | gentamicin | primary | yes | 367 | 0 | 0 | 48 | none | none | none | none | none |
| *Salmonella enterica* | kanamycin | primary | yes | 121 | 0 | 0 | 33 | none | none | none | none | none |
| *Salmonella enterica* | nalidixic acid | primary | yes | 27 | 316 | 106 | 51 | -0.963 | -1.000 | -0.760 | -0.957 | -0.946 |
| *Salmonella enterica* | spectinomycin | primary | yes | 0 | 0 | 0 | 208 | none | none | none | none | none |
| *Salmonella enterica* | streptomycin | primary | yes | 621 | 0 | 0 | 19 | none | none | none | none | none |
| *Salmonella enterica* | sulfamethoxazole | primary | yes | 86 | 0 | 0 | 6 | none | none | none | none | none |
| *Salmonella enterica* | sulfamethoxazole/trimethoprim | primary | yes | 272 | 0 | 0 | 24 | none | none | none | none | none |
| *Salmonella enterica* | sulfathiazole | primary | yes | 175 | 0 | 0 | 0 | none | none | none | none | none |
| *Salmonella enterica* | sulfisoxazole | primary | yes | 437 | 0 | 0 | 6 | none | none | none | none | none |
| *Salmonella enterica* | tetracycline | primary | yes | 1038 | 0 | 0 | 80 | none | none | none | none | none |
| *Salmonella enterica* | trimethoprim | primary | yes | 235 | 0 | 0 | 21 | none | none | none | none | none |
| *Shigella flexneri* | ampicillin | primary | yes | 92 | 4 | 222 | 0 | none | none | none | none | none |
| *Shigella flexneri* | azithromycin | primary | yes | 178 | 0 | 0 | 8 | none | none | none | none | none |
| *Shigella flexneri* | ciprofloxacin | primary | yes | 0 | 24 | 9 | 0 | none | none | none | none | none |
| *Shigella flexneri* | trimethoprim | primary | yes | 274 | 0 | 0 | 6 | none | none | none | none | none |
| *Shigella sonnei* | amikacin | primary | yes | 0 | 0 | 0 | 270 | none | none | none | none | none |
| *Shigella sonnei* | ampicillin | primary | yes | 770 | 1 | 3 | 74 | none | none | none | none | none |
| *Shigella sonnei* | azithromycin | primary | yes | 576 | 0 | 0 | 26 | none | none | none | none | none |
| *Shigella sonnei* | cefoperazone | primary | yes | 226 | 0 | 0 | 0 | none | none | none | none | none |
| *Shigella sonnei* | cefotaxime | primary | yes | 59 | 0 | 0 | 53 | none | none | none | none | none |
| *Shigella sonnei* | cefoxitin | primary | yes | 56 | 0 | 0 | 49 | none | none | none | none | none |
| *Shigella sonnei* | ceftazidime | primary | yes | 50 | 0 | 0 | 48 | none | none | none | none | none |
| *Shigella sonnei* | ceftriaxone | primary | yes | 429 | 0 | 3 | 7 | none | none | none | none | none |
| *Shigella sonnei* | ciprofloxacin | primary | yes | 24 | 561 | 19 | 2 | -0.607 | -0.312 | -0.276 | -0.227 | -0.529 |
| *Shigella sonnei* | colistin | primary | yes | 35 | 0 | 0 | 0 | none | none | none | none | none |
| *Shigella sonnei* | ertapenem | primary | yes | 1 | 2 | 1 | 97 | none | none | none | none | none |
| *Shigella sonnei* | gentamicin | primary | yes | 245 | 0 | 0 | 268 | none | none | none | none | none |
| *Shigella sonnei* | nalidixic acid | primary | yes | 24 | 508 | 9 | 3 | -0.250 | -0.278 | -0.651 | -0.700 | -0.956 |
| *Shigella sonnei* | streptomycin | primary | yes | 448 | 0 | 0 | 23 | none | none | none | none | none |
| *Shigella sonnei* | sulfamethoxazole/trimethoprim | primary | yes | 522 | 0 | 0 | 2 | none | none | none | none | none |
| *Shigella sonnei* | tetracycline | primary | yes | 670 | 0 | 0 | 34 | none | none | none | none | none |
| *Shigella sonnei* | tobramycin | primary | yes | 0 | 0 | 0 | 462 | none | none | none | none | none |
| *Shigella sonnei* | trimethoprim | primary | yes | 1030 | 0 | 0 | 6 | none | none | none | none | none |
| *Staphylococcus aureus* | cefoxitin | primary | yes | 680 | 0 | 2 | 48 | none | none | none | none | none |
| *Staphylococcus aureus* | chloramphenicol | primary | yes | 88 | 0 | 4 | 31 | none | none | none | none | none |
| *Staphylococcus aureus* | ciprofloxacin | primary | yes | 0 | 581 | 0 | 44 | none | none | none | none | none |
| *Staphylococcus aureus* | clindamycin | primary | yes | 636 | 0 | 0 | 57 | none | none | none | none | none |
| *Staphylococcus aureus* | daptomycin | primary | yes | 0 | 2 | 0 | 45 | none | none | none | none | none |
| *Staphylococcus aureus* | erythromycin | primary | yes | 867 | 0 | 0 | 81 | none | none | none | none | none |
| *Staphylococcus aureus* | fusidic acid | primary | yes | 166 | 0 | 0 | 117 | none | none | none | none | none |
| *Staphylococcus aureus* | gentamicin | primary | yes | 323 | 0 | 0 | 44 | none | none | none | none | none |
| *Staphylococcus aureus* | mupirocin | primary | yes | 42 | 8 | 1 | 0 | none | none | none | none | none |
| *Staphylococcus aureus* | oxacillin | M | yes | 527 | 0 | 0 | 41 | none | none | none | none | none |
| *Staphylococcus aureus* | penicillin | primary | yes | 1203 | 1 | 8 | 40 | none | none | none | none | none |
| *Staphylococcus aureus* | sulfamethoxazole/trimethoprim | primary | yes | 118 | 6 | 153 | 8 | none | none | none | none | none |
| *Staphylococcus aureus* | tetracycline | primary | yes | 308 | 0 | 72 | 28 | none | none | none | none | none |
| *Streptococcus agalactiae* | clindamycin | primary | yes | 115 | 0 | 0 | 0 | none | none | none | none | none |
| *Streptococcus agalactiae* | erythromycin | primary | yes | 116 | 0 | 0 | 0 | none | none | none | none | none |
| *Streptococcus agalactiae* | tetracycline | primary | yes | 104 | 0 | 0 | 0 | none | none | none | none | none |
| *Streptococcus pneumoniae* | azithromycin | primary | yes | 136 | 0 | 0 | 11 | none | none | none | none | none |
| *Streptococcus pneumoniae* | cefepime | M | yes | 0 | 0 | 0 | 49 | none | none | none | none | none |
| *Streptococcus pneumoniae* | cefotaxime | M | yes | 0 | 0 | 0 | 43 | none | none | none | none | none |
| *Streptococcus pneumoniae* | ceftriaxone | M | yes | 0 | 0 | 0 | 56 | none | none | none | none | none |
| *Streptococcus pneumoniae* | cefuroxime | primary | yes | 0 | 0 | 0 | 184 | none | none | none | none | none |
| *Streptococcus pneumoniae* | chloramphenicol | primary | yes | 92 | 0 | 0 | 42 | none | none | none | none | none |
| *Streptococcus pneumoniae* | clindamycin | primary | yes | 155 | 0 | 0 | 11 | none | none | none | none | none |
| *Streptococcus pneumoniae* | erythromycin | primary | yes | 436 | 0 | 0 | 24 | none | none | none | none | none |
| *Streptococcus pneumoniae* | meropenem | primary | yes | 0 | 0 | 0 | 72 | none | none | none | none | none |
| *Streptococcus pneumoniae* | penicillin | M | yes | 0 | 0 | 195 | 1 | none | none | none | none | none |
| *Streptococcus pneumoniae* | sulfamethoxazole/trimethoprim | primary | yes | 0 | 0 | 0 | 343 | none | none | none | none | none |
| *Streptococcus pneumoniae* | tetracycline | primary | yes | 294 | 0 | 0 | 20 | none | none | none | none | none |
| *Vibrio cholerae* | nalidixic acid | primary | yes | 0 | 49 | 0 | 0 | none | none | none | none | none |
| *Vibrio cholerae* | sulfamethoxazole | primary | yes | 37 | 0 | 0 | 0 | none | none | none | none | none |
| *Vibrio cholerae* | trimethoprim | primary | yes | 38 | 0 | 0 | 2 | none | none | none | none | none |
