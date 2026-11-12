# Outside-study test: results, species by species

This page reads the result of [outside-study-results.md](outside-study-results.md) one species at a time. It is the result of a retrospective test on newly deposited studies from the same database. That page gives the question, the counts over all species, the limits and the statements the test supports.

Each figure here is a new estimate set against the report's figure, as in "98.9% against 90.3%". For a learner the report's figure is the lineage-split figure. For the known-gene rule it is the fixed rule's figure on the whole cohort.

A learner holds on a dataset when neither its sensitivity nor its specificity reads below. That count is not in the plan and was added after the results were read. The datasets a learner holds on are part of the species' datasets. The two populations are separate counts and are never added together.

The rows named here are in the table "Rows cited in the text" below. Every dataset's row is in the species' `metrics.tsv`, and `reading.txt` gives the same rows in plain words.

## Primary population

These genomes come from studies the cohort does not hold.

***S. aureus*, 803 genomes, 7 datasets.** One study, PRJNA1050907, holds 74.5% to 100% of the resistant genomes on every dataset. Levofloxacin has no lineage-split figure in its report, so its rows are not read. For the set covering machine, sensitivity reads below on 1 of the 6 datasets read, and oxacillin has 98.9% sensitivity against 90.3%. The known-gene rule gives every scored genome the same call on 4 of the 7 datasets. On tetracycline it calls all 731 genomes resistant.

***S. pneumoniae*, 503 genomes, 7 datasets.** Sensitivity never reads below for a learner. Specificity reads below on 18 of the 21 learner rows. On chloramphenicol 124 of the set covering machine's 125 false positives are in two of the three studies. It calls 1 of 262 susceptible genomes resistant in PRJNA320861, 79 of 95 in PRJEB82317 and 45 of 55 in PRJNA1357085.

On every *S. pneumoniae* dataset, more than half of the resistant genomes come from one study and more than half of the susceptible genomes from another. There sensitivity and specificity are measured mostly on different studies.

***E. coli*, 601 genomes, 18 datasets.** The set covering machine holds on 1 of the 18 datasets, the decision tree on 0, and the lasso on 0 of the 3 it is read on. The set covering machine reads below on sensitivity on 10 datasets and on specificity on 13. Its amikacin sensitivity is 9.4% against 46.7%, and its ceftazidime sensitivity is 22.9% against 87.2%. On sulfamethoxazole/trimethoprim its specificity is 39.3% against 83.5%.

Two facts limit this batch. One study, PRJNA841811, holds 71.6% to 90.6% of the resistant genomes on every dataset. The species check marks 328 of the 601 genomes as borderline, and the plan keeps and flags a borderline genome, so they are scored.

***M. tuberculosis*, 1,152 genomes, 11 datasets.** For the set covering machine, sensitivity reads below on 1 of the 11 datasets and specificity on 3. Ethambutol has 31.0% specificity against 60.9%. Pyrazinamide goes the other way, with 21.7% sensitivity against 58.0%. On levofloxacin both figures read above, with 94.3% sensitivity against 59.8% and 94.0% specificity against 79.2%. One study, PRJEB48275, holds at least half of the resistant genomes on 9 of the 11 datasets.

***Salmonella enterica*, 2,107 genomes, 17 datasets.** This batch has the most studies, up to 39 on a dataset. On two drugs the set covering machine reads below where the known-gene rule does not. On sulfamethoxazole the set covering machine has 30.6% sensitivity against 87.0%, and the known-gene rule has 98.9%. On tetracycline the set covering machine's specificity is 43.0% against 96.1%, and the known-gene rule's is 98.5%.

The other learners split. The lasso holds on sulfamethoxazole, with 93.0% sensitivity, and the decision tree reads below there, with 11.4% against 82.6%. On tetracycline the decision tree holds, with 96.6% sensitivity and 99.0% specificity.

***K. pneumoniae*, 2,014 genomes, 21 datasets.** No single study holds this batch: PRJEB29742 holds the most resistant genomes on 11 of the 21 datasets, with 24.0% to 90.6% of them. The set covering machine reads below on sensitivity on 11 datasets: meropenem has 61.3% against 85.6%, while the known-gene rule has 90.5%.

The decision tree reads below on specificity on 15 datasets, and on ertapenem its specificity is 12.6% against 89.0%. On tigecycline the set covering machine's specificity is 3.7% against 66.4%, with sensitivity at 97.7% against 23.7%.

***S. flexneri*, 78 genomes, 1 dataset.** Every genome is from one study. For the set covering machine, sensitivity is consistent with the report's figure and specificity reads below it, on 31 resistant and 47 susceptible genomes. The decision tree holds.

## Same-study population

These genomes are new, and they come from studies the cohort already holds. The two populations are scored apart and are never added together.

Three sentences below set a same-study figure beside the primary figure of the same model. The plan does not compare the two populations, so that comparison was added after both results were read.

On every dataset here, one study holds at least 99.9% of the resistant genomes. *E. cloacae*, *E. coli* and *K. pneumoniae* share one study, PRJNA288601, so their 32 datasets are one study's.

***E. cloacae*, 71 genomes, 2 datasets.** This is the smallest batch, and the species check marks 25 of its 71 genomes as borderline. The set covering machine reads below on sensitivity on both datasets. On sulfamethoxazole/trimethoprim its sensitivity is 59.3% against 100.0%.

***E. faecium*, 368 genomes, 5 datasets.** The study is PRJNA771404. The set covering machine holds on 3 of the 5 datasets. Its gentamicin sensitivity is 13.3% against 84.6%, on 30 resistant genomes.

***E. coli*, 87 genomes, 13 datasets.** The set covering machine holds on 11 of the 13 datasets. Its ceftazidime sensitivity is 80.6% against 87.2%, where the primary population gave 22.9%. A dataset here has 28 to 56 resistant genomes, so its intervals are wide.

***K. pneumoniae*, 459 genomes scored of 470, 17 datasets.** The 11 genomes that failed the species check are each nearer the *K. quasipneumoniae* reference. The set covering machine holds on 14 of the 17 datasets. Its meropenem sensitivity is 92.9% against 85.6%, where the primary population gave 61.3%.

***Salmonella enterica*, 1,741 genomes, 8 datasets.** One study, PRJNA248792, holds 99.9% to 100% of the resistant genomes on every dataset. The set covering machine holds on 5 of the 8 datasets, and two of its other rows read below with a sensitivity under 5%.

On amoxicillin/clavulanate it calls none of the 301 resistant genomes resistant, and the known-gene rule's sensitivity is 15.9% against 93.6%. On trimethoprim its sensitivity is 4.0% against 82.0%, where the decision tree has 99.7%.

***S. pneumoniae*, 418 genomes, 4 datasets.** The study is PRJEB3084. The set covering machine holds on 3 of the 4 datasets. Its penicillin specificity is 97.9% against 92.4%, where the primary population gave 76.0%. The lasso's penicillin sensitivity is 29.8% against 95.4%.

## Rows

`scripts/outside-study-results/tables.py --species` prints these tables from the output files in `reference/outside-study/`. Each figure is a cell of those files or a cell shown as a percentage.

<!-- tables:begin -->
**Rows cited in the text**

Each estimate is followed by its 95% interval and its counts. The report's figure is the lineage-split figure for a learner, and the whole-cohort figure for the known-gene rule.

| Population | Species | Drug | Learner | Sensitivity | Report | Reading | Specificity | Report | Reading |
|---|---|---|---|---|---|---|---|---|---|
| primary | *S. aureus* | oxacillin | Set covering machine | 98.9% (97.8% to 99.6%), 652 of 659 | 90.3% | above | 94.6% (85.1% to 98.9%), 53 of 56 | 98.3% | consistent |
| primary | *S. aureus* | tetracycline | Known-gene rule | 100.0% (96.5% to 100.0%), 104 of 104 | 100.0% | consistent | 0.0% (0.0% to 0.6%), 0 of 627 | 0.0% | consistent |
| primary | *S. pneumoniae* | chloramphenicol | Set covering machine | 96.7% (90.6% to 99.3%), 87 of 90 | 68.7% | above | 69.7% (65.0% to 74.1%), 287 of 412 | 97.7% | below |
| primary | *S. pneumoniae* | penicillin | Set covering machine | 97.7% (88.0% to 99.9%), 43 of 44 | 95.9% | consistent | 76.0% (70.8% to 80.6%), 234 of 308 | 92.4% | below |
| primary | *E. coli* | amikacin | Set covering machine | 9.4% (3.5% to 19.3%), 6 of 64 | 46.7% | below | 93.9% (91.6% to 95.8%), 496 of 528 | 98.4% | below |
| primary | *E. coli* | ceftazidime | Set covering machine | 22.9% (17.0% to 29.8%), 41 of 179 | 87.2% | below | 94.8% (92.0% to 96.9%), 331 of 349 | 92.1% | consistent |
| primary | *E. coli* | sulfamethoxazole/trimethoprim | Set covering machine | 93.4% (90.8% to 95.5%), 428 of 458 | 89.8% | above | 39.3% (30.4% to 48.8%), 46 of 117 | 83.5% | below |
| primary | *M. tuberculosis* | ethambutol | Set covering machine | 86.0% (81.9% to 89.5%), 296 of 344 | 79.6% | above | 31.0% (26.7% to 35.5%), 138 of 445 | 60.9% | below |
| primary | *M. tuberculosis* | pyrazinamide | Set covering machine | 21.7% (18.0% to 25.9%), 96 of 442 | 58.0% | below | 68.4% (63.6% to 72.9%), 279 of 408 | 55.9% | above |
| primary | *M. tuberculosis* | levofloxacin | Set covering machine | 94.3% (88.5% to 97.7%), 115 of 122 | 59.8% | above | 94.0% (91.8% to 95.7%), 606 of 645 | 79.2% | above |
| primary | *S. enterica* | sulfamethoxazole | Set covering machine | 30.6% (25.2% to 36.5%), 83 of 271 | 87.0% | below | 97.5% (92.9% to 99.5%), 117 of 120 | 93.4% | consistent |
| primary | *S. enterica* | sulfamethoxazole | Decision tree | 11.4% (7.9% to 15.8%), 31 of 271 | 82.6% | below | 98.3% (94.1% to 99.8%), 118 of 120 | 95.6% | consistent |
| primary | *S. enterica* | sulfamethoxazole | Lasso | 93.0% (89.3% to 95.7%), 252 of 271 | 75.8% | above | 98.3% (94.1% to 99.8%), 118 of 120 | 88.2% | above |
| primary | *S. enterica* | sulfamethoxazole | Known-gene rule | 98.9% (96.8% to 99.8%), 268 of 271 | 93.5% | above | 100.0% (97.0% to 100.0%), 120 of 120 | 99.5% | consistent |
| primary | *S. enterica* | tetracycline | Set covering machine | 98.9% (97.4% to 99.6%), 435 of 440 | 77.0% | above | 43.0% (38.0% to 48.1%), 169 of 393 | 96.1% | below |
| primary | *S. enterica* | tetracycline | Decision tree | 96.6% (94.4% to 98.1%), 425 of 440 | 83.7% | above | 99.0% (97.4% to 99.7%), 389 of 393 | 93.3% | above |
| primary | *S. enterica* | tetracycline | Known-gene rule | 98.9% (97.4% to 99.6%), 435 of 440 | 92.8% | above | 98.5% (96.7% to 99.4%), 387 of 393 | 91.7% | above |
| primary | *K. pneumoniae* | meropenem | Set covering machine | 61.3% (58.3% to 64.2%), 658 of 1,074 | 85.6% | below | 91.5% (88.9% to 93.7%), 517 of 565 | 95.7% | below |
| primary | *K. pneumoniae* | meropenem | Known-gene rule | 90.5% (88.6% to 92.2%), 972 of 1,074 | 95.0% | below | 88.8% (86.0% to 91.3%), 502 of 565 | 91.5% | below |
| primary | *K. pneumoniae* | ertapenem | Decision tree | 93.5% (91.2% to 95.4%), 533 of 570 | 82.9% | above | 12.6% (9.6% to 16.1%), 54 of 428 | 89.0% | below |
| primary | *K. pneumoniae* | tigecycline | Set covering machine | 97.7% (94.7% to 99.3%), 213 of 218 | 23.7% | above | 3.7% (2.3% to 5.8%), 19 of 509 | 66.4% | below |
| primary | *S. flexneri* | ciprofloxacin | Set covering machine | 74.2% (55.4% to 88.1%), 23 of 31 | 69.7% | consistent | 83.0% (69.2% to 92.4%), 39 of 47 | 95.7% | below |
| same-study | *K. pneumoniae* | meropenem | Set covering machine | 92.9% (88.9% to 95.8%), 223 of 240 | 85.6% | above | 97.6% (94.4% to 99.2%), 201 of 206 | 95.7% | consistent |
| same-study | *E. coli* | ceftazidime | Set covering machine | 80.6% (64.0% to 91.8%), 29 of 36 | 87.2% | consistent | 95.8% (85.7% to 99.5%), 46 of 48 | 92.1% | consistent |
| same-study | *S. pneumoniae* | penicillin | Set covering machine | 98.8% (95.8% to 99.9%), 169 of 171 | 95.9% | consistent | 97.9% (95.2% to 99.3%), 234 of 239 | 92.4% | above |
| same-study | *S. pneumoniae* | penicillin | Lasso | 29.8% (23.1% to 37.3%), 51 of 171 | 95.4% | below | 99.2% (97.0% to 99.9%), 237 of 239 | 98.3% | consistent |
| same-study | *E. faecium* | gentamicin | Set covering machine | 13.3% (3.8% to 30.7%), 4 of 30 | 84.6% | below | 96.6% (94.0% to 98.3%), 313 of 324 | 89.1% | above |
| same-study | *E. cloacae* | sulfamethoxazole/trimethoprim | Set covering machine | 59.3% (38.8% to 77.6%), 16 of 27 | 100.0% | below | 90.9% (78.3% to 97.5%), 40 of 44 | 94.9% | consistent |
| same-study | *S. enterica* | amoxicillin/clavulanate | Set covering machine | 0.0% (0.0% to 1.2%), 0 of 301 | 91.7% | below | 100.0% (99.5% to 100.0%), 716 of 716 | 99.7% | consistent |
| same-study | *S. enterica* | amoxicillin/clavulanate | Known-gene rule | 15.9% (12.0% to 20.6%), 48 of 301 | 93.6% | below | 100.0% (99.5% to 100.0%), 716 of 716 | 91.2% | above |
| same-study | *S. enterica* | trimethoprim | Set covering machine | 4.0% (2.2% to 6.8%), 13 of 322 | 82.0% | below | 100.0% (99.7% to 100.0%), 1,382 of 1,382 | 78.0% | above |
| same-study | *S. enterica* | trimethoprim | Decision tree | 99.7% (98.3% to 100.0%), 321 of 322 | 87.9% | above | 99.9% (99.6% to 100.0%), 1,381 of 1,382 | 66.7% | above |

**Cited rows, study by study**

The same counts as the row above them in `metrics.tsv`, split by the genome's BioProject.

| Population | Species | Drug | Learner | Study | Resistant genomes | Of them called resistant | Susceptible genomes | Of them called resistant |
|---|---|---|---|---|---:|---:|---:|---:|
| primary | *S. pneumoniae* | chloramphenicol | Set covering machine | PRJEB82317 | 42 | 42 | 95 | 79 |
| primary | *S. pneumoniae* | chloramphenicol | Set covering machine | PRJNA1357085 | 48 | 45 | 55 | 45 |
| primary | *S. pneumoniae* | chloramphenicol | Set covering machine | PRJNA320861 | 0 | 0 | 262 | 1 |
| primary | *S. pneumoniae* | penicillin | Set covering machine | PRJEB82317 | 44 | 43 | 36 | 3 |
| primary | *S. pneumoniae* | penicillin | Set covering machine | PRJNA1357085 | 0 | 0 | 87 | 71 |
| primary | *S. pneumoniae* | penicillin | Set covering machine | PRJNA320861 | 0 | 0 | 185 | 0 |
<!-- tables:end -->

```sh
python3 -I scripts/outside-study-results/tables.py reference/outside-study --species --check docs/outside-study-results-by-species.md
```

The command exits 1 unless the tables above are exactly what the files give.
