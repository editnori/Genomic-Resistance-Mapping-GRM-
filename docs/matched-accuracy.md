# Matched accuracy

Not run yet. Stage 2's real reading has not run and is not authorised: it needs the merged script and the lead's word after the merge.

Stage 1 freezes which genome rows of revision 14's lineage split enter a matched comparison, opening no score.

## Question

The protocol asks how often the five learners' mean score ranks a resistant genome above a susceptible one in a group matched on dataset, test fold, lineage cluster, BioProject, lab method, standard and standard year. It is a retrospective read inside one database, BV-BRC.

## Steps

1. Check every input's sha256 first. Each report must match the release's checksum list.
2. Load each report through named keys, keeping only the count of each learner's scores.
3. Stop on a line count other than 1,285,111, or a lab record or genome table row with any insertion time blank, not a real date and time, or at or after 19:10:00 on 5 November 2026.
4. Stop on other than five folds, a lineage cluster or copy group in two folds, or a lab call other than the label table's among rows not labelled from an MIC.
5. Take each row out under the first reason that applies: a copy group, no single BioProject, no lab record, a blank method or standard, or more than one tuple of method, standard and year.
6. Group the rest. A group with one lab outcome leaves.
7. Scan the counts for every genome id and accession read, after a positive control, then write them.

A failed check prints one line, and a scan hit one line a file. Both exit 1 and write nothing, and a write that fails removes what it made.

## Readings

Reading A uses these rows. Reading B first takes out each genome and antibiotic pair with an intermediate lab record, then groups again. Row lists stay local. Only their sha256 and line counts enter the counts.

## Counts

Per dataset, species and overall: rows out by each reason, then the rows, genomes, groups and pairs left, and the same for exact MLST and MLST with HC20, beside the wider match on those rows. Overall alone: the largest BioProject's share of rows and of weight. Groups weigh equally inside a dataset, datasets inside a species and species overall.

Rows lacking a typing value are counted among rows entering the groups and among eligible rows; genomes missing from the genome table, over every report and over those that enter.

## Rules

- A copy group is a group the report lists: genomes of one dataset with identical canonical k-mer sets. The script compares no assembly.
- A BioProject is one part, PRJ, two letters and digits in any case. A genome whose one part has another shape has none and is counted apart. Blank genome table rows are skipped.
- The 5 November genome table gives the HC20 group and has no insertion time to check. Without it, no HC20 count is made.

## Stage one run

Stage 1 ran once on 27 November 2026 at main 00755885, ended with exit 0 and printed "scores opened: no". Ruling R59 gives these counts.

| Count | Reading A | Reading B |
| --- | --- | --- |
| Eligible rows | 49,621 | 48,156 |
| Genomes | 9,132 | 9,122 |
| Matched groups | 2,783 | 2,749 |
| Resistant-susceptible pairs | 283,118 | 272,840 |
| Datasets | 206 | 205 |
| Species | 16 | 16 |

Coverage, which every later figure carries: 49,621 of 215,319 genome rows and 16 of 23 species in reading A, and 48,156 of 215,319 rows and 16 of 23 species in reading B. Each row count of stage 2's export is within those 215,319 rows, and its inputs table gives both wholes.

- The BioProject with the most eligible rows: 6,308 of 49,621 rows in reading A, and 2771/48384 of the weight.
- The largest weight of any BioProject: 1/8 in both readings.

## Stage two steps

Before step 1, the run stops on an output path that is not a directory or already holds an export, a `local` directory or the line "scores opened: yes" of an earlier run, and on numpy other than the version the resampling is written for.

1. Check the sha256 of every report, both row lists, stage 1's committed count and typing tables and the 5 November genome table against values no argument changes. Stop on a missing file or a report whose id is not its file name.
2. Stop on a line repeated, out of stage 1's sort or giving a genome twice, a listed dataset or row no report holds, a fold or lineage cluster unlike its report's, or a report without scm, cart, rf, gbt and lasso in order.
3. Rebuild each list's groups. Stop on a group with one lab outcome, or on counts or typing support that differ from stage 1's tables at any level.
4. Stop on a report changed since step 1, a listed row moved, a catalogue call not true, false or null, fit sizes by fold not 5 finite numbers, or a dataset or species name with a full stop. Their loader keeps no value.
5. Write "scores opened: yes" to a file and print it. Only then does a second loader read, for listed rows, the five scores, each learner's fit size by fold and the catalogue call, and whether the report's top-level catalogue is present.
6. Measure both readings as below, scan the export for every genome id and accession read, after a positive control, and write it in integers, estimates in parts per million.

A file not as stage 1 writes it stops its step too. A committed table must have stage 1's columns, whole lines and no line twice. A row list must be UTF-8, hold a line, end with a line break and give 8 fields a line, none blank but the year, with whole-number folds and lineage clusters. The reports directory must hold the reports the release's checksum list names, and each report and genome table row must pass stage 1's loader.

A check that fails before the line prints one line, exits 1 and writes nothing. So does a line file that cannot be made or written. The run first removes a line file it made, so a rerun into that directory starts, and its stop line says when the file cannot be removed.

After the line, a run exits 0 with the export, or stops with one line and no export on one of four ends:

- A listed report that cannot be read, or whose bytes changed since step 1.
- A score null, true or false, not a number or not finite.
- A full stop in an export file, a failed positive control of the scan, or a scan hit, which prints one line a file.
- An export write that fails, which first removes what it made.

## Lead's readings

- A row's mean is its five scores added left to right as doubles, scm, cart, rf, gbt, lasso, then divided by 5. A score null, true or false, not a number or not finite stops the run.
- In a group each resistant row meets each susceptible row: a higher resistant mean counts 1, an exact tie one half, else 0. The group's value is that count over its pairs.
- Groups weigh equally in a dataset, datasets in a species, species overall, as exact fractions. With 6 or more species, 10,000 draws of species means, numpy's PCG64 at seed 42 for each estimate and the linear 2.5th and 97.5th percentiles give the bounds.
- Secondary reads keep the weights: four cells of catalogue calls, resistant member first; MLST, then MLST and HC20, each beside the wider match on the same rows; and the read without groups in a fold where a learner fit nothing.
- Unevaluable calls and reports without a drug rule are counted apart. Each cell is weighted over its own groups, datasets and species. A cell prints its estimate wherever it has pairs and gets no interval.
- A dataset or species under 10 pairs enters every mean. Its row gives counts and no estimate.

A row under 10 pairs carries no reading, yet its value still follows by arithmetic, to within the rounding, from the mean above it and the other rows, when that mean prints and no other row beneath it lacks an estimate.

## Ruling conditions

- The real reading runs once at the merged commit, after the lead's word, on the computer holding the 274 reports and both row lists. A run that stops before "scores opened: yes" may run again. One that printed it is the one reading.
- Both readings are reported, whatever their figures, each with its coverage.
- Each interval of the export carries one word from its bounds before rounding: `above` with the lower bound above one half, `below` with the upper bound below it, else `unresolved`. A claim above chance needs `above` in both readings.

## Limits

- Stage 1 never measures the 252,148 records dropped at the cut, a recorded figure.
- Species weigh equally, so the resampling of species means carries a thin species into the interval's width.
- No output names or ranks a feature. The export holds no row per genome, pair, group or BioProject.

Record: protocol `reference/matched-accuracy/PROTOCOL-1.0.md` sha256 e4c124466f7b85c649ea4684004831beb80002174e51317ad0ee7476a84401b0; amendment `reference/matched-accuracy/AMENDMENT-1.md` sha256 634ffb8db8587a34eba879465a414a518d4d1afdaf730bca69b87b5c782fe27e; scripts `scripts/matched-accuracy/stage1.py` and `stage2.py`; faults `scripts/matched-accuracy/planted_faults.py` and `stage2_faults.py`; release rev14-274-20261112 with its checksum list `reference/lineage-control/rev14-274-20261112.reports.sha256`; stage 1's export `reference/matched-accuracy/stage1-rev14-274-20261112/`.
