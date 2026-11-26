# Matched accuracy

Not run yet. Stage 1 counts which genome rows of revision 14's lineage split can enter a matched comparison and freezes that list; it opens no score, and stage 2 is not authorised.

## Question

The protocol asks how often the mean score of the five learners ranks a resistant genome above a susceptible one inside a matched group: the same dataset, test fold, lineage cluster, BioProject, lab method, standard and standard year. It is a retrospective read inside one database, BV-BRC, and never external validation.

## Steps

1. Check the sha256 of every input before any other read. Each report must equal its line of the release's own checksum list.
2. Load each report through a named list of keys. Of each learner's scores, the loader keeps the count alone.
3. Stop on a lab record inserted at or after 19:10:00 on 5 November 2026, or a line count other than 1,285,111.
4. Stop on a fold count other than five, a lineage cluster or copy group in two folds, or a lab call other than the label table's.
5. Take each row out under the first reason that applies: a copy group, no single BioProject, no lab record, a blank method or standard, or more than one tuple of method, standard and year.
6. Group the rest. A group with one lab outcome leaves.
7. Scan the counts for every genome id and BioProject accession read, after a positive control, then write them.

Each stop prints one line, writes nothing and exits 1.

## Readings

Reading A uses the rows above. Reading B first takes out every genome and antibiotic pair with an intermediate lab record, then forms the groups again. Both are always computed and written. Each frozen row list stays local, and only its sha256 enters the counts.

## Counts

Per dataset, per species and overall: rows out by each reason, rows, genomes, groups, pairs of one resistant and one susceptible row, datasets and species. Then the largest BioProject's share of rows and of weight, as whole numbers, and the same counts for exact MLST and for MLST with the HC20 group, beside the wider match on the same rows. Groups weigh equally inside a dataset, datasets inside a species and species overall.

## Limits

- Stage 2, the comparison of scores, is not authorised. A ruling on the counts of stage 1 comes first.
- The read is retrospective, inside one database, and is never external validation.
- The 252,148 records dropped at the cut are a recorded figure; stage 1 does not measure them.
- A copy group stands for samples duplicated within a dataset, since no id repeats inside a report.
- No output names a feature or ranks one.
- The tests and planted faults use made-up files alone.

Record: protocol `reference/matched-accuracy/PROTOCOL-1.0.md` sha256 e4c124466f7b85c649ea4684004831beb80002174e51317ad0ee7476a84401b0; amendment `reference/matched-accuracy/AMENDMENT-1.md` sha256 634ffb8db8587a34eba879465a414a518d4d1afdaf730bca69b87b5c782fe27e; script `scripts/matched-accuracy/stage1.py`; faults `scripts/matched-accuracy/planted_faults.py`; release rev14-274-20261112 with its checksum list `reference/lineage-control/rev14-274-20261112.reports.sha256`.
