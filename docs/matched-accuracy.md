# Matched accuracy

Not run yet. Stage 1 counts which genome rows of revision 14's lineage split can enter a matched comparison and freezes that list; it opens no score, and stage 2 is not authorised.

## Question

The protocol asks how often the five learners' mean score ranks a resistant genome above a susceptible one in a group matched on dataset, test fold, lineage cluster, BioProject, lab method, standard and standard year. It is a retrospective read inside one database, BV-BRC, and never external validation.

## Steps

1. Check every input's sha256 first; each report must match the release's checksum list.
2. Load each report through named keys, keeping only the count of each learner's scores.
3. Stop on a line count other than 1,285,111, or a lab record or genome table row with any insertion time blank, not a real date and time, or at or after 19:10:00 on 5 November 2026.
4. Stop on other than five folds, a lineage cluster or copy group in two folds, or a lab call other than the label table's among rows not labelled from an MIC.
5. Take each row out under the first reason that applies: a copy group, no single BioProject, no lab record, a blank method or standard, or more than one tuple of method, standard and year.
6. Group the rest; a group with one lab outcome leaves.
7. Scan the counts for every genome id and accession read, after a positive control, then write them.

A failed check prints one line, a scan hit one line a file; both exit 1 and write nothing, and a write that fails removes what it made.

## Readings

Reading A uses these rows. Reading B first takes out each genome and antibiotic pair with an intermediate lab record, then groups again. Row lists stay local; only their sha256 and line counts enter the counts.

## Counts

Per dataset, species and overall: rows out by each reason, then the rows, genomes, groups and pairs left, the largest BioProject's share of rows and of weight, and the same for exact MLST and MLST with HC20, beside the wider match on those rows. Groups weigh equally inside a dataset, datasets inside a species and species overall.

Rows lacking a typing value are counted among rows entering the groups and among eligible rows; genomes missing from the genome table, over every report and over those that enter.

## Rules

- A copy group is a group the report lists: genomes of one dataset with identical canonical k-mer sets. The script compares no assembly.
- A BioProject is one part, PRJ, two letters and digits in any case; a genome whose one part has another shape has none and is counted apart. Blank genome table rows are skipped.
- The 5 November genome table gives the HC20 group and has no insertion time to check; without it, no HC20 count is made.

## Limits

- Stage 2, the comparison of scores, is not authorised; a ruling on stage 1's counts comes first.
- The read is retrospective, inside one database, and never external validation.
- Stage 1 never measures the 252,148 records dropped at the cut, a recorded figure.
- No output names or ranks a feature.
- Beside made-up files and the repository's own pages and code, the tests read 3 committed files: the protocol, the amendment and the release's checksum list, on which one test runs the script.

Record: protocol `reference/matched-accuracy/PROTOCOL-1.0.md` sha256 e4c124466f7b85c649ea4684004831beb80002174e51317ad0ee7476a84401b0; amendment `reference/matched-accuracy/AMENDMENT-1.md` sha256 634ffb8db8587a34eba879465a414a518d4d1afdaf730bca69b87b5c782fe27e; script `scripts/matched-accuracy/stage1.py`; faults `scripts/matched-accuracy/planted_faults.py`; release rev14-274-20261112 with its checksum list `reference/lineage-control/rev14-274-20261112.reports.sha256`.
