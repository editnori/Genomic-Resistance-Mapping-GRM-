# Givens export

[Givens](https://github.com/editnori/givens) is a provisional shared browser of the definitions and assumptions of four projects. GRM writes its own part to one file, `reference/givens-export/grm.json`, with `scripts/givens-export/export.py`.

The file covers evaluation terms, the assumptions a report figure rests on and aggregate counts. It changes nothing in GRM's engine, and GRM reads nothing from Givens. [givens-terms.md](givens-terms.md) shows the same entries as a page to read.

## Contents

It holds 22 definitions, 24 links between them and 3 assumptions. Its 37 measured lines split into 29 read from a file and 8 quoted from a document. It marks 1 entry with what GRM has not counted.

The three assumptions:

- A genome in a test fold has no copy in the training folds.
- A cross-validated figure holds for genomes from a newly deposited study the cohort never saw.
- Every learner row of a dataset is scored on the same genomes.

## Fields

The file has the shape of a Givens module, so Givens can store it unchanged.

| Field | Holds |
|---|---|
| `nodes` | one definition each: `id`, `kind`, `label`, `definition`, and `quotes`, the sentences of this repository that back it |
| `edges` | one link each: `from`, `type`, `to` |
| `assumptions` | `statement`, `under`, `hidden_in` and `when_false`, with `quotes` |
| `measured` | on a definition or an assumption, the counts that bear on it |
| `unmeasured` | on a definition or an assumption, what GRM has not counted and why |
| `sources` | the script, its spec, each file read with its sha256, and each document quoted |
| `import_rule` | the rules for whoever shows the file |
| `left_out` | what the file never holds |

A measured line has `count`, `of`, `what`, `source` and `basis`. It also names its `evaluation_revision`, so a figure of revision 12 is never shown as revision 14's. A median has no count. Its line carries `values`, one for each revision, and a `note`.

## Sources

Each number has one of two bases.

- **File.** The script reads the number from a committed summary file. `read` says which field, and `sources` lists the file's sha256.
- **Document.** The number is quoted from a sentence of this repository's text, and `quote` holds the sentence. The script requires the sentence to be in the document and the number to be in the sentence.

The words come from `scripts/givens-export/spec.json`. A number written in the spec stops the script.

The files read:

| File | Read for |
|---|---|
| `reference/release-summary/rev12-274-20261105/index-counts.json` | learner rows that left genomes without a call |
| `reference/release-summary/rev12-274-20261105/guidance.json` | those rows, by learner |
| `reference/release-summary/rev12-274-20261105/summary.json` | two learners' median balanced accuracy on the lineage split |
| `reference/release-summary/rev14-274-20261112/index-counts.json` | learner rows, the rows that carry the flag, the error-rate counts compared, the datasets that hold labels inferred from a MIC |
| `reference/release-summary/rev14-274-20261112/guidance.json` | datasets that have each split |
| `reference/release-summary/rev14-274-20261112/summary.json` | datasets, lineage splits, the same two medians |
| `web/src/data/outside-study.json` | the outside-study readings, against revision 12 reports |

`tests/outside_study_results` holds the last file to `reference/outside-study`, and `tests/release_summary` holds the other six to their README.

The lines quoted from a document:

| Line | Document | A suite rebuilds it |
|---|---|---|
| 274 of 277 curated datasets have a report, and 3 of 277 have none | `README.md` | no |
| 58 of 277 datasets had a duplicate group at revision 9, and 55 of 277 had one that crossed the saved random folds | `docs/duplicate-isolation.md` | no |
| 73 of 274 reports differ from revision 12 in the rows of a learner with a null fit, and the best learner changed in 21 of 274 | `README.md` | no |
| the lasso holds on 8 of 25 and on 7 of 15 datasets against the revision 14 figure | `docs/outside-study-results.md` | `tests/outside_study_results` |

## Missing

The export marks what GRM has not measured.

- **Duplicate groups after revision 9.** The two counts of duplicate groups are of the revision 9 cohorts. No committed file counts them again on a later revision.
- **Six document lines.** The first three rows of the table above have no committed file behind them. The repository states them, and no test can rebuild them.
- **Three empty report fields.** A report field can be empty because a measure does not apply to a learner, because a denominator is zero, or because an older report lacks the field. GRM keeps no count of any of the three, and the export has no entry for them.

## Left out

The file names no dataset, species or antibiotic. It also leaves out:

- rows for one species or one drug class;
- GRM's features, models, panels and samples;
- the known-gene rule's rows and its catalog;
- MIC values and the MIC regression's figures;
- sequences, k-mers, variants and genes;
- the error-limit screen's counts.

The entry on labels inferred from a MIC counts datasets. It reads the datasets where the released index records at least one such label, of the datasets whose entry records the count, from the block `labels_from_mic_excluded` of `index-counts.json` ([release-summary.md](release-summary.md)). It holds no MIC value, and a dataset whose entry lacks the count is never read as zero.

The script enforces this with two lists and one screen. It reads numbers from the seven files above only. In each it reads only the named count fields, and it counts the datasets that have no lineage split without reading their ids.

The screen stops the script if any string it would write names a dataset, a species, an antibiotic or a drug class. It takes those names from the seven files and from GRM's antibiotic vocabulary in `src/drugs.rs`: every antibiotic, each spelling and brand the vocabulary folds, and every drug class. It refuses a class word in the singular and in the plural. The script writes none of these names.

## Use

Run from the repository root, with Python 3.12 and the standard library only:

```sh
python3 -I scripts/givens-export/export.py --write
python3 -I scripts/givens-export/export.py --check
```

The first command writes `reference/givens-export/grm.json`, `SHA256SUMS` beside it and the page `docs/givens-terms.md`. The list holds the file's sha256 on one line, as `sha256sum` prints it, so Givens can store it unchanged. The second command prints nothing and exits 1 unless all three are exactly what the script writes. Run the first again after a change to a quoted sentence, to a summary file or to the spec.

## Page

`docs/givens-terms.md` holds the file's terms, links, assumptions and counts in the file's order. The script builds it from the export alone. Its own sentences hold no number beyond the totals of the first paragraph, and it passes the same screen as the file.

The page leaves out the quoted sentences and each sha256, which stay in the file. The spec may not quote the page, because the script writes it.

## Checks

`tests/givens_export` runs in CI and fits nothing. It fails if:

- the committed file differs from what the script writes, or `SHA256SUMS` differs from the file;
- the page differs from what the script writes, lacks a term, an assumption, a row or a note of the file, or states a number the file lacks;
- the file breaks a rule that Givens' own check applies to a module;
- a count of the export differs from the sentence of `docs/outside-study-results.md` or `docs/release-summary.md` that states it;
- a line's revision differs from the revision its files name;
- an entry with no count lacks `unmeasured`.

It also plants faults in a copy of the spec, and each must stop the script: a species, an antibiotic, a dataset or a sequence in a string, a quote its document does not hold, a file or a field outside the lists, an unknown field, a count above its denominator, and a number written by hand. A copy of `index-counts.json` whose count of labels inferred from a MIC is missing or null must stop it too, so such a count is never written as zero.

The suite tries every name the screen has from the vocabulary alone, and each one stops the script. It also renames each table of a copy of the vocabulary. The script must then stop, because a table it cannot read would shrink the screen.

## Limits

- The file names no commit, because a committed file cannot name its own. Givens records the GRM commit it copied the file from.
- The export changes nothing in Givens by itself. Givens takes the file in with its own importer, which records the GRM commit.
- A quote shows that the repository states a sentence. It does not show that the sentence is right.
- The screen knows the names in the seven files and in GRM's antibiotic vocabulary. It knows no gene or variant, and it misses a name the vocabulary lacks, such as a new brand. A change to the spec's words still needs a read.
- The outside-study lines come from a retrospective test on newly deposited studies from the same database. It is no external validation. They read revision 12 reports, and the site serves revision 14.
- The export is a documentation and software change. It has no clinical effect.
