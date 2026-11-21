# TBProfiler comparison record

This directory holds the record of the comparison of TBProfiler 6.6.3, run with the database it ships on assemblies of the genomes of the revision 14 M. tuberculosis reports, with the phenotypes those reports record and with the five learners on the same genomes. The page that reads it is [`docs/tbprofiler-comparison.md`](../../docs/tbprofiler-comparison.md). `scripts/tbprofiler-comparison/record.py` writes every file here from the input set, and its `--check` mode fails on any byte it would not write. The input set is the call file of the TBProfiler run on assemblies of public BV-BRC genomes, with GRM's revision 14 M. tuberculosis reports.

## Files

| File | Rows | Contents |
|---|---|---|
| `reports.tsv` | 22 | Row file. One row per M. tuberculosis report of the release |
| `callers.tsv` | 880 | Row file. One row per way of reading, report and caller: the counts on the genomes with a call |
| `intervals.tsv` | 946 | Row file. The lineage interval's ends, for the primary way of reading |
| `disagreements.tsv` | 20 | Row file. One row per report: the primary way's counts by outcome and grade |
| `grades.tsv` | 18 | Row file. One row per tool drug: the shipped database's entries by source and grade |
| `order.tsv` | 5 | Row file. One row per event of the plan and the scoring run |
| `settings.tsv` | 14 | Row file. The settings the script applies and the sha256 of the plan, the addendum and the compact file |
| `totals.tsv` | 15 | Row file. Counts over all reports, and the catalog row's figures over the reports of the plan's set that have one |
| `set-aside.tsv` | 11 | Row file. Exploratory: the counts both ways for the one entry the tool's rule set aside |
| `checks.tsv` | 4 | The checks the build ran on the input set |
| `inputs.sha256` | 8 | The sha256 of each file of the input set |
| `measures.tsv` | 880 | Derived. One row per way of reading, report and caller |
| `endpoints.tsv` | 80 | Derived. One row per way of reading and report |
| `medians.tsv` | 264 | Derived. One row per way of reading, set, caller and measure |
| `comparisons.tsv` | 80 | Derived. One row per way of reading, set and learner row: each learner minus the tool |
| `screen.tsv` | 8 | Derived. One row per way of reading and set: the 0.02 screen |
| `bounds.tsv` | 32 | Derived. Exploratory: one row per way of reading, set and bound on every genome |
| `summary.tsv` | 29 | Derived. Single figures the page quotes |

`SHA256SUMS` gives the sha256 of every other file, and this README is what the script writes.

## Columns

Each line of `inputs.sha256` gives the sha256 of a file of the input set and its name. The other files' columns, in file order:

| File | Column | Gives |
|---|---|---|
| `reports.tsv` | `dataset` | the report, as species and antibiotic |
| `reports.tsv` | `antibiotic` | the report's antibiotic |
| `reports.tsv` | `kind` | primary (the tool names the drug), proxy (read from another drug by a shared name stem) or not named |
| `reports.tsv` | `tool_drugs` | the tool drug or drugs read, with `,` between them |
| `reports.tsv` | `graded` | yes when the shipped database grades at least one entry of the tool drug, no when it grades none, `-` for a report that is not primary |
| `reports.tsv` | `genomes` | genomes of the report |
| `reports.tsv` | `catalog_row` | yes when the report has a catalog row |
| `callers.tsv` | `way` | the way of reading the tool's output |
| `callers.tsv` | `dataset` | the report, as species and antibiotic |
| `callers.tsv` | `caller` | tool, or a learner and split |
| `callers.tsv` | `genomes` | genomes of the report |
| `callers.tsv` | `with_call` | genomes with a call of the tool |
| `callers.tsv` | `no_call_resistant` | genomes without a call, recorded resistant |
| `callers.tsv` | `no_call_susceptible` | genomes without a call, recorded susceptible |
| `callers.tsv` | `tp` | resistant genomes called resistant |
| `callers.tsv` | `fp` | susceptible genomes called resistant |
| `callers.tsv` | `tn` | susceptible genomes called susceptible |
| `callers.tsv` | `fn` | resistant genomes called susceptible |
| `callers.tsv` | `caller_only_right` | genomes this caller gets right and the tool does not |
| `callers.tsv` | `tool_only_right` | genomes the tool gets right and this caller does not |
| `intervals.tsv` | `kind` | endpoint, difference, median or median_difference |
| `intervals.tsv` | `set` | by_name (the plan's set), graded (the addendum's set) or `-` |
| `intervals.tsv` | `dataset` | the report, or `-` for a median |
| `intervals.tsv` | `caller` | tool, or a learner and split |
| `intervals.tsv` | `measure` | sensitivity, specificity, balanced_accuracy or minus_tool |
| `intervals.tsv` | `low` | the lower end of the lineage interval |
| `intervals.tsv` | `high` | the upper end of the lineage interval |
| `intervals.tsv` | `resamples_with_value` | resamples in which the figure has a value |
| `disagreements.tsv` | `dataset` | the report, as species and antibiotic |
| `disagreements.tsv` | `tool_r_label_s_assoc` | resistant call on grade Assoc w R, recorded susceptible |
| `disagreements.tsv` | `tool_r_label_s_interim` | resistant call on grade Assoc w R - Interim alone, recorded susceptible |
| `disagreements.tsv` | `tool_s_label_r_other` | susceptible call with an entry of another grade reported, recorded resistant |
| `disagreements.tsv` | `tool_s_label_r_none` | susceptible call with nothing reported for the drug, recorded resistant |
| `disagreements.tsv` | `both_r_assoc` | resistant call on grade Assoc w R, recorded resistant |
| `disagreements.tsv` | `both_r_interim` | resistant call on grade Assoc w R - Interim alone, recorded resistant |
| `disagreements.tsv` | `both_s_other` | susceptible call with an entry of another grade reported, recorded susceptible |
| `disagreements.tsv` | `both_s_none` | susceptible call with nothing reported for the drug, recorded susceptible |
| `disagreements.tsv` | `no_call_label_r` | no call, recorded resistant |
| `disagreements.tsv` | `no_call_label_s` | no call, recorded susceptible |
| `grades.tsv` | `tool_drug` | the tool drug |
| `grades.tsv` | `assoc_w_r` | entries of source WHO catalogue v2, grade Assoc w R |
| `grades.tsv` | `assoc_w_r_interim` | entries of source WHO catalogue v2, grade Assoc w R - Interim |
| `grades.tsv` | `own_uncertain_significance` | the database's own entries, grade Uncertain significance |
| `grades.tsv` | `own_not_assoc_w_r_interim` | the database's own entries, grade Not assoc w R - Interim |
| `grades.tsv` | `own_ungraded` | the database's own entries with no grade |
| `grades.tsv` | `counted_grade_entries` | entries of the two counted grades |
| `order.tsv` | `event` | the event |
| `order.tsv` | `time` | its UTC time |
| `order.tsv` | `source` | the file the time is read from |
| `settings.tsv` | `name` | the name of a setting |
| `settings.tsv` | `value` | the setting |
| `totals.tsv` | `name` | the name of a count over all reports or of a catalog row figure |
| `totals.tsv` | `value` | the figure |
| `set-aside.tsv` | `name` | the name of a count |
| `set-aside.tsv` | `value` | the count |
| `checks.tsv` | `check` | a check the build ran on the input set |
| `checks.tsv` | `count` | the count the check gives |
| `measures.tsv` | `way` | the way of reading the tool's output |
| `measures.tsv` | `dataset` | the report, as species and antibiotic |
| `measures.tsv` | `caller` | tool, or a learner and split |
| `measures.tsv` | `sensitivity` | `tp` over `tp` plus `fn` |
| `measures.tsv` | `specificity` | `tn` over `tn` plus `fp` |
| `measures.tsv` | `balanced_accuracy` | the mean of sensitivity and specificity |
| `measures.tsv` | `minus_tool` | this caller's balanced accuracy minus the tool's, `-` for the tool |
| `endpoints.tsv` | `way` | the way of reading the tool's output |
| `endpoints.tsv` | `dataset` | the report, as species and antibiotic |
| `endpoints.tsv` | `kind` | primary or proxy |
| `endpoints.tsv` | `scored` | yes when at least 25 resistant and 25 susceptible genomes have a call |
| `endpoints.tsv` | `by_name` | yes for a scored primary report |
| `endpoints.tsv` | `graded` | yes for a report of the graded set |
| `endpoints.tsv` | `tool_balanced_accuracy` | the tool's balanced accuracy |
| `endpoints.tsv` | `best_lineage_learner` | the learner whose lineage split row has the highest balanced accuracy |
| `endpoints.tsv` | `best_lineage_balanced_accuracy` | that learner's balanced accuracy |
| `endpoints.tsv` | `gap` | that balanced accuracy minus the tool's |
| `endpoints.tsv` | `within_margin` | yes when `gap` is at most 0.02 |
| `endpoints.tsv` | `sensitivity_lower` | sensitivity on every genome, no call counted as wrong |
| `endpoints.tsv` | `sensitivity_upper` | sensitivity on every genome, no call counted as right |
| `endpoints.tsv` | `specificity_lower` | specificity on every genome, no call counted as wrong |
| `endpoints.tsv` | `specificity_upper` | specificity on every genome, no call counted as right |
| `medians.tsv` | `way` | the way of reading the tool's output |
| `medians.tsv` | `set` | by_name (the plan's set), graded (the addendum's set) or `-` |
| `medians.tsv` | `caller` | tool, or a learner and split |
| `medians.tsv` | `measure` | the measure |
| `medians.tsv` | `endpoints` | reports in the set |
| `medians.tsv` | `median` | the median over the set's reports |
| `medians.tsv` | `min` | the smallest |
| `medians.tsv` | `max` | the largest |
| `medians.tsv` | `low` | the lineage interval's lower end, `-` for other ways of reading |
| `medians.tsv` | `high` | the lineage interval's upper end, `-` for other ways of reading |
| `comparisons.tsv` | `way` | the way of reading the tool's output |
| `comparisons.tsv` | `set` | by_name (the plan's set), graded (the addendum's set) or `-` |
| `comparisons.tsv` | `caller` | a learner and split |
| `comparisons.tsv` | `endpoints` | reports in the set |
| `comparisons.tsv` | `tool_higher` | reports where the tool's balanced accuracy is higher |
| `comparisons.tsv` | `tool_lower` | reports where it is lower |
| `comparisons.tsv` | `equal` | reports where the two are equal |
| `comparisons.tsv` | `tool_ahead_or_within` | reports where the learner minus the tool is at most 0.02 |
| `comparisons.tsv` | `median_minus_tool` | the median of the learner minus the tool |
| `comparisons.tsv` | `low` | the lineage interval's lower end, `-` for other ways of reading |
| `comparisons.tsv` | `high` | the lineage interval's upper end, `-` for other ways of reading |
| `comparisons.tsv` | `caller_only_right` | summed over the set |
| `comparisons.tsv` | `tool_only_right` | summed over the set |
| `screen.tsv` | `way` | the way of reading the tool's output |
| `screen.tsv` | `set` | by_name (the plan's set), graded (the addendum's set) or `-` |
| `screen.tsv` | `endpoints` | reports in the set |
| `screen.tsv` | `within_best` | reports with `within_margin` yes |
| `screen.tsv` | `needed` | half the set's reports, rounded up |
| `bounds.tsv` | `way` | the way of reading the tool's output |
| `bounds.tsv` | `set` | by_name (the plan's set), graded (the addendum's set) or `-` |
| `bounds.tsv` | `bound` | the bound on every genome |
| `bounds.tsv` | `endpoints` | reports in the set |
| `bounds.tsv` | `median` | the median over the set's reports |
| `bounds.tsv` | `min` | the smallest |
| `bounds.tsv` | `max` | the largest |
| `summary.tsv` | `name` | the name of a figure |
| `summary.tsv` | `value` | the figure |

## Rules

**Calls.** The four ways of reading are `who_assoc` (an entry graded Assoc w R or Assoc w R - Interim, the primary way), `who_assoc_only` (Assoc w R alone), `tool_default` (any entry the tool reports for the drug) and `who_assoc_small` (the primary way without large changes). A genome without a call is never read as susceptible: every caller is counted on the genomes where the tool has a call.

**Learners.** A learner calls resistant at a report score of 0.5 or more, on the lineage and the random split. The catalog row is GRM's own reader of the AMRFinderPlus catalog.

**Catalog row.** The catalog row appears only as six-report totals, in `totals.tsv` and again in `summary.tsv`, because on three reports a per-report count of it would describe one genome. Its median is given to 3 decimals for the same reason: in full, with the other totals, it would show whether one of those genomes had a call.

**Sets.** A primary report is scored when at least 25 resistant and 25 susceptible genomes have a call. `by_name` holds the scored primary reports, the plan's set; `graded` holds those whose tool drug the database grades, the addendum's set. Proxy reports enter no set.

**Interval.** 2,000 resamples of the lineage groups, the reports' lineage ids joined across all the reports and numbered in the text order of their smallest genome id. Each resample draws as many group numbers as there are groups with Python's `random.Random(20261114).randrange`, one call per draw, and takes every genome of each group drawn. A figure has no value in a resample where a class has no genome. `low` and `high` are the values at sorted index int(0.025 n) and n - 1 - int(0.025 n) of the n resamples with a value: the 51st and the 1,950th when n is 2,000. Medians skip a report without a value.

**Rounding.** Every figure here is the full value, apart from the catalog row's median. The page rounds each one once.

## Limits

No row is about one genome, and no cell names a genome. The tool's raw output and the assemblies no longer exist; the call file, kept outside the repository, is the lowest layer the counts are rebuilt from. `set-aside.tsv` and the counts behind it rest on the raw output alone and cannot be rebuilt. `--check` alone cannot hold a row file to the input set; that needs the input set.

## Guards

`--check` with the input set rebuilds every file and fails on any byte that differs. `--check` alone reads the row files, derives every derived file and this README again, and fails on any byte that differs or on a file that differs from `SHA256SUMS`. The input set is not in the repository: it holds the call file and the reports' scores, which are per genome.

`--out` writes only to a new path, an empty directory or an earlier record of this script whose files agree with its own `SHA256SUMS`, and refuses a symlink.

## Commands

```sh
python3 -I scripts/tbprofiler-comparison/record.py --inputs <dir> --out reference/tbprofiler-comparison
python3 -I scripts/tbprofiler-comparison/record.py --inputs <dir> --check reference/tbprofiler-comparison
python3 -I scripts/tbprofiler-comparison/record.py --check reference/tbprofiler-comparison
python3 -I -m unittest discover -s tests/tbprofiler_comparison -t .
```
