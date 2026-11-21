# AMRFinderPlus comparison record

This directory is the record behind [`docs/amrfinderplus-comparison.md`](../../docs/amrfinderplus-comparison.md): the AMRFinderPlus program, GRM's known-gene rule and GRM's five learners scored on the same genomes and recorded labels of public BV-BRC datasets. `scripts/amrfinderplus-comparison/record.py` writes every file of it, and its `--check` mode fails on any byte it would not write.

This is a retrospective comparison inside one public database, BV-BRC. It is not a prospective test and not external validation. Neither the AMRFinderPlus program nor GRM's known-gene rule reads a label, but the catalog both read is newer than most of these genomes and may have been curated from the same public isolates, which favours both against the learners.

## Files

| File | Rows | Contents |
|---|---|---|
| `README.md` |  | this file |
| `SHA256SUMS` |  | the sha256 of every other file |
| `callers.tsv` | 3775 | row file: counts and measures per dataset and caller, on the genomes and labels of the release |
| `datasets.tsv` | 253 | row file: one row per dataset with a published rule |
| `differences.tsv` | 3269 | row file: paired differences per dataset of the primary reading, with each one's own lineage interval |
| `figures.tsv` | 7 | derived: single figures the page reads that no summary row holds |
| `intervals.tsv` | 103 | row file: the species interval of each summary row, as the source modes redraw it |
| `order.tsv` | 5 | row file: the clock times of the plan, its three addenda and the final scoring |
| `sets.tsv` | 8 | derived: each set, its role and size |
| `settings.tsv` | 20 | row file: the program, database, release and interval settings |
| `summary.tsv` | 103 | derived: per set, contrast and measure, over the set's datasets |

Every file is plain text, and each TSV has a header row. A file holds counts and measures per dataset and caller, per set and contrast or per learner: no genome id, no sequence, no rule and no row per genome.

## Columns

- `datasets.tsv`: `dataset`, `species`, `drug_class`, `organism_option`, `lineage_split`, `read_before_the_plan`, `served_learner`.
- `callers.tsv`: `dataset`, `caller`, `genomes`, `scored`, `no_call_resistant`, `no_call_susceptible`, `resistant`, `susceptible`, `tp`, `fn`, `tn`, `fp`, `sensitivity`, `sens_exact_lo`, `sens_exact_hi`, `sens_lineage_lo`, `sens_lineage_hi`, `specificity`, `spec_exact_lo`, `spec_exact_hi`, `spec_lineage_lo`, `spec_lineage_hi`, `balanced_accuracy`, `ba_lineage_lo`, `ba_lineage_hi`.
- `differences.tsv`: `dataset`, `reading`, `contrast`, `measure`, `first`, `second`, `first_minus_second`, `lineage_lo`, `lineage_hi`, `first_only_right`, `second_only_right`, `mcnemar_p`, `calls_agree`, `scored`.
- `intervals.tsv`: `set`, `contrast`, `measure`, `median_species_lo`, `median_species_hi`, `mean_species_lo`, `mean_species_hi`, `share_species_lo`, `share_species_hi`.
- `order.tsv`: `utc`, `event`, `file`, `sha256`.
- `settings.tsv`: `key`, `value`.
- `sets.tsv`: `set`, `role`, `datasets`, `species`, `definition`.
- `summary.tsv`: `set`, `contrast`, `measure`, `role`, `datasets`, `species`, `ahead`, `level`, `behind`, `median`, `q1`, `q3`, `mean`, `sign_test_p`, `median_species_lo`, `median_species_hi`, `mean_species_lo`, `mean_species_hi`, `own_interval_above_zero`, `own_interval_spans_zero`, `own_interval_below_zero`, `learner_at_most_0.02_ahead`, `share_learner_at_most_0.02_ahead`, `share_species_lo`, `share_species_hi`, `within_0.02_either_way`.
- `figures.tsv`: `key`, `set`, `value`, `role`.

- `datasets.tsv`: `organism_option` is whether the program has an organism option for the species, `lineage_split` whether the release holds a lineage split for the dataset, `read_before_the_plan` whether its calls were read before the plan was hashed, and `served_learner` the learner the release serves, picked on the random split.
- `callers.tsv`: `caller` is `program:<reading>` for the AMRFinderPlus program, `grm_known_gene_rule` for GRM's known-gene rule and `learner:<learner>:<split>` for a learner. A resistant genome called resistant is `tp`. `*_exact_*` are Clopper-Pearson ends, which treat genomes as independent, and `*_lineage_*` ends from resampling the report's lineages.
- `differences.tsv`: `reading` is the program's reading; the record keeps the primary one. A dataset's difference is `first` minus `second`. `first_only_right` counts genomes only the first caller calls right, and `mcnemar_p` is the exact two-sided p from the two such counts.
- `summary.tsv`: `ahead`, `level` and `behind` count the datasets where the first named caller is ahead of, level with or behind the second. `own_interval_*` count datasets by their own lineage interval. `learner_at_most_0.02_ahead` counts datasets where the learner is at most 0.02 ahead of the program, or behind it, as the plan's rule M defines within 0.02; `within_0.02_either_way` counts those within 0.02 on both sides. Neither is given for the program against the known-gene rule.

## Rules

- **Three callers.** The AMRFinderPlus program is read under each report's rule. `tool_default` is the primary reading; `tool_default_intact` counts no row whose method the program marks as partial or as an internal stop, `genes_core` only gene rows of core scope, and `genes_all_scopes` gene rows of core or plus scope. GRM's known-gene rule is the release's own call from the catalog, and a learner calls resistant at a score of 0.5 or more.
- **Measures.** Sensitivity and specificity come from each caller's counts, with 0 for an empty class, and balanced accuracy is their mean. Differences are exact, so level means a difference of exactly 0.
- **Over datasets.** The median, quartiles and mean of the per-dataset differences, and the exact two-sided sign test on ahead against behind. The species interval resamples the set's species with replacement and takes every dataset of a drawn species. The plan fixed 4,000 resamples and the seed string, and named neither the generator nor the draw: the source modes use `random.Random` seeded with the string, `randrange` over the species sorted by name and `statistics.median`, which give the final scoring's ends. `intervals.tsv` holds the ends; the record does not hold the seed string.
- **Roles.** `primary` rows are the plan's headline: the program against the known-gene rule on primary-new-235 in each measure, and each learner and the served learner against the program in balanced accuracy on its lineage split. `secondary` rows are those the plan listed beside it. `exploratory` rows and figures were chosen after the results were read.

## Commands

```sh
python3 -I scripts/amrfinderplus-comparison/record.py --check reference/amrfinderplus-comparison
python3 -I scripts/amrfinderplus-comparison/record.py --source "$FINAL" --check reference/amrfinderplus-comparison
```

The first rebuilds the derived files, this README and `SHA256SUMS` from the row files. The second rebuilds every file from the final scoring, which is not public.
