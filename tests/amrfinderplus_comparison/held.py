"""The paragraphs and list items of docs/amrfinderplus-comparison.md that the suite holds word for word.

lines.py rebuilds the section Reading, the paragraphs that give the record's figures and every table from the
record. Every other paragraph and list item of the page stands here, in page order, and test_page holds each to
the page once. A deliberate rewording of one of them is a change here too; the figures in them are held to the
record by test_page as well."""

HELD = [
    ('This page compares three callers on the same genomes and recorded labels of public BV-BRC datasets: the '
     "AMRFinderPlus program read under each report's rule, GRM's known-gene rule and GRM's five learners. This is "
     'a retrospective comparison inside one public database, BV-BRC. It is not a prospective test and not external '
     "validation. Neither the AMRFinderPlus program nor GRM's known-gene rule reads a label, but the catalog both "
     'read is newer than most of these genomes and may have been curated from the same public isolates, which '
     'favours both against the learners.'),
    ('The program read the genomes of the 253 datasets of release `rev14-274-20261112` that have a published rule, '
     '23,751 genomes in all. The headline set, primary-new-235, is the 235 datasets with an organism option and a '
     'lineage split whose calls were not read before the plan was hashed. Every figure on this page comes from the '
     'final scoring of 16 Nov, and the section Review gives what an independent review narrowed.'),
    ('The plan (sha256 `5cae853c`) was hashed on 14 Nov, while the program ran and before any call of the datasets '
     'outside the first 2 was scored. Those 2 were read before the plan, so primary-new-235 leaves them out. The '
     'plan fixed the sets, the four readings of the program, the measures, the species interval and rule M.'),
    ('Its headline is the program against the known-gene rule on primary-new-235 in each measure, and each learner '
     'and the served learner against the program in balanced accuracy on the lineage split. Every other summary it '
     'lists is secondary, and a figure it did not list is exploratory.'),
    ('Three addenda (sha256 `d6eda978`, `c8b8500e` and `a6239776`) change where and in what order the program ran '
     'and which repeat checks stand in for lost files. Each says it changes no comparison, set, measure, interval, '
     'reading or rule of the plan. The second was hashed after part of the cohort had been scored, and it changes '
     'only the order of the genomes still to run.'),
    ('- **Datasets.** A dataset is one species and one antibiotic of the release with a published rule. Its '
     'genomes and labels are those the release reports for it.'),
    ('- **The AMRFinderPlus program.** Version `4.2.7` with database `2026-08-07.1`, given the organism option '
     "where the species has one. A genome is called resistant when every group of the report's rule is met by at "
     'least one counted row. `tool_default`, the primary reading, counts rows of core scope; `tool_default_intact` '
     'counts no row whose method the program marks as partial or as an internal stop; `genes_core` counts gene '
     'rows of core scope only; and `genes_all_scopes` counts gene rows of core or plus scope. The record keeps all '
     'four readings per dataset, and this page reads `tool_default` alone.'),
    ("- **GRM's known-gene rule.** The call the release itself makes for each genome from the catalog. It is not "
     "the program's call, and the page never gives it the program's name."),
    ("- **GRM's learners.** scm, cart, rf, gbt and lasso, from the scores the release's reports hold, calling "
     'resistant at a score of 0.5 or more, on the random split and on the lineage split. The learner the release '
     'serves is picked on the random split.'),
    ("- **Measures.** Sensitivity, specificity and balanced accuracy, their mean, from each caller's counts on "
     'each dataset. A difference is the first caller minus the second on one dataset, exact, so level means a '
     'difference of exactly 0.'),
    ('- **Over datasets.** The counts ahead, level and behind, the median and quartiles, the mean and the exact '
     "two-sided sign test. The species interval draws the set's species with replacement, as many as it has, takes "
     'every dataset of each drawn species, and repeats 4,000 times; its ends are the 101st and the 3,900th sorted '
     'value. A set with fewer than 6 species gets none.'),
    ('- **Generator and draw.** The plan fixed the seed string and named neither the generator nor the draw. The '
     "record's script uses Python's `random.Random` seeded with that string, `randrange` over the species sorted "
     "by name and `statistics.median`, and gives the final scoring's ends at 6 decimals."),
    ("Each row gives the program read under the report's rule minus GRM's known-gene rule over a set's datasets. "
     "primary-new-235 is the plan's headline, and primary-237 adds the 2 datasets read before the plan."),
    ('Each row gives one learner minus the program on the lineage split of primary-new-235, in balanced accuracy '
     'with the counts ahead, level and behind, and in sensitivity beside specificity. The served learner is the '
     'learner the release serves, picked on the random split, and the highest of the five is an oracle no one can '
     'pick in advance.'),
    ('The next table is the random split, secondary for the five learners and the served learner: training and '
     'test share lineages, so it is no test on unseen lineages, and the served learner is picked on it.'),
    ("These sets are drug classes of the datasets. The table gives each set's counts, median and mean in balanced "
     'accuracy, with the species interval where a set has at least 6 species.'),
    'Both sets have 2 species, too few for a species interval.',
    ("An independent review recomputed the figures from the final scoring's tables and the release before this "
     'record was written. It narrowed the reading in these ways:'),
    ('- The counts ahead, level and behind and the median lead. The mean follows with the reason it sits above the '
     'median, and no sentence joins ahead with level.'),
    ('- At 0.02 the page gives the two-sided count. The one-sided count, the learner at most 0.02 ahead of the '
     'program or behind it, appears only with that definition.'),
    ('- The sentence on the catalog travels with the figures that set a learner beside the program: in the first '
     "paragraph, the section Reading and the record's README."),
    ('- Each row and figure is labelled primary, secondary or exploratory as the plan fixed it, and the 26 '
     'quinolone datasets of primary-new-235 are named by drug class alone.'),
    ('- Every random-split figure carries the label that training and test share lineages and that the served '
     'learner is picked on that split.'),
    ('- The two disagreement tables of the final scoring stay out of the record, as some of their rows would give '
     "the program's call on one genome."),
    ("- Rule M appears only in the review's words, as a proposed screen no one adopted and never an applicable FDA "
     'validation criterion.'),
    ("- The plan's sensitivity run on the newest catalog database (2026-10-29.1) was not run, so these figures "
     'hold for database 2026-08-07.1 only.'),
    ("- The program's tables per genome are gone, so its figures can be checked only at the counts per dataset and "
     'per fold. The record holds the counts per dataset.'),
    ('- The species interval depends on a generator and a draw the plan did not name. The record holds its ends, '
     'and not the seed string the plan fixed.'),
    '- Every caller is scored against the same recorded labels of one database, BV-BRC, with one curation.',
    '- The page gives no figure for a single dataset, and it names no gene and no catalog row.',
    ('- `reference/amrfinderplus-comparison/` holds the record: 11 files that '
     '`scripts/amrfinderplus-comparison/record.py` writes from the final scoring. Its README lists each file and '
     "its columns, and `SHA256SUMS` gives each file's sha256."),
    ('- The record holds rows per dataset and caller, per set and contrast, and per learner, and single figures. '
     'It holds no genome id, no label or call per genome and no rule.'),
    ("- The final scoring's tables, the plan and the set lists are not public, so the repository holds the record "
     'alone. `--check` alone derives every derived file from the row files.'),
    ('- `tests/amrfinderplus_comparison` rebuilds the section Reading and every table from the record, holds every '
     'other paragraph and list item word for word, and checks that each other number on the page stands in the '
     'record.'),
]
