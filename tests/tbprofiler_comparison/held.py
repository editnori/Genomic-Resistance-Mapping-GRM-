"""The paragraphs and list items of docs/tbprofiler-comparison.md that the suite holds word for word.

lines.py rebuilds the section Reading, the paragraph on the reports the tool does not name, the catalog row's
paragraph, the set-aside paragraph and every table from the record, and test_page holds the review's sentence on the
catalogue's overlap in its four sections. Every other paragraph and list item of the page stands here, in page order,
and test_page holds each to the page once. A deliberate rewording of one of them is a change here too; the figures in
them are held to the record by test_page as well."""

HELD = [
    ('This page counts the calls of TBProfiler 6.6.3, a public catalogue-based caller, against the phenotypes '
     "recorded in GRM's 22 M. tuberculosis reports of evaluation revision 14, beside the five learners' calls "
     'on the same genomes. This is a retrospective comparison on assemblies of public BV-BRC genomes. It is not '
     'a prospective test and not external validation.'),
    ('This is a diagnostic comparison of a published catalogue-based caller: a call is a marker read against a '
     'phenotype and establishes no mechanism. TBProfiler is built for sequencing reads and was run here on '
     'assemblies, which cannot show a mixed population, and it gave no call in 819 of 21,737 genome-drug rows '
     'where an assembly did not cover a target region in full; these figures say nothing about the tool on '
     'reads.'),
    ('The plan (sha256 `a12a03b3`) was hashed at 16:42:57Z on 14 Nov and its addendum (sha256 `099af373`) at '
     "16:51:41Z, both before any recorded phenotype was read against a call, by the hash list's own note. The "
     'scoring run checked both hashes at 17:28:38Z, before the evaluator read any phenotype.'),
    ('The plan fixed `who_assoc` as the primary way of reading, the rule that scores a report, sensitivity '
     'beside specificity per report, and the medians over the 16 scored reports whose drug the tool names, with '
     'the same for each learner on the same genomes, each learner minus the tool, and the lineage interval of '
     'each.'),
    ('The shipped database grades no entry for cycloserine or para-aminosalicylic acid (`grades.tsv`), so the '
     "tool cannot call either resistant; this was written into the plan's addendum, hashed at 16:51:41Z, before "
     'the evaluator first read a recorded phenotype after 17:28:42Z. The addendum added the 14 scored reports '
     "whose drug the database grades as a second set beside the plan's 16, and changed no figure the plan "
     'defines.'),
    ('TBProfiler read each assembly with the database it ships, `tbdb` commit `3ed1cb99`, whose graded entries '
     'carry the source label "WHO catalogue v2". The tool names 18 drugs. A report is primary when the tool has a '
     'drug of the same name, with rifampin read as rifampicin. The four ways of reading its output:'),
    ('- **who_assoc**, the primary way: resistant when the tool reports an entry graded Assoc w R or Assoc w R '
     '- Interim for the drug.'),
    '- **who_assoc_only**: resistant on an entry graded Assoc w R alone.',
    ("- **tool_default**: resistant on any entry the tool reports for the drug, the database's own ungraded and "
     'Uncertain significance entries included, as its printed report does.'),
    ('- **who_assoc_small**: as who_assoc, leaving out an insertion or deletion of 50 bases or more or the '
     'loss of a whole coding region, which a contig end of an assembly can mimic.'),
    ('Each genome and report gets one call: resistant, susceptible or no call, the last where the way of reading '
     'finds no resistant entry and an assembly did not cover a target region in full; a resistant call stands '
     'whatever the coverage. A no call is never read as susceptible. Each learner and the catalog '
     'row are scored on the same genomes, those where the tool has a call.'),
    ('A learner calls resistant at a report score of 0.5 or more, on the lineage split and the random split. '
     'Each learner was scored only on genomes held out of its fit, and production selects its settings inside '
     'each fit.'),
    ('A report is scored when at least 25 resistant and 25 susceptible genomes have a call. Sensitivity is the '
     'share of resistant genomes called resistant, specificity the share of susceptible genomes called '
     'susceptible, and balanced accuracy their mean.'),
    ("Each of the 2,000 resamples draws 881 group numbers with Python's "
     '`random.Random(20261114).randrange(881)`, one call per draw, with the groups numbered in the text order '
     "of their smallest genome id. The interval's ends are the 51st and the 1,950th of the 2,000 sorted values."),
    'An interval needs at least 40 resamples with a value; every interval here had a value in all 2,000.',
    ("The groups are the reports' own lineage ids, joined across the 22 reports: two genomes share a group when "
     'any report puts them in one lineage. A resample takes every genome of each group drawn. A figure has no '
     'value in a resample where a class has no genome, and a median skips a report without one. Brackets on '
     'this page give this 95% lineage interval.'),
    'Every figure on this page is rounded once from the full value in the record.',
    "The shipped database's entries of the two counted grades, per tool drug:",
    ('TBProfiler under `who_assoc` on the 16 primary reports, with the lineage interval in brackets. Graded '
     "says whether the database grades an entry for the report's drug."),
    ('The plan read 4 reports from another drug with a shared name stem, and the call is resistant when any '
     'drug read is. They enter no set and no median.'),
    ("Each cell gives the median over the set's reports, the range in parentheses and the lineage interval in "
     'brackets. The plan fixed the set of 16, and the addendum added the 14. Each median is taken on its own, '
     'so a median sensitivity and a median specificity need not come from one report; the tables above give '
     'both per report.'),
    ("On the lineage split and the same genomes, TBProfiler's balanced accuracy was ahead of or within 0.02 of "
     'the best of the five learners on 12 of 16 endpoints and 12 of the 14 graded ones; the 0.02 margin is an '
     'adapted research screen the plan proposed and no one has adopted, and it is never an applicable FDA '
     'validation criterion.'),
    ("Per learner and split: on how many reports the tool's balanced accuracy is higher or lower, on how many "
     'the tool is ahead or the learner at most 0.02 ahead, the median of the learner minus the tool with its lineage '
     'interval, and the genome-drug rows only one of the two gets right, summed over the set.'),
    ("The best lineage learner on each report of the plan's set is picked after the fact, which favours the "
     "learners. Gap is its balanced accuracy minus the tool's."),
    ("The tool's medians and ranges under each way of reading over each set, and the screen's count. "
     "`who_assoc` is the primary way and the only one the plan's screen reads; the other three show how the "
     'calls move with the grades counted or with large changes left out; their counts in the last column are '
     'exploratory, and the plan named no screen for them. The screen and its label stand in the section '
     'Learners beside the tool.'),
    ('Under `who_assoc`, per primary report: the resistant calls on a susceptible genome by the grade they rest '
     'on, the susceptible calls on a resistant genome by what the tool reported for the drug, and the genomes '
     'with no call by label. An entry of another grade is one the tool reports that `who_assoc` does not count.'),
    '**Exploratory.** The figures in this section were read after the outcomes and are in no plan.',
    ('Counting each no call as wrong or as right gives bounds on every genome. Reading a no call as '
     "susceptible, as the tool's printed report does, gives the lower sensitivity bound and the upper "
     'specificity bound.'),
    "These counts rest on the tool's raw output alone and cannot be rebuilt from the call file.",
    ("An independent review recounted every row of the first report's table from the call file and the reports "
     'and found no blocker. It narrowed the first report in these ways, and this page takes each in the '
     "review's form."),
    ("- The plan's words alone do not fix the lineage interval. The review read the evaluator's two choices, "
     "which the section Method states, and this record's script gives the review's ends."),
    '- Some cells of the first report were rounded twice. Every figure here is rounded once.',
    ("- The word the first report used for the screen's outcome overstates a screen no one has adopted. The "
     "section Learners beside the tool gives the review's wording, and the screen's count has no interval here."),
    ("- The first report's account of how the two drugs without a graded entry were found rests on an unlogged "
     "run. This page gives the fact from the database's grade table."),
    ("- The catalogue's possible overlap with these genomes goes beside every figure that sets the tool beside "
     'a learner.'),
    ("- The tool's raw output and the assemblies it read no longer exist; the per-genome call file, kept "
     'outside this repository, is the lowest layer these counts can be rebuilt from.'),
    ('- The tool read assemblies. It is built for sequencing reads, and an assembly cannot show a mixed '
     'population.'),
    "- The overlap of the catalogue's source isolates with these genomes was not counted.",
    '- The best of the five learners is picked per report after the fact, which favours the learners.',
    ("- A no call leaves a genome out for every caller. The section Exploratory gives the tool's bounds on "
     'every genome.'),
    "- One cutoff: a learner's score of 0.5 or more calls resistant.",
    ('- The record leaves out these figures the plan named: the exact per-report interval, which treats related '
     'genomes as independent, the McNemar and sign tests, the quartiles of the differences, the per-fold '
     'measures, and sensitivity and specificity over all genome-drug rows taken together.'),
    '- One database release, `tbdb` commit `3ed1cb99`, and one tool version.',
    '- The page gives no figure for a single genome, and it names no genome.',
    ('- `reference/tbprofiler-comparison/` holds the record: 20 files that '
     '`scripts/tbprofiler-comparison/record.py` writes from the input set. Its README lists each file and its '
     "columns, and `SHA256SUMS` gives each file's sha256."),
    ('- The record holds rows per report and caller, per way of reading, per tool drug, per grade and per '
     'learner, and single counts. It holds no genome id, label, call or score per genome.'),
    "- `checks.tsv` gives 0 as the count of the input set's genome ids found in the record.",
    ("- The input set holds the call file and the reports' scores, which are per genome, so the repository does "
     'not hold it. `--check` alone derives every derived file from the row files.'),
    ('- `tests/tbprofiler_comparison` rebuilds the section Reading, three paragraphs and every table from the '
     'record, holds every other paragraph and list item word for word, and checks that each other number stands '
     'in the record.'),
]
