"""The paragraphs and list items of docs/time-split.md that the suite holds word for word.

lines.py rebuilds the section Reading, the two exploratory paragraphs under the study and time table and every
table from the record. Every other paragraph and list item of the page stands here, in page order, and
test_page holds each to the page once. A deliberate rewording of one of them is a change here too; the
figures in them are held to the record by test_page as well."""

HELD = [
    ('This page compares two ways of choosing training genomes. The past arm trains only on genomes public in '
     "BV-BRC by each dataset's cut day. The random arm trains on a random draw of the same size and class balance. "
     'Both are scored on the same later genomes of studies new after that day, in 24 datasets of 4 species at '
     'evaluation revision 14. This is a split by first-public day inside one database, BV-BRC. It is not a '
     'prospective test and not external validation.'),
    ('Five learners fit each arm: scm, cart, rf, gbt and lasso. As the test genomes are later genomes of new '
     'studies, this is a retrospective test on newly deposited studies from the same database. Every figure on '
     'this page comes from the rerun after the review, which the section Review describes, except where that '
     'section gives what the fix changed.'),
    ('The plan (sha256 `d415dbd4`) was hashed with its analysis script (sha256 `3264c0f6`) before any fold file, '
     "fit or result of the experiment existed, by its hash list's own note. It fixed the design, the rule that "
     'admits a dataset, recall and log loss per learner on three contrasts, the species interval and the rule for '
     'a cost.'),
    ('The plan named specificity and balanced accuracy as supporting figures and gave them no rule. The table '
     "layout script was hashed after the first result and, by its hash list's own note, before any score was read."),
    ("The rerun's plan (sha256 `5b886192`) was written after the results of 14 Nov and the review's counts had "
     'been read. It fixed what runs and what is compared, and adds no rule. The 14 Nov rule applies to the rerun '
     'as written.'),
    ('**Exploratory.** Every figure the plan did not fix is exploratory and enters no rule: the share of log loss '
     'from wrong calls at the clip and balanced accuracy beside recall and log loss in the main table, which the '
     "rerun's plan added, the comparison of the study and time contrasts, the share of past genomes in each arm's "
     'training genomes as a reading of the study contrast, and the interval over every ordered species draw.'),
    ('A dataset is one species and one antibiotic, and its cohort is the genomes `grm learn` evaluates for it. A '
     "genome's day is the day its record entered BV-BRC, so it is first public inside this database. A genome's "
     'study is its BV-BRC BioProject.'),
    ('- **Cut day.** The earliest day by which at least 60% of the cohort was public. Past genomes have a day on '
     'or before it, later genomes a day after it.'),
    ('- **Test studies.** A study is new when none of its genomes of the species is public by the cut day. A keyed '
     'hash of the accession picks about half of the new studies with a later genome; it depends on no dataset or '
     'label.'),
    ('- **Test set.** A keyed hash picks two in three of the later genomes of the test studies. The rest are the '
     'reserve. All three arms are scored on exactly this test set.'),
    ('- **Arms.** past trains on the past genomes. study-out draws from the past genomes and later genomes of '
     'other studies, with no genome of a test study. random draws from the study-out pool and the reserve, so it '
     "holds genomes of the test genomes' own studies and period."),
    '- **Size.** Each arm trains on as many resistant and as many susceptible genomes as the past arm.',
    ('- **Left out.** A later genome with no recorded study is in no pool and no test set. A copy of a test '
     'genome, with identical inputs, is in no pool.'),
    ('- **Entry.** A dataset enters when its test set holds at least 20 genomes of each class, each arm trains on '
     'at least 25 of each class, and at least 20 reserve genomes and 20 later genomes of other studies remain in '
     'the pools.'),
    ('- **Datasets.** 24 datasets in 4 species enter: 6 of *Escherichia coli*, 8 of *Klebsiella pneumoniae*, 2 of '
     '*Pseudomonas aeruginosa* and 8 of *Salmonella enterica*.'),
    ('- **Fits.** Each arm is one fold file for `grm learn --fold-file` ([fold-file.md](fold-file.md)), with the '
     'default settings of `grm learn` (seed 42, 4,000 screened features, 3 inner folds). Each fit selects '
     'candidate features and tunes its settings inside its own training genomes.'),
    ("Three contrasts compare the arms, each a paired difference on one dataset's test genomes. total is past "
     "minus random, the plan's question. study is study-out minus random: the test genomes' own studies absent "
     'from training against present. time is past minus study-out: later genomes of other studies absent against '
     'present. The plan gave the study and time contrasts no rule.'),
    ('A score of 0.5 or more calls resistant, and recall is the share of resistant genomes called resistant. Log '
     'loss, unlike recall, reads the score itself: it clips each score to [0.0001, 1 - 0.0001], so one wrong call '
     'at the clip adds 9.21 to a sum over the test genomes. scm scores only 0 or 1, so its log loss is about 9.21 '
     'times its share of wrong calls, the calls its recall and specificity count.'),
    ('Each contrast is read as its median over the 24 datasets, with a 95% interval: 10,000 draws of the 4 species '
     'with replacement (seed 20261114), each the median over the datasets of the drawn species, from the 251st to '
     'the 9,750th sorted draw.'),
    ('4 species give only 35 distinct draws (256 in order), so the interval excludes zero only when nearly every '
     'mix of species agrees in sign. Over all 256 ordered draws the same two quantiles give the same reading in 60 '
     'of 60 readings.'),
    ('A median reads worse when its interval lies wholly on the worse side of zero, better when it lies wholly on '
     'the other side, and not shown otherwise. Worse is lower recall, balanced accuracy or specificity, or a '
     'higher log loss.'),
    ('For one learner the total contrast shows a cost when recall or log loss reads worse and neither reads '
     'better. The experiment shows a cost when at least 3 of 5 learners show one and no learner reads better on '
     'recall or log loss, on at least 12 datasets in 3 species.'),
    'Past minus random. Each cell gives the median over 24 datasets, the interval and the reading.',
    'Study is study-out minus random, and time is past minus study-out.',
    ("The study contrast takes the test genomes' own studies out of training, and it also raises the share of past "
     'genomes in training from 76.19% in the random arm to 81.47% in the study-out arm. So it moves study '
     'membership and the mix of periods together. Later genomes can also differ from past ones in place, '
     'laboratory and method, and the design separates none of these.'),
    ('Test genomes of the 24 datasets added up per learner and arm, as the table layout script adds them; a genome '
     'tested in several datasets of a species counts once per dataset. These sums are a description, and the '
     "plan's summary is the median over datasets above."),
    ('The last column is exploratory: the share of the pooled log loss from wrong calls at the clip, to read '
     "beside recall. scm's share is 100.0% at one decimal, as a right call at the clip adds about 0.0001 to the "
     'sum.'),
    ('The test sets hold 3,390 genome and dataset pairs, 1,729 of them resistant, from 828 distinct genomes. A '
     'genome tested in several datasets of its species is one pair per dataset and one distinct genome. A test set '
     'holds 54 to 230 genomes (median 138.5), and each arm trains on 478 to 1,207 genomes (median 815.5), the same '
     'number in all three arms of a dataset.'),
    ('The tested genomes come from 22 distinct studies. The per-species counts below add up to 26, as 4 studies '
     'have tested genomes in two species. A test set holds genomes of 1 to 6 studies (1 to 7 drawn), and 3 '
     "datasets hold one. In 6 datasets a drawn study has no tested genome; the rerun's report, not the record, "
     'gives the reason: all its later genomes went to the reserve.'),
    ('An independent review of the run of 14 Nov found that genomes a fold file lists as unused reached the fits '
     'by two paths. The support floor, the fewest training genomes a candidate feature must be present in and '
     'absent from, came from the whole cohort, and a copy group of unused genomes alone moved tuning onto grouped '
     'draws. With `grm learn --fold-file` fixed to take both from the tested and training genomes alone, the 72 '
     'fold files ran again unchanged.'),
    ('47 of the 360 fits changed a score (3,224 of 50,850 scores), in 28 of the 72 result files and 20 of the 24 '
     "datasets. 6 of the 60 readings changed, none of them on the total contrast, and the rule's answer stayed "
     "cost shown. cart's and gbt's study recall, gbt's study specificity and cart's time log loss now read not "
     "shown; lasso's time log loss and balanced accuracy now read worse."),
    ('The review also found that the test sets hold genomes of 22 distinct studies where the first report said 26, '
     "which counted each species' studies apart, and 1 to 6 studies per dataset where it said 1 to 7. The tables "
     "on this page are the rerun's alone."),
    ('- One database, BV-BRC, with one curation of its labels and studies. Training and test genomes take their '
     'labels from the same kind of record.'),
    ('- The day is the day a genome entered BV-BRC, which can differ from the day its sample became public in the '
     'sequence archives. A split on the archive day could order genomes differently.'),
    '- Later genomes can differ from past ones in study, place, laboratory and method all at once.',
    ('- The study contrast takes the test studies out of training and also raises the share of past genomes, so '
     'the design does not separate study membership from the mix of periods, nor either from place, laboratory or '
     'method.'),
    '- Past genomes with no recorded study train in every arm, and later ones are left out.',
    ("- Two in three of each test study's later genomes are tested, so the random arm holds a smaller share of a "
     "test genome's study than a five-fold random split would, which trains on four in five."),
    '- One draw per arm, one cut per dataset and 4 species: read the direction as shown and the size as coarse.',
    '- The interval resamples species only. It holds nothing of the draw, the seed or other species and datasets.',
    '- Recall and specificity rest on one cutoff: a score of 0.5 or more calls resistant.',
    ("- scm's log loss follows only its count of wrong calls at the clip, so it says nothing more about scm's "
     'scores than recall and specificity do.'),
    '- The page gives no figure for a single dataset, and it names no feature and no study.',
    ('- `reference/time-split/` holds the record: 18 files that `scripts/time-split/record.py` writes from the '
     "input set. Its README lists each file and its columns, and `SHA256SUMS` gives each file's sha256."),
    ('- The record holds rows per dataset, species, arm and learner, and single counts. It holds no genome id, '
     'label or score per genome, no fold file and no list of test-set studies.'),
    "- `checks.tsv` gives 0 as the count of the input set's genome ids found in the record.",
    ('- The input set holds the result files, whose scores are per genome, so the repository does not hold it. '
     '`--check` alone derives every derived file from the row files.'),
    ('- `tests/time_split` rebuilds the section Reading, the study and time reading and every table from the '
     'record, holds every other paragraph and list item word for word, and checks that each other number in '
     'figures stands in the record.'),
]
