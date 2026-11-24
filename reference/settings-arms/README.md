# The settings reading: the record

This directory holds the tables behind [docs/settings-arms-result.md](../../docs/settings-arms-result.md): plan block C5, a nested search of each learner's settings among 31 trials under the default-anchored rule at 0.04 (`search`), and a call cutoff chosen on inner scores from a 19-point grid (`cutoff`), each against the release's rule, beside a twin and the plain best with no anchor (`plain-best`), read once on the first run's 55 datasets in 8 species, each with a lineage and a study split. Every row is an aggregate over a dataset, a species, an arm, a fold count or a learner. No file holds a genome, a sequence or a feature, and no file holds a trial's settings: a trial is named by its id alone.

`scripts/c5-settings/score.py --count` wrote the four tables under `count/`. `checks.txt` is the printout of `scripts/c5-settings/run_checks.py`, the checks of the run's copies made before any score was read. `scripts/c5-settings/page_tables.py write` checked the reading against the 16 gates of `gates.txt`, and wrote every other file here once they passed. `reading.txt` lists the reading's commands, each with the sha256 of what it wrote. The rules of the reading were fixed at 05:52Z on 23 Nov 2026, before any C5 score was read (`docs/assumption-plan-corrections.md`, entry "C5 design").

| File | Holds |
|---|---|
| `checks.txt` | the checks of the run's copies before the reading, counts only |
| `gates.txt` | the gates the table script passed before it wrote a table |
| `reading.txt` | the reading's commands and the sha256 of each one's output |
| `count/count.tsv`, `count/conditions.tsv` | the count that decides, and the plan's conditions per arm as `score.py` read them |
| `count/count-plain-best.tsv`, `count/count-no-constant.tsv` | the counts beside it, each in no Holm family and deciding nothing |
| `page-conditions.tsv` | the conditions per arm, `decides` 1 for ruling 15 and 0 for the stricter reading beside it; condition 6 is the maintainer's |
| `page-count.tsv` | every row of the count that decides, with the 15 rows of the Holm family marked and conditions 1 and 2 read on each; the cutoff arm's log loss rows marked equal by construction, deciding nothing |
| `page-beside.tsv` | every row of the two counts beside it, each marked as deciding nothing |
| `page-scm.tsv` | `scm`'s search and plain best against the default on balanced accuracy, beside the synthetic loss of 0.0126 (0.0088 to 0.0163) the plan names for the default-anchored rule |
| `page-choices.tsv` | per split and learner, the outer folds searched: those that kept trial 0, the plain best's place, the folds where trial 0 or every trial failed inside, and the cutoffs moved from 0.5 |
| `page-trials.tsv` | per split, learner and trial id, the folds the search picked it on and the folds it was the plain best on; ids alone |
| `page-cutoffs.tsv` | per split, learner and grid point, the folds whose cutoff it was |
| `page-species.tsv` | per species, named in full, split, pair and learner: better, worse and the same on the primary measure; it decides nothing |
| `page-datasets.tsv` | per dataset, species and antibiotic in full, split, pair and learner: the default's value and the arm's on balanced accuracy, recall and log loss, the differences, the draws scored, the genome count, the constant-model draws and a study run refused |
| `page-checks.tsv` | `score.py`'s checks of every result as checks passed of checks made, and the table script's own checks of the cutoff arm |
| `page-binaries.tsv` | each species' commit, binary and code sha256, and its restarts |

All 8 species ran one binary by its code sha256, `343172462952f291381d765a54530c7f913bf101b04b4f9ac5a6ad27a16f9ee8`, built from 1 commit: `0d987ba`.

A dataset's figure is the mean of its draws: seeds 42, 43 and 44 on the lineage split, the twin's 45, 46 and 47, and one study run at the release's seeds. The cutoff arm moves no score, so its log loss equals the default's on every dataset by construction, and its log loss rows decide nothing. `scm` scores 0 or 1, so its cutoff arm is the default on every fold. A split by lineage or study is a retrospective split inside one public database, never external validation.

`SHA256SUMS` gives the sha256 of every other file here. `scripts/c5-settings/page_tables.py print` prints the page's tables from these files.
