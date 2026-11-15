# The tie rule experiment of 15 November 2026: the record

This folder holds the record behind [docs/tie-rule.md](../../docs/tie-rule.md). The experiment fitted 169 datasets of release `rev14-274-20261112` (evaluation revision 14) twice, once with the present rule that the tree learners use to break ties between equal splits and once with an exact rule, and an independent review recomputed its figures on the same day. A change under the exact rule is a sensitivity of a figure to the tie rule. It does not show that either rule is more accurate. Every figure is cross-validated on public BV-BRC genomes.

`scripts/tie-rule/record.py` wrote every file here but this README, `PLAN.md` and `PLAN.sha256` (which it copies unchanged from the experiment's output and checks against each other) and `SHA256SUMS`, which it writes last, over every other file. It reads the five arms' reports and traces and the release's index, and fits nothing. The reports and traces are not here; `inputs.sha256` gives the sha256 of each.

## The arms

| Arm | What ran |
|---|---|
| A | `main`'s binary at `852ebc9` |
| B | the branch `exp/tie-rule-exact` binary at `88af130`, flag off: the present rule |
| C | B with `--exp-tie-trace`: the trace of the present trees |
| D | the branch's binary with `--exp-tie-rule exact`: the exact rule |
| E | D with `--exp-tie-trace`: the trace of the exact trees |

Every arm used the release's matrices, folds and seeds. `settings.tsv` gives the full commit ids, the binaries' sha256, and the settings of the paired test.

## Commands

From the repository root, with `RUNS` the folder holding `A/` to `E/`, `R` the unpacked release and `P` the experiment's folder holding `PLAN.md` and `PLAN.sha256`. Python 3.12, standard library only.

```sh
python3 -I scripts/tie-rule/record.py "$RUNS" "$R" "$P" --out reference/tie-rule
python3 -I scripts/tie-rule/record.py "$RUNS" "$R" "$P" --check reference/tie-rule
```

The second prints the name of each file that is not the text the script writes, and fails if there is one. CI has neither the reports nor the release, so it does not run it; `tests/tie_rule` runs the script on SYNTHETIC reports written in the test, and checks the committed files against each other. To check the files' bytes, run `sha256sum -c SHA256SUMS` from this folder.

## Files

The row files hold what the script reads from the reports and traces. Every other `.tsv` file is derived from the row files alone, so it can be written again from this folder.

| File | Kind | One row per |
|---|---|---|
| `datasets.tsv` | row | dataset of the release: species, genomes, resistant genomes, features, and whether it was run |
| `measures.tsv` | row | dataset, split and learner: recall, specificity, balanced accuracy and log loss under each rule and the change, correct calls, held-out calls that change, and whether the scores, the metrics and the fold balanced accuracy change |
| `reproduction.tsv` | row | dataset: A against the release, B against A, C against B and E against D, byte for byte with the `seconds` fields masked |
| `seconds.tsv` | row | dataset and learner: the learner's seconds summed over both splits, in arms A, B and D |
| `traces.tsv` | row | dataset and kind (`cart`, `rf`, `gbt` for the logistic models, `gbt-mic` for the MIC models): fits, the trace counts of the present trees (C) and of the exact trees (E), and the tree digests compared position by position |
| `inputs.sha256` | row | report of A, B and D and trace of C and E, in the form `sha256sum` writes |
| `splits.tsv` | derived | dataset and split: the best learner under each rule, the gap from the best learner to the second under the present rule, and the tie effect, the largest change in balanced accuracy over cart, rf and gbt |
| `order.tsv` | derived | split and pair of learners: the datasets where their order by balanced accuracy changes |
| `paired.tsv` | derived | tree learner, split and measure: the mean change, exact minus present, with its 95% interval from a species bootstrap |
| `summary.tsv` | derived | figure: every count, share, median and largest change the page quotes, by key |
| `settings.tsv` | fixed | setting: the commits, binaries and test settings, as the script holds them |
| `PLAN.md`, `PLAN.sha256` | copied | the experiment's plan, fixed before any run, and its sha256 |

A held-out call is a score at or above 0.5. Log loss is the mean over held-out genomes of the negative log likelihood of the score clamped to [1e-12, 1 - 1e-12]. The best learner is the index's rule (`src/report.rs`, `best_learner`): the highest balanced accuracy, a tie to the later learner of scm, cart, rf, gbt and lasso. The trace counts come from the branch's own counters, one line per fit: nodes searched, nodes with more than one candidate at the best exact impurity (`exact_tied`), and nodes where the exact pick differs from the present one, as an exact tie (`differ_tie`), a strictly better exact split (`differ_order`) or, for gbt, a lower candidate with the same two groups (`differ_partition`).

## Not here

The reports and traces `inputs.sha256` lists, the matrices, the binaries, and the experiment's and the review's reports and scripts stay on the machine that ran them. No string of a report other than dataset ids, species slugs and learner names is written here.

## Limits

- 169 of the 274 datasets were run. The 105 of the 4 largest species were not, and nothing here speaks for them.
- For gbt, two candidates are equal only when they divide the node's fitting samples into the same two groups. Candidates with equal gains and different groups are not counted.
- A trace counts nodes per fit. It holds no record per node.
- The seconds were measured on a machine shared with other work, and are rough.
