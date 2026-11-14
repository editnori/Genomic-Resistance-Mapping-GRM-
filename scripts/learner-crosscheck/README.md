# Learner cross-check

GRM runs its own code for each learner. This kit fits each learner on 16 small test panels, once with GRM and once with a standard version of the same method, and writes down where the two differ and why.

It is a software check on small test panels. It says nothing about clinical performance and nothing about how well a learner predicts resistance. It reads no genome data.

The result of the run of 14 November 2026 is in [`reference/learner-crosscheck`](../../reference/learner-crosscheck/README.md), with the plan that was fixed before any comparison was read. `tests/learner_crosscheck` checks that result's shape and counts and needs none of the libraries below.

## What is compared

| Learner | GRM | Reference |
|---|---|---|
| Set covering machine | `scm::fit`, both kinds, 7 values of p, at most 6 rules: 14 fits per panel | `scm_reference.py`, an independent implementation in exact arithmetic, and pyscm-ml 1.1.2 |
| Decision tree | `tree::fit_tree`, depth 8 | scikit-learn 1.7.2, one split at a time on each node's samples |
| Random forest | `tree::fit_forest`, 200 trees | scikit-learn 1.7.2, one split at a time on the samples and the features GRM drew |
| Boosted trees | `gbt::fit`, 400 rounds | XGBoost 3.0.5, one tree at a time from GRM's margin, and its own run end to end |
| L1 logistic regression | not fitted here | its own suite, `cargo test --release --lib learn::lasso` (docs/lasso-stopping.md) |

The panels are the 16 of `tests/data/lasso-stopping-fixtures.txt`: 2 to 128 presence features on 60 to 508 samples, each with a label per sample and a training mask. Eight are named for a public dataset and a fold, eight for the case they hold (`one_class_train`, `tied_counts`, and so on). The file does not say how they were made.

No library can be handed GRM's random draws, and scikit-learn breaks equal gains in a seeded random order. So the tree checks do not ask a library to grow the same forest. The harness writes the draws GRM made, the library searches one split on exactly those samples and features, and exact fractions give the set of features with the largest gain. `PLAN.md` in the result folder gives every rule, and what counts as a difference with a shown cause.

## Files

| File | What it does | Needs |
|---|---|---|
| `main.rs` | The GRM side, a Cargo example (`learner-crosscheck`). It calls the public functions of `grm::learn` on each panel and writes every fitted model, and every random draw, as JSON. It changes nothing in `src/`. | Rust 1.98.1 |
| `panels.py` | Decodes the panel file without the Rust code. | nothing |
| `scm_reference.py` | The independent set covering machine. | nothing |
| `exact.py` | Exact split gains of a class-weighted Gini tree. | nothing |
| `check_scm.py` | The set covering machine against `scm_reference.py` and pyscm-ml. | pyscm-ml |
| `check_trees.py` | The decision tree and the forest against scikit-learn. | scikit-learn |
| `check_gbt.py` | The boosted trees against XGBoost. | XGBoost |
| `compare.py` | Runs the three checks on every panel and writes the tables. | all of them |
| `table.py` | The shape of the tables and the counts that follow from them. | nothing |
| `exploratory.py` | Four counts written after the first comparison was read. No outcome depends on them. | all of them |
| `check_env.py` | Checks that the environment holds exactly the pinned versions. | nothing |
| `run.sh` | The whole check, from this checkout. | all of them |

## Run it

Python 3.12 on Linux x86-64. The reference libraries go into an environment of their own, each file named by its sha256 on pypi.org:

```sh
python3 -m venv ENV
ENV/bin/pip install --no-deps --require-hashes -r scripts/learner-crosscheck/requirements-build.txt
ENV/bin/pip install --no-deps --require-hashes -r scripts/learner-crosscheck/requirements.txt
ENV/bin/pip install --no-deps --require-hashes --no-build-isolation --no-binary pyscm-ml -r scripts/learner-crosscheck/requirements-scm.txt
```

pyscm-ml 1.1.2 has no Linux wheel, so the third command builds it from its source archive with the setuptools and wheel of the first, and needs a C++ compiler. pyscm-ml is under the GPL, version 3. The kit installs it and calls it; no code of it is in this repository.

Then, with Rust 1.98.1:

```sh
scripts/learner-crosscheck/run.sh ENV OUT
```

`OUT` must be a new folder. `run.sh` checks the environment, builds the harness, runs it with 1 and with 8 threads and requires the same bytes from both, runs `compare.py` and `exploratory.py`, and writes `OUT/run.json` with the commit, the tree of `src/`, the compiler and the times. On 14 November 2026 this took 18 seconds with the harness already built, on 8 threads.

XGBoost's wheel prints a warning about an old glibc each time it loads. That is the wheel's own note about how it was built, and it changes no result.

## What it writes

| File | Content |
|---|---|
| `summary.tsv` | One row per panel and check: the outcome, and the units compared, matched, explained and unexplained. |
| `differences.tsv` | One row per difference: panel, learner, check, unit, what differs, the cause, and the numbers that show the cause. An empty cause means none could be shown. |
| `totals.json` | Per check: the panels by outcome and the sums of the units. |
| `detail.json` | Each check's further counts per panel, and the settings production would select on the panel. |
| `versions.json` | Python, the libraries, the sha256 of the panel file and of each panel's GRM fits. |
| `exploratory.json` | The exploratory counts. |
| `run.json`, `environment.tsv` | The commit, the compiler, the times, and the pinned and installed versions. |
| `grm/`, `panels/` | Every GRM fit and everything each check recorded. These are large and are not committed. |

An outcome is `matched` (the plan's rule held), `explained` (there is a difference and each has a shown cause), `unexplained` (a difference without a shown cause, which is a finding) or `not applicable` (the comparison cannot run on the panel).

`python3 -I scripts/learner-crosscheck/table.py OUT` checks the tables of a run and exits 1 if their counts do not add up.

To replace the committed result, copy the eight files of the table above, without `grm/` and `panels/`, into `reference/learner-crosscheck`, update the sha256 lines and the numbers in its README, and update the expected totals in `tests/learner_crosscheck`. The test suite fails until the tables, the README and the Methods text of the site agree.

## Limits

- The panels are small and few. A learner that agrees here can still differ on other inputs.
- A match is the plan's rule, not always an identical model. The end-to-end boosted fit matches when both sides keep the same number of trees and their probabilities agree within 0.0001.
- Where several features have exactly equal gain, no library gives a single answer, so the reference there is the exact arithmetic of `exact.py`.
- The harness uses the production defaults. It does not run `grm evaluate`, so the settings production selects inside each fit are recorded but not compared.
- The L1 logistic regression is not fitted here.
