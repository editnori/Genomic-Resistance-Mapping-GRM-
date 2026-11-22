# MIC models on a fold plan

An experiment, off by default. `grm learn --exp-mic` fits MIC models on a fold plan the caller picks, in place of the MIC regression's own folds. Without the flag, `grm learn` is unchanged, byte for byte. No report is read or written, `EVALUATION_REVISION` stays 14 and the report schema stays 2.

A report's MIC regression (`evaluate_mic` in `src/learn/eval.rs`) draws stratified folds over the genomes with an MIC, so a test genome's relative or copy can train its fold. Here the folds come from a plan:

- the lineage split that `grm learn --seed S` draws (`grouped_folds` over the cohort's lineages), with `--exp-mic-out DIR`, written to `DIR/<dataset>.lineage-s<S>.mic.json`;
- a fold file (`docs/fold-file.md`), with `--fold-file F --fold-out D`, at the release's seeds, written to `D/<dataset>.<scheme>.mic.json`.

`--exp-mic` takes a comma list of arms:

- `released`: the report's censored MIC regression, fold for fold, on the plan's folds. On the folds `evaluate_mic` draws it gives the same predictions, bit for bit.
- `ordered`: one logistic boosted fit per dilution threshold of the training genomes. An exact MIC enters every threshold, one at most v only the thresholds at v and above, and one above v only those at v and below. Each threshold's probability is put back on its training share, made non-increasing by pool-adjacent-violators, and read off at one half.

Each result also holds the training median of each fold (the lower median of the training genomes' exact MICs) and the breakpoint `infer_breakpoint` fits on each fold's training genomes. Every fit and breakpoint reads only the fold's training genomes with an MIC. The release's gate holds: at least 100 genomes with an MIC and 50 exact.

A result records the evaluation revision, the input fingerprint the dataset's report records for the same inputs, the settings, the plan and the arms, then per fold the counts, breakpoint, median and fit counts, and one row per test genome with an MIC: its fold, MIC, censor, reported call and each arm's prediction, with no genome id. An existing result is refused unless `--force`. `--exp-mic` is refused beside `--learner-seed`, the `--exp-tie-*` options and `--exp-null-refit`, and `--exp-mic-out` beside `--fold-file`.

Any split here is a retrospective split inside one public database. A split by lineage or by study is not external validation, and an agreement line of 0.90 is an adapted research screen, not an FDA validation criterion.

Tests: `src/learn/micplan.rs`, `analysis::tests::the_mic_experiment_writes_one_result_per_plan_and_nothing_else` and `tests::learn_takes_the_mic_experiment_only_with_an_output_and_no_other_experiment` in `src/main.rs`.
