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

Any split here is a retrospective split inside one public database. A split by lineage or by study is not external validation, and an agreement line of 0.90 is an adapted research screen, never an applicable FDA validation criterion.

Tests: `src/learn/micplan.rs`, `analysis::tests::the_mic_experiment_writes_one_result_per_plan_and_nothing_else` and `tests::learn_takes_the_mic_experiment_only_with_an_output_and_no_other_experiment` in `src/main.rs`.

## Running one species

`scripts/c6-mic/` runs the experiment for one species on a box of its own and scores it:

- `run_species.sh KIT SPECIES LINEAGE STUDY OUT [JOBS]` takes the species' datasets on the lineage split and on the study split as comma lists (`-` for none). It builds `grm` from the checkout and records its commit and sha256, builds the species' data root from public files checked against a pinned list, writes the study plan for seed 47, fits, checks that every result exists, runs the probe, scores, and writes `OUT/DONE` with the sha256 of its list of output files. Any failure writes `OUT/FAILED` and exits 1. Run again, it skips finished steps and fits only missing results. It stops nothing but itself.
- `expect.py` lists the results a run must hold. The lineage split runs at seeds 42 to 47. The study split runs the five reviewed study plans (`study-whole`, `study-d1` to `study-d4`) at seeds 42 to 46 and `c6mic-s47` at 47. Seeds 42, 43 and 44 fit `released,ordered`; seeds 45, 46 and 47, the twin, fit `released` only. Every result holds the training median.
- `study_plans.py` writes `c6mic-s47`. A unit is the genomes all five reviewed study plans test in one fold, so a unit keeps each study whole. Units are placed by the reviewed plans' rule B, in the order of `sha256("c6mic-20261122|study-s47|" + the unit's lowest genome id)`. A plan that splits a unit or a k-mer copy group, or leaves a fold without both classes, is not written.
- `score.py` checks each result first: format, plan and arms, the settings and input fingerprint against the dataset's release report, the test counts of each fold, each fold's training median recounted from the other folds, and a prediction for every tested genome of an available arm. The fingerprint hashes every setting, the seed too, and the release holds reports at seed 42 only, so a lineage result at seeds 43 to 47 must hold the report's settings but the seed and the cohort of the dataset's seed-42 result (its counts, and each tested genome's MIC, censor and call); the study plans run at the release's seeds and match the report itself. A dataset whose curated breakpoint differs from the release's passes only through its probe (below). Then it writes one row per dataset, split, seed and arm (`released`, `ordered`, `median`): EA-exact, then category agreement, the very major error rate and the major error rate, each with its count and denominator. EA-exact reads MIC values only. The other three read the R and S labels, through the call and through each fold's training breakpoint. An unavailable arm, or a fold with MIC test genomes and no breakpoint, leaves the figures empty and is counted. No figure is taken over part of the folds. The scores hold no row per genome and no genome id.
- The probe. The input fingerprint hashes a dataset's breakpoint, which curation fits on every MIC pair of the species and drug. A table that passes one more genome through QC moves the breakpoint's `n` and `agreement` and so the fingerprint, while the dataset's genomes, calls and MICs stay the same, and no fit reads that breakpoint: each fold fits its own on its training genomes. `release_specs.py` lists the datasets whose antibiotic and class minimum are the release's and whose breakpoint is not (it stops on any other difference) and writes a copy of `datasets.json` with the release's breakpoints. The run refits seed 42 of those datasets on a copy of the data root holding it, into `results/probe`, and needs the release's fingerprint from each. The scorer then needs each probe to equal the run's seed-42 result but for its fingerprint and wall time, and every result of the dataset at the release's settings to hold the run's seed-42 key.

Tests: `tests/c6_mic`, on SYNTHETIC cohorts and results only.
