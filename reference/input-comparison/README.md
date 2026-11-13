# Input comparison: the second reading across species

`second-reading.json` is the JSON that `scripts/input-comparison/second_reading.py` writes for H-1's four
species reports. [docs/input-comparison-runner.md](../../docs/input-comparison-runner.md), "H-1's reading",
prints the same numbers as tables. The reading changes no report, and it fits and tests nothing again.

It was made from the root of the repository with the command below, then copied here unchanged:

```
python3 -I -B scripts/input-comparison/second_reading.py \
    SAL/report/ablation-report.json SA/report/ablation-report.json \
    NG/report/ablation-report.json MTB/report/ablation-report.json \
    --json /tmp/reading.json \
    --predict cap-4000=better --predict cap-1000=worse \
    --predict-on "Staphylococcus aureus" --predict-on "Mycobacterium tuberculosis"
```

SAL, SA, NG and MTB stand for the four `combine` WORK folders, given by full path. These are the reports it read, by sha256:

```
358c8066115f02d10d7bceb58944f12f4e7335ecfb108c8fc5e4e007ee2164a1  MTB/report/ablation-report.json
1a7a12a06cc5588c431f6dcb77246c2281ff128761f67c4a9b0ef2a4e73e4f80  NG/report/ablation-report.json
6ce7867da5a2ddb8dd2d58951540f82cc5c13fb0aabae171d419da2e113dd7be  SAL/report/ablation-report.json
e44b702dfb44b7ce916a76def768a89fb54f5bb84c2d8b666a61fe1f3ee58a1c  SA/report/ablation-report.json
```

The run used Python 3.12 and numpy 1.26.4. It wrote `second-reading.json` with sha256
`8facf0640741c9991f823fdcb4034f942fdfc17a68b1d35c9d09f85f5d8bd606`.

## The reading by learner

`learner-reading.json` is the JSON that `scripts/input-comparison/learner_reading.py` writes for the same run:
each learner's log loss beside its calls, per split, change and learner
([docs/input-comparison-runner.md](../../docs/input-comparison-runner.md), "A per-learner reading with calls").
It was designed after the results were read, and it decides nothing. It was made with:

```
python3 -I -B scripts/input-comparison/learner_reading.py \
    SAL SA NG MTB --json learner-reading.json --sensitivity
```

SAL, SA, NG and MTB stand for the four `combine` WORK folders themselves, given by full path, not their reports:
the reading pairs the prediction rows they hold. Those folders are not kept in this repository, so this file
cannot be rebuilt without them. It is copied here unchanged as `learner-reading.json` with sha256
`105c2d1cc4d707a38673683e79d594ef1b9a1dbe58afcfe8583320df4e14acdc`. The script that wrote it had sha256
`be5008be8c4e110a772631375f0c3192c805e4020aa542acf2e65aab5f0d4049`, the script in this repository at the time.

## The plan fixed before the results

`reading-plan.md` is the plan that fixed the second reading on 11 November 2026: what to count, the calibration
rule, the firm rule and the four predictions
([docs/input-comparison-runner.md](../../docs/input-comparison-runner.md), "H-1's reading", under "Timing").
It was a private working note until 13 November 2026. It is copied here unchanged as `reading-plan.md` with sha256
`57ba6f212a57e2516477f020d1c3ae225b2f8a6c37fc9941cedcfd7be1b8c24d`, the value the runner document has named since pull request #82.

Three things about it need saying:

- **Nothing public attests its hour.** The note gives its own time, 06:15 UTC on 11 November. The sha256 first
  became public in pull request #82, opened at 19:58 UTC that day, after the results. So the file shows what the plan
  said, and the claim that it came before the results rests on the author's record.
- **It names files that are not in this repository.** The rebuilt reports under `rehearsal/`, the ledger row and
  `analysis/h1_contrasts.py` belong to the author's working folder. `scripts/input-comparison/second_reading.py`
  gives the numbers that script first read.
- **H-1 is the working name** of the input comparison's run on four species.

## What the web app bundles

The web app bundles part of both readings, `web/src/data/input-comparison.json` (of the reading by learner,
the lineage split only).
`scripts/input-comparison/site_data.py` writes that part, and `tests/input_comparison/test_site_data.py`
holds the bundled file to the script's output:

```
python3 -I -B scripts/input-comparison/site_data.py reference/input-comparison/second-reading.json reference/input-comparison/learner-reading.json > web/src/data/input-comparison.json
python3 -I -B scripts/input-comparison/site_data.py reference/input-comparison/second-reading.json reference/input-comparison/learner-reading.json --check web/src/data/input-comparison.json
```
