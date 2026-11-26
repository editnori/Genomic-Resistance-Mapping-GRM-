# The run on the datasets of M: its tables

This directory holds the tables of plan block C3's run on its 35 datasets of M ([docs/input-arms.md](../../docs/input-arms.md), The run on the datasets of M). The run has not happened, so one table is here: the first run's step times, which the cost estimate reads.

| File | Holds |
|---|---|
| `first-run-times.tsv` | each step of the first run in the 9 species of M: species, step, whole seconds, and the sha256 of the species' `logs/timing.tsv` the row was read from |

`scripts/c3-inputs/cost.py times` wrote it from the 9 timing tables, each held to its own `SHA256SUMS` line, and `scripts/c3-inputs/cost.py measured` reads it by default. The first run's own tables are in [reference/input-arms](../input-arms/), whose `SHA256SUMS` does not list this file.
