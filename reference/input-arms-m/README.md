# The run on the datasets of M: its tables

This directory holds the tables of plan block C3's run on its 35 datasets of M ([docs/input-arms.md](../../docs/input-arms.md), The run on the datasets of M). The run has not happened, so one table is here: the first run's step times, which the cost estimate reads.

| File | Holds |
|---|---|
| `first-run-times.tsv` | each step of the first run in the 9 species of M: species, step, whole seconds, and the sha256 of the species' `logs/timing.tsv` the row was read from |

`scripts/c3-inputs/cost.py times` wrote it from the 9 timing tables, each held to its own `SHA256SUMS` line, and `scripts/c3-inputs/cost.py measured` reads it by default. The first run's own tables are in [reference/input-arms](../input-arms/), whose `SHA256SUMS` does not list this file.

## After the run

`scripts/c3-inputs/m_page.py` copies here the 7 tables `m_reading.py` writes on the count machine, as they are, writes `SHA256SUMS` beside them, and writes nothing else here.

| File | Copied by | Holds |
|---|---|---|
| `m-floor.tsv` | `m_page.py floor`, before any score is read | each dataset of M's genomes by class, those listed and those left, and its second reading: in it, under the floor or a short list |
| `m-left-out.tsv` | `m_page.py floor` | each dataset left out of the second reading, with its reason and its counts by class |
| `m-count.tsv` | `m_page.py record`, after the count | the three counts in each of the two readings |
| `m-conditions.tsv` | `m_page.py record` | the conditions each count meets, in each reading |
| `m-both.tsv` | `m_page.py record` | per arm, conditions 1 to 5 met in both readings |
| `m-size.tsv` | `m_page.py record` | per reading, arm and learner: the datasets entered, better, worse and the same |
| `m-datasets.tsv` | `m_page.py record` | per reading, dataset, split, pair and learner; the second reading only for a dataset in it |
| `SHA256SUMS` | `m_page.py floor`, then `record` again | the sha256 of every file here but itself and this README |

`m_page.py print` prints every table from these files alone. After the floor `m_page.py left` writes the left-out page, `docs/input-arms-m-left-out.md`, and after the count `m_page.py page` writes the count page, `docs/input-arms-m-count.md`, each whole from these files; every other table stays here. `tests/c3_inputs/test_m_record.py` checks this directory against both pages, byte for byte, and every other page under `docs/` in each of its three states: before the run, after the floor and after the count.
