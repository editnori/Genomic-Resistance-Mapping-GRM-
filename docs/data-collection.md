# Data collection: one checked copy of the project's files

`scripts/collect/collect.py` lists files with their size, sha256 and class, packs one class into archives, and checks an archive against the list it carries. It uploads nothing and needs only Python's standard library.

## Classes

| Class | Packed | Holds |
|---|---|---|
| `public-result` | yes | published release files |
| `internal-record` | yes | messages, notes, ledgers, reviews, logs |
| `public-input` | yes | public genomes and tables |
| `regenerable` | no, listed with its size | builds, caches, environments, an unpacked copy of a listed archive |
| `refused` | never | any path holding `sealed`, `held-out`, `ks30`, `ks50`, `vumc`, `clinical`, `private`, `ksd` or `sigdrug` |

A refused path stops the command: it is named, and nothing in it is opened. An archive holds one class, so internal records never share a file with public results. A `public-result` file is packed only if a published checksum file also names its sha256 (`--published`).

## Collect

Write the rules, one `CLASS<TAB>PATTERN` per line. The first match wins, and a file that matches no rule stops the list. A pattern ending in `/` names a `regenerable` directory. `node_modules`, `__pycache__` and a directory holding `CACHEDIR.TAG` or `pyvenv.cfg` need no rule.

```tsv
public-result	records/release/*.json
internal-record	*
```

```sh
C="python3 -I scripts/collect/collect.py"
$C list --rules rules.tsv --root records=/path/to/records --out list.tsv
$C pack --list list.tsv --class internal-record --name grm-20261114 --root records=/path/to/records --out archive
$C pack --list list.tsv --class public-result --name grm-20261114 --root records=/path/to/records --out archive \
    --published /path/to/records/release/SHA256SUMS
$C summary --place this-machine=list.tsv
```

`archive/` holds `grm-20261114.<class>.partNNN.tar.gz`, each at most 1.9 GB (`--max-bytes`), and `grm-20261114.<class>.parts.sha256`. Every part carries the class's whole list.

`summary` takes one `--place NAME=LIST` per machine. For each class it prints the bytes, the bytes after duplicates are removed, and what only one place holds.

## Restore check

Check the parts without unpacking:

```sh
$C check --parts archive/grm-20261114.internal-record.parts.sha256
```

Or unpack into a new, empty folder and check that:

```sh
mkdir restore
tar -xzf archive/grm-20261114.internal-record.part001.tar.gz -C restore
$C check restore/grm-20261114.internal-record
```

`check` exits 0 only if every listed file is there with its size and sha256, and nothing else is. Otherwise it exits 1 and names each missing, changed or added file.

## Tests

`tests/collect` runs in CI on a synthetic project. `python3 -I tests/collect/mutants.py` plants faults in the tool and prints the tests that catch each.
