# Data collection: one checked copy of the project's files

`scripts/collect/collect.py` lists files with their size, sha256 and class, packs each class into its own archives, and checks an archive against the list it carries. It uploads nothing and needs only Python's standard library.

## Classes

| Class | Packed | Holds |
|---|---|---|
| `public-result` | yes | published release files |
| `internal-record` | yes | notes, messages, ledgers, logs |
| `public-input` | yes | public genomes and tables |
| `regenerable` | no, listed with its size | builds, caches, environments |
| `refused` | never | any path holding `sealed`, `held-out`, `ks30`, `ks50`, `vumc`, `clinical`, `private`, `ksd` or `sigdrug` |

A path holding a barred word is never opened and stops the command, or with `list --record-refused` becomes a row in no archive: look at it by hand. The guard reads names only: a hard link with a clean name is read. A `public-result` file is packed only if `--published`, the checksum file of a release already published, lists it by name and content.

## Collect

Rules are `CLASS<TAB>PATTERN` lines: the first match wins, and a file no rule matches stops the list. A pattern ending in `/` names a `regenerable` directory. `node_modules`, `__pycache__` and a directory holding `CACHEDIR.TAG` or `pyvenv.cfg` need no rule.

```tsv
public-result	records/release/*.json
internal-record	*
```

```sh
C="python3 -I scripts/collect/collect.py"
$C list --rules rules.tsv --root records=/path/to/records --out list.tsv
$C pack --list list.tsv --class internal-record --name grm-20261114 --root records=/path/to/records --out archive/internal
$C pack --list list.tsv --class public-result --name grm-20261114 --root records=/path/to/records --out archive/public \
    --published /path/to/release/SHA256SUMS
$C summary --place this-machine=list.tsv
```

Each class has its own folder: parts, each at most 1.9 GB (`--max-bytes`) and carrying the whole list, and a `parts.sha256`. Send `archive/public` alone. Parts hold file times, so a second pack gives other bytes; the list repeats.

`summary` counts copies by sha256, one `--place` per machine.

## Restore check

Check the parts first:

```sh
$C check --parts archive/internal/grm-20261114.internal-record.parts.sha256
```

It trusts the `parts.sha256` beside them, so keep a copy apart. Then unpack every part:

```sh
mkdir restore
for p in archive/internal/grm-20261114.internal-record.part*.tar.gz; do tar -xzf "$p" -C restore; done
$C check restore/grm-20261114.internal-record
```

`check DIR` proves only that the folder matches the list inside it. `check` exits 0 only if every listed file is there with its size and sha256 and nothing else is; else 1, naming each difference.
