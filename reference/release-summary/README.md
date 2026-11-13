# Release summary reference outputs

These files are the outputs of the scripts in `scripts/release-summary/` on two public GRM releases. `docs/release-summary.md` says what each output holds and how each number is made. They were computed on 2026-11-13 at the commit that added them.

## Inputs

Each release is one archive. The scripts read two things from it: `index.json` and the folder `intervals/`, one file per dataset.

| Folder | Release asset | Asset sha256 | Evaluation revision | Datasets |
| --- | --- | --- | --- | --- |
| `rev12-274-20261105/` | `grm-release-rev12-274-20261105.tar.gz` | `cefaef3957f8c9e161b32e1989cd4561c9a314430fa7eeb75b97eefebbbbf5c4` | 12 | 274 |
| `rev14-274-20261112/` | `rev14-274-20261112.tar.zst` | `c67155f598dbb1ebe1a4ab908e793124ebc058430a33c1c67d7e2058910256de` | 14 | 274 |

The sha256 of the two inputs of each release:

| Folder | `index.json` | `intervals/` (274 files) |
| --- | --- | --- |
| `rev12-274-20261105/` | `3e467cf05c24ed5840a9c5b25e515cb75b7ab3262b0e6c67c5105ba4a4890bc6` | `545cbe2af908483d1c102c8707014a6c6a4d80bf6b4669580801d491ec253d0c` |
| `rev14-274-20261112/` | `1e549f017271b127633ef0f75c7eb9bd88c01c6b361a875f82fe0b4f73ff9fbb` | `cb8fb922e97892a4b3a839faf9abba9c38a3c5aab439c3bce0d0cc7904bbf59b` |

The `index.json` value is the line `./index.json` of the release's own checksum list. The `intervals/` value is one sha256 over the folder's `sha256  name` lines, sorted by name, as this command prints it from the release folder:

```sh
(cd intervals && sha256sum -- *.json | LC_ALL=C sort -k2 | sha256sum)
```

## Commands

Run from the repository root, with `R` the unpacked release folder (the one that holds `index.json`) and `OUT` the matching folder here. Python 3.12, standard library only.

```sh
python3 -I scripts/release-summary/uncertainty.py  "$R/index.json" "$R/intervals" > "$OUT/uncertainty.json"
python3 -I scripts/release-summary/summary.py      "$R/index.json" "$R/intervals" > "$OUT/summary.json"
python3 -I scripts/release-summary/catalog.py      "$R/index.json" "$R/intervals" > "$OUT/catalog.json"
python3 -I scripts/release-summary/index_counts.py "$R/index.json" "$R/intervals" > "$OUT/index-counts.json"
```

`uncertainty.py` uses its default of 4000 bootstrap repetitions and seed 20261111; both are in its output. To put the two releases side by side:

```sh
python3 -I scripts/release-summary/compare.py \
  reference/release-summary/rev12-274-20261105/uncertainty.json \
  reference/release-summary/rev14-274-20261112/uncertainty.json
```

## Outputs

The sha256 of each output, in the form `sha256sum -c` reads from this folder. `tests/release_summary` checks these lines against the files.

```
ed3fcce262ae58bc4446c79423a70aff27f2e19e1643bd2d5b8e06147ab406e5  rev12-274-20261105/uncertainty.json
c9de54d7987b891696ec5711bcaa4e677c11620ff729ffdb01a96cf1a9f0443d  rev12-274-20261105/summary.json
8982d32ea7a3028f734ca650d38cd92798b1441f2de9b891ae0e31398b51d55a  rev12-274-20261105/catalog.json
7ce1ea63fb4a89abc8d39ed1991e47c467bc14057e9341cf17e61f9612b6e5fe  rev12-274-20261105/index-counts.json
4f6b0d4181280c0a79ef5a18dde0350ee5d027ed40d0b9580cb935426af7cd03  rev14-274-20261112/uncertainty.json
7dfc2108b8e8853c438cd2e475703cd185c748bd7d03186c2aff71aced3282db  rev14-274-20261112/summary.json
2d74fab696e2c7a5bb0814ac86082116d2c787f0440cdf5ba872b302532da6bd  rev14-274-20261112/catalog.json
b0352b81ddf10dc3f5abaaea39baa1870c3bc54f93f8972603728e30b1d17d52  rev14-274-20261112/index-counts.json
```

## Against the earlier runs

The same scripts ran from outside the repository before they moved here. Against those runs:

- `summary.json` is byte-identical for both releases.
- `uncertainty.json` differs only in the lower end of some bootstrap intervals, each by 0.0001 to 0.001: 15 intervals of revision 12 and 14 of revision 14. The earlier runs took the lower end one order statistic too low (`docs/release-summary.md`, Changes against the earlier runs).
- `catalog.json` differs in two counts, the same in both releases: `rule.calls_nothing_resistant` 32 is now 29 and `rule.calls_everything_resistant` 31 is now 30. The earlier counts included four rule entries that are not constant calls.
- `index-counts.json` comes from one script where there were two, so its layout is new. Every count the earlier files hold is unchanged. Both files add `revision`, `rule` and `datasets`, and the revision-12 file adds the flag counts the earlier revision-12 counter did not print. Both releases have 0 mismatches.
