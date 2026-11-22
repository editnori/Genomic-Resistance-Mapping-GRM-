"""The results one species' C9 run must hold (docs/served-pick.md), as the list score.py reads.

  python3 -I expect.py <species_slug> <datasets> <study datasets> <seeds> <twin seeds> > expect.tsv

<datasets> and <study datasets> are comma lists (the second may be empty or
"-"); every study dataset must be a listed dataset. Each dataset has a
lineage result at every seed of <seeds> (role draw, 3 inner folds) and of
<twin seeds> (role twin, no inner fold); each study dataset has one study
result at seed 42 (role study, 3 inner folds), as grm learn --fold-file runs
at the release's seeds.

Columns: dataset, species_slug, split, scheme, seed, inner, role, file (the
result's name, as grm learn --exp-pick writes it).
"""

import sys

RELEASE_SEED = 42
INNER = 3
STUDY_SCHEME = "study-whole"
COLUMNS = ["dataset", "species_slug", "split", "scheme", "seed", "inner", "role", "file"]


def comma(text):
    return [] if text in ("", "-") else text.split(",")


def rows(species, datasets, study, seeds, twin):
    for d in datasets + study:
        if not d.startswith(species + "__"):
            raise ValueError(f"{d} is not a dataset of {species}")
    for name, items in (("datasets", datasets), ("study datasets", study), ("seeds", seeds), ("twin seeds", twin)):
        if len(set(items)) != len(items):
            raise ValueError(f"the {name} list names one twice")
    if not datasets:
        raise ValueError("no dataset listed")
    if not set(study) <= set(datasets):
        raise ValueError("a study dataset is not a listed dataset")
    if set(seeds) & set(twin):
        raise ValueError("a seed is both a draw's and the twin's")
    out = []
    for d in sorted(datasets):
        for role, inner, ss in (("draw", INNER, seeds), ("twin", 0, twin)):
            for s in ss:
                out.append([d, species, "lineage", "lineage", s, inner, role, f"{d}.lineage-s{s}.pick.json"])
        if d in study:
            out.append([d, species, "study", STUDY_SCHEME, RELEASE_SEED, INNER, "study",
                        f"{d}.{STUDY_SCHEME}-s{RELEASE_SEED}.pick.json"])
    return out


def main(argv):
    if len(argv) != 5:
        print(__doc__.strip().splitlines()[2], file=sys.stderr)
        return 2
    seeds = [int(s) for s in comma(argv[3])]
    twin = [int(s) for s in comma(argv[4])]
    print("\t".join(COLUMNS))
    for r in rows(argv[0], comma(argv[1]), comma(argv[2]), seeds, twin):
        print("\t".join(str(x) for x in r))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
