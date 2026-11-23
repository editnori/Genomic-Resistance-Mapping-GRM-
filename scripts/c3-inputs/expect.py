"""The results one species' C3 run must hold (docs/input-arms.md), as the list score.py reads.

  python3 -I expect.py <species_slug> <datasets> <A datasets> <study datasets> <seeds> <twin seeds> [<arms>]

<datasets>, <A datasets> and <study datasets> are comma lists ("-" or empty
for none); the A and study datasets must be listed datasets. <arms> is a
comma list of the arms the run fits beside the default and the twin
(default: program,k31-program, the arms built so far); the program arms run
on the A datasets only.

Each dataset has, for the default (k31) and each arm it runs, a lineage
result at every seed of <seeds> (role draw) and, if it is a study dataset,
one study result at seed 42 (role study), as grm learn --fold-file runs at
the release's seeds; the twin is the default at every seed of <twin seeds>
(role twin), lineage only. Every result is fitted with no inner fold.

Columns: dataset, species_slug, arm, split, scheme, seed, role, file (the
result under the picks folder, <arm>/<name> as grm learn --exp-pick writes
it) and inputs (the arm's column counts, <arm>/<dataset>.<tag>.inputs.json,
or "-" for the default and the twin, which write none). An arm's study
result that grm learn refused for no candidate stands as the refusal record
run_species.sh writes beside it, <arm>/<dataset>.study-whole-s42.refused.json,
which score.py scores as the constant model (ruling 16).
"""

import sys

RELEASE_SEED = 42
STUDY_SCHEME = "study-whole"
DEFAULT = "k31"
TWIN = "twin"
L_ARMS = ("k31-counts", "unitigs", "k31-aa")
A_ARMS = ("program", "k31-program")
BUILT = ("program", "k31-program")
COLUMNS = ["dataset", "species_slug", "arm", "split", "scheme", "seed", "role", "file", "inputs"]


def comma(text):
    return [] if text in ("", "-") else text.split(",")


def rows(species, datasets, a, study, seeds, twin, arms=BUILT):
    for d in datasets + a + study:
        if not d.startswith(species + "__"):
            raise ValueError(f"{d} is not a dataset of {species}")
    for name, items in (("datasets", datasets), ("A datasets", a), ("study datasets", study), ("seeds", seeds),
                        ("twin seeds", twin), ("arms", arms)):
        if len(set(items)) != len(items):
            raise ValueError(f"the {name} list names one twice")
    if not datasets:
        raise ValueError("no dataset listed")
    if not set(a) <= set(datasets) or not set(study) <= set(datasets):
        raise ValueError("an A or study dataset is not a listed dataset")
    if set(seeds) & set(twin):
        raise ValueError("a seed is both a draw's and the twin's")
    if RELEASE_SEED not in seeds:
        raise ValueError("the draws leave out the release's seed")
    unknown = [x for x in arms if x not in L_ARMS + A_ARMS]
    if unknown:
        raise ValueError(f"no such arm: {','.join(unknown)}")
    out = []
    for d in sorted(datasets):
        for arm in (DEFAULT,) + tuple(x for x in arms if x in L_ARMS or d in a):
            inputs = arm != DEFAULT
            for s in seeds:
                out.append([d, species, arm, "lineage", "lineage", s, "draw", f"{arm}/{d}.lineage-s{s}.pick.json",
                            f"{arm}/{d}.drawn.inputs.json" if inputs else "-"])
            if d in study:
                out.append([d, species, arm, "study", STUDY_SCHEME, RELEASE_SEED, "study",
                            f"{arm}/{d}.{STUDY_SCHEME}-s{RELEASE_SEED}.pick.json",
                            f"{arm}/{d}.{STUDY_SCHEME}.inputs.json" if inputs else "-"])
        for s in twin:
            out.append([d, species, TWIN, "lineage", "lineage", s, "twin", f"{TWIN}/{d}.lineage-s{s}.pick.json", "-"])
    return out


def main(argv):
    if len(argv) not in (6, 7):
        print(__doc__.strip().splitlines()[2], file=sys.stderr)
        return 2
    seeds = [int(s) for s in comma(argv[4])]
    twin = [int(s) for s in comma(argv[5])]
    arms = comma(argv[6]) if len(argv) == 7 else list(BUILT)
    print("\t".join(COLUMNS))
    for r in rows(argv[0], comma(argv[1]), comma(argv[2]), comma(argv[3]), seeds, twin, arms):
        print("\t".join(str(x) for x in r))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
