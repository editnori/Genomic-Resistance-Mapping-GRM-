"""The results one species' C7 run must hold (docs/outside-learners.md), as the list score.py reads.

  python3 -I expect.py <species_slug> <datasets> <study datasets> <seeds> <twin seeds> <arm seeds> [<arms>]

<datasets> and <study datasets> are comma lists ("-" or empty for none); the
study datasets must be listed datasets. <seeds> are the default's draws (they
must hold 42), <twin seeds> the twin's, <arm seeds> the library arms' (a
subset of <seeds> holding 42). <arms> is the library arms that run, a comma
list: all three (the default), or the two boosters when the elastic net is
left out, which then has no row.

Each dataset has, for the default (k31), a lineage result at every seed of
<seeds> (role draw), for each library arm (xgboost, lightgbm, elasticnet) a
lineage result at every seed of <arm seeds> (role draw), and, if it is a
study dataset, one study result at seed 42 for the default and for each arm
(role study); the twin is the default at every seed of <twin seeds> (role
twin), lineage only.

Columns: dataset, species_slug, arm, split, scheme, seed, role, file. The
file is relative to the run's folder: the default's and the twin's are the
C3 stage's results, c3/results/picks/<arm>/<dataset>.<scheme>-s<seed>.pick.json;
an arm's is results/outside/<arm>/<dataset>.<scheme>-s<seed>.outside.json.
"""

import sys

RELEASE_SEED = 42
STUDY_SCHEME = "study-whole"
DEFAULT = "k31"
TWIN = "twin"
ARMS = ("xgboost", "lightgbm", "elasticnet")
COLUMNS = ["dataset", "species_slug", "arm", "split", "scheme", "seed", "role", "file"]


def comma(text):
    return [] if text in ("", "-") else text.split(",")


def picked(arm, d, scheme, s):
    return f"c3/results/picks/{arm}/{d}.{scheme}-s{s}.pick.json"


def outside(arm, d, scheme, s):
    return f"results/outside/{arm}/{d}.{scheme}-s{s}.outside.json"


def rows(species, datasets, study, seeds, twin, arm_seeds, arms=ARMS):
    for d in datasets + study:
        if not d.startswith(species + "__"):
            raise ValueError(f"{d} is not a dataset of {species}")
    for name, items in (("datasets", datasets), ("study datasets", study), ("seeds", seeds), ("twin seeds", twin),
                        ("arm seeds", arm_seeds)):
        if len(set(items)) != len(items):
            raise ValueError(f"the {name} list names one twice")
    if not datasets:
        raise ValueError("no dataset listed")
    if not set(study) <= set(datasets):
        raise ValueError("a study dataset is not a listed dataset")
    if set(seeds) & set(twin):
        raise ValueError("a seed is both a draw's and the twin's")
    if RELEASE_SEED not in seeds or RELEASE_SEED not in arm_seeds:
        raise ValueError("the draws leave out the release's seed")
    if not set(arm_seeds) <= set(seeds):
        raise ValueError("an arm seed is not a seed of the default's draws")
    if tuple(arms) not in (ARMS, ARMS[:2]):
        raise ValueError("the arms are all three, or the two boosters with the elastic net left out")
    out = []
    for d in sorted(datasets):
        for s in seeds:
            out.append([d, species, DEFAULT, "lineage", "lineage", s, "draw", picked(DEFAULT, d, "lineage", s)])
        for arm in arms:
            for s in arm_seeds:
                out.append([d, species, arm, "lineage", "lineage", s, "draw", outside(arm, d, "lineage", s)])
        if d in study:
            out.append([d, species, DEFAULT, "study", STUDY_SCHEME, RELEASE_SEED, "study",
                        picked(DEFAULT, d, STUDY_SCHEME, RELEASE_SEED)])
            for arm in arms:
                out.append([d, species, arm, "study", STUDY_SCHEME, RELEASE_SEED, "study",
                            outside(arm, d, STUDY_SCHEME, RELEASE_SEED)])
        for s in twin:
            out.append([d, species, TWIN, "lineage", "lineage", s, "twin", picked(TWIN, d, "lineage", s)])
    return out


def main(argv):
    if len(argv) not in (6, 7):
        print(__doc__.strip().splitlines()[2], file=sys.stderr)
        return 2
    seeds, twin, arm_seeds = ([int(s) for s in comma(a)] for a in argv[3:6])
    arms = tuple(comma(argv[6])) if len(argv) == 7 else ARMS
    print("\t".join(COLUMNS))
    for r in rows(argv[0], comma(argv[1]), comma(argv[2]), seeds, twin, arm_seeds, arms):
        print("\t".join(str(x) for x in r))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
