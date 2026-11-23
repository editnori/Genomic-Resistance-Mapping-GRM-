"""The results one species' C4 run must hold (docs/screen-arms.md), as the list score.py reads.

  python3 -I expect.py <species_slug> <datasets> <study datasets> <seeds> <twin seeds> [<arms>]

<datasets> and <study datasets> are comma lists ("-" or empty for none);
the study datasets must be listed datasets. <arms> is a comma list of the
arms the run fits beside the default and the twin (default:
cmh-4000,cmh-16000,floor-3, every arm of the block).

Each dataset has, for the default (k31) and each arm, a lineage result at
every seed of <seeds> (role draw) and, if it is a study dataset, one study
result at seed 42 (role study), as grm learn --fold-file runs at the
release's seeds; the twin is the default at every seed of <twin seeds>
(role twin), lineage only. Every result is fitted with no inner fold.

Columns are scripts/c3-inputs/expect.py's: dataset, species_slug, arm,
split, scheme, seed, role, file (the result under the picks folder,
<arm>/<name> as grm learn --exp-pick writes it) and inputs (the arm's
counts, <arm>/<dataset>.<tag>-s<seed>.inputs.json with tag drawn or the
study scheme, or "-" for the default and the twin, which write none).
"""

import importlib.util
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent


def _c3():
    spec = importlib.util.spec_from_file_location("c3_expect_for_c4", HERE.parent / "c3-inputs" / "expect.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


C3 = _c3()
RELEASE_SEED = C3.RELEASE_SEED
STUDY_SCHEME = C3.STUDY_SCHEME
DEFAULT = C3.DEFAULT
TWIN = C3.TWIN
COLUMNS = C3.COLUMNS
comma = C3.comma
# Each arm's grm learn flags beside --exp-pick DIR --exp-pick-inner 0.
ARM_FLAGS = {
    "cmh-4000": ["--exp-screen", "cmh", "--force"],
    "cmh-16000": ["--exp-screen", "cmh", "--screened-features", "16000", "--force"],
    "floor-3": ["--exp-floor", "3", "--force"],
}
ARMS = tuple(ARM_FLAGS)


def inputs_name(arm, dataset, tag, seed):
    return f"{arm}/{dataset}.{tag}-s{seed}.inputs.json"


def rows(species, datasets, study, seeds, twin, arms=ARMS):
    for d in datasets + study:
        if not d.startswith(species + "__"):
            raise ValueError(f"{d} is not a dataset of {species}")
    for name, items in (("datasets", datasets), ("study datasets", study), ("seeds", seeds),
                        ("twin seeds", twin), ("arms", arms)):
        if len(set(items)) != len(items):
            raise ValueError(f"the {name} list names one twice")
    if not datasets:
        raise ValueError("no dataset listed")
    if not set(study) <= set(datasets):
        raise ValueError("a study dataset is not a listed dataset")
    if set(seeds) & set(twin):
        raise ValueError("a seed is both a draw's and the twin's")
    if RELEASE_SEED not in seeds:
        raise ValueError("the draws leave out the release's seed")
    unknown = [x for x in arms if x not in ARMS]
    if unknown:
        raise ValueError(f"no such arm: {','.join(unknown)}")
    out = []
    for d in sorted(datasets):
        for arm in (DEFAULT,) + tuple(arms):
            inputs = arm != DEFAULT
            for s in seeds:
                out.append([d, species, arm, "lineage", "lineage", s, "draw", f"{arm}/{d}.lineage-s{s}.pick.json",
                            inputs_name(arm, d, "drawn", s) if inputs else "-"])
            if d in study:
                out.append([d, species, arm, "study", STUDY_SCHEME, RELEASE_SEED, "study",
                            f"{arm}/{d}.{STUDY_SCHEME}-s{RELEASE_SEED}.pick.json",
                            inputs_name(arm, d, STUDY_SCHEME, RELEASE_SEED) if inputs else "-"])
        for s in twin:
            out.append([d, species, TWIN, "lineage", "lineage", s, "twin", f"{TWIN}/{d}.lineage-s{s}.pick.json", "-"])
    return out


def main(argv):
    if len(argv) not in (5, 6):
        print(__doc__.strip().splitlines()[2], file=sys.stderr)
        return 2
    seeds = [int(s) for s in comma(argv[3])]
    twin = [int(s) for s in comma(argv[4])]
    arms = comma(argv[5]) if len(argv) == 6 else list(ARMS)
    print("\t".join(COLUMNS))
    for r in rows(argv[0], comma(argv[1]), comma(argv[2]), seeds, twin, arms):
        print("\t".join(str(x) for x in r))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
