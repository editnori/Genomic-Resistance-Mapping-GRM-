"""The results one species' C5 run must hold (docs/settings-arms.md), as the list score.py reads.

  python3 -I expect.py <species_slug> <datasets> <study datasets> <seeds> <twin seeds>

<datasets> and <study datasets> are comma lists ("-" or empty for none); the
study datasets must be listed datasets. The default and the twin are block
C3's rows (scripts/c3-inputs/expect.py, with no arm beside them), with the
default named `default`: per dataset a lineage result at every seed of
<seeds> (role draw), one study result at seed 42 (role study) for a study
dataset, and the twin, the default at every seed of <twin seeds> (role twin,
lineage only). The search and the cutoff arm stand beside each draw and
study result of the default and of no twin: they share the driver's choices
files, one per outer fold, choices/<dataset>.<scheme>-s<seed>.f<k>.choices.json
for k from 0 to FOLDS - 1, and their file column is that name with <k>
standing for the fold.

Columns: dataset, species_slug, arm, split, scheme, seed, role, file.
"""

import importlib.util
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
DEFAULT = "default"
TWIN = "twin"
ARMS = ("search", "cutoff")
FOLDS = 5
COLUMNS = ["dataset", "species_slug", "arm", "split", "scheme", "seed", "role", "file"]


def _c3():
    spec = importlib.util.spec_from_file_location("c3_expect_for_c5", HERE.parent / "c3-inputs" / "expect.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


C3 = _c3()
comma = C3.comma


def choices_name(dataset, scheme, seed, fold="<k>"):
    return f"choices/{dataset}.{scheme}-s{seed}.f{fold}.choices.json"


def rows(species, datasets, study, seeds, twin):
    # C3's rows check the lists (species, no name twice, study datasets listed, the release's seed among the
    # draws, no seed both a draw's and the twin's) and give the default's and the twin's rows.
    base = C3.rows(species, datasets, [], study, seeds, twin, ())
    out = []
    for d, sp, arm, split, scheme, seed, role, name, _ in base:
        if arm == C3.DEFAULT:
            arm, name = DEFAULT, DEFAULT + name[len(C3.DEFAULT):]
        out.append([d, sp, arm, split, scheme, seed, role, name])
        if arm == DEFAULT and role in ("draw", "study"):
            for a in ARMS:
                out.append([d, sp, a, split, scheme, seed, role, choices_name(d, scheme, seed)])
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
