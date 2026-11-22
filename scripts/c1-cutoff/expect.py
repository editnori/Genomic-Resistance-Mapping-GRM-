"""The results one species' cohort run must hold (docs/cohort-draw.md), as the list score.py reads.

  python3 -I expect.py <species_slug> <draws.tsv> <twin> > expect.tsv

<draws.tsv> is draw.py's list. Every cohort file it wrote is fitted at seed
42, the release's; with <twin> 1, the cohorts at r 25 are fitted again at
seed 45. A refused draw has no result. Every dataset must begin with
<species_slug>__.

Columns: dataset, species_slug, r, draw, seed, study, cohort (the cohort
file's name), cohort_sha256, file (the result's name under the results
directory, as grm learn --exp-cohort writes it).
"""

import csv
import sys

SEED = 42
TWIN_SEED = 45
TWIN_R = 25
COLUMNS = ["dataset", "species_slug", "r", "draw", "seed", "study", "cohort", "cohort_sha256", "file"]


def rows(species, draws, twin):
    out = []
    for d in draws:
        if not d["dataset"].startswith(species + "__"):
            raise ValueError(f"{d['dataset']} is not a dataset of {species}")
        if d["written"] != "1":
            continue
        r = int(d["r"])
        seeds = [SEED, TWIN_SEED] if twin and r == TWIN_R else [SEED]
        for seed in seeds:
            out.append({
                "dataset": d["dataset"], "species_slug": species, "r": r, "draw": int(d["draw"]),
                "seed": seed, "study": d["study"], "cohort": d["file"], "cohort_sha256": d["sha256"],
                "file": f"{d['dataset']}.r{r}.d{d['draw']}.s{seed}.cohort.json",
            })
    names = [r["file"] for r in out]
    if len(set(names)) != len(names):
        raise ValueError("two results share a name")
    return out


def main(argv):
    if len(argv) != 3 or argv[2] not in ("0", "1"):
        print(__doc__.strip().splitlines()[2], file=sys.stderr)
        return 2
    with open(argv[1], newline="") as f:
        draws = list(csv.DictReader(f, delimiter="\t"))
    table = rows(argv[0], draws, argv[2] == "1")
    if not table:
        print("no cohort to fit", file=sys.stderr)
        return 2
    sys.stdout.write("\t".join(COLUMNS) + "\n")
    for r in table:
        sys.stdout.write("\t".join(str(r[c]) for c in COLUMNS) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
