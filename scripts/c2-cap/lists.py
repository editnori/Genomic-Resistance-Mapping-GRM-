"""The dataset lists of plan block C2 (docs/cap-arms.md, Datasets).

  python3 -I lists.py <c3 lists.tsv> <release dir> <out dir>

<c3 lists.tsv> is L as scripts/c3-inputs/lists.py wrote it: the release
reports with a lineage figure, with their species and antibiotic in full,
M, S, the first run (L less M) and the label rule each one's measures depend
on (the label count's mark). <release dir> is the release
(rev14-274-20261112), whose species/<slug>.json lists each species' cohort.
CAP is L in the species whose cohort reached the cap of 2,000 genomes.

Writes <out dir>/lists.tsv (CAP's rows: dataset, species in full,
antibiotic in full, the report's counts, in_M, in_S, first_run, the label
mark) and <out dir>/species.tsv (per species: CAP, M, the first run, and
the first run in S), and prints the totals. The first run is CAP less M;
the study split runs on the first run's datasets in S.
"""

import csv
import json
import sys
from pathlib import Path

CAP = 2000
COLUMNS = ["dataset", "species_slug", "species", "antibiotic", "samples", "resistant", "susceptible", "lineages",
           "in_M", "in_S", "first_run", "labels_depend_on"]
SPECIES_COLUMNS = ["species_slug", "species", "cohort", "CAP", "M", "first_run", "first_run_in_S"]


def read(path):
    with open(path, encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write(path, columns, rows):
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, columns, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def cohort_sizes(release, slugs):
    out = {}
    for slug in sorted(slugs):
        with open(Path(release) / "species" / f"{slug}.json", encoding="utf-8") as f:
            out[slug] = len(json.load(f)["genomes"])
    return out


def lists(l_rows, sizes):
    """CAP's rows and the per-species rows."""
    cap = [r for r in l_rows if sizes[r["species_slug"]] >= CAP]
    per = []
    for slug in sorted({r["species_slug"] for r in cap}):
        rows = [r for r in cap if r["species_slug"] == slug]
        per.append({"species_slug": slug, "species": rows[0]["species"], "cohort": sizes[slug], "CAP": len(rows),
                    "M": sum(r["in_M"] == "1" for r in rows),
                    "first_run": sum(r["first_run"] == "1" for r in rows),
                    "first_run_in_S": sum(r["first_run"] == "1" and r["in_S"] == "1" for r in rows)})
    for r in cap:
        if (r["first_run"] == "1") == (r["in_M"] == "1"):
            raise SystemExit(f"{r['dataset']}: the first run must be exactly the datasets not in M")
    return cap, per


def main(argv):
    if len(argv) != 4:
        raise SystemExit(__doc__)
    l_rows = read(argv[1])
    sizes = cohort_sizes(argv[2], {r["species_slug"] for r in l_rows})
    if any(n > CAP for n in sizes.values()):
        raise SystemExit("a release cohort is past the cap of 2,000")
    cap, per = lists(l_rows, sizes)
    out = Path(argv[3])
    out.mkdir(parents=True, exist_ok=True)
    write(out / "lists.tsv", COLUMNS, cap)
    write(out / "species.tsv", SPECIES_COLUMNS, per)
    print(f"CAP {len(cap)} in {len(per)} species, M {sum(p['M'] for p in per)}, "
          f"first run {sum(p['first_run'] for p in per)}, first run in S {sum(p['first_run_in_S'] for p in per)}")


if __name__ == "__main__":
    main(sys.argv)
