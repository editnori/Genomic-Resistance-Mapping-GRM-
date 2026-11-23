"""The dataset lists of plan block C7 (docs/outside-learners.md, Datasets).

  python3 -I lists.py <c3 lists.tsv> <out dir>

<c3 lists.tsv> is plan block C3's list (scripts/c3-inputs/lists.tsv): the
release reports with a lineage figure (L), with their species and antibiotic
in full, M, S, the first run and the label rules each one's measures depend
on. C7 runs on L, as C3's L arms do; its first run is L less M, and the
study split runs on the first run's datasets in S.

Writes <out dir>/lists.tsv (L's rows with the arms and the splits each
dataset runs) and <out dir>/species.tsv (per species, in full: L, L less M,
and the first run's datasets in S), and prints the totals.
"""

import csv
import sys
from pathlib import Path

ARMS = ["k31", "twin", "xgboost", "lightgbm", "elasticnet"]
COLUMNS = ["dataset", "species_slug", "species", "antibiotic", "samples", "resistant", "susceptible", "lineages",
           "in_M", "in_S", "first_run", "labels_depend_on", "arms", "splits"]
SPECIES_COLUMNS = ["species_slug", "species", "L", "L_not_M", "L_not_M_in_S"]


def read(path):
    with open(path, encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def lists(c3_rows):
    seen = set()
    out = []
    for r in c3_rows:
        if r["dataset"] in seen:
            raise SystemExit(f"{r['dataset']} is listed twice")
        seen.add(r["dataset"])
        if r["first_run"] == "1" and r["in_M"] == "1":
            raise SystemExit(f"{r['dataset']} is in M and in the first run")
        row = {k: r[k] for k in COLUMNS if k in r}
        row["arms"] = ",".join(ARMS)
        row["splits"] = "lineage,study" if r["first_run"] == "1" and r["in_S"] == "1" else "lineage"
        out.append(row)
    return out


def species_rows(rows):
    out = []
    for sp in sorted({r["species_slug"] for r in rows}):
        rs = [r for r in rows if r["species_slug"] == sp]
        first = [r for r in rs if r["first_run"] == "1"]
        out.append({"species_slug": sp, "species": rs[0]["species"], "L": len(rs), "L_not_M": len(first),
                    "L_not_M_in_S": sum(r["in_S"] == "1" for r in first)})
    return out


def write(path, rows, columns):
    with open(path, "w", encoding="utf-8", newline="") as f:
        f.write("\t".join(columns) + "\n")
        for r in rows:
            f.write("\t".join(str(r[c]) for c in columns) + "\n")


def main(argv):
    if len(argv) != 2:
        print(__doc__.strip().splitlines()[2], file=sys.stderr)
        return 2
    rows = lists(read(argv[0]))
    out = Path(argv[1])
    out.mkdir(parents=True, exist_ok=True)
    write(out / "lists.tsv", rows, COLUMNS)
    sp = species_rows(rows)
    write(out / "species.tsv", sp, SPECIES_COLUMNS)
    first_sp = sum(r["L_not_M"] > 0 for r in sp)
    print(f"L {sum(r['L'] for r in sp)} in {len(sp)} species; L less M {sum(r['L_not_M'] for r in sp)} in {first_sp} "
          f"species, {sum(r['L_not_M_in_S'] for r in sp)} in S")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
