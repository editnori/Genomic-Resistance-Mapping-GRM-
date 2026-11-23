"""The dataset lists of plan block C3 (docs/input-arms.md, Datasets).

  python3 -I lists.py <l.tsv> <amrfinderplus datasets.tsv> <out dir>

<l.tsv> is the list of L as scripts/c9-pick/lists.py wrote it
(scripts/c9-pick/l.tsv): the release reports with a lineage figure, with
their species and antibiotic in full, M, S, the first run and the label
rules each one's measures depend on. A is the rule
scripts/assumption-d0-c8/miss_sets.py applies to the AMRFinderPlus
comparison's record (reference/amrfinderplus-comparison/datasets.tsv): a
dataset with an organism option and a lineage split. Every A dataset must be
in L.

Writes <out dir>/lists.tsv (L's rows with in_A, first_run_A and the arms each
dataset runs) and <out dir>/species.tsv (per species, in full: L, L less M,
A, A less M, and the first runs in S), and prints the totals. The first run
of an arm is its set less M. The study split runs on the first run's
datasets in S.
"""

import csv
import sys
from pathlib import Path

L_ARMS = ["k31", "twin", "k31-counts", "unitigs", "k31-aa"]
A_ARMS = ["program", "k31-program"]
COLUMNS = ["dataset", "species_slug", "species", "antibiotic", "samples", "resistant", "susceptible", "lineages",
           "in_M", "in_S", "in_A", "first_run", "first_run_A", "pairs_S+R", "labels_depend_on", "arms"]
SPECIES_COLUMNS = ["species_slug", "species", "L", "L_not_M", "L_not_M_in_S", "A", "A_not_M", "A_not_M_in_S",
                   "L_not_M_with_S_and_R_pair", "A_not_M_with_S_and_R_pair"]


def read(path):
    with open(path, encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def a_set(rows):
    return {r["dataset"] for r in rows if r["organism_option"] == "yes" and r["lineage_split"] == "yes"}


def lists(l_rows, a):
    by = {r["dataset"]: r for r in l_rows}
    outside = sorted(a - set(by))
    if outside:
        raise SystemExit(f"{len(outside)} A datasets are not in L: {','.join(outside)}")
    out = []
    for r in l_rows:
        in_a = r["dataset"] in a
        row = {k: r[k] for k in COLUMNS if k in r}
        row["in_A"] = int(in_a)
        row["first_run_A"] = int(in_a and r["first_run"] == "1")
        arms = list(L_ARMS) + (A_ARMS if in_a else [])
        row["arms"] = ",".join(arms)
        out.append(row)
    return out


def species_rows(rows):
    out = []
    for sp in sorted({r["species_slug"] for r in rows}):
        rs = [r for r in rows if r["species_slug"] == sp]
        first = [r for r in rs if r["first_run"] == "1"]
        first_a = [r for r in rs if r["first_run_A"] == 1]
        out.append({
            "species_slug": sp, "species": rs[0]["species"], "L": len(rs), "L_not_M": len(first),
            "L_not_M_in_S": sum(r["in_S"] == "1" for r in first), "A": sum(r["in_A"] for r in rs),
            "A_not_M": len(first_a), "A_not_M_in_S": sum(r["in_S"] == "1" for r in first_a),
            "L_not_M_with_S_and_R_pair": sum(int(r["pairs_S+R"]) > 0 for r in first),
            "A_not_M_with_S_and_R_pair": sum(int(r["pairs_S+R"]) > 0 for r in first_a),
        })
    return out


def write(path, rows, columns):
    with open(path, "w", encoding="utf-8", newline="") as f:
        f.write("\t".join(columns) + "\n")
        for r in rows:
            f.write("\t".join(str(r[c]) for c in columns) + "\n")


def main(argv):
    if len(argv) != 3:
        print(__doc__.strip().splitlines()[2], file=sys.stderr)
        return 2
    rows = lists(read(argv[0]), a_set(read(argv[1])))
    out = Path(argv[2])
    out.mkdir(parents=True, exist_ok=True)
    write(out / "lists.tsv", rows, COLUMNS)
    sp = species_rows(rows)
    write(out / "species.tsv", sp, SPECIES_COLUMNS)
    tot = {c: sum(r[c] for r in sp) for c in SPECIES_COLUMNS[2:]}
    print(f"L {tot['L']} in {len(sp)} species; L less M {tot['L_not_M']}, {tot['L_not_M_in_S']} in S; "
          f"A {tot['A']} in {sum(r['A'] > 0 for r in sp)} species; A less M {tot['A_not_M']}, {tot['A_not_M_in_S']} in S; "
          f"with an S and R pair: {tot['L_not_M_with_S_and_R_pair']} of L less M, {tot['A_not_M_with_S_and_R_pair']} of A less M")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
