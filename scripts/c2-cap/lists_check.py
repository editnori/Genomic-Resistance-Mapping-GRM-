"""Holds a run's dataset list to the first run (docs/cap-arms.md, Datasets).

  python3 -I lists_check.py <lists.tsv> <marks.tsv> <species slug> <datasets>

<datasets> is a comma list. Each must be a dataset of the species, listed
once, in CAP and in the first run by lists.tsv, and not M by the label
count's marks.tsv (pairs_S+I or pairs_I+R above 0). Exits 1 naming the
count of each kind of stray, else prints the count.
"""

import csv
import sys


def read(path):
    with open(path, encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def strays(lists_rows, marks_rows, slug, listed):
    first = {r["dataset"] for r in lists_rows if r["species_slug"] == slug and r["first_run"] == "1"}
    marked = {r["dataset"]: int(r["pairs_S+I"]) > 0 or int(r["pairs_I+R"]) > 0 for r in marks_rows}
    out = {"twice": len(listed) - len(set(listed)),
           "other species": sum(not d.startswith(slug + "__") for d in listed),
           "not first run": sum(d not in first for d in listed),
           "M": sum(marked.get(d, True) for d in listed)}
    return {k: v for k, v in out.items() if v}


def main(argv):
    if len(argv) != 5:
        raise SystemExit(__doc__)
    listed = [d for d in argv[4].split(",") if d]
    bad = strays(read(argv[1]), read(argv[2]), argv[3], listed)
    if bad or not listed:
        print("refused: " + (", ".join(f"{k} {v}" for k, v in bad.items()) or "no dataset"))
        return 1
    print(f"datasets: {len(listed)}, each a first-run dataset of {argv[3]}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
