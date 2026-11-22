"""The plan's sets A and MISS from the repository's comparison record.

    python3 -I miss_sets.py [--species SLUG] [--list a|miss]

A: datasets with an organism option and a lineage split
(reference/amrfinderplus-comparison/datasets.tsv). MISS: datasets of A where
the program's balanced accuracy (`program:tool_default`) exceeds a learner's
on the lineage split (`learner:<name>:lineage`) by 0.10 or more, each
learner on its own set. The record keeps six decimals, so a gap within 1e-6
of 0.10 is reported apart; D0 recomputes every gap from the calls and the
scores. These sets are chosen by their result: they order comparisons and
test nothing.
"""

import argparse
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REF = os.path.join(HERE, "..", "..", "reference", "amrfinderplus-comparison")
LEARNERS = ("scm", "cart", "rf", "gbt", "lasso")


def sets(ref=REF):
    with open(os.path.join(ref, "datasets.tsv"), encoding="utf-8", newline="") as f:
        a = [r["dataset"] for r in csv.DictReader(f, delimiter="\t")
             if r["organism_option"] == "yes" and r["lineage_split"] == "yes"]
    with open(os.path.join(ref, "callers.tsv"), encoding="utf-8", newline="") as f:
        ba = {(r["dataset"], r["caller"]): float(r["balanced_accuracy"]) for r in csv.DictReader(f, delimiter="\t")}
    miss = {l: [] for l in LEARNERS}
    close = []
    for d in a:
        p = ba[(d, "program:tool_default")]
        for l in LEARNERS:
            gap = p - ba[(d, f"learner:{l}:lineage")]
            if abs(gap - 0.10) < 1e-6:
                close.append((d, l))
            if gap >= 0.10 - 1e-9:
                miss[l].append(d)
    return a, miss, close


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--species")
    ap.add_argument("--list", choices=("a", "miss"))
    ap.add_argument("--reference", default=REF)
    args = ap.parse_args(argv)
    a, miss, close = sets(args.reference)
    keep = (lambda d: d.split("__")[0] == args.species) if args.species else (lambda d: True)
    union = sorted({d for v in miss.values() for d in v})
    if args.list:
        # One id per line, and nothing at all for a species with none.
        for d in a if args.list == "a" else union:
            if keep(d):
                print(d)
        return
    species = sorted({d.split("__")[0] for d in union})
    print(f"A: {len(a)} datasets; MISS: {len(union)} datasets in {len(species)} species "
          f"(sets chosen by their result; they order comparisons and test nothing)")
    print("per learner: " + ", ".join(f"{l} {sum(keep(d) for d in miss[l])}" for l in LEARNERS))
    if close:
        print(f"within 1e-6 of the cut: {len(close)}", file=sys.stderr)


if __name__ == "__main__":
    main()
