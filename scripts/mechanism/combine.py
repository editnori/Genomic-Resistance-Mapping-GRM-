"""combine.py: the readings of the mechanism split, over species.

    python3 -I combine.py --run DIR [--run DIR ...] --contrasts contrasts.tsv --out DIR

Each --run is one species' OUT of run_species.sh and must hold DONE; its
results/ hold strata.tsv (M1) and m2-cells.tsv and timesplit.tsv (M2). This
fits nothing and reads no genome. It writes:

- m1-reading.tsv: per split, set and row, the paired difference recall(gene
  only) less recall(point only) over the datasets where both strata hold 10
  or more resistant genomes: datasets, species, median, the species interval
  (4,000 resamples of species, seeded by the row's name), gene ahead, point
  ahead, same (margin 0.02), the two-sided sign test, and for the five
  learners on the lineage split of the primary set, Holm over the five and
  the reading. Every other row is beside and decides nothing.
- m1-neither.tsv: per split, set and row, the median recall on the neither
  stratum over datasets with 10 or more resistant genomes in it.
- m2-reading.tsv and m2-summary.tsv: Holm over the cells that pass, each
  statistic apart, per cell kind, period pair and guard; a cell moves at a
  Holm p under 0.05. Only the BioProject vote on the primary periods over
  species and class cells decides; the rest is beside.
- timesplit-link.tsv: Spearman's rank correlation, over the time-split
  datasets, of the share of resistant later test genomes with an unseen
  determinant against each learner's past-less-random total-contrast
  difference in recall and in log loss (exploratory, coarse).

Python 3.12, standard library only; run with -I.
"""

import argparse
import collections
import os
import statistics
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (MARGIN, holm, r12, read_tsv, sign_p_two_sided, spearman,  # noqa: E402
                    species_interval, write_tsv)

LEARNERS = ("scm", "cart", "rf", "gbt", "lasso")
ROWS = LEARNERS + ("catalog", "program")
M1_HEADER = ["split", "set", "row", "datasets", "species", "median", "lo", "hi", "gene_ahead", "point_ahead",
             "same", "sign_p", "holm_p", "reading", "decides", "label"]
NEITHER_HEADER = ["split", "set", "row", "datasets", "median_recall", "label"]
M2_HEADER = ["species_slug", "species", "cell", "class", "dataset", "antibiotic", "periods", "guard", "statistic",
             "value", "p", "holm_p", "moves", "decides"]
SUMMARY_HEADER = ["cell", "periods", "guard", "statistic", "cells", "pass", "move", "decides"]
LINK_HEADER = ["learner", "reading", "measure", "datasets", "spearman", "label"]
RANDOM_LABEL = "training and test share lineages"


def fmt(x):
    return "" if x is None else f"{r12(x):.12f}"


def read_runs(runs):
    strata, cells, ts = [], [], []
    for run in runs:
        if not os.path.isfile(os.path.join(run, "DONE")):
            print(f"combine.py: {run} has no DONE; a run that did not finish is not read", file=sys.stderr)
            sys.exit(1)
        res = os.path.join(run, "results")
        strata += read_tsv(os.path.join(res, "m1", "strata.tsv"))
        cells += read_tsv(os.path.join(res, "m2", "m2-cells.tsv"))
        ts += read_tsv(os.path.join(res, "m2", "timesplit.tsv"))
    return strata, cells, ts


def m1(strata):
    recall = {}
    for r in strata:
        if r["kind"] == "stratum" and r["recall"]:
            recall[(r["split"], r["set"], r["row"], r["dataset"], r["stratum"])] = (float(r["recall"]),
                                                                                    r["species_slug"])
    datasets = sorted({(r["split"], r["set"], r["row"], r["dataset"]) for r in strata})
    diffs = collections.defaultdict(lambda: collections.defaultdict(list))
    for split, dset, row, d in datasets:
        g = recall.get((split, dset, row, d, "gene_only"))
        p = recall.get((split, dset, row, d, "point_only"))
        if g and p:
            diffs[(split, dset, row)][g[1]].append(r12(g[0] - p[0]))
    out = []
    for split in ("lineage", "random"):
        for dset in ("primary", "M"):
            for row in ROWS:
                by = diffs.get((split, dset, row), {})
                vals = [v for vs in by.values() for v in vs]
                name = f"m1|{split}|{dset}|{row}|gene_only-point_only"
                lo, hi = species_interval(by, name) if vals else (None, None)
                ga = sum(v > MARGIN for v in vals)
                pa = sum(v < -MARGIN for v in vals)
                out.append({"split": split, "set": dset, "row": row, "datasets": len(vals), "species": len(by),
                            "median": fmt(statistics.median(vals)) if vals else "", "lo": fmt(lo), "hi": fmt(hi),
                            "gene_ahead": ga, "point_ahead": pa, "same": len(vals) - ga - pa,
                            "sign_p": fmt(sign_p_two_sided(ga, pa)) if vals else "",
                            "decides": int(split == "lineage" and dset == "primary" and row in LEARNERS),
                            "label": RANDOM_LABEL if split == "random" else ""})
    main = [r for r in out if r["decides"]]
    ps = [float(r["sign_p"]) if r["sign_p"] else 1.0 for r in main]
    for r, h in zip(main, holm(ps)):
        r["holm_p"] = fmt(h)
        if h < 0.05 and r["gene_ahead"] > r["point_ahead"]:
            r["reading"] = "gene only higher"
        elif h < 0.05 and r["point_ahead"] > r["gene_ahead"]:
            r["reading"] = "point only higher"
        else:
            r["reading"] = "no difference shown"
    neither = []
    for split in ("lineage", "random"):
        for dset in ("primary", "M"):
            for row in ROWS:
                vals = [v[0] for k, v in recall.items() if k[:3] == (split, dset, row) and k[4] == "neither"]
                neither.append({"split": split, "set": dset, "row": row, "datasets": len(vals),
                                "median_recall": fmt(statistics.median(vals)) if vals else "",
                                "label": RANDOM_LABEL if split == "random" else ""})
    return out, neither


def m2(cells):
    out, summary = [], []
    groups = collections.defaultdict(list)
    for r in cells:
        for stat, value, p in (("point", r.get("d_point", ""), r.get("p_point", "")),
                               ("tvd", r.get("tvd", ""), r.get("p_tvd", ""))):
            key = (r["cell"], r["periods"], r["guard"], stat)
            row = {k: r.get(k, "") for k in M2_HEADER[:8]}
            row.update(statistic=stat, value=value, p=p,
                       decides=int(key[:3] == ("species_class", "primary", "bioproject")))
            row["_pass"] = r.get("passes") == "1"
            groups[key].append(row)
    for key in sorted(groups):
        rows = groups[key]
        passing = [r for r in rows if r["_pass"] and r["p"]]
        for r, h in zip(passing, holm([float(r["p"]) for r in passing])):
            r["holm_p"] = fmt(h)
            r["moves"] = int(h < 0.05)
        for r in rows:
            r.pop("_pass")
            out.append(r)
        summary.append({"cell": key[0], "periods": key[1], "guard": key[2], "statistic": key[3],
                        "cells": len(rows), "pass": len(passing),
                        "move": sum(1 for r in passing if r.get("moves") == 1),
                        "decides": int(key[:3] == ("species_class", "primary", "bioproject"))})
    return out, summary


def link(ts, contrasts):
    share = {r["dataset"]: float(r["share_unseen"]) for r in ts if r["reading"] == "point_as_one" and r["share_unseen"]}
    out = []
    by = collections.defaultdict(dict)
    for r in contrasts:
        if r["contrast"] == "total":
            by[r["learner"]][r["dataset"]] = r
    for learner in sorted(by):
        for measure in ("d_recall", "d_log_loss"):
            ds = sorted(d for d in by[learner] if d in share and by[learner][d][measure] not in ("", "NA"))
            rho = spearman([share[d] for d in ds], [float(by[learner][d][measure]) for d in ds])
            out.append({"learner": learner, "reading": "point_as_one", "measure": measure, "datasets": len(ds),
                        "spearman": fmt(rho), "label": "exploratory; a retrospective split inside one database"})
    return out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", action="append", required=True)
    ap.add_argument("--contrasts", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args(argv)
    strata, cells, ts = read_runs(args.run)
    m1_rows, neither = m1(strata)
    m2_rows, summary = m2(cells)
    link_rows = link(ts, read_tsv(args.contrasts))
    os.makedirs(args.out, exist_ok=True)
    write_tsv(os.path.join(args.out, "m1-reading.tsv"), M1_HEADER, m1_rows)
    write_tsv(os.path.join(args.out, "m1-neither.tsv"), NEITHER_HEADER, neither)
    write_tsv(os.path.join(args.out, "m2-reading.tsv"), M2_HEADER, m2_rows)
    write_tsv(os.path.join(args.out, "m2-summary.tsv"), SUMMARY_HEADER, summary)
    write_tsv(os.path.join(args.out, "timesplit-link.tsv"), LINK_HEADER, link_rows)
    print(f"combine: {len(args.run)} species runs; m1 rows {len(m1_rows)}, m2 rows {len(m2_rows)}, "
          f"time-split link rows {len(link_rows)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
