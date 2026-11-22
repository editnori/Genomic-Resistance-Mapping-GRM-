"""The dataset lists of plan block C1 part B (docs/cohort-draw.md).

  python3 -I lists.py <release dir> <label count table> <study plan list> <out dir>

BIG is the plan's set: the release's reports whose lineage split ran (its
`unavailable_reason` is empty) and whose smaller class holds 200 or more
genomes. Each is read from <release dir>/datasets/<dataset>.json: the counts
and the lineage split's reason, nothing else.

M is a dataset whose label count row (<label count table>, the `datasets.tsv`
of the label count: columns `pairs_S+I` and `pairs_I+R`) holds any genome
pair of S and I calls or of I and R calls. Its labels depend on the
maintainer's label rule, so it waits. The column `pairs_S+R` counts the
genomes with both an S and an R call; the most resistant call relabels them.

S is a dataset with a reviewed study plan: a line of <study plan list> (the
sha256 list of the reviewed fold files) that names
folds/<dataset>.study-whole.json.

Writes <out dir>/big.tsv, one row per BIG dataset, and
<out dir>/counts.tsv, every count with its denominator. The first run takes
the BIG datasets not in M. No genome id is read or written.
"""

import csv
import json
import re
import sys
from pathlib import Path

SMALLER_CLASS = 200
PLAN = re.compile(r"(?:^|[ /])folds/([a-z0-9_]+__[a-z0-9_]+)\.study-whole\.json$")
COLUMNS = [
    "dataset", "species_slug", "species", "antibiotic", "samples", "resistant", "susceptible",
    "smaller_class", "lineages", "in_M", "in_S", "first_run", "pairs_S+R",
]


def read_report(path):
    """The fields of one release report the lists read."""
    with open(path) as f:
        rep = json.load(f)
    lineage = [s for s in rep["splits"] if s["scheme"] == "lineage"]
    if len(lineage) != 1:
        raise ValueError(f"{path.name}: not one lineage split")
    c = rep["counts"]
    if c["resistant"] + c["susceptible"] != c["samples"]:
        raise ValueError(f"{path.name}: the class counts do not sum to the samples")
    return {
        "dataset": rep["id"], "species_slug": rep["species_slug"], "species": rep["species"],
        "antibiotic": rep["antibiotic"], "samples": c["samples"], "resistant": c["resistant"],
        "susceptible": c["susceptible"], "smaller_class": min(c["resistant"], c["susceptible"]),
        "lineages": lineage[0]["n_groups"], "lineage_ran": lineage[0].get("unavailable_reason") is None,
    }


def read_labels(path):
    """dataset -> (in M, pairs_S+R) from the label count table."""
    out = {}
    with open(path, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            m = int(row["pairs_S+I"]) > 0 or int(row["pairs_I+R"]) > 0
            out[row["dataset"]] = (m, int(row["pairs_S+R"]))
    return out


def read_study(path):
    out = set()
    with open(path) as f:
        for line in f:
            hit = PLAN.search(line.rstrip("\n"))
            if hit:
                out.add(hit.group(1))
    return out


def lists(reports, labels, study):
    """The BIG rows and the count table from the parsed inputs."""
    ids = {r["dataset"] for r in reports}
    if len(ids) != len(reports):
        raise ValueError("a dataset has two reports")
    if set(labels) != ids:
        raise ValueError("the label count table and the release name other datasets")
    if not study <= ids:
        raise ValueError("a study plan names a dataset outside the release")
    big = []
    for r in sorted(reports, key=lambda r: r["dataset"]):
        if r["lineage_ran"] and r["smaller_class"] >= SMALLER_CLASS:
            m, sr = labels[r["dataset"]]
            row = {k: r[k] for k in COLUMNS if k in r}
            row.update({"in_M": int(m), "in_S": int(r["dataset"] in study), "first_run": int(not m), "pairs_S+R": sr})
            big.append(row)
    n = len(reports)
    m_all = sum(1 for d in ids if labels[d][0])
    first = [r for r in big if r["first_run"]]
    counts = [
        ("reports in the release", n, n),
        ("reports whose lineage split ran", sum(r["lineage_ran"] for r in reports), n),
        (f"reports whose smaller class holds {SMALLER_CLASS} or more", sum(r["smaller_class"] >= SMALLER_CLASS for r in reports), n),
        ("BIG datasets", len(big), n),
        ("BIG species", len({r["species_slug"] for r in big}), len({r["species_slug"] for r in reports})),
        ("M datasets", m_all, n),
        ("S datasets", len(study), n),
        ("BIG datasets in M", sum(r["in_M"] for r in big), len(big)),
        ("BIG datasets not in M (the first run)", len(first), len(big)),
        ("BIG datasets in S", sum(r["in_S"] for r in big), len(big)),
        ("first-run datasets in S", sum(r["in_S"] for r in first), len(first)),
        ("BIG datasets in M and S", sum(r["in_M"] and r["in_S"] for r in big), len(big)),
        ("first-run species", len({r["species_slug"] for r in first}), len({r["species_slug"] for r in big})),
        ("first-run datasets with an S and R pair", sum(r["pairs_S+R"] > 0 for r in first), len(first)),
        ("not-M datasets with an S and R pair", sum(1 for d in ids if not labels[d][0] and labels[d][1] > 0), n - m_all),
    ]
    for sp in sorted({r["species_slug"] for r in big}):
        rows = [r for r in big if r["species_slug"] == sp]
        fr = [r for r in rows if r["first_run"]]
        counts.append((f"{sp}: BIG, first run, first run in S", f"{len(rows)}, {len(fr)}, {sum(r['in_S'] for r in fr)}", len(rows)))
    return big, counts


def main(argv):
    if len(argv) != 4:
        print(__doc__.strip().splitlines()[2], file=sys.stderr)
        return 2
    release, labels, study, out = Path(argv[0]), Path(argv[1]), Path(argv[2]), Path(argv[3])
    reports = [read_report(p) for p in sorted((release / "datasets").glob("*.json"))]
    big, counts = lists(reports, read_labels(labels), read_study(study))
    out.mkdir(parents=True, exist_ok=True)
    with open(out / "big.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, COLUMNS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(big)
    with open(out / "counts.tsv", "w") as f:
        f.write("count\tvalue\tof\n")
        for name, value, of in counts:
            f.write(f"{name}\t{value}\t{of}\n")
    for name, value, of in counts:
        print(f"{name}: {value} of {of}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
