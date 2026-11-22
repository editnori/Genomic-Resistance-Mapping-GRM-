"""The dataset lists of plan block C9 (docs/served-pick.md).

  python3 -I lists.py <release dir> <label count table> <study plan list> <out dir>

L is the plan's set: the release's reports whose lineage split ran (its
`unavailable_reason` is empty). Each is read from
<release dir>/datasets/<dataset>.json: the species, the antibiotic, the
counts and the lineage split's reason, nothing else.

M is a dataset whose label count row (<label count table>, the label
count's `datasets.tsv`: columns `pairs_S+I` and `pairs_I+R`) holds any
genome pair of S and I calls or of I and R calls. Its labels depend on the
maintainer's label rule, so it waits. The column `pairs_S+R` counts genome
pairs with both an S and an R call; the most resistant call labels them R.

S is a dataset with a reviewed study plan: a line of <study plan list> (the
sha256 list of the reviewed fold files) that names
folds/<dataset>.study-whole.json.

Writes <out dir>/l.tsv, one row per L dataset with its species and
antibiotic in full and the label rules its measures depend on, and
<out dir>/counts.tsv, every count with its denominator, per species too.
The first run takes L less M. No genome id is read or written.
"""

import csv
import json
import re
import sys
from pathlib import Path

PLAN = re.compile(r"(?:^|[ /])folds/([a-z0-9_]+__[a-z0-9_]+)\.study-whole\.json$")
COLUMNS = [
    "dataset", "species_slug", "species", "antibiotic", "samples", "resistant", "susceptible",
    "lineages", "in_M", "in_S", "first_run", "pairs_S+R", "labels_depend_on",
]
# Which of the four label rules (r-or-s-wins, the release's; intermediate-aside;
# most-resistant; intermediate-class) a dataset's labels, and so every measure
# on it, depend on. Not in M, three of the four agree; the most resistant call
# labels an S and R pair R, so a dataset with one depends on that rule against
# the other three. In M, the labels wait for the maintainer's rule.
DEPENDS = {"none": "none of the four", "most": "most-resistant against the other three", "M": "the maintainer's rule"}


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
        "susceptible": c["susceptible"], "lineages": lineage[0]["n_groups"],
        "lineage_ran": lineage[0].get("unavailable_reason") is None,
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
    """The L rows and the count table from the parsed inputs."""
    ids = {r["dataset"] for r in reports}
    if len(ids) != len(reports):
        raise ValueError("a dataset has two reports")
    if set(labels) != ids:
        raise ValueError("the label count table and the release name other datasets")
    if not study <= ids:
        raise ValueError("a study plan names a dataset outside the release")
    rows = []
    for r in sorted(reports, key=lambda r: r["dataset"]):
        if r["lineage_ran"]:
            m, sr = labels[r["dataset"]]
            row = {k: r[k] for k in COLUMNS if k in r}
            row.update({"in_M": int(m), "in_S": int(r["dataset"] in study), "first_run": int(not m), "pairs_S+R": sr,
                        "labels_depend_on": DEPENDS["M" if m else "most" if sr > 0 else "none"]})
            rows.append(row)
    n = len(reports)
    m_all = sum(1 for d in ids if labels[d][0])
    first = [r for r in rows if r["first_run"]]
    in_m = [r for r in rows if r["in_M"]]
    species = len({r["species_slug"] for r in reports})
    counts = [
        ("reports in the release", n, n),
        ("L: reports whose lineage split ran", len(rows), n),
        ("L species", len({r["species_slug"] for r in rows}), species),
        ("M datasets in the release", m_all, n),
        ("S datasets in the release", len(study), n),
        ("S datasets in L", sum(r["in_S"] for r in rows), len(study)),
        ("L datasets in M", len(in_m), len(rows)),
        ("L datasets not in M (the first run)", len(first), len(rows)),
        ("first-run species", len({r["species_slug"] for r in first}), len({r["species_slug"] for r in rows})),
        ("first-run datasets in S (the first study run)", sum(r["in_S"] for r in first), len(first)),
        ("L datasets in M and S", sum(r["in_S"] for r in in_m), len(in_m)),
        ("first-run datasets with an S and R pair", sum(r["pairs_S+R"] > 0 for r in first), len(first)),
        ("first-run datasets whose labels depend on none of the four rules",
         sum(r["labels_depend_on"] == DEPENDS["none"] for r in first), len(first)),
    ]
    for sp in sorted({r["species_slug"] for r in rows}):
        sr = [r for r in rows if r["species_slug"] == sp]
        fr = [r for r in sr if r["first_run"]]
        counts.append((f"{sp}: L, not M, M, not M in S, M in S",
                       f"{len(sr)}, {len(fr)}, {len(sr) - len(fr)}, {sum(r['in_S'] for r in fr)}, "
                       f"{sum(r['in_S'] for r in sr if r['in_M'])}", len(sr)))
    return rows, counts


def main(argv):
    if len(argv) != 4:
        print(__doc__.strip().splitlines()[2], file=sys.stderr)
        return 2
    release, labels, study, out = Path(argv[0]), Path(argv[1]), Path(argv[2]), Path(argv[3])
    reports = [read_report(p) for p in sorted((release / "datasets").glob("*.json"))]
    rows, counts = lists(reports, read_labels(labels), read_study(study))
    out.mkdir(parents=True, exist_ok=True)
    with open(out / "l.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, COLUMNS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    with open(out / "counts.tsv", "w") as f:
        f.write("count\tvalue\tof\n")
        for name, value, of in counts:
            f.write(f"{name}\t{value}\t{of}\n")
    for name, value, of in counts:
        print(f"{name}: {value} of {of}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
