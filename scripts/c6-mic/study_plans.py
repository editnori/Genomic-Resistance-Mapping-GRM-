"""The sixth study plan of the MIC fold-plan experiment (docs/mic-plan.md), `c6mic-s47`.

  python3 -I study_plans.py <release dir> <plan dir> <dataset list> <out dir>

<plan dir> holds the five reviewed study plans of each dataset as fold files
(docs/fold-file.md): <dataset>.study-whole.json and <dataset>.study-d1.json
to <dataset>.study-d4.json. <dataset list> names one dataset per line.

From the release's report <release dir>/datasets/<dataset>.json this reads
the samples' ids and reported calls, and the k-mer copy groups
(`kset_duplicates`), and nothing else: no score, metric or model. The study
of each genome is not read. The plan is placed on units instead: a unit is
the genomes that every one of the five plans tests in one fold. Each reviewed
plan keeps every study whole, so a unit is a union of whole studies and the
new plan keeps every study whole too.

Placement is rule B of the reviewed plans, the evaluator's lineage step on
units in a keyed order: units sorted by
sha256("c6mic-20261122|study-s47|" + the unit's lowest genome id), each into
the fold that least raises the squared class imbalance summed over folds,
ties to the lowest fold. Each plan is checked: every genome in one test fold,
both classes in every fold's test and training genomes, no copy group split,
no unit split. A dataset whose plan fails a check gets no file.

Writes <out dir>/<dataset>.c6mic-s47.json (grm-fold-file-v1) and
<out dir>/study-plans.tsv: one row per dataset with its genome and unit
counts, the check and the file's sha256. Exits 1 when any plan fails a check.
The fold files list genome ids, so they stay on the machine that fits them.
"""

import hashlib
import json
import sys
from pathlib import Path

K = 5
KEY = "c6mic-20261122|study-s47"
SCHEME = "c6mic-s47"
REVIEWED = ["study-whole", "study-d1", "study-d2", "study-d3", "study-d4"]


def gnum(g):
    """A BV-BRC genome id as numbers, so 9.10 sorts after 9.9."""
    a, b = g.split(".", 1)
    return (int(a), int(b))


def order_key(lowest):
    return hashlib.sha256(f"{KEY}|{lowest}".encode()).hexdigest()


def read_plan(doc, ids, name="plan"):
    """The test fold of each genome in `ids` order, from a fold file's
    document, or a ValueError unless it tests every genome once in K folds."""
    idx = {g: i for i, g in enumerate(ids)}
    fold = [None] * len(ids)
    if doc.get("format") != "grm-fold-file-v1" or len(doc.get("test", [])) != K:
        raise ValueError(f"{name}: not a fold file of {K} test folds")
    if doc.get("train") or doc.get("unused"):
        raise ValueError(f"{name}: lists genomes to train or leave unused")
    for f, members in enumerate(doc["test"]):
        for g in members:
            if g not in idx:
                raise ValueError(f"{name}: a genome outside the cohort")
            if fold[idx[g]] is not None:
                raise ValueError(f"{name}: a genome listed twice")
            fold[idx[g]] = f
    if None in fold:
        raise ValueError(f"{name}: a cohort genome in no test fold")
    return fold


def units_of(plans):
    """Genomes every plan tests in the same fold, as index lists in sample
    order, ordered by first member."""
    cells = {}
    for i in range(len(plans[0])):
        cells.setdefault(tuple(p[i] for p in plans), []).append(i)
    return sorted(cells.values(), key=lambda m: m[0])


def place(y, order, k=K):
    """The evaluator's placement step (`grouped_folds` after its order)."""
    total_r = sum(y)
    total_s = len(y) - total_r
    load = [[0, 0] for _ in range(k)]
    fold = [None] * len(y)
    for m in order:
        r = sum(1 for i in m if y[i])
        s = len(m) - r

        def cost(f):
            lr, ls = load[f]
            return ((lr + r) ** 2 - lr**2) / max(total_r, 1) ** 2 + ((ls + s) ** 2 - ls**2) / max(total_s, 1) ** 2

        best = min(range(k), key=lambda f: (cost(f), f))
        load[best][0] += r
        load[best][1] += s
        for i in m:
            fold[i] = best
    return fold


def plan(ids, y, units):
    order = sorted(units, key=lambda m: order_key(min((ids[i] for i in m), key=gnum)))
    return place(y, order)


def check(y, fold, dups, units, k=K):
    """Problems with a plan; empty when it passes."""
    n = len(y)
    if len(fold) != n or any(f is None or not 0 <= f < k for f in fold):
        return ["a genome without a test fold"]
    out = []
    for f in range(k):
        for side, part in (
            ("test", [y[i] for i in range(n) if fold[i] == f]),
            ("training", [y[i] for i in range(n) if fold[i] != f]),
        ):
            if not part or all(part) or not any(part):
                out.append(f"fold {f}'s {side} genomes lack a class")
    if any(len({fold[i] for i in g}) > 1 for g in dups):
        out.append("a copy group is split")
    if any(len({fold[i] for i in u}) > 1 for u in units):
        out.append("a unit is split")
    return out


def fold_file(dataset, ids, fold, k=K):
    test = [[ids[i] for i in range(len(ids)) if fold[i] == f] for f in range(k)]
    doc = {
        "format": "grm-fold-file-v1",
        "dataset": dataset,
        "scheme": SCHEME,
        "test": test,
        "train": [],
        "unused": [],
    }
    return (json.dumps(doc, indent=1) + "\n").encode()


def make(dataset, report, reviewed):
    """The plan for one dataset: (fold file bytes or None, row). `report` is
    the release report's document and `reviewed` the five plans' documents
    in REVIEWED order."""
    ids = [s["id"] for s in report["samples"]]
    y = [bool(s["resistant"]) for s in report["samples"]]
    dups = [sorted(g) for g in (report.get("kset_duplicates") or [])]
    plans = [read_plan(doc, ids, f"{dataset}.{name}") for doc, name in zip(reviewed, REVIEWED)]
    units = units_of(plans)
    fold = plan(ids, y, units)
    problems = check(y, fold, dups, units)
    data = None if problems else fold_file(dataset, ids, fold)
    row = {
        "dataset": dataset,
        "genomes": len(ids),
        "units": len(units),
        "check": "; ".join(problems) or "ok",
        "sha256": hashlib.sha256(data).hexdigest() if data else "",
    }
    return data, row


def main(argv):
    if len(argv) != 4:
        print(__doc__.strip().splitlines()[2], file=sys.stderr)
        return 2
    release, plans_dir, listing, out = map(Path, argv)
    out.mkdir(parents=True, exist_ok=True)
    datasets = [line.strip() for line in listing.read_text().splitlines() if line.strip()]
    rows = []
    for did in datasets:
        report = json.loads((release / "datasets" / f"{did}.json").read_text())
        reviewed = [json.loads((plans_dir / f"{did}.{p}.json").read_text()) for p in REVIEWED]
        data, row = make(did, report, reviewed)
        if data is not None:
            (out / f"{did}.{SCHEME}.json").write_bytes(data)
        rows.append(row)
    cols = ["dataset", "genomes", "units", "check", "sha256"]
    with open(out / "study-plans.tsv", "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")
    failed = sum(1 for r in rows if r["check"] != "ok")
    print(f"study plans: {len(rows) - failed} written, {failed} failed a check")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
