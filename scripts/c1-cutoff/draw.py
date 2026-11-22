"""The cohort draws of plan block C1 part B (docs/cohort-draw.md, the draw).

  python3 -I draw.py <release dir> <datasets> <lines> <draws> <out dir> [<study plan dir>]

<datasets> is a comma list of dataset ids, <lines> a comma list of r (genomes
per class) and <draws> the number of draws per dataset and line, numbered
from 0. From the release's report <release dir>/datasets/<dataset>.json this
reads each sample's genome id, reported call (`resistant`) and lineage
number, and nothing else: no score, metric or model.

For draw j of dataset d a lineage l gets the key
sha256("c1b-20261122|" + d + "|draw" + j + "|lineage|" + l) and a genome g
the key sha256("c1b-20261122|" + d + "|draw" + j + "|genome|" + g), compared
as hex text. A class's order: the lineages holding a genome of the class by
lineage key, and inside each its genomes of the class by genome key; pass 1
takes the first genome of each lineage, pass 2 the second of each lineage
that has one, and so on. The cohort at r is the first r genomes of each
class's order. The keys do not hold r, so a draw's cohort at a smaller r
lies inside its cohort at a larger one.

A class with r or fewer genomes refuses the draw at r: no cohort file, and
the reason in draws.tsv. No class is filled from another.

With <study plan dir> (the reviewed fold files, <dataset>.study-whole.json),
each cohort genome of a dataset that has one gets its study-whole fold.

Writes <out dir>/<dataset>.r<r>.d<j>.cohort.json (grm-cohort-draw-v1; it
lists genome ids, so it stays on the machine that fits it) and
<out dir>/draws.tsv: one row per dataset, line and draw with the cohort
file's sha256 and lineage counts, or the refusal. draws.tsv holds no genome id.
"""

import hashlib
import json
import sys
from pathlib import Path

KEY = "c1b-20261122"
FORMAT = "grm-cohort-draw-v1"
STUDY_FOLDS = 5
COLUMNS = [
    "dataset", "r", "draw", "written", "reason", "resistant_in_dataset", "susceptible_in_dataset",
    "cohort_lineages", "resistant_lineages", "susceptible_lineages", "study", "file", "sha256",
]


def key(dataset, draw, kind, name):
    return hashlib.sha256(f"{KEY}|{dataset}|draw{draw}|{kind}|{name}".encode()).hexdigest()


def class_order(samples, dataset, draw, resistant):
    """The draw's order of one class's genomes, as sample indices."""
    by_lineage = {}
    for i, s in enumerate(samples):
        if s["resistant"] == resistant:
            by_lineage.setdefault(s["lineage"], []).append(i)
    lineages = sorted(by_lineage, key=lambda l: key(dataset, draw, "lineage", l))
    lists = [sorted(by_lineage[l], key=lambda i: key(dataset, draw, "genome", samples[i]["id"])) for l in lineages]
    order = []
    depth = 0
    while True:
        taken = [m[depth] for m in lists if depth < len(m)]
        if not taken:
            return order
        order.extend(taken)
        depth += 1


def cohort(samples, dataset, draw, r):
    """(sample indices of the cohort, None) or (None, the refusal)."""
    out = []
    for resistant in (True, False):
        n = sum(1 for s in samples if s["resistant"] == resistant)
        if n <= r:
            name = "resistant" if resistant else "susceptible"
            return None, f"class short of r: {n} {name} for r {r}"
        out.extend(class_order(samples, dataset, draw, resistant)[:r])
    return out, None


def read_samples(release, dataset):
    with open(Path(release) / "datasets" / f"{dataset}.json") as f:
        rep = json.load(f)
    if rep["id"] != dataset:
        raise ValueError(f"{dataset}: the report names another dataset")
    samples = [{"id": s["id"], "resistant": bool(s["resistant"]), "lineage": int(s["lineage"])} for s in rep["samples"]]
    if len({s["id"] for s in samples}) != len(samples):
        raise ValueError(f"{dataset}: a genome listed twice")
    return samples


def read_study(path, ids):
    """genome id -> study-whole fold, checked to test every cohort genome once."""
    with open(path) as f:
        doc = json.load(f)
    if doc.get("format") != "grm-fold-file-v1" or len(doc.get("test", [])) != STUDY_FOLDS:
        raise ValueError(f"{path.name}: not a fold file of {STUDY_FOLDS} test folds")
    fold = {}
    for k, members in enumerate(doc["test"]):
        for g in members:
            if g in fold:
                raise ValueError(f"{path.name}: a genome listed twice")
            fold[g] = k
    if set(fold) != set(ids):
        raise ValueError(f"{path.name}: does not test exactly the dataset's genomes")
    return fold


def draw_dataset(samples, dataset, lines, draws, study=None):
    """[(row, document or None)] for every line and draw of one dataset."""
    n_r = sum(1 for s in samples if s["resistant"])
    out = []
    for r in lines:
        for j in range(draws):
            row = {
                "dataset": dataset, "r": r, "draw": j, "written": 0, "reason": "",
                "resistant_in_dataset": n_r, "susceptible_in_dataset": len(samples) - n_r,
                "cohort_lineages": "", "resistant_lineages": "", "susceptible_lineages": "",
                "study": int(study is not None), "file": "", "sha256": "",
            }
            rows, why = cohort(samples, dataset, j, r)
            if rows is None:
                row["reason"] = why
                out.append((row, None))
                continue
            doc = {
                "format": FORMAT, "dataset": dataset, "r": r, "draw": j,
                "genomes": [samples[i]["id"] for i in rows],
                "resistant": [samples[i]["resistant"] for i in rows],
                "lineage": [samples[i]["lineage"] for i in rows],
            }
            if study is not None:
                doc["study_fold"] = [study[samples[i]["id"]] for i in rows]
            row["written"] = 1
            row["cohort_lineages"] = len({samples[i]["lineage"] for i in rows})
            row["resistant_lineages"] = len({samples[i]["lineage"] for i in rows if samples[i]["resistant"]})
            row["susceptible_lineages"] = len({samples[i]["lineage"] for i in rows if not samples[i]["resistant"]})
            row["file"] = f"{dataset}.r{r}.d{j}.cohort.json"
            out.append((row, doc))
    return out


def encode(doc):
    return (json.dumps(doc, indent=1, sort_keys=True) + "\n").encode()


def main(argv):
    if len(argv) not in (5, 6):
        print(__doc__.strip().splitlines()[2], file=sys.stderr)
        return 2
    release, names, lines, draws, out = argv[0], argv[1].split(","), argv[2].split(","), int(argv[3]), Path(argv[4])
    lines = [int(r) for r in lines]
    if len(set(names)) != len(names) or len(set(lines)) != len(lines) or min(lines) < 1 or draws < 1:
        print("a dataset or line twice, a line below 1 or no draw", file=sys.stderr)
        return 2
    plans = Path(argv[5]) if len(argv) == 6 else None
    out.mkdir(parents=True, exist_ok=True)
    table = []
    for d in names:
        samples = read_samples(release, d)
        study = None
        if plans is not None and (plans / f"{d}.study-whole.json").is_file():
            study = read_study(plans / f"{d}.study-whole.json", [s["id"] for s in samples])
        for row, doc in draw_dataset(samples, d, lines, draws, study):
            if doc is not None:
                data = encode(doc)
                part = out / (row["file"] + ".part")
                part.write_bytes(data)
                part.replace(out / row["file"])
                row["sha256"] = hashlib.sha256(data).hexdigest()
            table.append(row)
    with open(out / "draws.tsv.part", "w") as f:
        f.write("\t".join(COLUMNS) + "\n")
        for row in table:
            f.write("\t".join(str(row[c]) for c in COLUMNS) + "\n")
    (out / "draws.tsv.part").replace(out / "draws.tsv")
    written = sum(r["written"] for r in table)
    print(f"draws: {written} cohort files of {len(table)} draws in {len(names)} datasets, {len(table) - written} refused")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
