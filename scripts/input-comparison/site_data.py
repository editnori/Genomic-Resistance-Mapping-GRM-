"""site_data.py READING LEARNERS [--check FILE]: the part of the input comparison's readings that the web app bundles.

READING is reference/input-comparison/second-reading.json, the JSON second_reading.py wrote for H-1's four
species. LEARNERS is reference/input-comparison/learner-reading.json, the JSON learner_reading.py wrote for the
same run. The script writes, as JSON, what the site's "Input settings against the default" section counts:

- each species with its datasets and the datasets with a lineage answer;
- per arm, in the runner document's order: the datasets better and worse over the four species with the sign
  test's p and the median of the dataset means, the same counts per species, and the family C contrasts made,
  inside the margin, beyond it and exactly zero, summed over the species;
- the panel family's size and tested contrasts, the firm decisions by direction and by learner, and the
  predictions;
- the log-loss margin, which the four species specs under scripts/input-comparison/h1 share;
- from the reading by learner, the lineage split only: per arm and learner, the datasets better, worse and
  unchanged with the sign test's p, the summed calls' recall and specificity under the default and under the
  arm, and the learner's share of the five learners' pooled change; and the rows the page names (NAMED).

It recomputes nothing the reading decided: every number is the reading's, or a sum or a count of its records.

  python3 -I -B scripts/input-comparison/site_data.py reference/input-comparison/second-reading.json \
      reference/input-comparison/learner-reading.json > web/src/data/input-comparison.json
  python3 -I -B scripts/input-comparison/site_data.py reference/input-comparison/second-reading.json \
      reference/input-comparison/learner-reading.json --check web/src/data/input-comparison.json

With --check it prints nothing and exits 1 unless FILE is exactly the text the script prints.
"""
import json
import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
# The arms in the order the runner document's tables list them; the two twins last.
ARMS = ["cap-1000", "cap-4000", "grr-v1", "k21", "k25", "k-c", "k-h", "k-p", "k-c-h-p", "budget-1000", "budget-16000",
        "sentinel-seed-twin", "sentinel-tuning-twin"]
SPECS = ["salmonella_enterica", "staphylococcus_aureus", "neisseria_gonorrhoeae", "mycobacterium_tuberculosis"]
COUNTS = ("better", "worse", "unchanged", "datasets")
# The learners in the site's order, and the rows of the reading by learner the page names: (dataset, arm, learner).
LEARNERS = ["scm", "cart", "rf", "gbt", "lasso"]
NAMED = [("mycobacterium_tuberculosis__clofazimine", "cap-4000", "scm")]


def margin():
    """The log-loss margin of the four species specs, which must agree."""
    values = {json.loads((HERE / "h1" / f"h1-{s}.json").read_text(encoding="utf-8"))["margin"]["log_loss"] for s in SPECS}
    if len(values) != 1:
        raise SystemExit(f"the species specs disagree on the log-loss margin: {sorted(values)}")
    return values.pop()


def learner(contrast_id):
    """'<dataset>/C:rf:lineage:default->arm:log_loss' -> 'rf'."""
    return contrast_id.split("/", 1)[1].split(":")[1]


def by_learner(path, species):
    """The lineage split of the reading by learner, for the arms in ARMS."""
    reading = json.loads(pathlib.Path(path).read_text(encoding="utf-8"))
    mine = [(s["species"], s["datasets"], s["lineage_answer"]) for s in reading["species"]]
    theirs = [(s["species"], s["datasets"], s["datasets"] - len(s["unanswered"])) for s in species]
    if mine != theirs:
        raise SystemExit(f"the two readings disagree on the species: {mine} {theirs}")
    lineage = reading["tables"]["lineage"]
    if sorted(lineage) != sorted(ARMS) or any(sorted(lineage[a]) != sorted(LEARNERS) for a in ARMS):
        raise SystemExit("the reading by learner's arms or learners are not the document's")
    rows = []
    for dataset, arm, learner_ in NAMED:
        found = [r for r in reading["rows"] if r["split"] == "lineage" and (r["dataset"], r["arm"], r["learner"]) == (dataset, arm, learner_)]
        if len(found) != 1:
            raise SystemExit(f"{dataset} {arm} {learner_}: {len(found)} lineage rows")
        r = found[0]
        rows.append({"dataset": dataset, "arm": arm, "learner": learner_, "estimate": r["estimate"], "repeats": r["repeats"], "default": r["calls_default"], "arm_calls": r["calls_arm"]})

    def entry(t):
        return {
            **{k: t[k] for k in COUNTS},
            "p": t["p"],
            "recall": [t["rates_default"]["recall"], t["rates_arm"]["recall"]],
            "specificity": [t["rates_default"]["specificity"], t["rates_arm"]["specificity"]],
            "share": t["share"],
        }

    return {
        "exploratory": reading["exploratory"],
        "arms": {a: {l: entry(lineage[a][l]) for l in LEARNERS} for a in ARMS},
        "named": rows,
    }


def build(path, learner_path):
    reading = json.loads(pathlib.Path(path).read_text(encoding="utf-8"))
    species, panel = reading["species"], reading["panel"]
    for s in species:
        if sorted(s["arms"]) != sorted(ARMS) or sorted(s["by_dataset"]) != sorted(ARMS):
            raise SystemExit(f"{s['species']}: the reading's arms are not the document's")
    if sorted(panel["by_dataset"]) != sorted(ARMS):
        raise SystemExit("the panel's arms are not the document's")
    firm = [d for d in panel["decisions"] if d["firm"]]
    by_side = {k: sum(1 for d in firm if d["in_species"] == k) for k in ("improves", "worse")}
    if by_side != panel["firm"]:
        raise SystemExit(f"the firm decisions {by_side} do not add up to the reading's {panel['firm']}")
    learners = {}
    for d in firm:
        learners[learner(d["id"])] = learners.get(learner(d["id"]), 0) + 1
    p = panel["by_dataset"]
    split = by_learner(learner_path, species)
    out = {
        "source": "reference/input-comparison/second-reading.json",
        "margin": margin(),
        "species": [{"name": s["species"], "datasets": s["datasets"], "answered": s["datasets"] - len(s["unanswered"])} for s in species],
        "family": {"size": panel["family"]["size"], "tested": panel["family"]["tested"]},
        "firm": {"better": panel["firm"]["improves"], "worse": panel["firm"]["worse"], "learners": dict(sorted(learners.items()))},
        "learner_reading": {"source": "reference/input-comparison/learner-reading.json", "split": "lineage", "exploratory": split["exploratory"], "named": split["named"]},
        "predictions": {
            "species": panel["predictions"]["species"],
            "tests": [{k: t[k] for k in ("arm", "predicted", "better", "worse", "unchanged", "p_one_sided")} for t in panel["predictions"]["tests"]],
        },
        "arms": [
            {
                "id": a,
                "twin": species[0]["arms"][a]["twin"],
                "all": {**{k: p[a][k] for k in COUNTS}, "p": p[a]["p"], "median": p[a]["median"]},
                "species": [{k: s["by_dataset"][a][k] for k in COUNTS} for s in species],
                "contrasts": {k: sum(s["arms"][a][k] for s in species) for k in ("made", "inside", "beyond", "zero")},
                "learners": split["arms"][a],
            }
            for a in ARMS
        ],
    }
    # One line per species and arm, so a change in the reading reads as a change of those lines.
    lines = ["{"]
    for k in ("source", "margin", "species", "family", "firm", "learner_reading", "predictions"):
        lines.append(f"{json.dumps(k)}:{json.dumps(out[k], separators=(',', ':'))},")
    lines.append('"arms":[')
    lines += [json.dumps(a, separators=(",", ":")) + ("," if i < len(out["arms"]) - 1 else "") for i, a in enumerate(out["arms"])]
    lines.append("]}")
    text = "\n".join(lines) + "\n"
    assert json.loads(text) == out
    return text


def main(argv):
    args = argv[1:]
    if len(args) not in (2, 4) or (len(args) == 4 and args[2] != "--check"):
        raise SystemExit(__doc__)
    text = build(args[0], args[1])
    if len(args) == 4:
        if pathlib.Path(args[3]).read_text(encoding="utf-8") == text:
            return 0
        print("the bundled file differs from what the readings give", file=sys.stderr)
        return 1
    sys.stdout.write(text)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
