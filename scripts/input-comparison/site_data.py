"""site_data.py READING [--check FILE]: the part of the input comparison's second reading that the web app bundles.

READING is reference/input-comparison/second-reading.json, the JSON second_reading.py wrote for H-1's four
species. The script writes, as JSON, what the site's "Input settings against the default" section counts:

- each species with its datasets and the datasets with a lineage answer;
- per arm, in the runner document's order: the datasets better and worse over the four species with the sign
  test's p and the median of the dataset means, the same counts per species, and the family C contrasts made,
  inside the margin, beyond it and exactly zero, summed over the species;
- the panel family's size and tested contrasts, the firm decisions by direction and by learner, and the
  predictions;
- the log-loss margin, which the four species specs under scripts/input-comparison/h1 share.

It recomputes nothing the reading decided: every number is the reading's, or a sum or a count of its records.

  python3 -I -B scripts/input-comparison/site_data.py reference/input-comparison/second-reading.json > web/src/data/input-comparison.json
  python3 -I -B scripts/input-comparison/site_data.py reference/input-comparison/second-reading.json --check web/src/data/input-comparison.json

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


def margin():
    """The log-loss margin of the four species specs, which must agree."""
    values = {json.loads((HERE / "h1" / f"h1-{s}.json").read_text(encoding="utf-8"))["margin"]["log_loss"] for s in SPECS}
    if len(values) != 1:
        raise SystemExit(f"the species specs disagree on the log-loss margin: {sorted(values)}")
    return values.pop()


def learner(contrast_id):
    """'<dataset>/C:rf:lineage:default->arm:log_loss' -> 'rf'."""
    return contrast_id.split("/", 1)[1].split(":")[1]


def build(path):
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
    out = {
        "source": "reference/input-comparison/second-reading.json",
        "margin": margin(),
        "species": [{"name": s["species"], "datasets": s["datasets"], "answered": s["datasets"] - len(s["unanswered"])} for s in species],
        "family": {"size": panel["family"]["size"], "tested": panel["family"]["tested"]},
        "firm": {"better": panel["firm"]["improves"], "worse": panel["firm"]["worse"], "learners": dict(sorted(learners.items()))},
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
            }
            for a in ARMS
        ],
    }
    # One line per species and arm, so a change in the reading reads as a change of those lines.
    lines = ["{"]
    for k in ("source", "margin", "species", "family", "firm", "predictions"):
        lines.append(f"{json.dumps(k)}:{json.dumps(out[k], separators=(',', ':'))},")
    lines.append('"arms":[')
    lines += [json.dumps(a, separators=(",", ":")) + ("," if i < len(out["arms"]) - 1 else "") for i, a in enumerate(out["arms"])]
    lines.append("]}")
    text = "\n".join(lines) + "\n"
    assert json.loads(text) == out
    return text


def main(argv):
    args = argv[1:]
    if len(args) not in (1, 3) or (len(args) == 3 and args[1] != "--check"):
        raise SystemExit(__doc__)
    text = build(args[0])
    if len(args) == 3:
        if pathlib.Path(args[2]).read_text(encoding="utf-8") == text:
            return 0
        print("the bundled file differs from what the reading gives", file=sys.stderr)
        return 1
    sys.stdout.write(text)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
