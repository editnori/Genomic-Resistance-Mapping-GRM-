"""site_data.py ROOT [--check FILE]: the outside-study rows the web app bundles.

ROOT is reference/outside-study. The script reads it with tables.py's own
reading rules (load, holds, frozen_and_scored) and writes, as JSON, what the
site's "Newly deposited genomes" section counts and draws:

- for each population, each run's frozen and scored genomes; its datasets with
  their species and drug (tables.py's drug names), whether they have estimates
  and the top study's share of the resistant genomes; and every row with
  estimates for the learners the documents read: its dataset, its learner and
  its two readings;
- for the set covering machine's rows, the report's lineage-split figures and
  the new estimates with their 95% intervals, as the files write them.

It estimates nothing and adds no population to another. The site computes
every count from these rows at render time.

  python3 -I scripts/outside-study-results/site_data.py reference/outside-study > web/src/data/outside-study.json
  python3 -I scripts/outside-study-results/site_data.py reference/outside-study --check web/src/data/outside-study.json

With --check it prints nothing and exits 1 unless FILE is exactly the text
the script prints.
"""
import gzip
import importlib.util
import json
import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("outside_study_tables", HERE / "tables.py")
tables = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(tables)

# The learner the site draws, whose figures the rows carry.
DRAWN = "scm"


def number(text):
    """A cell as a number, or None for an empty cell."""
    return None if text == "" else float(text)


def species_name(slug):
    """'staphylococcus_aureus' -> 'Staphylococcus aureus', as the site's index names a species."""
    genus, _, rest = slug.partition("_")
    return f"{genus.capitalize()} {rest.replace('_', ' ')}"


def run_header(folder):
    """The run's evaluation revision and interval level, from its results.json.gz."""
    with gzip.open(folder / "results.json.gz", "rt", encoding="utf-8") as f:
        r = json.load(f)
    return r["evaluation_revision"], r["confidence"]


def build(root):
    root = pathlib.Path(root)
    data = tables.load(root)
    revisions, levels, populations = set(), set(), []
    for pop, _ in tables.POPULATIONS:
        if pop not in data:
            continue
        species, datasets, rows = [], [], []
        for slug, (reading, ds, metrics, _) in data[pop].items():
            revision, level = run_header(root / pop / slug)
            revisions.add(revision)
            levels.add(level)
            frozen, scored, _ = tables.frozen_and_scored(reading)
            species.append({"slug": slug, "frozen": frozen, "scored": scored})
            for d in ds:
                datasets.append(
                    {
                        "id": d["dataset"],
                        "species": species_name(slug),
                        "antibiotic": tables.drug(d["dataset"]),
                        "estimates": d["reportable"] == "yes",
                        "top_study_resistant_share": float(d["top_study_resistant_share"]),
                    }
                )
            for r in metrics:
                if r["learner"] not in tables.LEARNERS or r["reportable"] != "yes":
                    continue
                row = {"dataset": r["dataset"], "learner": r["learner"], "sensitivity": r["reading_sensitivity"], "specificity": r["reading_specificity"]}
                if r["learner"] == DRAWN:
                    row["figures"] = {
                        k: {
                            "report": number(r["lineage_" + k]),
                            "estimate": float(r[k]),
                            "lo": float(r[k + "_lo"]),
                            "hi": float(r[k + "_hi"]),
                        }
                        for k in ("sensitivity", "specificity")
                    }
                rows.append(row)
        populations.append({"id": pop, "species": species, "datasets": datasets, "rows": rows})
    if len(revisions) != 1 or len(levels) != 1:
        raise SystemExit(f"the runs disagree on the evaluation revision or the interval level: {sorted(revisions)} {sorted(levels)}")
    out = {
        "source": "reference/outside-study",
        "evaluation_revision": revisions.pop(),
        "confidence": levels.pop(),
        "populations": populations,
    }
    # One line per run, dataset and row, so a change in the files reads as a
    # change of those lines.
    lines = ["{"]
    for k in ("source", "evaluation_revision", "confidence"):
        lines.append(f"{json.dumps(k)}:{json.dumps(out[k], separators=(',', ':'))},")
    lines.append('"populations":[')
    for i, p in enumerate(out["populations"]):
        lines.append("{" + f'"id":{json.dumps(p["id"])},')
        for key in ("species", "datasets", "rows"):
            lines.append(f'"{key}":[')
            items = p[key]
            lines += [json.dumps(x, separators=(",", ":")) + ("," if j < len(items) - 1 else "") for j, x in enumerate(items)]
            lines.append("]" + ("," if key != "rows" else ""))
        lines.append("}" + ("," if i < len(out["populations"]) - 1 else ""))
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
        print("the bundled file differs from what the output files give", file=sys.stderr)
        return 1
    sys.stdout.write(text)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
