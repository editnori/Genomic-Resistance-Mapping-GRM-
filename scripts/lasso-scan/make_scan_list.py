"""Write the lasso scan list (scan-list.json) from the revision-7 reports.

usage: python3 -I -B make_scan_list.py CANDIDATE_ROOT > scan-list.json

CANDIDATE_ROOT holds species/<slug>/results/<id>.json for the 275
revision-7 candidate reports, such as
data/backups/revision7-all23-candidate-20261027. The scan set is chosen by
measured revision-7 lasso cost: the five species with the most lasso
seconds, and in each its three costliest datasets plus its datasets with the
most samples and the most features. That is 18 datasets, listed most
expensive first, the order the scan runs them in.
"""
import json
import sys
from pathlib import Path

SPECIES = 5
COSTLIEST = 3


def lasso_timing(report):
    seconds, fits = 0.0, 0
    for split in report["splits"]:
        if split.get("unavailable_reason"):
            continue
        (row,) = [r for r in split["learners"] if r["learner"] == "lasso"]
        seconds += row["seconds"]
        fits += split["folds"]
    return seconds, fits


def main(root):
    rows = []
    for path in sorted(Path(root).glob("species/*/results/*.json")):
        report = json.loads(path.read_text())
        assert report["evaluation"]["revision"] == 7, path
        seconds, fits = lasso_timing(report)
        rows.append({
            "species": report["species_slug"],
            "id": report["id"],
            "samples": report["counts"]["samples"],
            "features": report["counts"]["features"],
            "lasso_seconds": seconds,
            "outer_fits": fits,
            "report_seconds": report["seconds"],
            "input_fingerprint": report["input_fingerprint"],
        })
    assert len(rows) == 275, len(rows)
    for rank, row in enumerate(sorted(rows, key=lambda r: -r["lasso_seconds"]), 1):
        row["rank"] = rank
    by_species = {}
    for row in rows:
        by_species.setdefault(row["species"], []).append(row)
    heaviest = sorted(by_species, key=lambda s: -sum(r["lasso_seconds"] for r in by_species[s]))[:SPECIES]
    chosen = {}
    for slug in heaviest:
        ranked = sorted(by_species[slug], key=lambda r: -r["lasso_seconds"])
        picks = [(r, f"cost rank {i + 1} in species") for i, r in enumerate(ranked[:COSTLIEST])]
        picks.append((max(ranked, key=lambda r: r["samples"]), "most samples in species"))
        picks.append((max(ranked, key=lambda r: r["features"]), "most features in species"))
        for row, why in picks:
            chosen.setdefault(row["id"], dict(row, chosen_for=[]))["chosen_for"].append(why)
    datasets = sorted(chosen.values(), key=lambda r: r["rank"])
    return {
        "description": "The large-species lasso scan list; docs/lasso-stopping.md registers its stop rule.",
        "selection": __doc__.split("The scan set is ")[1].split("\n\n")[0].replace("\n", " ").strip(),
        "revision7": {
            "source": "the 275 revision-7 candidate reports in data/backups/revision7-all23-candidate-20261027",
            "reports": len(rows),
            "lasso_seconds": sum(r["lasso_seconds"] for r in rows),
            "report_seconds": sum(r["report_seconds"] for r in rows),
            "species_by_lasso_seconds": heaviest,
        },
        "datasets": datasets,
    }


if __name__ == "__main__":
    print(json.dumps(main(sys.argv[1]), indent=1))
