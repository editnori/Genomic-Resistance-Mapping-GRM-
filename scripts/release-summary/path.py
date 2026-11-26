"""path.py INDEX.json INTERVALS_DIR

Reads the path from lab records to datasets that a release publishes in its index: the counts of the index's `funnel`
block, sums over its `species` rows and sums over its `datasets` rows. Prints one JSON document on stdout. The script
reads the index only. It takes the intervals directory so that every script has the same command line.

Revisions understood: 12 and 14, read from the index's evaluation_revision. Any other revision stops the script with a
message and exit status 1.

The output holds:

  revision, datasets (the index's dataset rows) and species_rows (its species rows);
  funnel: the `funnel` block's counts (src/curate.rs:474-506) under plain keys. Records dropped by name and genomes
      failing a check each give their parts and their total; a reason a genome fails a check is written as the index
      writes it. With no `funnel` block in the index the block is null, never zeros;
  species: sums over the species rows (src/report.rs:984-995) of datasets_curated, datasets_reported, genomes_labeled
      and genomes_analyzed, the highest genomes_analyzed of any species, and how many species sit at it;
  dataset_rows: sums over the dataset rows of n and n_resistant;
  checks: two checks of the funnel block, each with both sides, or null with no funnel block. Genomes with a record
      equal those that pass the checks plus those without metadata plus those that fail a check, and the five kinds
      of pair add up to the pairs. A failed check stops the script with exit status 1 and names both sides. A release
      whose label rule sets pairs aside for an intermediate call (pairs_intermediate_mixed, src/curate.rs:490-495)
      fails the second check, since the script reads the five kinds only.

A value that is not a count (a whole number, 0 or more) stops the script with a message and exit status 1."""
import json
import sys
from pathlib import Path

REVISIONS = {12, 14}

index = json.loads(Path(sys.argv[1]).read_text())
revision = index.get("evaluation_revision")
if revision not in REVISIONS:
    sys.exit(f"path.py: {sys.argv[1]} is evaluation revision {revision!r}; "
             f"this script understands revisions {', '.join(map(str, sorted(REVISIONS)))}")


def count(block, key, where):
    v = block.get(key) if isinstance(block, dict) else None
    if isinstance(v, bool) or not isinstance(v, int) or v < 0:
        sys.exit(f"path.py: {sys.argv[1]}: {where} has {key} {v!r}, which is not a count")
    return v


def rows_of(key):
    rows = index.get(key)
    if not isinstance(rows, list):
        sys.exit(f"path.py: {sys.argv[1]}: the index's {key} is {rows!r}, which is not a list of rows")
    return rows


datasets, species = rows_of("datasets"), rows_of("species")
out = {"revision": revision, "datasets": len(datasets), "species_rows": len(species), "funnel": None}

f = index.get("funnel")
if f is not None:
    c = lambda key: count(f, key, "funnel")
    by_name = {"drug_class": c("records_class_label"), "not_an_antibiotic": c("records_not_antibiotic"),
               "name_not_resolved": c("records_unresolved_name")}
    failed = f.get("genomes_failed_qc") if isinstance(f, dict) else None
    if not isinstance(failed, dict):
        sys.exit(f"path.py: {sys.argv[1]}: funnel has genomes_failed_qc {failed!r}, which is not a block of counts")
    reasons = {why: count(failed, why, "funnel genomes_failed_qc") for why in failed}
    kinds = {"reported_resistant_or_susceptible": c("pairs_reported"), "label_from_mic": c("pairs_labeled_from_mic"),
             "no_label": c("pairs_unlabeled"), "intermediate_only": c("pairs_intermediate_only"),
             "conflict": c("pairs_conflict")}
    out["funnel"] = {
        "records": c("records"),
        "records_dropped_by_name": {**by_name, "total": sum(by_name.values())},
        "records_no_call_no_mic": c("records_no_measure"),
        "genomes_with_a_record": c("genomes_with_records"),
        "genomes_without_metadata": c("genomes_without_metadata"),
        "genomes_failing_a_check": {"reasons": reasons, "total": sum(reasons.values())},
        "genomes_passing": c("genomes_passing_qc"),
        "pairs": c("pairs"),
        "pair_kinds": kinds,
    }

names = {"datasets_curated": "datasets_curated", "datasets_with_report": "datasets_reported",
         "genomes_with_any_label": "genomes_labeled", "genomes_analyzed": "genomes_analyzed"}
per_row = []
for i, r in enumerate(species):
    where = f"species row {(r.get('slug') if isinstance(r, dict) else None) or i}"
    per_row.append({plain: count(r, key, where) for plain, key in names.items()})
analyzed = [r["genomes_analyzed"] for r in per_row]
top = max(analyzed) if analyzed else None
out["species"] = {**{plain: sum(r[plain] for r in per_row) for plain in names},
                  "genomes_analyzed_highest": top, "species_at_highest": analyzed.count(top) if analyzed else None}

out["dataset_rows"] = {key: sum(count(d, key, f"dataset {d.get('id') if isinstance(d, dict) else None}")
                                for d in datasets) for key in ("n", "n_resistant")}

out["checks"] = None
if out["funnel"] is not None:
    p = out["funnel"]
    out["checks"] = {
        "genomes": {"with_a_record": p["genomes_with_a_record"],
                    "passing_plus_without_metadata_plus_failing":
                        p["genomes_passing"] + p["genomes_without_metadata"] + p["genomes_failing_a_check"]["total"]},
        "pairs": {"pairs": p["pairs"], "sum_of_five_kinds": sum(p["pair_kinds"].values())},
    }
    for name, sides in out["checks"].items():
        (ka, a), (kb, b) = sides.items()
        if a != b:
            sys.exit(f"path.py: {sys.argv[1]}: the {name} check fails: {ka} is {a} and {kb} is {b}")
json.dump(out, sys.stdout, indent=1)
