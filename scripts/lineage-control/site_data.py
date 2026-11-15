#!/usr/bin/env python3
"""site_data.py TABLE [--check FILE]: the within-lineage verdicts the web app bundles.

TABLE is a file of reference/lineage-control that verdicts.py wrote, named
after its release. The script writes, as JSON, what the dataset inspector
needs to attach a verdict to a rule term: for each dataset its report's input
fingerprint, and for each rule term, in order, its feature, its presence flag
and its verdict. The inspector shows a verdict only on a report of the same
evaluation revision and fingerprint whose rule at that position is the same,
so a later release shows none until its own table is written.

It tests nothing: the counts and verdicts are the table's.

  python3 -I scripts/lineage-control/site_data.py reference/lineage-control/rev14-274-20261112.tsv > web/src/data/lineage-control.json
  python3 -I scripts/lineage-control/site_data.py reference/lineage-control/rev14-274-20261112.tsv --check web/src/data/lineage-control.json

With --check it prints nothing and exits 1 unless FILE is exactly the text
the script prints.
"""
import csv
import importlib.util
import json
import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("lineage_control_verdicts", HERE / "verdicts.py")
verdicts = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(verdicts)


def build(table):
    """The bundle's text for one verdict table."""
    table = pathlib.Path(table)
    with table.open(encoding="utf-8", newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    if not rows:
        raise SystemExit(f"{table}: no rows")
    revisions = {r["revision"] for r in rows}
    if len(revisions) != 1:
        raise SystemExit(f"{table}: reports of revisions {sorted(revisions)}; a table holds one release")
    datasets = {}
    counts = dict.fromkeys(verdicts.VERDICTS, 0)
    for r in rows:
        if r["verdict"] not in counts:
            raise SystemExit(f"{table}: {r['dataset']} term {r['position']} has the verdict {r['verdict']!r}")
        entry = datasets.setdefault(r["dataset"], {"fingerprint": r["input_fingerprint"], "terms": []})
        if entry["fingerprint"] != r["input_fingerprint"]:
            raise SystemExit(f"{table}: {r['dataset']} has two input fingerprints")
        if int(r["position"]) != len(entry["terms"]) + 1:
            raise SystemExit(f"{table}: {r['dataset']} term {r['position']} is out of order")
        entry["terms"].append([int(r["feature"]), r["presence"] == "1", r["verdict"]])
        counts[r["verdict"]] += 1
    data = {
        "release": table.stem,
        "evaluation_revision": int(revisions.pop()),
        "min_strata": verdicts.MIN_STRATA,
        "alpha": verdicts.ALPHA,
        "terms": len(rows),
        "counts": counts,
        "datasets": dict(sorted(datasets.items())),
    }
    return json.dumps(data, separators=(",", ":")) + "\n"


def main(argv):
    if len(argv) == 2:
        sys.stdout.write(build(argv[1]))
        return 0
    if len(argv) == 4 and argv[2] == "--check":
        return 0 if pathlib.Path(argv[3]).read_text(encoding="utf-8") == build(argv[1]) else 1
    raise SystemExit(__doc__)


if __name__ == "__main__":
    sys.exit(main(sys.argv))
