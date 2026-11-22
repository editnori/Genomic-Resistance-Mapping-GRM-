"""site_data.py ROOT [--check FILE]: the other-database summaries the web app bundles.

ROOT is reference/other-database. The script reads aggregate TSV files only
(species.tsv, totals.tsv, learners.tsv) and writes, as JSON, the species counts,
learner totals, and per-species rule results for the site's other-database
section. It never reads genomes.tsv and holds no genome or BioSample row.

  python3 -I scripts/other-database/site_data.py reference/other-database > web/src/data/other-database.json
  python3 -I scripts/other-database/site_data.py reference/other-database --check web/src/data/other-database.json
"""

import csv
import json
import pathlib
import sys

LEARNERS = ("scm", "cart", "lasso")


def load_tsv(path):
    with pathlib.Path(path).open(encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def species_name(slug):
    genus, _, rest = slug.partition("_")
    return f"{genus.capitalize()} {rest.replace('_', ' ')}"


def opt_float(text):
    return None if text == "" else float(text)


def opt_int(text):
    return None if text == "" else int(text)


def build(root):
    root = pathlib.Path(root)
    sp_rows = load_tsv(root / "species.tsv")
    tot_rows = load_tsv(root / "totals.tsv")
    lr_rows = load_tsv(root / "learners.tsv")

    species = [
        {
            "slug": r["species"],
            "name": species_name(r["species"]),
            "part": r["part"],
            "scored": int(r["scored"]),
            "reports_scored": int(r["reports_scored"]),
            "shared_bioproject": int(r["shared_bioproject"]),
            "largest_bioproject": r["largest_bioproject"],
            "largest_bioproject_share": float(r["largest_bioproject_share"]),
            "largest_bioproject_training": int(r["largest_bioproject_training"]),
        }
        for r in sp_rows
    ]

    tot_by_learner = {r["learner"]: r for r in tot_rows}
    totals = {}
    for l in LEARNERS:
        r = tot_by_learner[l]
        totals[l] = {
            "species": int(r["species"]),
            "plan_meets": int(r["plan_meets"]),
            "plan_fails": int(r["plan_fails"]),
            "plan_not_assessable": int(r["plan_not_assessable"]),
            "estimable": int(r["estimable"]),
            "holds": int(r["holds"]),
            "recall_below": int(r["recall_below"]),
            "log_loss_worse": int(r["log_loss_worse"]),
            "specificity_under_half": int(r["specificity_under_half"]),
        }

    by_species = [
        {
            "species": r["species"],
            "name": species_name(r["species"]),
            "part": r["part"],
            "learner": r["learner"],
            "estimable": int(r["estimable"]),
            "holds": int(r["holds"]),
            "recall_below": int(r["recall_below"]),
            "log_loss_worse": int(r["log_loss_worse"]),
            "specificity_under_half": int(r["specificity_under_half"]),
            "lowest_specificity": opt_float(r["lowest_specificity"]),
            "lowest_specificity_dataset": r["lowest_specificity_dataset"] or None,
            "rule": r["rule"],
            "margin": opt_int(r["margin"]),
        }
        for r in lr_rows
        if r["learner"] in LEARNERS
    ]

    out = {
        "source": "reference/other-database",
        "evaluation_revision": 14,
        "species": species,
        "totals": totals,
        "by_species": by_species,
    }
    text = json.dumps(out, separators=(",", ":"), ensure_ascii=True) + "\n"
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
        print("the bundled file differs from what the record gives", file=sys.stderr)
        return 1
    sys.stdout.write(text)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
