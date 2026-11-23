"""drug_rules.py: the drug map the AMRFinderPlus comparison used, as a file.

    python3 -I drug_rules.py --reports DIR --lists lists.tsv --out drug_rules.tsv
    python3 -I drug_rules.py --reports DIR --lists lists.tsv --check drug_rules.tsv

Each release report holds its published rule, `catalog.rule`, the output of
catalog::drug_rule (src/catalog.rs). This writes one row per dataset of A
(scripts/c3-inputs/lists.tsv, in_A): dataset, species, antibiotic, the
report's drug_class and the rule, groups joined by " & " and tokens of a
group by " | ", in the catalogue's own class and subclass names. With
--check it stops unless the file equals what the reports give.

Python 3.12, standard library only; run with -I.
"""

import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import amr_rule as ar  # noqa: E402
from common import read_tsv, write_tsv  # noqa: E402

HEADER = ["dataset", "species_slug", "species", "antibiotic", "drug_class", "rule"]


def rows_from(reports, lists, species=None):
    out = []
    for r in sorted(read_tsv(lists), key=lambda r: r["dataset"]):
        if r["in_A"] != "1" or (species and r["species_slug"] != species):
            continue
        with open(os.path.join(reports, r["dataset"] + ".json"), encoding="utf-8") as f:
            rep = json.load(f)
        if rep.get("id") != r["dataset"]:
            ar.refuse(f"{r['dataset']}: the report file holds another dataset")
        if not rep.get("catalog") or not rep["catalog"].get("rule"):
            ar.refuse(f"{r['dataset']}: a dataset of A has no published rule")
        rule = rep["catalog"]["rule"]
        if ar.parse_rule(ar.format_rule(rule)) != rule:
            ar.refuse(f"{r['dataset']}: the rule does not survive the file's form")
        out.append({"dataset": r["dataset"], "species_slug": rep["species_slug"], "species": rep.get("species", ""),
                    "antibiotic": rep.get("antibiotic", ""), "drug_class": rep.get("drug_class", ""),
                    "rule": ar.format_rule(rule)})
    return out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--reports", required=True)
    ap.add_argument("--lists", required=True)
    ap.add_argument("--species")
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--out")
    g.add_argument("--check")
    args = ap.parse_args(argv)
    try:
        rows = rows_from(args.reports, args.lists, args.species)
        if args.out:
            write_tsv(args.out, HEADER, rows)
            print(f"drug_rules: {len(rows)} datasets written")
            return 0
        have = {r["dataset"]: r for r in read_tsv(args.check)}
        bad = [r["dataset"] for r in rows if {k: have.get(r["dataset"], {}).get(k) for k in HEADER} != r]
        if bad:
            ar.refuse(f"{len(bad)} datasets differ from {os.path.basename(args.check)}: {','.join(bad)}")
        print(f"drug_rules: {len(rows)} datasets equal the file")
        return 0
    except ar.Refusal as e:
        print(f"drug_rules.py: refused: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
