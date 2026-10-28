#!/usr/bin/env python3
"""Field-level diff of the prefix cohort's revision-8 reports against the revision-7 2k candidate.

usage: python3 -I diff_revision7_prefix.py --revision8 PREFIX_DATA_ROOT --revision7 CANDIDATE_DATA_ROOT
           --out OUT.json [--slug mycobacterium_tuberculosis]

The prefix cohort's matrix is byte-identical to the frozen 2k matrix, and revision 8 changed only
the lasso (docs/lasso-stopping.md, step 4). A whole-file comparison cannot show that, because
`models.lasso_fit` sits before `rf_importance`, so this compares field by field. For each of the
datasets in both roots:

- identical: input_fingerprint (the fingerprint excludes code, and its function is unchanged since
  revision 7), eval_config, identity, counts, lineage threshold, catalog groups, samples, catalog,
  panels, unsupported assays and mic; each split's scheme, fold count, group count, test_fold and
  unavailable reason; every learner row other than lasso (scm, cart, rf, gbt, catalog), apart from
  its wall-clock `seconds`; and the final scm, cart, cart_depth, rf_importance and gbt_importance;
- allowed to differ: the lasso rows, models.lasso and models.lasso_fit (revision 8 adds the solver
  evidence); revision-8 lasso rows must carry lasso_fits;
- not compared: features (feature cards follow the lasso coefficients too) and seconds.

The revision-8 side must be evaluation revision 8 and the other revision 7, both schema 2. Any
difference in a field that must be identical exits 1; a missing dataset or wrong revision exits 2.
"""
import argparse, json, sys
from pathlib import Path

SAME_TOP = ["input_fingerprint", "eval_config", "id", "species", "species_slug", "antibiotic", "drug_class", "breakpoint", "counts",
            "lineage_threshold", "catalog_groups", "samples", "catalog", "panels", "unsupported_assays", "mic"]
SAME_SPLIT = ["scheme", "folds", "n_groups", "test_fold", "unavailable_reason"]
SAME_MODELS = ["scm", "cart", "cart_depth", "rf_importance", "gbt_importance"]
NOT_COMPARED = ["features", "seconds"]


class Stop(Exception):
    pass


def without_seconds(row):
    return {k: v for k, v in row.items() if k != "seconds"}


def first_difference(a, b, where):
    """Path of the first difference between two JSON values, or None."""
    if type(a) is not type(b):
        return where
    if isinstance(a, dict):
        for key in sorted(set(a) | set(b)):
            if key not in a or key not in b:
                return f"{where}.{key}"
            d = first_difference(a[key], b[key], f"{where}.{key}")
            if d:
                return d
        return None
    if isinstance(a, list):
        if len(a) != len(b):
            return f"{where} (length {len(a)} vs {len(b)})"
        for i, (x, y) in enumerate(zip(a, b)):
            d = first_difference(x, y, f"{where}[{i}]")
            if d:
                return d
        return None
    return None if a == b else where


def diff_report(r8, r7):
    if r8.get("schema") != 2 or r7.get("schema") != 2:
        raise Stop(f"{r8.get('id')}: both reports must be schema 2")
    if r8["evaluation"]["revision"] != 8 or r7["evaluation"]["revision"] != 7:
        raise Stop(f"{r8.get('id')}: expected revision 8 against revision 7, got {r8['evaluation']['revision']} and {r7['evaluation']['revision']}")
    differences = []
    for key in SAME_TOP:
        if (key in r8) != (key in r7):
            differences.append(key)
        elif key in r8:
            d = first_difference(r8[key], r7[key], key)
            if d:
                differences.append(d)
    s8 = {s["scheme"]: s for s in r8["splits"]}
    s7 = {s["scheme"]: s for s in r7["splits"]}
    if set(s8) != set(s7):
        differences.append("splits (schemes)")
    lasso_rows = 0
    for scheme in sorted(set(s8) & set(s7)):
        a, b = s8[scheme], s7[scheme]
        for key in SAME_SPLIT:
            if a.get(key) != b.get(key):
                differences.append(f"splits.{scheme}.{key}")
        la = {row["learner"]: row for row in a["learners"]}
        lb = {row["learner"]: row for row in b["learners"]}
        if list(la) != list(lb):
            differences.append(f"splits.{scheme}.learners (names)")
        for name in la:
            if name == "lasso":
                lasso_rows += 1
                if "lasso_fits" not in la[name]:
                    differences.append(f"splits.{scheme}.learners.lasso.lasso_fits (missing at revision 8)")
                continue
            if name in lb:
                d = first_difference(without_seconds(la[name]), without_seconds(lb[name]), f"splits.{scheme}.learners.{name}")
                if d:
                    differences.append(d)
    for key in SAME_MODELS:
        d = first_difference(r8["models"].get(key), r7["models"].get(key), f"models.{key}")
        if d:
            differences.append(d)
    if "lasso_fit" not in r8["models"]:
        differences.append("models.lasso_fit (missing at revision 8)")
    return {"id": r8["id"], "identical_required_fields": not differences, "differences": differences,
            "lasso_rows_allowed_to_differ": lasso_rows, "not_compared": NOT_COMPARED}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--revision8", required=True, type=Path)
    ap.add_argument("--revision7", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--slug", default="mycobacterium_tuberculosis")
    args = ap.parse_args(argv)
    try:
        d8 = args.revision8 / "species" / args.slug / "results"
        d7 = args.revision7 / "species" / args.slug / "results"
        names8 = sorted(p.name for p in d8.glob("*.json"))
        names7 = sorted(p.name for p in d7.glob("*.json"))
        if not names8 or names8 != names7:
            raise Stop(f"dataset sets differ: only revision 8 {sorted(set(names8) - set(names7))}, only revision 7 {sorted(set(names7) - set(names8))}")
        rows = [diff_report(json.loads((d8 / n).read_text()), json.loads((d7 / n).read_text())) for n in names8]
    except (Stop, OSError, KeyError, ValueError) as err:
        print(f"diff_revision7_prefix: STOP: {err}", file=sys.stderr)
        return 2
    out = {"revision8_root": str(args.revision8), "revision7_root": str(args.revision7), "datasets": len(rows),
           "identical": sum(r["identical_required_fields"] for r in rows), "rows": rows}
    with args.out.open("x") as f:
        json.dump(out, f, indent=1)
        f.write("\n")
    print(json.dumps({k: v for k, v in out.items() if k != "rows"}, indent=1))
    return 0 if out["identical"] == len(rows) else 1


if __name__ == "__main__":
    sys.exit(main())
