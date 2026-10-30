#!/usr/bin/env python3
"""Field-level diff of the rerun's prefix reports against revision N's MTB 2k candidate, every field.

usage: python3 -I diff_reference_prefix.py --prefix PREFIX_DATA_ROOT --reference REFERENCE_DATA_ROOT
           --out OUT.json [--slug mycobacterium_tuberculosis]

Plan revision 6, section 7b (pending reference). The cleaned prefix is revision N's own 2k
selection (prepare_execution_roots.py --cohorts cleaned), its matrix is byte-identical to revision
N's candidate's (section 7b's cmp), and both sides ran the same solver, EvalConfig and seeds from
the same inputs (section 7b's fingerprint view). So,
unlike diff_revision7_prefix.py, which had to let the lasso differ across a solver change, nothing
may differ here but wall-clock time:

- not compared: the report's top-level `seconds` and each learner row's `seconds`, which time the
  run and not its result;
- identical: everything else, recursively, including evaluation (revision, validity, reasons,
  scope), input_fingerprint, eval_config, samples, catalog, panels, mic, every split's folds and
  every learner row with its scores, metrics and lasso_fits, every model (lasso and
  models.lasso_fit included) and the feature cards. A field present on one side only is a
  difference.

Both roots must hold the same species/<slug>/genomes.txt, byte for byte (the reference must be
revision N's candidate for this very prefix), and every report on both sides must be schema 2 and
evaluation revision REVISION_N from revision_n.py, which stops with exit 2 while that is a
placeholder. Any difference exits 1; a different genome list, a missing or extra dataset, a wrong
revision or an unreadable report exits 2.

Lasso is included on purpose. Revision N's lasso can land on a different point of a degenerate
optimum (exact complement columns on a fold's fitting samples) with held-out scores up to about
0.1 apart at equal objective (the revision-9 note of 2026-10-29). Within one binary and one input
that should not happen, so a lasso-only difference here is a stop for review, not noise.
"""
import argparse
import importlib
import json
import sys
from pathlib import Path

# M1: python3 -I (3.11+ implies -P) leaves this script's directory off sys.path. Put back that one
# directory, first, and require revision_n to load from it.
LANE = Path(__file__).resolve().parent
if sys.path[:1] != [str(LANE)]:
    sys.path.insert(0, str(LANE))


def _sibling(name):
    module = importlib.import_module(name)
    if Path(module.__file__).resolve().parent != LANE:
        raise ImportError(f"{name} loaded from {module.__file__}, not from {LANE}")
    return module


try:
    rev_n = _sibling("revision_n")
except ImportError as _exc:
    print(f"diff_reference_prefix: STOP: cannot load revision_n: {_exc}", file=sys.stderr)
    sys.exit(2)

WALL_CLOCK = ["seconds", "splits[*].learners[*].seconds"]


class Stop(Exception):
    pass


def without_wall_clock(report):
    """The report without the fields WALL_CLOCK names; a copy, the input is not changed."""
    out = {k: v for k, v in report.items() if k != "seconds"}
    if isinstance(out.get("splits"), list):
        out["splits"] = [{**s, "learners": [{k: v for k, v in row.items() if k != "seconds"} if isinstance(row, dict) else row
                                            for row in s["learners"]]}
                         if isinstance(s, dict) and isinstance(s.get("learners"), list) else s for s in out["splits"]]
    return out


def differences(a, b, where):
    """Every path at which two JSON values differ, outermost first."""
    if type(a) is not type(b):
        return [where]
    if isinstance(a, dict):
        out = []
        for key in sorted(set(a) | set(b)):
            if key not in a or key not in b:
                out.append(f"{where}.{key}" if where else key)
            else:
                out.extend(differences(a[key], b[key], f"{where}.{key}" if where else key))
        return out
    if isinstance(a, list):
        if len(a) != len(b):
            return [f"{where} (length {len(a)} vs {len(b)})"]
        return [d for i, (x, y) in enumerate(zip(a, b)) for d in differences(x, y, f"{where}[{i}]")]
    return [] if a == b else [where]


def diff_report(prefix, reference, revision):
    for side, r in (("prefix", prefix), ("reference", reference)):
        if not isinstance(r, dict) or r.get("schema") != 2 or not isinstance(r.get("evaluation"), dict):
            raise Stop(f"{side} report {r.get('id') if isinstance(r, dict) else r!r}: not a schema-2 report")
        if r["evaluation"].get("revision") != revision:
            raise Stop(f"{side} report {r.get('id')}: evaluation revision {r['evaluation'].get('revision')!r}, not {revision}")
    found = differences(without_wall_clock(prefix), without_wall_clock(reference), "")
    return {"id": prefix.get("id"), "identical": not found, "differences": found[:200], "difference_count": len(found),
            "not_compared": WALL_CLOCK}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0], allow_abbrev=False)
    ap.add_argument("--prefix", required=True, type=Path)
    ap.add_argument("--reference", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--slug", default="mycobacterium_tuberculosis")
    args = ap.parse_args(argv)
    try:
        revision = rev_n.require(["REVISION_N"])["REVISION_N"]
        genomes = [root / "species" / args.slug / "genomes.txt" for root in (args.prefix, args.reference)]
        if genomes[0].read_bytes() != genomes[1].read_bytes():
            raise Stop(f"{genomes[1]} is not the prefix's genome list {genomes[0]}: not revision N's candidate for this prefix")
        dp = args.prefix / "species" / args.slug / "results"
        dr = args.reference / "species" / args.slug / "results"
        names_p = sorted(p.name for p in dp.glob("*.json"))
        names_r = sorted(p.name for p in dr.glob("*.json"))
        if not names_p or names_p != names_r:
            raise Stop(f"dataset sets differ: only the prefix {sorted(set(names_p) - set(names_r))}, "
                       f"only the reference {sorted(set(names_r) - set(names_p))}")
        rows = [diff_report(json.loads((dp / n).read_text()), json.loads((dr / n).read_text()), revision) for n in names_p]
    except (Stop, rev_n.Unfilled, rev_n.Malformed, OSError, KeyError, TypeError, ValueError) as err:
        print(f"diff_reference_prefix: STOP: {err}", file=sys.stderr)
        return 2
    out = {"prefix_root": str(args.prefix), "reference_root": str(args.reference), "evaluation_revision": revision,
           "datasets": len(rows), "identical": sum(r["identical"] for r in rows), "rows": rows}
    with args.out.open("x") as f:
        json.dump(out, f, indent=1)
        f.write("\n")
    print(json.dumps({k: v for k, v in out.items() if k != "rows"}, indent=1))
    return 0 if out["identical"] == len(rows) else 1


if __name__ == "__main__":
    sys.exit(main())
