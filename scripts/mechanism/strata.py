"""strata.py: block M1 of the mechanism split, recall by mechanism, for one species.

    python3 -I strata.py --species SLUG --reports DIR --lists lists.tsv \
        --drug-rules drug_rules.tsv --tables DIR --calls calls.tsv \
        --callers callers.tsv --out DIR

Reads the release's reports of the species' datasets in A
(scripts/c3-inputs/lists.tsv), the AMRFinderPlus program's tables
(DIR/<genome id>.tsv), the kit adapter's calls.tsv and the comparison's
record callers.tsv. Fits nothing. Writes aggregate tables only: counts and
recalls per dataset, split, row and stratum, with species, antibiotic and
drug class in full, and no genome id.

Checks, before any stratum is counted (PLAN.md, M1):
- each report's catalog.rule equals its row of drug_rules.tsv;
- the hard check: for each split and each row with metrics (the five
  learners and the catalog row, GRM's known-gene rule), tp, fn, tn and fp
  recomputed from `scores` at 0.5 or more equal the report's `metrics`, and
  the catalog row's calls equal samples[].catalog_resistant;
- the program's call made here equals the adapter's tool_default call on
  every genome, and its counts equal the record's program:tool_default row.
A dataset that fails any check is named in datasets.tsv and left out.

Strata, from the counted rows that meet the dataset's rule, never from the
label: point only, gene only, both, neither; beside them, partial only and
with a POINT_DISRUPT row, and all genomes.

Python 3.12, standard library only; run with -I.
"""

import argparse
import collections
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import amr_rule as ar  # noqa: E402
from common import RECALL_FLOOR, read_tsv, recall_cell, shown, write_tsv  # noqa: E402

LEARNERS = ("scm", "cart", "rf", "gbt", "lasso")
ROWS = LEARNERS + ("catalog", "program")
SPLITS = ("lineage", "random")
CUTOFF = 0.5  # predict::CALL_CUTOFF
STRATA = ("point_only", "gene_only", "both", "neither")
BESIDE = ("partial_only", "with_point_disrupt", "all")
ROW_NAMES = {"catalog": "GRM's known-gene rule (catalog row)", "program": "the AMRFinderPlus program"}

STRATA_HEADER = ["dataset", "species_slug", "species", "antibiotic", "drug_class", "set", "labels_depend_on",
                 "split", "row", "stratum", "kind", "resistant", "resistant_called", "susceptible",
                 "susceptible_called", "recall"]
POOLED_HEADER = ["species_slug", "species", "drug_class", "set", "split", "row", "stratum", "datasets",
                 "resistant_pairs", "resistant_called", "susceptible_pairs", "susceptible_called", "recall"]
DATASETS_HEADER = ["dataset", "species_slug", "species", "antibiotic", "drug_class", "set", "labels_depend_on",
                   "genomes", "with_table", "rule_check", "hard_check", "learner_rows_equal", "learner_rows",
                   "program_check", "record_check", "entered", "why"]
NOCALL_HEADER = ["dataset", "split", "learner", "no_call_resistant", "no_call_susceptible"]


def confusion(truth, calls):
    c = {"tp": 0, "fn": 0, "tn": 0, "fp": 0}
    for t, p in zip(truth, calls):
        c[("tp" if p else "fn") if t else ("fp" if p else "tn")] += 1
    return c


def calls_of(scores):
    out = []
    for s in scores:
        if not isinstance(s, (int, float)) or s != s:
            return None
        out.append(s >= CUTOFF)
    return out


def split_rows(report):
    """{scheme: {learner: row}} for the splits the release holds."""
    out = {}
    for s in report["splits"]:
        out[s["scheme"]] = {lr["learner"]: lr for lr in s["learners"]}
    return out


def hard_check(report):
    """(passed, equal rows, rows, why, row calls {(split, row): [bool]}, no_call rows)."""
    truth = [s["resistant"] for s in report["samples"]]
    n = len(truth)
    by = split_rows(report)
    calls, nocall, why = {}, [], []
    equal = total = 0
    for split in SPLITS:
        if split not in by:
            why.append(f"no {split} split")
            continue
        for name in LEARNERS + ("catalog",):
            lr = by[split].get(name)
            total += 1
            if lr is None:
                why.append(f"{split}: no {name} row")
                continue
            sc = lr.get("scores") or []
            c = calls_of(sc) if len(sc) == n else None
            m = lr.get("metrics")
            if c is None or m is None:
                why.append(f"{split} {name}: scores or metrics missing")
                continue
            got = confusion(truth, c)
            if any(got[k] != m.get(k) for k in ("tp", "fn", "tn", "fp")):
                why.append(f"{split} {name}: counts from scores differ from metrics")
                continue
            if name == "catalog" and c != [s.get("catalog_resistant") for s in report["samples"]]:
                why.append(f"{split} catalog: calls differ from samples[].catalog_resistant")
                continue
            equal += 1
            calls[(split, name)] = c
            if name in LEARNERS:
                nc = lr.get("no_call") or {"resistant": 0, "susceptible": 0}
                nocall.append({"split": split, "learner": name, "no_call_resistant": int(nc["resistant"]),
                               "no_call_susceptible": int(nc["susceptible"])})
    return not why and equal == total, equal, total, "; ".join(why), calls, nocall


def stratum_of(rows, rule):
    """(stratum, partial only, with POINT_DISRUPT) of one genome in one
    dataset, from its table rows and the dataset's rule alone."""
    met = ar.rule_rows(rows, rule)
    point = any(ar.is_point(r) for r in met)
    gene = any(not ar.is_point(r) for r in met)
    stratum = "both" if point and gene else "point_only" if point else "gene_only" if gene else "neither"
    partial_only = bool(met) and all(ar.is_partial(r) for r in met)
    disrupt = any(r["subtype"] == "POINT_DISRUPT" for r in met)
    return stratum, partial_only, disrupt


def genome_strata(samples, tables, rule):
    """Per sample, (stratum, partial only, with POINT_DISRUPT), or None
    where the genome has no table (the program failed on it)."""
    return [None if tables.get(s["id"]) is None else stratum_of(tables[s["id"]], rule) for s in samples]


def dataset_tables(report, tables, rule, row_calls, program, keep):
    """Count rows per split, row and stratum. `keep` marks the samples with a
    table; every count is over those samples."""
    truth = [s["resistant"] for s in report["samples"]]
    strata = genome_strata(report["samples"], tables, rule)
    out = []
    for split in SPLITS:
        for row in ROWS:
            c = program if row == "program" else row_calls.get((split, row))
            if c is None:
                continue
            cells = collections.defaultdict(lambda: [0, 0, 0, 0])
            for i, k in enumerate(keep):
                if not k:
                    continue
                st, partial, disrupt = strata[i]
                names = [(st, "stratum"), ("all", "beside")]
                if partial:
                    names.append(("partial_only", "beside"))
                if disrupt:
                    names.append(("with_point_disrupt", "beside"))
                for name in names:
                    cell = cells[name]
                    if truth[i]:
                        cell[0] += 1
                        cell[1] += c[i]
                    else:
                        cell[2] += 1
                        cell[3] += c[i]
            for name, kind in [(s, "stratum") for s in STRATA] + [(b, "beside") for b in BESIDE]:
                r, rc, s, sc = cells.get((name, kind), [0, 0, 0, 0])
                out.append({"split": split, "row": row, "stratum": name, "kind": kind,
                            "_r": r, "_rc": rc, "_s": s, "_sc": sc})
    return out


def shown_row(base, row):
    d = dict(base)
    d.update({"split": row["split"], "row": row["row"], "stratum": row["stratum"], "kind": row["kind"],
              "resistant": shown(row["_r"]), "resistant_called": shown(row["_rc"]),
              "susceptible": shown(row["_s"]), "susceptible_called": shown(row["_sc"]),
              "recall": recall_cell(row["_rc"], row["_r"]) if row["_r"] else ""})
    return d


def pooled_rows(entered):
    """Strata under the recall floor pooled per species, drug class and set,
    over genome and dataset pairs."""
    pool = collections.defaultdict(lambda: [set(), 0, 0, 0, 0])
    for base, rows in entered:
        for row in rows:
            if row["kind"] != "stratum" or row["_r"] >= RECALL_FLOOR:
                continue
            key = (base["species_slug"], base["species"], base["drug_class"], base["set"], row["split"], row["row"],
                   row["stratum"])
            p = pool[key]
            p[0].add(base["dataset"])
            p[1] += row["_r"]
            p[2] += row["_rc"]
            p[3] += row["_s"]
            p[4] += row["_sc"]
    out = []
    for key in sorted(pool):
        ds, r, rc, s, sc = pool[key]
        out.append(dict(zip(POOLED_HEADER[:7], key), datasets=len(ds), resistant_pairs=shown(r),
                        resistant_called=shown(rc), susceptible_pairs=shown(s), susceptible_called=shown(sc),
                        recall=recall_cell(rc, r) if r else ""))
    return out


def read_calls(path, wanted):
    """The adapter's tool_default call per (dataset, genome) of `wanted`."""
    out = {}
    with open(path, encoding="utf-8", newline="") as f:
        head = f.readline().rstrip("\n").split("\t")
        ix = {h: i for i, h in enumerate(head)}
        for line in f:
            p = line.rstrip("\n").split("\t")
            if p[ix["variant"]] != "tool_default" or p[ix["dataset_id"]] not in wanted:
                continue
            key = (p[ix["dataset_id"]], p[ix["genome_id"]])
            if key in out:
                ar.refuse(f"{key[0]}: a genome has two tool_default calls in the adapter's file")
            out[key] = p[ix["call"]]
    return out


def run(args):
    lists = [r for r in read_tsv(args.lists) if r["species_slug"] == args.species and r["in_A"] == "1"]
    if not lists:
        ar.refuse(f"{args.species}: no dataset of A in the list")
    rules = {r["dataset"]: r for r in read_tsv(args.drug_rules)}
    record = {(r["dataset"], r["caller"]): r for r in read_tsv(args.callers)}
    wanted = {r["dataset"] for r in lists}
    adapter = read_calls(args.calls, wanted)
    tables = {}
    dataset_rows, strata_rows, nocall_rows, entered = [], [], [], []
    hard_pass = rows_equal = rows_total = 0
    for meta in sorted(lists, key=lambda r: r["dataset"]):
        d = meta["dataset"]
        path = os.path.join(args.reports, d + ".json")
        with open(path, encoding="utf-8") as f:
            report = json.load(f)
        if report.get("id") != d or report.get("species_slug") != args.species:
            ar.refuse(f"{d}: the report file holds another dataset")
        dset = "M" if meta["in_M"] == "1" else "primary"
        base = {"dataset": d, "species_slug": args.species, "species": meta["species"],
                "antibiotic": meta["antibiotic"], "drug_class": report.get("drug_class", ""), "set": dset,
                "labels_depend_on": meta["labels_depend_on"] if dset == "M" else ""}
        why = []
        rule = report["catalog"]["rule"] if report.get("catalog") else None
        listed = rules.get(d)
        rule_ok = rule is not None and listed is not None and ar.parse_rule(listed["rule"]) == rule
        if not rule_ok:
            why.append("catalog.rule differs from drug_rules.tsv")
        ok, eq, tot, hwhy, row_calls, nocall = hard_check(report)
        rows_equal += eq
        rows_total += tot
        hard_pass += ok
        if hwhy:
            why.append(hwhy)
        samples = report["samples"]
        keep, program = [], []
        prog_ok = True
        for s in samples:
            a = adapter.get((d, s["id"]))
            if a is None:
                ar.refuse(f"{d}: the adapter's file has no tool_default call for a genome of the report")
            if a == "no_call":
                keep.append(False)
                program.append(False)
                continue
            if s["id"] not in tables:
                tables[s["id"]] = ar.read_table(os.path.join(args.tables, s["id"] + ".tsv"))
            call = ar.program_call(tables[s["id"]], rule) if rule else False
            keep.append(True)
            program.append(call)
            prog_ok &= call == (a == "resistant")
        if not prog_ok:
            why.append("the program's call differs from the adapter's")
        truth_k = [s["resistant"] for s, k in zip(samples, keep) if k]
        prog_k = [c for c, k in zip(program, keep) if k]
        pc = confusion(truth_k, prog_k)
        rec = record.get((d, "program:tool_default"))
        rec_ok = rec is not None and all(int(rec[k]) == pc[k] for k in ("tp", "fn", "tn", "fp")) \
            and int(rec["scored"]) == len(truth_k)
        if not rec_ok:
            why.append("the program's counts differ from the record's program:tool_default row")
        entered_ok = rule_ok and ok and prog_ok and rec_ok
        tab = {s["id"]: tables.get(s["id"]) if k else None for s, k in zip(samples, keep)}
        dataset_rows.append(dict(base, genomes=shown(len(samples)), with_table=shown(sum(keep)),
                                 rule_check=int(rule_ok), hard_check=int(ok), learner_rows_equal=eq,
                                 learner_rows=tot, program_check=int(prog_ok), record_check=int(rec_ok),
                                 entered=int(entered_ok), why="; ".join(why)))
        for r in nocall:
            nocall_rows.append({"dataset": d, "split": r["split"], "learner": r["learner"],
                                "no_call_resistant": shown(r["no_call_resistant"]),
                                "no_call_susceptible": shown(r["no_call_susceptible"])})
        if not entered_ok:
            continue
        rows = dataset_tables(report, tab, rule, row_calls, program, keep)
        entered.append((base, rows))
        strata_rows.extend(shown_row(base, r) for r in rows)
    os.makedirs(args.out, exist_ok=True)
    write_tsv(os.path.join(args.out, "datasets.tsv"), DATASETS_HEADER, dataset_rows)
    write_tsv(os.path.join(args.out, "strata.tsv"), STRATA_HEADER, strata_rows)
    write_tsv(os.path.join(args.out, "pooled.tsv"), POOLED_HEADER, pooled_rows(entered))
    write_tsv(os.path.join(args.out, "nocall.tsv"), NOCALL_HEADER, nocall_rows)
    summary = {"species": args.species, "datasets": len(lists), "hard_check_pass": hard_pass,
               "learner_rows_equal": rows_equal, "learner_rows": rows_total, "entered": len(entered),
               "strata_rows": len(strata_rows)}
    with open(os.path.join(args.out, "summary.json"), "w", encoding="utf-8") as f:
        f.write(json.dumps(summary, indent=1, sort_keys=True) + "\n")
    print(f"strata: {args.species}: hard check {hard_pass} of {len(lists)} datasets pass "
          f"({rows_equal} of {rows_total} rows equal); {len(entered)} entered; {len(strata_rows)} strata rows")
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser()
    for a in ("--species", "--reports", "--lists", "--drug-rules", "--tables", "--calls", "--callers", "--out"):
        ap.add_argument(a, required=True)
    args = ap.parse_args(argv)
    try:
        return run(args)
    except ar.Refusal as e:
        print(f"strata.py: refused: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
