#!/usr/bin/env python3
"""Build the record of the TBProfiler comparison on the revision 14 M. tuberculosis reports.

    python3 -I scripts/tbprofiler-comparison/record.py --inputs DIR --out reference/tbprofiler-comparison
    python3 -I scripts/tbprofiler-comparison/record.py --inputs DIR --check reference/tbprofiler-comparison
    python3 -I scripts/tbprofiler-comparison/record.py --check reference/tbprofiler-comparison

With --inputs the script reads the call file and the reports' compact file,
counts every caller again on every report and way of reading, redraws the
lineage interval, and writes the row files; it then derives every other file
from the row files. --out writes the record, --check rebuilds it and fails on
any byte that differs. --check alone reads the record's row files, derives the
other files again and fails on any byte that differs from them or from
SHA256SUMS.

The input set DIR holds:

    calls.tsv                     the tool's call per genome, report and way of reading
    rev14_compact.jsonl.gz        the revision 14 reports: genome ids, labels, lineage ids, learner scores
    database-grades-by-drug.tsv   the shipped database's entries per tool drug, by source and grade
    PLAN.md, PLAN-ADDENDUM.md     the plan of 14 Nov and its addendum
    PLAN.sha256                   their hash list, with the time each hash was taken
    finish.log                    the log of the scoring run
    exploratory.json              the counts both ways for the entry the tool's rule set aside

The record holds no genome id, label, call or score per genome. Its rows are
per report and caller, per way of reading, per tool drug, per grade, per
learner, per event or a single count. The script reads nothing else and opens
no network connection.
"""
import argparse
import gzip
import hashlib
import json
import math
import random
import re
import sys
from collections import Counter
from pathlib import Path

SPECIES = "mycobacterium_tuberculosis"
WAYS = ("who_assoc", "who_assoc_only", "tool_default", "who_assoc_small")
PRIMARY_WAY = "who_assoc"
LEARNERS = ("scm", "cart", "rf", "gbt", "lasso")
SPLITS = ("lineage", "random")
LEARNER_CALLERS = tuple(f"{l}:{s}" for s in SPLITS for l in LEARNERS)
BOOT_CALLERS = ("tool",) + LEARNER_CALLERS
SETS = ("by_name", "graded")
MEASURES = ("sensitivity", "specificity", "balanced_accuracy")
BOUNDS = ("sensitivity_lower", "sensitivity_upper", "specificity_lower", "specificity_upper")
CALLS = {"resistant": True, "susceptible": False, "no_call": None}
GRADES = {"Assoc w R": "assoc", "Assoc w R - Interim": "interim"}
MIN_CLASS = 25
CUTOFF = 0.5
MARGIN = 0.02
RESAMPLES = 2000
SEED = 20261114
LOW_RANK, HIGH_RANK = 51, 1950  # the 51st and the 1,950th of 2,000 sorted values
MIN_VALUES = 40
CATALOG_COUNTS = ("catalog_endpoints", "catalog_rows", "catalog_tool_only_right", "catalog_only_right")
TOOL_VERSION = "6.6.3"
DATABASE_COMMIT = "3ed1cb99"
FIELD = 32  # bits per count in the packed bootstrap vectors

ROW_FILES = ("reports.tsv", "callers.tsv", "intervals.tsv", "disagreements.tsv", "grades.tsv", "order.tsv",
             "settings.tsv", "totals.tsv", "set-aside.tsv", "checks.tsv", "inputs.sha256")
DERIVED_FILES = ("measures.tsv", "endpoints.tsv", "medians.tsv", "comparisons.tsv", "screen.tsv", "bounds.tsv",
                 "summary.tsv")
FILES = ROW_FILES + DERIVED_FILES + ("README.md", "SHA256SUMS")

COUNT_COLS = ["genomes", "with_call", "no_call_resistant", "no_call_susceptible", "tp", "fp", "tn", "fn"]
DIS_COLS = ["tool_r_label_s_assoc", "tool_r_label_s_interim", "tool_s_label_r_other", "tool_s_label_r_none",
            "both_r_assoc", "both_r_interim", "both_s_other", "both_s_none", "no_call_label_r", "no_call_label_s"]
GRADE_COLS = ["assoc_w_r", "assoc_w_r_interim", "own_uncertain_significance", "own_not_assoc_w_r_interim",
              "own_ungraded", "counted_grade_entries"]
COLUMNS = {
    "reports.tsv": ["dataset", "antibiotic", "kind", "tool_drugs", "graded", "genomes", "catalog_row"],
    "callers.tsv": ["way", "dataset", "caller"] + COUNT_COLS + ["caller_only_right", "tool_only_right"],
    "intervals.tsv": ["kind", "set", "dataset", "caller", "measure", "low", "high", "resamples_with_value"],
    "disagreements.tsv": ["dataset"] + DIS_COLS,
    "grades.tsv": ["tool_drug"] + GRADE_COLS,
    "order.tsv": ["event", "time", "source"],
    "settings.tsv": ["name", "value"],
    "totals.tsv": ["name", "value"],
    "set-aside.tsv": ["name", "value"],
    "checks.tsv": ["check", "count"],
    "measures.tsv": ["way", "dataset", "caller", "sensitivity", "specificity", "balanced_accuracy",
                     "minus_tool"],
    "endpoints.tsv": ["way", "dataset", "kind", "scored", "by_name", "graded", "tool_balanced_accuracy",
                      "best_lineage_learner", "best_lineage_balanced_accuracy", "gap", "within_margin"]
    + [f"{b}" for b in BOUNDS],
    "medians.tsv": ["way", "set", "caller", "measure", "endpoints", "median", "min", "max", "low", "high"],
    "comparisons.tsv": ["way", "set", "caller", "endpoints", "tool_higher", "tool_lower", "equal",
                        "tool_ahead_or_within", "median_minus_tool", "low", "high", "caller_only_right",
                        "tool_only_right"],
    "screen.tsv": ["way", "set", "endpoints", "within_best", "needed"],
    "bounds.tsv": ["way", "set", "bound", "endpoints", "median", "min", "max"],
    "summary.tsv": ["name", "value"],
}
INTS = set(COUNT_COLS) | set(DIS_COLS) | set(GRADE_COLS) | {
    "caller_only_right", "tool_only_right", "resamples_with_value", "count", "endpoints", "tool_higher",
    "tool_lower", "equal", "tool_ahead_or_within", "within_best", "needed"}
FLOATS = {"low", "high"}
FIXED_SETTINGS = {"tool_version": TOOL_VERSION, "database_commit": DATABASE_COMMIT, "min_class": str(MIN_CLASS),
                  "cutoff": repr(CUTOFF), "margin": repr(MARGIN), "interval_level": "0.95",
                  "resamples": str(RESAMPLES), "seed": str(SEED), "low_rank": str(LOW_RANK),
                  "high_rank": str(HIGH_RANK), "min_resamples_with_value": str(MIN_VALUES)}
ORDER_EVENTS = ("plan_hashed", "addendum_hashed", "hashes_checked", "last_step_before_scoring", "scoring_done")
# Shapes of what must never enter the record: genome ids of this input set's shape, accessions, base runs.
BANNED_SHAPES = (r"(?<![\w.])\d{4,7}\.\d{1,5}(?![\w.])", r"\bPRJ[A-Z]{1,2}\d+", r"\bSAM[NED][A-Z]?\d+",
                 r"\bGC[AF]_\d+", r"\b[SE]R[RSX]\d{5,}", r"(?i)[acgt]{20,}")


class Stop(Exception):
    pass


def need(cond, msg):
    if not cond:
        raise Stop(msg)


def sha256_bytes(b):
    return hashlib.sha256(b).hexdigest()


def fmt(v):
    """A cell: ints as written, floats as Python's shortest round-trip form, None as '-'."""
    if v is None:
        return "-"
    if isinstance(v, bool):
        return "yes" if v else "no"
    if isinstance(v, float):
        need(math.isfinite(v), "a cell is not a finite number")
        return repr(v)
    return str(v)


def tsv(name, rows):
    cols = COLUMNS[name]
    out = ["\t".join(cols)]
    for r in rows:
        need(set(r) == set(cols), f"{name}: a row has other cells than its columns")
        cells = [fmt(r[c]) for c in cols]
        for c in cells:
            need(c != "" and "\t" not in c and "\n" not in c, f"{name}: a cell is empty or holds a tab or a line break")
        out.append("\t".join(cells))
    return "\n".join(out) + "\n"


def parse_tsv(name, text):
    """A file of the record as dicts, held to its columns; counts and interval ends parsed."""
    lines = text.split("\n")
    need(len(lines) >= 2 and lines[-1] == "", f"{name}: no final newline")
    head = lines[0].split("\t")
    need(head == COLUMNS[name], f"{name}: the columns are not {COLUMNS[name]}")
    rows = []
    for i, line in enumerate(lines[1:-1], 2):
        cells = line.split("\t")
        need(len(cells) == len(head), f"{name}:{i}: {len(cells)} cells, the header has {len(head)}")
        row = {}
        for c, v in zip(head, cells):
            if c in INTS and name not in ("settings.tsv", "totals.tsv", "set-aside.tsv", "summary.tsv"):
                need(re.fullmatch(r"0|[1-9]\d*", v) is not None, f"{name}:{i}: {c} is not a count")
                row[c] = int(v)
            elif c in FLOATS and v == "-":
                row[c] = None
            elif c in FLOATS:
                try:
                    row[c] = float(v)
                except ValueError:
                    raise Stop(f"{name}:{i}: {c} is not a number")
                need(fmt(row[c]) == v, f"{name}:{i}: {c} is not in its shortest form")
            else:
                row[c] = v
        rows.append(row)
    return rows


def pairs(rows):
    return {r["name"]: r["value"] for r in rows}


# ----------------------------------------------------------------- measures

def ratio(a, b):
    return 0.0 if b == 0 else a / b


def rates(tp, fp, tn, fn, strict=False):
    """(sensitivity, specificity, balanced accuracy). strict: None when a class has no genome (a resample);
    otherwise a ratio over zero is 0.0, the reports' rule."""
    if strict and (tp + fn == 0 or tn + fp == 0):
        return None, None, None
    sens, spec = ratio(tp, tp + fn), ratio(tn, tn + fp)
    return sens, spec, (sens + spec) / 2.0


def median(values):
    v = sorted(values)
    need(v, "a median over no value")
    mid = len(v) // 2
    return v[mid] if len(v) % 2 else (v[mid - 1] + v[mid]) / 2.0


def interval(values):
    """The two ends: of the sorted values with a value, the one at int(0.025 n) and the one as far from the top."""
    v = sorted(x for x in values if x is not None)
    if len(v) < MIN_VALUES:
        return None, None, len(v)
    k = int(0.025 * len(v))
    return v[k], v[len(v) - 1 - k], len(v)


def confusion(truth, calls):
    tp = fp = tn = fn = 0
    for t, c in zip(truth, calls):
        if t:
            tp += c
            fn += not c
        else:
            fp += c
            tn += not c
    return tp, fp, tn, fn


# ----------------------------------------------------------------- inputs

def read_calls(path):
    """{dataset: {"kind", "tool_drugs", "rows": {genome: row}}} and the row count."""
    text = path.read_text(encoding="utf-8")
    lines = text.split("\n")
    need(lines[-1] == "", "calls.tsv: no final newline")
    head = lines[0].split("\t")
    want = ["dataset_id", "genome_id", "kind", "tool_drugs"] + list(WAYS) + ["top_grade", "other_grade",
                                                                           "targets_covered", "large_only"]
    need(head == want, f"calls.tsv: the columns are not {want}")
    out = {}
    for n, line in enumerate(lines[1:-1], 2):
        row = dict(zip(head, line.split("\t")))
        need(len(row) == len(head) and line.count("\t") == len(head) - 1, f"calls.tsv:{n}: not {len(head)} cells")
        for w in WAYS:
            need(row[w] in CALLS, f"calls.tsv:{n}: a call that is not resistant, susceptible or no_call")
        need(row["kind"] in ("primary", "proxy"), f"calls.tsv:{n}: a kind that is not primary or proxy")
        need(row["top_grade"] in ("-",) + tuple(GRADES), f"calls.tsv:{n}: a grade the record does not count")
        need(row["other_grade"] in ("yes", "no"), f"calls.tsv:{n}: other_grade is not yes or no")
        need(row["dataset_id"].startswith(SPECIES + "__"), f"calls.tsv:{n}: a report of another species")
        e = out.setdefault(row["dataset_id"], {"kind": row["kind"], "tool_drugs": row["tool_drugs"], "rows": {}})
        need((e["kind"], e["tool_drugs"]) == (row["kind"], row["tool_drugs"]),
             f"calls.tsv:{n}: the kind or the tool drugs change inside a report")
        need(row["genome_id"] not in e["rows"], f"calls.tsv:{n}: a genome twice in one report")
        e["rows"][row["genome_id"]] = row
    need(out, "calls.tsv: no call rows")
    return out, len(lines) - 2


def read_compact(path):
    """{dataset: report} for the species' reports in the compact file."""
    out = {}
    with gzip.open(path, "rt", encoding="utf-8") as f:
        for line in f:
            if not line.strip():
                continue
            d = json.loads(line)
            if d.get("species_slug") != SPECIES:
                continue
            need(d["id"] not in out, f"compact file: {d['id']} twice")
            n = len(d["sid"])
            need(len(set(d["sid"])) == n, f"compact file: {d['id']} lists a genome twice")
            need(len(d["y"]) == n and set(d["y"]) <= {0, 1}, f"compact file: {d['id']}: labels are not one 0 or 1 per genome")
            need(len(d["lineage"]) == n, f"compact file: {d['id']}: not one lineage id per genome")
            for s in SPLITS:
                for l in LEARNERS:
                    sc = d["splits"][s]["scores"][l]
                    need(len(sc) == n and all(isinstance(x, (int, float)) and not isinstance(x, bool) and 0 <= x <= 1
                                              for x in sc), f"compact file: {d['id']}: the {l} {s} scores")
            cat = d.get("catalog_resistant")
            need(isinstance(cat, list) and len(cat) == n, f"compact file: {d['id']}: catalog_resistant")
            need(all(isinstance(x, bool) for x in cat) or all(x is None for x in cat),
                 f"compact file: {d['id']}: catalog_resistant mixes calls and gaps")
            out[d["id"]] = d
    need(out, "compact file: no report of the species")
    return out


def read_grades(path):
    lines = path.read_text(encoding="utf-8").split("\n")
    want = ["tool_drug", "WHO catalogue v2: Assoc w R", "WHO catalogue v2: Assoc w R - Interim",
            "tbdb: Uncertain significance", "tbdb: Not assoc w R - Interim", "tbdb: none", "counted_grade_entries"]
    need(lines[0].split("\t") == want and lines[-1] == "", "database-grades-by-drug.tsv: not the grade table")
    rows = []
    for line in lines[1:-1]:
        cells = line.split("\t")
        need(len(cells) == len(want) and all(re.fullmatch(r"0|[1-9]\d*", c) for c in cells[1:]),
             "database-grades-by-drug.tsv: a row that is not a drug and six counts")
        r = {"tool_drug": cells[0]}
        r.update({c: int(v) for c, v in zip(GRADE_COLS, cells[1:])})
        need(r["counted_grade_entries"] == r["assoc_w_r"] + r["assoc_w_r_interim"],
             f"database-grades-by-drug.tsv: {cells[0]}: the counted entries are not the two counted grades")
        rows.append(r)
    need(len({r["tool_drug"] for r in rows}) == len(rows), "database-grades-by-drug.tsv: a drug twice")
    return rows


def read_plans(I):
    """The plan's and the addendum's sha256 and the time the hash list gives each."""
    listed, times, last = {}, {}, None
    for line in (I / "PLAN.sha256").read_text(encoding="utf-8").split("\n"):
        m = re.fullmatch(r"([0-9a-f]{64})  (\S+)", line)
        if m:
            listed[m.group(2)] = m.group(1)
            last = m.group(2)
            continue
        t = re.search(r"\b(\d{4}-\d\d-\d\dT\d\d:\d\d:\d\dZ)\b", line)
        if line.startswith("#") and t and last and last not in times:
            times[last] = t.group(1)
    for name, what in (("PLAN.md", "plan"), ("PLAN-ADDENDUM.md", "addendum")):
        need(name in listed and name in times, f"PLAN.sha256: no hash and time for the {what}")
        need(sha256_bytes((I / name).read_bytes()) == listed[name], f"not the hashed {what}")
    need(times["PLAN.md"] < times["PLAN-ADDENDUM.md"], "the addendum was hashed before the plan")
    plan = (I / "PLAN.md").read_text(encoding="utf-8")
    need(f"TBProfiler {TOOL_VERSION}" in plan and f"`{DATABASE_COMMIT}`" in plan,
         "the plan names another tool version or database commit")
    return listed, times


def read_log(path):
    """The times of the finish log's steps."""
    found = {}
    for line in path.read_text(encoding="utf-8").split("\n"):
        m = re.fullmatch(r"(\d{4}-\d\d-\d\dT\d\d:\d\d:\d\dZ) (.+)", line)
        if m:
            found.setdefault(m.group(2), m.group(1))
    out = {}
    for key, words in (("hashes_checked", "plan and addendum match PLAN.sha256"),
                       ("last_step_before_scoring", "input limit done"), ("scoring_done", "evaluator done")):
        need(words in found, f"finish.log: no line '{words}'")
        out[key] = found[words]
    return out


def read_set_aside(path):
    body = json.loads(path.read_text(encoding="utf-8"))
    need(str(body.get("label", "")).startswith("EXPLORATORY"), "exploratory.json: not labelled exploratory")
    rule = body["rule_that_sets_an_entry_aside"]
    touched = rule["primary_endpoints_touched"]
    need(len(touched) == 1, "exploratory.json: the rule touched other than one report")
    (dataset, t), = touched.items()
    out = [("result_files_with_a_note", rule["result_files_with_a_note"]),
           ("entries_set_aside", sum(rule["entries_set_aside_by_drug_and_grade"].values())),
           ("dataset", dataset.split("__", 1)[1])]
    for name, key in (("as_run", "counts_as_run"), ("counted", "counts_with_the_entry_counted")):
        c = t[key]
        need(len(c) == 4 and all(isinstance(x, int) for x in c), f"exploratory.json: {key} is not four counts")
        out += [(f"{k}_{name}", v) for k, v in zip(("tp", "fp", "tn", "fn"), c)]
    return out


def lineage_groups(reports):
    """{genome: group number}: lineage ids joined across the reports, groups numbered in the text order of their
    smallest genome id."""
    parent = {}

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for d in reports:
        first = {}
        for gid, lin in zip(d["sid"], d["lineage"]):
            parent.setdefault(gid, gid)
            if lin in first:
                ra, rb = find(first[lin]), find(gid)
                if ra != rb:
                    parent[max(ra, rb)] = min(ra, rb)
            else:
                first[lin] = gid
    roots = sorted({find(a) for a in parent})
    index = {r: i for i, r in enumerate(roots)}
    return {a: index[find(a)] for a in parent}, len(roots)


# ----------------------------------------------------------------- build

def callers_of(d):
    """{caller: calls on every genome of the report}: each learner's scores at the cutoff. The catalog row is read
    apart, by catalog_figures, and has no row per report."""
    out = {}
    for s in SPLITS:
        for l in LEARNERS:
            out[f"{l}:{s}"] = [x >= CUTOFF for x in d["splits"][s]["scores"][l]]
    return out


def catalog_figures(reports, calls, members):
    """The catalog row over the reports of the plan's set that have one, under the primary way, on the genomes with a
    call: reports, rows, the rows only the tool and only the catalog row get right, and the median of the tool's
    balanced accuracy minus the catalog row's, to 3 decimals. On three reports a count of the catalog row per report
    would describe one genome, and the median in full would show whether one of them had a call."""
    diffs, rows, only_t, only_c = [], 0, 0, 0
    for dataset in sorted(members):
        d = reports[dataset]
        cat = d["catalog_resistant"]
        if not cat or cat[0] is None:
            continue
        entry = calls[dataset]
        truth, tool, catalog = [], [], []
        for g, y, c in zip(d["sid"], d["y"], cat):
            t = CALLS[entry["rows"][g][PRIMARY_WAY]]
            if t is not None:
                truth.append(bool(y))
                tool.append(t)
                catalog.append(c)
        rows += len(truth)
        only_t += sum(1 for y, a, b in zip(truth, tool, catalog) if a == y and b != y)
        only_c += sum(1 for y, a, b in zip(truth, tool, catalog) if b == y and a != y)
        diffs.append(rates(*confusion(truth, tool))[2] - rates(*confusion(truth, catalog))[2])
    med = "%.3f" % median(diffs) if diffs else "-"
    return [("catalog_endpoints", len(diffs)), ("catalog_rows", rows), ("catalog_tool_only_right", only_t),
            ("catalog_only_right", only_c), ("tool_minus_catalog_median", med)]


def score(dataset, d, entry, way, others):
    """The rows of callers.tsv for one report and way, and (for the primary way) its packed group vectors."""
    ids = d["sid"]
    need(set(entry["rows"]) == set(ids), f"{dataset}: the call file's genomes are not the report's")
    truth_all = [bool(y) for y in d["y"]]
    tool_all = [CALLS[entry["rows"][g][way]] for g in ids]
    keep = [i for i, c in enumerate(tool_all) if c is not None]
    truth = [truth_all[i] for i in keep]
    tool = [tool_all[i] for i in keep]
    base = {"genomes": len(ids), "with_call": len(keep),
            "no_call_resistant": sum(1 for i, c in enumerate(tool_all) if c is None and truth_all[i]),
            "no_call_susceptible": sum(1 for i, c in enumerate(tool_all) if c is None and not truth_all[i])}
    rows = []
    callers = {"tool": tool}
    for name in LEARNER_CALLERS:
        if name in others:
            callers[name] = [others[name][i] for i in keep]
    for name, calls in callers.items():
        tp, fp, tn, fn = confusion(truth, calls)
        only_c = sum(1 for t, a, b in zip(truth, calls, tool) if a == t and b != t) if name != "tool" else 0
        only_t = sum(1 for t, a, b in zip(truth, calls, tool) if a != t and b == t) if name != "tool" else 0
        rows.append(dict(base, way=way, dataset=dataset, caller=name, tp=tp, fp=fp, tn=tn, fn=fn,
                         caller_only_right=only_c, tool_only_right=only_t))
    return rows, keep, truth, callers


def pack(group_of, ids, keep, truth, callers):
    """{group: one integer holding tp, fp, tn, fn of every bootstrap caller on the group's genomes}."""
    table = {}
    for pos, i in enumerate(keep):
        g = group_of[ids[i]]
        t = truth[pos]
        v = 0
        for j, name in enumerate(BOOT_CALLERS):
            p = callers[name][pos]
            field = (0 if p else 3) if t else (1 if p else 2)
            v += 1 << (FIELD * (4 * j + field))
        table[g] = table.get(g, 0) + v
    return table


def bootstrap(tables, n_groups, sets):
    """The lineage interval's rows: per report and caller, per learner minus the tool, and the medians over each
    set. One generator seeded 20261114 draws every resample, one randrange call per group drawn."""
    rng = random.Random(SEED)
    mask = (1 << FIELD) - 1
    width = 4 * len(BOOT_CALLERS)
    datasets = list(tables)
    per = {d: [] for d in datasets}  # per resample: [(sens, spec, ba) per caller]
    for _ in range(RESAMPLES):
        drawn = Counter(rng.randrange(n_groups) for _ in range(n_groups))
        for d in datasets:
            t = tables[d]
            total = sum(m * t[g] for g, m in drawn.items() if g in t)
            f = [(total >> (FIELD * k)) & mask for k in range(width)]
            per[d].append([rates(*f[4 * j:4 * j + 4], strict=True) for j in range(len(BOOT_CALLERS))])
    rows = []
    for d in datasets:
        for j, name in enumerate(BOOT_CALLERS):
            for mi, m in enumerate(MEASURES):
                lo, hi, n = interval(r[j][mi] for r in per[d])
                rows.append({"kind": "endpoint", "set": "-", "dataset": d, "caller": name, "measure": m,
                             "low": lo, "high": hi, "resamples_with_value": n})
        for j, name in enumerate(BOOT_CALLERS[1:], 1):
            diffs = [None if r[j][2] is None or r[0][2] is None else r[j][2] - r[0][2] for r in per[d]]
            lo, hi, n = interval(diffs)
            rows.append({"kind": "difference", "set": "-", "dataset": d, "caller": name, "measure": "minus_tool",
                         "low": lo, "high": hi, "resamples_with_value": n})
    for s in SETS:
        members = [d for d in datasets if d in sets[s]]
        for j, name in enumerate(BOOT_CALLERS):
            for mi, m in enumerate(MEASURES):
                meds = []
                for b in range(RESAMPLES):
                    v = [per[d][b][j][mi] for d in members if per[d][b][j][mi] is not None]
                    meds.append(median(v) if v else None)
                lo, hi, n = interval(meds)
                rows.append({"kind": "median", "set": s, "dataset": "-", "caller": name, "measure": m,
                             "low": lo, "high": hi, "resamples_with_value": n})
        for j, name in enumerate(BOOT_CALLERS[1:], 1):
            meds = []
            for b in range(RESAMPLES):
                v = [per[d][b][j][2] - per[d][b][0][2] for d in members
                     if per[d][b][j][2] is not None and per[d][b][0][2] is not None]
                meds.append(median(v) if v else None)
            lo, hi, n = interval(meds)
            rows.append({"kind": "median_difference", "set": s, "dataset": "-", "caller": name,
                         "measure": "minus_tool", "low": lo, "high": hi, "resamples_with_value": n})
    return rows


def disagreement(dataset, d, entry):
    ids, truth = d["sid"], [bool(y) for y in d["y"]]
    c = dict.fromkeys(DIS_COLS, 0)
    for g, t in zip(ids, truth):
        row = entry["rows"][g]
        call = CALLS[row[PRIMARY_WAY]]
        if call is None:
            c["no_call_label_r" if t else "no_call_label_s"] += 1
        elif call:
            need(row["top_grade"] in GRADES, f"{dataset}: a resistant call with no counted grade")
            c[("both_r_" if t else "tool_r_label_s_") + GRADES[row["top_grade"]]] += 1
        else:
            kind = "other" if row["other_grade"] == "yes" else "none"
            c[("tool_s_label_r_" if t else "both_s_") + kind] += 1
    return dict(c, dataset=dataset)


def build(I):
    I = Path(I)
    names = ("calls.tsv", "rev14_compact.jsonl.gz", "database-grades-by-drug.tsv", "PLAN.md", "PLAN-ADDENDUM.md",
             "PLAN.sha256", "finish.log", "exploratory.json")
    for n in names:
        need((I / n).is_file(), f"{n}: not in the input set")
    listed, plan_times = read_plans(I)
    log = read_log(I / "finish.log")
    need(plan_times["PLAN-ADDENDUM.md"] < log["hashes_checked"] <= log["last_step_before_scoring"]
         <= log["scoring_done"], "a plan was hashed after scoring began, or the log's steps are out of order")
    calls, call_rows = read_calls(I / "calls.tsv")
    reports = read_compact(I / "rev14_compact.jsonl.gz")
    grades = read_grades(I / "database-grades-by-drug.tsv")
    graded_drugs = {r["tool_drug"] for r in grades if r["counted_grade_entries"] > 0}
    need(set(calls) <= set(reports), "calls.tsv names a report the compact file does not hold")
    group_of, n_groups = lineage_groups([reports[k] for k in sorted(reports)])
    report_rows, caller_rows, dis_rows, tables, scored = [], [], [], {}, {}
    joined = 0
    for dataset in sorted(reports):
        d = reports[dataset]
        e = calls.get(dataset)
        has_cat = d["catalog_resistant"][0] is not None if d["catalog_resistant"] else False
        if e is None:
            report_rows.append({"dataset": dataset, "antibiotic": d["antibiotic"], "kind": "not named",
                                "tool_drugs": "-", "graded": "-", "genomes": len(d["sid"]),
                                "catalog_row": "yes" if has_cat else "no"})
            continue
        drugs = e["tool_drugs"].split(",")
        need(all(x in {r["tool_drug"] for r in grades} for x in drugs), f"{dataset}: a tool drug the grade table lacks")
        need(e["kind"] == "proxy" or len(drugs) == 1, f"{dataset}: a primary report with more than one tool drug")
        graded = ("yes" if drugs[0] in graded_drugs else "no") if e["kind"] == "primary" else "-"
        report_rows.append({"dataset": dataset, "antibiotic": d["antibiotic"], "kind": e["kind"],
                            "tool_drugs": e["tool_drugs"], "graded": graded, "genomes": len(d["sid"]),
                            "catalog_row": "yes" if has_cat else "no"})
        others = callers_of(d)
        for way in WAYS:
            rows, keep, truth, callers = score(dataset, d, e, way, others)
            caller_rows += rows
            if way == PRIMARY_WAY:
                tables[dataset] = pack(group_of, d["sid"], keep, truth, callers)
                tool = rows[0]
                scored[dataset] = (e["kind"] == "primary" and tool["tp"] + tool["fn"] >= MIN_CLASS
                                   and tool["tn"] + tool["fp"] >= MIN_CLASS, graded == "yes")
        dis_rows.append(disagreement(dataset, d, e))
        joined += 1
    sets = {"by_name": {k for k, (s, g) in scored.items() if s}, "graded": {k for k, (s, g) in scored.items() if s and g}}
    interval_rows = bootstrap(tables, n_groups, sets)
    ids = {g for d in reports.values() for g in d["sid"]}
    call_ids = {g for e in calls.values() for g in e["rows"]}
    need(call_ids == ids, "the call file's genomes are not the reports' genomes")
    order = [{"event": "plan_hashed", "time": plan_times["PLAN.md"], "source": "the plan's hash list"},
             {"event": "addendum_hashed", "time": plan_times["PLAN-ADDENDUM.md"], "source": "the plan's hash list"}]
    order += [{"event": k, "time": log[k], "source": "the log of the scoring run"} for k in ORDER_EVENTS[2:]]
    settings = dict(FIXED_SETTINGS, plan_sha256=listed["PLAN.md"], addendum_sha256=listed["PLAN-ADDENDUM.md"],
                    compact_sha256=sha256_bytes((I / "rev14_compact.jsonl.gz").read_bytes()))
    kinds = Counter(r["kind"] for r in report_rows)
    totals = [("reports_of_the_species", len(reports)), ("reports_with_calls", len(calls)),
              ("primary_reports", kinds["primary"]), ("proxy_reports", kinds["proxy"]),
              ("not_named_reports", kinds["not named"]), ("call_rows", call_rows), ("distinct_genomes", len(ids)),
              ("lineage_groups", n_groups), ("tool_drugs", len(grades)), ("graded_tool_drugs", len(graded_drugs))]
    totals += catalog_figures(reports, calls, sets["by_name"])
    rows_text = {
        "reports.tsv": tsv("reports.tsv", report_rows),
        "callers.tsv": tsv("callers.tsv", sorted(caller_rows, key=lambda r: (WAYS.index(r["way"]), r["dataset"]))),
        "intervals.tsv": tsv("intervals.tsv", interval_rows),
        "disagreements.tsv": tsv("disagreements.tsv", dis_rows),
        "grades.tsv": tsv("grades.tsv", grades),
        "order.tsv": tsv("order.tsv", order),
        "settings.tsv": tsv("settings.tsv", [{"name": k, "value": v} for k, v in settings.items()]),
        "totals.tsv": tsv("totals.tsv", [{"name": k, "value": v} for k, v in totals]),
        "set-aside.tsv": tsv("set-aside.tsv", [{"name": k, "value": v} for k, v in read_set_aside(I / "exploratory.json")]),
    }
    inputs = [I / n for n in names]
    rows_text["inputs.sha256"] = "".join(f"{sha256_bytes(p.read_bytes())}  {p.name}\n" for p in inputs)
    checks = [("reports_whose_genomes_equal_the_call_file", joined), ("plans_match_their_hash", 2),
              ("call_rows_read", call_rows)]
    return finish(rows_text, checks, ids)


# ----------------------------------------------------------------- derive

def check_rows(rows):
    """The row files agree with one another and hold one row per report, caller and way."""
    reps = rows["reports.tsv"]
    need([r["dataset"] for r in reps] == sorted(r["dataset"] for r in reps), "reports.tsv: not in dataset order")
    need(len({r["dataset"] for r in reps}) == len(reps), "reports.tsv: a report twice")
    for r in reps:
        need(r["kind"] in ("primary", "proxy", "not named"), f"reports.tsv: {r['dataset']}: kind")
        need(r["graded"] in (("yes", "no") if r["kind"] == "primary" else ("-",)), f"reports.tsv: {r['dataset']}: graded")
        need(r["genomes"] > 0, f"reports.tsv: {r['dataset']}: genomes")
    named = [r for r in reps if r["kind"] != "not named"]
    want = []
    for w in WAYS:
        for r in named:
            want += [(w, r["dataset"], c) for c in ("tool",) + LEARNER_CALLERS]
    have = [(r["way"], r["dataset"], r["caller"]) for r in rows["callers.tsv"]]
    need(have == want, "callers.tsv: not one row per way, report and caller")
    genomes = {r["dataset"]: int(r["genomes"]) for r in reps}
    first = {}
    for r in rows["callers.tsv"]:
        need(r["with_call"] + r["no_call_resistant"] + r["no_call_susceptible"] == r["genomes"] == genomes[r["dataset"]],
             f"callers.tsv: {r['way']} {r['dataset']} {r['caller']}: the genomes do not add up")
        need(r["tp"] + r["fp"] + r["tn"] + r["fn"] == r["with_call"],
             f"callers.tsv: {r['way']} {r['dataset']} {r['caller']}: the confusion counts do not add up")
        k = (r["way"], r["dataset"])
        if r["caller"] == "tool":
            first[k] = r
            need(r["caller_only_right"] == r["tool_only_right"] == 0, "callers.tsv: the tool's own row counts pairs")
        else:
            t = first[k]
            need(r["tp"] + r["fn"] == t["tp"] + t["fn"] and r["no_call_resistant"] == t["no_call_resistant"]
                 and r["no_call_susceptible"] == t["no_call_susceptible"],
                 f"callers.tsv: {r['way']} {r['dataset']} {r['caller']}: other genomes than the tool's")
            right_c, right_t = r["tp"] + r["tn"], t["tp"] + t["tn"]
            need(right_c - right_t == r["caller_only_right"] - r["tool_only_right"]
                 and r["caller_only_right"] + r["tool_only_right"] <= r["with_call"],
                 f"callers.tsv: {r['way']} {r['dataset']} {r['caller']}: the pair counts do not follow")
    dis = rows["disagreements.tsv"]
    need([r["dataset"] for r in dis] == [r["dataset"] for r in named], "disagreements.tsv: not one row per report")
    for r in dis:
        t = first[(PRIMARY_WAY, r["dataset"])]
        need((r["both_r_assoc"] + r["both_r_interim"], r["tool_r_label_s_assoc"] + r["tool_r_label_s_interim"],
              r["both_s_other"] + r["both_s_none"], r["tool_s_label_r_other"] + r["tool_s_label_r_none"],
              r["no_call_label_r"], r["no_call_label_s"])
             == (t["tp"], t["fp"], t["tn"], t["fn"], t["no_call_resistant"], t["no_call_susceptible"]),
             f"disagreements.tsv: {r['dataset']}: other counts than callers.tsv")
    drugs = {r["tool_drug"] for r in rows["grades.tsv"]}
    for r in named:
        need(all(x in drugs for x in r["tool_drugs"].split(",")), f"reports.tsv: {r['dataset']}: a drug grades.tsv lacks")
    graded = {r["tool_drug"] for r in rows["grades.tsv"] if r["counted_grade_entries"] > 0}
    for r in rows["grades.tsv"]:
        need(r["counted_grade_entries"] == r["assoc_w_r"] + r["assoc_w_r_interim"], "grades.tsv: counted entries")
    for r in named:
        if r["kind"] == "primary":
            need((r["graded"] == "yes") == (r["tool_drugs"] in graded), f"reports.tsv: {r['dataset']}: graded")
    ev = rows["order.tsv"]
    need([r["event"] for r in ev] == list(ORDER_EVENTS), "order.tsv: not the five events")
    need([r["time"] for r in ev] == sorted(r["time"] for r in ev), "order.tsv: the times are out of order")
    st = pairs(rows["settings.tsv"])
    for k, v in FIXED_SETTINGS.items():
        need(st.get(k) == v, f"settings.tsv: {k} is not {v}")
    tot = pairs(rows["totals.tsv"])
    need(int(tot["primary_reports"]) + int(tot["proxy_reports"]) + int(tot["not_named_reports"])
         == int(tot["reports_of_the_species"]) == len(reps), "totals.tsv: the reports do not add up")
    need(int(tot["reports_with_calls"]) == len(named), "totals.tsv: reports with calls")
    need(int(tot["call_rows"]) == sum(int(r["genomes"]) for r in named), "totals.tsv: the call rows are not the genomes")
    for k in CATALOG_COUNTS:
        need(re.fullmatch(r"0|[1-9]\d*", tot.get(k, "")) is not None, f"totals.tsv: {k} is not a count")
    need(re.fullmatch(r"-|-?\d\.\d{3}", tot.get("tool_minus_catalog_median", "")) is not None,
         "totals.tsv: the catalog row's median is not given to 3 decimals")
    need(int(tot["catalog_tool_only_right"]) + int(tot["catalog_only_right"]) <= int(tot["catalog_rows"]),
         "totals.tsv: the catalog row's pair counts exceed its rows")
    sa = pairs(rows["set-aside.tsv"])
    t = first[(PRIMARY_WAY, SPECIES + "__" + sa["dataset"])]
    need([int(sa[f"{k}_as_run"]) for k in ("tp", "fp", "tn", "fn")] == [t["tp"], t["fp"], t["tn"], t["fn"]],
         "set-aside.tsv: the counts as run are not callers.tsv's")
    need(sum(int(sa[f"{k}_counted"]) for k in ("tp", "fp", "tn", "fn")) == t["with_call"],
         "set-aside.tsv: the counts with the entry counted do not add up")
    iv = rows["intervals.tsv"]
    want = []
    for r in named:
        want += [("endpoint", "-", r["dataset"], c, m) for c in BOOT_CALLERS for m in MEASURES]
        want += [("difference", "-", r["dataset"], c, "minus_tool") for c in LEARNER_CALLERS]
    for s in SETS:
        want += [("median", s, "-", c, m) for c in BOOT_CALLERS for m in MEASURES]
        want += [("median_difference", s, "-", c, "minus_tool") for c in LEARNER_CALLERS]
    need([(r["kind"], r["set"], r["dataset"], r["caller"], r["measure"]) for r in iv] == want,
         "intervals.tsv: not one row per report, caller and measure and per set")
    for r in iv:
        need(r["resamples_with_value"] <= int(st["resamples"]), "intervals.tsv: more resamples than were drawn")
        need((r["low"] is None) == (r["high"] is None) == (r["resamples_with_value"] < MIN_VALUES)
             and (r["low"] is None or r["low"] <= r["high"]), "intervals.tsv: the ends do not follow")
    return first


def derive(rows):
    first = check_rows(rows)
    reps = {r["dataset"]: r for r in rows["reports.tsv"]}
    named = [r["dataset"] for r in rows["reports.tsv"] if r["kind"] != "not named"]
    by = {(r["way"], r["dataset"], r["caller"]): r for r in rows["callers.tsv"]}
    iv = {(r["kind"], r["set"], r["dataset"], r["caller"], r["measure"]): r for r in rows["intervals.tsv"]}
    margin = float(pairs(rows["settings.tsv"])["margin"])
    min_class = int(pairs(rows["settings.tsv"])["min_class"])
    meas, ends, members = [], [], {}
    rate = {}
    for r in rows["callers.tsv"]:
        rate[(r["way"], r["dataset"], r["caller"])] = rates(r["tp"], r["fp"], r["tn"], r["fn"])
    for r in rows["callers.tsv"]:
        k = (r["way"], r["dataset"], r["caller"])
        s, p, b = rate[k]
        minus = None if r["caller"] == "tool" else b - rate[(r["way"], r["dataset"], "tool")][2]
        meas.append({"way": r["way"], "dataset": r["dataset"], "caller": r["caller"], "sensitivity": s,
                     "specificity": p, "balanced_accuracy": b, "minus_tool": minus})
    for w in WAYS:
        members[(w, "by_name")], members[(w, "graded")] = [], []
        for d in named:
            t = by[(w, d, "tool")]
            kind = reps[d]["kind"]
            scored = t["tp"] + t["fn"] >= min_class and t["tn"] + t["fp"] >= min_class
            in_name = kind == "primary" and scored
            in_graded = in_name and reps[d]["graded"] == "yes"
            if in_name:
                members[(w, "by_name")].append(d)
            if in_graded:
                members[(w, "graded")].append(d)
            tb = rate[(w, d, "tool")][2]
            best = max(LEARNERS, key=lambda l: (rate[(w, d, f"{l}:lineage")][2], -LEARNERS.index(l)))
            bb = rate[(w, d, f"{best}:lineage")][2]
            ncr, ncs = t["no_call_resistant"], t["no_call_susceptible"]
            ends.append({"way": w, "dataset": d, "kind": kind, "scored": "yes" if scored else "no",
                         "by_name": "yes" if in_name else "no", "graded": "yes" if in_graded else "no",
                         "tool_balanced_accuracy": tb, "best_lineage_learner": best,
                         "best_lineage_balanced_accuracy": bb, "gap": bb - tb,
                         "within_margin": "yes" if bb - tb <= margin else "no",
                         "sensitivity_lower": ratio(t["tp"], t["tp"] + t["fn"] + ncr),
                         "sensitivity_upper": ratio(t["tp"] + ncr, t["tp"] + t["fn"] + ncr),
                         "specificity_lower": ratio(t["tn"], t["tn"] + t["fp"] + ncs),
                         "specificity_upper": ratio(t["tn"] + ncs, t["tn"] + t["fp"] + ncs)})
    med_rows, cmp_rows, scr_rows, bnd_rows = [], [], [], []
    endrow = {(e["way"], e["dataset"]): e for e in ends}
    for w in WAYS:
        for s in SETS:
            ds = members[(w, s)]
            need(ds, f"no endpoint in the set {s} under {w}")
            for c in BOOT_CALLERS:
                for mi, m in enumerate(MEASURES):
                    v = [rate[(w, d, c)][mi] for d in ds]
                    lo = hi = None
                    if w == PRIMARY_WAY:
                        x = iv[("median", s, "-", c, m)]
                        lo, hi = x["low"], x["high"]
                    med_rows.append({"way": w, "set": s, "caller": c, "measure": m, "endpoints": len(v),
                                     "median": median(v), "min": min(v), "max": max(v), "low": lo, "high": hi})
            for c in LEARNER_CALLERS:
                diffs = [rate[(w, d, c)][2] - rate[(w, d, "tool")][2] for d in ds]
                lo = hi = None
                if w == PRIMARY_WAY:
                    x = iv[("median_difference", s, "-", c, "minus_tool")]
                    lo, hi = x["low"], x["high"]
                cmp_rows.append({"way": w, "set": s, "caller": c, "endpoints": len(ds),
                                 "tool_higher": sum(x < 0 for x in diffs), "tool_lower": sum(x > 0 for x in diffs),
                                 "equal": sum(x == 0 for x in diffs),
                                 "tool_ahead_or_within": sum(x <= margin for x in diffs),
                                 "median_minus_tool": median(diffs), "low": lo, "high": hi,
                                 "caller_only_right": sum(by[(w, d, c)]["caller_only_right"] for d in ds),
                                 "tool_only_right": sum(by[(w, d, c)]["tool_only_right"] for d in ds)})
            scr_rows.append({"way": w, "set": s, "endpoints": len(ds),
                             "within_best": sum(endrow[(w, d)]["within_margin"] == "yes" for d in ds),
                             "needed": math.ceil(len(ds) / 2)})
            for bname in BOUNDS:
                v = [endrow[(w, d)][bname] for d in ds]
                bnd_rows.append({"way": w, "set": s, "bound": bname, "endpoints": len(v), "median": median(v),
                                 "min": min(v), "max": max(v)})
    summary = summarize(rows, by, members, rate, margin, endrow)
    return {"measures.tsv": tsv("measures.tsv", meas), "endpoints.tsv": tsv("endpoints.tsv", ends),
            "medians.tsv": tsv("medians.tsv", med_rows), "comparisons.tsv": tsv("comparisons.tsv", cmp_rows),
            "screen.tsv": tsv("screen.tsv", scr_rows), "bounds.tsv": tsv("bounds.tsv", bnd_rows),
            "summary.tsv": tsv("summary.tsv", [{"name": k, "value": v} for k, v in summary])}


def summarize(rows, by, members, rate, margin, endrow):
    """Single figures the page quotes."""
    w = PRIMARY_WAY
    out = []
    for s in SETS:
        ds = members[(w, s)]
        t = [by[(w, d, "tool")] for d in ds]
        out += [(f"{s}_endpoints", len(ds)), (f"{s}_rows", sum(r["genomes"] for r in t)),
                (f"{s}_rows_with_call", sum(r["with_call"] for r in t)),
                (f"{s}_rows_without_call", sum(r["no_call_resistant"] + r["no_call_susceptible"] for r in t)),
                (f"{s}_rows_without_call_label_resistant", sum(r["no_call_resistant"] for r in t)),
                (f"{s}_sensitivity_under_half", sum(rate[(w, d, "tool")][0] < 0.5 for d in ds))]
    reps = rows["reports.tsv"]
    out += [("proxy_endpoints", sum(r["kind"] == "proxy" for r in reps)),
            ("not_named_endpoints", sum(r["kind"] == "not named" for r in reps)),
            ("graded_tool_drugs_in_by_name", len(members[(w, "graded")])),
            ("ungraded_in_by_name", len(members[(w, "by_name")]) - len(members[(w, "graded")]))]
    tot = pairs(rows["totals.tsv"])
    reps = {r["dataset"]: r for r in rows["reports.tsv"]}
    cat = [d for d in members[(w, "by_name")] if reps[d]["catalog_row"] == "yes"]
    need(int(tot["catalog_endpoints"]) == len(cat), "totals.tsv: the catalog row's reports are not those of the plan's "
                                                   "set with a catalog row")
    need(int(tot["catalog_rows"]) == sum(by[(w, d, "tool")]["with_call"] for d in cat),
         "totals.tsv: the catalog row's rows are not the tool's rows with a call on its reports")
    out += [("catalog_endpoints", tot["catalog_endpoints"]), ("tool_minus_catalog_median", tot["tool_minus_catalog_median"]),
            ("catalog_rows", tot["catalog_rows"]), ("catalog_tool_only_right", tot["catalog_tool_only_right"]),
            ("catalog_only_right", tot["catalog_only_right"])]
    out.append(("catalog_calls_agree", int(tot["catalog_rows"]) - int(tot["catalog_tool_only_right"])
                - int(tot["catalog_only_right"])))
    sa = pairs(rows["set-aside.tsv"])
    d = SPECIES + "__" + sa["dataset"]
    best = endrow[(w, d)]["best_lineage_balanced_accuracy"]
    for name in ("as_run", "counted"):
        tp, fp, tn, fn = (int(sa[f"{k}_{name}"]) for k in ("tp", "fp", "tn", "fn"))
        s, p, b = rates(tp, fp, tn, fn)
        out += [(f"set_aside_sensitivity_{name}", s), (f"set_aside_specificity_{name}", p),
                (f"set_aside_gap_{name}", best - b)]
    out.append(("record_files", len(FILES)))
    return out


def finish(rows_text, checks, ids):
    """Derive, scan for what must not enter, and add checks.tsv, README.md and SHA256SUMS."""
    rows = {n: parse_tsv(n, rows_text[n]) for n in ROW_FILES if n not in ("checks.tsv", "inputs.sha256")}
    derived = derive(rows)
    out = dict(rows_text)
    out.update(derived)
    text = "".join(out.values())
    tokens = set(re.findall(r"[0-9A-Za-z_.]+", text)) | set(re.findall(r"[0-9.]+", text))
    found = len(ids & tokens)
    need(found == 0, f"{found} genome ids of the input set stand in the record")
    scan(out)
    checks = list(checks) + [("genome_ids_of_the_inputs_in_the_record", found)]
    out["checks.tsv"] = tsv("checks.tsv", [{"check": k, "count": v} for k, v in checks])
    out["README.md"] = readme(out)
    scan({"README.md": out["README.md"]})
    out["SHA256SUMS"] = sums(out)
    return {name: out[name] for name in FILES}


def scan(files):
    for shape in BANNED_SHAPES:
        for name, body in files.items():
            need(re.search(shape, body) is None, f"{name}: text of the shape {shape}")


def sums(out):
    return "".join(f"{sha256_bytes(out[n].encode('utf-8'))}  {n}\n" for n in sorted(out) if n != "SHA256SUMS")


# ----------------------------------------------------------------- README

DESCRIBE = {
    "reports.tsv": "Row file. One row per M. tuberculosis report of the release",
    "callers.tsv": "Row file. One row per way of reading, report and caller: the counts on the genomes with a call",
    "intervals.tsv": "Row file. The lineage interval's ends, for the primary way of reading",
    "disagreements.tsv": "Row file. One row per report: the primary way's counts by outcome and grade",
    "grades.tsv": "Row file. One row per tool drug: the shipped database's entries by source and grade",
    "order.tsv": "Row file. One row per event of the plan and the scoring run",
    "settings.tsv": "Row file. The settings the script applies and the sha256 of the plan, the addendum and the "
                    "compact file",
    "totals.tsv": "Row file. Counts over all reports, and the catalog row's figures over the reports of the plan's set "
                  "that have one",
    "set-aside.tsv": "Row file. Exploratory: the counts both ways for the one entry the tool's rule set aside",
    "checks.tsv": "The checks the build ran on the input set",
    "inputs.sha256": "The sha256 of each file of the input set",
    "measures.tsv": "Derived. One row per way of reading, report and caller",
    "endpoints.tsv": "Derived. One row per way of reading and report",
    "medians.tsv": "Derived. One row per way of reading, set, caller and measure",
    "comparisons.tsv": "Derived. One row per way of reading, set and learner row: each learner minus the tool",
    "screen.tsv": "Derived. One row per way of reading and set: the 0.02 screen",
    "bounds.tsv": "Derived. Exploratory: one row per way of reading, set and bound on every genome",
    "summary.tsv": "Derived. Single figures the page quotes",
}
GLOSS = {
    "dataset": "the report, as species and antibiotic", "antibiotic": "the report's antibiotic",
    "kind": "primary (the tool names the drug), proxy (read from another drug by a shared name stem) or not named",
    "tool_drugs": "the tool drug or drugs read, with `,` between them", "graded": "yes when the shipped database "
    "grades at least one entry of the tool drug, no when it grades none, `-` for a report that is not primary",
    "genomes": "genomes of the report", "catalog_row": "yes when the report has a catalog row",
    "way": "the way of reading the tool's output", "caller": "tool, or a learner and split",
    "with_call": "genomes with a call of the tool", "no_call_resistant": "genomes without a call, recorded resistant",
    "no_call_susceptible": "genomes without a call, recorded susceptible",
    "tp": "resistant genomes called resistant", "fp": "susceptible genomes called resistant",
    "tn": "susceptible genomes called susceptible", "fn": "resistant genomes called susceptible",
    "caller_only_right": "genomes this caller gets right and the tool does not",
    "tool_only_right": "genomes the tool gets right and this caller does not",
    "set": "by_name (the plan's set), graded (the addendum's set) or `-`", "measure": "the measure",
    "low": "the lower end of the lineage interval", "high": "the upper end of the lineage interval",
    "resamples_with_value": "resamples in which the figure has a value",
    "tool_r_label_s_assoc": "resistant call on grade Assoc w R, recorded susceptible",
    "tool_r_label_s_interim": "resistant call on grade Assoc w R - Interim alone, recorded susceptible",
    "tool_s_label_r_other": "susceptible call with an entry of another grade reported, recorded resistant",
    "tool_s_label_r_none": "susceptible call with nothing reported for the drug, recorded resistant",
    "both_r_assoc": "resistant call on grade Assoc w R, recorded resistant",
    "both_r_interim": "resistant call on grade Assoc w R - Interim alone, recorded resistant",
    "both_s_other": "susceptible call with an entry of another grade reported, recorded susceptible",
    "both_s_none": "susceptible call with nothing reported for the drug, recorded susceptible",
    "no_call_label_r": "no call, recorded resistant", "no_call_label_s": "no call, recorded susceptible",
    "tool_drug": "the tool drug", "assoc_w_r": "entries of source WHO catalogue v2, grade Assoc w R",
    "assoc_w_r_interim": "entries of source WHO catalogue v2, grade Assoc w R - Interim",
    "own_uncertain_significance": "the database's own entries, grade Uncertain significance",
    "own_not_assoc_w_r_interim": "the database's own entries, grade Not assoc w R - Interim",
    "own_ungraded": "the database's own entries with no grade",
    "counted_grade_entries": "entries of the two counted grades",
    "event": "the event", "time": "its UTC time", "source": "the file the time is read from",
    "name": "the name of the figure", "value": "the figure", "check": "a check the build ran on the input set",
    "count": "the count the check gives", "sensitivity": "`tp` over `tp` plus `fn`",
    "specificity": "`tn` over `tn` plus `fp`", "balanced_accuracy": "the mean of sensitivity and specificity",
    "minus_tool": "this caller's balanced accuracy minus the tool's, `-` for the tool",
    "scored": "yes when at least 25 resistant and 25 susceptible genomes have a call",
    "by_name": "yes for a scored primary report", "tool_balanced_accuracy": "the tool's balanced accuracy",
    "best_lineage_learner": "the learner whose lineage split row has the highest balanced accuracy",
    "best_lineage_balanced_accuracy": "that learner's balanced accuracy",
    "gap": "that balanced accuracy minus the tool's", "within_margin": "yes when `gap` is at most 0.02",
    "sensitivity_lower": "sensitivity on every genome, no call counted as wrong",
    "sensitivity_upper": "sensitivity on every genome, no call counted as right",
    "specificity_lower": "specificity on every genome, no call counted as wrong",
    "specificity_upper": "specificity on every genome, no call counted as right",
    "endpoints": "reports in the set", "median": "the median over the set's reports", "min": "the smallest",
    "max": "the largest", "tool_higher": "reports where the tool's balanced accuracy is higher",
    "tool_lower": "reports where it is lower", "equal": "reports where the two are equal",
    "tool_ahead_or_within": "reports where the learner minus the tool is at most 0.02",
    "median_minus_tool": "the median of the learner minus the tool", "within_best": "reports with `within_margin` yes",
    "needed": "half the set's reports, rounded up", "bound": "the bound on every genome",
}
GLOSS_BY_FILE = {
    ("endpoints.tsv", "kind"): "primary or proxy", ("endpoints.tsv", "graded"): "yes for a report of the graded set",
    ("medians.tsv", "low"): "the lineage interval's lower end, `-` for other ways of reading",
    ("medians.tsv", "high"): "the lineage interval's upper end, `-` for other ways of reading",
    ("comparisons.tsv", "low"): "the lineage interval's lower end, `-` for other ways of reading",
    ("comparisons.tsv", "high"): "the lineage interval's upper end, `-` for other ways of reading",
    ("comparisons.tsv", "caller"): "a learner and split",
    ("comparisons.tsv", "caller_only_right"): "summed over the set", ("comparisons.tsv", "tool_only_right"): "summed over the set",
    ("intervals.tsv", "kind"): "endpoint, difference, median or median_difference",
    ("intervals.tsv", "dataset"): "the report, or `-` for a median",
    ("intervals.tsv", "measure"): "sensitivity, specificity, balanced_accuracy or minus_tool",
    ("bounds.tsv", "median"): "the median over the set's reports",
    ("settings.tsv", "name"): "the name of a setting", ("settings.tsv", "value"): "the setting",
    ("totals.tsv", "name"): "the name of a count over all reports or of a catalog row figure",
    ("summary.tsv", "name"): "the name of a figure",
    ("set-aside.tsv", "name"): "the name of a count", ("set-aside.tsv", "value"): "the count",
}


def readme(out):
    def count(n):
        return out[n].count("\n") - (0 if n == "inputs.sha256" else 1)

    lines = ["# TBProfiler comparison record", "",
             "This directory holds the record of the comparison of TBProfiler 6.6.3, run with the database "
             "it ships on assemblies of the genomes of the revision 14 M. tuberculosis reports, with the phenotypes "
             "those reports record and with the five learners on the same genomes. The page that reads it is "
             "[`docs/tbprofiler-comparison.md`](../../docs/tbprofiler-comparison.md). "
             "`scripts/tbprofiler-comparison/record.py` writes every file here from the input set, and its `--check` "
             "mode fails on any byte it would not write.", "",
             "The input set is the call file of the TBProfiler run on assemblies of public BV-BRC genomes, with GRM's "
             "revision 14 M. tuberculosis reports.", "",
             "## Files", "", "| File | Rows | Contents |", "|---|---|---|"]
    for n in ROW_FILES + DERIVED_FILES:
        lines.append(f"| `{n}` | {count(n):,} | {DESCRIBE[n]} |")
    lines += ["", "`SHA256SUMS` gives the sha256 of every other file, and this README is what the script writes.", "",
              "## Columns", "",
              "Each line of `inputs.sha256` gives the sha256 of a file of the input set and its name. The other files' "
              "columns, in file order:", "", "| File | Column | Gives |", "|---|---|---|"]
    for n in ROW_FILES + DERIVED_FILES:
        for c in COLUMNS.get(n, []):
            gloss = GLOSS_BY_FILE.get((n, c), GLOSS.get(c))
            need(gloss is not None, f"README: no words for the column {c} of {n}")
            lines.append(f"| `{n}` | `{c}` | {gloss} |")
    lines += ["", "## Rules", "",
              "**Calls.** The four ways of reading are `who_assoc` (an entry graded Assoc w R or Assoc w R - Interim, "
              "the primary way), `who_assoc_only` (Assoc w R alone), `tool_default` (any entry the tool reports for "
              "the drug) and `who_assoc_small` (the primary way without large changes). A genome without a call is "
              "never read as susceptible: every caller is counted on the genomes where the tool has a call.", "",
              "**Learners.** A learner calls resistant at a report score of 0.5 or more, on the lineage and the "
              "random split. The catalog row is GRM's own reader of the AMRFinderPlus catalog.", "",
              "**Catalog row.** The catalog row appears only as six-report totals, in `totals.tsv` and again in "
              "`summary.tsv`, because on three reports a per-report count of it would describe one genome. Its median "
              "is given to 3 decimals for the same reason: in full, with the other totals, it would show whether one of "
              "those genomes had a call.", "",
              "**Sets.** A primary report is scored when at least 25 resistant and 25 susceptible genomes have a "
              "call. `by_name` holds the scored primary reports, the plan's set; `graded` holds those whose tool "
              "drug the database grades, the addendum's set. Proxy reports enter no set.", "",
              "**Interval.** 2,000 resamples of the lineage groups, the reports' lineage ids joined across all the "
              "reports and numbered in the text order of their smallest genome id. Each resample draws as many group "
              "numbers as there are groups with Python's `random.Random(20261114).randrange`, one call per draw, "
              "and takes every genome of each group drawn. A figure has no value in a resample where a class has no "
              "genome. `low` and `high` are the values at sorted index int(0.025 n) and n - 1 - int(0.025 n) of the "
              "n resamples with a value: the 51st and the 1,950th when n is 2,000. Medians skip a report without "
              "a value.", "",
              "**Rounding.** Every figure here is the full value, apart from the catalog row's median. The page rounds each "
              "one once.", "",
              "## Limits", "",
              "No row is about one genome, and no cell names a genome. The tool's raw output and the assemblies "
              "no longer exist; the call file, kept outside the repository, is the lowest layer the counts are "
              "rebuilt from. `set-aside.tsv` and the counts behind it rest on the raw output alone and cannot be "
              "rebuilt. `--check` alone cannot hold a row file to the input set; that needs the input set.", "",
              "## Guards", "",
              "`--check` with the input set rebuilds every file and fails on any byte that differs. `--check` alone "
              "reads the row files, derives every derived file and this README again, and fails on any byte that "
              "differs or on a file that differs from `SHA256SUMS`. The input set is not in the repository: it holds "
              "the call file and the reports' scores, which are per genome.", "",
              "`--out` writes only to a new path, an empty directory or an earlier record of this script whose files "
              "agree with its own `SHA256SUMS`, and refuses a symlink.", "",
              "## Commands", "", "```sh",
              "python3 -I scripts/tbprofiler-comparison/record.py --inputs <dir> --out reference/tbprofiler-comparison",
              "python3 -I scripts/tbprofiler-comparison/record.py --inputs <dir> --check reference/tbprofiler-comparison",
              "python3 -I scripts/tbprofiler-comparison/record.py --check reference/tbprofiler-comparison",
              "python3 -I -m unittest discover -s tests/tbprofiler_comparison -t .", "```", ""]
    return "\n".join(lines)


# ----------------------------------------------------------------- record on disk

def read_record(R):
    R = Path(R)
    need(R.is_dir(), f"{R}: no record")
    present = sorted(p.name for p in R.iterdir())
    need(present == sorted(FILES), f"{R}: files {sorted(set(present) ^ set(FILES))} differ from the record's list")
    return {n: (R / n).read_text(encoding="utf-8") for n in FILES}


def check_alone(R):
    """Derive from the record's own row files and compare every derived byte."""
    have = read_record(R)
    rows = {n: parse_tsv(n, have[n]) for n in ROW_FILES if n not in ("checks.tsv", "inputs.sha256")}
    parse_tsv("checks.tsv", have["checks.tsv"])
    need(all(re.fullmatch(r"[0-9a-f]{64}  \S+", x) for x in have["inputs.sha256"].rstrip("\n").split("\n")),
         "inputs.sha256: a line that is not a sha256 and a name")
    want = dict(have)
    want.update(derive(rows))
    scan({n: want[n] for n in FILES if n != "SHA256SUMS"})
    want["README.md"] = readme(want)
    want["SHA256SUMS"] = sums(want)
    bad = [n for n in FILES if want[n] != have[n]]
    need(not bad, "differs from what the row files give: " + ", ".join(bad))
    return len(FILES)


def check_out(out):
    """--out takes a new path, an empty directory, or an earlier record of this script: its own files, none a
    symlink, each with the sha256 its own SHA256SUMS gives. Anything else is left as it is."""
    need(not out.is_symlink(), f"{out}: a symlink")
    if not out.exists():
        return
    need(out.is_dir(), f"{out}: not a directory")
    present = sorted(out.iterdir())
    if not present:
        return
    need(sorted(p.name for p in present) == sorted(FILES), f"{out}: holds files this script did not write")
    for p in present:
        need(not p.is_symlink() and p.is_file(), f"{out}/{p.name}: a symlink or not a plain file")
    listed = {}
    for line in (out / "SHA256SUMS").read_bytes().decode("utf-8", "replace").split("\n")[:-1]:
        m = re.fullmatch(r"([0-9a-f]{64})  (\S+)", line)
        need(m is not None and m.group(2) not in listed, f"{out}/SHA256SUMS: not one sha256 and name per line")
        listed[m.group(2)] = m.group(1)
    need(sorted(listed) == sorted(n for n in FILES if n != "SHA256SUMS"),
         f"{out}/SHA256SUMS: does not list the record's files")
    bad = [n for n, h in sorted(listed.items()) if sha256_bytes((out / n).read_bytes()) != h]
    need(not bad, f"{out}: holds files this script did not write, or changed since: " + ", ".join(bad))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--inputs")
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--out")
    g.add_argument("--check")
    a = ap.parse_args(argv)
    try:
        if a.out:
            need(a.inputs is not None, "--out needs --inputs")
            out = Path(a.out)
            check_out(out)
            record = build(a.inputs)
            out.mkdir(parents=True, exist_ok=True)
            for n, text in record.items():
                (out / n).write_text(text, encoding="utf-8")
            print(f"wrote {len(record)} files to {out}")
        elif a.inputs:
            record = build(a.inputs)
            have = read_record(a.check)
            bad = [n for n in FILES if have[n] != record[n]]
            need(not bad, "differs from the input set's record: " + ", ".join(bad))
            print(f"{len(FILES)} files rebuilt from the input set, byte for byte")
        else:
            n = check_alone(a.check)
            print(f"{n} files agree with the row files and SHA256SUMS")
    except Stop as e:
        print(f"record.py: {e}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
