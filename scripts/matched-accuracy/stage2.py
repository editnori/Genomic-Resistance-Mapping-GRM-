"""Stage 2 of the matched accuracy protocol (docs/matched-accuracy.md): the five learners' mean score, compared
inside each matched group of stage 1's frozen rows.

  python3 -I stage2.py --reports DIR --rows-a FILE --rows-b FILE --hc20 FILE --out DIR

Inputs. DIR holds the release's reports, each checked against its line of the release's sums list as stage 1
checks it. --rows-a and --rows-b are stage 1's two row lists, each checked against the sha256 stage 1 printed in
its one real run (ruling R59), written below as constants that no argument changes. --hc20 is the 5 Nov genome
table, checked against the constant that equals its line of stage 1's committed inputs table. Stage 1's count
table and typing table are read from their committed copies under
reference/matched-accuracy/stage1-rev14-274-20261112/, each checked against a constant.

Steps.
1. Checks before any score. Each failure prints one line and exits 1 with no export: a sha256 that differs (a
   report, a row list, the committed count or typing table, the genome table); numpy other than 1.26.4; a repeated
   line or lines out of stage 1's sort; a report whose id is not its file name; a listed report without scm,
   cart, rf, gbt and lasso in that order under the lineage split; a listed row its report does not hold, or whose
   fold or lineage cluster differs from the report's; a group with one lab outcome; groups, pairs, rows, genomes,
   datasets or species, rebuilt from a list, that differ from stage 1's count table at any level; the support of
   the MLST read or the MLST and HC20 read, rebuilt from a list, that differs from stage 1's typing table at any
   level; in a listed report, a dataset or species name with a full stop, a listed row's catalogue call other
   than true, false or null, or a fit size by fold (fold_complexity) of the five other than 5 finite numbers. Up
   to here reports load through stage 1's loader, which drops every score, and through count_bad, which returns
   counts alone.
2. The line. Write "scores opened: yes" to <out>/scores-opened.txt, then print it. After it, a run stops only on
   a report that cannot be read or whose bytes changed, a score null, true or false, not a number or not finite,
   a hit of the export's scan, or a write that fails.
3. The second loader. For listed rows alone: the five scores, the five fold_complexity lists, the row's
   catalogue call, and whether the report's top-level catalogue is present or null.
4. Primary measure. A row's mean is its five scores added left to right as IEEE doubles, in the order scm,
   cart, rf, gbt, lasso, then divided by 5. Inside a group each resistant row meets each susceptible row: a
   higher resistant mean counts 1, an exact tie one half. A group's value is that count over its pairs; groups
   weigh equally inside a dataset, datasets inside a species, species overall, as exact fractions. With 6 or
   more species, the species means are resampled 10,000 times (numpy PCG64, seed 42) and the linear 2.5th and
   97.5th percentiles bound the estimate.
5. Secondary measures, each with the primary weights: the four cells of catalogue calls, the typing reads
   (MLST, then MLST and HC20, each beside the wider match on the same surviving rows) and the read without
   groups in a fold with a null fit.
6. Export. Numbers as integers alone, estimates in parts per million; a dataset or species under 10 pairs gives
   its counts and no estimate. Beside each interval stands one word read from its two bounds before they are
   rounded: above when the lower bound is above one half, below when the upper bound is below one half, else
   unresolved. Before any file is written, the export's scan: no file may hold a full stop, then a scan for
   every genome id and BioProject accession the run read, after a positive control.

Python 3.12, standard library, with numpy for the resampling alone; run with -I.
"""

import argparse
import collections
import fractions
import functools
import hashlib
import importlib.util
import json
import math
import operator
import os
import re
import sys
from pathlib import Path

import numpy

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def _stage1():
    spec = importlib.util.spec_from_file_location("matched_accuracy_stage1_of_stage2", HERE / "stage1.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


stage1 = _stage1()
Stop = stage1.Stop
FIVE = stage1.FIVE

# The sha256 of each row list as stage 1 printed it in its one real run (ruling R59). No argument changes them.
ROW_LIST_SHA256 = {"A": "b135a8ccd24f3c9b87113498c58bb76181fc8c143135260df01eafbd8ad7dd26",
                   "B": "e3558b1b5a7a4a9f8ff23c79aca502c25586fcd6f889e4dac8d788cdd0b45d66"}
# Stage 1's export of its one real run, committed byte for byte; stage 2 reads two of its files.
STAGE1_EXPORT = ROOT / "reference" / "matched-accuracy" / "stage1-rev14-274-20261112"
COUNTS = STAGE1_EXPORT / "counts.tsv"
COUNTS_SHA256 = "d686f07a3d3d34985b15a7be86e61c13aeac6c026fb925eba25c53f1496ced3e"
TYPING_TABLE = STAGE1_EXPORT / "typing.tsv"
TYPING_TABLE_SHA256 = "d25d1fedfbda087f40f91d4ca026705fa8d1d00808de309e06280629f90982fb"
# The 5 Nov genome table's sha256, equal to its line of the committed inputs.tsv of stage 1. No argument changes it.
HC20_SHA256 = "d60e2b15f2cc0ed0466a1016da1761869fb42cc7e9fdd817c1ff7051e759a80e"
NUMPY_VERSION = "1.26.4"
SEED = 42
DRAWS = 10000
PERCENTILES = (2.5, 97.5)
MIN_SPECIES = 6
MIN_PAIRS = 10
PPM = 1000000
OPENED = "scores opened: yes"
OPENED_FILE = "scores-opened.txt"
READINGS = ("A", "B")
LEVELS = ("overall", "species", "dataset")
SUPPORT = ("rows", "genomes", "groups", "pairs", "datasets", "species")
ROW_FIELDS = ("dataset", "genome", "fold", "lineage", "bioproject", "method", "standard", "year")
CALLS = {True: "resistant", False: "susceptible"}
CELLS = ((True, True), (True, False), (False, True), (False, False))
APART = ("unevaluable", "no drug rule")
TYPING = (("mlst", ("mlst",), (("mlst", ("mlst",)), ("primary", ()))),
          ("mlst_hc20", ("mlst", "hc20"), (("mlst_hc20", ("mlst", "hc20")), ("mlst", ("mlst",)), ("primary", ()))))
# The checks that run before the line, by name: check() alone calls them, and run() calls check() before
# open_scores(). A test reads this list.
CHECKS = ("check_out", "check_numpy", "check_count_table", "check_typing_table", "check_row_list", "check_hc20",
          "report_files", "parse_rows", "load_reports", "read_hc20", "match_rows", "check_one_outcome",
          "check_counts", "check_typing", "check_values")
# The counts count_bad returns for one listed report, in the order check_values reads them, each with its stop.
BAD = (("moved", "listed rows that moved after they were checked: {n}"),
       ("name", "its dataset or species name holds a full stop, which no export file may hold"),
       ("catalog_resistant", "listed rows whose catalog_resistant is not true, false or null: {n}")) + \
    tuple((f"fold_complexity of {n}", f"the fit size by fold (fold_complexity) of {n} is not a list of 5 finite "
                                      f"numbers") for n in FIVE)
# The columns of stage 1's typing table a row list rebuilds: all but the rows lacking a value among the rows that
# enter the groups, which counts rows outside every eligible group, and no row list holds those.
TYPING_REBUILT = ("lacking_value_of_eligible_rows", "rows", "genomes", "groups", "pairs", "datasets", "species",
                  "primary_groups", "primary_pairs", "mlst_groups", "mlst_pairs", "rows_lost_from_primary")
MEASURE = ("reading", "level", "name", "estimate_ppm", "rows", "groups", "pairs", "datasets", "species")
LOST = ("lost_rows", "lost_groups", "lost_pairs", "lost_datasets", "lost_species")
INTERVAL = ("reading", "read", "match", "lower_ppm", "upper_ppm", "chance", "species", "draws")
CELL_COLUMNS = ("reading", "resistant_member_call", "susceptible_member_call", "estimate_ppm", "pairs", "groups",
                "datasets", "species")
TYPING_COLUMNS = ("reading", "step", "match", "level", "name", "estimate_ppm", "rows", "groups", "pairs",
                  "datasets", "species") + LOST
NULL_COLUMNS = MEASURE + LOST


# 1. The checks before any score.

def check_out(out):
    """The output directory: absent, or a directory without export, local or the line of an earlier run."""
    out = Path(out)
    if out.exists() and not out.is_dir():
        raise Stop("the output path is not a directory")
    if any((out / name).exists() for name in ("export", "local", OPENED_FILE)):
        raise Stop("the output directory already holds export, local or the line of an earlier run")


def check_numpy():
    if numpy.__version__ != NUMPY_VERSION:
        raise Stop("numpy is not the version the resampling is written for")


def check_count_table():
    """{(level, name, reading): {count: int}} of stage 1's committed count table, checked against its sha256."""
    path = Path(COUNTS)
    if not path.is_file():
        raise Stop("stage 1's committed count table is missing")
    if stage1.sha256(path) != COUNTS_SHA256:
        raise Stop("stage 1's committed count table differs from its sha256")
    lines = path.read_bytes().decode("utf-8").split("\n")
    if lines[-1] != "" or lines[0].split("\t") != list(stage1.COUNT_COLUMNS):
        raise Stop("stage 1's committed count table does not have stage 1's columns")
    table = {}
    for line in lines[1:-1]:
        cells = dict(zip(stage1.COUNT_COLUMNS, line.split("\t")))
        if len(line.split("\t")) != len(stage1.COUNT_COLUMNS) or cells["level"] not in LEVELS \
                or cells["reading"] not in READINGS or not all(cells[c].isdigit() for c in SUPPORT):
            raise Stop("a line of stage 1's committed count table is not a level, a name, a reading and counts")
        k = (cells["level"], cells["name"], cells["reading"])
        if k in table:
            raise Stop("stage 1's committed count table gives one level, name and reading twice")
        table[k] = {c: int(cells[c]) for c in SUPPORT}
    return table


def check_typing_table():
    """{(level, name, reading, step): {column: int or ""}} of stage 1's committed typing table, checked against its
    sha256. The MLST read gives mlst_groups and mlst_pairs blank; every other cell is an integer."""
    path = Path(TYPING_TABLE)
    if not path.is_file():
        raise Stop("stage 1's committed typing table is missing")
    if stage1.sha256(path) != TYPING_TABLE_SHA256:
        raise Stop("stage 1's committed typing table differs from its sha256")
    lines = path.read_bytes().decode("utf-8").split("\n")
    if lines[-1] != "" or lines[0].split("\t") != list(stage1.TYPING_COLUMNS):
        raise Stop("stage 1's committed typing table does not have stage 1's columns")
    table = {}
    for line in lines[1:-1]:
        fields = line.split("\t")
        cells = dict(zip(stage1.TYPING_COLUMNS, fields))
        numbers = stage1.TYPING_COLUMNS[4:]
        blank = ("mlst_groups", "mlst_pairs") if cells.get("step") == "mlst" else ()
        if len(fields) != len(stage1.TYPING_COLUMNS) or cells["level"] not in LEVELS \
                or cells["reading"] not in READINGS or cells["step"] not in [t[0] for t in TYPING] \
                or not all(cells[c] == "" if c in blank else cells[c].isdigit() for c in numbers):
            raise Stop("a line of stage 1's committed typing table is not a level, a name, a reading, a step and "
                       "counts")
        k = (cells["level"], cells["name"], cells["reading"], cells["step"])
        if k in table:
            raise Stop("stage 1's committed typing table gives one level, name, reading and step twice")
        table[k] = {c: "" if c in blank else int(cells[c]) for c in numbers}
    return table


def check_row_list(path, reading):
    """The row list's bytes, checked against the sha256 stage 1 printed for the reading."""
    path = Path(path)
    if not path.is_file():
        raise Stop(f"the row list of reading {reading} is missing")
    data = path.read_bytes()
    if hashlib.sha256(data).hexdigest() != ROW_LIST_SHA256[reading]:
        raise Stop(f"the row list of reading {reading} differs from the sha256 stage 1 printed")
    return data


def check_hc20(path):
    path = Path(path)
    if not path.is_file():
        raise Stop("the genome table, 5 Nov, is missing")
    if stage1.sha256(path) != HC20_SHA256:
        raise Stop("the genome table, 5 Nov, differs from its sha256")


def report_files(directory):
    return stage1.report_files(directory)


def load_reports(files):
    """(reports that enter, reports left out), through stage 1's loader one file at a time: a report whose id is
    not its file name stops the run."""
    entering, left_out = [], []
    for name, path in files:
        e, o = stage1.load_reports([(name, path)])
        for rep in e + o:
            if rep["id"] != name:
                raise Stop(f"{name}.json gives the id {rep['id']}, not its file name")
        entering += e
        left_out += o
    return entering, left_out


def parse_rows(data, reading):
    """[row] of one row list, as stage 1's row_list writes it: 8 tab-separated fields a line, sorted, UTF-8, LF."""
    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError:
        raise Stop(f"the row list of reading {reading} is not UTF-8") from None
    if not text:
        raise Stop(f"the row list of reading {reading} holds no line")
    if not text.endswith("\n"):
        raise Stop(f"the row list of reading {reading} does not end with a line break")
    lines = [line + "\n" for line in text[:-1].split("\n")]
    rows, seen = [], set()
    for i, line in enumerate(lines):
        if i and line == lines[i - 1]:
            raise Stop(f"a line of the row list of reading {reading} is repeated")
        if i and line < lines[i - 1]:
            raise Stop(f"the lines of the row list of reading {reading} are out of stage 1's sort")
        fields = line[:-1].split("\t")
        if len(fields) != len(ROW_FIELDS) or "\r" in line:
            raise Stop(f"a line of the row list of reading {reading} does not hold {len(ROW_FIELDS)} fields")
        r = dict(zip(ROW_FIELDS, fields))
        if not all(r[k] for k in ROW_FIELDS[:-1] if k not in ("fold", "lineage")) \
                or not all(re.fullmatch(r"-?[0-9]+", r[k], re.ASCII) and str(int(r[k])) == r[k]
                           for k in ("fold", "lineage")):
            raise Stop(f"a line of the row list of reading {reading} has a blank field or a fold or lineage "
                       f"cluster that is not an integer")
        r["fold"], r["lineage"], r["year"] = int(r["fold"]), int(r["lineage"]), r["year"] or None
        if (r["dataset"], r["genome"]) in seen:
            raise Stop(f"the row list of reading {reading} gives one genome of one dataset twice")
        seen.add((r["dataset"], r["genome"]))
        rows.append(r)
    return rows


def match_rows(rows, reports, reading):
    """Each listed row beside its report's row: the report must hold the five learners in order under the lineage
    split and the genome, with the same fold and lineage cluster. Adds species, lab call, MLST and row index."""
    out = []
    for r in rows:
        rep = reports.get(r["dataset"])
        if rep is None:
            raise Stop(f"reading {reading} lists a dataset no report gives")
        names = [e["learner"] for e in rep["lineage_split"]["learners"]]
        if [n for n in names if n in FIVE] != list(FIVE):
            raise Stop(f"reading {reading} lists a report without scm, cart, rf, gbt, lasso in that order under "
                       f"the lineage split")
        i = rep["index"].get(r["genome"])
        if i is None:
            raise Stop(f"reading {reading} lists a row its report does not hold")
        s = rep["samples"][i]
        if rep["lineage_split"]["test_fold"][i] != r["fold"]:
            raise Stop(f"reading {reading} lists a row whose fold differs from its report's")
        if s["lineage"] != r["lineage"]:
            raise Stop(f"reading {reading} lists a row whose lineage cluster differs from its report's")
        out.append(dict(r, species=rep["species_slug"], resistant=s["resistant"],
                        mlst=(s["mlst"] or "").strip() or None, index=i))
    return out


def check_one_outcome(rows, reading):
    """{group key: rows} of the groups rebuilt from a list; a group with one lab outcome stops the run."""
    groups, rest = stage1.groups_of(rows)
    if rest:
        raise Stop(f"reading {reading}: a group rebuilt from the row list holds one lab outcome")
    return groups


def check_counts(groups, table, reading):
    """Rows, genomes, groups, pairs, datasets and species of the rebuilt groups, at each level, against stage 1's
    count table."""
    sup = stage1.support(groups)
    zero = dict.fromkeys(SUPPORT, 0)
    for level in LEVELS:
        mine = sup.get(level, {})
        theirs = {name: c for (lv, name, rd), c in table.items() if lv == level and rd == reading}
        for name in sorted(set(mine) | set(theirs)):
            if name not in theirs or mine.get(name, zero) != theirs[name]:
                raise Stop(f"reading {reading}: the counts rebuilt from the row list differ from stage 1's count "
                           f"table at the {level} level")


def check_typing(rows, groups, hc20, table, reading):
    """The support of the MLST read and of the MLST and HC20 read, rebuilt from a list with stage 1's typing(), at
    each level, against stage 1's typing table. A name either side lacks counts as zero on that side, and a name
    stage 1's table lacks stops the run."""
    theirs = {k: c for k, c in table.items() if k[2] == reading}
    names = {level: {k[1] for k in theirs if k[0] == level} | set(stage1.support(groups).get(level, {}))
             for level in LEVELS}
    mine = stage1.typing(rows, groups, hc20, names)
    for step, _, _ in TYPING:
        for level in LEVELS:
            for name in sorted(names[level]):
                k = (level, name, reading, step)
                if k not in theirs or {c: mine[step][level][name][c] for c in TYPING_REBUILT} \
                        != {c: theirs[k][c] for c in TYPING_REBUILT}:
                    raise Stop(f"reading {reading}: the {step} support rebuilt from the row list differs from stage "
                               f"1's typing table at the {level} level")


def read_hc20(path, wanted):
    return stage1.read_hc20(path, wanted)


def finite(x):
    """True for a number, not true or false, that is finite as a double."""
    if isinstance(x, bool) or not isinstance(x, (int, float)):
        return False
    try:
        return math.isfinite(float(x))
    except OverflowError:
        return False


def count_bad(path, digest, listed):
    """The loader before the line: {kind: count} of the bad values of one listed report, the kinds of BAD. It reads
    the report's bytes checked against the sha256 checked before, the listed rows' ids and catalogue calls, the
    names and each of the five learners' fit sizes by fold, and returns integers alone: no score, no catalogue call
    and no fit size."""
    data = Path(path).read_bytes()
    if hashlib.sha256(data).hexdigest() != digest:
        raise Stop(f"{Path(path).name} changed after it was checked")
    doc = json.loads(data)
    try:
        split = next(s for s in doc["splits"] if isinstance(s, dict) and s.get("scheme") == "lineage")
        entries = {e["learner"]: e for e in split["learners"] if e.get("learner") in FIVE}
        rows = [doc["samples"][i] if i < len(doc["samples"]) else {} for _, i in listed]
        out = {"moved": sum(1 for (genome, _), s in zip(listed, rows) if s.get("id") != genome),
               "name": sum(1 for x in (doc.get("id"), doc.get("species_slug")) if "." in x),
               "catalog_resistant": sum(1 for s in rows if not (s.get("catalog_resistant") is None
                                                                or isinstance(s.get("catalog_resistant"), bool)))}
        for n in FIVE:
            c = entries[n].get("fold_complexity")
            out[f"fold_complexity of {n}"] = int(not (isinstance(c, list) and len(c) == stage1.FOLDS
                                                      and all(finite(x) for x in c)))
        return out
    finally:
        del doc


def check_values(listed):
    """Each listed report's counts from count_bad: any count above 0 stops the run, naming the report and what."""
    for d, (path, digest, rows) in listed.items():
        bad = count_bad(path, digest, rows)
        for kind, words in BAD:
            if bad[kind]:
                raise Stop(f"{d}: " + words.format(n=bad[kind]))


def check(args):
    """Every check before any score. Returns what the measure needs, and no score."""
    check_out(args.out)
    check_numpy()
    table = check_count_table()
    typing_table = check_typing_table()
    data = {name: check_row_list(path, name) for name, path in (("A", args.rows_a), ("B", args.rows_b))}
    check_hc20(args.hc20)
    files = report_files(args.reports)
    report_sums = {name: digest for digest, _, name in
                   (line.partition("  ") for line in Path(stage1.REPORT_SUMS).read_text(encoding="utf-8").splitlines())}
    lists = {name: parse_rows(data[name], name) for name in READINGS}
    entering, left_out = load_reports(files)
    reports = {}
    for rep in entering + left_out:
        rep["index"] = {s["id"]: i for i, s in enumerate(rep["samples"])}
        reports[rep["id"]] = rep
    wanted = {r["genome"] for name in READINGS for r in lists[name]}
    hc20, hc20_ids, hc20_rows, hc20_missing = read_hc20(args.hc20, wanted)
    rows, groups = {}, {}
    for name in READINGS:
        rows[name] = match_rows(lists[name], reports, name)
        groups[name] = check_one_outcome(rows[name], name)
        check_counts(groups[name], table, name)
        check_typing(rows[name], groups[name], hc20, typing_table, name)
    paths = dict(files)
    listed = {d: (paths[d], report_sums["./datasets/" + d + ".json"],
                  sorted({(r["genome"], r["index"]) for name in READINGS for r in rows[name] if r["dataset"] == d},
                         key=lambda x: x[1]))
              for d in sorted({r["dataset"] for name in READINGS for r in rows[name]})}
    check_values(listed)
    return {
        "rows": rows, "groups": groups, "hc20": hc20, "listed": listed,
        "sums": {"reports": stage1.REPORT_SUMS_SHA256, "counts": COUNTS_SHA256, "typing": TYPING_TABLE_SHA256,
                 "A": hashlib.sha256(data["A"]).hexdigest(), "B": hashlib.sha256(data["B"]).hexdigest(),
                 "hc20": HC20_SHA256},
        "facts": {"reports": len(files), "entering": len(entering), "table_lines": len(table),
                  "typing_lines": len(typing_table),
                  "entering_rows": sum(len(rep["samples"]) for rep in entering),
                  "entering_species": len({rep["species_slug"] for rep in entering}),
                  "lines": {name: len(lists[name]) for name in READINGS}, "hc20_rows": hc20_rows,
                  "hc20_missing": hc20_missing, "listed": len(listed)},
        "ids": {s["id"] for rep in entering + left_out for s in rep["samples"]} | wanted | hc20_ids,
        "accessions": {p for name in READINGS for r in rows[name] for p in stage1.bioproject_parts(r["bioproject"])
                       if stage1.is_accession(p)},
    }


# 2. The line.

class Opened:
    """Made by open_scores once the line is written and printed. The second loader takes nothing else."""

    def __init__(self, path):
        self.path = path


def open_scores(out):
    """Write the line to <out>/scores-opened.txt, then print it. A write that fails stops the run before it, and
    first removes a line file this run made, so a rerun into <out> can start."""
    out = Path(out)
    made = None
    try:
        out.mkdir(parents=True, exist_ok=True)
        fd = os.open(out / OPENED_FILE, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
        made = out / OPENED_FILE
        with os.fdopen(fd, "w", encoding="utf-8") as f:
            f.write(OPENED + "\n")
    except OSError as e:
        kept = ""
        if made is not None:
            try:
                made.unlink()
            except OSError:
                kept = ", and its file could not be removed"
        raise Stop(f"the line could not be written ({type(e).__name__}){kept}") from None
    print(OPENED, flush=True)
    return Opened(out / OPENED_FILE)


# 3. The second loader.

def component(x):
    """A score as a double. Null, true or false, a value that is not a number, or one not finite stops the run."""
    if isinstance(x, bool) or not isinstance(x, (int, float)):
        raise Stop("a learner's score is null, true or false, or not a number")
    try:
        x = float(x)
    except OverflowError:
        raise Stop("a learner's score is not finite") from None
    if not math.isfinite(x):
        raise Stop("a learner's score is not finite")
    return x


def read_scores(path, digest, listed, opened):
    """The second loader: for the listed rows of one report alone, the five scores and the catalogue call; the
    five fold_complexity lists; the top-level catalogue as present or None. It takes the token open_scores makes
    and does not read the line file again. A report that cannot be read stops the run, and its bytes are checked
    again against the sha256 checked before the line."""
    if not isinstance(opened, Opened):
        raise Stop("a score was asked for before the line")
    try:
        data = Path(path).read_bytes()
    except OSError as e:
        raise Stop(f"{Path(path).name} could not be read after the line ({type(e).__name__})") from None
    if hashlib.sha256(data).hexdigest() != digest:
        raise Stop(f"{Path(path).name} changed after it was checked")
    doc = json.loads(data)
    try:
        return pick(doc, listed)
    finally:
        del doc


def pick(doc, listed):
    """The values of the listed rows. Each row's id, catalogue call and the fit sizes were checked before the line
    by count_bad, on the bytes read_scores checks again, so the scores alone can stop the run here."""
    split = next(s for s in doc["splits"] if isinstance(s, dict) and s.get("scheme") == "lineage")
    entries = {e["learner"]: e for e in split["learners"] if e.get("learner") in FIVE}
    rows = {}
    for genome, i in listed:
        call = doc["samples"][i].get("catalog_resistant")
        rows[genome] = (tuple(component(entries[n]["scores"][i]) for n in FIVE), call)
    complexity = {n: tuple(float(x) for x in entries[n]["fold_complexity"]) for n in FIVE}
    return {"rows": rows, "fold_complexity": complexity, "catalog": "present" if doc.get("catalog") is not None
            else None}


# 4. The measure.

def row_mean(five):
    """The five scores added left to right as IEEE doubles, in the order scm, cart, rf, gbt, lasso, then / 5."""
    return functools.reduce(operator.add, [component(x) for x in five]) / 5


def points(r, s):
    """Half points of one pair, a resistant mean r against a susceptible mean s: higher 2, an exact tie 1, else 0."""
    return 2 if r > s else 1 if r == s else 0


def group_count(resistant, susceptible):
    """(half points, pairs): every resistant mean meets every susceptible mean."""
    half = 0
    for r in resistant:
        for s in susceptible:
            half += points(r, s)
    return half, len(resistant) * len(susceptible)


def compare(groups):
    """{group key: (species, dataset, half points, pairs, rows)} of groups of rows that carry a mean."""
    out = {}
    for k, g in groups.items():
        half, pairs = group_count([r["mean"] for r in g if r["resistant"]],
                                  [r["mean"] for r in g if not r["resistant"]])
        if pairs:
            out[k] = (g[0]["species"], g[0]["dataset"], half, pairs, len(g))
    return out


def mean(values):
    return sum(values, fractions.Fraction(0)) / len(values)


def levels(compared):
    """{level: {name: summary}}: groups equal inside a dataset, datasets inside a species, species overall. A
    summary holds the exact value and the rows, groups, pairs, datasets and species beneath it."""
    by_dataset, datasets_of = collections.defaultdict(list), collections.defaultdict(set)
    for species, dataset, half, pairs, rows in compared.values():
        by_dataset[dataset].append((half, pairs, rows))
        datasets_of[species].add(dataset)
    dataset_value = {}
    for d, gs in by_dataset.items():
        dataset_value[d] = mean([fractions.Fraction(h, 2 * p) for h, p, _ in gs])
    species_value = {sp: mean([dataset_value[d] for d in ds]) for sp, ds in datasets_of.items()}
    overall = mean([species_value[sp] for sp in species_value]) if species_value else None

    def count(datasets):
        gs = [g for d in datasets for g in by_dataset[d]]
        return {"rows": sum(g[2] for g in gs), "groups": len(gs), "pairs": sum(g[1] for g in gs),
                "datasets": len(datasets), "species": len({sp for sp, ds in datasets_of.items() if ds & set(datasets)})}

    out = {"overall": {"all": dict(count(set(by_dataset)), value=overall)}, "species": {}, "dataset": {}}
    for sp, ds in datasets_of.items():
        out["species"][sp] = dict(count(ds), value=species_value[sp])
    for d in by_dataset:
        out["dataset"][d] = dict(count({d}), value=dataset_value[d])
    return out


def interval(values):
    """(lower, upper) of the species means, given in sorted order of the species' names: 10,000 draws with
    replacement at the observed count, numpy PCG64 seed 42, a new generator for each estimate, and numpy's linear
    2.5th and 97.5th percentiles of the draws' means."""
    means = numpy.array([float(v) for v in values], dtype=numpy.float64)
    rng = numpy.random.Generator(numpy.random.PCG64(SEED))
    n = len(means)
    draws = means[rng.integers(0, n, size=(DRAWS, n))].mean(axis=1)
    lower, upper = numpy.percentile(draws, PERCENTILES)
    return float(lower), float(upper)


def summarise(compared):
    """The levels of one comparison, and its interval when 6 or more species enter."""
    lv = levels(compared)
    names = sorted(lv["species"])
    bounds = interval([lv["species"][sp]["value"] for sp in names]) if len(names) >= MIN_SPECIES else None
    return {"levels": lv, "interval": bounds, "species": len(names)}


def with_means(rows, scored):
    return [dict(r, mean=scored[(r["dataset"], r["genome"])][0], call=scored[(r["dataset"], r["genome"])][1])
            for r in rows]


def cells(groups, rule):
    """{cell: compared} for the four cells of catalogue calls, resistant member first, and the pairs counted
    apart: a member's unevaluable call, or a report with no drug rule."""
    out = {c: {} for c in CELLS + APART}
    for k, g in groups.items():
        rs = [r for r in g if r["resistant"]]
        ss = [r for r in g if not r["resistant"]]
        tally = collections.defaultdict(lambda: [0, 0, set()])
        for r in rs:
            for s in ss:
                if rule[r["dataset"]] is None:
                    c = "no drug rule"
                elif r["call"] is None or s["call"] is None:
                    c = "unevaluable"
                else:
                    c = (r["call"], s["call"])
                t = tally[c]
                t[0] += points(r["mean"], s["mean"])
                t[1] += 1
                t[2] |= {r["genome"], s["genome"]}
        for c, (half, pairs, members) in tally.items():
            out[c][k] = (g[0]["species"], g[0]["dataset"], half, pairs, len(members))
    return out


def typing(rows, hc20):
    """{step: {match: summary}}: each narrower match on the rows with its typing values, and the wider matches on
    the same surviving rows."""
    out = {}
    for step, extra, matches in TYPING:
        typed = [dict(r, hc20=hc20.get(r["genome"])) for r in rows]
        typed = [r for r in typed if all(r[k] is not None for k in extra)]
        narrow, _ = stage1.groups_of(typed, extra)
        surviving = [r for g in narrow.values() for r in g]
        out[step] = {match: summarise(compare(stage1.groups_of(surviving, keys)[0])) for match, keys in matches}
    return out


def null_fits(groups, complexity):
    """The primary comparison without every group whose fold has a null fit (fold_complexity 0) in any of the
    five learners of its report. Entry f of a learner's fold_complexity is the fit of fold f
    (src/learn/eval.rs, pool_folds)."""
    kept = {k: g for k, g in groups.items()
            if not any(complexity[g[0]["dataset"]][n][k[1]] == 0 for n in FIVE)}
    return summarise(compare(kept))


def measure(checked, opened):
    """Both readings, after the line: the primary measure and each secondary read."""
    scored, complexity, rule = {}, {}, {}
    for d, (path, digest, listed) in checked["listed"].items():
        got = read_scores(path, digest, listed, opened)
        for genome, (five, call) in got["rows"].items():
            scored[(d, genome)] = (row_mean(five), call)
        complexity[d] = got["fold_complexity"]
        rule[d] = got["catalog"]
    out = {}
    for name in READINGS:
        rows = with_means(checked["rows"][name], scored)
        groups = stage1.groups_of(rows)[0]
        primary = summarise(compare(groups))
        out[name] = {"primary": primary, "cells": {c: levels(v) for c, v in cells(groups, rule).items()},
                     "typing": typing(rows, checked["hc20"]), "null_fits": null_fits(groups, complexity)}
    return out


# 5. The export.

def ppm(value):
    """An exact fraction as an integer of parts per million, rounded half to even."""
    return round(fractions.Fraction(value) * PPM)


def estimate(level, s):
    """The estimate a row gives: none without pairs, and none for a dataset or species under 10 pairs."""
    if s is None or s["value"] is None or (level != "overall" and s["pairs"] < MIN_PAIRS):
        return ""
    return ppm(s["value"])


ZERO = {"rows": 0, "groups": 0, "pairs": 0, "datasets": 0, "species": 0, "value": None}


def level_rows(lv, names, lost_from=None):
    """One row a level and name, each value beside its counts; with lost_from, what the primary loses."""
    out = []
    for level in LEVELS:
        for n in sorted(names[level]):
            s = lv[level].get(n, ZERO)
            row = {"level": level, "name": n, "estimate_ppm": estimate(level, s)}
            row.update({c: s[c] for c in ("rows", "groups", "pairs", "datasets", "species")})
            if lost_from is not None:
                p = lost_from[level].get(n, ZERO)
                row.update({"lost_" + c: p[c] - s[c] for c in ("rows", "groups", "pairs", "datasets", "species")})
            out.append(row)
    return out


def chance(bounds):
    """The word beside an interval, read from its two bounds before they are rounded to parts per million: above
    when the lower bound is above one half, below when the upper bound is below one half, else unresolved. An
    estimate with no interval has no word."""
    if bounds is None:
        return ""
    lower, upper = bounds
    return "above" if lower > 0.5 else "below" if upper < 0.5 else "unresolved"


def bounds_row(name, read, match, summary):
    b = summary["interval"]
    return {"reading": name, "read": read, "match": match, "lower_ppm": "" if b is None else ppm(b[0]),
            "upper_ppm": "" if b is None else ppm(b[1]), "chance": chance(b), "species": summary["species"],
            "draws": 0 if b is None else DRAWS}


def numpy_parts():
    return [int(x) for x in numpy.__version__.split(".")]


def export_files(checked, result):
    """{file name: bytes}: integers and names alone, no decimal point, no genome id, pair, group or BioProject."""
    measure_rows, interval_rows, cell_rows, typing_rows, null_rows = [], [], [], [], []
    for name in READINGS:
        r = result[name]
        primary = r["primary"]["levels"]
        names = {level: set(primary[level]) for level in LEVELS}
        measure_rows += [dict(x, reading=name) for x in level_rows(primary, names)]
        interval_rows.append(bounds_row(name, "primary", "primary", r["primary"]))
        for c in CELLS + APART:
            s = r["cells"][c]["overall"]["all"] if r["cells"][c]["overall"]["all"]["groups"] else ZERO
            first, second = (CALLS[c[0]], CALLS[c[1]]) if c in CELLS else (c, c)
            cell_rows.append({"reading": name, "resistant_member_call": first, "susceptible_member_call": second,
                              "estimate_ppm": "" if c in APART else estimate("overall", s),
                              **{k: s[k] for k in ("pairs", "groups", "datasets", "species")}})
        for step, matches in r["typing"].items():
            for match, summary in matches.items():
                rows = level_rows(summary["levels"], names, primary if match == "primary" else None)
                typing_rows += [dict(x, reading=name, step=step, match=match) for x in rows]
                interval_rows.append(bounds_row(name, step, match, summary))
        null_rows += [dict(x, reading=name) for x in level_rows(r["null_fits"]["levels"], names, primary)]
        interval_rows.append(bounds_row(name, "null_fits", "primary", r["null_fits"]))
    sums, facts = checked["sums"], checked["facts"]
    major, minor, micro = numpy_parts()
    inputs = [
        {"input": "reports", "sha256": sums["reports"], "count": facts["reports"],
         "note": f"each report equals its line of the release's sums list; the sha256 is that list's; listed in "
                 f"a row list: {facts['listed']}"},
        {"input": "genome rows of the reports that enter", "sha256": "-", "count": facts["entering_rows"],
         "note": f"rows of the {facts['entering']} reports holding the five learners under the lineage split, counted "
                 f"before the line; each figure's rows are within them"},
        {"input": "species of the reports that enter", "sha256": "-", "count": facts["entering_species"],
         "note": "counted before the line; each figure's species are within them"},
        {"input": "row list, reading A", "sha256": sums["A"], "count": facts["lines"]["A"],
         "note": "lines; equal to the sha256 stage 1 printed"},
        {"input": "row list, reading B", "sha256": sums["B"], "count": facts["lines"]["B"],
         "note": "lines; equal to the sha256 stage 1 printed"},
        {"input": "count table of stage 1", "sha256": sums["counts"], "count": facts["table_lines"],
         "note": "lines of the committed copy; the groups rebuilt from each row list equal it at every level"},
        {"input": "typing table of stage 1", "sha256": sums["typing"], "count": facts["typing_lines"],
         "note": "lines of the committed copy; the typing support rebuilt from each row list equals it at every "
                 "level, but for the rows lacking a value among the rows that enter the groups"},
        {"input": "genome table, 5 Nov", "sha256": sums["hc20"], "count": facts["hc20_rows"],
         "note": f"rows; it gives the HC20 group; listed genomes it lacks: {facts['hc20_missing']}"},
        {"input": "numpy version, major", "sha256": "-", "count": major, "note": "the resampling's library"},
        {"input": "numpy version, minor", "sha256": "-", "count": minor, "note": "the resampling's library"},
        {"input": "numpy version, micro", "sha256": "-", "count": micro, "note": "the resampling's library"},
    ]
    out = {
        "measure.tsv": stage1.tsv(MEASURE, measure_rows),
        "interval.tsv": stage1.tsv(INTERVAL, interval_rows),
        "cells.tsv": stage1.tsv(CELL_COLUMNS, cell_rows),
        "typing.tsv": stage1.tsv(TYPING_COLUMNS, typing_rows),
        "null_fits.tsv": stage1.tsv(NULL_COLUMNS, null_rows),
        "inputs.tsv": stage1.tsv(("input", "sha256", "count", "note"), inputs) + (OPENED + "\n").encode("utf-8"),
    }
    return out


def scan_export(files, ids, accessions):
    """The export's scan, before any file is written: a full stop in any file stops the run, then [(file name,
    hits)] of stage 1's scan for every genome id and accession read, after its positive control. Names were
    checked before the line, so only a fault of this script can put a full stop here."""
    for file, data in files.items():
        if b"." in data:
            raise Stop(f"{file} holds a decimal point")
    return stage1.scan(files, ids, accessions)


def write_export(out, files):
    """Write <out>/export. On an OSError, remove every file it made, then stop."""
    made = []
    try:
        export = Path(out) / "export"
        export.mkdir()
        made.append(export)
        for name, data in files.items():
            fd = os.open(export / name, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
            made.append(export / name)
            with os.fdopen(fd, "wb") as f:
                f.write(data)
    except OSError as e:
        left = 0
        for m in reversed(made):
            try:
                m.rmdir() if m.is_dir() else m.unlink()
            except OSError:
                left += 1
        raise Stop(f"the export could not be written ({type(e).__name__}); of the {len(made)} paths it made, "
                   f"{left} could not be removed") from None


def said(result):
    lines = []
    for name in READINGS:
        p = result[name]["primary"]
        o = p["levels"]["overall"]["all"]
        b = p["interval"]
        tail = f"interval {ppm(b[0])} to {ppm(b[1])} ppm over {p['species']} species and {DRAWS} draws" \
            if b else f"no interval: {p['species']} species, under {MIN_SPECIES}"
        lines.append(f"reading {name}: estimate {ppm(o['value'])} ppm over {o['groups']} groups, {o['pairs']} pairs, "
                     f"{o['datasets']} datasets and {o['species']} species; {tail}")
    return lines


def run(args):
    """(export files, ids and accessions read, printout). The line comes after check() and before any score."""
    checked = check(args)
    opened = open_scores(args.out)
    result = measure(checked, opened)
    return export_files(checked, result), checked["ids"], checked["accessions"], said(result)


def parse(argv):
    p = argparse.ArgumentParser(prog="stage2.py", description="Stage 2 of the matched accuracy protocol.")
    p.add_argument("--reports", required=True)
    p.add_argument("--rows-a", required=True)
    p.add_argument("--rows-b", required=True)
    p.add_argument("--hc20", required=True)
    p.add_argument("--out", required=True)
    return p.parse_args(argv)


def main(argv):
    args = parse(argv)
    try:
        files, ids, accessions, lines = run(args)
        found = scan_export(files, ids, accessions)
        if not found:
            write_export(args.out, files)
    except Stop as e:
        print(f"stage2: stopped: {e}")
        return 1
    if found:
        for name, n in found:
            print(f"{name} {n}")
        return 1
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
