"""year_mix.py: the mix of catalogue strata by collection year, from the
release's reports alone (docs/mechanism-year-result.md).

    python3 -I scripts/mechanism/year_mix.py plan
    python3 -I scripts/mechanism/year_mix.py read --reports DIR --lists LISTS --published STRATA --out READING
    python3 -I scripts/mechanism/year_mix.py write --reading READING --reports DIR --out OUT
    python3 -I scripts/mechanism/year_mix.py scan --reports DIR PATH ...

`read` reads every report under DIR once. It checks guards 1, 2, 3 and 5 of
the plan, rebuilds M1's four strata from GRM's known-gene rule (each genome's
`groups` against `catalog.groups[].counts_for_rule`) and holds the rebuild to
the eight published counts of STRATA (`page-m1-strata.tsv`). If all eight
are equal the year read runs on the four strata; if not, on the reader's call
and no call (`samples[].catalog_resistant`). It writes READING, counts per
dataset, year and stratum, which stays where the reports are. `write` reads
READING, checks guard 4 and writes the tables of OUT. Both end with guard 6:
no genome id and no catalogue group name of the reports stands in a written
file. Every message prints counts, species and antibiotics only; no genome id
and no group name is printed, hashed or written.

It fits nothing, scores no genome and runs no program: every figure comes
from the scores, labels and calls the reports hold.

Python 3.12, standard library only; run with -I.
"""

import argparse
import csv
import glob
import hashlib
import json
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from common import r12, recall_cell  # noqa: E402

LEARNERS = ("scm", "cart", "rf", "gbt", "lasso")
SPLIT = "lineage"
CUTOFF = 0.5
THIN = 100            # a year with fewer genomes is thin
SPECIES_FLOOR = 10    # species.tsv keeps a species, antibiotic and stratum with 10 or more resistant genomes
UNKNOWN = "year_unknown"
YEAR_RE = re.compile(r"[0-9]{4}")
FOUR = ("point_only", "gene_only", "both", "neither")
TWO = ("reader_call", "reader_no_call", "reader_no_rule")
ALL = "all_strata"
READING_SCHEMA = "mechanism-year-reading/1"

STEP1_HEADER = ["stratum", "published_resistant", "published_susceptible", "rebuilt_resistant",
                "rebuilt_susceptible", "equal"]
COUNTS = ["genomes", "resistant", "susceptible", "resistant_called_mean5", "susceptible_called_mean5"]
YEARS_HEADER = ["year", "stratum"] + COUNTS + ["recall_mean5", "thin"]
TOTALS_HEADER = ["stratum"] + COUNTS + ["recall_mean5"]
SHARES_HEADER = ["year", "stratum", "resistant", "year_resistant", "resistant_share", "susceptible",
                 "year_susceptible", "susceptible_share", "year_genomes", "thin"]
SPECIES_HEADER = ["species", "antibiotic", "stratum"] + COUNTS + ["recall_mean5", "year_unknown"]
TABLES = ("step1.tsv", "years.tsv", "totals.tsv", "shares.tsv", "species.tsv", "guards.txt")
NOT_USED = ("in_M", "not_in_A", "not_in_list", "no_lineage_split")


class Stop(Exception):
    """A guard failed; nothing more is written."""


# ---------------------------------------------------------------- small rules

def year_key(text):
    """A year from samples[].year: four ASCII digits and nothing else, not even
    a blank; anything else, a missing year included, is year_unknown."""
    if isinstance(text, str) and YEAR_RE.fullmatch(text):
        return text
    return UNKNOWN


def mean_called(scores):
    """The five learners' mean vote: resistant at a mean of 0.5 or more. The
    sum is rounded once (math.fsum), so scores whose decimals add to 2.5,
    such as 0, 0.6, 0.6, 0.6 and 0.7, are called whatever their order."""
    return math.fsum(scores) >= CUTOFF * len(LEARNERS)


def genome_scores(rows, i):
    """Genome i's held-out scores, one per learner, in LEARNERS' order."""
    return [rows[name]["scores"][i] for name in LEARNERS]


def set_of(lrow):
    """'primary', or the reason a dataset is not used."""
    if lrow is None:
        return "not_in_list"
    if lrow["in_A"] != "1":
        return "not_in_A"
    if lrow["in_M"] != "0":
        return "in_M"
    return "primary"


def is_thin(genomes):
    return genomes < THIN


def stratum_of(point, gene):
    if point and gene:
        return "both"
    if point:
        return "point_only"
    if gene:
        return "gene_only"
    return "neither"


def point_flag(sample, report):
    """No field of a report marks a catalogue group as a point mutation: a
    catalog.groups row holds class, subclass, counts_for_rule, group and two
    counts, and GRM's reader holds acquired genes only (src/catalog.rs)."""
    return False


def counting_groups(report):
    """Indices into catalog_groups of the groups whose catalog.groups row
    counts for the dataset's rule."""
    catalog = report.get("catalog") or {}
    counts = {g.get("group"): g.get("counts_for_rule") is True for g in catalog.get("groups") or []}
    return {i for i, name in enumerate(report.get("catalog_groups") or []) if counts.get(name, False)}


def rebuilt_stratum(sample, report, counting):
    gene = any(g in counting for g in sample.get("groups") or [])
    return stratum_of(point_flag(sample, report), gene)


def reader_group(sample):
    call = sample.get("catalog_resistant")
    if call is True:
        return "reader_call"
    if call is False:
        return "reader_no_call"
    return "reader_no_rule"


def label_cells(samples):
    """tp, fp, tn, fn of the reader's call against the label, and the genomes
    whose call is null."""
    tp = fp = tn = fn = null = 0
    for s in samples:
        call = s.get("catalog_resistant")
        if call is None:
            null += 1
        elif s["resistant"]:
            tp, fn = (tp + 1, fn) if call else (tp, fn + 1)
        else:
            fp, tn = (fp + 1, tn) if call else (fp, tn + 1)
    return tp, fp, tn, fn, null


def cell_disagreement(report):
    """None when the four cells counted from samples equal catalog.metrics, or
    the report has no catalogue row and every call is null; else the reason."""
    tp, fp, tn, fn, null = label_cells(report["samples"])
    catalog = report.get("catalog")
    if catalog is None:
        return None if null == len(report["samples"]) else "calls but no catalogue row"
    metrics = catalog.get("metrics") or {}
    want = tuple(metrics.get(k) for k in ("tp", "fp", "tn", "fn"))
    if null or (tp, fp, tn, fn) != want:
        return f"cells {tp} {fp} {tn} {fn} null {null}, metrics {' '.join(str(x) for x in want)}"
    return None


def lineage_rows(report):
    """learner -> row of the lineage split, or None with no lineage split."""
    found = [s for s in report.get("splits") or [] if s.get("scheme") == SPLIT]
    if not found:
        return None
    if len(found) > 1:
        raise Stop(f"guard 3: {report['species']} with {report['antibiotic']} has {len(found)} lineage splits")
    return {row.get("learner"): row for row in found[0].get("learners") or []}


def a_score(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x) and 0 <= x <= 1


def check_scores(report, rows):
    """Guard 3: each learner has a row, scores as long as samples, each a number from 0 to 1."""
    name = f"{report['species']} with {report['antibiotic']}"
    n = len(report["samples"])
    for learner in LEARNERS:
        if learner not in rows:
            raise Stop(f"guard 3: {name} has no {learner} row on the lineage split")
        scores = rows[learner].get("scores")
        if not isinstance(scores, list) or len(scores) != n:
            got = len(scores) if isinstance(scores, list) else "no"
            raise Stop(f"guard 3: {name}, {learner}: {got} scores for {n} samples")
        bad = sum(not a_score(x) for x in scores)
        if bad:
            raise Stop(f"guard 3: {name}, {learner}: {bad} scores are not a number from 0 to 1")


def strata_match(published, rebuilt):
    return all(tuple(published[s]) == tuple(rebuilt[s]) for s in FOUR)


def read_tsv(path):
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def published_strata(path):
    """The eight counts of page-m1-strata.tsv's primary lineage stratum rows,
    after checking they are the same for every row (learner) of the file."""
    by_row = {}
    for r in read_tsv(path):
        if (r["set"], r["split"], r["kind"]) == ("primary", SPLIT, "stratum") and r["stratum"] in FOUR:
            by_row.setdefault(r["row"], {})[r["stratum"]] = (int(r["resistant"]), int(r["susceptible"]))
    if not by_row:
        raise Stop("step 1: the published file holds no primary lineage stratum row")
    forms = {tuple(sorted(v.items())) for v in by_row.values()}
    if len(forms) != 1:
        raise Stop(f"step 1: the published stratum rows differ between the file's {len(by_row)} rows")
    one = next(iter(by_row.values()))
    if set(one) != set(FOUR):
        raise Stop(f"step 1: the published file gives {len(one)} of the 4 strata")
    return one, len(by_row)


# ---------------------------------------------------------------- read

def blank():
    return [0, 0, 0, 0, 0]


def add(cell, resistant, called):
    cell[0] += 1
    if resistant:
        cell[1] += 1
        cell[3] += called
    else:
        cell[2] += 1
        cell[4] += called


def read(reports_dir, lists_path, published_path):
    files = sorted(glob.glob(os.path.join(reports_dir, "*.json")))
    if not files:
        raise Stop(f"guard 1: no report under {reports_dir}")
    lists = {r["dataset"]: r for r in read_tsv(lists_path)}
    published, published_rows = published_strata(published_path)
    failed, disagree, names_differ = [], [], []
    not_used = {k: 0 for k in NOT_USED}
    checked = no_catalog = used = learner_rows = 0
    seen = set()
    rebuilt = {s: [0, 0] for s in FOUR}
    cells = {"four": {}, "two": {}}
    totals = {"four": {s: blank() for s in FOUR}, "two": {s: blank() for s in TWO}}
    ids, names = set(), set()
    for path in files:
        stem = os.path.basename(path)[:-len(".json")]
        try:
            with open(path, encoding="utf-8") as handle:
                report = json.load(handle)
            samples = report["samples"]
        except (OSError, ValueError, KeyError, TypeError):
            failed.append(stem)
            continue
        collect_names(report, ids, names)
        if report.get("catalog") is None:
            no_catalog += 1
        else:
            checked += 1
        why = cell_disagreement(report)
        if why is not None:
            disagree.append(f"{report.get('species')} with {report.get('antibiotic')}")
        lrow = lists.get(stem)
        kind = set_of(lrow)
        if kind != "primary":
            not_used[kind] += 1
            continue
        seen.add(stem)
        if (report.get("species"), report.get("antibiotic")) != (lrow["species"], lrow["antibiotic"]):
            names_differ.append(stem)
        rows = lineage_rows(report)
        if rows is None:
            not_used["no_lineage_split"] += 1
            continue
        check_scores(report, rows)
        learner_rows += len(LEARNERS)
        used += 1
        counting = counting_groups(report)
        species, antibiotic = report["species"], report["antibiotic"]
        for i, s in enumerate(samples):
            resistant = s["resistant"] is True
            called = int(mean_called(genome_scores(rows, i)))
            year = year_key(s.get("year"))
            four = rebuilt_stratum(s, report, counting)
            two = reader_group(s)
            rebuilt[four][0 if resistant else 1] += 1
            for mode, stratum in (("four", four), ("two", two)):
                add(cells[mode].setdefault((species, antibiotic, year, stratum), blank()), resistant, called)
                add(totals[mode][stratum], resistant, called)
    if failed:
        raise Stop(f"guard 1: {len(failed)} of {len(files)} reports do not parse: {', '.join(failed)}")
    if disagree:
        raise Stop(f"guard 2: {len(disagree)} reports' cells differ from catalog.metrics: {'; '.join(sorted(disagree))}")
    missing = sorted(k for k, r in lists.items() if set_of(r) == "primary" and k not in seen)
    if missing:
        raise Stop(f"guard 5: {len(missing)} primary datasets of the list have no report: {', '.join(missing)}")
    if names_differ:
        raise Stop(f"guard 5: {len(names_differ)} reports name another species or antibiotic than the list")
    if used + sum(not_used.values()) != len(files):
        raise Stop(f"guard 5: used {used} plus not used {sum(not_used.values())} is not {len(files)}")
    rebuilt = {s: tuple(v) for s, v in rebuilt.items()}
    mode = "four" if strata_match(published, rebuilt) else "two"
    step1 = []
    for s in FOUR + (ALL,):
        p = published[s] if s != ALL else tuple(sum(published[x][j] for x in FOUR) for j in (0, 1))
        b = rebuilt[s] if s != ALL else tuple(sum(rebuilt[x][j] for x in FOUR) for j in (0, 1))
        step1.append({"stratum": s, "published_resistant": p[0], "published_susceptible": p[1],
                      "rebuilt_resistant": b[0], "rebuilt_susceptible": b[1], "equal": "yes" if p == b else "no"})
    strata = FOUR if mode == "four" else TWO
    reading = {
        "schema": READING_SCHEMA,
        "mode": mode,
        "strata": list(strata),
        "reports_found": len(files),
        "reports_read": len(files) - len(failed),
        "used": used,
        "not_used": not_used,
        "catalogue_rows_checked": checked,
        "no_catalogue_row": no_catalog,
        "learner_rows": learner_rows,
        "published_rows": published_rows,
        "step1": step1,
        "totals": {s: dict(zip(COUNTS, totals[mode][s])) for s in strata},
        "cells": [dict(zip(["species", "antibiotic", "year", "stratum"], k), **dict(zip(COUNTS, v)))
                  for k, v in sorted(cells[mode].items())],
    }
    return reading, ids, names


# ---------------------------------------------------------------- guard 6

def collect_names(report, ids, names):
    """The genome ids and catalogue group names of one report, into two sets
    held in memory only; never printed or written."""
    for s in report.get("samples") or []:
        if isinstance(s.get("id"), str) and s["id"]:
            ids.add(s["id"])
    for name in report.get("catalog_groups") or []:
        if isinstance(name, str) and name:
            names.add(name)
    for g in (report.get("catalog") or {}).get("groups") or []:
        if isinstance(g.get("group"), str) and g["group"]:
            names.add(g["group"])


def names_from(reports_dir):
    ids, names = set(), set()
    files = sorted(glob.glob(os.path.join(reports_dir, "*.json")))
    for path in files:
        with open(path, encoding="utf-8") as handle:
            collect_names(json.load(handle), ids, names)
    return ids, names, len(files)


TOKEN = re.compile(r"[A-Za-z0-9_.\-]+")


def scan_hits(paths, ids, names):
    """Files holding a genome id (a whole token) and files holding a group name
    (with no letter or digit beside it). Counts only."""
    pattern = None
    if names:
        alts = "|".join(re.escape(n) for n in sorted(names, key=lambda n: (-len(n), n)))
        pattern = re.compile(rf"(?<![A-Za-z0-9])(?:{alts})(?![A-Za-z0-9])")
    id_files = name_files = 0
    for path in paths:
        with open(path, encoding="utf-8") as handle:
            text = handle.read()
        if ids & set(TOKEN.findall(text)):
            id_files += 1
        if pattern is not None and pattern.search(text):
            name_files += 1
    return id_files, name_files


def guard6(paths, ids, names, where):
    id_files, name_files = scan_hits(paths, ids, names)
    line = (f"guard 6 ({where}): {len(paths)} files, {id_files} with a genome id, {name_files} with a catalogue "
            f"group name, against {len(ids)} ids and {len(names)} group names")
    return line, id_files + name_files == 0


# ---------------------------------------------------------------- write

def year_order(year):
    return (1, "") if year == UNKNOWN else (0, year)


def share(n, d):
    return f"{r12(n / d):.12f}" if d else ""


def recall_of(c):
    return recall_cell(c["resistant_called_mean5"], c["resistant"])


def check_sums(by_year, totals, strata):
    """Guard 4: each stratum's counts summed over the years equal its total."""
    for s in strata:
        for k in COUNTS:
            summed = sum(c[k] for (y, x), c in by_year.items() if x == s)
            if summed != totals[s][k]:
                raise Stop(f"guard 4: {s} {k} sums to {summed} over the years, its total is {totals[s][k]}")


def check_reading(reading):
    if reading.get("schema") != READING_SCHEMA:
        raise Stop("write: the reading is not year_mix.py read's")
    step = {r["stratum"]: r for r in reading["step1"]}
    equal = all(step[s]["equal"] == "yes" for s in FOUR)
    if reading["mode"] != ("four" if equal else "two") or reading["strata"] != list(FOUR if equal else TWO):
        raise Stop("write: the reading's strata do not follow from step 1")
    for r in reading["step1"]:
        same = (r["published_resistant"], r["published_susceptible"]) == (r["rebuilt_resistant"], r["rebuilt_susceptible"])
        if r["equal"] != ("yes" if same else "no"):
            raise Stop(f"write: step 1's {r['stratum']} row says equal {r['equal']}")


def tables(reading):
    check_reading(reading)
    strata = reading["strata"]
    totals = reading["totals"]
    by_year = {}
    by_species = {}
    for c in reading["cells"]:
        if c["stratum"] not in strata:
            raise Stop(f"write: a cell of stratum {c['stratum']} outside the reading's strata")
        y = by_year.setdefault((c["year"], c["stratum"]), dict.fromkeys(COUNTS, 0))
        sp = by_species.setdefault((c["species"], c["antibiotic"], c["stratum"]),
                                   dict(dict.fromkeys(COUNTS, 0), year_unknown=0))
        for k in COUNTS:
            y[k] += c[k]
            sp[k] += c[k]
        if c["year"] == UNKNOWN:
            sp["year_unknown"] += c["genomes"]
    years = sorted({y for y, _ in by_year}, key=year_order)
    for y in years:
        for s in strata:
            by_year.setdefault((y, s), dict.fromkeys(COUNTS, 0))
    check_sums(by_year, totals, strata)
    out = {"step1.tsv": [STEP1_HEADER] + [[str(r[h]) for h in STEP1_HEADER] for r in reading["step1"]]}
    year_rows, share_rows = [YEARS_HEADER], [SHARES_HEADER]
    for y in years:
        yr = sum(by_year[(y, s)]["resistant"] for s in strata)
        ys = sum(by_year[(y, s)]["susceptible"] for s in strata)
        yg = sum(by_year[(y, s)]["genomes"] for s in strata)
        thin = "yes" if is_thin(yg) else "no"
        for s in strata:
            c = by_year[(y, s)]
            year_rows.append([y, s] + [str(c[k]) for k in COUNTS] + [recall_of(c), thin])
            share_rows.append([y, s, str(c["resistant"]), str(yr), share(c["resistant"], yr), str(c["susceptible"]),
                               str(ys), share(c["susceptible"], ys), str(yg), thin])
    out["years.tsv"] = year_rows
    out["totals.tsv"] = [TOTALS_HEADER] + [[s] + [str(totals[s][k]) for k in COUNTS] + [recall_of(totals[s])]
                                           for s in strata]
    out["shares.tsv"] = share_rows
    order = {s: i for i, s in enumerate(strata)}
    out["species.tsv"] = [SPECIES_HEADER] + [
        [sp, ab, s] + [str(c[k]) for k in COUNTS] + [recall_of(c), str(c["year_unknown"])]
        for (sp, ab, s), c in sorted(by_species.items(), key=lambda kv: (kv[0][0], kv[0][1], order[kv[0][2]]))
        if c["resistant"] >= SPECIES_FLOOR]
    return out, years


def guard_lines(reading, years):
    nu = reading["not_used"]
    thin = sum(1 for y in years if is_thin(sum(c["genomes"] for c in reading["cells"] if c["year"] == y)))
    return [
        f"guard 1: {reading['reports_found']} reports found, {reading['reports_read']} read",
        f"guard 2: {reading['catalogue_rows_checked']} reports with a catalogue row, their cells equal catalog.metrics;"
        f" {reading['no_catalogue_row']} with no catalogue row and no call",
        f"guard 3: {reading['used']} datasets, {reading['learner_rows']} learner rows on the lineage split,"
        f" each with scores of the samples' length, each a number from 0 to 1",
        f"guard 4: {len(reading['strata'])} strata, each summed over {len(years)} year rows"
        f" ({thin} thin, {UNKNOWN} included where present) equals its total",
        f"guard 5: {reading['used']} used + " + " + ".join(f"{nu[k]} {k}" for k in NOT_USED)
        + f" = {reading['reports_found']} reports",
        f"step 1: {reading['published_rows']} published rows the same; strata {reading['mode']}: "
        + ", ".join(reading["strata"]),
    ]


def sha256_of(path):
    with open(path, "rb") as handle:
        return hashlib.sha256(handle.read()).hexdigest()


def write_tsv(path, rows):
    with open(path, "w", encoding="utf-8", newline="") as handle:
        handle.write("".join("\t".join(r) + "\n" for r in rows))


def write(reading, out_dir, ids, names):
    if os.path.isdir(out_dir) and os.listdir(out_dir):
        raise Stop(f"write: {out_dir} is not empty")
    out, years = tables(reading)
    lines = guard_lines(reading, years)
    os.makedirs(out_dir, exist_ok=True)
    paths = []
    for name in TABLES[:-1]:
        paths.append(os.path.join(out_dir, name))
        write_tsv(paths[-1], out[name])
    paths.append(os.path.join(out_dir, "guards.txt"))
    with open(paths[-1], "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines) + "\n")
    line6, clean = guard6(paths, ids, names, "tables")
    if not clean:
        for p in paths:
            os.remove(p)
        raise Stop(line6)
    with open(paths[-1], "a", encoding="utf-8") as handle:
        handle.write(line6 + "\n")
    sums = "".join(f"{sha256_of(p)}  {os.path.basename(p)}\n" for p in paths)
    with open(os.path.join(out_dir, "SHA256SUMS"), "w", encoding="utf-8") as handle:
        handle.write(sums)
    return lines + [line6]


# ---------------------------------------------------------------- the commands

def plan():
    return {
        "learners": list(LEARNERS), "split": SPLIT, "cutoff": CUTOFF, "call": "mean of the five scores >= cutoff",
        "year": "four ASCII digits, else " + UNKNOWN, "thin_under_genomes": THIN,
        "species_floor_resistant": SPECIES_FLOOR, "recall_floor_resistant": 10,
        "primary": "lists.tsv in_A 1 and in_M 0", "strata_four": list(FOUR), "strata_two": list(TWO),
        "columns": {"step1.tsv": STEP1_HEADER, "years.tsv": YEARS_HEADER, "totals.tsv": TOTALS_HEADER,
                    "shares.tsv": SHARES_HEADER, "species.tsv": SPECIES_HEADER},
        "not_used": list(NOT_USED),
    }


def step1_lines(reading):
    lines = [f"step 1: {r['stratum']} published {r['published_resistant']} resistant {r['published_susceptible']}"
             f" susceptible; rebuilt {r['rebuilt_resistant']} resistant {r['rebuilt_susceptible']} susceptible;"
             f" equal {r['equal']}" for r in reading["step1"]]
    if reading["mode"] == "two":
        lines.append("step 1: the rebuild differs from the published counts; the year read runs on the reader's call"
                     " and no call")
    return lines


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("plan")
    r = sub.add_parser("read")
    r.add_argument("--reports", required=True)
    r.add_argument("--lists", required=True)
    r.add_argument("--published", required=True)
    r.add_argument("--out", required=True)
    w = sub.add_parser("write")
    w.add_argument("--reading", required=True)
    w.add_argument("--reports", required=True)
    w.add_argument("--out", required=True)
    s = sub.add_parser("scan")
    s.add_argument("--reports", required=True)
    s.add_argument("paths", nargs="+")
    args = ap.parse_args(argv)
    try:
        if args.cmd == "plan":
            print(json.dumps(plan(), indent=1))
            return 0
        if args.cmd == "read":
            if os.path.exists(args.out):
                raise Stop(f"read: {args.out} exists")
            reading, ids, names = read(args.reports, args.lists, args.published)
            with open(args.out, "w", encoding="utf-8") as handle:
                json.dump(reading, handle, indent=1, sort_keys=True)
                handle.write("\n")
            line6, clean = guard6([args.out], ids, names, "reading")
            if not clean:
                os.remove(args.out)
                raise Stop(line6)
            for line in guard_lines_read(reading) + step1_lines(reading) + [line6]:
                print(line)
            print(f"read: {len(reading['cells'])} cells, strata {reading['mode']}")
            return 0
        if args.cmd == "write":
            with open(args.reading, encoding="utf-8") as handle:
                reading = json.load(handle)
            ids, names, _ = names_from(args.reports)
            for line in write(reading, args.out, ids, names):
                print(line)
            print(f"write: {len(TABLES)} tables and SHA256SUMS in {args.out}")
            return 0
        ids, names, n = names_from(args.reports)
        line6, clean = guard6(args.paths, ids, names, f"scan over {n} reports")
        print(line6)
        return 0 if clean else 1
    except Stop as e:
        print(f"stop: {e}")
        return 2


def guard_lines_read(reading):
    nu = reading["not_used"]
    return [
        f"guard 1: {reading['reports_found']} reports found, {reading['reports_read']} read",
        f"guard 2: {reading['catalogue_rows_checked']} reports with a catalogue row, cells equal catalog.metrics;"
        f" {reading['no_catalogue_row']} with no catalogue row and no call",
        f"guard 3: {reading['used']} datasets, {reading['learner_rows']} learner rows",
        f"guard 5: {reading['used']} used + " + " + ".join(f"{nu[k]} {k}" for k in NOT_USED)
        + f" = {reading['reports_found']} reports",
    ]


if __name__ == "__main__":
    sys.exit(main())
