"""curve_tables.py: the tables of docs/lineage-curve-result.md, block M4 of the
mechanism split (the lineage curve for Mycobacterium tuberculosis) read once,
as files and as the page prints them, and the checks of the tool run.

    python3 -I curve_tables.py check RUN_DIR CHECKS
    python3 -I curve_tables.py write READING LISTS BED OUT
    python3 -I curve_tables.py print OUT

`check` reads run_tbprofiler.sh's output folder (plan.txt, failed.tsv,
results/<genome>.results.json, run.start, run.end) and writes CHECKS, counts
only: the genomes planned, the result files, the failed lines and their
reasons, the result files that name TBProfiler 6.6.3 and `tbdb` commit
3ed1cb99, and the wall time. Of a result file it reads `id` and the two
version fields alone. It prints no genome id and exits 1 when a check fails.

READING is the reading's directory:

    reading.txt   the reading's commands, each with the sha256 of what it wrote
    checks.txt    what `check` wrote over the tool run
    m4/           what tb_curve.py wrote: m4-curve.tsv, m4-coverage.tsv,
                  m4-convergence.tsv, m4-genes.tsv, m4-summary.json

LISTS is scripts/c3-inputs/lists.tsv; BED is the tool's tbdb.bed, whose
fifth column gives the gene names the WHO catalogue uses in that database.

tb_curve.py's own tables are the reading; this script does not compute them
again in another way. `write` passes every gate of GATES, stopping at the
first that fails and naming it, then copies m4/ under OUT as it is and writes
beside it page-endpoints.tsv (one row per report of the species, species and
antibiotic in full, set M marked), then gates.txt, README.md and SHA256SUMS.
OUT must be absent or empty.

`print` prints the page's tables from OUT alone, figures rounded to 3
decimals. Every row is an aggregate over an antibiotic, a lineage group or a
gene: no file holds a genome id, a BioProject accession, a variant's own
name, a sequence or a position.

Python 3.12, standard library only; run with -I.
"""

import csv
import hashlib
import json
import os
import shutil
import sys
from datetime import datetime
from decimal import ROUND_HALF_EVEN, Decimal
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import scan  # noqa: E402

WANT_TOOL = "6.6.3"
WANT_DB = "3ed1cb99"
FLOOR = 5
FOLDED = "other genes"
M4_FILES = ("m4-convergence.tsv", "m4-coverage.tsv", "m4-curve.tsv", "m4-genes.tsv", "m4-summary.json")
GROUPINGS = ("main_lineage", "sub_lineage", "cluster")
CHECK_KEYS = ("planned", "results", "failed_lines", "failed_genomes", "accounted", "neither", "both", "unplanned",
              "versions_named", "versions_other", "failed_by_reason", "start", "end", "wall_seconds")
GATES = [
    (1, "layout", "the reading holds reading.txt, checks.txt and tb_curve.py's five files under m4/, and nothing more "
                  "under m4/"),
    (2, "versions", "checks.txt counts every result file naming TBProfiler 6.6.3 and tbdb 3ed1cb99 and none naming "
                    "another, and m4-summary.json names the same two"),
    (3, "accounting", "result files plus failed lines equal the genomes planned, with no planned genome left with "
                      "neither, none with both and no result outside the plan; m4-summary.json gives the same counts"),
    (4, "datasets", "the summary's endpoints and reports not scored are exactly the species' reports in lists.tsv, "
                    "named as lists.tsv names them, and every row of every table carries in_M as lists.tsv gives it"),
    (5, "floors", "a share stands only over 5 or more resistant genomes, and no curve over under 5"),
    (6, "names", "no file names a variant, a position or a sequence: the gene column holds only gene names of "
                 "tbdb.bed or `other genes`"),
    (7, "genes", "a gene is named only where 5 or more resistant genomes of the antibiotic carry it"),
]
ENDPOINT_COLUMNS = ["dataset", "species", "antibiotic", "in_M", "labels_depend_on", "kind", "tool_drugs",
                    "resistant_with_result", "resistant_no_result", "catalogue_variants", "genes_named",
                    "without_determinant", "share_without", "groups_main", "groups_sub", "groups_cluster",
                    "no_main", "no_sub", "no_cluster", "one_lineage_main", "one_lineage_sub"]


class Stop(Exception):
    pass


# ---------- small pieces ----------

def load(path):
    with open(path, encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t", quoting=csv.QUOTE_NONE))


def write_tsv(path, rows, header):
    with open(path, "w", encoding="utf-8", newline="") as handle:
        handle.write("\t".join(header) + "\n")
        for r in rows:
            cells = ["" if r.get(h) is None else str(r.get(h)) for h in header]
            if any("\t" in c or "\n" in c for c in cells):
                raise Stop(f"{Path(path).name}: a cell holds a tab or a newline")
            handle.write("\t".join(cells) + "\n")


def sha256_of(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def is_int(text):
    return isinstance(text, str) and text.isascii() and text.isdigit()


def species_rows(lrows, slug):
    return [r for r in lrows if r["species_slug"] == slug]


def in_m(lrow):
    return lrow["in_M"] == "1"


def bed_genes(path):
    genes = set()
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            row = line.rstrip("\n").split("\t")
            if len(row) >= 6 and row[4]:
                genes.add(row[4])
    if not genes:
        raise Stop("tbdb.bed names no gene")
    return genes


# ---------- check: the tool run, counts only ----------

def names_pinned(d, gid):
    """A result file of this genome that names TBProfiler 6.6.3 and tbdb 3ed1cb99."""
    sw = d.get("pipeline", {}).get("software_version")
    db = (d.get("pipeline", {}).get("db_version") or {}).get("commit")
    return d.get("id") == gid and sw == WANT_TOOL and db == WANT_DB


def run_counts(run_dir):
    """The checks of a run_tbprofiler.sh folder as a dict of CHECK_KEYS. Genome ids are compared and never
    returned."""
    run = Path(run_dir)
    plan = [x.strip() for x in (run / "plan.txt").read_text(encoding="utf-8").splitlines() if x.strip()]
    planned = set(plan)
    failed_rows = []
    if (run / "failed.tsv").is_file():
        failed_rows = [x.split("\t") for x in (run / "failed.tsv").read_text(encoding="utf-8").splitlines() if x]
    failed = {r[0] for r in failed_rows}
    reasons = {}
    for r in failed_rows:
        reason = r[1] if len(r) > 1 and r[1] else "none given"
        reasons[reason] = reasons.get(reason, 0) + 1
    suffix = ".results.json"
    results = {p.name[: -len(suffix)]: p for p in (run / "results").iterdir() if p.name.endswith(suffix)}
    named = other = 0
    for gid, path in results.items():
        with open(path, encoding="utf-8") as handle:
            d = json.load(handle)
        if names_pinned(d, gid):
            named += 1
        else:
            other += 1
    start = (run / "run.start").read_text().strip() if (run / "run.start").is_file() else ""
    end = (run / "run.end").read_text().strip() if (run / "run.end").is_file() else ""
    wall = ""
    if start and end:
        f = "%Y-%m-%dT%H:%M:%SZ"
        wall = str(int((datetime.strptime(end, f) - datetime.strptime(start, f)).total_seconds()))
    in_plan = set(results) & planned
    return {"planned": len(plan), "results": len(results), "failed_lines": len(failed_rows),
            "failed_genomes": len(failed), "accounted": len(in_plan) + len(failed & planned - set(results)),
            "neither": len(planned - set(results) - failed), "both": len(set(results) & failed),
            "unplanned": len((set(results) | failed) - planned), "versions_named": named, "versions_other": other,
            "failed_by_reason": ",".join(f"{k}={reasons[k]}" for k in sorted(reasons)) or "none",
            "start": start, "end": end, "wall_seconds": wall}


def checks_pass(c):
    """The run's checks as (name, passes) pairs, from a dict of CHECK_KEYS (ints or their text)."""
    n = {k: int(c[k]) for k in ("planned", "results", "failed_lines", "failed_genomes", "accounted", "neither",
                                "both", "unplanned", "versions_named", "versions_other")}
    return [
        ("result files plus failed lines equal the genomes planned", n["results"] + n["failed_lines"] == n["planned"]),
        ("every planned genome has a result file or a failed line", n["neither"] == 0 and n["accounted"] == n["planned"]),
        ("no genome has both a result file and a failed line", n["both"] == 0),
        ("no result file or failed line outside the plan", n["unplanned"] == 0),
        ("no genome has two failed lines", n["failed_lines"] == n["failed_genomes"]),
        (f"every result file names TBProfiler {WANT_TOOL} and tbdb {WANT_DB}",
         n["versions_other"] == 0 and n["versions_named"] == n["results"]),
    ]


def write_checks(counts, path):
    lines = [f"{k}: {counts[k]}" for k in CHECK_KEYS]
    for name, ok in checks_pass(counts):
        lines.append(f"check: {'pass' if ok else 'FAIL'}: {name}")
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")
    return all(ok for _, ok in checks_pass(counts))


def read_checks(path):
    out = {}
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        key, _, value = line.partition(": ")
        if key in CHECK_KEYS:
            out[key] = value
    missing = [k for k in CHECK_KEYS if k not in out]
    if missing:
        raise Stop(f"checks.txt lacks {', '.join(missing)}")
    return out


# ---------- gates ----------

def gate_layout(reading):
    for need in ("reading.txt", "checks.txt"):
        if not (reading / need).is_file():
            raise Stop(f"gate 1 layout: no {need}")
    if not (reading / "m4").is_dir():
        raise Stop("gate 1 layout: no m4/")
    got = sorted(p.name for p in (reading / "m4").iterdir())
    if got != sorted(M4_FILES):
        raise Stop(f"gate 1 layout: m4/ holds {got}, not tb_curve.py's five files")


def read_m4(reading):
    m4 = reading / "m4"
    tables = {name: load(m4 / name) for name in M4_FILES if name.endswith(".tsv")}
    with open(m4 / "m4-summary.json", encoding="utf-8") as handle:
        summary = json.load(handle)
    return tables, summary


def gate_versions(checks, summary):
    if int(checks["versions_other"]) != 0 or checks["versions_named"] != checks["results"]:
        raise Stop(f"gate 2 versions: {checks['versions_other']} result files name another tool or database")
    if summary.get("tool") != WANT_TOOL or summary.get("database_commit") != WANT_DB:
        raise Stop("gate 2 versions: m4-summary.json names another tool or database")


def gate_accounting(checks, summary):
    for name, ok in checks_pass(checks):
        if not ok and "names TBProfiler" not in name:
            raise Stop(f"gate 3 accounting: {name}: fails")
    planned, results, failed = (int(summary[k]) for k in ("genomes_planned", "genomes_with_result", "genomes_failed"))
    if planned != results + failed:
        raise Stop(f"gate 3 accounting: the summary gives {planned} planned, {results} with a result and {failed} "
                   f"failed")
    if (str(planned), str(results), str(failed)) != (checks["planned"], checks["results"], checks["failed_genomes"]):
        raise Stop("gate 3 accounting: m4-summary.json and checks.txt give other counts")


def gate_datasets(tables, summary, lrows):
    rows = species_rows(lrows, summary["species"])
    by_ab = {r["antibiotic"]: r for r in rows}
    names = {r["species"] for r in rows}
    if len(names) != 1:
        raise Stop(f"gate 4 datasets: lists.tsv names {summary['species']} {len(names)} ways")
    seen = [e["antibiotic"] for e in summary["endpoints"]] + [e["antibiotic"] for e in summary["not_scored_reports"]]
    if sorted(seen) != sorted(by_ab) or len(seen) != len(set(seen)):
        raise Stop("gate 4 datasets: the summary's reports are not exactly the species' reports in lists.tsv")
    for e in summary["endpoints"]:
        if bool(e["in_M"]) != in_m(by_ab[e["antibiotic"]]):
            raise Stop(f"gate 4 datasets: {e['antibiotic']}: the summary's set M mark is not lists.tsv's")
    for name, trows in tables.items():
        for r in trows:
            lr = by_ab.get(r["antibiotic"])
            if lr is None:
                raise Stop(f"gate 4 datasets: {name} names an antibiotic lists.tsv does not")
            if r["species"] not in names:
                raise Stop(f"gate 4 datasets: {name} names the species otherwise than lists.tsv")
            if r["in_M"] != ("1" if in_m(lr) else "0"):
                raise Stop(f"gate 4 datasets: {name}: {r['antibiotic']} carries in_M {r['in_M']}, lists.tsv "
                           f"{lr['in_M']}")


def coverage_hits(tables):
    """Coverage rows that write a share over under 5 resistant genomes."""
    return sum(1 for r in tables["m4-coverage.tsv"]
               if r["share_without"] != "" and not (is_int(r["resistant"]) and int(r["resistant"]) >= FLOOR))


def curve_hits(tables):
    """Curve rows over under 5 resistant genomes."""
    return sum(1 for r in tables["m4-curve.tsv"] if not (is_int(r["resistant"]) and int(r["resistant"]) >= FLOOR))


def convergence_hits(tables, summary):
    """Convergence rows that write a share where the antibiotic has under 5 resistant genomes with a group of that
    grouping, or where the summary has no such antibiotic."""
    hits = 0
    ends = {e["antibiotic"]: e for e in summary["endpoints"]}
    for r in tables["m4-convergence.tsv"]:
        e = ends.get(r["antibiotic"])
        if r["share_covered"] == "":
            continue
        if e is None:
            hits += 1
            continue
        total = int(e["resistant_with_result"]) - int(e["resistant_with_no_group"][r["grouping"]])
        if total < FLOOR:
            hits += 1
    return hits


def share_hits(tables, summary):
    """Rows that write a share over under 5 resistant genomes, or a curve over under 5."""
    return coverage_hits(tables) + curve_hits(tables) + convergence_hits(tables, summary)


def gate_floors(tables, summary):
    hits = share_hits(tables, summary)
    if hits:
        raise Stop(f"gate 5 floors: {hits} rows write a share or a curve over under 5 resistant genomes")


def name_hits(reading, genes):
    """Files under m4/ holding a variant's own name, a position or a sequence, and gene cells that are not a gene
    name of tbdb.bed."""
    hits = 0
    for name in M4_FILES:
        seqs, positions = scan.scan_text(str(reading / "m4" / name))
        hits += seqs + positions
    for r in load(reading / "m4" / "m4-genes.tsv"):
        if r["gene"] != FOLDED and r["gene"] not in genes:
            hits += 1
    return hits


def gate_names(reading, genes):
    hits = name_hits(reading, genes)
    if hits:
        raise Stop(f"gate 6 names: {hits} cells or words name a variant, a position, a sequence or no gene of "
                   f"tbdb.bed")


def gene_hits(tables):
    """Gene rows naming a gene carried by under 5 resistant genomes."""
    return sum(1 for r in tables["m4-genes.tsv"]
               if r["gene"] != FOLDED and not (is_int(r["resistant_with_gene"]) and int(r["resistant_with_gene"]) >= FLOOR))


def gate_genes(tables):
    hits = gene_hits(tables)
    if hits:
        raise Stop(f"gate 7 genes: {hits} gene rows name a gene carried by under 5 resistant genomes")


def check_reading(reading, lrows, genes):
    reading = Path(reading)
    gate_layout(reading)
    tables, summary = read_m4(reading)
    checks = read_checks(reading / "checks.txt")
    gate_versions(checks, summary)
    gate_accounting(checks, summary)
    gate_datasets(tables, summary, lrows)
    gate_floors(tables, summary)
    gate_names(reading, genes)
    gate_genes(tables)
    return tables, summary, checks


# ---------- the page's tables ----------

def page_endpoints(tables, summary, lrows):
    """One row per report of the species, species and antibiotic in full, set M marked as lists.tsv marks it."""
    rows = species_rows(lrows, summary["species"])
    ends = {e["antibiotic"]: e for e in summary["endpoints"]}
    out = []
    for lr in sorted(rows, key=lambda r: r["antibiotic"]):
        ab = lr["antibiotic"]
        o = {"dataset": lr["dataset"], "species": lr["species"], "antibiotic": ab,
             "in_M": "1" if in_m(lr) else "0", "labels_depend_on": lr["labels_depend_on"]}
        e = ends.get(ab)
        if e is None:
            o["kind"] = "none"
            out.append(o)
            continue
        genes = [r for r in tables["m4-genes.tsv"] if r["antibiotic"] == ab]
        overall = [r for r in tables["m4-coverage.tsv"] if r["antibiotic"] == ab and r["grouping"] == "all"]
        if len(overall) != 1:
            raise Stop(f"m4-coverage.tsv holds {len(overall)} overall rows for {ab}")
        o.update({"kind": e["kind"], "tool_drugs": ",".join(e["tool_drugs"]),
                  "resistant_with_result": e["resistant_with_result"], "resistant_no_result": e["resistant_no_result"],
                  "catalogue_variants": sum(int(r["catalogue_variants"]) for r in genes),
                  "genes_named": sum(1 for r in genes if r["gene"] != FOLDED),
                  "without_determinant": overall[0]["without_determinant"], "share_without": overall[0]["share_without"]})
        for g, short in zip(GROUPINGS, ("main", "sub", "cluster")):
            groups = {r["group"] for r in tables["m4-coverage.tsv"] if r["antibiotic"] == ab and r["grouping"] == g}
            o[f"groups_{short}"] = len(groups)
            o[f"no_{short}"] = e["resistant_with_no_group"][g]
        for g, short in (("main_lineage", "main"), ("sub_lineage", "sub")):
            one = [r for r in tables["m4-convergence.tsv"]
                   if r["antibiotic"] == ab and r["grouping"] == g and r["lineages_seen"] == "1"]
            o[f"one_lineage_{short}"] = one[0]["determinants"] if one else ""
        out.append(o)
    return out


README = """# The lineage curve's reading of {date}: the record

This directory holds the tables behind [docs/lineage-curve-result.md](../../docs/lineage-curve-result.md): block M4 of the mechanism split, the lineage curve for Mycobacterium tuberculosis, read once over TBProfiler {tool} with `tbdb` commit {db} run again on the {planned} genomes of the species' {reports} revision 14 reports. Every row is an aggregate over an antibiotic, a lineage group or a gene. No file holds a genome id, a BioProject accession, a variant's own name, a sequence or a position.

`scripts/mechanism/tb_curve.py` wrote the five files under `m4/`; they are copied here as it wrote them. `scripts/mechanism/curve_tables.py write` checked the reading against the seven gates of `gates.txt` and wrote `page-endpoints.tsv` from those files and `scripts/c3-inputs/lists.tsv`, species and antibiotic in full, with no new measure. `reading.txt` lists the reading's commands with the sha256 of what each wrote, and `checks.txt` the checks of the tool run, counts only. The rules of the reading are `docs/mechanism-split.md`'s, section "M3 and M4", fixed before any run.

| File | Holds |
|---|---|
| `gates.txt` | the seven gates the reading passed before a table was written |
| `reading.txt` | the reading's commands and the sha256 of each one's output |
| `checks.txt` | the checks of the tool run, counts only |
| `m4/m4-curve.tsv` | per antibiotic, grouping and step: distinct graded variants as lineage groups are added, the median and the 2.5 and 97.5 percentiles over 1,000 seeded random orders, and the count in the order by group size, largest first |
| `m4/m4-coverage.tsv` | per antibiotic, grouping and group, and over all: resistant genomes, those with no graded variant, and the share where 5 or more |
| `m4/m4-convergence.tsv` | per antibiotic and grouping: graded variants seen in 1, 2, 3 and 4 or more lineages, and the resistant genomes each count covers |
| `m4/m4-genes.tsv` | per antibiotic and gene, as the WHO catalogue names it: graded variants seen and resistant genomes with one; genes in 1 to 4 resistant genomes share the row `other genes` |
| `m4/m4-summary.json` | the tool, the database, the genomes planned, with a result and failed, and per report its kind, tool drugs, set M mark and the resistant genomes with no result or no group |
| `page-endpoints.tsv` | one row per report, species and antibiotic in full, set M marked: the counts the page reads |

"The WHO catalogue" here is the catalogue as TBProfiler {tool} ships it (`tbdb` {db}). A curve describes the deposits of one public database. A share is no prevalence: the genomes are BV-BRC deposits chosen for having lab results. A split by lineage is a retrospective split inside one database, never external validation. A table names no cause and no new marker.

`SHA256SUMS` gives the sha256 of every other file here. `scripts/mechanism/curve_tables.py print` prints the page's tables from these files.
"""


def write(reading, lists, bed, out):
    reading, out = Path(reading), Path(out)
    lrows = load(lists)
    tables, summary, checks = check_reading(reading, lrows, bed_genes(bed))
    if out.exists() and any(out.iterdir()):
        raise Stop("OUT is not empty")
    (out / "m4").mkdir(parents=True, exist_ok=True)
    for name in M4_FILES:
        shutil.copyfile(reading / "m4" / name, out / "m4" / name)
    for name in ("reading.txt", "checks.txt"):
        shutil.copyfile(reading / name, out / name)
    write_tsv(out / "page-endpoints.tsv", page_endpoints(tables, summary, lrows), ENDPOINT_COLUMNS)
    with open(out / "gates.txt", "w", encoding="utf-8") as handle:
        for n, name, words in GATES:
            handle.write(f"{n} {name}: {words}\n")
    date = "an unknown date"
    if checks["end"]:
        day = datetime.strptime(checks["end"], "%Y-%m-%dT%H:%M:%SZ")
        date = f"{day.day} {day:%B %Y}"
    (out / "README.md").write_text(README.format(date=date, tool=WANT_TOOL, db=WANT_DB,
                                                 planned=f"{int(checks['planned']):,}",
                                                 reports=len(species_rows(lrows, summary["species"]))),
                                   encoding="utf-8")
    files = sorted(p.relative_to(out).as_posix() for p in out.rglob("*") if p.is_file())
    with open(out / "SHA256SUMS", "w", encoding="utf-8") as handle:
        for name in files:
            handle.write(f"{sha256_of(out / name)}  {name}\n")
    return files


# ---------- printing ----------

def d3(x):
    if x in ("", None):
        return "none"
    value = Decimal(str(x)).quantize(Decimal("0.001"), rounding=ROUND_HALF_EVEN)
    if value == 0:
        value = Decimal("0.000")
    return str(value)


def count_words(x):
    return "none" if x in ("", None) else f"{int(x):,}"


def num(x):
    """A curve figure as written: an int, or a median of two that may end in .5."""
    v = Decimal(str(x))
    return str(int(v)) if v == v.to_integral_value() else str(v.normalize())


def band(r):
    return f"{num(r['determinants_median'])} ({num(r['determinants_lo'])} to {num(r['determinants_hi'])})"


def table(header, rows):
    lines = ["| " + " | ".join(header) + " |", "|" + "|".join("---" for _ in header) + "|"]
    lines += ["| " + " | ".join(str(c) for c in row) + " |" for row in rows]
    return lines


def set_word(r):
    return "M" if r["in_M"] == "1" else "not M"


def scored(endpoints):
    return [r for r in endpoints if r["kind"] != "none"]


def curve_of(curve, ab, grouping):
    return sorted((r for r in curve if r["antibiotic"] == ab and r["grouping"] == grouping), key=lambda r: int(r["step"]))


def print_reading(endpoints, curve):
    """Per scored report: distinct graded variants, the main-lineage curve's first and last steps, the share with no
    graded variant and the variants seen in one main lineage only."""
    out = []
    for r in scored(endpoints):
        steps = curve_of(curve, r["antibiotic"], "main_lineage")
        first = band(steps[0]) if steps else "none"
        last = num(steps[-1]["determinants_median"]) if steps else "none"
        out.append([r["antibiotic"], r["kind"], set_word(r), count_words(r["resistant_with_result"]),
                    r["catalogue_variants"], r["groups_main"], first, last,
                    f"{count_words(r['without_determinant'])} ({d3(r['share_without'])})",
                    r["one_lineage_main"] or "none"])
    return table(["Antibiotic", "Kind", "Set", "Resistant genomes", "Graded variants", "Main lineages",
                  "First lineage, median (band)", "All lineages", "No graded variant (share)",
                  "In one main lineage only"], out)


def print_curves(endpoints, curve, grouping):
    """The curve of each scored report under one grouping: every step, median (band) and largest first."""
    out = []
    for r in scored(endpoints):
        for s in curve_of(curve, r["antibiotic"], grouping):
            out.append([r["antibiotic"], f"{s['step']} of {s['groups']}", band(s), s["determinants_largest_first"]])
    return table(["Antibiotic", "Groups added", "Median (2.5 to 97.5 percentile)", "Largest first"], out)


def print_curve_ends(endpoints, curve):
    """Under the sublineage and the cluster groupings: groups, the first step and the last."""
    out = []
    for r in scored(endpoints):
        row = [r["antibiotic"]]
        for grouping in ("sub_lineage", "cluster"):
            steps = curve_of(curve, r["antibiotic"], grouping)
            row += ([steps[0]["groups"], band(steps[0]), steps[0]["determinants_largest_first"],
                     num(steps[-1]["determinants_median"])] if steps else ["none"] * 4)
        out.append(row)
    return table(["Antibiotic", "Sublineages", "First sublineage, median (band)", "Largest sublineage",
                  "All sublineages", "Clusters", "First cluster, median (band)", "Largest cluster", "All clusters"], out)


def print_without(endpoints, coverage):
    """No graded variant: over all, then per main lineage."""
    groups = sorted({r["group"] for r in coverage if r["grouping"] == "main_lineage"})
    out = []
    for r in scored(endpoints):
        cells = {c["group"]: c for c in coverage if c["antibiotic"] == r["antibiotic"] and c["grouping"] == "main_lineage"}
        row = [r["antibiotic"], f"{count_words(r['without_determinant'])} of {count_words(r['resistant_with_result'])}",
               d3(r["share_without"])]
        for g in groups:
            c = cells.get(g)
            row.append("none" if c is None else f"{c['without_determinant']} of {c['resistant']}"
                       + (f" ({d3(c['share_without'])})" if c["share_without"] else ""))
        out.append(row)
    return table(["Antibiotic", "No graded variant", "Share"] + groups, out)


def print_spread(endpoints, convergence):
    out = []
    for r in scored(endpoints):
        for grouping in ("main_lineage", "sub_lineage"):
            rows = {c["lineages_seen"]: c for c in convergence
                    if c["antibiotic"] == r["antibiotic"] and c["grouping"] == grouping}
            if not rows:
                continue
            row = [r["antibiotic"], "main lineage" if grouping == "main_lineage" else "sublineage"]
            for b in ("1", "2", "3", "4+"):
                c = rows[b]
                row.append(f"{c['determinants']} ({count_words(c['resistant_covered'])}"
                           + (f", {d3(c['share_covered'])}" if c["share_covered"] else "") + ")")
            out.append(row)
    return table(["Antibiotic", "Grouping", "In 1 (resistant covered, share)", "In 2", "In 3", "In 4 or more"], out)


def print_genes(endpoints, genes):
    out = []
    for r in scored(endpoints):
        rows = [g for g in genes if g["antibiotic"] == r["antibiotic"]]
        named = sorted((g for g in rows if g["gene"] != FOLDED), key=lambda g: (-int(g["resistant_with_gene"]), g["gene"]))
        for g in named + [g for g in rows if g["gene"] == FOLDED]:
            gene = g["gene"] if g["gene"] == FOLDED else f"*{g['gene']}*"
            out.append([r["antibiotic"], gene, g["catalogue_variants"], count_words(g["resistant_with_gene"])])
    return table(["Antibiotic", "Gene", "Graded variants seen", "Resistant genomes with one"], out)


def print_datasets(endpoints):
    out = []
    for r in endpoints:
        out.append([f"*{r['species']}*", r["antibiotic"], set_word(r), r["kind"],
                    r["tool_drugs"].replace(",", ", ").replace("_", " ") or "none",
                    count_words(r.get("resistant_with_result")), count_words(r.get("resistant_no_result")),
                    count_words(r.get("no_main")), count_words(r.get("no_cluster"))])
    return table(["Species", "Antibiotic", "Set", "Kind", "Tool drugs", "Resistant with a result",
                  "Resistant with no result", "No main lineage", "No cluster"], out)


def tables(folder):
    folder = Path(folder)
    endpoints = load(folder / "page-endpoints.tsv")
    curve = load(folder / "m4" / "m4-curve.tsv")
    return {
        "reading": print_reading(endpoints, curve),
        "curves_main": print_curves(endpoints, curve, "main_lineage"),
        "curve_ends": print_curve_ends(endpoints, curve),
        "without": print_without(endpoints, load(folder / "m4" / "m4-coverage.tsv")),
        "spread": print_spread(endpoints, load(folder / "m4" / "m4-convergence.tsv")),
        "genes": print_genes(endpoints, load(folder / "m4" / "m4-genes.tsv")),
        "datasets": print_datasets(endpoints),
    }


def main(argv):
    if len(argv) == 3 and argv[0] == "check":
        counts = run_counts(argv[1])
        ok = write_checks(counts, argv[2])
        print(Path(argv[2]).read_text(encoding="utf-8"), end="")
        return 0 if ok else 1
    if len(argv) == 5 and argv[0] == "write":
        try:
            files = write(*argv[1:])
        except Stop as e:
            print(f"stopped: {e}", file=sys.stderr)
            return 1
        print(f"gates: {len(GATES)} of {len(GATES)} pass; wrote {len(files) + 1} files")
        return 0
    if len(argv) == 2 and argv[0] == "print":
        for name, lines in tables(argv[1]).items():
            print(f"[{name}]")
            print("\n".join(lines))
            print()
        return 0
    sys.stderr.write(__doc__)
    return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
