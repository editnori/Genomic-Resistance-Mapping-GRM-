"""catalogue_size.py: does the number of graded WHO catalogue variants of an
antibiotic track its resistant share, or only how many resistant genomes there
were to study? docs/catalogue-size-result.md, read once by the rule its
section "Rule" fixed before any count was read.

    python3 -I catalogue_size.py plan
    python3 -I catalogue_size.py read REPORTS ENDPOINTS READING
    python3 -I catalogue_size.py write READING ENDPOINTS SPECIES_DIR OUT
    python3 -I catalogue_size.py print OUT

`plan` prints the rule as the page states it.

`read` matches each row of ENDPOINTS (the lineage curve's
reference/lineage-curve/page-endpoints.tsv) to the release report
REPORTS/<dataset>.json, takes from each report `counts.samples`,
`counts.resistant` and `counts.susceptible` and nothing else, and writes
READING/release-counts.tsv (one row per report) and READING/read.txt (guards
1, 2 and 6). A report's `samples` list, whose entries carry the genome ids, is
dropped as the file is parsed. READING must be absent or empty.

`write` proves the rank correlation on three hand-made cases (guard 4), joins
READING/release-counts.tsv to ENDPOINTS, and writes under OUT antibiotics.tsv
(one row per report, species and antibiotic in full), correlations.tsv (one
row per set and pair, guard 3), README.md, reading.txt (copied from READING,
each sha256 it lists checked first), guards.txt and SHA256SUMS. Before
guards.txt takes guard 5's line, every file under OUT and READING is scanned
for every genome id of SPECIES_DIR/*.json (the release's species files, their
`genomes[].id`, held in memory only) and for a BioProject accession, a
sequence or a position; any hit stops the run. OUT must be absent or empty.

`print` prints the page's tables from OUT alone, figures rounded to 3
decimals.

Each guard prints counts only and stops the run with exit 1 when it fails. No
p-value, no test, no fit, no fitted line and no model: the correlation is
Pearson's over average ranks, written here.

Python 3.12, standard library only; run with -I.
"""

import csv
import hashlib
import json
import math
import os
import re
import sys
from decimal import ROUND_HALF_EVEN, Decimal
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import scan  # noqa: E402

MARGIN = 0.15
TOL = 1e-9
MIN_SET = 3
Y = "catalogue_variants"
# Guard 5's second net, a genome id's shape: digits, a dot and digits, the first part 3 to 7 digits not led by 0, with
# no digit on either side. A letter or an underscore beside it does not hide it.
ID_SHAPE = re.compile(rb"(?<!\d)[1-9]\d{2,6}\.\d+(?!\d)")
PAIRS = ("resistant_share", "resistant", "genomes")
SETS = ("rule", "primary only")
ROLES = {"rule": "the rule's set", "primary only": "beside, decides nothing"}
BRANCHES = ("tracks resistant genomes", "tracks share", "neither by the margin")
REASONS = ("not read by the tool", "no graded variant", "set M")
COUNTS = ("samples", "resistant", "susceptible")
ENDPOINT_NEEDS = ("dataset", "species", "antibiotic", "in_M", "kind", "resistant_with_result", "resistant_no_result",
                  "catalogue_variants", "genes_named", "without_determinant", "share_without")
CARRIED = ("in_M", "kind", "resistant_with_result", "catalogue_variants", "genes_named", "without_determinant",
           "share_without")
RELEASE_COLUMNS = ["dataset", "genomes", "resistant", "susceptible"]
TABLE_COLUMNS = ["species", "antibiotic", "in_M", "kind", "genomes", "resistant", "susceptible", "resistant_share",
                 "resistant_with_result", "catalogue_variants", "genes_named", "without_determinant", "share_without",
                 "in_correlation", "left_out_because"]
CORRELATION_COLUMNS = ["set", "role", "x", "y", "n", "rho", "loo_min", "loo_min_without", "loo_max", "loo_max_without",
                       "decides", "comparison_gives", "antibiotics"]
# The three hand-made cases of guard 4: x, y and the answer worked out by hand.
# 1: ranks of y are 1, 3, 2, 5, 4, so the squared rank differences sum to 4 and rho = 1 - 6 * 4 / (5 * 24) = 0.8.
# 2: y falls as x rises: rho = -1.
# 3: x ties at 2: ranks 1, 2.5, 2.5, 4 against 1, 2, 3, 4; rho = 4.5 / sqrt(4.5 * 5) = sqrt(0.9).
CASES = (
    ("no ties", (1, 2, 3, 4, 100), (1, 30, 20, 500, 400), 0.8),
    ("reversed", (1, 2, 3, 4, 5), (50, 40, 30, 20, 10), -1.0),
    ("ties", (1, 2, 2, 4), (1, 2, 3, 4), math.sqrt(0.9)),
)
RULE = """This rule was written and committed before the release's counts were read, and it is not changed after.

- **Antibiotics.** The rows of `reference/lineage-curve/page-endpoints.tsv` whose `catalogue_variants` is 1 or more and whose `in_M` is 0. A row with `in_M` 1 waits on the maintainer's label rule: it stands in the table and is left out of every correlation. A row the tool does not read, or with no graded variant, stands in the table and is left out too. Every correlation of a set uses the same antibiotics.
- **Pairs.** The Spearman rank correlation of `catalogue_variants` against each of: the resistant share (`counts.resistant` over `counts.samples` of the antibiotic's release report), the resistant genomes (`counts.resistant`) and the genomes (`counts.samples`). Tied values take the average of their ranks, and the correlation is Pearson's over the ranks.
- **Leave one out.** For each pair, the correlation over the set with one antibiotic dropped, for each antibiotic in turn: the smallest and the largest, each with the antibiotic dropped.
- **Decision.** Let R be the correlation with the resistant genomes and S the one with the resistant share, each in absolute value and unrounded. If R is as large as S or larger, the answer is that catalogue size tracks how many resistant genomes there were to study, and target size is not separable here. If S is larger than R by 0.15 or more, the answer is that the number of graded variants tracks the resistant share, in one species, one catalogue and one database, and no causal reading follows. If S is larger than R by less than 0.15, the answer is that these numbers cannot separate the two. A difference within 0.000000001 of a bound counts as meeting it.
- **Beside.** The correlation with the genomes is reported and decides nothing. The same three pairs over the antibiotics of the set that the tool reads by their own drug (kind `primary`) are reported and decide nothing: a proxy counts another drug's graded variants.
- **Not done.** No p-value, no test, no fit, no fitted line and no model."""
GUARDS = [
    (1, "match", "every antibiotic of page-endpoints.tsv is matched to exactly one release report by its dataset name; "
                 "the unmatched on either side are listed and counted"),
    (2, "sums", "each report's counts.resistant and counts.susceptible add to its counts.samples"),
    (3, "sets", "the antibiotics of every correlation of a set are the same, listed, with n printed beside each "
                "correlation"),
    (4, "rank correlation", "the rank correlation gives the known answers of three hand-made cases, one with ties, "
                            "before it touches a real number"),
    (5, "ids", "every genome id of the release's species files is searched for in full in every written file and "
               "stands in none; no token of a genome id's shape stands there either, and the kit's BioProject, "
               "sequence and position scans find nothing"),
    (6, "resistant", "for each report the tool reads, counts.resistant equals the page's resistant_with_result plus "
                     "resistant_no_result"),
]


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


def r12(x):
    return round(float(x) + 0.0, 12)


def f12(x):
    return f"{r12(x):.12f}"


def is_count(text):
    return isinstance(text, str) and text.isascii() and text.isdigit()


def empty_dir(path):
    path = Path(path)
    if path.exists() and any(path.iterdir()):
        raise Stop(f"{path.name} is not empty")
    path.mkdir(parents=True, exist_ok=True)
    return path


# ---------- the rank correlation ----------

def ranks(values):
    """Ranks from 1, tied values each taking the average of the ranks they span."""
    order = sorted(range(len(values)), key=lambda i: values[i])
    out = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
            j += 1
        for k in range(i, j + 1):
            out[order[k]] = (i + j) / 2 + 1
        i = j + 1
    return out


def pearson(x, y):
    if len(x) != len(y) or len(x) < 2:
        raise Stop("a correlation needs two columns of one length, 2 or more")
    mx, my = sum(x) / len(x), sum(y) / len(y)
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    sxx = sum((a - mx) ** 2 for a in x)
    syy = sum((b - my) ** 2 for b in y)
    if sxx == 0 or syy == 0:
        raise Stop("a column whose values are all tied has no rank correlation")
    return sxy / math.sqrt(sxx * syy)


def spearman(x, y):
    return pearson(ranks(x), ranks(y))


def leave_one_out(x, y, names):
    """(smallest, its antibiotic dropped, largest, its antibiotic dropped) over the set less one, each in turn."""
    values = [(r12(spearman(x[:i] + x[i + 1:], y[:i] + y[i + 1:])), names[i]) for i in range(len(names))]
    low = min(values, key=lambda v: (v[0], v[1]))
    high = min(values, key=lambda v: (-v[0], v[1]))
    return low[0], low[1], high[0], high[1]


def prove_spearman():
    """Guard 4: the three hand-made cases, before any real number."""
    for name, x, y, want in CASES:
        got = spearman(list(x), list(y))
        if abs(got - want) > 1e-12:
            raise Stop(f"guard 4 rank correlation: the case {name!r} gives {got!r}, not {want!r}")
    return f"4 rank correlation: {len(CASES)} of {len(CASES)} hand-made cases give their known answers"


def branch(rho_resistant, rho_share):
    """The rule's decision, on the unrounded correlations."""
    r, s = abs(rho_resistant), abs(rho_share)
    if r >= s - TOL:
        return BRANCHES[0]
    if s - r >= MARGIN - TOL:
        return BRANCHES[1]
    return BRANCHES[2]


# ---------- read: the release's three counts ----------

def drop_samples(pairs):
    # a report's `samples` list (each entry carries a genome id) is dropped as the file is parsed
    return {k: v for k, v in pairs if not (k == "samples" and isinstance(v, list))}


def counts_of(path):
    with open(path, encoding="utf-8") as handle:
        report = json.load(handle, object_pairs_hook=drop_samples)
    counts = report.get("counts") if isinstance(report, dict) else None
    if not isinstance(counts, dict):
        raise Stop(f"{Path(path).name}: no counts")
    out = {}
    for key in COUNTS:
        v = counts.get(key)
        if not isinstance(v, int) or isinstance(v, bool) or v < 0:
            raise Stop(f"{Path(path).name}: counts.{key} is not a count")
        out[key] = v
    return out


def endpoints_of(path):
    rows = load(path)
    if not rows:
        raise Stop("page-endpoints.tsv holds no row")
    missing = [c for c in ENDPOINT_NEEDS if c not in rows[0]]
    if missing:
        raise Stop(f"page-endpoints.tsv lacks {', '.join(missing)}")
    for key in ("dataset", "antibiotic"):
        names = [r[key] for r in rows]
        if len(set(names)) != len(names):
            raise Stop(f"page-endpoints.tsv names a {key} twice")
    slugs = {r["dataset"].split("__", 1)[0] for r in rows}
    if len(slugs) != 1 or any("__" not in r["dataset"] for r in rows):
        raise Stop("page-endpoints.tsv holds more than one species")
    return rows, slugs.pop()


def report_names(reports, slug):
    return sorted(p.name[:-len(".json")] for p in Path(reports).iterdir()
                  if p.is_file() and p.name.startswith(slug + "__") and p.name.endswith(".json"))


def gate_match(page_datasets, report_datasets):
    """Guard 1: each page row to exactly one report, and back."""
    page, reports = list(page_datasets), list(report_datasets)
    no_report = sorted(set(page) - set(reports))
    no_row = sorted(set(reports) - set(page))
    if no_report or no_row:
        raise Stop(f"guard 1 match: {len(no_report)} page rows with no report ({', '.join(no_report) or 'none'}), "
                   f"{len(no_row)} reports with no page row ({', '.join(no_row) or 'none'})")
    return (f"1 match: {len(page)} page rows, {len(reports)} reports of the species, {len(page)} matched, "
            f"0 page rows with no report, 0 reports with no page row")


def gate_sums(counts):
    """Guard 2: resistant plus susceptible is the report's genomes."""
    bad = sorted(d for d, c in counts.items() if c["resistant"] + c["susceptible"] != c["samples"])
    if bad:
        raise Stop(f"guard 2 sums: {len(bad)} reports do not add up ({', '.join(bad)})")
    return f"2 sums: {len(counts)} of {len(counts)} reports add up"


def gate_page_agrees(rows, counts):
    """Guard 6: the report's resistant genomes are the page's, where the tool reads the report."""
    read = [r for r in rows if r["kind"] != "none"]
    bad = []
    for r in read:
        if not (is_count(r["resistant_with_result"]) and is_count(r["resistant_no_result"])):
            bad.append(r["dataset"])
        elif counts[r["dataset"]]["resistant"] != int(r["resistant_with_result"]) + int(r["resistant_no_result"]):
            bad.append(r["dataset"])
    if bad:
        raise Stop(f"guard 6 resistant: {len(bad)} reports differ from the page ({', '.join(sorted(bad))})")
    return f"6 resistant: {len(read)} of {len(read)} reports the tool reads agree with the page"


def read(reports, endpoints, reading):
    rows, slug = endpoints_of(endpoints)
    lines = [gate_match([r["dataset"] for r in rows], report_names(reports, slug))]
    counts = {r["dataset"]: counts_of(Path(reports) / f"{r['dataset']}.json") for r in rows}
    lines.append(gate_sums(counts))
    lines.append(gate_page_agrees(rows, counts))
    reading = empty_dir(reading)
    write_tsv(reading / "release-counts.tsv",
              [{"dataset": d, "genomes": c["samples"], "resistant": c["resistant"], "susceptible": c["susceptible"]}
               for d, c in counts.items()], RELEASE_COLUMNS)
    (reading / "read.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    return lines


# ---------- write: the table and the correlations ----------

def left_out_because(r):
    if r["kind"] == "none":
        return REASONS[0]
    if r["catalogue_variants"] == "0":
        return REASONS[1]
    if r["in_M"] == "1":
        return REASONS[2]
    return ""


def share_of(resistant, genomes):
    return Fraction(resistant, genomes)


def antibiotic_rows(endpoints, release):
    by = {r["dataset"]: r for r in release}
    out = []
    for e in endpoints:
        c = by[e["dataset"]]
        for key in RELEASE_COLUMNS[1:]:
            if not is_count(c[key]):
                raise Stop(f"release-counts.tsv: {e['dataset']} {key} is not a count")
        genomes, resistant = int(c["genomes"]), int(c["resistant"])
        if genomes == 0:
            raise Stop(f"release-counts.tsv: {e['dataset']} has no genome")
        reason = left_out_because(e)
        if reason == "" and not is_count(e[Y]):
            raise Stop(f"page-endpoints.tsv: {e['antibiotic']} {Y} is not a count")
        row = {"species": e["species"], "antibiotic": e["antibiotic"], "genomes": genomes, "resistant": resistant,
               "susceptible": int(c["susceptible"]), "resistant_share": f12(share_of(resistant, genomes)),
               "in_correlation": "0" if reason else "1", "left_out_because": reason,
               "_share": share_of(resistant, genomes)}
        row.update({k: e[k] for k in CARRIED})
        out.append(row)
    return out


def set_rows(rows, name):
    chosen = [r for r in rows if r["in_correlation"] == "1"]
    return chosen if name == "rule" else [r for r in chosen if r["kind"] == "primary"]


def rows_for(rows, name, x):
    return set_rows(rows, name)


def value(r, x):
    return r["_share"] if x == "resistant_share" else int(r[x])


def correlate(rows, name, x):
    chosen = rows_for(rows, name, x)
    if len(chosen) < MIN_SET:
        raise Stop(f"the set {name!r} holds {len(chosen)} antibiotics, under {MIN_SET}")
    names = [r["antibiotic"] for r in chosen]
    xs, ys = [value(r, x) for r in chosen], [int(r[Y]) for r in chosen]
    rho = spearman(xs, ys)
    low, low_without, high, high_without = leave_one_out(xs, ys, names)
    return {"set": name, "role": ROLES[name], "x": x, "y": Y, "n": len(chosen), "_rho": rho, "rho": f12(rho), "loo_min": f12(low),
            "loo_min_without": low_without, "loo_max": f12(high), "loo_max_without": high_without,
            "antibiotics": "; ".join(names)}


def gate_sets(correlations):
    """Guard 3: one list of antibiotics per set, n beside each correlation."""
    lines = []
    for name in SETS:
        mine = [c for c in correlations if c["set"] == name]
        if len({c["antibiotics"] for c in mine}) != 1 or len(mine) != len(PAIRS):
            raise Stop(f"guard 3 sets: the pairs of the set {name!r} use different antibiotics")
        lines.append(f"3 sets: {name}: " + ", ".join(f"{c['x']} n={c['n']}" for c in mine)
                     + f"; antibiotics: {mine[0]['antibiotics']}")
    return lines


def decide(correlations):
    for name in SETS:
        mine = {c["x"]: c for c in correlations if c["set"] == name}
        gives = branch(mine["resistant"]["_rho"], mine["resistant_share"]["_rho"])
        for c in mine.values():
            c["comparison_gives"] = gives
            c["decides"] = "yes" if name == "rule" and c["x"] in ("resistant", "resistant_share") else "no"
    return correlations


def ids_of(species_dir):
    ids = set()
    files = sorted(Path(species_dir).glob("*.json"))
    if not files:
        raise Stop("guard 5 ids: no species file")
    for path in files:
        ids |= scan.load_ids(path)
    return ids, len(files)


def id_hits(paths, ids):
    """Guard 5's nets over the written files, counts only: every id of the list searched for as a plain substring of
    each file, wherever it stands; tokens of a genome id's shape, the second net; and the kit's BioProject, sequence
    and position scans."""
    listed = shaped = other = 0
    for path in paths:
        data = Path(path).read_bytes()
        listed += sum(data.count(i) for i in ids)
        shaped += len(ID_SHAPE.findall(data))
        seq, pos = scan.scan_text(path)
        other += len(scan.BIOPROJECT.findall(data)) + seq + pos
    return listed, shaped, other


def gate_ids(out, reading, species_dir):
    """Guard 5: every written file searched for every genome id of the release, with the shape of an id counted
    beside; counts only."""
    ids, species = ids_of(species_dir)
    paths = sorted(p for d in (out, reading) for p in Path(d).rglob("*") if p.is_file())
    listed, shaped, other = id_hits(paths, ids)
    if listed or shaped or other:
        raise Stop(f"guard 5 ids: {listed} hits of the id list, {shaped} tokens of a genome id's shape and {other} "
                   f"BioProject, sequence or position hits in {len(paths)} files")
    return (f"5 ids: {len(paths)} files; 0 hits of {len(ids)} genome ids of {species} species files, each searched "
            f"for in full; 0 tokens of a genome id's shape; 0 BioProject, sequence or position hits")


def check_reading_sums(reading):
    lines = (Path(reading) / "reading.txt").read_text(encoding="utf-8").splitlines()
    listed = 0
    for line in lines:
        digest, _, name = line.partition("  ")
        if len(digest) == 64 and all(ch in "0123456789abcdef" for ch in digest) and name:
            if not (Path(reading) / name).is_file() or sha256_of(Path(reading) / name) != digest:
                raise Stop(f"reading.txt: {name} does not have the sha256 it lists")
            listed += 1
    if not listed:
        raise Stop("reading.txt lists no sha256")
    return listed


def write(reading, endpoints, species_dir, out):
    reading = Path(reading)
    first = prove_spearman()
    rows, _ = endpoints_of(endpoints)
    release = load(reading / "release-counts.tsv")
    check_reading_sums(reading)
    lines = (reading / "read.txt").read_text(encoding="utf-8").splitlines()
    gate_match([r["dataset"] for r in rows], [r["dataset"] for r in release])
    table = antibiotic_rows(rows, release)
    correlations = decide([correlate(table, name, x) for name in SETS for x in PAIRS])
    lines += gate_sets(correlations) + [first]
    out = empty_dir(out)
    write_tsv(out / "antibiotics.tsv", table, TABLE_COLUMNS)
    write_tsv(out / "correlations.tsv", correlations, CORRELATION_COLUMNS)
    (out / "README.md").write_text(README, encoding="utf-8")
    (out / "reading.txt").write_bytes((reading / "reading.txt").read_bytes())
    (out / "guards.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    lines.append(gate_ids(out, reading, species_dir))
    (out / "guards.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    files = sorted(p.relative_to(out).as_posix() for p in out.rglob("*") if p.is_file())
    with open(out / "SHA256SUMS", "w", encoding="utf-8") as handle:
        for name in files:
            handle.write(f"{sha256_of(out / name)}  {name}\n")
    return lines, correlations


README = """# The catalogue size reading: the record

This directory holds the tables behind [docs/catalogue-size-result.md](../../docs/catalogue-size-result.md): whether the number of graded variants of the WHO catalogue among an antibiotic's resistant genomes tracks the resistant share of its genomes or only how many resistant genomes there were to study, for the antibiotics of Mycobacterium tuberculosis. It was read once, by the rule the page's section "Rule" fixed before any count was read. Every row is an aggregate over an antibiotic or over a set of antibiotics. No file holds a genome id, a BioProject accession, a variant's own name, a gene's name, a sequence or a position.

`scripts/mechanism/catalogue_size.py read` took three counts from each of the species' release reports, `counts.samples`, `counts.resistant` and `counts.susceptible`, and nothing else. `scripts/mechanism/catalogue_size.py write` joined them to the lineage curve's `reference/lineage-curve/page-endpoints.tsv` (page `docs/lineage-curve-result.md`) and wrote the two tables. `in_M`, `kind`, `resistant_with_result`, `catalogue_variants`, `genes_named`, `without_determinant` and `share_without` are carried from that file unchanged.

| File | Holds |
|---|---|
| `antibiotics.tsv` | one row per report of the species, species and antibiotic in full: its genomes, resistant and susceptible genomes and resistant share from the release report, the lineage curve's columns carried unchanged, set M marked, and whether the antibiotic enters the correlations and, if not, why |
| `correlations.tsv` | one row per set and pair: the set and its role, the rule's set or beside it and deciding nothing; the Spearman rank correlation of `catalogue_variants` against the resistant share, the resistant genomes or the genomes, its n, the smallest and largest value with one antibiotic dropped and which one, whether the row decides, what the rule's comparison gives over the set, and the set's antibiotics |
| `guards.txt` | the six guards the reading passed, counts only |
| `reading.txt` | the reading's commands and the sha256 of each one's output |

A share is no prevalence: the genomes are BV-BRC deposits chosen for having lab results. A correlation over antibiotics of one species, one catalogue and one database names no cause. These are public bacterial genomes with public accession numbers, and the record holds no field about a patient.

`SHA256SUMS` gives the sha256 of every other file here. `scripts/mechanism/catalogue_size.py print` prints the page's tables from these files.
"""


# ---------- printing ----------

def d3(x):
    if x in ("", None):
        return "none"
    v = Decimal(str(x)).quantize(Decimal("0.001"), rounding=ROUND_HALF_EVEN)
    return str(Decimal("0.000") if v == 0 else v)


def count_words(x):
    return "none" if x in ("", None) else f"{int(x):,}"


def table(header, rows):
    return ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)] + \
        ["| " + " | ".join(r) + " |" for r in rows]


PAIR_WORDS = {"resistant_share": "resistant share", "resistant": "resistant genomes", "genomes": "genomes"}


def print_antibiotics(rows):
    return table(["Antibiotic", "Set", "Kind", "Genomes", "Resistant", "Susceptible", "Resistant share",
                  "Graded variants", "Genes named", "No graded variant (share)", "In the correlations"],
                 [[r["antibiotic"], "M" if r["in_M"] == "1" else "not M", r["kind"], count_words(r["genomes"]),
                   count_words(r["resistant"]), count_words(r["susceptible"]), d3(r["resistant_share"]),
                   count_words(r["catalogue_variants"]), count_words(r["genes_named"]),
                   "none" if r["without_determinant"] == "" else
                   f"{count_words(r['without_determinant'])} ({d3(r['share_without'])})",
                   "yes" if r["in_correlation"] == "1" else f"no: {r['left_out_because']}"] for r in rows])


def print_correlations(rows):
    return table(["Set", "Role", "Graded variants against", "n", "Spearman", "One dropped, smallest (dropped)",
                  "One dropped, largest (dropped)", "Decides", "The comparison gives"],
                 [[c["set"], c["role"], PAIR_WORDS[c["x"]], c["n"], d3(c["rho"]), f"{d3(c['loo_min'])} ({c['loo_min_without']})",
                   f"{d3(c['loo_max'])} ({c['loo_max_without']})", c["decides"], c["comparison_gives"]]
                  for c in rows])


def print_sets(rows):
    seen = {}
    for c in rows:
        seen.setdefault(c["set"], (c["n"], c["antibiotics"]))
    return table(["Set", "n", "Antibiotics"], [[name, n, names.replace("; ", ", ")] for name, (n, names) in seen.items()])


def print_guards(lines):
    words = {str(n): name for n, name, _ in GUARDS}
    return table(["Guard", "Printed"], [[f"{line.split(' ', 1)[0]}. {words[line.split(' ', 1)[0]]}",
                                         line.split(": ", 1)[1].replace("|", "/")] for line in lines])


def tables(folder):
    folder = Path(folder)
    rows, correlations = load(folder / "antibiotics.tsv"), load(folder / "correlations.tsv")
    guards = (folder / "guards.txt").read_text(encoding="utf-8").splitlines()
    return {"antibiotics": print_antibiotics(rows), "correlations": print_correlations(correlations),
            "sets": print_sets(correlations), "guards": print_guards(guards)}


def main(argv):
    try:
        if argv == ["plan"]:
            print(RULE)
            return 0
        if len(argv) == 4 and argv[0] == "read":
            print("\n".join(read(*argv[1:])))
            return 0
        if len(argv) == 5 and argv[0] == "write":
            lines, correlations = write(*argv[1:])
            print("\n".join(lines))
            rule = [c for c in correlations if c["set"] == "rule"][0]
            print(f"guards: {len(GUARDS)} of {len(GUARDS)} pass; the rule's comparison gives: {rule['comparison_gives']}")
            return 0
        if len(argv) == 2 and argv[0] == "print":
            for name, lines in tables(argv[1]).items():
                print(f"[{name}]")
                print("\n".join(lines))
                print()
            return 0
    except Stop as e:
        print(f"stopped: {e}", file=sys.stderr)
        return 1
    sys.stderr.write(__doc__)
    return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
