"""The record and the page tables of plan block C3's run on its datasets of M (docs/input-arms.md, Page tables).

  python3 -I m_page.py floor <tables dir> <record dir>
  python3 -I m_page.py record <tables dir> <record dir>
  python3 -I m_page.py print <record dir>
  python3 -I m_page.py onpage <record dir> <page file>

<tables dir> is where m_reading.py wrote its tables on the count machine.
<record dir> is reference/input-arms-m, which holds README.md and
first-run-times.tsv before the run.

floor runs after m_reading.py floor and before any score is read. It reads
m-floor.tsv and m-left-out.tsv and nothing else of <tables dir>, not even
its listing, holds them to the gates below, copies them into <record dir>,
and writes SHA256SUMS over every file there but itself and README.md. It
runs once, on a record that holds README.md and first-run-times.tsv alone.

record runs after m_reading.py count. It stops unless <tables dir> holds
exactly m_reading.py's 7 tables, its two floor tables are byte for byte the
ones already in the record, and every gate passes on the 7. Then it copies
the other 5 in and writes SHA256SUMS again.

print prints each page table from <record dir> alone, as Markdown, each
under its name: the datasets left out of the second reading, the datasets
by species, the count rows (lineage balanced accuracy) and the rows beside
them, conditions 1 to 5 by ruling 15, the conditions met in both readings,
the size table and the datasets, in page_tables.py's number forms. With
the floor's two tables alone it prints the first two and says the count
has not been taken.

onpage exits 0 only when <page file> holds the left-out table as print
gives it, line after line. Otherwise it exits 1 with a line that names the
first row missing and no figure.

record, print and onpage stop on any gate, naming the file and the row and
no figure:

- a file missing, a file not on the record's list, a column list that is
  not m_reading.py's, or a file that differs from its SHA256SUMS line;
- a dataset of m-floor.tsv in none of the three states, or the datasets of
  m-floor.tsv not exactly those lists.tsv marks in_M, each once;
- a dataset in m-left-out.tsv whose reason is not its state in
  m-floor.tsv, or whose counts are not its row's there, or a dataset left
  out that m-left-out.tsv lacks;
- a second-reading row in m-datasets.tsv for a dataset whose state in
  m-floor.tsv is not in, or a dataset in the second reading whose rows in
  the two readings do not pair;
- a second-reading row of m-count.tsv or m-size.tsv that counts more
  datasets of M than m-floor.tsv gives the state in;
- a token of pack.py's GENOME_ID shape with 2 or more digits before the
  point, anywhere in a table. This is a shape check; the exact check
  against the species' ids runs on the run's machine.

Standard library only; it reads no score and no genome.
"""

import csv
import hashlib
import importlib.util
import sys
from decimal import ROUND_HALF_EVEN, Decimal
from pathlib import Path

HERE = Path(__file__).resolve().parent
README, SUMS, TIMES = "README.md", "SHA256SUMS", "first-run-times.tsv"
FLOOR_FILES = ("m-floor.tsv", "m-left-out.tsv")
COUNT_FILES = ("m-count.tsv", "m-conditions.tsv", "m-both.tsv", "m-size.tsv", "m-datasets.tsv")
TABLES = FLOOR_FILES + COUNT_FILES
LISTED = (README, SUMS, TIMES) + TABLES
READ = "m_reading"
READINGS = ("first", "second")
COUNTS = ("first run", "M", "all datasets")
IN, UNDER_FLOOR, SHORT_LIST = "in", "under the floor", "short list"
STATES = (IN, UNDER_FLOOR, SHORT_LIST)
LEFT_WORDS = {UNDER_FLOOR: "left out, under the floor", SHORT_LIST: "left out, short list"}
DEFAULT, PRIMARY, DECIDING = "k31", "balanced_accuracy", "ruling 15"
SHOWN_CONDITIONS = ("1 and 2", "3", "4", "5")

# m_reading.py's column lists, written out: tests/c3_inputs/test_m_page.py holds each to m_reading.py's own.
COLUMNS = {
    "m-floor.tsv": ["dataset", "species_slug", "species", "antibiotic", "cohort", "resistant", "susceptible", "listed",
                    "listed_resistant", "listed_susceptible", "second_resistant", "second_susceptible", "floor",
                    "enters", "second_reading"],
    "m-left-out.tsv": ["dataset", "species", "antibiotic", "reason", "cohort", "resistant", "susceptible", "listed",
                       "listed_resistant", "listed_susceptible", "second_resistant", "second_susceptible", "floor"],
    "m-count.tsv": ["count", READ, "split", "pair", "learner", "measure", "orientation", "withdrawn", "datasets",
                    "entered", "not_scored", "constant", "better", "worse", "same", "median_difference", "species",
                    "interval_lo", "interval_hi", "sign_p", "holm_p"],
    "m-conditions.tsv": ["count", READ, "arm", "reading", "condition", "reads", "met", "learners"],
    "m-both.tsv": ["arm", "condition", "reads", "met_first", "met_second", "met_both"],
    "m-size.tsv": [READ, "arm", "learner", "datasets", "entered", "better", "worse", "same"],
    "m-datasets.tsv": [READ, "dataset", "species", "antibiotic", "split", "pair", "learner", "draws",
                       "arm_draws_scored", "default_draws_scored", "constant_draws", "study_refused", "arm_value",
                       "default_value", "difference", "arm_recall", "default_recall", "recall_difference",
                       "listed_resistant", "listed_susceptible", "second_reading"],
}
# page_tables.py's words, written out: the test holds them to page_tables.py's own.
ARM_WORDS = {"k31": "K31, the default", "unitigs": "unitigs", "k31-counts": "K31 plus k-mer counts",
             "k31-aa": "K31 plus amino-acid k-mers", "program": "the AMRFinderPlus program's calls alone",
             "k31-program": "K31 plus the AMRFinderPlus program's calls",
             "twin": "the twin, the default at seeds 45, 46 and 47"}
MEASURE_WORDS = {"balanced_accuracy": "balanced accuracy", "recall": "recall", "specificity": "specificity",
                 "vme": "very major error rate", "me": "major error rate", "log_loss": "log loss"}
COUNT_WORDS = {"first run": "the first run as read", "M": "M alone", "all datasets": "all datasets"}


class Stop(Exception):
    pass


def _pack():
    spec = importlib.util.spec_from_file_location("c3_pack_for_m_page", HERE / "pack.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# pack.py's own pattern: digits, a point and digits, not inside a longer number and not after a letter.
GENOME_ID = _pack().GENOME_ID


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def m_datasets(lists_path):
    """The datasets lists.tsv marks in_M."""
    with open(lists_path, newline="", encoding="utf-8") as handle:
        return {r["dataset"] for r in csv.DictReader(handle, delimiter="\t") if r["in_M"] == "1"}


# ---------- reading and gates ----------

class Table:
    """One table as m_reading.py writes it: a header line, then one line per row, each cell between tabs."""

    def __init__(self, name, data):
        self.name, self.data = name, data
        try:
            self.text = data.decode("utf-8")
        except UnicodeDecodeError:
            raise Stop(f"{name}: not UTF-8 text") from None
        lines = self.text.split("\n")
        if lines and lines[-1] == "":
            lines.pop()
        self.lines = lines
        self.header = lines[0].split("\t") if lines else []
        self.rows = []
        for i, line in enumerate(lines[1:], 2):
            cells = line.split("\t")
            if len(cells) != len(self.header):
                raise Stop(f"{name} line {i}: its cells do not fit the header")
            self.rows.append((i, dict(zip(self.header, cells))))


def read_tables(directory, names):
    """{name: Table} of these files of <directory>, opened by name alone: the directory is not listed."""
    out = {}
    for name in names:
        path = Path(directory) / name
        if not path.is_file():
            raise Stop(f"{name} is missing from {directory}")
        out[name] = Table(name, path.read_bytes())
    return out


def gate_columns(tables):
    for name, t in tables.items():
        if t.header != COLUMNS[name]:
            raise Stop(f"{name} line 1: its columns are not m_reading.py's")


def gate_ids(tables):
    """Stop on a token of pack.py's GENOME_ID shape with 2 or more digits before the point, naming no token."""
    for name, t in tables.items():
        for i, line in enumerate(t.lines, 1):
            for m in GENOME_ID.finditer(line):
                if len(m.group(1).split(".")[0]) >= 2:
                    raise Stop(f"{name} line {i}: a token of a genome id's shape, digits before the point")


def gate_floor(floor, m_set):
    """m-floor.tsv: each row in one of the three states, and exactly the datasets lists.tsv marks in_M, each once."""
    seen = set()
    for i, r in floor.rows:
        d = r["dataset"]
        if r["second_reading"] not in STATES:
            raise Stop(f"m-floor.tsv line {i}: {d}: its second reading is none of in, under the floor, short list")
        if d in seen:
            raise Stop(f"m-floor.tsv line {i}: {d} a second time")
        if d not in m_set:
            raise Stop(f"m-floor.tsv line {i}: {d} is not a dataset lists.tsv marks in_M")
        seen.add(d)
    for d in sorted(m_set - seen):
        raise Stop(f"m-floor.tsv lacks {d}, a dataset lists.tsv marks in_M")


def gate_left_out(floor, left):
    """m-left-out.tsv: each dataset left out in m-floor.tsv, once, with its state as reason and its row's counts."""
    by = {r["dataset"]: r for _, r in floor.rows}
    seen = set()
    for i, r in left.rows:
        d = r["dataset"]
        f = by.get(d)
        if f is None:
            raise Stop(f"m-left-out.tsv line {i}: {d} is not in m-floor.tsv")
        if d in seen:
            raise Stop(f"m-left-out.tsv line {i}: {d} a second time")
        seen.add(d)
        if f["second_reading"] == IN:
            raise Stop(f"m-left-out.tsv line {i}: {d} is in the second reading by m-floor.tsv")
        if r["reason"] != f["second_reading"]:
            raise Stop(f"m-left-out.tsv line {i}: {d}: its reason is not its state in m-floor.tsv")
        if any(r[c] != f[c] for c in COLUMNS["m-left-out.tsv"] if c != "reason"):
            raise Stop(f"m-left-out.tsv line {i}: {d}: its counts are not its row's in m-floor.tsv")
    for i, f in floor.rows:
        if f["second_reading"] != IN and f["dataset"] not in seen:
            raise Stop(f"m-floor.tsv line {i}: {f['dataset']} is left out and m-left-out.tsv lacks it")


def gate_datasets(floor, datasets):
    """m-datasets.tsv: a second-reading row only for a dataset in that reading, and its rows paired across the two."""
    state = {r["dataset"]: r["second_reading"] for _, r in floor.rows}
    at = {g: {} for g in READINGS}
    for i, r in datasets.rows:
        d, g = r["dataset"], r[READ]
        if d not in state:
            raise Stop(f"m-datasets.tsv line {i}: {d} is not in m-floor.tsv")
        if g not in READINGS:
            raise Stop(f"m-datasets.tsv line {i}: {d}: its reading is neither first nor second")
        if g == "second" and state[d] != IN:
            raise Stop(f"m-datasets.tsv line {i}: a second-reading row of {d}, whose state in m-floor.tsv is not in")
        key = (d, r["split"], r["pair"], r["learner"])
        if key in at[g]:
            raise Stop(f"m-datasets.tsv line {i}: {d} {r['split']} {r['pair']} {r['learner']} a second time")
        at[g][key] = i
    for g, other in (("first", "second"), ("second", "first")):
        for key, i in at[g].items():
            if state[key[0]] == IN and key not in at[other]:
                raise Stop(f"m-datasets.tsv line {i}: {' '.join(key)} is in the second reading and has no "
                           f"{other}-reading row")


def whole(name, i, x):
    if not x.isdigit():
        raise Stop(f"{name} line {i}: a count of datasets that is not a whole number")
    return int(x)


def gate_counts(floor, count, size):
    """No second-reading row of m-count.tsv or m-size.tsv counts more datasets of M than m-floor.tsv gives the state
    in: the M count's datasets, and the others' less the first run's count of the first reading on the same row."""
    n_in = sum(1 for _, r in floor.rows if r["second_reading"] == IN)
    base = {(r["split"], r["pair"], r["learner"], r["measure"]): whole("m-count.tsv", i, r["datasets"])
            for i, r in count.rows if (r["count"], r[READ]) == ("first run", "first")}

    def m_part(name, i, r, key, of_m):
        if of_m:
            return max(whole(name, i, r["datasets"]), whole(name, i, r["entered"]))
        if key not in base:
            raise Stop(f"{name} line {i}: no row of the first run's count in the first reading beside it")
        return max(whole(name, i, r["datasets"]), whole(name, i, r["entered"])) - base[key]

    for i, r in count.rows:
        if r[READ] != "second":
            continue
        key = (r["split"], r["pair"], r["learner"], r["measure"])
        if m_part("m-count.tsv", i, r, key, r["count"] == "M") > n_in:
            raise Stop(f"m-count.tsv line {i}: the second reading's {r['count']} row for {' '.join(key)} counts more "
                       "datasets of M than m-floor.tsv gives the state in")
    for i, r in size.rows:
        if r[READ] != "second":
            continue
        key = ("lineage", f"{r['arm']}-{DEFAULT}", r["learner"], PRIMARY)
        if m_part("m-size.tsv", i, r, key, False) > n_in:
            raise Stop(f"m-size.tsv line {i}: the second reading's row for {r['arm']} {r['learner']} counts more "
                       "datasets of M than m-floor.tsv gives the state in")


def gate_tables(tables, m_set):
    """Every gate on the tables alone: the floor's two, and the count's five when they are there."""
    gate_columns(tables)
    gate_ids(tables)
    gate_floor(tables["m-floor.tsv"], m_set)
    gate_left_out(tables["m-floor.tsv"], tables["m-left-out.tsv"])
    if "m-count.tsv" in tables:
        gate_datasets(tables["m-floor.tsv"], tables["m-datasets.tsv"])
        gate_counts(tables["m-floor.tsv"], tables["m-count.tsv"], tables["m-size.tsv"])


def state(record):
    """The record's state by the tables it holds: none, floor (the floor's two alone) or count."""
    have = {p.name for p in Path(record).iterdir()}
    if not have & set(TABLES):
        return "none"
    return "count" if have & set(COUNT_FILES) else "floor"


def sums_text(record):
    files = sorted(p.name for p in Path(record).iterdir() if p.name not in (SUMS, README))
    return "".join(f"{sha256(Path(record) / name)}  {name}\n" for name in files)


def gate_files(record):
    """The record holds its list's files for its state, nothing else, each as its SHA256SUMS line gives it."""
    record = Path(record)
    for p in sorted(record.iterdir()):
        if p.name not in LISTED or not p.is_file():
            raise Stop(f"{p.name} is not on the record's list")
    want = (README, TIMES) + FLOOR_FILES + (COUNT_FILES if state(record) == "count" else ()) + (SUMS,)
    for name in want:
        if not (record / name).is_file():
            raise Stop(f"{name} is missing from the record")
    listed = {}
    for i, line in enumerate((record / SUMS).read_text(encoding="utf-8").splitlines(), 1):
        digest, _, name = line.partition("  ")
        if name in listed or name not in want or name in (SUMS, README):
            raise Stop(f"SHA256SUMS line {i}: {name} is no file of the record to sum")
        listed[name] = digest
    for name in want:
        if name in (SUMS, README):
            continue
        if name not in listed:
            raise Stop(f"{name} has no SHA256SUMS line")
        if sha256(record / name) != listed[name]:
            raise Stop(f"{name} differs from its SHA256SUMS line")


def gate(record, m_set):
    """Every gate on the record. Returns its tables; stops in the state none, where m-floor.tsv is missing."""
    gate_files(record)
    names = TABLES if state(record) == "count" else FLOOR_FILES
    tables = read_tables(record, names)
    gate_tables(tables, m_set)
    return tables


# ---------- the commands ----------

def floor(tables_dir, record, m_set):
    """Copy the floor's two tables into a record that holds README.md and first-run-times.tsv alone."""
    record = Path(record)
    for p in sorted(record.iterdir()):
        if p.name not in (README, TIMES):
            raise Stop(f"the record holds {p.name}: the floor's tables are copied once, into a record that holds "
                       "README.md and first-run-times.tsv alone")
    for name in (README, TIMES):
        if not (record / name).is_file():
            raise Stop(f"{name} is missing from the record")
    tables = read_tables(tables_dir, FLOOR_FILES)
    gate_tables(tables, m_set)
    for name in FLOOR_FILES:
        (record / name).write_bytes(tables[name].data)
    (record / SUMS).write_text(sums_text(record), encoding="utf-8")
    gate(record, m_set)
    return f"floor: {len(FLOOR_FILES)} tables copied; SHA256SUMS lists {len(sums_text(record).splitlines())} files"


def record_count(tables_dir, record, m_set):
    """After the count: the other 5 tables into a record that holds the floor's two, all 7 gated first."""
    tables_dir, record = Path(tables_dir), Path(record)
    have = sorted(p.name for p in tables_dir.iterdir())
    for name in have:
        if name not in TABLES:
            raise Stop(f"{name} in {tables_dir} is not one of m_reading.py's 7 tables")
    for name in TABLES:
        if name not in have:
            raise Stop(f"{name} is missing from {tables_dir}")
    if state(record) != "floor":
        raise Stop("the record does not hold the floor's two tables alone: floor runs first, and record once")
    recorded = gate(record, m_set)
    tables = read_tables(tables_dir, TABLES)
    for name in FLOOR_FILES:
        if tables[name].data != recorded[name].data:
            raise Stop(f"{name} in {tables_dir} differs from the record's")
    gate_tables(tables, m_set)
    for name in COUNT_FILES:
        (record / name).write_bytes(tables[name].data)
    (record / SUMS).write_text(sums_text(record), encoding="utf-8")
    gate(record, m_set)
    return f"record: {len(COUNT_FILES)} tables copied; SHA256SUMS lists {len(sums_text(record).splitlines())} files"


# ---------- page_tables.py's number forms ----------

def d3(x):
    if x in ("", None):
        return "none"
    value = Decimal(str(x)).quantize(Decimal("0.001"), rounding=ROUND_HALF_EVEN)
    if value == 0:
        value = Decimal("0.000")
    return str(value)


def pval(x):
    if x in ("", None):
        return "none"
    p = float(x)
    return "< 0.001" if p < 0.001 else d3(p)


def count_words(x):
    return f"{int(x):,}"


def of(part, whole_):
    return f"{count_words(part)} of {count_words(whole_)}"


def interval(lo, hi):
    return "none" if lo in ("", None) else f"{d3(lo)} to {d3(hi)}"


def yes(x):
    return {"1": "met", "0": "not met", "": "the maintainer's"}[str(x)]


def table(header, rows):
    lines = ["| " + " | ".join(header) + " |", "|" + "|".join("---" for _ in header) + "|"]
    lines += ["| " + " | ".join(str(c) for c in row) + " |" for row in rows]
    return lines


# ---------- the page's tables ----------

LEFT_OUT_HEADER = ["Species", "Antibiotic", "Second reading", "Cohort", "Listed", "Left after the list", "Floor"]
SPECIES_HEADER = ["Species", "Datasets of M", "In the second reading", "Left out, under the floor",
                  "Left out, short list"]
COUNT_HEADER = ["Count", "Reading", "Arm", "Learner", "Better", "Worse", "Same", "Median", "Interval", "Sign p",
                "Holm p", "Constant"]
BESIDE_HEADER = ["Count", "Reading", "Split", "Arm", "Learner", "Measure", "Better is", "Entered", "Better", "Worse",
                 "Same", "Median", "Interval"]
CONDITIONS_HEADER = ["Count", "Reading", "Arm", "Condition", "Read", "Learners"]
BOTH_HEADER = ["Arm", "Condition", "First reading", "Second reading", "Both"]
SIZE_HEADER = ["Reading", "Arm", "Learner", "Entered", "Better", "Worse", "Same"]
DATASETS_HEADER = ["Species", "Antibiotic", "Arm", "Learner", "First reading: arm", "First reading: K31",
                   "First reading: difference", "Second reading: arm", "Second reading: K31",
                   "Second reading: difference"]


def classes(r, s):
    return f"{count_words(r)} R, {count_words(s)} S"


def pair_words(pair):
    return ARM_WORDS[pair[: -len(f"-{DEFAULT}")]]


def print_left_out(floor, left):
    """Each dataset left out of the second reading, named in full, with its reason and its counts by class."""
    if not left.rows:
        return [f"0 of the {len(floor.rows)} datasets of M are left out of the second reading."]
    rows = sorted((r for _, r in left.rows), key=lambda r: (r["species"], r["antibiotic"], r["dataset"]))
    return table(LEFT_OUT_HEADER, [[f"*{r['species']}*", r["antibiotic"], LEFT_WORDS[r["reason"]],
                                    classes(r["resistant"], r["susceptible"]),
                                    classes(r["listed_resistant"], r["listed_susceptible"]),
                                    classes(r["second_resistant"], r["second_susceptible"]),
                                    f"{count_words(r['floor'])} of each class"] for r in rows])


def print_species(floor):
    """The datasets of M by species: in the second reading, left out under the floor, left out for a short list."""
    by = {}
    for _, r in floor.rows:
        by.setdefault(r["species"], dict.fromkeys(STATES, 0))[r["second_reading"]] += 1
    rows = []
    for species in sorted(by):
        n = sum(by[species].values())
        rows.append([f"*{species}*", count_words(n)] + [of(by[species][s], n) for s in STATES])
    n = len(floor.rows)
    rows.append([f"All {len(by)} species", count_words(n)]
                + [of(sum(v[s] for v in by.values()), n) for s in STATES])
    return table(SPECIES_HEADER, rows)


def reads_line(floor):
    """The line above each table that holds the second reading: how many of the datasets of M each reading reads."""
    n = len(floor.rows)
    k = {s: sum(1 for _, r in floor.rows if r["second_reading"] == s) for s in STATES}
    return (f"Datasets of M read: {of(n, n)} in the first reading, {of(k[IN], n)} in the second "
            f"({count_words(k[UNDER_FLOOR])} left out under the floor, {count_words(k[SHORT_LIST])} for a short "
            "list).")


def print_count(count):
    """The count rows on the primary measure, lineage split, for each of the three counts in each reading."""
    rows = [[COUNT_WORDS[r["count"]], r[READ], pair_words(r["pair"]), f"`{r['learner']}`",
             of(r["better"], r["entered"]), of(r["worse"], r["entered"]), of(r["same"], r["entered"]),
             d3(r["median_difference"]), interval(r["interval_lo"], r["interval_hi"]), pval(r["sign_p"]),
             pval(r["holm_p"]), r["constant"]]
            for _, r in count.rows if r["split"] == "lineage" and r["measure"] == PRIMARY]
    return table(COUNT_HEADER, rows)


def print_beside(count):
    """Every other count row but the study split's other measures, as page_tables.py's table beside the count."""
    rows = []
    for _, r in count.rows:
        if (r["split"], r["measure"]) == ("lineage", PRIMARY) or (r["split"] == "study" and r["measure"] != PRIMARY):
            continue
        rows.append([COUNT_WORDS[r["count"]], r[READ], r["split"], pair_words(r["pair"]), f"`{r['learner']}`",
                     MEASURE_WORDS[r["measure"]], r["orientation"], of(r["entered"], r["datasets"]), r["better"],
                     r["worse"], r["same"], d3(r["median_difference"]), interval(r["interval_lo"], r["interval_hi"])])
    return table(BESIDE_HEADER, rows)


def print_conditions(conditions):
    """Conditions 1 to 5 by ruling 15, which decides, for each count and reading."""
    rows = [[COUNT_WORDS[r["count"]], r[READ], ARM_WORDS[r["arm"]], r["condition"], yes(r["met"]), r["learners"]]
            for _, r in conditions.rows if r["reading"] == DECIDING and r["condition"] in SHOWN_CONDITIONS]
    return table(CONDITIONS_HEADER, rows)


def print_both(both):
    rows = [[ARM_WORDS[r["arm"]], r["condition"], yes(r["met_first"]), yes(r["met_second"]), yes(r["met_both"])]
            for _, r in both.rows]
    return table(BOTH_HEADER, rows)


def print_size(size):
    rows = [[r[READ], ARM_WORDS[r["arm"]], f"`{r['learner']}`", of(r["entered"], r["datasets"]),
             of(r["better"], r["entered"]), of(r["worse"], r["entered"]), of(r["same"], r["entered"])]
            for _, r in size.rows]
    return table(SIZE_HEADER, rows)


def print_datasets(floor, datasets):
    """One row per dataset of M, pair and learner on the lineage split: the first reading always, the second only
    for a dataset in it; a dataset left out reads its reason in each of the second reading's cells."""
    state = {r["dataset"]: r["second_reading"] for _, r in floor.rows}
    second = {(r["dataset"], r["pair"], r["learner"]): r for _, r in datasets.rows
              if r[READ] == "second" and r["split"] == "lineage"}
    rows = []
    for _, r in datasets.rows:
        if r[READ] != "first" or r["split"] != "lineage":
            continue
        cells = [f"*{r['species']}*", r["antibiotic"], pair_words(r["pair"]), f"`{r['learner']}`",
                 d3(r["arm_value"]), d3(r["default_value"]), d3(r["difference"])]
        if state[r["dataset"]] == IN:
            s = second[(r["dataset"], r["pair"], r["learner"])]
            cells += [d3(s["arm_value"]), d3(s["default_value"]), d3(s["difference"])]
        else:
            cells += [LEFT_WORDS[state[r["dataset"]]]] * 3
        rows.append(cells)
    return table(DATASETS_HEADER, rows)


NOT_TAKEN = "The count has not been taken: the record holds the floor's two tables alone."


def tables(record, m_set):
    """{name: lines} of the page's tables, from the record alone, after every gate."""
    t = gate(record, m_set)
    floor_t = t["m-floor.tsv"]
    out = {"left-out": print_left_out(floor_t, t["m-left-out.tsv"]), "species": print_species(floor_t)}
    if "m-count.tsv" not in t:
        out["count"] = [NOT_TAKEN]
        return out
    line = [reads_line(floor_t), ""]
    out["count"] = line + print_count(t["m-count.tsv"])
    out["beside"] = line + print_beside(t["m-count.tsv"])
    out["conditions"] = line + print_conditions(t["m-conditions.tsv"])
    out["both"] = line + print_both(t["m-both.tsv"])
    out["size"] = line + print_size(t["m-size.tsv"])
    out["datasets"] = line + print_datasets(floor_t, t["m-datasets.tsv"])
    return out


def count_headers():
    """The header lines of the tables print gives only after the count: none stands on a page before it."""
    return [table(h, [])[0] for h in (COUNT_HEADER, BESIDE_HEADER, CONDITIONS_HEADER, BOTH_HEADER, SIZE_HEADER,
                                      DATASETS_HEADER)]


def run_length(text, lines):
    """How many of `lines`, from the first, stand in `text` as whole lines one after another, at their best start."""
    page = text.split("\n")
    best = 0
    for start in range(len(page)):
        k = 0
        while k < len(lines) and start + k < len(page) and page[start + k] == lines[k]:
            k += 1
        best = max(best, k)
    return best


def onpage(record, page, m_set):
    """(True, line) when <page> holds the left-out table as print gives it; else (False, a line naming the first
    row missing, with no figure)."""
    t = gate(record, m_set)
    lines = print_left_out(t["m-floor.tsv"], t["m-left-out.tsv"])
    k = run_length(Path(page).read_text(encoding="utf-8"), lines)
    if k == len(lines):
        return True, f"onpage: {page} holds the left-out table as print gives it"
    if not t["m-left-out.tsv"].rows:
        what = "one line, that no dataset of M is left out"
    elif k == 0:
        what = "header row"
    elif k == 1:
        what = "rule row under its header"
    else:
        rows = sorted((r for _, r in t["m-left-out.tsv"].rows),
                      key=lambda r: (r["species"], r["antibiotic"], r["dataset"]))
        r = rows[k - 2]
        what = f"row of {r['species']} with {r['antibiotic']}"
    return False, f"onpage: {page} lacks the left-out table's {what}"


def main(argv, lists=HERE / "lists.tsv"):
    try:
        if len(argv) == 3 and argv[0] == "floor":
            print(floor(argv[1], argv[2], m_datasets(lists)))
            return 0
        if len(argv) == 3 and argv[0] == "record":
            print(record_count(argv[1], argv[2], m_datasets(lists)))
            return 0
        if len(argv) == 2 and argv[0] == "print":
            for name, lines in tables(argv[1], m_datasets(lists)).items():
                print(f"[{name}]")
                print("\n".join(lines))
                print()
            return 0
        if len(argv) == 3 and argv[0] == "onpage":
            ok, line = onpage(argv[1], argv[2], m_datasets(lists))
            print(line)
            return 0 if ok else 1
    except Stop as e:
        print(f"m_page: stopped: {e}")
        return 1
    print("\n".join(__doc__.strip().splitlines()[2:6]), file=sys.stderr)
    return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
