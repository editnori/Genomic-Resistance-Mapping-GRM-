"""The counts and page tables of plan block C3's run on its datasets of M (docs/input-arms.md, The run on the datasets of M).

  python3 -I m_reading.py check <reference dir>
  python3 -I m_reading.py floor <out dir> <species folder> [<species folder> ...]
  python3 -I m_reading.py count <reference dir> <out dir> <species folder> [<species folder> ...]

check rebuilds the first run's dataset rows from its committed record,
<reference dir>/page-datasets.tsv (reference/input-arms), held to the
folder's SHA256SUMS: per dataset, split, arm and learner, balanced accuracy
and recall, the only measures the record holds. From them score.py's own
count, on those two measures, must give back every row of the published
count (page-count.tsv's count.tsv rows of those measures) and every row of
conditions 1 to 5 (page-conditions.tsv), at 12 decimals. It stops on any
difference.

floor reads each species folder's results/m/floor.tsv (m_list.py) and no
result and no score. Every dataset of M on lists.tsv must be there once,
and no first-run dataset. It writes to <out dir>:

- m-floor.tsv: every dataset of M with its counts by class in both
  readings, and whether its second reading is shown, counted and not
  shown, or under the floor;
- m-left-out.tsv: each dataset of M left out of the second reading, named in
  full, with its counts by class.

count runs check, then holds <out dir>'s m-floor.tsv and m-left-out.tsv to
the ones floor writes from the same folders, byte for byte, before it reads
any score. Then it reads each species folder of the run on M, as the run's
pack holds it: score/datasets.tsv (score.py, the first reading) and
score/second/datasets.tsv (second.py, the second reading). It takes three
counts in each reading, each its own call of score.py's count, so that
condition 2's Holm family of 30 is taken apart for each count and each
reading and nothing is corrected across them:

- the first run as read: the record's rows, the same in both readings, as
  a dataset outside M holds no genome with an intermediate call;
- M alone: the datasets of M, on score.py's six measures; it decides
  nothing (the entry Pick rule);
- all datasets: the record's rows and those of M, on balanced accuracy and
  recall only, the measures the record holds.

No count leaves out the constant-model datasets: the entry M run reads the
counts once, with no other subset. A dataset of M under the floor of 25 a
class is not in the second reading's counts. Every other dataset of M is in
every count of both readings (ruling R27).

Small lists (the entry M run, Small lists). In a dataset that lists 1
genome, the first reading less the second gives back that genome's class,
its calls and its scores. So a dataset's second-reading figures are written
to a page table only where each class lists no genome or at least
SMALLEST_LISTED (floor.tsv's listed counts by class). Every other dataset of
M that enters the second reading is counted and not shown: m-datasets.tsv
gives its first reading, its listed counts by class and that it is counted
and not shown, and no table gives its second-reading figures. A count row
of the second reading in which that dataset's figure would be the only one
on no page (the record's and the shown datasets' are on the page) gives its
datasets and its counts entered, and no better, worse, same, median,
interval or p value, so its figure cannot be taken back out by
subtraction; m-count.tsv and m-size.tsv say so in a column `figures`. The
count itself and its conditions are taken on every row.

The tables, written to <out dir>, hold counts per dataset at most and no
genome:

- m-count.tsv and m-conditions.tsv: each count and reading, in score.py's
  columns; m-count.tsv adds `datasets_not_shown` and `figures`;
- m-both.tsv: per arm, conditions 1 to 5 of ruling 15 on the count over all
  datasets in each reading, and whether they are met in both;
- m-size.tsv: per reading, arm (the twin's too) and learner, the datasets entered, better,
  worse and the same on the count over all datasets (lineage, balanced
  accuracy, by more than 0.02 at 12 decimals), the counts the entry Pick
  rule's size rule reads; it decides no pick;
- m-datasets.tsv: per reading, dataset of M, split, pair and learner, as
  page_tables.py's page-datasets.tsv, with the listed counts by class and
  the second reading's state; a second-reading row only for a dataset shown.

It prints the row counts only.
"""

import importlib.util
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
READINGS = ("first", "second")
# The column that names the M run's reading; a condition row's own `reading` is ruling 15 or the stricter one.
READ = "m_reading"
COUNTS = ("first run", "M", "all datasets")
RECORD_MEASURES = ("balanced_accuracy", "recall")
DECIDING = "ruling 15"
# A dataset's second-reading figures go in a page table only where each class lists no genome or at least this many
# (floor.tsv's listed counts by class). Five is a convention for a smallest cell (the entry M run, Small lists).
SMALLEST_LISTED = 5
SHOWN, NOT_SHOWN, UNDER_FLOOR = "shown", "counted and not shown", "under the floor"
# The cells of a count row that a dataset's figure moves: withheld where that figure is the only one on no page.
WITHHELD = ("better", "worse", "same", "median_difference", "interval_lo", "interval_hi", "sign_p", "holm_p")


class Stop(Exception):
    pass


def _load(name, file):
    spec = importlib.util.spec_from_file_location(name, HERE / file)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


SCORE = _load("c3_score_for_m_reading", "score.py")
PAGE = _load("c3_page_tables_for_m_reading", "page_tables.py")
CHECKS = _load("c3_checks_for_m_reading", "m_checks.py")
LIST = _load("c3_list_for_m_reading", "m_list.py")
read_tsv, write_tsv, fmt, r12 = SCORE.read_tsv, SCORE.write_tsv, SCORE.C9.fmt, SCORE.r12
DEFAULT, FAMILY, LEARNERS, PRIMARY, MARGIN = SCORE.DEFAULT, SCORE.FAMILY, SCORE.LEARNERS, SCORE.PRIMARY, SCORE.MARGIN


def scorer(measures):
    """score.py loaded afresh with its count over `measures`: count reads the module's MEASURES."""
    module = _load("c3_score_over_" + "_".join(measures), "score.py")
    module.MEASURES = list(measures)
    return module


def held(reference, names=("page-datasets.tsv", "page-count.tsv", "page-conditions.tsv")):
    try:
        sums = CHECKS.held_sums(reference)
    except CHECKS.Stop as e:
        raise Stop(str(e))
    for name in names:
        if sums.get(name) != CHECKS.sha256(Path(reference) / name):
            raise Stop(f"{name} differs from SHA256SUMS")


def record_rows(reference):
    """datasets.tsv rows, per dataset, split, arm and learner, rebuilt from page-datasets.tsv."""
    cells = {}

    def put(row):
        key = (row["dataset"], row["split"], row["arm"], row["learner"])
        if key in cells and cells[key] != row:
            raise Stop(f"page-datasets.tsv gives {key} two ways")
        cells[key] = row

    for r in read_tsv(Path(reference) / "page-datasets.tsv"):
        arm = r["pair"][: -len(f"-{DEFAULT}")]
        base = {"dataset": r["dataset"], "species_slug": r["dataset"].split("__")[0], "species": r["species"],
                "antibiotic": r["antibiotic"], "split": r["split"], "learner": r["learner"], "draws": r["draws"]}
        put(dict(base, arm=arm, scored_draws=r["arm_draws_scored"], constant_draws=r["constant_draws"],
                 scored=str(int(r["arm_draws_scored"] == r["draws"])), balanced_accuracy=r["arm_value"],
                 recall=r["arm_recall"]))
        put(dict(base, arm=DEFAULT, scored_draws=r["default_draws_scored"], constant_draws="0",
                 scored=str(int(r["default_draws_scored"] == r["draws"])), balanced_accuracy=r["default_value"],
                 recall=r["default_recall"]))
    return list(cells.values())


def same_cell(a, b):
    a, b = fmt(a), fmt(b)
    if a == b:
        return True
    try:
        return r12(float(a)) == r12(float(b))
    except ValueError:
        return False


def record_count(reference):
    """(the record's rows, score.py's count of them on balanced accuracy and recall, its conditions)."""
    held(reference)
    rows = record_rows(reference)
    two = scorer(RECORD_MEASURES)
    counted = two.count(rows)
    return rows, counted, two.conditions(counted)


def check(reference):
    """The first run's count and conditions given back from the record. Returns the record's rows."""
    rows, counted, conds = record_count(reference)
    return (rows,) + given_back(reference, counted, conds)


def given_back(reference, counted, conds):
    """(count rows, condition rows) compared: every published row given back at 12 decimals, or Stop."""
    key = lambda r: (r["split"], r["pair"], r["learner"], r["measure"])
    want = {key(r): r for r in read_tsv(Path(reference) / "page-count.tsv")
            if r["file"] == "count.tsv" and r["measure"] in RECORD_MEASURES}
    have = {key(r): r for r in counted}
    if set(have) != set(want):
        raise Stop(f"check: {len(set(want) - set(have))} published count rows not given back, "
                   f"{len(set(have) - set(want))} rows not published")
    for k in sorted(want):
        for c in SCORE.COUNT_COLUMNS:
            if not same_cell(have[k][c], want[k][c]):
                raise Stop(f"check: count row {k}, {c}: published {want[k][c]}, from the record {fmt(have[k][c])}")
    ckey = lambda r: (r["arm"], r["reading"], r["condition"])
    cwant = {ckey(r): r for r in read_tsv(Path(reference) / "page-conditions.tsv")}
    chave = {ckey(r): r for r in conds}
    if set(chave) != set(cwant):
        raise Stop("check: the conditions given back are not the published conditions' rows")
    for k in sorted(cwant):
        for c in SCORE.CONDITION_COLUMNS:
            if fmt(chave[k][c]) != cwant[k][c]:
                raise Stop(f"check: condition {k}, {c}: published {cwant[k][c]}, from the record {fmt(chave[k][c])}")
    return len(want), len(cwant)


def floor_tables(folders, lists_rows):
    """[(species folder, its floor.tsv rows)], from results/m/floor.tsv alone: no result and no score is read.
    Each dataset of M on lists.tsv is there once, and nothing outside M."""
    m = {r["dataset"] for r in lists_rows if r["in_M"] == "1"}
    out, seen = [], set()
    for top in folders:
        top = Path(top)
        path = top / "results" / "m" / "floor.tsv"
        if not path.is_file():
            raise Stop(f"{top.name}: no results/m/floor.tsv")
        f_rows = read_tsv(path)
        named = {r["dataset"] for r in f_rows}
        if named & seen or len(named) != len(f_rows):
            raise Stop(f"{top.name}: a dataset of M is read twice")
        seen |= named
        out.append((top, f_rows))
    if seen != m:
        raise Stop(f"the species folders hold {len(seen & m)} of the {len(m)} datasets of M and "
                   f"{len(seen - m)} datasets outside M")
    return out


def species_rows(folders, lists_rows):
    """({reading: M dataset rows}, floor rows), each dataset of M on lists.tsv there once."""
    out, floor = {k: [] for k in READINGS}, []
    for top, f_rows in floor_tables(folders, lists_rows):
        first = read_tsv(top / "score" / "datasets.tsv")
        second = read_tsv(top / "score" / "second" / "datasets.tsv")
        if {r["dataset"] for r in first} != {r["dataset"] for r in f_rows}:
            raise Stop(f"{top.name}: the first reading's datasets are not floor.tsv's")
        enters = {r["dataset"] for r in f_rows if r["enters"] == "1"}
        if {r["dataset"] for r in second} != enters:
            raise Stop(f"{top.name}: the second reading's datasets are not those floor.tsv enters")
        out["first"] += first
        out["second"] += second
        floor += f_rows
    return out, floor


def shown(row):
    """True when a dataset's second-reading figures may go in a page table: each class lists no genome or at least
    SMALLEST_LISTED, by floor.tsv's listed counts by class."""
    return all(int(row[c]) == 0 or int(row[c]) >= SMALLEST_LISTED for c in ("listed_resistant", "listed_susceptible"))


def second_state(row):
    """A floor row's second reading: shown, counted and not shown, or under the floor."""
    if row["enters"] != "1":
        return UNDER_FLOOR
    return SHOWN if shown(row) else NOT_SHOWN


def counts(record, m_rows):
    """{(count, reading): (count rows, condition rows)}, each its own count with its own Holm family."""
    two, six = scorer(RECORD_MEASURES), scorer(SCORE.MEASURES)
    clash = {(r["dataset"], r["split"], r["arm"], r["learner"]) for r in record} & {
        (r["dataset"], r["split"], r["arm"], r["learner"]) for rs in m_rows.values() for r in rs}
    if clash:
        raise Stop("a dataset of M is in the first run's record")
    out = {}
    for reading in READINGS:
        for name, module, rows in (("first run", two, record), ("M", six, m_rows[reading]),
                                   ("all datasets", two, record + m_rows[reading])):
            try:
                counted = module.count(rows)
            except SystemExit as e:
                raise Stop(f"the count {name} in the {reading} reading: {e}")
            out[(name, reading)] = (counted, module.conditions(counted))
    return out


def both_rows(taken):
    out = []
    first = {(r["arm"], r["condition"]): r for r in taken[("all datasets", "first")][1] if r["reading"] == DECIDING}
    second = {(r["arm"], r["condition"]): r for r in taken[("all datasets", "second")][1] if r["reading"] == DECIDING}
    for k in first:
        a, b = first[k]["met"], second[k]["met"]
        if a == "":
            continue
        out.append({"arm": k[0], "condition": k[1], "reads": first[k]["reads"], "met_first": a, "met_second": b,
                    "met_both": int(bool(a) and bool(b))})
    return out


def left_out_rows(floor):
    return [{c: r[c] for c in LEFT_OUT_COLUMNS} for r in floor if r["enters"] != "1"]


def floor_rows(floor):
    return [dict(r, second_reading=second_state(r)) for r in floor]


def m_dataset_rows(m_rows, floor):
    """page-datasets.tsv's rows for the datasets of M in each reading, with the listed counts by class and the second
    reading's state. A second-reading row only for a dataset whose second reading is shown."""
    by = {r["dataset"]: r for r in floor}
    out = []
    for reading in READINGS:
        rows = [r for r in m_rows[reading] if reading == "first" or second_state(by[r["dataset"]]) == SHOWN]
        slugs = sorted({r["species_slug"] for r in rows})
        data = type("Data", (), {"species": slugs,
                                 "score": {s: {"datasets.tsv": [r for r in rows if r["species_slug"] == s]}
                                           for s in slugs}})
        for r in PAGE.page_datasets(data):
            f = by[r["dataset"]]
            out.append(dict(r, listed_resistant=f["listed_resistant"], listed_susceptible=f["listed_susceptible"],
                            second_reading=second_state(f), **{READ: reading}))
    return out


def entered(rows):
    """{(split, pair, learner, measure): the datasets score.py's count enters}, by its own rule: the arm's row and the
    default's both scored, with both values."""
    cell = {(r["dataset"], r["split"], r["arm"], r["learner"]): r for r in rows}
    out = {}
    for (d, split, arm, learner), a in cell.items():
        b = cell.get((d, split, DEFAULT, learner))
        if arm == DEFAULT or b is None or not int(a["scored"]) or not int(b["scored"]):
            continue
        for m in SCORE.MEASURES:
            if a.get(m) not in (None, "") and b.get(m) not in (None, ""):
                out.setdefault((split, f"{arm}-{DEFAULT}", learner, m), set()).add(d)
    return out


def not_shown(record, m_rows, floor):
    """{(count, split, pair, learner, measure): (datasets counted and not shown, figures)} for the second reading.
    A row's figures are withheld where exactly one dataset's value in it is on no page and that dataset is counted
    and not shown: the record's datasets and the shown datasets of M are on the page on balanced accuracy and
    recall."""
    hidden = {r["dataset"] for r in floor if second_state(r) == NOT_SHOWN}
    on_page = {r["dataset"] for r in record} | {r["dataset"] for r in floor if second_state(r) == SHOWN}
    out = {}
    for name, rows in (("M", m_rows["second"]), ("all datasets", record + m_rows["second"])):
        for key, ds in entered(rows).items():
            off = ds - on_page if key[3] in RECORD_MEASURES else ds
            one = len(off) == 1 and off <= hidden
            out[(name,) + key] = (len(ds & hidden), NOT_SHOWN if one else SHOWN)
    return out


def count_rows(taken, hide):
    out = []
    for (c, g), (rows, _) in taken.items():
        for r in rows:
            n, figures = hide.get((c, r["split"], r["pair"], r["learner"], r["measure"]), (0, SHOWN)) \
                if g == "second" else (0, SHOWN)
            row = dict(r, count=c, datasets_not_shown=n, figures=figures, **{READ: g})
            if figures != SHOWN:
                row.update({k: None for k in WITHHELD})
            out.append(row)
    return out


def size_rows(counted):
    """From m-count.tsv's rows: the count over all datasets, lineage, the primary measure."""
    out = []
    for reading in READINGS:
        for r in (x for x in counted if x["count"] == "all datasets" and x[READ] == reading):
            if r["split"] == "lineage" and r["measure"] == PRIMARY:
                out.append({READ: reading, "arm": r["pair"][: -len(f"-{DEFAULT}")], "learner": r["learner"],
                            "datasets": r["datasets"], "entered": r["entered"], "better": r["better"],
                            "worse": r["worse"], "same": r["same"], "figures": r["figures"]})
    return out


COUNT_COLUMNS = ["count", READ] + SCORE.COUNT_COLUMNS + ["datasets_not_shown", "figures"]
CONDITION_COLUMNS = ["count", READ] + SCORE.CONDITION_COLUMNS
BOTH_COLUMNS = ["arm", "condition", "reads", "met_first", "met_second", "met_both"]
SIZE_COLUMNS = [READ, "arm", "learner", "datasets", "entered", "better", "worse", "same", "figures"]
LEFT_OUT_COLUMNS = ["dataset", "species", "antibiotic", "cohort", "resistant", "susceptible", "listed",
                    "listed_resistant", "listed_susceptible", "second_resistant", "second_susceptible", "floor"]
FLOOR_PAGE_COLUMNS = LIST.FLOOR_COLUMNS + ["second_reading"]
M_DATASET_COLUMNS = [READ] + PAGE.DATASET_TABLE_COLUMNS + ["listed_resistant", "listed_susceptible", "second_reading"]
FLOOR_FILES = ("m-floor.tsv", "m-left-out.tsv")


def write_floor(out, floor):
    """The floor form's two tables, from floor.tsv's rows alone."""
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    write_tsv(out / "m-floor.tsv", floor_rows(floor), FLOOR_PAGE_COLUMNS)
    write_tsv(out / "m-left-out.tsv", left_out_rows(floor), LEFT_OUT_COLUMNS)


def same_floor(out, floor):
    """Stop unless <out> holds m-floor.tsv and m-left-out.tsv as the floor form writes them from these rows."""
    with tempfile.TemporaryDirectory() as d:
        write_floor(d, floor)
        for name in FLOOR_FILES:
            have = Path(out) / name
            if not have.is_file():
                raise Stop(f"{name} is not in {out}: the floor form runs first")
            if have.read_bytes() != (Path(d) / name).read_bytes():
                raise Stop(f"{name} differs from the one the floor form writes from these species folders")


def write(out, taken, floor, m_rows, record):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    counted = count_rows(taken, not_shown(record, m_rows, floor))
    write_tsv(out / "m-count.tsv", counted, COUNT_COLUMNS)
    write_tsv(out / "m-conditions.tsv", [dict(r, count=c, **{READ: g}) for (c, g), (_, rows) in taken.items()
                                         for r in rows], CONDITION_COLUMNS)
    write_tsv(out / "m-both.tsv", both_rows(taken), BOTH_COLUMNS)
    write_tsv(out / "m-size.tsv", size_rows(counted), SIZE_COLUMNS)
    write_tsv(out / "m-datasets.tsv", m_dataset_rows(m_rows, floor), M_DATASET_COLUMNS)
    return counted


def states(floor):
    return "; ".join(f"{sum(1 for r in floor if second_state(r) == s)} {s}" for s in (SHOWN, NOT_SHOWN, UNDER_FLOOR))


def floor_form(out, folders, lists_rows):
    """The floor form: m-floor.tsv and m-left-out.tsv from the species' floor.tsv files alone. Returns its line."""
    floor = [r for _, rows in floor_tables(folders, lists_rows) for r in rows]
    write_floor(out, floor)
    return f"floor: datasets of M {len(floor)}; second reading: {states(floor)}"


def count_form(record, out, folders, lists_rows):
    """The count form after the first run's check gave back its record: the floor form's tables held to these
    folders' floor.tsv files before any score is read, then the three counts in both readings and their tables.
    Yields its lines."""
    same_floor(out, [r for _, rows in floor_tables(folders, lists_rows) for r in rows])
    yield "floor: m-floor.tsv and m-left-out.tsv are the floor form's"
    m_rows, floor = species_rows(folders, lists_rows)
    taken = counts(record, m_rows)
    counted = write(out, taken, floor, m_rows, record)
    withheld = sum(1 for r in counted if r["figures"] != SHOWN)
    yield (f"count: datasets of M {len(floor)}; second reading: {states(floor)}; "
           + "; ".join(f"{c}, {g} reading: {len(rows)} count rows" for (c, g), (rows, _) in taken.items())
           + f"; count rows counted and not shown: {withheld}")


def main(argv):
    try:
        if len(argv) == 2 and argv[0] == "check":
            _, n_count, n_cond = check(argv[1])
            print(f"check: {n_count} of {n_count} count rows and {n_cond} of {n_cond} condition rows of the first run "
                  "given back from page-datasets.tsv")
            return 0
        if len(argv) >= 3 and argv[0] == "floor":
            print(floor_form(argv[1], argv[2:], read_tsv(HERE / "lists.tsv")))
            return 0
        if len(argv) >= 4 and argv[0] == "count":
            record, n_count, n_cond = check(argv[1])
            print(f"check: {n_count} of {n_count} count rows and {n_cond} of {n_cond} condition rows given back")
            for line in count_form(record, argv[2], argv[3:], read_tsv(HERE / "lists.tsv")):
                print(line, flush=True)
            return 0
    except Stop as e:
        print(f"m_reading: stopped: {e}")
        return 1
    print("\n".join(__doc__.strip().splitlines()[2:5]), file=sys.stderr)
    return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
