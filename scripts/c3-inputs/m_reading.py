"""The counts and page tables of plan block C3's run on its datasets of M (docs/input-arms.md, The run on the datasets of M).

  python3 -I m_reading.py check <reference dir>
  python3 -I m_reading.py count <reference dir> <out dir> <species folder> [<species folder> ...]

check rebuilds the first run's dataset rows from its committed record,
<reference dir>/page-datasets.tsv (reference/input-arms), held to the
folder's SHA256SUMS: per dataset, split, arm and learner, balanced accuracy
and recall, the only measures the record holds. From them score.py's own
count, on those two measures, must give back every row of the published
count (page-count.tsv's count.tsv rows of those measures) and every row of
conditions 1 to 5 (page-conditions.tsv), at 12 decimals. It stops on any
difference.

count runs check, then reads each species folder of the run on M, as the
run's pack holds it: score/datasets.tsv (score.py, the first reading),
score/second/datasets.tsv (second.py, the second reading) and
results/m/floor.tsv (m_list.py). Every dataset of M on lists.tsv must be
there once, and no first-run dataset. It then takes three counts in each
reading, each its own call of score.py's count, so that condition 2's Holm
family of 30 is taken apart for each count and each reading and nothing is
corrected across them:

- the first run as read: the record's rows, the same in both readings, as
  a dataset outside M holds no genome with an intermediate call;
- M alone: the datasets of M, on score.py's six measures; it decides
  nothing (the entry Pick rule);
- all datasets: the record's rows and those of M, on balanced accuracy and
  recall only, the measures the record holds.

No count leaves out the constant-model datasets: the entry M run reads the
counts once, with no other subset. A dataset of M under the floor of 25 a
class is not in the second reading's counts. The tables, written to
<out dir>, hold counts per dataset at most and no genome:

- m-count.tsv and m-conditions.tsv: each count and reading, in score.py's
  columns;
- m-both.tsv: per arm, conditions 1 to 5 of ruling 15 on the count over all
  datasets in each reading, and whether they are met in both;
- m-size.tsv: per reading, arm (the twin's too) and learner, the datasets entered, better,
  worse and the same on the count over all datasets (lineage, balanced
  accuracy, by more than 0.02 at 12 decimals), the counts the entry Pick
  rule's size rule reads; it decides no pick;
- m-floor.tsv: every dataset of M with its counts by class in both readings;
- m-left-out.tsv: each dataset of M left out of the second reading, named in
  full, with its counts by class;
- m-datasets.tsv: per reading, dataset of M, split, pair and learner, as
  page_tables.py's page-datasets.tsv.

It prints the row counts only.
"""

import importlib.util
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
READINGS = ("first", "second")
# The column that names the M run's reading; a condition row's own `reading` is ruling 15 or the stricter one.
READ = "m_reading"
COUNTS = ("first run", "M", "all datasets")
RECORD_MEASURES = ("balanced_accuracy", "recall")
DECIDING = "ruling 15"


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


def species_rows(folders, lists_rows):
    """({reading: M dataset rows}, floor rows), each dataset of M on lists.tsv there once."""
    m = {r["dataset"] for r in lists_rows if r["in_M"] == "1"}
    out, floor, seen = {k: [] for k in READINGS}, [], set()
    for top in folders:
        top = Path(top)
        f_rows = read_tsv(top / "results" / "m" / "floor.tsv")
        first = read_tsv(top / "score" / "datasets.tsv")
        second = read_tsv(top / "score" / "second" / "datasets.tsv")
        named = {r["dataset"] for r in f_rows}
        if {r["dataset"] for r in first} != named:
            raise Stop(f"{top.name}: the first reading's datasets are not floor.tsv's")
        enters = {r["dataset"] for r in f_rows if r["enters"] == "1"}
        if {r["dataset"] for r in second} != enters:
            raise Stop(f"{top.name}: the second reading's datasets are not those floor.tsv enters")
        if named & seen:
            raise Stop(f"{top.name}: a dataset of M is read twice")
        seen |= named
        out["first"] += first
        out["second"] += second
        floor += f_rows
    if seen != m:
        raise Stop(f"the species folders hold {len(seen & m)} of the {len(m)} datasets of M and "
                   f"{len(seen - m)} datasets outside M")
    return out, floor


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


def size_rows(taken):
    out = []
    for reading in READINGS:
        for r in taken[("all datasets", reading)][0]:
            if r["split"] == "lineage" and r["measure"] == PRIMARY:
                out.append({READ: reading, "arm": r["pair"][: -len(f"-{DEFAULT}")], "learner": r["learner"],
                            "datasets": r["datasets"], "entered": r["entered"], "better": r["better"],
                            "worse": r["worse"], "same": r["same"]})
    return out


def left_out_rows(floor):
    return [{c: r[c] for c in LEFT_OUT_COLUMNS} for r in floor if r["enters"] != "1"]


def m_dataset_rows(m_rows):
    out = []
    for reading in READINGS:
        slugs = sorted({r["species_slug"] for r in m_rows[reading]})
        data = type("Data", (), {"species": slugs,
                                 "score": {s: {"datasets.tsv": [r for r in m_rows[reading] if r["species_slug"] == s]}
                                           for s in slugs}})
        out += [dict(r, **{READ: reading}) for r in PAGE.page_datasets(data)]
    return out


COUNT_COLUMNS = ["count", READ] + SCORE.COUNT_COLUMNS
CONDITION_COLUMNS = ["count", READ] + SCORE.CONDITION_COLUMNS
BOTH_COLUMNS = ["arm", "condition", "reads", "met_first", "met_second", "met_both"]
SIZE_COLUMNS = [READ, "arm", "learner", "datasets", "entered", "better", "worse", "same"]
LEFT_OUT_COLUMNS = ["dataset", "species", "antibiotic", "cohort", "resistant", "susceptible", "listed",
                    "listed_resistant", "listed_susceptible", "second_resistant", "second_susceptible", "floor"]
M_DATASET_COLUMNS = [READ] + PAGE.DATASET_TABLE_COLUMNS


def write(out, taken, floor, m_rows):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    write_tsv(out / "m-count.tsv", [dict(r, count=c, **{READ: g}) for (c, g), (rows, _) in taken.items() for r in rows],
              COUNT_COLUMNS)
    write_tsv(out / "m-conditions.tsv", [dict(r, count=c, **{READ: g}) for (c, g), (_, rows) in taken.items()
                                         for r in rows], CONDITION_COLUMNS)
    write_tsv(out / "m-both.tsv", both_rows(taken), BOTH_COLUMNS)
    write_tsv(out / "m-size.tsv", size_rows(taken), SIZE_COLUMNS)
    write_tsv(out / "m-floor.tsv", floor, LIST.FLOOR_COLUMNS)
    write_tsv(out / "m-left-out.tsv", left_out_rows(floor), LEFT_OUT_COLUMNS)
    write_tsv(out / "m-datasets.tsv", m_dataset_rows(m_rows), M_DATASET_COLUMNS)


def main(argv):
    try:
        if len(argv) == 2 and argv[0] == "check":
            _, n_count, n_cond = check(argv[1])
            print(f"check: {n_count} of {n_count} count rows and {n_cond} of {n_cond} condition rows of the first run "
                  "given back from page-datasets.tsv")
            return 0
        if len(argv) >= 4 and argv[0] == "count":
            record, n_count, n_cond = check(argv[1])
            print(f"check: {n_count} of {n_count} count rows and {n_cond} of {n_cond} condition rows given back")
            m_rows, floor = species_rows(argv[3:], read_tsv(HERE / "lists.tsv"))
            taken = counts(record, m_rows)
            write(argv[2], taken, floor, m_rows)
            enters = sum(1 for r in floor if r["enters"] == "1")
            print(f"count: datasets of M {len(floor)}, {enters} in the second reading; "
                  + "; ".join(f"{c}, {g} reading: {len(rows)} count rows" for (c, g), (rows, _) in taken.items()))
            return 0
    except Stop as e:
        print(f"m_reading: stopped: {e}")
        return 1
    print("\n".join(__doc__.strip().splitlines()[2:4]), file=sys.stderr)
    return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
