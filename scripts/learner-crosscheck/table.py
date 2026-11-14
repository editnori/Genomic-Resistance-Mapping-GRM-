"""table.py DIRECTORY

The result table of the learner cross-check: its shape, and the counts that
follow from it. It needs no library, so the test suite reads the committed
table (reference/learner-crosscheck) with it.

`summary.tsv` has one row per panel and check. `differences.tsv` has one row
per difference found. `totals.json` holds, per check, the count of panels by
outcome and the sum of the units compared. Run as a script it prints the
problems of the tables in DIRECTORY and exits 1 if there is one.

An outcome is one of:

- matched: the check's rule in the plan held on the panel;
- explained: there is a difference and each has a cause that is shown;
- unexplained: there is a difference without a shown cause, a finding;
- not applicable: the comparison cannot run on the panel.
"""

import csv
import json
import pathlib
import sys

# (learner, check, the unit that `compared` counts)
CHECKS = (
    ("scm", "fits", "fits"),
    ("cart", "nodes", "nodes"),
    ("rf", "trees", "trees"),
    ("gbt", "tree by tree", "trees"),
    ("gbt", "end to end", "panels"),
)
OUTCOMES = ("matched", "explained", "unexplained", "not applicable")
SUMMARY_COLUMNS = ("panel", "samples", "features", "fit_samples", "learner", "check", "outcome", "units", "compared", "matched", "explained", "unexplained")
DIFFERENCE_COLUMNS = ("panel", "learner", "check", "unit", "what", "cause", "evidence")
COUNTS = ("samples", "features", "fit_samples", "compared", "matched", "explained", "unexplained")
# One difference row is one unit here; the tree checks also write rows for scores and free trees.
ROW_PER_UNIT = (("scm", "fits"), ("gbt", "tree by tree"), ("gbt", "end to end"))


def write(path, columns, rows):
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, columns, delimiter="\t", lineterminator="\n", quoting=csv.QUOTE_NONE, escapechar=None)
        writer.writeheader()
        for row in rows:
            for column in columns:
                if any(c in str(row[column]) for c in "\t\n\r"):
                    raise ValueError(f"a tab or line break in {column}: {row[column]!r}")
            writer.writerow({column: row[column] for column in columns})


def read(path, columns):
    with open(path, encoding="utf-8", newline="") as handle:
        lines = handle.read().split("\n")
    if lines[-1] != "":
        raise ValueError(f"{path}: no line end at the end")
    if tuple(lines[0].split("\t")) != tuple(columns):
        raise ValueError(f"{path}: the header is not {' '.join(columns)}")
    rows = []
    for number, line in enumerate(lines[1:-1], 2):
        fields = line.split("\t")
        if len(fields) != len(columns):
            raise ValueError(f"{path} line {number}: {len(fields)} fields")
        rows.append(dict(zip(columns, fields)))
    return rows


def read_summary(path):
    rows = read(path, SUMMARY_COLUMNS)
    for row in rows:
        for column in COUNTS:
            if not row[column].isdigit():
                raise ValueError(f"{path}: {column} of {row['panel']} is no count")
            row[column] = int(row[column])
    return rows


def read_differences(path):
    return read(path, DIFFERENCE_COLUMNS)


def totals(summary):
    """Per check: panels by outcome, and the sums of the unit counts."""
    out = []
    for learner, check, units in CHECKS:
        rows = [row for row in summary if (row["learner"], row["check"]) == (learner, check)]
        out.append(
            {
                "learner": learner,
                "check": check,
                "units": units,
                "panels": {outcome: sum(row["outcome"] == outcome for row in rows) for outcome in OUTCOMES},
                "compared": sum(row["compared"] for row in rows),
                "units_matched": sum(row["matched"] for row in rows),
                "units_explained": sum(row["explained"] for row in rows),
                "units_unexplained": sum(row["unexplained"] for row in rows),
            }
        )
    return out


def problems(summary, differences, panels=None):
    """What is wrong with the two tables' shape and counts; an empty list if nothing is."""
    found = []
    names = []
    for row in summary:
        if row["panel"] not in names:
            names.append(row["panel"])
    if panels is not None and names != list(panels):
        found.append("the panels are not the fixture file's, in its order")
    wanted = [(name, learner, check) for name in names for learner, check, _ in CHECKS]
    if [(row["panel"], row["learner"], row["check"]) for row in summary] != wanted:
        found.append("summary.tsv has not one row per panel and check, in order")
    units = {(learner, check): unit for learner, check, unit in CHECKS}
    by_key = {}
    for row in differences:
        by_key.setdefault((row["panel"], row["learner"], row["check"]), []).append(row)
        if not row["unit"] or not row["what"]:
            found.append(f"a difference of {row['panel']} names no unit or says not what differs")
    known = {(row["panel"], row["learner"], row["check"]) for row in summary}
    for key in by_key:
        if key not in known:
            found.append(f"a difference names {key}, which summary.tsv has not")
    for row in summary:
        key = (row["panel"], row["learner"], row["check"])
        where = f"{row['panel']} {row['learner']} {row['check']}"
        if row["outcome"] not in OUTCOMES:
            found.append(f"{where}: unknown outcome {row['outcome']!r}")
            continue
        if units.get(key[1:]) != row["units"]:
            found.append(f"{where}: the unit is not {units.get(key[1:])}")
        if row["fit_samples"] > row["samples"] or row["samples"] == 0 or row["features"] == 0:
            found.append(f"{where}: sample or feature counts")
        if row["matched"] + row["explained"] + row["unexplained"] != row["compared"]:
            found.append(f"{where}: matched, explained and unexplained do not add up to compared")
        mine = by_key.get(key, [])
        with_cause = sum(bool(d["cause"]) for d in mine)
        if row["outcome"] == "matched" and (mine or row["matched"] != row["compared"]):
            found.append(f"{where}: matched, with a difference")
        if row["outcome"] == "explained" and (not mine or with_cause != len(mine) or row["unexplained"]):
            found.append(f"{where}: explained, without a cause for every difference")
        if row["outcome"] == "unexplained" and with_cause == len(mine):
            found.append(f"{where}: unexplained, and no difference lacks a cause")
        if row["outcome"] == "not applicable" and (mine or row["compared"]):
            found.append(f"{where}: not applicable, with units or differences")
        if key[1:] in ROW_PER_UNIT and (row["explained"], row["unexplained"]) != (with_cause, len(mine) - with_cause):
            found.append(f"{where}: the unit counts are not the counts of its difference rows")
    return found


def main(argv):
    if len(argv) != 2:
        sys.exit(__doc__.splitlines()[0])
    directory = pathlib.Path(argv[1])
    summary = read_summary(directory / "summary.tsv")
    found = problems(summary, read_differences(directory / "differences.tsv"))
    if json.loads((directory / "totals.json").read_text(encoding="utf-8")) != totals(summary):
        found.append("totals.json is not what summary.tsv adds up to")
    for line in found:
        print(line)
    sys.exit(1 if found else 0)


if __name__ == "__main__":
    main(sys.argv)
