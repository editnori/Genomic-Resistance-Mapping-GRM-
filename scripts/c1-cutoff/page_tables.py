#!/usr/bin/env python3
"""Print the tables of docs/cohort-result.md from the reading's files.

The files are the ones `read.py` wrote on its one run, kept under
reference/cohort-result/, with part A's line counts beside them. This program
counts nothing: every cell it prints is a cell of those files, with a value
rounded to 3 decimals for the page. tests/c1_cutoff checks that each line it
prints stands in the page.

    python3 -I scripts/c1-cutoff/page_tables.py reference/cohort-result
"""
import csv
import sys
from decimal import ROUND_HALF_EVEN, Decimal
from pathlib import Path

LEARNERS = ["scm", "cart", "rf", "gbt", "lasso"]
LOWER = ["10", "15", "20"]
TWIN = "25-s45"
SPLITS = [("lineage", "the lineage split"), ("random", "the random split"), ("study", "the study split")]
TRUTHS = [("all", "all scored genomes"), ("new_lineage", "new lineages")]
INTERVALS = [("levelled", "levelled"), ("fold_range", "fold range")]


def load(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def d3(text):
    value = Decimal(text).quantize(Decimal("0.001"), rounding=ROUND_HALF_EVEN)
    if value == 0:
        value = Decimal("0.000")
    return str(value)


def count(text):
    return f"{int(text):,}"


def of(part, whole):
    return f"{count(part)} of {count(whole)}"


def part_failed(reading):
    head = "does not meet: "
    return reading[len(head):] if reading.startswith(head) else reading


def table(header, rows):
    lines = ["| " + " | ".join(header) + " |", "|" + "|".join("---" for _ in header) + "|"]
    lines += ["| " + " | ".join(row) + " |" for row in rows]
    return lines


def learner_rows(rows, line):
    by = {r["learner"]: r for r in rows if r["line"] == line and r["row"] == "learner"}
    return [by[name] for name in LEARNERS]


def deciding(rows):
    out = []
    for line in LOWER:
        for r in learner_rows(rows, line):
            out.append([
                line, f"`{r['learner']}`", r["better"], r["worse"], r["same"],
                d3(r["median_difference"]), f"{d3(r['interval_lo'])} to {d3(r['interval_hi'])}",
                of(r["covered"], r["coverage_n"]), of(r["base_covered"], r["base_coverage_n"]),
                d3(r["coverage_minus_25"]), part_failed(r["reading"]),
            ])
    return table(["Line", "Learner", "Better", "Worse", "Same", "Median", "Interval", "At line", "At 25", "Gap", "Row"], out)


def twin(rows):
    out = []
    for r in learner_rows(rows, TWIN):
        out.append([
            f"`{r['learner']}`", r["better"], r["worse"], r["same"],
            d3(r["median_difference"]), f"{d3(r['interval_lo'])} to {d3(r['interval_hi'])}",
            of(r["covered"], r["coverage_n"]), of(r["base_covered"], r["base_coverage_n"]),
            d3(r["coverage_minus_25"]), part_failed(r["reading"]),
        ])
    return table(["Learner", "Seed better", "Seed worse", "Seed same", "Seed median", "Seed interval", "At seed 45", "At seed 42", "Seed gap", "Twin row"], out)


def fifty(rows):
    pick = [r for r in rows if r["split"] == "lineage" and r["truth_set"] == "all" and r["coverage"] == "levelled"]
    out = []
    for r in learner_rows(pick, "50"):
        out.append([
            f"`{r['learner']}`", r["better"], r["worse"], r["same"],
            d3(r["median_difference"]), f"{d3(r['interval_lo'])} to {d3(r['interval_hi'])}",
            of(r["covered"], r["coverage_n"]), of(r["base_covered"], r["base_coverage_n"]),
        ])
    return table(["Learner", "Better at 50", "Worse at 50", "Same at 50", "Median at 50", "Interval at 50", "Covered at 50", "Covered at 25"], out)


def widths(deciding_rows, twin_rows, other_rows):
    fifty_rows = [r for r in other_rows if r["split"] == "lineage" and r["truth_set"] == "all" and r["coverage"] == "levelled"]
    out = []
    for name in LEARNERS:
        cells = [f"`{name}`"]
        for line in LOWER:
            r = [x for x in learner_rows(deciding_rows, line) if x["learner"] == name][0]
            cells.append(d3(r["median_width"]))
        base = [x for x in learner_rows(deciding_rows, "20") if x["learner"] == name][0]
        cells.append(d3(base["base_median_width"]))
        r50 = [x for x in learner_rows(fifty_rows, "50") if x["learner"] == name][0]
        cells.append(d3(r50["median_width"]))
        out.append(cells)
    return table(["Learner", "Width at 10", "Width at 15", "Width at 20", "Width at 25", "Width at 50"], out)


def others(deciding_rows, twin_rows, other_rows):
    """One row per split, truth set and interval: learners meeting of 5 at each line, with the datasets entered."""
    def cell(rows, line):
        learners = [r for r in rows if r["line"] == line and r["row"] == "learner"]
        total = [r for r in rows if r["line"] == line and r["row"] != "learner"]
        entered = sorted({int(r["entered"]) for r in learners})
        datasets = sorted({int(r["datasets"]) for r in learners})
        span = str(entered[0]) if len(entered) == 1 else f"{entered[0]} to {entered[-1]}"
        return f"{total[0]['learners_meeting']} of 5 ({span} of {datasets[0]})"

    out = []
    for split, split_words in SPLITS:
        for truth, truth_words in TRUTHS:
            for interval, interval_words in INTERVALS:
                rows = [r for r in other_rows if r["split"] == split and r["truth_set"] == truth and r["coverage"] == interval]
                if (split, truth, interval) == ("lineage", "all", "levelled"):
                    rows = rows + deciding_rows + twin_rows
                out.append([split_words, truth_words, interval_words] + [cell(rows, line) for line in LOWER + [TWIN]])
    return table(["Split", "Truth", "Interval", "Line 10", "Line 15", "Line 20", "Twin"], out)


def measure(measures, name, header):
    out = []
    pick = [r for r in measures if r["split"] == "lineage" and r["truth_set"] == "all" and r["measure"] == name]
    for line in LOWER + ["25"]:
        for learner in LEARNERS:
            if line == "25":
                rows = [r for r in pick if r["line"] == "20" and r["side"] == "25" and r["learner"] == learner]
            else:
                rows = [r for r in pick if r["line"] == line and r["side"] == "line" and r["learner"] == learner]
            by = {r["who"]: r for r in rows}
            est, tru = by["estimate"], by["truth"]
            if name == "log_loss":
                out.append([
                    line, f"`{learner}`",
                    d3(est["pooled_value"]), count(est["pooled_denominator"]), d3(est["median_over_datasets"]),
                    d3(tru["pooled_value"]), count(tru["pooled_denominator"]), d3(tru["median_over_datasets"]),
                ])
            else:
                out.append([
                    line, f"`{learner}`",
                    f"{d3(est['pooled_value'])} ({of(est['pooled_count'], est['pooled_denominator'])})", d3(est["median_over_datasets"]),
                    f"{d3(tru['pooled_value'])} ({of(tru['pooled_count'], tru['pooled_denominator'])})", d3(tru["median_over_datasets"]),
                    of(tru["truth_inside_estimate_wilson"], tru["of_results"]),
                ])
    return table(header, out)


def part_a(lines):
    out = []
    for r in lines:
        out.append([
            r["min_class"], r["datasets"], r["added_to_the_277"], r["removed_from_the_277"],
            r["species_with_added"], r["species_new_at_this_line"], r["added_in_species_with_no_dataset_at_25"],
        ])
    return table(["Per class", "Datasets", "Added", "Removed", "Species", "New species", "Their datasets"], out)


def tables(folder):
    folder = Path(folder)
    deciding_rows = load(folder / "page-deciding.tsv")
    twin_rows = load(folder / "page-twin.tsv")
    other_rows = load(folder / "page-other.tsv")
    measures = load(folder / "page-measures.tsv")
    lines = load(folder / "part-a-lines.tsv")
    return {
        "deciding": deciding(deciding_rows),
        "twin": twin(twin_rows),
        "fifty": fifty(other_rows),
        "widths": widths(deciding_rows, twin_rows, other_rows),
        "others": others(deciding_rows, twin_rows, other_rows),
        "recall": measure(measures, "recall", ["Line", "Learner", "Estimate recall", "Estimate median", "Truth recall", "Truth median", "Inside"]),
        "specificity": measure(measures, "specificity", ["Line", "Learner", "Estimate specificity", "Median estimate", "Truth specificity", "Median truth", "Truth inside"]),
        "log_loss": measure(measures, "log_loss", ["Line", "Learner", "Estimate loss", "Estimate genomes", "Estimate's median", "Truth loss", "Truth genomes", "Truth's median"]),
        "part_a": part_a(lines),
    }


def main(argv):
    if len(argv) != 2:
        sys.stderr.write(__doc__)
        return 2
    for name, lines in tables(argv[1]).items():
        print(f"[{name}]")
        print("\n".join(lines))
        print()
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
