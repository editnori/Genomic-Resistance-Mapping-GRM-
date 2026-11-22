#!/usr/bin/env python3
"""Print the tables of docs/refit-result.md from the reading's files.

The files are the ones `read_refit.py` wrote on its one run, kept under
reference/refit-result/. Every cell this program prints is a cell of those
files, with a share or an overlap rounded to 3 decimals for the page. One
table is counted here: D0 by species, the rows of page-d0-datasets.tsv per
species and cause, where a pair whose folds tie counts in each tied cause.
The table of tied pairs lists rows of that file and counts nothing.
tests/assumption_d0_c8 checks that each line it prints stands in the page.

    python3 -I scripts/assumption-d0-c8/page_tables.py reference/refit-result
"""
import csv
import sys
from decimal import ROUND_HALF_EVEN, Decimal
from pathlib import Path

LEARNERS = ["scm", "cart", "rf", "gbt", "lasso"]
CAUSES = ["absent", "removed", "lineage", "chosen", "unused"]
# The plan's row D0, Ordering, and the entry "D0 causes" for the fifth cause.
BLOCKS = {
    "absent": "C3",
    "removed": "C4's floor arm",
    "lineage": "C4's adjusted arms",
    "chosen": "none",
    "unused": "C5",
}
SCHEMES = [("lineage", "lineage"), ("random", "random")]
OVER_FIVE = "over_five"
NO_FOLD = "no missed fold"


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


def table(header, rows):
    lines = ["| " + " | ".join(header) + " |", "|" + "|".join("---" for _ in header) + "|"]
    lines += ["| " + " | ".join(row) + " |" for row in rows]
    return lines


def order(rows):
    by = {r["cause"]: r for r in rows}
    out = []
    for cause in CAUSES:
        r = by[cause]
        out.append([
            cause.capitalize(), of(r["miss_pairs_in_cause"], r["miss_pairs"]), d3(r["share_of_miss_pairs"]),
            r["largest_of_the_four"], BLOCKS[cause],
        ])
    return table(["Cause", "Miss pairs", "Share", "Largest", "Block"], out)


def learner_name(r):
    return f"`{r['learner']}`" if r["learner"] in LEARNERS else "All five"


def causes(rows):
    out = []
    for r in rows:
        out.append(
            [learner_name(r), r["miss_pairs"]]
            + [r[f"datasets_{cause}"] for cause in CAUSES]
            + [r["ties"], r["pairs_with_no_missed_fold"]]
        )
    return table(["Learner", "Pairs", "Absent", "Removed", "Lineage", "Chosen", "Unused", "Tied", "No fold"], out)


def shares(rows):
    out = []
    for r in rows:
        if r["learner"] not in LEARNERS:
            continue
        out.append(
            [learner_name(r)] + [d3(r[f"share_of_gap_to_the_amrfinderplus_program_{cause}"]) for cause in CAUSES]
        )
    return table(["Learner", "Absent gap", "Removed gap", "Lineage gap", "Chosen gap", "Unused gap"], out)


def by_species(rows):
    """Counted here: datasets, pairs and pairs per cause of each species in D0."""
    names = []
    for r in rows:
        if r["species"] not in names:
            names.append(r["species"])
    out = []
    for name in names:
        mine = [r for r in rows if r["species"] == name]
        cells = [f"*{name}*", str(len({r["dataset"] for r in mine})), str(len(mine))]
        for cause in CAUSES:
            cells.append(str(sum(1 for r in mine if r["cause"] != NO_FOLD and cause in r["cause"].split("+"))))
        out.append(cells)
    return table(
        ["Species", "D0 datasets", "D0 pairs", "Absent pairs", "Removed pairs", "Lineage pairs", "Chosen pairs", "Unused pairs"],
        out,
    )


def tied(rows):
    """The pairs whose missed folds tie between causes, and the pair with no missed fold."""
    out = []
    for r in rows:
        if "+" in r["cause"] or r["cause"] == NO_FOLD:
            words = NO_FOLD if r["cause"] == NO_FOLD else ", ".join(c.capitalize() for c in r["cause"].split("+"))
            out.append([f"*{r['species']}*", r["antibiotic"], f"`{r['learner']}`", words])
    return table(["Species", "Antibiotic", "Learner", "Causes"], out)


def scheme_rows(rows, scheme):
    by = {r["learner"]: r for r in rows if r["scheme"] == scheme}
    return [by[name] for name in LEARNERS + [OVER_FIVE]]


def overlap(rows):
    out = []
    for scheme, words in SCHEMES:
        for r in scheme_rows(rows, scheme):
            pattern = "none" if r["pattern_median"] == "NA" else d3(r["pattern_median"])
            out.append([
                words, learner_name(r), of(r["defined"], r["units"]), d3(r["q1"]), d3(r["median"]), d3(r["q3"]),
                r["units_fewer_than_two_folds_chose"], pattern,
            ])
    return table(["Split", "Learner", "Defined", "First quartile", "Median", "Third quartile", "Under two", "Pattern median"], out)


def family(rows):
    out = []
    for scheme, words in SCHEMES:
        for r in scheme_rows(rows, scheme):
            out.append([
                words, learner_name(r), count(r["distinct_summed"]),
                of(r["with_card_summed"], r["distinct_summed"]),
                of(r["with_amrfinderplus_catalog_family_summed"], r["with_card_summed"]),
                d3(r["catalog_family_share_of_carded"]),
            ])
    return table(["Folds' split", "Learner", "Distinct", "With card", "With family", "Family share"], out)


def checks(rows):
    out = []
    for r in rows:
        out.append([
            f"*{r['species']}*", r["datasets"], r["d0_datasets"], r["check_dataset_rows"], count(r["check_learner_rows"]),
            r["breakpoints_set_to_the_release"], r["feature_lists_fetched"],
            count(r["samples_with_other_descriptive_fields"]), r["cards_with_other_locus"],
        ])
    return table(["Species", "Datasets", "In D0", "Dataset rows", "Learner rows", "Breakpoints", "Lists", "Samples", "Cards"], out)


def names(rows):
    """One paragraph a species: its antibiotics in full, and those D0 reads."""
    out = []
    for r in rows:
        all_names = r["antibiotics"].split("; ")
        d0_names = [] if r["d0_antibiotics"] == "none" else r["d0_antibiotics"].split("; ")
        if not d0_names:
            tail = "D0 reads none of them."
        elif len(d0_names) == 1:
            tail = f"D0 reads 1 dataset: {d0_names[0]}."
        else:
            tail = f"D0 reads {len(d0_names)} datasets: {', '.join(d0_names)}."
        out.append(f"*{r['species']}* has {len(all_names)} datasets: {', '.join(all_names)}. {tail}")
    return out


def tables(folder):
    folder = Path(folder)
    d0_order = load(folder / "page-d0-order.tsv")
    d0_learners = load(folder / "page-d0-learners.tsv")
    d0_datasets = load(folder / "page-d0-datasets.tsv")
    c8 = load(folder / "page-c8-summary.tsv")
    species = load(folder / "page-species.tsv")
    return {
        "order": order(d0_order),
        "causes": causes(d0_learners),
        "shares": shares(d0_learners),
        "by_species": by_species(d0_datasets),
        "tied": tied(d0_datasets),
        "overlap": overlap(c8),
        "family": family(c8),
        "checks": checks(species),
        "names": names(species),
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
