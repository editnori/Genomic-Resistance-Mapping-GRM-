"""overlap_faults.py: the faults of the learner overlap's table script
(learner_overlap.py, O1 to O17 and O29 to O31) and of its second pass
(overlap_second.py, O18 to O28), each planted in turn, and the SYNTHETIC test
that must catch it.

    python3 -I scripts/mechanism/overlap_faults.py [--out faults.tsv]

Run from the repository root. For each fault the named test runs once as the
script stands, where it must pass, and once with the fault planted, where it
must fail. Each fault changes one line of the script's logic. Exit 1 unless
every fault is caught. The faults are planted in memory only.

Python 3.12, standard library only; run with -I.
"""

import argparse
import fractions
import os
import sys
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
for p in (HERE, ROOT):
    if p not in sys.path:
        sys.path.insert(0, p)

from tests.mechanism import common as fixtures  # noqa: E402

planted = fixtures.load("planted_faults")
lo = fixtures.load("learner_overlap")
os2 = fixtures.load("overlap_second")

T = "tests.mechanism.test_learner_overlap."
T2 = "tests.mechanism.test_overlap_second."


def weighted_ratios(rows):
    # fault O7: each report's ratio counted once per genome
    return [x for r in rows for x in [lo.ratio_of(r)[1]] if x is not None for _ in range(int(r["genomes"]))]


def weighted_nonull_ratios(rows):
    # fault O25: each report's ratio outside null folds counted once per genome
    return [x for r in rows for x in [os2.nonull_of(r)[1]] if x is not None for _ in range(int(r["genomes"]))]


def yardstick_of_cells(rows):
    # fault O29: the pooled yardstick summed from the rounded cells in place of the counts
    return sum(fractions.Fraction(r["expected_all5"]) for r in rows)


def dec_of_cell(x, places):
    # fault O30: a page figure printed from the record's 3-place rounding of it, rounded a second time
    return "none" if x is None or x == "" else f"{float(lo.rec(x)):,.{places}f}"


def every_cell_skipped(path, decimal_columns):
    # fault O31: check 6 leaving every cell of a record table to the id list, not only the computed decimal columns
    with open(path, encoding="utf-8") as f:
        text = f.read()
    return len(lo.BARRED_KEY.findall(text)), 0


def majority(scores):
    # fault O19: three or more of the five calls in place of their mean
    return sum(1 for x in scores if x >= lo.CUTOFF) >= 3


def every_first_row(first_record):
    # fault O26: every report of the first record pooled, whatever its resistant count
    return set(os2.first_rows(first_record))


def every_fold_null(lr):
    # fault O16: every fold of a learner with one null fit read as null
    return set(range(len(lr.get("fold_complexity") or []))) if 0 in (lr.get("fold_complexity") or []) else set()


FAULTS = [
    ("O1", "a score of exactly 0.5 read as a miss", T + "Counts.test_a_score_of_exactly_the_cutoff_is_a_call",
     lambda: mock.patch.object(lo, "is_call", lambda score: score > lo.CUTOFF)),
    ("O2", "the random split read in place of the lineage split", T + "Counts.test_only_the_lineage_split_is_read",
     lambda: mock.patch.object(lo, "SCHEME", "random")),
    ("O3", "the catalogue reader counted among the learners",
     T + "Counts.test_the_catalogue_reader_is_not_a_learner",
     lambda: mock.patch.object(lo, "SEPARATE", ())),
    ("O4", "a score list of another length let through",
     T + "Stops.test_a_score_list_of_another_length_stops_the_read",
     lambda: mock.patch.object(lo, "scores_fit", lambda scores, n: True)),
    ("O5", "a report with one learner kept", T + "Skips.test_no_lineage_split_and_one_learner_are_skipped_and_counted",
     lambda: mock.patch.object(lo, "MIN_LEARNERS", 1)),
    ("O6", "the yardstick from the mean miss rate", T + "Counts.test_all_five_miss_the_same_genomes",
     lambda: mock.patch.object(lo, "expected_all5",
                               lambda resistant, misses: resistant * (sum(misses) / len(misses) / resistant) ** 5)),
    ("O7", "the median weighted by genome count", T + "Pooled.test_the_median_weights_each_report_once",
     lambda: mock.patch.object(lo, "ratio_values", weighted_ratios)),
    ("O8", "a report with a yardstick of 0 dropped from the pool",
     T + "Pooled.test_a_report_with_a_yardstick_of_0_is_reported_apart",
     lambda: mock.patch.object(lo, "pooled_in", lambda row: row["learner_count"] == 5 and row["miss_all5"] != ""
                               and lo.ratio_of(row)[1] is not None and row["resistant"] >= lo.MIN_RESISTANT)),
    ("O9", "an all-five miss the catalogue reader calls counted as uncalled",
     T + "Counts.test_misses_with_no_catalogue_call",
     lambda: mock.patch.object(lo, "uncalled", lambda catalog_call: catalog_call)),
    ("O10", "a genome id kept as a report is parsed", T + "Ids.test_no_id_is_held_after_parsing",
     lambda: mock.patch.object(lo, "drop_ids", lambda pairs: dict(pairs))),
    ("O11", "check 6 blind to an id-shaped token", T + "Ids.test_an_id_shaped_token_in_the_reading_stops_the_write",
     lambda: mock.patch.object(lo, "id_hits", lambda paths: 0)),
    ("O12", "counts that differ from the samples let through",
     T + "Stops.test_counts_that_differ_from_the_samples_stop_the_read",
     lambda: mock.patch.object(lo, "counts_agree", lambda counts, labels: True)),
    ("O13", "a report under 10 resistant genomes pooled",
     T + "Pooled.test_a_report_under_10_resistant_stays_out_of_the_pool",
     lambda: mock.patch.object(lo, "MIN_RESISTANT", 1)),
    ("O14", "stored metrics that differ from the recount let through",
     T + "Stops.test_stored_metrics_that_differ_from_the_recount_stop_the_read",
     lambda: mock.patch.object(lo, "metrics_agree", lambda labels, lr: True)),
    ("O15", "catalogue calls that differ from catalog_resistant let through",
     T + "Stops.test_catalogue_calls_that_differ_from_catalog_resistant_stop_the_read",
     lambda: mock.patch.object(lo, "catalog_agrees", lambda lr, catalog: True)),
    ("O16", "every fold of a learner with a null fit read as null", T + "Counts.test_misses_in_a_null_fit_fold",
     lambda: mock.patch.object(lo, "null_folds", every_fold_null)),
    ("O17", "a null catalogue call read as a call",
     T + "Counts.test_a_report_with_no_catalogue_rule_is_counted_and_marked",
     lambda: mock.patch.object(lo, "uncalled", lambda catalog_call: catalog_call is False)),
    ("O18", "a score sum of exactly 2.5 read as a miss by the mean", T2 + "Vote.test_a_sum_of_exactly_2_5_is_a_call",
     lambda: mock.patch.object(os2, "vote", lambda scores: os2.math.fsum(scores) > os2.VOTE_SUM)),
    ("O19", "the majority of the five read in place of their mean", T2 + "Vote.test_the_mean_is_not_the_majority",
     lambda: mock.patch.object(os2, "vote", majority)),
    ("O20", "a plain sum over 5 read in place of math.fsum",
     T2 + "Vote.test_the_vote_edge_is_counted_and_the_read_goes_on",
     lambda: mock.patch.object(os2, "vote", os2.plain_vote)),
    ("O21", "the null folds of one learner cut in place of all five's",
     T2 + "NullCut.test_null_folds_of_two_learners_are_both_cut",
     lambda: mock.patch.object(os2, "null_union", lambda learners: lo.null_folds(learners["gbt"]))),
    ("O22", "first-pass columns that differ from the record let through",
     T2 + "SecondStops.test_a_first_record_that_differs_stops_the_read",
     lambda: mock.patch.object(os2, "first_agrees", lambda row, first: True)),
    ("O23", "a vote that parts from the five let through",
     T2 + "Vote.test_a_vote_that_parts_from_the_five_stops_the_read",
     lambda: mock.patch.object(os2, "arithmetic_holds", lambda counts, row: True)),
    ("O24", "a null-fit cut that does not add up let through",
     T2 + "NullCut.test_a_cut_that_does_not_add_up_stops_the_read",
     lambda: mock.patch.object(os2, "null_cut_holds", lambda counts, row: True)),
    ("O25", "the median outside null folds weighted by genome count",
     T2 + "SecondPooled.test_the_median_outside_null_folds_weights_each_report_once",
     lambda: mock.patch.object(os2, "nonull_ratios", weighted_nonull_ratios)),
    ("O26", "a report under 10 resistant genomes pooled in the second pass",
     T2 + "SecondPooled.test_a_report_under_10_resistant_stays_out_of_the_pool",
     lambda: mock.patch.object(os2, "pooled_keys", every_first_row)),
    ("O27", "the second pass's shape scan blind to an id-shaped token",
     T2 + "SecondStops.test_an_id_shaped_token_in_the_reading_stops_the_write",
     lambda: mock.patch.object(os2, "scan", lambda paths: 0)),
    ("O28", "a genome 3 of the five miss counted with those all five miss",
     T2 + "Vote.test_the_mean_misses_what_three_or_more_of_the_five_miss_here",
     lambda: mock.patch.object(os2, "bucket", lambda k: "all5" if k >= 3 else ("none" if k == 0 else "1to4"))),
    ("O29", "the pooled yardstick summed from the rounded cells", T + "Places.test_the_pooled_yardstick_sums_the_counts",
     lambda: mock.patch.object(lo, "yardstick_sum", yardstick_of_cells)),
    ("O30", "a page figure rounded twice, through the record's 3 places",
     T + "Places.test_a_page_figure_is_rounded_once_from_the_counts",
     lambda: mock.patch.object(lo, "dec", dec_of_cell)),
    ("O31", "check 6 blind to an id-shaped token in a column that holds no computed decimal",
     T + "Places.test_check_6_scans_every_column_but_the_computed_decimals",
     lambda: mock.patch.object(lo, "record_hits", every_cell_skipped)),
]


def run_all():
    rows = []
    for n, name, test_id, plant in FAULTS:
        clean, with_fault = planted.outcome(test_id), planted.outcome(test_id, plant)
        rows.append({"fault": n, "planted": name, "test": test_id.rsplit(".", 2)[-2] + "." + test_id.rsplit(".", 1)[-1],
                     "as_written": clean, "with_fault": with_fault,
                     "caught": "yes" if clean == "pass" and with_fault == "fail" else "no"})
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--out")
    args = ap.parse_args(argv)
    rows = run_all()
    header = ["fault", "planted", "test", "as_written", "with_fault", "caught"]
    lines = ["\t".join(header)] + ["\t".join(str(r[h]) for h in header) for r in rows]
    if args.out:
        with open(args.out, "w", encoding="utf-8") as f:
            f.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    caught = sum(r["caught"] == "yes" for r in rows)
    print(f"planted faults: {caught} of {len(rows)} caught")
    return 0 if caught == len(rows) else 1


if __name__ == "__main__":
    sys.exit(main())
