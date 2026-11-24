"""overlap_faults.py: the faults of the learner overlap's table script
(learner_overlap.py), each planted in turn, and the SYNTHETIC test that must
catch it.

    python3 -I scripts/mechanism/overlap_faults.py [--out faults.tsv]

Run from the repository root. For each fault the named test runs once as the
script stands, where it must pass, and once with the fault planted, where it
must fail. Each fault changes one line of the script's logic. Exit 1 unless
every fault is caught. The faults are planted in memory only.

Python 3.12, standard library only; run with -I.
"""

import argparse
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

T = "tests.mechanism.test_learner_overlap."


def weighted_ratios(rows):
    # fault O7: each report's ratio counted once per genome
    return [float(r["ratio_all5"]) for r in rows if r["ratio_all5"] != "" for _ in range(int(r["genomes"]))]


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
     lambda: mock.patch.object(lo, "pooled_in", lambda row: row["learner_count"] == 5 and row["ratio_all5"] != ""
                               and row["resistant"] >= lo.MIN_RESISTANT)),
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
