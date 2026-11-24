"""year_faults.py: the faults of the mechanism year read's script
(year_mix.py), each planted in turn, and the SYNTHETIC test that must catch
it (docs/mechanism-year-result.md).

    python3 -I scripts/mechanism/year_faults.py [--out faults.tsv]

Run from the repository root. Each fault replaces one rule of year_mix.py
with a one-line change, in memory only. The named test runs once as the
script stands, where it must pass, and once with the fault planted, where it
must fail. Exit 1 unless every fault is caught.

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
ym = fixtures.load("year_mix")

T = "tests.mechanism.test_year_mix."

FAULTS = [
    ("Y1", "a mean of exactly 0.5 not called", T + "Rules.test_a_mean_of_exactly_half_is_called",
     lambda: mock.patch.object(ym, "mean_called", lambda s: sum(s) / 5 > 0.5)),
    ("Y2", "the five scores added in order, so a mean of 0.5 can fall under it",
     T + "Rules.test_a_mean_of_exactly_half_is_called",
     lambda: mock.patch.object(ym, "mean_called", lambda s: sum(s) / 5 >= 0.5)),
    ("Y3", "one learner's score read as the vote", T + "Read.test_the_call_reads_all_five_learners",
     lambda: mock.patch.object(ym, "genome_scores", lambda rows, i: [rows["scm"]["scores"][i]] * 5)),
    ("Y4", "a year read from the first four characters", T + "Rules.test_a_year_not_four_digits_is_unknown",
     lambda: mock.patch.object(ym, "year_key", lambda t: t[:4] if isinstance(t, str) and len(t) >= 4 else ym.UNKNOWN)),
    ("Y5", "cells that differ from catalog.metrics let through",
     T + "Read.test_cells_that_disagree_with_metrics_stop_the_read",
     lambda: mock.patch.object(ym, "cell_disagreement", lambda report: None)),
    ("Y6", "a scores list of another length let through", T + "Read.test_a_scores_length_mismatch_stops",
     lambda: mock.patch.object(ym, "check_scores", lambda report, rows: None)),
    ("Y7", "year counts that do not add up to the total let through",
     T + "Read.test_year_counts_that_do_not_add_up_stop_the_write",
     lambda: mock.patch.object(ym, "check_sums", lambda by_year, totals, strata: None)),
    ("Y8", "a year of exactly 100 genomes flagged thin", T + "Read.test_a_year_under_100_genomes_is_thin",
     lambda: mock.patch.object(ym, "is_thin", lambda genomes: genomes <= 100)),
    ("Y9", "the rebuild held to gene only and neither alone, so a mismatch reads as a match",
     T + "Read.test_a_rebuild_that_does_not_match_falls_back_to_two_groups",
     lambda: mock.patch.object(ym, "strata_match",
                               lambda p, b: all(tuple(p[s]) == tuple(b[s]) for s in ("gene_only", "neither")))),
    ("Y10", "a dataset of M read as primary", T + "Read.test_datasets_outside_the_primary_set_are_skipped_and_counted",
     lambda: mock.patch.object(ym, "set_of", lambda r: "not_in_list" if r is None else
                               ("not_in_A" if r["in_A"] != "1" else "primary"))),
    ("Y11", "a stratum read from the label", T + "Read.test_the_rebuild_reads_groups_that_count_for_the_rule",
     lambda: mock.patch.object(ym, "rebuilt_stratum", lambda s, r, c: "gene_only" if s["resistant"] else "neither")),
    ("Y12", "every group counted for the rule", T + "Read.test_the_rebuild_reads_groups_that_count_for_the_rule",
     lambda: mock.patch.object(ym, "counting_groups", lambda r: set(range(len(r.get("catalog_groups") or []))))),
    ("Y13", "a null call put with no call", T + "Read.test_a_null_call_is_its_own_group",
     lambda: mock.patch.object(ym, "reader_group",
                               lambda s: "reader_call" if s.get("catalog_resistant") else "reader_no_call")),
    ("Y14", "a recall written over under 10 resistant genomes", T + "Read.test_recall_only_over_ten_resistant",
     lambda: mock.patch.object(ym, "recall_of", lambda c: ym.share(c["resistant_called_mean5"], c["resistant"]))),
    ("Y15", "the scan of the written files skipped", T + "Read.test_an_id_in_a_written_file_stops_the_write",
     lambda: mock.patch.object(ym, "scan_hits", lambda paths, ids, names: (0, 0))),
    ("Y16", "the species table's floor at 9 resistant genomes",
     T + "Read.test_the_species_table_keeps_strata_with_ten_or_more_resistant",
     lambda: mock.patch.object(ym, "SPECIES_FLOOR", 9)),
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
