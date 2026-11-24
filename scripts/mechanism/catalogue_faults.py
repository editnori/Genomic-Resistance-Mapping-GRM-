"""catalogue_faults.py: the faults of the catalogue size reading's script
(catalogue_size.py), each planted in turn, and the SYNTHETIC test that must
catch it.

    python3 -I scripts/mechanism/catalogue_faults.py [--out faults.tsv]

Run from the repository root. For each fault the named test runs once as the
script stands, where it must pass, and once with the fault planted, where it
must fail. Exit 1 unless every fault is caught. Each fault changes one line of
the script's logic, in memory only.

Python 3.12, standard library only; run with -I.
"""

import argparse
import os
import sys
from fractions import Fraction
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
for p in (HERE, ROOT):
    if p not in sys.path:
        sys.path.insert(0, p)

from tests.mechanism import common as fixtures  # noqa: E402

planted = fixtures.load("planted_faults")
cs = fixtures.load("catalogue_size")

T = "tests.mechanism.test_catalogue_size."


def ordinal_ranks(values):
    # fault K1: tied values take the ranks of their order, not their average
    out = [0.0] * len(values)
    for rank, i in enumerate(sorted(range(len(values)), key=lambda i: values[i]), 1):
        out[i] = float(rank)
    return out


def signed_branch(rho_resistant, rho_share):
    # fault K3: the rule compares the signed correlations
    r, s = rho_resistant, rho_share
    return cs.BRANCHES[0] if r >= s - cs.TOL else cs.BRANCHES[1] if s - r >= cs.MARGIN - cs.TOL else cs.BRANCHES[2]


def no_margin(rho_resistant, rho_share):
    # fault K4: the share branch whenever the share's correlation is the stronger
    r, s = abs(rho_resistant), abs(rho_share)
    return cs.BRANCHES[0] if r >= s - cs.TOL else cs.BRANCHES[1]


def tie_to_share(rho_resistant, rho_share):
    # fault K5: equal strength leaves the resistant-genome branch
    r, s = abs(rho_resistant), abs(rho_share)
    return cs.BRANCHES[0] if r > s else cs.BRANCHES[1] if s - r >= cs.MARGIN - cs.TOL else cs.BRANCHES[2]


def set_m_kept(r):
    # fault K6: a row of set M enters the correlations
    return cs.REASONS[0] if r["kind"] == "none" else cs.REASONS[1] if r["catalogue_variants"] == "0" else ""


def ungraded_kept(r):
    # fault K7: a row with no graded variant enters the correlations
    return cs.REASONS[0] if r["kind"] == "none" else cs.REASONS[2] if r["in_M"] == "1" else ""


def largest_is_smallest(x, y, names):
    # fault K14: the leave-one-out largest is taken as the smallest
    values = [(cs.r12(cs.spearman(x[:i] + x[i + 1:], y[:i] + y[i + 1:])), names[i]) for i in range(len(names))]
    low = min(values)
    return low[0], low[1], low[0], low[1]


FAULTS = [
    ("K1", "tied values ranked by their order, not by the average of their ranks",
     T + "Spearman.test_a_case_with_ties_gives_its_known_answer",
     lambda: mock.patch.object(cs, "ranks", ordinal_ranks)),
    ("K2", "Pearson's correlation over the values in place of their ranks",
     T + "Spearman.test_a_case_with_no_ties_gives_its_known_answer",
     lambda: mock.patch.object(cs, "spearman", cs.pearson)),
    ("K3", "the rule compares the signed correlations", T + "Rule.test_the_rule_reads_absolute_values",
     lambda: mock.patch.object(cs, "branch", signed_branch)),
    ("K4", "the margin of 0.15 dropped", T + "Rule.test_a_share_stronger_by_less_than_the_margin_selects_neither",
     lambda: mock.patch.object(cs, "branch", no_margin)),
    ("K5", "equal strength read as the share's",
     T + "Rule.test_a_resistant_correlation_as_strong_as_the_share_selects_resistant_genomes",
     lambda: mock.patch.object(cs, "branch", tie_to_share)),
    ("K6", "a row of set M in the correlations", T + "Write.test_set_m_is_in_the_table_and_out_of_the_correlations",
     lambda: mock.patch.object(cs, "left_out_because", set_m_kept)),
    ("K7", "a row with no graded variant in the correlations",
     T + "Write.test_an_antibiotic_with_no_graded_variant_is_left_out",
     lambda: mock.patch.object(cs, "left_out_because", ungraded_kept)),
    ("K8", "a page antibiotic with no release report let through",
     T + "Read.test_a_page_antibiotic_with_no_report_stops",
     lambda: mock.patch.object(cs, "gate_match", lambda page, reports: "1 match")),
    ("K9", "counts that do not add up let through", T + "Read.test_counts_that_do_not_add_up_stop",
     lambda: mock.patch.object(cs, "gate_sums", lambda counts: "2 sums")),
    ("K10", "a resistant count unlike the lineage curve page's let through",
     T + "Read.test_a_resistant_count_unlike_the_page_stops",
     lambda: mock.patch.object(cs, "gate_page_agrees", lambda rows, counts: "6 resistant")),
    ("K11", "a pair over another set of antibiotics let through", T + "Write.test_a_pair_over_another_set_stops",
     lambda: mock.patch.object(cs, "gate_sets", lambda correlations: [])),
    ("K12", "the hand-made cases skipped", T + "Write.test_a_wrong_rank_correlation_stops_before_any_correlation",
     lambda: mock.patch.object(cs, "prove_spearman", lambda: "4 rank correlation")),
    ("K13", "the scan for genome ids skipped", T + "Write.test_a_genome_id_in_an_output_stops",
     lambda: mock.patch.object(cs, "id_hits", lambda paths, ids: 0)),
    ("K14", "the leave-one-out largest taken as the smallest",
     T + "Spearman.test_leave_one_out_gives_the_smallest_and_largest_with_the_antibiotic_dropped",
     lambda: mock.patch.object(cs, "leave_one_out", largest_is_smallest)),
    ("K15", "a report's samples list kept when the report is parsed", T + "Read.test_the_samples_list_is_never_kept",
     lambda: mock.patch.object(cs, "drop_samples", dict)),
    ("K16", "the resistant share taken over the susceptible genomes",
     T + "Write.test_the_resistant_share_is_over_all_genomes",
     lambda: mock.patch.object(cs, "share_of", lambda resistant, genomes: Fraction(resistant, genomes - resistant))),
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
