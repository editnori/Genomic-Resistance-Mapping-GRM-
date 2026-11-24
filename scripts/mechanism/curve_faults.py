"""curve_faults.py: the faults of the lineage curve's table script
(curve_tables.py), each planted in turn, and the SYNTHETIC test that must
catch it.

    python3 -I scripts/mechanism/curve_faults.py [--out faults.tsv]

Run from the repository root. For each fault the named test runs once as the
script stands, where it must pass, and once with the fault planted, where it
must fail. Exit 1 unless every fault is caught. The faults are planted in
memory only.

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
ct = fixtures.load("curve_tables")

T = "tests.mechanism.test_curve_tables."


def versions_only(c):
    # fault C9: the run check keeps its version check alone
    n = {k: int(c[k]) for k in ("results", "versions_named", "versions_other")}
    return [("every result file names the tool and database", n["versions_other"] == 0
             and n["versions_named"] == n["results"])]


FAULTS = [
    ("C1", "a share written over under 5 resistant genomes",
     T + "Gates.test_a_share_under_5_resistant_genomes_is_refused",
     lambda: mock.patch.object(ct, "share_hits", lambda tables, summary: 0)),
    ("C2", "a variant's own name in an output", T + "Gates.test_a_variants_own_name_in_an_output_is_refused",
     lambda: mock.patch.object(ct, "name_hits", lambda reading, genes: 0)),
    ("C3", "a gene seen in 1 to 4 resistant genomes named",
     T + "Gates.test_a_gene_in_1_to_4_resistant_genomes_named_is_refused",
     lambda: mock.patch.object(ct, "gene_hits", lambda tables: 0)),
    ("C4", "a dataset in M left unmarked", T + "Gates.test_a_dataset_of_m_left_unmarked_is_refused",
     lambda: mock.patch.object(ct, "in_m", lambda lrow: False)),
    ("C5", "a result file of another tool or database version read",
     T + "Gates.test_a_result_of_another_tool_version_is_refused",
     lambda: mock.patch.object(ct, "gate_versions", lambda checks, summary: None)),
    ("C6", "a planned genome with neither a result nor a line in failed.tsv passed",
     T + "Gates.test_a_planned_genome_with_neither_a_result_nor_a_failed_line_is_refused",
     lambda: mock.patch.object(ct, "gate_accounting", lambda checks, summary: None)),
    ("C7", "a curve written over under 5 resistant genomes", T + "Gates.test_a_curve_over_under_5_is_refused",
     lambda: mock.patch.object(ct, "share_hits", lambda tables, summary: 0)),
    ("C8", "the run check reads a result file of another tool version as the pinned one",
     T + "Check.test_a_result_of_another_tool_version_fails_the_check",
     lambda: mock.patch.object(ct, "names_pinned", lambda d, gid: True)),
    ("C9", "the run check passes a planned genome with neither a result nor a failed line",
     T + "Check.test_a_planned_genome_with_neither_fails_the_check",
     lambda: mock.patch.object(ct, "checks_pass", versions_only)),
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
