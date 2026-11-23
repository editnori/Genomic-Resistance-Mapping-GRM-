"""page_faults.py: the faults of the result page's table script
(page_tables.py), each planted in turn, and the SYNTHETIC test that must
catch it.

    python3 -I scripts/mechanism/page_faults.py [--out faults.tsv]

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
page = fixtures.load("page_tables")

T = "tests.mechanism.test_page_tables."


def found_species(reading, lrows):
    # fault P2: the species are those found, not those of A
    return sorted(p.name for p in (reading / "runs").iterdir() if p.is_dir())


def lineage_reads_bioproject(r):
    # fault P3: a lineage guard row looks up the BioProject (study) row's figures
    guard = "bioproject" if r["guard"] == "lineage" else r["guard"]
    return (r["species_slug"], r["cell"], r["class"], r["dataset"], r["periods"], guard)


FAULTS = [
    ("P1", "a dataset of M counted in the primary reading", T + "Gates.test_a_dataset_of_m_marked_primary_is_refused",
     lambda: mock.patch.object(page, "in_m", lambda lrow: False)),
    ("P2", "a species dropped", T + "Gates.test_a_species_missing_from_the_reading_is_refused",
     lambda: mock.patch.object(page, "gate_species", found_species)),
    ("P3", "a study row read as a lineage row", T + "Tables.test_a_lineage_guard_row_keeps_its_own_figures",
     lambda: mock.patch.object(page, "m2_key", lineage_reads_bioproject)),
    ("P4", "a cell under its floor written as a share", T + "Gates.test_a_share_in_a_cell_under_its_floor_is_refused",
     lambda: mock.patch.object(page, "floor_hits", lambda path: {})),
    ("P5", "a guard read as the deciding figure", T + "Tables.test_the_deciding_count_reads_the_bioproject_vote_alone",
     lambda: mock.patch.object(page, "DECIDING", ("species_class", "primary", "largest_out"))),
    ("P6", "a first-run table mixed in", T + "Gates.test_a_first_run_copy_is_refused",
     lambda: mock.patch.object(page, "gate_runs", lambda reading, species, commit: None)),
    ("P7", "a species named by its slug, not in full", T + "Tables.test_species_and_datasets_are_named_in_full",
     lambda: mock.patch.object(page, "species_name", lambda lrows, slug: slug)),
]


def outcome(test_id, plant=None):
    """As planted_faults.outcome; a class fixture is built again under the plant."""
    return planted.outcome(test_id, plant)


def run_all():
    rows = []
    for n, name, test_id, plant in FAULTS:
        clean, with_fault = outcome(test_id), outcome(test_id, plant)
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
