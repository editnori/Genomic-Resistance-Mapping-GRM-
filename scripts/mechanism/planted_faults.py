"""planted_faults.py: the faults of docs/mechanism-split.md, each planted in
the scripts in turn, and the SYNTHETIC test that must catch it: the ten of
the design, and the first run's (23 Nov 2026), from the database fetched
before the program is installed to an unknown symbol let through.

    python3 -I scripts/mechanism/planted_faults.py [--out faults.tsv]

Run from the repository root. For each fault the named test runs once as the
scripts stand, where it must pass, and once with the fault planted, where it
must fail. The table gives each fault, its test and both outcomes. Exit 1
unless every fault is caught.

Python 3.12, standard library only; run with -I.
"""

import argparse
import contextlib
import io
import os
import sys
import unittest
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
for p in (HERE, ROOT):
    if p not in sys.path:
        sys.path.insert(0, p)

from tests.mechanism import common as fixtures  # noqa: E402
from tests.mechanism import test_tools  # noqa: E402

ar = fixtures.load("amr_rule")
strata = fixtures.load("strata")
history = fixtures.load("history")


def label_strata(samples, tables, rule):
    # fault 1: a susceptible genome is put in "neither" whatever its table holds
    return [None if tables.get(s["id"]) is None else
            (strata.stratum_of(tables[s["id"]], rule) if s["resistant"] else ("neither", False, False, False))
            for s in samples]


def first_split(report):
    # fault 5: every split reads the first split's learner rows
    rows = {lr["learner"]: lr for lr in report["splits"][0]["learners"]}
    return {s["scheme"]: rows for s in report["splits"]}


def imputed_year(text, _orig=history.year_of):
    # fault 8: a genome with no year gets 2015
    return _orig(text) or 2015


def genome_permutation(units, observed, seed_name, statistic, _orig=history.permutation_p):
    # fault 9: period labels shuffled over genomes, not over BioProjects
    return _orig([(u[0], (u[1], i), [v]) for u in units for i, v in enumerate(u[2])], observed, seed_name, statistic)


def install_after_db(text):
    # fault 11: the program's environment installed after the database fetch,
    # the order of the first run (defect A)
    start, end = text.index("# The program's environment, once,"), text.index("db() {\n")
    block, text = text[start:end], text[:start] + text[end:]
    at = text.index("step db db\n") + len("step db db\n")
    return text[:at] + "\n" + block + text[at:]


def first_child_family(self, symbol, subtype="AMR", _orig=ar.Database.family):
    # fault 12: a node above the gene families given its first child's family
    fam = _orig(self, symbol, subtype)
    if fam == symbol and symbol not in self.families:
        children = [n for n, p in self.parent.items() if p == symbol and n in self.families]
        return children[0] if children else fam
    return fam


def let_through(self, symbol, subtype="AMR", _orig=ar.Database.family):
    # fault 13: a symbol in neither file let through under its own name
    try:
        return _orig(self, symbol, subtype)
    except ar.Refusal:
        if subtype != "AMR":
            raise
        return symbol


FAULTS = [
    (1, "a stratum read from the label", "tests.mechanism.test_strata.Strata.test_stratum_comes_from_the_table_not_the_label",
     lambda: mock.patch.object(strata, "genome_strata", label_strata)),
    (2, "a POINT row counted as a gene row", "tests.mechanism.test_strata.Strata.test_point_rows_make_point_only",
     lambda: mock.patch.object(ar, "is_point", lambda row: False)),
    (3, "a partial call counted as a full one",
     "tests.mechanism.test_strata.Strata.test_partial_calls_and_point_disrupt_are_counted_apart",
     lambda: mock.patch.object(ar, "is_partial", lambda row: False)),
    (4, "a recall written under the cell floor",
     "tests.mechanism.test_strata.Strata.test_recall_only_over_ten_resistant_and_small_counts_hidden",
     lambda: mock.patch.object(strata, "recall_cell", lambda called, res: f"{called / res:.12f}" if res else "")),
    (5, "scores read from the wrong split or learner",
     "tests.mechanism.test_strata.Strata.test_scores_are_read_from_their_own_split_and_learner",
     lambda: mock.patch.object(strata, "split_rows", first_split)),
    (6, "a BioProject given one vote per genome", "tests.mechanism.test_history.Vote.test_one_vote_per_bioproject",
     lambda: mock.patch.object(history, "bioproject_unit", lambda genome: id(genome))),
    (7, "a period boundary off by one year", "tests.mechanism.test_history.Years.test_period_ends_are_inclusive_and_exact",
     lambda: mock.patch.object(history, "PERIODS", {"primary": ((2010, 2015), (2016, 2019)),
                                                    "secondary": ((2000, 2010), (2011, 2019))})),
    (8, "an imputed year", "tests.mechanism.test_history.Vote.test_no_year_is_imputed",
     lambda: mock.patch.object(history, "year_of", imputed_year)),
    (9, "a permutation that moves genomes and not BioProjects",
     "tests.mechanism.test_history.Vote.test_permutation_moves_bioprojects_not_genomes",
     lambda: mock.patch.object(history, "permutation_p", genome_permutation)),
    (10, "a family named other than as the database names it",
     "tests.mechanism.test_history.Vote.test_family_mix_names_the_databases_families",
     lambda: mock.patch.object(ar.Database, "family", lambda self, symbol, subtype="AMR": symbol)),
    (11, "the database fetched before the program is installed",
     "tests.mechanism.test_tools.InstallOrder.test_with_db_dash_the_install_ends_before_the_fetch",
     lambda: mock.patch.object(test_tools.InstallOrder, "text", install_after_db(test_tools.InstallOrder.text))),
    (12, "a node above the gene families given its first child's family",
     "tests.mechanism.test_history.NodeFamilies.test_a_node_above_the_families_is_its_own_family",
     lambda: mock.patch.object(ar.Database, "family", first_child_family)),
    (13, "a symbol in neither database file let through",
     "tests.mechanism.test_history.NodeFamilies.test_a_symbol_in_neither_file_still_refuses",
     lambda: mock.patch.object(ar.Database, "family", let_through)),
]


def outcome(test_id, plant=None):
    suite = unittest.defaultTestLoader.loadTestsFromName(test_id)
    result = unittest.TestResult()
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
        if plant is None:
            suite.run(result)
        else:
            with plant():
                suite.run(result)
    if result.testsRun != 1 or result.skipped:
        return "not run"
    return "pass" if result.wasSuccessful() else "fail"


def run_all():
    rows = []
    for n, name, test_id, plant in FAULTS:
        clean, planted = outcome(test_id), outcome(test_id, plant)
        rows.append({"fault": n, "planted": name, "test": test_id.rsplit(".", 2)[-2] + "." + test_id.rsplit(".", 1)[-1],
                     "as_written": clean, "with_fault": planted,
                     "caught": "yes" if clean == "pass" and planted == "fail" else "no"})
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
