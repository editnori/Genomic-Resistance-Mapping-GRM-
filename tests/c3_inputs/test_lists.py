"""scripts/c3-inputs/lists.py and expect.py: the sets L and A and the results a run must hold (docs/input-arms.md).

The list tests read the committed public lists (scripts/c9-pick/l.tsv and
reference/amrfinderplus-comparison/datasets.tsv); every other dataset name
is SYNTHETIC.
"""

import contextlib
import csv
import io
import tempfile
import unittest
from pathlib import Path

from tests.c3_inputs import common

lists = common.load("lists")
expect = common.load("expect")
L_TSV = common.ROOT / "scripts" / "c9-pick" / "l.tsv"
A_TSV = common.ROOT / "reference" / "amrfinderplus-comparison" / "datasets.tsv"


def read(path):
    with open(path, encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


class Lists(unittest.TestCase):
    def test_a_is_the_rows_with_an_organism_option_and_a_lineage_split(self):
        rows = [{"dataset": "s__a", "organism_option": "yes", "lineage_split": "yes"},
                {"dataset": "s__b", "organism_option": "no", "lineage_split": "yes"},
                {"dataset": "s__c", "organism_option": "yes", "lineage_split": "no"}]
        self.assertEqual(lists.a_set(rows), {"s__a"})

    def test_an_a_dataset_outside_l_stops_the_lists(self):
        with self.assertRaises(SystemExit):
            lists.lists([{"dataset": "s__a", "first_run": "1"}], {"s__a", "s__z"})

    def test_the_committed_lists_are_what_the_script_writes(self):
        with tempfile.TemporaryDirectory() as tmp:
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                self.assertEqual(lists.main([str(L_TSV), str(A_TSV), tmp]), 0)
            self.assertEqual((Path(tmp) / "lists.tsv").read_text(), (common.KIT / "lists.tsv").read_text())
        self.assertEqual(buf.getvalue().strip(),
                         "L 268 in 23 species; L less M 233, 55 in S; A 237 in 21 species; A less M 204, 54 in S; "
                         "with an S and R pair: 32 of L less M, 21 of A less M")

    def test_the_program_arms_run_on_a_alone(self):
        rows = read(common.KIT / "lists.tsv")
        for r in rows:
            self.assertEqual("program" in r["arms"].split(","), r["in_A"] == "1", r["dataset"])
        self.assertEqual(sum(r["in_A"] == "1" for r in rows), 237)
        self.assertEqual({r["species"] for r in rows if r["in_A"] == "0" and r["dataset"].split("__")[0] in
                          ("mycobacterium_tuberculosis", "staphylococcus_haemolyticus")},
                         {"Mycobacterium tuberculosis", "Staphylococcus haemolyticus"})


class Expect(unittest.TestCase):
    def rows(self, **kw):
        args = dict(species=common.SPECIES, datasets=[common.DATASET, common.OTHER], a=[common.DATASET],
                    study=[common.DATASET], seeds=[42, 43, 44], twin=[45, 46, 47])
        args.update(kw)
        return expect.rows(**args)

    def test_each_arm_runs_each_seed_and_the_study_split_once(self):
        rows = self.rows()
        by = {}
        for r in rows:
            by.setdefault((r[0], r[2], r[3]), []).append(r[5])
        self.assertEqual(by[(common.DATASET, "k31", "lineage")], [42, 43, 44])
        self.assertEqual(by[(common.DATASET, "program", "study")], [42])
        self.assertEqual(by[(common.DATASET, "twin", "lineage")], [45, 46, 47])
        self.assertNotIn((common.OTHER, "program", "lineage"), by)
        self.assertNotIn((common.DATASET, "twin", "study"), by)
        self.assertEqual(len(rows), 2 * (3 + 3) + 2 * 3 + 3 * 1)

    def test_only_an_arm_names_an_inputs_file(self):
        for r in self.rows():
            want = "-" if r[2] in ("k31", "twin") else f"{r[2]}/{r[0]}.{'drawn' if r[3] == 'lineage' else r[4]}.inputs.json"
            self.assertEqual(r[8], want)
            self.assertTrue(r[7].startswith(r[2] + "/"))

    def test_a_twin_at_a_draw_seed_is_refused(self):
        with self.assertRaises(ValueError):
            self.rows(twin=[44, 45])

    def test_draws_without_the_release_seed_are_refused(self):
        with self.assertRaises(ValueError):
            self.rows(seeds=[43, 44])

    def test_an_unknown_arm_or_an_a_dataset_outside_the_list_is_refused(self):
        with self.assertRaises(ValueError):
            self.rows(arms=["k21"])
        with self.assertRaises(ValueError):
            self.rows(a=[common.SPECIES + "__drug_z"])


if __name__ == "__main__":
    unittest.main()
