"""Tests of scripts/c5-settings/lists.py and expect.py. The list test reads the committed lists
scripts/c3-inputs/lists.tsv and scripts/c5-settings/lists.tsv; the expect tests use made-up names."""

import contextlib
import csv
import io
import tempfile
import unittest
from pathlib import Path

from tests.c5_settings import common

LISTS = common.load("lists")
E = common.load("expect")


class Lists(unittest.TestCase):
    def test_the_list_is_s_with_s_less_m_first_and_is_the_committed_one(self):
        with tempfile.TemporaryDirectory() as tmp:
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                self.assertEqual(LISTS.main([str(common.ROOT / "scripts" / "c3-inputs" / "lists.tsv"), tmp]), 0)
            self.assertEqual(buf.getvalue(), "S 69 in 9 species; S less M 55 in 8 species\n")
            made = (Path(tmp) / "lists.tsv").read_bytes()
        self.assertEqual(made, (common.KIT / "lists.tsv").read_bytes())
        with open(common.KIT / "lists.tsv", newline="") as f:
            rows = list(csv.DictReader(f, delimiter="\t"))
        self.assertEqual(len(rows), 69)
        first = [r["first_run"] for r in rows]
        self.assertEqual(first, ["1"] * 55 + ["0"] * 14)
        self.assertTrue(all(r["in_M"] == ("0" if r["first_run"] == "1" else "1") for r in rows))
        self.assertEqual(len({r["species_slug"] for r in rows if r["first_run"] == "1"}), 8)


class Expect(unittest.TestCase):
    def rows(self, **kw):
        args = dict(species=common.SPECIES, datasets=[common.DATASET, common.OTHER], study=[common.DATASET],
                    seeds=[42, 43, 44], twin=[45, 46, 47])
        args.update(kw)
        return E.rows(**args)

    def test_each_draw_holds_the_default_and_both_arms_and_the_twin_its_own_seeds(self):
        rows = self.rows()
        # Per dataset: 3 draws and, for the first, 1 study run, each with 4 rows (the default, the search, the
        # cutoff arm and the plain best beside them); 3 twin rows.
        self.assertEqual(len(rows), (4 * 4 + 3) + (3 * 4 + 3))
        self.assertIn([common.DATASET, common.SPECIES, "plain-best", "lineage", "lineage", 43, "draw",
                       f"choices/{common.DATASET}.lineage-s43.f<k>.choices.json"], rows)
        self.assertIn([common.DATASET, common.SPECIES, "search", "study", "study-whole", 42, "study",
                       f"choices/{common.DATASET}.study-whole-s42.f<k>.choices.json"], rows)
        self.assertIn([common.OTHER, common.SPECIES, "default", "lineage", "lineage", 44, "draw",
                       f"default/{common.OTHER}.lineage-s44.pick.json"], rows)

    def test_no_arm_runs_at_a_twin_seed(self):
        rows = self.rows()
        self.assertEqual({r[2] for r in rows if r[6] == "twin"}, {"twin"})
        self.assertEqual({r[5] for r in rows if r[2] in ("search", "cutoff")}, {42, 43, 44})

    def test_a_twin_at_a_draw_seed_is_refused(self):
        with self.assertRaises(ValueError):
            self.rows(twin=[44, 45])
        with self.assertRaises(ValueError):
            self.rows(study=["other_species__drug"])


if __name__ == "__main__":
    unittest.main()
