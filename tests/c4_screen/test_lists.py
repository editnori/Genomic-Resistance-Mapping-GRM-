"""scripts/c4-screen/lists.py, expect.py and fit_times.py: the lists, the results a run must hold and the fit
times (docs/screen-arms.md).

The list tests read the committed public lists (scripts/c3-inputs/lists.tsv);
every other dataset name, result and time is SYNTHETIC.
"""

import contextlib
import csv
import io
import json
import tempfile
import unittest
from pathlib import Path

from tests.c4_screen import common

lists = common.load("lists")
expect = common.load("expect")
fit_times = common.load("fit_times")
C3_LISTS = common.ROOT / "scripts" / "c3-inputs" / "lists.tsv"


def read(path):
    with open(path, encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


class Lists(unittest.TestCase):
    def test_the_committed_lists_are_what_the_script_writes(self):
        with tempfile.TemporaryDirectory() as tmp:
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                self.assertEqual(lists.main([str(C3_LISTS), tmp]), 0)
            self.assertEqual((Path(tmp) / "lists.tsv").read_text(), (common.KIT / "lists.tsv").read_text())
            species = read(Path(tmp) / "species.tsv")
        self.assertEqual(buf.getvalue().strip(),
                         "L 268 in 23 species; L less M 233 in 23 species, 55 in S; the floor arm moves 92 in 11 "
                         "species, 24 in S")
        self.assertEqual(len(species), 23)
        self.assertEqual(sum(int(r["L_not_M_floor_moves"]) for r in species), 92)

    def test_the_first_run_comes_first_and_runs_every_arm(self):
        rows = read(common.KIT / "lists.tsv")
        first = [r["first_run"] for r in rows]
        self.assertEqual(first, sorted(first, reverse=True))
        for r in rows:
            self.assertEqual(r["arms"], "k31,cmh-4000,cmh-16000,floor-3,twin" if r["first_run"] == "1" else "-")

    def test_the_floor_arm_moves_a_dataset_of_800_genomes_or_more(self):
        self.assertEqual([lists.default_floor(n) for n in (100, 799, 800, 1599, 1600)], [3, 3, 4, 7, 8])
        for r in read(common.KIT / "lists.tsv"):
            self.assertEqual(r["floor_moves"] == "1", int(r["samples"]) >= 800, r["dataset"])

    def test_enterobacter_cloacae_has_nine_datasets_eight_in_s_and_none_moves(self):
        rows = [r for r in read(common.KIT / "lists.tsv") if r["species_slug"] == "enterobacter_cloacae"]
        self.assertEqual((len(rows), sum(r["in_S"] == "1" for r in rows), sum(r["floor_moves"] == "1" for r in rows)),
                         (9, 8, 0))
        self.assertEqual({r["species"] for r in rows}, {"Enterobacter cloacae"})


class Expect(unittest.TestCase):
    def rows(self, **kw):
        args = dict(species=common.SPECIES, datasets=[common.DATASET, common.OTHER], study=[common.DATASET],
                    seeds=[42, 43, 44], twin=[45, 46, 47])
        args.update(kw)
        return expect.rows(**args)

    def test_each_arm_runs_each_seed_and_the_study_split_once(self):
        rows = self.rows()
        by = {}
        for r in rows:
            by.setdefault((r[0], r[2], r[3]), []).append(r[5])
        for arm in ("k31",) + common.ARMS:
            self.assertEqual(by[(common.DATASET, arm, "lineage")], [42, 43, 44])
            self.assertEqual(by[(common.DATASET, arm, "study")], [42])
            self.assertNotIn((common.OTHER, arm, "study"), by)
        self.assertEqual(by[(common.DATASET, "twin", "lineage")], [45, 46, 47])
        self.assertNotIn((common.DATASET, "twin", "study"), by)
        self.assertEqual(len(rows), 2 * 4 * 3 + 4 + 2 * 3)

    def test_each_result_is_on_its_own_split(self):
        for r in self.rows():
            self.assertIn((r[3], r[4]), (("lineage", "lineage"), ("study", "study-whole")), r)
            self.assertTrue(r[7].endswith(f".{r[4]}-s{r[5]}.pick.json"), r)

    def test_only_an_arm_names_an_inputs_file(self):
        for r in self.rows():
            tag = "drawn" if r[3] == "lineage" else r[4]
            want = "-" if r[2] in ("k31", "twin") else f"{r[2]}/{r[0]}.{tag}-s{r[5]}.inputs.json"
            self.assertEqual(r[8], want)
            self.assertTrue(r[7].startswith(r[2] + "/"))

    def test_each_arms_flags_are_its_own(self):
        self.assertEqual(expect.ARM_FLAGS["cmh-4000"], ["--exp-screen", "cmh", "--force"])
        self.assertEqual(expect.ARM_FLAGS["cmh-16000"],
                         ["--exp-screen", "cmh", "--screened-features", "16000", "--force"])
        self.assertEqual(expect.ARM_FLAGS["floor-3"], ["--exp-floor", "3", "--force"])

    def test_a_twin_at_a_draw_seed_or_draws_without_seed_42_are_refused(self):
        with self.assertRaises(ValueError):
            self.rows(twin=[44, 45])
        with self.assertRaises(ValueError):
            self.rows(seeds=[43, 44])

    def test_an_unknown_arm_or_a_study_dataset_outside_the_list_is_refused(self):
        with self.assertRaises(ValueError):
            self.rows(arms=["cmh-8000"])
        with self.assertRaises(ValueError):
            self.rows(study=[common.SPECIES + "__drug_z"])

    def test_the_command_line_prints_the_rows(self):
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            self.assertEqual(expect.main([common.SPECIES, common.DATASET, "-", "42", "45", "floor-3"]), 0)
        lines = buf.getvalue().splitlines()
        self.assertEqual(lines[0].split("\t"), common.c3.COLUMNS)
        self.assertEqual(len(lines), 1 + 2 + 1)


class FitTimes(unittest.TestCase):
    def test_every_fit_time_of_every_inputs_file_is_a_row(self):
        d = common.pick(scheme="lineage")
        with tempfile.TemporaryDirectory() as tmp:
            for arm in common.ARMS:
                (Path(tmp) / arm).mkdir()
                with open(Path(tmp) / arm / f"{common.DATASET}.drawn-s42.inputs.json", "w") as f:
                    json.dump(common.inputs(d, arm=arm), f)
            out = Path(tmp) / "fits.tsv"
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                self.assertEqual(fit_times.main([tmp, str(out)]), 0)
            rows = read(out)
        # Two splits of three folds, five learners each, per arm.
        self.assertEqual(buf.getvalue().strip(), f"fit times: {3 * 2 * 3 * 5} rows")
        self.assertEqual(list(rows[0]), fit_times.COLUMNS)
        self.assertEqual({(r["arm"], r["rule"], r["budget"]) for r in rows},
                         {("cmh-4000", "cmh", "4000"), ("cmh-16000", "cmh", "16000"), ("floor-3", "chi-square", "4000")})


if __name__ == "__main__":
    unittest.main()
