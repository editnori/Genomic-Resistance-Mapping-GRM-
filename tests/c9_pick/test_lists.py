"""scripts/c9-pick/lists.py and expect.py: the dataset lists and the expected results (docs/served-pick.md).

Every report, label count row and dataset name is SYNTHETIC.
"""

import contextlib
import io
import json
import tempfile
import unittest
from pathlib import Path

from tests.c9_pick import common

lists = common.load("lists")
expect = common.load("expect")


def rep(d, r, s, ran=True, lineages=30):
    return {"dataset": d, "species_slug": d.split("__")[0], "species": "Synthia " + d.split("__")[0][-1],
            "antibiotic": d.split("__")[1], "samples": r + s, "resistant": r, "susceptible": s, "lineages": lineages,
            "lineage_ran": ran}


class Lists(unittest.TestCase):
    def test_l_is_a_lineage_figure_split_by_m_with_s_marked(self):
        reports = [rep("a_x__d1", 50, 90), rep("a_x__d2", 40, 90), rep("a_x__d3", 30, 30, ran=False),
                   rep("b_y__d1", 50, 25), rep("b_y__d2", 25, 25)]
        labels = {"a_x__d1": (True, 0), "a_x__d2": (False, 3), "a_x__d3": (False, 0), "b_y__d1": (False, 0),
                  "b_y__d2": (True, 2)}
        rows, counts = lists.lists(reports, labels, {"a_x__d1", "b_y__d1", "a_x__d3"})
        self.assertEqual([r["dataset"] for r in rows], ["a_x__d1", "a_x__d2", "b_y__d1", "b_y__d2"])
        self.assertEqual([(r["in_M"], r["in_S"], r["first_run"]) for r in rows], [(1, 1, 0), (0, 0, 1), (0, 1, 1), (1, 0, 0)])
        got = {name: (value, of) for name, value, of in counts}
        self.assertEqual(got["L: reports whose lineage split ran"], (4, 5))
        self.assertEqual(got["S datasets in L"], (2, 3))
        self.assertEqual(got["L datasets not in M (the first run)"], (2, 4))
        self.assertEqual(got["first-run datasets in S (the first study run)"], (1, 2))
        self.assertEqual(got["L datasets in M and S"], (1, 2))
        self.assertEqual(got["first-run datasets with an S and R pair"], (1, 2))
        self.assertEqual(got["a_x: L, not M, M, not M in S, M in S"], ("2, 1, 1, 0, 1", 2))
        self.assertEqual(rows[0]["species"], "Synthia x")
        self.assertEqual(rows[0]["antibiotic"], "d1")

    def test_a_label_table_naming_other_datasets_or_a_plan_outside_the_release_is_refused(self):
        reports = [rep("a_x__d1", 50, 90)]
        with self.assertRaises(ValueError):
            lists.lists(reports, {"a_x__d1": (False, 0), "a_x__d9": (False, 0)}, set())
        with self.assertRaises(ValueError):
            lists.lists(reports, {"a_x__d1": (False, 0)}, {"a_x__d9"})

    def test_it_reads_reports_labels_and_plans_from_files(self):
        with tempfile.TemporaryDirectory() as tmp:
            top = Path(tmp)
            (top / "release" / "datasets").mkdir(parents=True)
            for d, ran in (("a_x__d1", True), ("a_x__d2", False)):
                split = {"scheme": "lineage", "n_groups": 12}
                if not ran:
                    split["unavailable_reason"] = "too few lineages"
                (top / "release" / "datasets" / f"{d}.json").write_text(json.dumps(
                    {"id": d, "species_slug": "a_x", "species": "Synthia x", "antibiotic": d[-2:],
                     "counts": {"samples": 30, "resistant": 10, "susceptible": 20},
                     "splits": [{"scheme": "random", "n_groups": 30}, split]}))
            (top / "labels.tsv").write_text("dataset\tpairs_S+I\tpairs_I+R\tpairs_S+R\na_x__d1\t0\t0\t1\na_x__d2\t0\t1\t0\n")
            (top / "plans.sha256").write_text("0" * 64 + "  folds/a_x__d1.study-whole.json\n" + "1" * 64
                                              + "  folds/a_x__d2.study-other.json\n")
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                rc = lists.main([str(top / "release"), str(top / "labels.tsv"), str(top / "plans.sha256"), str(top / "out")])
            self.assertEqual(rc, 0)
            lines = (top / "out" / "l.tsv").read_text().splitlines()
            self.assertEqual(lines[0].split("\t"), lists.COLUMNS)
            self.assertEqual(lines[1].split("\t"), ["a_x__d1", "a_x", "Synthia x", "d1", "30", "10", "20", "12", "0", "1", "1", "1"])
            self.assertEqual(len(lines), 2)


class Expect(unittest.TestCase):
    def test_each_dataset_has_its_draws_and_twins_and_each_study_dataset_one_study_result(self):
        rows = expect.rows("a_x", ["a_x__d2", "a_x__d1"], ["a_x__d2"], [42, 43, 44], [45, 46, 47])
        self.assertEqual(len(rows), 2 * 6 + 1)
        self.assertEqual(rows[0], ["a_x__d1", "a_x", "lineage", "lineage", 42, 3, "draw", "a_x__d1.lineage-s42.pick.json"])
        self.assertEqual([r[4:7] for r in rows if r[0] == "a_x__d1"],
                         [[42, 3, "draw"], [43, 3, "draw"], [44, 3, "draw"], [45, 0, "twin"], [46, 0, "twin"], [47, 0, "twin"]])
        self.assertEqual(rows[-1], ["a_x__d2", "a_x", "study", "study-whole", 42, 3, "study",
                                    "a_x__d2.study-whole-s42.pick.json"])

    def test_the_names_match_what_grm_learn_writes(self):
        # src/learn/inner_pick.rs: file_name(dataset, scheme, seed).
        self.assertEqual(expect.rows("a_x", ["a_x__d1"], [], [42], [])[0][7], "a_x__d1.lineage-s42.pick.json")

    def test_a_list_it_cannot_run_is_refused(self):
        for args in (("a_x", ["b_y__d1"], [], [42], []), ("a_x", ["a_x__d1", "a_x__d1"], [], [42], []),
                     ("a_x", ["a_x__d1"], ["a_x__d2"], [42], []), ("a_x", ["a_x__d1"], [], [42], [42]),
                     ("a_x", [], [], [42], [])):
            with self.assertRaises(ValueError, msg=str(args)):
                expect.rows(*args)

    def test_the_command_takes_dash_for_no_twin_and_no_study(self):
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            self.assertEqual(expect.main(["a_x", "a_x__d1", "-", "42", "-"]), 0)
        lines = buf.getvalue().splitlines()
        self.assertEqual(lines[0].split("\t"), expect.COLUMNS)
        self.assertEqual(len(lines), 2)


if __name__ == "__main__":
    unittest.main()
