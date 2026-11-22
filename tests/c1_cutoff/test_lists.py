"""scripts/c1-cutoff/lists.py and expect.py: the dataset lists and the expected results (docs/cohort-draw.md).

Every report, label count row and dataset name is SYNTHETIC.
"""

import contextlib
import io
import tempfile
import unittest
from pathlib import Path

from tests.c1_cutoff import common

lists = common.load("lists")
expect = common.load("expect")


def rep(d, r, s, ran=True):
    return {"dataset": d, "species_slug": d.split("__")[0], "species": "x", "antibiotic": "y", "samples": r + s,
            "resistant": r, "susceptible": s, "smaller_class": min(r, s), "lineages": 30, "lineage_ran": ran}


class Lists(unittest.TestCase):
    def test_big_is_a_lineage_figure_and_200_in_the_smaller_class_split_by_m_and_s(self):
        reports = [rep("a_x__d1", 200, 900), rep("a_x__d2", 199, 900), rep("a_x__d3", 300, 300, ran=False),
                   rep("b_y__d1", 500, 250), rep("b_y__d2", 250, 250)]
        labels = {"a_x__d1": (True, 0), "a_x__d2": (False, 0), "a_x__d3": (False, 0), "b_y__d1": (False, 4), "b_y__d2": (False, 0)}
        big, counts = lists.lists(reports, labels, {"a_x__d1", "b_y__d2"})
        self.assertEqual([r["dataset"] for r in big], ["a_x__d1", "b_y__d1", "b_y__d2"])
        self.assertEqual([(r["in_M"], r["in_S"], r["first_run"]) for r in big], [(1, 1, 0), (0, 0, 1), (0, 1, 1)])
        c = {name: (value, of) for name, value, of in counts}
        self.assertEqual(c["BIG datasets"], (3, 5))
        self.assertEqual(c["BIG datasets not in M (the first run)"], (2, 3))
        self.assertEqual(c["first-run datasets in S"], (1, 2))
        self.assertEqual(c["first-run datasets with an S and R pair"], (1, 2))
        self.assertEqual(c["BIG species"], (2, 2))

    def test_the_study_list_and_the_label_table_must_name_the_releases_datasets(self):
        reports = [rep("a_x__d1", 200, 900)]
        with self.assertRaisesRegex(ValueError, "label count table"):
            lists.lists(reports, {}, set())
        with self.assertRaisesRegex(ValueError, "study plan"):
            lists.lists(reports, {"a_x__d1": (False, 0)}, {"a_x__d9"})

    def test_a_study_plan_line_names_its_dataset(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / "plans.sha256"
            p.write_text("# header\n" + "0" * 64 + "  PLAN.md\n" + "1" * 64 + "  work/folds/a_x__d1.study-whole.json\n"
                         + "2" * 64 + "  work/folds/a_x__d1.study-d1.json\n" + "3" * 64 + "  folds/b_y__d2.study-whole.json\n")
            self.assertEqual(lists.read_study(p), {"a_x__d1", "b_y__d2"})


class Expect(unittest.TestCase):
    def draws(self):
        out = []
        for r in (10, 25):
            for j in range(2):
                out.append({"dataset": common.DATASET, "r": str(r), "draw": str(j), "written": "1", "study": "0",
                            "file": f"{common.DATASET}.r{r}.d{j}.cohort.json", "sha256": f"{r}{j}"})
        out.append({"dataset": "synthia_coli__drug_b", "r": "10", "draw": "0", "written": "0", "study": "0", "file": "", "sha256": ""})
        return out

    def test_every_cohort_at_seed_42_and_the_r_25_cohorts_again_at_45(self):
        rows = expect.rows("synthia_coli", self.draws(), True)
        self.assertEqual([(r["r"], r["draw"], r["seed"]) for r in rows],
                         [(10, 0, 42), (10, 1, 42), (25, 0, 42), (25, 0, 45), (25, 1, 42), (25, 1, 45)])
        self.assertEqual(rows[3]["file"], f"{common.DATASET}.r25.d0.s45.cohort.json")
        self.assertEqual(len(expect.rows("synthia_coli", self.draws(), False)), 4)
        with self.assertRaisesRegex(ValueError, "not a dataset of"):
            expect.rows("other_coli", self.draws(), True)

    def test_main_refuses_a_twin_flag_other_than_0_or_1(self):
        with contextlib.redirect_stderr(io.StringIO()):
            self.assertEqual(expect.main(["synthia_coli", "draws.tsv", "2"]), 2)


if __name__ == "__main__":
    unittest.main()
