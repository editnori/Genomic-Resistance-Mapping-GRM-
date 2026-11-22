"""scripts/c6-mic/expect.py: the results one species' MIC run must hold (docs/mic-plan.md).

The dataset names are SYNTHETIC.
"""

import contextlib
import io
import unittest

from tests.c6_mic import common

expect = common.load("expect")


class Expect(unittest.TestCase):
    def test_six_plans_per_dataset_and_split_with_the_twin_on_the_released_model(self):
        rows = expect.rows("synthia_coli", ["synthia_coli__drug_a", "synthia_coli__drug_b"], ["synthia_coli__drug_a"])
        self.assertEqual(len(rows), 18)
        lineage = [r for r in rows if r["dataset"] == "synthia_coli__drug_a" and r["split"] == "lineage"]
        self.assertEqual([r["plan"] for r in lineage], [f"lineage-s{s}" for s in range(42, 48)])
        study = [r for r in rows if r["split"] == "study"]
        self.assertEqual([r["plan"] for r in study],
                         ["study-whole", "study-d1", "study-d2", "study-d3", "study-d4", "c6mic-s47"])
        for r in rows:
            want = "released,ordered" if r["seed"] in (42, 43, 44) else "released"
            self.assertEqual(r["arms"], want, r)
            self.assertEqual(r["file"], f"{r['split']}/{r['dataset']}.{r['plan']}.mic.json")

    def test_a_dataset_of_another_species_or_listed_twice_is_refused(self):
        with self.assertRaisesRegex(ValueError, "not a dataset of"):
            expect.datasets("synthia_coli__drug_a,other_coli__drug_a", "synthia_coli")
        with self.assertRaisesRegex(ValueError, "twice"):
            expect.datasets("synthia_coli__drug_a,synthia_coli__drug_a", "synthia_coli")
        self.assertEqual(expect.datasets("-", "synthia_coli"), [])

    def test_main_writes_the_columns_score_reads(self):
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            code = expect.main(["synthia_coli", "synthia_coli__drug_a", "-"])
        self.assertEqual(code, 0)
        lines = out.getvalue().splitlines()
        self.assertEqual(lines[0].split("\t"), common.load("score").EXPECT_COLUMNS)
        self.assertEqual(len(lines), 7)
        with contextlib.redirect_stderr(io.StringIO()):
            self.assertEqual(expect.main(["synthia_coli", "-", "-"]), 2)


if __name__ == "__main__":
    unittest.main()
