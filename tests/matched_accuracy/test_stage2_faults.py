"""scripts/matched-accuracy/stage2_faults.py: each fault's text occurs once in stage2.py, plants code that compiles,
and names one test of this suite that runs alone. These tests plant nothing and run no fault; the runner does."""

import ast
import importlib.util
import unittest

from tests.matched_accuracy import stage2_common as c
from tests.matched_accuracy.test_faults import test_ids

# The rules of stage 2, each with at least one fault.
RULES = ("row mean", "tie", "weights", "interval", "chance", "the line", "ten pairs", "export", "typing", "counts",
         "groups", "row list", "rows", "learners", "null fits")
# The faults the brief and the lead's amendments of 27 Nov ask for, by what each plants.
ASKED = {
    "S01": "sum", "S02": "math.fsum", "S03": "tie counted 1", "S04": "tie counted 0",
    "S05": "pairs pooled across the groups", "S06": "groups pooled across the datasets",
    "S07": "datasets pooled across species", "S08": "seed other than 42", "S09": "percentile",
    "S12": "score read before the line", "S14": "10-pair rule taken off", "S15": "decimal point",
    "S17": "typing support check taken off", "S25": "taken from the rounded bounds",
}


def load_faults():
    spec = importlib.util.spec_from_file_location("matched_accuracy_stage2_faults", c.FAULTS)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


faults = load_faults()


class Faults(unittest.TestCase):
    def test_each_faults_text_occurs_once_in_its_file(self):
        for n, _, _, path, old, new, _ in faults.FAULTS:
            text = (c.ROOT / path).read_text(encoding="utf-8")
            self.assertEqual(text.count(old), 1, f"fault {n}: {path}")
            self.assertNotEqual(old, new, f"fault {n}")

    def test_each_planted_file_still_compiles(self):
        for n, _, _, path, old, new, _ in faults.FAULTS:
            compile((c.ROOT / path).read_text(encoding="utf-8").replace(old, new), path, "exec")

    def test_each_faults_test_is_one_test_of_this_suite(self):
        ids = test_ids()
        for n, *_, test in faults.FAULTS:
            self.assertIn(test, ids, f"fault {n}")
            self.assertEqual([i for i in ids if test in i], [test], f"fault {n}")

    def test_the_faults_are_named_once_cover_every_rule_and_plant_in_stage_2_alone(self):
        names = [f[0] for f in faults.FAULTS]
        self.assertEqual(names, [f"S{i:02d}" for i in range(1, len(names) + 1)])
        self.assertEqual({f[1] for f in faults.FAULTS}, set(RULES))
        self.assertEqual({f[3] for f in faults.FAULTS}, {"scripts/matched-accuracy/stage2.py"})

    def test_each_asked_fault_is_in_the_table(self):
        by = {f[0]: f[2] for f in faults.FAULTS}
        for n, words in ASKED.items():
            self.assertIn(words, by[n], n)

    def test_the_sum_faults_change_the_row_mean_of_the_briefs_values(self):
        # The values 0, 0.12, 0.99, 0.99, 0.4: left to right 2.4999999999999996; sum and math.fsum 2.5.
        values = (0, 0.12, 0.99, 0.99, 0.4)
        self.assertNotEqual(sum(values), c.stage2.row_mean(values) * 5)
        self.assertEqual(sum(values), 2.5)

    def test_the_runner_is_stage_1s_given_this_table(self):
        tree = ast.parse(c.FAULTS.read_text(encoding="utf-8"))
        main = next(f for f in tree.body if isinstance(f, ast.FunctionDef) and f.name == "main")
        self.assertIn("runner.FAULTS = FAULTS", ast.unparse(main))
        self.assertIn('HERE / "planted_faults.py"', c.FAULTS.read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
