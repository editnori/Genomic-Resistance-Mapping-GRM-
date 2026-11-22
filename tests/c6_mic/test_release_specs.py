"""scripts/c6-mic/release_specs.py on SYNTHETIC dataset settings and reports (docs/mic-plan.md)."""

import contextlib
import io
import json
import os
import tempfile
import unittest

from tests.c6_mic import common

specs = common.load("release_specs")


def spec(drug, n, min_class=25):
    return {"id": f"synthia_coli__{drug}", "antibiotic": drug, "min_class": min_class,
            "breakpoint": None if n is None else {"r_above": 1, "agreement": 0.75, "n": n}}


def report(drug, n, min_class=25):
    return {"antibiotic": drug, "counts": {"min_class": min_class},
            "breakpoint": None if n is None else {"r_above": 1, "agreement": 0.5, "n": n}}


def run(root, reports, listed):
    """release_specs.main on SYNTHETIC files: (exit code, out dir files or None)."""
    with tempfile.TemporaryDirectory() as tmp:
        os.makedirs(os.path.join(tmp, "release", "datasets"))
        for drug, doc in reports.items():
            with open(os.path.join(tmp, "release", "datasets", f"synthia_coli__{drug}.json"), "w") as f:
                json.dump(doc, f)
        with open(os.path.join(tmp, "datasets.json"), "w") as f:
            json.dump(root, f)
        out = os.path.join(tmp, "out")
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            code = specs.main([os.path.join(tmp, "datasets.json"), os.path.join(tmp, "release"),
                               ",".join(f"synthia_coli__{d}" for d in listed), out])
        if not os.path.isdir(out):
            return code, None
        with open(os.path.join(out, "datasets.json")) as f:
            patched = json.load(f)
        with open(os.path.join(out, "differ.txt")) as f:
            differ = f.read()
        head, rows = common.read_tsv(os.path.join(out, "release-specs.tsv"))
        return code, (patched, differ, head, rows)


class ReleaseSpecs(unittest.TestCase):
    def test_only_a_listed_dataset_with_another_breakpoint_takes_the_release_one(self):
        # drug_a: no breakpoint on either side. drug_b: the root's breakpoint
        # read one more pair (n 11 against 10). drug_c differs too but is not
        # listed, so it stays as the root has it.
        root = [spec("drug_a", None), spec("drug_b", 11), spec("drug_c", 7)]
        reports = {"drug_a": report("drug_a", None), "drug_b": report("drug_b", 10),
                   "drug_c": report("drug_c", 6)}
        code, (patched, differ, head, rows) = run(root, reports, ["drug_a", "drug_b"])
        self.assertEqual(code, 0)
        self.assertEqual(differ, "synthia_coli__drug_b")
        self.assertEqual(patched[0], root[0])
        self.assertEqual(patched[1], dict(root[1], breakpoint=reports["drug_b"]["breakpoint"]))
        self.assertEqual(patched[2], root[2])
        self.assertEqual(head, specs.COLUMNS)
        self.assertEqual([(r["dataset"], r["breakpoint"], r["root_n"], r["release_n"]) for r in rows], [
            ("synthia_coli__drug_a", "equal", "-", "-"), ("synthia_coli__drug_b", "differs", "11", "10")])

    def test_another_class_minimum_or_a_missing_dataset_writes_nothing(self):
        root = [spec("drug_a", 10, min_class=25), spec("drug_b", 10)]
        reports = {"drug_a": report("drug_a", 10, min_class=20), "drug_b": report("drug_b", 10)}
        self.assertEqual(run(root, reports, ["drug_a"]), (1, None))
        self.assertEqual(run(root, reports, ["drug_b", "drug_z"]), (1, None))
        self.assertEqual(run(root, {"drug_a": reports["drug_a"]}, ["drug_b"]), (1, None))


if __name__ == "__main__":
    unittest.main()
