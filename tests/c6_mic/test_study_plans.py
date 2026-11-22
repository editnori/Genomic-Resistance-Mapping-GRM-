"""scripts/c6-mic/study_plans.py on a SYNTHETIC cohort (docs/mic-plan.md, the study plan for seed 47).

The cohort has 40 made-up genomes, 9000.1 to 9000.40, in eight made-up
studies of five (genome i in study i // 5), resistant when i is even. The
five SYNTHETIC reviewed plans keep each study whole and test study s in fold
(s + p) mod 5 in plan p, so studies s and s + 5 share a fold in every plan.
The units are then {0, 5}, {1, 6}, {2, 7}, {3} and {4}, as studies.
"""

import contextlib
import hashlib
import io
import json
import pathlib
import tempfile
import unittest

from tests.c6_mic import common

plans = common.load("study_plans")

N = 40
IDS = [f"9000.{i + 1}" for i in range(N)]
STUDY = [i // 5 for i in range(N)]


def report(dups=None, resistant=None):
    y = resistant or [i % 2 == 0 for i in range(N)]
    doc = {"id": "synthia_coli__drug_a", "samples": [{"id": g, "resistant": r} for g, r in zip(IDS, y)]}
    if dups is not None:
        doc["kset_duplicates"] = dups
    return doc


def reviewed():
    docs = []
    for p, name in enumerate(plans.REVIEWED):
        test = [[IDS[i] for i in range(N) if (STUDY[i] + p) % 5 == f] for f in range(5)]
        docs.append({"format": "grm-fold-file-v1", "dataset": "synthia_coli__drug_a", "scheme": name,
                     "test": test, "train": [], "unused": []})
    return docs


def folds_of(data):
    doc = json.loads(data)
    fold = {}
    for f, members in enumerate(doc["test"]):
        for g in members:
            assert g not in fold, g
            fold[g] = f
    return doc, [fold[g] for g in IDS]


class StudyPlan(unittest.TestCase):
    def test_the_plan_keeps_every_unit_and_study_whole_and_tests_each_genome_once(self):
        data, row = plans.make("synthia_coli__drug_a", report(), reviewed())
        self.assertIsNotNone(data, row)
        doc, fold = folds_of(data)
        self.assertEqual((doc["format"], doc["scheme"], doc["train"], doc["unused"]),
                         ("grm-fold-file-v1", "c6mic-s47", [], []))
        self.assertEqual(sorted(g for t in doc["test"] for g in t), sorted(IDS))
        units = [[i for i in range(N) if STUDY[i] in s] for s in ({0, 5}, {1, 6}, {2, 7}, {3}, {4})]
        for u in units:
            self.assertEqual(len({fold[i] for i in u}), 1, u)
        for s in range(8):
            self.assertEqual(len({fold[i] for i in range(N) if STUDY[i] == s}), 1, s)
        self.assertEqual((row["genomes"], row["units"], row["check"]), (40, 5, "ok"))
        self.assertEqual(row["sha256"], hashlib.sha256(data).hexdigest())

    def test_units_go_in_their_keyed_order_each_to_the_cheapest_fold(self):
        # Rule B: units sorted by sha256("c6mic-20261122|study-s47|" + lowest id).
        # Each unit holds both classes, so an empty fold is always the cheapest
        # and the k-th unit in that order goes to fold k.
        units = [[i for i in range(N) if STUDY[i] in s] for s in ({0, 5}, {1, 6}, {2, 7}, {3}, {4})]
        key = lambda u: hashlib.sha256(f"c6mic-20261122|study-s47|{IDS[u[0]]}".encode()).hexdigest()
        data, _ = plans.make("synthia_coli__drug_a", report(), reviewed())
        _, fold = folds_of(data)
        for k, u in enumerate(sorted(units, key=key)):
            self.assertEqual(fold[u[0]], k)
        again, _ = plans.make("synthia_coli__drug_a", report(), reviewed())
        self.assertEqual(again, data)

    def test_the_lowest_id_is_read_as_numbers(self):
        self.assertLess(plans.gnum("9000.9"), plans.gnum("9000.10"))
        self.assertEqual(min(["9000.10", "9000.9"], key=plans.gnum), "9000.9")

    def test_a_copy_group_across_two_units_gets_no_plan(self):
        # Genomes 0 (study 0) and 15 (study 3) given one k-mer set: the units put
        # them in two folds, so the plan fails its check and no file is written.
        data, row = plans.make("synthia_coli__drug_a", report(dups=[[0, 15]]), reviewed())
        self.assertIsNone(data)
        self.assertIn("a copy group is split", row["check"])
        data, row = plans.make("synthia_coli__drug_a", report(dups=[[0, 25]]), reviewed())
        self.assertIsNotNone(data, row)  # studies 0 and 5 share a unit

    def test_a_fold_without_both_classes_gets_no_plan(self):
        y = [STUDY[i] in (0, 5) for i in range(N)]  # only the first unit is resistant
        data, row = plans.make("synthia_coli__drug_a", report(resistant=y), reviewed())
        self.assertIsNone(data)
        self.assertIn("lack a class", row["check"])

    def test_check_names_a_split_unit(self):
        fold = [i % 5 for i in range(N)]
        y = [i % 2 == 0 for i in range(N)]
        self.assertIn("a unit is split", plans.check(y, fold, [], [[0, 1]]))
        self.assertEqual(plans.check(y, fold, [], [[0, 5]]), [])

    def test_a_reviewed_plan_must_test_the_cohort_once(self):
        cases = []
        doc = reviewed()[0]
        doc["test"][0].append("9000.99")
        cases.append((doc, "outside the cohort"))
        doc = reviewed()[0]
        doc["test"][1].append(doc["test"][0][0])
        cases.append((doc, "listed twice"))
        doc = reviewed()[0]
        doc["test"][0].pop()
        cases.append((doc, "in no test fold"))
        doc = reviewed()[0]
        doc["train"] = [doc["test"][0].pop()]
        cases.append((doc, "train or leave unused"))
        doc = reviewed()[0]
        doc["test"] = doc["test"][:4]
        cases.append((doc, "5 test folds"))
        for doc, words in cases:
            with self.assertRaisesRegex(ValueError, words):
                plans.read_plan(doc, IDS, "synthetic")

    def test_main_writes_the_plan_and_a_list_with_no_genome_id(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = pathlib.Path(tmp)
            (tmp / "release" / "datasets").mkdir(parents=True)
            (tmp / "release" / "datasets" / "synthia_coli__drug_a.json").write_text(json.dumps(report()))
            (tmp / "plans").mkdir()
            for doc in reviewed():
                (tmp / "plans" / f"synthia_coli__drug_a.{doc['scheme']}.json").write_text(json.dumps(doc))
            (tmp / "list.txt").write_text("synthia_coli__drug_a\n")
            out = tmp / "out"
            with contextlib.redirect_stdout(io.StringIO()):
                code = plans.main([str(tmp / "release"), str(tmp / "plans"), str(tmp / "list.txt"), str(out)])
            self.assertEqual(code, 0)
            self.assertTrue((out / "synthia_coli__drug_a.c6mic-s47.json").is_file())
            listing = (out / "study-plans.tsv").read_text()
            self.assertNotIn("9000.", listing)
            head, rows = common.read_tsv(out / "study-plans.tsv")
            self.assertEqual(head, ["dataset", "genomes", "units", "check", "sha256"])
            self.assertEqual(rows[0]["check"], "ok")
            y = [STUDY[i] in (0, 5) for i in range(N)]
            (tmp / "release" / "datasets" / "synthia_coli__drug_a.json").write_text(json.dumps(report(resistant=y)))
            with contextlib.redirect_stdout(io.StringIO()):
                code = plans.main([str(tmp / "release"), str(tmp / "plans"), str(tmp / "list.txt"), str(tmp / "out2")])
            self.assertEqual(code, 1)
            self.assertFalse((tmp / "out2" / "synthia_coli__drug_a.c6mic-s47.json").exists())


if __name__ == "__main__":
    unittest.main()
