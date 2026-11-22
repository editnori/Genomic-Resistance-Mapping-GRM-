"""scripts/c1-cutoff/draw.py: the keyed round-robin draw (docs/cohort-draw.md).

Every sample is SYNTHETIC (tests/c1_cutoff/common.py). Three tests plant a
fault in a copy of the draw and check that the test that guards it fails.
"""

import contextlib
import hashlib
import io
import json
import tempfile
import unittest
from pathlib import Path

from tests.c1_cutoff import common

draw = common.load("draw")
D = common.DATASET


def ids(smp, rows):
    return sorted(smp[i]["id"] for i in rows)


class Draw(unittest.TestCase):
    def setUp(self):
        self.smp = common.samples()

    # ---- the properties a test guards, written once so a planted fault can be run through them ----

    def seeded(self, cohort):
        """Same inputs, same cohort; another draw number or dataset, another cohort."""
        a, _ = cohort(self.smp, D, 0, 10)
        self.assertEqual(a, cohort(self.smp, D, 0, 10)[0])
        self.assertNotEqual(ids(self.smp, a), ids(self.smp, cohort(self.smp, D, 1, 10)[0]))
        self.assertNotEqual(ids(self.smp, a), ids(self.smp, cohort(self.smp, "synthia_coli__drug_b", 0, 10)[0]))

    def refuses_short(self, cohort):
        small = common.samples(resistant=10, susceptible=40)
        rows, why = cohort(small, D, 0, 10)
        self.assertIsNone(rows)
        self.assertEqual(why, "class short of r: 10 resistant for r 10")
        rows, why = cohort(small, D, 0, 9)
        self.assertEqual(len(rows), 18)

    def keeps_lineages(self, encode_doc):
        doc = encode_doc(self.smp)
        by_id = {s["id"]: s["lineage"] for s in self.smp}
        self.assertEqual(doc["lineage"], [by_id[g] for g in doc["genomes"]])

    # ---- the tests ----

    def test_a_draw_is_its_keys_and_moves_with_the_draw_number_and_the_dataset(self):
        self.seeded(draw.cohort)
        # The keys are the design's, compared as hex text.
        self.assertEqual(draw.key(D, 3, "genome", "9000.7"),
                         hashlib.sha256(f"c1b-20261122|{D}|draw3|genome|9000.7".encode()).hexdigest())

    def test_each_class_walks_its_lineages_round_robin_in_key_order(self):
        for resistant in (True, False):
            order = draw.class_order(self.smp, D, 2, resistant)
            members = [i for i, s in enumerate(self.smp) if s["resistant"] == resistant]
            self.assertEqual(sorted(order), members)
            lin = sorted({self.smp[i]["lineage"] for i in members}, key=lambda l: draw.key(D, 2, "lineage", l))
            # Pass 1: one genome of each lineage, in lineage key order.
            self.assertEqual([self.smp[i]["lineage"] for i in order[:len(lin)]], lin)
            # Inside a lineage, genomes go by genome key.
            for l in lin:
                got = [i for i in order if self.smp[i]["lineage"] == l]
                self.assertEqual(got, sorted(got, key=lambda i: draw.key(D, 2, "genome", self.smp[i]["id"])))

    def test_the_cohort_at_a_smaller_r_lies_inside_the_cohort_at_a_larger_r(self):
        for j in range(5):
            prev = set()
            for r in (10, 15, 20, 25, 50):
                rows, why = draw.cohort(self.smp, D, j, r)
                self.assertIsNone(why)
                self.assertEqual(sum(self.smp[i]["resistant"] for i in rows), r)
                self.assertEqual(len(rows), 2 * r)
                self.assertTrue(prev <= set(rows))
                prev = set(rows)

    def test_r_genomes_of_a_class_span_as_many_lineages_as_it_has_up_to_r(self):
        for r in (10, 25, 50):
            rows, _ = draw.cohort(self.smp, D, 0, r)
            for resistant in (True, False):
                have = {s["lineage"] for s in self.smp if s["resistant"] == resistant}
                got = {self.smp[i]["lineage"] for i in rows if self.smp[i]["resistant"] == resistant}
                self.assertEqual(len(got), min(r, len(have)))

    def test_a_class_of_r_or_fewer_genomes_refuses_the_draw(self):
        self.refuses_short(draw.cohort)

    def test_each_cohort_genome_keeps_its_report_lineage_and_study_fold(self):
        study = {s["id"]: i % 5 for i, s in enumerate(self.smp)}
        for row, doc in draw.draw_dataset(self.smp, D, [10, 25], 2, study):
            self.assertEqual(row["written"], 1)
            self.keeps_lineages(lambda smp, doc=doc: doc)
            self.assertEqual(doc["study_fold"], [study[g] for g in doc["genomes"]])
            self.assertEqual(set(doc), {"format", "dataset", "r", "draw", "genomes", "resistant", "lineage", "study_fold"})

    def test_main_writes_cohort_files_and_a_list_with_no_genome_id(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            small = common.samples(resistant=12, susceptible=70, seed=5)
            common.write_release(tmp / "rel", [common.report(self.smp), common.report(small, dataset="synthia_coli__drug_b")])
            with contextlib.redirect_stdout(io.StringIO()) as said:
                self.assertEqual(draw.main([str(tmp / "rel"), f"{D},synthia_coli__drug_b", "10,25", "2", str(tmp / "out")]), 0)
            self.assertEqual(said.getvalue(), "draws: 6 cohort files of 8 draws in 2 datasets, 2 refused\n")
            table = (tmp / "out" / "draws.tsv").read_text()
            self.assertNotIn("9000.", table)
            lines = [l.split("\t") for l in table.splitlines()[1:]]
            self.assertEqual(len(lines), 8)
            refused = [l for l in lines if l[3] == "0"]
            self.assertEqual([(l[0], l[1]) for l in refused], [("synthia_coli__drug_b", "25")] * 2)
            for l in lines:
                if l[3] == "1":
                    data = (tmp / "out" / l[11]).read_bytes()
                    self.assertEqual(hashlib.sha256(data).hexdigest(), l[12])
                    self.assertEqual(json.loads(data)["format"], "grm-cohort-draw-v1")

    # ---- planted faults ----

    def test_planted_a_draw_that_ignores_its_seed_fails_the_seed_test(self):
        real = draw.key

        def fixed(dataset, j, kind, name):  # the draw number and the dataset left out of the key
            return real("x", 0, kind, name)

        try:
            draw.key = fixed
            with self.assertRaises(AssertionError):
                self.seeded(draw.cohort)
        finally:
            draw.key = real

    def test_planted_a_short_class_accepted_silently_fails_the_refusal_test(self):
        def silent(samples, dataset, j, r):
            out = []
            for resistant in (True, False):
                out.extend(draw.class_order(samples, dataset, j, resistant)[:r])
            return out, None

        with self.assertRaises(AssertionError):
            self.refuses_short(silent)

    def test_planted_a_draw_that_splits_a_lineage_group_fails_the_lineage_test(self):
        def split(smp):
            rows, _ = draw.cohort(smp, D, 0, 25)
            # Each genome gets a lineage of its own, so a lineage group splits.
            return {"genomes": [smp[i]["id"] for i in rows], "lineage": list(range(len(rows)))}

        with self.assertRaises(AssertionError):
            self.keeps_lineages(split)


if __name__ == "__main__":
    unittest.main()
