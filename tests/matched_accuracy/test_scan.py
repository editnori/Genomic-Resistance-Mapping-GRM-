"""The scan of scripts/matched-accuracy/stage1.py before any file is written, on SYNTHETIC bytes and releases."""

import unittest
from unittest import mock

from tests.matched_accuracy import common
from tests.matched_accuracy.common import outcome

s1 = common.stage1


class Hits(unittest.TestCase):
    def test_a_genome_id_inside_a_file_is_a_hit(self):
        self.assertEqual(s1.hits(b"name\tgz-1\n", {"gz-1"}, set()), 1)
        self.assertEqual(s1.hits(b"gz-1 gz-1 gz-2", {"gz-1", "gz-2"}, set()), 3)
        self.assertEqual(s1.hits(b"counts\t3\n", {"gz-1"}, set()), 0)

    def test_an_accession_inside_a_file_is_a_hit(self):
        self.assertEqual(s1.hits(b"x pz-1 y", set(), {"pz-1"}), 1)
        self.assertEqual(s1.hits(b"x pz-1 y", {"gz-1"}, {"pz-9"}), 0)

    def test_a_token_with_the_shape_of_a_real_genome_id_is_a_hit(self):
        # This test needs the shape of a real genome id (digits, a dot, digits), and of a real accession; each is
        # put together here from pieces and is in no list the scan was given.
        shaped = "9" * 3 + "." + "9" * 5
        self.assertEqual(s1.hits(f"x {shaped} y".encode(), set(), set()), 1)
        self.assertEqual(s1.hits(b"a share of 0.25", set(), set()), 1)
        accession = "PRJ" + "NA" + "9" * 6
        self.assertEqual(s1.hits(f"x {accession}".encode(), set(), set()), 1)


class Scan(unittest.TestCase):
    FILES = {"a.tsv": b"rows\t3\n", "b.tsv": b"groups\t1\n"}

    def test_a_clean_export_gives_no_hit(self):
        self.assertEqual(outcome(s1.scan, self.FILES, {"gz-1"}, {"pz-1"}), ("ok", []))

    def test_each_file_with_a_hit_is_named_with_its_count(self):
        files = dict(self.FILES, **{"b.tsv": b"gz-1 pz-1 gz-1\n"})
        self.assertEqual(outcome(s1.scan, files, {"gz-1"}, {"pz-1"}), ("ok", [("b.tsv", 3)]))

    def test_the_positive_control_stops_a_scan_that_finds_nothing(self):
        with mock.patch.object(s1, "hits", return_value=0):
            got = outcome(s1.scan, self.FILES, {"gz-1"}, {"pz-1"})
        self.assertEqual(got, ("stopped", "the scan's positive control found no genome id in a file that holds one"))

    def test_the_positive_control_plants_a_real_id_in_a_copy_of_an_export_file(self):
        seen = []

        def spy(data, ids, accessions):
            seen.append(data)
            return 1
        with mock.patch.object(s1, "hits", side_effect=spy):
            s1.scan(self.FILES, {"gz-2", "gz-1"}, {"pz-1"})
        self.assertEqual(seen[:2], [b"rows\t3\ngz-1", b"rows\t3\npz-1"])
        self.assertEqual(self.FILES["a.tsv"], b"rows\t3\n")


class RunWithAHit(common.Case, unittest.TestCase):
    def test_a_hit_prints_each_file_and_its_count_and_writes_nothing(self):
        # A made-up genome id that is also part of a species name: counts.tsv holds synthia_gamma in 6 dataset rows
        # and 2 species rows; typing.tsv in 3 x 2 x 2 dataset rows and 2 x 2 species rows.
        rel = common.release()
        rel["genomes"].append({"genome_id": "synthia_gamma", "bioproject_accession": "pz-q-1", "mlst": "made-up",
                               "date_inserted": "2025-06-22T08:30:00Z"})
        run = self.run_stage1(rel)
        self.assertEqual((run.code, run.lines), (1, ["counts.tsv 8", "typing.tsv 16"]))
        self.assertFalse(run.wrote())


if __name__ == "__main__":
    unittest.main()
