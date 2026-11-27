"""The scan of scripts/matched-accuracy/stage1.py before any file is written, on SYNTHETIC bytes and releases.

Some tests need the shape of a real genome id or accession; each such token is put together here from pieces,
made of zeros or nines, and says so."""

import re
import unittest
from unittest import mock

from tests.matched_accuracy import common
from tests.matched_accuracy.common import outcome

s1 = common.stage1
SHAPES = {
    # Each made up from pieces: the shape of a real accession, never a real one.
    "a BioProject in lowercase": ("PRJ" + "ZZ").lower() + "9" * 6,
    "a BioProject after an underscore": "x_" + "PRJ" + "ZZ" + "9" * 6,
    "a BioSample after a letter": "x" + "SAM" + "ZZ" + "9" * 6,
    "an assembly in lowercase": ("GC" + "F_").lower() + "9" * 9,
    "a run": "S" + "RR" + "9" * 7,
    "a run of another registry": "E" + "RR" + "9" * 7,
    "a sample": "E" + "RS" + "9" * 7,
    "an experiment": "S" + "RX" + "9" * 7,
}


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

    def test_an_accession_without_case_inside_a_word_or_of_a_run_or_sample_is_a_hit(self):
        for name, token in SHAPES.items():
            with self.subTest(name=name):
                self.assertEqual(s1.hits(f"name\t{token}\n".encode(), set(), set()), 1)

    def test_a_biosample_accession_is_a_hit(self):
        # Made up from pieces: the shape of a real BioSample accession, never a real one.
        self.assertEqual(s1.hits(f"name\t{'SAM' + 'ZZ' + '9' * 6}\n".encode(), set(), set()), 1)

    def test_no_dataset_name_of_the_release_has_a_shape(self):
        # The 274 dataset names of the committed sums list reach counts.tsv and typing.tsv; a shape in one would
        # stop the real run.
        names = [line.split("./datasets/", 1)[1][:-len(".json")]
                 for line in common.REPORT_SUMS.read_text(encoding="utf-8").splitlines()]
        self.assertEqual(len(names), 274)
        self.assertEqual([n for n in names if s1.shaped(n)], [])
        self.assertEqual([n for n in names if re.search(r"[0-9]", n)], [])


class Scan(unittest.TestCase):
    FILES = {"a.tsv": b"rows\t3\n", "b.tsv": b"groups\t1\n", "c.tsv": b"pairs\t2\n"}
    # A shaped read id (digits, a dot, digits, made of nines here) beside two ids without a shape.
    IDS = {"9" * 3 + "." + "9" * 5, "gz-2", "gz-1"}

    def test_a_clean_export_gives_no_hit(self):
        self.assertEqual(outcome(s1.scan, self.FILES, self.IDS, {"pz-1"}), ("ok", []))

    def test_each_file_with_a_hit_is_named_with_its_count(self):
        files = dict(self.FILES, **{"b.tsv": b"gz-1 pz-1 gz-1\n"})
        self.assertEqual(outcome(s1.scan, files, self.IDS, {"pz-1"}), ("ok", [("b.tsv", 3)]))

    def test_the_positive_control_stops_a_scan_that_finds_nothing(self):
        with mock.patch.object(s1, "hits", return_value=0):
            got = outcome(s1.scan, self.FILES, self.IDS, {"pz-1"})
        self.assertEqual(got, ("stopped", "the scan's positive control found no read genome id without an id's "
                                          "shape in a copy of an export file that holds one"))

    def test_the_positive_control_plants_each_token_apart_in_a_copy_of_every_export_file(self):
        seen = []
        real = s1.hits

        def spy(data, ids, accessions):
            seen.append(data)
            return real(data, ids, accessions)
        with mock.patch.object(s1, "hits", side_effect=spy):
            self.assertEqual(s1.scan(self.FILES, self.IDS, {"pz-1"}), [])
        tokens = ["gz-1", "pz-1"] + [token for _, token in s1.CONTROL_TOKENS]
        want = list(self.FILES.values())
        for token in tokens:
            want += [data + b"\n" + token.encode() + b"\n" for _, data in sorted(self.FILES.items())]
        self.assertEqual(seen, want)
        self.assertEqual(self.FILES["a.tsv"], b"rows\t3\n")

    def test_the_control_plants_one_token_of_each_shape(self):
        tokens = [token for _, token in s1.CONTROL_TOKENS]
        self.assertEqual([t for t in tokens if s1.ID_SHAPED.fullmatch(t)], tokens[:1])
        for token in tokens[1:]:
            with self.subTest(token=token):
                self.assertIsNotNone(s1.ACCESSION_SHAPED.fullmatch(token))
                self.assertIsNone(s1.ID_SHAPED.search(token))
        for alternative in (r"PRJ[A-Z]{1,2}\d+", r"SAM[A-Z]{1,2}\d+", r"GC[AF]_\d+", r"[SED]R[RSXP]\d+"):
            with self.subTest(alternative=alternative):
                self.assertTrue(any(re.fullmatch(alternative, t, re.I) for t in tokens[1:]))
        self.assertTrue(any(t != t.upper() for t in tokens[1:]))

    def test_when_every_read_id_has_a_shape_a_made_up_id_stands_in(self):
        plants = s1.control_plants({"9" * 3 + "." + "9" * 5}, set())
        self.assertEqual(plants[:1], [("made-up genome id without an id's shape", s1.CONTROL_ID,
                                       {"9" * 3 + "." + "9" * 5, s1.CONTROL_ID})])
        self.assertFalse(s1.shaped(s1.CONTROL_ID))

    def test_a_scan_cut_to_the_first_file_fails_the_control(self):
        def first_file(files, ids, accessions):
            name = sorted(files)[0]
            n = s1.hits(files[name], ids, accessions)
            return {name: n} if n else {}
        with mock.patch.object(s1, "find", side_effect=first_file):
            got = outcome(s1.scan, self.FILES, self.IDS, {"pz-1"})
        self.assertEqual(got[0], "stopped")
        self.assertIn("positive control found no read genome id without an id's shape", got[1])

    def test_a_scan_cut_to_one_shape_fails_the_control(self):
        cuts = {"the id shape alone": ("ACCESSION_SHAPED", re.compile(r"(?!)")),
                "the BioProject shape alone": ("ACCESSION_SHAPED", re.compile(r"PRJ[A-Z]{1,2}\d+", re.I)),
                "no id shape": ("ID_SHAPED", re.compile(r"(?!)"))}
        for name, (attr, cut) in cuts.items():
            with self.subTest(cut=name):
                with mock.patch.object(s1, attr, cut):
                    got = outcome(s1.scan, self.FILES, self.IDS, {"pz-1"})
                self.assertEqual(got[0], "stopped")
                self.assertIn("-shaped token in a copy of an export file", got[1])

    def test_a_scan_that_looks_up_no_read_id_fails_the_control(self):
        def shapes_alone(data, ids, accessions):
            text = data.decode("utf-8")
            return len(s1.ID_SHAPED.findall(text)) + len(s1.ACCESSION_SHAPED.findall(text))
        with mock.patch.object(s1, "hits", side_effect=shapes_alone):
            got = outcome(s1.scan, self.FILES, self.IDS, {"pz-1"})
        self.assertEqual(got, ("stopped", "the scan's positive control found no read genome id without an id's "
                                          "shape in a copy of an export file that holds one"))


class RunWithAHit(common.Case, unittest.TestCase):
    def test_a_hit_prints_each_file_and_its_count_and_writes_nothing(self):
        # A made-up genome id that is also part of a species name: counts.tsv holds synthia_gamma in 6 dataset rows
        # and 2 species rows; typing.tsv in 3 x 2 x 2 dataset rows and 2 x 2 species rows.
        rel = common.release()
        rel["genomes"].append({"genome_id": "synthia_gamma", "bioproject_accession": "-", "mlst": "made-up",
                               "date_inserted": "2025-06-22T08:30:00Z"})
        run = self.run_stage1(rel)
        self.assertEqual((run.code, run.lines), (1, ["counts.tsv 8", "typing.tsv 16"]))
        self.assertFalse(run.wrote())

    def test_a_genome_id_read_from_the_5_nov_table_alone_is_a_hit(self):
        # The same made-up id, in a row of the 5 Nov table alone: no other input names it.
        rel = common.release()
        rel["hc20"].append({"genome_id": "synthia_gamma", "genome_name": "made-up", "mlst": "made-up",
                            "cgmlst_hc20": "hc-1"})
        run = self.run_stage1(rel)
        self.assertEqual((run.code, run.lines), (1, ["counts.tsv 8", "typing.tsv 16"]))
        self.assertFalse(run.wrote())


if __name__ == "__main__":
    unittest.main()
