"""scripts/c2-cap/bioprojects.py: the added genomes' BioProjects, by id, from BV-BRC's genome table only (docs/cap-arms.md).

Every genome id, BioProject and answer here is SYNTHETIC, and the transport
is a stand-in: no test opens a network connection.
"""

import gzip
import json
import tempfile
import unittest
from pathlib import Path

from tests.c2_cap import common

bioprojects = common.load("bioprojects")

ENDPOINT = "https://www.bv-brc.org/api/genome/"


class Transport:
    """A stand-in for BV-BRC: records each request and answers from `known` ({genome: BioProject or None})."""

    def __init__(self, known, extra=None, fail_first=0):
        self.known, self.extra, self.fail_first = known, extra or [], fail_first
        self.requests = []

    def __call__(self, url, body):
        self.requests.append((url, body))
        if self.fail_first:
            self.fail_first -= 1
            raise OSError("synthetic connection reset")
        ids = body[len("in(genome_id,("):body.index("))")].split(",")
        rows = []
        for g in ids:
            if g in self.known:
                rows.append({"genome_id": g} if self.known[g] is None
                            else {"genome_id": g, "bioproject_accession": self.known[g]})
        return rows + self.extra


def selected(body):
    part = [p for p in body.split("&") if p.startswith("select(")]
    return part[0][len("select("):-1].split(",") if len(part) == 1 else None


class Fetch(unittest.TestCase):
    def setUp(self):
        self.ids = [f"6{i:05d}.1" for i in range(1203)]
        self.known = {g: (f"PRJNA{i % 7}" if i % 5 else None) for i, g in enumerate(self.ids) if i % 11}

    def test_the_fetch_asks_the_genome_table_for_genome_id_and_bioproject_accession_only(self):
        send = Transport(self.known)
        answers, batches = bioprojects.fetch(self.ids, send, sleep=lambda s: None)
        self.assertEqual(batches, 3)
        self.assertEqual(len(send.requests), 3)
        asked = []
        for url, body in send.requests:
            self.assertEqual(url, ENDPOINT)
            self.assertEqual(selected(body), ["genome_id", "bioproject_accession"])
            self.assertTrue(body.startswith("in(genome_id,("))
            asked += body[len("in(genome_id,("):body.index("))")].split(",")
        self.assertEqual(asked, self.ids)
        self.assertEqual(answers, self.known)

    def test_a_query_of_another_field_or_another_table_is_refused(self):
        url, body = bioprojects.query(self.ids[:3])
        bioprojects.check_query(url, body)
        with self.assertRaisesRegex(ValueError, "other fields"):
            bioprojects.check_query(url, body.replace("bioproject_accession", "bioproject_accession,collection_year"))
        with self.assertRaisesRegex(ValueError, "other fields"):
            bioprojects.check_query(url, body.replace("&select(genome_id,bioproject_accession)", ""))
        with self.assertRaisesRegex(ValueError, "another table"):
            bioprojects.check_query("https://www.bv-brc.org/api/genome_amr/", body)

    def test_an_answer_with_another_field_an_unasked_genome_or_a_repeat_stops_the_fetch(self):
        for extra, why in (([{"genome_id": self.ids[1], "collection_year": "2001"}], "not asked"),
                           ([{"genome_id": "5.1", "bioproject_accession": "PRJNA1"}], "not asked"),
                           ([{"genome_id": self.ids[1], "bioproject_accession": "PRJNA1"}], "twice")):
            with self.subTest(why=why):
                send = Transport(self.known, extra=extra)
                with self.assertRaises(ValueError):
                    bioprojects.fetch(self.ids, send, sleep=lambda s: None)

    def test_a_dropped_connection_is_tried_again_and_then_stops_the_fetch(self):
        waits = []
        send = Transport(self.known, fail_first=2)
        answers, _ = bioprojects.fetch(self.ids[:10], send, sleep=waits.append)
        self.assertEqual((len(send.requests), waits), (3, [1, 2]))
        self.assertEqual(answers, {g: v for g, v in self.known.items() if g in self.ids[:10]})
        with self.assertRaises(OSError):
            bioprojects.fetch(self.ids[:10], Transport(self.known, fail_first=9), sleep=lambda s: None)


class Main(unittest.TestCase):
    def test_answers_stay_with_the_genome_rows_and_the_counts_hold_no_id(self):
        base = [f"9{i:03d}.1" for i in range(20)]
        added = [f"9{i:03d}.2" for i in range(12)]
        known = {g: (f"PRJ{i % 3}" if i < 9 else None) for i, g in enumerate(added) if i != 11}
        master = {added[0]: "PRJ0", added[1]: "PRJ9", base[0]: "PRJ5"}
        with tempfile.TemporaryDirectory() as tmp:
            t = Path(tmp)
            common.write_ids(t / "cohorts" / "cap-2000.txt", base)
            common.write_ids(t / "cohorts" / "cap-all.txt", sorted(base + added))
            common.meta(t / "meta.jsonl.gz", master)
            send = Transport(known)
            rc, text = common.quiet(bioprojects.main, ["bioprojects.py", str(t / "cohorts"), str(t / "meta.jsonl.gz"),
                                                       str(t / "answers.jsonl.gz"), str(t / "counts.json")], send)
            doc = json.loads((t / "counts.json").read_text(encoding="utf-8"))
            with gzip.open(t / "answers.jsonl.gz", "rt", encoding="utf-8") as f:
                rows = [json.loads(line) for line in f]
            counts_text = (t / "counts.json").read_text(encoding="utf-8")
        self.assertEqual(rc, 0)
        asked = [g for _, body in send.requests for g in body[len("in(genome_id,("):body.index("))")].split(",")]
        self.assertEqual(asked, sorted(added))
        self.assertEqual({k: doc[k] for k in ("asked", "answered", "with_bioproject", "without_bioproject",
                                              "not_answered", "master_same", "master_other", "batches")},
                         {"asked": 12, "answered": 11, "with_bioproject": 9, "without_bioproject": 2,
                          "not_answered": 1, "master_same": 1, "master_other": 1, "batches": 1})
        self.assertEqual((doc["endpoint"], doc["table"], doc["fields"]),
                         (ENDPOINT, "genome", ["genome_id", "bioproject_accession"]))
        self.assertRegex(doc["fetched_utc"], r"^\d{4}-\d\d-\d\dT\d\d:\d\d:\d\dZ$")
        self.assertEqual(rows, [{"genome_id": g, "bioproject_accession": known[g]} for g in sorted(known) if known[g]])
        self.assertFalse([g for g in base + added if g in counts_text or g in text])


if __name__ == "__main__":
    unittest.main()
