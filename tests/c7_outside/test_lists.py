"""scripts/c7-outside/lists.py and scan.py on SYNTHETIC lists and files, and the tracked lists against C3's."""

import contextlib
import csv
import io
import json
import tempfile
import unittest
from pathlib import Path

from tests.c7_outside import common as C

lists = C.load("lists")
scan = C.load("scan")
C3_LIST = C.ROOT / "scripts" / "c3-inputs" / "lists.tsv"


def read(path):
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def c3_row(d, first_run="1", in_m="0", in_s="0"):
    sp = d.split("__")[0]
    return {"dataset": d, "species_slug": sp, "species": "Synthia coli", "antibiotic": d.split("__")[1], "samples": "90",
            "resistant": "40", "susceptible": "50", "lineages": "12", "in_M": in_m, "in_S": in_s, "in_A": "0",
            "first_run": first_run, "first_run_A": "0", "pairs_S+R": "0", "labels_depend_on": "none of the four",
            "arms": "k31,twin"}


class Lists(unittest.TestCase):
    def test_the_first_run_in_s_gets_the_study_split(self):
        rows = lists.lists([c3_row(C.DATASET, in_s="1"), c3_row(C.OTHER), c3_row("synthia_coli__drug_c", "0", "1", "1")])
        self.assertEqual([r["splits"] for r in rows], ["lineage,study", "lineage", "lineage"])
        self.assertTrue(all(r["arms"] == "k31,twin,xgboost,lightgbm,elasticnet" for r in rows))
        sp = lists.species_rows(rows)
        self.assertEqual((sp[0]["L"], sp[0]["L_not_M"], sp[0]["L_not_M_in_S"]), (3, 2, 1))

    def test_a_dataset_twice_or_in_m_and_the_first_run_is_refused(self):
        with self.assertRaises(SystemExit):
            lists.lists([c3_row(C.DATASET), c3_row(C.DATASET)])
        with self.assertRaises(SystemExit):
            lists.lists([c3_row(C.DATASET, "1", "1")])

    def test_the_tracked_lists_are_the_ones_c3s_list_gives(self):
        with tempfile.TemporaryDirectory() as tmp:
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(lists.main([str(C3_LIST), tmp]), 0)
            for name in ("lists.tsv", "species.tsv"):
                self.assertEqual((Path(tmp) / name).read_bytes(), (C.KIT / name).read_bytes(), name)
        rows = read(C.KIT / "lists.tsv")
        self.assertEqual(len(rows), len(read(C3_LIST)))
        self.assertTrue(all(r["species"] and r["antibiotic"] for r in rows))


class Scan(unittest.TestCase):
    def test_an_id_of_the_species_is_found_and_counted_never_printed(self):
        with tempfile.TemporaryDirectory() as tmp:
            sp = Path(tmp) / "species.json"
            sp.write_text(json.dumps({"genomes": [{"id": "1234.56"}, {"id": "78.9"}]}))
            d = Path(tmp) / "out"
            d.mkdir()
            (d / "a.log").write_text("fold 1234.56 done; 1234.567 is another number; v1234.56 is a word\n")
            (d / "b.json").write_text('{"x": 0.5, "n": 78.9}')
            self.assertEqual(scan.scan(scan.ids_of(sp), [str(d)]), (2, 2))
            (d / "a.log").write_text("clean\n")
            (d / "b.json").write_text("{}")
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                self.assertEqual(scan.main([str(sp), str(d)]), 0)
            self.assertEqual(buf.getvalue(), "genome ids: 0 in 2 files\n")


if __name__ == "__main__":
    unittest.main()
