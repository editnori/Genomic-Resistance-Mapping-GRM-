"""docs/cohort-result.md against the committed tables of reference/cohort-result/.

This is the one test file of the suite that reads committed result tables:
the aggregates the reading program wrote on its one run, by species,
antibiotic, learner and line. It fits nothing and reads no genome.
"""

import csv
import hashlib
import unittest

from tests.c1_cutoff import common

page_tables = common.load("page_tables")

RECORD = common.ROOT / "reference" / "cohort-result"
PAGE = common.ROOT / "docs" / "cohort-result.md"


def rows(name):
    with open(RECORD / name, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class Page(unittest.TestCase):
    def setUp(self):
        self.lines = PAGE.read_text().split("\n")

    def test_every_printed_table_row_stands_in_the_page(self):
        tables = page_tables.tables(RECORD)
        self.assertEqual(
            sorted(tables),
            ["deciding", "fifty", "log_loss", "others", "part_a", "recall", "specificity", "twin", "widths"],
        )
        for name, printed in tables.items():
            self.assertGreater(len(printed), 2, name)
            self.assertIn(printed[1], self.lines, name)
            for line in printed[:1] + printed[2:]:
                self.assertEqual(self.lines.count(line), 1, f"{name}: {line}")

    def test_the_reading_is_the_one_the_summary_gives(self):
        summary = (RECORD / "summary.txt").read_text().split("\n")
        meeting = {}
        for r in rows("page-deciding.tsv"):
            if r["row"] != "learner":
                meeting[r["line"]] = (r["learners_meeting"], r["meeting"], r["candidate"])
        self.assertEqual(sorted(meeting), ["10", "15", "20"])
        for line, (count, names, candidate) in meeting.items():
            self.assertEqual(candidate, "0")
            self.assertIn(f"line {line}: {count} of 5 learners meet ({names or 'none'}): not a candidate", summary)
        self.assertEqual([meeting[line][0] for line in ("10", "15", "20")], ["0", "0", "1"])
        self.assertEqual(meeting["20"][1], "rf")
        text = "\n".join(self.lines)
        self.assertIn("Line 10 has 0, line 15 has 0 and line 20 has 1 (`rf`).", text)
        twin = [r for r in rows("page-twin.tsv") if r["row"] != "learner"]
        self.assertEqual(len(twin), 1)
        self.assertEqual((twin[0]["learners_meeting"], twin[0]["meeting"]), ("3", "cart,rf,gbt"))
        self.assertIn("meet the row on 3 of 5 learners (`cart`, `rf` and `gbt`)", text)
        self.assertIn("the row holds against a change of fold seed alone", "\n".join(summary))

    def test_the_page_names_every_species_and_antibiotic_of_the_run(self):
        species = rows("page-species.tsv")
        datasets = rows("page-datasets.tsv")
        self.assertEqual((len(species), len(datasets)), (10, 89))
        for s in species:
            mine = [d for d in datasets if d["species"] == s["species"]]
            self.assertEqual(len(mine), int(s["datasets"]))
            self.assertEqual(", ".join(d["antibiotic"] for d in mine), s["antibiotics"])
            plans = sum(1 for d in mine if d["study_plan"] == "yes")
            self.assertEqual(plans, int(s["study_datasets"]))
            line = f"| *{s['species']}* | {s['antibiotics']} | {s['datasets']} | {s['study_datasets']} |"
            self.assertEqual(self.lines.count(line), 1, line)
        self.assertEqual(self.lines.count("| All | | 89 | 27 |"), 1)
        self.assertEqual(sum(int(s["results"]) for s in species), 2670)

    def test_the_record_holds_the_files_its_checksums_give(self):
        listed = {}
        for line in (RECORD / "SHA256SUMS").read_text().splitlines():
            digest, name = line.split("  ", 1)
            listed[name] = digest
        found = sorted(
            str(p.relative_to(RECORD)) for p in RECORD.rglob("*") if p.is_file() and p.name != "SHA256SUMS"
        )
        self.assertEqual(sorted(listed), found)
        for name, digest in listed.items():
            self.assertEqual(hashlib.sha256((RECORD / name).read_bytes()).hexdigest(), digest, name)


if __name__ == "__main__":
    unittest.main()
