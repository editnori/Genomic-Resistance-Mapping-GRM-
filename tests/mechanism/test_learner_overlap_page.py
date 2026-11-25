"""docs/learner-overlap-result.md against the committed record of reference/learner-overlap/.

This test file reads committed result tables: the aggregates the table script
wrote from the one reading of the learner overlap, per report, per species and
pooled. It fits nothing and reads no genome.
"""

import hashlib
import re
import unittest

from tests.mechanism import common

lo = common.load("learner_overlap")
faults = common.load("overlap_faults")

RECORD = common.ROOT / "reference" / "learner-overlap"
PAGE = common.ROOT / "docs" / "learner-overlap-result.md"
INDEX = common.ROOT / "docs" / "INDEX.md"
SECTIONS = ["Reading", "Overlap", "Yardstick", "Species", "Limits", "Method", "Plan", "Datasets"]
TABLES = {"Overlap": "overlap", "Yardstick": "yardstick", "Species": "species", "Datasets": "datasets"}
FIVE = ("scm", "cart", "rf", "gbt", "lasso")
NUMBER = r"\d[\d,]*(?:\.\d+)?"


def sections(text):
    parts = re.split(r"^## (.+)$", text, flags=re.M)
    return parts[0], dict(zip(parts[1::2], parts[2::2]))


def prose(text):
    return "\n".join(line for line in text.split("\n") if not line.startswith("|"))


class Page(unittest.TestCase):
    def setUp(self):
        self.text = PAGE.read_text(encoding="utf-8")
        self.lines = self.text.split("\n")
        self.head, self.sections = sections(self.text)
        self.rows = lo.read_tsv(RECORD / "datasets.tsv")
        self.pooled = lo.read_tsv(RECORD / "pooled.tsv")
        self.gates = (RECORD / "gates.txt").read_text(encoding="utf-8").splitlines()

    def test_the_sections_stand_in_order(self):
        self.assertEqual(list(self.sections), SECTIONS)

    def test_every_printed_table_stands_whole_in_its_section(self):
        tables = lo.tables(str(RECORD))
        self.assertEqual(sorted(tables), sorted(TABLES.values()))
        for section, name in TABLES.items():
            printed = tables[name]
            self.assertGreater(len(printed), 2, name)
            self.assertIn("\n".join(printed) + "\n", self.sections[section], name)
            for line in printed[:1] + printed[2:]:
                self.assertEqual(self.lines.count(line), 1, f"{name}: {line}")

    def test_the_record_holds_the_files_its_checksums_give(self):
        listed = {}
        for line in (RECORD / "SHA256SUMS").read_text(encoding="utf-8").splitlines():
            digest, name = line.split("  ", 1)
            listed[name] = digest
        self.assertEqual(sorted(listed), sorted(lo.RECORD_FILES))
        for name, digest in listed.items():
            self.assertEqual(hashlib.sha256((RECORD / name).read_bytes()).hexdigest(), digest, name)
        found = sorted(p.name for p in RECORD.iterdir())
        self.assertEqual(found, sorted(list(listed) + ["SHA256SUMS", "README.md", "reading.txt"]))

    def test_the_counts_of_each_report_add_up(self):
        self.assertEqual(len(self.rows), len({(r["species"], r["antibiotic"]) for r in self.rows}))
        for r in self.rows:
            with self.subTest(report=(r["species"], r["antibiotic"])):
                n = {k: int(v) for k, v in r.items() if k not in ("species", "antibiotic", "learners", "catalog_rule",
                                                                   "expected_all5", "ratio_all5")}
                self.assertEqual((r["learners"], n["learner_count"]), (",".join(FIVE), 5))
                self.assertEqual(n["genomes"], n["resistant"] + n["susceptible"])
                self.assertEqual(sum(n[f"miss_{k}"] for k in range(6)), n["resistant"])
                self.assertEqual(sum(n[f"wrong_{k}"] for k in range(6)), n["susceptible"])
                self.assertEqual(n["miss_5"], n["miss_all5"])
                self.assertLessEqual(n["miss_all5"], min(n[f"miss_{x}"] for x in FIVE))
                self.assertLessEqual(n["miss_all5_no_catalog_call"], n["miss_all5"])
                self.assertLessEqual(n["miss_all5_in_null_fold"], n["miss_all5"])
                # each learner's misses are the genomes it misses among those missed by k of the five
                self.assertEqual(sum(n[f"miss_{x}"] for x in FIVE), sum(k * n[f"miss_{k}"] for k in range(6)))
                self.assertIn(r["catalog_rule"], ("yes", "no"))
                if r["catalog_rule"] == "no":
                    self.assertEqual(n["miss_all5_no_catalog_call"], n["miss_all5"])

    def test_the_yardstick_and_ratio_follow_from_the_miss_columns(self):
        for r in self.rows:
            exp = lo.expected_all5(int(r["resistant"]), [int(r[f"miss_{x}"]) for x in FIVE])
            self.assertEqual(r["expected_all5"], lo.d12(exp), r["antibiotic"])
            self.assertEqual(r["ratio_all5"], lo.d12(int(r["miss_all5"]) / exp) if exp > 0 else "", r["antibiotic"])

    def test_the_pooled_rows_are_what_the_script_pools(self):
        rows = [dict(r, learner_count=int(r["learner_count"]), resistant=int(r["resistant"])) for r in self.rows]
        again = [{k: str(v) for k, v in p.items()} for p in lo.pooled_rows(rows)]
        self.assertEqual(again, self.pooled)
        self.assertEqual(self.pooled[0]["datasets"], str(len(self.rows)))

    def test_the_gates_all_pass_and_stand_in_the_method(self):
        checks = [line for line in self.gates if line.startswith("check ")]
        self.assertEqual([line.split(" ", 2)[1] for line in checks], [str(n) for n in range(1, 11)])
        self.assertFalse(any("STOP" in line or "failed" in line for line in self.gates))
        for line in checks + [self.gates[-1]]:
            self.assertEqual(self.sections["Method"].count(line + "\n"), 1, line)
        self.assertEqual(self.gates[-1], "reading rule: stack")
        self.assertEqual(lo.reading_rule(float(self.pooled[0]["ratio_median"])), "stack")

    def test_the_reading_commands_stand_as_reading_txt_gives(self):
        reading = (RECORD / "reading.txt").read_text(encoding="utf-8").splitlines()
        commands = [line for line in reading if line.startswith("python3 -I ")]
        self.assertEqual(len(commands), 4)
        self.assertIn("```\n" + "\n".join(commands) + "\n```", self.sections["Method"])
        listed = {line.split("  ", 1)[1]: line.split("  ", 1)[0] for line in reading if "  " in line}
        for name in lo.RECORD_FILES + ("SHA256SUMS",):
            self.assertEqual(listed[f"record/{name}"], hashlib.sha256((RECORD / name).read_bytes()).hexdigest(), name)

    def test_the_fault_table_names_each_planted_faults_test(self):
        method = self.sections["Method"]
        for fault, words, test, _ in faults.FAULTS:
            self.assertEqual(self.lines.count(f"| {fault}. {words[0].upper()}{words[1:]} | "
                                              f"`{'.'.join(test.rsplit('.', 2)[-2:])}` |"), 1, fault)
        self.assertIn(f"All {len(faults.FAULTS)} were caught.", method)

    def test_the_page_gives_the_suites_own_test_counts(self):
        loader = unittest.TestLoader()
        here = loader.loadTestsFromName("tests.mechanism.test_learner_overlap_page").countTestCases()
        script = loader.loadTestsFromName("tests.mechanism.test_learner_overlap").countTestCases()
        second = loader.loadTestsFromName("tests.mechanism.test_overlap_second").countTestCases()
        method = self.sections["Method"]
        self.assertIn(f"`tests/mechanism/test_learner_overlap.py` holds {script} SYNTHETIC tests", method)
        self.assertIn(f"`tests/mechanism/test_overlap_second.py` holds {second} SYNTHETIC tests", method)
        self.assertIn(f"`tests/mechanism/test_learner_overlap_page.py` holds {here} tests", method)

    def test_every_figure_of_the_prose_stands_in_the_record(self):
        # The plan's constants (the 0.5 cutoff, 10, five, 1 to 4, fewer than 2), the release's revision and its
        # 31-mers; every other figure must stand in a printed table or in gates.txt.
        known = {"0", "0.5", "1", "2", "4", "5", "10", "14", "31"}
        for lines in lo.tables(str(RECORD)).values():
            for line in lines:
                known |= {n.replace(",", "") for n in re.findall(NUMBER, line)}
        for line in self.gates:
            known |= {n.replace(",", "") for n in re.findall(NUMBER, line)}
        head = "\n".join(line for line in self.head.split("\n") if not line.startswith("# "))
        texts = {"head": head}
        texts.update({s: self.sections[s] for s in ("Reading", "Overlap", "Yardstick", "Species", "Limits",
                                                    "Datasets")})
        for name, text in texts.items():
            for n in re.findall(NUMBER, prose(text)):
                self.assertIn(n.replace(",", ""), known, f"{name}: {n}")

    def test_the_words_keep_the_limits(self):
        for sentence in re.split(r"(?<=[.:])\s", self.text):
            for words in ("p-value", "null hypothesis", "would fail", "external validation"):
                if words in sentence:
                    self.assertRegex(sentence, r"\b(no|not|never|nothing|neither)\b", sentence)
            if "no difference" in sentence:
                self.assertNotIn("equal", sentence)
        self.assertNotIn("—", self.text)
        self.assertNotIn("<!--", self.text)
        self.assertNotRegex(self.text, r"\b(GCA|GCF)_\d|\bPRJ[A-Z]{2}\d|\b[cnpg]\.[A-Z-]?\d")
        # no gene or variant name: three lowercase letters and a capital (mecA, gyrA, blaKPC), or a point change
        self.assertNotRegex(self.text, r"\b(?:bla[A-Z]|[a-z]{3}[A-Z]\d*\b|[A-Z]\d{2,4}[A-Z]\b)")

    def test_the_skipped_reports_are_named(self):
        skipped = [line.split(": ", 1)[1] for line in self.gates if line.startswith("  skipped, ")]
        self.assertEqual(len(skipped), 6)
        sentence = self.sections["Datasets"].split("Skipped, with fewer than 2 learners on the lineage split: ")[1]
        for who in skipped:
            species, antibiotic = who.split(" with ")
            self.assertIn(f"*{species}* with {antibiotic}", sentence)

    def test_the_index_gives_the_page_one_row(self):
        index = INDEX.read_text(encoding="utf-8")
        self.assertEqual(index.count("| [learner-overlap-result.md](learner-overlap-result.md) |"), 1)

    def test_no_record_file_holds_an_id_shaped_token(self):
        self.assertEqual(lo.id_hits(sorted(str(p) for p in RECORD.iterdir())), 0)


if __name__ == "__main__":
    unittest.main()
