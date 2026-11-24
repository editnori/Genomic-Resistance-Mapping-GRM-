"""docs/lineage-curve-result.md against the committed tables of reference/lineage-curve/.

This test file reads committed result tables: the aggregates the table script
wrote from the one reading of block M4, by antibiotic, lineage group and
gene. It fits nothing and reads no genome.
"""

import csv
import hashlib
import json
import re
import unittest

from tests.mechanism import common

curve_tables = common.load("curve_tables")
curve_faults = common.load("curve_faults")

RECORD = common.ROOT / "reference" / "lineage-curve"
PAGE = common.ROOT / "docs" / "lineage-curve-result.md"
LISTS = common.ROOT / "scripts" / "c3-inputs" / "lists.tsv"
SECTIONS = ["Reading", "Data", "Curves", "No graded variant", "Spread over lineages", "Genes", "Seen after", "Limits",
            "Method", "Datasets"]
TABLES = {"Reading": ["reading"], "Curves": ["curves_main", "curve_ends"], "No graded variant": ["without"],
          "Spread over lineages": ["spread"], "Genes": ["genes"], "Datasets": ["datasets"]}
d3, count = curve_tables.d3, curve_tables.count_words


def rows(name):
    with open(RECORD / name, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sections(text):
    parts = re.split(r"^## (.+)$", text, flags=re.M)
    return dict(zip(parts[1::2], parts[2::2]))


def blocks(text):
    return [b for b in re.split(r"\n\s*\n", text) if b.strip()]


def checks():
    return dict(line.split(": ", 1) for line in (RECORD / "checks.txt").read_text(encoding="utf-8").splitlines()
                if ": " in line and not line.startswith("check: "))


class Page(unittest.TestCase):
    def setUp(self):
        self.text = PAGE.read_text(encoding="utf-8")
        self.lines = self.text.split("\n")
        self.sections = sections(self.text)
        self.endpoints = rows("page-endpoints.tsv")
        self.scored = [r for r in self.endpoints if r["kind"] != "none"]

    def test_the_sections_stand_in_order(self):
        self.assertEqual(list(self.sections), SECTIONS)

    def test_every_printed_table_stands_whole_in_its_section(self):
        tables = curve_tables.tables(RECORD)
        self.assertEqual(sorted(tables), sorted(t for names in TABLES.values() for t in names))
        for section, names in TABLES.items():
            for name in names:
                printed = tables[name]
                self.assertGreater(len(printed), 2, name)
                self.assertIn("\n".join(printed) + "\n", self.sections[section], name)
                for line in printed[:1] + printed[2:]:
                    self.assertEqual(self.lines.count(line), 1, f"{name}: {line}")

    def test_the_reports_are_the_species_reports_named_in_full(self):
        want = sorted((r["species"], r["antibiotic"]) for r in curve_tables.load(LISTS)
                      if r["species_slug"] == "mycobacterium_tuberculosis")
        self.assertEqual(sorted((r["species"], r["antibiotic"]) for r in self.endpoints), want)
        datasets = self.sections["Datasets"]
        for r in self.endpoints:
            self.assertIn(f"| *{r['species']}* | {r['antibiotic']} |", datasets)

    def test_the_run_checks_stand_as_checks_txt_gives(self):
        c = checks()
        text = (RECORD / "checks.txt").read_text(encoding="utf-8")
        self.assertNotIn("FAIL", text)
        self.assertEqual(text.count("check: pass: "), len(curve_tables.checks_pass(c)))
        method = self.sections["Method"]
        for key in ("planned", "results", "failed_lines"):
            self.assertIn(f"| {key.replace('_', ' ')} | {count(c[key])} |", method)
        self.assertIn(f"| wall time | {count(c['wall_seconds'])} s, from {c['start']} to {c['end']} |", method)
        self.assertIn(f"| failed, by reason | {c['failed_by_reason']} |", method)

    def test_the_gates_stand_as_gates_txt_gives(self):
        gates = (RECORD / "gates.txt").read_text(encoding="utf-8").splitlines()
        self.assertEqual(len(gates), len(curve_tables.GATES))
        for n, line in enumerate(gates, 1):
            name, words = line.split(" ", 1)[1].split(": ", 1)
            self.assertEqual(self.lines.count(f"| {n}. {name} | {words} |"), 1, line)
        self.assertIn("checked the reading against seven gates and wrote `reference/lineage-curve/`. All seven passed:",
                      self.sections["Method"])

    def test_the_reading_commands_stand_as_reading_txt_gives(self):
        reading = (RECORD / "reading.txt").read_text(encoding="utf-8").splitlines()
        method = self.sections["Method"]
        for line in reading:
            if line.startswith("python3 -I "):
                self.assertIn(line, method)
        listed = {line.split("  ", 1)[1]: line.split("  ", 1)[0] for line in reading if "  " in line}
        for name in curve_tables.M4_FILES:
            self.assertEqual(listed[name], hashlib.sha256((RECORD / "m4" / name).read_bytes()).hexdigest(), name)

    def test_the_fault_table_names_each_planted_faults_test(self):
        method = self.sections["Method"]
        for fault, words, test, _ in curve_faults.FAULTS:
            self.assertEqual(self.lines.count(f"| {fault}. {words[0].upper()}{words[1:]} | "
                                              f"`{test.split('test_curve_tables.')[1]}` |"), 1, fault)
        self.assertIn(f"All {len(curve_faults.FAULTS)} were caught.", method)

    def test_the_page_gives_the_suites_own_test_counts(self):
        loader = unittest.TestLoader()
        here = loader.loadTestsFromName("tests.mechanism.test_lineage_page").countTestCases()
        tables = loader.loadTestsFromName("tests.mechanism.test_curve_tables").countTestCases()
        method = self.sections["Method"]
        self.assertIn(f"`tests/mechanism/test_curve_tables.py` holds {tables} SYNTHETIC tests", method)
        self.assertIn(f"`tests/mechanism/test_lineage_page.py` holds {here} tests of this page", method)

    def test_the_summary_and_the_endpoints_agree(self):
        with open(RECORD / "m4" / "m4-summary.json", encoding="utf-8") as handle:
            summary = json.load(handle)
        c = checks()
        self.assertEqual((summary["tool"], summary["database_commit"]), ("6.6.3", "3ed1cb99"))
        self.assertEqual((summary["genomes_planned"], summary["genomes_with_result"], summary["genomes_failed"]),
                         (c["planned"], c["results"], c["failed_genomes"]))
        ends = {e["antibiotic"]: e for e in summary["endpoints"]}
        for r in self.scored:
            self.assertEqual(r["resistant_with_result"], ends[r["antibiotic"]]["resistant_with_result"])
            self.assertEqual(r["in_M"], "1" if ends[r["antibiotic"]]["in_M"] else "0")

    def test_every_figure_of_the_prose_stands_in_the_record(self):
        # The rules' own constants, the run's dates and the release's revision; every other figure must stand in a
        # printed table or count the record's reports, groups or genomes.
        fixed = {"2.5", "97.5", "5", "0.1", "0.8", "14", "2026", "24"}
        with open(RECORD / "m4" / "m4-summary.json", encoding="utf-8") as handle:
            summary = json.load(handle)
        coverage = rows("m4/m4-coverage.tsv")
        known = set(fixed) | {str(summary["orders"]), summary["genomes_planned"], str(len(self.endpoints)),
                              str(len(self.scored)), str(len({r["group"] for r in coverage
                                                              if r["grouping"] == "main_lineage"}))}
        known |= {str(sum(r["kind"] == k for r in self.scored)) for k in ("primary", "proxy")}
        for lines in curve_tables.tables(RECORD).values():
            for line in lines:
                known |= {n.replace(",", "") for n in re.findall(r"\d[\d,]*(?:\.\d+)?", line)}
        for section in ("Reading", "Data", "Curves", "No graded variant", "Spread over lineages", "Genes", "Seen after"):
            prose = "\n".join(line for line in self.sections[section].split("\n") if not line.startswith("|"))
            prose = re.sub(r"lineage\d|La\d", "", prose.replace(summary["tool"], "").replace(summary["database_commit"], ""))
            for n in re.findall(r"\d[\d,]*(?:\.\d+)?", prose):
                self.assertIn(n.replace(",", ""), known, f"{section}: {n}")

    def test_the_words_keep_the_limits(self):
        for sentence in re.split(r"(?<=[.:])\s", self.text):
            if "external validation" in sentence:
                self.assertRegex(sentence, r"\b(never|not)\b", sentence)
            if "prevalence" in sentence:
                self.assertRegex(sentence, r"\b(no|never|not)\b", sentence)
            self.assertNotRegex(sentence, r"(?<!no )\b(causes?|caused|explains?|explained|drives?|driven)\b"
                                          r"(?! resistance in a)", sentence)
        self.assertNotIn("—", self.text)
        self.assertNotIn("<!--", self.text)
        self.assertNotRegex(self.text, r"\b(GCA|GCF)_\d|\bPRJ[A-Z]{2}\d|\b[cnpg]\.[A-Z-]?\d")

    def test_the_record_holds_the_files_its_checksums_give(self):
        listed = {}
        for line in (RECORD / "SHA256SUMS").read_text(encoding="utf-8").splitlines():
            digest, name = line.split("  ", 1)
            listed[name] = digest
        found = sorted(p.relative_to(RECORD).as_posix() for p in RECORD.rglob("*")
                       if p.is_file() and p.name != "SHA256SUMS")
        self.assertEqual(sorted(listed), found)
        for name, digest in listed.items():
            self.assertEqual(hashlib.sha256((RECORD / name).read_bytes()).hexdigest(), digest, name)


if __name__ == "__main__":
    unittest.main()
