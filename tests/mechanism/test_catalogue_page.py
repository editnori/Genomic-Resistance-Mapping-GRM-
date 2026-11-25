"""docs/catalogue-size-result.md against the committed tables of reference/catalogue-size/.

This test file reads committed result tables: the aggregates the catalogue
size reading wrote, one row per antibiotic and one per set and pair, and the
lineage curve's reference/lineage-curve/page-endpoints.tsv they were joined
to. It fits nothing and reads no genome.
"""

import csv
import hashlib
import re
import unittest
from decimal import Decimal
from fractions import Fraction

from tests.mechanism import common

cs = common.load("catalogue_size")
faults = common.load("catalogue_faults")

RECORD = common.ROOT / "reference" / "catalogue-size"
CURVE = common.ROOT / "reference" / "lineage-curve"
PAGE = common.ROOT / "docs" / "catalogue-size-result.md"
SECTIONS = ["Rule", "Reading", "Antibiotics", "Correlations", "Limits", "Method"]
TABLES = {"Antibiotics": ["antibiotics"], "Correlations": ["correlations", "sets"], "Method": ["guards"]}
ANSWERS = {"tracks resistant genomes": "catalogue size tracks how many resistant genomes there were to study",
           "tracks share": "the number of graded variants tracks the resistant share",
           "neither by the margin": "these numbers cannot separate the two"}
OUTCOMES = "The rule declared three outcomes before the read."
ORDINALS = ("first", "second", "third")


def rows(folder, name):
    with open(folder / name, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t", quoting=csv.QUOTE_NONE))


def sections(text):
    parts = re.split(r"^## (.+)$", text, flags=re.M)
    return dict(zip(parts[1::2], parts[2::2]))


class Page(unittest.TestCase):
    def setUp(self):
        self.text = PAGE.read_text(encoding="utf-8")
        self.lines = self.text.split("\n")
        self.sections = sections(self.text)
        self.table = rows(RECORD, "antibiotics.tsv")
        self.correlations = rows(RECORD, "correlations.tsv")
        self.rule = {c["x"]: c for c in self.correlations if c["set"] == "rule"}

    def gap(self):
        """By how much the share's correlation passes the resistant genomes', in absolute value, to 3 decimals."""
        return cs.d3(abs(Decimal(self.rule["resistant_share"]["rho"])) - abs(Decimal(self.rule["resistant"]["rho"])))

    def test_the_sections_stand_in_order(self):
        self.assertEqual(list(self.sections), SECTIONS)

    def test_every_printed_table_stands_whole_in_its_section(self):
        tables = cs.tables(RECORD)
        self.assertEqual(sorted(tables), sorted(t for names in TABLES.values() for t in names))
        for section, names in TABLES.items():
            for name in names:
                printed = tables[name]
                self.assertGreater(len(printed), 2, name)
                self.assertIn("\n".join(printed) + "\n", self.sections[section], name)
                for line in printed[:1] + printed[2:]:
                    self.assertEqual(self.lines.count(line), 1, f"{name}: {line}")

    def test_the_table_carries_the_lineage_curves_columns_unchanged(self):
        page = rows(CURVE, "page-endpoints.tsv")
        self.assertEqual([(r["species"], r["antibiotic"]) for r in self.table],
                         [(r["species"], r["antibiotic"]) for r in page])
        for mine, theirs in zip(self.table, page):
            for key in cs.CARRIED:
                self.assertEqual(mine[key], theirs[key], (mine["antibiotic"], key))
            self.assertEqual(int(mine["resistant"]) + int(mine["susceptible"]), int(mine["genomes"]))
            self.assertEqual(mine["resistant_share"], cs.f12(Fraction(int(mine["resistant"]), int(mine["genomes"]))))
            if theirs["kind"] != "none":
                self.assertEqual(int(mine["resistant"]),
                                 int(theirs["resistant_with_result"]) + int(theirs["resistant_no_result"]))
            self.assertEqual(mine["left_out_because"], cs.left_out_because(theirs))

    def test_the_correlations_follow_from_the_table(self):
        for c in self.correlations:
            chosen = cs.set_rows([dict(r, _share=Fraction(int(r["resistant"]), int(r["genomes"]))) for r in self.table],
                                 c["set"])
            self.assertEqual(c["antibiotics"], "; ".join(r["antibiotic"] for r in chosen))
            self.assertEqual(int(c["n"]), len(chosen))
            xs, ys = [cs.value(r, c["x"]) for r in chosen], [int(r["catalogue_variants"]) for r in chosen]
            self.assertEqual(c["rho"], cs.f12(cs.spearman(xs, ys)), (c["set"], c["x"]))
            low, low_without, high, high_without = cs.leave_one_out(xs, ys, [r["antibiotic"] for r in chosen])
            self.assertEqual((c["loo_min"], c["loo_min_without"], c["loo_max"], c["loo_max_without"]),
                             (cs.f12(low), low_without, cs.f12(high), high_without), (c["set"], c["x"]))
        self.assertEqual({(c["set"], c["role"]) for c in self.correlations},
                         {("rule", "the rule's set"), ("primary only", "beside, decides nothing")})
        self.assertEqual({c["role"] for c in self.correlations if c["decides"] == "yes"}, {"the rule's set"})
        for name in cs.SETS:
            mine = {c["x"]: c for c in self.correlations if c["set"] == name}
            gives = cs.branch(float(mine["resistant"]["rho"]), float(mine["resistant_share"]["rho"]))
            self.assertEqual({c["comparison_gives"] for c in mine.values()}, {gives}, name)

    def test_the_reading_gives_the_answer_the_rule_selected(self):
        selected = self.rule["resistant"]["comparison_gives"]
        self.assertEqual(self.rule["resistant_share"]["comparison_gives"], selected)
        reading = self.sections["Reading"]
        self.assertIn(ANSWERS[selected], reading)
        self.assertIn(ANSWERS[selected], cs.RULE)
        outcomes = [p for p in reading.split("\n\n") if p.strip().startswith(OUTCOMES)]
        self.assertEqual(len(outcomes), 1)
        places = [outcomes[0].find(ANSWERS[b]) for b in cs.BRANCHES]
        self.assertTrue(all(i >= 0 for i in places) and places == sorted(places), places)
        self.assertIn(f"The numbers selected the {ORDINALS[cs.BRANCHES.index(selected)]}", outcomes[0])
        self.assertIn(f"{self.gap()} against the {cs.MARGIN}", outcomes[0])
        rest = reading.replace(outcomes[0], "")
        for branch, words in ANSWERS.items():
            if branch != selected:
                self.assertNotIn(words, rest, branch)
        first = reading.strip().split("\n")[0]
        self.assertIn(f"Over the {self.rule['resistant']['n']} antibiotics the rule reads", first)
        self.assertIn(f"{cs.d3(self.rule['resistant_share']['rho'])} with the resistant share", first)
        self.assertIn(f"{cs.d3(self.rule['resistant']['rho'])} with the resistant genomes", first)

    def test_every_figure_of_the_prose_stands_in_the_record(self):
        # the rule's margin, the floor of 5 the lineage curve names genes by and the species' cap of 2,000 deposits;
        # the tool's version and database commit are taken out; every other figure must stand in a printed table or
        # count the record's reports or left-out rows
        fixed = {"0.15", "5", "2000"}
        known = set(fixed) | {self.gap(), str(len(self.table)), str(sum(r["in_correlation"] == "1" for r in self.table))}
        known |= {str(sum(r["left_out_because"] == why for r in self.table)) for why in cs.REASONS}
        for lines in cs.tables(RECORD).values():
            for line in lines:
                known |= {n.replace(",", "") for n in re.findall(r"\d[\d,]*(?:\.\d+)?", line)}
        for section in ("Reading", "Antibiotics", "Correlations", "Limits"):
            prose = "\n".join(line for line in self.sections[section].split("\n") if not line.startswith("|"))
            prose = prose.replace("6.6.3", "").replace("3ed1cb99", "")
            for n in re.findall(r"\d[\d,]*(?:\.\d+)?", prose):
                self.assertIn(n.replace(",", ""), known, f"{section}: {n}")

    def test_the_left_out_rows_are_named_with_their_counts(self):
        text = self.sections["Antibiotics"]
        for why, words in zip(cs.REASONS, ("are not read by the tool", "have no graded variant", "are in set M")):
            names = [r["antibiotic"] for r in self.table if r["left_out_because"] == why]
            self.assertIn(f"{len(names)} {words} ({' and '.join(names)})", text)
        self.assertIn(f"Of the {len(self.table)} reports", text)
        self.assertIn(f"which read the other {sum(r['in_correlation'] == '1' for r in self.table)}.", text)

    def test_the_reading_commands_stand_as_reading_txt_gives(self):
        reading = (RECORD / "reading.txt").read_text(encoding="utf-8").splitlines()
        method = self.sections["Method"]
        commands = [line for line in reading if line.startswith("python3 -I ")]
        self.assertEqual(len(commands), 1)
        self.assertIn(commands[0], method)
        self.assertEqual(sorted(line.split("  ", 1)[1] for line in reading if "  " in line),
                         ["read.txt", "release-counts.tsv"])

    def test_the_fault_table_names_each_planted_faults_test(self):
        method = self.sections["Method"]
        for fault, words, test, _ in faults.FAULTS:
            self.assertEqual(self.lines.count(f"| {fault}. {words[0].upper()}{words[1:]} | "
                                              f"`{test.split('test_catalogue_size.')[1]}` |"), 1, fault)
        self.assertIn(f"All {len(faults.FAULTS)} were caught.", method)

    def test_the_page_gives_the_suites_own_test_counts(self):
        loader = unittest.TestLoader()
        here = loader.loadTestsFromName("tests.mechanism.test_catalogue_page").countTestCases()
        script = loader.loadTestsFromName("tests.mechanism.test_catalogue_size").countTestCases()
        method = self.sections["Method"]
        self.assertIn(f"`tests/mechanism/test_catalogue_size.py` holds {script} SYNTHETIC tests", method)
        self.assertIn(f"`tests/mechanism/test_catalogue_page.py` holds {here} tests of this page", method)

    def test_the_words_keep_the_limits(self):
        genes = {r["gene"] for r in rows(CURVE / "m4", "m4-genes.tsv")} - {"other genes"}
        self.assertTrue(genes)
        for gene in genes:
            self.assertNotRegex(self.text, rf"(?<![A-Za-z]){re.escape(gene)}(?![A-Za-z])", gene)
        for sentence in re.split(r"(?<=[.:])\s", self.text):
            if any(word in sentence for word in ("prevalence", "safer", "causal")):
                self.assertRegex(sentence, r"\b(no|never|not|Nothing)\b", sentence)
            self.assertNotRegex(sentence, r"\b(causes?|caused|explains?|explained|drives?|driven)\b", sentence)
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
        self.assertEqual((RECORD / "README.md").read_text(encoding="utf-8"), cs.README)


if __name__ == "__main__":
    unittest.main()
