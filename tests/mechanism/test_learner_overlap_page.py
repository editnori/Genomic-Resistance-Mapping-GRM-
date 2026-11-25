"""docs/learner-overlap-result.md against the committed record of reference/learner-overlap/.

This test file reads committed result tables: the aggregates the table script
wrote from the one reading of the learner overlap, and those the second pass's
script wrote under second-pass/, per report, per species and pooled. It fits
nothing and reads no genome.
"""

import hashlib
import pathlib
import re
import unittest

from tests.mechanism import common

lo = common.load("learner_overlap")
os2 = common.load("overlap_second")
faults = common.load("overlap_faults")

RECORD = common.ROOT / "reference" / "learner-overlap"
SECOND = RECORD / "second-pass"
PAGE = common.ROOT / "docs" / "learner-overlap-result.md"
INDEX = common.ROOT / "docs" / "INDEX.md"
SECTIONS = ["Reading", "Overlap", "Yardstick", "Species", "Mean vote", "Null-fit cut", "Limits", "Method", "Plan",
            "Datasets"]
TABLES = {"Overlap": "overlap", "Yardstick": "yardstick", "Species": "species", "Datasets": "datasets"}
SECOND_TABLES = {"mean vote": "Mean vote", "mean vote by species": "Mean vote", "null-fit cut": "Null-fit cut"}
COUNTS = ["genomes", "resistant", "susceptible"] + os2.COUNT_COLUMNS
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
        self.votes = lo.read_tsv(SECOND / "vote.tsv")
        self.vote_pooled = lo.read_tsv(SECOND / "vote-pooled.tsv")
        self.second_gates = (SECOND / "gates.txt").read_text(encoding="utf-8").splitlines()
        self.method = self.sections["Method"].split("\n")

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

    def test_every_second_pass_table_stands_whole_in_its_section(self):
        tables = os2.tables(str(SECOND), str(RECORD))
        self.assertEqual(sorted(tables), sorted(SECOND_TABLES))
        for name, section in SECOND_TABLES.items():
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
        self.assertEqual(found, sorted(list(listed) + ["SHA256SUMS", "README.md", "reading.txt", "second-pass"]))
        listed = {}
        for line in (SECOND / "SHA256SUMS").read_text(encoding="utf-8").splitlines():
            digest, name = line.split("  ", 1)
            listed[name] = digest
        self.assertEqual(sorted(listed), sorted(os2.RECORD_FILES))
        for name, digest in listed.items():
            self.assertEqual(hashlib.sha256((SECOND / name).read_bytes()).hexdigest(), digest, name)
        found = sorted(p.name for p in SECOND.iterdir())
        self.assertEqual(found, sorted(list(listed) + ["SHA256SUMS", "README.md", "reading.txt"]))

    def test_the_counts_of_each_report_add_up(self):
        self.assertEqual(len(self.rows), len({(r["species"], r["antibiotic"]) for r in self.rows}))
        for r in self.rows:
            with self.subTest(report=(r["species"], r["antibiotic"])):
                n = {k: int(v) for k, v in r.items() if k not in ("species", "antibiotic", "learners", "catalog_rule",
                                                                   "expected_all5")}
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

    def test_the_yardstick_follows_from_the_miss_columns_and_no_ratio_is_a_column(self):
        # a report's ratio is given as the counts it comes from, never as a cell (the re-write of 25 November 2026)
        self.assertEqual(list(self.rows[0]), lo.DATASET_COLUMNS)
        self.assertEqual(list(self.votes[0]), os2.VOTE_COLUMNS)
        self.assertNotIn("ratio_all5", lo.DATASET_COLUMNS)
        self.assertNotIn("nonull_ratio_all5", os2.VOTE_COLUMNS)
        for r in self.rows:
            exp = lo.expected_all5(int(r["resistant"]), [int(r[f"miss_{x}"]) for x in FIVE])
            self.assertEqual(r["expected_all5"], lo.rec(exp), r["antibiotic"])

    def test_the_pooled_rows_are_what_the_script_pools(self):
        rows = [dict(r, learner_count=int(r["learner_count"]), resistant=int(r["resistant"])) for r in self.rows]
        again = [{k: str(v) for k, v in p.items()} for p in lo.pooled_rows(rows)]
        self.assertEqual(again, self.pooled)
        self.assertEqual(self.pooled[0]["datasets"], str(len(self.rows)))

    def test_the_mean_vote_columns_add_up_and_agree_with_the_first_pass(self):
        self.assertEqual([(v["species"], v["antibiotic"]) for v in self.votes],
                         [(r["species"], r["antibiotic"]) for r in self.rows])
        for v, r in zip(self.votes, self.rows):
            with self.subTest(report=(v["species"], v["antibiotic"])):
                n = {c: int(v[c]) for c in COUNTS}
                f = {k: int(x) for k, x in r.items() if k not in ("species", "antibiotic", "learners", "catalog_rule",
                                                                   "expected_all5")}
                for c in ("genomes", "resistant", "susceptible"):
                    self.assertEqual(n[c], f[c], c)
                for kind, first in (("miss", "miss"), ("wrong", "wrong")):
                    self.assertEqual(n[f"mean_{kind}"], sum(n[f"mean_{kind}_by_{b}"] for b in ("all5", "1to4", "none")))
                    self.assertEqual(n[f"mean_{kind}_by_none"], 0)
                    self.assertEqual(n[f"mean_{kind}_by_all5"], f[f"{first}_5"])
                    self.assertLessEqual(n[f"mean_{kind}_by_1to4"], sum(f[f"{first}_{k}"] for k in range(1, 5)))
                self.assertEqual(n["nonull_miss_all5"] + f["miss_all5_in_null_fold"], f["miss_all5"])
                self.assertLessEqual(n["nonull_resistant"], n["resistant"])
                for x in FIVE:
                    self.assertLessEqual(n["nonull_miss_all5"], n[f"nonull_miss_{x}"])
                    self.assertLessEqual(n[f"nonull_miss_{x}"], f[f"miss_{x}"])
                exp = lo.expected_all5(n["nonull_resistant"], [n[f"nonull_miss_{x}"] for x in FIVE])
                self.assertEqual(v["nonull_expected_all5"], lo.rec(exp))

    def test_the_second_pooled_rows_are_what_the_script_pools(self):
        again = [{k: str(x) for k, x in p.items()} for p in os2.pooled_rows(self.votes, str(RECORD))]
        self.assertEqual(again, self.vote_pooled)
        self.assertEqual(self.vote_pooled[0]["datasets"], self.pooled[0]["datasets"])

    def test_the_gates_all_pass_and_stand_in_the_method(self):
        checks = [line for line in self.gates if line.startswith("check ")]
        self.assertEqual([line.split(" ", 2)[1] for line in checks], [str(n) for n in range(1, 11)])
        self.assertFalse(any("STOP" in line or "failed" in line for line in self.gates))
        for line in checks + [self.gates[-1]]:
            self.assertEqual(self.method.count(line), 1, line)
        self.assertEqual(self.gates[-1], "reading rule: stack")
        # the rule reads the exact median, worked from the count columns, and the record's 3 places agree with it
        exact = lo.exact_pool(lo.scopes(self.rows)[0][2])["ratio_median"]
        self.assertEqual((lo.reading_rule(exact), lo.rec(exact)), ("stack", self.pooled[0]["ratio_median"]))

    def test_the_second_gates_all_pass_and_stand_in_the_method(self):
        second = [line for line in self.second_gates if line.startswith("second check ")]
        self.assertEqual([line.split(" ", 3)[2] for line in second], [str(n) for n in range(1, 7)])
        self.assertFalse(any("STOP" in line or "failed" in line for line in self.second_gates))
        # the first pass's checks, run again, print the lines they printed; its check 6 runs at the write
        self.assertEqual([line for line in self.second_gates if line.startswith("check ")],
                         [line for line in self.gates if line.startswith("check ") and not line.startswith("check 6 ")])
        for line in second + [self.second_gates[-1]]:
            self.assertEqual(self.method.count(line), 1, line)
        self.assertEqual(self.second_gates[-1], "reading rule outside null-fit folds: stack")
        groups = os2.scopes(self.votes, str(RECORD))
        exact = os2.exact_pool(groups[0][2], os2.first_rows(str(RECORD)))["nonull_ratio_median"]
        self.assertEqual((lo.reading_rule(exact), lo.rec(exact)), ("stack", self.vote_pooled[0]["nonull_ratio_median"]))

    def reading_and_rewrite(self, directory, files):
        """reading.txt: the reading's commands and outputs, then the re-write's, which ran write and print again
        on the kept reading, unchanged, and no read; the committed record is the re-write's."""
        first, mark, rewrite = (directory / "reading.txt").read_text(encoding="utf-8").partition("\n# the re-write, ")
        self.assertTrue(mark)
        first, rewrite = first.splitlines(), rewrite.splitlines()
        commands = [line for line in first if line.startswith("python3 -I ")]
        self.assertEqual(len(commands), 4)
        self.assertIn("```\n" + "\n".join(commands) + "\n```", self.sections["Method"])
        self.assertEqual([line for line in rewrite if line.startswith("python3 -I ")],
                         [c for c in commands if " write " in c or " print " in c])
        listed = {line.split("  ", 1)[1]: line.split("  ", 1)[0] for line in first if "  " in line}
        relisted = {line.split("  ", 1)[1]: line.split("  ", 1)[0] for line in rewrite if "  " in line}
        self.assertEqual(relisted["rewrite/reading/reading.json"], listed["reading/reading.json"])
        for name in files + ("SHA256SUMS",):
            self.assertEqual(relisted[f"rewrite/record/{name}"],
                             hashlib.sha256((directory / name).read_bytes()).hexdigest(), name)
        return commands

    def test_the_reading_commands_stand_as_reading_txt_gives(self):
        self.reading_and_rewrite(RECORD, lo.RECORD_FILES)

    def test_the_second_reading_commands_stand_as_its_reading_txt_gives(self):
        commands = self.reading_and_rewrite(SECOND, os2.RECORD_FILES)
        self.assertTrue(all(" scripts/mechanism/overlap_second.py " in line for line in commands))

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
        # 31-mers; every other figure must stand in a printed table or in a gates.txt.
        known = {"0", "0.5", "1", "2", "4", "5", "10", "14", "31"}
        for lines in lo.tables(str(RECORD)).values():
            for line in lines:
                known |= {n.replace(",", "") for n in re.findall(NUMBER, line)}
        for lines in os2.tables(str(SECOND), str(RECORD)).values():
            for line in lines:
                known |= {n.replace(",", "") for n in re.findall(NUMBER, line)}
        for line in self.gates + self.second_gates:
            known |= {n.replace(",", "") for n in re.findall(NUMBER, line)}
        head = "\n".join(line for line in self.head.split("\n") if not line.startswith("# "))
        texts = {"head": head}
        texts.update({s: self.sections[s] for s in ("Reading", "Overlap", "Yardstick", "Species", "Mean vote",
                                                    "Null-fit cut", "Limits", "Datasets")})
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
        # every cell of the four tables but those of the computed decimal columns, which have 3 places and so may
        # have the shape, and the other eight files whole; the decimal columns are left to the id list on the box
        files = lo.files_under(str(RECORD))
        self.assertEqual(len(files), 12)
        tables = {RECORD / "datasets.tsv": lo.DECIMAL_COLUMNS, RECORD / "pooled.tsv": lo.DECIMAL_COLUMNS,
                  SECOND / "vote.tsv": os2.DECIMAL_COLUMNS, SECOND / "vote-pooled.tsv": os2.DECIMAL_COLUMNS}
        self.assertEqual(lo.id_hits([p for p in files if pathlib.Path(p) not in tables]), 0)
        for path, columns in tables.items():
            hits, skipped = lo.record_hits(str(path), columns)
            self.assertEqual(hits, 0, path.name)
            self.assertGreater(skipped, 0, path.name)

    def test_every_decimal_of_the_record_has_3_places(self):
        for path, columns in ((RECORD / "datasets.tsv", lo.DECIMAL_COLUMNS), (RECORD / "pooled.tsv", lo.DECIMAL_COLUMNS),
                              (SECOND / "vote.tsv", os2.DECIMAL_COLUMNS),
                              (SECOND / "vote-pooled.tsv", os2.DECIMAL_COLUMNS)):
            for r in lo.read_tsv(str(path)):
                for c in columns:
                    if c in r and r[c] != "":
                        self.assertRegex(r[c], r"^\d+\.\d{3}$", f"{path.name} {c}")


if __name__ == "__main__":
    unittest.main()
