"""docs/mechanism-year-result.md against the committed tables of reference/mechanism-year/.

This test file reads committed result tables: the aggregates year_mix.py wrote
from the one reading, by year, group, species and antibiotic, and M1's
committed per-dataset table, which no-call-strata.tsv sums. It fits nothing
and reads no genome.
"""

import csv
import hashlib
import re
import unittest

from tests.mechanism import common

ym = common.load("year_mix")
year_faults = common.load("year_faults")

RECORD = common.ROOT / "reference" / "mechanism-year"
PAGE = common.ROOT / "docs" / "mechanism-year-result.md"
LISTS = common.ROOT / "scripts" / "c3-inputs" / "lists.tsv"
M1 = common.ROOT / "reference" / "mechanism-split" / "page-m1-strata.tsv"
M1_DATASETS = common.ROOT / "reference" / "mechanism-split" / "page-m1-datasets.tsv"
SECTIONS = ["Reading", "Plan", "Step 1", "No call by stratum", "Groups", "Years", "Species", "Seen after", "Limits",
            "Method", "Datasets"]
TABLES = {"Step 1": ["step1"], "No call by stratum": ["no_call"], "Groups": ["totals"], "Years": ["years", "ranges"],
          "Datasets": ["datasets"]}
STRATA_COMMAND = ("python3 -I scripts/mechanism/year_mix.py strata --m1 reference/mechanism-split/page-m1-datasets.tsv "
                  "--lists scripts/c3-inputs/lists.tsv --record reference/mechanism-year")
FOUR = ["point_only", "gene_only", "both", "neither"]
d3, count = ym.d3, ym.count_words


def rows(path):
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sections(text):
    parts = re.split(r"^## (.+)$", text, flags=re.M)
    return dict(zip(parts[1::2], parts[2::2]))


def m1_rows():
    return [r for r in rows(M1) if (r["set"], r["split"], r["kind"]) == ("primary", "lineage", "stratum")]


def m1_no_call():
    """Recomputed here from M1's committed per-dataset table, apart from the script: per stratum, the
    resistant genomes and those the catalog row calls, on the lineage split of the primary set."""
    out = {s: [0, 0] for s in FOUR}
    for r in rows(M1_DATASETS):
        if r["set"] == "primary" and r["split"] == "lineage" and r["row"] == "catalog":
            for s in FOUR:
                out[s][0] += int(r[f"{s}_resistant"])
                out[s][1] += int(r[f"{s}_called"])
    return out


def guards():
    return (RECORD / "guards.txt").read_text(encoding="utf-8")


def prose(section):
    """A section's sentences: no table line, no code block and nothing in backticks."""
    text = re.sub(r"```.*?```", "", section, flags=re.S)
    text = "\n".join(line for line in text.split("\n") if not line.startswith("|"))
    return re.sub(r"`[^`]*`", "", text)


class Page(unittest.TestCase):
    def setUp(self):
        self.text = PAGE.read_text(encoding="utf-8")
        self.lines = self.text.split("\n")
        self.sections = sections(self.text)
        self.summary = ym.summary(str(RECORD))
        self.totals = {r["stratum"]: r for r in rows(RECORD / "totals.tsv")}
        self.step1 = {r["stratum"]: r for r in rows(RECORD / "step1.tsv")}

    def test_the_sections_stand_in_order(self):
        self.assertEqual(list(self.sections), SECTIONS)

    def test_every_printed_table_stands_whole_in_its_section(self):
        tables = ym.page_tables(str(RECORD))
        self.assertEqual(sorted(tables), sorted(t for names in TABLES.values() for t in names))
        for section, names in TABLES.items():
            for name in names:
                printed = tables[name]
                self.assertGreater(len(printed), 2, name)
                self.assertIn("\n".join(printed) + "\n", self.sections[section], name)
                for line in printed[:1] + printed[2:]:
                    self.assertEqual(self.lines.count(line), 1, f"{name}: {line}")

    def test_the_page_says_the_four_strata_were_not_readable(self):
        self.assertIn("no", [self.step1[s]["equal"] for s in FOUR])
        self.assertEqual(self.summary["strata"], ["reader_call", "reader_no_call"])
        self.assertIn("The four strata of block M1 were not readable from the release reports", self.text)
        self.assertIn("The four strata of M1 were not readable from the release reports.", self.sections["Reading"])
        g = self.step1["gene_only"]
        self.assertIn(f"{count(g['rebuilt_resistant'])} resistant genomes in gene only against M1's "
                      f"{count(g['published_resistant'])}.", self.sections["Reading"])

    def test_the_reading_gives_the_figures_the_record_gives(self):
        s, reading = self.summary, self.sections["Reading"]
        nc, call = self.totals["reader_no_call"], self.totals["reader_call"]
        self.assertEqual(self.totals["reader_no_rule"]["genomes"], "0")
        everyone = self.step1["all_strata"]
        self.assertEqual(int(nc["resistant"]) + int(call["resistant"]), int(everyone["published_resistant"]))
        self.assertIn(f"calls {count(nc['resistant_called_mean5'])} of the reader's {count(nc['resistant'])} no-call "
                      f"resistant genomes (recall {d3(nc['recall_mean5'])}), and {count(call['resistant_called_mean5'])} "
                      f"of the {count(call['resistant'])} resistant genomes the reader calls ({d3(call['recall_mean5'])})"
                      ". This reading counts genomes and reports recall: it fits nothing and reports no loss.", reading)
        self.assertIn(f"Among the {len(s['kept_years'])} years that are not thin, {ym.year_runs(s['kept_years'])}, ",
                      reading)
        g = s["groups"]["reader_no_call"]
        low, high = g["resistant_share"]
        self.assertIn(f"runs from {ym.share_words(low, 'resistant')} to {ym.share_words(high, 'resistant')}.", reading)
        low, high = g["recall"]
        self.assertIn(f"recall runs from {ym.recall_words(low)} to {ym.recall_words(high)}, over the "
                      f"{g['years_with_recall']} of those years", reading)
        g = s["groups"]["reader_call"]
        low, high = g["recall"]
        self.assertEqual(g["years_with_recall"], len(s["kept_years"]))
        self.assertIn(f"it runs from {ym.recall_words(low)} to {ym.recall_words(high)}, over all "
                      f"{len(s['kept_years'])}.", reading)
        self.assertIn("compares the AMRFinderPlus program's calls between genomes collected in 2010 to 2014 and in "
                      "2015 to 2019, one vote per BioProject. This read differs at each step: every year stands on its "
                      "own, with no vote and no BioProject grouping, and the groups are the reader's call rather than "
                      "the four strata.", reading)

    def test_the_rows_with_no_year_come_first_and_are_not_assumed_alike(self):
        s, reading = self.summary, self.sections["Reading"]
        everyone = self.step1["all_strata"]
        genomes = int(everyone["published_resistant"]) + int(everyone["published_susceptible"])
        self.assertEqual(s["genomes"], genomes)
        self.assertEqual(sum(int(r["genomes"]) for r in self.totals.values()), genomes)
        self.assertGreater(3 * s["unknown_genomes"], genomes)
        used = re.search(r"guard 5: (\d+) used", guards()).group(1)
        paragraph = (f"**Rows with no year.** Of the {count(genomes)} genome rows of the {used} primary datasets, a "
                     f"genome once per dataset, {count(s['unknown_genomes'])} carry no four-digit year: over a third. "
                     "Every year row below describes the genomes that carry a year, and the rows without one are not "
                     "assumed to be like them.")
        self.assertIn(paragraph, reading)
        self.assertLess(reading.index(paragraph), reading.index("Among the "))

    def test_the_two_no_call_counts_are_related_stratum_by_stratum(self):
        m1 = m1_no_call()
        table = {r["stratum"]: r for r in rows(RECORD / ym.NO_CALL_TABLE)}
        self.assertEqual(list(table), FOUR + ["program_row", "all_strata"])
        expect = {s: m1[s] for s in FOUR}
        expect["program_row"] = [sum(m1[s][i] for s in FOUR[:3]) for i in (0, 1)]
        expect["all_strata"] = [sum(m1[s][i] for s in FOUR) for i in (0, 1)]
        for name, (resistant, called) in expect.items():
            r = table[name]
            self.assertEqual((int(r["resistant"]), int(r["called"]), int(r["no_call"])),
                             (resistant, called, resistant - called), name)
            self.assertAlmostEqual(float(r["no_call_share"]), (resistant - called) / resistant, places=12, msg=name)
        for s in FOUR:
            self.assertEqual(str(m1[s][0]), self.step1[s]["published_resistant"], s)
        no = {s: count(m1[s][0] - m1[s][1]) for s in FOUR}
        yes = {s: count(m1[s][1]) for s in FOUR}
        held = {s: count(m1[s][0]) for s in FOUR}
        everyone = count(expect["all_strata"][0])
        no_call = count(self.totals["reader_no_call"]["resistant"])
        gap = count(expect["program_row"][0] - expect["program_row"][1])
        self.assertIn(
            f"**Two counts with no call, related stratum by stratum.** At the program's stage, the AMRFinderPlus "
            f"program reports no row that meets the antibiotic's rule, neither a point mutation nor an acquired gene, "
            f"for {held['neither']} of the {everyone} resistant genomes: M1's `neither` stratum. At the reader's stage, "
            f"GRM's known-gene rule makes no call on {no_call} of the same {everyone}: `reader_no_call` here, and the "
            "count a user of the release meets", self.sections["Reading"])
        self.assertIn(
            f"The {no_call} are not the {held['neither']} plus a remainder: {no['neither']} of the {held['neither']} are "
            f"among them, and the rule calls the other {yes['neither']} resistant, so the two counts cannot be "
            f"subtracted. The gap between the stages, resistant genomes with a program row that meets the rule on which "
            f"the rule makes no call, is {gap}, and {no['point_only']} of them are point only. The rule calls "
            f"{yes['point_only']} of the {held['point_only']} point-only resistant genomes, against {yes['gene_only']} "
            f"of the {held['gene_only']} gene-only and {yes['both']} of the {held['both']} with both: GRM's known-gene "
            "rule is in effect a gene-presence rule, and it almost never fires on a point mutation.",
            self.sections["Reading"])
        method = self.sections["Method"]
        self.assertIn("```\n" + STRATA_COMMAND + "\n```", method)
        self.assertIn("Its only new input is `page-m1-datasets.tsv`, a public table committed with M1. It opens no "
                      "release report and no score", method)
        self.assertIn("Added after the one reading.", self.sections["No call by stratum"])
        self.assertIn(STRATA_COMMAND, (RECORD / "README.md").read_text(encoding="utf-8"))
        for name in SECTIONS:
            for sentence in re.split(r"(?<=[.:;])\s", self.sections[name]):
                if re.search(r"\bmiss(es|ed)?\b", sentence):
                    self.assertRegex(sentence, r"reader|program|known-gene rule", sentence)

    def test_the_four_no_call_counts_sum_to_the_readings_no_call(self):
        table = {r["stratum"]: r for r in rows(RECORD / ym.NO_CALL_TABLE)}
        for column in ("resistant", "called", "no_call"):
            self.assertEqual(sum(int(table[s][column]) for s in FOUR), int(table["all_strata"][column]), column)
            self.assertEqual(sum(int(table[s][column]) for s in FOUR[:3]), int(table["program_row"][column]), column)
        self.assertEqual(int(table["all_strata"]["no_call"]), int(self.totals["reader_no_call"]["resistant"]))
        self.assertEqual(int(table["all_strata"]["called"]), int(self.totals["reader_call"]["resistant"]))
        self.assertEqual(table["all_strata"]["resistant"], self.step1["all_strata"]["published_resistant"])
        used = re.search(r"guard 5: (\d+) used", guards()).group(1)
        self.assertIn(f"over its {used} datasets, one row each, per M1 stratum", self.sections["No call by stratum"])
        self.assertIn(f"Its sums equal the reading's own: {count(table['all_strata']['resistant'])} resistant, "
                      f"{count(table['all_strata']['called'])} called and {count(table['all_strata']['no_call'])} "
                      "with no call", self.sections["No call by stratum"])

    def test_step_1_stands_as_m1_published_it(self):
        m1 = m1_rows()
        names = sorted({r["row"] for r in m1})
        for stratum in FOUR:
            these = {(r["resistant"], r["susceptible"]) for r in m1 if r["stratum"] == stratum}
            self.assertEqual(these, {(self.step1[stratum]["published_resistant"],
                                      self.step1[stratum]["published_susceptible"])}, stratum)
        self.assertIn(f"step 1: {len(names)} published rows the same;", guards())
        self.assertIn(f"The file's {len(names)} rows of each stratum give the same counts.", self.sections["Step 1"])
        everyone = self.step1["all_strata"]
        for side in ("resistant", "susceptible"):
            self.assertEqual(sum(int(self.step1[s][f"published_{side}"]) for s in FOUR),
                             int(everyone[f"published_{side}"]))
            self.assertEqual(everyone[f"published_{side}"], everyone[f"rebuilt_{side}"])

    def test_the_readers_call_equals_m1s_catalogue_row(self):
        catalog = {r["stratum"]: r for r in m1_rows() if r["row"] == "catalog"}
        call = self.totals["reader_call"]
        self.assertEqual(sum(int(catalog[s]["resistant_called"]) for s in FOUR), int(call["resistant"]))
        self.assertEqual(sum(int(catalog[s]["susceptible_called"]) for s in FOUR), int(call["susceptible"]))
        words = (", ".join(count(catalog[s]["resistant_called"]) for s in FOUR[:3])
                 + f" and {count(catalog['neither']['resistant_called'])} resistant; "
                 + ", ".join(count(catalog[s]["susceptible_called"]) for s in FOUR[:3])
                 + f" and {count(catalog['neither']['susceptible_called'])} susceptible")
        self.assertIn(f"the {count(call['resistant'])} resistant and {count(call['susceptible'])} susceptible genomes "
                      f"of `reader_call` equal", self.sections["Step 1"])
        self.assertIn(f"({words})", self.sections["Step 1"])

    def test_the_datasets_are_the_primary_set_named_in_full(self):
        primary = sorted((r["species"], r["antibiotic"]) for r in rows(LISTS) if (r["in_A"], r["in_M"]) == ("1", "0"))
        species = sorted({(r["species"], r["antibiotic"]) for r in rows(RECORD / "species.tsv")})
        self.assertEqual(species, primary)
        self.assertIn(f"guard 5: {len(primary)} used + ", guards())
        self.assertIn(f"The {len(primary)} datasets of the primary set", self.sections["Datasets"])
        used, in_m, not_a, not_list, no_split = re.search(
            r"guard 5: (\d+) used \+ (\d+) in_M \+ (\d+) not_in_A \+ (\d+) not_in_list \+ (\d+) no_lineage_split",
            guards()).groups()
        self.assertEqual(no_split, "0")
        self.assertIn(f"Not used: {in_m} reports in set M, {not_a} outside A and {not_list} not in",
                      self.sections["Datasets"])

    def test_the_guards_stand_as_guards_txt_gives(self):
        text = guards()
        self.assertIn("```\n" + text + "```\n", self.sections["Method"])
        self.assertEqual([line.split(":")[0] for line in text.splitlines()],
                         ["guard 1", "guard 2", "guard 3", "guard 4", "guard 5", "step 1", "guard 6 (tables)"])
        self.assertRegex(text, r"guard 6 \(tables\): 6 files, 0 with a genome id, 0 with a catalogue group name, "
                               r"against [1-9]\d* ids and [1-9]\d* group names")

    def test_the_reading_commands_and_checksums_stand_as_reading_txt_gives(self):
        reading = (RECORD / "reading.txt").read_text(encoding="utf-8").splitlines()
        commands = [line for line in reading if line.startswith("python3 -I ")]
        self.assertEqual(len(commands), 3)
        self.assertIn("```\n" + "\n".join(commands) + "\n```", self.sections["Method"])
        listed = {line.split("  ", 1)[1]: line.split("  ", 1)[0] for line in reading if "  " in line
                  and not line.startswith(("#", "python3"))}
        self.assertEqual(listed["SHA256SUMS"], hashlib.sha256((RECORD / "SHA256SUMS").read_bytes()).hexdigest())
        for name in ("scripts/c3-inputs/lists.tsv", "reference/mechanism-split/page-m1-strata.tsv"):
            self.assertEqual(listed[name], hashlib.sha256((common.ROOT / name).read_bytes()).hexdigest(), name)

    def test_the_record_holds_the_files_its_checksums_give(self):
        listed = {}
        for line in (RECORD / "SHA256SUMS").read_text(encoding="utf-8").splitlines():
            digest, name = line.split("  ", 1)
            listed[name] = digest
        self.assertEqual(sorted(listed), sorted(ym.TABLES))
        found = sorted(p.relative_to(RECORD).as_posix() for p in RECORD.rglob("*") if p.is_file())
        self.assertEqual(found, sorted(list(ym.TABLES) + ["README.md", "SHA256SUMS", "reading.txt", ym.NO_CALL_TABLE]))
        for name, digest in listed.items():
            self.assertEqual(hashlib.sha256((RECORD / name).read_bytes()).hexdigest(), digest, name)

    def test_the_fault_table_names_each_planted_faults_test(self):
        for fault, words, test, _ in year_faults.FAULTS:
            line = f"| {fault}. {words[0].upper()}{words[1:]} | `{test.split('test_year_mix.')[1]}` |"
            self.assertEqual(self.lines.count(line), 1, fault)
        self.assertIn(f"All {len(year_faults.FAULTS)} were caught.", self.sections["Method"])

    def test_the_page_gives_the_suites_own_test_counts(self):
        loader = unittest.TestLoader()
        here = loader.loadTestsFromName("tests.mechanism.test_year_page").countTestCases()
        script = loader.loadTestsFromName("tests.mechanism.test_year_mix").countTestCases()
        method = self.sections["Method"]
        self.assertIn(f"`tests/mechanism/test_year_mix.py` holds {script} SYNTHETIC tests", method)
        self.assertIn(f"`tests/mechanism/test_year_page.py` holds {here} tests of this page", method)

    def test_every_figure_of_the_prose_stands_in_the_record(self):
        # The rules' constants (0.5, 10, 100, the cap of 2,000 a species, 3 decimals), the release's revision, the
        # dates of the read, of the table added after it and of the maintainer's ruling, and M2's periods; every
        # other figure must stand in a printed table, in guards.txt, in M1's catalogue row or count the record's
        # years, rows or tests.
        fixed = {"0.5", "1", "3", "10", "100", "2000", "14", "23", "24", "25", "2026", "2010", "2014", "2015", "2019"}
        s = self.summary
        known = set(fixed) | {str(len(s["kept_years"])), str(len(s["thin_years"])), str(s["unknown_genomes"]),
                              str(s["genomes"]),
                              str(s["species_pairs"]), str(s["species_names"]), str(len(year_faults.FAULTS))}
        known |= {str(n) for n in s["species_rows"].values()}
        known |= {str(g["years_with_recall"]) for g in s["groups"].values()}
        loader = unittest.TestLoader()
        known |= {str(loader.loadTestsFromName(f"tests.mechanism.{name}").countTestCases())
                  for name in ("test_year_page", "test_year_mix")}
        known |= set(re.findall(r"\d+", guards()))
        for r in m1_rows():
            if r["row"] == "catalog":
                known |= {r["resistant_called"], r["susceptible_called"]}
        for lines in ym.page_tables(str(RECORD)).values():
            for line in lines:
                known |= {n.replace(",", "") for n in re.findall(r"\d[\d,]*(?:\.\d+)?", line)}
        head = self.text.split("\n## ", 1)[0]
        for name, section in [("head", head)] + [(n, self.sections[n]) for n in SECTIONS if n != "Plan"]:
            words = re.sub(r"\b[MY]\d+\b|\bsha256\b", "", prose(section))
            for n in re.findall(r"\d[\d,]*(?:\.\d+)?", words):
                self.assertIn(n.rstrip(",").replace(",", ""), known, f"{name}: {n}")

    def test_the_words_keep_the_limits(self):
        for sentence in re.split(r"(?<=[.:;])\s", self.text):
            for word in ("external validation", "prevalence", "trend", "slope", "loss"):
                if word in sentence:
                    self.assertRegex(sentence, r"\b(no|never|not)\b", sentence)
            if re.search(r"\b(grows?|shrinks?|rising|falling|increas\w*|decreas\w*)\b", sentence):
                self.assertRegex(sentence, r"\b(no|never|not)\b", sentence)
        # The plan's own words rule names both; no other sentence gives a difference as equal.
        for name in SECTIONS[:1] + SECTIONS[2:]:
            for sentence in re.split(r"(?<=[.:;])\s", self.sections[name]):
                if "difference" in sentence:
                    self.assertNotIn("equal", sentence)
        self.assertNotIn("—", self.text)
        self.assertNotIn("<!--", self.text)
        self.assertNotRegex(prose(self.text), r"\bI\b|\b[Mm]y\b|\b[Ww]e\b")
        self.assertNotRegex(self.text, r"\b(GCA|GCF)_\d|\bPRJ[A-Z]{2}\d|\b[cnpg]\.[A-Z-]?\d")


if __name__ == "__main__":
    unittest.main()
