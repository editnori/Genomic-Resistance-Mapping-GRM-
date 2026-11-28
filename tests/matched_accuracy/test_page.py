"""docs/matched-accuracy.md and the words of the files of both stages: "Not run yet." of stage 2's real reading,
bare headers, no result, stage 1's run and counts read from the committed copy of its export, stage 2's steps and
ends word for word, the words the page never uses, its row in docs/INDEX.md and the workflow's comment."""

import csv
import io
import re
import unittest

from tests.matched_accuracy import common, stage2_common
from tests.matched_accuracy.test_reference import AMENDMENT_SHA256, PROTOCOL_SHA256

# Words the page never uses, and words sent back by the reader it is written for. Written so that this file
# holds none of them.
BANNED = (r"\b(validated|generalis|generaliz|robust|reliab|real-world|real world|provenance|receipts?|manifest|"
          r"fold[e]rs?|safe|help|explain|becau[s]e|cau[s]e|due to|clinic|patient|sc[o]pe|threshold|mechanis|confer|"
          r"proves?|pro[o]f|guarantee|h[e]ld|spot[- ]checked|la[n]es?|census|landed|lands|cleanly|traps?|snags?|"
          r"breakthroughs?|synced|frames?|keepers?|guarded|own[e]rs?)\b")
ROLE = r"(?i)\b(mach[i]nes?|b[o]x(es)?|contr[o]ller|threa[d]s?|work[e]rs?|assist[a]nt|mod[e]ls?|vend[o]rs?|agen[t]s?)\b"
MAKER = r"(?i)\b(cl[a]ude|anthr[o]pic|co-author[e]d-by|generated[ ]with|open[a]i|chatg[p]t|g[p]t-?\d|gem[i]ni)\b"
MARKER = r"(?i)\bk-?mers?\b|\b(genes?|alleles?|mutations?|variants?|markers?|determinants?|snps?)\b"
# The numbers the page may give: the revision, the checks' figures, the steps and the stages; stage 1's run, its
# counts and coverage and the BioProject facts, each read from the committed export below; the protocol's constants
# and the reports' count. None is a result of stage 2.
NUMBERS = {"14", "1", "2", "3", "4", "5", "6", "7", "19:10:00", "2026", "1,285,111", "252,148",
           "27", "00755885", "0", "49,621", "48,156", "9,132", "9,122", "2,783", "2,749", "283,118", "272,840", "206",
           "205", "16", "215,319", "23", "6,308", "2771", "48384", "8", "10", "10,000", "42", "274"}
# The protocol's two percentiles, the one place the page writes a decimal point.
PERCENTILES = "the linear 2.5th and 97.5th percentiles"
# The lead's sentence of 27 Nov on a row under 10 pairs, with its two conditions; the review's check 16 cut its
# first 7 words, which the list item above it says.
UNDER_10 = ("A row under 10 pairs carries no reading, yet its value still follows by arithmetic, "
            "to within the rounding, from the mean above it and the other rows, when that mean prints and no other row "
            "beneath it lacks an estimate.")
# What a copy group is, in the words of its rule. The one place the page may name a k-mer.
COPY_GROUP_WORDS = "genomes of one dataset with identical canonical k-mer sets"
# Stage 2's part, line by line, as the code has it: the stops before step 1, the 6 steps, the stops on format, the
# line itself and every end after it. Stage 2's second review gave 6 steps that a cut left unpinned, and 4 stops, the
# top-level catalogue read and the line file's stop the page left out (its findings 11 and 12). Its findings 3 to 5,
# changed after the one reading, took out the two ends that were not allowed and the empty line file.
STAGE_TWO = (
    "## Stage two steps",
    "Before step 1, the run stops on an output path that is not a directory or already holds an export, a `local` "
    'directory or the line "scores opened: yes" of an earlier run, and on numpy other than the version the '
    "resampling is written for.",
    "1. Check the sha256 of every report, both row lists, stage 1's committed count and typing tables and the 5 "
    "November genome table against values no argument changes. Stop on a missing file or a report whose id is not "
    "its file name.",
    "2. Stop on a line repeated, out of stage 1's sort or giving a genome twice, a listed dataset or row no report "
    "holds, a fold or lineage cluster unlike its report's, or a report without scm, cart, rf, gbt and lasso in "
    "order.",
    "3. Rebuild each list's groups. Stop on a group with one lab outcome, or on counts or typing support that differ "
    "from stage 1's tables at any level.",
    "4. Stop on a report changed since step 1, a listed row moved, a catalogue call not true, false or null, fit "
    "sizes by fold not 5 finite numbers, or a dataset or species name with a full stop. Their loader keeps no value.",
    '5. Write "scores opened: yes" to a file and print it. Only then does a second loader read, for listed rows, the '
    "five scores, each learner's fit size by fold and the catalogue call, and whether the report's top-level "
    "catalogue is present.",
    "6. Measure both readings as below, scan the export for every genome id and accession read, after a positive "
    "control, and write it in integers, estimates in parts per million.",
    "A file not as stage 1 writes it stops its step too. A committed table must have stage 1's columns, whole lines "
    "and no line twice. A row list must be UTF-8, hold a line, end with a line break and give 8 fields a line, none "
    "blank but the year, with whole-number folds and lineage clusters. The reports directory must hold the reports "
    "the release's checksum list names, and each report and genome table row must pass stage 1's loader.",
    "A check that fails before the line prints one line, exits 1 and writes nothing. So does a line file that cannot "
    "be made or written. The run first removes a line file it made, so a rerun into that directory starts, and its "
    "stop line says when the file cannot be removed.",
    "After the line, a run exits 0 with the export, or stops with one line and no export on one of four ends:",
    "- A listed report that cannot be read, or whose bytes changed since step 1.",
    "- A score null, true or false, not a number or not finite.",
    "- A full stop in an export file, a failed positive control of the scan, or a scan hit, which prints one line a "
    "file.",
    "- An export write that fails, which first removes what it made.",
)


def prose(text):
    text = re.sub(r"```.*?```", " ", text, flags=re.S)
    return re.sub(r"`[^`]*`", "CODE", text).replace("*", "")


def sentences(text):
    out = []
    for line in prose(text).split("\n"):
        line = line.strip()
        if not line or line.startswith(("|", "#", "Record: ")):
            continue
        line = re.sub(r"^(- |\d+\. )", "", line)
        out += [s for s in re.split(r"(?<=[^A-Z][.!?])\s+", line) if s]
    return out


class Page(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.page = common.PAGE.read_text(encoding="utf-8")
        cls.lines = cls.page.rstrip("\n").split("\n")
        cls.body = "\n".join(cls.lines[:-1])

    def test_not_run_yet_speaks_of_stage_2s_real_reading_alone(self):
        self.assertEqual(self.page.count("Not run yet."), 1)
        self.assertEqual(self.lines[2], "Not run yet. Stage 2's real reading has not run and is not authorised: it "
                                        "needs the merged script and the lead's word after the merge.")

    def test_under_1685_words(self):
        # 500 until the fixes of the review added the two counts of two kinds, the BioProject shape and the 5
        # November table; 600 until stage 2's part added its steps, the lead's readings, ruling R59's conditions,
        # stage 1's run with R59's counts, the coverage and the two BioProject facts; 1,350 until ruling R62 cut 7
        # passages that repeat another line, under the 1,345 words of stage 2's first page; 1,345 until ruling R65
        # had the page name every stop before the line and every end after it.
        self.assertLess(len(self.page.split()), 1685)

    def test_the_record_line_gives_the_protocols_and_the_amendments_sha256(self):
        self.assertTrue(self.lines[-1].startswith("Record: "))
        for sha in (PROTOCOL_SHA256, AMENDMENT_SHA256):
            self.assertIn(f"sha256 {sha}", self.lines[-1])
            self.assertEqual(self.page.count(sha), 1)

    def test_paths_and_hashes_only_in_the_record_line(self):
        # A slash between digits is a fraction of weight, not a path.
        for line in self.lines[:-1]:
            self.assertIsNone(re.search(r"\.(py|md|tsv|json|jsonl|txt|gz)\b|\b[0-9a-f]{16,}\b|(?<!\d)/|/(?!\d)", line),
                              line)

    def test_headers_are_bare_nouns(self):
        headers = re.findall(r"^#+ (.+)$", self.page, re.M)
        self.assertEqual(headers[0], "Matched accuracy")
        self.assertLessEqual(len(headers), 12)
        for h in headers:
            self.assertIsNone(re.search(r"[:?!.,;()0-9]", h), h)
            self.assertIsNone(re.match(r"(?i)(the|a|an|each|every|how|why|what) ", h), h)
            self.assertLessEqual(len(h.split()), 4, h)

    def test_lengths(self):
        items = 0
        for block in re.split(r"\n\s*\n", self.body):
            for item in re.split(r"\n(?=- |\d+\. )", block):
                if re.match(r"(- |\d+\. )", item):
                    items += 1
                    self.assertLess(len(item.split()), 45, item)
                elif not item.startswith(("#", "|")):
                    self.assertLess(len(item.split()), 85, item)
        # 30 until ruling R65 had the page list the four allowed ends after the line.
        self.assertLessEqual(items, 34)

    def test_no_result(self):
        self.assertEqual(self.body.count(PERCENTILES), 1)
        body = self.body.replace(PERCENTILES, "")
        found = set(re.findall(r"(?<![A-Za-z\d])(?:\d[\d,:]*\d|\d)(?![A-Za-z\d])", body))
        self.assertLessEqual(found, NUMBERS)
        self.assertIsNone(re.search(r"\d\.\d|%", body))
        for s in sentences(self.page):
            self.assertIsNone(re.search(r"(?i)\b(was|were|found|came out|turned out|shows?|showed)\b", s), s)

    def test_stage_2s_real_reading_has_not_run_and_waits_for_the_leads_word(self):
        for s in sentences(self.page):
            if re.search(r"(?i)authoris", s):
                self.assertIn("is not authorised", s)
        said = [s for s in sentences(self.page) if "real reading" in s]
        self.assertEqual(said, ["Stage 2's real reading has not run and is not authorised: it needs the merged script "
                                "and the lead's word after the merge.",
                                "The real reading runs once at the merged commit, after the lead's word, on the "
                                "computer holding the 274 reports and both row lists."])

    def test_a_retrospective_read_inside_one_database_and_never_external_validation(self):
        self.assertIn("a retrospective read inside one database", self.page)
        self.assertIsNone(re.search(r"(?i)external\s+validation", self.page))
        self.assertIsNone(re.search(r"(?i)prospective", self.page))

    def test_words_it_never_uses(self):
        page = self.page.replace(COPY_GROUP_WORDS, "")
        for pattern in (BANNED, ROLE, MAKER, MARKER):
            self.assertIsNone(re.search(pattern, page, re.I), pattern)
        self.assertNotIn(chr(0x2014), self.page)
        self.assertIsNone(re.search(r"(?i)\brather than\b", self.page))

    def test_a_copy_group_is_a_listed_group_of_identical_canonical_k_mer_sets(self):
        self.assertEqual(self.page.count(COPY_GROUP_WORDS), 1)
        self.assertIn("A copy group is a group the report lists: " + COPY_GROUP_WORDS + ". The script compares no "
                      "assembly.", self.page)

    def test_a_stop_prints_one_line_and_a_scan_hit_one_line_a_file(self):
        self.assertIn("A failed check prints one line, and a scan hit one line a file. Both exit 1 and write "
                      "nothing, and a write that fails removes what it made.", self.page)
        self.assertNotIn("Each stop prints one line", self.page)

    def test_the_committed_files_the_tests_read_are_named_by_the_workflow_and_not_the_page(self):
        # Ruling R62 cut the page's line on them, which repeated the workflow's comment and docs/ci.md.
        text = common.WORKFLOW.read_text(encoding="utf-8")
        block = text.split("\n          - suite: tests/matched_accuracy\n", 1)[0].rsplit("\n          - suite:", 1)[1]
        for path in ("reference/matched-accuracy", "reference/lineage-control/rev14-274-20261112.reports.sha256",
                     "reference/matched-accuracy/stage1-rev14-274-20261112"):
            self.assertIn(path, block)
        self.assertEqual(len(list(stage2_common.STAGE1_COPY.iterdir())), 5)
        self.assertNotIn("committed files", self.page)
        self.assertNotIn("made-up files alone", self.page)

    def test_a_row_list_gives_the_counts_its_sha256_and_line_count(self):
        self.assertIn("Row lists stay local. Only their sha256 and line counts enter the counts.", self.page)
        self.assertNotIn("only its sha256 enters", self.page)

    def test_the_lab_call_is_checked_among_rows_not_labelled_from_an_mic(self):
        self.assertIn("a lab call other than the label table's among rows not labelled from an MIC.", self.page)

    def test_the_two_counts_of_rows_lacking_a_typing_value_and_of_genomes_missing(self):
        self.assertIn("Rows lacking a typing value are counted among rows entering the groups and among eligible "
                      "rows; genomes missing from the genome table, over every report and over those that enter.",
                      self.page)

    def test_a_bioproject_is_one_part_of_an_accessions_shape(self):
        self.assertIn("A BioProject is one part, PRJ, two letters and digits in any case. A genome whose one part has "
                      "another shape has none and is counted apart. Blank genome table rows are skipped.", self.page)

    def test_the_5_november_table_has_no_insertion_time_and_gives_the_hc20_counts(self):
        self.assertIn("The 5 November genome table gives the HC20 group and has no insertion time to check. Without "
                      "it, no HC20 count is made.", self.page)

    def test_a_semicolon_joins_no_two_full_sentences(self):
        # The two left join a clause and an elliptical one, or the items of a list.
        self.assertEqual([line[:30] for line in self.body.split("\n") if ";" in line],
                         ["Rows lacking a typing value ar", "- Secondary reads keep the wei"])

    def test_the_export_holds_no_row_per_genome_pair_group_or_bioproject(self):
        # Ruling R59 item 6's words.
        self.assertIn("- No output names or ranks a feature. The export holds no row per genome, pair, group or "
                      "BioProject.\n", self.page)
        self.assertNotIn("stands for one genome", self.page)

    def test_in_the_index(self):
        index = common.INDEX.read_text(encoding="utf-8")
        self.assertEqual(index.count("[matched-accuracy.md](matched-accuracy.md)"), 1)
        row = next(line for line in index.split("\n") if line.startswith("| [matched-accuracy.md]"))
        self.assertIn("stage 2's real reading is not authorised. Not run yet.", row)

    def test_the_workflows_comment_names_the_page(self):
        text = common.WORKFLOW.read_text(encoding="utf-8")
        block = text.split("\n          - suite: tests/matched_accuracy\n", 1)[0].rsplit("\n          - suite:", 1)[1]
        comments = [line.strip() for line in block.split("\n") if line.strip().startswith("#")]
        self.assertTrue(comments[0].startswith("# SYNTHETIC"), comments)
        self.assertIn("docs/matched-accuracy.md", comments[0])


def committed(name):
    """The rows of one file of stage 1's committed export."""
    text = (stage2_common.STAGE1_COPY / name).read_text(encoding="utf-8")
    return list(csv.DictReader(io.StringIO(text), delimiter="\t"))


def n(value):
    return f"{int(value):,}"


class StageOneRun(unittest.TestCase):
    """Stage 1's run, R59's counts, the coverage and the BioProject facts, each read from the committed export."""

    @classmethod
    def setUpClass(cls):
        cls.page = common.PAGE.read_text(encoding="utf-8")

    def test_stage_1_ran_once_at_main_00755885(self):
        self.assertIn('Stage 1 ran once on 27 November 2026 at main 00755885, ended with exit 0 and printed "scores '
                      'opened: no". Ruling R59 gives these counts.\n', self.page)

    def test_r59s_counts_are_the_committed_count_tables(self):
        overall = {r["reading"]: r for r in committed("counts.tsv") if r["level"] == "overall"}
        labels = (("Eligible rows", "rows"), ("Genomes", "genomes"), ("Matched groups", "groups"),
                  ("Resistant-susceptible pairs", "pairs"), ("Datasets", "datasets"), ("Species", "species"))
        table = [line for line in self.page.split("\n") if line.startswith("| ")]
        self.assertEqual(table[0], "| Count | Reading A | Reading B |")
        self.assertEqual(table[2:], [f"| {label} | {n(overall['A'][c])} | {n(overall['B'][c])} |"
                                     for label, c in labels])
        self.assertEqual(n(overall["A"]["rows"]), "49,621")

    def test_the_coverage_of_both_readings_is_the_committed_count_tables(self):
        rows = committed("counts.tsv")
        cover = {}
        for reading in ("A", "B"):
            entering = sum(int(r["genome_rows"]) for r in rows
                           if (r["level"], r["reading"], r["reports_left_out"]) == ("dataset", reading, "0"))
            species = [r for r in rows if (r["level"], r["reading"]) == ("species", reading)]
            entered = [r for r in species if int(r["rows"]) > 0]
            overall = next(r for r in rows if (r["level"], r["reading"]) == ("overall", reading))
            cover[reading] = (n(overall["rows"]), n(entering), len(entered), len(species))
        a, b = cover["A"], cover["B"]
        self.assertIn(f"Coverage, which every later figure carries: {a[0]} of {a[1]} genome rows and {a[2]} of {a[3]} "
                      f"species in reading A, and {b[0]} of {b[1]} rows and {b[2]} of {b[3]} species in reading B. Each "
                      f"row count of stage 2's export is within those {a[1]} rows, and its inputs table gives both "
                      f"wholes.", self.page)
        self.assertEqual(cover, {"A": ("49,621", "215,319", 16, 23), "B": ("48,156", "215,319", 16, 23)})

    def test_the_two_bioproject_facts_are_the_committed_shares(self):
        shares = {r["reading"]: r for r in committed("shares.tsv")}
        a = shares["A"]
        self.assertIn(f"- The BioProject with the most eligible rows: {n(a['largest_bioproject_rows'])} of "
                      f"{n(a['eligible_rows'])} rows in reading A, and {a['its_weight_numerator']}/"
                      f"{a['its_weight_denominator']} of the weight.", self.page)
        self.assertEqual({(r["largest_weight_numerator"], r["largest_weight_denominator"]) for r in shares.values()},
                         {("1", "8")})
        self.assertIn("- The largest weight of any BioProject: 1/8 in both readings.", self.page)
        self.assertEqual((a["largest_bioproject_rows"], a["its_weight_numerator"], a["its_weight_denominator"]),
                         ("6308", "2771", "48384"))


class LeadsReadings(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.page = common.PAGE.read_text(encoding="utf-8")

    def test_a_row_under_10_pairs_prints_no_estimate_and_its_value_follows_by_arithmetic(self):
        self.assertEqual(self.page.count(UNDER_10), 1)
        self.assertIn("- A dataset or species under 10 pairs enters every mean. Its row gives counts and no estimate.",
                      self.page)
        self.assertIn("A cell prints its estimate wherever it has pairs and gets no interval.", self.page)

    def test_each_cell_is_weighted_over_its_own_groups_datasets_and_species(self):
        self.assertIn("Each cell is weighted over its own groups, datasets and species. A cell prints", self.page)

    def test_the_tie_and_the_reading_of_an_interval(self):
        self.assertIn("a higher resistant mean counts 1, an exact tie one half, else 0", self.page)
        self.assertEqual(self.page.count("A claim above chance needs `above` in both readings."), 1)
        self.assertIn("`above` with the lower bound above one half, `below` with the upper bound below it, else "
                      "`unresolved`", self.page)
        self.assertIn("from its bounds before rounding", self.page)
        self.assertIn("added left to right as doubles", self.page)


class StageTwoSteps(unittest.TestCase):
    """Stage 2's part word for word, its names and messages read from the script, and ruling R59's one reading."""

    @classmethod
    def setUpClass(cls):
        cls.page = common.PAGE.read_text(encoding="utf-8")
        cls.script = stage2_common.SCRIPT.read_text(encoding="utf-8")

    def test_the_part_is_word_for_word_the_steps_and_ends_of_the_code_today(self):
        part = self.page[self.page.index("## Stage two steps\n"):self.page.index("## Lead's readings\n")]
        self.assertEqual(tuple(line for line in part.split("\n") if line.strip()), STAGE_TWO)

    def test_its_names_and_messages_are_the_scripts(self):
        stage2 = stage2_common.stage2
        self.assertEqual(stage2.OPENED, "scores opened: yes")
        self.assertEqual(stage2.FIVE, ("scm", "cart", "rf", "gbt", "lasso"))
        self.assertEqual((len(stage2.ROW_FIELDS), stage2.ROW_FIELDS[-1], stage2.stage1.FOLDS), (8, "year", 5))
        for message in ("a score was asked for before the line",
                        "the output directory already holds export, local or the line of an earlier run",
                        "numpy is not the version the resampling is written for"):
            self.assertEqual(self.script.count(f'raise Stop("{message}")'), 1, message)
        self.assertEqual(self.script.count('"catalog": "present" if doc.get("catalog") is not None'), 1)

    def test_a_run_that_printed_the_line_is_the_one_reading(self):
        # Ruling R59 item 4's rule.
        self.assertIn('- The real reading runs once at the merged commit, after the lead\'s word, on the computer '
                      'holding the 274 reports and both row lists. A run that stops before "scores opened: yes" may '
                      'run again. One that printed it is the one reading.\n', self.page)


class FileWords(unittest.TestCase):
    """The words of every file of both stages: no maker, no long dash, and no real id, accession or sequence."""

    @classmethod
    def setUpClass(cls):
        paths = [common.SCRIPT, common.FAULTS, stage2_common.SCRIPT, stage2_common.FAULTS, common.PAGE] + \
            sorted((common.ROOT / common.SUITE).glob("*.py"))
        cls.texts = {p.relative_to(common.ROOT).as_posix(): p.read_text(encoding="utf-8") for p in paths}

    def test_no_maker_and_no_long_dash(self):
        for name, text in self.texts.items():
            with self.subTest(file=name):
                self.assertIsNone(re.search(MAKER, text))
                self.assertNotIn(chr(0x2014), text)
                self.assertIsNone(re.search(r"(?i)\brather than\b", text))

    def test_no_id_accession_or_sequence(self):
        for name, text in self.texts.items():
            with self.subTest(file=name):
                self.assertIsNone(re.search(r"(?i)(?:PRJ[A-Z]{1,2}|SAM[A-Z]{1,2}|GC[AF]_|[SED]R[RSXP])\d+", text))
                self.assertIsNone(re.search(r"\b\d{3,7}\.\d{3,6}\b", text))
                self.assertIsNone(re.search(r"\b[ACGT]{8,}\b", text))

    def test_no_name_of_a_marker_in_the_script_and_its_runner(self):
        for name in ("scripts/matched-accuracy/stage1.py", "scripts/matched-accuracy/planted_faults.py",
                     "scripts/matched-accuracy/stage2.py", "scripts/matched-accuracy/stage2_faults.py"):
            self.assertIsNone(re.search(MARKER, self.texts[name]), name)


if __name__ == "__main__":
    unittest.main()
