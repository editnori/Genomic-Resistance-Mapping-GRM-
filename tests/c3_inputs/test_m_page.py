"""scripts/c3-inputs/m_page.py: the record and the page tables of the run on the datasets of M.

Every table here is SYNTHETIC: m_reading.py's own floor and count forms
write the 7 tables from test_m_reading.py's made-up species directories
(tests/c3_inputs/m_tables.py), 6 made-up datasets of M in two made-up
species, no genome. Only the column lists, words and number forms are read
from the committed m_reading.py and page_tables.py, and the count of the
datasets of M from the committed lists.tsv. A test of a rule names the
lines and words it expects itself.
"""

import contextlib
import io
import os
import shutil
import sys
import tempfile
import unittest
from pathlib import Path

from tests.c3_inputs import common, m_tables
from tests.c3_inputs import test_m_reading as mr

m_page = m_tables.m_page
reading = mr.reading
page_tables = common.load("page_tables")
m_list = common.load("m_list")
M_SET = m_tables.M_SET
SEVEN = ("m-floor.tsv", "m-left-out.tsv", "m-count.tsv", "m-conditions.tsv", "m-both.tsv", "m-size.tsv",
         "m-datasets.tsv")
LEFT_OUT = [
    "| Species | Antibiotic | Second reading | Cohort | Listed | Left after the list | Floor |",
    "|---|---|---|---|---|---|---|",
    "| *Synth c* | drug 5 | left out, short list | 40 R, 40 S | 1 R, 0 S | 39 R, 40 S | 25 of each class |",
    "| *Synth c* | drug 6 | left out, short list | 40 R, 40 S | 4 R, 0 S | 36 R, 40 S | 25 of each class |",
    "| *Synth d* | drug 10 | left out, under the floor | 30 R, 40 S | 6 R, 2 S | 24 R, 38 S | 25 of each class |"]
SPECIES = [
    "| Species | Datasets of M | In the second reading | Left out, under the floor | Left out, short list |",
    "|---|---|---|---|---|",
    "| *Synth c* | 3 | 1 of 3 | 0 of 3 | 2 of 3 |",
    "| *Synth d* | 3 | 2 of 3 | 1 of 3 | 0 of 3 |",
    "| All 2 species | 6 | 3 of 6 | 1 of 6 | 2 of 6 |"]
READS = ("Datasets of M read: 6 of 6 in the first reading, 3 of 6 in the second (1 left out under the floor, 2 for a "
         "short list).")
NOT_TAKEN = "The count has not been taken: the record holds the floor's two tables alone."
NONE_LEFT = "0 of the 6 datasets of M are left out of the second reading."
LEFT_AT = "The design page, `docs/input-arms.md`, names each dataset left out."
BUILT = {}
WATCH = {"on": None, "hooked": False}


def _audit(event, args):
    seen = WATCH["on"]
    if seen is not None and event in ("open", "os.listdir", "os.scandir") and args:
        if isinstance(args[0], (str, bytes, os.PathLike)):
            seen.append((event, os.fsdecode(args[0])))


@contextlib.contextmanager
def watched():
    """Each file this process opens and each directory it lists while the block runs, as (event, path). The audit
    hook that records them is added once and records nothing outside such a block."""
    if not WATCH["hooked"]:
        sys.addaudithook(_audit)
        WATCH["hooked"] = True
    WATCH["on"] = []
    try:
        yield WATCH["on"]
    finally:
        WATCH["on"] = None


def setUpModule():
    BUILT["tmp"] = tempfile.TemporaryDirectory()
    top = Path(BUILT["tmp"].name)
    BUILT["seven"] = m_tables.seven(top / "default")
    BUILT["all in"] = m_tables.seven(top / "all-in", m_tables.ALL_IN)


def tearDownModule():
    BUILT.pop("tmp").cleanup()


def lines_of(path):
    return Path(path).read_text(encoding="utf-8").split("\n")


def set_cell(directory, name, n, column, value, sums=True):
    """Line n of one table of <directory> with one cell set, by m_tables.edit; a tables directory takes sums=False."""
    def change(text):
        lines = text.split("\n")
        cells = lines[n - 1].split("\t")
        cells[lines[0].split("\t").index(column)] = value
        lines[n - 1] = "\t".join(cells)
        return "\n".join(lines)
    m_tables.edit(directory, name, change, sums)


def line_of(path, *words):
    """The 1-based line of <path> that holds every one of words, which must be one line."""
    found = [i for i, line in enumerate(lines_of(path), 1) if all(w in line.split("\t") for w in words)]
    assert len(found) == 1, (path, words, found)
    return found[0]


class Made(unittest.TestCase):
    """Each test in a directory of its own, with the SYNTHETIC tables, a record and a lists.tsv."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.top = Path(self.tmp.name)
        self.made_n = 0

    def tearDown(self):
        self.tmp.cleanup()

    def made(self, upto="count", which="seven"):
        """(tables dir, record dir), new on each call: the record as it stands before the run (upto none), after
        floor (upto floor) or after floor and record (upto count)."""
        self.made_n += 1
        base = self.top / f"made-{self.made_n}"
        tables = base / "tables"
        shutil.copytree(BUILT[which], tables)
        record = m_tables.record(base)
        if upto in ("floor", "count"):
            m_page.floor(tables, record, M_SET)
        if upto == "count":
            m_page.record_count(tables, record, M_SET)
        return tables, record

    def run_main(self, *argv):
        lists = m_tables.lists_file(self.top / "lists.tsv")
        with contextlib.redirect_stdout(io.StringIO()) as out:
            rc = m_page.main([str(a) for a in argv], lists=lists)
        return rc, out.getvalue()

    def stops(self, record, pattern):
        with self.assertRaisesRegex(m_page.Stop, pattern):
            m_page.tables(record, M_SET)


class Columns(unittest.TestCase):
    def test_the_column_lists_and_states_are_m_readings_own(self):
        self.assertEqual(m_page.COLUMNS, {
            "m-floor.tsv": reading.FLOOR_PAGE_COLUMNS, "m-left-out.tsv": reading.LEFT_OUT_COLUMNS,
            "m-count.tsv": reading.COUNT_COLUMNS, "m-conditions.tsv": reading.CONDITION_COLUMNS,
            "m-both.tsv": reading.BOTH_COLUMNS, "m-size.tsv": reading.SIZE_COLUMNS,
            "m-datasets.tsv": reading.M_DATASET_COLUMNS})
        self.assertEqual(sorted(m_page.TABLES), sorted(SEVEN))
        self.assertEqual(m_page.FLOOR_FILES, reading.FLOOR_FILES)
        self.assertEqual((m_page.STATES, m_page.READINGS, m_page.COUNTS, m_page.READ),
                         (reading.STATES, reading.READINGS, reading.COUNTS, reading.READ))
        self.assertEqual((m_page.DEFAULT, m_page.PRIMARY, m_page.DECIDING),
                         (page_tables.DEFAULT, page_tables.PRIMARY, page_tables.DECIDING))
        self.assertEqual((m_page.FLOOR, m_page.SMALLEST_LISTED), (m_list.FLOOR, reading.SMALLEST_LISTED))
        self.assertEqual(m_page.COLUMNS["m-floor.tsv"][:-1], m_list.FLOOR_COLUMNS)
        for name, cells in list(m_page.COUNT_CELLS.items()) + list(m_page.MET_CELLS.items()):
            self.assertLess(set(cells), set(m_page.COLUMNS[name]), name)

    def test_the_words_and_number_forms_are_page_tables_own(self):
        self.assertEqual((m_page.ARM_WORDS, m_page.MEASURE_WORDS), (page_tables.ARM_WORDS, page_tables.MEASURE_WORDS))
        for x in ("", None, "0", "0.0005", "-0.0005", "0.001", "0.0015", "0.0025", "0.123456789012", "-0.02", "1",
                  "0.9995", "1e-05", 0.0004, 0.001, -0.1234):
            self.assertEqual(m_page.d3(x), page_tables.d3(x), x)
            self.assertEqual(m_page.pval(x), page_tables.pval(x), x)
            self.assertEqual(m_page.interval(x, "0.5"), page_tables.interval(x, "0.5"), x)
        self.assertEqual((m_page.pval("0.001"), m_page.pval("0.000999")), ("0.001", "< 0.001"))
        for a, b in (("0", "0"), ("7", "232"), ("1285", "1285111")):
            self.assertEqual(m_page.of(a, b), page_tables.of(a, b))
        for x in ("1", "0", ""):
            self.assertEqual(m_page.yes(x), page_tables.yes(x))
        self.assertEqual(m_page.table(["A", "B"], [[1, "x"]]), page_tables.table(["A", "B"], [[1, "x"]]))

    def test_the_committed_lists_tsv_marks_35_datasets_of_m(self):
        self.assertEqual(len(m_page.m_datasets(common.KIT / "lists.tsv")), 35)


class Floor(Made):
    def test_floor_copies_the_two_tables_and_sums_every_file_but_itself_and_the_readme(self):
        tables, record = self.made("none")
        rc, out = self.run_main("floor", tables, record)
        self.assertEqual((rc, out), (0, "floor: 2 tables copied; SHA256SUMS lists 3 files\n"))
        self.assertEqual(sorted(p.name for p in record.iterdir()),
                         ["README.md", "SHA256SUMS", "first-run-times.tsv", "m-floor.tsv", "m-left-out.tsv"])
        sums = [line.split("  ") for line in (record / "SHA256SUMS").read_text().splitlines()]
        self.assertEqual([n for _, n in sums], ["first-run-times.tsv", "m-floor.tsv", "m-left-out.tsv"])
        for digest, name in sums:
            self.assertEqual(digest, m_page.sha256(record / name), name)
        for name in ("m-floor.tsv", "m-left-out.tsv"):
            self.assertEqual((record / name).read_bytes(), (tables / name).read_bytes(), name)

    def test_floor_reads_no_other_table_and_does_not_list_the_tables_directory(self):
        # All 7 tables stand in the tables directory as files. An audit hook records each open and each listing while
        # floor runs: under that directory it opens its two tables alone and lists nothing.
        tables, record = self.made("none")
        self.assertEqual(sorted(p.name for p in tables.iterdir()), sorted(SEVEN))
        self.assertTrue(all((tables / name).is_file() for name in SEVEN))
        with watched() as seen:
            rc, _ = self.run_main("floor", tables, record)
        self.assertEqual(rc, 0)
        under = sorted((event, Path(p).name) for event, p in seen if Path(p) == tables or tables in Path(p).parents)
        self.assertEqual(under, [("open", "m-floor.tsv"), ("open", "m-left-out.tsv")])
        self.assertTrue(any(Path(p) == record / "m-floor.tsv" for _, p in seen))
        self.assertEqual(m_page.state(record), "floor")

    def test_floor_runs_once_on_a_record_that_holds_the_readme_and_the_times_alone(self):
        tables, record = self.made("floor")
        before = {p.name: p.read_bytes() for p in record.iterdir()}
        rc, out = self.run_main("floor", tables, record)
        self.assertEqual(rc, 1)
        self.assertIn("m_page: stopped: the record holds SHA256SUMS: the floor's tables are copied once", out)
        self.assertEqual({p.name: p.read_bytes() for p in record.iterdir()}, before)
        other = m_tables.record(self.top / "other")
        (other / "notes.txt").write_text("x\n")
        with self.assertRaisesRegex(m_page.Stop, "the record holds notes.txt"):
            m_page.floor(tables, other, M_SET)

    def test_floor_stops_on_floor_tables_that_break_a_gate_and_writes_nothing(self):
        tables, record = self.made("none")
        with self.assertRaisesRegex(m_page.Stop, "m-floor.tsv lacks synth_x__drug_1, a dataset lists.tsv marks in_M"):
            m_page.floor(tables, record, M_SET | {"synth_x__drug_1"})
        self.assertEqual(sorted(p.name for p in record.iterdir()), ["README.md", "first-run-times.tsv"])


class Record(Made):
    def test_record_copies_the_other_five_and_sums_every_file_again(self):
        tables, record = self.made("floor")
        rc, out = self.run_main("record", tables, record)
        self.assertEqual((rc, out), (0, "record: 5 tables copied; SHA256SUMS lists 8 files\n"))
        self.assertEqual(sorted(p.name for p in record.iterdir()),
                         sorted(("README.md", "SHA256SUMS", "first-run-times.tsv") + SEVEN))
        sums = [line.split("  ") for line in (record / "SHA256SUMS").read_text().splitlines()]
        self.assertEqual([n for _, n in sums], sorted(("first-run-times.tsv",) + SEVEN))
        for digest, name in sums:
            self.assertEqual(digest, m_page.sha256(record / name), name)
        for name in SEVEN:
            self.assertEqual((record / name).read_bytes(), (tables / name).read_bytes(), name)

    def test_record_stops_unless_the_tables_directory_holds_exactly_the_seven(self):
        tables, record = self.made("floor")
        (tables / "notes.txt").write_text("x\n")
        with self.assertRaisesRegex(m_page.Stop, "notes.txt in .* is not one of m_reading.py's 7 tables"):
            m_page.record_count(tables, record, M_SET)
        (tables / "notes.txt").unlink()
        kept = (tables / "m-size.tsv").read_bytes()
        (tables / "m-size.tsv").unlink()
        with self.assertRaisesRegex(m_page.Stop, "m-size.tsv is missing from"):
            m_page.record_count(tables, record, M_SET)
        (tables / "m-size.tsv").write_bytes(kept)
        self.assertEqual(m_page.record_count(tables, record, M_SET),
                         "record: 5 tables copied; SHA256SUMS lists 8 files")

    def test_record_stops_on_floor_tables_unlike_the_records(self):
        # One susceptible genome more in synth_c__drug_6's cohort, in both floor tables and every count it moves, so
        # every gate passes on the new tables and only the bytes differ.
        tables, record = self.made("floor")
        old, new = "\t80\t40\t40\t4\t4\t0\t36\t40\t25", "\t81\t40\t41\t4\t4\t0\t36\t41\t25"
        for name in ("m-floor.tsv", "m-left-out.tsv"):
            text = (tables / name).read_text()
            self.assertEqual(text.count(old), 1, name)
            (tables / name).write_text(text.replace(old, new))
        m_page.gate_tables(m_page.read_tables(tables, SEVEN), M_SET)
        with self.assertRaisesRegex(m_page.Stop, "^m-floor.tsv in .* differs from the record's$"):
            m_page.record_count(tables, record, M_SET)
        self.assertEqual(m_page.state(record), "floor")

    def test_record_gates_the_seven_before_it_copies_any(self):
        # A second-reading row for a dataset left out for a short list: the record keeps the floor's two alone.
        tables, record = self.made("floor")
        path = tables / "m-datasets.tsv"
        lines = lines_of(path)
        row = next(line for line in lines if line.startswith("first\tsynth_c__drug_5\t"))
        path.write_text("\n".join(lines[:-1] + ["second" + row[len("first"):], ""]))
        before = {p.name: p.read_bytes() for p in record.iterdir()}
        with self.assertRaisesRegex(m_page.Stop, "a second-reading row of synth_c__drug_5"):
            m_page.record_count(tables, record, M_SET)
        self.assertEqual({p.name: p.read_bytes() for p in record.iterdir()}, before)

    def test_record_runs_after_the_floor_and_once(self):
        tables, record = self.made("none")
        with self.assertRaisesRegex(m_page.Stop, "floor runs first, and record once"):
            m_page.record_count(tables, record, M_SET)
        m_page.floor(tables, record, M_SET)
        m_page.record_count(tables, record, M_SET)
        with self.assertRaisesRegex(m_page.Stop, "floor runs first, and record once"):
            m_page.record_count(tables, record, M_SET)


class Print(Made):
    def test_print_gives_the_left_out_and_species_tables_as_written_here(self):
        _, record = self.made()
        printed = m_page.tables(record, M_SET)
        self.assertEqual(list(printed), ["left-out", "species", "count", "beside", "conditions", "both", "size",
                                         "datasets"])
        self.assertEqual((printed["left-out"], printed["species"]), (LEFT_OUT, SPECIES))

    def test_the_left_out_table_goes_by_species_then_antibiotic_and_not_by_dataset_id(self):
        # synth_d__drug_10 put in a made-up species, Synth a, that sorts before Synth c while its id sorts after
        # synth_c's; and synth_c__drug_5's antibiotic named drug 9, which sorts after drug 6 while its id sorts before
        # synth_c__drug_6's.
        tables, record = self.made("none")
        for name in ("m-floor.tsv", "m-left-out.tsv"):
            text = (tables / name).read_text()
            for old, new in (("synth_d__drug_10\tsynth_d\tSynth d\t", "synth_d__drug_10\tsynth_d\tSynth a\t"),
                             ("synth_d__drug_10\tSynth d\t", "synth_d__drug_10\tSynth a\t"),
                             ("\tdrug 5\t", "\tdrug 9\t")):
                text = text.replace(old, new)
            (tables / name).write_text(text)
        m_page.floor(tables, record, M_SET)
        left = m_page.tables(record, M_SET)["left-out"]
        self.assertEqual(left, LEFT_OUT[:2] + [LEFT_OUT[4].replace("*Synth d*", "*Synth a*"), LEFT_OUT[3],
                                               LEFT_OUT[2].replace("drug 5", "drug 9")])
        page = self.top / "sorted.md"
        page.write_text("\n".join(left[:3] + left[4:]) + "\n")
        self.assertEqual(self.run_main("onpage", record, page),
                         (1, f"onpage: {page} lacks the left-out table's row of Synth c with drug 6\n"))

    def test_print_with_the_floor_tables_alone_says_the_count_has_not_been_taken(self):
        _, record = self.made("floor")
        self.assertEqual(m_page.tables(record, M_SET), {"left-out": LEFT_OUT, "species": SPECIES,
                                                        "count": [NOT_TAKEN]})
        rc, out = self.run_main("print", record)
        self.assertEqual(rc, 0)
        self.assertEqual(out, "[left-out]\n" + "\n".join(LEFT_OUT) + "\n\n[species]\n" + "\n".join(SPECIES)
                         + "\n\n[count]\n" + NOT_TAKEN + "\n\n")

    def test_page_gives_the_result_pages_six_tables_from_the_record(self):
        tables, record = self.made()
        printed, blocks = m_page.tables(record, M_SET), m_page.page(record, M_SET)
        self.assertEqual(tuple(blocks), ("reads", "both", "species", "count-all", "recall-all", "study-all", "count-m",
                                         "files"))
        self.assertEqual(m_page.ON_PAGE, tuple(blocks))
        self.assertEqual(blocks["reads"], [READS + " " + LEFT_AT])
        self.assertEqual(blocks["both"], [
            "**Both readings.** Per arm, conditions 1 to 5 of ruling 15 on the count over all datasets.", ""]
            + printed["both"][2:])
        self.assertEqual(len(blocks["both"]), 2 + 2 + 20)
        self.assertEqual(blocks["species"], [
            "**Datasets read.** The datasets of M by species: in the second reading, left out under the floor, or left "
            "out for a short list.", ""] + SPECIES)
        rows = mr.score.read_tsv(tables / "m-count.tsv")
        pd = page_tables
        for name, count, line in (
                ("count-all", "all datasets",
                 "**All datasets.** The first run's datasets and M's, on the lineage split and balanced accuracy."),
                ("count-m", "M", "**M alone.** The datasets of M, on the lineage split and balanced accuracy. "
                                 "This count decides nothing.")):
            want = page_tables.table(
                ["Reading", "Arm", "Learner", "Better", "Worse", "Same", "Median", "Interval", "Sign p", "Holm p",
                 "Constant"],
                [[r["m_reading"], pd.pair_words(r["pair"]), f"`{r['learner']}`", pd.of(r["better"], r["entered"]),
                  pd.of(r["worse"], r["entered"]), pd.of(r["same"], r["entered"]), pd.d3(r["median_difference"]),
                  pd.interval(r["interval_lo"], r["interval_hi"]), pd.pval(r["sign_p"]), pd.pval(r["holm_p"]),
                  r["constant"]] for r in rows
                 if (r["count"], r["split"], r["measure"]) == (count, "lineage", "balanced_accuracy")])
            # 6 arm rows by 5 learners in each of the 2 readings, under the header and rule rows.
            self.assertEqual(len(want), 2 + 60)
            self.assertEqual(blocks[name], [line, ""] + want)
        for name, split, measure, line, n in (
                ("recall-all", "lineage", "recall",
                 "**All datasets, recall.** The first run's datasets and M's, on the lineage split and recall.", 60),
                ("study-all", "study", "balanced_accuracy",
                 "**All datasets, study split.** The first run's datasets and M's, on the study split and balanced "
                 "accuracy.", 50)):
            # The lineage split's 6 arm rows by 5 learners, or the study split's 5 by 5, in each reading.
            want = page_tables.table(
                ["Reading", "Arm", "Learner", "Better", "Worse", "Same", "Median", "Interval"],
                [[r["m_reading"], pd.pair_words(r["pair"]), f"`{r['learner']}`", pd.of(r["better"], r["entered"]),
                  pd.of(r["worse"], r["entered"]), pd.of(r["same"], r["entered"]), pd.d3(r["median_difference"]),
                  pd.interval(r["interval_lo"], r["interval_hi"])] for r in rows
                 if (r["count"], r["split"], r["measure"]) == ("all datasets", split, measure)])
            self.assertEqual(len(want), 2 + n)
            self.assertEqual(blocks[name], [line, ""] + want)
        self.assertNotEqual(blocks["count-all"][2:], blocks["count-m"][2:])
        self.assertNotIn(m_page.table(m_page.BESIDE_HEADER, [])[0], {b[2] for b in blocks.values() if len(b) > 2})
        # 20 + 3 + 60 + 60 + 50 + 60 rows: the species table's 2 made-up species and its total row.
        self.assertEqual(blocks["files"], [
            "**Files.** The 6 tables above, 253 rows for the 2 species of M, come from the record, "
            "`reference/input-arms-m/`, which holds every other table: each count on every measure and split "
            "(`m-count.tsv`), each count's conditions (`m-conditions.tsv`), the size table (`m-size.tsv`) and each "
            "dataset's figures (`m-datasets.tsv`). `scripts/c3-inputs/m_page.py print` prints them."])
        rc, out = self.run_main("page", record)
        self.assertEqual((rc, out), (0, "\n\n".join("\n".join(lines) for lines in blocks.values()) + "\n"))
        # print still gives its 8 tables, and the record keeps its 5.
        self.assertEqual(list(printed), ["left-out", "species", "count", "beside", "conditions", "both", "size",
                                         "datasets"])

    def test_the_result_page_gives_the_reads_line_once_above_its_first_table(self):
        _, record = self.made()
        blocks = m_page.page(record, M_SET)
        lines = [line for b in blocks.values() for line in b]
        self.assertEqual([line for line in lines if line.startswith("Datasets of M read:")], [READS + " " + LEFT_AT])
        self.assertEqual(list(blocks)[0], "reads")
        self.assertEqual(sum(READS in line for line in lines), 1)

    def test_page_stops_before_the_count(self):
        _, record = self.made("floor")
        self.assertEqual(self.run_main("page", record), (
            1, "m_page: stopped: the count has not been taken, and the result page takes its tables after it\n"))
        _, record = self.made("none")
        self.assertEqual(self.run_main("page", record), (1, "m_page: stopped: m-floor.tsv is missing from the record\n"))

    def test_print_stops_on_a_record_with_no_table_of_m(self):
        _, record = self.made("none")
        rc, out = self.run_main("print", record)
        self.assertEqual((rc, out), (1, "m_page: stopped: m-floor.tsv is missing from the record\n"))

    def test_print_gives_the_count_rows_in_page_tables_number_forms(self):
        tables, record = self.made()
        printed = m_page.tables(record, M_SET)
        rows = mr.score.read_tsv(tables / "m-count.tsv")
        lin = [r for r in rows if r["split"] == "lineage" and r["measure"] == "balanced_accuracy"]
        rest = [r for r in rows if (r["split"], r["measure"]) != ("lineage", "balanced_accuracy")
                and not (r["split"] == "study" and r["measure"] != "balanced_accuracy")]
        # 6 arm rows by 5 learners for each of the 3 counts in each of the 2 readings; beside them, per reading, the
        # first run's and all datasets' 30 recall rows and 25 study rows, and M alone's 150 other rows and 25 study rows.
        self.assertEqual((len(lin), len(rest)), (180, 570))
        words = {"first run": "the first run as read", "M": "M alone", "all datasets": "all datasets"}
        pd = page_tables
        want = page_tables.table(
            ["Count", "Reading", "Arm", "Learner", "Better", "Worse", "Same", "Median", "Interval", "Sign p", "Holm p",
             "Constant"],
            [[words[r["count"]], r["m_reading"], pd.pair_words(r["pair"]), f"`{r['learner']}`",
              pd.of(r["better"], r["entered"]), pd.of(r["worse"], r["entered"]), pd.of(r["same"], r["entered"]),
              pd.d3(r["median_difference"]), pd.interval(r["interval_lo"], r["interval_hi"]), pd.pval(r["sign_p"]),
              pd.pval(r["holm_p"]), r["constant"]] for r in lin])
        self.assertEqual(printed["count"], [READS, ""] + want)
        beside = page_tables.table(
            ["Count", "Reading", "Split", "Arm", "Learner", "Measure", "Better is", "Entered", "Better", "Worse",
             "Same", "Median", "Interval"],
            [[words[r["count"]], r["m_reading"], r["split"], pd.pair_words(r["pair"]), f"`{r['learner']}`",
              pd.MEASURE_WORDS[r["measure"]], r["orientation"], pd.of(r["entered"], r["datasets"]), r["better"],
              r["worse"], r["same"], pd.d3(r["median_difference"]), pd.interval(r["interval_lo"], r["interval_hi"])]
             for r in rest])
        self.assertEqual(printed["beside"], [READS, ""] + beside)
        self.assertEqual([line.split(" | ")[0][2:] for line in printed["count"][4::30]],
                         ["the first run as read", "M alone", "all datasets"] * 2)

    def test_print_gives_conditions_1_to_5_by_ruling_15_for_each_count_and_reading(self):
        tables, record = self.made()
        printed = m_page.tables(record, M_SET)["conditions"]
        self.assertEqual(printed[:4], [READS, "", "| Count | Reading | Arm | Condition | Read | Learners |",
                                       "|---|---|---|---|---|---|"])
        rows = [r for r in mr.score.read_tsv(tables / "m-conditions.tsv")
                if r["reading"] == "ruling 15" and r["condition"] != "6"]
        # 3 counts by 2 readings by 5 arms by conditions 1 and 2, 3, 4 and 5.
        self.assertEqual(len(rows), 120)
        self.assertEqual({r["condition"] for r in rows}, {"1 and 2", "3", "4", "5"})
        self.assertEqual(len(printed), 4 + 120)
        words = {"first run": "the first run as read", "M": "M alone", "all datasets": "all datasets"}
        for line, r in zip(printed[4:], rows):
            self.assertEqual(line, f"| {words[r['count']]} | {r['m_reading']} | {page_tables.ARM_WORDS[r['arm']]} | "
                                   f"{r['condition']} | {page_tables.yes(r['met'])} | {r['learners']} |")

    def test_print_gives_the_both_table_and_the_size_table(self):
        tables, record = self.made()
        printed = m_page.tables(record, M_SET)
        both = mr.score.read_tsv(tables / "m-both.tsv")
        self.assertEqual(printed["both"], [READS, "", "| Arm | Condition | First reading | Second reading | Both |",
                                           "|---|---|---|---|---|"] + [
            f"| {page_tables.ARM_WORDS[r['arm']]} | {r['condition']} | {page_tables.yes(r['met_first'])} | "
            f"{page_tables.yes(r['met_second'])} | {page_tables.yes(r['met_both'])} |" for r in both])
        self.assertEqual(len(both), 20)
        size = mr.score.read_tsv(tables / "m-size.tsv")
        self.assertEqual(printed["size"], [READS, "", "| Reading | Arm | Learner | Entered | Better | Worse | Same |",
                                           "|---|---|---|---|---|---|---|"] + [
            f"| {r['m_reading']} | {page_tables.ARM_WORDS[r['arm']]} | `{r['learner']}` | "
            f"{r['entered']} of {r['datasets']} | {r['better']} of {r['entered']} | {r['worse']} of {r['entered']} | "
            f"{r['same']} of {r['entered']} |" for r in size])
        self.assertEqual(len(size), 60)

    def test_print_gives_entered_of_the_datasets_and_each_count_of_those_entered(self):
        # Every SYNTHETIC row enters all its datasets; one first-reading row of each table here enters one fewer.
        _, record = self.made()
        size_head = ["first", "unitigs", "rf"]
        count_head = ["M", "first", "lineage", "unitigs-k31", "rf", "recall"]
        for name, head in (("m-size.tsv", size_head), ("m-count.tsv", count_head)):
            lines = lines_of(record / name)
            columns = lines[0].split("\t")
            n = next(i for i, line in enumerate(lines) if line.split("\t")[:len(head)] == head)
            cells = lines[n].split("\t")
            self.assertEqual(cells[columns.index("entered")], cells[columns.index("datasets")])
            cells[columns.index("entered")] = str(int(cells[columns.index("datasets")]) - 1)
            lines[n] = "\t".join(cells)
            m_tables.edit(record, name, lambda t: "\n".join(lines))
        printed = m_page.tables(record, M_SET)
        r = [r for r in mr.score.read_tsv(record / "m-size.tsv")
             if [r["m_reading"], r["arm"], r["learner"]] == size_head][0]
        self.assertEqual(int(r["entered"]), int(r["datasets"]) - 1)
        self.assertIn(f"| first | unitigs | `rf` | {r['entered']} of {r['datasets']} | {r['better']} of {r['entered']} | "
                      f"{r['worse']} of {r['entered']} | {r['same']} of {r['entered']} |", printed["size"])
        c = [r for r in mr.score.read_tsv(record / "m-count.tsv")
             if [r["count"], r["m_reading"], r["split"], r["pair"], r["learner"], r["measure"]] == count_head][0]
        self.assertIn(f"| M alone | first | lineage | unitigs | `rf` | recall | {c['orientation']} | "
                      f"{c['entered']} of {c['datasets']} | {c['better']} | {c['worse']} | {c['same']} | "
                      f"{page_tables.d3(c['median_difference'])} | "
                      f"{page_tables.interval(c['interval_lo'], c['interval_hi'])} |", printed["beside"])

    def test_a_dataset_left_out_reads_its_reason_in_each_of_its_second_reading_cells(self):
        tables, record = self.made()
        printed = m_page.tables(record, M_SET)["datasets"]
        self.assertEqual(printed[2], "| Species | Antibiotic | Arm | Learner | First reading: arm | First reading: K31 | "
                                     "First reading: difference | Second reading: arm | Second reading: K31 | "
                                     "Second reading: difference |")
        rows = [line[2:-2].split(" | ") for line in printed[4:]]
        # 6 datasets by 6 arm rows by 5 learners on the lineage split.
        self.assertEqual(len(rows), 180)
        reason = {"drug 5": "left out, short list", "drug 6": "left out, short list",
                  "drug 10": "left out, under the floor"}
        both = mr.score.read_tsv(tables / "m-datasets.tsv")
        first = [r for r in both if r["m_reading"] == "first" and r["split"] == "lineage"]
        second = {(r["dataset"], r["pair"], r["learner"]): r for r in both
                  if r["m_reading"] == "second" and r["split"] == "lineage"}
        self.assertEqual({k[0] for k in second}, {"synth_c__drug_7", "synth_d__drug_8", "synth_d__drug_9"})
        self.assertEqual(len(first), 180)
        d3 = page_tables.d3
        for cells, r in zip(rows, first):
            want = [f"*{r['species']}*", r["antibiotic"], page_tables.pair_words(r["pair"]), f"`{r['learner']}`",
                    d3(r["arm_value"]), d3(r["default_value"]), d3(r["difference"])]
            s = second.get((r["dataset"], r["pair"], r["learner"]))
            if r["antibiotic"] in reason:
                self.assertIsNone(s)
                want += [reason[r["antibiotic"]]] * 3
            else:
                want += [d3(s["arm_value"]), d3(s["default_value"]), d3(s["difference"])]
            self.assertEqual(cells, want)
            self.assertFalse({"", "-", "0", "none"} & set(cells[7:]) and r["antibiotic"] in reason, cells)
        self.assertEqual(sum(cells[1] in reason for cells in rows), 90)
        self.assertEqual(sum(cells[1] not in reason for cells in rows), len(second))

    def test_each_table_with_the_second_reading_says_how_many_datasets_of_m_it_reads(self):
        _, record = self.made()
        printed = m_page.tables(record, M_SET)
        for name in ("count", "beside", "conditions", "both", "size", "datasets"):
            self.assertEqual(printed[name][:2], [READS, ""], name)
        for name in ("left-out", "species"):
            self.assertNotIn(READS, printed[name], name)

    def test_no_dataset_left_out_is_one_line_and_every_second_reading_cell_a_figure(self):
        _, record = self.made("count", "all in")
        printed = m_page.tables(record, M_SET)
        self.assertEqual(printed["left-out"], [NONE_LEFT])
        self.assertEqual(printed["species"][-1], "| All 2 species | 6 | 6 of 6 | 0 of 6 | 0 of 6 |")
        self.assertEqual(printed["size"][0], "Datasets of M read: 6 of 6 in the first reading, 6 of 6 in the second "
                                             "(0 left out under the floor, 0 for a short list).")
        for line in printed["datasets"][4:]:
            for c in line[2:-2].split(" | ")[4:]:
                float(c)

    def test_floor_then_record_then_print_twice_gives_the_same_bytes(self):
        tables, record = self.made("none")
        self.assertEqual(self.run_main("floor", tables, record)[0], 0)
        self.assertEqual(self.run_main("record", tables, record)[0], 0)
        one, two = self.run_main("print", record), self.run_main("print", record)
        self.assertEqual(one[0], 0)
        self.assertEqual(one[1].encode("utf-8"), two[1].encode("utf-8"))
        self.assertEqual(one[1].count("\n[") + 1, 8)
        for name in SEVEN:
            self.assertEqual((record / name).read_bytes(), (tables / name).read_bytes(), name)


class OnPage(Made):
    def page(self, lines, name="page.md"):
        path = self.top / name
        path.write_text("# A page\n\n**Left-out table.** It goes on the page as `reference/input-arms-m/m-left-out.tsv`."
                        "\n\n" + "\n".join(lines) + "\n\nMore words.\n", encoding="utf-8")
        return path

    def test_onpage_passes_a_page_holding_the_table_as_print_gives_it(self):
        _, record = self.made("floor")
        page = self.page(LEFT_OUT)
        self.assertEqual(self.run_main("onpage", record, page), (0, f"onpage: {page} holds the left-out table as "
                                                                    "print gives it\n"))
        # After the count, the same page passes the same check.
        _, full = self.made("count")
        self.assertEqual(m_page.onpage(full, page, M_SET), (True, f"onpage: {page} holds the left-out table as print "
                                                                  "gives it"))

    def test_onpage_refuses_a_page_that_names_the_path_alone(self):
        _, record = self.made("floor")
        page = self.page([])
        self.assertEqual(self.run_main("onpage", record, page),
                         (1, f"onpage: {page} lacks the left-out table's header row\n"))

    def test_onpage_names_the_first_row_missing_when_one_row_changed(self):
        _, record = self.made("floor")
        changed = LEFT_OUT[:3] + [LEFT_OUT[3].replace("| 4 R, 0 S |", "| 3 R, 0 S |")] + LEFT_OUT[4:]
        page = self.page(changed)
        rc, out = self.run_main("onpage", record, page)
        self.assertEqual((rc, out), (1, f"onpage: {page} lacks the left-out table's row of Synth c with drug 6\n"))
        self.assertNotIn("4 R", out)
        self.assertEqual(self.run_main("onpage", record, self.page(LEFT_OUT[:4], "short.md"))[1],
                         f"onpage: {self.top / 'short.md'} lacks the left-out table's row of Synth d with drug 10\n")

    def test_onpage_refuses_two_rows_in_the_other_order(self):
        _, record = self.made("floor")
        page = self.page(LEFT_OUT[:2] + [LEFT_OUT[3], LEFT_OUT[2]] + LEFT_OUT[4:])
        self.assertEqual(self.run_main("onpage", record, page),
                         (1, f"onpage: {page} lacks the left-out table's row of Synth c with drug 5\n"))

    def test_onpage_with_no_dataset_left_out_reads_its_one_line(self):
        _, record = self.made("floor", "all in")
        self.assertEqual(self.run_main("onpage", record, self.page([NONE_LEFT]))[0], 0)
        page = self.page(LEFT_OUT, "other.md")
        self.assertEqual(self.run_main("onpage", record, page),
                         (1, f"onpage: {page} lacks the left-out table's one line, that no dataset of M is left out\n"))

    def test_onpage_stops_on_a_record_that_breaks_a_gate(self):
        _, record = self.made("floor")
        (record / "notes.txt").write_text("x\n")
        self.assertEqual(self.run_main("onpage", record, self.page(LEFT_OUT)),
                         (1, "m_page: stopped: notes.txt is not on the record's list\n"))

    def test_onpage_refuses_a_row_after_the_last_row_or_above_the_header_with_no_blank_line(self):
        _, record = self.made("floor")
        made_up = "| *Synth z* | drug 99 | left out, short list | 40 R, 40 S | 1 R, 0 S | 39 R, 40 S | 25 of each class |"
        for name, lines, side in (("twice.md", LEFT_OUT + [LEFT_OUT[-1]], "after the left-out table's last row"),
                                  ("added.md", LEFT_OUT + [made_up], "after the left-out table's last row"),
                                  ("words.md", LEFT_OUT + ["More words."], "after the left-out table's last row"),
                                  ("above.md", [made_up] + LEFT_OUT, "before the left-out table's header row")):
            with self.subTest(name=name):
                page = self.page(lines, name)
                rc, out = self.run_main("onpage", record, page)
                self.assertEqual((rc, out), (1, f"onpage: {page} lacks a blank line or the file's edge {side}\n"))
                self.assertNotIn("99", out.replace(str(page), ""))

    def test_onpage_refuses_the_table_inside_a_code_block_a_comment_or_a_details_block(self):
        _, record = self.made("floor")
        for name, head, tail in (("fence.md", ["```"], ["```"]), ("fence-blank.md", ["```", ""], ["", "```"]),
                                 ("tilde.md", ["~~~~ text"], ["~~~~"]),
                                 ("indented.md", ["   ```"], ["```"]), ("open.md", ["```"], []),
                                 ("comment.md", ["<!--"], ["-->"]), ("comment-blank.md", ["<!--", ""], ["", "-->"]),
                                 ("details.md", ["<details>", ""], ["", "</details>"]),
                                 ("details-summary.md", ["<DETAILS open><summary>Rows</summary>", ""], ["", "</details>"])):
            with self.subTest(name=name):
                page = self.page(head + LEFT_OUT + tail, name)
                self.assertEqual(self.run_main("onpage", record, page),
                                 (1, f"onpage: {page} lacks the left-out table's header row\n"))
        _, none_left = self.made("floor", "all in")
        page = self.page(["<!--", NONE_LEFT, "-->"], "zero-comment.md")
        self.assertEqual(self.run_main("onpage", none_left, page),
                         (1, f"onpage: {page} lacks the left-out table's one line, that no dataset of M is left out\n"))
        page = self.page([NONE_LEFT, "More words."], "zero-run-on.md")
        self.assertEqual(self.run_main("onpage", none_left, page),
                         (1, f"onpage: {page} lacks a blank line or the file's edge after the left-out table's one "
                             "line\n"))

    def test_onpage_passes_the_table_after_closed_blocks_at_the_files_edge_and_with_crlf_line_ends(self):
        _, record = self.made("floor")
        ok = "holds the left-out table as print gives it\n"
        before = ["```", "| code |", "```", "", "<!-- a note -->", "", "<details>", "", "folded", "", "</details>", "",
                  "A `<!--` and a `<details>` in code spans hide nothing.", ""]
        page = self.page(before + LEFT_OUT, "after-blocks.md")
        self.assertEqual(self.run_main("onpage", record, page), (0, f"onpage: {page} {ok}"))
        edge = self.top / "edge.md"
        edge.write_text("\n".join(LEFT_OUT), encoding="utf-8")
        self.assertEqual(self.run_main("onpage", record, edge), (0, f"onpage: {edge} {ok}"))
        crlf = self.page(LEFT_OUT, "crlf.md")
        crlf.write_bytes(crlf.read_bytes().replace(b"\n", b"\r\n"))
        self.assertEqual(self.run_main("onpage", record, crlf), (0, f"onpage: {crlf} {ok}"))
        # The table twice, once in a comment and once as its own block, passes on the second.
        twice = self.page(["<!--"] + LEFT_OUT + ["-->", ""] + LEFT_OUT, "twice-shown.md")
        self.assertEqual(self.run_main("onpage", record, twice), (0, f"onpage: {twice} {ok}"))

    def test_onpage_stops_on_a_page_that_does_not_exist(self):
        _, record = self.made("floor")
        gone = self.top / "gone.md"
        self.assertEqual(self.run_main("onpage", record, gone), (1, f"m_page: stopped: {gone} does not exist\n"))
        (self.top / "dir.md").mkdir()
        self.assertEqual(self.run_main("onpage", record, self.top / "dir.md"),
                         (1, f"m_page: stopped: {self.top / 'dir.md'} is not a file\n"))


class Gates(Made):
    def test_gate_a_file_missing(self):
        _, record = self.made()
        (record / "m-size.tsv").unlink()
        (record / "SHA256SUMS").write_text(m_page.sums_text(record))
        self.stops(record, "^m-size.tsv is missing from the record$")
        other = m_tables.record(self.top / "two")
        shutil.copy(record / "m-floor.tsv", other / "m-floor.tsv")
        shutil.copy(record / "m-left-out.tsv", other / "m-left-out.tsv")
        self.stops(other, "^SHA256SUMS is missing from the record$")
        (other / "SHA256SUMS").write_text(m_page.sums_text(other))
        (other / "first-run-times.tsv").unlink()
        self.stops(other, "^first-run-times.tsv is missing from the record$")

    def test_gate_a_file_not_on_the_list(self):
        _, record = self.made()
        (record / "notes.txt").write_text("x\n")
        self.stops(record, "^notes.txt is not on the record's list$")
        (record / "notes.txt").unlink()
        (record / "m-extra").mkdir()
        self.stops(record, "^m-extra is not on the record's list$")

    def test_gate_a_column_list_unlike_m_readings(self):
        _, record = self.made()
        m_tables.edit(record, "m-size.tsv", lambda t: t.replace("better\tworse", "worse\tbetter", 1))
        self.stops(record, "^m-size.tsv line 1: its columns are not m_reading.py's$")

    def test_gate_a_file_unlike_its_sha256sums_line(self):
        _, record = self.made()
        m_tables.edit(record, "m-count.tsv", lambda t: t + "\n", sums=False)
        self.stops(record, "^m-count.tsv differs from its SHA256SUMS line$")
        sums = (record / "SHA256SUMS").read_text()
        (record / "SHA256SUMS").write_text(m_page.sums_text(record).replace("m-count.tsv\n", "m-other.tsv\n"))
        self.stops(record, "^SHA256SUMS line 4: m-other.tsv is no file of the record to sum$")
        (record / "SHA256SUMS").write_text("".join(l + "\n" for l in m_page.sums_text(record).splitlines()
                                                   if not l.endswith("  m-size.tsv")))
        self.stops(record, "^m-size.tsv has no SHA256SUMS line$")
        self.assertIn("m-count.tsv", sums)

    def test_gate_a_second_reading_row_of_a_dataset_not_in_it(self):
        # One dataset left out for each reason: synth_c__drug_5 for a short list, synth_d__drug_10 under the floor.
        for dataset in ("synth_c__drug_5", "synth_d__drug_10"):
            with self.subTest(dataset=dataset):
                _, record = self.made()

                def add(text):
                    row = next(line for line in text.split("\n") if line.startswith(f"first\t{dataset}\t"))
                    return text + "second" + row[len("first"):] + "\n"

                m_tables.edit(record, "m-datasets.tsv", add)
                n = len(lines_of(record / "m-datasets.tsv")) - 1
                self.stops(record, f"^m-datasets.tsv line {n}: a second-reading row of {dataset}, whose state in "
                                   "m-floor.tsv is not in$")

    def test_gate_rows_of_a_dataset_in_it_that_do_not_pair(self):
        _, record = self.made()
        gone = "\t".join(["second", "synth_d__drug_8", "Synth d", "drug 8", "lineage", "unitigs-k31", "rf"]) + "\t"
        m_tables.edit(record, "m-datasets.tsv",
                      lambda t: "".join(l + "\n" for l in t.split("\n")[:-1] if not l.startswith(gone)))
        n = line_of(record / "m-datasets.tsv", "first", "synth_d__drug_8", "lineage", "unitigs-k31", "rf")
        self.stops(record, f"^m-datasets.tsv line {n}: synth_d__drug_8 lineage unitigs-k31 rf is in the second reading "
                           "and has no second-reading row$")

    def test_gate_a_second_reading_row_counting_more_datasets_than_are_in_it(self):
        # 3 datasets of M are in the second reading: the M count's lineage rows count 3, all datasets' 11 (8 of the
        # first run and 3 of M), and the size table's 11. One more on any of them stops.
        _, record = self.made()
        self.assertEqual(len(m_page.tables(record, M_SET)), 8)
        for name, head, cells in (("m-count.tsv", ["M", "second", "lineage", "k31-aa-k31", "gbt", "balanced_accuracy"],
                                   ("\t3\t3\t", "\t4\t3\t")),
                                  ("m-count.tsv", ["all datasets", "second", "study", "program-k31", "rf", "recall"],
                                   ("\t11\t11\t", "\t11\t12\t")),
                                  ("m-size.tsv", ["second", "twin", "lasso"], ("\t11\t11\t", "\t12\t11\t"))):
            with self.subTest(name=name, row=head):
                tables, record = self.made()
                path = record / name
                lines = lines_of(path)
                n = next(i for i, l in enumerate(lines, 1) if l.split("\t")[:len(head)] == head)
                self.assertIn(cells[0], lines[n - 1])
                lines[n - 1] = lines[n - 1].replace(cells[0], cells[1], 1)
                m_tables.edit(record, name, lambda t: "\n".join(lines))
                self.stops(record, f"^{name} line {n}: the second reading's .* counts more datasets of M than "
                                   "m-floor.tsv gives the state in$")
        # The first reading's rows count every dataset of M, 6, and pass.
        _, record = self.made()
        self.assertEqual({r["datasets"] for r in mr.score.read_tsv(record / "m-count.tsv")
                          if r["count"] == "M" and r["m_reading"] == "first"}, {"6"})

    def test_gate_a_left_out_reason_unlike_its_state(self):
        _, record = self.made()
        m_tables.edit(record, "m-left-out.tsv", lambda t: t.replace("drug 5\tshort list", "drug 5\tunder the floor"))
        self.stops(record, "^m-left-out.tsv line 2: synth_c__drug_5: its reason is not its state in m-floor.tsv$")

    def test_gate_left_out_counts_unlike_the_floor_row(self):
        _, record = self.made()
        m_tables.edit(record, "m-left-out.tsv", lambda t: t.replace("\t36\t40\t25", "\t36\t41\t25"))
        self.stops(record, "^m-left-out.tsv line 3: synth_c__drug_6: its counts are not its row's in m-floor.tsv$")

    def test_gate_a_dataset_left_out_that_the_left_out_table_lacks(self):
        _, record = self.made()
        m_tables.edit(record, "m-left-out.tsv",
                      lambda t: "".join(l + "\n" for l in t.split("\n")[:-1] if not l.startswith(mr.UNDER)))
        self.stops(record, f"^m-floor.tsv line 7: {mr.UNDER} is left out and m-left-out.tsv lacks it$")
        floor_line = next(l for l in lines_of(record / "m-floor.tsv") if l.startswith("synth_c__drug_7\t"))
        cells = floor_line.split("\t")
        extra = "\t".join([cells[0], cells[2], cells[3], "in"] + cells[4:13])
        m_tables.edit(record, "m-left-out.tsv", lambda t: t + extra + "\n")
        self.stops(record, "^m-left-out.tsv line 4: synth_c__drug_7 is in the second reading by m-floor.tsv$")

    def test_gate_floor_datasets_unlike_those_lists_tsv_marks_in_m(self):
        _, record = self.made()
        with self.assertRaisesRegex(m_page.Stop, "^m-floor.tsv line 7: synth_d__drug_10 is not a dataset lists.tsv "
                                                 "marks in_M$"):
            m_page.tables(record, M_SET - {"synth_d__drug_10"})
        with self.assertRaisesRegex(m_page.Stop, "^m-floor.tsv lacks synth_x__drug_1, a dataset lists.tsv marks "
                                                 "in_M$"):
            m_page.tables(record, M_SET | {"synth_x__drug_1"})
        m_tables.edit(record, "m-floor.tsv", lambda t: t + t.split("\n")[3] + "\n")
        self.stops(record, "^m-floor.tsv line 8: synth_c__drug_7 a second time$")

    def test_gate_a_floor_state_none_of_the_three(self):
        for state in ("shown", ""):
            with self.subTest(state=state):
                _, record = self.made()
                m_tables.edit(record, "m-floor.tsv", lambda t: t.replace("\t1\tin\n", f"\t1\t{state}\n", 1))
                self.stops(record, "^m-floor.tsv line 4: synth_c__drug_7: its second reading is none of in, under the "
                                   "floor, short list$")

    def test_gate_a_written_state_unlike_the_one_enters_and_the_listed_counts_give(self):
        # synth_d__drug_10 (enters 0) and synth_c__drug_5 (1 resistant genome listed) moved to in, in both floor tables,
        # their left-out rows dropped and the sums written again; synth_c__drug_7 (5 listed, SMALLEST_LISTED) moved to a
        # short list. m_reading.py works each state out from enters and the listed counts by class.
        for dataset, n, state in (("synth_d__drug_10", 7, "in"), ("synth_c__drug_5", 2, "in"),
                                  ("synth_c__drug_7", 4, "short list")):
            with self.subTest(dataset=dataset):
                _, record = self.made("floor")
                set_cell(record, "m-floor.tsv", n, "second_reading", state)
                if state == "in":
                    m_tables.edit(record, "m-left-out.tsv",
                                  lambda t: "".join(l + "\n" for l in t.split("\n")[:-1] if not l.startswith(dataset)))
                why = (f"m-floor.tsv line {n}: {dataset}: its second reading is not the state enters and the listed "
                       "counts give")
                self.stops(record, f"^{why}$")
                page = self.top / f"{dataset}.md"
                page.write_text("\n".join(LEFT_OUT) + "\n")
                self.assertEqual(self.run_main("onpage", record, page), (1, f"m_page: stopped: {why}\n"))

    def test_gate_floor_counts_unlike_m_lists(self):
        # synth_c__drug_7, line 4: cohort 80 of 40 R and 40 S, 5 R listed, 35 R and 40 S left, floor 25, enters 1;
        # synth_d__drug_10, line 7: 24 R left, enters 0.
        for n, column, value, why in (
                (4, "cohort", "81", "its cohort is not resistant plus susceptible"),
                (4, "listed", "6", "its listed genomes are not its listed resistant plus susceptible"),
                (4, "second_resistant", "36", "its second-reading counts are not the cohort's less the list"),
                (4, "second_susceptible", "39", "its second-reading counts are not the cohort's less the list"),
                (4, "floor", "20", "its floor is not m_list.py's"),
                (7, "enters", "1", "its enters is not 1 exactly when each class keeps the floor"),
                (4, "enters", "0", "its enters is not 1 exactly when each class keeps the floor")):
            with self.subTest(column=column, value=value):
                _, record = self.made()
                set_cell(record, "m-floor.tsv", n, column, value)
                dataset = lines_of(record / "m-floor.tsv")[n - 1].split("\t")[0]
                with self.assertRaises(m_page.Stop) as caught:
                    m_page.tables(record, M_SET)
                self.assertEqual(str(caught.exception), f"m-floor.tsv line {n}: {dataset}: {why}")

    def test_gate_a_met_cell_that_is_not_0_1_or_empty(self):
        for name, n, column, value in (("m-conditions.tsv", 2, "met", "0.5"), ("m-conditions.tsv", 2, "met", "yes"),
                                       ("m-both.tsv", 3, "met_both", "2"), ("m-both.tsv", 3, "met_first", "")):
            with self.subTest(name=name, column=column, value=value):
                _, record = self.made()
                set_cell(record, name, n, column, value)
                if value == "":
                    # An empty met is a met cell's form; the both gate stops it, since m_reading.py writes no such row.
                    self.stops(record, f"^{name} line {n}: its row is not m-conditions.tsv's$")
                    continue
                with self.assertRaises(m_page.Stop) as caught:
                    m_page.tables(record, M_SET)
                self.assertEqual(str(caught.exception), f"{name} line {n}: its {column} is not 0, 1 or empty")

    def test_gate_a_count_cell_that_is_not_a_whole_number(self):
        _, base = self.made()
        count = lines_of(base / "m-count.tsv")
        better = count[1].split("\t")[count[0].split("\t").index("better")]
        for name, n, column, value in (("m-count.tsv", 2, "better", better + ".0"), ("m-count.tsv", 2, "better", "1.5"),
                                       ("m-count.tsv", 2, "species", "-1"), ("m-size.tsv", 2, "entered", "2e1"),
                                       ("m-datasets.tsv", 2, "draws", "3.0"), ("m-floor.tsv", 2, "listed", "1.0"),
                                       ("m-left-out.tsv", 2, "floor", " 25"), ("m-conditions.tsv", 2, "learners", "5.0"),
                                       ("m-count.tsv", 2, "entered", "")):
            with self.subTest(name=name, column=column, value=value):
                _, record = self.made()
                set_cell(record, name, n, column, value)
                with self.assertRaises(m_page.Stop) as caught:
                    m_page.tables(record, M_SET)
                self.assertEqual(str(caught.exception), f"{name} line {n}: its {column} is not a whole number")
        # The stricter reading's condition 3 has no count of learners, and passes.
        self.assertIn("\tstricter\t3\t", (base / "m-conditions.tsv").read_text())
        self.assertEqual(len(m_page.tables(base, M_SET)), 8)

    def test_gate_a_dataset_of_m_with_no_first_reading_row_on_the_lineage_split(self):
        # Every m-datasets.tsv row of synth_d__drug_8 (in), the first-reading rows of synth_c__drug_5 (a short list),
        # or its lineage rows alone.
        for dataset, n, drop in (("synth_d__drug_8", 5, lambda l: l.split("\t")[1] == "synth_d__drug_8"),
                                 ("synth_c__drug_5", 2, lambda l: l.startswith("first\tsynth_c__drug_5\t")),
                                 ("synth_c__drug_5", 2, lambda l: l.startswith("first\tsynth_c__drug_5\t")
                                  and l.split("\t")[4] == "lineage")):
            with self.subTest(dataset=dataset):
                _, record = self.made()
                m_tables.edit(record, "m-datasets.tsv",
                              lambda t: "".join(l + "\n" for l in t.split("\n")[:-1] if not drop(l)))
                self.stops(record, f"^m-floor.tsv line {n}: {dataset} has no first-reading row on the lineage split "
                                   "in m-datasets.tsv$")

    def test_gate_a_first_run_row_of_the_second_reading_unlike_its_first(self):
        _, record = self.made()
        count = lines_of(record / "m-count.tsv")
        head = count[0].split("\t")
        n = next(i for i, l in enumerate(count, 1)
                 if l.split("\t")[:6] == ["first run", "second", "lineage", "unitigs-k31", "rf", "balanced_accuracy"])
        key = "lineage unitigs-k31 rf balanced_accuracy"
        for column in ("datasets", "median_difference", "holm_p"):
            with self.subTest(column=column):
                _, record = self.made()
                cell = count[n - 1].split("\t")[head.index(column)]
                new = str(int(cell) + 1) if column == "datasets" else ("0.5" if cell != "0.5" else "0.25")
                set_cell(record, "m-count.tsv", n, column, new)
                self.stops(record, f"^m-count.tsv line {n}: the first run's row for {key} is not the same in both "
                                   "readings$")
        _, record = self.made()
        m_tables.edit(record, "m-count.tsv", lambda t: "\n".join(l for i, l in enumerate(t.split("\n"), 1) if i != n))
        first = next(i for i, l in enumerate(lines_of(record / "m-count.tsv"), 1)
                     if l.split("\t")[:6] == ["first run", "first", "lineage", "unitigs-k31", "rf", "balanced_accuracy"])
        self.stops(record, f"^m-count.tsv line {first}: the first run's row for {key} has no second-reading row$")
        _, record = self.made()
        cond = lines_of(record / "m-conditions.tsv")
        k = next(i for i, l in enumerate(cond, 1)
                 if l.split("\t")[:5] == ["first run", "second", "unitigs", "ruling 15", "4"])
        met = cond[k - 1].split("\t")[6]
        set_cell(record, "m-conditions.tsv", k, "met", "0" if met == "1" else "1")
        self.stops(record, f"^m-conditions.tsv line {k}: the first run's row for unitigs ruling 15 4 is not the same in "
                           "both readings$")

    def test_gate_a_both_row_unlike_both_readings_or_unlike_the_conditions_table(self):
        # m-both.tsv line 3: unitigs, condition 3, met in the second reading alone; line 2: unitigs, 1 and 2, met in both.
        _, record = self.made()
        both = lines_of(record / "m-both.tsv")
        self.assertEqual(both[2].split("\t")[3:], ["0", "1", "0"])
        self.assertEqual(both[1].split("\t")[3:], ["1", "1", "1"])
        set_cell(record, "m-both.tsv", 3, "met_both", "1")
        self.stops(record, "^m-both.tsv line 3: its Both is not met exactly when both readings are$")
        _, record = self.made()
        for column in ("met_first", "met_second", "met_both"):
            set_cell(record, "m-both.tsv", 2, column, "0")
        self.stops(record, "^m-both.tsv line 2: its row is not m-conditions.tsv's$")
        _, record = self.made()
        m_tables.edit(record, "m-both.tsv", lambda t: t.replace("; at least 3 of 5 carry\t", "; at least 2 of 5 carry\t", 1))
        self.stops(record, "^m-both.tsv line 2: its row is not m-conditions.tsv's$")
        _, record = self.made()
        m_tables.edit(record, "m-both.tsv", lambda t: t + t.split("\n")[1] + "\n")
        self.stops(record, f"^m-both.tsv line {len(both)}: unitigs 1 and 2 a second time$")
        _, record = self.made()
        m_tables.edit(record, "m-both.tsv", lambda t: "\n".join(l for i, l in enumerate(t.split("\n"), 1) if i != 2))
        n = next(i for i, l in enumerate(lines_of(record / "m-conditions.tsv"), 1)
                 if l.split("\t")[:5] == ["all datasets", "first", "unitigs", "ruling 15", "1 and 2"])
        self.stops(record, f"^m-conditions.tsv line {n}: m-both.tsv lacks its row for unitigs 1 and 2$")

    def test_gate_a_token_of_a_genome_ids_shape(self):
        # pack.py's shape with 2 or more digits before the point stops, naming the line and not the token; one digit
        # before it, or the token after a letter, passes. Each made-up token starts 00, as no taxon id does, so no real
        # genome id can equal it. The token goes in the first row of the count on M alone, which no other row repeats.
        for token, stops in (("00.000", True), ("0.5", False), ("x00.000", False), ("_00.0", True),
                             ("0.123456789012", False), ("000.00.0", True)):
            with self.subTest(token=token):
                _, record = self.made()
                n = next(i for i, l in enumerate(lines_of(record / "m-conditions.tsv"), 1)
                         if l.startswith("M\tfirst\t") and "after Holm over 30" in l)

                def put(text):
                    lines = text.split("\n")
                    lines[n - 1] = lines[n - 1].replace("after Holm over 30", f"after {token}")
                    return "\n".join(lines)

                m_tables.edit(record, "m-conditions.tsv", put)
                if stops:
                    with self.assertRaises(m_page.Stop) as caught:
                        m_page.tables(record, M_SET)
                    self.assertEqual(str(caught.exception), f"m-conditions.tsv line {n}: a token of a genome id's shape, "
                                                            "digits before the point")
                    self.assertNotIn(token.strip("x_"), str(caught.exception))
                else:
                    self.assertEqual(len(m_page.tables(record, M_SET)), 8)

    def test_each_command_stops_on_a_new_gate_where_it_applies(self):
        # floor on a state unlike the one its counts give, in the tables directory; record, print, onpage and page on a
        # both row met in both where the second reading is not met.
        tables, record = self.made("none")
        set_cell(tables, "m-floor.tsv", 7, "second_reading", "in", sums=False)
        m_tables.edit(tables, "m-left-out.tsv", lambda t: "".join(l + "\n" for l in t.split("\n")[:-1]
                                                                  if not l.startswith(mr.UNDER)), sums=False)
        self.assertEqual(self.run_main("floor", tables, record), (
            1, f"m_page: stopped: m-floor.tsv line 7: {mr.UNDER}: its second reading is not the state enters and the "
               "listed counts give\n"))
        self.assertEqual(sorted(p.name for p in record.iterdir()), ["README.md", "first-run-times.tsv"])
        why = "m_page: stopped: m-both.tsv line 3: its Both is not met exactly when both readings are\n"
        tables, record = self.made("floor")
        set_cell(tables, "m-both.tsv", 3, "met_both", "1", sums=False)
        self.assertEqual(self.run_main("record", tables, record), (1, why))
        self.assertEqual(m_page.state(record), "floor")
        _, record = self.made()
        set_cell(record, "m-both.tsv", 3, "met_both", "1")
        page = self.top / "page.md"
        page.write_text("\n".join(LEFT_OUT) + "\n")
        for argv in (("print", record), ("onpage", record, page), ("page", record)):
            self.assertEqual(self.run_main(*argv), (1, why), argv[0])

    def test_record_print_and_onpage_each_stop_on_a_gate(self):
        tables, record = self.made("floor")
        page = self.top / "page.md"
        page.write_text("\n".join(LEFT_OUT) + "\n")
        m_tables.edit(record, "m-left-out.tsv", lambda t: t.replace("drug 6\tshort list", "drug 6\tunder the floor"))
        why = "m_page: stopped: m-left-out.tsv line 3: synth_c__drug_6: its reason is not its state in m-floor.tsv\n"
        self.assertEqual(self.run_main("print", record), (1, why))
        self.assertEqual(self.run_main("onpage", record, page), (1, why))
        self.assertEqual(self.run_main("record", tables, record)[0], 1)


class States(Made):
    def test_the_state_of_a_record_follows_the_tables_it_holds(self):
        tables, record = self.made("none")
        self.assertEqual(m_page.state(record), "none")
        m_page.floor(tables, record, M_SET)
        self.assertEqual(m_page.state(record), "floor")
        m_page.record_count(tables, record, M_SET)
        self.assertEqual(m_page.state(record), "count")
        self.assertEqual([line.split(" | ")[0] for line in m_page.count_headers()],
                         ["| Count", "| Count", "| Count", "| Arm", "| Reading", "| Species", "| Reading", "| Reading"])

    def test_five_tables_stay_in_the_record_and_no_block_of_the_page_holds_their_header(self):
        _, record = self.made()
        printed, blocks = m_page.tables(record, M_SET), m_page.page(record, M_SET)
        heads = m_page.record_headers()
        self.assertEqual(heads, [printed[name][2] for name in ("count", "beside", "conditions", "size",
                                                               "datasets")])
        self.assertEqual(len(set(heads)), 5)
        self.assertLess(set(heads), set(m_page.count_headers()))
        on = {line for lines in blocks.values() for line in lines}
        self.assertEqual(on & set(heads), set())
        # Each table row of a block is a header of the count, so no block stands on a page before the count.
        self.assertEqual({lines[2] for name, lines in blocks.items() if name not in ("reads", "species", "files")}
                         - set(m_page.count_headers()), set())
        # The both and species tables' own headers; the two lineage counts share one, and recall and the study split
        # theirs.
        self.assertEqual(len({lines[2] for name, lines in blocks.items() if name not in ("reads", "files")}), 4)

    def test_no_row_of_the_five_tables_that_stay_in_the_record_is_a_line_of_the_page(self):
        # The rows print gives the five, and each count's rows on every split and measure in the page's own forms, but
        # those the page gives: the record test fails on any of them on a page after the count.
        _, record = self.made()
        printed, blocks = m_page.tables(record, M_SET), m_page.page(record, M_SET)
        on = {line for lines in blocks.values() for line in lines}
        rows = {line for name in ("count", "beside", "conditions", "size", "datasets") for line in printed[name][4:]}
        self.assertEqual(len(rows), 180 + 570 + 120 + 60 + 180)
        self.assertEqual(rows & on, set())
        count = m_page.read_tables(record, ["m-count.tsv"])["m-count.tsv"]
        forms = set()
        for c in m_page.COUNTS:
            forms |= set(m_page.print_count(count, c)[2:])
            for split, measure in {(r["split"], r["measure"]) for _, r in count.rows}:
                forms |= set(m_page.print_other(count, c, split, measure)[2:])
        self.assertLess({line for name in ("count-all", "recall-all", "study-all", "count-m")
                         for line in blocks[name][4:]}, forms)


if __name__ == "__main__":
    unittest.main()
