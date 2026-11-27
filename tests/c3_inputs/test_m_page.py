"""scripts/c3-inputs/m_page.py: the record of the run on the datasets of M and the two pages written from it.

Every table here is SYNTHETIC: m_reading.py's own floor and count forms
write the 7 tables from test_m_reading.py's made-up species directories
(tests/c3_inputs/m_tables.py), 6 made-up datasets of M in two made-up
species, no genome. Only the column lists, words, lists and number forms are
read from the committed m_reading.py, score.py and page_tables.py, and the
count of the datasets of M from the committed lists.tsv. A test of a rule
names the lines and words it expects itself.
"""

import contextlib
import io
import os
import re
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from tests.c3_inputs import common, m_tables
from tests.c3_inputs import test_m_reading as mr

m_page = m_tables.m_page
reading = mr.reading
score = mr.score
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
LEFT_AT = "The left-out page, `docs/input-arms-m-left-out.md`, names each dataset left out."
WRITTEN = ("The script `scripts/c3-inputs/m_page.py {}` writes this page whole from the record, "
           "`reference/input-arms-m/`, and no one edits it by hand.")
# The left-out page as left gives it on the SYNTHETIC record: 11 lines, the table's rows on lines 9 to 11.
LEFT = ["# Left-out datasets of M", "", WRITTEN.format("left"), "", READS, ""] + LEFT_OUT
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
    # A copy of m_page.py and pack.py beside a SYNTHETIC lists.tsv, for runs as the lead runs the script.
    BUILT["kit"] = top / "kit"
    BUILT["kit"].mkdir()
    for name in ("m_page.py", "pack.py"):
        shutil.copy(common.KIT / name, BUILT["kit"] / name)
    m_tables.lists_file(BUILT["kit"] / "lists.tsv")


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


def change_rows(directory, name, change, sums=True):
    """One table of <directory> rewritten through change(header, rows) -> rows, rows as {column: cell}."""
    def put(text_):
        lines = text_.split("\n")
        header = lines[0].split("\t")
        rows = [dict(zip(header, line.split("\t"))) for line in lines[1:-1]]
        return "\n".join(["\t".join(header)] + ["\t".join(r[c] for c in header) for r in change(header, rows)]) + "\n"
    m_tables.edit(directory, name, put, sums)


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

    def run_bytes(self, *argv):
        """(exit, the bytes main writes to standard output)."""
        lists = m_tables.lists_file(self.top / "lists.tsv")
        out = io.TextIOWrapper(io.BytesIO(), encoding="utf-8", newline="")
        with contextlib.redirect_stdout(out):
            rc = m_page.main([str(a) for a in argv], lists=lists)
        out.flush()
        return rc, out.buffer.getvalue()

    def run_main(self, *argv):
        rc, out = self.run_bytes(*argv)
        return rc, out.decode("utf-8")

    def cli(self, *argv):
        """(exit, standard output, standard error) of the kit's copy of m_page.py run as the lead runs it."""
        done = subprocess.run([sys.executable, "-I", "-B", str(BUILT["kit"] / "m_page.py"), *map(str, argv)],
                              capture_output=True, timeout=300)
        return done.returncode, done.stdout, done.stderr

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
        for name, cells in (list(m_page.COUNT_CELLS.items()) + list(m_page.MET_CELLS.items())
                            + list(m_page.FIGURE_CELLS.items())):
            self.assertLess(set(cells), set(m_page.COLUMNS[name]), name)

    def test_the_lists_are_score_pys_and_m_readings_own(self):
        self.assertEqual((list(m_page.LEARNERS), list(m_page.FAMILY), list(m_page.MEASURES), m_page.TWIN),
                         (list(score.LEARNERS), list(score.FAMILY), list(score.MEASURES), score.TWIN))
        self.assertEqual((set(m_page.LOWER_IS_BETTER), m_page.HOLM_ALPHA), (set(score.LOWER_IS_BETTER),
                                                                           score.HOLM_ALPHA))
        self.assertEqual((m_page.RECORD_MEASURES, m_page.SPLITS), (reading.RECORD_MEASURES, page_tables.SPLITS))
        self.assertEqual(m_page.COUNT_MEASURES, {"first run": reading.RECORD_MEASURES, "M": tuple(score.MEASURES),
                                                 "all datasets": reading.RECORD_MEASURES})

    def test_the_rows_the_lists_give_are_the_rows_m_reading_writes_in_its_order(self):
        # m_reading.py's tables from the SYNTHETIC species directories, row for row: 1,100 count rows, 330 condition
        # rows and 60 size rows; and the conditions by reading are score.py's own.
        tables = m_page.read_tables(m_tables.seven(Path(self.enterContext(tempfile.TemporaryDirectory()))), SEVEN)
        for name, n in (("m-count.tsv", 1100), ("m-conditions.tsv", 330), ("m-size.tsv", 60)):
            have = [tuple(r[c] for c in m_page.ROW_KEYS[name]) for _, r in tables[name].rows]
            self.assertEqual(have, m_page.rows_given(name), name)
            self.assertEqual(len(have), n, name)
        pairs = {(r["reading"], r["condition"]) for _, r in tables["m-conditions.tsv"].rows}
        self.assertEqual(pairs, {(k, c) for k, cs in m_page.CONDITIONS.items() for c in cs})
        self.assertEqual({(r["orientation"], r["measure"] in score.LOWER_IS_BETTER)
                          for _, r in tables["m-count.tsv"].rows}, {("higher", False), ("lower", True)})
        self.assertEqual({r["condition"] for _, r in tables["m-both.tsv"].rows}, set(m_page.SHOWN_CONDITIONS))

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

    def test_m_page_runs_where_numpy_cannot_be_imported(self):
        # numpy stands as a module that fails to import, as on a machine that lacks it; m_page.py still loads and
        # gives its use line, so it imports numpy through no script.
        code = ("import runpy, sys; sys.modules['numpy'] = None; sys.argv = ['m_page.py']; "
                f"runpy.run_path({str(common.KIT / 'm_page.py')!r}, run_name='__main__')")
        done = subprocess.run([sys.executable, "-I", "-B", "-c", code], capture_output=True, text=True, timeout=300)
        self.assertEqual(done.returncode, 2, done.stderr)
        self.assertEqual(done.stderr.split("\n")[0], "  python3 -I m_page.py floor <tables dir> <record dir>")
        self.assertEqual(len(done.stderr.strip().split("\n")), 6)


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
        # The left-out page takes the rows in the same order, on its lines 9 to 11.
        self.assertEqual(m_page.left_page(record, M_SET).decode("utf-8").split("\n")[6:], left + [""])

    def test_print_with_the_floor_tables_alone_says_the_count_has_not_been_taken(self):
        _, record = self.made("floor")
        self.assertEqual(m_page.tables(record, M_SET), {"left-out": LEFT_OUT, "species": SPECIES,
                                                        "count": [NOT_TAKEN]})
        rc, out = self.run_main("print", record)
        self.assertEqual(rc, 0)
        self.assertEqual(out, "[left-out]\n" + "\n".join(LEFT_OUT) + "\n\n[species]\n" + "\n".join(SPECIES)
                         + "\n\n[count]\n" + NOT_TAKEN + "\n\n")

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
        # Every SYNTHETIC row enters all its datasets. Here one first-reading row of the count on M alone (recall) and
        # one of the count over all datasets (lineage, balanced accuracy) enter one fewer, with one fewer better, so
        # better, worse and same still sum to entered and better still beats worse; the size table's row copies the
        # second, as m_reading.py copies it.
        _, record = self.made()
        size_head = ["first", "unitigs", "rf"]
        count_head = ["M", "first", "lineage", "unitigs-k31", "rf", "recall"]
        for name, head in (("m-size.tsv", size_head), ("m-count.tsv", count_head),
                           ("m-count.tsv", ["all datasets", "first", "lineage", "unitigs-k31", "rf",
                                            "balanced_accuracy"])):
            def change(header, rows, head=head):
                r = next(r for r in rows if [r[c] for c in header[:len(head)]] == head)
                self.assertEqual(r["entered"], r["datasets"])
                self.assertGreater(int(r["better"]), int(r["worse"]) + 1)
                r["entered"], r["better"] = str(int(r["entered"]) - 1), str(int(r["better"]) - 1)
                return rows
            change_rows(record, name, change)
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


# The count page's run-in lines and its Files line, as this test writes them.
BOTH_LINE = "**Both readings.** Per arm, conditions 1 to 5 of ruling 15 on the count over all datasets."
SPECIES_LINE = ("**Datasets read.** The datasets of M by species: in the second reading, left out under the floor, or "
                "left out for a short list.")
ALL_LINE = ("**All datasets.** The first run's datasets and M's, on the lineage split and balanced accuracy. "
            "Condition 1 is 1 where the interval lies wholly above 0. Condition 2 is 1 where better beats worse and "
            "the Holm p, the sign test's p after Holm's correction, is under 0.05. Each mark comes from the record's "
            "full figures, not the 3 decimals shown.")
RECALL_LINE = "**All datasets, recall.** The first run's datasets and M's, on the lineage split and recall."
STUDY_LINE = "**All datasets, study split.** The first run's datasets and M's, on the study split and balanced accuracy."
M_LINE = "**M alone.** The datasets of M, on the lineage split and balanced accuracy. This count decides nothing."
M_RECALL_LINE = "**M alone, recall.** The datasets of M, on the lineage split and recall."
M_STUDY_LINE = "**M alone, study split.** The datasets of M, on the study split and balanced accuracy."
FILES = ("**Files.** The 8 tables above, 363 rows for the 2 species of M, come from the record, "
         "`reference/input-arms-m/`, which holds every other table: each count on every measure and split "
         "(`m-count.tsv`), each count's conditions (`m-conditions.tsv`), the size table (`m-size.tsv`) and each "
         "dataset's figures (`m-datasets.tsv`). `scripts/c3-inputs/m_page.py print` prints them.")
ALL_READS = ("Datasets of M read: 6 of 6 in the first reading, 6 of 6 in the second (0 left out under the floor, 0 for a "
             "short list).")
FULL = ["Reading", "Arm", "Learner", "Better", "Worse", "Same", "Median", "Interval", "Sign p", "Holm p", "Constant"]


def text(lines):
    return ("\n".join(lines) + "\n").encode("utf-8")


def marked(r):
    """Conditions 1 and 2 of one count row as this test works them out from the row's full figures."""
    one = r["interval_lo"] != "" and float(r["interval_lo"]) > 0
    two = r["holm_p"] != "" and float(r["holm_p"]) < 0.05 and int(r["better"]) > int(r["worse"])
    return [str(int(one)), str(int(two))]


def other_cells(r):
    pd = page_tables
    return [r["m_reading"], pd.pair_words(r["pair"]), f"`{r['learner']}`", pd.of(r["better"], r["entered"]),
            pd.of(r["worse"], r["entered"]), pd.of(r["same"], r["entered"]), pd.d3(r["median_difference"]),
            pd.interval(r["interval_lo"], r["interval_hi"])]


def full_cells(r):
    return other_cells(r) + [page_tables.pval(r["sign_p"]), page_tables.pval(r["holm_p"]), r["constant"]]


def want_page(record):
    """(the count page's lines as this test writes them from the record's tables, the rows of each of its 8 tables)."""
    pd = page_tables
    count = mr.score.read_tsv(record / "m-count.tsv")

    def rows(only, split, measure):
        return [r for r in count if (r["count"], r["split"], r["measure"]) == (only, split, measure)]

    both = [[pd.ARM_WORDS[r["arm"]], r["condition"], pd.yes(r["met_first"]), pd.yes(r["met_second"]),
             pd.yes(r["met_both"])] for r in mr.score.read_tsv(record / "m-both.tsv")]
    blocks = [
        (BOTH_LINE, pd.table(["Arm", "Condition", "First reading", "Second reading", "Both"], both)),
        (SPECIES_LINE, SPECIES),
        (ALL_LINE, pd.table(FULL + ["Condition 1", "Condition 2"],
                            [full_cells(r) + marked(r) for r in rows("all datasets", "lineage", "balanced_accuracy")])),
        (RECALL_LINE, pd.table(FULL[:8], [other_cells(r) for r in rows("all datasets", "lineage", "recall")])),
        (STUDY_LINE, pd.table(FULL[:8], [other_cells(r) for r in rows("all datasets", "study", "balanced_accuracy")])),
        (M_LINE, pd.table(FULL, [full_cells(r) for r in rows("M", "lineage", "balanced_accuracy")])),
        (M_RECALL_LINE, pd.table(FULL[:8], [other_cells(r) for r in rows("M", "lineage", "recall")])),
        (M_STUDY_LINE, pd.table(FULL[:8], [other_cells(r) for r in rows("M", "study", "balanced_accuracy")]))]
    lines = ["# Count tables of M", "", WRITTEN.format("page"), "", f"{READS} {LEFT_AT}"]
    for line, block in blocks:
        lines += ["", line, ""] + block
    return lines + ["", FILES], [len(block) - 2 for _, block in blocks]


def kinds(data):
    """Each line of a written page by kind: title, table, blank or paragraph."""
    out = []
    for line in data.decode("utf-8").split("\n")[:-1]:
        out.append(("title" if line.startswith("#") else "table" if line.startswith("|") else
                    "blank" if line == "" else "paragraph", line))
    return out


class Written(Made):
    def test_left_gives_the_left_out_page_as_written_here_after_the_floor_and_after_the_count(self):
        for upto in ("floor", "count"):
            with self.subTest(upto=upto):
                _, record = self.made(upto)
                self.assertEqual(m_page.left_page(record, M_SET), text(LEFT))
                self.assertEqual(self.run_bytes("left", record), (0, text(LEFT)))

    def test_left_with_no_dataset_left_out_gives_its_one_line(self):
        _, record = self.made("floor", "all in")
        self.assertEqual(m_page.left_page(record, M_SET), text(LEFT[:4] + [ALL_READS, "", NONE_LEFT]))

    def test_left_stops_before_the_floor(self):
        _, record = self.made("none")
        self.assertEqual(self.run_bytes("left", record), (1, b"m_page: stopped: m-floor.tsv is missing from the record\n"))

    def test_page_gives_the_count_pages_eight_tables_from_the_record(self):
        # The tables in their order with their rows: conditions 1 to 5 in both readings, 20; the datasets by species,
        # 3 for the 2 made-up species; the count over all datasets with its marks, 60, on recall, 60, and on the study
        # split, 50, which holds no twin; and M alone on the same three.
        _, record = self.made()
        lines, rows = want_page(record)
        self.assertEqual(rows, [20, 3, 60, 60, 50, 60, 60, 50])
        self.assertEqual(sum(rows), 363)
        self.assertEqual(m_page.count_page(record, M_SET), text(lines))
        self.assertEqual(self.run_bytes("page", record), (0, text(lines)))
        self.assertEqual([line for kind, line in kinds(text(lines)) if line.startswith("**")],
                         [BOTH_LINE, SPECIES_LINE, ALL_LINE, RECALL_LINE, STUDY_LINE, M_LINE, M_RECALL_LINE,
                          M_STUDY_LINE, FILES])

    def test_the_count_page_gives_the_reads_line_once_above_its_first_table(self):
        _, record = self.made()
        lines = m_page.count_page(record, M_SET).decode("utf-8").split("\n")
        self.assertEqual(lines[4:7], [f"{READS} {LEFT_AT}", "", BOTH_LINE])
        self.assertEqual(sum(READS in line for line in lines), 1)
        self.assertEqual(sum(LEFT_AT in line for line in lines), 1)

    def test_page_stops_before_the_count(self):
        for upto, why in (("none", "m-floor.tsv is missing from the record"),
                          ("floor", "the count has not been taken, and the count page takes its tables after it")):
            with self.subTest(upto=upto):
                _, record = self.made(upto)
                self.assertEqual(self.run_bytes("page", record), (1, f"m_page: stopped: {why}\n".encode("utf-8")))

    def test_the_marks_come_from_the_records_full_figures_and_not_the_three_decimals_shown(self):
        # A low end of 0.0004 prints 0.000 and lies above 0; a Holm p of 0.0495 prints 0.050 and is under 0.05, and one
        # of 0.0505 prints 0.050 too and is not.
        self.assertEqual([m_page.d3("0.0004"), m_page.pval("0.0495"), m_page.pval("0.0505")], ["0.000", "0.050", "0.050"])
        for cells, want in (({"interval_lo": "0.0004", "holm_p": "0.0495", "better": "3", "worse": "1"}, (1, 1)),
                            ({"interval_lo": "-0.0004", "holm_p": "0.0505", "better": "3", "worse": "1"}, (0, 0)),
                            ({"interval_lo": "1e-05", "holm_p": "1e-05", "better": "2", "worse": "0"}, (1, 1)),
                            ({"interval_lo": "0", "holm_p": "0.01", "better": "1", "worse": "1"}, (0, 0)),
                            ({"interval_lo": "0.2", "holm_p": "0.01", "better": "1", "worse": "2"}, (1, 0)),
                            ({"interval_lo": "", "holm_p": "", "better": "0", "worse": "0"}, (0, 0))):
            self.assertEqual(m_page.marks(cells), want, cells)
        # On the page: one row of all datasets given those figures keeps both marks, so every gate passes, and the
        # page marks it 1 and 1 beside 0.000 and 0.050.
        _, record = self.made()
        head = ["all datasets", "first", "lineage", "unitigs-k31", "scm", "balanced_accuracy"]

        def change(header, rows):
            r = next(r for r in rows if [r[c] for c in header[:6]] == head)
            self.assertEqual(marked(r), ["1", "1"])
            r["interval_lo"], r["holm_p"] = "0.0004", "0.0495"
            return rows

        change_rows(record, "m-count.tsv", change)
        lines, _ = want_page(record)
        self.assertEqual(m_page.count_page(record, M_SET), text(lines))
        row = next(line for line in lines if line.startswith("| first | unitigs | `scm` |") and line.count(" | ") == 12)
        cells = row[2:-2].split(" | ")
        self.assertEqual((cells[7][:9], cells[9], cells[11:]), ("0.000 to ", "0.050", ["1", "1"]))

    def test_each_written_page_ends_in_one_newline_and_holds_no_html_and_no_code_block(self):
        _, record = self.made()
        _, all_in = self.made("count", "all in")
        for data in (m_page.left_page(record, M_SET), m_page.left_page(all_in, M_SET),
                     m_page.count_page(record, M_SET), m_page.count_page(all_in, M_SET)):
            self.assertTrue(data.endswith(b"\n") and not data.endswith(b"\n\n"))
            self.assertNotIn(b"\r", data)
            self.assertIsNone(re.search(rb"<[A-Za-z!/?]", data))
            for line in data.split(b"\n"):
                self.assertFalse(line.startswith((b"```", b"~~~", b"    ", b"\t", b">", b"- ", b"* ")), line)
            self.assertNotIn(b"\n\n\n", data)

    def test_each_written_page_keeps_the_guides_limits(self):
        # Under 85 words a paragraph and 30 a cell, at most 12 titles, and no list.
        _, record = self.made()
        for data in (m_page.left_page(record, M_SET), m_page.count_page(record, M_SET)):
            lines = kinds(data)
            self.assertLessEqual(sum(kind == "title" for kind, _ in lines), 12)
            for kind, line in lines:
                if kind == "paragraph":
                    self.assertLess(len(line.split()), 85, line)
                if kind == "table":
                    for cell in line[2:-2].split(" | "):
                        self.assertLess(len(cell.split()), 30, cell)

    def test_left_and_page_from_the_command_line_give_the_functions_bytes(self):
        _, record = self.made()
        for command, want in (("left", m_page.left_page(record, M_SET)), ("page", m_page.count_page(record, M_SET))):
            with self.subTest(command=command):
                self.assertEqual(self.cli(command, record), (0, want, b""))

    def test_left_and_page_give_the_same_bytes_each_time(self):
        _, record = self.made()
        for command in ("left", "page"):
            self.assertEqual(self.run_bytes(command, record), self.run_bytes(command, record))


def gives(page, words):
    return f"onpage: {page} {words} the bytes `m_page.py left` gives\n"


class OnPage(Made):
    def page(self, lines, name="page.md", data=None):
        path = self.top / name
        path.write_bytes(text(lines) if data is None else data)
        return path

    def test_onpage_passes_the_bytes_left_gives_and_no_other(self):
        for upto in ("floor", "count"):
            with self.subTest(upto=upto):
                _, record = self.made(upto)
                page = self.page(LEFT)
                self.assertEqual(self.run_main("onpage", record, page),
                                 (0, f"onpage: {page} is the bytes `m_page.py left` gives\n"))
                # The left-out table alone, as the design page held it, is not the page.
                page = self.page(LEFT_OUT)
                self.assertEqual(self.run_main("onpage", record, page), (1, gives(page, "line 1 differs from")))

    def test_onpage_names_the_first_line_that_differs_and_no_text_or_figure_of_the_page(self):
        _, record = self.made("floor")
        for n, old, new in ((9, "| 40 R, 40 S |", "| 41 R, 40 S |"), (11, "| 24 R, 38 S |", "| 23 R, 38 S |"),
                            (10, "drug 6", "drug 7"), (5, "6 of 6", "7 of 6"), (1, "# Left", "# left")):
            with self.subTest(line=n, old=old):
                lines = list(LEFT)
                self.assertEqual(lines[n - 1].count(old), 1)
                lines[n - 1] = lines[n - 1].replace(old, new)
                page = self.page(lines)
                rc, out = self.run_main("onpage", record, page)
                self.assertEqual((rc, out), (1, gives(page, f"line {n} differs from")))
                self.assertEqual(re.findall(r"[0-9]+", out[len(f"onpage: {page} "):]), [str(n)])

    def test_onpage_refuses_a_page_respaced_or_with_other_line_ends(self):
        _, record = self.made("floor")
        for n, line in ((9, LEFT[8] + " "), (9, " " + LEFT[8]), (10, LEFT[9].replace(" | ", "|")),
                        (8, "| --- | --- | --- | --- | --- | --- | --- |"), (3, LEFT[2] + "\t"), (6, " "),
                        (7, LEFT[6].replace("| Species |", "|Species |"))):
            with self.subTest(line=n, text=line):
                lines = list(LEFT)
                lines[n - 1] = line
                page = self.page(lines)
                self.assertEqual(self.run_main("onpage", record, page), (1, gives(page, f"line {n} differs from")))
        for data, n in ((text(LEFT).replace(b"\n", b"\r\n"), 1), (text(LEFT)[:-1], 11),
                        (b"\xef\xbb\xbf" + text(LEFT), 1), (text(LEFT).replace(b" R, ", b"\xc2\xa0R, ", 1), 9)):
            with self.subTest(data=data[:12]):
                page = self.page(None, data=data)
                self.assertEqual(self.run_main("onpage", record, page), (1, gives(page, f"line {n} differs from")))

    def test_onpage_refuses_a_page_that_ends_early_or_runs_on(self):
        _, record = self.made("floor")
        for lines, words in ((LEFT[:10], "ends before line 11 of"), (LEFT[:6], "ends before line 7 of"),
                             (LEFT + [""], "runs on at line 12, past the end of"),
                             (LEFT + LEFT_OUT[2:3], "runs on at line 12, past the end of"),
                             (LEFT + LEFT, "runs on at line 12, past the end of")):
            with self.subTest(lines=len(lines), words=words):
                page = self.page(lines)
                self.assertEqual(self.run_main("onpage", record, page), (1, gives(page, words)))
        page = self.page(None, data=b"")
        self.assertEqual(self.run_main("onpage", record, page), (1, gives(page, "ends before line 1 of")))

    def test_onpage_refuses_rows_in_another_order_twice_or_with_one_gone(self):
        _, record = self.made("floor")
        for lines, n in ((LEFT[:8] + [LEFT[9], LEFT[8], LEFT[10]], 9), (LEFT[:9] + [LEFT[8]] + LEFT[9:], 10),
                         (LEFT[:9] + LEFT[10:], 10), (LEFT[:6] + LEFT[6:][::-1], 7),
                         (LEFT[:9] + [LEFT[8], LEFT[9]], 10)):
            with self.subTest(lines=lines[6:]):
                page = self.page(lines)
                self.assertEqual(self.run_main("onpage", record, page), (1, gives(page, f"line {n} differs from")))

    def test_onpage_refuses_the_table_inside_a_code_block_a_comment_or_a_details_block(self):
        _, record = self.made("floor")
        for opens, closes in ((["```"], ["```"]), (["<!--"], ["-->"]), (["<details>", ""], ["", "</details>"])):
            with self.subTest(opens=opens):
                page = self.page(LEFT[:6] + opens + LEFT[6:] + closes)
                self.assertEqual(self.run_main("onpage", record, page), (1, gives(page, "line 7 differs from")))
        page = self.page(LEFT[:6] + ["    " + line for line in LEFT[6:]])
        self.assertEqual(self.run_main("onpage", record, page), (1, gives(page, "line 7 differs from")))

    def test_onpage_with_no_dataset_left_out_reads_its_one_line(self):
        _, record = self.made("floor", "all in")
        lines = LEFT[:4] + [ALL_READS, "", NONE_LEFT]
        page = self.page(lines)
        self.assertEqual(self.run_main("onpage", record, page),
                         (0, f"onpage: {page} is the bytes `m_page.py left` gives\n"))
        page = self.page(LEFT)
        self.assertEqual(self.run_main("onpage", record, page), (1, gives(page, "line 5 differs from")))
        page = self.page(lines[:-1])
        self.assertEqual(self.run_main("onpage", record, page), (1, gives(page, "ends before line 7 of")))

    def test_onpage_runs_the_gates_before_it_opens_the_page(self):
        _, record = self.made("floor")
        m_tables.edit(record, "m-left-out.tsv", lambda t: t.replace("drug 6\tshort list", "drug 6\tunder the floor"))
        why = "m_page: stopped: m-left-out.tsv line 3: synth_c__drug_6: its reason is not its state in m-floor.tsv\n"
        self.assertEqual(self.run_main("onpage", record, self.top / "missing.md"), (1, why))
        page = self.page(LEFT)
        with watched() as seen:
            self.assertEqual(self.run_main("onpage", record, page), (1, why))
        self.assertFalse(any(Path(p) == page for _, p in seen))
        _, none = self.made("none")
        self.assertEqual(self.run_main("onpage", none, page),
                         (1, "m_page: stopped: m-floor.tsv is missing from the record\n"))

    def test_onpage_stops_on_a_page_that_does_not_exist(self):
        _, record = self.made("floor")
        missing = self.top / "missing.md"
        self.assertEqual(self.run_main("onpage", record, missing), (1, f"m_page: stopped: {missing} does not exist\n"))
        folder = self.top / "page.md"
        folder.mkdir()
        self.assertEqual(self.run_main("onpage", record, folder), (1, f"m_page: stopped: {folder} is not a file\n"))


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
        # The row of all datasets enters one more, with one more the same, so its better, worse and same still sum to
        # entered.
        for name, head, cells in (("m-count.tsv", ["M", "second", "lineage", "k31-aa-k31", "gbt", "balanced_accuracy"],
                                   {"datasets": ("3", "4")}),
                                  ("m-count.tsv", ["all datasets", "second", "study", "program-k31", "rf", "recall"],
                                   {"entered": ("11", "12"), "same": None}),
                                  ("m-size.tsv", ["second", "twin", "lasso"], {"datasets": ("11", "12")})):
            with self.subTest(name=name, row=head):
                _, record = self.made()
                at = []

                def change(header, rows, head=head, cells=cells):
                    i = next(i for i, r in enumerate(rows) if [r[c] for c in header[:len(head)]] == head)
                    at.append(i + 2)
                    for c, old_new in cells.items():
                        if old_new is None:
                            rows[i][c] = str(int(rows[i][c]) + 1)
                        else:
                            self.assertEqual(rows[i][c], old_new[0])
                            rows[i][c] = old_new[1]
                    return rows

                change_rows(record, name, change)
                self.stops(record, f"^{name} line {at[0]}: the second reading's .* counts more datasets of M than "
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
                (4, "floor", "26", "its floor is not m_list.py's"),
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
                                       ("m-count.tsv", 2, "entered", ""), ("m-count.tsv", 2, "better", ""),
                                       ("m-count.tsv", 2, "better", "1" * 13)):
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
        # With the row gone, the rows gate names the row the lists give that the table lacks.
        _, record = self.made()
        m_tables.edit(record, "m-count.tsv", lambda t: "\n".join(l for i, l in enumerate(t.split("\n"), 1) if i != n))
        self.stops(record, f"^m-count.tsv lacks its row for first run second {key}$")
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
        for argv in (("print", record), ("left", record), ("onpage", record, page), ("page", record)):
            self.assertEqual(self.run_main(*argv), (1, why), argv[0])

    def test_record_print_and_onpage_each_stop_on_a_gate(self):
        tables, record = self.made("floor")
        page = self.top / "page.md"
        page.write_text("\n".join(LEFT_OUT) + "\n")
        m_tables.edit(record, "m-left-out.tsv", lambda t: t.replace("drug 6\tshort list", "drug 6\tunder the floor"))
        why = "m_page: stopped: m-left-out.tsv line 3: synth_c__drug_6: its reason is not its state in m-floor.tsv\n"
        for argv in (("print", record), ("left", record), ("page", record), ("onpage", record, page)):
            self.assertEqual(self.run_main(*argv), (1, why), argv[0])
        self.assertEqual(self.run_main("record", tables, record)[0], 1)


def where(**cells):
    return lambda r: all(r[c] == v for c, v in cells.items())


def break_cell(pred, column, value, why):
    """A break of m-count.tsv: the first row that meets pred takes value(cell) in column; the gate's line names it."""
    def change(rows):
        i = next(i for i, r in enumerate(rows) if pred(r))
        rows[i][column] = value(rows[i][column])
        return rows, f"m-count.tsv line {i + 2}: {why}"
    return change


def twice(rows):
    i = next(i for i, r in enumerate(rows) if where(count="M", split="lineage", measure="balanced_accuracy")(r))
    return rows + [dict(rows[i])], f"m-count.tsv line {len(rows) + 2}: its row stands a second time"


def one_gone(rows):
    i = next(i for i, r in enumerate(rows)
             if where(count="M", m_reading="second", split="lineage", measure="balanced_accuracy")(r))
    assert (rows[i]["pair"], rows[i]["learner"]) == ("unitigs-k31", "scm")
    return rows[:i] + rows[i + 1:], "m-count.tsv lacks its row for M second lineage unitigs-k31 scm balanced_accuracy"


SUM = "its better, worse and same do not sum to entered"
M_ROW = where(count="M", split="lineage", measure="balanced_accuracy")
# The second reviewer's 12 breaks of m-count.tsv's columns, B1 to B12, on the SYNTHETIC record.
BREAKS = [
    ("B1", break_cell(M_ROW, "better", lambda v: str(int(v) + 1), SUM)),
    ("B2", break_cell(where(count="all datasets", m_reading="first", split="lineage", measure="recall"), "better",
                      lambda v: "99", SUM)),
    ("B3", break_cell(lambda r: r["count"] == "all datasets" and r["holm_p"] != "", "holm_p", lambda v: "1.5",
                      "its holm_p is outside 0 to 1")),
    ("B4", break_cell(lambda r: r["count"] == "all datasets" and r["interval_lo"] != "", "interval_lo",
                      lambda v: "0.9", "its interval's low end is above its high end")),
    ("B5", break_cell(where(count="all datasets"), "median_difference", lambda v: "x",
                      "its median_difference is not a number")),
    ("B6", break_cell(M_ROW, "measure", lambda v: "balanced accuracy", "its measure is not on its list")),
    ("B7", break_cell(where(count="all datasets", split="lineage", measure="recall"), "split", lambda v: "Lineage",
                      "its split is not on its list")),
    ("B8", twice),
    ("B9", break_cell(M_ROW, "m_reading", lambda v: "third", "its m_reading is not on its list")),
    ("B10", break_cell(M_ROW, "learner", lambda v: "zz", "its learner is not on its list")),
    ("B11", break_cell(where(count="M"), "pair", lambda v: "zz-k31", "its pair is not on its list")),
    ("B12", one_gone)]


def break_count(directory, brk, sums=True):
    """m-count.tsv of <directory> through one break; returns the gate's line."""
    why = []

    def change(header, rows):
        rows, line = brk(rows)
        why.append(line)
        return rows

    change_rows(directory, "m-count.tsv", change, sums)
    return why[0]


class CountGates(Made):
    def at(self, record, name, *head):
        """The 1-based line of one table of <record> whose first cells are head, which must be one line."""
        found = [i for i, line in enumerate(lines_of(record / name), 1) if line.split("\t")[:len(head)] == list(head)]
        self.assertEqual(len(found), 1, (name, head))
        return found[0]

    def test_gate_a_figure_cell_that_is_not_a_number(self):
        for name, column, value in (("m-count.tsv", "median_difference", "x"), ("m-count.tsv", "interval_lo", "0.5.1"),
                                    ("m-count.tsv", "sign_p", "nan"), ("m-count.tsv", "holm_p", "1e5"),
                                    ("m-count.tsv", "interval_hi", " 0.5"), ("m-count.tsv", "median_difference", "1,5"),
                                    ("m-count.tsv", "median_difference", "1" * 13),
                                    ("m-datasets.tsv", "arm_value", "none"),
                                    ("m-datasets.tsv", "recall_difference", "--0.1")):
            with self.subTest(name=name, column=column, value=value):
                _, record = self.made()
                set_cell(record, name, 2, column, value)
                with self.assertRaises(m_page.Stop) as caught:
                    m_page.tables(record, M_SET)
                self.assertEqual(str(caught.exception), f"{name} line 2: its {column} is not a number")
        # C9's number forms pass, on a row of the count on M alone that no other gate reads by its median.
        for value in ("1e-05", "-0.0", "1.0", "0.123456789012", "-0.5", "0", ""):
            with self.subTest(value=value):
                _, record = self.made()
                n = self.at(record, "m-count.tsv", "M", "first", "lineage", "unitigs-k31", "scm", "balanced_accuracy")
                set_cell(record, "m-count.tsv", n, "median_difference", value)
                self.assertEqual(len(m_page.tables(record, M_SET)), 8)

    def test_gate_a_name_cell_that_is_not_on_its_list(self):
        for name, column, value in (
                ("m-count.tsv", "count", "all"), ("m-count.tsv", "m_reading", "third"),
                ("m-count.tsv", "split", "Lineage"), ("m-count.tsv", "pair", "zz-k31"), ("m-count.tsv", "learner", "zz"),
                ("m-count.tsv", "measure", "balanced accuracy"), ("m-count.tsv", "orientation", "up"),
                ("m-conditions.tsv", "count", "M alone"), ("m-conditions.tsv", "m_reading", "both"),
                ("m-conditions.tsv", "arm", "twin"), ("m-conditions.tsv", "reading", "ruling 16"),
                ("m-conditions.tsv", "condition", "7"), ("m-both.tsv", "arm", "k31"), ("m-both.tsv", "condition", "6"),
                ("m-size.tsv", "m_reading", "third"), ("m-size.tsv", "arm", "k31"), ("m-size.tsv", "learner", "svm"),
                ("m-datasets.tsv", "m_reading", "third"), ("m-datasets.tsv", "split", "random"),
                ("m-datasets.tsv", "pair", "k31-k31"), ("m-datasets.tsv", "learner", "zz"),
                ("m-datasets.tsv", "second_reading", "out")):
            with self.subTest(name=name, column=column, value=value):
                _, record = self.made()
                set_cell(record, name, 2, column, value)
                with self.assertRaises(m_page.Stop) as caught:
                    m_page.tables(record, M_SET)
                self.assertEqual(str(caught.exception), f"{name} line 2: its {column} is not on its list")

    def test_gate_a_row_twice_a_row_no_list_gives_and_a_row_the_table_lacks(self):
        for name in ("m-count.tsv", "m-conditions.tsv", "m-size.tsv"):
            with self.subTest(name=name, row="twice"):
                _, record = self.made()
                m_tables.edit(record, name, lambda t: t + t.split("\n")[1] + "\n")
                self.stops(record, f"^{name} line {len(lines_of(record / name)) - 1}: its row stands a second time$")
        # The first run's count on specificity, which reads balanced accuracy and recall alone; and ruling 15's
        # condition 1, which it takes with condition 2.
        for name, column, value in (("m-count.tsv", "measure", "specificity"), ("m-conditions.tsv", "condition", "1")):
            with self.subTest(name=name, row="none"):
                _, record = self.made()
                set_cell(record, name, 2, column, value)
                self.stops(record, f"^{name} line 2: its row is none the lists give$")
        for name, key in (("m-count.tsv", "first run first lineage unitigs-k31 scm balanced_accuracy"),
                          ("m-conditions.tsv", "first run first unitigs ruling 15 1 and 2"),
                          ("m-size.tsv", "first unitigs scm")):
            with self.subTest(name=name, row="gone"):
                _, record = self.made()
                m_tables.edit(record, name, lambda t: "\n".join(l for i, l in enumerate(t.split("\n"), 1) if i != 2))
                self.stops(record, f"^{name} lacks its row for {key}$")

    def test_gate_better_worse_and_same_that_do_not_sum_to_entered(self):
        n = None
        for cells in ({"better": lambda v: str(int(v) + 1)}, {"same": lambda v: str(int(v) + 1)},
                      {"entered": lambda v: str(int(v) - 1)}, {"worse": lambda v: str(int(v) + 1),
                                                               "better": lambda v: str(int(v) + 1)}):
            with self.subTest(cells=sorted(cells)):
                _, record = self.made()
                n = self.at(record, "m-count.tsv", "M", "first", "lineage", "unitigs-k31", "scm", "balanced_accuracy")
                for column, value in cells.items():
                    old = lines_of(record / "m-count.tsv")[n - 1].split("\t")[m_page.COLUMNS["m-count.tsv"].index(column)]
                    set_cell(record, "m-count.tsv", n, column, value(old))
                self.stops(record, f"^m-count.tsv line {n}: its better, worse and same do not sum to entered$")
        # One better fewer and one more the same still sum, and pass.
        _, record = self.made()
        lines = lines_of(record / "m-count.tsv")
        cells = dict(zip(lines[0].split("\t"), lines[n - 1].split("\t")))
        set_cell(record, "m-count.tsv", n, "better", str(int(cells["better"]) - 1))
        set_cell(record, "m-count.tsv", n, "same", str(int(cells["same"]) + 1))
        self.assertEqual(len(m_page.tables(record, M_SET)), 8)

    def test_gate_a_p_value_outside_0_to_1(self):
        for column, value in (("sign_p", "1.5"), ("holm_p", "1.0000001"), ("sign_p", "-0.01"), ("holm_p", "-0.0001")):
            with self.subTest(column=column, value=value):
                _, record = self.made()
                n = self.at(record, "m-count.tsv", "all datasets", "first", "lineage", "unitigs-k31", "scm",
                            "balanced_accuracy")
                set_cell(record, "m-count.tsv", n, column, value)
                self.stops(record, f"^m-count.tsv line {n}: its {column} is outside 0 to 1$")
        for value in ("0", "1", "1.0", "0.0"):
            with self.subTest(value=value):
                _, record = self.made()
                n = self.at(record, "m-count.tsv", "M", "first", "lineage", "unitigs-k31", "scm", "balanced_accuracy")
                set_cell(record, "m-count.tsv", n, "sign_p", value)
                self.assertEqual(len(m_page.tables(record, M_SET)), 8)

    def test_gate_an_interval_with_its_low_end_above_its_high_end_or_one_end_alone(self):
        for column, value, why in (("interval_lo", "0.9", "its interval's low end is above its high end"),
                                   ("interval_hi", "-0.9", "its interval's low end is above its high end"),
                                   ("interval_lo", "", "its interval has one end alone"),
                                   ("interval_hi", "", "its interval has one end alone")):
            with self.subTest(column=column, value=value):
                _, record = self.made()
                n = self.at(record, "m-count.tsv", "M", "first", "lineage", "unitigs-k31", "scm", "balanced_accuracy")
                set_cell(record, "m-count.tsv", n, column, value)
                self.stops(record, f"^m-count.tsv line {n}: {why}$")
        # A low end equal to its high end passes.
        _, record = self.made()
        n = self.at(record, "m-count.tsv", "M", "first", "lineage", "unitigs-k31", "scm", "balanced_accuracy")
        lines = lines_of(record / "m-count.tsv")
        set_cell(record, "m-count.tsv", n, "interval_lo",
                 dict(zip(lines[0].split("\t"), lines[n - 1].split("\t")))["interval_hi"])
        self.assertEqual(len(m_page.tables(record, M_SET)), 8)

    def test_gate_a_size_row_unlike_its_row_of_the_count_over_all_datasets(self):
        for column, value in (("better", lambda v: str(int(v) - 1)), ("datasets", lambda v: str(int(v) + 1)),
                              ("same", lambda v: str(int(v) + 1))):
            with self.subTest(column=column):
                _, record = self.made()
                old = lines_of(record / "m-size.tsv")[1].split("\t")[m_page.COLUMNS["m-size.tsv"].index(column)]
                set_cell(record, "m-size.tsv", 2, column, value(old))
                self.stops(record, "^m-size.tsv line 2: its row is not its row of m-count.tsv$")

    def test_gate_a_count_of_learners_unlike_the_count_of_the_marks(self):
        # Condition 1 alone and condition 2 alone by the stricter reading, both by ruling 15, and the larger of the
        # twin's two counts for ruling 15's condition 3, for the count over all datasets and the count on M alone.
        for count, reading_, condition in (("all datasets", "stricter", "1"), ("M", "stricter", "2"),
                                           ("all datasets", "ruling 15", "1 and 2"), ("M", "ruling 15", "3")):
            with self.subTest(count=count, reading=reading_, condition=condition):
                _, record = self.made()
                n = self.at(record, "m-conditions.tsv", count, "first", "unitigs", reading_, condition)
                old = lines_of(record / "m-conditions.tsv")[n - 1].split("\t")[-1]
                set_cell(record, "m-conditions.tsv", n, "learners", str((int(old) + 1) % 6))
                self.stops(record, f"^m-conditions.tsv line {n}: its learners is not the count of the marks "
                                   "m-count.tsv gives$")
        # One count row's Holm p made 0.5: its condition 2 mark goes to 0, and ruling 15's count of the learners that
        # meet both conditions no longer sums its marks.
        _, record = self.made()
        head = ["all datasets", "first", "lineage", "unitigs-k31", "scm", "balanced_accuracy"]

        def change(header, rows):
            r = next(r for r in rows if [r[c] for c in header[:6]] == head)
            self.assertEqual(marked(r), ["1", "1"])
            r["holm_p"] = "0.5"
            return rows

        change_rows(record, "m-count.tsv", change)
        n = self.at(record, "m-conditions.tsv", "all datasets", "first", "unitigs", "ruling 15", "1 and 2")
        self.stops(record, f"^m-conditions.tsv line {n}: its learners is not the count of the marks m-count.tsv gives$")


class Breaks(Made):
    def test_twelve_column_breaks_stop_record_print_left_page_and_onpage_with_one_gate_line(self):
        page = self.top / "page.md"
        page.write_bytes(text(LEFT))
        for label, brk in BREAKS:
            with self.subTest(label):
                tables, record = self.made("floor")
                why = break_count(tables, brk, sums=False)
                self.assertEqual(self.run_main("record", tables, record), (1, f"m_page: stopped: {why}\n"))
                self.assertEqual(m_page.state(record), "floor")
                _, record = self.made()
                self.assertEqual(break_count(record, brk), why)
                for argv in (("print", record), ("left", record), ("page", record), ("onpage", record, page)):
                    self.assertEqual(self.run_main(*argv), (1, f"m_page: stopped: {why}\n"), argv[0])

    def test_the_three_breaks_that_ended_in_a_traceback_give_one_line_from_the_command_line(self):
        page = self.top / "page.md"
        page.write_bytes(text(LEFT))
        for label, brk in BREAKS:
            if label not in ("B5", "B6", "B11"):
                continue
            with self.subTest(label):
                tables, record = self.made("floor")
                why = break_count(tables, brk, sums=False)
                line = f"m_page: stopped: {why}\n".encode("utf-8")
                self.assertEqual(self.cli("record", tables, record), (1, line, b""))
                _, record = self.made()
                break_count(record, brk)
                for argv in (("print", record), ("left", record), ("page", record), ("onpage", record, page)):
                    self.assertEqual(self.cli(*argv), (1, line, b""), argv[0])

    def test_a_record_or_tables_directory_that_is_not_a_directory_stops_with_one_line(self):
        nope = self.top / "nope"
        page = self.top / "page.md"
        page.write_bytes(text(LEFT))
        for argv in (("print", nope), ("left", nope), ("page", nope), ("onpage", nope, page)):
            self.assertEqual(self.run_main(*argv), (1, f"m_page: stopped: the record {nope} is not a directory\n"),
                             argv[0])
        tables, record = self.made("floor")
        self.assertEqual(self.run_main("record", nope, record),
                         (1, f"m_page: stopped: the tables directory {nope} is not a directory\n"))
        self.assertEqual(self.run_main("record", tables, page),
                         (1, f"m_page: stopped: the record {page} is not a directory\n"))
        self.assertEqual(self.run_main("floor", tables, nope),
                         (1, f"m_page: stopped: the record {nope} is not a directory\n"))

    def test_an_error_no_gate_names_stops_with_its_kind_alone_and_no_traceback(self):
        def broken(*args):
            raise KeyError("synth_c__drug_5")

        page = self.top / "page.md"
        page.write_bytes(text(LEFT))
        tables, record = self.made("floor")
        _, full = self.made()
        _, bare = self.made("none")
        why = (1, "m_page: stopped: m_page.py met a KeyError that no gate names\n")
        with mock.patch.object(m_page, "gate_tables", broken):
            for argv in (("floor", tables, bare), ("record", tables, record), ("print", full), ("left", full),
                         ("page", full), ("onpage", full, page)):
                self.assertEqual(self.run_main(*argv), why, argv[0])
        self.assertEqual(sorted(p.name for p in bare.iterdir()), ["README.md", "first-run-times.tsv"])


class States(Made):
    def test_the_headers_of_this_run_are_eleven_and_are_the_headers_print_and_page_give(self):
        self.assertEqual(m_page.headers(), [
            ["Species", "Antibiotic", "Second reading", "Cohort", "Listed", "Left after the list", "Floor"],
            ["Species", "Datasets of M", "In the second reading", "Left out, under the floor", "Left out, short list"],
            ["Count"] + FULL,
            ["Count", "Reading", "Split", "Arm", "Learner", "Measure", "Better is", "Entered", "Better", "Worse", "Same",
             "Median", "Interval"],
            ["Count", "Reading", "Arm", "Condition", "Read", "Learners"],
            ["Arm", "Condition", "First reading", "Second reading", "Both"],
            ["Reading", "Arm", "Learner", "Entered", "Better", "Worse", "Same"],
            ["Species", "Antibiotic", "Arm", "Learner", "First reading: arm", "First reading: K31",
             "First reading: difference", "Second reading: arm", "Second reading: K31", "Second reading: difference"],
            FULL + ["Condition 1", "Condition 2"], FULL, FULL[:8]])
        self.assertEqual(len({tuple(h) for h in m_page.headers()}), 11)
        _, record = self.made()
        given = set()
        for lines in list(m_page.tables(record, M_SET).values()) + [
                m_page.count_page(record, M_SET).decode("utf-8").split("\n")]:
            for i, line in enumerate(lines):
                if re.fullmatch(r"\|(---\|)+", line):
                    given.add(tuple(lines[i - 1][2:-2].split(" | ")))
        self.assertEqual(given, {tuple(h) for h in m_page.headers()})
