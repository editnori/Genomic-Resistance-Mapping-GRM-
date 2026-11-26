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
import tempfile
import unittest
from pathlib import Path

from tests.c3_inputs import common, m_tables
from tests.c3_inputs import test_m_reading as mr

m_page = m_tables.m_page
reading = mr.reading
page_tables = common.load("page_tables")
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
BUILT = {}


def setUpModule():
    BUILT["tmp"] = tempfile.TemporaryDirectory()
    top = Path(BUILT["tmp"].name)
    BUILT["seven"] = m_tables.seven(top / "default")
    BUILT["all in"] = m_tables.seven(top / "all-in", m_tables.ALL_IN)


def tearDownModule():
    BUILT.pop("tmp").cleanup()


def lines_of(path):
    return Path(path).read_text(encoding="utf-8").split("\n")


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

    def test_the_words_and_number_forms_are_page_tables_own(self):
        self.assertEqual((m_page.ARM_WORDS, m_page.MEASURE_WORDS), (page_tables.ARM_WORDS, page_tables.MEASURE_WORDS))
        for x in ("", None, "0", "0.0005", "-0.0005", "0.0015", "0.0025", "0.123456789012", "-0.02", "1", "0.9995",
                  "1e-05", 0.0004, -0.1234):
            self.assertEqual(m_page.d3(x), page_tables.d3(x), x)
            self.assertEqual(m_page.pval(x), page_tables.pval(x), x)
            self.assertEqual(m_page.interval(x, "0.5"), page_tables.interval(x, "0.5"), x)
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
        # Each count table is a directory, so opening one fails; the directory can be passed and not listed.
        tables, record = self.made("none")
        for name in SEVEN[2:]:
            (tables / name).unlink()
            (tables / name).mkdir()
        os.chmod(tables, 0o300)
        try:
            self.assertEqual(self.run_main("floor", tables, record)[0], 0)
        finally:
            os.chmod(tables, 0o700)
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
        # The same count changed in both floor tables, so every other gate passes and only the bytes differ.
        tables, record = self.made("floor")
        for name in ("m-floor.tsv", "m-left-out.tsv"):
            text = (tables / name).read_text()
            self.assertEqual(text.count("\t36\t40\t25"), 1, name)
            (tables / name).write_text(text.replace("\t36\t40\t25", "\t36\t41\t25"))
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
        _, record = self.made()
        m_tables.edit(record, "m-floor.tsv", lambda t: t.replace("\t1\tin\n", "\t1\tshown\n", 1))
        self.stops(record, "^m-floor.tsv line 4: synth_c__drug_7: its second reading is none of in, under the floor, "
                           "short list$")

    def test_gate_a_token_of_a_genome_ids_shape(self):
        # pack.py's shape with 2 or more digits before the point stops, naming the line and not the token; one digit
        # before it, or the token after a letter, passes. Each made-up token starts 00, as no taxon id does, so no real
        # genome id can equal it.
        for token, stops in (("00.000", True), ("0.5", False), ("x00.000", False), ("_00.0", True),
                             ("0.123456789012", False), ("000.00.0", True)):
            with self.subTest(token=token):
                _, record = self.made()
                m_tables.edit(record, "m-conditions.tsv", lambda t: t.replace("after Holm over 30", f"after {token}", 1))
                if stops:
                    with self.assertRaises(m_page.Stop) as caught:
                        m_page.tables(record, M_SET)
                    self.assertRegex(str(caught.exception), "^m-conditions.tsv line 2: a token of a genome id's shape, "
                                                            "digits before the point$")
                    self.assertNotIn(token.strip("x_"), str(caught.exception))
                else:
                    self.assertEqual(len(m_page.tables(record, M_SET)), 8)

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
                         ["| Count", "| Count", "| Count", "| Arm", "| Reading", "| Species"])


if __name__ == "__main__":
    unittest.main()
