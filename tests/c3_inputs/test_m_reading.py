"""scripts/c3-inputs/m_reading.py: the three counts in both readings, and the page tables of the run on M.

The check of the first run runs for real on the committed record,
reference/input-arms/ (page-datasets.tsv, page-count.tsv,
page-conditions.tsv, SHA256SUMS). The counts and tables run on SYNTHETIC
datasets.tsv rows: made-up species, datasets and figures, no genome.
"""

import contextlib
import io
import shutil
import tempfile
import unittest
from pathlib import Path

from tests.c3_inputs import common, m_common

reading = common.load("m_reading")
score = reading.SCORE
REFERENCE = common.ROOT / "reference" / "input-arms"
FAMILY = score.FAMILY
ARMS = [score.DEFAULT] + FAMILY
LEARNERS = score.LEARNERS
MEASURES = score.MEASURES
RECORD = [f"synth_{s}__drug_{i}" for s in "ab" for i in range(1, 5)]
M_SETS = {"synth_c": [f"synth_c__drug_{i}" for i in (5, 6, 7)], "synth_d": [f"synth_d__drug_{i}" for i in (8, 9, 10)]}
UNDER = "synth_d__drug_10"
# Listed genomes by class (R, S): 1 genome, 4 and 0, 5 and 0, none, 6 and 7, and the dataset under the floor.
LISTED = {"synth_c__drug_5": (1, 0), "synth_c__drug_6": (4, 0), "synth_c__drug_7": (5, 0), "synth_d__drug_8": (0, 0),
          "synth_d__drug_9": (6, 7), UNDER: (6, 2)}
NOT_SHOWN = {"synth_c__drug_5", "synth_c__drug_6"}
SHOWN = {"synth_c__drug_7", "synth_d__drug_8", "synth_d__drug_9"}


def figure(i, arm, learner, measure, shift=0):
    """A made-up figure in [0, 1]: the default's, or the default's plus a difference. The difference is one
    value across datasets except for the first learner's balanced accuracy in unitigs and the twin, and the
    second reading's shift moves that learner's balanced accuracy alone, so few rows draw resamples."""
    base = 0.70 + 0.01 * (i % 5) + 0.001 * MEASURES.index(measure)
    if arm == score.DEFAULT:
        return base
    j, k = ARMS.index(arm), LEARNERS.index(learner)
    lift = 0.05 if arm == "unitigs" else 0.0
    first = k == 0 and measure == score.PRIMARY
    spread = ((i * 7) % 5 - 2) * 0.01 if first and arm in ("unitigs", score.TWIN) else 0.0
    return base + lift + spread + ((j * 3 + k * 5 + (shift if first else 0)) % 11 - 3) * 0.01


def rows_of(datasets, measures, start, shift=0):
    """datasets.tsv rows as text, as read back from a file: every arm, both splits (the twin on lineage)."""
    out = []
    for n, d in enumerate(datasets):
        for split in ("lineage", "study"):
            for arm in ARMS:
                if split == "study" and arm == score.TWIN:
                    continue
                for learner in LEARNERS:
                    r = {"dataset": d, "species_slug": d.split("__")[0], "species": "Synth " + d.split("__")[0][6:],
                         "antibiotic": d.split("__")[1].replace("_", " "), "split": split, "arm": arm,
                         "learner": learner, "draws": "3" if split == "lineage" else "1",
                         "scored_draws": "3" if split == "lineage" else "1", "constant_draws": "0", "scored": "1",
                         "n": "80", "tp": "30", "fp": "10", "tn": "30", "fn": "10"}
                    for m in measures:
                        r[m] = repr(round(figure(start + n, arm, learner, m, shift), 12))
                    out.append(r)
    return out


def floor_rows(slug, listed=LISTED):
    out = []
    for d in M_SETS[slug]:
        lr, ls = listed[d]
        r, s = (30, 40) if d == UNDER else (40, 40)
        out.append({"dataset": d, "species_slug": slug, "species": "Synth " + slug[6:],
                    "antibiotic": d.split("__")[1].replace("_", " "), "cohort": r + s, "resistant": r,
                    "susceptible": s, "listed": lr + ls, "listed_resistant": lr, "listed_susceptible": ls,
                    "second_resistant": r - lr, "second_susceptible": s - ls, "floor": 25,
                    "enters": int(min(r - lr, s - ls) >= 25)})
    return out


def folders(top, listed=LISTED, scores=True):
    """One species folder of the run on M per synthetic species, as its pack holds it. With scores False, only
    floor.tsv: the floor form reads nothing else."""
    out = []
    start = len(RECORD)
    for slug, datasets in M_SETS.items():
        f = Path(top) / slug
        if scores:
            first = rows_of(datasets, MEASURES, start)
            second = rows_of([d for d in datasets if d != UNDER], MEASURES, start, shift=4)
            m_common.write_tsv(f / "score" / "datasets.tsv", score.DATASET_COLUMNS, first)
            m_common.write_tsv(f / "score" / "second" / "datasets.tsv", score.DATASET_COLUMNS, second)
        m_common.write_tsv(f / "results" / "m" / "floor.tsv", reading.LIST.FLOOR_COLUMNS, floor_rows(slug, listed))
        start += len(datasets)
        out.append(f)
    return out


def tables(top, fs, record, listed=LISTED):
    """The floor form's tables, then the count's, written to <top>/tables. Returns (taken, the count rows)."""
    out = Path(top) / "tables"
    m_rows, floor = reading.species_rows(fs, LISTS)
    reading.write_floor(out, floor)
    taken = reading.counts(record, m_rows)
    return taken, reading.write(out, taken, floor, m_rows, record)


LISTS = [{"dataset": d, "in_M": "1"} for ds in M_SETS.values() for d in ds] + [
    {"dataset": d, "in_M": "0"} for d in RECORD]


class FirstRun(unittest.TestCase):
    """The first run's count and conditions, given back from the committed record."""

    @classmethod
    def setUpClass(cls):
        cls.rows, cls.counted, cls.conds = reading.record_count(REFERENCE)

    def copy(self, d, name, edit):
        ref = Path(d) / "ref"
        shutil.copytree(REFERENCE, ref)
        lines = (ref / name).read_text().splitlines(keepends=True)
        lines[1] = edit(lines[1])
        (ref / name).write_text("".join(lines))
        return ref

    def test_the_committed_record_gives_back_110_count_rows_and_55_condition_rows(self):
        self.assertEqual(reading.given_back(REFERENCE, self.counted, self.conds), (110, 55))
        self.assertEqual({r["arm"] for r in self.rows}, set(ARMS))
        self.assertEqual({r["measure"] for r in self.counted}, set(reading.RECORD_MEASURES))

    def test_a_published_count_cell_unlike_the_records_count_stops(self):
        with tempfile.TemporaryDirectory() as d:
            ref = self.copy(d, "page-count.tsv", lambda l: l.replace("\t19\t29\t184\t", "\t20\t28\t184\t", 1))
            with self.assertRaisesRegex(reading.Stop, "count row .*better"):
                reading.given_back(ref, self.counted, self.conds)

    def test_a_published_condition_unlike_the_records_conditions_stops(self):
        with tempfile.TemporaryDirectory() as d:
            ref = self.copy(d, "page-conditions.tsv", lambda l: l[: l.rindex("\t0\t0\n")] + "\t1\t0\n")
            with self.assertRaisesRegex(reading.Stop, "condition .*met"):
                reading.given_back(ref, self.counted, self.conds)

    def test_a_record_file_unlike_its_sha256sums_stops_before_any_count(self):
        for name in ("page-datasets.tsv", "page-count.tsv", "page-conditions.tsv"):
            with tempfile.TemporaryDirectory() as d:
                ref = self.copy(d, name, lambda l: l.replace("\t", "\t ", 1))
                with self.assertRaisesRegex(reading.Stop, f"{name} differs from SHA256SUMS"):
                    reading.check(ref)
                with contextlib.redirect_stdout(io.StringIO()) as out:
                    self.assertEqual(reading.main(["count", str(ref), str(Path(d) / "out"), str(Path(d) / "none")]), 1)
                self.assertIn("stopped", out.getvalue())


class Counts(unittest.TestCase):
    """Three counts in each reading on synthetic rows, each its own call of score.py's count."""

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.record = rows_of(RECORD, reading.RECORD_MEASURES, 0)
        cls.m_rows, cls.floor = reading.species_rows(folders(cls.tmp.name), LISTS)
        cls.taken = reading.counts(cls.record, cls.m_rows)
        cls.out = Path(cls.tmp.name) / "tables"
        reading.write_floor(cls.out, cls.floor)
        reading.write(cls.out, cls.taken, cls.floor, cls.m_rows, cls.record)

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def table(self, name):
        return score.read_tsv(self.out / name)

    def test_three_counts_in_each_of_two_readings(self):
        self.assertEqual(sorted(self.taken), sorted((c, g) for c in reading.COUNTS for g in reading.READINGS))
        for (c, g), (rows, conds) in self.taken.items():
            want = MEASURES if c == "M" else list(reading.RECORD_MEASURES)
            self.assertEqual(len(rows), 11 * len(LEARNERS) * len(want), (c, g))
            self.assertEqual({r["measure"] for r in rows}, set(want), (c, g))
            self.assertEqual(len(conds), 5 * 11, (c, g))
        lineage = lambda c, g: {r["datasets"] for r in self.taken[(c, g)][0] if r["split"] == "lineage"}
        self.assertEqual((lineage("first run", "first"), lineage("M", "first"), lineage("all datasets", "first")),
                         ({8}, {6}, {14}))
        self.assertEqual((lineage("M", "second"), lineage("all datasets", "second")), ({5}, {13}))
        self.assertEqual(self.taken[("first run", "first")], self.taken[("first run", "second")])

    def test_holm_is_taken_apart_for_each_count_and_each_reading(self):
        fam = lambda c, g: [r for r in self.taken[(c, g)][0] if r["split"] == "lineage" and r["measure"] == score.PRIMARY]
        for key in self.taken:
            rows = fam(*key)
            self.assertEqual(len(rows), 30, key)
            self.assertEqual([r["holm_p"] for r in rows], score.C9.holm([r["sign_p"] for r in rows]), key)
        a, b = fam("all datasets", "first"), fam("all datasets", "second")
        across = score.C9.holm([r["sign_p"] for r in a + b])
        self.assertNotEqual([r["holm_p"] for r in a], across[:30])

    def test_the_m_count_reads_six_measures_and_the_record_counts_two(self):
        m = {r["measure"] for r in self.table("m-count.tsv") if r["count"] == "M"}
        rest = {r["measure"] for r in self.table("m-count.tsv") if r["count"] != "M"}
        self.assertEqual((m, rest), (set(MEASURES), set(reading.RECORD_MEASURES)))

    def test_the_left_out_table_names_each_dataset_under_the_floor_with_its_counts(self):
        rows = self.table("m-left-out.tsv")
        self.assertEqual([r["dataset"] for r in rows], [UNDER])
        r = rows[0]
        self.assertEqual((r["species"], r["antibiotic"], r["resistant"], r["susceptible"], r["listed_resistant"],
                          r["listed_susceptible"], r["second_resistant"], r["second_susceptible"], r["floor"]),
                         ("Synth d", "drug 10", "30", "40", "6", "2", "24", "38", "25"))
        self.assertEqual(len(self.table("m-floor.tsv")), 6)

    def test_the_size_table_gives_entered_better_and_worse_per_reading_arm_and_learner(self):
        rows = self.table("m-size.tsv")
        self.assertEqual(len(rows), 2 * len(FAMILY) * len(LEARNERS))
        self.assertEqual({r["arm"] for r in rows}, set(FAMILY))
        for r in rows:
            c = next(x for x in self.taken[("all datasets", r[reading.READ])][0] if x["split"] == "lineage"
                     and x["measure"] == score.PRIMARY and x["pair"] == f"{r['arm']}-k31" and x["learner"] == r["learner"])
            self.assertEqual((r["entered"], r["better"], r["worse"], r["same"]),
                             (str(c["entered"]), str(c["better"]), str(c["worse"]), str(c["same"])))
            self.assertEqual(int(r["better"]) + int(r["worse"]) + int(r["same"]), int(r["entered"]))
        self.assertNotIn("pick", " ".join(reading.SIZE_COLUMNS))

    def test_the_both_table_meets_a_condition_only_when_both_readings_meet_it(self):
        rows = self.table("m-both.tsv")
        self.assertEqual(len(rows), 4 * 5)
        for r in rows:
            self.assertEqual(r["met_both"], str(int(r["met_first"] == "1" and r["met_second"] == "1")))

    def test_the_dataset_table_holds_the_first_reading_of_each_and_the_second_of_those_shown(self):
        rows = self.table("m-datasets.tsv")
        first = {r["dataset"] for r in rows if r[reading.READ] == "first"}
        second = {r["dataset"] for r in rows if r[reading.READ] == "second"}
        self.assertEqual((len(first), second), (6, SHOWN))
        self.assertEqual(list(rows[0]), reading.M_DATASET_COLUMNS)
        state = {r["dataset"]: r["second_reading"] for r in rows}
        self.assertEqual(state, {**{d: "shown" for d in SHOWN}, **{d: "counted and not shown" for d in NOT_SHOWN},
                                 UNDER: "under the floor"})
        for r in rows:
            self.assertEqual((int(r["listed_resistant"]), int(r["listed_susceptible"])), LISTED[r["dataset"]])

    def test_a_dataset_listing_1_or_4_and_0_is_counted_and_its_second_reading_not_shown(self):
        # In each count of the second reading the datasets of M not shown are counted: 5 enter, 3 of them shown.
        self.assertEqual({r["datasets"] for r in self.taken[("M", "second")][0] if r["split"] == "lineage"}, {5})
        counted = self.table("m-count.tsv")
        second = [r for r in counted if r[reading.READ] == "second" and r["count"] != "first run"]
        self.assertEqual({r["datasets_not_shown"] for r in second if r["split"] == "lineage"}, {"2"})
        self.assertEqual({r["figures"] for r in counted}, {"shown"})
        floor = {r["dataset"]: r["second_reading"] for r in self.table("m-floor.tsv")}
        self.assertEqual({d for d, s in floor.items() if s == "counted and not shown"}, NOT_SHOWN)

    def test_no_file_holds_a_second_reading_figure_of_a_dataset_not_shown(self):
        # Each dataset's second-reading rows carry a figure no other row has; a dataset not shown leaves none behind.
        with tempfile.TemporaryDirectory() as d:
            fs = folders(d)
            marks = {}
            for f in fs:
                path = f / "score" / "second" / "datasets.tsv"
                rows = score.read_tsv(path)
                for r in rows:
                    if r["arm"] != score.DEFAULT and r["learner"] == LEARNERS[1]:
                        n = sorted(LISTED).index(r["dataset"])
                        r["balanced_accuracy"] = repr(0.512345678901 + n * 1e-9)
                        marks.setdefault(r["dataset"], set()).add(r["balanced_accuracy"][:12])
                m_common.write_tsv(path, score.DATASET_COLUMNS, rows)
            tables(d, fs, self.record)
            text = "".join(p.read_text() for p in sorted((Path(d) / "tables").iterdir()))
            for ds, want in ((SHOWN, True), (NOT_SHOWN, False)):
                for name in ds:
                    self.assertEqual(any(m in text for m in marks[name]), want, name)

    def test_an_m_dataset_in_the_record_stops(self):
        clash = self.record + [dict(r) for r in self.m_rows["first"] if r["dataset"] == "synth_c__drug_5"]
        with self.assertRaisesRegex(reading.Stop, "in the first run's record"):
            reading.counts(clash, self.m_rows)


class SmallLists(unittest.TestCase):
    def test_each_class_lists_none_or_at_least_5_to_be_shown(self):
        row = lambda r, s: {"listed_resistant": str(r), "listed_susceptible": str(s), "enters": "1"}
        for (r, s), want in (((1, 0), "counted and not shown"), ((4, 0), "counted and not shown"),
                             ((0, 4), "counted and not shown"), ((5, 0), "shown"), ((0, 0), "shown"),
                             ((5, 1), "counted and not shown"), ((5, 9), "shown")):
            self.assertEqual(reading.second_state(row(r, s)), want, (r, s))
        self.assertEqual(reading.second_state(dict(row(0, 0), enters="0")), "under the floor")
        self.assertEqual(reading.SMALLEST_LISTED, 5)

    def test_a_count_row_whose_one_figure_off_the_page_is_a_dataset_not_shown_gives_no_figures(self):
        # With 4 and 0 listing none, the dataset listing 1 genome is the one value of balanced accuracy and recall
        # on no page in each count row of the second reading: those rows give counts of datasets alone.
        listed = dict(LISTED, synth_c__drug_6=(0, 0))
        record = rows_of(RECORD, reading.RECORD_MEASURES, 0)
        with tempfile.TemporaryDirectory() as d:
            taken, counted = tables(d, folders(d, listed), record, listed)
            out = Path(d) / "tables"
            for r in counted:
                off = r[reading.READ] == "second" and r["count"] != "first run" \
                    and r["measure"] in reading.RECORD_MEASURES
                self.assertEqual(r["figures"], "counted and not shown" if off else "shown", r)
                if off:
                    self.assertEqual({r[k] for k in reading.WITHHELD}, {None})
                    self.assertEqual((r["datasets_not_shown"], r["entered"] > 0), (1, True))
            size = score.read_tsv(out / "m-size.tsv")
            self.assertEqual({(r[reading.READ], r["figures"], r["better"]) for r in size if r[reading.READ] == "second"},
                             {("second", "counted and not shown", "")})
            self.assertTrue(all(r["better"] != "" for r in size if r[reading.READ] == "first"))
            # The count itself is taken on every row: the conditions read the dataset not shown.
            self.assertEqual(len(score.read_tsv(out / "m-conditions.tsv")), 6 * 5 * 11)


class Floor(unittest.TestCase):
    def test_the_floor_form_writes_the_left_out_and_floor_tables_from_floor_tsv_alone(self):
        with tempfile.TemporaryDirectory() as d:
            fs = folders(d, scores=False)
            out = Path(d) / "tables"
            said = reading.floor_form(out, fs, LISTS)
            self.assertEqual(sorted(p.name for p in out.iterdir()), ["m-floor.tsv", "m-left-out.tsv"])
            self.assertEqual([r["dataset"] for r in score.read_tsv(out / "m-left-out.tsv")], [UNDER])
            self.assertEqual(len(score.read_tsv(out / "m-floor.tsv")), 6)
            self.assertEqual(said, "floor: datasets of M 6; second reading: 3 shown; 2 counted and not shown; "
                                   "1 under the floor")
            self.assertFalse(any((f / "score").exists() for f in fs))
            # The tables the count would write from the full folders are these, byte for byte.
            full = folders(Path(d) / "full")
            reading.same_floor(out, reading.species_rows(full, LISTS)[1])
            with contextlib.redirect_stdout(io.StringIO()) as printed:
                self.assertEqual(reading.main(["floor", str(out), str(Path(d) / "none")]), 1)
            self.assertIn("m_reading: stopped: none: no results/m/floor.tsv", printed.getvalue())

    def test_the_count_stops_unless_the_floor_forms_tables_are_there_and_the_same(self):
        with tempfile.TemporaryDirectory() as d:
            fs = folders(d)
            m_rows, floor = reading.species_rows(fs, LISTS)
            out = Path(d) / "tables"
            with self.assertRaisesRegex(reading.Stop, "m-floor.tsv is not in .*the floor form runs first"):
                reading.same_floor(out, floor)
            reading.write_floor(out, floor)
            reading.same_floor(out, floor)
            moved = [dict(r, enters="1") if r["dataset"] == UNDER else r for r in floor]
            with self.assertRaisesRegex(reading.Stop, "m-floor.tsv differs from the one the floor form writes"):
                reading.same_floor(out, moved)
            left = out / "m-left-out.tsv"
            left.write_text(left.read_text().replace("drug 10", "drug 11"))
            with self.assertRaisesRegex(reading.Stop, "m-left-out.tsv differs"):
                reading.same_floor(out, floor)

    def test_the_count_form_reads_no_score_until_the_floor_forms_tables_are_held_to_floor_tsv(self):
        # With no floor form run, the count stops on the floor tables before a score file is read: here the first
        # species' score file is no TSV at all. After the floor form, the count runs.
        record = rows_of(RECORD, reading.RECORD_MEASURES, 0)
        with tempfile.TemporaryDirectory() as d:
            fs = folders(d)
            out = Path(d) / "tables"
            kept = (fs[0] / "score" / "datasets.tsv").read_bytes()
            (fs[0] / "score" / "datasets.tsv").write_bytes(b"\xff not a table\n")
            with self.assertRaisesRegex(reading.Stop, "m-floor.tsv is not in .*the floor form runs first"):
                list(reading.count_form(record, out, fs, LISTS))
            (fs[0] / "score" / "datasets.tsv").write_bytes(kept)
            reading.floor_form(out, fs, LISTS)
            said = list(reading.count_form(record, out, fs, LISTS))
            self.assertEqual(said[0], "floor: m-floor.tsv and m-left-out.tsv are the floor form's")
            self.assertTrue(said[1].startswith("count: datasets of M 6; second reading: 3 shown; "
                                               "2 counted and not shown; 1 under the floor; "), said[1])


class SpeciesFolders(unittest.TestCase):
    def test_each_dataset_of_m_is_read_once_and_nothing_else(self):
        with tempfile.TemporaryDirectory() as d:
            fs = folders(d)
            m_rows, floor = reading.species_rows(fs, LISTS)
            self.assertEqual(len(floor), 6)
            with self.assertRaisesRegex(reading.Stop, "5 of the 5 datasets of M and 1 datasets outside M"):
                reading.species_rows(fs, LISTS[1:] + [{"dataset": "synth_c__drug_5", "in_M": "0"}])
            with self.assertRaisesRegex(reading.Stop, "read twice"):
                reading.species_rows(fs + fs[:1], LISTS)
            path = fs[0] / "results" / "m" / "floor.tsv"
            rows = score.read_tsv(path)
            m_common.write_tsv(path, reading.LIST.FLOOR_COLUMNS, rows + rows[:1])
            with self.assertRaisesRegex(reading.Stop, "synth_c: a dataset of M is read twice"):
                reading.species_rows(fs, LISTS)
            m_common.write_tsv(path, reading.LIST.FLOOR_COLUMNS, rows)
            with self.assertRaisesRegex(reading.Stop, "3 of the 6"):
                reading.species_rows(fs[:1], LISTS)

    def test_a_second_reading_holding_a_dataset_under_the_floor_stops(self):
        with tempfile.TemporaryDirectory() as d:
            fs = folders(d)
            path = fs[1] / "score" / "second" / "datasets.tsv"
            rows = score.read_tsv(path) + rows_of([UNDER], MEASURES, 20)
            m_common.write_tsv(path, score.DATASET_COLUMNS, rows)
            with self.assertRaisesRegex(reading.Stop, "not those floor.tsv enters"):
                reading.species_rows(fs, LISTS)


if __name__ == "__main__":
    unittest.main()
