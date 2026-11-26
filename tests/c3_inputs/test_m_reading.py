"""scripts/c3-inputs/m_reading.py: the three counts in both readings, and the page tables of the run on M.

The check of the first run runs for real on the committed record,
reference/input-arms/ (page-datasets.tsv, page-count.tsv,
page-conditions.tsv, SHA256SUMS). The counts and tables run on SYNTHETIC
species folders: made-up species, datasets, listed counts and figures, no
genome. A test of a rule names the cells and values it expects itself.
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
TWO = ("balanced_accuracy", "recall")
READ = "m_reading"
RECORD = [f"synth_{s}__drug_{i}" for s in "ab" for i in range(1, 5)]
M_SETS = {"synth_c": [f"synth_c__drug_{i}" for i in (5, 6, 7)], "synth_d": [f"synth_d__drug_{i}" for i in (8, 9, 10)]}
UNDER = "synth_d__drug_10"
# Listed genomes by class (R, S): 1 genome, 4 and 0, 5 and 0, none, 6 and 7, and the dataset under the floor, whose
# susceptible class lists 2: the floor is read first.
LISTED = {"synth_c__drug_5": (1, 0), "synth_c__drug_6": (4, 0), "synth_c__drug_7": (5, 0), "synth_d__drug_8": (0, 0),
          "synth_d__drug_9": (6, 7), UNDER: (6, 2)}
SHORT = {"synth_c__drug_5", "synth_c__drug_6"}
IN = {"synth_c__drug_7", "synth_d__drug_8", "synth_d__drug_9"}
STATE = {**{d: "short list" for d in SHORT}, **{d: "in" for d in IN}, UNDER: "under the floor"}
FLOOR_LINE = ("floor: datasets of M 6; second reading: 3 in; 1 left out under the floor; "
              "2 left out for a short list")
SEVEN = ("m-both.tsv", "m-conditions.tsv", "m-count.tsv", "m-datasets.tsv", "m-floor.tsv", "m-left-out.tsv",
         "m-size.tsv")


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


def floor_row(d, slug, r, s, lr, ls):
    """One floor.tsv row as m_list.py writes it: the floor of 25 read on the classes less the list."""
    return {"dataset": d, "species_slug": slug, "species": "Synth " + slug[6:],
            "antibiotic": d.split("__")[1].replace("_", " "), "cohort": r + s, "resistant": r, "susceptible": s,
            "listed": lr + ls, "listed_resistant": lr, "listed_susceptible": ls, "second_resistant": r - lr,
            "second_susceptible": s - ls, "floor": 25, "enters": int(min(r - lr, s - ls) >= 25)}


def floor_rows(slug, listed=LISTED):
    return [floor_row(d, slug, *((30, 40) if d == UNDER else (40, 40)), *listed[d]) for d in M_SETS[slug]]


def folders(top, listed=LISTED, scores=True):
    """One species folder of the run on M per synthetic species, as its pack holds it. With scores False, only
    floor.tsv: the floor form reads nothing else. The second reading's file holds every dataset that passes the
    floor, as second.py writes it, those with a short list too."""
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


def edit_second(fs, change):
    """Rewrite each species folder's second-reading rows through change(row), which edits a row in place."""
    for f in fs:
        path = f / "score" / "second" / "datasets.tsv"
        rows = score.read_tsv(path)
        for r in rows:
            change(r)
        m_common.write_tsv(path, score.DATASET_COLUMNS, rows)


def tables(top, fs, record):
    """The floor form's tables, then the count form's, written to <top>/tables. Returns the count form's lines."""
    out = Path(top) / "tables"
    reading.floor_form(out, fs, LISTS)
    return list(reading.count_form(record, out, fs, LISTS))


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
        self.assertEqual({r["measure"] for r in self.counted}, set(TWO))

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
    """Three counts in each reading on synthetic species folders, each its own call of score.py's count."""

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.record = rows_of(RECORD, TWO, 0)
        fs = folders(cls.tmp.name)
        cls.m_rows, cls.floor = reading.species_rows(fs, LISTS)
        cls.taken = reading.counts(cls.record, cls.m_rows)
        cls.out = Path(cls.tmp.name) / "tables"
        tables(cls.tmp.name, fs, cls.record)

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def table(self, name):
        return score.read_tsv(self.out / name)

    def header(self, name):
        return (self.out / name).read_text().splitlines()[0].split("\t")

    def test_three_counts_in_each_of_two_readings(self):
        self.assertEqual(sorted(self.taken), sorted((c, g) for c in ("first run", "M", "all datasets")
                                                    for g in ("first", "second")))
        for (c, g), (rows, conds) in self.taken.items():
            want = MEASURES if c == "M" else list(TWO)
            self.assertEqual(len(rows), 11 * len(LEARNERS) * len(want), (c, g))
            self.assertEqual({r["measure"] for r in rows}, set(want), (c, g))
            self.assertEqual(len(conds), 5 * 11, (c, g))
        lineage = lambda c, g: {r["datasets"] for r in self.taken[(c, g)][0] if r["split"] == "lineage"}
        self.assertEqual((lineage("first run", "first"), lineage("M", "first"), lineage("all datasets", "first")),
                         ({8}, {6}, {14}))
        self.assertEqual((lineage("M", "second"), lineage("all datasets", "second")), ({3}, {11}))
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
        self.assertEqual((m, rest), ({"balanced_accuracy", "recall", "specificity", "vme", "me", "log_loss"},
                                     {"balanced_accuracy", "recall"}))

    def test_the_left_out_table_names_each_dataset_left_out_with_its_reason_and_counts(self):
        # Both reasons, in floor.tsv's order: two lists of 1 and 4 resistant genomes, and the dataset under the floor.
        self.assertEqual(self.header("m-left-out.tsv"),
                         ["dataset", "species", "antibiotic", "reason", "cohort", "resistant", "susceptible", "listed",
                          "listed_resistant", "listed_susceptible", "second_resistant", "second_susceptible", "floor"])
        self.assertEqual([list(r.values()) for r in self.table("m-left-out.tsv")], [
            ["synth_c__drug_5", "Synth c", "drug 5", "short list", "80", "40", "40", "1", "1", "0", "39", "40", "25"],
            ["synth_c__drug_6", "Synth c", "drug 6", "short list", "80", "40", "40", "4", "4", "0", "36", "40", "25"],
            [UNDER, "Synth d", "drug 10", "under the floor", "70", "30", "40", "8", "6", "2", "24", "38", "25"]])

    def test_the_floor_table_gives_each_datasets_second_reading_in_one_of_three_states(self):
        rows = self.table("m-floor.tsv")
        self.assertEqual(self.header("m-floor.tsv")[-2:], ["enters", "second_reading"])
        self.assertEqual({r["dataset"]: r["second_reading"] for r in rows}, STATE)
        self.assertEqual({r["dataset"]: r["enters"] for r in rows},
                         {**{d: "1" for d in IN | SHORT}, UNDER: "0"})

    def test_the_size_table_gives_entered_better_and_worse_per_reading_arm_and_learner(self):
        rows = self.table("m-size.tsv")
        self.assertEqual(len(rows), 2 * len(FAMILY) * len(LEARNERS))
        self.assertEqual({r["arm"] for r in rows}, set(FAMILY))
        for r in rows:
            c = next(x for x in self.taken[("all datasets", r[READ])][0] if x["split"] == "lineage"
                     and x["measure"] == score.PRIMARY and x["pair"] == f"{r['arm']}-k31" and x["learner"] == r["learner"])
            self.assertEqual((r["entered"], r["better"], r["worse"], r["same"]),
                             (str(c["entered"]), str(c["better"]), str(c["worse"]), str(c["same"])))
            self.assertEqual(int(r["better"]) + int(r["worse"]) + int(r["same"]), int(r["entered"]))
        self.assertNotIn("pick", " ".join(self.header("m-size.tsv")))

    def test_the_both_table_meets_a_condition_only_when_both_readings_meet_it(self):
        rows = self.table("m-both.tsv")
        self.assertEqual(len(rows), 4 * 5)
        for r in rows:
            self.assertEqual(r["met_both"], str(int(r["met_first"] == "1" and r["met_second"] == "1")))

    def test_the_dataset_table_holds_the_first_reading_of_each_and_the_second_of_those_in_it(self):
        rows = self.table("m-datasets.tsv")
        first = {r["dataset"] for r in rows if r[READ] == "first"}
        second = {r["dataset"] for r in rows if r[READ] == "second"}
        self.assertEqual((first, second), (IN | SHORT | {UNDER}, IN))
        head = self.header("m-datasets.tsv")
        self.assertEqual((head[0], head[-3:]), ("m_reading", ["listed_resistant", "listed_susceptible",
                                                              "second_reading"]))
        self.assertEqual({r["dataset"]: r["second_reading"] for r in rows}, STATE)
        for r in rows:
            self.assertEqual((int(r["listed_resistant"]), int(r["listed_susceptible"])), LISTED[r["dataset"]])

    def test_each_second_reading_count_row_leaves_out_the_datasets_left_out(self):
        # 8 datasets in the record; 6 of M in the first reading, the 3 in it in the second: on every row, both splits.
        want = {("first run", "first"): 8, ("first run", "second"): 8, ("M", "first"): 6, ("M", "second"): 3,
                ("all datasets", "first"): 14, ("all datasets", "second"): 11}
        got = {}
        for r in self.table("m-count.tsv"):
            got.setdefault((r["count"], r[READ]), set()).add((r["datasets"], r["entered"], r["split"]))
        self.assertEqual(got, {k: {(str(n), str(n), "lineage"), (str(n), str(n), "study")} for k, n in want.items()})

    def test_no_table_holds_a_second_reading_figure_of_a_dataset_left_out(self):
        # Each dataset's second-reading rows carry a figure no other row has; a dataset left out leaves none behind.
        with tempfile.TemporaryDirectory() as d:
            fs = folders(d)
            marks = {}

            def mark(r):
                if r["arm"] != score.DEFAULT and r["learner"] == LEARNERS[1]:
                    r["balanced_accuracy"] = repr(0.512345678901 + sorted(LISTED).index(r["dataset"]) * 1e-9)
                    marks.setdefault(r["dataset"], set()).add(r["balanced_accuracy"][:12])

            edit_second(fs, mark)
            tables(d, fs, self.record)
            text = "".join((Path(d) / "tables" / name).read_text() for name in SEVEN)
            for ds, want in ((IN, True), (SHORT, False)):
                for name in ds:
                    self.assertEqual(any(m in text for m in marks[name]), want, name)

    def test_no_table_has_a_column_or_a_cell_for_a_dataset_counted_and_not_shown(self):
        self.assertEqual(self.header("m-count.tsv"),
                         ["count", "m_reading", "split", "pair", "learner", "measure", "orientation", "withdrawn",
                          "datasets", "entered", "not_scored", "constant", "better", "worse", "same",
                          "median_difference", "species", "interval_lo", "interval_hi", "sign_p", "holm_p"])
        self.assertEqual(self.header("m-size.tsv"),
                         ["m_reading", "arm", "learner", "datasets", "entered", "better", "worse", "same"])
        text = "".join((self.out / name).read_text() for name in SEVEN)
        for word in ("not shown", "not_shown", "figures", "shown"):
            self.assertNotIn(word, text)
        # Every second-reading count row gives its figures: none is blank where its first-reading row is not.
        first = {(r["count"], r["split"], r["pair"], r["learner"], r["measure"]): r
                 for r in self.table("m-count.tsv") if r[READ] == "first"}
        for r in self.table("m-count.tsv"):
            if r[READ] == "second":
                a = first[(r["count"], r["split"], r["pair"], r["learner"], r["measure"])]
                for c in ("better", "worse", "same", "median_difference", "interval_lo", "interval_hi", "sign_p",
                          "holm_p"):
                    self.assertEqual(r[c] == "", a[c] == "", (r["count"], r["pair"], r["measure"], c))

    def test_an_m_dataset_in_the_record_stops(self):
        clash = self.record + [dict(r) for r in self.m_rows["first"] if r["dataset"] == "synth_c__drug_5"]
        with self.assertRaisesRegex(reading.Stop, "in the first run's record"):
            reading.counts(clash, self.m_rows)


class SmallLists(unittest.TestCase):
    def test_a_class_listing_1_to_4_genomes_leaves_the_second_reading_and_0_or_5_does_not(self):
        row = lambda r, s, enters="1": {"listed_resistant": str(r), "listed_susceptible": str(s), "enters": enters}
        cases = {(1, 0): "short list", (2, 0): "short list", (3, 0): "short list", (4, 0): "short list",
                 (0, 1): "short list", (0, 2): "short list", (0, 3): "short list", (0, 4): "short list",
                 (5, 4): "short list", (4, 5): "short list", (5, 1): "short list", (1, 1): "short list",
                 (0, 0): "in", (5, 0): "in", (0, 5): "in", (5, 5): "in", (6, 7): "in", (40, 0): "in"}
        for (r, s), want in cases.items():
            self.assertEqual(reading.second_state(row(r, s)), want, (r, s))
        # The floor is read first: under it, a dataset is left out under the floor, whatever it lists.
        for r, s in ((0, 0), (5, 0), (1, 0), (4, 5)):
            self.assertEqual(reading.second_state(row(r, s, "0")), "under the floor", (r, s))
        self.assertEqual(reading.SMALLEST_LISTED, 5)

    def test_the_floor_form_gives_each_lists_state_from_floor_tsv_alone(self):
        slug = "synth_e"
        cases = [((40, 40, 1, 0), "short list"), ((40, 40, 2, 0), "short list"), ((40, 40, 0, 3), "short list"),
                 ((40, 40, 4, 0), "short list"), ((40, 40, 0, 0), "in"), ((40, 40, 5, 0), "in"),
                 ((40, 40, 0, 5), "in"), ((40, 40, 5, 4), "short list"), ((40, 40, 4, 5), "short list"),
                 ((40, 40, 5, 5), "in"), ((26, 40, 2, 0), "under the floor"), ((24, 40, 0, 0), "under the floor")]
        names = [f"{slug}__drug_{i}" for i in range(1, len(cases) + 1)]
        with tempfile.TemporaryDirectory() as d:
            top = Path(d) / slug
            m_common.write_tsv(top / "results" / "m" / "floor.tsv", reading.LIST.FLOOR_COLUMNS,
                               [floor_row(n, slug, *c) for n, (c, _) in zip(names, cases)])
            out = Path(d) / "tables"
            said = reading.floor_form(out, [top], [{"dataset": n, "in_M": "1"} for n in names])
            self.assertEqual(said, "floor: datasets of M 12; second reading: 4 in; 2 left out under the floor; "
                                   "6 left out for a short list")
            self.assertEqual([r["second_reading"] for r in score.read_tsv(out / "m-floor.tsv")],
                             [want for _, want in cases])
            self.assertEqual([(r["dataset"], r["reason"]) for r in score.read_tsv(out / "m-left-out.tsv")],
                             [(n, want) for n, (_, want) in zip(names, cases) if want != "in"])

    def test_changing_the_second_reading_rows_of_a_dataset_left_out_for_a_short_list_changes_no_table(self):
        # Every second-reading row of the two datasets with a short list is changed, three ways and by any amount:
        # the 7 tables and the lines come out byte for byte as they do with those rows as second.py wrote them.
        record = rows_of(RECORD, TWO, 0)
        fig = lambda x: float(x) if x not in ("", None) else 0.5

        def ways(k):
            def change(r):
                if r["dataset"] not in SHORT:
                    return
                for m in MEASURES:
                    r[m] = repr(round((1 - fig(r[m]), 0.0, (fig(r[m]) + 0.123456789012) % 1)[k], 12))
                for c, v in (("n", "7"), ("tp", "1"), ("fp", "2"), ("tn", "3"), ("fn", "1")):
                    r[c] = v if k else str(int(r[c]) * 3)
                if k == 2:
                    r["scored"], r["scored_draws"] = "0", "0"
            return change

        with tempfile.TemporaryDirectory() as d:
            base = Path(d) / "base"
            lines = tables(base, folders(base), record)
            want = {name: (base / "tables" / name).read_bytes() for name in SEVEN}
            self.assertEqual(sorted(p.name for p in (base / "tables").iterdir()), list(SEVEN))
            for k in range(3):
                top = Path(d) / f"way{k}"
                fs = folders(top)
                before = (fs[0] / "score" / "second" / "datasets.tsv").read_bytes()
                edit_second(fs, ways(k))
                self.assertNotEqual((fs[0] / "score" / "second" / "datasets.tsv").read_bytes(), before, k)
                self.assertEqual(tables(top, fs, record), lines, k)
                for name in SEVEN:
                    self.assertEqual((top / "tables" / name).read_bytes(), want[name], (k, name))


class Floor(unittest.TestCase):
    def test_the_floor_form_writes_the_left_out_and_floor_tables_from_floor_tsv_alone(self):
        with tempfile.TemporaryDirectory() as d:
            fs = folders(d, scores=False)
            out = Path(d) / "tables"
            said = reading.floor_form(out, fs, LISTS)
            self.assertEqual(sorted(p.name for p in out.iterdir()), ["m-floor.tsv", "m-left-out.tsv"])
            self.assertEqual([r["dataset"] for r in score.read_tsv(out / "m-left-out.tsv")],
                             ["synth_c__drug_5", "synth_c__drug_6", UNDER])
            self.assertEqual(len(score.read_tsv(out / "m-floor.tsv")), 6)
            self.assertEqual(said, FLOOR_LINE)
            self.assertFalse(any((f / "score").exists() for f in fs))
            with contextlib.redirect_stdout(io.StringIO()) as printed:
                self.assertEqual(reading.main(["floor", str(out), str(Path(d) / "none")]), 1)
            self.assertIn("m_reading: stopped: none: no results/m/floor.tsv", printed.getvalue())

    def test_the_floor_forms_tables_are_the_count_forms_and_it_opens_no_score_file(self):
        # In full species folders each score file is made a directory, so opening one fails; the floor form runs.
        record = rows_of(RECORD, TWO, 0)
        with tempfile.TemporaryDirectory() as d:
            fs = folders(d)
            out = Path(d) / "tables"
            kept = {}
            for f in fs:
                for p in (f / "score" / "datasets.tsv", f / "score" / "second" / "datasets.tsv"):
                    kept[p] = p.read_bytes()
                    p.unlink()
                    p.mkdir()
            try:
                line = reading.floor_form(out, fs, LISTS)
            except OSError as e:
                self.fail(f"the floor form opened a score file: {e}")
            self.assertEqual(line, FLOOR_LINE)
            with self.assertRaises(IsADirectoryError):
                reading.species_rows(fs, LISTS)
            for p, data in kept.items():
                p.rmdir()
                p.write_bytes(data)
            floor_made = {name: (out / name).read_bytes() for name in ("m-floor.tsv", "m-left-out.tsv")}
            # The count form holds them to its own, written from the same floor.tsv files, and leaves them as they are.
            said = list(reading.count_form(record, out, fs, LISTS))
            self.assertEqual(said[0], "floor: m-floor.tsv and m-left-out.tsv are the floor form's")
            own = Path(d) / "own"
            reading.write_floor(own, reading.species_rows(fs, LISTS)[1])
            for name, data in floor_made.items():
                self.assertEqual((own / name).read_bytes(), data, name)
                self.assertEqual((out / name).read_bytes(), data, name)

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
        record = rows_of(RECORD, TWO, 0)
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
            self.assertTrue(said[1].startswith("count: datasets of M 6; "), said[1])


class Lines(unittest.TestCase):
    def test_both_forms_lines_give_the_three_states_counts(self):
        # A driver reads each line's start: "floor: datasets of M " from the floor form; "floor: " and
        # "count: datasets of M " from the count form, after main's "check: ".
        record = rows_of(RECORD, TWO, 0)
        with tempfile.TemporaryDirectory() as d:
            fs = folders(d)
            out = Path(d) / "tables"
            self.assertEqual(reading.floor_form(out, fs, LISTS), FLOOR_LINE)
            self.assertEqual(list(reading.count_form(record, out, fs, LISTS)), [
                "floor: m-floor.tsv and m-left-out.tsv are the floor form's",
                "count: datasets of M 6; second reading: 3 in; 1 left out under the floor; 2 left out for a short list; "
                "first run, first reading: 110 count rows; M, first reading: 330 count rows; "
                "all datasets, first reading: 110 count rows; first run, second reading: 110 count rows; "
                "M, second reading: 330 count rows; all datasets, second reading: 110 count rows"])


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

    def test_a_second_reading_lacking_a_dataset_that_passes_the_floor_stops(self):
        # A dataset left out for a short list passes the floor, so second.py writes its rows; lacking them, or a
        # dataset in the second reading, the count form stops after the floor form's tables and before any count.
        record = rows_of(RECORD, TWO, 0)
        for gone in ("synth_c__drug_5", "synth_c__drug_6", "synth_c__drug_7"):
            with tempfile.TemporaryDirectory() as d:
                fs = folders(d)
                path = fs[0] / "score" / "second" / "datasets.tsv"
                m_common.write_tsv(path, score.DATASET_COLUMNS,
                                   [r for r in score.read_tsv(path) if r["dataset"] != gone])
                with self.assertRaisesRegex(reading.Stop, "synth_c: the second reading's datasets are not those "
                                                          "floor.tsv enters"):
                    reading.species_rows(fs, LISTS)
                out = Path(d) / "tables"
                reading.floor_form(out, fs, LISTS)
                with self.assertRaisesRegex(reading.Stop, "not those floor.tsv enters"):
                    list(reading.count_form(record, out, fs, LISTS))
                self.assertEqual(sorted(p.name for p in out.iterdir()), ["m-floor.tsv", "m-left-out.tsv"], gone)


if __name__ == "__main__":
    unittest.main()
