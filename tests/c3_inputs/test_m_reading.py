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


def floor_rows(slug):
    out = []
    for d in M_SETS[slug]:
        lr, ls = (6, 2) if d == UNDER else (1, 1)
        r, s = (30, 40) if d == UNDER else (40, 40)
        out.append({"dataset": d, "species_slug": slug, "species": "Synth " + slug[6:],
                    "antibiotic": d.split("__")[1].replace("_", " "), "cohort": r + s, "resistant": r,
                    "susceptible": s, "listed": lr + ls, "listed_resistant": lr, "listed_susceptible": ls,
                    "second_resistant": r - lr, "second_susceptible": s - ls, "floor": 25,
                    "enters": int(min(r - lr, s - ls) >= 25)})
    return out


def folders(top):
    """One species folder of the run on M per synthetic species, as its pack holds it."""
    out = []
    start = len(RECORD)
    for slug, datasets in M_SETS.items():
        f = Path(top) / slug
        first = rows_of(datasets, MEASURES, start)
        second = rows_of([d for d in datasets if d != UNDER], MEASURES, start, shift=4)
        m_common.write_tsv(f / "score" / "datasets.tsv", score.DATASET_COLUMNS, first)
        m_common.write_tsv(f / "score" / "second" / "datasets.tsv", score.DATASET_COLUMNS, second)
        m_common.write_tsv(f / "results" / "m" / "floor.tsv", reading.LIST.FLOOR_COLUMNS, floor_rows(slug))
        start += len(datasets)
        out.append(f)
    return out


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
        reading.write(cls.out, cls.taken, cls.floor, cls.m_rows)

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

    def test_the_dataset_table_holds_both_readings_of_each_dataset_of_m(self):
        rows = self.table("m-datasets.tsv")
        first = {r["dataset"] for r in rows if r[reading.READ] == "first"}
        second = {r["dataset"] for r in rows if r[reading.READ] == "second"}
        self.assertEqual((len(first), len(second), first - second), (6, 5, {UNDER}))
        self.assertEqual(list(rows[0]), reading.M_DATASET_COLUMNS)

    def test_an_m_dataset_in_the_record_stops(self):
        clash = self.record + [dict(r) for r in self.m_rows["first"] if r["dataset"] == "synth_c__drug_5"]
        with self.assertRaisesRegex(reading.Stop, "in the first run's record"):
            reading.counts(clash, self.m_rows)


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
