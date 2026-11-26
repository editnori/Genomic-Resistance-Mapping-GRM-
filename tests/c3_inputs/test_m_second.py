"""scripts/c3-inputs/second.py: the second reading, with the listed genomes out of every outer test fold.

Every fit, score file, list and floor table here is SYNTHETIC (m_common.py):
the ids are numbers of the form a genome id takes, made up in the test.
"""

import tempfile
import unittest
from pathlib import Path

from tests.c3_inputs import common, m_common

second = common.load("second")
score = common.load("score")
expect = common.load("expect")
checks = m_common.checks
D = m_common.M_SET
UNDER = "synthia_coli__drug_n"
NAMES = {D: ("Synthia coli", "drug m"), UNDER: ("Synthia coli", "drug n")}
ARMS = ["program", "k31-counts"]


def floor_row(dataset, gs, listed):
    r = sum(1 for g in gs if g[1])
    s = len(gs) - r
    lr = sum(1 for g in gs if g[0] in listed and g[1])
    ls = len(listed) - lr
    second_r, second_s = r - lr, s - ls
    return {"dataset": dataset, "resistant": str(r), "susceptible": str(s), "listed": str(len(listed)),
            "listed_resistant": str(lr), "listed_susceptible": str(ls), "second_resistant": str(second_r),
            "second_susceptible": str(second_s), "floor": "25", "enters": str(int(min(second_r, second_s) >= 25))}


class Fixture:
    """One dataset that enters (80 genomes, 40 R) and one under the floor (50 genomes, 25 R, two R listed)."""

    def __init__(self, top, listed=None, under=True, refused=(), change=None):
        self.gs = {D: m_common.genomes(80), UNDER: m_common.genomes(50)}
        ids = [g[0] for g in self.gs[D]]
        self.listed = {D: set(ids[3:7]) | {ids[40]} if listed is None else set(listed),
                       UNDER: {self.gs[UNDER][0][0], self.gs[UNDER][2][0]}}
        datasets = [D, UNDER] if under else [D]
        rows = expect.rows(m_common.SPECIES, datasets, datasets, datasets, [42, 43], [45], ARMS)
        self.expect, self.picks, self.scores, _, self.study = m_common.lay_out(
            top, rows, lambda e: self.gs[e["dataset"]], refused=refused)
        self.floor = [floor_row(d, self.gs[d], self.listed[d]) for d in datasets]
        if change:
            change(self)

    def run(self):
        return second.second(self.expect, self.picks, self.scores, self.study, NAMES, self.floor,
                             {r["dataset"]: self.listed[r["dataset"]] for r in self.floor})


def first_reading(fx, datasets=(D,)):
    """score.py's draws.tsv rows on the results as fitted, for the rows of the datasets given."""
    rows = [e for e in fx.expect if e["dataset"] in datasets]
    docs = {e["file"]: checks.load(fx.picks / e["file"]) for e in rows}
    return score.score_rows(rows, docs, NAMES)


class Again(unittest.TestCase):
    def test_only_the_listed_genomes_leave_and_training_stays_as_fitted(self):
        gs = m_common.genomes(80)
        doc = m_common.pick_from(gs)
        left = {gs[0][0], gs[5][0], gs[10][0], gs[2][0]}  # folds 0, 0, 0 and 2; three R
        got = second.again(doc, gs, left)
        kept = m_common.pick_from([g for g in gs if g[0] not in left])
        for f, (a, b, k) in enumerate(zip(got["folds"], doc["folds"], kept["folds"])):
            self.assertEqual((a["train"], a["train_resistant"]), (b["train"], b["train_resistant"]), f)
            self.assertEqual((a["test"], a["test_resistant"], a["outer"]), (k["test"], k["test_resistant"], k["outer"]))
            if f in (1, 3, 4):
                self.assertIs(a, b)
        self.assertEqual(got["folds"][0]["test"], doc["folds"][0]["test"] - 3)
        self.assertEqual((got["tested"], got["tested_resistant"]), (76, 37))
        self.assertEqual(doc["tested"], 80)

    def test_an_unavailable_split_stays_unavailable(self):
        doc = dict(m_common.pick_from([]), unavailable_reason="too few lineages")
        self.assertIs(second.again(doc, [], {"90007.1"}), doc)


class Reading(unittest.TestCase):
    def test_the_rows_are_the_scores_of_the_results_less_the_listed_genomes(self):
        with tempfile.TemporaryDirectory() as d:
            fx = Fixture(d)
            draws, enters, under = fx.run()
            self.assertEqual((enters, under), ([D], [UNDER]))
            self.assertEqual({r["dataset"] for r in draws}, {D})
            kept = [g for g in fx.gs[D] if g[0] not in fx.listed[D]]
            want = {e["file"]: m_common.pick_from(kept, e["dataset"], e["scheme"], int(e["seed"]))
                    for e in fx.expect if e["dataset"] == D}
            self.assertEqual(draws, score.score_rows([e for e in fx.expect if e["dataset"] == D], want, NAMES))
            self.assertTrue(all(r["n"] == 75 for r in draws))

    def test_with_no_listed_genome_the_rows_are_the_first_readings(self):
        with tempfile.TemporaryDirectory() as d:
            fx = Fixture(d, listed=set(), under=False)
            draws, enters, under = fx.run()
            self.assertEqual((enters, under), ([D], []))
            self.assertEqual(draws, first_reading(fx))

    def test_a_dataset_under_the_floor_has_no_row(self):
        with tempfile.TemporaryDirectory() as d:
            fx = Fixture(d)
            self.assertEqual(fx.floor[1]["enters"], "0")
            self.assertEqual((fx.floor[1]["second_resistant"], fx.floor[1]["second_susceptible"]), ("23", "25"))
            draws, _, under = fx.run()
            self.assertEqual(under, [UNDER])
            self.assertFalse(any(r["dataset"] == UNDER for r in draws))

    def test_a_refused_study_fit_is_the_constant_model_on_the_defaults_folds_less_the_list(self):
        name = f"program/{D}.study-whole-s42.pick.json"
        with tempfile.TemporaryDirectory() as d:
            fx = Fixture(d, refused={name}, under=False)
            draws, _, _ = fx.run()
            got = [r for r in draws if r["arm"] == "program" and r["split"] == "study"]
            self.assertEqual(len(got), len(m_common.LEARNERS))
            base = second.again(checks.load(fx.picks / f"k31/{D}.study-whole-s42.pick.json"), fx.gs[D], fx.listed[D])
            want = score.measures(score.constant_calls(base))
            for r in got:
                self.assertEqual((r["constant"], r["n"]), (1, 75))
                self.assertEqual({k: r[k] for k in want}, want)


class Stops(unittest.TestCase):
    def stops(self, change, pattern):
        with tempfile.TemporaryDirectory() as d:
            fx = Fixture(d, change=change)
            with self.assertRaisesRegex(second.Stop, pattern):
                fx.run()

    def test_counts_that_do_not_come_back_from_the_scores_stop(self):
        def change(fx):
            path = fx.picks / "k31-counts" / f"{D}.lineage-s43.pick.json"
            doc = checks.load(path)
            doc["folds"][1]["outer"][0]["fp"] += 1
            doc["folds"][1]["outer"][0]["tn"] -= 1
            m_common.write_json(path, doc)
        self.stops(change, "do not come back from the scores")

    def test_a_listed_genome_outside_the_cohort_stops(self):
        def change(fx):
            fx.listed[D] = fx.listed[D] | {"90007.999"}
            fx.floor[0] = floor_row(D, fx.gs[D], fx.listed[D])
        self.stops(change, "not in the fit's cohort")

    def test_kept_counts_unlike_the_floor_tables_stop(self):
        def change(fx):
            fx.floor[0].update(listed_resistant=str(int(fx.floor[0]["listed_resistant"]) + 1),
                               listed_susceptible=str(int(fx.floor[0]["listed_susceptible"]) - 1),
                               second_resistant=str(int(fx.floor[0]["second_resistant"]) - 1),
                               second_susceptible=str(int(fx.floor[0]["second_susceptible"]) + 1))
        self.stops(change, "the fit keeps")

    def test_a_floor_table_unlike_its_own_rule_stops(self):
        self.stops(lambda fx: fx.floor[1].update(enters="1"), "floor of 25 a class")
        self.stops(lambda fx: fx.floor[0].update(second_resistant="1"), "cohort's less the list")
        self.stops(lambda fx: fx.listed[D].pop(), "the list file")
        self.stops(lambda fx: fx.floor.pop(), "not expect.tsv's")

    def test_a_missing_result_or_score_file_stops(self):
        self.stops(lambda fx: (fx.picks / f"k31/{D}.lineage-s43.pick.json").unlink(), "no result")
        self.stops(lambda fx: (fx.scores / f"program/{D}.lineage-s42.scores.json").unlink(), "were not kept")


class Lists(unittest.TestCase):
    def test_a_list_file_naming_a_genome_twice_is_refused(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "x.txt"
            p.write_text("90007.1\n90007.2\n\n")
            self.assertEqual(second.read_list(p), {"90007.1", "90007.2"})
            p.write_text("90007.1\n90007.1\n")
            with self.assertRaises(second.Stop):
                second.read_list(p)


if __name__ == "__main__":
    unittest.main()
