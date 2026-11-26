"""scripts/c3-inputs/cost.py: the cost of the run on the datasets of M, in the plan's units, before it starts.

model reads the committed lists.tsv, page-builds.tsv and the entry C3 arms'
cost model (docs/assumption-plan-corrections.md). measured runs on
SYNTHETIC copies: made-up timing tables, one per species of M.
"""

import contextlib
import hashlib
import io
import tempfile
import unittest
from pathlib import Path

from tests.c3_inputs import common

cost = common.load("cost")
LISTS = cost.read_tsv(common.KIT / "lists.tsv")
NAMES = {"k31": "K31", "twin": "the twin", "k31-counts": "K31 plus counts", "k31-aa": "K31 plus amino-acid k-mers",
         "unitigs": "unitigs", "program": "the program alone", "k31-program": "K31 plus the program"}


def entry_cost():
    """The entry C3 arms' Cost paragraph, as committed."""
    text = (common.ROOT / "docs" / "assumption-plan-corrections.md").read_text()
    paragraphs = [p for p in text.split("\n\n") if p.startswith("**Cost.** Each arm's fit costs its ratio on the pilot")]
    assert len(paragraphs) == 1
    return paragraphs[0]


class Model(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Read here, not at import: page_faults.py runs this suite in a copy of scripts/ and tests/ alone.
        cls.builds = cost.read_tsv(common.ROOT / "reference" / "input-arms" / "page-builds.tsv")

    def test_the_figures_are_the_entry_c3_arms_cost_models(self):
        text = entry_cost()
        for arm, (c, u) in cost.ARM_COST.items():
            self.assertIn(f"{NAMES[arm]} {c:.2f} and {u:.2f}", text, arm)
        self.assertIn(f"the program on every species of A {cost.PROGRAM_ON_SPECIES[0]:.2f} and "
                      f"{cost.PROGRAM_ON_SPECIES[1]:.2f}", text)
        self.assertIn(f"{cost.RUN_TOTAL[0]:.2f} xlarge-hours central and {cost.RUN_TOTAL[1]:.2f} upper", text)
        self.assertIn(f"{cost.FIRST_RUN_BOXES} boxes of {cost.XLARGE_CPUS} CPUs and 32 GB", text)
        self.assertIn(f"${cost.DOLLARS_PER_XLARGE_HOUR:.2f} an hour", text)

    def test_the_data_steps_are_the_runs_total_less_its_arm_rows(self):
        parts = cost.model(LISTS, self.builds)
        data = next(p for p in parts if p[0].startswith("data steps"))
        boxes = len({r["species_slug"] for r in LISTS if r["in_M"] == "1"})
        self.assertEqual(boxes, 9)
        self.assertEqual((round(data[1] * 23 / boxes, 2), round(data[2] * 23 / boxes, 2)), (3.89, 5.65))
        self.assertEqual(cost.units(1.0), "1.00 xlarge-hours, 16.0 core-hours, $0.20")

    def test_a_matrix_arm_takes_the_larger_of_its_samples_share_and_the_genomes_share(self):
        parts = {p[0].split(":")[0]: p for p in cost.model(LISTS, self.builds)}
        share = lambda arm: float(parts[arm][0].rsplit(" ", 1)[1])
        g = cost.genomes(self.builds)
        m = {r["species_slug"] for r in LISTS if r["in_M"] == "1"}
        by_genomes = round(sum(g[s] for s in m) / sum(g.values()), 4)
        for arm in cost.MATRIX_ARMS:
            self.assertEqual(share(arm), max(by_genomes, round(cost.samples(LISTS, arm, False) /
                                                                cost.samples(LISTS, arm, True), 4)), arm)
        self.assertLess(share("k31"), by_genomes)

    def test_the_printed_totals_on_the_committed_files(self):
        with contextlib.redirect_stdout(io.StringIO()) as out:
            self.assertEqual(cost.main(["model"]), 0)
        lines = out.getvalue().splitlines()
        self.assertIn("total, central: 44.76 xlarge-hours, 716.1 core-hours, $8.95", lines)
        self.assertIn("total, upper: 102.77 xlarge-hours, 1644.4 core-hours, $20.55", lines)
        self.assertTrue(lines[-1].startswith("with the program's tables reused: central 27.77 xlarge-hours"), lines[-1])

    def test_a_first_run_unlike_page_builds_stops(self):
        with self.assertRaisesRegex(cost.Stop, "not page-builds.tsv's"):
            cost.model(LISTS, [r for r in self.builds if r["species_slug"] != self.builds[0]["species_slug"]])


STEPS = [("release", 10), ("build", 45), ("curated", 20), ("inputs_counts", 30), ("prodigal", 15), ("inputs_aa", 5),
         ("bcalm", 40), ("inputs_unitigs", 10), ("program_run", 100), ("fit-k31-s42", 360), ("fit-twin-s45", 300),
         ("study-k31-program-x", 60)]


def copies(top, steps=STEPS, change=None):
    for slug in sorted({r["species_slug"] for r in LISTS if r["in_M"] == "1"}):
        f = Path(top) / slug
        (f / "logs").mkdir(parents=True)
        rows = [s for s in steps if s[0] != "study-k31-program-x" or first_samples(slug, "k31-program")]
        text = "step\twall_s\tpeak_kb\n" + "".join(f"{s}\t{w:.2f}\t1000\n" for s, w in rows)
        (f / "logs" / "timing.tsv").write_text(text)
        (f / "SHA256SUMS").write_text(f"{hashlib.sha256(text.encode()).hexdigest()}  logs/timing.tsv\n")
        if change:
            change(f)


def first_samples(slug, arm):
    col = "first_run_A" if arm in ("program", "k31-program") else "first_run"
    return sum(int(r["samples"]) for r in LISTS if r["species_slug"] == slug and r[col] == "1")


def m_samples(slug, arm):
    a = arm in ("program", "k31-program")
    return sum(int(r["samples"]) for r in LISTS if r["species_slug"] == slug and r["in_M"] == "1"
               and (not a or r["in_A"] == "1"))


class Measured(unittest.TestCase):
    def test_each_species_costs_its_first_run_box_with_its_fits_scaled_to_its_datasets_of_m(self):
        with tempfile.TemporaryDirectory() as d:
            copies(d)
            parts = dict((p[0].split(":")[0], p[1]) for p in cost.measured(d, LISTS))
        self.assertEqual(len(parts), 9)
        for slug, xh in parts.items():
            fits = sum(w * m_samples(slug, "k31") / first_samples(slug, "k31") for _, w in (STEPS[9], STEPS[10]))
            if first_samples(slug, "k31-program"):
                fits += 60 * m_samples(slug, "k31-program") / first_samples(slug, "k31-program")
            again = 360 + 300 + (60 if first_samples(slug, "k31-program") else 0) if slug == cost.CHECK_SPECIES else 0
            want = fits + 100 + 100 + 75 + 45 + again
            self.assertAlmostEqual(xh, want / 3600, places=9, msg=slug)

    def test_a_timing_table_unlike_its_sha256sums_stops(self):
        def change(f):
            p = f / "logs" / "timing.tsv"
            p.write_text(p.read_text().replace("\t360.00\t", "\t36.00\t"))
        with tempfile.TemporaryDirectory() as d:
            copies(d, change=change)
            with self.assertRaisesRegex(cost.Stop, "differs from its SHA256SUMS"):
                cost.measured(d, LISTS)

    def test_a_table_with_no_build_or_a_fit_of_no_arm_stops(self):
        for steps, why in (([s for s in STEPS if s[0] != "build"], "no build step"),
                           (STEPS + [("fit-zzz-s42", 1)], "no arm")):
            with tempfile.TemporaryDirectory() as d:
                copies(d, steps=steps)
                with self.assertRaisesRegex(cost.Stop, why):
                    cost.measured(d, LISTS)

    def test_the_steps_fall_to_their_kind(self):
        self.assertEqual([cost.kind(s) for s in ("fit-k31-program-s42-m", "study-k31-counts-x", "fit-twin-s45",
                                                 "inputs_aa", "program_tables", "m_list", "same_build")],
                         [("fit", "k31-program"), ("fit", "k31-counts"), ("fit", "twin"), ("build", "k31-aa"),
                          ("program", None), ("data", None), ("data", None)])


if __name__ == "__main__":
    unittest.main()
