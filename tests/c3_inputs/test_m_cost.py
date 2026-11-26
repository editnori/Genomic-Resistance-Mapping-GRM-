"""scripts/c3-inputs/cost.py: the cost of the run on the datasets of M, in the plan's units, before it starts.

model reads the committed lists.tsv, page-builds.tsv and the entry C3 arms'
cost model (docs/assumption-plan-corrections.md). measured runs on
SYNTHETIC step tables and copies: made-up timing tables, one per species of
M. Page reads the committed reference/input-arms/first-run-times.tsv and
holds docs/input-arms.md's Cost to what cost.py prints from it.
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


DATA = [("release", 10), ("build", 45), ("curated", 20), ("inputs_counts", 30), ("prodigal", 15), ("inputs_aa", 5),
        ("bcalm", 40), ("inputs_unitigs", 10)]
FITS = {"k31": 360, "unitigs": 200, "k31-counts": 150, "k31-aa": 170, "program": 40, "k31-program": 300, "twin": 300}
PROGRAM_ARMS = ("program", "k31-program")
M_SPECIES = sorted({r["species_slug"] for r in LISTS if r["in_M"] == "1"})
RAN = cost.CHECK_SPECIES   # in these synthetic tables the check species is the one whose first run ran the program


def rows_where(slug, arm, first, study):
    a = arm in PROGRAM_ARMS
    col = ("first_run_A" if a else "first_run") if first else "in_M"
    return [r for r in LISTS if r["species_slug"] == slug and r[col] == "1" and (first or not a or r["in_A"] == "1")
            and (not study or r["in_S"] == "1")]


def n(slug, arm, first, study=False):
    return sum(int(r["samples"]) for r in rows_where(slug, "k31" if arm == "twin" else arm, first, study))


def steps(slug):
    """A made-up first-run box: data steps, a seed-42 fit of each arm (the twin at 45), a study fit of each arm with
    first-run datasets in S, and the program: run where the first run ran it, its tables copied elsewhere."""
    out = list(DATA) + [(f"fit-{a}-s{45 if a == 'twin' else 42}", w) for a, w in FITS.items()]
    out += [(f"study-{a}-synth_x", 60) for a in cost.STUDY_ARMS if n(slug, a, True, study=True)]
    if any(r["in_A"] == "1" for r in rows_where(slug, "k31", False, False)):
        out += [("program_env", 30), ("program_run", 100)] if slug == RAN else [("program_tables", 4 if slug < "m" else 2)]
    return out


def table(change=None):
    rows = [{"species_slug": s, "step": step, "seconds": str(w), "source_sha256": "a" * 64}
            for s in M_SPECIES for step, w in steps(s)]
    return change(rows) if change else rows


def want(slug, program_again):
    lineage = sum(w * n(slug, a, False) / n(slug, a, True) for a, w in FITS.items() if n(slug, a, True))
    study = 0.0
    for a in cost.STUDY_ARMS:
        m_s, f_s = n(slug, a, False, True), n(slug, a, True, True)
        if m_s:
            study += 60 * m_s / f_s if f_s else FITS[a] * m_s / n(slug, a, True)
    data = sum(w for s, w in DATA if s not in ("build", "inputs_counts", "prodigal", "inputs_aa", "bcalm", "inputs_unitigs"))
    builds = 30 + 15 + 5 + 40 + 10
    program = 0.0
    if any(r["in_A"] == "1" for r in rows_where(slug, "k31", False, False)):
        if program_again:
            program = 30 + 100 * (15 / 15)
        else:
            program = 4 if slug == RAN else (4 if slug < "m" else 2)
    again = sum(FITS.values()) + 60 * sum(1 for a in cost.STUDY_ARMS if n(slug, a, True, True)) if slug == RAN else 0
    return (lineage + study + data + builds + 2 * 45 + program + again) / 3600


class Measured(unittest.TestCase):
    def test_each_species_costs_its_box_with_its_lineage_fits_scaled_to_its_datasets_of_m(self):
        parts, _ = cost.measured(table(), LISTS)
        got = dict((p[0].split(":")[0], p[1]) for p in parts)
        self.assertEqual(sorted(got), M_SPECIES)
        for slug, xh in got.items():
            self.assertAlmostEqual(xh, want(slug, False), places=9, msg=slug)

    def test_a_study_fit_is_scaled_by_the_samples_in_s(self):
        # Acinetobacter baumannii: its datasets of M in S over its first-run datasets in S, 4,897 over 576.
        slug = "acinetobacter_baumannii"
        self.assertEqual((n(slug, "k31", False, True), n(slug, "k31", True, True)), (4897, 576))
        parts, _ = cost.measured(table(), LISTS)
        line = next(p[0] for p in parts if p[0].startswith(slug))
        study = sum(60 * n(slug, a, False, True) / n(slug, a, True, True) for a in cost.STUDY_ARMS
                    if n(slug, a, True, True))
        self.assertIn(f"study fits {study / 3600:.2f},", line)

    def test_a_species_with_m_in_s_and_no_first_run_study_fit_is_costed_and_named(self):
        slug = "enterococcus_faecium"
        self.assertEqual((n(slug, "k31", False, True), n(slug, "k31", True, True)), (333, 0))
        parts, notes = cost.measured(table(), LISTS)
        self.assertTrue(any(note.startswith(f"{slug}: no first-run study fit") and "k31 (333 samples)" in note
                            for note in notes), notes)
        line = next(p[0] for p in parts if p[0].startswith(slug))
        self.assertNotIn("study fits 0.00,", line)

    def test_the_program_with_its_tables_reused_or_run_again(self):
        again, notes = cost.measured(table(), LISTS, program_again=True)
        for p in again:
            slug = p[0].split(":")[0]
            self.assertAlmostEqual(p[1], want(slug, True), places=9, msg=slug)
        reused, notes = cost.measured(table(), LISTS)
        self.assertIn(f"{RAN}: the program, the largest program_tables step of the others, 4 s", notes)
        tb = "mycobacterium_tuberculosis"
        self.assertFalse(any(r["in_A"] == "1" for r in rows_where(tb, "k31", False, False)))
        self.assertIn("program 0.00,", next(p[0] for p in again if p[0].startswith(tb)))

    def test_the_grm_build_is_counted_twice_and_the_left_out_steps_are_named(self):
        parts, _ = cost.measured(table(), LISTS)
        self.assertTrue(all(f"two grm builds {90 / 3600:.2f}," in p[0] for p in parts))
        with contextlib.redirect_stdout(io.StringIO()) as out:
            self.assertEqual(cost.main(["measured", "--tables-reused"]), 0)
        self.assertEqual(out.getvalue().splitlines()[-1], "left out: m_list, Same counts, score.py, second.py, the pack, "
                                                          "the git archive copies, box start-up and idle time")

    def test_a_table_with_no_build_or_a_fit_of_no_arm_stops(self):
        for change, why in ((lambda rows: [r for r in rows if r["step"] != "build"], "no build step"),
                            (lambda rows: rows + [dict(rows[0], step="fit-zzz-s42")], "no arm"),
                            (lambda rows: rows[1:] + [dict(rows[0], source_sha256="b" * 64)], "more than one source"),
                            (lambda rows: [r for r in rows if r["species_slug"] != RAN], "not the 9 species of M")):
            with self.assertRaisesRegex(cost.Stop, why):
                cost.measured(table(change), LISTS)

    def test_the_steps_fall_to_their_kind(self):
        self.assertEqual([cost.kind(s) for s in ("fit-k31-program-s42-m", "study-k31-counts-x", "fit-twin-s45",
                                                 "inputs_aa", "program_tables", "m_list", "same_build", "build")],
                         [("lineage", "k31-program"), ("study", "k31-counts"), ("lineage", "twin"), ("build", "k31-aa"),
                          ("program", None), ("data", None), ("data", None), ("grm", None)])


def copies(top, change=None):
    """SYNTHETIC copies of the first-run boxes: each species' logs/timing.tsv and SHA256SUMS."""
    for slug in M_SPECIES:
        f = Path(top) / slug
        (f / "logs").mkdir(parents=True)
        text = "step\twall_s\tpeak_kb\n" + "".join(f"{s}\t{w - 0.5:.2f}\t1000\n" for s, w in steps(slug))
        (f / "logs" / "timing.tsv").write_text(text)
        (f / "SHA256SUMS").write_text(f"{hashlib.sha256(text.encode()).hexdigest()}  logs/timing.tsv\n")
        if change:
            change(f)


class Times(unittest.TestCase):
    def test_the_times_form_writes_whole_seconds_with_each_tables_sha256(self):
        with tempfile.TemporaryDirectory() as d:
            copies(d)
            with contextlib.redirect_stdout(io.StringIO()) as out:
                self.assertEqual(cost.main(["times", d, str(Path(d) / "times.tsv")]), 0)
            rows = cost.read_tsv(Path(d) / "times.tsv")
            self.assertEqual(list(rows[0]), cost.TIMES_COLUMNS)
            # Half a second up: each made-up wall time is a whole second less 0.5.
            self.assertEqual([(r["species_slug"], r["step"], int(r["seconds"])) for r in rows],
                             [(s, step, w) for s in M_SPECIES for step, w in steps(s)])
            for slug in M_SPECIES:
                digest = hashlib.sha256((Path(d) / slug / "logs" / "timing.tsv").read_bytes()).hexdigest()
                self.assertEqual({r["source_sha256"] for r in rows if r["species_slug"] == slug}, {digest})
            self.assertEqual(cost.measured(rows, LISTS)[0], cost.measured(table(), LISTS)[0])

    def test_a_timing_table_unlike_its_sha256sums_stops(self):
        def change(f):
            p = f / "logs" / "timing.tsv"
            p.write_text(p.read_text().replace("\t359.50\t", "\t35.50\t"))
        with tempfile.TemporaryDirectory() as d:
            copies(d, change=change)
            with self.assertRaisesRegex(cost.Stop, "differs from its SHA256SUMS"):
                cost.times_rows(d, LISTS)

    def test_the_committed_table_holds_whole_seconds_of_the_9_species_of_m(self):
        rows = cost.read_tsv(cost.TIMES)
        self.assertEqual(sorted({r["species_slug"] for r in rows}), M_SPECIES)
        self.assertTrue(all(r["seconds"].isdigit() for r in rows))
        self.assertEqual(len({(r["species_slug"], r["source_sha256"]) for r in rows}), 9)


def printed(option):
    with contextlib.redirect_stdout(io.StringIO()) as out:
        assert cost.main(["measured", option]) == 0
    return next(l for l in out.getvalue().splitlines() if l.startswith("total, ")).split(": ", 1)[1]


class Page(unittest.TestCase):
    """docs/input-arms.md gives the measured figures as cost.py prints them from the committed table."""

    def test_the_pages_cost_figures_are_cost_pys_on_the_committed_table(self):
        page = (common.ROOT / "docs" / "input-arms.md").read_text()
        para = page[page.index("**Cost.**"):page.index("**Scaling.**")]
        reused, again = printed("--tables-reused"), printed("--program-again")
        self.assertIn(reused, para)
        self.assertIn(again, para)
        self.assertLess(para.index(reused), para.index(again))
        with contextlib.redirect_stdout(io.StringIO()) as out:
            self.assertEqual(cost.main(["model"]), 0)
        central = next(l for l in out.getvalue().splitlines() if l.startswith("total, central: ")).split(": ", 1)[1]
        upper = next(l for l in out.getvalue().splitlines() if l.startswith("total, upper: ")).split(": ", 1)[1]
        for figure in (central.split(",")[0], upper.split(",")[0]):
            self.assertIn(figure, para)
        self.assertLess(para.index(again), para.index(central.split(",")[0]))
        self.assertIn("upper bound", para)
        self.assertIn("45 xlarge-hours", para)


if __name__ == "__main__":
    unittest.main()
