"""scripts/mechanism/combine.py: the readings over species (docs/mechanism-split.md).

Every species, dataset, table row and figure is SYNTHETIC.
"""

import contextlib
import hashlib
import io
import random
import statistics
import tempfile
import unittest
from pathlib import Path

from tests.mechanism import common as c

combine = c.load("combine")
C = c.load("common")

STRATA = ["dataset", "species_slug", "species", "antibiotic", "drug_class", "set", "labels_depend_on", "split", "row",
          "stratum", "kind", "resistant", "resistant_called", "susceptible", "susceptible_called", "recall"]
CELLS = ["species_slug", "species", "cell", "class", "dataset", "antibiotic", "periods", "guard", "genomes",
         "with_year_share", "earlier_genomes", "later_genomes", "earlier_units", "later_units", "passes",
         "earlier_point", "later_point", "d_point", "p_point", "earlier_gene", "later_gene", "d_gene", "tvd", "p_tvd",
         "note"]
TS = ["dataset", "species_slug", "species", "antibiotic", "reading", "resistant_later_test", "with_determinant",
      "with_unseen_determinant", "no_determinant", "share_unseen"]
CONTRASTS = ["dataset", "species", "learner", "contrast", "d_recall", "d_log_loss", "d_specificity",
             "d_balanced_accuracy"]


def strata_rows(sp, diffs, split="lineage", dset="primary"):
    """diffs: {row: gene_only less point_only} for one dataset of species sp."""
    out = []
    d = f"{sp}__drug_a"
    for row, diff in diffs.items():
        for stratum, rec in (("point_only", 0.5), ("gene_only", 0.5 + diff), ("neither", 0.2)):
            out.append({"dataset": d, "species_slug": sp, "species": sp, "set": dset, "split": split, "row": row,
                        "stratum": stratum, "kind": "stratum", "resistant": "10", "recall": f"{rec:.12f}"})
    return out


def lay(top, species_rows, cells=(), ts=(), done=True):
    runs = []
    for sp, rows in species_rows.items():
        run = Path(top) / sp
        c.write_tsv(run / "results" / "m1" / "strata.tsv", STRATA, rows)
        c.write_tsv(run / "results" / "m2" / "m2-cells.tsv", CELLS, [r for r in cells if r["species_slug"] == sp])
        c.write_tsv(run / "results" / "m2" / "timesplit.tsv", TS, [r for r in ts if r["species_slug"] == sp])
        if done:
            (run / "DONE").write_text("DONE\n")
        runs.append(run)
    return runs


def go(top, runs, contrasts=()):
    c.write_tsv(Path(top) / "contrasts.tsv", CONTRASTS, contrasts)
    argv = ["--contrasts", str(Path(top) / "contrasts.tsv"), "--out", str(Path(top) / "out")]
    for r in runs:
        argv += ["--run", str(r)]
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
        rc = combine.main(argv)
    return rc, Path(top) / "out"


class M1(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)

    def test_reading_per_learner_with_holm_over_five(self):
        sp = {f"synthia_s{i:02d}": strata_rows(f"synthia_s{i:02d}", {"scm": 0.3, "cart": -0.3, "rf": 0.01})
              for i in range(10)}
        rc, out = go(self.tmp.name, lay(self.tmp.name, sp))
        self.assertEqual(rc, 0)
        got = {(r["split"], r["set"], r["row"]): r for r in c.read_tsv(out / "m1-reading.tsv")}
        scm, cart, rf = (got[("lineage", "primary", x)] for x in ("scm", "cart", "rf"))
        self.assertEqual((scm["datasets"], scm["species"], scm["gene_ahead"], scm["point_ahead"]), ("10", "10", "10", "0"))
        self.assertEqual((scm["median"], scm["lo"], scm["hi"]), ("0.300000000000",) * 3)
        self.assertEqual(scm["sign_p"], f"{2 / 1024:.12f}")
        self.assertEqual(scm["holm_p"], f"{5 * 2 / 1024:.12f}")
        self.assertEqual(scm["reading"], "gene only higher")
        self.assertEqual(cart["reading"], "point only higher")
        self.assertEqual((rf["same"], rf["reading"]), ("10", "no difference shown"))
        self.assertEqual(got[("lineage", "primary", "lasso")]["reading"], "no difference shown")
        self.assertEqual(got[("lineage", "primary", "program")]["decides"], "0")

    def test_a_set_under_six_species_has_no_interval(self):
        sp = {f"synthia_s{i}": strata_rows(f"synthia_s{i}", {"scm": 0.3}) for i in range(5)}
        rc, out = go(self.tmp.name, lay(self.tmp.name, sp))
        scm = {(r["split"], r["set"], r["row"]): r for r in c.read_tsv(out / "m1-reading.tsv")}[("lineage", "primary", "scm")]
        self.assertEqual((scm["median"], scm["lo"], scm["hi"]), ("0.300000000000", "", ""))

    def test_random_split_is_labelled_and_decides_nothing(self):
        sp = {"synthia_s1": strata_rows("synthia_s1", {"scm": 0.3}, split="random")}
        rc, out = go(self.tmp.name, lay(self.tmp.name, sp))
        r = {(r["split"], r["set"], r["row"]): r for r in c.read_tsv(out / "m1-reading.tsv")}[("random", "primary", "scm")]
        self.assertEqual((r["decides"], r["label"], r["reading"]), ("0", "training and test share lineages", ""))

    def test_neither_stratum_median(self):
        sp = {f"synthia_s{i}": strata_rows(f"synthia_s{i}", {"scm": 0.3}) for i in range(3)}
        rc, out = go(self.tmp.name, lay(self.tmp.name, sp))
        r = {(r["split"], r["set"], r["row"]): r for r in c.read_tsv(out / "m1-neither.tsv")}[("lineage", "primary", "scm")]
        self.assertEqual((r["datasets"], r["median_recall"]), ("3", "0.200000000000"))

    def test_a_run_with_no_done_is_not_read(self):
        runs = lay(self.tmp.name, {"synthia_s1": strata_rows("synthia_s1", {"scm": 0.3})}, done=False)
        with self.assertRaises(SystemExit):
            go(self.tmp.name, runs)


class M2(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)

    def cell(self, cls, passes, p, guard="bioproject"):
        return {"species_slug": "synthia_s1", "species": "Synthia s1", "cell": "species_class", "class": cls,
                "periods": "primary", "guard": guard, "passes": passes, "d_point": "0.1", "p_point": p, "tvd": "0.1",
                "p_tvd": p}

    def test_holm_over_passing_cells_only(self):
        cells = [self.cell("A", "1", "0.001"), self.cell("B", "1", "0.01"), self.cell("C", "1", "0.04"),
                 self.cell("D", "1", "0.5"), self.cell("E", "0", "")]
        runs = lay(self.tmp.name, {"synthia_s1": []}, cells=cells)
        rc, out = go(self.tmp.name, runs)
        got = {(r["class"], r["statistic"]): r for r in c.read_tsv(out / "m2-reading.tsv")}
        self.assertEqual([got[(k, "point")]["holm_p"] for k in "ABCDE"],
                         [f"{x:.12f}" for x in (0.004, 0.03, 0.08, 0.5)] + [""])
        self.assertEqual([got[(k, "point")]["moves"] for k in "ABCDE"], ["1", "1", "0", "0", ""])
        s = {(r["guard"], r["statistic"]): r for r in c.read_tsv(out / "m2-summary.tsv")}[("bioproject", "point")]
        self.assertEqual((s["cells"], s["pass"], s["move"], s["decides"]), ("5", "4", "2", "1"))

    def test_guards_are_beside(self):
        runs = lay(self.tmp.name, {"synthia_s1": []}, cells=[self.cell("A", "1", "0.001", guard="largest_out")])
        rc, out = go(self.tmp.name, runs)
        r = c.read_tsv(out / "m2-reading.tsv")[0]
        self.assertEqual((r["decides"], r["moves"]), ("0", "1"))


class Link(unittest.TestCase):
    def test_spearman_over_the_time_split_datasets(self):
        with tempfile.TemporaryDirectory() as tmp:
            ts = [{"dataset": f"synthia_s1__d{i}", "species_slug": "synthia_s1", "reading": "point_as_one",
                   "share_unseen": f"{i / 10:.12f}"} for i in range(5)]
            con = [{"dataset": f"synthia_s1__d{i}", "learner": "scm", "contrast": "total", "d_recall": str(-i),
                    "d_log_loss": str(i * i)} for i in range(5)]
            con.append({"dataset": "synthia_s1__d0", "learner": "scm", "contrast": "first", "d_recall": "9",
                        "d_log_loss": "9"})
            rc, out = go(tmp, lay(tmp, {"synthia_s1": []}, ts=ts), con)
            got = {r["measure"]: r for r in c.read_tsv(out / "timesplit-link.tsv")}
            self.assertEqual((got["d_recall"]["datasets"], got["d_recall"]["spearman"]), ("5", "-1.000000000000"))
            self.assertEqual(got["d_log_loss"]["spearman"], "1.000000000000")


class Shared(unittest.TestCase):
    def test_sign_test_holm_and_spearman_by_hand(self):
        self.assertEqual(C.sign_p_two_sided(0, 0), 1.0)
        self.assertAlmostEqual(C.sign_p_two_sided(5, 0), 2 / 32)
        self.assertAlmostEqual(C.sign_p_two_sided(3, 3), 1.0)
        self.assertEqual([round(x, 12) for x in C.holm([0.04, 0.01, 0.03])], [0.06, 0.03, 0.06])
        self.assertIsNone(C.spearman([1, 2], [1, 2]))
        self.assertIsNone(C.spearman([1, 1, 1], [1, 2, 3]))
        self.assertAlmostEqual(C.spearman([1, 2, 3, 4], [1, 3, 2, 4]), 0.8)

    def test_counts_are_written_as_they_are_and_recall_needs_ten(self):
        self.assertEqual([C.shown(n) for n in (0, 1, 4, 5)], ["0", "1", "4", "5"])
        self.assertEqual(C.recall_cell(9, 9), "")
        self.assertEqual(C.recall_cell(5, 10), "0.500000000000")
        with self.assertRaises(TypeError):
            C.shown(5.0)
        with self.assertRaises(TypeError):
            C.shown(True)
        with self.assertRaises(ValueError):
            C.shown(-1)

    def test_species_interval_draws_species_from_the_seed(self):
        by = {f"s{i}": [float(i)] for i in range(8)}
        a = C.species_interval(by, "m1|lineage|primary|scm|gene_only-point_only")
        self.assertEqual(a, C.species_interval(by, "m1|lineage|primary|scm|gene_only-point_only"))
        self.assertLessEqual(a[0], a[1])
        self.assertEqual(C.species_interval({f"s{i}": [1.0] for i in range(5)}, "x"), (None, None))


    def test_species_interval_ends_are_the_101st_and_3900th_of_4000(self):
        # the plan's words drawn again by a plain loop: resamples of the sorted
        # species, seeded by the row's name, medians to 12 decimals, the ends
        # at the 101st and 3,900th of 4,000 (the 2nd and 39th of 40)
        by = {f"s{i:02d}": [(i * 37 % 101) / 101 + i / 10 ** 6] for i in range(60)}
        name = "m1|lineage|primary|scm|gene_only-point_only"

        def medians(draws):
            rng = random.Random(int(hashlib.sha256(name.encode("utf-8")).hexdigest(), 16))
            keys = sorted(by)
            out = []
            for _ in range(draws):
                vals = [v for _ in keys for v in by[keys[rng.randrange(len(keys))]]]
                out.append(round(statistics.median(vals), 12))
            return sorted(out)

        meds = medians(4000)
        self.assertEqual(C.species_interval(by, name), (meds[100], meds[3899]))
        few = medians(40)
        self.assertNotEqual((few[0], few[39]), (few[1], few[38]))
        self.assertEqual(C.species_interval(by, name, resamples=40), (few[1], few[38]))

if __name__ == "__main__":
    unittest.main()
