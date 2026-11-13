"""guidance.py (scripts/release-summary) on two synthetic releases.

The first is built here: 14 datasets whose rows sit on the edges of the guidance's wording. The second is the synthetic
release of test_synthetic.py, on which guidance.py must read the rows and give the screen classes uncertainty.py gives.
Each expected value below is worked out by hand from the tables in this file; the comment beside it shows how.
"""

import json
import pathlib
import tempfile
import unittest

from tests.release_summary import test_synthetic as base

TB = "mycobacterium_tuberculosis"
# (x, n, lo, hi) of VME and of ME. The bounds are set by hand: the script reads the served bounds and recomputes none.
VME_OK = (2, 200, 0.001, 0.04)
ME_OK = (2, 200, 0.001, 0.036)  # 1%
VME_AT_UPPER = (5, 200, 0.01, 0.075)

# id, species, resistant genomes, susceptible genomes, then the VME and ME rows of the set covering machine. Every
# dataset's best learner is the set covering machine, except LASSO_BEST's.
DATASETS = [
    ("alpha__upper_at", "alpha", 200, 200, VME_AT_UPPER, ME_OK),  # upper bound exactly 7.5%: does not meet
    ("alpha__upper_under", "alpha", 200, 200, (5, 200, 0.01, 0.0749999), ME_OK),  # just under 7.5%: meets
    ("alpha__lower_at", "alpha", 200, 200, (8, 200, 0.015, 0.05), ME_OK),  # lower bound exactly 1.5%: does not meet
    ("alpha__lower_under", "alpha", 200, 200, (8, 200, 0.0149999, 0.05), ME_OK),  # just under 1.5%: meets
    ("beta__me_at", "beta", 200, 1000, VME_OK, (30, 1000, 0.02, 0.043)),  # ME exactly 3%: does not meet
    ("beta__me_under", "beta", 200, 1000, VME_OK, (29, 1000, 0.019, 0.041)),  # 2.9%, one genome fewer: meets
    ("beta__no_vme", "beta", 200, 200, None, ME_OK),  # no VME rate: not read
    ("beta__no_me", "beta", 200, 200, VME_OK, None),  # no ME rate: not read
    (f"{TB}__drug_a", TB, 200, 200, VME_OK, ME_OK),  # meets, and is left out without M. tuberculosis
    (f"{TB}__drug_b", TB, 200, 200, VME_AT_UPPER, ME_OK),  # does not meet
    # zero misses: the exact upper bound 1 - 0.025 ** (1 / n) is 0.0755 with 47 resistant genomes, 0.0740 with 48
    ("gamma__resistant_47", "gamma", 47, 200, (0, 47, 0.0, 1 - 0.025 ** (1 / 47)), ME_OK),  # does not meet
    ("gamma__resistant_48", "gamma", 48, 200, (0, 48, 0.0, 1 - 0.025 ** (1 / 48)), ME_OK),  # meets
    ("gamma__no_call", "gamma", 200, 200, (2, 190, 0.001, 0.042), ME_OK),  # 190 of 200 resistant counted: meets
    ("delta__lasso_best", "delta", 200, 200, VME_AT_UPPER, ME_OK),  # does not meet; its best learner's row does
]
LASSO_BEST = "delta__lasso_best"  # the lasso is best here; it has rows on the random split only, and only here
NO_LINEAGE = "delta__lasso_best"  # no lineage split; every other dataset has the same rows on both splits


def rate(row):
    return None if row is None else dict(zip(("x", "n", "lo", "hi"), row))


def build_edges(folder, revision=14, method="clopper-pearson", confidence=0.95):
    """Write index.json and intervals/ for the edge release into folder."""
    folder = pathlib.Path(folder)
    (folder / "intervals").mkdir(parents=True)
    datasets = []
    for ident, species, n_res, n_sus, vme, me in DATASETS:
        rows = {"random": [{"learner": "scm", "vme": rate(vme), "me": rate(me)}]}
        if ident == LASSO_BEST:
            rows["random"].append({"learner": "lasso", "vme": rate(VME_OK), "me": rate(ME_OK)})
        if ident != NO_LINEAGE:
            rows["lineage"] = rows["random"]
        datasets.append({"id": ident, "species_slug": species, "n": n_res + n_sus, "n_resistant": n_res,
                         "best_learner": "lasso" if ident == LASSO_BEST else "scm"})
        intervals = {"dataset": ident, "evaluation_revision": revision, "method": method, "confidence": confidence,
                     "splits": [{"scheme": s, "learners": r} for s, r in rows.items()]}
        (folder / "intervals" / f"{ident}.json").write_text(json.dumps(intervals))
    index = {"evaluation_revision": revision, "generated_at": "2026-01-01T00:00:00Z", "datasets": datasets}
    (folder / "index.json").write_text(json.dumps(index))
    return folder


def parts(row):
    return {k: row[k] for k in ("n", "vme_upper", "vme_lower", "vme", "me", "meets_wording")}


class Edges(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.edges = build_edges(pathlib.Path(cls.tmp.name) / "edges")
        cls.out = base.run_json("guidance.py", cls.edges)
        cls.rows = cls.out["rows"]

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_the_criteria_are_stated_as_numbers(self):
        c = self.out["criteria"]
        self.assertEqual(c["wording"]["vme_upper_below"], 0.075)
        self.assertEqual(c["wording"]["vme_lower_below"], 0.015)
        self.assertEqual(c["wording"]["me_rate_below"], 0.03)
        self.assertIs(c["wording"]["strict"], True)
        self.assertEqual((c["screen"]["vme_upper_at_most"], c["screen"]["me_upper_at_most"]), (0.015, 0.03))
        self.assertEqual(c["without_tuberculosis"], TB)
        # 1 - 0.025 ** (1 / 47) is 0.0755 and 1 - 0.025 ** (1 / 48) is 0.0740: Table 8 starts at 48
        self.assertEqual(c["table_8_first_resistant"], 48)
        self.assertEqual((self.out["revision"], self.out["datasets"]), (14, 14))

    def test_each_part_at_its_edge(self):
        # 14 rows, less beta__no_vme and beta__no_me: 12 read.
        # Upper bound not under 7.5%: alpha__upper_at, drug_b, gamma__resistant_47, delta (4), so 8 meet that part.
        # Lower bound not under 1.5%: alpha__lower_at alone, so 11. Both VME parts: 12 - 4 - 1 = 7.
        # ME not under 3%: beta__me_at alone, so 11. All three: the 7 less beta__me_at = 6.
        self.assertEqual(parts(self.rows["scm:random"]),
                         {"n": 12, "vme_upper": 8, "vme_lower": 11, "vme": 7, "me": 11, "meets_wording": 6})

    def test_the_screen_beside_it(self):
        # No VME upper bound is at or under 1.5%, and no lower bound is over 1.5% (VME) or 3% (ME): none is decided
        self.assertEqual(self.rows["scm:random"]["screen"], {"shown_to_meet_both": 0, "shown_to_miss_one": 0, "not_decided": 12})

    def test_a_row_without_both_rates_is_not_read(self):
        listed = {e["dataset"] for split in ("random", "lineage") for e in self.out["best_meets_wording"][split]}
        self.assertNotIn("beta__no_vme", listed)
        self.assertNotIn("beta__no_me", listed)
        # learners with no row at all are read on no dataset
        for who in ("cart", "rf", "gbt", "catalog"):
            for split in ("random", "lineage"):
                self.assertEqual(parts(self.rows[f"{who}:{split}"]),
                                 {"n": 0, "vme_upper": 0, "vme_lower": 0, "vme": 0, "me": 0, "meets_wording": 0})
        self.assertEqual(list(self.rows), [f"{w}:{s}" for s in ("random", "lineage")
                                           for w in ("best", "scm", "cart", "rf", "gbt", "lasso", "catalog")])

    def test_tuberculosis_is_left_out_of_the_second_total(self):
        # drug_a meets and drug_b does not: 12 - 2 read, 6 - 1 meet
        self.assertEqual(self.rows["scm:random"]["without_tuberculosis"], {"n": 10, "meets_wording": 5})
        self.assertEqual(self.rows["best:random"]["without_tuberculosis"], {"n": 10, "meets_wording": 6})

    def test_table_8_and_the_rows_short_of_every_genome(self):
        r = self.rows["scm:random"]
        self.assertEqual(r["vme_n_under_table_8"], 1)  # gamma__resistant_47; 48 is not under
        self.assertEqual(r["fewer_than_every_genome"], 1)  # gamma__no_call counts 190 of 200

    def test_the_best_learner_and_its_list(self):
        # as the set covering machine, but on delta the lasso's row, which meets every part
        self.assertEqual(parts(self.rows["best:random"]),
                         {"n": 12, "vme_upper": 9, "vme_lower": 11, "vme": 8, "me": 11, "meets_wording": 7})
        self.assertEqual(parts(self.rows["lasso:random"]),
                         {"n": 1, "vme_upper": 1, "vme_lower": 1, "vme": 1, "me": 1, "meets_wording": 1})
        listed = self.out["best_meets_wording"]["random"]
        self.assertEqual([(e["dataset"], e["learner"]) for e in listed], [
            ("alpha__upper_under", "scm"), ("alpha__lower_under", "scm"), ("beta__me_under", "scm"),
            (f"{TB}__drug_a", "scm"), ("gamma__resistant_48", "scm"), ("gamma__no_call", "scm"),
            (LASSO_BEST, "lasso")])
        # the served bounds and the rate are not rounded: 0.0749999 does not print as 0.075
        self.assertEqual(listed[0]["vme"], {"x": 5, "n": 200, "lo": 0.01, "hi": 0.0749999})
        self.assertEqual(listed[2]["me"], {"x": 29, "n": 1000, "rate": 29 / 1000})
        self.assertEqual(listed[4]["vme"], {"x": 0, "n": 48, "lo": 0.0, "hi": 1 - 0.025 ** (1 / 48)})

    def test_the_lineage_split(self):
        # the random rows without delta, which has no lineage split. The set covering machine's row on delta fails the
        # upper bound and passes the lower bound and ME, so n, vme_lower and me each drop by one.
        self.assertEqual(parts(self.rows["scm:lineage"]),
                         {"n": 11, "vme_upper": 8, "vme_lower": 10, "vme": 7, "me": 10, "meets_wording": 6})
        self.assertEqual(parts(self.rows["best:lineage"]), parts(self.rows["scm:lineage"]))
        self.assertEqual(self.rows["lasso:lineage"]["n"], 0)
        self.assertEqual([e["dataset"] for e in self.out["best_meets_wording"]["lineage"]],
                         [e["dataset"] for e in self.out["best_meets_wording"]["random"][:6]])

    def test_revision_12_is_read_the_same_way(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = base.run_json("guidance.py", build_edges(pathlib.Path(tmp) / "rev12", 12))
        self.assertEqual(out["revision"], 12)
        self.assertEqual(out["rows"], self.rows)
        self.assertEqual(out["best_meets_wording"], self.out["best_meets_wording"])

    def test_unknown_revision_is_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = build_edges(pathlib.Path(tmp) / "rev13", 13)
            p = base.run("guidance.py", folder / "index.json", folder / "intervals", check=False)
        self.assertEqual((p.returncode, p.stdout), (1, ""))
        self.assertIn("evaluation revision 13", p.stderr)
        self.assertIn("understands revisions 12, 14", p.stderr)

    def test_other_intervals_are_refused(self):
        for name, kwargs in (("wilson", {"method": "wilson"}), ("ninety", {"confidence": 0.9})):
            with self.subTest(name), tempfile.TemporaryDirectory() as tmp:
                folder = build_edges(pathlib.Path(tmp) / name, **kwargs)
                p = base.run("guidance.py", folder / "index.json", folder / "intervals", check=False)
                self.assertEqual((p.returncode, p.stdout), (1, ""))
                self.assertIn("clopper-pearson 95% intervals only", p.stderr)


class OnTheSyntheticRelease(unittest.TestCase):
    """The release of test_synthetic.py, in revisions 12 and 14."""

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        root = pathlib.Path(cls.tmp.name)
        cls.out, cls.unc = {}, {}
        for revision in (12, 14):
            folder = base.build(root / f"rev{revision}", revision)
            cls.out[revision] = base.run_json("guidance.py", folder)
            cls.unc[revision] = base.run_json("uncertainty.py", folder, "--reps", "20")

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_the_rows_and_screen_classes_of_uncertainty(self):
        for revision in (12, 14):
            bars = self.unc[revision]["error_bars"]["rows"]
            self.assertEqual(list(self.out[revision]["rows"]), list(bars))
            for key, r in self.out[revision]["rows"].items():
                with self.subTest(revision=revision, row=key):
                    self.assertEqual(r["n"], bars[key]["n"])
                    self.assertEqual(r["screen"], {k: bars[key][k] for k in r["screen"]})

    def test_the_wording_reads_looser_than_the_screen(self):
        r = self.out[14]["rows"]["scm:random"]
        # drug1 (VME lower bound 0.0150001), drug2 (ME 10/200) and drug4 (lower bound 0.015, ME 7/200) do not meet;
        # drug0, drug3 and the 7 default rows (VME 2/200 in 0.001 to 0.04, ME 4/200) do. One is shown to meet the screen.
        self.assertEqual(parts(r), {"n": 12, "vme_upper": 12, "vme_lower": 10, "vme": 10, "me": 10, "meets_wording": 9})
        self.assertEqual(r["screen"]["shown_to_meet_both"], 1)
        self.assertEqual(r["without_tuberculosis"], {"n": 12, "meets_wording": 9})

    def test_rows_short_of_every_genome_by_revision(self):
        # cart, lineage, gamma__drug11 leaves out 2 resistant no-calls in both revisions. The lasso's lineage row on
        # gamma__drug10 leaves out its 4 + 6 no-calls in revision 12 and is flagged and scored on every genome in 14.
        for revision, lasso in ((12, 1), (14, 0)):
            rows = self.out[revision]["rows"]
            with self.subTest(revision=revision):
                self.assertEqual(rows["cart:lineage"]["fewer_than_every_genome"], 1)
                self.assertEqual(rows["lasso:lineage"]["fewer_than_every_genome"], lasso)
                self.assertEqual(sum(r["fewer_than_every_genome"] for r in rows.values()), 1 + lasso)


if __name__ == "__main__":
    unittest.main()
