"""The designed release of made_up2.py, 8 made-up species, through stage 1 and stage 2, against figures worked by
hand: one exact tie, one fold with a null fit, one dataset under 10 pairs, one report without a catalogue rule
and one row lacking an MLST type. Then two releases of one species for the typing reads, the designed release
with a ninth species and a report left out for the wholes, and its first 5 species, which give an estimate and
no interval.

Every value here is made up: the species, the rows and their scores (made_up2.designed())."""

import csv
import fractions
import io
import unittest

import numpy

from tests.matched_accuracy import made_up2, stage2_common as c

F = fractions.Fraction


def mean(values):
    return sum(values, F(0)) / len(values)


# Each group's value worked by hand from made_up2.designed(): half points over 2 x pairs. {species: [[group
# values of one dataset]]}, reading A.
A = {
    "synthia_kappa": [[F(5, 6), F(7, 8)]],         # group 2: R 0.6 ties S 0.6, R 0.7 beats it; both beat S 0.3
    "synthia_lambda": [[F(1)], [F(2, 3), F(3, 4)]],
    "synthia_mu": [[F(3, 4), F(0)]],               # 6 pairs, under 10
    "synthia_nu": [[F(5, 6)]],                     # the report without a catalogue rule
    "synthia_xi": [[F(5, 6), F(1, 2)]],            # group 2 sits in fold 4, where gbt has a null fit
    "synthia_omicron": [[F(1), F(1, 2)]],          # the S row at 0.7 lacks an MLST type
    "synthia_rho": [[F(5, 6), F(1, 2)], [F(3, 4)]],
    "synthia_tau": [[F(0), F(1)]],
}
# Reading B leaves out kappa's S row at 0.3, which carries an I record: group 2 keeps the tie and one win.
B = dict(A, synthia_kappa=[[F(5, 6), F(3, 4)]])


def species_values(plan):
    return {sp: mean([mean(d) for d in ds]) for sp, ds in plan.items()}


def overall(plan):
    return mean(list(species_values(plan).values()))


def ppm(x):
    return round(x * 1000000)


def interval_by_hand(plan):
    """numpy's PCG64 at seed 42 over the hand species values in sorted name order, linear percentiles; the bounds in
    parts per million and the word read from the bounds before they are rounded."""
    values = numpy.array([float(v) for _, v in sorted(species_values(plan).items())])
    rng = numpy.random.Generator(numpy.random.PCG64(42))
    draws = values[rng.integers(0, len(values), size=(10000, len(values)))].mean(axis=1)
    lower, upper = (float(b) for b in numpy.percentile(draws, [2.5, 97.5]))
    word = "above" if lower > 0.5 else "below" if upper < 0.5 else "unresolved"
    return ppm(F(lower)), ppm(F(upper)), word


class Designed(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.made = c.Run()

    @classmethod
    def tearDownClass(cls):
        cls.made.close()

    def row(self, file, reading, level, name, **match):
        return next(r for r in self.made.table(file) if (r["reading"], r["level"], r["name"]) == (reading, level, name)
                    and all(r[k] == v for k, v in match.items()))

    def test_the_run_passes_and_prints_both_readings(self):
        self.assertEqual(self.made.code, 0, self.made.lines)
        self.assertEqual(self.made.lines[0], "scores opened: yes")
        self.assertEqual(self.made.lines[1], "reading A: estimate 692708 ppm over 17 groups, 102 pairs, 10 datasets "
                                             "and 8 species; interval 572917 to 796875 ppm over 8 species and 10000 "
                                             "draws")
        self.assertTrue(self.made.lines[2].startswith("reading B: estimate 684896 ppm over 17 groups, 100 pairs"))

    def test_the_overall_figure_worked_by_hand(self):
        # kappa 41/48, lambda (1 + 17/24)/2 = 41/48, mu 3/8, nu 5/6, xi 2/3, omicron 3/4, rho (2/3 + 3/4)/2 =
        # 17/24, tau 1/2. Over 48: 41 + 41 + 18 + 40 + 32 + 36 + 34 + 24 = 266, so 266/48 / 8 = 133/192.
        self.assertEqual(overall(A), F(133, 192))
        self.assertEqual(ppm(F(133, 192)), 692708)
        self.assertEqual(self.made.measure("A", "overall", "all")["estimate_ppm"], "692708")
        # Reading B: kappa (5/6 + 3/4)/2 = 19/24, 3/48 lower: 263/384.
        self.assertEqual(overall(B), F(263, 384))
        self.assertEqual(self.made.measure("B", "overall", "all")["estimate_ppm"], str(ppm(F(263, 384))))

    def test_each_species_and_dataset_value_is_the_hand_ones(self):
        for name, plan in (("A", A), ("B", B)):
            for sp, v in species_values(plan).items():
                with self.subTest(reading=name, species=sp):
                    got = self.made.measure(name, "species", sp)
                    want = "" if int(got["pairs"]) < 10 else str(ppm(v))
                    self.assertEqual(got["estimate_ppm"], want)

    def test_the_one_tie_counts_one_half(self):
        # Counted 1, kappa's group 2 would be 1; counted 0, 3/4. One half gives 7/8 and kappa 41/48.
        self.assertEqual(self.made.measure("A", "dataset", "synthia_kappa__drug_a")["estimate_ppm"],
                         str(ppm(mean([F(5, 6), F(7, 8)]))))
        self.assertEqual(ppm(F(41, 48)), 854167)

    def test_a_dataset_under_10_pairs_gives_counts_and_no_estimate_and_enters_every_mean(self):
        for level, name in (("dataset", "synthia_mu__drug_a"), ("species", "synthia_mu")):
            got = self.made.measure("A", level, name)
            self.assertEqual((got["estimate_ppm"], got["pairs"], got["groups"], got["rows"]), ("", "6", "2", "7"))
        # The overall figure of the test above counts mu's 3/8; without it, the mean of 7 species differs.
        without = mean([v for sp, v in species_values(A).items() if sp != "synthia_mu"])
        self.assertNotEqual(ppm(without), 692708)
        # In reading B kappa falls to 8 pairs: counts, and no estimate.
        got = self.made.measure("B", "dataset", "synthia_kappa__drug_a")
        self.assertEqual((got["estimate_ppm"], got["pairs"]), ("", "8"))

    def test_the_null_fit_read_leaves_out_xis_fold_4_group(self):
        # xi falls to 5/6: 266 - 32 + 40 = 274, over 384, 137/192.
        plan = dict(A, synthia_xi=[[F(5, 6)]])
        self.assertEqual(overall(plan), F(137, 192))
        got = self.row("null_fits.tsv", "A", "overall", "all")
        self.assertEqual(got["estimate_ppm"], str(ppm(F(137, 192))))
        self.assertEqual([got[k] for k in ("lost_rows", "lost_groups", "lost_pairs", "lost_datasets", "lost_species")],
                         ["4", "1", "4", "0", "0"])
        xi = self.row("null_fits.tsv", "A", "species", "synthia_xi")
        self.assertEqual((xi["estimate_ppm"], xi["pairs"]), ("", "6"))

    def test_the_four_cells_and_the_pairs_counted_apart(self):
        # Resistant member called resistant, susceptible member susceptible: every pair but rho's 8 and nu's 12.
        # By hand: kappa 41/48, lambda 41/48, mu 3/8, xi 2/3, omicron 3/4, rho 1 (its R row called susceptible,
        # its S row called resistant and its unevaluable S row each leave the cell), tau 1/2; nu has no drug rule.
        cell = {"synthia_kappa": F(41, 48), "synthia_lambda": F(41, 48), "synthia_mu": F(3, 8), "synthia_xi": F(2, 3),
                "synthia_omicron": F(3, 4), "synthia_rho": F(1), "synthia_tau": F(1, 2)}
        self.assertEqual(mean(list(cell.values())), F(5, 7))
        rows = {(r["reading"], r["resistant_member_call"], r["susceptible_member_call"]):
                (r["estimate_ppm"], r["pairs"], r["groups"], r["datasets"], r["species"])
                for r in self.made.table("cells.tsv")}
        self.assertEqual(rows[("A", "resistant", "susceptible")], (str(ppm(F(5, 7))), "82", "16", "9", "7"))
        # rho's S row at 0.9, called resistant, beats both R rows: 0 over 2 pairs, which a cell prints.
        self.assertEqual(rows[("A", "resistant", "resistant")], ("0", "2", "1", "1", "1"))
        # rho's R row at 0.4, called susceptible: it loses to 0.5 and beats 0.3 and 0.2, 2/3.
        self.assertEqual(rows[("A", "susceptible", "susceptible")], (str(ppm(F(2, 3))), "3", "1", "1", "1"))
        self.assertEqual(rows[("A", "susceptible", "resistant")], ("", "0", "0", "0", "0"))
        self.assertEqual(rows[("A", "unevaluable", "unevaluable")], ("", "3", "1", "1", "1"))
        self.assertEqual(rows[("A", "no drug rule", "no drug rule")], ("", "12", "1", "1", "1"))
        self.assertEqual(sum(int(v[1]) for k, v in rows.items() if k[0] == "A"), 102)

    def test_the_mlst_read_leaves_out_the_row_lacking_a_type(self):
        # MLST match: kappa (1 + 1 + 7/8)/3 = 23/24, ST-2 splitting group 1; lambda's dataset 1 splits in ST-1 (1)
        # and ST-2 (1), its ST-3 rows holding no R row; omicron's group 1 loses its untyped S row at 0.7 and stays
        # 1, and its group 2 loses its one R row of ST-2, leaving R 0.5 against S 0.4 and 0.3, 1.
        mlst = dict(A, synthia_kappa=[[F(1), F(1), F(7, 8)]], synthia_lambda=[[F(1), F(1)], [F(2, 3), F(3, 4)]],
                    synthia_omicron=[[F(1), F(1)]])
        self.assertEqual(overall(mlst), F(283, 384))
        got = next(r for r in self.made.table("typing.tsv") if (r["reading"], r["step"], r["match"], r["level"]) ==
                   ("A", "mlst", "mlst", "overall"))
        self.assertEqual(got["estimate_ppm"], str(ppm(F(283, 384))))
        wider = next(r for r in self.made.table("typing.tsv") if (r["reading"], r["step"], r["match"], r["level"]) ==
                     ("A", "mlst", "primary", "overall"))
        # The same 78 rows by the primary match: 4 rows lost, the untyped S row, omicron's lone ST-2 R row and
        # lambda's two ST-3 S rows.
        self.assertEqual((wider["rows"], wider["lost_rows"], wider["lost_pairs"]), ("78", "4", "10"))
        primary = dict(A, synthia_omicron=[[F(1), F(1)]],
                       synthia_lambda=[[F(1)], [F(2, 3), F(3, 4)]])
        self.assertEqual(wider["estimate_ppm"], str(ppm(overall(primary))))

    def test_the_lost_columns_stand_on_the_primary_match_rows_alone(self):
        # On a row of the narrower match the lost columns are empty; on the wider primary match, integers.
        rows = self.made.table("typing.tsv")
        self.assertEqual({r["match"] for r in rows}, {"mlst", "mlst_hc20", "primary"})
        for r in rows:
            lost = [r[k] for k in c.stage2.LOST]
            with self.subTest(step=r["step"], match=r["match"], level=r["level"], name=r["name"]):
                if r["match"] == "primary":
                    self.assertTrue(all(x.isdigit() for x in lost), lost)
                else:
                    self.assertEqual(lost, [""] * 5)
        mlst = next(r for r in rows if (r["reading"], r["step"], r["match"], r["level"]) == ("A", "mlst", "mlst",
                                                                                            "overall"))
        self.assertEqual((mlst["rows"], mlst["lost_rows"], mlst["lost_pairs"]), ("78", "", ""))

    def test_each_interval_is_numpys_over_the_hand_species_values(self):
        got = {(r["reading"], r["read"], r["match"]):
               (int(r["lower_ppm"]), int(r["upper_ppm"]), r["chance"], r["species"], r["draws"])
               for r in self.made.table("interval.tsv")}
        self.assertEqual(got[("A", "primary", "primary")], interval_by_hand(A) + ("8", "10000"))
        self.assertEqual(got[("B", "primary", "primary")], interval_by_hand(B) + ("8", "10000"))
        self.assertEqual(got[("A", "null_fits", "primary")][:3], interval_by_hand(dict(A, synthia_xi=[[F(5, 6)]])))
        # Reading A's lower bound, 572917 parts per million, is above one half.
        self.assertEqual(got[("A", "primary", "primary")][:3], (572917, 796875, "above"))


def typing_row(made, step, match, level="overall"):
    return next(r for r in made.table("typing.tsv") if (r["reading"], r["step"], r["match"], r["level"]) ==
                ("A", step, match, level))


COLUMNS = ("estimate_ppm", "rows", "groups", "pairs", "datasets", "species") + c.stage2.LOST


class TypingRows(unittest.TestCase):
    """Two made-up releases of one species: a group with MLST, no HC20 and both outcomes, and MLST values with
    blanks around them."""

    def test_a_group_with_mlst_no_hc20_and_both_outcomes_leaves_the_mlst_and_hc20_read(self):
        plan = [[made_up2.dataset("drug_a", [
            made_up2.group(0, 10, [made_up2.R(0.9), made_up2.R(0.6), made_up2.S(0.5), made_up2.S(0.2)]),  # 1
            made_up2.group(1, 11, [made_up2.R(0.8, hc20=None), made_up2.R(0.2, hc20=None),
                                   made_up2.S(0.3, hc20=None), made_up2.S(0.7, hc20=None)]),             # 1/2
        ])]]
        made = c.Run(plan)
        self.addCleanup(made.close)
        self.assertEqual(made.code, 0, made.lines)
        # The MLST read holds both groups, (1 + 1/2) / 2 = 3/4; the MLST and HC20 read holds the first alone.
        self.assertEqual([typing_row(made, "mlst", "mlst")[k] for k in COLUMNS],
                         ["750000", "8", "2", "8", "1", "1", "", "", "", "", ""])
        self.assertEqual([typing_row(made, "mlst_hc20", "mlst_hc20")[k] for k in COLUMNS],
                         ["1000000", "4", "1", "4", "1", "1", "", "", "", "", ""])
        # Its wider primary match runs on the same 4 rows and loses the other 4 rows, 1 group and 4 pairs.
        self.assertEqual([typing_row(made, "mlst_hc20", "primary")[k] for k in COLUMNS],
                         ["1000000", "4", "1", "4", "1", "1", "4", "1", "4", "0", "0"])

    def test_an_mlst_value_with_blanks_around_it_is_read_trimmed_and_blanks_alone_lack_a_value(self):
        plan = [[made_up2.dataset("drug_a", [
            made_up2.group(0, 10, [made_up2.R(0.9, mlst=" ST-2 "), made_up2.R(0.4, mlst="ST-2"),
                                   made_up2.S(0.5, mlst="ST-2"), made_up2.S(0.1, mlst="   ")]),
        ])]]
        made = c.Run(plan)
        self.addCleanup(made.close)
        self.assertEqual(made.code, 0, made.lines)
        # " ST-2 " joins ST-2: R 0.9 beats S 0.5 and R 0.4 does not, 1/2 over 2 pairs; the S row of blanks leaves.
        self.assertEqual([typing_row(made, "mlst", "mlst")[k] for k in COLUMNS[:6]], ["500000", "3", "1", "2", "1", "1"])
        self.assertEqual(made.measure("A", "overall", "all")["estimate_ppm"], "750000")


class Wholes(unittest.TestCase):
    """inputs.tsv gives the genome rows and the species of the reports that enter, counted before the line. The
    designed release gains a ninth species whose report enters with no eligible row, and a report of a tenth
    species that holds none of the five learners and is left out."""

    def test_the_rows_and_species_of_the_reports_that_enter(self):
        def before(release):
            rows = [{"id": made_up2.made_up.gid("ch", n), "resistant": n % 2 == 0, "fold": n % 5, "lineage": 30 + n,
                     "mlst": "ST-1", "from_mic": False} for n in range(1, 8)]
            release["reports"].append(made_up2.made_up.report("synthia_chi__drug_a", "Synthia chi", "synthia_chi",
                                                              "drug_a", rows, None, learners=()))
        plan = made_up2.designed() + [[made_up2.dataset("drug_a", [
            made_up2.group(0, 10, [made_up2.R(0.5)]), made_up2.group(1, 11, [made_up2.S(0.5)])])]]
        made = c.Run(plan, before=before)
        self.addCleanup(made.close)
        self.assertEqual(made.code, 0, made.lines)
        entering = [d for d in made.release["reports"]
                    if next(s for s in d["splits"] if s["scheme"] == "lineage")["learners"]]
        self.assertEqual(len(entering), 11)
        want = sum(len(d["samples"]) for d in entering)
        self.assertEqual(sum(len(d["samples"]) for d in made.release["reports"]), want + 7)
        # Stage 1's count table gives the same whole, as its dataset rows of the reports that enter, which on the
        # real release add to 215,319.
        counts = list(csv.DictReader(io.StringIO(made.counts().read_text(encoding="utf-8")), delimiter="\t"))
        self.assertEqual(sum(int(r["genome_rows"]) for r in counts
                             if (r["level"], r["reading"], r["reports_left_out"]) == ("dataset", "A", "0")), want)
        rows = {r["input"]: r for r in made.table("inputs.tsv")}
        self.assertEqual(rows["genome rows of the reports that enter"]["count"], str(want))
        self.assertEqual(rows["species of the reports that enter"]["count"], "9")
        # 8 species enter the measure; the ninth enters with no eligible row.
        self.assertEqual(made.measure("A", "overall", "all")["species"], "8")


class FiveSpecies(unittest.TestCase):
    """The first 5 species of the designed release: an estimate, and no interval."""

    @classmethod
    def setUpClass(cls):
        cls.made = c.Run(made_up2.designed()[:5])

    @classmethod
    def tearDownClass(cls):
        cls.made.close()

    def test_5_species_give_the_estimate_and_no_interval(self):
        self.assertEqual(self.made.code, 0, self.made.lines)
        five = {sp: v for sp, v in A.items() if sp in ("synthia_kappa", "synthia_lambda", "synthia_mu", "synthia_nu",
                                                       "synthia_xi")}
        # Over 48: 41 + 41 + 18 + 40 + 32 = 172, so 172/240 = 43/60.
        self.assertEqual(overall(five), F(43, 60))
        self.assertEqual(self.made.measure("A", "overall", "all")["estimate_ppm"], str(ppm(F(43, 60))))
        # Groups 2 + 3 + 2 + 1 + 2 = 10; pairs 10 + 22 + 6 + 12 + 10 = 60.
        self.assertEqual(self.made.lines[1], f"reading A: estimate {ppm(F(43, 60))} ppm over 10 groups, 60 pairs, 6 "
                                             f"datasets and 5 species; no interval: 5 species, under 6")
        for r in self.made.table("interval.tsv"):
            self.assertEqual((r["lower_ppm"], r["upper_ppm"], r["chance"], r["species"], r["draws"]),
                             ("", "", "", "5", "0"))


if __name__ == "__main__":
    unittest.main()
