"""scripts/time-split/stagger_design.py on SYNTHETIC metadata: each rule of the frozen plan's Design items 1 to 8
and of R33 points 1 to 6, by values the test names itself, the variant rule step by step, both clocks, the 60% check,
the barred paths and the fold files."""

import datetime
import hashlib
import json
import os
import shutil
import tempfile
import unittest
from pathlib import Path

from tests.time_split import stagger_fixture as fx

sd = fx.design_script()

# A hand world of 13 cohort genomes: id, date, study, label, kset digest. Dates are whole numbers, as years are.
HAND = [
    ("g01", 1, "pa", "R", "d01"),
    ("g02", 2, "pa", "S", "d07"),    # the digest of g07, a test genome of window 2
    ("g03", 3, "pb", "R", "d03"),
    ("g04", 4, "pb", "S", "d04"),
    ("g05", 5, "pc", "R", "d11"),    # the digest of g11, a test genome of window 3
    ("g06", 6, "pd", "S", "d06"),
    ("g07", 7, "pe", "R", "d07"),
    ("g08", 8, "pb", "S", "d08"),
    ("g09", 9, None, "R", "d09"),    # no recorded study
    ("g10", 10, "pf", "S", "d10"),
    ("g11", 11, "pg", "R", "d11"),
    ("g12", 12, "pe", "S", "d12"),
    ("g13", None, "ph", "R", "d13"),  # undated
]
# Genomes of the species list outside the cohort: id, date, study.
LISTED_ONLY = [("x01", 9, "pf"), ("x02", None, "pg")]


def hand_plan(shares=(25, 50, 75), changes=None, clock="collection-year"):
    rows = {g: (d, s, l, k) for g, d, s, l, k in HAND}
    extra = {g: (d, s) for g, d, s in LISTED_ONLY}
    for g, row in (changes or {}).items():
        if g in rows:
            rows[g] = row
        else:
            extra[g] = row
    cohort = list(rows)
    date = {g: r[0] for g, r in rows.items()}
    date.update({g: r[0] for g, r in extra.items()})
    study = {g: r[1] for g, r in rows.items()}
    study.update({g: r[1] for g, r in extra.items()})
    labels = {g: r[2] for g, r in rows.items()}
    digests = {g: r[3] for g, r in rows.items()}
    listed = cohort + list(extra)
    return sd.plan_dataset("hand__drug", clock, shares, cohort, labels, date, study, listed, digests), labels


def window_with(test, arms):
    """A plan's window from (n R, n S) per set, for pairs_of: test and arms {(arm, k): (n R, n S)}."""
    labels = {}

    def make(name, counts):
        out = []
        for cls, n in zip("RS", counts):
            for i in range(n):
                g = f"syn-{name}-{cls}-{i:03d}"
                labels[g] = cls
                out.append(g)
        return out

    w = {"test": make("t", test), "copies": set(), "arms": {}}
    for (arm, k), counts in arms.items():
        w["arms"][(arm, k)] = make(f"{arm}{k}", counts)
    return w, labels


class CutPoints(unittest.TestCase):
    def test_ten_dates_at_forty_sixty_eighty(self):
        # ceil(0.4 x 10) = 4, ceil(0.6 x 10) = 6, ceil(0.8 x 10) = 8.
        self.assertEqual(sd.cut_points([10, 9, 8, 7, 6, 5, 4, 3, 2, 1], (40, 60, 80)), [4, 6, 8])

    def test_seven_dates_round_the_place_up(self):
        # ceil(2.8) = 3, ceil(4.2) = 5, ceil(5.6) = 6.
        self.assertEqual(sd.cut_points([1, 2, 3, 4, 5, 6, 7], (40, 60, 80)), [3, 5, 6])

    def test_days_as_text(self):
        days = ["2020-01-0%d" % i for i in range(1, 10)] + ["2020-01-10"]
        self.assertEqual(sd.cut_points(days, (50, 75)), ["2020-01-05", "2020-01-08"])

    def test_equal_cut_points_leave_empty_windows(self):
        dates = [1, 1, 1, 1, 1, 1, 1, 1, 2, 3]
        cuts = sd.cut_points(dates, (40, 60, 80))
        self.assertEqual(cuts, [1, 1, 1])
        self.assertEqual([sd.window_of(d, cuts) for d in (1, 2, 3)], [None, 3, 3])

    def test_a_variant_keeps_its_origins_whatever_its_cut_points(self):
        plan, _ = hand_plan(shares=(40, 60, 80), changes={g: (1, s, l, k) for g, _, s, l, k in HAND[:12]})
        self.assertEqual(plan["cuts"], [1, 1, 1])
        self.assertEqual(sorted(plan["windows"]), [2, 3])


class Windows(unittest.TestCase):
    def test_after_one_cut_up_to_the_next_the_last_to_the_end(self):
        cuts = [3, 6, 8]
        got = [sd.window_of(d, cuts) for d in (None, 2, 3, 4, 6, 7, 8, 9, 100)]
        self.assertEqual(got, [None, None, None, 1, 1, 2, 2, 3, 3])


class HandWorld(unittest.TestCase):
    """The hand world at shares 25, 50 and 75: 12 dated genomes give the cut points 3, 6 and 9."""

    def setUp(self):
        self.plan, self.labels = hand_plan()

    def test_cut_points(self):
        self.assertEqual(self.plan["cuts"], [3, 6, 9])
        self.assertEqual((self.plan["dated"], self.plan["undated"]), (12, 1))

    def test_test_set_of_each_window(self):
        self.assertEqual(self.plan["windows"][2]["test"], ["g07"])
        self.assertEqual(self.plan["windows"][3]["test"], ["g11"])

    def test_a_study_known_by_an_origin(self):
        # pf is carried by x01, dated 9, the cut point of window 3: known, so g10 is not tested.
        self.assertNotIn("g10", self.plan["windows"][3]["test"])
        # Dated after the cut point, x01 makes pf known no more, and g10 is tested.
        later, _ = hand_plan(changes={"x01": (10, "pf")})
        self.assertEqual(later["windows"][3]["test"], ["g10", "g11"])

    def test_a_genome_with_no_study_is_tested_nowhere(self):
        for w in self.plan["windows"].values():
            self.assertNotIn("g09", w["test"])

    def test_an_undated_genome_trains_nowhere_and_is_tested_nowhere(self):
        for w in self.plan["windows"].values():
            self.assertNotIn("g13", w["test"])
            for members in w["arms"].values():
                self.assertNotIn("g13", members)

    def test_an_undated_genome_makes_no_study_known(self):
        # pg is carried by x02, undated, alone: g11 of pg stays new and is tested.
        self.assertIn("g11", self.plan["windows"][3]["test"])

    def test_copies_per_window(self):
        w2, w3 = self.plan["windows"][2], self.plan["windows"][3]
        self.assertEqual(w2["copies"], {"g02"})
        self.assertEqual(w3["copies"], {"g05"})
        # g02 trains nowhere on window 2 and trains on window 3; g05 the other way round.
        self.assertNotIn("g02", w2["arms"][("all", 2)])
        self.assertIn("g02", w3["arms"][("all", 1)])
        self.assertIn("g05", w2["arms"][("all", 2)])
        self.assertNotIn("g05", w3["arms"][("all", 3)])

    def test_all_k(self):
        w2, w3 = self.plan["windows"][2], self.plan["windows"][3]
        self.assertEqual(w2["arms"][("all", 1)], ["g01", "g03"])
        self.assertEqual(w2["arms"][("all", 2)], ["g01", "g03", "g04", "g05", "g06"])
        self.assertEqual(w3["arms"][("all", 2)], ["g01", "g02", "g03", "g04", "g06"])
        # A test genome of window 2, g07, trains origin 3.
        self.assertEqual(w3["arms"][("all", 3)], ["g01", "g02", "g03", "g04", "g06", "g07", "g08", "g09"])

    def test_recent_k_takes_the_latest(self):
        w2, w3 = self.plan["windows"][2], self.plan["windows"][3]
        # Window 2: all(1) holds 2 R and 0 S.
        self.assertEqual(w2["arms"][("recent", 2)], ["g03", "g05"])
        # Window 3: all(1) holds 2 R and 1 S.
        self.assertEqual(w3["arms"][("recent", 2)], ["g01", "g03", "g06"])
        self.assertEqual(w3["arms"][("recent", 3)], ["g07", "g08", "g09"])

    def test_recent_k_holds_as_many_of_each_class_as_all_1(self):
        w2, w3 = self.plan["windows"][2], self.plan["windows"][3]
        self.assertEqual(sd.counts(w2["arms"][("recent", 2)], self.labels), (2, 0))
        self.assertEqual(sd.counts(w3["arms"][("recent", 2)], self.labels), (2, 1))
        self.assertEqual(sd.counts(w3["arms"][("recent", 3)], self.labels), (2, 1))

    def test_recent_1_is_all_1(self):
        for w in self.plan["windows"].values():
            self.assertEqual(w["arms"][("recent", 1)], w["arms"][("all", 1)])


class Ties(unittest.TestCase):
    def test_the_plans_keyed_hash_to_the_letter(self):
        self.assertEqual(sd.tie_key("public-day", "alpha_species__drug_one", "syn-a-0001"),
                         "0c33b9ecb219f4aaa44ec4c5686e81075d073dc85edbb59efe21592b642daf03")
        self.assertEqual(sd.tie_key("collection-year", "alpha_species__drug_one", "syn-a-0001"),
                         "46460bd99ed4543e88fbef1fbe924cf7dd5b7305c053f2cd6ba6501524034cd7")

    def test_a_tie_goes_to_the_smaller_hash_of_its_clock(self):
        pool = ["syn-t-1", "syn-t-2", "syn-t-3"]
        labels = {g: "R" for g in pool}
        date = {"syn-t-1": 5, "syn-t-2": 5, "syn-t-3": 3}
        for clock, winner in (("public-day", "syn-t-1"), ("collection-year", "syn-t-2")):
            got = sd.recent(pool, {"R": 1, "S": 0}, labels, date, lambda g: sd.tie_key(clock, "hand__drug", g))
            self.assertEqual(got, [winner], clock)


class Pairs(unittest.TestCase):
    def arms(self, older=(25, 25), fresh=(30, 30), arm="all"):
        return {(arm, 1): older, (arm, 2): fresh}

    def test_a_pair_at_the_minimums(self):
        w, labels = window_with((20, 20), self.arms())
        self.assertEqual(sd.pairs_of({"windows": {2: w}}, labels, 20, "stale-all"), [(1, 2)])

    def test_a_test_set_under_the_minimum_in_either_class(self):
        for test in ((19, 20), (20, 19)):
            w, labels = window_with(test, self.arms())
            self.assertEqual(sd.pairs_of({"windows": {2: w}}, labels, 20, "stale-all"), [], test)

    def test_a_training_set_under_25_in_either_class(self):
        for older, fresh in (((24, 25), (30, 30)), ((25, 24), (30, 30)), ((25, 25), (30, 24))):
            w, labels = window_with((20, 20), self.arms(older, fresh))
            self.assertEqual(sd.pairs_of({"windows": {2: w}}, labels, 20, "stale-all"), [], (older, fresh))

    def test_each_contrast_reads_its_own_arms(self):
        arms = {("all", 1): (25, 25), ("all", 2): (30, 30), ("recent", 1): (25, 25), ("recent", 2): (25, 24)}
        w, labels = window_with((20, 20), arms)
        plan = {"windows": {2: w}}
        self.assertEqual(sd.pairs_of(plan, labels, 20, "stale-all"), [(1, 2)])
        self.assertEqual(sd.pairs_of(plan, labels, 20, "stale-same-size"), [])

    def test_two_training_sets_of_the_same_genomes_make_no_pair(self):
        # Shares 25 and 50 give one cut point when the first 6 dates are equal: all(1) is all(2).
        plan, _ = hand_plan(shares=(25, 50), changes={g: (1, s, l, k) for g, _, s, l, k in HAND[:6]})
        self.assertEqual(plan["cuts"], [1, 1])
        w = plan["windows"][2]
        self.assertEqual(w["arms"][("all", 1)], w["arms"][("all", 2)])
        # The same at the sizes of a pair: 25 of each class in both, the same genomes.
        w, labels = window_with((20, 20), self.arms())
        w["arms"][("all", 2)] = list(w["arms"][("all", 1)])
        self.assertEqual(sd.pairs_of({"windows": {2: w}}, labels, 20, "stale-all"), [])


def count_world(base):
    """One made-up species of 200 genomes, each public on its own day from 2020-01-01, each of a study of its own,
    none with a collection year, R and S in turn. drug_a labels all 200, drug_b the first 160."""
    start = datetime.date(2020, 1, 1)
    genomes = [{"id": f"syn-n-{i:03d}", "slug": "count_species",
                "day": (start + datetime.timedelta(days=i)).isoformat() + "T00:00:00Z", "year": None,
                "study": f"syn-study-n-{i:03d}", "digest": f"syn-digest-n-{i}", "in_matrix": True} for i in range(200)]
    datasets = [{"id": f"count_species__{drug}", "species_slug": "count_species", "antibiotic": drug,
                 "labels": {g["id"]: "RS"[i % 2] for i, g in enumerate(genomes[:n])}}
                for drug, n in (("drug_a", 200), ("drug_b", 160))]
    return fx.write_custom(base, genomes, datasets)


class Counts(unittest.TestCase):
    """count on the made-up species of count_world, shares 40, 60 and 80 on the public day.

    drug_a, 200 genomes: cut points at the 80th, 120th and 160th; windows 2 and 3 each test 20 R and 20 S, all of new
    studies, and all(1) holds 40 of each class. So stale all has the pairs (1, 2) and (2, 3) at horizon 1 and (1, 3)
    at horizon 2, at every test minimum. drug_b, 160 genomes: windows of 16 R and 16 S, so it has the same pairs at
    test minimums 15 and 10, and none at 20."""

    @classmethod
    def setUpClass(cls):
        cls.tmp = Path(tempfile.mkdtemp(prefix="stagger-count-"))
        rows = sd.count(sd.load_inputs(*count_world(cls.tmp / "w")))[0]
        cls.rows = {(r["clock"], r["shares"], r["test_minimum"], r["contrast"], r["horizon"]): r for r in rows}

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp)

    def row(self, minimum, horizon, contrast="stale-all"):
        r = self.rows[("public-day", "40-60-80", minimum, contrast, horizon)]
        return r["datasets"], r["species"], r["pairs"]

    def test_a_dataset_with_two_pairs_at_a_horizon_counts_once(self):
        # Test minimum 20: drug_a alone, 2 pairs at horizon 1 and 1 at horizon 2, on both contrasts.
        self.assertEqual(self.row(20, 1), (1, 1, 2))
        self.assertEqual(self.row(20, 2), (1, 1, 1))
        self.assertEqual(self.row(20, 1, "stale-same-size"), (1, 1, 2))

    def test_two_datasets_of_one_species_count_one_species(self):
        # Test minimum 15: drug_a and drug_b, one species, 2 pairs each at horizon 1 and 1 each at horizon 2.
        self.assertEqual(self.row(15, 1), (2, 1, 4))
        self.assertEqual(self.row(10, 2), (2, 1, 2))


def entry(shares, minimum, datasets, species):
    return {"shares": shares, "minimum": minimum, "datasets": datasets, "species": species}


def taken(entries):
    t, descriptive, steps = sd.choose(entries)
    return sd.variant_name(t["shares"], t["minimum"]), descriptive, [(s, len(k)) for s, _, k in steps]


class VariantRule(unittest.TestCase):
    def test_variants_in_the_plans_order(self):
        names = [sd.variant_name(s, m) for s, m in sd.variants()]
        self.assertEqual(names[:4], ["40-60-80/20", "40-60-80/15", "40-60-80/10", "50-70/20"])
        self.assertEqual(len(names), 12)

    def test_step_1_needs_12_datasets_and_3_species(self):
        got = taken([entry((40, 60, 80), 20, 11, 5), entry((40, 60, 80), 15, 30, 2), entry((50, 70), 20, 12, 3)])
        self.assertEqual(got[0], "50-70/20")
        self.assertFalse(got[1])

    def test_step_2_the_most_origins(self):
        got = taken([entry((50, 70), 20, 40, 6), entry((40, 60, 80), 10, 12, 3)])
        self.assertEqual(got[0], "40-60-80/10")

    def test_step_3_the_highest_test_minimum(self):
        got = taken([entry((40, 60, 80), 10, 40, 6), entry((40, 60, 80), 15, 12, 3), entry((50, 70), 20, 50, 6)])
        self.assertEqual(got[0], "40-60-80/15")

    def test_step_4_the_most_datasets(self):
        # 50-75/15 holds the most datasets and loses at step 3; step 4 then takes 40-60/20 over 50-70/20.
        got = taken([entry((40, 60, 80), 20, 11, 3), entry((50, 70), 20, 13, 3), entry((40, 60), 20, 14, 3),
                     entry((50, 75), 15, 30, 3)])
        self.assertEqual(got[0], "40-60/20")
        self.assertEqual(got[2][:4], [("1", 3), ("2", 3), ("3", 2), ("4", 1)])

    def test_a_tie_goes_to_the_first_in_the_plans_order(self):
        got = taken([entry((50, 70), 20, 14, 3), entry((40, 60), 20, 14, 3), entry((50, 75), 20, 14, 3)])
        self.assertEqual(got[0], "50-70/20")
        self.assertEqual(got[2], [("1", 3), ("2", 3), ("3", 3), ("4", 3), ("5", 1)])

    def test_none_reaches_12_in_3_the_most_datasets(self):
        got = taken([entry((40, 60, 80), 20, 11, 4), entry((40, 60, 80), 15, 30, 2), entry((50, 70), 20, 30, 1)])
        self.assertEqual(got[0], "40-60-80/15")
        self.assertEqual(got[2], [("1", 0), ("2", 2), ("3", 1)])

    def test_none_reaches_12_in_3_the_clock_is_descriptive_only(self):
        self.assertTrue(taken([entry((40, 60, 80), 20, 11, 4), entry((50, 70), 20, 30, 2)])[1])
        self.assertFalse(taken([entry((40, 60, 80), 20, 12, 3), entry((50, 70), 20, 30, 2)])[1])

    def test_the_rule_reads_stale_all_at_horizon_1(self):
        rows = []
        for shares, minimum in sd.variants():
            for contrast in sd.CONTRASTS:
                for h in range(1, len(shares)):
                    # Stale all holds 12 in 3 at horizon 1 for 50-70/10 alone; same size for 40-60-80/20 alone.
                    big = (contrast == "stale-all" and h == 1 and (shares, minimum) == ((50, 70), 10)) or \
                          (contrast == "stale-same-size" and h == 1 and (shares, minimum) == ((40, 60, 80), 20)) or \
                          (contrast == "stale-all" and h == 2 and (shares, minimum) == ((40, 60, 80), 20))
                    rows.append({"clock": "public-day", "shares": "-".join(map(str, shares)), "test_minimum": minimum,
                                 "contrast": contrast, "horizon": h, "datasets": 20 if big else 3,
                                 "species": 4 if big else 1})
        got = taken(sd.choice_entries(rows, "public-day"))
        self.assertEqual(got[:2], ("50-70/10", False))


class Clocks(unittest.TestCase):
    def test_the_collection_year_is_four_digits(self):
        cases = [("2015", 2015), (2015, 2015), (" 2015 ", 2015), ("2015-2016", None), ("15", None), (2015.0, None),
                 ("", None), (None, None), ("20155", None)]
        for value, want in cases:
            self.assertEqual(sd.year_of(value), want, value)

    def test_each_clock_its_own_dates_and_hash(self):
        with tempfile.TemporaryDirectory(prefix="stagger-") as tmp:
            genomes = [{"id": f"syn-c-{i:03d}", "slug": "clock_species", "day": f"2021-02-{i + 1:02d}T00:00:00Z",
                        "year": None if i == 3 else 2000 + (9 - i), "study": f"syn-study-{i}",
                        "digest": f"syn-digest-{i}", "in_matrix": True} for i in range(10)]
            labels = {g["id"]: "RS"[i % 2] for i, g in enumerate(genomes)}
            root, fetch, built = fx.write_custom(tmp, genomes, [
                {"id": "clock_species__drug", "species_slug": "clock_species", "antibiotic": "drug",
                 "labels": labels}])
            inputs = sd.load_inputs(root, fetch, built)
            day, year = inputs["dates"]["public-day"], inputs["dates"]["collection-year"]
            self.assertEqual(day["syn-c-000"], "2021-02-01")
            self.assertEqual((year["syn-c-000"], year["syn-c-003"]), (2009, None))
            spec = inputs["specs"][0]
            cohort, lab, sp = sd.cohort_of(inputs, spec)
            by_day = sd.plan_dataset(spec["id"], "public-day", (50, 70), cohort, lab, day, inputs["study"],
                                     sp["listed"], sp["digests"])
            by_year = sd.plan_dataset(spec["id"], "collection-year", (50, 70), cohort, lab, year, inputs["study"],
                                      sp["listed"], sp["digests"])
            # Days run with the ids, years against them; syn-c-003 has no year.
            self.assertEqual(by_day["cuts"], ["2021-02-05", "2021-02-07"])
            self.assertEqual(by_day["windows"][2]["test"], ["syn-c-007", "syn-c-008", "syn-c-009"])
            # 9 years: places ceil(4.5) = 5 and ceil(6.3) = 7 of 2000, 2001, 2002, 2003, 2004, 2005, 2007, 2008, 2009.
            self.assertEqual(by_year["cuts"], [2004, 2007])
            self.assertEqual(by_year["windows"][2]["test"], ["syn-c-000", "syn-c-001"])
            self.assertEqual((by_year["dated"], by_year["undated"]), (9, 1))
            for members in by_year["windows"][2]["arms"].values():
                self.assertNotIn("syn-c-003", members)


CHECK_HEAD = "dataset\tspecies\tantibiotic\tcohort\tcohort_resistant\tcut_day\tshare_public_by_cut\n"


def check_world(base, undated=()):
    """24 made-up datasets of one species in the 14 Nov record's shape. Genome i is public on 2020-01-(i + 1);
    syn-k-009 is out of the matrix; dataset k labels genomes 0 to 7 + k. Genomes named in undated have no day.
    Returns (paths, {dataset: the cut day at 60%, worked out here}, {dataset: cohort})."""
    genomes = [{"id": f"syn-k-{i:03d}", "slug": "check_species",
                "day": None if i in undated else f"2020-01-{i + 1:02d}T10:00:00Z", "year": 2020, "study": None,
                "digest": f"syn-digest-{i}", "in_matrix": i != 9} for i in range(31)]
    datasets, cuts, cohorts = [], {}, {}
    for k in range(24):
        name = f"check_species__drug_{k:02d}"
        labelled = list(range(8 + k))
        datasets.append({"id": name, "species_slug": "check_species", "antibiotic": f"drug_{k:02d}",
                         "labels": {genomes[i]["id"]: "RS"[i % 2] for i in labelled}})
        cohort = [i for i in labelled if i != 9]
        days = sorted(i + 1 for i in cohort if i not in undated)
        # The place ceil(0.6 x n) of the n dated genomes, in whole numbers.
        cuts[name] = f"2020-01-{days[-(-6 * len(days) // 10) - 1]:02d}"
        cohorts[name] = len(cohort)
    return fx.write_custom(base, genomes, datasets), cuts, cohorts


class Check60(unittest.TestCase):
    """A made-up dataset table of 24 rows in the 14 Nov record's columns, right, and broken one way at a time."""

    @classmethod
    def setUpClass(cls):
        cls.tmp = Path(tempfile.mkdtemp(prefix="stagger-check-"))
        cls.paths, cls.cuts, cls.cohorts = check_world(cls.tmp / "w")
        cls.names = sorted(cls.cuts)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp)

    def table(self, name, names=None, cuts=None):
        """The table of the given datasets in order, each with its cut day here or as given in cuts."""
        path = self.tmp / name
        rows = "".join(f"{n}\tcheck_species\t{n.split('__')[1]}\t{self.cohorts.get(n, 9)}\t1\t"
                       f"{(cuts or {}).get(n, self.cuts.get(n, '2020-01-05'))}\t0.6\n"
                       for n in (self.names if names is None else names))
        path.write_text(CHECK_HEAD + rows)
        return path

    def run_check(self, table, out, paths=None):
        return fx.run(sd, "check60", *(paths or self.paths), table, self.tmp / out)

    def stopped(self, table, message):
        code, text, err = self.run_check(table, table.stem + "-out")
        self.assertEqual((code, text), (1, ""))
        self.assertEqual(err, f"stagger_design: stopped: {message}\n")
        # The stop names no cut day, and nothing is written.
        self.assertNotIn("2020-", err)
        self.assertFalse((self.tmp / (table.stem + "-out")).exists())

    def test_the_cut_days_of_the_table(self):
        # check_species__drug_00: genomes 0 to 7, 8 dated, ceil(4.8) = 5: the 5th day, 2020-01-05.
        self.assertEqual(self.cuts["check_species__drug_00"], "2020-01-05")
        # check_species__drug_23: genomes 0 to 30 less syn-k-009, 30 dated, ceil(18) = 18: day 19, 2020-01-19.
        self.assertEqual(self.cuts["check_species__drug_23"], "2020-01-19")
        code, text, _ = self.run_check(self.table("right.tsv"), "right")
        self.assertEqual(code, 0)
        self.assertIn("check60: 24 of 24 cut days are the record's", text)
        rows = (self.tmp / "right" / "check60.tsv").read_text().splitlines()
        self.assertEqual(len(rows), 25)
        self.assertEqual(rows[24].split("\t")[:6], ["check_species__drug_23", "30", "30", "0", "2020-01-19",
                                                    "2020-01-19"])

    def test_a_cut_day_other_than_the_tables(self):
        code, text, _ = self.run_check(self.table("wrong.tsv", cuts={"check_species__drug_05": "2020-01-06"}),
                                       "wrong")
        self.assertEqual(code, 1)
        self.assertIn("check60: 23 of 24 cut days are the record's", text)

    def test_a_table_of_23_datasets_stops(self):
        need = "the 60% check needs 24 dataset rows, each a different dataset of the metadata; the reference table has"
        self.stopped(self.table("short.tsv", self.names[:11] + self.names[12:]), need + " 23")
        self.stopped(self.table("one.tsv", self.names[:1]), need + " 1")
        self.stopped(self.table("more.tsv", self.names + self.names[:1]), need + " 25")

    def test_the_count_is_24_whatever_the_table(self):
        self.assertEqual(sd.CHECK_DATASETS, 24)
        source = fx.DESIGN.read_text(encoding="utf-8")
        self.assertEqual(source.count("CHECK_DATASETS = "), 1)
        self.assertNotIn("environ", source)
        # A right table of 2 rows passes when a test passes its own count, and stops the script, which takes 24.
        two = self.table("two.tsv", self.names[:2])
        rows, ok = sd.check60(sd.load_inputs(*self.paths), two, 2)
        self.assertEqual((len(rows), ok), (2, True))
        self.stopped(two, "the 60% check needs 24 dataset rows, each a different dataset of the metadata; the "
                          "reference table has 2")

    def test_one_dataset_twice_and_one_left_out_stops(self):
        self.stopped(self.table("twice.tsv", self.names[:23] + self.names[:1]),
                     "the 60% check needs 24 different datasets; the reference table's 24 rows name 23, and these "
                     "more than once: check_species__drug_00")

    def test_a_dataset_the_metadata_lacks_stops(self):
        self.stopped(self.table("lacks.tsv", self.names[:23] + ["check_species__drug_99"]),
                     "the 60% check needs 24 datasets of the metadata; the metadata lacks 1 of the reference table's "
                     "datasets: check_species__drug_99")

    def test_a_dataset_name_in_upper_case_stops(self):
        # Names are matched letter for letter: the metadata has no dataset of that name.
        names = self.names[:5] + [self.names[5].upper()] + self.names[6:]
        self.stopped(self.table("upper.tsv", names),
                     "the 60% check needs 24 datasets of the metadata; the metadata lacks 1 of the reference table's "
                     "datasets: CHECK_SPECIES__DRUG_05")

    def test_a_dataset_name_with_a_space_after_it_stops(self):
        names = self.names[:5] + [self.names[5] + " "] + self.names[6:]
        self.stopped(self.table("spaced.tsv", names),
                     "the 60% check needs 24 datasets of the metadata; the metadata lacks 1 of the reference table's "
                     "datasets: check_species__drug_05 ")

    def test_a_row_of_one_cell_stops_with_one_line(self):
        # 24 rows, the 6th its dataset's name alone: one stop line where a traceback was (R41).
        path = self.table("cell.tsv")
        lines = path.read_text().split("\n")
        lines[6] = self.names[5]
        path.write_text("\n".join(lines))
        self.stopped(path, "1 of the reference table's 24 rows hold another count of cells than its header's 7: row 6")

    def test_a_row_one_cell_short_stops(self):
        # A table with a last column the check does not read, as the 14 Nov record's has, and the 6th row without
        # that cell. Every cut day is still the record's, so only the count of cells shows it.
        path = self.table("short-row.tsv")
        lines = path.read_text().split("\n")
        lines[:25] = [lines[0] + "\ttest"] + [line + "\t3" for line in lines[1:25]]
        lines[6] = lines[6].rsplit("\t", 1)[0]
        path.write_text("\n".join(lines))
        self.stopped(path, "1 of the reference table's 24 rows hold another count of cells than its header's 8: row 6")

    def test_a_row_with_more_cells_than_the_header_stops(self):
        # Rows 6 and 7 on one line, row 7's cut day two days early. Read by the header's seven names, that line is
        # row 6 alone, and row 7's wrong day is never checked (R45). The stop names the row before the count of rows.
        self.assertEqual(self.cuts[self.names[6]], "2020-01-08")
        path = self.table("long-row.tsv", cuts={self.names[6]: "2020-01-06"})
        lines = path.read_text().split("\n")
        lines[6:8] = [lines[6] + "\t" + lines[7]]
        path.write_text("\n".join(lines))
        self.stopped(path, "1 of the reference table's 23 rows hold another count of cells than its header's 7: row 6")
        # Handed its own count of 23, the check stops too, where it would read every cut day as the record's.
        with self.assertRaises(sd.Stop):
            sd.check60(sd.load_inputs(*self.paths), path, 23)

    def test_an_undated_cohort_genome_fails_its_dataset(self):
        # syn-k-030 is labelled in check_species__drug_23 alone. Undated, it leaves 29 dated genomes there:
        # ceil(17.4) = 18, the 18th of days 1 to 9 and 11 to 30, day 19, the table's. The cut day agrees, and the
        # dataset still fails, as the 14 Nov rule needs every cohort genome dated.
        paths, cuts, _ = check_world(self.tmp / "w-undated", undated=(30,))
        self.assertEqual(cuts["check_species__drug_23"], "2020-01-19")
        code, text, _ = self.run_check(self.table("undated.tsv"), "undated", paths)
        self.assertEqual(code, 1)
        self.assertIn("check60: 23 of 24 cut days are the record's", text)
        rows = [r.split("\t") for r in (self.tmp / "undated" / "check60.tsv").read_text().splitlines()]
        self.assertEqual(rows[24][:6] + rows[24][-1:], ["check_species__drug_23", "30", "30", "1", "2020-01-19",
                                                        "2020-01-19", "0"])
        self.assertEqual({r[-1] for r in rows[1:24]}, {"1"})


class BarredPaths(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = Path(tempfile.mkdtemp(prefix="stagger-paths-"))
        cls.root, cls.fetch, cls.built, _ = fx.write_world(cls.tmp / "w", n=60, species=fx.SPECIES[:1],
                                                         drugs=fx.DRUGS[:1])

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp)

    def link(self, name, target):
        path = self.tmp / name
        if not path.exists():
            os.symlink(target, path)
        return path

    def stopped(self, argv, role, word):
        code, out, err = fx.run(sd, *argv)
        self.assertEqual(code, 1)
        self.assertEqual(out, "")
        self.assertEqual(err, f"stagger_design: stopped: the {role} path holds '{word}'; nothing under it is opened\n")
        self.assertNotIn(str(self.tmp), err)

    def test_each_word_and_role(self):
        out = self.tmp / "out"
        self.stopped(["count", self.link("latetest-root", self.root), self.fetch, self.built, out], "data root",
                     "latetest")
        self.stopped(["count", self.root, self.link("held-out-fetch", self.fetch), self.built, out], "fetch dir",
                     "held-out")
        self.stopped(["count", self.root, self.fetch, self.link("sealed-built", self.built), out], "built root",
                     "sealed")
        self.stopped(["folds", self.root, self.fetch, self.built, self.tmp / "Sealed-out"], "out dir", "sealed")
        self.stopped(["check60", self.root, self.fetch, self.built, self.tmp / "latetest.tsv", out],
                     "reference table", "latetest")
        self.assertFalse(out.exists())

    def test_a_link_that_leads_into_a_sealed_place(self):
        target = self.tmp / "sealed-set"
        if not target.exists():
            shutil.copytree(self.fetch, target)
        self.stopped(["count", self.root, self.link("plain-name", target), self.built, self.tmp / "o2"], "fetch dir",
                     "sealed")


class Folds(unittest.TestCase):
    """count and folds on the made-up world of 3 species and 5 drugs."""

    @classmethod
    def setUpClass(cls):
        cls.tmp = Path(tempfile.mkdtemp(prefix="stagger-folds-"))
        cls.root, cls.fetch, cls.built, cls.labels = fx.write_world(cls.tmp / "w")
        cls.out = cls.tmp / "design"
        cls.code, cls.text, cls.err = fx.run(sd, "folds", cls.root, cls.fetch, cls.built, cls.out)
        cls.design = json.loads((cls.out / "design.json").read_text())
        cls.cohorts = {}
        inputs = sd.load_inputs(cls.root, cls.fetch, cls.built)
        for spec in inputs["specs"]:
            cls.cohorts[spec["id"]] = sd.cohort_of(inputs, spec)[0]

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp)

    def fits(self):
        for d in self.design["datasets"]:
            for w in d["windows"]:
                for f in w["fits"]:
                    yield d, w, f

    def test_it_ran(self):
        self.assertEqual(self.code, 0, self.err)
        self.assertIn("count: 60 rows; it read no score and no result", self.text)

    def test_a_count_row_per_clock_share_set_minimum_contrast_and_horizon(self):
        rows = (self.out / "counts.tsv").read_text().splitlines()
        self.assertEqual(rows[0].split("\t"), sd.COUNT_COLUMNS)
        # Per clock: 3 minimums x 2 contrasts x (2 + 1 + 1 + 1) horizons.
        self.assertEqual(len(rows) - 1, 2 * 3 * 2 * 5)

    def test_the_choice_names_the_variant_of_each_clock(self):
        rows = [r.split("\t") for r in (self.out / "choice.tsv").read_text().splitlines()[1:]]
        for clock in sd.CLOCKS:
            last = [r for r in rows if r[0] == clock][-1]
            got = self.design["clocks"][clock]
            self.assertEqual(last[4], sd.variant_name(got["shares"], got["test_minimum"]))
            self.assertEqual(last[2], "a tie goes to the first in the plan's order")

    def test_every_cohort_genome_once_in_each_fold_file(self):
        seen = set()
        for d, w, f in self.fits():
            if f["file"] in seen:
                continue
            seen.add(f["file"])
            body = json.loads((self.out / "fold-files" / f["file"]).read_text())
            self.assertEqual(set(body), {"format", "dataset", "scheme", "test", "train", "unused"})
            self.assertEqual((body["format"], body["dataset"], body["scheme"]), ("grm-fold-file-v1", d["id"],
                                                                                  f["scheme"]))
            listed = body["test"][0] + body["train"] + body["unused"]
            self.assertEqual(sorted(listed), sorted(self.cohorts[d["id"]]), f["file"])
            self.assertEqual((len(body["test"][0]), len(body["train"])), (w["test"], f["train"]))
        self.assertGreater(len(seen), 100)

    def test_each_file_named_with_the_sha256_of_its_bytes(self):
        on_disk = sorted(p.name for p in (self.out / "fold-files").iterdir())
        named = sorted({f["file"] for _, _, f in self.fits()})
        self.assertEqual(on_disk, named)
        for _, _, f in self.fits():
            data = (self.out / "fold-files" / f["file"]).read_bytes()
            self.assertEqual(f["sha256"], hashlib.sha256(data).hexdigest(), f["file"])

    def test_both_fits_of_a_pair_test_the_same_genomes(self):
        for d in self.design["datasets"]:
            for w in d["windows"]:
                tests = {json.dumps(json.loads((self.out / "fold-files" / f["file"]).read_text())["test"])
                         for f in w["fits"]}
                self.assertEqual(len(tests), 1, (d["id"], w["window"]))

    def test_recent_1_and_all_1_share_one_file(self):
        shared = 0
        for d, w, f in self.fits():
            if f["origin"] == 1 and f["arm"] == "recent":
                shared += 1
                self.assertEqual(f["same_file_as"], "all")
                self.assertTrue(f["scheme"].endswith("-o1-all"))
        self.assertGreater(shared, 0)
        names = [p.name for p in (self.out / "fold-files").iterdir()]
        self.assertFalse([n for n in names if "-o1-recent" in n])
        self.assertIn("files serve both all(1) and recent(1)", self.text)

    def test_the_design_holds_no_genome_id(self):
        text = (self.out / "design.json").read_text() + (self.out / "counts.tsv").read_text()
        self.assertNotIn("syn-a-", text)
        self.assertNotIn("syn-study-", text)

    def test_a_species_not_built_gets_no_fold_file(self):
        out = self.tmp / "unbuilt"
        root, fetch, built, _ = fx.write_world(self.tmp / "w-unbuilt", built=False)
        code, _, err = fx.run(sd, "folds", root, fetch, built, out)
        self.assertEqual(code, 0, err)
        design = json.loads((out / "design.json").read_text())
        self.assertFalse(list((out / "fold-files").iterdir()) if (out / "fold-files").exists() else [])
        self.assertTrue(all(d["why_not"] for d in design["datasets"]))

    def test_an_out_dir_that_holds_a_file_is_refused(self):
        code, _, err = fx.run(sd, "count", self.root, self.fetch, self.built, self.out)
        self.assertEqual(code, 1)
        self.assertIn("the out dir holds a file", err)


if __name__ == "__main__":
    unittest.main()
