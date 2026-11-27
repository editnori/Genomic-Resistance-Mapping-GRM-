"""The measure of scripts/matched-accuracy/stage2.py on made-up numbers: the row mean, the pair count, the weights,
the parts per million, the interval, the word beside it and the 10-pair rule. Every number here is made up for the
test."""

import fractions
import math
import random
import unittest

import numpy

from tests.matched_accuracy.stage2_common import stage2

F = fractions.Fraction
# The brief's five values: added left to right they give 2.4999999999999996, where Python 3.12's compensated sum
# and math.fsum give 2.5.
FIVE_VALUES = (0, 0.12, 0.99, 0.99, 0.4)


def by_hand_left_to_right(values):
    total = values[0]
    total = total + values[1]
    total = total + values[2]
    total = total + values[3]
    total = total + values[4]
    return total / 5


def compared(groups):
    """{key: (species, dataset, half points, pairs, rows)} from [(species, dataset, half points, pairs)]."""
    return {(sp, d, i): (sp, d, h, p, 2) for i, (sp, d, h, p) in enumerate(groups)}


def interval_by_hand(values):
    """The interval as the protocol writes it, computed here apart from the script."""
    means = numpy.array(values, dtype=numpy.float64)
    rng = numpy.random.Generator(numpy.random.PCG64(42))
    idx = rng.integers(0, len(means), size=(10000, len(means)))
    lo, hi = numpy.percentile(means[idx].mean(axis=1), [2.5, 97.5])
    return float(lo), float(hi)


class RowMean(unittest.TestCase):
    def test_the_five_scores_are_added_left_to_right_then_divided_by_5(self):
        got = stage2.row_mean(FIVE_VALUES)
        self.assertEqual(got, by_hand_left_to_right(FIVE_VALUES))
        self.assertEqual(repr(got), "0" + "." + "4" + "9" * 15)
        self.assertNotEqual(got, sum(FIVE_VALUES) / 5)

    def test_a_row_of_the_five_values_loses_to_a_row_of_five_halves(self):
        # Left to right the row's mean is just under one half, so it loses; a compensated sum would tie it.
        half = stage2.row_mean((0.5,) * 5)
        self.assertEqual(stage2.group_count([stage2.row_mean(FIVE_VALUES)], [half]), (0, 1))

    def test_the_order_is_scm_cart_rf_gbt_lasso(self):
        self.assertNotEqual(stage2.row_mean(FIVE_VALUES), stage2.row_mean(FIVE_VALUES[::-1]))
        self.assertEqual(stage2.row_mean(FIVE_VALUES[::-1]), 0.5)

    def test_a_component_null_true_text_or_not_finite_stops(self):
        for bad in (None, True, False, "0", float("nan"), float("inf"), 10 ** 400):
            with self.subTest(value=repr(bad)[:12]):
                with self.assertRaises(stage2.Stop):
                    stage2.row_mean((0.5, 0.5, bad, 0.5, 0.5))


class Pairs(unittest.TestCase):
    def test_a_higher_resistant_mean_counts_1_a_tie_one_half_and_lower_0(self):
        self.assertEqual(stage2.group_count([0.6], [0.5]), (2, 1))
        self.assertEqual(stage2.group_count([0.5], [0.5]), (1, 1))
        self.assertEqual(stage2.group_count([0.4], [0.5]), (0, 1))

    def test_every_resistant_row_meets_every_susceptible_row(self):
        # Two ties, two wins: 6 half points over 4 pairs, 3/4.
        self.assertEqual(stage2.group_count([0.6, 0.6], [0.6, 0.3]), (6, 4))
        self.assertEqual(stage2.group_count([0.9, 0.7], [0.8, 0.2, 0.1]), (10, 6))


class Weights(unittest.TestCase):
    def test_pairs_are_not_pooled_across_groups(self):
        # One dataset: a group of value 1 over 1 pair, and one of value 0 over 3 pairs. Equal groups give 1/2.
        lv = stage2.levels(compared([("s", "d", 2, 1), ("s", "d", 0, 3)]))
        self.assertEqual(lv["dataset"]["d"]["value"], F(1, 2))

    def test_groups_are_not_pooled_across_datasets(self):
        # One species: a dataset of one group of value 1, and one of three groups of value 0. Equal datasets: 1/2.
        lv = stage2.levels(compared([("s", "d1", 2, 1), ("s", "d2", 0, 1), ("s", "d2", 0, 1), ("s", "d2", 0, 1)]))
        self.assertEqual(lv["species"]["s"]["value"], F(1, 2))

    def test_datasets_are_not_pooled_across_species(self):
        # A species of one dataset of value 1, and one of three datasets of value 0. Equal species: 1/2.
        lv = stage2.levels(compared([("s1", "a", 2, 1), ("s2", "b", 0, 1), ("s2", "c", 0, 1), ("s2", "e", 0, 1)]))
        self.assertEqual(lv["overall"]["all"]["value"], F(1, 2))

    def test_values_are_exact_fractions(self):
        lv = stage2.levels(compared([("s1", "a", 1, 3), ("s2", "b", 2, 7), ("s2", "c", 5, 9)]))
        for level in stage2.LEVELS:
            for s in lv[level].values():
                self.assertIsInstance(s["value"], fractions.Fraction)
        self.assertEqual(lv["overall"]["all"]["value"], (F(1, 6) + (F(1, 7) + F(5, 18)) / 2) / 2)

    def test_the_counts_beneath_each_value(self):
        lv = stage2.levels(compared([("s1", "a", 1, 3), ("s2", "b", 2, 7), ("s2", "c", 5, 9)]))
        self.assertEqual({k: lv["overall"]["all"][k] for k in ("rows", "groups", "pairs", "datasets", "species")},
                         {"rows": 6, "groups": 3, "pairs": 19, "datasets": 3, "species": 2})
        self.assertEqual(lv["species"]["s2"]["pairs"], 16)


class Ppm(unittest.TestCase):
    def test_rounded_half_to_even_from_the_exact_fraction(self):
        self.assertEqual(stage2.ppm(F(1, 2000000)), 0)
        self.assertEqual(stage2.ppm(F(3, 2000000)), 2)
        self.assertEqual(stage2.ppm(F(263, 384)), 684896)
        self.assertEqual(stage2.ppm(F(1, 3)), 333333)
        self.assertEqual(stage2.ppm(F(2, 3)), 666667)
        self.assertEqual(stage2.ppm(1), 1000000)

    def test_a_bound_is_rounded_from_the_percentiles_double(self):
        # The double nearest 2.5 parts per million lies just above it, so its exact value rounds up.
        x = 0.0000025
        self.assertGreater(F(x), F(25, 10 ** 7))
        self.assertEqual(stage2.ppm(x), 3)
        self.assertEqual(stage2.ppm(0.5), 500000)
        self.assertEqual(stage2.ppm(0.375), 375000)


class Interval(unittest.TestCase):
    def values(self, n, seed=7):
        rng = random.Random(seed)
        return [F(rng.randint(1, 10 ** 6), 10 ** 6) for _ in range(n)]

    def test_the_bounds_are_those_of_pcg64_seed_42_and_the_linear_percentiles(self):
        for n in (6, 9, 16):
            with self.subTest(species=n):
                v = self.values(n)
                self.assertEqual(stage2.interval(v), interval_by_hand([float(x) for x in v]))

    def test_another_percentile_method_gives_other_bounds(self):
        # The check above can tell the methods apart on these values.
        v = [float(x) for x in self.values(9)]
        rng = numpy.random.Generator(numpy.random.PCG64(42))
        draws = numpy.array(v)[rng.integers(0, 9, size=(10000, 9))].mean(axis=1)
        for method in ("lower", "higher", "nearest", "midpoint"):
            self.assertNotEqual(tuple(numpy.percentile(draws, [2.5, 97.5], method=method)), interval_by_hand(v),
                                method)

    def test_each_estimate_starts_a_new_generator_at_seed_42(self):
        v = self.values(8)
        first, second = stage2.interval(v), stage2.interval(v)
        self.assertEqual(first, second)
        self.assertEqual(first, interval_by_hand([float(x) for x in v]))

    def test_species_are_drawn_in_sorted_order_of_their_names(self):
        v = self.values(7, seed=3)
        names = ["sp_g", "sp_a", "sp_e", "sp_b", "sp_f", "sp_c", "sp_d"]
        groups = [(names[i], f"d{i}", 2 * x.numerator, x.denominator) for i, x in enumerate(v)]
        s = stage2.summarise(compared(groups))
        ordered = [float(x) for _, x in sorted(zip(names, v))]
        self.assertEqual(s["interval"], interval_by_hand(ordered))
        self.assertNotEqual(interval_by_hand(ordered), interval_by_hand([float(x) for x in v]))

    def test_fewer_than_6_species_give_the_estimate_and_no_interval(self):
        v = self.values(5)
        s = stage2.summarise(compared([(f"s{i}", f"d{i}", 2 * x.numerator, x.denominator) for i, x in enumerate(v)]))
        self.assertIsNone(s["interval"])
        self.assertEqual(s["species"], 5)
        self.assertEqual(s["levels"]["overall"]["all"]["value"], sum(v, F(0)) / 5)
        v = self.values(6)
        s = stage2.summarise(compared([(f"s{i}", f"d{i}", 2 * x.numerator, x.denominator) for i, x in enumerate(v)]))
        self.assertIsNotNone(s["interval"])


class Chance(unittest.TestCase):
    """The word beside an interval, read from its bounds before they are rounded to parts per million."""

    def test_a_lower_bound_just_above_one_half_that_rounds_to_500000_gives_above(self):
        for lower in (math.nextafter(0.5, 1), 0.5000004):
            with self.subTest(lower=lower):
                self.assertGreater(lower, 0.5)
                self.assertEqual(stage2.ppm(lower), 500000)
                self.assertEqual(stage2.chance((lower, 0.7)), "above")

    def test_a_lower_bound_just_below_one_half_that_rounds_to_500000_gives_unresolved(self):
        for lower in (math.nextafter(0.5, 0), 0.4999996):
            with self.subTest(lower=lower):
                self.assertLess(lower, 0.5)
                self.assertEqual(stage2.ppm(lower), 500000)
                self.assertEqual(stage2.chance((lower, 0.7)), "unresolved")

    def test_an_upper_bound_below_one_half_gives_below_and_a_bound_at_one_half_unresolved(self):
        self.assertEqual(stage2.chance((0.3, math.nextafter(0.5, 0))), "below")
        self.assertEqual(stage2.chance((0.3, 0.4999996)), "below")
        self.assertEqual(stage2.chance((0.3, 0.5)), "unresolved")
        self.assertEqual(stage2.chance((0.5, 0.7)), "unresolved")
        self.assertEqual(stage2.chance((0.3, 0.7)), "unresolved")

    def test_an_estimate_with_no_interval_has_no_word(self):
        self.assertEqual(stage2.chance(None), "")
        row = stage2.bounds_row("A", "primary", "primary", {"interval": None, "species": 5})
        self.assertEqual((row["lower_ppm"], row["upper_ppm"], row["chance"], row["draws"]), ("", "", "", 0))

    def test_the_word_stands_beside_the_rounded_bounds(self):
        row = stage2.bounds_row("A", "primary", "primary", {"interval": (0.5000004, 0.6), "species": 6})
        self.assertEqual((row["lower_ppm"], row["upper_ppm"], row["chance"]), (500000, 600000, "above"))
        self.assertEqual(stage2.INTERVAL.index("chance"), stage2.INTERVAL.index("upper_ppm") + 1)


class TenPairs(unittest.TestCase):
    def test_a_dataset_or_species_under_10_pairs_gives_no_estimate(self):
        for level in ("dataset", "species"):
            self.assertEqual(stage2.estimate(level, {"value": F(1, 2), "pairs": 9}), "")
            self.assertEqual(stage2.estimate(level, {"value": F(1, 2), "pairs": 10}), 500000)

    def test_the_overall_figure_gives_its_estimate_and_no_pairs_none(self):
        self.assertEqual(stage2.estimate("overall", {"value": F(1, 4), "pairs": 3}), 250000)
        self.assertEqual(stage2.estimate("overall", {"value": None, "pairs": 0}), "")

    def test_a_thin_dataset_enters_its_species_mean(self):
        lv = stage2.levels(compared([("s", "thin", 2, 1), ("s", "wide", 0, 20)]))
        self.assertEqual(lv["species"]["s"]["value"], F(1, 2))
        self.assertEqual(stage2.estimate("dataset", lv["dataset"]["thin"]), "")
        self.assertEqual(stage2.estimate("species", lv["species"]["s"]), 500000)


if __name__ == "__main__":
    unittest.main()
