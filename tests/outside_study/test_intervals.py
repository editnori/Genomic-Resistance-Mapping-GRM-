"""Clopper-Pearson intervals: the values src/intervals.rs tests against
(scipy's beta.ppf) and the closed forms at the boundaries."""

import unittest

from tests.outside_study.fixtures import S

ALPHA = 0.05
# src/intervals.rs, `SCIPY`: (x, n, lower, upper) from scipy.stats.beta.ppf.
SCIPY = [
    (0, 6, 0.0, 0.45925812643990044),
    (0, 7, 0.0, 0.40961639722500337),
    (0, 21, 0.0, 0.1610976152190796),
    (0, 1886, 0.0, 0.00195401602039081),
    (3, 111, 0.005608683216675948, 0.0769555750675793),
    (4, 87, 0.012667375820641409, 0.11355344796666053),
    (5, 65, 0.025448467019043023, 0.17045629048892497),
    (5, 83, 0.019847321241740542, 0.1350082396842335),
    (5, 209, 0.0078124846568847685, 0.05494333888926265),
    (6, 110, 0.020277720969548958, 0.11494662512228843),
    (6, 280, 0.007903552986266925, 0.04605508562773042),
    (7, 1319, 0.0021362919591686753, 0.010903806367681847),
    (8, 101, 0.03481691039879477, 0.15011620100508763),
    (8, 317, 0.010956946291487912, 0.04911870180113821),
    (8, 1886, 0.0018330248422895048, 0.008340805876968928),
    (10, 10, 0.6915028921812392, 1.0),
    (11, 11, 0.7150858470818455, 1.0),
    (11, 120, 0.04665226965524853, 0.15809854119583328),
    (14, 57, 0.14126857773596005, 0.37761030480194585),
    (22, 175, 0.08048892787899473, 0.18412749654984673),
    (23, 451, 0.032599285708675686, 0.07554156137766942),
    (24, 1886, 0.008169922719802096, 0.01887546766118562),
    (26, 80, 0.22445958517745335, 0.4389370456044484),
    (41, 304, 0.09855207074221606, 0.1784867586358306),
    (53, 408, 0.0988407923817262, 0.16645003066433706),
    (63, 848, 0.05755743818575499, 0.09405440618587402),
    (75, 619, 0.09650626678664745, 0.14949904178259882),
    (78, 618, 0.10105919143253055, 0.15500146718542052),
    (83, 161, 0.43556769905190557, 0.5949061282319127),
    (208, 530, 0.3506387002748109, 0.43547319573686855),
]


class ClopperPearson(unittest.TestCase):
    def test_matches_the_scipy_values_src_intervals_rs_uses(self):
        for x, n, lo, hi in SCIPY:
            got_lo, got_hi = S.clopper_pearson(x, n)
            self.assertAlmostEqual(got_lo, lo, delta=1e-12, msg=f"{x}/{n} lower")
            self.assertAlmostEqual(got_hi, hi, delta=1e-12, msg=f"{x}/{n} upper")

    def test_closed_forms_at_the_boundaries(self):
        for n in (1, 2, 6, 25, 400, 1886, 10000):
            lo, hi = S.clopper_pearson(0, n)
            self.assertEqual(lo, 0.0)
            self.assertAlmostEqual(hi, 1.0 - (ALPHA / 2.0) ** (1.0 / n), delta=1e-12)
            lo, hi = S.clopper_pearson(n, n)
            self.assertEqual(hi, 1.0)
            self.assertAlmostEqual(lo, (ALPHA / 2.0) ** (1.0 / n), delta=1e-12)
        lo, hi = S.clopper_pearson(1, 2)
        root = (1.0 - ALPHA / 2.0) ** 0.5
        self.assertAlmostEqual(lo, 1.0 - root, delta=1e-12)
        self.assertAlmostEqual(hi, root, delta=1e-12)

    def test_no_interval_without_genomes(self):
        self.assertIsNone(S.clopper_pearson(0, 0))
        self.assertIsNone(S.clopper_pearson(3, 2))

    def test_tails_hold_alpha_over_two(self):
        for x, n in ((1, 9), (5, 83), (40, 300), (48, 60)):
            lo, hi = S.clopper_pearson(x, n)
            self.assertAlmostEqual(S.beta_cdf(lo, x, n - x + 1), ALPHA / 2, delta=1e-12)
            self.assertAlmostEqual(1 - S.beta_cdf(hi, x + 1, n - x), ALPHA / 2, delta=1e-12)
            self.assertLess(lo, x / n)
            self.assertLess(x / n, hi)


if __name__ == "__main__":
    unittest.main()
