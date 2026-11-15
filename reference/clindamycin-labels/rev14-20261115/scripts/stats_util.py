"""Small exact statistics with the standard library only (no numpy or scipy on this machine).

clopper_pearson(k, n)      exact two-sided 95% interval for a binomial proportion
fisher_two_sided(a,b,c,d)  Fisher's exact two-sided p for the 2x2 table [[a, b], [c, d]]
"""
import math


def _log_binom_pmf(k, n, p):
    if p <= 0.0:
        return 0.0 if k == 0 else float("-inf")
    if p >= 1.0:
        return 0.0 if k == n else float("-inf")
    return (math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)
            + k * math.log(p) + (n - k) * math.log1p(-p))


def _binom_cdf(k, n, p):
    """P(X <= k)."""
    if k < 0:
        return 0.0
    if k >= n:
        return 1.0
    return min(1.0, math.fsum(math.exp(_log_binom_pmf(i, n, p)) for i in range(0, k + 1)))


def _binom_sf(k, n, p):
    """P(X >= k)."""
    if k <= 0:
        return 1.0
    if k > n:
        return 0.0
    return min(1.0, math.fsum(math.exp(_log_binom_pmf(i, n, p)) for i in range(k, n + 1)))


def clopper_pearson(k, n, conf=0.95):
    """Exact interval: lower solves P(X >= k | p) = alpha/2, upper solves P(X <= k | p) = alpha/2."""
    if n == 0:
        return (float("nan"), float("nan"))
    alpha = 1.0 - conf
    if k == 0:
        lo = 0.0
    else:
        a, b = 0.0, 1.0
        for _ in range(200):
            m = (a + b) / 2
            if _binom_sf(k, n, m) < alpha / 2:
                a = m
            else:
                b = m
        lo = (a + b) / 2
    if k == n:
        hi = 1.0
    else:
        a, b = 0.0, 1.0
        for _ in range(200):
            m = (a + b) / 2
            if _binom_cdf(k, n, m) > alpha / 2:
                a = m
            else:
                b = m
        hi = (a + b) / 2
    return (lo, hi)


def _log_hyper(a, r1, r2, c1):
    """log P(top-left = a) with row sums r1, r2 and first column sum c1."""
    return (math.lgamma(r1 + 1) - math.lgamma(a + 1) - math.lgamma(r1 - a + 1)
            + math.lgamma(r2 + 1) - math.lgamma(c1 - a + 1) - math.lgamma(r2 - c1 + a + 1)
            - (math.lgamma(r1 + r2 + 1) - math.lgamma(c1 + 1) - math.lgamma(r1 + r2 - c1 + 1)))


def fisher_two_sided(a, b, c, d):
    """Two-sided p: the sum of the probabilities of every table no more likely than the observed one."""
    r1, r2, c1 = a + b, c + d, a + c
    if r1 == 0 or r2 == 0 or c1 == 0 or (b + d) == 0:
        return float("nan")
    lo, hi = max(0, c1 - r2), min(r1, c1)
    obs = _log_hyper(a, r1, r2, c1)
    return min(1.0, math.fsum(math.exp(lp) for lp in (_log_hyper(x, r1, r2, c1) for x in range(lo, hi + 1))
                              if lp <= obs + 1e-9))


def odds_ratio(a, b, c, d):
    """Sample odds ratio (a*d)/(b*c); inf or nan when a cell is empty."""
    if b * c == 0:
        return float("inf") if a * d > 0 else float("nan")
    return (a * d) / (b * c)


if __name__ == "__main__":
    # Self-check against textbook values.
    lo, hi = clopper_pearson(0, 10)
    assert lo == 0.0 and abs(hi - 0.3084971) < 1e-6, (lo, hi)
    lo, hi = clopper_pearson(5, 10)
    assert abs(lo - 0.1870860) < 1e-6 and abs(hi - 0.8129140) < 1e-6, (lo, hi)
    lo, hi = clopper_pearson(10, 10)
    assert abs(lo - 0.6915029) < 1e-6 and hi == 1.0, (lo, hi)
    # Fisher's tea-tasting table [[3,1],[1,3]]: two-sided p = 34/70.
    assert abs(fisher_two_sided(3, 1, 1, 3) - 34 / 70) < 1e-12
    # [[1,9],[11,3]]: two-sided p = 0.002759456 (R fisher.test).
    assert abs(fisher_two_sided(1, 9, 11, 3) - 0.002759456) < 1e-8, fisher_two_sided(1, 9, 11, 3)
    print("stats_util self-check passed")
