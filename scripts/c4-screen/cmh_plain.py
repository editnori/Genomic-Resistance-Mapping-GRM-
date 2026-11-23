"""A second, plain implementation of the C4 adjusted screen (docs/screen-arms.md).

  python3 -I cmh_plain.py <out.json>

The Cochran-Mantel-Haenszel statistic of each row of a presence table on a
fit's fitting samples, over strata, with no continuity correction, in exact
rationals:

  D = sum over informative strata of (a - m r / n)
  V = sum over informative strata of m (n - m) r (n - r) / (n^2 (n - 1))
  statistic = D^2 / V, or 0 when V is 0

where a stratum is informative when it holds two or more fitting samples and
both classes; n and r are its fitting samples and resistant ones; m and a
the row's present ones and present resistant ones. Only fitting samples are
read. The screen keeps the top k by the statistic, ties by the lower row.

It writes SYNTHETIC cases, made up from fixed seeds, with each case's
statistics and top k, to <out.json>: tests/c4_screen/cmh_cases.json, which
the Rust tests in src/learn/screen_arm.rs read and
tests/c4_screen/test_cmh_plain.py holds to this script.
"""

import json
import random
import sys
from fractions import Fraction

FORMAT = "grm-c4-cmh-cases-1"


def strata(labels, fit, lineage):
    """(n, r) of each informative stratum among the fitting samples in ascending
    lineage id, each fitting sample's informative stratum index or None, and the
    count of strata among the fitting samples."""
    by = {}
    for i, f in enumerate(fit):
        if f:
            n, r = by.get(lineage[i], (0, 0))
            by[lineage[i]] = (n + 1, r + labels[i])
    nr, index = [], {}
    for sid in sorted(by):
        n, r = by[sid]
        if n >= 2 and 0 < r < n:
            index[sid] = len(nr)
            nr.append((n, r))
    slot = [index.get(lineage[i]) if f else None for i, f in enumerate(fit)]
    return nr, slot, len(by)


def statistic(row, labels, fit, lineage):
    nr, slot, _ = strata(labels, fit, lineage)
    m = [0] * len(nr)
    a = [0] * len(nr)
    for i, s in enumerate(slot):
        if s is not None and row[i]:
            m[s] += 1
            a[s] += labels[i]
    d = Fraction(0)
    v = Fraction(0)
    for (n, r), mk, ak in zip(nr, m, a):
        d += ak - Fraction(mk * r, n)
        v += Fraction(mk * (n - mk) * r * (n - r), n * n * (n - 1))
    return d * d / v if v else Fraction(0)


def chi_square(row, labels, fit):
    """The release's 2x2 Pearson chi-square on the fitting samples, exactly."""
    a = sum(1 for i, f in enumerate(fit) if f and row[i] and labels[i])
    b = sum(1 for i, f in enumerate(fit) if f and row[i] and not labels[i])
    np_ = sum(1 for i, f in enumerate(fit) if f and labels[i])
    nn = sum(1 for i, f in enumerate(fit) if f and not labels[i])
    c, d = np_ - a, nn - b
    denom = (a + b) * (c + d) * np_ * nn
    return Fraction((np_ + nn) * (a * d - b * c) ** 2, denom) if denom else Fraction(0)


def top(stats, k):
    order = sorted(range(len(stats)), key=lambda f: (-stats[f], f))
    return order[:k]


def counts(row, labels, fit, lineage):
    """The row's (m, a) per informative stratum: rows with equal counts get equal
    statistics in any arithmetic, so only they may tie."""
    nr, slot, _ = strata(labels, fit, lineage)
    m = [0] * len(nr)
    a = [0] * len(nr)
    for i, s in enumerate(slot):
        if s is not None and row[i]:
            m[s] += 1
            a[s] += labels[i]
    return tuple(m), tuple(a)


def make(name, seed, n, n_strata, rows, fit_share, k, one_class=0, lone=0, duplicate=0, ids=None):
    rng = random.Random(seed)
    sids = ids or list(range(n_strata))
    lineage = [rng.choice(sids) for _ in range(n)]
    labels = [int(rng.random() < 0.45) for _ in range(n)]
    fit = [int(rng.random() < fit_share) for _ in range(n)]
    # Strata of one class: every fitting member resistant, or every one susceptible.
    for j in range(one_class):
        sid = 10_000 + j
        for i in rng.sample(range(n), 4):
            lineage[i] = sid
            labels[i] = j % 2
    # Strata of one fitting sample.
    for j in range(lone):
        i = rng.randrange(n)
        lineage[i] = 20_000 + j
        fit[i] = 1
    table = []
    seen = {}
    while len(table) < rows:
        p = rng.uniform(0.05, 0.95)
        row = [int(rng.random() < p) for _ in range(n)]
        s = statistic(row, labels, fit, lineage)
        key = counts(row, labels, fit, lineage)
        # Two rows whose exact statistics tie must have equal counts.
        if any(s == t and key != other for t, other in seen.items()):
            continue
        seen[s] = key
        table.append(row)
    for _ in range(duplicate):
        table.insert(rng.randrange(len(table) + 1), list(table[rng.randrange(len(table))]))
    stats = [statistic(r, labels, fit, lineage) for r in table]
    order = top(stats, k)
    nr, _, all_strata = strata(labels, fit, lineage)
    return {
        "name": name,
        "labels": "".join(map(str, labels)),
        "fit": "".join(map(str, fit)),
        "strata": lineage,
        "rows": ["".join(map(str, r)) for r in table],
        "k": k,
        "strata_in_fit": all_strata,
        "informative": len(nr),
        "statistics": [float(s) for s in stats],
        "top": order,
        "top_at_zero": sum(1 for f in order if stats[f] == 0),
    }


def no_informative(seed):
    """Every stratum among the fitting samples has one class: every statistic is
    0 and the top k are the lowest rows."""
    rng = random.Random(seed)
    n = 24
    lineage = [i % 6 for i in range(n)]
    labels = [int(lineage[i] % 2 == 0) for i in range(n)]
    fit = [1] * n
    table = [[int(rng.random() < 0.5) for _ in range(n)] for _ in range(8)]
    stats = [statistic(r, labels, fit, lineage) for r in table]
    nr, _, all_strata = strata(labels, fit, lineage)
    order = top(stats, 5)
    return {"name": "no-informative-stratum", "labels": "".join(map(str, labels)), "fit": "".join(map(str, fit)),
            "strata": lineage, "rows": ["".join(map(str, r)) for r in table], "k": 5, "strata_in_fit": all_strata,
            "informative": len(nr), "statistics": [float(s) for s in stats], "top": order,
            "top_at_zero": sum(1 for f in order if stats[f] == 0)}


def cases():
    return [
        make("one-stratum", 1, 30, 1, 12, 1.0, 5),
        make("six-strata", 2, 60, 6, 20, 1.0, 8),
        make("fit-subset", 3, 80, 7, 24, 0.75, 10),
        make("one-class-and-lone-strata", 4, 70, 5, 20, 0.8, 9, one_class=2, lone=2),
        make("duplicate-rows", 5, 50, 4, 15, 0.9, 12, duplicate=3),
        make("budget-above-rows", 6, 40, 3, 6, 1.0, 50),
        make("large-ids", 7, 64, 5, 16, 0.85, 7, ids=[4_000_000_000, 7, 123_456, 2, 999_999]),
        no_informative(8),
    ]


def main(argv):
    if len(argv) != 1:
        print(__doc__.strip().splitlines()[2], file=sys.stderr)
        return 2
    doc = {"format": FORMAT, "note": "SYNTHETIC: made-up labels, strata and rows from fixed seeds",
           "cases": cases()}
    with open(argv[0], "w", encoding="utf-8") as f:
        json.dump(doc, f, indent=1)
        f.write("\n")
    print(f"{len(doc['cases'])} SYNTHETIC cases written")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
