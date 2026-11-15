"""An independent set covering machine for the learner cross-check.

Written from the method's definition (Marchand and Shawe-Taylor 2002), not
from GRM's code, and in exact arithmetic. It needs no library.

A conjunction calls a sample resistant when every rule holds, a disjunction
when any rule holds. A rule is the presence or the absence of one feature.
The fit is greedy. Each step takes the rule with the largest utility,
`covered - p * lost`: for a conjunction, `covered` counts the susceptible
training samples the rule newly rules out and `lost` the resistant training
samples it wrongly rules out. A disjunction is the same search with the
classes swapped and every chosen rule negated.

What the definition leaves open is taken from the documentation of
`src/learn/scm.rs`:

- equal utilities: fewer lost, then the lower feature, then presence before
  absence;
- the fit stops at `max_rules`, when no sample of the class to rule out
  remains, or when the best rule rules out none.

Utilities are compared as fractions of the float `p`, so no rounding can
decide a step here. Sets of samples are ints used as bitsets (panels.py).
"""

from fractions import Fraction

KINDS = ("conjunction", "disjunction")


def _count(bits):
    return bin(bits).count("1")


def fit(rows, y, train, n, kind, p, max_rules):
    """The greedy path on the samples of `train`.

    Returns (rules, steps): each rule is (feature, presence) as the model
    stores it, each step is (covered, lost, utility), the utility a Fraction.
    """
    if kind not in KINDS:
        raise ValueError(f"unknown kind {kind!r}")
    everyone = (1 << n) - 1
    keep, rule_out = train & y, train & ~y & everyone
    if kind == "disjunction":
        keep, rule_out = rule_out, keep
    penalty = Fraction(p)
    rules, steps = [], []
    while len(rules) < max_rules and rule_out:
        best = None
        for feature, row in enumerate(rows):
            for presence in (True, False):
                holds = row if presence else ~row & everyone
                covered = _count(rule_out & ~holds)
                lost = _count(keep & ~holds)
                # Larger is better in every field: True sorts above False.
                rank = (covered - penalty * lost, -lost, -feature, presence)
                if best is None or rank > best[0]:
                    best = (rank, covered, lost, feature, presence, holds)
        if best is None:
            break
        rank, covered, lost, feature, presence, holds = best
        if covered == 0:
            break
        rule_out &= holds
        keep &= holds
        rules.append((feature, presence if kind == "conjunction" else not presence))
        steps.append((covered, lost, rank[0]))
    return rules, steps


def predict(rows, n, kind, rules):
    """The samples the model calls resistant."""
    everyone = (1 << n) - 1
    called = everyone if kind == "conjunction" else 0
    for feature, presence in rules:
        holds = rows[feature] if presence else ~rows[feature] & everyone
        called = called & holds if kind == "conjunction" else called | holds
    return called
