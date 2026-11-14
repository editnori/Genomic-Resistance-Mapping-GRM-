"""Exact split gains of a class-weighted Gini tree, for the learner cross-check.

scikit-learn breaks equal gains in a seeded random order and GRM takes the
lowest feature, so a tied node has no single reference answer. These
functions compute every candidate's gain as a fraction, which gives the set
of features with the largest gain. It needs no library.

A tree trained on P resistant and N susceptible samples weighs a resistant
sample T / (2 P) and a susceptible one T / (2 N), T = P + N ("balanced" class
weights; a class with no sample counts as 1, as in `src/learn/tree.rs`). A
node with p resistant and n susceptible samples then has the weighted counts
a = p T / (2 P) and b = n T / (2 N), and its Gini impurity times its weight is

    (a + b) * 2 * a / (a + b) * b / (a + b) = 2 a b / (a + b).

The gain of a split is that quantity at the node minus its sum over the two
children. Sets of samples are ints used as bitsets (panels.py).
"""

from fractions import Fraction


def _count(bits):
    return bin(bits).count("1")


class Weights:
    """The class weights of one tree, from its own training samples."""

    def __init__(self, n_resistant, n_susceptible):
        total = n_resistant + n_susceptible
        self.resistant = Fraction(total, 2 * max(n_resistant, 1))
        self.susceptible = Fraction(total, 2 * max(n_susceptible, 1))
        self.n_resistant, self.n_susceptible = n_resistant, n_susceptible

    def mass(self, p, n):
        """Gini impurity times weight of a node with p resistant and n susceptible samples."""
        a, b = p * self.resistant, n * self.susceptible
        return 2 * a * b / (a + b) if a + b else Fraction(0)

    def value(self, p, n):
        """The node's score: its class-weighted resistant share, 1/2 for an empty node."""
        a, b = p * max(self.n_susceptible, 1), n * max(self.n_resistant, 1)
        return Fraction(a, a + b) if a + b else Fraction(1, 2)


def gain(weights, resistant, susceptible, row):
    """The gain of splitting a node on a feature, or None if a side is empty.

    `resistant` and `susceptible` are the node's training samples by class,
    `row` the feature's presence set.
    """
    p, n = _count(resistant), _count(susceptible)
    pp, pn = _count(resistant & row), _count(susceptible & row)
    if pp + pn == 0 or pp + pn == p + n:
        return None
    return weights.mass(p, n) - weights.mass(pp, pn) - weights.mass(p - pp, n - pn)


def best(weights, resistant, susceptible, rows, candidates):
    """(largest gain, the candidates that reach it, ascending), or (None, []).

    A split needs a positive gain, so a result of zero or less means a leaf.
    """
    top, at = None, []
    for feature in sorted(set(candidates)):
        g = gain(weights, resistant, susceptible, rows[feature])
        if g is None:
            continue
        if top is None or g > top:
            top, at = g, [feature]
        elif g == top:
            at.append(feature)
    return top, at
