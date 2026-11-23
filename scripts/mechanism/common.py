"""Small pieces the mechanism scripts share: the cell floor, rounding, seeded
draws, the species interval, the sign test and Holm.

Python 3.12, standard library only; run scripts with -I.
"""

import csv
import hashlib
import math
import os
import random
import statistics

FLOOR = 5            # a catalogue name needs 5 or more genomes in its row's cell; M4's share floor
RECALL_FLOOR = 10    # a recall is written only over 10 or more resistant genomes
MARGIN = 0.02
RESAMPLES = 4000
LO_INDEX, HI_INDEX = 100, 3899   # the 101st and 3,900th sorted values
MIN_SPECIES = 6      # a set under 6 species gets no species interval


def r12(x):
    return None if x is None else round(x + 0.0, 12)


def shown(n):
    """A genome count as a table writes it: as it is, 0 to 4 included (the
    lead's ruling on question 12, 23 Nov). A rate keeps its floor."""
    if not isinstance(n, int) or isinstance(n, bool):
        raise TypeError(f"a count must be an int, not {n!r}")
    if n < 0:
        raise ValueError("a count is never negative")
    return str(n)


def recall_cell(called, resistant):
    """The recall a table writes: only over 10 or more resistant genomes."""
    return f"{r12(called / resistant):.12f}" if resistant >= RECALL_FLOOR else ""


def seed_of(name):
    return int(hashlib.sha256(name.encode("utf-8")).hexdigest(), 16)


def species_interval(by_species, name, resamples=RESAMPLES):
    """95% interval of the median from resamples of species: `by_species` maps
    a species to its datasets' values. None under MIN_SPECIES species."""
    keys = sorted(by_species)
    if len(keys) < MIN_SPECIES:
        return None, None
    rng = random.Random(seed_of(name))
    meds = []
    for _ in range(resamples):
        vals = []
        for _ in keys:
            vals.extend(by_species[keys[rng.randrange(len(keys))]])
        meds.append(r12(statistics.median(vals)))
    meds.sort()
    return meds[LO_INDEX * resamples // RESAMPLES], meds[HI_INDEX * resamples // RESAMPLES]


def sign_p_two_sided(a, b):
    """Exact two-sided sign test of a against b by a fair coin."""
    n = a + b
    if n == 0:
        return 1.0
    k = min(a, b)
    tail = sum(math.comb(n, i) for i in range(0, k + 1)) / 2 ** n
    return min(1.0, 2 * tail)


def holm(ps):
    """Holm-adjusted p values, in the input order."""
    order = sorted(range(len(ps)), key=lambda i: ps[i])
    out, run = [None] * len(ps), 0.0
    for rank, i in enumerate(order):
        run = max(run, min(1.0, (len(ps) - rank) * ps[i]))
        out[i] = run
    return out


def read_tsv(path):
    with open(path, encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t", quoting=csv.QUOTE_NONE))


def write_tsv(path, header, rows):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            cells = []
            for h in header:
                v = r.get(h, "")
                cells.append("" if v is None else str(v))
            if any("\t" in c or "\n" in c for c in cells):
                raise ValueError(f"{path}: a cell holds a tab or a newline")
            f.write("\t".join(cells) + "\n")


def spearman(xs, ys):
    """Spearman's rank correlation with average ranks for ties; None when
    either side is constant or there are under 3 pairs."""
    if len(xs) != len(ys) or len(xs) < 3:
        return None

    def ranks(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        out = [0.0] * len(v)
        i = 0
        while i < len(v):
            j = i
            while j + 1 < len(v) and v[order[j + 1]] == v[order[i]]:
                j += 1
            for k in range(i, j + 1):
                out[order[k]] = (i + j) / 2 + 1
            i = j + 1
        return out

    rx, ry = ranks(xs), ranks(ys)
    mx, my = sum(rx) / len(rx), sum(ry) / len(ry)
    sxy = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    sxx = sum((a - mx) ** 2 for a in rx)
    syy = sum((b - my) ** 2 for b in ry)
    if sxx == 0 or syy == 0:
        return None
    return sxy / math.sqrt(sxx * syy)
