"""A plain computation of stage 2's figures straight from a made-up plan, apart from both scripts, beside the
export of scripts/matched-accuracy/stage2.py at 3 seeds of made_up2.randomised().

The plan gives each group's rows as made_up2 writes them, so this file reads no report, no row list and nothing
stage 1 wrote; it calls no function of stage 1 or stage 2. Every value here is made up."""

import fractions
import unittest

import numpy

from tests.matched_accuracy import made_up2, stage2_common as c

F = fractions.Fraction
SEEDS = (1, 2, 3)


def row_mean(five):
    total = 0.0
    for x in five:          # scm, cart, rf, gbt, lasso, left to right
        total = total + x
    return total / 5


def plain(plan, reading, keep=lambda ds, g: True, cell=None):
    """{"overall": value, "species": {name: (value, pairs)}, "dataset": {name: (value, pairs)}}.

    reading B leaves out rows with an I record; a group needs an R row and an S row; keep(ds, g) picks groups;
    cell, a pair of catalogue calls, keeps only the pairs of that cell."""
    species = {}
    for k, datasets in enumerate(plan):
        sp = f"synthia_{made_up2.NAMES[k]}"
        for ds in datasets:
            values, pairs = [], 0
            for g in ds["groups"]:
                if not keep(ds, g):
                    continue
                rows = [r for r in g["rows"] if not (reading == "B" and r["intermediate"])]
                rs = [r for r in rows if r["resistant"]]
                ss = [r for r in rows if not r["resistant"]]
                if not rs or not ss:
                    continue
                half = n = 0
                for r in rs:
                    for s in ss:
                        if cell is not None and (not ds["rule"] or (r["catalog"], s["catalog"]) != cell):
                            continue
                        a, b = row_mean(r["five"]), row_mean(s["five"])
                        half += 2 if a > b else 1 if a == b else 0
                        n += 1
                if n:
                    values.append(F(half, 2 * n))
                    pairs += n
            if values:
                species.setdefault(sp, {})[f"{sp}__{ds['drug']}"] = (sum(values, F(0)) / len(values), pairs)
    out = {"species": {}, "dataset": {}}
    for sp, ds in species.items():
        out["dataset"].update(ds)
        out["species"][sp] = (sum((v for v, _ in ds.values()), F(0)) / len(ds), sum(p for _, p in ds.values()))
    out["overall"] = sum((v for v, _ in out["species"].values()), F(0)) / len(out["species"])
    return out


def no_null_fit(ds, g):
    return not any(g["fold"] in folds for folds in ds["null_folds"].values())


def bounds(values):
    means = numpy.array([float(v) for v in values])
    rng = numpy.random.Generator(numpy.random.PCG64(42))
    draws = means[rng.integers(0, len(means), size=(10000, len(means)))].mean(axis=1)
    lower, upper = (float(b) for b in numpy.percentile(draws, [2.5, 97.5]))
    word = "above" if lower > 0.5 else "below" if upper < 0.5 else "unresolved"
    return [str(round(F(lower) * 1000000)), str(round(F(upper) * 1000000)), word]


def ppm(value, pairs=10):
    return "" if pairs < 10 else str(round(value * 1000000))


class Plain(c.Case, unittest.TestCase):
    def test_the_script_agrees_with_the_plain_computation_at_3_seeds(self):
        for seed in SEEDS:
            plan = made_up2.randomised(seed)
            run = self.run_stage2(plan)
            self.assert_passes(run)
            measure = {(r["reading"], r["level"], r["name"]): r["estimate_ppm"] for r in run.table("measure.tsv")}
            nulls = {(r["reading"], r["level"], r["name"]): r["estimate_ppm"] for r in run.table("null_fits.tsv")}
            cells = {(r["reading"], r["resistant_member_call"], r["susceptible_member_call"]): r["estimate_ppm"]
                     for r in run.table("cells.tsv")}
            interval = {(r["reading"], r["read"]): [r["lower_ppm"], r["upper_ppm"], r["chance"]]
                        for r in run.table("interval.tsv") if r["match"] == "primary"}
            for reading in ("A", "B"):
                with self.subTest(seed=seed, reading=reading):
                    p = plain(plan, reading)
                    self.assertEqual(measure[(reading, "overall", "all")], ppm(p["overall"]))
                    for level in ("species", "dataset"):
                        got = {n: v for (rd, lv, n), v in measure.items() if (rd, lv) == (reading, level)}
                        self.assertEqual(got, {n: ppm(v, pairs) for n, (v, pairs) in p[level].items()})
                    self.assertEqual(interval[(reading, "primary")],
                                     bounds(v for _, (v, _) in sorted(p["species"].items())))
                    q = plain(plan, reading, keep=no_null_fit)
                    self.assertEqual(nulls[(reading, "overall", "all")], ppm(q["overall"]))
                    self.assertEqual(interval[(reading, "null_fits")],
                                     bounds(v for _, (v, _) in sorted(q["species"].items())))
                    for first, second in ((True, False), (False, False), (True, True)):
                        name = (reading, c.stage2.CALLS[first], c.stage2.CALLS[second])
                        try:
                            want = ppm(plain(plan, reading, cell=(first, second))["overall"])
                        except ZeroDivisionError:   # no pair of the cell: the cell prints no estimate
                            want = ""
                        self.assertEqual(cells[name], want, name)


if __name__ == "__main__":
    unittest.main()
