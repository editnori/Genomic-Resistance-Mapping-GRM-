"""Reads the MIC comparison of plan block C6 part B from the scores of score.py (docs/mic-plan.md).

  python3 -I compare.py --out DIR [--results N] [--lineage N] [--study N] [--species N] SCORE_DIR...

Each SCORE_DIR is one species' `score/` folder from score.py, holding
scores.tsv and checks.tsv, with the run's expect.tsv in the folder above it.
The reading's rules were fixed before any score was read. This script
applies rules 1 to 10 of them:

  1. A dataset's value for a split and arm is the mean of EA-exact
     (ea_exact_count / ea_exact_n of each draw) over seeds 42, 43 and 44.
     The twin is the released model over seeds 45, 46 and 47. Values are
     exact fractions.
  2. A draw is defined for an arm when its row is not unavailable, leaves no
     fold undefined and has ea_exact_n above 0. A dataset enters a pair when
     both sides are defined at all three of their draws. A dataset left out
     is named with its reason. No value is filled in.
  3. Pairs on each split: P1 released minus median, P2 ordered minus
     released, P3 twin minus released.
  4. Per pair: compared, better (difference above 0.02), same (within 0.02
     either way) and worse (below -0.02).
  5. The median difference over the datasets compared, with a 95% interval
     from 4,000 species resamples read at the 101st and 3,900th sorted
     values, by the method and generator of
     scripts/amrfinderplus-comparison/record.py: random.Random seeded with a
     string, randrange over the species sorted by name, every dataset of a
     drawn species, statistics.median.
  6. A one-sided exact sign test that better exceeds worse, the same left
     out, with Holm over P1, P2 and P3 of the lineage split. The study
     split's figures decide nothing.
  7. Part B's release rule: on the lineage split, the count of P1 datasets
     whose difference is under 0.02, over the P1 datasets compared.
  8. Beside it: category agreement and the very major and major error rates,
     per arm and split, as the median over datasets of the three-draw mean
     with the number of datasets where it is defined, and the pooled count
     over the pooled denominator of the three draws.
  9. Label marks: EA-exact reads MIC values only; the other three read the R
     and S labels and each fold's training breakpoint.
 10. Outputs are aggregates: one row per dataset and split, one row per pair
     and split, and a summary. No genome id, sequence or marker.

The reading stops, and writes nothing, when an expect.tsv, checks.tsv or
scores.tsv lacks a row, holds a row twice or holds a row it should not, when
a checks.tsv row is not ok, or when a total given on the command line
differs. Standard library only; run it as python3 -I.
"""

import argparse
import csv
import math
import os
import random
import statistics
import sys
from fractions import Fraction

SAME = Fraction(1, 50)
RESAMPLES = 4000
LOW_RANK = round(0.025 * (RESAMPLES - 1))  # 100: the 101st and the 3,900th of the sorted draws
ALPHA = Fraction(1, 20)
# Fixed at 14:48Z on 22 Nov 2026, before any score was read. The reading's rules name the
# generator of record.py, which seeds it with a string, and fix no string.
SEED_STRING = "c6read-20261122|<split>|<pair>|ea_exact"
SPLITS = ["lineage", "study"]
DRAW_SEEDS = ["42", "43", "44"]
TWIN_SEEDS = ["45", "46", "47"]
SEEDS = DRAW_SEEDS + TWIN_SEEDS
# The rows each dataset and split must hold: (seed, arm).
WANT_ROWS = [(s, a) for s in DRAW_SEEDS for a in ("released", "ordered", "median")] + \
            [(s, a) for s in TWIN_SEEDS for a in ("released", "median")]
EXPECT_ARMS = {s: "released,ordered" for s in DRAW_SEEDS} | {s: "released" for s in TWIN_SEEDS}
# Each value of rule 1: (name, score arm, seeds).
VALUES = [("released", "released", DRAW_SEEDS), ("ordered", "ordered", DRAW_SEEDS),
          ("median", "median", DRAW_SEEDS), ("twin", "released", TWIN_SEEDS)]
PAIRS = [("P1", "released", "median"), ("P2", "ordered", "released"), ("P3", "twin", "released")]
# Rule 8: (name, count column, denominator column).
BESIDE = [("category_agreement", "category_count", "category_n"),
          ("very_major_rate", "very_major_count", "very_major_n"),
          ("major_rate", "major_count", "major_n")]
EXPECT_COLUMNS = ["dataset", "species_slug", "split", "plan", "seed", "arms", "file"]
CHECK_COLUMNS = [
    "dataset", "split", "plan", "seed", "file", "exists", "format", "fingerprint",
    "plan_and_arms", "fold_counts", "median_check", "predictions", "ok",
]
SCORE_COLUMNS = [
    "dataset", "species_slug", "antibiotic", "split", "plan", "seed", "draw", "arm",
    "unavailable", "folds", "folds_scored", "folds_undefined", "folds_without_breakpoint",
    "ea_exact", "ea_exact_count", "ea_exact_n",
    "category_agreement", "category_count", "category_n",
    "very_major_rate", "very_major_count", "very_major_n",
    "major_rate", "major_count", "major_n",
    "ea_exact_reads", "category_reads",
]
EA_READS = "MIC values only"
BESIDE_READS = "R and S labels and the training breakpoint"
LABELS = ("not M: the same under label rules (a), (b) and (d); under (c) a dataset holding a pair "
          "with an R and an S record is checked again")


class Stop(Exception):
    """An input the reading does not read on from."""


# ---------------------------------------------------------------- inputs

def read_tsv(path, columns):
    if not os.path.isfile(path):
        raise Stop(f"{path}: missing")
    with open(path, newline="") as f:
        reader = csv.reader(f, delimiter="\t", quoting=csv.QUOTE_NONE)
        rows = list(reader)
    if not rows or rows[0] != columns:
        raise Stop(f"{path}: columns {rows[0] if rows else []}, want {columns}")
    out = []
    for i, cells in enumerate(rows[1:], start=2):
        if len(cells) != len(columns):
            raise Stop(f"{path}: line {i} has {len(cells)} cells, want {len(columns)}")
        out.append(dict(zip(columns, cells)))
    return out


def keyed(rows, key, path):
    out = {}
    for r in rows:
        k = tuple(r[c] for c in key)
        if k in out:
            raise Stop(f"{path}: row {k} twice")
        out[k] = r
    return out


def read_species(score_dir):
    """One species' rows: (species, {(dataset, split): {(seed, arm): score row}}, results)."""
    expect_path = os.path.join(os.path.dirname(os.path.normpath(score_dir)), "expect.tsv")
    checks_path = os.path.join(score_dir, "checks.tsv")
    scores_path = os.path.join(score_dir, "scores.tsv")
    expect = read_tsv(expect_path, EXPECT_COLUMNS)
    checks = read_tsv(checks_path, CHECK_COLUMNS)
    scores = read_tsv(scores_path, SCORE_COLUMNS)
    if not expect:
        raise Stop(f"{expect_path}: no result listed")
    species = {r["species_slug"] for r in expect}
    if len(species) != 1:
        raise Stop(f"{expect_path}: species {sorted(species)}, want one")
    (slug,) = species

    exp = keyed(expect, ["dataset", "split", "seed"], expect_path)
    plan_of = {}
    for (dataset, split, seed), r in exp.items():
        if split not in SPLITS or seed not in SEEDS:
            raise Stop(f"{expect_path}: {dataset} {split} seed {seed} is not in the reading")
        if r["arms"] != EXPECT_ARMS[seed]:
            raise Stop(f"{expect_path}: {dataset} {split} seed {seed} fits {r['arms']}, want {EXPECT_ARMS[seed]}")
        plan_of[(dataset, split, seed)] = r["plan"]
    sets = sorted({(d, s) for d, s, _ in exp})
    for d, s in sets:
        for seed in SEEDS:
            if (d, s, seed) not in exp:
                raise Stop(f"{expect_path}: {d} {s} lacks seed {seed}")

    chk = keyed(checks, ["dataset", "split", "plan", "seed", "file"], checks_path)
    want_checks = {tuple(r[c] for c in ["dataset", "split", "plan", "seed", "file"]) for r in expect}
    for k in sorted(want_checks - set(chk)):
        raise Stop(f"{checks_path}: lacks the row of {k}")
    for k in sorted(set(chk) - want_checks):
        raise Stop(f"{checks_path}: holds {k}, which expect.tsv does not list")
    for k, r in sorted(chk.items()):
        if r["ok"] != "yes":
            raise Stop(f"{checks_path}: {k} is not ok")

    sc = keyed(scores, ["dataset", "split", "seed", "arm"], scores_path)
    want_scores = {(d, s, seed, arm) for d, s in sets for seed, arm in WANT_ROWS}
    for k in sorted(want_scores - set(sc)):
        raise Stop(f"{scores_path}: lacks the row of {k}")
    for k in sorted(set(sc) - want_scores):
        raise Stop(f"{scores_path}: holds {k}, which the reading does not list")
    by_set = {}
    for (d, s, seed, arm), r in sc.items():
        if r["species_slug"] != slug:
            raise Stop(f"{scores_path}: {d} {s} seed {seed} {arm} is species {r['species_slug']}")
        if r["plan"] != plan_of[(d, s, seed)]:
            raise Stop(f"{scores_path}: {d} {s} seed {seed} {arm} has plan {r['plan']}, "
                       f"expect.tsv {plan_of[(d, s, seed)]}")
        if r["draw"] != ("draw" if seed in DRAW_SEEDS else "twin"):
            raise Stop(f"{scores_path}: {d} {s} seed {seed} {arm} is marked {r['draw']}")
        by_set.setdefault((d, s), {})[(seed, arm)] = r
    for (d, s), rows in by_set.items():
        names = {r["antibiotic"] for r in rows.values()}
        if len(names) != 1 or "" in names:
            raise Stop(f"{scores_path}: {d} {s} has antibiotic {sorted(names)}")
    return slug, by_set, len(expect)


def read_all(score_dirs):
    """{(dataset, split): rows} over every species, the species, and the result count."""
    by_set, species, results = {}, [], 0
    for score_dir in score_dirs:
        slug, rows, n = read_species(score_dir)
        if slug in species:
            raise Stop(f"{score_dir}: species {slug} twice")
        species.append(slug)
        for key in rows:
            if key in by_set:
                raise Stop(f"{score_dir}: {key} twice")
        by_set.update(rows)
        results += n
    return by_set, sorted(species), results


# ---------------------------------------------------------------- rules

def integer(cell, what):
    try:
        return int(cell)
    except ValueError:
        raise Stop(f"{what}: {cell!r} is not a count") from None


def draw_ea(row):
    """Rule 2 for one draw: (EA-exact as a fraction, None) or (None, why not)."""
    where = f"{row['dataset']} {row['split']} seed {row['seed']} {row['arm']}"
    if row["unavailable"]:
        return None, f"unavailable at seed {row['seed']}"
    if integer(row["folds_undefined"], where) > 0:
        return None, f"{row['folds_undefined']} fold(s) undefined at seed {row['seed']}"
    n = integer(row["ea_exact_n"], where)
    if n <= 0:
        return None, f"no exact MIC tested at seed {row['seed']}"
    return Fraction(integer(row["ea_exact_count"], where), n), None


def arm_value(rows, arm, seeds):
    """Rule 1: (three-draw mean, count sum, denominator sum, None) or (None, ..., why not)."""
    values, why = [], []
    count = n = 0
    for seed in seeds:
        row = rows[(seed, arm)]
        v, reason = draw_ea(row)
        if v is None:
            why.append(reason)
        else:
            values.append(v)
            count += int(row["ea_exact_count"])
            n += int(row["ea_exact_n"])
    if why:
        return None, None, None, "; ".join(why)
    return sum(values) / len(values), count, n, None


def beside_value(rows, arm, seeds, count_col, n_col):
    """Rule 8 for a dataset: (three-draw mean, count sum, denominator sum) or Nones.
    A draw's measure is defined when the row is not unavailable, leaves no fold
    undefined and has the measure's count over a denominator above 0."""
    values = []
    count = n = 0
    for seed in seeds:
        row = rows[(seed, arm)]
        where = f"{row['dataset']} {row['split']} seed {seed} {arm}"
        if row["unavailable"] or integer(row["folds_undefined"], where) > 0:
            return None, None, None
        if row[count_col] == "" or integer(row[n_col], where) <= 0:
            return None, None, None
        c, d = integer(row[count_col], where), integer(row[n_col], where)
        values.append(Fraction(c, d))
        count, n = count + c, n + d
    return sum(values) / len(values), count, n


def count_class(diff):
    """Rule 4: better above 0.02, worse below -0.02, else same."""
    if diff > SAME:
        return "better"
    if diff < -SAME:
        return "worse"
    return "same"


def sign_p(better, worse):
    """Rule 6: one-sided exact binomial P(X >= better), X ~ Bin(better + worse, 1/2)."""
    n = better + worse
    if n == 0:
        return Fraction(1)
    return Fraction(sum(math.comb(n, i) for i in range(better, n + 1)), 2 ** n)


def holm(ps):
    """Holm's adjusted p values, in the order given."""
    m = len(ps)
    order = sorted(range(m), key=lambda i: (ps[i], i))
    out, running = [None] * m, Fraction(0)
    for rank, i in enumerate(order):
        running = max(running, min(Fraction(1), (m - rank) * ps[i]))
        out[i] = running
    return out


def seed_string(split, pair):
    return SEED_STRING.replace("<split>", split).replace("<pair>", pair)


def species_interval(by_species, seed):
    """Rule 5, as record.py's redraw: the medians of 4,000 species resamples,
    sorted and read at the 101st and the 3,900th. None with no species."""
    sp = sorted(by_species)
    if not sp:
        return None
    flt = {s: [float(x) for x in by_species[s]] for s in sp}
    rng = random.Random(seed)
    meds = []
    for _ in range(RESAMPLES):
        xs = []
        for _ in range(len(sp)):
            s = sp[rng.randrange(len(sp))]
            xs.extend(flt[s])
        meds.append(statistics.median(xs))
    meds.sort()
    return meds[LOW_RANK], meds[RESAMPLES - 1 - LOW_RANK]


def sign_word(x):
    return "above 0" if x > 0 else "below 0" if x < 0 else "0"


# ---------------------------------------------------------------- reading

def species_of(dataset):
    return dataset.split("__", 1)[0]


def one_count(rows, col):
    vals = sorted({integer(r[col], f"{r['dataset']} {r['split']} {col}") for r in rows.values()})
    return str(vals[0]) if len(vals) == 1 else f"{vals[0]}-{vals[-1]}"


def read_sets(by_set):
    """Rules 1, 2, 4 and 8 for every dataset and split."""
    out = []
    for (dataset, split), rows in sorted(by_set.items(), key=lambda kv: (SPLITS.index(kv[0][1]), kv[0][0])):
        one = next(iter(rows.values()))
        rec = {"dataset": dataset, "species_slug": one["species_slug"], "antibiotic": one["antibiotic"],
               "split": split, "mic_genomes": one_count(rows, "category_n"),
               "exact_mics": one_count(rows, "ea_exact_n"), "values": {}, "pairs": {}, "beside": {}}
        for name, arm, seeds in VALUES:
            rec["values"][name] = arm_value(rows, arm, seeds)
            for m, count_col, n_col in BESIDE:
                rec["beside"][(name, m)] = beside_value(rows, arm, seeds, count_col, n_col)
        for pair, first, second in PAIRS:
            a, b = rec["values"][first], rec["values"][second]
            why = []
            if a[0] is None:
                why.append(f"{first}: {a[3]}")
            if b[0] is None:
                why.append(f"{second}: {b[3]}")
            if why:
                rec["pairs"][pair] = (None, "left out", "; ".join(why))
            else:
                diff = a[0] - b[0]
                rec["pairs"][pair] = (diff, count_class(diff), "")
        out.append(rec)
    return out


def read_pairs(sets):
    """Rules 3 to 7 for each pair and split."""
    out = []
    for split in SPLITS:
        mine = [r for r in sets if r["split"] == split]
        for pair, first, second in PAIRS:
            diffs = [(r, r["pairs"][pair][0]) for r in mine if r["pairs"][pair][0] is not None]
            left = [r for r in mine if r["pairs"][pair][0] is None]
            better = sum(1 for _, d in diffs if d > SAME)
            worse = sum(1 for _, d in diffs if d < -SAME)
            row = {"split": split, "pair": pair, "first": first, "second": second,
                   "datasets": len(mine), "compared": len(diffs), "left_out": len(left),
                   "better": better, "same": len(diffs) - better - worse, "worse": worse,
                   "species": len({r["species_slug"] for r, _ in diffs}),
                   "median": None, "low": None, "high": None,
                   "sign_p": sign_p(better, worse), "holm_p": None,
                   "under": None, "seed_string": seed_string(split, pair),
                   "left": [(r["dataset"], r["pairs"][pair][2]) for r in left]}
            if diffs:
                row["median"] = statistics.median([d for _, d in diffs])
                by_species = {}
                for r, d in diffs:
                    by_species.setdefault(r["species_slug"], []).append(d)
                row["low"], row["high"] = species_interval(by_species, row["seed_string"])
            if split == "lineage" and pair == "P1":
                row["under"] = sum(1 for _, d in diffs if d < SAME)
            out.append(row)
    lineage = [r for r in out if r["split"] == "lineage"]
    for r, p in zip(lineage, holm([r["sign_p"] for r in lineage])):
        r["holm_p"] = p
    return out


def read_beside(sets):
    """Rule 8 per value and split."""
    out = []
    for split in SPLITS:
        mine = [r for r in sets if r["split"] == split]
        for name, _, _ in VALUES:
            for m, _, _ in BESIDE:
                got = [r["beside"][(name, m)] for r in mine if r["beside"][(name, m)][0] is not None]
                out.append({"split": split, "arm": name, "measure": m, "datasets": len(mine),
                            "defined": len(got),
                            "median": statistics.median([g[0] for g in got]) if got else None,
                            "count": sum(g[1] for g in got), "n": sum(g[2] for g in got)})
    return out


# ---------------------------------------------------------------- outputs

def f6(x):
    if x is None:
        return ""
    s = "%.6f" % float(x)
    return "0.000000" if s == "-0.000000" else s


def p6(x):
    return "" if x is None else "%.6g" % float(x)


def set_columns():
    cols = ["dataset", "species_slug", "antibiotic", "split", "mic_genomes", "exact_mics"]
    for name, _, _ in VALUES:
        cols += [f"{name}_ea_exact", f"{name}_ea_count", f"{name}_ea_n", f"{name}_defined"]
    for pair, _, _ in PAIRS:
        cols += [f"{pair}_difference", f"{pair}_class", f"{pair}_left_out"]
    for name, _, _ in VALUES:
        for m, _, _ in BESIDE:
            cols += [f"{name}_{m}", f"{name}_{m}_count", f"{name}_{m}_n"]
    return cols + ["ea_exact_reads", "beside_reads", "labels"]


def set_row(r):
    out = {k: r[k] for k in ["dataset", "species_slug", "antibiotic", "split", "mic_genomes", "exact_mics"]}
    for name, _, _ in VALUES:
        v, c, n, why = r["values"][name]
        out[f"{name}_ea_exact"] = f6(v)
        out[f"{name}_ea_count"] = "" if c is None else c
        out[f"{name}_ea_n"] = "" if n is None else n
        out[f"{name}_defined"] = "yes" if why is None else why
    for pair, _, _ in PAIRS:
        d, cls, why = r["pairs"][pair]
        out[f"{pair}_difference"], out[f"{pair}_class"], out[f"{pair}_left_out"] = f6(d), cls, why
    for name, _, _ in VALUES:
        for m, _, _ in BESIDE:
            v, c, n = r["beside"][(name, m)]
            out[f"{name}_{m}"] = f6(v)
            out[f"{name}_{m}_count"] = "" if c is None else c
            out[f"{name}_{m}_n"] = "" if n is None else n
    out.update({"ea_exact_reads": EA_READS, "beside_reads": BESIDE_READS, "labels": LABELS})
    return out


PAIR_COLUMNS = [
    "split", "pair", "first", "second", "datasets", "compared", "left_out", "better", "same", "worse",
    "species", "median_difference", "median_exact", "interval_low", "interval_high", "median_sign",
    "sign_p_one_sided", "holm_p", "holm_below_0.05", "under_0.02", "under_0.02_over_half",
    "decides", "seed_string", "reads",
]


def pair_row(r):
    lineage = r["split"] == "lineage"
    return {
        "split": r["split"], "pair": r["pair"], "first": r["first"], "second": r["second"],
        "datasets": r["datasets"], "compared": r["compared"], "left_out": r["left_out"],
        "better": r["better"], "same": r["same"], "worse": r["worse"], "species": r["species"],
        "median_difference": f6(r["median"]), "median_exact": "" if r["median"] is None else str(r["median"]),
        "interval_low": f6(r["low"]), "interval_high": f6(r["high"]),
        "median_sign": "" if r["median"] is None else sign_word(r["median"]),
        "sign_p_one_sided": p6(r["sign_p"]), "holm_p": p6(r["holm_p"]),
        "holm_below_0.05": "" if r["holm_p"] is None else ("yes" if r["holm_p"] < ALPHA else "no"),
        "under_0.02": "" if r["under"] is None else r["under"],
        "under_0.02_over_half": "" if r["under"] is None else ("yes" if 2 * r["under"] > r["compared"] else "no"),
        "decides": "Holm family" if lineage else "no, shown beside",
        "seed_string": r["seed_string"], "reads": EA_READS,
    }


def write_tsv(path, columns, rows):
    with open(path, "w") as f:
        f.write("\t".join(columns) + "\n")
        for r in rows:
            f.write("\t".join(str(r[c]) for c in columns) + "\n")


def summary(species, results, by_set, sets, pairs, beside):
    n_split = {s: sum(1 for r in sets if r["split"] == s) for s in SPLITS}
    lines = [
        "MIC comparison, plan block C6 part B: the reading of rules 1 to 10",
        f"results: {results} in {len(species)} species: {', '.join(species)}",
        f"score rows: {sum(len(v) for v in by_set.values())}",
        f"datasets: lineage {n_split['lineage']}, study {n_split['study']}",
        "",
        "pairs (better above 0.02, same within 0.02, worse below -0.02; difference = first minus second)",
    ]
    for r in pairs:
        lines.append(
            f"{r['split']} {r['pair']} {r['first']} - {r['second']}: compared {r['compared']} of {r['datasets']}, "
            f"left out {r['left_out']}, better {r['better']}, same {r['same']}, worse {r['worse']}; "
            f"median {f6(r['median'])} ({f6(r['low'])} to {f6(r['high'])}) over {r['species']} species, "
            f"sign {sign_word(r['median']) if r['median'] is not None else 'none'}; "
            f"one-sided sign p {p6(r['sign_p'])}"
            + (f", Holm p {p6(r['holm_p'])}" if r["holm_p"] is not None else ", decides nothing"))
    p1 = next(r for r in pairs if r["split"] == "lineage" and r["pair"] == "P1")
    lines += [
        "",
        f"release rule (lineage P1 under 0.02): {p1['under']} of {p1['compared']}; "
        f"more than half: {'yes' if 2 * p1['under'] > p1['compared'] else 'no'}",
        "",
        "beside it (median over datasets of the three-draw mean; datasets defined; pooled count / pooled n)",
    ]
    for b in beside:
        lines.append(f"{b['split']} {b['arm']} {b['measure']}: median {f6(b['median'])}, "
                     f"defined {b['defined']} of {b['datasets']}, pooled {b['count']}/{b['n']}")
    lines += ["", "left out"]
    left = [(r["split"], r["pair"], d, why) for r in pairs for d, why in r["left"]]
    lines += [f"{s} {p} {d}: {why}" for s, p, d, why in left] or ["none"]
    lines += ["", "labels",
              f"EA-exact reads {EA_READS}.",
              f"Category agreement and both error rates read the {BESIDE_READS}.",
              f"Every dataset: {LABELS}."]
    return "\n".join(lines) + "\n"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    ap.add_argument("--out", required=True)
    ap.add_argument("--results", type=int)
    ap.add_argument("--lineage", type=int)
    ap.add_argument("--study", type=int)
    ap.add_argument("--species", type=int)
    ap.add_argument("score_dirs", nargs="+")
    args = ap.parse_args(argv)
    names = ["datasets.tsv", "pairs.tsv", "summary.txt"]
    try:
        for name in names:
            if os.path.exists(os.path.join(args.out, name)):
                raise Stop(f"{os.path.join(args.out, name)} exists")
        by_set, species, results = read_all(args.score_dirs)
        got = {"results": results, "species": len(species),
               "lineage": sum(1 for _, s in by_set if s == "lineage"),
               "study": sum(1 for _, s in by_set if s == "study")}
        for k in ["results", "lineage", "study", "species"]:
            want = getattr(args, k)
            if want is not None and got[k] != want:
                raise Stop(f"{k}: {got[k]}, want {want}")
        sets = read_sets(by_set)
        pairs = read_pairs(sets)
        beside = read_beside(sets)
    except Stop as e:
        print(f"stop: {e}", file=sys.stderr)
        return 2
    os.makedirs(args.out, exist_ok=True)
    write_tsv(os.path.join(args.out, "datasets.tsv"), set_columns(), [set_row(r) for r in sets])
    write_tsv(os.path.join(args.out, "pairs.tsv"), PAIR_COLUMNS, [pair_row(r) for r in pairs])
    text = summary(species, results, by_set, sets, pairs, beside)
    with open(os.path.join(args.out, "summary.txt"), "w") as f:
        f.write(text)
    print(text, end="")
    return 0


if __name__ == "__main__":
    sys.exit(main())
