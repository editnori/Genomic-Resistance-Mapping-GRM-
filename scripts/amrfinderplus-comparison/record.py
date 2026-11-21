#!/usr/bin/env python3
"""Write and check reference/amrfinderplus-comparison, the record of the AMRFinderPlus comparison.

From 14 to 16 Nov 2026 the AMRFinderPlus program (4.2.7, database 2026-08-07.1) was run on the
genomes of the 253 datasets of release rev14-274-20261112 that have a published rule, and its
rows were read under each report's rule. Its calls were scored beside GRM's known-gene rule and
GRM's five learners on the same genomes and recorded labels. A review then recomputed the figures
and narrowed the reading. This script writes the record of that comparison from the final
scoring's tables per dataset, which hold counts and measures only.

Modes:

    record.py --source FINAL --out DIR     write the record from the final scoring
    record.py --source FINAL --check DIR   fail unless DIR is what the final scoring gives
    record.py --check DIR                  fail unless DIR's derived files, README and SHA256SUMS
                                           are what its row files give

FINAL is the final scoring's directory. The plan's set lists sit beside it in ../plan-lists, and
the plan, its addenda, their hash files and the run logs in ../.. and ../logs. None of them is
public, so only the last mode runs without them. The source modes redraw every species interval
from the plan's seed string and hold it, and every summary cell, to the final scoring's own
summary. Standard library only; run it as python3 -I.
"""

import argparse
import gzip
import hashlib
import json
import math
import pathlib
import random
import re
import statistics
import sys
from fractions import Fraction

LEARNERS = ["scm", "cart", "rf", "gbt", "lasso"]
PICKS = LEARNERS + ["served_best", "highest_of_five"]
SPLITS = ["lineage", "random"]
MEASURES = ["ba", "sens", "spec"]
PR = "program_minus_known_gene_rule"
PROGRAM = "program:tool_default"
RULE = "grm_known_gene_rule"
READINGS = ["tool_default", "tool_default_intact", "genes_core", "genes_all_scopes"]
WITHIN = Fraction(1, 50)
LARGE = Fraction(1, 10)
RESAMPLES = 4000
FEWEST_SPECIES = 6
LOW_RANK = round(0.025 * (RESAMPLES - 1))  # 100: the 101st and the 3900th of the sorted draws

# Shapes that never belong in the record: a genome id and a run of 20 or more bases.
GENOME_ID = re.compile(r"(?<![\d.])\d{2,7}\.\d{3,6}(?![\d.])")
BASES = re.compile(r"(?i)[acgt]{20,}")

COLUMNS = {
    "datasets.tsv": ["dataset", "species", "drug_class", "organism_option", "lineage_split",
                     "read_before_the_plan", "served_learner"],
    "callers.tsv": ["dataset", "caller", "genomes", "scored", "no_call_resistant", "no_call_susceptible",
                    "resistant", "susceptible", "tp", "fn", "tn", "fp", "sensitivity", "sens_exact_lo",
                    "sens_exact_hi", "sens_lineage_lo", "sens_lineage_hi", "specificity", "spec_exact_lo",
                    "spec_exact_hi", "spec_lineage_lo", "spec_lineage_hi", "balanced_accuracy", "ba_lineage_lo",
                    "ba_lineage_hi"],
    "differences.tsv": ["dataset", "reading", "contrast", "measure", "first", "second", "first_minus_second",
                        "lineage_lo", "lineage_hi", "first_only_right", "second_only_right", "mcnemar_p",
                        "calls_agree", "scored"],
    "intervals.tsv": ["set", "contrast", "measure", "median_species_lo", "median_species_hi", "mean_species_lo",
                      "mean_species_hi", "share_species_lo", "share_species_hi"],
    "order.tsv": ["utc", "event", "file", "sha256"],
    "settings.tsv": ["key", "value"],
    "sets.tsv": ["set", "role", "datasets", "species", "definition"],
    "summary.tsv": ["set", "contrast", "measure", "role", "datasets", "species", "ahead", "level", "behind",
                    "median", "q1", "q3", "mean", "sign_test_p", "median_species_lo", "median_species_hi",
                    "mean_species_lo", "mean_species_hi", "own_interval_above_zero", "own_interval_spans_zero",
                    "own_interval_below_zero", "learner_at_most_0.02_ahead", "share_learner_at_most_0.02_ahead",
                    "share_species_lo", "share_species_hi", "within_0.02_either_way"],
    "figures.tsv": ["key", "set", "value", "role"],
}
ROW_FILES = ["datasets.tsv", "callers.tsv", "differences.tsv", "intervals.tsv", "order.tsv", "settings.tsv"]
DERIVED = ["sets.tsv", "summary.tsv", "figures.tsv"]

# The sets, in the order the record gives them: name, role, definition over datasets.tsv.
SETS = [
    ("primary-new-235", "primary", "organism option, lineage split, not read before the plan"),
    ("primary-237", "secondary", "organism option and lineage split"),
    ("rule-M-quinolone-28", "secondary", "quinolone datasets of primary-237"),
    ("rule-M-quinolone-new-26", "secondary", "quinolone datasets of primary-new-235"),
    ("quinolone-organism-option-30", "secondary", "quinolone datasets with an organism option"),
    ("carbapenem-20", "secondary", "carbapenem datasets"),
    ("no-organism-option-12", "secondary", "datasets in a species without an organism option"),
    ("organism-option-no-lineage-split-4", "secondary", "organism option and no lineage split"),
]
SET_NAMES = [s[0] for s in SETS]
PLAN_LISTS = {"primary-new-235", "primary-237", "rule-M-quinolone-28", "quinolone-organism-option-30",
              "carbapenem-20", "no-organism-option-12", "organism-option-no-lineage-split-4"}


def members(name, d):
    """Whether dataset row d belongs to set name."""
    org, lin = d["organism_option"] == "yes", d["lineage_split"] == "yes"
    new, quin = d["read_before_the_plan"] == "no", d["drug_class"] == "quinolone"
    return {
        "primary-new-235": org and lin and new,
        "primary-237": org and lin,
        "rule-M-quinolone-28": org and lin and quin,
        "rule-M-quinolone-new-26": org and lin and new and quin,
        "quinolone-organism-option-30": org and quin,
        "carbapenem-20": d["drug_class"] == "carbapenem",
        "no-organism-option-12": not org,
        "organism-option-no-lineage-split-4": org and not lin,
    }[name]


def contrast_name(split, pick):
    return "learner_minus_program:%s:%s" % (split, pick)


def targets():
    """Every (set, contrast, measure) the summary gives, in its order."""
    out = []
    for name in SET_NAMES:
        if name == "primary-new-235":
            out += [(name, PR, m) for m in MEASURES]
            out += [(name, contrast_name(s, p), m) for s in SPLITS for p in PICKS for m in MEASURES]
        elif name == "organism-option-no-lineage-split-4":
            out += [(name, PR, "ba"), (name, contrast_name("random", "served_best"), "ba")]
        else:
            out += [(name, PR, m) for m in (MEASURES if name == "primary-237" else ["ba"])]
            out += [(name, contrast_name("lineage", p), "ba") for p in PICKS]
            out.append((name, contrast_name("random", "served_best"), "ba"))
    return out


def role(name, con, m):
    """primary, secondary or exploratory, as the plan fixed them before any result of stages 2 to 4."""
    if con == PR:
        if name == "primary-new-235":
            return "primary"
        return "exploratory" if name == "rule-M-quinolone-new-26" else "secondary"
    _, split, pick = con.split(":")
    if name == "primary-new-235":
        if m != "ba" or pick == "highest_of_five":
            return "exploratory"
        return "primary" if split == "lineage" else "secondary"
    if name == "rule-M-quinolone-28":
        return "secondary"
    if name == "rule-M-quinolone-new-26":
        return "secondary" if con == contrast_name("lineage", "served_best") else "exploratory"
    return "exploratory" if pick == "highest_of_five" else "secondary"


# The settings, each checked in the source modes against the text named beside it.
SETTINGS = [
    ("program", "AMRFinderPlus"),
    ("program_version", "4.2.7"),
    ("database_version", "2026-08-07.1"),
    ("newest_database_not_run", "2026-10-29.1"),
    ("release", "rev14-274-20261112"),
    ("datasets_scored", "253"),
    ("genomes_read", "23751"),
    ("score_cutoff", "0.5"),
    ("resamples", str(RESAMPLES)),
    ("fewest_species", str(FEWEST_SPECIES)),
    ("interval_level", "0.95"),
    ("interval_low_rank", str(LOW_RANK + 1)),
    ("interval_high_rank", str(RESAMPLES - LOW_RANK)),
    ("generator", "random.Random(seed string)"),
    ("draw", "randrange over the species of the set sorted by name"),
    ("median", "statistics.median"),
    ("mean", "statistics.fmean"),
    ("within", "0.02"),
    ("large", "0.10"),
    ("rule_m_needs", "14"),
]
SETTINGS_IN_PLAN = {
    "newest_database_not_run": "2026-10-29.1",
    "release": "Release `rev14-274-20261112`",
    "score_cutoff": "at 0.5 with `>=`",
    "resamples": "4,000 resamples",
    "fewest_species": "fewer than 6 species",
    "within": "<= 0.02 in balanced accuracy",
    "rule_m_needs": "at least 14 of the 28",
}
# Not named by the plan: the generator, the draw and the two statistics the review found.
ORDER = [
    ("plan hashed", "PLAN.md", "PLAN.sha256", "results/logs/plan-hash-time.txt"),
    ("first addendum hashed", "PLAN-ADDENDUM-20261115.md", "PLAN-ADDENDUM-20261115.sha256",
     "results/logs/plan-addendum-hash-time.txt"),
    ("second addendum hashed", "PLAN-ADDENDUM-20261115-2.md", "PLAN-ADDENDUM-20261115-2.sha256",
     "results/logs/plan-addendum-2-hash-time.txt"),
    ("third addendum hashed", "PLAN-ADDENDUM-20261116.md", "PLAN-ADDENDUM-20261116.sha256",
     "results/logs/plan-addendum-3-hash-time.txt"),
]
FINAL_LOG = "results/logs/analyze-final.log"


class Stop(Exception):
    """An input or a record the script will not write from."""


# ---------------------------------------------------------------- tables

def parse(text):
    lines = text.rstrip("\n").split("\n")
    head = lines[0].split("\t")
    rows = []
    for line in lines[1:]:
        cells = line.split("\t")
        if len(cells) != len(head):
            raise Stop("a row has %d cells where the header has %d" % (len(cells), len(head)))
        rows.append(dict(zip(head, cells)))
    return head, rows


def tsv(name, rows):
    cols = COLUMNS[name]
    out = ["\t".join(cols)]
    for r in rows:
        out.append("\t".join(str(r[c]) for c in cols))
    return "\n".join(out) + "\n"


def f6(x):
    s = "%.6f" % x
    return "0.000000" if s == "-0.000000" else s


def opt6(x):
    return "" if x is None else f6(x)


# ---------------------------------------------------------------- measures

def rates(r):
    tp, fn, tn, fp = (int(r[k]) for k in ("tp", "fn", "tn", "fp"))
    sens = Fraction(tp, tp + fn) if tp + fn else Fraction(0)
    spec = Fraction(tn, tn + fp) if tn + fp else Fraction(0)
    return {"sens": sens, "spec": spec, "ba": (sens + spec) / 2}


def sign_p(ahead, behind):
    n, k = ahead + behind, min(ahead, behind)
    if not n:
        return 1.0
    return min(1.0, 2 * sum(math.comb(n, i) for i in range(k + 1)) / 2 ** n)


def describe(values):
    fl = [float(v) for v in values]
    q = statistics.quantiles(fl, n=4, method="inclusive") if len(fl) > 1 else [fl[0]] * 3
    ahead = sum(1 for v in values if v > 0)
    behind = sum(1 for v in values if v < 0)
    return {"ahead": ahead, "level": len(values) - ahead - behind, "behind": behind,
            "median": statistics.median(fl), "q1": q[0], "q3": q[2], "mean": statistics.fmean(fl),
            "sign_test_p": "%.6g" % sign_p(ahead, behind)}


def redraw(by_species, seed):
    """The species interval: median, mean and share at most 0.02 ahead, or None under 6 species."""
    sp = sorted(by_species)
    if len(sp) < FEWEST_SPECIES:
        return None
    flt = {s: [float(x) for x in by_species[s]] for s in sp}
    near = {s: [x <= WITHIN for x in by_species[s]] for s in sp}
    rng = random.Random(seed)
    meds, means, shares = [], [], []
    for _ in range(RESAMPLES):
        xs, ws = [], []
        for _ in range(len(sp)):
            s = sp[rng.randrange(len(sp))]
            xs.extend(flt[s])
            ws.extend(near[s])
        meds.append(statistics.median(xs))
        means.append(statistics.fmean(xs))
        shares.append(sum(ws) / len(ws))
    ends = []
    for a in (meds, means, shares):
        a.sort()
        ends += [a[LOW_RANK], a[RESAMPLES - 1 - LOW_RANK]]
    return ends


class Cohort:
    """The row files read as datasets, callers and per-dataset differences."""

    def __init__(self, tables):
        self.datasets = tables["datasets.tsv"]
        self.by_id = {d["dataset"]: d for d in self.datasets}
        self.callers = {(r["dataset"], r["caller"]): r for r in tables["callers.tsv"]}
        self.diff = {(r["dataset"], r["contrast"], r["measure"], r["first"]): r for r in tables["differences.tsv"]}
        self.sets = {name: [d["dataset"] for d in self.datasets if members(name, d)] for name in SET_NAMES}

    def rates(self, ds, caller):
        return rates(self.callers[(ds, caller)])

    def pick(self, ds, split, pick):
        """The learner a pick names on one dataset."""
        if pick == "served_best":
            return self.by_id[ds]["served_learner"]
        if pick == "highest_of_five":
            best, top = None, None
            for l in LEARNERS:
                ba = self.rates(ds, "learner:%s:%s" % (l, split))["ba"]
                if top is None or ba > top:
                    best, top = l, ba
            return best
        return pick

    def units(self, name, con):
        """The datasets of a set that a contrast reads: a lineage contrast needs a lineage split."""
        sel = self.sets[name]
        if con.startswith("learner_minus_program:lineage:"):
            sel = [ds for ds in sel if self.by_id[ds]["lineage_split"] == "yes"]
        return sel

    def value(self, ds, con, m):
        program = self.rates(ds, PROGRAM)[m]
        if con == PR:
            return program - self.rates(ds, RULE)[m]
        _, split, pick = con.split(":")
        return self.rates(ds, "learner:%s:%s" % (self.pick(ds, split, pick), split))[m] - program

    def own_row(self, ds, con, m):
        if con == PR:
            return self.diff[(ds, PR, m, PROGRAM)]
        _, split, pick = con.split(":")
        return self.diff[(ds, "learner_minus_program:" + split, m,
                          "learner:%s:%s" % (self.pick(ds, split, pick), split))]

    def by_species(self, name, con, m):
        out = {}
        for ds in self.units(name, con):
            out.setdefault(self.by_id[ds]["species"], []).append(self.value(ds, con, m))
        return out


# ---------------------------------------------------------------- the row files' rules

def validate(tables):
    for name in ROW_FILES:
        for r in tables[name]:
            for cell in r.values():
                if GENOME_ID.search(cell) or BASES.search(cell):
                    raise Stop("%s: a cell has the shape of a genome id or of a base run" % name)
    ds_rows = tables["datasets.tsv"]
    ids = [d["dataset"] for d in ds_rows]
    if not ids or ids != sorted(set(ids)):
        raise Stop("datasets.tsv: datasets not unique and in order")
    for d in ds_rows:
        if not re.fullmatch(r"[a-z]+_[a-z]+__[a-z0-9_]+", d["dataset"]) or \
                not d["dataset"].startswith(d["species"] + "__"):
            raise Stop("datasets.tsv: %s is not species__antibiotic" % d["dataset"])
        for col in ("organism_option", "lineage_split", "read_before_the_plan"):
            if d[col] not in ("yes", "no"):
                raise Stop("datasets.tsv: %s %s is %r" % (d["dataset"], col, d[col]))
        if d["served_learner"] not in LEARNERS:
            raise Stop("datasets.tsv: %s serves %r" % (d["dataset"], d["served_learner"]))
    want = []
    for d in ds_rows:
        callers = ["program:" + r for r in READINGS] + [RULE]
        splits = SPLITS if d["lineage_split"] == "yes" else ["random"]
        callers += ["learner:%s:%s" % (l, s) for s in splits for l in LEARNERS]
        want += [(d["dataset"], c) for c in sorted(callers)]
    got = [(r["dataset"], r["caller"]) for r in tables["callers.tsv"]]
    if got != want:
        raise Stop("callers.tsv: rows are not every caller of every dataset, in order")
    for r in tables["callers.tsv"]:
        where = "callers.tsv %s %s" % (r["dataset"], r["caller"])
        n = {k: int(r[k]) for k in ("genomes", "scored", "no_call_resistant", "no_call_susceptible", "resistant",
                                     "susceptible", "tp", "fn", "tn", "fp")}
        if n["tp"] + n["fn"] != n["resistant"] or n["tn"] + n["fp"] != n["susceptible"] or \
                n["resistant"] + n["susceptible"] != n["scored"] or \
                n["scored"] + n["no_call_resistant"] + n["no_call_susceptible"] != n["genomes"]:
            raise Stop("%s: the counts do not add up" % where)
        rr = rates(r)
        for col, key in (("sensitivity", "sens"), ("specificity", "spec"), ("balanced_accuracy", "ba")):
            if r[col] != f6(rr[key]):
                raise Stop("%s: %s does not follow from the counts" % (where, col))
    cohort = Cohort(tables)
    seen = set()
    for r in tables["differences.tsv"]:
        where = "differences.tsv %s %s %s %s" % (r["dataset"], r["contrast"], r["measure"], r["first"])
        key = (r["dataset"], r["contrast"], r["measure"], r["first"])
        if key in seen or r["reading"] != "tool_default" or r["dataset"] not in cohort.by_id:
            raise Stop("%s: repeated, of another reading or of no dataset" % where)
        seen.add(key)
        first, second = cohort.callers.get((r["dataset"], r["first"])), cohort.callers.get((r["dataset"], r["second"]))
        if first is None or second is None:
            raise Stop("%s: names a caller callers.tsv does not hold" % where)
        exact = rates(first)[r["measure"]] - rates(second)[r["measure"]]
        if abs(float(r["first_minus_second"]) - float(exact)) > 5e-7:
            raise Stop("%s: first_minus_second does not follow from callers.tsv" % where)
        if int(r["first_only_right"]) + int(r["second_only_right"]) + int(r["calls_agree"]) != int(r["scored"]) or \
                int(r["scored"]) != int(first["scored"]):
            raise Stop("%s: the call counts do not add up" % where)
        if not float(r["lineage_lo"]) <= float(r["lineage_hi"]):
            raise Stop("%s: the lineage interval is reversed" % where)
    for name in SET_NAMES:
        if not cohort.sets[name]:
            raise Stop("set %s holds no dataset" % name)
    keys = [(r["set"], r["contrast"], r["measure"]) for r in tables["intervals.tsv"]]
    if keys != targets():
        raise Stop("intervals.tsv: rows are not the summary's rows, in order")
    for r in tables["intervals.tsv"]:
        name, con, m = r["set"], r["contrast"], r["measure"]
        few = len(cohort.by_species(name, con, m)) < FEWEST_SPECIES
        cols = ["median_species_lo", "median_species_hi", "mean_species_lo", "mean_species_hi"]
        share = ["share_species_lo", "share_species_hi"]
        if con == PR and any(r[c] for c in share):
            raise Stop("intervals.tsv: %s %s %s gives a share for the program against the rule" % (name, con, m))
        cols += [] if con == PR else share
        for lo, hi in zip(cols[::2], cols[1::2]):
            if few != (r[lo] == "" and r[hi] == "") or (not few and float(r[lo]) > float(r[hi])):
                raise Stop("intervals.tsv: %s %s %s %s" % (name, con, m, lo))
    utc = [r["utc"] for r in tables["order.tsv"]]
    if utc != sorted(utc) or len(utc) != len(ORDER) + 1:
        raise Stop("order.tsv: not %d events in clock order" % (len(ORDER) + 1))
    settings = [(r["key"], r["value"]) for r in tables["settings.tsv"]]
    if settings != SETTINGS:
        raise Stop("settings.tsv: the settings differ from the script's")


# ---------------------------------------------------------------- the derived files

def summarize(cohort, intervals):
    rows = []
    for name, con, m in targets():
        units = cohort.units(name, con)
        vals = [cohort.value(ds, con, m) for ds in units]
        dsc = describe(vals)
        row = {"set": name, "contrast": con, "measure": m, "role": role(name, con, m), "datasets": len(units),
               "species": len({cohort.by_id[ds]["species"] for ds in units})}
        row.update({k: dsc[k] for k in ("ahead", "level", "behind", "sign_test_p")})
        row.update({k: f6(dsc[k]) for k in ("median", "q1", "q3", "mean")})
        iv = intervals[(name, con, m)]
        for k in ("median_species_lo", "median_species_hi", "mean_species_lo", "mean_species_hi",
                  "share_species_lo", "share_species_hi"):
            row[k] = iv[k]
        if con == PR or m == "ba":
            own = [cohort.own_row(ds, con, m) for ds in units]
            above = sum(1 for r in own if float(r["lineage_lo"]) > 0)
            below = sum(1 for r in own if float(r["lineage_hi"]) < 0)
            row.update({"own_interval_above_zero": above, "own_interval_spans_zero": len(own) - above - below,
                        "own_interval_below_zero": below})
        else:
            row.update({"own_interval_above_zero": "", "own_interval_spans_zero": "", "own_interval_below_zero": ""})
        if con == PR:
            row.update({"learner_at_most_0.02_ahead": "", "share_learner_at_most_0.02_ahead": "",
                        "within_0.02_either_way": ""})
        else:
            near = sum(1 for v in vals if v <= WITHIN)
            row.update({"learner_at_most_0.02_ahead": near, "share_learner_at_most_0.02_ahead": f6(near / len(vals)),
                        "within_0.02_either_way": sum(1 for v in vals if abs(v) <= WITHIN)})
        rows.append(row)
    return rows


def figures(cohort):
    """The single figures the page reads that no summary row holds."""
    s235 = cohort.sets["primary-new-235"]
    pr = {ds: cohort.value(ds, PR, "ba") for ds in s235}
    large = [ds for ds in s235 if pr[ds] >= LARGE]
    level = [ds for ds in s235 if pr[ds] == 0]
    same = sum(1 for ds in level if cohort.own_row(ds, PR, "ba")["calls_agree"] == cohort.own_row(ds, PR, "ba")["scored"])
    served = contrast_name("lineage", "served_best")
    top = sum(1 for ds in s235
              if cohort.value(ds, served, "ba") == cohort.value(ds, contrast_name("lineage", "highest_of_five"), "ba"))
    m28 = cohort.sets["rule-M-quinolone-28"]
    edge = []
    for ds in m28:
        program = cohort.rates(ds, PROGRAM)["ba"]
        learner = program + cohort.value(ds, served, "ba")
        if program == Fraction(1, 2) and learner < program - WITHIN:
            edge.append(learner)
    edge.sort(reverse=True)
    rows = [
        ("program_ahead_by_0.10_or_more", "primary-new-235", len(large), "exploratory"),
        ("of_them_quinolone", "primary-new-235",
         sum(1 for ds in large if cohort.by_id[ds]["drug_class"] == "quinolone"), "exploratory"),
        ("level_with_the_same_call_on_every_genome", "primary-new-235", same, "exploratory"),
        ("served_learner_highest_of_five_on_the_lineage_split", "primary-new-235", top, "exploratory"),
        ("served_learner_behind_by_more_than_0.02_where_the_program_is_at_0.5", "rule-M-quinolone-28", len(edge),
         "exploratory"),
    ]
    rows += [("served_learner_balanced_accuracy_%d_of_them" % (i + 1), "rule-M-quinolone-28", f6(v), "exploratory")
             for i, v in enumerate(edge)]
    return [dict(zip(COLUMNS["figures.tsv"], r)) for r in rows]


def derive(tables):
    cohort = Cohort(tables)
    intervals = {(r["set"], r["contrast"], r["measure"]): r for r in tables["intervals.tsv"]}
    sets = []
    for name, rl, definition in SETS:
        sel = cohort.sets[name]
        sets.append({"set": name, "role": rl, "datasets": len(sel),
                     "species": len({cohort.by_id[ds]["species"] for ds in sel}), "definition": definition})
    return {
        "sets.tsv": tsv("sets.tsv", sets),
        "summary.tsv": tsv("summary.tsv", summarize(cohort, intervals)),
        "figures.tsv": tsv("figures.tsv", figures(cohort)),
    }


README_ROWS = {
    "README.md": "this file",
    "SHA256SUMS": "the sha256 of every other file",
    "datasets.tsv": "row file: one row per dataset with a published rule",
    "callers.tsv": "row file: counts and measures per dataset and caller, on the genomes and labels of the release",
    "differences.tsv": "row file: paired differences per dataset of the primary reading, with each one's own lineage interval",
    "intervals.tsv": "row file: the species interval of each summary row, as the source modes redraw it",
    "order.tsv": "row file: the clock times of the plan, its three addenda and the final scoring",
    "settings.tsv": "row file: the program, database, release and interval settings",
    "sets.tsv": "derived: each set, its role and size",
    "summary.tsv": "derived: per set, contrast and measure, over the set's datasets",
    "figures.tsv": "derived: single figures the page reads that no summary row holds",
}
COLUMN_NOTES = [
    ("datasets.tsv", "`organism_option` is whether the program has an organism option for the species, "
                     "`lineage_split` whether the release holds a lineage split for the dataset, "
                     "`read_before_the_plan` whether its calls were read before the plan was hashed, and "
                     "`served_learner` the learner the release serves, picked on the random split."),
    ("callers.tsv", "`caller` is `program:<reading>` for the AMRFinderPlus program, `grm_known_gene_rule` for "
                    "GRM's known-gene rule and `learner:<learner>:<split>` for a learner. A resistant genome "
                    "called resistant is `tp`. `*_exact_*` are Clopper-Pearson ends, which treat genomes as "
                    "independent, and `*_lineage_*` ends from resampling the report's lineages."),
    ("differences.tsv", "`reading` is the program's reading; the record keeps the primary one. A dataset's "
                        "difference is `first` minus `second`. `first_only_right` counts genomes only the first "
                        "caller calls right, and `mcnemar_p` is the exact two-sided p from the two such counts."),
    ("summary.tsv", "`ahead`, `level` and `behind` count the datasets where the first named caller is ahead of, "
                    "level with or behind the second. `own_interval_*` count datasets by their own lineage "
                    "interval. `learner_at_most_0.02_ahead` counts datasets where the learner is at most 0.02 "
                    "ahead of the program, or behind it, as the plan's rule M defines within 0.02; "
                    "`within_0.02_either_way` counts those within 0.02 on both sides. Neither is given for the "
                    "program against the known-gene rule."),
]


def readme(files):
    lines = [
        "# AMRFinderPlus comparison record",
        "",
        "This directory is the record behind [`docs/amrfinderplus-comparison.md`](../../docs/amrfinderplus-comparison.md): "
        "the AMRFinderPlus program, GRM's known-gene rule and GRM's five learners scored on the same genomes and "
        "recorded labels of public BV-BRC datasets. `scripts/amrfinderplus-comparison/record.py` writes every "
        "file of it, and its `--check` mode fails on any byte it would not write.",
        "",
        "This is a retrospective comparison inside one public database, BV-BRC. It is not a prospective test "
        "and not external validation. Neither the AMRFinderPlus program nor GRM's known-gene rule reads a label, "
        "but the catalog both read is newer than most of these genomes and may have been curated from the same "
        "public isolates, which favours both against the learners.",
        "",
        "## Files",
        "",
        "| File | Rows | Contents |",
        "|---|---|---|",
    ]
    for name in sorted(files, key=lambda x: (x != "README.md", x != "SHA256SUMS", x)):
        rows = "" if not name.endswith(".tsv") else str(files[name].count("\n") - 1)
        lines.append("| `%s` | %s | %s |" % (name, rows, README_ROWS[name]))
    lines += [
        "",
        "Every file is plain text, and each TSV has a header row. A file holds counts and measures per dataset "
        "and caller, per set and contrast or per learner: no genome id, no sequence, no rule and no row per genome.",
        "",
        "## Columns",
        "",
    ]
    for name in COLUMNS:
        lines.append("- `%s`: %s." % (name, ", ".join("`%s`" % c for c in COLUMNS[name])))
    lines += [""] + ["- `%s`: %s" % (name, note) for name, note in COLUMN_NOTES] + [
        "",
        "## Rules",
        "",
        "- **Three callers.** The AMRFinderPlus program is read under each report's rule. `tool_default` is the "
        "primary reading; `tool_default_intact` counts no row whose method the program marks as partial or as an "
        "internal stop, `genes_core` only gene rows of core scope, and `genes_all_scopes` gene rows of core or plus scope. GRM's known-gene "
        "rule is the release's own call from the catalog, and a learner calls resistant at a score of 0.5 or more.",
        "- **Measures.** Sensitivity and specificity come from each caller's counts, with 0 for an empty class, "
        "and balanced accuracy is their mean. Differences are exact, so level means a difference of exactly 0.",
        "- **Over datasets.** The median, quartiles and mean of the per-dataset differences, and the exact "
        "two-sided sign test on ahead against behind. The species interval resamples the set's species with "
        "replacement and takes every dataset of a drawn species. The plan fixed 4,000 resamples and the seed "
        "string, and named neither the generator nor the draw: the source modes use `random.Random` seeded with "
        "the string, `randrange` over the species sorted by name and `statistics.median`, which give the final "
        "scoring's ends. `intervals.tsv` holds the ends; the record does not hold the seed string.",
        "- **Roles.** `primary` rows are the plan's headline: the program against the known-gene rule on "
        "primary-new-235 in each measure, and each learner and the served learner against the program in "
        "balanced accuracy on its lineage split. `secondary` rows are those the plan listed beside it. "
        "`exploratory` rows and figures were chosen after the results were read.",
        "",
        "## Commands",
        "",
        "```sh",
        "python3 -I scripts/amrfinderplus-comparison/record.py --check reference/amrfinderplus-comparison",
        "python3 -I scripts/amrfinderplus-comparison/record.py --source \"$FINAL\" --check reference/amrfinderplus-comparison",
        "```",
        "",
        "The first rebuilds the derived files, this README and `SHA256SUMS` from the row files. The second "
        "rebuilds every file from the final scoring, which is not public.",
        "",
    ]
    return "\n".join(lines)


def assemble(row_texts):
    """Every file of the record from the row files' texts."""
    tables = {}
    for name in ROW_FILES:
        head, rows = parse(row_texts[name])
        if head != COLUMNS[name]:
            raise Stop("%s: columns %s" % (name, head))
        tables[name] = rows
    validate(tables)
    files = dict(row_texts)
    files.update(derive(tables))
    files["README.md"] = ""
    files["SHA256SUMS"] = ""
    files["README.md"] = readme(files)
    del files["SHA256SUMS"]
    sums = "".join("%s  %s\n" % (hashlib.sha256(files[k].encode("utf-8")).hexdigest(), k) for k in sorted(files))
    files["SHA256SUMS"] = sums
    return files


# ---------------------------------------------------------------- the sources

def read(path):
    with open(path, encoding="utf-8") as fh:
        return fh.read()


def read_gz(path):
    with gzip.open(path, "rt", encoding="utf-8", newline="") as fh:
        return fh.read()


def from_sources(source):
    final = pathlib.Path(source)
    results = final.parent
    top = results.parent
    texts = {}

    # Datasets, from the per-dataset headline, the readings table and the callers.
    _, head = parse(read(final / "per-dataset-headline.tsv"))
    _, by_dataset = parse(read(final / "by-dataset.tsv"))
    point = {r["dataset_id"]: r["point_screened"] for r in by_dataset if r["variant"] == "tool_default"}
    _, callers = parse(read_gz(final / "by-dataset-callers.tsv.gz"))
    lineage = {r["dataset_id"] for r in callers if r["caller"] == "learner:scm:lineage"}
    datasets = []
    for h in sorted(head, key=lambda r: r["dataset_id"]):
        ds = h["dataset_id"]
        row = {"dataset": ds, "species": h["species_slug"], "drug_class": h["drug_class"],
               "organism_option": point[ds], "lineage_split": "yes" if ds in lineage else "no",
               "read_before_the_plan": h["read_before_this_plan"], "served_learner": h["served_best"]}
        if h["primary_set"] != ("yes" if row["organism_option"] == "yes" and row["lineage_split"] == "yes" else "no"):
            raise Stop("per-dataset-headline.tsv: %s primary_set is not organism option and lineage split" % ds)
        datasets.append(row)
    texts["datasets.tsv"] = tsv("datasets.tsv", datasets)
    rows = [dict(r, dataset=r["dataset_id"]) for r in callers]
    rows.sort(key=lambda r: (r["dataset"], r["caller"]))
    texts["callers.tsv"] = tsv("callers.tsv", rows)

    # Differences of the primary reading, with variant renamed reading.
    _, diffs = parse(read_gz(final / "differences.tsv.gz"))
    rows = [dict(r, dataset=r["dataset_id"], reading=r["variant"]) for r in diffs if r["variant"] == "tool_default"]
    texts["differences.tsv"] = tsv("differences.tsv", rows)

    # The sets, held to the plan's lists.
    tables = {name: parse(texts[name])[1] for name in ("datasets.tsv", "callers.tsv", "differences.tsv")}
    cohort = Cohort(tables)
    for name in sorted(PLAN_LISTS):
        listed = sorted(x.strip() for x in read(results / "plan-lists" / (name + ".txt")).split("\n") if x.strip())
        if listed != sorted(cohort.sets[name]):
            raise Stop("set %s differs from the plan's list" % name)

    # The species intervals, redrawn from the plan's seed string, and every summary cell held to the
    # final scoring's own summary.
    plan = read(top / "PLAN.md")
    seeds = [x for x in re.findall(r"seed string `([^`]+)`", plan) if "<set>" in x]
    if len(seeds) != 1 or not all(x in seeds[0] for x in ("<variant>", "<contrast>", "<measure>")):
        raise Stop("PLAN.md: no single seed string")
    _, cohort_rows = parse(read_gz(final / "cohort-summary.tsv.gz"))
    cs = {(r["set"], r["contrast"], r["measure"]): r for r in cohort_rows if r["variant"] == "tool_default"}
    intervals = []
    for name, con, m in targets():
        seed = seeds[0].replace("<set>", name).replace("<variant>", "tool_default")
        seed = seed.replace("<contrast>", con).replace("<measure>", m)
        ends = redraw(cohort.by_species(name, con, m), seed)
        row = {"set": name, "contrast": con, "measure": m}
        keys = ["median_species_lo", "median_species_hi", "mean_species_lo", "mean_species_hi", "share_species_lo",
                "share_species_hi"]
        for i, k in enumerate(keys):
            row[k] = "" if ends is None or (con == PR and k.startswith("share")) else f6(ends[i])
        theirs = cs.get((name, con, m))
        if theirs is None:
            raise Stop("cohort-summary.tsv.gz: no row for %s %s %s" % (name, con, m))
        for k in keys:
            t = theirs[k.replace("share_species", "share_within_species")]
            if row[k] != (t if t == "" else f6(float(t))) and not (con == PR and k.startswith("share")):
                raise Stop("cohort-summary.tsv.gz: %s %s %s %s is %s, redrawn %s" % (name, con, m, k, t, row[k]))
        intervals.append(row)
    texts["intervals.tsv"] = tsv("intervals.tsv", intervals)

    # Order and settings.
    order = []
    for event, name, sums, when in ORDER:
        listed = [ln.split() for ln in read(top / sums).split("\n") if ln.strip()]
        digest = hashlib.sha256((top / name).read_bytes()).hexdigest()
        if [digest, name] not in listed:
            raise Stop("%s does not give the sha256 of %s" % (sums, name))
        utc = read(top / when).strip()
        if not re.fullmatch(r"2026-11-1[4-6]T\d\d:\d\d:\d\dZ", utc):
            raise Stop("%s holds no clock time" % when)
        order.append({"utc": utc, "event": event, "file": name, "sha256": digest})
    first = read(top / FINAL_LOG).split("\n", 1)[0]
    found = re.match(r"(2026-11-\d\dT\d\d:\d\d:\d\dZ) scope final: (\d+) datasets", first)
    if not found or found.group(2) != str(len(datasets)):
        raise Stop("%s: no start of the final scoring" % FINAL_LOG)
    order.append({"utc": found.group(1), "event": "final scoring started", "file": "", "sha256": ""})
    texts["order.tsv"] = tsv("order.tsv", order)
    run_file = json.loads(read(final / "run-manifest.json"))
    adapter = json.loads(read(final / "adapter-summary.json"))
    given = dict(SETTINGS)
    for key, value in (("program", run_file.get("tool")), ("program_version", run_file.get("software_version")),
                       ("database_version", run_file.get("database_version")),
                       ("datasets_scored", str(adapter.get("datasets_scored"))),
                       ("genomes_read", str(adapter.get("genomes_read")))):
        if given[key] != value:
            raise Stop("settings: %s is %r in the final scoring" % (key, value))
    for key, text in SETTINGS_IN_PLAN.items():
        if text not in plan:
            raise Stop("PLAN.md does not give %s as %r" % (key, text))
    texts["settings.tsv"] = tsv("settings.tsv", [{"key": a, "value": b} for a, b in SETTINGS])

    # Every summary cell of the record, held to the final scoring's summary.
    files = assemble(texts)
    _, summary = parse(files["summary.tsv"])
    same = {"ahead": "first_ahead", "behind": "second_ahead", "level": "equal", "datasets": "units",
            "species": "species", "median": "median", "q1": "q1", "q3": "q3", "mean": "mean",
            "own_interval_above_zero": "lineage_interval_above_zero",
            "own_interval_spans_zero": "lineage_interval_spans_zero",
            "own_interval_below_zero": "lineage_interval_below_zero",
            "learner_at_most_0.02_ahead": "within_0.02", "share_learner_at_most_0.02_ahead": "share_within_0.02"}
    for r in summary:
        theirs = cs[(r["set"], r["contrast"], r["measure"])]
        for mine, col in same.items():
            if r[mine] == "" and r["contrast"] == PR and mine.endswith("ahead"):
                continue
            if r[mine] != theirs[col]:
                raise Stop("cohort-summary.tsv.gz: %s %s %s %s is %s, here %s"
                           % (r["set"], r["contrast"], r["measure"], col, theirs[col], r[mine]))
        a, b = float(r["sign_test_p"]), float(theirs["sign_test_p"])
        if abs(a - b) > 1e-5 * max(abs(b), 1e-300):
            raise Stop("cohort-summary.tsv.gz: %s %s %s sign_test_p" % (r["set"], r["contrast"], r["measure"]))
    return texts


def from_record(directory):
    d = pathlib.Path(directory)
    texts = {}
    for name in ROW_FILES:
        path = d / name
        if not path.is_file():
            raise Stop("%s is missing" % name)
        texts[name] = read(path)
    return texts


# ---------------------------------------------------------------- main

def check(files, directory):
    d = pathlib.Path(directory)
    on_disk = {p.relative_to(d).as_posix() for p in d.rglob("*") if p.is_file()}
    bad = sorted(on_disk ^ set(files))
    for name in sorted(files):
        path = d / name
        if path.is_file() and path.read_bytes() != files[name].encode("utf-8"):
            bad.append(name)
    return sorted(set(bad))


def earlier_record(out, files):
    """Whether a directory holds nothing, or an earlier record of this script and nothing else."""
    found = list(out.iterdir())
    if not found:
        return True
    if any(not p.is_file() for p in found) or {p.name for p in found} - set(files):
        return False
    sums, readme_file = out / "SHA256SUMS", out / "README.md"
    if not (sums.is_file() and readme_file.is_file()):
        return False
    listed = {ln.split("  ", 1)[-1] for ln in sums.read_text(encoding="utf-8", errors="replace").splitlines()}
    first = readme_file.read_text(encoding="utf-8", errors="replace").split("\n", 1)[0]
    return "README.md" in listed and first == files["README.md"].split("\n", 1)[0]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--source", help="the final scoring's directory")
    group = ap.add_mutually_exclusive_group(required=True)
    group.add_argument("--out", help="write the record here")
    group.add_argument("--check", help="fail unless this directory is the record")
    args = ap.parse_args(argv)
    if args.out and not args.source:
        ap.error("--out needs --source")
    try:
        texts = from_sources(args.source) if args.source else from_record(args.check)
        files = assemble(texts)
    except (Stop, OSError, KeyError, ValueError, IndexError, ZeroDivisionError, json.JSONDecodeError) as exc:
        print("record.py: stopped: %s" % exc, file=sys.stderr)
        return 1
    if args.out:
        out = pathlib.Path(args.out)
        if out.exists() and not (out.is_dir() and earlier_record(out, files)):
            print("record.py: stopped: %s holds files that are not an earlier record of this script" % out,
                  file=sys.stderr)
            return 1
        out.mkdir(parents=True, exist_ok=True)
        for name, text in files.items():
            (out / name).write_bytes(text.encode("utf-8"))
        print("record.py: wrote %d files to %s" % (len(files), out))
        return 0
    bad = check(files, args.check)
    for name in bad:
        print("record.py: differs: %s" % name, file=sys.stderr)
    if bad:
        return 1
    print("record.py: %d files as written" % len(files))
    return 0


if __name__ == "__main__":
    sys.exit(main())
