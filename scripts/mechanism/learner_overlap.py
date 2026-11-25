"""learner_overlap.py: do the five learners of the revision 14 release miss
the same resistant genomes on the lineage split? (docs/learner-overlap-result.md)

    python3 -I learner_overlap.py plan
    python3 -I learner_overlap.py read REPORTS_DIR --out READING_DIR [--expect 274]
    python3 -I learner_overlap.py write READING_DIR --out RECORD_DIR
    python3 -I learner_overlap.py print RECORD_DIR

`plan` prints the fixed rules as JSON. `read` reads every report of
REPORTS_DIR, runs the page's checks 1 to 5 and 7 to 10 and writes one file,
reading.json: per report, counts only. `write` turns a reading into
datasets.tsv, pooled.tsv, gates.txt and SHA256SUMS, after check 6 finds no
genome id in any file of the reading or the record. `print` prints the page's
tables from a record.

Every decimal of the record and of the page is worked from the counts as an
exact fraction and rounded once: the record's to 3 places (the lead's ruling of
25 November 2026, after a review, in place of the plan's 12), the page's to the
places its tables show. No pooled figure and no page figure is worked from a
rounded cell. A report's ratio is not a column of the record (the same ruling):
it is given as the counts it comes from, `miss_all5` over the yardstick worked
from `resistant` and the five learners' miss columns.

Nothing is fitted, trained or scored again: every count comes from the scores
and labels the reports store. Every key named `id` is dropped as a report is
parsed, so no genome id is held, printed or written. A failed stop check
writes nothing and exits 1.

Python 3.12, standard library only; run with -I.
"""

import argparse
import collections
import fractions
import hashlib
import json
import math
import os
import re
import shutil
import statistics
import sys
import tempfile

CUTOFF = 0.5                 # the release's own: a score of 0.5 or more is a resistant call
SCHEME = "lineage"           # the only split read
SEPARATE = ("catalog",)      # GRM's reader of the AMRFinderPlus catalogue: not a learner
FIVE = ("scm", "cart", "rf", "gbt", "lasso")
EXPECT_FILES = 274
MIN_RESISTANT = 10           # a report enters the pooled rows with 10 or more resistant genomes
MIN_LEARNERS = 2             # a report with fewer learners is skipped and counted
STACK = 10                   # a median ratio of 10 or more: the misses stack
KS = range(6)
# The shape of a BV-BRC genome id (a taxon number, a point, a genome number), standing alone. A computed decimal of
# the record has 3 places and may have this shape, so check 6 scans every cell of a record table but those of its
# computed decimal columns, which the write fills from exact fractions (ruled 25 Nov 2026).
ID_SHAPE = re.compile(r"(?<![\w.])\d{2,7}\.\d{1,9}(?![\w.])")
PLACES = 3                   # every decimal of the record: worked exactly, rounded once to 3 places
BARRED_KEY = re.compile(r'"(?:id|samples)"\s*:')
READING = "reading.json"
RECORD_FILES = ("datasets.tsv", "pooled.tsv", "gates.txt")
DATASET_COLUMNS = (["species", "antibiotic", "genomes", "resistant", "susceptible", "learners", "learner_count",
                    "catalog_rule"]
                   + [f"miss_{k}" for k in KS] + [f"wrong_{k}" for k in KS] + [f"miss_{n}" for n in FIVE]
                   + ["miss_all5", "miss_all5_no_catalog_call", "miss_all5_in_null_fold", "null_fold_pairs",
                      "expected_all5"])
SUMMED = (["genomes", "resistant", "susceptible"] + [f"miss_{k}" for k in KS] + [f"wrong_{k}" for k in KS]
          + [f"miss_{n}" for n in FIVE] + ["miss_all5", "miss_all5_no_catalog_call", "miss_all5_in_null_fold",
                                           "null_fold_pairs"])
POOLED_COLUMNS = (["scope", "species", "datasets"] + SUMMED
                  + ["datasets_no_catalog_rule", "miss_all5_no_catalog_rule", "expected_all5", "datasets_with_ratio",
                     "datasets_expected_zero", "ratio_median", "ratio_q1", "ratio_q3", "datasets_ratio_over_1",
                     "share_all5_median", "share_some_median"])
DECIMAL_COLUMNS = ("expected_all5", "ratio_median", "ratio_q1", "ratio_q3", "share_all5_median",
                   "share_some_median")


class Stop(Exception):
    """A stop check failed: nothing is written."""


# ---- reading one report -------------------------------------------------------------------------------------------

def drop_ids(pairs):
    """Every key named `id` goes as the file is parsed: no genome id is held."""
    return {k: v for k, v in pairs if k != "id"}


def load_report(path):
    with open(path, encoding="utf-8") as f:
        return json.load(f, object_pairs_hook=drop_ids)


def is_call(score):
    return score >= CUTOFF


def lineage_splits(report):
    return [s for s in report.get("splits") or [] if s.get("scheme") == SCHEME]


def learner_entries(split):
    return [lr for lr in split.get("learners") or [] if lr.get("learner") not in SEPARATE]


def separate_entries(split):
    return {lr.get("learner"): lr for lr in split.get("learners") or [] if lr.get("learner") in SEPARATE}


def good_score(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x) and 0 <= x <= 1


def scores_fit(scores, n):
    return isinstance(scores, list) and len(scores) == n and all(good_score(x) for x in scores)


def counts_agree(counts, labels):
    return counts.get("resistant") == sum(labels) and counts.get("samples") == len(labels)


def metrics_agree(labels, lr):
    stored = lr.get("metrics") or {}
    return recount(labels, lr["scores"]) == {k: stored.get(k) for k in ("tp", "fn", "fp", "tn")}


def catalog_agrees(lr, catalog):
    return [is_call(x) for x in lr["scores"]] == catalog


def uncalled(catalog_call):
    """The catalogue reader leaves the genome without a resistant call: false, or null where the report holds no
    catalogue rule for its antibiotic."""
    return not catalog_call


def catalog_shape(catalog):
    """'rule' when every call is true or false, 'no rule' when every call is null, None for a mix."""
    if all(isinstance(v, bool) for v in catalog):
        return "rule"
    if all(v is None for v in catalog):
        return "no rule"
    return None


def recount(labels, scores):
    c = {"tp": 0, "fn": 0, "fp": 0, "tn": 0}
    for t, x in zip(labels, scores):
        c[("tp" if is_call(x) else "fn") if t else ("fp" if is_call(x) else "tn")] += 1
    return c


def null_folds(lr):
    return {f for f, c in enumerate(lr.get("fold_complexity") or []) if c == 0}


def null_genomes(labels, folds, nulls):
    out = {"resistant": 0, "susceptible": 0}
    for t, f in zip(labels, folds):
        if f in nulls:
            out["resistant" if t else "susceptible"] += 1
    return out


def misses_of(labels, scores):
    return sum(1 for t, x in zip(labels, scores) if t and not is_call(x))


def tally(labels, catalog, folds, learners):
    """Counts of one report: `learners` maps a name to its scores and null folds."""
    names = list(learners)
    n_learners = len(names)
    miss_k, wrong_k = [0] * (n_learners + 1), [0] * (n_learners + 1)
    nulls = set().union(*(learners[n]["nulls"] for n in names)) if names else set()
    miss_all = miss_all_no_cat = miss_all_null = 0
    for i, t in enumerate(labels):
        called = sum(1 for n in names if is_call(learners[n]["scores"][i]))
        if t:
            k = n_learners - called
            miss_k[k] += 1
            if k == n_learners:
                miss_all += 1
                miss_all_no_cat += 1 if uncalled(catalog[i]) else 0
                miss_all_null += 1 if folds[i] in nulls else 0
        else:
            wrong_k[called] += 1
    return {"learners": names, "miss_k": miss_k, "wrong_k": wrong_k,
            "miss_by_learner": {n: misses_of(labels, learners[n]["scores"]) for n in names},
            "miss_all": miss_all, "miss_all_no_catalog_call": miss_all_no_cat, "miss_all_in_null_fold": miss_all_null,
            "null_fold_pairs": sum(len(learners[n]["nulls"]) for n in names)}


def name_of(report):
    return f"{report.get('species')} with {report.get('antibiotic')}"


# ---- the read -----------------------------------------------------------------------------------------------------

def report_files(directory):
    return sorted(e.name for e in os.scandir(directory) if e.is_file() and e.name.endswith(".json"))


def read_reports(directory, expect=EXPECT_FILES):
    """Read every report; return the reading, or raise Stop after printing every check's counts."""
    files = report_files(directory)
    fails = collections.defaultdict(list)       # stop check -> what failed (report names, never genomes)
    skipped = {"no_lineage": [], "under_two_learners": []}
    datasets, seen, sets, mismatch_rows, no_rule = [], set(), collections.Counter(), [], []
    parsed = 0
    for name in files:
        try:
            report = load_report(os.path.join(directory, name))
            if not isinstance(report, dict):
                raise ValueError("not an object")
        except (OSError, ValueError) as e:
            fails[1].append(f"{name}: {type(e).__name__}")
            continue
        parsed += 1
        who = name_of(report)
        key = (report.get("species"), report.get("antibiotic"))
        if key in seen:
            fails[8].append(who)
        seen.add(key)
        splits = lineage_splits(report)
        if not splits:
            skipped["no_lineage"].append(who)
            continue
        if len(splits) > 1:
            fails[5].append(f"{who}: {len(splits)} lineage splits")
            continue
        split = splits[0]
        entries = learner_entries(split)
        if len(entries) < MIN_LEARNERS:
            skipped["under_two_learners"].append(who)
            continue
        samples = report.get("samples") or []
        n = len(samples)
        labels = [s.get("resistant") for s in samples]
        catalog = [s.get("catalog_resistant") for s in samples]
        folds = split.get("test_fold")
        shape = catalog_shape(catalog)
        if not all(isinstance(v, bool) for v in labels) or shape is None:
            fails[3].append(f"{who}: a label is not true or false, or catalogue calls mix null with true or false")
            continue
        bad = [lr.get("learner") for lr in split.get("learners") or [] if not scores_fit(lr.get("scores"), n)]
        if bad or not (isinstance(folds, list) and len(folds) == n):
            fails[2].append(f"{who}: " + ", ".join(str(b) for b in bad) if bad else f"{who}: test_fold")
            continue
        counts = report.get("counts") or {}
        if not counts_agree(counts, labels):
            fails[3].append(who)
            continue
        learners = {}
        for lr in entries:
            if not metrics_agree(labels, lr):
                fails[7].append(f"{who}: {lr.get('learner')}")
            nulls = null_folds(lr)
            stored = lr.get("no_call") or {"resistant": 0, "susceptible": 0}
            if null_genomes(labels, folds, nulls) != {k: stored.get(k) for k in ("resistant", "susceptible")}:
                mismatch_rows.append(f"{who}: {lr.get('learner')}")
            learners[lr.get("learner")] = {"scores": lr["scores"], "nulls": nulls}
        for sname, lr in separate_entries(split).items():
            if sname == "catalog" and not catalog_agrees(lr, catalog):
                fails[9].append(who)
        sets[tuple(sorted(learners))] += 1
        if shape == "no rule":
            no_rule.append(who)
        row = {"species": report.get("species"), "antibiotic": report.get("antibiotic"), "genomes": n,
               "resistant": sum(labels), "susceptible": n - sum(labels), "catalog_rule": shape == "rule"}
        row.update(tally(labels, catalog, folds, learners))
        datasets.append(row)
    release = max(sets, key=lambda s: (sets[s], s)) if sets else ()
    differ = sorted(f"{d['species']} with {d['antibiotic']}: {','.join(d['learners'])}" for d in datasets
                    if tuple(sorted(d["learners"])) != release)
    n_skipped = sum(len(v) for v in skipped.values())
    lines = [
        f"check 1 parse: {len(files)} files found, {parsed} read, {expect} expected: "
        + ("pass" if parsed == len(files) == expect and not fails[1] else "STOP"),
        f"check 2 lengths: {len(fails[2])} reports with a score list or test_fold of another length or a bad score: "
        + ("pass" if not fails[2] else "STOP"),
        f"check 3 counts: {len(fails[3])} reports whose counts or labels differ from their samples; "
        f"{len(no_rule)} used reports with no catalogue rule (catalog_resistant null on every genome): "
        + ("pass" if not fails[3] else "STOP"),
        f"check 4 learner sets: {len(sets)} sets; the release's is {','.join(release)} on {sets.get(release, 0)} "
        f"reports; {len(differ)} reports differ",
        f"check 5 skipped: {len(skipped['no_lineage'])} with no lineage split, {len(skipped['under_two_learners'])} "
        f"with fewer than 2 learners; {n_skipped} skipped and {len(datasets)} used of {len(files)} found: "
        + ("pass" if n_skipped + len(datasets) == len(files) and not fails[5] and not fails[1] else "STOP"),
        f"check 7 stored metrics: {len(fails[7])} learner rows whose counts at the cutoff differ from metrics: "
        + ("pass" if not fails[7] else "STOP"),
        f"check 8 pairs: {len(fails[8])} species and antibiotic pairs on two reports: "
        + ("pass" if not fails[8] else "STOP"),
        f"check 9 catalogue: {len(fails[9])} reports whose catalog calls differ from catalog_resistant: "
        + ("pass" if not fails[9] else "STOP"),
        f"check 10 null folds: {len(mismatch_rows)} learner rows whose null-fold genomes differ from no_call",
    ]
    lines += [f"  differs: {d}" for d in differ]
    lines += [f"  skipped, no lineage split: {d}" for d in sorted(skipped["no_lineage"])]
    lines += [f"  skipped, fewer than 2 learners: {d}" for d in sorted(skipped["under_two_learners"])]
    lines += [f"  null folds differ: {d}" for d in sorted(mismatch_rows)]
    lines += [f"  no catalogue rule: {d}" for d in sorted(no_rule)]
    for check in sorted(fails):
        lines += [f"  check {check} failed: {d}" for d in sorted(fails[check])]
    stop = any("STOP" in line for line in lines[:9])
    if not stop and release != tuple(sorted(FIVE)):
        lines.append(f"  the release's learner set {','.join(release)} is not the plan's five")
        stop = True
    if stop:
        raise Stop("\n".join(lines))
    return {"cutoff": CUTOFF, "scheme": SCHEME, "five": list(FIVE), "found": len(files), "read": parsed,
            "used": len(datasets), "skipped": {k: len(v) for k, v in skipped.items()},
            "learner_sets": {",".join(s): c for s, c in sorted(sets.items())}, "differ": differ,
            "null_fold_mismatch": len(mismatch_rows), "no_catalog_rule": len(no_rule), "checks": lines,
            "datasets": sorted(datasets, key=lambda d: (d["species"], d["antibiotic"]))}


# ---- check 6 ------------------------------------------------------------------------------------------------------

def id_hits(paths):
    """Id-shaped tokens and `id` or `samples` keys in the files: a count, never the hit."""
    hits = 0
    for path in paths:
        with open(path, encoding="utf-8") as f:
            text = f.read()
        hits += len(ID_SHAPE.findall(text)) + len(BARRED_KEY.findall(text))
    return hits


def record_hits(path, decimal_columns):
    """Id-shaped tokens in every cell of a record table but those of its computed decimal columns, and `id` or
    `samples` keys anywhere in it: (hits, cells of the decimal columns left to the id list). Counts only."""
    with open(path, encoding="utf-8") as f:
        text = f.read()
    hits, skipped = len(BARRED_KEY.findall(text)), 0
    lines = text.split("\n")
    head = lines[0].split("\t")
    hits += sum(len(ID_SHAPE.findall(cell)) for cell in head)
    for line in lines[1:]:
        for name, cell in zip(head, line.split("\t")):
            if name in decimal_columns and cell != "":
                skipped += 1
            else:
                hits += len(ID_SHAPE.findall(cell))
        hits += sum(len(ID_SHAPE.findall(cell)) for cell in line.split("\t")[len(head):])
    return hits, skipped


def files_under(directory):
    return sorted(os.path.join(root, f) for root, _, names in os.walk(directory) for f in names)


def fresh_out(out):
    if os.path.lexists(out) and (not os.path.isdir(out) or os.listdir(out)):
        raise Stop(f"{out} exists and is not an empty directory")


def publish(tmp, out):
    if os.path.isdir(out):
        os.rmdir(out)
    os.replace(tmp, out)


# ---- the record ---------------------------------------------------------------------------------------------------

def rec(x):
    """A decimal of the record: an exact value rounded once, half to even, to PLACES places."""
    q = round(fractions.Fraction(x) * 10 ** PLACES)
    if q < 0:
        raise Stop("a decimal of the record is under 0")
    whole, part = divmod(q, 10 ** PLACES)
    return f"{whole}.{part:0{PLACES}d}"


def expected_all5(resistant, misses):
    """`resistant` times the product of the miss rates, as an exact fraction."""
    if resistant == 0:
        return fractions.Fraction(0)
    out = fractions.Fraction(resistant)
    for m in misses:
        out *= fractions.Fraction(m, resistant)
    return out


def is_five(d):
    return tuple(sorted(d["learners"])) == tuple(sorted(FIVE))


def dataset_row(d):
    n_learners = len(d["learners"])
    row = {k: d[k] for k in ("species", "antibiotic", "genomes", "resistant", "susceptible")}
    row["learners"] = ",".join(d["learners"])
    row["learner_count"] = n_learners
    row["catalog_rule"] = "yes" if d["catalog_rule"] else "no"
    if n_learners <= 5:
        for k in KS:
            row[f"miss_{k}"] = d["miss_k"][k] if k <= n_learners else ""
            row[f"wrong_{k}"] = d["wrong_k"][k] if k <= n_learners else ""
    for n in FIVE:
        row[f"miss_{n}"] = d["miss_by_learner"].get(n, "")
    if is_five(d):
        exp = expected_all5(d["resistant"], [d["miss_by_learner"][n] for n in FIVE])
        row.update({"miss_all5": d["miss_all"], "miss_all5_no_catalog_call": d["miss_all_no_catalog_call"],
                    "miss_all5_in_null_fold": d["miss_all_in_null_fold"], "null_fold_pairs": d["null_fold_pairs"],
                    "expected_all5": rec(exp)})
    return row


def ratio_of(r):
    """A report row's yardstick and ratio, worked again from its count columns as exact fractions (None: no ratio)."""
    exp = expected_all5(int(r["resistant"]), [int(r[f"miss_{n}"]) for n in FIVE])
    return exp, (int(r["miss_all5"]) / exp if exp > 0 else None)


def pooled_in(row):
    return row["learner_count"] == 5 and row["miss_all5"] != "" and row["resistant"] >= MIN_RESISTANT


def median_of(values):
    return statistics.median(values)


def quartiles(values):
    if len(values) == 1:
        return values[0], values[0]
    q = statistics.quantiles(values, n=4, method="inclusive")
    return q[0], q[2]


def ratio_values(rows):
    """The reports' ratios, exact, from their count columns: never from a rounded cell."""
    return [x for x in (ratio_of(r)[1] for r in rows) if x is not None]


def yardstick_sum(rows):
    return sum((ratio_of(r)[0] for r in rows), fractions.Fraction(0))


def exact_pool(rows):
    """A scope's decimals as exact fractions, worked from the count columns of its reports."""
    ratios = sorted(ratio_values(rows))
    v = {"expected_all5": yardstick_sum(rows), "ratios": ratios,
         "share_all5_median": median_of([fractions.Fraction(int(r["miss_all5"]), int(r["resistant"])) for r in rows]),
         "share_some_median": median_of([fractions.Fraction(sum(int(r[f"miss_{k}"]) for k in (1, 2, 3, 4)),
                                                            int(r["resistant"])) for r in rows])}
    if ratios:
        v["ratio_median"] = median_of(ratios)
        v["ratio_q1"], v["ratio_q3"] = quartiles(ratios)
    return v


def pool(rows, scope, species=""):
    out = {"scope": scope, "species": species, "datasets": len(rows)}
    for c in SUMMED:
        out[c] = sum(int(r[c]) for r in rows)
    out["datasets_no_catalog_rule"] = sum(1 for r in rows if r["catalog_rule"] == "no")
    out["miss_all5_no_catalog_rule"] = sum(int(r["miss_all5"]) for r in rows if r["catalog_rule"] == "no")
    v = exact_pool(rows)
    out["expected_all5"] = rec(v["expected_all5"])
    out["datasets_with_ratio"] = len(v["ratios"])
    out["datasets_expected_zero"] = len(rows) - len(v["ratios"])
    if v["ratios"]:
        out.update({"ratio_median": rec(v["ratio_median"]), "ratio_q1": rec(v["ratio_q1"]),
                    "ratio_q3": rec(v["ratio_q3"]), "datasets_ratio_over_1": sum(1 for x in v["ratios"] if x > 1)})
    else:
        out.update({"ratio_median": "", "ratio_q1": "", "ratio_q3": "", "datasets_ratio_over_1": 0})
    out["share_all5_median"] = rec(v["share_all5_median"])
    out["share_some_median"] = rec(v["share_some_median"])
    return out


def scopes(rows):
    """The pooled reports, all and per species: [(scope, species, rows)]."""
    rows = [dict(r, learner_count=int(r["learner_count"]), resistant=int(r["resistant"])) for r in rows]
    kept = [r for r in rows if pooled_in(r)]
    if not kept:
        raise Stop("no report enters the pooled rows")
    return [("all", "", kept)] + [("species", sp, [r for r in kept if r["species"] == sp])
                                  for sp in sorted({r["species"] for r in kept})]


def pooled_rows(rows):
    return [pool(kept, scope, sp) for scope, sp, kept in scopes(rows)]


def reading_rule(median):
    if median >= STACK:
        return "stack"
    if median > 1:
        return "above"
    return "none shown"


def write_tsv(path, columns, rows):
    with open(path, "w", encoding="utf-8", newline="") as f:
        f.write("\t".join(columns) + "\n")
        for r in rows:
            cells = ["" if r.get(c) is None else str(r.get(c, "")) for c in columns]
            if any("\t" in c or "\n" in c for c in cells):
                raise Stop(f"{os.path.basename(path)}: a cell holds a tab or a newline")
            f.write("\t".join(cells) + "\n")


def read_tsv(path):
    with open(path, encoding="utf-8", newline="") as f:
        lines = f.read().split("\n")
    head = lines[0].split("\t")
    return [dict(zip(head, line.split("\t"))) for line in lines[1:] if line]


def sha256(path):
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def write_record(reading_dir, out):
    with open(os.path.join(reading_dir, READING), encoding="utf-8") as f:
        reading = json.load(f)
    if (reading.get("cutoff"), reading.get("scheme"), reading.get("five")) != (CUTOFF, SCHEME, list(FIVE)):
        raise Stop("the reading was made under other rules")
    fresh_out(out)
    rows = [dataset_row(d) for d in reading["datasets"]]
    pooled = pooled_rows(rows)
    parent = os.path.dirname(os.path.abspath(out))
    tmp = tempfile.mkdtemp(prefix=".overlap-", dir=parent)
    try:
        write_tsv(os.path.join(tmp, "datasets.tsv"), DATASET_COLUMNS, rows)
        write_tsv(os.path.join(tmp, "pooled.tsv"), POOLED_COLUMNS, pooled)
        reading_files = files_under(reading_dir)
        hits, skipped = id_hits(reading_files), 0
        for name in ("datasets.tsv", "pooled.tsv"):
            h, s = record_hits(os.path.join(tmp, name), DECIMAL_COLUMNS)
            hits, skipped = hits + h, skipped + s
        line6 = (f"check 6 ids: {len(reading_files) + 4} files scanned, {hits} id-shaped tokens or barred keys "
                 f"outside the {skipped:,} cells of the computed decimal columns: " + ("pass" if hits == 0 else "STOP"))
        exact = exact_pool(scopes(rows)[0][2])
        rule = reading_rule(exact["ratio_median"]) if exact["ratios"] else "none shown"
        with open(os.path.join(tmp, "gates.txt"), "w", encoding="utf-8") as f:
            f.write("\n".join(reading["checks"][:5] + [line6] + reading["checks"][5:]) + "\n")
            f.write(f"reading rule: {rule}\n")
        with open(os.path.join(tmp, "SHA256SUMS"), "w", encoding="utf-8") as f:
            for name in RECORD_FILES:
                f.write(f"{sha256(os.path.join(tmp, name))}  {name}\n")
        hits += id_hits([os.path.join(tmp, "gates.txt"), os.path.join(tmp, "SHA256SUMS")])
        if hits:
            raise Stop(f"check 6 ids: {hits} id-shaped tokens or barred keys: STOP")
        publish(tmp, out)
    finally:
        if os.path.isdir(tmp):
            shutil.rmtree(tmp)
    return rows, pooled, line6, rule


def write_reading(reading, out):
    fresh_out(out)
    parent = os.path.dirname(os.path.abspath(out))
    tmp = tempfile.mkdtemp(prefix=".overlap-", dir=parent)
    try:
        with open(os.path.join(tmp, READING), "w", encoding="utf-8") as f:
            json.dump(reading, f, indent=1, sort_keys=True)
            f.write("\n")
        hits = id_hits([os.path.join(tmp, READING)])
        if hits:
            raise Stop(f"check 6 ids: {hits} id-shaped tokens or barred keys in the reading: STOP")
        publish(tmp, out)
    finally:
        if os.path.isdir(tmp):
            shutil.rmtree(tmp)


# ---- the page's tables --------------------------------------------------------------------------------------------

def count(x):
    return f"{int(x):,}"


def dec(x, places):
    """A page figure: an exact value rounded once, half to even, to the places the table shows."""
    if x is None or x == "":
        return "none"
    q = round(fractions.Fraction(x) * 10 ** places)
    if q < 0:
        raise Stop("a page figure is under 0")
    whole, part = divmod(q, 10 ** places)
    return f"{whole:,}.{part:0{places}d}" if places else f"{whole:,}"


def tables(record):
    """The page's tables: counts from the record, decimals worked from its counts, never from a rounded cell."""
    pooled = read_tsv(os.path.join(record, "pooled.tsv"))
    rows = read_tsv(os.path.join(record, "datasets.tsv"))
    exact = [exact_pool(kept) for _, _, kept in scopes(rows)]
    if [p["species"] for p in pooled] != [sp for _, sp, _ in scopes(rows)]:
        raise Stop("pooled.tsv and datasets.tsv give other scopes")
    a, e = pooled[0], exact[0]
    out = {}
    out["overlap"] = (["| Of the five | Resistant genomes missed | Susceptible genomes called resistant |",
                       "|---|---:|---:|"]
                      + [f"| {k} | {count(a[f'miss_{k}'])} | {count(a[f'wrong_{k}'])} |" for k in KS]
                      + [f"| 1 to 4 | {count(sum(int(a[f'miss_{k}']) for k in (1, 2, 3, 4)))} "
                         f"| {count(sum(int(a[f'wrong_{k}']) for k in (1, 2, 3, 4)))} |",
                         f"| total | {count(a['resistant'])} | {count(a['susceptible'])} |"])
    out["yardstick"] = ["| Measure | Value |", "|---|---:|",
                        f"| pooled reports | {count(a['datasets'])} |",
                        f"| with a ratio | {count(a['datasets_with_ratio'])} |",
                        f"| yardstick of 0, reported apart | {count(a['datasets_expected_zero'])} |",
                        f"| missed by all five, summed | {count(a['miss_all5'])} |",
                        f"| of those, in a fold with a null fit | {count(a['miss_all5_in_null_fold'])} |",
                        f"| of those, a catalogue call | "
                        f"{count(int(a['miss_all5']) - int(a['miss_all5_no_catalog_call']))} |",
                        f"| of those, no catalogue call | {count(a['miss_all5_no_catalog_call'])} |",
                        f"| of those, on a report with no catalogue rule | {count(a['miss_all5_no_catalog_rule'])} |",
                        f"| reports with no catalogue rule | {count(a['datasets_no_catalog_rule'])} |",
                        f"| yardstick, summed | {dec(e['expected_all5'], 2)} |",
                        f"| median ratio | {dec(e.get('ratio_median'), 1)} |",
                        f"| quartiles of the ratio | {dec(e.get('ratio_q1'), 1)} to {dec(e.get('ratio_q3'), 1)} |",
                        f"| reports with a ratio over 1 | {count(a['datasets_ratio_over_1'])} |",
                        f"| median share missed by all five | {dec(e['share_all5_median'], 3)} |",
                        f"| median share missed by 1 to 4 | {dec(e['share_some_median'], 3)} |"]
    out["species"] = (["| Species | Reports | Resistant | Missed by all five | Yardstick | Median ratio "
                       "| Share by all five | Share by 1 to 4 |", "|---|---:|---:|---:|---:|---:|---:|---:|"]
                      + [f"| *{p['species']}* | {count(p['datasets'])} | {count(p['resistant'])} "
                         f"| {count(p['miss_all5'])} | {dec(x['expected_all5'], 2)} "
                         f"| {dec(x.get('ratio_median'), 1)} "
                         f"| {dec(x['share_all5_median'], 3)} | {dec(x['share_some_median'], 3)} |"
                         for p, x in zip(pooled[1:], exact[1:])])
    out["datasets"] = (["| Species | Antibiotic | Genomes | Resistant | Missed by all five | No catalogue call "
                        "| In a null-fit fold | Yardstick | Ratio |", "|---|---|---:|---:|---:|---:|---:|---:|---:|"]
                       + [f"| *{r['species']}* | {r['antibiotic']} | {count(r['genomes'])} | {count(r['resistant'])} "
                          + (f"| {count(r['miss_all5'])} "
                             f"| {count(r['miss_all5_no_catalog_call']) if r['catalog_rule'] == 'yes' else 'no rule'} "
                             f"| {count(r['miss_all5_in_null_fold'])} | {dec(ratio_of(r)[0], 2)} "
                             f"| {dec(ratio_of(r)[1], 1)} |"
                             if r["miss_all5"] != "" else "| | | | | |")
                          for r in rows])
    return out


# ---- the command line ---------------------------------------------------------------------------------------------

def plan():
    return {"cutoff": CUTOFF, "call": "score >= cutoff", "scheme": SCHEME, "separate": list(SEPARATE),
            "five": list(FIVE), "expect_files": EXPECT_FILES, "min_learners": MIN_LEARNERS,
            "min_resistant": MIN_RESISTANT, "stack": STACK,
            "id_shape": ID_SHAPE.pattern, "datasets_columns": DATASET_COLUMNS, "pooled_columns": POOLED_COLUMNS}


def main(argv=None):
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("plan")
    r = sub.add_parser("read")
    r.add_argument("reports")
    r.add_argument("--out", required=True)
    r.add_argument("--expect", type=int, default=EXPECT_FILES)
    w = sub.add_parser("write")
    w.add_argument("reading")
    w.add_argument("--out", required=True)
    p = sub.add_parser("print")
    p.add_argument("record")
    args = ap.parse_args(argv)
    try:
        if args.cmd == "plan":
            print(json.dumps(plan(), indent=1))
        elif args.cmd == "read":
            fresh_out(args.out)
            reading = read_reports(args.reports, args.expect)
            write_reading(reading, args.out)
            print("\n".join(reading["checks"]))
            print(f"read: {reading['used']} reports used of {reading['found']}; wrote {READING}")
        elif args.cmd == "write":
            rows, pooled, line6, rule = write_record(args.reading, args.out)
            print(line6)
            print(f"wrote datasets.tsv ({len(rows)} rows), pooled.tsv ({len(pooled)} rows), gates.txt, SHA256SUMS")
            print(f"pooled: {pooled[0]['datasets']} reports, median ratio {pooled[0]['ratio_median'] or 'none'}; "
                  f"reading rule: {rule}")
        else:
            for name, lines in tables(args.record).items():
                print(f"## {name}\n" + "\n".join(lines) + "\n")
    except Stop as e:
        print(str(e))
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
