"""C8 of the assumption plan: do the learners choose the same inputs across
folds, and how many of the chosen inputs carry a catalogue family?

    python3 -I c8.py --reports DIR --marks DATASETS_TSV --out DIR

It reads each report's `fold_models` and feature cards and fits nothing.
Per dataset and learner it writes, for the lineage and the random split:

- the median pairwise Jaccard of the five folds' chosen feature sets, with a
  feature counted by its matrix row (the strict reading,
  docs/report-rule-fields.md "Feature numbers");
- for scm, cart and lasso, the same median with two features counted as one
  when their patterns agree on the samples both folds train on;
- the distinct features chosen over the split's folds, how many of them have
  a feature card, and how many of those cards name a catalogue family.

Only final-model features have a card, so a fold feature without one has no
note to read: it is counted apart and never as "no family". No feature, row,
k-mer or family is written: the tables hold counts and overlaps only.

The per-learner summary gives the median and quartiles of the lineage and
random overlaps over datasets, the datasets under 0.2 and over 0.8, and the
pooled share of chosen features whose card names a family, for all
datasets, for those `intcount` marks M and for the rest.

    python3 -I c8.py combine --out DIR --one-fold zero|undefined \
        --over-five pairs|dataset-median --bands none|LOW,HIGH [--names TSV] C8_DIR...

`combine` reads the c8-datasets.tsv of each species (one C8_DIR each) and
writes, per learner and over the five, the plan's measure: the median over
datasets of the per-dataset median pairwise Jaccard, with quartiles, beside
it the distinct features chosen and those whose card names a catalogue
family, the units counted and those left undefined, and why. The rules
that are open (docs/assumption-plan.md leaves them to the lead) are its
arguments, and it has no default for them. It refuses any cell that is not
a count, an overlap, a mark, a learner or a dataset name, and never prints
such a cell.
"""

import argparse
import csv
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import foldread as fr  # noqa: E402

LOW, HIGH = 0.2, 0.8

DATASET_HEADER = (
    ["dataset", "species", "m_mark", "learner"]
    + [
        f"{s}_{k}"
        for s in fr.SCHEMES
        for k in (
            "folds",
            "pairs",
            "median_jaccard",
            "median_jaccard_pattern",
            "distinct",
            "with_card",
            "with_family",
        )
    ]
)

SUMMARY_HEADER = [
    "subset",
    "learner",
    "scheme",
    "datasets",
    "with_overlap",
    "q1",
    "median",
    "q3",
    "under_0.2",
    "over_0.8",
    "pattern_median",
    "distinct_features",
    "with_card",
    "with_family",
    "family_share_of_carded",
]


FINAL_HEADER = ["dataset", "species", "m_mark", "learner", "features", "with_card", "with_family"]
FINAL_SUMMARY_HEADER = ["subset", "learner", "datasets", "with_a_feature", "with_a_family_feature",
                        "features", "with_card", "with_family", "family_share_of_carded"]


def final_rows(report, mark):
    """Per learner, the final model's features, how many have a card and how
    many cards name a catalogue family. A released report has these without
    fold models."""
    fams = fr.card_families(report)
    out = []
    for learner in fr.LEARNERS:
        feats = fr.final_set(report, learner)
        carded = [f for f in feats if f in fams]
        out.append([report["id"], report["species_slug"], mark, learner, len(feats), len(carded),
                    sum(1 for f in carded if fams[f])])
    return out


def final_summary(rows):
    col = {h: i for i, h in enumerate(FINAL_HEADER)}
    out = []
    for name, keep in (("all", lambda r: True), ("not_M", lambda r: r[col["m_mark"]] is False),
                       ("M", lambda r: r[col["m_mark"]] is True), ("unmarked", lambda r: r[col["m_mark"]] is None)):
        chosen = [r for r in rows if keep(r)]
        if name == "unmarked" and not chosen:
            continue
        # The last row, any_summed: its dataset columns count a dataset once
        # when any learner qualifies it; its feature columns are summed over
        # the five learners, so a feature two learners chose counts twice.
        for learner in list(fr.LEARNERS) + ["any_summed"]:
            lr = [r for r in chosen if learner == "any_summed" or r[col["learner"]] == learner]
            if learner == "any_summed":
                # One row per dataset: is any learner's final feature carded with a family?
                by = {}
                for r in lr:
                    d = by.setdefault(r[col["dataset"]], [0, 0])
                    d[0] += r[col["features"]] > 0
                    d[1] += r[col["with_family"]] > 0
                datasets, with_feat, with_fam = len(by), sum(v[0] > 0 for v in by.values()), \
                    sum(v[1] > 0 for v in by.values())
            else:
                datasets = len(lr)
                with_feat = sum(r[col["features"]] > 0 for r in lr)
                with_fam = sum(r[col["with_family"]] > 0 for r in lr)
            feats = sum(r[col["features"]] for r in lr)
            carded = sum(r[col["with_card"]] for r in lr)
            family = sum(r[col["with_family"]] for r in lr)
            out.append([name, learner, datasets, with_feat, with_fam, feats, carded, family,
                        family / carded if carded else None])
    return out


def dataset_rows(report, mark):
    """One row per learner for one report."""
    pats = fr.patterns(report)
    fams = fr.card_families(report)
    out = []
    for learner in fr.LEARNERS:
        row = [report["id"], report["species_slug"], mark, learner]
        for scheme in fr.SCHEMES:
            sets = fr.fold_sets(report, scheme, learner)
            if sets is None:
                row += [0, 0, None, None, None, None, None]
                continue
            s = fr.split(report, scheme)
            pairs = fr.pairwise(sets)
            pat = None
            if learner in fr.PATTERN_LEARNERS:
                pat = fr.median(fr.pairwise_by_pattern(sets, pats, s["test_fold"], s["folds"]))
            union = frozenset().union(*sets)
            carded = [f for f in union if f in fams]
            family = [f for f in carded if fams[f]]
            row += [len(sets), len(pairs), fr.median(pairs), pat, len(union), len(carded), len(family)]
        out.append(row)
    return out


def summarize(rows):
    """Per subset, learner and scheme: the spread of overlaps over datasets."""
    col = {h: i for i, h in enumerate(DATASET_HEADER)}
    subsets = (
        ("all", lambda r: True),
        ("not_M", lambda r: r[col["m_mark"]] is False),
        ("M", lambda r: r[col["m_mark"]] is True),
        ("unmarked", lambda r: r[col["m_mark"]] is None),
    )
    out = []
    for name, keep in subsets:
        chosen = [r for r in rows if keep(r)]
        if name == "unmarked" and not chosen:
            continue
        for learner in fr.LEARNERS:
            lr = [r for r in chosen if r[col["learner"]] == learner]
            for scheme in fr.SCHEMES:
                ran = [r for r in lr if r[col[f"{scheme}_folds"]]]
                vals = [r[col[f"{scheme}_median_jaccard"]] for r in ran]
                vals = [v for v in vals if v is not None]
                pv = [r[col[f"{scheme}_median_jaccard_pattern"]] for r in ran]
                pv = [v for v in pv if v is not None]
                q = fr.quartiles(vals) or (None, None, None)
                distinct = sum(r[col[f"{scheme}_distinct"]] for r in ran)
                carded = sum(r[col[f"{scheme}_with_card"]] for r in ran)
                family = sum(r[col[f"{scheme}_with_family"]] for r in ran)
                out.append(
                    [
                        name,
                        learner,
                        scheme,
                        len(ran),
                        len(vals),
                        q[0],
                        q[1],
                        q[2],
                        sum(v < LOW for v in vals),
                        sum(v > HIGH for v in vals),
                        fr.median(pv) if learner in fr.PATTERN_LEARNERS else None,
                        distinct,
                        carded,
                        family,
                        family / carded if carded else None,
                    ]
                )
    return out


# ---------------------------------------------------------------- combine

NAMES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "outside-study", "frozen", "datasets.tsv")
DATASET_RE = re.compile(r"[a-z][a-z0-9_]*__[a-z0-9][a-z0-9_]*\Z")
SPECIES_RE = re.compile(r"[a-z][a-z0-9_]*\Z")
FIELDS = ("folds", "pairs", "median_jaccard", "median_jaccard_pattern", "distinct", "with_card", "with_family")
SPLIT_NOT_RUN = "the split did not run"
NO_FOLD = "no fold chose a feature"
ONE_FOLD = "one fold chose a feature"
NO_LEARNER = "no learner defined"
OVER_FIVE = "over_five"
READS = {
    "lineage": "median pairwise Jaccard of the five lineage folds' chosen features; a pair where both folds "
               "chose nothing is left out; decides nothing alone",
    "random": "the random split: shown, decides nothing",
}
COMBINED_HEADER = ["species", "antibiotic", "dataset", "learner", "scheme", "folds", "folds_choosing", "pairs",
                   "median_jaccard", "median_jaccard_pattern", "distinct", "with_card",
                   "with_amrfinderplus_catalog_family", "defined"]
COMBINED_SUMMARY_HEADER = [
    "learner", "scheme", "unit", "units", "defined", "undefined_split_not_run", "undefined_no_fold_chose",
    "undefined_one_fold_chose", "undefined_no_learner_defined", "q1", "median", "q3",
    "units_fewer_than_two_folds_chose", "pattern_defined",
    "pattern_median", "band_low", "band_high", "under_band_low", "over_band_high", "distinct_summed",
    "with_card_summed", "with_amrfinderplus_catalog_family_summed", "catalog_family_share_of_carded", "reads"]


class CombineError(ValueError):
    """A C8 table the combine does not read on from."""


def species_name(slug):
    """A species slug in full, as the data names it: `klebsiella_aerogenes`
    is Klebsiella aerogenes."""
    return slug.replace("_", " ").capitalize()


def read_names(path):
    """Dataset id to (species in full, antibiotic in full) from the frozen
    table of the release's 274 datasets."""
    with open(path, encoding="utf-8", newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t", quoting=csv.QUOTE_NONE))
    if not rows or not {"dataset_id", "species_slug", "antibiotic"} <= set(rows[0]):
        raise CombineError(f"{path}: needs the columns dataset_id, species_slug and antibiotic")
    out = {}
    for r in rows:
        d, sp, ab = r["dataset_id"], r["species_slug"], r["antibiotic"].strip()
        if d in out or not ab or not d.startswith(sp + "__"):
            raise CombineError(f"{path}: the row of {d} is twice, has no antibiotic or another species")
        out[d] = (species_name(sp), ab)
    return out


LABEL = {"dataset": "a dataset name", "species": "a species name", "mark": "a mark", "learner": "a learner",
         "int": "a count", "int_na": "a count", "overlap": "an overlap"}


def _cell(path, line, col, text, kind):
    """One typed cell; the error names the place and never the text."""
    where = f"{path}: line {line}, column {col}"
    if kind == "dataset":
        if DATASET_RE.match(text):
            return text
    elif kind == "species":
        if SPECIES_RE.match(text):
            return text
    elif kind == "mark":
        if text in ("yes", "no", "NA"):
            return {"yes": True, "no": False, "NA": None}[text]
    elif kind == "learner":
        if text in fr.LEARNERS:
            return text
    elif kind in ("int", "int_na"):
        if kind == "int_na" and text == "NA":
            return None
        if text.isdigit() and text.isascii():
            return int(text)
    elif kind == "overlap":
        if text == "NA":
            return None
        if re.fullmatch(r"[01]\.\d{6}", text):
            v = float(text)
            if 0.0 <= v <= 1.0:
                return v
    raise CombineError(f"{where}: not {LABEL[kind]}; a C8 table holds counts, overlaps, marks, learners and "
                       "dataset names only, and this cell is not printed")


def _kind(col):
    if col in ("dataset", "species", "learner"):
        return col
    if col == "m_mark":
        return "mark"
    field = col.split("_", 1)[1]
    if field in ("folds", "pairs"):
        return "int"
    if field.startswith("median_jaccard"):
        return "overlap"
    return "int_na"


def choosing(folds, pairs):
    """The counts of folds that chose a feature that give `pairs` kept pairs:
    a pair is left out only when both its folds chose nothing, so k folds
    that chose give C(folds, 2) - C(folds - k, 2) pairs."""
    def c2(n):
        return n * (n - 1) // 2
    return [k for k in range(folds + 1) if c2(folds) - c2(folds - k) == pairs]


def check_scheme(where, learner, s):
    """The cells of one split must fit together as c8.py writes them."""
    if s["folds"] == 0:
        if s["pairs"] != 0 or any(s[k] is not None for k in FIELDS[2:]):
            raise CombineError(f"{where}: a split that did not run holds a value")
        return
    if any(s[k] is None for k in ("distinct", "with_card", "with_family")):
        raise CombineError(f"{where}: a split that ran lacks a count")
    ks = choosing(s["folds"], s["pairs"])
    if not ks:
        raise CombineError(f"{where}: {s['pairs']} pairs cannot come from {s['folds']} folds")
    if (s["pairs"] == 0) != (s["median_jaccard"] is None):
        raise CombineError(f"{where}: the overlap is undefined exactly when no pair is kept, and it is not")
    pattern = learner in fr.PATTERN_LEARNERS
    if pattern and (s["pairs"] == 0) != (s["median_jaccard_pattern"] is None):
        raise CombineError(f"{where}: the pattern overlap is undefined exactly when no pair is kept, and it is not")
    if not pattern and s["median_jaccard_pattern"] is not None:
        raise CombineError(f"{where}: {learner} has no pattern overlap")
    if not s["with_family"] <= s["with_card"] <= s["distinct"]:
        raise CombineError(f"{where}: with a family, with a card and distinct are not in that order")
    if (s["distinct"] == 0) != (ks == [0]):
        raise CombineError(f"{where}: features chosen and folds that chose disagree")
    if ks == [1] and (s["median_jaccard"] != 0.0 or (pattern and s["median_jaccard_pattern"] != 0.0)):
        raise CombineError(f"{where}: one fold chose, so every kept pair is 0, and the overlap is not")


def read_table(path):
    """One species' c8-datasets.tsv: typed rows, every cell checked."""
    with open(path, encoding="utf-8", newline="") as f:
        lines = list(csv.reader(f, delimiter="\t", quoting=csv.QUOTE_NONE))
    if not lines or lines[0] != DATASET_HEADER:
        raise CombineError(f"{path}: its columns are not c8.py's")
    kinds = [_kind(c) for c in DATASET_HEADER]
    rows = []
    for n, cells in enumerate(lines[1:], start=2):
        if len(cells) != len(DATASET_HEADER):
            raise CombineError(f"{path}: line {n} has {len(cells)} cells, want {len(DATASET_HEADER)}")
        v = [_cell(path, n, c, t, k) for c, t, k in zip(DATASET_HEADER, cells, kinds)]
        r = dict(zip(DATASET_HEADER, v))
        if not r["dataset"].startswith(r["species"] + "__"):
            raise CombineError(f"{path}: line {n}: {r['dataset']} is not of {r['species']}")
        row = {"dataset": r["dataset"], "species": r["species"], "m_mark": r["m_mark"], "learner": r["learner"]}
        for s in fr.SCHEMES:
            row[s] = {k: r[f"{s}_{k}"] for k in FIELDS}
            check_scheme(f"{path}: line {n}, {s}", r["learner"], row[s])
        rows.append(row)
    return rows


def define(s, one_fold):
    """(value, pattern value, reason left undefined or None, folds that chose)."""
    if s["folds"] == 0:
        return None, None, SPLIT_NOT_RUN, "NA"
    ks = choosing(s["folds"], s["pairs"])
    text = str(ks[0]) if len(ks) == 1 else f"{ks[0]}-{ks[-1]}"
    if s["pairs"] == 0:
        return None, None, NO_FOLD, text
    if min(ks) < 2 <= max(ks):
        raise CombineError(f"{s['pairs']} pairs of {s['folds']} folds cannot tell one fold that chose from two")
    if ks == [1] and one_fold == "undefined":
        return None, None, ONE_FOLD, text
    return s["median_jaccard"], s["median_jaccard_pattern"], None, text


def combine_tables(tables, one_fold, over_five, bands, names):
    """tables: (path, rows) per species. Returns the long dataset rows and
    the summary rows. `one_fold`, `over_five` and `bands` are the open
    rules; see the module text."""
    if one_fold not in ("zero", "undefined"):
        raise CombineError("--one-fold must be zero or undefined")
    if over_five not in ("pairs", "dataset-median"):
        raise CombineError("--over-five must be pairs or dataset-median")
    seen_species, seen = {}, {}
    for path, rows in tables:
        species = {r["species"] for r in rows}
        if len(species) != 1:
            raise CombineError(f"{path}: holds {len(species)} species, want one")
        (sp,) = species
        if sp in seen_species:
            raise CombineError(f"{path}: species {sp} again, first in {seen_species[sp]}")
        seen_species[sp] = path
        for r in rows:
            key = (r["dataset"], r["learner"])
            if key in seen:
                raise CombineError(f"{path}: {r['dataset']} {r['learner']} twice")
            seen[key] = r
    datasets = sorted({d for d, _ in seen})
    for d in datasets:
        have = [l for l in fr.LEARNERS if (d, l) in seen]
        if len(have) != len(fr.LEARNERS):
            raise CombineError(f"{d}: rows for {len(have)} of the five learners")
        if d not in names:
            raise CombineError(f"{d}: no species and antibiotic name in the names table")
    long_rows, cells = [], {}
    for d in datasets:
        sp, ab = names[d]
        for l in fr.LEARNERS:
            r = seen[(d, l)]
            for s in fr.SCHEMES:
                v, p, why, k = define(r[s], one_fold)
                cells[(d, l, s)] = (v, p, why, r[s])
                long_rows.append([sp, ab, d, l, s, r[s]["folds"], k, r[s]["pairs"], v, p, r[s]["distinct"],
                                  r[s]["with_card"], r[s]["with_family"], "yes" if why is None else why])
    low = high = None
    if bands != "none":
        try:
            low, high = (float(x) for x in bands.split(","))
        except ValueError:
            raise CombineError("--bands must be none or LOW,HIGH") from None
        if not 0.0 <= low < high <= 1.0:
            raise CombineError("--bands must hold 0 <= LOW < HIGH <= 1")
    summary = []
    for s in fr.SCHEMES:
        for l in list(fr.LEARNERS) + [OVER_FIVE]:
            why = {SPLIT_NOT_RUN: 0, NO_FOLD: 0, ONE_FOLD: 0, NO_LEARNER: 0}
            learners = fr.LEARNERS if l == OVER_FIVE else (l,)
            if l == OVER_FIVE and over_five == "dataset-median":
                unit, vals, pvals = "dataset", [], []
                for d in datasets:
                    mine = [cells[(d, x, s)] for x in learners]
                    got = [c[0] for c in mine if c[2] is None]
                    pg = [c[1] for x, c in zip(learners, mine) if c[2] is None and x in fr.PATTERN_LEARNERS]
                    if got:
                        vals.append(fr.median(got))
                    else:
                        why[NO_LEARNER] += 1
                    if pg:
                        pvals.append(fr.median(pg))
                units = len(datasets)
            else:
                unit = "dataset and learner"
                mine = [cells[(d, x, s)] for d in datasets for x in learners]
                vals = [c[0] for c in mine if c[2] is None]
                pvals = [c[1] for d in datasets for x in learners
                         for c in [cells[(d, x, s)]] if c[2] is None and x in fr.PATTERN_LEARNERS]
                for c in mine:
                    if c[2] is not None:
                        why[c[2]] += 1
                units = len(mine)
            ran = [cells[(d, x, s)][3] for d in datasets for x in learners if cells[(d, x, s)][3]["folds"]]
            # Ruling 1: the count of units where fewer than two folds chose (none, or one alone) is said beside
            # the median, whether one fold alone counts 0 or is left undefined.
            fewer = None if unit == "dataset" else sum(max(choosing(c["folds"], c["pairs"])) < 2 for c in ran)
            q = fr.quartiles(vals) or (None, None, None)
            carded = sum(c["with_card"] for c in ran)
            family = sum(c["with_family"] for c in ran)
            has_pattern = l == OVER_FIVE or l in fr.PATTERN_LEARNERS
            summary.append([
                l, s, unit, units, len(vals), why[SPLIT_NOT_RUN], why[NO_FOLD], why[ONE_FOLD],
                why[NO_LEARNER] if unit == "dataset" else None, q[0], q[1], q[2], fewer,
                len(pvals) if has_pattern else None, fr.median(pvals) if has_pattern else None,
                low, high, None if low is None else sum(v < low for v in vals),
                None if high is None else sum(v > high for v in vals),
                sum(c["distinct"] for c in ran), carded, family, family / carded if carded else None, READS[s]])
    return long_rows, summary


def combine_main(argv):
    ap = argparse.ArgumentParser(prog="c8.py combine", description="C8 over species")
    ap.add_argument("--out", required=True)
    ap.add_argument("--one-fold", required=True, choices=("zero", "undefined"),
                    help="a dataset and learner where one fold alone chose a feature: its kept pairs are all 0, "
                         "so its overlap is 0 (zero), or it is left undefined (undefined)")
    ap.add_argument("--over-five", required=True, choices=("pairs", "dataset-median"),
                    help="the row over the five learners: every dataset and learner value pooled (pairs), or per "
                         "dataset the median of its defined learners (dataset-median)")
    ap.add_argument("--bands", required=True, help="none, or LOW,HIGH: also count the values under LOW and over HIGH")
    ap.add_argument("--names", default=NAMES, help="dataset_id, species_slug and antibiotic per dataset")
    ap.add_argument("--expect-datasets", type=int, help="datasets the tables must hold")
    ap.add_argument("dirs", nargs="+", help="directories holding one species' c8-datasets.tsv")
    a = ap.parse_args(argv)
    try:
        tables = [(p, read_table(os.path.join(p, "c8-datasets.tsv"))) for p in a.dirs]
        long_rows, summary = combine_tables(tables, a.one_fold, a.over_five, a.bands, read_names(a.names))
        n = len({r[2] for r in long_rows})
        if a.expect_datasets is not None and n != a.expect_datasets:
            raise CombineError(f"{n} datasets, expected {a.expect_datasets}")
    except CombineError as e:
        print(f"stop: {e}", file=sys.stderr)
        return 2
    os.makedirs(a.out, exist_ok=True)
    fr.write_tsv(os.path.join(a.out, "c8-combined-datasets.tsv"), COMBINED_HEADER, long_rows)
    fr.write_tsv(os.path.join(a.out, "c8-combined-summary.tsv"), COMBINED_SUMMARY_HEADER, summary)
    print(f"c8 combine: {len(tables)} species, {n} datasets, {len(long_rows)} rows")
    return 0


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if argv[:1] == ["combine"]:
        return combine_main(argv[1:])
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--reports", required=True, help="a directory of dataset reports with fold_models")
    ap.add_argument("--marks", help="intcount's cut table (datasets.tsv); without it no dataset is marked")
    ap.add_argument("--only", help="a file of dataset ids, one per line; default every report in --reports")
    ap.add_argument("--out", required=True)
    ap.add_argument("--final", action="store_true",
                    help="count the final models' features only; for reports without fold_models")
    a = ap.parse_args(argv)
    marks = fr.read_marks(a.marks) if a.marks else {}
    names = sorted(n for n in os.listdir(a.reports) if n.endswith(".json"))
    if a.only:
        with open(a.only, encoding="utf-8") as f:
            want = [line.strip() for line in f if line.strip()]
        have = {n[:-5] for n in names}
        missing = [d for d in want if d not in have]
        if missing:
            raise SystemExit(f"{len(missing)} listed datasets have no report, the first {missing[0]}")
        names = [d + ".json" for d in sorted(want)]
    if not names:
        raise SystemExit(f"no report in {a.reports}")
    rows = []
    for n in names:
        report = fr.load_report(os.path.join(a.reports, n), need_models=not a.final)
        if marks and report["id"] not in marks:
            raise SystemExit(f"{report['id']} is not in the marks table")
        rows += (final_rows if a.final else dataset_rows)(report, marks.get(report["id"]))
    os.makedirs(a.out, exist_ok=True)
    if a.final:
        fr.write_tsv(os.path.join(a.out, "c8-final-datasets.tsv"), FINAL_HEADER, rows)
        fr.write_tsv(os.path.join(a.out, "c8-final-summary.tsv"), FINAL_SUMMARY_HEADER, final_summary(rows))
        print(f"c8 final models: {len(names)} reports, {len(rows)} rows")
        return
    fr.write_tsv(os.path.join(a.out, "c8-datasets.tsv"), DATASET_HEADER, rows)
    fr.write_tsv(os.path.join(a.out, "c8-summary.tsv"), SUMMARY_HEADER, summarize(rows))
    print(f"c8: {len(names)} reports, {len(rows)} rows")


if __name__ == "__main__":
    sys.exit(main())
