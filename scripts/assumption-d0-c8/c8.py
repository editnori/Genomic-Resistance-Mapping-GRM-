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
"""

import argparse
import os
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


def main(argv=None):
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
    main()
