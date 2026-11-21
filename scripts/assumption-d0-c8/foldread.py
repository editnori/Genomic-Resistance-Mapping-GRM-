"""Readers shared by C8 and D0 of the assumption plan.

They read what a dataset report keeps about its outer folds (`fold_models`,
docs/report-rule-fields.md), the feature cards' catalogue families, and the
test fold of each sample. Nothing here fits, draws or ranks anything.

A feature is a row of the dataset matrix, the numbering `models`,
`features` and `fold_models` share. A catalogue family is named only as the
report's card gives it, which is the AMRFinderPlus catalog's own name.
"""

import base64
import csv
import json
import os
import statistics

LEARNERS = ("scm", "cart", "rf", "gbt", "lasso")
SCHEMES = ("lineage", "random")
# The learners whose fold model is its whole feature list and whose
# features have a presence pattern in `fold_models.features`.
PATTERN_LEARNERS = ("scm", "cart", "lasso")


class ReportError(ValueError):
    """A report lacks what a reader needs, or holds a value it cannot read."""


def load_report(path, need_models=True):
    with open(path, encoding="utf-8") as f:
        report = json.load(f)
    if need_models and "fold_models" not in report:
        raise ReportError(
            f"{os.path.basename(path)}: no fold_models; the reports of release "
            "rev14-274-20261112 hold none, so rebuild the report with grm learn first"
        )
    return report


def split(report, scheme):
    """The report's split of `scheme`, or None when it did not run."""
    for s in report["splits"]:
        if s["scheme"] == scheme:
            return None if s.get("unavailable_reason") else s
    return None


def fold_split(report, scheme):
    """The `fold_models` entry of `scheme`, or None."""
    for s in report["fold_models"]["splits"]:
        if s["scheme"] == scheme:
            return s
    return None


def chosen(model, learner):
    """The features one fold model of `learner` chooses, as a set of rows.

    scm: its rules' features. cart: the features its splits test. rf and
    gbt: the up to 50 heaviest features the report keeps for the fold.
    lasso: the features of its nonzero weights. A null lasso (`None`) and a
    boosting list with no entry choose nothing.
    """
    if model is None:
        return frozenset()
    if learner == "scm":
        return frozenset(int(r["feature"]) for r in model["rules"])
    if learner == "cart":
        return frozenset(
            int(n["split"]["feature"]) for n in model["nodes"] if n.get("split")
        )
    if learner in ("rf", "gbt"):
        return frozenset(int(w["feature"]) for w in model)
    if learner == "lasso":
        return frozenset(int(c[0]) for c in model.get("coef", []))
    raise ReportError(f"unknown learner {learner}")


def fold_sets(report, scheme, learner):
    """Per outer fold, the set of features `learner` chose; None when the
    split did not run. Entry f is the model fitted with fold f held out."""
    entry = fold_split(report, scheme)
    if entry is None:
        return None
    models = entry[learner]
    s = split(report, scheme)
    if s is None:
        raise ReportError(f"{report['id']}: fold_models has {scheme} but the split did not run")
    if len(models) != s["folds"]:
        raise ReportError(
            f"{report['id']}: {scheme} {learner} has {len(models)} fold models for {s['folds']} folds"
        )
    return [chosen(m, learner) for m in models]


FINAL_KEYS = {"scm": "scm", "cart": "cart", "rf": "rf_importance", "gbt": "gbt_importance", "lasso": "lasso"}


def final_set(report, learner):
    """The features the final model of `learner` chooses (`models`, the same
    shapes as the fold models)."""
    return chosen(report["models"].get(FINAL_KEYS[learner]), learner)


def decode_bits(text, n):
    """A presence pattern as an int: sample i is bit i. `text` is base64 of
    little-endian 64-bit words, as on a feature card."""
    raw = base64.b64decode(text, validate=True)
    if len(raw) % 8:
        raise ReportError("a pattern is not a whole number of 64-bit words")
    value = int.from_bytes(raw, "little")
    if value >> n:
        raise ReportError("a pattern sets a bit past the last sample")
    return value


def patterns(report):
    """Feature row to presence pattern (int), for every feature that has one:
    the fold models' patterns, then the cards'."""
    n = len(report["samples"])
    out = {}
    for f in report["fold_models"]["features"]:
        out[int(f["feature"])] = decode_bits(f["bits"], n)
    for c in report["features"]:
        p = decode_bits(c["bits"], n)
        if out.setdefault(int(c["feature"]), p) != p:
            raise ReportError(f"{report['id']}: feature {c['feature']} has two patterns")
    return out


def card_families(report):
    """Feature row to its card's catalogue families; a feature with no card
    is absent from the map."""
    return {int(c["feature"]): list(c["catalog_families"]) for c in report["features"]}


def jaccard(a, b):
    """|a & b| / |a | b|; None when both are empty."""
    union = len(a | b)
    if union == 0:
        return None
    return len(a & b) / union


def pairwise(sets):
    """Jaccard of every pair of folds, in fold order, leaving out a pair
    whose two folds both chose nothing."""
    out = []
    for i in range(len(sets)):
        for j in range(i + 1, len(sets)):
            v = jaccard(sets[i], sets[j])
            if v is not None:
                out.append(v)
    return out


def train_mask(test_fold, folds):
    """Per pair of folds (f, g), the samples that train in both, as an int."""
    masks = {}
    for f in range(folds):
        for g in range(f + 1, folds):
            m = 0
            for i, t in enumerate(test_fold):
                if t is not None and t != f and t != g:
                    m |= 1 << i
            masks[(f, g)] = m
    return masks


def pairwise_by_pattern(sets, pats, test_fold, folds):
    """As `pairwise`, but two features count as one when their patterns agree
    on every sample that trains in both folds of the pair. Every feature
    must have a pattern."""
    masks = train_mask(test_fold, folds)
    out = []
    for i in range(len(sets)):
        for j in range(i + 1, len(sets)):
            m = masks[(i, j)]
            try:
                a = {pats[x] & m for x in sets[i]}
                b = {pats[x] & m for x in sets[j]}
            except KeyError as e:
                raise ReportError(f"feature {e.args[0]} has no pattern") from None
            v = jaccard(a, b)
            if v is not None:
                out.append(v)
    return out


def median(values):
    return statistics.median(values) if values else None


def quartiles(values):
    """(q1, median, q3) by the inclusive method; None for an empty list."""
    if not values:
        return None
    if len(values) == 1:
        return (values[0], values[0], values[0])
    q = statistics.quantiles(values, n=4, method="inclusive")
    return (q[0], q[1], q[2])


def read_marks(path):
    """Datasets marked M: `intcount`'s cut table shows an S+I or I+R pair."""
    with open(path, encoding="utf-8", newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    need = {"dataset", "pairs_S+I", "pairs_I+R"}
    if not rows or not need <= set(rows[0]):
        raise ReportError(f"{path}: needs the columns {sorted(need)}")
    return {r["dataset"]: int(r["pairs_S+I"]) > 0 or int(r["pairs_I+R"]) > 0 for r in rows}


def fmt(v, places=6):
    if v is None:
        return "NA"
    if isinstance(v, bool):
        return "yes" if v else "no"
    if isinstance(v, float):
        return f"{v:.{places}f}"
    return str(v)


def write_tsv(path, header, rows):
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8", newline="") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            if len(r) != len(header):
                raise ValueError(f"{path}: a row has {len(r)} cells for {len(header)} columns")
            f.write("\t".join(fmt(v) for v in r) + "\n")
    os.replace(tmp, path)
