"""Is a refit report the release's rule? Equality counts only.

    python3 -I refit_check.py REFIT_RESULTS RELEASE_DIR DATASETS OUT_DIR

REFIT_RESULTS holds the refit's <dataset>.json; RELEASE_DIR the unpacked
release (datasets/<dataset>.json); DATASETS a comma list. For every dataset
it compares the refit report with the release's report:

  report fields  every top-level key of the release's report except
                 `seconds` must be equal; `samples` is compared on each
                 sample's fitted fields (id, resistant, lineage, from_mic,
                 mic, groups, catalog_resistant) and, apart, on the four
                 descriptive fields (country, year, host, mlst), which come
                 from genome metadata the fingerprint does not hash;
  learners       for each split and learner, every field but `seconds`
                 must be equal; a learner that differs is named with the
                 count of differing fields and the largest absolute
                 difference over its numbers;
  added fields   `fold_models` (one entry per split, each learner with one
                 model per fold), `samples[].rule_groups`, and `selection`
                 are allowed; any other key the release lacks is counted;
  feature cards  the top-level `features` list is compared card by card, in
                 order, one to one: a card's `kmer` and `locus` may differ,
                 since both are written after the fit from genome files and
                 CDS lists fetched again; every other card field, the number
                 and the order of the cards must be equal. The cards whose
                 `kmer` or `locus` differs are counted per dataset. A `kmer`
                 or `locus` key anywhere else is compared as any field.

A column compares the two files' text: with `seconds`, the three allowed
fields and the four descriptive fields taken out, and the values of each
top-level feature card's `kmer` and `locus` set to one mark (the keys stay),
both are written again in their own key order and must be the same string,
which also checks the key order.

Writes OUT_DIR/datasets.tsv (one row per dataset: field names that differ,
no value; its last column, cards_with_other_locus, counts the feature cards
whose `kmer` or `locus` differs) and OUT_DIR/learners.tsv (one row per
dataset, split and learner), and prints the totals. Exit 1 when any fitted
field or learner differs. No genome id, value or feature is written.
"""

import json
import math
from collections import OrderedDict
import os
import sys

FITTED = ("id", "resistant", "lineage", "from_mic", "mic", "groups", "catalog_resistant")
DESCRIPTIVE = ("country", "year", "host", "mlst")
ALLOWED_NEW = {"fold_models", "selection"}
LEARNERS = ("scm", "cart", "rf", "gbt", "lasso")
DROP = {"seconds", "fold_models", "rule_groups", "selection"} | set(DESCRIPTIVE)
# Written after the fit from refetched genome files and CDS lists: allowed to
# differ in a top-level feature card, and nowhere else.
CARD_ALLOWED = ("kmer", "locus")
MARK = "[allowed to differ]"


def masked(card):
    """A feature card with the allowed fields' values set to one mark."""
    if not isinstance(card, dict):
        return card
    return OrderedDict((k, MARK if k in CARD_ALLOWED else v) for k, v in card.items())


def cards(a, b):
    """(equal but for the allowance, cards whose kmer or locus differs)."""
    if not isinstance(a, list) or not isinstance(b, list):
        return a == b, 0
    other = sum(isinstance(x, dict) and isinstance(y, dict) and any(x.get(k) != y.get(k) for k in CARD_ALLOWED)
                for x, y in zip(a, b))
    return len(a) == len(b) and all(masked(x) == masked(y) for x, y in zip(a, b)), other


def stripped(x):
    if isinstance(x, dict):
        return OrderedDict((k, stripped(v)) for k, v in x.items() if k not in DROP)
    if isinstance(x, list):
        return [stripped(v) for v in x]
    return x


def numbers(x):
    """Every number in x, in a fixed walk order, with the walk's shape."""
    if isinstance(x, bool) or x is None or isinstance(x, str):
        yield ("v", x)
    elif isinstance(x, (int, float)):
        yield ("n", x)
    elif isinstance(x, list):
        yield ("l", len(x))
        for v in x:
            yield from numbers(v)
    elif isinstance(x, dict):
        yield ("d", tuple(sorted(x)))
        for k in sorted(x):
            yield from numbers(x[k])


def gap(a, b):
    """(same shape, largest absolute difference over the numbers)."""
    wa, wb = list(numbers(a)), list(numbers(b))
    if len(wa) != len(wb):
        return False, math.inf
    worst = 0.0
    for (ka, va), (kb, vb) in zip(wa, wb):
        if ka != kb:
            return False, math.inf
        if ka == "n":
            if isinstance(va, float) and isinstance(vb, float) and math.isnan(va) and math.isnan(vb):
                continue
            worst = max(worst, abs(va - vb))
        elif va != vb:
            return False, math.inf
    return True, worst


def text_form(text):
    """The report as the text check compares it."""
    doc = stripped(json.loads(text, object_pairs_hook=OrderedDict))
    if isinstance(doc, dict) and isinstance(doc.get("features"), list):
        doc["features"] = [masked(c) for c in doc["features"]]
    return json.dumps(doc)


def check(refit, release):
    row = {"report_fields_differ": [], "descriptive_samples_differ": 0, "unexpected_new_keys": [],
           "fold_models": "no", "rule_groups_samples": 0, "cards_with_other_locus": 0}
    for k in release:
        if k == "seconds":
            continue
        if k == "samples":
            a, b = refit.get("samples", []), release["samples"]
            if len(a) != len(b) or any({f: s.get(f) for f in FITTED} != {f: t.get(f) for f in FITTED}
                                       for s, t in zip(a, b)):
                row["report_fields_differ"].append("samples")
            row["descriptive_samples_differ"] = sum(any(s.get(f) != t.get(f) for f in DESCRIPTIVE)
                                                    for s, t in zip(a, b))
            row["rule_groups_samples"] = sum("rule_groups" in s for s in a)
            continue
        if k == "splits":
            a, b = refit.get("splits", []), release["splits"]
            strip = lambda s: {x: y for x, y in s.items() if x != "learners"}
            if len(a) != len(b) or any(strip(s) != strip(t) for s, t in zip(a, b)):
                row["report_fields_differ"].append("splits")
            continue
        if k == "features":
            same, row["cards_with_other_locus"] = cards(refit.get("features"), release["features"])
            if not same:
                row["report_fields_differ"].append("features")
            continue
        if refit.get(k) != release[k]:
            row["report_fields_differ"].append(k)
    row["unexpected_new_keys"] = sorted(set(refit) - set(release) - ALLOWED_NEW)
    fm = refit.get("fold_models")
    if isinstance(fm, dict):
        ran = [s for s in release["splits"] if s.get("learners")]
        ok = len(fm.get("splits", [])) == len(ran) and all(
            fs.get("scheme") == s["scheme"] and all(len(fs.get(l) or []) == s["folds"] for l in LEARNERS)
            for fs, s in zip(fm["splits"], ran))
        row["fold_models"] = "complete" if ok else "incomplete"
    learners = []
    for s, t in zip(refit.get("splits", []), release["splits"]):
        by = {x.get("learner"): x for x in s.get("learners", [])}
        for lt in t.get("learners", []):
            name = lt.get("learner")
            ls = by.get(name)
            if ls is None:
                learners.append((t["scheme"], name, "missing", 0, "NA"))
                continue
            keys = sorted((set(lt) | set(ls)) - {"seconds"})
            nd, worst = 0, 0.0
            for k in keys:
                if lt.get(k) != ls.get(k):
                    nd += 1
                    same, g = gap(ls.get(k), lt.get(k))
                    worst = max(worst, g if same else math.inf)
            learners.append((t["scheme"], name, "equal" if nd == 0 else "differs", nd,
                             "0" if nd == 0 else ("shape" if math.isinf(worst) else f"{worst:.3g}")))
    return row, learners


def main():
    refit_dir, release_dir, listed, out = sys.argv[1:5]
    os.makedirs(out, exist_ok=True)
    ds_rows, ln_rows, bad = [], [], 0
    for d in [x for x in listed.split(",") if x]:
        p = os.path.join(refit_dir, d + ".json")
        if not os.path.isfile(p):
            ds_rows.append([d, "no report", "", "", "", "", "", "", ""])
            bad += 1
            continue
        with open(p, encoding="utf-8") as f:
            refit_text = f.read()
        with open(os.path.join(release_dir, "datasets", d + ".json"), encoding="utf-8") as f:
            release_text = f.read()
        # Plain dicts for the field checks (key order aside), ordered ones for the text.
        refit, release = json.loads(refit_text), json.loads(release_text)
        text_equal = text_form(refit_text) == text_form(release_text)
        row, learners = check(refit, release)
        differ = row["report_fields_differ"]
        ds_rows.append([d, "equal" if not differ else "differs", ",".join(differ) or "-",
                        str(row["descriptive_samples_differ"]), row["fold_models"],
                        str(row["rule_groups_samples"]), ",".join(row["unexpected_new_keys"]) or "-",
                        "yes" if text_equal else "no", str(row["cards_with_other_locus"])])
        bad += bool(differ) or row["fold_models"] != "complete" or bool(row["unexpected_new_keys"]) or not text_equal
        for scheme, name, verdict, nd, worst in learners:
            ln_rows.append([d, scheme, name, verdict, str(nd), worst])
            bad += verdict != "equal"
    with open(os.path.join(out, "datasets.tsv"), "w", encoding="utf-8") as f:
        f.write("dataset\tfitted_fields\tfields_that_differ\tsamples_with_other_descriptive_fields"
                "\tfold_models\tsamples_with_rule_groups\tunexpected_new_keys\ttext_equal_less_allowed"
                "\tcards_with_other_locus\n")
        for r in ds_rows:
            f.write("\t".join(r) + "\n")
    with open(os.path.join(out, "learners.tsv"), "w", encoding="utf-8") as f:
        f.write("dataset\tsplit\tlearner\tverdict\tfields_that_differ\tlargest_abs_difference\n")
        for r in ln_rows:
            f.write("\t".join(r) + "\n")
    eq_d = sum(r[1] == "equal" for r in ds_rows)
    eq_l = sum(r[3] == "equal" for r in ln_rows)
    other = [int(r[8]) for r in ds_rows if r[8]]
    print(f"datasets: {len(ds_rows)}; fitted fields equal {eq_d}, not equal {len(ds_rows) - eq_d}; "
          f"fold_models complete {sum(r[4] == 'complete' for r in ds_rows)}; "
          f"text equal less the allowed fields {sum(r[7] == 'yes' for r in ds_rows)}; "
          f"datasets with a card whose kmer or locus differs {sum(n > 0 for n in other)}, such cards {sum(other)}")
    print(f"learners (dataset x split x learner): {len(ln_rows)}; equal {eq_l}, not equal {len(ln_rows) - eq_l}")
    for r in ln_rows:
        if r[3] != "equal":
            print(f"  differs: {r[0]} {r[1]} {r[2]}, {r[4]} fields, largest absolute difference {r[5]}")
    for r in ds_rows:
        if r[1] != "equal":
            print(f"  report fields differ: {r[0]}: {r[2]}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
