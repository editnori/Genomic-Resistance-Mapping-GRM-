#!/usr/bin/env python3
"""Size of the null-fold effect in one revision's reports.

For each learner row of the audited revision with a null-model fold, put
the pooled metrics with those folds' isolates scored 0.5 (a resistant call,
as schema-1 reports stored them) beside the metrics with them left out (no
call, as schema-2 reports store them). A row whose every fold is null has
no metrics without them. Also checks that a zero fold_complexity for gbt or
lasso always means a null-model fold, and the reverse.

usage: python3 -I nullimpact.py OUT.json --revision N [--root DIR] [--reshaped]
"""
import glob
import json
import os
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import audit  # noqa: E402

p = audit.options("Null-fold effect on pooled metrics.")
p.add_argument("out")
args = p.parse_args()
audit.configure(args, need=())
fg = audit.FoldGen()
rows, mismatch = [], []
null_by = Counter()
for path in sorted(glob.glob(audit.ROOT + "/data/species/*/results/*.json")):
    rep = json.load(open(path))
    if rep["evaluation"]["revision"] != audit.REVISION:
        continue
    schema = rep.get("schema", 1)
    truth = [bool(s["resistant"]) for s in rep["samples"]]
    folds = audit.fold_partition(rep, fg, truth)
    k = audit.fold_settings(rep)[0]
    n = len(truth)
    for split in rep["splits"]:
        if split.get("unavailable_reason"):
            continue
        f = folds[split["scheme"]]
        for lr in split["learners"]:
            name = lr["learner"]
            if name not in audit.LEARNERS:
                continue
            null = []
            for fold in range(k):
                idx = [i for i in range(n) if f[i] == fold]
                is_null = audit.is_null_fold(schema, [lr["scores"][i] for i in idx])
                zero = lr["fold_complexity"][fold] == 0.0
                if name in ("gbt", "lasso") and is_null != zero:
                    mismatch.append([rep["id"], split["scheme"], name, fold, is_null, lr["fold_complexity"][fold]])
                if is_null:
                    null.append(fold)
                    null_by["%s/%s" % (split["scheme"], name)] += 1
            if not null:
                continue
            ranked = name in audit.RANKED
            half = audit.metrics(truth, [0.5 if f[i] in null else lr["scores"][i] for i in range(n)], ranked)
            keep = [i for i in range(n) if f[i] not in null]
            kept = audit.metrics([truth[i] for i in keep], [lr["scores"][i] for i in keep], ranked)

            def without(key):
                return None if kept is None else kept[key]

            rows.append({
                "id": rep["id"], "scheme": split["scheme"], "learner": name, "null_folds": null,
                "stored": "without_null" if schema >= 2 else "with_null",
                "n": half["n"], "n_without_null": without("n"),
                "vme": half["vme"], "vme_without_null": without("vme"),
                "me": half["me"], "me_without_null": without("me"),
                "ba": half["balanced_accuracy"], "ba_without_null": without("balanced_accuracy"),
                "auroc": half["auroc"], "auroc_without_null": without("auroc"),
            })
measured = [r for r in rows if r["vme_without_null"] is not None]
out = {"revision": audit.REVISION, "null_folds": dict(sorted(null_by.items())), "rows": rows, "complexity_mismatch": mismatch,
       "rows_without_any_call": sum(1 for r in rows if r["n_without_null"] is None),
       "rows_vme_drop_ge_5pts": sum(1 for r in measured if r["vme_without_null"] - r["vme"] >= 0.05),
       "max_vme_understatement": max((r["vme_without_null"] - r["vme"] for r in measured), default=None)}
json.dump(out, open(args.out, "w"), indent=1, sort_keys=True)
print(json.dumps({k: v for k, v in out.items() if k != "rows"}, indent=1, sort_keys=True))


def shown(v):
    return None if v is None else round(v, 3)


for r in sorted(measured, key=lambda r: r["vme"] - r["vme_without_null"])[:8]:
    print("%-48s %-7s %-5s folds=%s n=%d->%d VME %.3f->%.3f ME %.3f->%.3f AUROC %s->%s" % (
        r["id"], r["scheme"], r["learner"], r["null_folds"], r["n"], r["n_without_null"],
        r["vme"], r["vme_without_null"], r["me"], r["me_without_null"],
        shown(r["auroc"]), shown(r["auroc_without_null"])))
