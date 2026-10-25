#!/usr/bin/env python3
"""Call changes at the 0.5 threshold from the backup to the audited revision.

For every report of the audited revision whose backup copy has the same
ordered samples and labels, rebuild the outer folds and compare, per
held-out fold, the null-model and constant 0.5 folds and the samples scored
exactly 0.5 in other folds. A call is resistant when score >= 0.5. A
schema-2 null fold makes no call; a schema-1 null fold scored 0.5, which is
a resistant call.

usage: python3 -I transitions.py OUT.json --revision N --backup DIR [--root DIR]
"""
import glob
import json
import os
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import audit  # noqa: E402

NO_CALL = "no call"


def fold_state(schema, held):
    """(call, value) of one held-out fold: NO_CALL for a schema-2 null fold,
    R or S and the score for a constant fold, (None, None) otherwise."""
    if schema >= 2 and audit.is_null_fold(schema, held):
        return NO_CALL, None
    values = set(held)
    if len(values) == 1 and None not in values:
        v = next(iter(values))
        return ("R" if v >= 0.5 else "S"), v
    return None, None


p = audit.options("Call changes between a run and its backup.")
p.add_argument("out")
args = p.parse_args()
audit.configure(args, need=("backup",))
fg = audit.FoldGen()
null_rows, tie_rows, skipped = [], [], []
pooled = []
for path in sorted(glob.glob(audit.ROOT + "/data/species/*/results/*.json")):
    new = json.load(open(path))
    if new["evaluation"]["revision"] != audit.REVISION:
        continue
    bpath = "%s/data/species/%s/results/%s" % (audit.BACKUP, new["species_slug"], os.path.basename(path))
    if not os.path.exists(bpath):
        skipped.append([new["id"], "no backup"])
        continue
    old = json.load(open(bpath))
    ids = [s["id"] for s in new["samples"]]
    truth = [bool(s["resistant"]) for s in new["samples"]]
    if ids != [s["id"] for s in old["samples"]] or truth != [bool(s["resistant"]) for s in old["samples"]]:
        skipped.append([new["id"], "cohort changed"])
        continue
    schemas = old.get("schema", 1), new.get("schema", 1)
    folds = audit.fold_partition(new, fg, truth)
    n = len(truth)
    for sn, so in zip(new["splits"], old["splits"]):
        assert sn["scheme"] == so["scheme"]
        if sn.get("unavailable_reason") or so.get("unavailable_reason"):
            continue
        f = folds[sn["scheme"]]
        old_by = {lr["learner"]: lr for lr in so["learners"]}
        for lr in sn["learners"]:
            name = lr["learner"]
            if name not in audit.LEARNERS:
                continue
            a, b = old_by[name]["scores"], lr["scores"]
            ma, mb = old_by[name]["metrics"], lr["metrics"]
            keys = ("tp", "fp", "tn", "fn", "vme", "me", "balanced_accuracy")
            pooled.append({"id": new["id"], "scheme": sn["scheme"], "learner": name,
                           "old": None if ma is None else {k: ma[k] for k in keys},
                           "new": None if mb is None else {k: mb[k] for k in keys}})
            for k in range(audit.fold_settings(new)[0]):
                idx = [i for i in range(n) if f[i] == k]
                ca, va = fold_state(schemas[0], [a[i] for i in idx])
                cb, vb = fold_state(schemas[1], [b[i] for i in idx])
                if ca == NO_CALL or cb == NO_CALL or (va is not None and abs(va - 0.5) < 1e-9) or vb == 0.5:
                    r = sum(truth[i] for i in idx)
                    null_rows.append({"id": new["id"], "scheme": sn["scheme"], "learner": name, "fold": k,
                                      "n": len(idx), "r": r, "s": len(idx) - r,
                                      "old_value": va, "new_value": vb, "old_call": ca, "new_call": cb,
                                      "new_complexity": lr["fold_complexity"][k],
                                      "old_complexity": old_by[name]["fold_complexity"][k]})
                    continue
                for i in idx:
                    if a[i] is None or b[i] is None:
                        continue
                    if b[i] == 0.5 or (abs(a[i] - 0.5) < 1e-9 and a[i] != b[i]):
                        tie_rows.append({"id": new["id"], "scheme": sn["scheme"], "learner": name, "fold": k,
                                         "truth": "R" if truth[i] else "S", "old": a[i], "new": b[i],
                                         "old_call": "R" if a[i] >= 0.5 else "S", "new_call": "R" if b[i] >= 0.5 else "S"})


def tally(rows, key):
    return dict(sorted(Counter(key(r) for r in rows).items()))


def flipped(old_call, new_call, field):
    return sum(r[field] for r in null_rows if r["old_call"] == old_call and r["new_call"] == new_call)


out = {
    "revision": audit.REVISION,
    "skipped": skipped,
    "null_fold_transitions": tally(null_rows, lambda r: "%s/%s %s->%s" % (r["scheme"], r["learner"], r["old_call"], r["new_call"])),
    "null_fold_samples_flipped": {
        "R_now_called_R": flipped("S", "R", "r"),
        "S_now_called_R": flipped("S", "R", "s"),
        "R_now_called_S": flipped("R", "S", "r"),
        "S_now_called_S": flipped("R", "S", "s"),
        "R_now_uncalled": sum(r["r"] for r in null_rows if r["old_call"] in ("R", "S") and r["new_call"] == NO_CALL),
        "S_now_uncalled": sum(r["s"] for r in null_rows if r["old_call"] in ("R", "S") and r["new_call"] == NO_CALL),
        "R_uncalled_before_R_now": flipped(NO_CALL, "R", "r"),
        "S_uncalled_before_R_now": flipped(NO_CALL, "R", "s"),
    },
    "null_datasets": sorted({r["id"] for r in null_rows}),
    "tie_transitions": tally(tie_rows, lambda r: "%s/%s truth %s: %s->%s" % (r["scheme"], r["learner"], r["truth"], r["old_call"], r["new_call"])),
    "tie_datasets": len({r["id"] for r in tie_rows}),
    "null_rows": null_rows,
    "tie_rows": tie_rows,
    "pooled_changes": [p for p in pooled if p["old"] != p["new"]],
}
json.dump(out, open(args.out, "w"), indent=1, sort_keys=True)
print(json.dumps({k: v for k, v in out.items() if k not in ("null_rows", "tie_rows", "pooled_changes")}, indent=1, sort_keys=True))
print("pooled learner rows changed:", len(out["pooled_changes"]), "of", len(pooled))
