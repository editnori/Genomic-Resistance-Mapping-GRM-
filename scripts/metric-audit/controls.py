#!/usr/bin/env python3
"""Negative controls: show audit.py flags wrong folds, scores, and counts.

Each control writes a mutated copy of a report of the audited revision
under controls/ beside OUT.json and audits it. A control
passes when the expected issue appears; the unmutated copy must give no
issue. The schema-2 controls run when the source report is schema 2.

usage: python3 -I controls.py SOURCE_REPORT.json OUT.json --revision N [--root DIR]
           (--backup DIR --log FILE --launch UTC | --reshaped)
"""
import copy
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import audit  # noqa: E402

p = audit.options("Negative controls for audit.py.")
p.add_argument("source")
p.add_argument("out")
args = p.parse_args()
audit.configure(args)
base = json.load(open(args.source))
work = os.path.join(os.path.dirname(os.path.abspath(args.out)), "controls")
os.makedirs(work, exist_ok=True)
fg = audit.FoldGen()
log_done, log_recomputed = audit.done_datasets()
manifest = audit.load_manifest(audit.ROOT)


def run(name, rep, expect):
    path = "%s/%s.json" % (work, name)
    json.dump(rep, open(path, "w"))
    rec = audit.audit_file(path, fg, {}, log_done, log_recomputed, manifest)
    text = json.dumps(rec["issues"])
    hit = all(e in text for e in expect) if expect else rec["issues"] == []
    os.remove(path)
    return {"control": name, "expect": expect, "passed": bool(hit),
            "issues": [i.get("where", "") + " " + str(i.get("what", i.get("field", ""))) for i in rec["issues"]][:8]}


results = []
results.append(run("unmutated", copy.deepcopy(base), []))

rep = copy.deepcopy(base)
lr = rep["splits"][0]["learners"][2]  # random/rf
i = next(j for j, s in enumerate(rep["samples"]) if s["resistant"] and lr["scores"][j] is not None and lr["scores"][j] >= 0.5)
lr["scores"][i] = 0.1
results.append(run("one_rf_score_flipped", rep, ["random/rf/pooled", "random/rf/fold_ba"]))

rep = copy.deepcopy(base)
saved_seed = audit.SEED
if rep.get("eval_config"):
    rep["eval_config"]["seed"] += 1  # folds are rebuilt from the recorded seed
else:
    audit.SEED = saved_seed + 1
results.append(run("wrong_fold_seed", rep, ["/fold_ba"]))
audit.SEED = saved_seed

rep = copy.deepcopy(base)
rep["splits"][1]["learners"][3]["fold_balanced_accuracy"].reverse()  # lineage/gbt
results.append(run("lineage_fold_order_reversed", rep, ["lineage/gbt/fold_ba"]))

rep = copy.deepcopy(base)
rep["counts"]["labels_from_mic"] = 1
rep["samples"][0]["from_mic"] = True
results.append(run("inferred_label", rep, ["inferred labels present"]))

rep = copy.deepcopy(base)
rep["samples"][0]["resistant"] = not rep["samples"][0]["resistant"]
results.append(run("label_flipped", rep, ["counts.resistant", "/pooled"]))

rep = copy.deepcopy(base)
r_total = sum(1 for s in rep["samples"] if s["resistant"])
s_total = len(rep["samples"]) - r_total
good = {"id": "gc_gyra", "name": "fixture", "mirrors": "fixture", "targets": [],
        "metrics": dict(audit.from_counts(r_total - 7, 3, s_total - 5, 4), auroc=None, auprc=None),
        "detected": [["mutation detected", r_total - 7, 3]],
        "no_call": {"resistant": 3, "susceptible": 2}}
rep["panels"].append(good)
results.append(run("probe_consistent", rep, []))

rep = copy.deepcopy(base)
bad = copy.deepcopy(good)
bad["no_call"]["resistant"] = 2  # one resistant sample unaccounted for
rep["panels"].append(bad)
results.append(run("probe_no_call_short", rep, ["probe R partition"]))

rep = copy.deepcopy(base)
rep["panels"].append(copy.deepcopy(good))
rep["unsupported_assays"] = [{"id": "gc_gyra", "name": "f", "mirrors": "f", "reason": "fixture"}]
results.append(run("unsupported_also_scored", rep, ["unsupported assay also scored"]))

if base.get("mic"):
    rep = copy.deepcopy(base)
    # An exact MIC within one dilution: a 3-dilution shift must break EA.
    j = next(k for k, v in enumerate(rep["mic"]["predicted"])
             if v is not None and rep["samples"][k]["mic"]["censor"] == "Exact"
             and abs(audit.rust_round(v) - rep["samples"][k]["mic"]["log2"]) <= 1)
    rep["mic"]["predicted"][j] += 3.0
    results.append(run("mic_prediction_shifted", rep, ["mic/"]))

if base["panels"]:
    rep = copy.deepcopy(base)
    rep["panels"][0]["metrics"]["tp"] += 1
    rep["panels"][0]["metrics"]["fn"] -= 1
    results.append(run("catalog_panel_count_moved", rep, ["panel/" + base["panels"][0]["id"]]))

if not audit.RESHAPED:
    rep = copy.deepcopy(base)
    rep["evaluation"]["revision"] = audit.REVISION - 1
    results.append(run("earlier_revision_written_after_launch", rep, ["non-revision-%d file written after launch" % audit.REVISION]))

if base.get("schema", 1) >= 2:
    rep = copy.deepcopy(base)
    lr = rep["splits"][0]["learners"][3]  # random/gbt
    lr["no_call"]["resistant"] += 1
    results.append(run("no_call_one_short", rep, ["random/gbt/no_call"]))

    rep = copy.deepcopy(base)
    folds = rep["splits"][0]["test_fold"]
    a = folds.index(0)
    b = folds.index(1)
    folds[a], folds[b] = folds[b], folds[a]
    results.append(run("test_fold_swapped", rep, ["random/test_fold"]))

    if base.get("eval_config"):
        rep = copy.deepcopy(base)
        rep["eval_config"]["forest"]["seed"] = 17
        results.append(run("forest_seed_recorded", rep, ["eval_config.forest fields"]))

json.dump(results, open(args.out, "w"), indent=1)
for r in results:
    print("%-38s %s  %s" % (r["control"], "PASS" if r["passed"] else "FAIL", "; ".join(r["issues"][:3])))
