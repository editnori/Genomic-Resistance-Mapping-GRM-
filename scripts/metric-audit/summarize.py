#!/usr/bin/env python3
"""Aggregate one audit.py run into the receipt's evidence tables.

Null stays null: a learner row without calls has null rates, a schema-2
null fold has null balanced accuracy, and value lists sort null first
instead of failing or counting it as zero.

usage: python3 -I summarize.py RUN.json OUT.json
"""
import json
import sys
from collections import Counter, defaultdict

run = json.load(open(sys.argv[1]))
recs = [r for r in run["records"] if "error" not in r]
# The run's own revision; runs from the revision-6 audit did not record it.
current = run["summary"].get("revision", 6)
out = {"run": sys.argv[1], "summary": run["summary"]}


def null_first(values):
    return sorted(values, key=lambda v: (v is not None, v if v is not None else 0))


def folds_of(r):
    return (r.get("eval_config") or {}).get("folds", 5)


by_rev = defaultdict(list)
for r in recs:
    by_rev[r["revision"]].append(r)

cov = {}
for rev, rs in sorted(by_rev.items()):
    species = Counter(r["species"] for r in rs)
    learners = Counter((x["scheme"], x["learner"]) for r in rs for x in r["learners"])
    rows = [x for r in rs for x in r["learners"] if x["learner"] != "catalog"]
    uncalled = Counter()
    for x in rows:
        if x.get("no_call"):
            uncalled["%s/%s" % (x["scheme"], x["learner"])] += x["no_call"]["resistant"] + x["no_call"]["susceptible"]
    cov[rev] = {
        "files": len(rs),
        "schemas": dict(Counter(str(r.get("schema", 1)) for r in rs)),
        "species": dict(sorted(species.items())),
        "learner_rows": sum(learners.values()),
        "learner_rows_by_scheme": dict(Counter(s for (s, _) in learners.elements())),
        "fold_ba_values": sum(folds_of(r) for r in rs for x in r["learners"] if x["learner"] != "catalog"),
        "fold_ba_null": sum(1 for r in rs for f in r["null_folds"] if f["fold_ba"] is None),
        "rows_with_no_call": sum(1 for x in rows if x.get("no_call") and x["no_call"]["resistant"] + x["no_call"]["susceptible"]),
        "rows_without_calls": sum(1 for x in rows if x.get("called") == 0),
        "no_call_isolates": {key: v for key, v in sorted(uncalled.items()) if v},
        "best_learner": dict(Counter(str((r.get("best") or {}).get("learner")) for r in rs)),
        "files_bit_exact": sum(1 for r in rs if r["metrics_bit_exact"]),
        "catalog_panels": sum(1 for r in rs for p in r["panels"]),
        "catalog_panels_from_hits": sum(1 for r in rs for p in r["panels"] if p.get("recomputed_from_hits")),
        "probe_panels": sum(len(r["probe_panels"]) for r in rs),
        "rf_cards_checked": sum(r.get("rf_cards_checked", 0) for r in rs),
        "mic_reports": sum(1 for r in rs if r["mic"]),
        "mic_categorical": sum(1 for r in rs if r["mic"] and r["mic"]["categorical"]),
        "mic_absent_with_enough_mics": [r["id"] for r in rs if r.get("mic_absent_with_enough_mics")],
        "lineage_unavailable": [[r["id"], u["reason"]] for r in rs for u in r.get("unavailable", [])],
        "issues": sum(len(r["issues"]) for r in rs),
        "counted_hits_match": sum(1 for r in rs if r["counted_hits_match_samples"]),
        "samples_missing_from_counted": sum(r["samples_missing_from_counted"] for r in rs),
        "labels_from_mic_total": sum(r["counts"]["labels_from_mic"] for r in rs),
        "written_after_launch": sum(1 for r in rs if r["written_after_launch"]),
        "identical_to_backup": sum(1 for r in rs if r["identical_to_backup"]),
        "log_done": sum(1 for r in rs if r["log_done"]),
        "manifest_revision": dict(Counter(str(r["manifest_revision"]) for r in rs)),
        "manifest_compared": sum(1 for r in rs if r.get("manifest_compared")),
        "published_revision": dict(Counter(str(r["published_revision"]) for r in rs)),
        "published_identical": sum(1 for r in rs if r["published_identical"]),
        "has_unsupported_field": sum(1 for r in rs if r.get("has_unsupported_field")),
    }
out["coverage"] = cov

# Null-model folds (all held-out scores null in schema 2, all exactly 0.5 in
# schema 1), fitted folds that call every isolate at exactly 0.5, and other
# constant folds.
nulls = {}
for rev, rs in sorted(by_rev.items()):
    c = Counter((f["scheme"], f["learner"]) for r in rs for f in r["null_folds"])
    ds = Counter((f["scheme"], f["learner"]) for r in rs for f in {(g["scheme"], g["learner"]): g for g in r["null_folds"]}.values())
    nulls[rev] = {
        "folds": {"%s/%s" % k: v for k, v in sorted(c.items())},
        "datasets": {"%s/%s" % k: v for k, v in sorted(ds.items())},
        "fold_ba_values": null_first({f["fold_ba"] for r in rs for f in r["null_folds"]}),
        "complexities": null_first({f["complexity"] for r in rs for f in r["null_folds"]}),
        "half_called": [[r["id"], f] for r in rs for f in r.get("half_called_folds", [])],
        "constant_non_half": [[r["id"], f] for r in rs for f in r["constant_non_half_folds"]],
        "half_ties": {k: sum(r["half_ties_in_non_null_folds"].get(k, 0) for r in rs)
                      for k in sorted({k for r in rs for k in r["half_ties_in_non_null_folds"]})},
        "half_tie_datasets": sum(1 for r in rs if r["half_ties_in_non_null_folds"]),
        "folds_with_class_below_5": sum(len(r["folds_with_class_below_5"]) for r in rs),
    }
out["null_models"] = nulls

# Null folds before and after for the same datasets of the run's revision.
pairs = [r for r in by_rev.get(current, []) if r.get("vs_backup")]
out["vs_backup"] = {
    "revision": current,
    "datasets": len(pairs),
    "same_ordered_ids": sum(1 for r in pairs if r["vs_backup"]["same_ordered_ids"]),
    "same_labels": sum(1 for r in pairs if r["vs_backup"]["same_labels"]),
    "changed_cohorts": [[r["id"], r["vs_backup"]["n_before"], r["vs_backup"]["n_after"], r["vs_backup"]["r_before"],
                         r["vs_backup"]["r_after"], r["vs_backup"]["ids_added"], r["vs_backup"]["ids_removed"]]
                        for r in pairs if not r["vs_backup"]["same_labels"]],
    "null_folds_after": sum(len(r["null_folds"]) for r in pairs),
    "null_folds_before": sum(r["vs_backup"].get("old_null_folds", 0) for r in pairs if "old_null_folds" in r["vs_backup"]),
    "learners_with_changed_scores": dict(Counter(x for r in pairs for x in r["vs_backup"].get("learners_with_changed_scores", []))),
    "unchanged_learners": dict(Counter(
        "%s/%s" % (s, l) for r in pairs if r["vs_backup"].get("same_ordered_ids")
        for s in ("random", "lineage") for l in ("scm", "cart", "rf", "gbt", "lasso", "catalog")
        if "%s/%s" % (s, l) not in r["vs_backup"].get("learners_with_changed_scores", [])
        and any(x["scheme"] == s and x["learner"] == l for x in r["learners"]))),
    "panel_changes": [[r["id"], r["vs_backup"]["old_panels"], [[p["id"], p["tp"], p["fp"], p["tn"], p["fn"]] for p in r["panels"]] + [[p["id"], p["tp"], p["fp"], p["tn"], p["fn"], p["no_call_r"], p["no_call_s"]] for p in r["probe_panels"]]]
                      for r in pairs if r["vs_backup"].get("old_panels") is not None and
                      sorted(map(str, [p[:5] for p in r["vs_backup"]["old_panels"]])) != sorted(map(str, [[p["id"], p["tp"], p["fp"], p["tn"], p["fn"]] for p in r["panels"] + r["probe_panels"]]))],
    "old_unsupported": sorted({u for r in pairs for u in r["vs_backup"].get("old_unsupported", [])}),
}

# Probe no-calls and unsupported assays.
probes = []
for r in recs:
    for p in r["probe_panels"]:
        probes.append([r["id"], r["revision"], p["id"], p["tp"], p["fp"], p["tn"], p["fn"], p["no_call_r"], p["no_call_s"], p["no_call_field"],
                       r["counts"]["resistant"], r["counts"]["susceptible"]])
out["probe_panels"] = probes
out["probe_totals"] = {
    "revision": current,
    "panels": sum(1 for p in probes if p[1] == current),
    "no_call_r": sum(p[7] for p in probes if p[1] == current),
    "no_call_s": sum(p[8] for p in probes if p[1] == current),
}
uns = [[r["id"], r["revision"], u["id"], u["reason"]] for r in recs for u in r["unsupported_assays"]]
out["unsupported"] = uns
out["unsupported_totals"] = dict(Counter("rev%s %s" % (u[1], u[2]) for u in uns))

# Kpn meropenem.
kpn = next((r for r in recs if r["id"] == "klebsiella_pneumoniae__meropenem"), None)
if kpn:
    out["kpn_meropenem"] = {k: kpn.get(k) for k in ("revision", "schema", "counts", "issues", "metrics_bit_exact", "identical_to_backup",
                                                    "written_after_launch", "log_done", "manifest_revision", "published_revision",
                                                    "published_identical", "vs_backup", "sha256", "mic", "panels", "best")}

json.dump(out, open(sys.argv[2], "w"), indent=1, sort_keys=True)
print(json.dumps({"coverage": cov, "null_models": nulls, "vs_backup": {k: v for k, v in out["vs_backup"].items() if k != "panel_changes"},
                  "panel_changes": len(out["vs_backup"]["panel_changes"]),
                  "probe_totals": out["probe_totals"], "unsupported_totals": out["unsupported_totals"]}, indent=1, sort_keys=True))
