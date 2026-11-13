"""summary.py INDEX.json [INTERVALS_DIR]

Paper-facing summary of a GRM release's index.json (a public file). Every number is recomputed from the index rows;
nothing is fitted. Prints one JSON document on stdout. INTERVALS_DIR is accepted so that every script of this folder
takes the same command line; this script does not read it.

Revisions understood: 12 and 14, read from the index's evaluation_revision. Any other revision stops the script
with a message and exit status 1. The script reads the balanced accuracies, VME and ME the index serves, so the
revision-14 scoring rule (docs/equal-scoring.md) reaches it only through those served numbers."""
import json, statistics as st, sys
from pathlib import Path

REVISIONS = {12, 14}
m = json.loads(Path(sys.argv[1]).read_text())
if m.get("evaluation_revision") not in REVISIONS:
    sys.exit(f"summary.py: {sys.argv[1]} is evaluation revision {m.get('evaluation_revision')!r}; "
             f"this script understands revisions {', '.join(map(str, sorted(REVISIONS)))}")
ds = m["datasets"]
L = ["scm", "cart", "rf", "gbt", "lasso"]
def ba(x): return None if x is None else x.get("balanced_accuracy")
def med(v):
    v = [x for x in v if x is not None]
    return (round(st.median(v), 4), len(v)) if v else (None, 0)
out = {"revision": m["evaluation_revision"], "generated_at": m["generated_at"], "datasets": len(ds),
       "species": len({d["species_slug"] for d in ds}), "genome_rows": sum(d["n"] for d in ds)}
# per learner medians on both splits, over datasets where the learner has both
per = {}
for l in L:
    r, g, gap = [], [], []
    for d in ds:
        e = next((x for x in d["learners"] if x["learner"] == l), None)
        if not e: continue
        a, b = ba(e.get("random")), ba(e.get("lineage"))
        if a is not None and b is not None:
            r.append(a); g.append(b); gap.append(a - b)
    per[l] = {"n_both": len(r), "median_random": med(r)[0], "median_lineage": med(g)[0], "median_gap": med(gap)[0],
              "share_gap_over_5pt": round(sum(x > 0.05 for x in gap) / len(gap), 3) if gap else None}
out["per_learner"] = per
# headline (best learner as the index names it)
br = [d["best_random"] for d in ds if d["best_random"] is not None]
bl = [d["best_lineage"] for d in ds if d["best_lineage"] is not None]
both = [(d["best_random"], d["best_lineage"]) for d in ds if d["best_random"] is not None and d["best_lineage"] is not None]
out["best"] = {"median_random": med(br), "median_lineage": med(bl), "n_both": len(both),
               "median_gap": round(st.median(a - b for a, b in both), 4),
               "drop_over_5pt": sum(a - b > 0.05 for a, b in both), "drop_over_10pt": sum(a - b > 0.10 for a, b in both),
               "lineage_at_least_random": sum(b >= a for a, b in both)}
out["best_learner_counts"] = {l: sum(d["best_learner"] == l for d in ds) for l in L}
# catalog comparison on each split: best learner vs known-gene rule
cmp = {}
for split, key in (("random", "best_random"), ("lineage", "best_lineage")):
    rows = [(d[key], ba(d["catalog"])) for d in ds if d.get("catalog") and ba(d["catalog"]) is not None and d[key] is not None]
    diff = [a - c for a, c in rows]
    cmp[split] = {"n": len(rows), "learner_ahead": sum(x > 0 for x in diff), "catalog_ahead": sum(x < 0 for x in diff),
                  "catalog_within_2pt_or_ahead": sum(x <= 0.02 for x in diff), "learner_ahead_over_2pt": sum(x > 0.02 for x in diff),
                  "median_diff": round(st.median(diff), 4)}
out["catalog_vs_best"] = cmp
out["catalog"] = {"with_rule": sum(1 for d in ds if d.get("catalog") and ba(d["catalog"]) is not None),
                  "median_ba": med([ba(d["catalog"]) for d in ds if d.get("catalog")])[0]}
# error-rate bars on the lineage split for the best learner: VME <= 1.5%, ME <= 3%
def best_entry(d, split):
    e = next((x for x in d["learners"] if x["learner"] == d["best_learner"]), None)
    return None if not e else e.get(split)
bars = {}
for split in ("random", "lineage"):
    es = [best_entry(d, split) for d in ds]
    es = [e for e in es if e and e.get("vme") is not None and e.get("me") is not None]
    bars[split] = {"n": len(es), "vme_le_1.5": sum(e["vme"] <= 0.015 for e in es), "me_le_3": sum(e["me"] <= 0.03 for e in es),
                   "both": sum(e["vme"] <= 0.015 and e["me"] <= 0.03 for e in es),
                   "median_vme": round(st.median(e["vme"] for e in es), 4), "median_me": round(st.median(e["me"] for e in es), 4)}
out["error_bars_best_learner"] = bars
# per species: datasets, median best random/lineage, gap
sp = {}
for d in ds: sp.setdefault(d["species"], []).append(d)
rows = []
for s, v in sp.items():
    b = [(d["best_random"], d["best_lineage"]) for d in v if d["best_random"] is not None and d["best_lineage"] is not None]
    rows.append({"species": s, "datasets": len(v), "with_lineage": len(b),
                 "median_random": med([a for a, _ in b])[0], "median_lineage": med([x for _, x in b])[0],
                 "median_gap": round(st.median(a - x for a, x in b), 4) if b else None})
# a species with no lineage score sorts last; a gap of 0.0 sorts by its value
out["per_species"] = sorted(rows, key=lambda r: (r["median_gap"] is None, -(r["median_gap"] or 0)))
# per drug class
dc = {}
for d in ds: dc.setdefault(d.get("drug_class") or "unclassified", []).append(d)
rows = []
for c, v in dc.items():
    b = [(d["best_random"], d["best_lineage"]) for d in v if d["best_random"] is not None and d["best_lineage"] is not None]
    if len(b) >= 5:
        rows.append({"class": c, "datasets": len(v), "median_random": med([a for a, _ in b])[0],
                     "median_lineage": med([x for _, x in b])[0], "median_gap": round(st.median(a - x for a, x in b), 4)})
out["per_drug_class_min5"] = sorted(rows, key=lambda r: -r["median_gap"])
# size effect: lineage BA by cohort size tercile
ns = sorted(d["n"] for d in ds); t1, t2 = ns[len(ns)//3], ns[2*len(ns)//3]
def tier(n): return "small" if n < t1 else ("mid" if n < t2 else "large")
tiers = {}
for d in ds:
    if d["best_lineage"] is not None: tiers.setdefault(tier(d["n"]), []).append(d["best_lineage"])
out["lineage_by_size"] = {"cuts": [t1, t2], **{k: med(v) for k, v in tiers.items()}}
out["mic"] = {"with_model": sum(d.get("mic_essential_agreement") is not None for d in ds),
              "median_ea": med([d.get("mic_essential_agreement") for d in ds])[0]}
out["no_lineage"] = [d["id"] for d in ds if d["best_lineage"] is None]
json.dump(out, sys.stdout, indent=1)
