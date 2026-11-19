"""catalog.py INDEX.json [INTERVALS_DIR]

The best learner against the known-gene rule, split by how much the rule itself can do on the dataset. Reads a GRM
release's index.json only and prints one JSON document on stdout. INTERVALS_DIR is accepted so that every script of
this folder takes the same command line; this script does not read it.

Revisions understood: 12 and 14, read from the index's evaluation_revision. Any other revision stops the script
with a message and exit status 1. The known-gene rule is fitted on nothing, so revision 14 changes none of its rows.

The rule is not fitted, so its one score answers for both splits. On a dataset where it calls nothing resistant, or
everything, its balanced accuracy is exactly 0.5 and "the learner is ahead" says little. The bands below sort the
datasets by the rule's own balanced accuracy. That is a split on an outcome, so a rule that scored high by chance on a
small dataset is sorted into a high band and looks better than it is; read the bands as a description, not as a test.

  rule            datasets with a rule entry; those where it calls nothing resistant, everything resistant, or sits at
                  exactly 0.5
  bands           per band of the rule's balanced accuracy: datasets, learner ahead on each split, rule ahead or within
                  2 points on the lineage split, and the median lineage difference (best learner minus rule)
  at_chance_by    where the rule sits under 0.55, by species and by drug class, with the way it fails: it finds under
                  20% of the resistant genomes ("misses_resistance": sensitivity under 0.2, zero included), or it calls
                  over 80% of the susceptible genomes resistant ("calls_all_resistant": specificity under 0.2, zero
                  included), or neither
"""
import collections
import json
import statistics as st
import sys
from pathlib import Path

REVISIONS = {12, 14}

index = json.loads(Path(sys.argv[1]).read_text())
if index.get("evaluation_revision") not in REVISIONS:
    sys.exit(f"catalog.py: {sys.argv[1]} is evaluation revision {index.get('evaluation_revision')!r}; "
             f"this script understands revisions {', '.join(map(str, sorted(REVISIONS)))}")
ds = index["datasets"]
cat = [d for d in ds if d.get("catalog")]
ba = lambda d: d["catalog"]["balanced_accuracy"]
out = {"datasets": len(ds), "rule": {
    # A rule that calls nothing resistant has sensitivity 0 and specificity 1; one that calls everything resistant has
    # sensitivity 1 and specificity 0. Sensitivity 0 alone also holds for a rule that misses every resistant genome
    # and calls some susceptible ones resistant, which is not a constant call.
    "with_an_entry": len(cat),
    "calls_nothing_resistant": sum(d["catalog"]["sensitivity"] == 0 and d["catalog"]["specificity"] == 1 for d in cat),
    "calls_everything_resistant": sum(d["catalog"]["sensitivity"] == 1 and d["catalog"]["specificity"] == 0 for d in cat),
    "balanced_accuracy_exactly_half": sum(ba(d) == 0.5 for d in cat)}}


def band(g):
    gl = [d for d in g if d["best_lineage"] is not None]
    diff = [d["best_lineage"] - ba(d) for d in gl]
    return {"datasets": len(g), "learner_ahead_random": sum(d["best_random"] > ba(d) for d in g),
            "with_a_lineage_score": len(gl), "learner_ahead_lineage": sum(x > 0 for x in diff),
            "rule_ahead_or_within_2pt_lineage": sum(x <= 0.02 for x in diff),
            "median_lineage_difference": round(st.median(diff), 4) if diff else None,
            "median_random_difference": round(st.median(d["best_random"] - ba(d) for d in g), 4) if g else None}


out["bands"] = {name: band([d for d in cat if lo <= ba(d) < hi]) for name, lo, hi in (
    ("under 0.55", 0, 0.55), ("0.55 to 0.70", 0.55, 0.70), ("0.70 to 0.85", 0.70, 0.85), ("0.85 and up", 0.85, 1.01),
    ("0.55 and up", 0.55, 1.01), ("all", 0, 1.01))}


def fails(d):
    c = d["catalog"]
    return "misses_resistance" if c["sensitivity"] < 0.2 else "calls_all_resistant" if c["specificity"] < 0.2 else "neither"


low = [d for d in cat if ba(d) < 0.55]
out["at_chance_by"] = {"failure": dict(collections.Counter(fails(d) for d in low))}
for field in ("species", "drug_class"):
    c = collections.defaultdict(lambda: [0, 0, collections.Counter()])
    for d in cat:
        c[d[field]][1] += 1
        if ba(d) < 0.55:
            c[d[field]][0] += 1
            c[d[field]][2][fails(d)] += 1
    out["at_chance_by"][field] = [{"name": k, "under_0.55": v[0], "of": v[1], **dict(v[2])}
                                  for k, v in sorted(c.items(), key=lambda kv: (-kv[1][0], kv[0])) if v[0]]
json.dump(out, sys.stdout, indent=1)
