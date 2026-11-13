"""index_counts.py INDEX.json INTERVALS_DIR

Can the counts behind every served error-rate interval be rebuilt from the index alone? For each row of each
intervals file (a learner on a split, or the known-gene rule) this rebuilds the two denominators and the two event
counts from index fields and compares them with the file's x and n. Prints one JSON document on stdout.

Revisions understood: 12 and 14, read from the index's evaluation_revision. Any other revision stops the script
with a message and exit status 1. The denominators follow the release's scoring rule:

  revision 12, every learner row, and revision 14, a row without the flag:
      very major error: n = n_resistant - no_call.resistant;           x = round(vme * n)
      major error:      n = (n - n_resistant) - no_call.susceptible;   x = round(me * n)
  revision 14, a learner row whose split entry carries null_fits_scored: true (docs/equal-scoring.md):
      very major error: n = n_resistant;                               x = round(vme * n)
      major error:      n = n - n_resistant;                           x = round(me * n)

From revision 14, no_call on a flagged row still counts the genomes of null-fit folds, but they have a score and a
call, so the row is scored on every genome. The known-gene rule's index entry has no no_call, so it counts as zero.

The output counts the rows that carry the flag, how many of them have a non-zero no_call, and where each
null_fits_scored key was seen ("split:true" is a split entry with the value true). A flag anywhere else, or in a
revision-12 index, is counted there and does not change the rule."""
import json
import sys
from pathlib import Path

REVISIONS = {12, 14}
EQUAL_SCORING = 14  # the first revision that scores a flagged row on every genome

index, idir = json.loads(Path(sys.argv[1]).read_text()), Path(sys.argv[2])
revision = index.get("evaluation_revision")
if revision not in REVISIONS:
    sys.exit(f"index_counts.py: {sys.argv[1]} is evaluation revision {revision!r}; "
             f"this script understands revisions {', '.join(map(str, sorted(REVISIONS)))}")
out = {"revision": revision,
       "rule": "flagged rows scored on every genome" if revision >= EQUAL_SCORING else "no_call genomes left out",
       "datasets": len(index["datasets"]), "files": 0, "rows": 0, "learner_rows": 0, "rule_rows": 0,
       "flagged_rows": 0, "flagged_with_no_call": 0, "unflagged_with_no_call": 0, "rates_compared": 0,
       "rates_null_in_file": 0, "mismatches": 0, "examples": [], "flag_keys_seen": {}}
for d in index["datasets"]:
    f = json.loads((idir / f"{d['id']}.json").read_text())
    out["files"] += 1
    for s in f["splits"]:
        for r in s["learners"]:
            flagged = False
            if r["learner"] == "catalog":
                e, nc = d.get("catalog"), {}
                out["rule_rows"] += 1
            else:
                row = next((x for x in d["learners"] if x["learner"] == r["learner"]), {})
                e = row.get(s["scheme"])
                nc = (e or {}).get("no_call") or {}
                for where, obj in (("learner", row), ("split", e or {})):
                    if "null_fits_scored" in obj:
                        k = f"{where}:{json.dumps(obj['null_fits_scored'])}"
                        out["flag_keys_seen"][k] = out["flag_keys_seen"].get(k, 0) + 1
                flagged = revision >= EQUAL_SCORING and (e or {}).get("null_fits_scored") is True
                out["learner_rows"] += 1
                has_nc = any(nc.get(k, 0) for k in ("resistant", "susceptible"))
                if flagged:
                    out["flagged_rows"] += 1
                    out["flagged_with_no_call"] += int(has_nc)
                else:
                    out["unflagged_with_no_call"] += int(has_nc)
            out["rows"] += 1
            for key, total, skipped in (("vme", d["n_resistant"], nc.get("resistant", 0)),
                                        ("me", d["n"] - d["n_resistant"], nc.get("susceptible", 0))):
                if r[key] is None:
                    out["rates_null_in_file"] += 1
                    continue
                out["rates_compared"] += 1
                n = total if flagged else total - skipped
                x = None if e is None or e.get(key) is None else round(e[key] * n)
                if (x, n) != (r[key]["x"], r[key]["n"]):
                    out["mismatches"] += 1
                    if len(out["examples"]) < 10:
                        out["examples"].append(f"{d['id']} {s['scheme']} {r['learner']} {key}: index {x}/{n}, "
                                               f"file {r[key]['x']}/{r[key]['n']}")
json.dump(out, sys.stdout, indent=1)
