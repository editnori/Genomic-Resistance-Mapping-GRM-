"""compare.py A/uncertainty.json B/uncertainty.json

The headline numbers of two releases side by side. Reads the two JSON documents uncertainty.py printed and prints one
line per statistic. Nothing is recomputed.

Revisions understood: 12 and 14, read from each document's "revision". Any other revision stops the script with a
message and exit status 1."""
import json
import sys
from pathlib import Path

REVISIONS = {12, 14}

a, b = (json.loads(Path(p).read_text()) for p in sys.argv[1:3])
for path, doc in zip(sys.argv[1:3], (a, b)):
    if doc.get("revision") not in REVISIONS:
        sys.exit(f"compare.py: {path} is revision {doc.get('revision')!r}; "
                 f"this script understands revisions {', '.join(map(str, sorted(REVISIONS)))}")


def cell(r):
    if r is None:
        return "-"
    ci = r.get("ci95_species")
    return f"{r['value']} n={r['n']}" + (f" [{ci[0]}, {ci[1]}]" if ci else "")


print(f"revisions {a['revision']} and {b['revision']}; datasets {a['datasets']} and {b['datasets']}; reps {a['reps']} and {b['reps']}")
print("checks:", a["checks"]["problem_count"], b["checks"]["problem_count"], "rows", a["checks"]["rows"], b["checks"]["rows"])
for block in ("medians", "gaps", "paired_lineage", "shares", "catalog"):
    print(f"== {block}")
    for k in a[block]:
        ra, rb = a[block][k], b[block].get(k)
        extra = ""
        if "first_at_least_second" in ra:
            extra = f" | first>=second {ra['first_at_least_second']} / {rb['first_at_least_second']}"
        print(f"{k:38s} {cell(ra):40s} | {cell(rb)}{extra}")
print("== selection")
for k in ("median_lineage_of_best", "median_lineage_of_oracle", "median_oracle_minus_best", "median_random_best_of_five"):
    print(f"{k:38s} {cell(a['selection'][k]):40s} | {cell(b['selection'][k])}")
print("best_is_also_lineage_best", a["selection"]["best_is_also_lineage_best"], "of", a["selection"]["of"], "|",
      b["selection"]["best_is_also_lineage_best"], "of", b["selection"]["of"])
print("median_random_single_learners", a["selection"]["median_random_single_learners"], "|", b["selection"]["median_random_single_learners"])
for l, ra in a["selection"]["median_best_minus_fixed_learner"].items():
    rb = b["selection"]["median_best_minus_fixed_learner"][l]
    print(f"best-{l:33s} {cell(ra):40s} | {cell(rb)} | best>=fixed {ra['best_at_least_fixed']} / {rb['best_at_least_fixed']}")
print("== error bars (n, point meets both, shown to meet, shown to miss, not decided)")
for k, ra in a["error_bars"]["rows"].items():
    rb = b["error_bars"]["rows"][k]
    f = lambda r: f"{r['n']:4d} {r['point_meets_both']:3d} {r['shown_to_meet_both']:3d} {r['shown_to_miss_one']:4d} {r['not_decided']:4d}"
    print(f"{k:18s} {f(ra)} | {f(rb)}")
print("resistant needed", a["error_bars"]["resistant_needed_for_zero_misses_to_show_vme"], b["error_bars"]["resistant_needed_for_zero_misses_to_show_vme"],
      "datasets with that many", a["error_bars"]["datasets_with_that_many_resistant"], b["error_bars"]["datasets_with_that_many_resistant"])
