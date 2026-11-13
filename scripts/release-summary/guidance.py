"""guidance.py INDEX.json INTERVALS_DIR

Reads a GRM release's error rates by the acceptance wording of FDA's 2009 guidance for antimicrobial susceptibility
test systems, next to the release's own stricter screen. Reads two public files only: the release's index.json and its
per-dataset intervals files (intervals/<id>.json: counts and Clopper-Pearson 95% bounds on the very major error VME and
the major error ME). Nothing is fitted and no genome is read. Prints one JSON document on stdout.

Revisions understood: 12 and 14, read from the index's evaluation_revision. Any other revision stops the script
with a message and exit status 1, and so does an intervals file whose method is not clopper-pearson at 95%. The script
reads the counts and bounds the release serves, so the revision-14 scoring rule (docs/equal-scoring.md) reaches it only
through those served numbers.

The guidance ("Class II Special Controls Guidance Document: Antimicrobial Susceptibility Test (AST) Systems", issued
August 28, 2009) lists as acceptable performance "A maj rate of < 3% based on the number of susceptible organisms
tested" and a vmj rate with "an upper 95% confidence limit for the true vmj rate of < 7.5% and the lower 95% confidence
limit for the true vmj < 1.5%". Here a row meets that wording when all three parts hold, each strict as written:

    vme_upper   the VME interval's upper bound hi < 0.075
    vme_lower   the VME interval's lower bound lo < 0.015
    me          the ME rate x / n < 0.03

The guidance does not apply to devices that detect genomic resistance features, nor to anti-mycobacterial agents. This
is a reading of the release by that wording: nothing here is validated, and nothing passes or fails a regulatory test.

For each split, the best learner (each dataset's best_learner), the five learners and the known-gene rule, over the rows
uncertainty.py reads for its error_bars (an intervals row with both rates):

  n                     datasets read
  meets_wording         the three parts together; vme_upper, vme_lower and me count each part alone, vme the two VME
                        parts together
  without_tuberculosis  n and meets_wording without the Mycobacterium tuberculosis datasets, which test
                        anti-mycobacterial drugs
  screen                the release's own screen as uncertainty.py counts it (VME <= 1.5% and ME <= 3%, by the
                        bounds): shown to meet both, shown to miss one, not decided
  vme_n_under_table_8   rows whose VME counts fewer resistant genomes than the first row of the guidance's Table 8
  fewer_than_every_genome
                        rows whose VME or ME counts fewer genomes than the dataset's class holds (the release's no-call
                        rule; the guidance counts every organism tested)

Table 8's first row is computed, not typed: the smallest number of resistant genomes for which zero misses give an
exact upper bound under 7.5%, 1 - 0.025 ** (1 / n) < 0.075. For the best learner on each split, best_meets_wording
lists the datasets that meet the wording, with the learner, the VME count and interval and the ME count and rate. In that
list the bounds are the served values and the rate is x / n, unrounded, so that a value just under a limit does not
print as the limit. Every other number of the output is a count or a limit."""
import json
import sys
from pathlib import Path

REVISIONS = {12, 14}
VME_UPPER, VME_LOWER, ME_RATE = 0.075, 0.015, 0.03  # the guidance's wording, each strict
VME_BAR, ME_BAR = 0.015, 0.03  # the release's own screen, on the bounds, at or under
TUBERCULOSIS = "mycobacterium_tuberculosis"

index, idir = json.loads(Path(sys.argv[1]).read_text()), Path(sys.argv[2])
if index.get("evaluation_revision") not in REVISIONS:
    sys.exit(f"guidance.py: {sys.argv[1]} is evaluation revision {index.get('evaluation_revision')!r}; "
             f"this script understands revisions {', '.join(map(str, sorted(REVISIONS)))}")
ds = index["datasets"]
L = ["scm", "cart", "rf", "gbt", "lasso"]

# the smallest number of resistant genomes at which zero misses give an exact 95% upper bound under 7.5%
table_8 = 1
while 1 - 0.025 ** (1 / table_8) >= VME_UPPER:
    table_8 += 1

out = {"revision": index["evaluation_revision"], "generated_at": index["generated_at"], "datasets": len(ds),
       "criteria": {
           "wording": {"vme_upper_below": VME_UPPER, "vme_lower_below": VME_LOWER, "me_rate_below": ME_RATE,
                       "strict": True, "vme_interval": "the served Clopper-Pearson 95% interval of VME",
                       "me_rate": "x / n of the served ME counts"},
           "screen": {"vme_upper_at_most": VME_BAR, "me_upper_at_most": ME_BAR,
                      "shown_to_miss": "VME lower bound over 1.5% or ME lower bound over 3%"},
           "without_tuberculosis": TUBERCULOSIS,
           "table_8_first_resistant": table_8}}

iv = {}
for d in ds:
    f = json.loads((idir / f"{d['id']}.json").read_text())
    if f.get("method") != "clopper-pearson" or f.get("confidence", 0.95) != 0.95:
        sys.exit(f"guidance.py: {d['id']}: the intervals are {f.get('method')!r} at {f.get('confidence')!r}; "
                 f"this script reads clopper-pearson 95% intervals only")
    iv[d["id"]] = {(s["scheme"], r["learner"]): r for s in f["splits"] for r in s["learners"]}


def wording(v, m):
    """The three parts of the guidance's wording on one row: (vme_upper, vme_lower, me)."""
    return v["hi"] < VME_UPPER, v["lo"] < VME_LOWER, m["x"] / m["n"] < ME_RATE


rows, listed = {}, {}
for split in ("random", "lineage"):
    for who in ["best"] + L + ["catalog"]:
        c = {"n": 0, "meets_wording": 0, "vme_upper": 0, "vme_lower": 0, "vme": 0, "me": 0,
             "without_tuberculosis": {"n": 0, "meets_wording": 0},
             "screen": {"shown_to_meet_both": 0, "shown_to_miss_one": 0, "not_decided": 0},
             "vme_n_under_table_8": 0, "fewer_than_every_genome": 0}
        for d in ds:
            learner = d["best_learner"] if who == "best" else who
            r = iv[d["id"]].get((split, learner))
            if r is None or r["vme"] is None or r["me"] is None:
                continue
            v, m = r["vme"], r["me"]
            upper, lower, me = wording(v, m)
            met = upper and lower and me
            c["n"] += 1
            c["vme_upper"] += upper
            c["vme_lower"] += lower
            c["vme"] += upper and lower
            c["me"] += me
            c["meets_wording"] += met
            if d["species_slug"] != TUBERCULOSIS:
                c["without_tuberculosis"]["n"] += 1
                c["without_tuberculosis"]["meets_wording"] += met
            if v["hi"] <= VME_BAR and m["hi"] <= ME_BAR:
                c["screen"]["shown_to_meet_both"] += 1
            elif v["lo"] > VME_BAR or m["lo"] > ME_BAR:
                c["screen"]["shown_to_miss_one"] += 1
            else:
                c["screen"]["not_decided"] += 1
            c["vme_n_under_table_8"] += v["n"] < table_8
            c["fewer_than_every_genome"] += v["n"] < d["n_resistant"] or m["n"] < d["n"] - d["n_resistant"]
            if who == "best" and met:
                listed.setdefault(split, []).append({
                    "dataset": d["id"], "learner": learner,
                    "vme": {"x": v["x"], "n": v["n"], "lo": v["lo"], "hi": v["hi"]},
                    "me": {"x": m["x"], "n": m["n"], "rate": m["x"] / m["n"]}})
        rows[f"{who}:{split}"] = c
out["rows"] = rows
out["best_meets_wording"] = {split: listed.get(split, []) for split in ("random", "lineage")}
json.dump(out, sys.stdout, indent=1)
