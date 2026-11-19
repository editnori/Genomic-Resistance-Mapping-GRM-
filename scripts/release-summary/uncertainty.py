"""uncertainty.py INDEX.json INTERVALS_DIR [--reps N]

Uncertainty on the main numbers of a GRM release, from two public files only: the release's index.json and its
per-dataset intervals files (intervals/<id>.json: confusion counts and Clopper-Pearson 95% bounds on the very major
error VME and the major error ME). Nothing is fitted and no genome is read. Prints one JSON document on stdout.

Revisions understood: 12 and 14, read from the index's evaluation_revision. Any other revision stops the script
with a message and exit status 1. Nothing here depends on the scoring rule that changed in revision 14: the script
reads the rates and bounds the release serves.

1. checks     every intervals file has the index's revision, names its dataset, and reproduces the index's VME and ME
2. medians    median balanced accuracy per learner and split, with two bootstrap 95% intervals:
              "datasets" resamples datasets as if independent; "species" resamples whole species (23 clusters), since
              datasets of one species share genomes and lineages. The species interval is the one to quote.
3. paired     median paired difference in lineage balanced accuracy between the set covering machine and each other
              learner, same intervals, and the random-to-lineage gap per learner
4. selection  the site's "best learner" is chosen on the random split. Its lineage score against the best lineage
              score any of the five reaches (an oracle no one can use) and against each fixed learner
5. bars       per split and learner: datasets whose point estimates meet VME <= 1.5% and ME <= 3%, whose upper bounds
              meet both (shown to meet), whose lower bound breaks one (shown to miss), and the rest (not decided)
6. catalog    best learner against the known-gene rule, with the species interval

The bootstrap seed and the number of repetitions are in the output, so a run is deterministic for a given input.
The percentile bootstrap with 23 clusters is rough; read the species intervals to one decimal of a point at most."""
import json
import random
import statistics as st
import sys
from pathlib import Path

REVISIONS = {12, 14}

a = sys.argv[1:]
reps = 4000
if "--reps" in a:
    i = a.index("--reps")
    reps = int(a[i + 1])
    del a[i:i + 2]
index, idir = json.loads(Path(a[0]).read_text()), Path(a[1])
if index.get("evaluation_revision") not in REVISIONS:
    sys.exit(f"uncertainty.py: {a[0]} is evaluation revision {index.get('evaluation_revision')!r}; "
             f"this script understands revisions {', '.join(map(str, sorted(REVISIONS)))}")
ds = index["datasets"]
L = ["scm", "cart", "rf", "gbt", "lasso"]
VME_BAR, ME_BAR = 0.015, 0.03
out = {"revision": index["evaluation_revision"], "generated_at": index["generated_at"], "datasets": len(ds),
       "species": len({d["species_slug"] for d in ds}), "reps": reps, "seed": 20261111}


def entry(d, learner, split):
    e = next((x for x in d["learners"] if x["learner"] == learner), None)
    return None if e is None else e.get(split)


def ba(d, learner, split):
    e = entry(d, learner, split)
    return None if e is None else e.get("balanced_accuracy")


# ---- 1. checks
iv, problems = {}, []
for d in ds:
    f = json.loads((idir / f"{d['id']}.json").read_text())
    if f["dataset"] != d["id"] or f["evaluation_revision"] != index["evaluation_revision"] or f["method"] != "clopper-pearson":
        problems.append(f"{d['id']}: header")
    rows = {}
    for s in f["splits"]:
        for r in s["learners"]:
            rows[(s["scheme"], r["learner"])] = r
    iv[d["id"]] = rows
    for (scheme, learner), r in rows.items():
        # the known-gene rule is fitted on nothing, so its one index entry answers for both splits
        e = d.get("catalog") if learner == "catalog" else entry(d, learner, scheme)
        if e is None:
            problems.append(f"{d['id']}: {scheme} {learner} has an intervals row and no index entry")
            continue
        for k in ("vme", "me"):
            if r[k] is not None and e.get(k) is not None and abs(r[k]["x"] / r[k]["n"] - e[k]) > 1e-9:
                problems.append(f"{d['id']}: {scheme} {learner} {k} {r[k]['x']}/{r[k]['n']} is not the index's {e[k]}")
out["checks"] = {"intervals_files": len(iv), "schemes": sorted({k[0] for r in iv.values() for k in r}),
                 "rows": sum(len(r) for r in iv.values()), "problems": problems[:20], "problem_count": len(problems)}

# ---- the series each statistic is taken over: one value or None per dataset
series = {}
for l in L:
    for s in ("random", "lineage"):
        series[f"{l}:{s}"] = [ba(d, l, s) for d in ds]
    series[f"{l}:gap"] = [None if ba(d, l, "random") is None or ba(d, l, "lineage") is None
                          else ba(d, l, "random") - ba(d, l, "lineage") for d in ds]
series["best:random"] = [d["best_random"] for d in ds]
series["best:lineage"] = [d["best_lineage"] for d in ds]
series["best:gap"] = [None if d["best_random"] is None or d["best_lineage"] is None else d["best_random"] - d["best_lineage"] for d in ds]
for l in L[1:]:
    series[f"scm-{l}:lineage"] = [None if ba(d, "scm", "lineage") is None or ba(d, l, "lineage") is None
                                  else ba(d, "scm", "lineage") - ba(d, l, "lineage") for d in ds]
oracle = []
for d in ds:
    v = [ba(d, l, "lineage") for l in L]
    v = [x for x in v if x is not None]
    oracle.append(max(v) if v and d["best_lineage"] is not None else None)
series["oracle:lineage"] = oracle
series["oracle-best:lineage"] = [None if o is None else o - d["best_lineage"] for o, d in zip(oracle, ds)]
for l in L:
    series[f"best-{l}:lineage"] = [None if d["best_lineage"] is None or ba(d, l, "lineage") is None
                                   else d["best_lineage"] - ba(d, l, "lineage") for d in ds]
cat = lambda d: (d.get("catalog") or {}).get("balanced_accuracy")
for s, key in (("random", "best_random"), ("lineage", "best_lineage")):
    series[f"best-catalog:{s}"] = [None if cat(d) is None or d[key] is None else d[key] - cat(d) for d in ds]
shares = {  # share statistics: (series, predicate)
    "best_ahead_of_catalog:random": ("best-catalog:random", lambda x: x > 0),
    "best_ahead_of_catalog:lineage": ("best-catalog:lineage", lambda x: x > 0),
    "catalog_ahead_or_within_2pt:lineage": ("best-catalog:lineage", lambda x: x <= 0.02),
    "gap_over_5pt:best": ("best:gap", lambda x: x > 0.05),
}


def med(idx, v):
    x = sorted(v[i] for i in idx if v[i] is not None)
    n = len(x)
    return None if n == 0 else (x[n // 2] if n % 2 else (x[n // 2 - 1] + x[n // 2]) / 2)


def share(idx, v, pred):
    x = [v[i] for i in idx if v[i] is not None]
    return None if not x else sum(pred(y) for y in x) / len(x)


def stats(idx):
    r = {k: med(idx, v) for k, v in series.items()}
    r.update({k: share(idx, series[s], p) for k, (s, p) in shares.items()})
    return r


everything = list(range(len(ds)))
point = stats(everything)
by_species = {}
for i, d in enumerate(ds):
    by_species.setdefault(d["species_slug"], []).append(i)
clusters = list(by_species.values())
rng = random.Random(out["seed"])
boot = {"datasets": {k: [] for k in point}, "species": {k: [] for k in point}}
for _ in range(reps):
    for scheme, idx in (("datasets", rng.choices(everything, k=len(everything))),
                        ("species", [i for c in rng.choices(clusters, k=len(clusters)) for i in c])):
        for k, v in stats(idx).items():
            if v is not None:
                boot[scheme][k].append(v)


def ci(v):
    """The central 95% of the bootstrap values: the order statistics k and len - 1 - k, k = round(0.025 (len - 1)),
    so the same number of values lies outside each end. With 4000 values these are the 101st and the 3900th."""
    if not v:
        return None
    v = sorted(v)
    k = round(0.025 * (len(v) - 1))
    return [round(v[k], 4), round(v[len(v) - 1 - k], 4)]


res = {}
for k, p in point.items():
    n = sum(x is not None for x in series[shares[k][0] if k in shares else k])
    res[k] = {"n": n, "value": None if p is None else round(p, 4), "ci95_datasets": ci(boot["datasets"][k]),
              "ci95_species": ci(boot["species"][k])}
out["medians"] = {k: res[k] for k in res if k.split(":")[0] in L + ["best"] and k.split(":")[1] in ("random", "lineage")}
out["gaps"] = {k: res[k] for k in res if k.endswith(":gap")}
out["paired_lineage"] = {k: {**res[k], "first_at_least_second": sum(x is not None and x >= 0 for x in series[k])}
                         for k in res if k.startswith("scm-")}
out["selection"] = {
    "note": "best = the learner with the highest random-split score; oracle = the highest lineage score of the five",
    "median_lineage_of_best": res["best:lineage"], "median_lineage_of_oracle": res["oracle:lineage"],
    "median_oracle_minus_best": res["oracle-best:lineage"],
    "best_is_also_lineage_best": sum(o is not None and abs(o - d["best_lineage"]) < 1e-12 for o, d in zip(oracle, ds)),
    "of": sum(o is not None for o in oracle),
    "median_best_minus_fixed_learner": {l: {**res[f"best-{l}:lineage"],
                                            "best_at_least_fixed": sum(x is not None and x >= 0 for x in series[f"best-{l}:lineage"])}
                                        for l in L},
    "median_random_best_of_five": res["best:random"],
    "median_random_single_learners": {l: res[f"{l}:random"]["value"] for l in L},
}
out["shares"] = {k: res[k] for k in shares}
out["catalog"] = {k: res[k] for k in res if k.startswith("best-catalog")}

# ---- 5. error bars, with the served bounds
bars = {}
for split in ("random", "lineage"):
    for who in ["best"] + L + ["catalog"]:
        c = {"n": 0, "point_meets_both": 0, "shown_to_meet_both": 0, "shown_to_miss_one": 0, "not_decided": 0,
             "shown_to_miss_vme": 0, "shown_to_miss_me": 0, "point_meets_vme": 0, "point_meets_me": 0}
        rn = []
        for d in ds:
            r = iv[d["id"]].get((split, d["best_learner"] if who == "best" else who))
            if r is None or r["vme"] is None or r["me"] is None:
                continue
            v, m = r["vme"], r["me"]
            c["n"] += 1
            rn.append(v["n"])
            c["point_meets_vme"] += v["x"] / v["n"] <= VME_BAR
            c["point_meets_me"] += m["x"] / m["n"] <= ME_BAR
            c["point_meets_both"] += v["x"] / v["n"] <= VME_BAR and m["x"] / m["n"] <= ME_BAR
            c["shown_to_miss_vme"] += v["lo"] > VME_BAR
            c["shown_to_miss_me"] += m["lo"] > ME_BAR
            if v["hi"] <= VME_BAR and m["hi"] <= ME_BAR:
                c["shown_to_meet_both"] += 1
            elif v["lo"] > VME_BAR or m["lo"] > ME_BAR:
                c["shown_to_miss_one"] += 1
            else:
                c["not_decided"] += 1
        c["median_resistant_scored"] = st.median(rn) if rn else None
        bars[f"{who}:{split}"] = c
out["error_bars"] = {"vme_bar": VME_BAR, "me_bar": ME_BAR, "rows": bars}
# how many resistant genomes a dataset needs before zero misses can show VME <= 1.5% at 95%: 1 - 0.025 ** (1 / n) <= 0.015
n = 1
while 1 - 0.025 ** (1 / n) > VME_BAR:
    n += 1
out["error_bars"]["resistant_needed_for_zero_misses_to_show_vme"] = n
out["error_bars"]["datasets_with_that_many_resistant"] = sum(d["n_resistant"] >= n for d in ds)
json.dump(out, sys.stdout, indent=1)
