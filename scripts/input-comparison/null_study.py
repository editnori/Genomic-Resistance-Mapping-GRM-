#!/usr/bin/env python3
"""How often the cohort arms come out `improves` with no k-mer signal (rulings after lane INPUT).

`run` generates one synthetic data root per cohort seed (synthetic_root.py: the phenotype is
carriage of a blaKPC family K cannot see, so every arm that changes only the K model has no
signal to find), runs the decided comparison on it (run.py, spec.json) and keeps each seed's
summary.json and J contrasts; the data roots are deleted after each seed. `report` pools the
kept seeds over the decision-bearing contrasts only: the lineage split's (5 learners per seed),
as J decided them (summary.json `raw_decisions`, before K's sentinel mask; a summary written
before the mask existed has only `decisions`, which are then J's). Random-split contrasts are
descriptive since main's ruling after the calibration study and never enter a rate; a contrast
not produced is counted apart and printed. Per arm: the share of its lineage contrasts that are
`improves`, with a Clopper-Pearson 95% interval, and the share of seeds with any `improves`.
Contrasts within a seed are correlated (one cohort, one test set), so the per-contrast interval
is optimistic; the per-seed share is the conservative unit. Each seed's sentinel firings (K's
`unreliable` list) are reported apart (PR #69's review). Summaries written before the round-12
arm rename (`K+C`, `K+H`, `K+P`, `K+C+H+P`) are read under today's names.

Every arm also gets cluster-robust rates with the seed as the cluster: each seed's share of
`improves` (and of `worse`) over its lineage contrasts, averaged over seeds, with a 95% Student t
interval on the per-seed shares (seeds are independent). The null study passes only when every arm
it judges had a lineage contrast in at least 20 seeds (REVIEW-NULL16 L1), and refuses when one of
its judged arms is absent. Its figures before 7 Nov 2026 came from the cluster test and are
superseded (docs/input-comparison-runner.md).
A seed is the one its cohort.json records, whatever its folder is called, and counts once (REVIEW-INPUT17B L2): its
copies must be byte-identical, and a summary without a readable seed is refused. A summary recorded under two seeds
is refused (REVIEW-INPUT18 F1). `run` writes each seed's summary.json last, so a seed killed before it reruns.

`report --calibration` is the calibration H-1 must pass (main's ruling 6 after VERIFY-H1STATS), on
seeds run with `--spec calibration.json` and H-1's exact decision code: the (repeat, fold) unit t
at 3 lineage repeats. Its twins are exchangeable with the default, so every decision on one is a
false positive: the runner's two sentinels (the learner-seed twin and the retuned twin), judged
on their own p as H-1 judges them, and the cohort-path arm. Each (twin, learner) row counts its
lineage contrasts, one per seed (seeds are independent), with improves and worse apart, each with
a Clopper-Pearson 95% interval. A row whose every contrast is exactly zero is zero by construction
(the twin doesn't reach that learner) and is reported, not judged. Every other row passes with at
least 100 contrasts, improves + worse at most 5% and the interval's upper bound at most 11%; every
learner needs a judged row; and every planted arm (k-h, k-p, k-c-h-p) must be `improves` for every
learner on every seed, as the run decided it (BH over the run's lineage contrasts). A seed whose
run made no lineage decision (run.py's exit 3) keeps its summary, is listed and counts toward no
row; a seed counts once however many folders hold it (REVIEW-INPUT17 L2).

`report --power` is the power curve (pre-registered in docs/input-comparison-runner.md, round 2), on seeds
run with `--spec power.json --visible V` for each size V: per arm, learner and V, the share of seeds whose
lineage contrast is `improves` (Clopper-Pearson 95%) and the mean estimate, and per arm and learner the
minimum detectable effect, the smallest V detected on at least 80% of seeds (None: above every tested V).

usage: python3 -B null_study.py run --seeds 101-110 --out DIR --grm GRM [--spec SPEC] [--jobs N] [--parallel P]
                                    [--visible V]
       python3 -B null_study.py report DIR [DIR ...] [--json OUT] [--calibration | --power]
"""
import argparse
import json
import hashlib
import math
import os
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "scripts" / "ablation"))
import stats  # noqa: E402  J, for its decision rule
COHORT_ARMS = ("cap-1000", "cap-4000", "grr-v1")
PLANTED_ARMS = ("k-h", "k-p", "k-c-h-p")
POWER_ARMS = ("k-h", "k-p")
# The block arms' names before the round-12 rename; the stored null-study seeds (5 Nov 2026) use them.
RENAMED = {"K+C": "k-c", "K+H": "k-h", "K+P": "k-p", "K+C+H+P": "k-c-h-p"}
ALPHA = 0.05
MIN_SEEDS = 20
# The calibration (main's ruling 6 after VERIFY-H1STATS): its twins, by arm, and its bars.
TWINS = {"sentinel-seed-twin": "learner-seed twin", "sentinel-tuning-twin": "retuned twin",
         "cohort-path": "cohort path"}
LEARNERS = ("scm", "cart", "rf", "gbt", "lasso")
MIN_CONTRASTS = 100
MAX_RATE = 0.05
MAX_UPPER = 0.11
# The power curve's detection bar.
DETECTED = 0.8
NO_DECISIONS_EXIT = 3


def seeds(text):
    out = []
    for part in text.split(","):
        lo, _, hi = part.partition("-")
        out += list(range(int(lo), int(hi or lo) + 1))
    return out


def synthetic_spec(spec, out):
    """A copy of `spec` for synthetic roots: "synthetic": true, no revision-12 gate (it pins the real master
    root), and a 1 MiB kset size (the spec's 27.5 MiB is a real genome's) for the disk check."""
    s = json.loads(Path(spec).read_text())
    s["synthetic"] = True
    s.pop("revision12_gate", None)
    s.setdefault("pool_budget", {})["default_kset_bytes"] = 2 ** 20
    path = out / f"synthetic-{Path(spec).name}"
    path.write_text(json.dumps(s, indent=1))
    return path, s


def one_seed(seed, out, grm, jobs, spec=HERE / "spec.json", visible=1.0):
    work = out / f"seed-{seed}"
    if (work / "summary.json").exists():
        return seed, 0.0
    started = time.monotonic()
    raw, run = work / "raw", work / "work"
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    with open(work / "log.txt", "wb") as log:
        subprocess.run([sys.executable, "-B", str(HERE / "synthetic_root.py"), str(raw), "--seed", str(seed),
                        *(["--visible", str(visible)] if visible != 1 else [])],
                       check=True, stdout=log, stderr=subprocess.STDOUT)
        path, s = synthetic_spec(spec, work)
        vocabulary = (["--sigdrug-vocabulary", str(raw / "vocabulary")]
                      if s.get("curate_inputs", {}).get("sigdrug_vocabulary") else [])
        rc = subprocess.run([sys.executable, "-B", str(HERE / "run.py"), str(path), "--raw", str(raw),
                             "--work", str(run), "--grm", grm, "--source-commit", "0" * 40, "--jobs", str(jobs),
                             *vocabulary], stdout=log, stderr=subprocess.STDOUT).returncode
    # A run with a dataset that made no lineage decision exits 3 after writing its summary; the seed keeps it.
    if rc not in (0, NO_DECISIONS_EXIT):
        raise subprocess.CalledProcessError(rc, "run.py")
    # summary.json last and whole (REVIEW-INPUT18 L2): a rerun skips a folder with one, so a seed killed before it
    # landed reruns, and one that has it also has the cohort.json that names its seed.
    truth = json.loads((raw / "truth.json").read_text())
    (work / "cohort.json").write_text(json.dumps({"seed": seed, "genomes": len(truth["genomes"]),
                                                  "resistant": sum(g["resistant"] for g in truth["genomes"]),
                                                  "visible": truth["config"].get("visible", 1.0)}))
    shutil.copy(run / "j" / "contrasts.json", work / "contrasts.json")
    shutil.copy(run / "summary.json", work / "summary.json.part")
    os.replace(work / "summary.json.part", work / "summary.json")
    for d in (raw, run):
        shutil.rmtree(d)
    return seed, time.monotonic() - started


def binom_cdf(k, n, p):
    return sum(math.comb(n, i) * p ** i * (1 - p) ** (n - i) for i in range(k + 1))


def clopper_pearson(k, n, level=0.95):
    """Exact two-sided interval for a binomial share k/n."""
    a = (1 - level) / 2

    def solve(f):
        lo, hi = 0.0, 1.0
        for _ in range(100):
            mid = (lo + hi) / 2
            lo, hi = (mid, hi) if f(mid) else (lo, mid)
        return (lo + hi) / 2

    lower = 0.0 if k == 0 else solve(lambda p: 1 - binom_cdf(k - 1, n, p) < a)
    upper = 1.0 if k == n else solve(lambda p: binom_cdf(k, n, p) > a)
    return [lower, upper]


def t975(df):
    """The 0.975 quantile of Student's t with df degrees of freedom, by bisection on its CDF."""
    def cdf(t):  # integrate the density with Simpson's rule from 0
        c = math.exp(math.lgamma((df + 1) / 2) - math.lgamma(df / 2)) / math.sqrt(df * math.pi)
        n = 2000
        h = t / n
        f = [c * (1 + (i * h) ** 2 / df) ** (-(df + 1) / 2) for i in range(n + 1)]
        return 0.5 + h / 3 * (f[0] + f[-1] + 4 * sum(f[1:-1:2]) + 2 * sum(f[2:-1:2]))
    lo, hi = 0.0, 50.0
    for _ in range(80):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if cdf(mid) < 0.975 else (lo, mid)
    return (lo + hi) / 2


def per_seed(shares):
    """Mean of per-seed shares with a 95% t interval (seed as the cluster)."""
    s = len(shares)
    mean = sum(shares) / s
    if s < 2:
        return {"mean": mean, "ci95": None, "seeds": s}
    sd = math.sqrt(sum((x - mean) ** 2 for x in shares) / (s - 1))
    half = t975(s - 1) * sd / math.sqrt(s)
    return {"mean": mean, "ci95": [mean - half, mean + half], "seeds": s}


def lineage_decisions(arm):
    """One arm's decision-bearing contrasts in one seed's summary.json: the lineage split's, as J decided them
    (`raw_decisions` where present, so K's sentinel mask hides nothing), and how many were not produced."""
    raw = arm.get("raw_decisions", {})
    made, missing = [], 0
    for key, shown in arm["decisions"].items():
        if key.rsplit(":", 1)[-1] != "lineage":
            continue
        if shown.startswith("not produced"):
            missing += 1
        else:
            made.append(raw.get(key, shown))
    return made, missing


def seed_of(summary):
    """The seed of a seed folder's summary.json: the integer `seed` in the cohort.json beside it (`one_seed` writes
    both). A summary whose seed can't be read is refused."""
    cohort = summary.parent / "cohort.json"
    try:
        seed = json.loads(cohort.read_bytes())["seed"]
    except (OSError, ValueError, KeyError, TypeError) as e:
        raise SystemExit(f"{summary}: no seed in {cohort} ({type(e).__name__}: {e})")
    if type(seed) is not int:
        raise SystemExit(f"{summary}: the seed in {cohort} is {seed!r}, not an integer")
    return seed



def visible_of(seed_dir):
    """The visible share a seed's root was written with (cohort.json), 1 when it records none."""
    try:
        return json.loads((seed_dir / "cohort.json").read_text()).get("visible", 1.0)
    except OSError:
        return 1.0


def report(dirs, calibration=False, power=False):
    for d in dirs:
        if not Path(d).is_dir():
            raise SystemExit(f"{d}: no such folder")
    # One summary per seed (REVIEW-INPUT17 L2), keyed by the seed cohort.json records, not the folder's name
    # (REVIEW-INPUT17B L2): a folder given twice, a seed in two folders or a copy under another name counts once when
    # the summaries are byte-identical and is refused when they differ. One summary under two seeds is refused too
    # (REVIEW-INPUT18 F1): every run records its own `seconds`, so only a copy with an edited cohort.json matches.
    seen, digests = {}, {}
    for path in sorted(p for d in dirs for p in Path(d).glob("seed-*/summary.json")):
        seed = seed_of(path)
        # A seed is its cohort.json's seed and its root's visible share: the power curve runs one seed at
        # several shares.
        first = seen.setdefault((seed, visible_of(path.parent)), path)
        if first is not path and first.read_bytes() != path.read_bytes():
            raise SystemExit(f"seed {seed} has two different summaries: {first} and {path}")
        other, at = digests.setdefault(hashlib.sha256(path.read_bytes()).digest(), (seed, path))
        if other != seed:
            raise SystemExit(f"seeds {other} and {seed} have byte-identical summaries: {at} and {path}")
    summaries = sorted(seen.values())
    if calibration:
        return calibrate(summaries)
    if power:
        return power_curve(summaries)
    arms = {}
    firings = {"seeds_with_a_firing": 0, "firings": 0, "learners": {}}
    for path in summaries:
        summary = json.loads(path.read_text())
        fired = summary.get("unreliable", [])
        firings["seeds_with_a_firing"] += bool(fired)
        firings["firings"] += len(fired)
        for u in fired:
            firings["learners"][u["learner"]] = firings["learners"].get(u["learner"], 0) + 1
        for arm, s in summary["arms"].items():
            a = arms.setdefault(RENAMED.get(arm, arm), {"contrasts": 0, "improves": 0, "worse": 0, "seeds": 0,
                                      "seeds_any_improves": 0, "not_produced": 0,
                                      "_improves": [], "_worse": []})
            decisions, missing = lineage_decisions(s)
            a["not_produced"] += missing
            if not decisions:
                continue
            a["_improves"].append(decisions.count("improves") / len(decisions))
            a["_worse"].append(decisions.count("worse") / len(decisions))
            a["contrasts"] += len(decisions)
            a["improves"] += decisions.count("improves")
            a["worse"] += decisions.count("worse")
            a["seeds"] += 1
            a["seeds_any_improves"] += "improves" in decisions
    for name, a in arms.items():
        if not a["contrasts"]:
            raise SystemExit(f"arm {name}: no lineage contrast was produced in any seed")
        a["per_seed_improves"] = per_seed(a.pop("_improves"))
        a["per_seed_worse"] = per_seed(a.pop("_worse"))
        a["improves_rate"] = a["improves"] / a["contrasts"]
        a["improves_ci95"] = clopper_pearson(a["improves"], a["contrasts"])
        a["seed_any_improves_ci95"] = clopper_pearson(a["seeds_any_improves"], a["seeds"])
    absent = [n for n in COHORT_ARMS + PLANTED_ARMS if n not in arms]
    if absent:
        raise SystemExit(f"null study: arm {', '.join(absent)} is in no seed's summary.json "
                         f"({len(summaries)} summaries); the study judges {', '.join(COHORT_ARMS + PLANTED_ARMS)}")
    cohort = {k: sum(arms[a][k] for a in COHORT_ARMS) for k in ("contrasts", "improves", "worse", "seeds",
                                                                  "seeds_any_improves")}
    cohort["improves_rate"] = cohort["improves"] / cohort["contrasts"]
    cohort["improves_ci95"] = clopper_pearson(cohort["improves"], cohort["contrasts"])
    planted = {a: arms[a]["improves"] / arms[a]["contrasts"] for a in PLANTED_ARMS}
    return {"seeds": len(summaries), "nominal_alpha": ALPHA, "arms": arms, "cohort_arms_pooled": cohort,
            "planted_improves_share": planted, "split": "lineage", "sentinel_firings": firings,
            "min_seeds": MIN_SEEDS,
            "pass": (all(arms[a]["seeds"] >= MIN_SEEDS for a in COHORT_ARMS + PLANTED_ARMS)
                     and cohort["improves_rate"] <= ALPHA
                     and all(arms[a]["improves_rate"] <= ALPHA for a in COHORT_ARMS)
                     and all(v == 1.0 for v in planted.values()))}


def calibrate(summaries):
    """The calibration H-1 must pass (module docstring), from each seed's contrasts.json and summary.json."""
    rows = {(arm, learner): {"twin": TWINS[arm], "arm": arm, "learner": learner, "contrasts": 0, "improves": 0,
                             "worse": 0, "zero": 0, "not_produced": 0}
            for arm in TWINS for learner in LEARNERS}
    planted = {(arm, learner): {"arm": arm, "learner": learner, "seeds": 0, "improves": 0, "missed": []}
               for arm in PLANTED_ARMS for learner in LEARNERS}
    no_decision, absent = [], []
    for path in summaries:
        seed = path.parent.name
        summary = json.loads(path.read_text())
        if summary.get("no_decision_datasets"):
            no_decision.append(seed)
            continue
        try:
            records = json.loads((path.parent / "contrasts.json").read_text())
        except (OSError, ValueError) as e:
            raise SystemExit(f"{seed}: no readable contrasts.json beside its summary ({e.__class__.__name__})")
        lineage = {(x["treatment"], x["learner"]): x for x in records if x["split"] == "lineage" and x["family"] == "C"}
        for arm in list(TWINS) + list(PLANTED_ARMS):
            if arm not in summary["arms"]:
                absent.append(f"{seed}: {arm}")
        for (arm, learner), row in rows.items():
            x = lineage.get((arm, learner))
            if x is None:
                row["not_produced"] += arm in summary["arms"]
                continue
            if TWINS[arm] != "cohort path" and x["adjust"] != "none":
                raise SystemExit(f"{seed}: {x['id']} was adjusted ({x['adjust']}); a sentinel is judged on its own p")
            row["contrasts"] += 1
            row["zero"] += x["estimate"] == 0 and x["ci"] == [0.0, 0.0]
            # J's own rule on the contrast's own p, as H-1 judges a sentinel (the cohort path's record is BH's).
            decision = stats.decide(x["ci"], x["p"])
            if decision != "inconclusive":
                row[decision] += 1
        for (arm, learner), row in planted.items():
            x = lineage.get((arm, learner))
            row["seeds"] += 1
            if x is not None and x["decision"] == "improves":
                row["improves"] += 1
            else:
                row["missed"].append(seed)
    out_rows = {}
    for (arm, learner), row in rows.items():
        n, hits = row["contrasts"], row["improves"] + row["worse"]
        row["zero_by_construction"] = n > 0 and row["zero"] == n
        row["rate"] = hits / n if n else None
        row["ci95"] = clopper_pearson(hits, n) if n else None
        row["directions"] = {d: {"rate": row[d] / n if n else None, "ci95": clopper_pearson(row[d], n) if n else None}
                             for d in ("improves", "worse")}
        row["judged"] = n > 0 and not row["zero_by_construction"]
        row["pass"] = (not row["judged"] or n >= MIN_CONTRASTS and row["rate"] <= MAX_RATE
                       and row["ci95"][1] <= MAX_UPPER)
        out_rows[f"{arm}:{learner}"] = row
    judged = {r["learner"] for r in out_rows.values() if r["judged"]}
    unjudged = [learner for learner in LEARNERS if learner not in judged]
    failed = sorted(name for name, r in out_rows.items() if not r["pass"])
    missed = sorted(f"{r['arm']}:{r['learner']}" for r in planted.values() if r["missed"])
    seeds = len(summaries) - len(no_decision)
    return {"seeds": len(summaries), "seeds_with_decisions": seeds, "no_decision_seeds": no_decision,
            "absent_arms": absent, "split": "lineage", "test": "repeat_fold_corrected_t",
            "bars": {"min_contrasts": MIN_CONTRASTS, "max_rate": MAX_RATE, "max_upper": MAX_UPPER},
            "rows": out_rows, "planted": {f"{r['arm']}:{r['learner']}": r for r in planted.values()},
            "unjudged_learners": unjudged, "failed_rows": failed, "planted_missed": missed,
            "pass": bool(seeds) and not absent and not unjudged and not failed and not missed}


def power_curve(summaries):
    """The power curve (module docstring) from each seed's cohort.json, contrasts.json and summary.json."""
    cells = {}
    for path in summaries:
        seed = path.parent.name
        visible = json.loads((path.parent / "cohort.json").read_text()).get("visible", 1.0)
        summary = json.loads(path.read_text())
        records = json.loads((path.parent / "contrasts.json").read_text())
        lineage = {(x["treatment"], x["learner"]): x for x in records if x["split"] == "lineage" and x["family"] == "C"}
        for arm in POWER_ARMS:
            if arm not in summary["arms"]:
                raise SystemExit(f"{seed}: no arm {arm} in its summary.json")
            for learner in LEARNERS:
                c = cells.setdefault((arm, learner, visible), {"arm": arm, "learner": learner, "visible": visible,
                                                                 "seeds": 0, "improves": 0, "estimates": []})
                x = lineage.get((arm, learner))
                c["seeds"] += 1
                c["improves"] += x is not None and x["decision"] == "improves"
                if x is not None and x["estimate"] is not None:
                    c["estimates"].append(x["estimate"])
    out = []
    for c in sorted(cells.values(), key=lambda c: (c["arm"], LEARNERS.index(c["learner"]), c["visible"])):
        estimates = c.pop("estimates")
        out.append({**c, "share": c["improves"] / c["seeds"], "ci95": clopper_pearson(c["improves"], c["seeds"]),
                    "mean_estimate": sum(estimates) / len(estimates) if estimates else None})
    mde = {}
    for arm in POWER_ARMS:
        for learner in LEARNERS:
            rows = [c for c in out if c["arm"] == arm and c["learner"] == learner]
            hit = [c for c in rows if c["share"] >= DETECTED]
            mde[f"{arm}:{learner}"] = ({"visible": hit[0]["visible"], "share": hit[0]["share"], "ci95": hit[0]["ci95"],
                                        "mean_estimate": hit[0]["mean_estimate"]} if hit else None)
    return {"seeds": len(summaries), "sizes": sorted({c["visible"] for c in out}), "bar": DETECTED, "cells": out,
            "minimum_detectable": mde}


def print_power(result):
    for c in result["cells"]:
        lo, hi = c["ci95"]
        mean = "n/a" if c["mean_estimate"] is None else f"{c['mean_estimate']:.4f}"
        print(f"{c['arm']:<5} {c['learner']:<6} visible {c['visible']:<6} improves {c['improves']:>2}/{c['seeds']:<2} = "
              f"{c['share']:.2f} [{lo:.2f}, {hi:.2f}]  mean log-loss gain {mean}")
    for key, m in result["minimum_detectable"].items():
        print(f"minimum detectable {key}: " + (f"visible {m['visible']} ({m['share']:.2f}, mean gain "
                                               f"{m['mean_estimate']:.4f})" if m else
                                               f"above every tested size {result['sizes']}"))


def print_calibration(result):
    for name, r in result["rows"].items():
        if not r["contrasts"]:
            print(f"{name:<28} no contrasts (not produced {r['not_produced']})")
            continue
        lo, hi = r["ci95"]
        state = "zero by construction" if r["zero_by_construction"] else ("pass" if r["pass"] else "FAIL")
        print(f"{name:<28} {r['improves'] + r['worse']:>3}/{r['contrasts']:<3} = {r['rate']:.3f} [{lo:.3f}, {hi:.3f}]  "
              f"improves {r['improves']}, worse {r['worse']}; not produced {r['not_produced']}; {state}")
    for name, r in result["planted"].items():
        print(f"planted {name:<20} improves {r['improves']}/{r['seeds']}" +
              (f"; missed {', '.join(r['missed'])}" if r["missed"] else ""))
    why = ([f"no seed made a lineage decision"] * (not result["seeds_with_decisions"])
           + [f"absent {a}" for a in result["absent_arms"]]
           + [f"{l} has no judged row" for l in result["unjudged_learners"]]
           + [f"row {n} fails" for n in result["failed_rows"]]
           + [f"planted {n} missed a seed" for n in result["planted_missed"]])
    print(f"calibration: seeds {result['seeds']} ({result['seeds_with_decisions']} with decisions; no decision: "
          f"{', '.join(result['no_decision_seeds']) or 'none'}); pass {result['pass']}"
          + (f" ({'; '.join(why)})" if why else ""))


def main(argv):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = p.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--seeds", required=True)
    r.add_argument("--out", required=True)
    r.add_argument("--grm", required=True)
    r.add_argument("--jobs", type=int, default=3)
    r.add_argument("--parallel", type=int, default=1)
    r.add_argument("--spec", default=str(HERE / "spec.json"))
    r.add_argument("--visible", type=float, default=1.0, help="synthetic_root.py --visible, for the power curve")
    q = sub.add_parser("report")
    q.add_argument("dirs", nargs="+")
    q.add_argument("--json")
    q.add_argument("--calibration", action="store_true")
    q.add_argument("--power", action="store_true")
    a = p.parse_args(argv)
    if a.cmd == "run":
        out = Path(a.out)
        out.mkdir(parents=True, exist_ok=True)
        with ThreadPoolExecutor(max_workers=a.parallel) as pool:
            for seed, seconds in pool.map(lambda s: one_seed(s, out, a.grm, a.jobs, a.spec, a.visible),
                                          seeds(a.seeds)):
                print(f"seed {seed}: {seconds:.0f} s", flush=True)
        return
    if a.calibration and a.power:
        raise SystemExit("report takes --calibration or --power, not both")
    result = report(a.dirs, a.calibration, a.power)
    if a.json:
        Path(a.json).write_text(json.dumps(result, indent=1, sort_keys=True))
    if a.calibration:
        print_calibration(result)
        return
    if a.power:
        print_power(result)
        return
    fmt = lambda r: f"{r['mean']:.3f} [{r['ci95'][0]:.3f}, {r['ci95'][1]:.3f}]" if r["ci95"] else f"{r['mean']:.3f}"
    for name, x in sorted(result["arms"].items()):
        print(f"{name:<14} per-seed improves {fmt(x['per_seed_improves'])}  worse {fmt(x['per_seed_worse'])}  "
              f"(contrasts: improves {x['improves']}/{x['contrasts']}, worse {x['worse']}; not produced "
              f"{x['not_produced']}; seeds with a lineage contrast {x['seeds']}/{result['seeds']})")
    f = result["sentinel_firings"]
    print(f"lineage split only, as J decided; the sentinel fired {f['firings']} times in {f['seeds_with_a_firing']} of "
          f"{result['seeds']} seeds {f['learners'] or ''}")
    print(f"lineage contrasts not produced: {sum(x['not_produced'] for x in result['arms'].values())}; "
          f"pass needs every judged arm in at least {result['min_seeds']} seeds")
    for name in COHORT_ARMS + PLANTED_ARMS:
        x = result["arms"][name]
        lo, hi = x["improves_ci95"]
        print(f"{name:<9} improves {x['improves']:>3}/{x['contrasts']:<3} = {x['improves_rate']:.3f} "
              f"[{lo:.3f}, {hi:.3f}]  worse {x['worse']:>3}  seeds with any improves "
              f"{x['seeds_any_improves']}/{x['seeds']}")
    c = result["cohort_arms_pooled"]
    print(f"cohort arms pooled: improves {c['improves']}/{c['contrasts']} = {c['improves_rate']:.3f} "
          f"[{c['improves_ci95'][0]:.3f}, {c['improves_ci95'][1]:.3f}] vs nominal {ALPHA}; "
          f"seeds {result['seeds']}; pass {result['pass']}")


if __name__ == "__main__":
    main(sys.argv[1:])
