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
`unreliable` list) are reported apart (PR #69's review).

Every arm also gets cluster-robust rates with the seed as the cluster: each seed's share of
`improves` (and of `worse`) over its lineage contrasts, averaged over seeds, with a 95% Student t
interval on the per-seed shares (seeds are independent). `--spec calibration.json` runs the
calibration arms, which are exchangeable with the default: there every `improves` or `worse`
is a false positive, and the procedure passes when each arm's mean per-seed rate of each is at
most 5% (main's ruling on the null study, 5 Nov 2026). Either study passes only when every arm
it judges had a lineage contrast in at least 20 seeds (REVIEW-NULL16 L1): a seed that gave an
arm none doesn't count toward that arm, however many summaries there are.

usage: python3 -B null_study.py run --seeds 101-110 --out DIR --grm GRM [--spec SPEC] [--jobs N] [--parallel P]
       python3 -B null_study.py report DIR [DIR ...] [--json OUT] [--calibration]
"""
import argparse
import json
import math
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
COHORT_ARMS = ("cap-1000", "cap-4000", "grr-v1")
PLANTED_ARMS = ("k-h", "k-p", "k-c-h-p")
ALPHA = 0.05
MIN_SEEDS = 20


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


def one_seed(seed, out, grm, jobs, spec=HERE / "spec.json"):
    work = out / f"seed-{seed}"
    if (work / "summary.json").exists():
        return seed, 0.0
    started = time.monotonic()
    raw, run = work / "raw", work / "work"
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    with open(work / "log.txt", "wb") as log:
        subprocess.run([sys.executable, "-B", str(HERE / "synthetic_root.py"), str(raw), "--seed", str(seed)],
                       check=True, stdout=log, stderr=subprocess.STDOUT)
        path, s = synthetic_spec(spec, work)
        vocabulary = (["--sigdrug-vocabulary", str(raw / "vocabulary")]
                      if s.get("curate_inputs", {}).get("sigdrug_vocabulary") else [])
        subprocess.run([sys.executable, "-B", str(HERE / "run.py"), str(path), "--raw", str(raw),
                        "--work", str(run), "--grm", grm, "--source-commit", "0" * 40, "--jobs", str(jobs),
                        *vocabulary], check=True, stdout=log, stderr=subprocess.STDOUT)
    shutil.copy(run / "summary.json", work / "summary.json")
    shutil.copy(run / "j" / "contrasts.json", work / "contrasts.json")
    truth = json.loads((raw / "truth.json").read_text())
    (work / "cohort.json").write_text(json.dumps({"seed": seed, "genomes": len(truth["genomes"]),
                                                  "resistant": sum(g["resistant"] for g in truth["genomes"])}))
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


def report(dirs, calibration=False):
    summaries = sorted(p for d in dirs for p in Path(d).glob("seed-*/summary.json"))
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
            a = arms.setdefault(arm, {"contrasts": 0, "improves": 0, "worse": 0, "seeds": 0,
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
    if calibration:
        ok = all(a["per_seed_improves"]["mean"] <= ALPHA and a["per_seed_worse"]["mean"] <= ALPHA
                 for a in arms.values())
        return {"seeds": len(summaries), "nominal_alpha": ALPHA, "cluster": "seed", "arms": arms,
                "split": "lineage", "sentinel_firings": firings, "min_seeds": MIN_SEEDS,
                "pass": all(a["seeds"] >= MIN_SEEDS for a in arms.values()) and ok}
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
    q = sub.add_parser("report")
    q.add_argument("dirs", nargs="+")
    q.add_argument("--json")
    q.add_argument("--calibration", action="store_true")
    a = p.parse_args(argv)
    if a.cmd == "run":
        out = Path(a.out)
        out.mkdir(parents=True, exist_ok=True)
        with ThreadPoolExecutor(max_workers=a.parallel) as pool:
            for seed, seconds in pool.map(lambda s: one_seed(s, out, a.grm, a.jobs, a.spec), seeds(a.seeds)):
                print(f"seed {seed}: {seconds:.0f} s", flush=True)
        return
    result = report(a.dirs, a.calibration)
    if a.json:
        Path(a.json).write_text(json.dumps(result, indent=1, sort_keys=True))
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
    if a.calibration:
        print(f"calibration: seeds {result['seeds']}; pass {result['pass']}")
        return
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
