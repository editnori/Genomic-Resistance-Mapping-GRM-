#!/usr/bin/env python3
"""Stage 2's entry rule and arms (DECISIONS D3, amended by main's ruling 5 after VERIFY-H1STATS).

Stage 2 is exploratory. A stage-1 setting (one arm) enters stage 2 on a dataset when, on the lineage
split's (repeat, fold) unit t, at least one reliable learner (one K didn't mark unreliable on the
dataset's species) is `improves` with its log-loss interval's lower end above the spec's margin
(0.02) and its AUROC-delta interval above 0, and no reliable learner is `worse`. VERIFY-H1STATS
showed a log-loss `improves` can be a real change with no resistance signal (budget-1000's rf), so
the decision alone doesn't let a setting in.

The rest of D3 stands. Block arms that enter combine by union. When two values of one other input
both enter, the one with the larger mean log-loss gain (the mean estimate over the reliable
learners) wins; ties go to the value nearer the default, then the smaller budget. k settings join
only when the k add-on ran. Each stage-2 arm is the block union alone, or the union plus one
winning value of one other input; a stage-2 spec says "stage": 2, and run.py's load_spec takes
such two-factor arms (never k with a cohort change: that is two inputs besides the blocks).

usage: python3 -B stage2.py STAGE1_SPEC REPORT_JSON [--k-addon]   (prints each dataset's entries and arms)
"""
import argparse
import importlib.util
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("input_comparison_run", HERE / "run.py")
run = sys.modules.get(_spec.name)
if run is None:
    run = importlib.util.module_from_spec(_spec)
    sys.modules[_spec.name] = run
    _spec.loader.exec_module(run)

# The order D3's tie-break walks for the inputs other than the blocks.
INPUTS = ("cap", "selection", "budget", "k")


def entry(arm, dataset, contrasts, unreliable, species, margin):
    """(enters, why) for one stage-1 arm on one dataset, from K's report records (`contrasts`, `unreliable`) and the
    dataset -> species map."""
    flagged = {u["learner"] for u in unreliable if u["species"] == species[dataset]}
    lineage = [x for x in contrasts if x["dataset"] == dataset and x["treatment"] == arm
               and x["split"] == "lineage" and x["family"] == "C" and x["learner"] not in flagged]
    if not lineage:
        return False, "no reliable learner has a lineage decision"
    worse = sorted(x["learner"] for x in lineage if x["decision"] == "worse")
    if worse:
        return False, f"worse for {', '.join(worse)}"
    clears = sorted(x["learner"] for x in lineage if x["decision"] == "improves" and x["ci"][0] > margin
                    and x["auroc_ci"] is not None and x["auroc_ci"][0] > 0)
    if not clears:
        return False, (f"no reliable learner improves with its log-loss interval above {margin:g} and its AUROC-delta "
                       "interval above 0")
    return True, f"enters on {', '.join(clears)}"


def gain(arm, dataset, contrasts, unreliable, species):
    """The arm's mean log-loss gain on the dataset: its lineage estimates averaged over the reliable learners."""
    flagged = {u["learner"] for u in unreliable if u["species"] == species[dataset]}
    values = [x["estimate"] for x in contrasts if x["dataset"] == dataset and x["treatment"] == arm
              and x["split"] == "lineage" and x["family"] == "C" and x["learner"] not in flagged
              and x["estimate"] is not None]
    return sum(values) / len(values) if values else float("-inf")


def distance(factor, arm, default):
    """How far the arm's value of `factor` is from the default's, for D3's tie-break."""
    if factor == "cap":
        return abs(arm["cap"] - default["cap"])
    if factor == "budget":
        return abs(arm["screened_features"] - default["screened_features"])
    if factor == "k":
        return abs(arm["k"] - default["k"])
    return 1  # another selection rule: every one is one step from the default's


def winner(factor, candidates, default):
    """D3's pick among the entered arms changing one input: (gain, arm) pairs."""
    return min(candidates, key=lambda c: (-c[0], distance(factor, c[1], default),
                                          c[1]["screened_features"] if factor == "budget" else 0, c[1]["name"]))[1]


def stage2_arms(spec_path, report, k_addon=False):
    """Per confirmatory dataset: each stage-1 arm's entry and reason, and the stage-2 arms (run.py arm entries for
    a "stage": 2 spec) built from the arms that entered."""
    spec, default, arms = run.load_spec(spec_path)
    species = {d["id"]: d["species"] for d in report["manifest"]["datasets"]}
    margin = spec["margin"]["log_loss"]
    out = {}
    for d in spec["datasets"]:
        if d["partition"] != "confirmatory":
            continue
        ds = d["id"]
        entered, why = [], {}
        for arm in arms:
            if arm.get("factor") in (None, "calibration") or len(arm["factors"]) != 1:
                continue
            ok, reason = entry(arm["name"], ds, report["contrasts"], report["unreliable"], species, margin)
            why[arm["name"]] = reason
            if ok and (arm["factors"] != ["k"] or k_addon):
                entered.append(arm)
            elif ok:
                why[arm["name"]] += "; k joins only when the k add-on ran"
        blocks = sorted({b for a in entered if a["factors"] == ["blocks"] for b in a["blocks"]},
                        key="KCHPNX".index)
        union = {"blocks": blocks} if blocks else {}
        picks = []
        for factor in INPUTS:
            candidates = [(gain(a["name"], ds, report["contrasts"], report["unreliable"], species), a)
                          for a in entered if a["factors"] == [factor]]
            if candidates:
                picks.append(winner(factor, candidates, default))
        stage2 = []
        if union:
            stage2.append({"name": "s2-blocks", **union})
        for pick in picks:
            keys = [k for f in pick["factors"] for k in run.FACTORS[f] if k in pick and pick[k] != default.get(k)]
            stage2.append({"name": f"s2-{pick['name']}", **union, **{k: pick[k] for k in keys}})
        out[ds] = {"entries": why, "entered": [a["name"] for a in entered], "arms": stage2}
    return out


def main(argv):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("spec")
    p.add_argument("report")
    p.add_argument("--k-addon", action="store_true", help="the k add-on ran, so k settings may join")
    a = p.parse_args(argv)
    report = json.loads(Path(a.report).read_text(encoding="utf-8"))
    print(json.dumps(stage2_arms(a.spec, report, a.k_addon), indent=1, sort_keys=True))


if __name__ == "__main__":
    main(sys.argv[1:])
