#!/usr/bin/env python3
"""A second reading of the input comparison's reports, one report per species, printed beside K's.

K's report decides inside one species: BH over that run's lineage contrasts, and rule 10, which marks every decision
of a learner on the species unreliable when one twin contrast of that learner rejects or isn't produced. This script
reads the reports of several species together and changes none of them. Nothing is fitted or tested again. Every
number is a count, a median, a binomial tail or a BH adjustment over the records J wrote
(docs/input-comparison-runner.md, "A second reading across species").

REPORT is a run's or a combine's `report/ablation-report.json`, loaded through K's validation (report.load_report).
Each REPORT holds one species. A species given twice is refused, and so are reports that differ in source commit,
grm binary or the synthetic flag.

A contrast is default -> arm for one dataset and learner, and an estimate above zero means the arm's log loss is
lower. The contrasts that decide are family C: the lineage split of the confirmatory datasets. A twin is a sentinel
arm, the default again with other seeds, so every true difference on a twin is zero.

Per species:
  arms         Per arm, the family C contrasts planned and made, and improves and worse as the report decided them.
               inside counts the made contrasts whose 90% interval lies wholly inside +/- the margin, beyond those
               whose estimate is at least the margin either way, and zero those whose estimate is exactly 0. Then
               the quartiles of the estimate and the median AUROC change. The twins are in the table, so each arm
               can be read against them.
  by dataset   Per arm, the datasets where the mean estimate over the learners is above zero (better) and below it
               (worse), with a two-sided exact sign test. A species' datasets share genomes, so this p guides and
               doesn't decide. The JSON holds each mean (dataset_means).
  calibration  Per learner, k rejections out of n contrasts of the twins that judge it: the family C contrasts
               among those rule 10 reads (report.judging_contrasts). A twin contrast rejects when its own p is
               under 0.05 and its 95% interval is off zero. One with no p is left out when the default has no
               lineage result for that dataset and learner, because the dataset then has no decision to protect.
               With the default's result present it counts as a rejection. flagged = P(k or more | n, 0.05) < 0.05.
               With no twin contrast to count, flagged is null and no decision of that learner is firm. The number
               of contrasts rule 10 marks unreliable prints beside the count.

Across the species (the panel):
  family       BH over every report's family C as one family, as stats.finalize_sets adjusts several sets: every
               planned contrast outside the ids the report judged on their own p, one without a p entering at
               p = 1.
  decisions    Every contrast that is improves or worse inside its species or in the panel family. firm = the same
               decision under both, on a learner the species' calibration count leaves unflagged.
  calibration  The calibration count per learner over all the species, with its exact 95% interval.
  by dataset   The sign test per arm over the datasets of all the species.
  predictions  Each --predict ARM=better|worse, read on the --predict-on species (every species without it): the
               datasets on the predicted side of zero against the other side, ties dropped, in a one-sided exact
               sign test at 0.05. pass needs every prediction.

usage: python3 -B second_reading.py REPORT [REPORT ...] [--json OUT] [--predict ARM=better|worse ...]
                                    [--predict-on SPECIES ...]
"""
import argparse
import importlib.util
import json
import statistics
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "ablation"))
import stats  # noqa: E402  J, for its BH and its decision rule
sys.path.insert(0, str(ROOT))
from scripts.ablation import report as k_report  # noqa: E402  K
from scripts.ablation import schema as k_schema  # noqa: E402

_spec = importlib.util.spec_from_file_location("input_comparison_null_study", HERE / "null_study.py")
null_study = sys.modules.get(_spec.name)
if null_study is None:
    null_study = importlib.util.module_from_spec(_spec)
    sys.modules[_spec.name] = null_study
    _spec.loader.exec_module(null_study)

# What must be the same in every report for them to be one comparison.
SHARED = ("source_commit", "binary_sha256", "synthetic")
SIDES = ("better", "worse")


def upper_tail(k, n, p):
    """P(X >= k) for X binomial(n, p), as the lower tail of n - X: a small tail keeps its digits, where 1 - cdf
    would lose them."""
    return min(1.0, null_study.binom_cdf(n - k, n, 1 - p))


def sign_test(values):
    """The values above zero (better) and below it (worse), with the two-sided exact sign test over those two."""
    better, worse = sum(v > 0 for v in values), sum(v < 0 for v in values)
    n = better + worse
    return {"datasets": len(values), "better": better, "worse": worse, "unchanged": len(values) - n,
            "median": statistics.median(values) if values else None,
            "p": min(1.0, 2 * upper_tail(max(better, worse), n, 0.5)) if n else None}


def calibration_count(k, n):
    """k rejections of n twin contrasts against a test of exact size ALPHA."""
    tail = upper_tail(k, n, stats.ALPHA) if n else None
    return {"rejections": k, "of": n, "tail": tail, "flagged": None if tail is None else tail < stats.ALPHA}


def parts_of(contrast_id):
    return k_schema.parse_key(k_schema.contrast_key(contrast_id))


def arm_row(planned, made, twin):
    estimates = sorted(x["estimate"] for x in made)
    q1, median, q3 = ((None,) * 3 if not made else (estimates[0],) * 3 if len(made) == 1
                      else statistics.quantiles(estimates, n=4, method="inclusive"))
    auroc = [x["auroc_delta"] for x in made if x["auroc_delta"] is not None]
    return {"twin": twin, "planned": planned, "made": len(made),
            "improves": sum(x["decision"] == "improves" for x in made),
            "worse": sum(x["decision"] == "worse" for x in made),
            "inside": sum(-x["margin"] < x["ci_tost"][0] and x["ci_tost"][1] < x["margin"] for x in made),
            "beyond": sum(abs(x["estimate"]) >= x["margin"] for x in made),
            "zero": sum(x["estimate"] == 0 for x in made),
            "q1": q1, "median": median, "q3": q3, "auroc_median": statistics.median(auroc) if auroc else None}


def read_species(report):
    """One species' reading from its loaded report: the module docstring's "Per species" block."""
    manifest = report["manifest"]
    names = sorted({d["species"] for d in manifest["datasets"]})
    if len(names) != 1:
        raise SystemExit(f"a report must hold one species, and this one holds {', '.join(names)}")
    sentinels = set(manifest.get("sentinel_arms", []))
    judging = {pid for pid, _ in k_report.judging_contrasts(manifest, report["planned_contrasts"])}
    ran = {(u["arm"], u["dataset"], u["learner"], u["split"]) for u in manifest["runs"] if u["status"] == "ok"}
    records = {x["id"]: x for x in report["contrasts"]}
    planned, made, counts = defaultdict(int), defaultdict(list), {}
    for pid in report["planned_contrasts"]:
        parts = parts_of(pid)
        if parts["family"] != "C":
            continue
        arm, learner, x = parts["treatment"], parts["learner"], records.get(pid)
        planned[arm] += 1
        if x is not None and x["estimate"] is not None:
            made[arm].append(x)
        count = counts.setdefault(learner, {"twins": set(), "k": 0, "n": 0, "left_out": 0, "no_test": 0})
        if pid not in judging:
            continue
        count["twins"].add(arm)
        if x is not None and x["p"] is not None:
            count["n"] += 1
            count["k"] += stats.decide(x["ci"], x["p"]) != "inconclusive"
        elif (parts["baseline"], k_schema.contrast_dataset(pid), learner, parts["split"]) in ran:
            count["n"] += 1
            count["k"] += 1
            count["no_test"] += 1
        else:
            count["left_out"] += 1
    means = {arm: {d: statistics.fmean(x["estimate"] for x in made[arm] if x["dataset"] == d)
                   for d in sorted({x["dataset"] for x in made[arm]})} for arm in planned}
    rule10 = {u["learner"]: len(u["contrasts"]) for u in report["unreliable"]}
    return {
        "species": names[0], "datasets": len(manifest["datasets"]),
        "unanswered": [u["dataset"] for u in k_report.unanswered_datasets(report)],
        "planned": sum(planned.values()), "made": sum(len(v) for v in made.values()),
        "arms": {arm: arm_row(planned[arm], made[arm], arm in sentinels) for arm in planned},
        "dataset_means": means,
        "by_dataset": {arm: sign_test(list(means[arm].values())) for arm in planned},
        "calibration": {learner: {"twins": sorted(c["twins"]), **calibration_count(c["k"], c["n"]),
                                  "left_out": c["left_out"], "no_test": c["no_test"],
                                  "rule10_contrasts": rule10.get(learner, 0)}
                        for learner, c in counts.items()},
    }


def bh_family(report):
    """The report's BH family as stats.finalize takes it: the planned family C contrasts outside the ids judged on
    their own p. Returns its records with a p, and how many planned contrasts have none and enter at p = 1."""
    outside = k_report.unadjusted_ids(report["manifest"], report["planned_contrasts"])
    family = {pid for pid in report["planned_contrasts"] if pid not in outside and parts_of(pid)["family"] == "C"}
    tested = [x for x in report["contrasts"] if x["id"] in family and x["p"] is not None]
    return tested, len(family) - len(tested)


def read_predictions(readings, predictions, predict_on):
    """Each (arm, side) prediction as a one-sided exact sign test over the datasets of the `predict_on` species."""
    on = [r for r in readings if not predict_on or r["species"] in predict_on]
    tests = []
    for arm, side in predictions:
        if not any(arm in r["arms"] for r in on):
            raise SystemExit(f"--predict names arm {arm}, which no report of {', '.join(r['species'] for r in on)} "
                             "plans")
        s = sign_test([v for r in on for v in r["dataset_means"].get(arm, {}).values()])
        n = s["better"] + s["worse"]
        tests.append({"arm": arm, "predicted": side, "better": s["better"], "worse": s["worse"],
                      "unchanged": s["unchanged"], "p_one_sided": upper_tail(s[side], n, 0.5) if n else None})
    return {"species": [r["species"] for r in on], "tests": tests,
            "pass": all(t["p_one_sided"] is not None and t["p_one_sided"] < stats.ALPHA for t in tests)}


def read_panel(reports, readings, predictions=(), predict_on=()):
    """What spans the species: the module docstring's "Across the species" block. `readings` holds read_species
    of each report, in the reports' order."""
    families = [bh_family(report) for report in reports]
    tested = [(reading, x) for reading, (records, _) in zip(readings, families) for x in records]
    untested = sum(n for _, n in families)
    adjusted = stats.bh([x["p"] for _, x in tested] + [1.0] * untested) if tested else []
    decisions = []
    for (reading, x), p_adj in zip(tested, adjusted):
        in_panel = stats.decide(x["ci"], p_adj)
        if x["decision"] == in_panel == "inconclusive":
            continue
        flagged = reading["calibration"][x["learner"]]["flagged"]
        decisions.append({"id": x["id"], "species": reading["species"], "in_species": x["decision"],
                          "in_panel": in_panel, "learner_flagged": flagged,
                          "firm": x["decision"] == in_panel and flagged is False, "estimate": x["estimate"],
                          "ci": x["ci"], "p": x["p"], "p_adj_species": x["p_adj"], "p_adj_panel": p_adj,
                          "auroc_delta": x["auroc_delta"]})
    calibration = {}
    for learner in sorted({learner for r in readings for learner in r["calibration"]}):
        rows = [r["calibration"][learner] for r in readings if learner in r["calibration"]]
        k, n = sum(c["rejections"] for c in rows), sum(c["of"] for c in rows)
        calibration[learner] = {**calibration_count(k, n), "ci95": null_study.clopper_pearson(k, n) if n else None}
    arms = list(dict.fromkeys(arm for r in readings for arm in r["arms"]))
    return {
        "family": {"tested": len(tested), "entered_at_p_1": untested, "size": len(tested) + untested},
        "decisions": decisions,
        "firm": {side: sum(d["firm"] and d["in_species"] == side for d in decisions)
                 for side in ("improves", "worse")},
        "one_family_only": sum(d["in_species"] != d["in_panel"] for d in decisions),
        "flagged_or_unjudged": sum(d["learner_flagged"] is not False for d in decisions),
        "calibration": calibration,
        "by_dataset": {arm: sign_test([v for r in readings for v in r["dataset_means"].get(arm, {}).values()])
                       for arm in arms},
        "predictions": read_predictions(readings, predictions, predict_on) if predictions else None,
    }


def read(reports, predictions=(), predict_on=()):
    """The whole reading of loaded reports: {"species": [one per report], "panel": {...}}."""
    readings = [read_species(report) for report in reports]
    names = [r["species"] for r in readings]
    twice = sorted({name for name in names if names.count(name) > 1})
    if twice:
        raise SystemExit(f"a species is read from one report, and these have more: {', '.join(twice)}")
    for field in SHARED:
        values = sorted({str(report["manifest"].get(field)) for report in reports})
        if len(values) > 1:
            raise SystemExit(f"the reports differ in {field} ({', '.join(values)}), so they aren't one comparison")
    unknown = sorted(set(predict_on) - set(names))
    if unknown:
        raise SystemExit(f"--predict-on names a species with no report: {', '.join(unknown)} (the reports hold "
                         f"{', '.join(names)})")
    return {"species": readings, "panel": read_panel(reports, readings, predictions, predict_on)}


def _n(value, spec=".4f"):
    return "-" if value is None else format(value, spec)


def _count(c):
    return f"{c['rejections']} of {c['of']} (tail {_n(c['tail'])}, flagged {c['flagged']})"


def render(result):
    """The reading as text, numbers rounded; the JSON keeps them whole."""
    lines = []
    for r in result["species"]:
        lines += [f"{r['species']}: {r['datasets']} datasets, {len(r['unanswered'])} the test can't answer; "
                  f"{r['made']} of {r['planned']} family C contrasts made",
                  f"{'arm':22}{'plan':>5}{'made':>5}{'impr':>5}{'worse':>6}{'inside':>7}{'beyond':>7}{'zero':>5}"
                  f"{'q1':>9}{'median':>9}{'q3':>9}{'auroc':>9} | {'datasets':>8}{'better':>7}{'worse':>6}"
                  f"{'sign p':>8}"]
        for arm, a in r["arms"].items():
            s = r["by_dataset"][arm]
            lines.append(f"{arm:22}{a['planned']:5}{a['made']:5}{a['improves']:5}{a['worse']:6}{a['inside']:7}"
                         f"{a['beyond']:7}{a['zero']:5}{_n(a['q1']):>9}{_n(a['median']):>9}{_n(a['q3']):>9}"
                         f"{_n(a['auroc_median']):>9} | {s['datasets']:8}{s['better']:7}{s['worse']:6}"
                         f"{_n(s['p']):>8}")
        for learner, c in r["calibration"].items():
            lines.append(f"calibration {learner:6} {', '.join(c['twins']) or 'no twin'}: {_count(c)}; "
                         f"{c['left_out']} left out, {c['no_test']} counted with no test; rule 10 marks "
                         f"{c['rule10_contrasts']} contrasts")
        lines.append("")
    p = result["panel"]
    lines += [f"panel: one BH family of {p['family']['size']} ({p['family']['tested']} tested, "
              f"{p['family']['entered_at_p_1']} at p = 1)",
              f"firm decisions: {p['firm']['improves']} improves, {p['firm']['worse']} worse; "
              f"{p['one_family_only']} decided under one family only; {p['flagged_or_unjudged']} on a flagged or "
              "unjudged learner"]
    for d in p["decisions"]:
        lines.append(f"  {'firm' if d['firm'] else '    '} species {d['in_species']:12} panel {d['in_panel']:12} "
                     f"{d['id']} estimate {_n(d['estimate'])} p_adj {_n(d['p_adj_species'], '.4g')} / "
                     f"{_n(d['p_adj_panel'], '.4g')} auroc {_n(d['auroc_delta'])} learner flagged "
                     f"{d['learner_flagged']}")
    for learner, c in p["calibration"].items():
        ci = c["ci95"]
        lines.append(f"calibration {learner:6} over the species: {_count(c)}"
                     + (f", 95% interval {ci[0]:.4f} to {ci[1]:.4f}" if ci else ""))
    for arm, s in p["by_dataset"].items():
        lines.append(f"by dataset {arm:22} {s['better']:3} better {s['worse']:3} worse of {s['datasets']:3}, "
                     f"median {_n(s['median'])}, sign p {_n(s['p'])}")
    if p["predictions"]:
        lines.append(f"predictions on {', '.join(p['predictions']['species'])}: pass {p['predictions']['pass']}")
        for t in p["predictions"]["tests"]:
            lines.append(f"  {t['arm']} predicted {t['predicted']}: {t['better']} better, {t['worse']} worse, "
                         f"{t['unchanged']} unchanged, one-sided p {_n(t['p_one_sided'])}")
    return "\n".join(lines)


def prediction(text):
    arm, _, side = text.rpartition("=")
    if not arm or side not in SIDES:
        raise argparse.ArgumentTypeError(f"{text!r} is not ARM=better or ARM=worse")
    return arm, side


def main(argv):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("reports", nargs="+", metavar="REPORT")
    p.add_argument("--json")
    p.add_argument("--predict", action="append", default=[], type=prediction, metavar="ARM=better|worse")
    p.add_argument("--predict-on", action="append", default=[], metavar="SPECIES")
    a = p.parse_args(argv)
    if a.predict_on and not a.predict:
        p.error("--predict-on needs a --predict")
    result = read([k_report.load_report(path) for path in a.reports], a.predict, a.predict_on)
    if a.json:
        Path(a.json).write_text(json.dumps(result, indent=1, sort_keys=True))
    print(render(result))


if __name__ == "__main__":
    main(sys.argv[1:])
