#!/usr/bin/env python3
"""A third reading of the input comparison: per learner, the run's calls beside its log loss.

Exploratory. It was designed on 11 Nov 2026, after H-1's results had been read, and it decides nothing: log loss stays
the registered endpoint, K's reports stand as written, and second_reading.py is unchanged
(docs/input-comparison-runner.md, "A per-learner reading with calls").

COMBINE is a run's or a combine's WORK folder, one species each: k/manifest.json, j/contrasts.json, j/planned.json,
j/missing.json, and arms/<arm>/<dataset>/{runs.json,predictions.jsonl,arm.json,plans/}. A species given twice is
refused, and so are folders that differ in source commit, grm binary or the synthetic flag.

What it takes from where:
  log loss   each contrast's estimate as J recorded it in j/contrasts.json (above zero: the arm's log loss is lower).
  pairing    the record's own pairing rule through stats.pair_units, the function J paired with: per repeat, the
             genomes the arm scored, the default restricted to them. A contrast J didn't make (j/missing.json) is
             left out. A record whose rows no longer pair, or whose estimate or unit count stats.compare_units doesn't
             give back from the rows under the learner's registered test, refuses the folder: the predictions are
             not the ones J read.
  calls      the row's `called`, the run's own call: `grm arm` writes called = score >= 0.5, and on a no-call row
             (no model, NaN score) called = false with the training prevalence as score (src/ablate/arm.rs at
             31b7edb, where it writes each PredictionRow). A no-call is then a susceptible call, as in J's thresholded
             endpoints. Each (repeat, genome) pair is one call, so three repeats count a genome three times.

Per split (lineage, which decides, and random, kept apart), arm and learner, over the folders' datasets:
  log loss   the datasets with an answer (a made contrast with an estimate); better, worse and zero by the sign of
             the estimate, with the two-sided exact sign test; median and mean.
  calls      TP, FN, FP, TN of the default and of the arm, summed over datasets and repeats; recall, specificity,
             balanced accuracy, precision and F1 from those sums, and the arm minus the default. A rate whose
             denominator is zero is null.
  disagree   the datasets where log loss improves (estimate > 0) while recall falls, and those where log loss is
             worse (< 0) while recall rises; recall per dataset from its own sums, ties and undefined recall left
             out. The twin arms change no input, so their counts are the reference for no change.
  share      each learner's part of the five learners' pooled change for that arm: its sum over datasets divided by
             the sum over every learner and dataset. Two weightings, each labelled: "units", the run's estimate,
             where every (repeat, outer fold) unit weighs the same; "genomes", each repeat's mean over its paired
             genomes, then the mean over repeats.
  cohort     for an arm under a subset pairing (another cohort: cap or selection), per dataset the default cohort's
             class counts from its plan files' labels and both arms' genome counts (arm.json). The arm's own labels
             are not in the folder, so its class mix prints as "not recorded".

--sensitivity adds a sensitivity reading. It loads each folder's report/ablation-report.json through K's validation,
recomputes every family C record from the rows with one test for every learner (the paired t, then the corrected t),
finalizes each species as J does, and counts second_reading.py's firm decisions on the result. It first recomputes
the records under the registered tests and refuses unless they are J's. A different formula doesn't validate a test,
and the datasets share genomes and species, so the count decides nothing.

usage: python3 -B learner_reading.py COMBINE [COMBINE ...] [--json OUT] [--sensitivity]
"""
import argparse
import importlib.util
import json
import math
import statistics
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "ablation"))
import stats  # noqa: E402  J, for its pairing, its unit test and its finalize
sys.path.insert(0, str(ROOT))
from scripts.ablation import report as k_report  # noqa: E402  K

_spec = importlib.util.spec_from_file_location("input_comparison_second_reading", HERE / "second_reading.py")
second_reading = sys.modules.get(_spec.name)
if second_reading is None:
    second_reading = importlib.util.module_from_spec(_spec)
    sys.modules[_spec.name] = second_reading
    _spec.loader.exec_module(second_reading)

DEFAULT = "default"
LEARNERS = ("scm", "cart", "rf", "gbt", "lasso")
SPLITS = ("lineage", "random")
COUNTS = ("tp", "fn", "fp", "tn")
RATES = ("recall", "specificity", "balanced_accuracy", "precision", "f1")
WEIGHTINGS = {"units": "the run's estimate: every (repeat, outer fold) unit weighs the same",
              "genomes": "each repeat's mean over its paired genomes, then the mean over repeats"}
WEIGHTED = {"units": "estimate", "genomes": "genome_weighted"}  # each weighting's field in a row
UNIFORM = {"paired_t": stats.UNIT_T, "corrected_t": stats.UNIT_CORRECTED_T}
CAVEAT = ("sensitivity reading only: a different formula doesn't validate a test, and the datasets share genomes and "
          "species, so this count decides nothing")
# J's recorded estimate against the one compare_units gives back from the rows: the same code on the same floats.
TOLERANCE = 1e-12


def refuse(message):
    raise SystemExit(f"learner_reading: {message}")


def _json(path):
    try:
        return json.loads(Path(path).read_text())
    except FileNotFoundError:
        refuse(f"{path} is missing")
    except json.JSONDecodeError as e:
        refuse(f"{path} is not JSON: {e}")


def read_rows(path):
    """{(learner, split): [row, ...]} of one predictions.jsonl."""
    out = defaultdict(list)
    try:
        with open(path, encoding="utf-8") as f:
            for line in f:
                if line.strip():
                    row = json.loads(line)
                    out[(row["learner"], row["split"])].append(row)
    except FileNotFoundError:
        refuse(f"{path} is missing")
    return out


def rates(c):
    """The call rates of one confusion count, None where a denominator is zero."""
    def ratio(a, b):
        return a / b if b else None
    recall, specificity = ratio(c["tp"], c["tp"] + c["fn"]), ratio(c["tn"], c["tn"] + c["fp"])
    return {"recall": recall, "specificity": specificity,
            "balanced_accuracy": None if recall is None or specificity is None else (recall + specificity) / 2,
            "precision": ratio(c["tp"], c["tp"] + c["fp"]),
            "f1": ratio(2 * c["tp"], 2 * c["tp"] + c["fp"] + c["fn"])}


def calls(rows, paired):
    """Confusion counts of both arms and the genome-weighted log-loss change on the genomes `paired` holds.

    The rows are rebuilt in pair_units' own order (repeats in order, the arm's genomes of each repeat sorted), and
    their labels and scores must equal its arrays, so the counts stand on exactly the genomes J paired."""
    tables = {DEFAULT: {}, paired.treatment: {}}
    for row in rows:
        if row["arm"] in tables:
            tables[row["arm"]][(row["repeat"], row["sample_id"])] = row
    arm_ids = defaultdict(set)
    for repeat, sample in tables[paired.treatment]:
        arm_ids[repeat].add(sample)
    repeats = sorted({repeat for repeat, _ in tables[DEFAULT]}, key=str)
    base_count, arm_count = dict.fromkeys(COUNTS, 0), dict.fromkeys(COUNTS, 0)
    y, base, treat, per_repeat = [], [], [], []
    for repeat in repeats:
        start = len(y)
        for sample in sorted(arm_ids[repeat], key=str):
            b, t = tables[DEFAULT][(repeat, sample)], tables[paired.treatment][(repeat, sample)]
            for count, row in ((base_count, b), (arm_count, t)):
                count[("tp" if row["called"] else "fn") if b["y"] else ("fp" if row["called"] else "tn")] += 1
            y.append(b["y"])
            base.append(float(b["score"]))
            treat.append(float(t["score"]))
        per_repeat.append((start, len(y)))
    y, base, treat = np.array(y, dtype=np.int64), np.array(base), np.array(treat)
    if not (np.array_equal(y, paired.y) and np.array_equal(base, paired.base) and np.array_equal(treat, paired.treat)):
        refuse(f"{paired.dataset} {paired.learner} {paired.split} {paired.treatment}: the call rows don't line up "
               "with the genomes pair_units paired")
    delta = stats.sample_log_loss(y, base) - stats.sample_log_loss(y, treat)
    genomes = statistics.fmean(float(delta[a:b].mean()) for a, b in per_repeat if b > a)
    return {"calls_default": base_count, "calls_arm": arm_count, "repeats": sum(b > a for a, b in per_repeat),
            "paired_calls": len(y), "genome_weighted": genomes}


def _same(a, b):
    if a is None or b is None:
        return a is None and b is None
    return math.isclose(a, b, rel_tol=TOLERANCE, abs_tol=TOLERANCE)


def recompute(paired, record, test):
    """compare_units' pending record for `paired` under `test`, in the record's family and margin."""
    (pending,) = stats.compare_units(paired, family=record["family"], margins={"log_loss": record["margin"]},
                                     development=[], test=test)
    return pending


def read_folder(folder, sensitivity=False):
    """One species' WORK folder: its per-contrast rows, cohort lines and, with `sensitivity`, J's pending records
    under each uniform test."""
    folder = Path(folder)
    manifest = _json(folder / "k" / "manifest.json")
    names = sorted({d["species"] for d in manifest["datasets"]})
    if len(names) != 1:
        refuse(f"{folder}: a folder must hold one species, and this one holds {', '.join(names)}")
    records = _json(folder / "j" / "contrasts.json")
    planned = _json(folder / "j" / "planned.json")
    missing = _json(folder / "j" / "missing.json")
    by_id = {x["id"]: x for x in records}
    if len(by_id) != len(records):
        refuse(f"{folder}: j/contrasts.json holds an id twice")
    unplanned = sorted(set(by_id) - set(planned))
    unaccounted = sorted(set(planned) - set(by_id) - set(missing))
    if unplanned or unaccounted or set(by_id) & set(missing):
        refuse(f"{folder}: j/ doesn't account for every planned contrast once (unplanned {unplanned[:3]}, neither "
               f"made nor missing {unaccounted[:3]}, both {sorted(set(by_id) & set(missing))[:3]})")
    sentinels = set(manifest.get("sentinel_arms", []))
    rows_out, cohort, pending = [], [], {name: [] for name in ("registered", *UNIFORM)}
    for d in manifest["datasets"]:
        ds = d["id"]
        mine = [x for x in records if x["dataset"] == ds]
        if not mine:
            continue
        default_rows = read_rows(folder / "arms" / DEFAULT / ds / "predictions.jsonl")
        for arm in sorted({x["treatment"] for x in mine}):
            arm_rows = read_rows(folder / "arms" / arm / ds / "predictions.jsonl")
            for x in (x for x in mine if x["treatment"] == arm):
                rows = default_rows[(x["learner"], x["split"])] + arm_rows[(x["learner"], x["split"])]
                try:
                    paired = stats.pair_units(rows, dataset=ds, learner=x["learner"], split=x["split"],
                                              baseline=DEFAULT, treatment=arm, pairing=x["pairing"])
                except stats.PairingError as e:
                    refuse(f"{folder}: {x['id']} is a record in j/contrasts.json, and its rows don't pair: {e}")
                registered = recompute(paired, x, stats.LEARNER_TESTS[x["learner"]])
                if not (_same(registered["estimate"], x["estimate"]) and registered["n_units"] == x["n_units"]):
                    refuse(f"{folder}: {x['id']}: the rows give estimate {registered['estimate']} over "
                           f"{registered['n_units']} units, and J recorded {x['estimate']} over {x['n_units']}")
                if sensitivity:
                    pending["registered"].append(registered)
                    for name, test in UNIFORM.items():
                        pending[name].append(registered if test == stats.LEARNER_TESTS[x["learner"]]
                                             else recompute(paired, x, test))
                if x["estimate"] is None:
                    continue
                rows_out.append({"species": names[0], "dataset": ds, "arm": arm, "learner": x["learner"],
                                 "split": x["split"], "twin": arm in sentinels, "pairing": x["pairing"],
                                 "estimate": x["estimate"], "decision": x["decision"], **calls(rows, paired)})
            if any(x["pairing"] in stats.SUBSET_PAIRINGS for x in mine if x["treatment"] == arm):
                cohort.append(cohort_line(folder, names[0], ds, arm))
    answered = sorted({x["dataset"] for x in records
                       if x["family"] == "C" and x["p"] is not None and x["treatment"] not in sentinels})
    return {"folder": str(folder), "species": names[0], "manifest": manifest, "datasets": len(manifest["datasets"]),
            "lineage_answer": len(answered), "planned": planned, "missing": missing, "rows": rows_out, "cohort": cohort, "pending": pending if sensitivity else None}


def cohort_line(folder, species, ds, arm):
    """The default cohort's class counts from its plan files' labels, and both arms' genome counts. The arm's
    labels are not in the folder."""
    plans = folder / "arms" / DEFAULT / ds / "plans"
    labels = {json.dumps([p["sample_ids"], p["y"]]) for p in
              (_json(f) for f in sorted(plans.glob("*.json"))
               if not f.name.startswith("refused-") and not f.name.endswith(".cohort.json"))}
    if len(labels) > 1:
        refuse(f"{plans}: the default's plan files disagree on its genomes or labels")
    y = json.loads(labels.pop())[1] if labels else None
    return {"species": species, "dataset": ds, "arm": arm,
            "default_genomes": _json(folder / "arms" / DEFAULT / ds / "arm.json")["samples"],
            "default_resistant": None if y is None else sum(y),
            "default_susceptible": None if y is None else len(y) - sum(y),
            "arm_genomes": _json(folder / "arms" / arm / ds / "arm.json")["samples"],
            "arm_class_mix": "not recorded"}


def summarize(rows):
    """The per-(split, arm, learner) table over the rows of every folder (the module docstring's block), and per
    (split, arm) the contrasts and the pooled sums each share divides by."""
    groups = defaultdict(list)
    for r in rows:
        groups[(r["split"], r["arm"], r["learner"])].append(r)
    pooled = {}
    for r in rows:
        total = pooled.setdefault(r["split"], {}).setdefault(r["arm"],
                                                             {"contrasts": 0, **dict.fromkeys(WEIGHTINGS, 0.0)})
        total["contrasts"] += 1
        for w, field in WEIGHTED.items():
            total[w] += r[field]
    out = {}
    for (split, arm, learner), rs in sorted(groups.items(), key=lambda kv: (SPLITS.index(kv[0][0]), kv[0][1],
                                                                           LEARNERS.index(kv[0][2]))):
        estimates = [r["estimate"] for r in rs]
        sign = second_reading.sign_test(estimates)
        sums = {side: {k: sum(r[f"calls_{side}"][k] for r in rs) for k in COUNTS} for side in (DEFAULT, "arm")}
        got = {side: rates(sums[side]) for side in sums}
        disagree = {"improves_recall_falls": 0, "worse_recall_rises": 0, "recall_defined": 0}
        for r in rs:
            before, after = rates(r["calls_default"])["recall"], rates(r["calls_arm"])["recall"]
            if before is None or after is None:
                continue
            disagree["recall_defined"] += 1
            disagree["improves_recall_falls"] += r["estimate"] > 0 and after < before
            disagree["worse_recall_rises"] += r["estimate"] < 0 and after > before
        total = pooled[split][arm]
        share = {w: sum(r[field] for r in rs) / total[w] if total[w] else None for w, field in WEIGHTED.items()}
        out.setdefault(split, {}).setdefault(arm, {})[learner] = {
            "twin": rs[0]["twin"], **sign, "mean": statistics.fmean(estimates), "default": sums[DEFAULT],
            "arm": sums["arm"], "rates_default": got[DEFAULT], "rates_arm": got["arm"],
            "rates_change": {k: None if got["arm"][k] is None or got[DEFAULT][k] is None
                             else got["arm"][k] - got[DEFAULT][k] for k in RATES},
            **disagree, "share": share}
    return out, pooled


def sensitivity(folders, reports):
    """Firm decisions (second_reading.read) on each species' report with every family C record recomputed under
    one test for every learner. The registered tests must give J's records back first."""
    def finalized(f, name):
        return stats.finalize(f["pending"][name], f["planned"], f["missing"],
                              k_report.unadjusted_ids(f["manifest"], f["planned"]))

    for f, report in zip(folders, reports):
        recorded = {x["id"]: x for x in report["contrasts"]}
        final = finalized(f, "registered")
        if len(final) != len(recorded) or any(
                x["id"] not in recorded or x["decision"] != recorded[x["id"]]["decision"]
                or not all(_same(x[k], recorded[x["id"]][k]) for k in ("estimate", "p", "p_adj")) for x in final):
            refuse(f"{f['folder']}: its records recomputed under the registered tests aren't its report's")
    readings = {"registered": second_reading.read(reports)}
    for name in UNIFORM:
        readings[name] = second_reading.read([{**report, "contrasts": finalized(f, name)}
                                              for f, report in zip(folders, reports)])
    out = {"caveat": CAVEAT, "tests": {"registered": dict(stats.LEARNER_TESTS),
                                       **{name: dict.fromkeys(LEARNERS, test) for name, test in UNIFORM.items()}}}
    for name, reading in readings.items():
        firm = [d for d in reading["panel"]["decisions"] if d["firm"]]
        by_learner = {learner: sum(d["id"].split("/")[1].split(":")[1] == learner for d in firm)
                      for learner in LEARNERS}
        out[name] = {"firm": reading["panel"]["firm"], "firm_by_learner": by_learner,
                     "flagged": {learner: c["flagged"] for learner, c in reading["panel"]["calibration"].items()},
                     "decisions": reading["panel"]["decisions"]}
    return out


def read(folders, reports=None):
    """The whole reading of read_folder() outputs; with `reports` (each folder's loaded K report, in order) the
    sensitivity reading too."""
    names = [f["species"] for f in folders]
    twice = sorted({name for name in names if names.count(name) > 1})
    if twice:
        refuse(f"a species is read from one folder, and these have more: {', '.join(twice)}")
    for field in second_reading.SHARED:
        values = sorted({str(f["manifest"].get(field)) for f in folders})
        if len(values) > 1:
            refuse(f"the folders differ in {field} ({', '.join(values)}), so they aren't one comparison")
    rows = [r for f in folders for r in f["rows"]]
    tables, pooled = summarize(rows)
    return {"exploratory": "designed on 11 Nov 2026, after the results were read; it decides nothing",
            "call_rule": "the row's called: score >= 0.5, a no-call row false (src/ablate/arm.rs at 31b7edb)",
            "weightings": WEIGHTINGS,
            "species": [{k: f[k] for k in ("species", "folder", "datasets", "lineage_answer")} for f in folders],
            "tables": tables, "pooled": pooled, "rows": rows, "cohort": [c for f in folders for c in f["cohort"]],
            "sensitivity": sensitivity(folders, reports) if reports is not None else None}


def _n(value, spec=".4f"):
    return "-" if value is None else format(value, spec)


def _move(rates_default, rates_arm, key):
    return f"{_n(rates_default[key], '.3f')}>{_n(rates_arm[key], '.3f')}"


def render(result):
    """The reading as text, numbers rounded; the JSON keeps them whole."""
    lines = [f"Per-learner reading with calls: exploratory, {result['exploratory']}.",
             f"Call: {result['call_rule']}; each (repeat, genome) is one call.", ""]
    for s in result["species"]:
        lines.append(f"{s['species']}: {s['datasets']} datasets, {s['lineage_answer']} with a lineage answer")
    for split, arms in result["tables"].items():
        lines += ["", f"== {split} split{' (decides)' if split == 'lineage' else ' (kept apart; decides nothing)'} ==",
                  "log loss: datasets, better/worse/zero, two-sided sign p, median, mean; share of the five learners' "
                  "pooled change, units = " + WEIGHTINGS["units"] + ", genomes = " + WEIGHTINGS["genomes"],
                  f"{'arm':22}{'learner':8}{'n':>4}{'better':>7}{'worse':>6}{'zero':>5}{'sign p':>8}{'median':>9}"
                  f"{'mean':>9}{'share units':>12}{'share genomes':>14}"]
        for arm, learners in arms.items():
            for learner, t in learners.items():
                lines.append(f"{arm:22}{learner:8}{t['datasets']:4}{t['better']:7}{t['worse']:6}{t['unchanged']:5}"
                             f"{_n(t['p']):>8}{_n(t['median']):>9}{_n(t['mean']):>9}"
                             f"{_n(t['share']['units'], '.3f'):>12}{_n(t['share']['genomes'], '.3f'):>14}")
            total = result["pooled"][split][arm]
            lines.append(f"{arm:22}{'pooled':8}{total['contrasts']:4} contrasts, sum units {total['units']:+.4f}, "
                         f"sum genomes {total['genomes']:+.4f}")
        lines += ["", "calls, summed over datasets and repeats: default TP FN FP TN | arm TP FN FP TN | default > arm "
                      "recall, specificity, balanced accuracy, precision, F1 | datasets where log loss improves and "
                      "recall falls, worse and recall rises, of those with recall on both"]
        for arm, learners in arms.items():
            for learner, t in learners.items():
                d, a = t["default"], t["arm"]
                lines.append(f"{arm:22}{learner:6} " + " ".join(f"{d[k]:6}" for k in COUNTS) + " |"
                             + " ".join(f"{a[k]:6}" for k in COUNTS) + " | "
                             + " ".join(_move(t["rates_default"], t["rates_arm"], k) for k in RATES)
                             + f" | {t['improves_recall_falls']} {t['worse_recall_rises']} of {t['recall_defined']}")
    if result["cohort"]:
        lines += ["", "cohort arms: the default cohort's class counts (its plan files' labels) and genomes; the arm's "
                      "genomes (arm.json); the arm's class mix"]
        for c in result["cohort"]:
            lines.append(f"{c['arm']:10}{c['dataset']:52} default {c['default_genomes']:5} genomes, "
                         f"{_n(c['default_resistant'], 'd')} R {_n(c['default_susceptible'], 'd')} S | arm "
                         f"{c['arm_genomes']:5} genomes, class mix {c['arm_class_mix']}")
    s = result["sensitivity"]
    if s:
        lines += ["", f"SENSITIVITY READING ({s['caveat']})"]
        for name in ("registered", *UNIFORM):
            tests = ", ".join(f"{learner} {test}" for learner, test in s["tests"][name].items())
            f = s[name]
            flagged = ", ".join(k for k, v in f["flagged"].items() if v) or "none"
            lines.append(f"  {name:12} firm {f['firm']['improves']} improves, {f['firm']['worse']} worse; by learner "
                         + ", ".join(f"{k} {v}" for k, v in f["firm_by_learner"].items())
                         + f"; calibration over the species flags {flagged} ({tests})")
    return "\n".join(lines)


def main(argv):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("folders", nargs="+", metavar="COMBINE")
    p.add_argument("--json")
    p.add_argument("--sensitivity", action="store_true")
    a = p.parse_args(argv)
    folders = [read_folder(path, a.sensitivity) for path in a.folders]
    reports = ([k_report.load_report(Path(path) / "report" / "ablation-report.json") for path in a.folders]
               if a.sensitivity else None)
    result = read(folders, reports)
    if a.json:
        Path(a.json).write_text(json.dumps(result, indent=1, sort_keys=True))
    print(render(result))


if __name__ == "__main__":
    main(sys.argv[1:])
