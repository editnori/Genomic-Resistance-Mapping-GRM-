#!/usr/bin/env python3
"""Reads the staggered split's result files against its design, and writes its tables, rule and predictions.

  python3 -I stagger_analyze.py <design dir> <result dir> <out dir> [--listed <list dir>]

<design dir> holds the design.json that `stagger_design.py folds` wrote; <result dir> the files `grm learn
--fold-file` wrote, <dataset>.<scheme>.json. It stops, writing nothing, unless each result names its fold file's
sha256, both fits of a pair scored the same genomes with the same labels, the training counts equal the design's
and recall worked out again equals the evaluator's: the checks of 14 Nov, unchanged.

Per clock, contrast, horizon, learner and measure it gives the median over datasets of the paired difference,
older fit minus the fresh fit of the window's own origin, and its 95% interval from 10,000 draws of species with
replacement, seed 20261125, the 251st and the 9,750th of the sorted draws, read as worse, better or not shown by
the rule of 14 Nov. A dataset with several pairs at one horizon gives their mean. Then the plan's rule for a
learner and for a contrast, descriptive only under 12 datasets or 3 species.

With --listed (one file a dataset of M, <dataset>.txt, one genome id a line, as m_list.py writes them) each
contrast is read twice, by R33 points 8 to 10. The second reading takes the listed genomes out of a pair's test
genomes for both fits, with no refit. A pair stays only when its test set lists no genome or at least 5 of each
class, and each class keeps the test minimum once they are out. lacking.tsv names each dataset the second reading
lacks at a horizon, with its reason. A pair left out for a short list gives no figure anywhere, and no count over
datasets takes it in. A verdict stands only when both readings give it.

What it writes holds rows per dataset, species, clock, contrast, horizon, arm and learner, and no genome id, study
or row per genome: it stops before writing if a genome id of the results or the lists stands in any of its text.

Its draws: one generator a clock and reading, seeded 20261125, drawn learner by learner, then by contrast, horizon
and measure, as the 14 Nov script drew learner, contrast and measure.

The predictions, each met or not met by a rule fixed here before any result exists (PREDICTIONS). A prediction is
read in each reading made and is met only when every reading made meets it:

- P1, public-day clock, stale same size, horizon 1: met when the contrast there is not descriptive only and recall
  reads not shown for at least 4 of 5 learners.
- P2, public-day clock, stale all, horizon 1: met when the contrast there is not descriptive only and at least 3 of
  5 learners show a cost with log loss reading worse.
- P3, collection-year clock: met when both contrasts at horizon 1 read not shown; descriptive only is not not
  shown.
"""

import importlib.util
import json
import math
import random
import re
import statistics
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("stagger_design_for_analyze", HERE / "stagger_design.py")
design_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(design_mod)

Stop, need, guard = design_mod.Stop, design_mod.need, design_mod.guard
CLOCKS, CONTRASTS, ARM_OF = design_mod.CLOCKS, design_mod.CONTRASTS, design_mod.ARM_OF
LEARNERS = ("scm", "cart", "rf", "gbt", "lasso")
MEASURES = ("recall", "log_loss", "specificity", "balanced_accuracy")
RULE_MEASURES = ("recall", "log_loss")
WORSE = {"recall": -1, "log_loss": 1, "specificity": -1, "balanced_accuracy": -1}
CUTOFF, CLIP = 0.5, 1e-4
RESAMPLES, SEED = 10000, 20261125
LOW_DRAW, HIGH_DRAW = 251, 9750
MIN_DATASETS, MIN_SPECIES, MIN_LEARNERS_WITH_COST = 12, 3, 3
LIST_MINIMUM = 5
REVISION = 14
RESULT_FORMAT = "grm-fold-file-result-v1"
SETTINGS = {"seed": 42, "screened_features": 4000, "inner_folds": 3}
READINGS = ("first", "second")
SHORT = "its test set lists genomes, under 5 of a class"
UNDER = "a class keeps under the test minimum once the listed genomes are out"
SHAPES = (r"\b\d{3,7}\.\d{3,6}\b", r"\bPRJ[A-Z]{1,2}\d+", r"\bSAM[NED][A-Z]?\d+", r"\bGC[AF]_\d+",
          r"\b[SE]R[RSX]\d{5,}", r"[ACGT]{20,}")

PREDICTIONS = (
    ("P1", "public-day", "stale-same-size", 1,
     "stale same size, horizon 1, recall: not shown for at least 4 of 5 learners"),
    ("P2", "public-day", "stale-all", 1, "stale all, horizon 1: a cost through log loss for at least 3 of 5"),
    ("P3", "collection-year", None, 1, "collection-year clock: not shown on every contrast, read at horizon 1"),
)

ARM_COLUMNS = ["clock", "reading", "dataset", "species", "contrast", "horizon", "arm", "learner", "pairs", "test",
               "test_resistant"] + list(MEASURES)
CONTRAST_COLUMNS = ["clock", "reading", "dataset", "species", "contrast", "horizon", "learner", "pairs", "test",
                    "test_resistant"] + ["d_" + m for m in MEASURES]
READING_COLUMNS = ["clock", "reading", "contrast", "horizon", "learner", "measure", "datasets", "species", "median",
                   "low", "high", "read_as", "datasets_worse", "datasets_better", "datasets_equal"]
RULE_COLUMNS = ["clock", "reading", "contrast", "horizon", "learner", "recall", "log_loss", "rule"]
VERDICT_COLUMNS = ["clock", "reading", "contrast", "horizon", "datasets", "species", "learners_with_cost",
                   "learners_better", "verdict"]
STAND_COLUMNS = ["clock", "contrast", "horizon", "first", "second", "stands"]
LACK_COLUMNS = ["clock", "contrast", "horizon", "dataset", "species", "reason"]
PREDICTION_COLUMNS = ["prediction", "rule", "first", "second", "outcome"]


# ----------------------------------------------------------------- measures

def measure(samples, scores, keep=None):
    """Recall, specificity, balanced accuracy and mean clipped log loss at the 0.5 cutoff, in sample order, over
    the samples whose index is in keep (all when keep is None)."""
    tp = fn = tn = fp = 0
    loss = 0.0
    for i, ((_, t), s) in enumerate(zip(samples, scores)):
        if keep is not None and i not in keep:
            continue
        called = s >= CUTOFF
        tp += t and called
        fn += t and not called
        tn += (not t) and not called
        fp += (not t) and called
        p = min(max(s, CLIP), 1 - CLIP)
        loss -= math.log(p if t else 1 - p)
    recall, specificity = tp / (tp + fn), tn / (tn + fp)
    return {"tp": tp, "fn": fn, "tn": tn, "fp": fp, "recall": recall, "specificity": specificity,
            "balanced_accuracy": (recall + specificity) / 2, "log_loss": loss / (tp + fn + tn + fp)}


def interval(by_species, rng):
    """The median over datasets, and the 251st and 9,750th of 10,000 medians over species drawn with replacement."""
    species = sorted(by_species)
    point = statistics.median([v for s in species for v in by_species[s]])
    draws = []
    for _ in range(RESAMPLES):
        values = [v for s in rng.choices(species, k=len(species)) for v in by_species[s]]
        draws.append(statistics.median(values))
    draws.sort()
    return point, draws[LOW_DRAW - 1], draws[HIGH_DRAW - 1]


def read(lo, hi, m):
    """Worse when the interval lies wholly on the worse side of zero, better on the other, else not shown."""
    if lo > 0 or hi < 0:
        return "worse" if (lo > 0) == (WORSE[m] > 0) else "better"
    return "not shown"


def learner_rule(recall, log_loss):
    """The plan's rule for one learner: a cost when recall or log loss reads worse and neither reads better."""
    got = (recall, log_loss)
    if "better" in got:
        return "mixed" if "worse" in got else "benefit"
    return "cost" if "worse" in got else "not shown"


def verdict(datasets, species, rules, reads, descriptive_clock):
    """The plan's rule for a contrast at a horizon."""
    if descriptive_clock or datasets < MIN_DATASETS or species < MIN_SPECIES:
        return "descriptive only"
    if any(r == "better" for l in rules for r in reads[l]):
        return "mixed"
    if sum(rules[l] == "cost" for l in rules) >= MIN_LEARNERS_WITH_COST:
        return "cost shown"
    return "not shown"


def second_case(samples, listed, test_minimum):
    """R33 point 9 for one pair: (kept indices or None, reason). None as kept means the pair leaves the reading."""
    on_list = [i for i, (g, _) in enumerate(samples) if g in listed]
    if not on_list:
        return set(range(len(samples))), ""
    lr = sum(1 for i in on_list if samples[i][1])
    ls = len(on_list) - lr
    if min(lr, ls) < LIST_MINIMUM:
        return None, SHORT
    r = sum(1 for _, t in samples if t) - lr
    s = len(samples) - sum(1 for _, t in samples if t) - ls
    if min(r, s) < test_minimum:
        return None, UNDER
    out = set(on_list)
    return {i for i in range(len(samples)) if i not in out}, ""


def prediction_met(pid, verdicts, rules, reads):
    """The rule of one prediction on one reading's verdicts, learner rules and readings."""
    _, clock, contrast, h, _ = next(p for p in PREDICTIONS if p[0] == pid)
    if pid == "P3":
        return all(verdicts.get((clock, c, h)) == "not shown" for c in CONTRASTS)
    v = verdicts.get((clock, contrast, h))
    if v is None or v == "descriptive only":
        return False
    if pid == "P1":
        return sum(reads[(clock, contrast, h, l)][0] == "not shown" for l in LEARNERS) >= 4
    return sum(rules[(clock, contrast, h, l)] == "cost" and reads[(clock, contrast, h, l)][1] == "worse"
               for l in LEARNERS) >= 3


# ----------------------------------------------------------------- results

def load_result(results, d, fit, window):
    """One result file, checked against the design. (samples, {learner: scores})."""
    name = f"{d['id']}.{fit['scheme']}.json"
    path = guard(results / name, "result file")
    need(path.exists(), f"{name}: no result file")
    try:
        r = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        raise Stop(f"{name}: not JSON")
    need(r.get("format") == RESULT_FORMAT, f"{name}: not {RESULT_FORMAT}")
    need(r.get("evaluation_revision") == REVISION, f"{name}: not revision {REVISION}")
    need(r.get("id") == d["id"] and r.get("scheme") == fit["scheme"], f"{name}: another dataset or scheme")
    need(r.get("fold_file_sha256") == fit["sha256"], f"{name}: does not name its fold file's sha256")
    try:
        return checked_result(name, r, d, fit, window)
    except (KeyError, TypeError, IndexError):
        raise Stop(f"{name}: a field the evaluator writes is missing")


def checked_result(name, r, d, fit, window):
    samples = []
    for s in r["samples"]:
        need(isinstance(s.get("genome_id"), str) and isinstance(s.get("resistant"), bool),
             f"{name}: a sample without an id or a label")
        samples.append((s["genome_id"], s["resistant"]))
    truth = [t for _, t in samples]
    need(len(truth) == window["test"] and sum(truth) == window["test_resistant"], f"{name}: test counts")
    need(len(r["folds"]) == 1, f"{name}: not one fold")
    fold = r["folds"][0]
    need(fold["test"] == len(truth) and fold["test_resistant"] == sum(truth), f"{name}: the fold's test counts")
    c = r["eval_config"]
    need({k: c.get(k) for k in SETTINGS} == SETTINGS, f"{name}: not the default settings")
    learners = r["split"]["learners"]
    need([x["learner"] for x in learners] == list(LEARNERS), f"{name}: not the five learners in order")
    scores = {}
    for x in learners:
        sc = x["scores"]
        need(len(sc) == len(samples), f"{name} {x['learner']}: {len(sc)} scores for {len(samples)} genomes")
        need(all(isinstance(v, (int, float)) and not isinstance(v, bool) and 0 <= v <= 1 for v in sc),
             f"{name} {x['learner']}: a score outside 0 to 1")
        m, e = measure(samples, sc), x["metrics"]
        need((m["tp"], m["fn"], m["tn"], m["fp"]) == (e["tp"], e["fn"], e["tn"], e["fp"]),
             f"{name} {x['learner']}: confusion counts other than the evaluator's")
        need(abs(m["recall"] - e["sensitivity"]) < 1e-12 and abs(m["specificity"] - e["specificity"]) < 1e-12,
             f"{name} {x['learner']}: recall worked out again is not the evaluator's")
        scores[x["learner"]] = sc
    return samples, scores, fold


def read_lists(folder):
    """{dataset: set of listed genome ids} from <dataset>.txt files."""
    lists = {}
    root = guard(folder, "list dir")
    need(root.is_dir(), "the list dir is not a directory")
    for path in sorted(root.glob("*.txt")):
        text = guard(path, "list file").read_text(encoding="utf-8")
        lists[path.name[:-4]] = {x.strip() for x in text.splitlines() if x.strip()}
    return lists


def collect(design, results, lists):
    """Per clock and reading, the dataset values of each contrast, horizon and learner, the lacking rows, and every
    genome id seen. Stops on any failed check."""
    values = {}   # (clock, reading, contrast, h) -> list of dataset rows
    lacking, ids = [], set()
    for lst in (lists or {}).values():
        ids |= lst
    for d in design["datasets"]:
        if d["why_not"]:
            continue
        clock = d["clock"]
        minimum = design["clocks"][clock]["test_minimum"]
        listed = (lists or {}).get(d["id"], set())
        loaded = {}
        windows = {w["window"]: w for w in d["windows"]}
        fits = {(w["window"], f["origin"], f["arm"]): f for w in d["windows"] for f in w["fits"]}
        for (j, k, arm), fit in sorted(fits.items()):
            if fit["file"] not in loaded:
                loaded[fit["file"]] = load_result(results, d, fit, windows[j])
            samples, _, fold = loaded[fit["file"]]
            need(fold["train"] == fit["train"] and fold["train_resistant"] == fit["train_resistant"],
                 f"{fit['file']}: training counts other than the design's for {arm}({k}) on window {j}")
            ids.update(g for g, _ in samples)
        for contrast in CONTRASTS:
            arm = ARM_OF[contrast]
            by_h = {}
            for p in d["pairs"][contrast]:
                j, k = p["window"], p["older"]
                older, fresh = loaded[fits[(j, k, arm)]["file"]], loaded[fits[(j, j, arm)]["file"]]
                need([g for g, _ in older[0]] == [g for g, _ in fresh[0]],
                     f"{d['id']} {clock} {contrast} {k} and {j}: the two fits scored other genomes")
                need(older[0] == fresh[0], f"{d['id']} {clock} {contrast} {k} and {j}: the same genomes, other labels")
                by_h.setdefault(p["horizon"], []).append((older, fresh))
            for h, pairs in sorted(by_h.items()):
                for reading in READINGS:
                    if reading == "second" and lists is None:
                        continue
                    kept, reasons = [], []
                    for older, fresh in pairs:
                        keep, why = (set(range(len(fresh[0]))), "") if reading == "first" else \
                            second_case(fresh[0], listed, minimum)
                        if keep is None:
                            reasons.append(why)
                        else:
                            kept.append((older, fresh, keep))
                    if not kept:
                        lacking.append({"clock": clock, "contrast": contrast, "horizon": h, "dataset": d["id"],
                                        "species": d["species_slug"],
                                        "reason": "; ".join(r for r in (SHORT, UNDER) if r in reasons)})
                        continue
                    values.setdefault((clock, reading, contrast, h), []).append(
                        dataset_rows(d, kept))
    return values, lacking, ids


def dataset_rows(d, kept):
    """One dataset's mean over its kept pairs at a horizon: per learner the older and fresh measures and the
    difference, with the test genomes the pairs scored."""
    out = {"dataset": d["id"], "species": d["species_slug"], "pairs": len(kept),
           "test": sum(len(keep) for _, _, keep in kept),
           "test_resistant": sum(sum(1 for i in keep if fresh[0][i][1]) for _, fresh, keep in kept), "learners": {}}
    for l in LEARNERS:
        older = [measure(o[0], o[1][l], keep) for o, _, keep in kept]
        fresh = [measure(f[0], f[1][l], keep) for _, f, keep in kept]
        mean = lambda xs: math.fsum(xs) / len(xs)
        out["learners"][l] = {
            "older": {m: mean([x[m] for x in older]) for m in MEASURES},
            "fresh": {m: mean([x[m] for x in fresh]) for m in MEASURES},
            "d": {m: mean([o[m] - f[m] for o, f in zip(older, fresh)]) for m in MEASURES}}
    return out


def analyze(design, values, readings_made):
    """Readings, learner rules and verdicts per clock and reading, drawn in the documented order."""
    reading_rows, rule_rows, verdict_rows = [], [], []
    verdicts, rules, reads = {}, {}, {}
    for clock in CLOCKS:
        descriptive_clock = design["clocks"][clock]["descriptive_only"]
        horizons = range(1, len(design["clocks"][clock]["shares"]))
        for reading in readings_made:
            rng = random.Random(SEED)
            got = {}
            for l in LEARNERS:
                for contrast in CONTRASTS:
                    for h in horizons:
                        rows = values.get((clock, reading, contrast, h), [])
                        if not rows:
                            continue
                        for m in MEASURES:
                            by_species = {}
                            for r in rows:
                                by_species.setdefault(r["species"], []).append(r["learners"][l]["d"][m])
                            point, lo, hi = interval(by_species, rng)
                            vs = [v for s in by_species.values() for v in s]
                            got[(contrast, h, l, m)] = read(lo, hi, m)
                            reading_rows.append({
                                "clock": clock, "reading": reading, "contrast": contrast, "horizon": h, "learner": l,
                                "measure": m, "datasets": len(rows), "species": len(by_species),
                                "median": f"{point:.6f}", "low": f"{lo:.6f}", "high": f"{hi:.6f}",
                                "read_as": got[(contrast, h, l, m)],
                                "datasets_worse": sum(1 for v in vs if v * WORSE[m] > 0),
                                "datasets_better": sum(1 for v in vs if v * WORSE[m] < 0),
                                "datasets_equal": sum(1 for v in vs if v == 0)})
            for contrast in CONTRASTS:
                for h in horizons:
                    rows = values.get((clock, reading, contrast, h), [])
                    if not rows:
                        continue
                    lr, lreads = {}, {}
                    for l in LEARNERS:
                        lreads[l] = tuple(got[(contrast, h, l, m)] for m in RULE_MEASURES)
                        lr[l] = learner_rule(*lreads[l])
                        rule_rows.append({"clock": clock, "reading": reading, "contrast": contrast, "horizon": h,
                                          "learner": l, "recall": lreads[l][0], "log_loss": lreads[l][1],
                                          "rule": lr[l]})
                        rules[(reading, clock, contrast, h, l)] = lr[l]
                        reads[(reading, clock, contrast, h, l)] = lreads[l]
                    n_species = len({r["species"] for r in rows})
                    v = verdict(len(rows), n_species, lr, lreads, descriptive_clock)
                    verdicts[(reading, clock, contrast, h)] = v
                    verdict_rows.append({"clock": clock, "reading": reading, "contrast": contrast, "horizon": h,
                                         "datasets": len(rows), "species": n_species,
                                         "learners_with_cost": sum(x == "cost" for x in lr.values()),
                                         "learners_better": sum("better" in lreads[x] for x in LEARNERS),
                                         "verdict": v})
    return reading_rows, rule_rows, verdict_rows, verdicts, rules, reads


def stands(design, verdicts, readings_made):
    rows = []
    for clock in CLOCKS:
        for contrast in CONTRASTS:
            for h in range(1, len(design["clocks"][clock]["shares"])):
                first = verdicts.get(("first", clock, contrast, h))
                if first is None:
                    continue
                second = verdicts.get(("second", clock, contrast, h), "no dataset") if "second" in readings_made \
                    else "not read"
                rows.append({"clock": clock, "contrast": contrast, "horizon": h, "first": first, "second": second,
                             "stands": first if first == second else (
                                 "second reading not made" if second == "not read" else "no verdict stands")})
    return rows


def predictions(verdicts, rules, reads, readings_made):
    rows = []
    for pid, _, _, _, text in PREDICTIONS:
        got = {}
        for reading in readings_made:
            v = {k[1:]: x for k, x in verdicts.items() if k[0] == reading}
            ru = {k[1:]: x for k, x in rules.items() if k[0] == reading}
            re_ = {k[1:]: x for k, x in reads.items() if k[0] == reading}
            got[reading] = "met" if prediction_met(pid, v, ru, re_) else "not met"
        rows.append({"prediction": pid, "rule": text, "first": got["first"],
                     "second": got.get("second", "not read"),
                     "outcome": "met" if all(x == "met" for x in got.values()) else "not met"})
    return rows


def tsv(columns, rows):
    return "\t".join(columns) + "\n" + "".join("\t".join(str(r[c]) for c in columns) + "\n" for r in rows)


def arm_and_contrast_rows(values):
    arms, contrasts = [], []
    for (clock, reading, contrast, h), rows in sorted(values.items(), key=lambda kv: (
            CLOCKS.index(kv[0][0]), READINGS.index(kv[0][1]), CONTRASTS.index(kv[0][2]), kv[0][3])):
        for r in rows:
            for l in LEARNERS:
                base = {"clock": clock, "reading": reading, "dataset": r["dataset"], "species": r["species"],
                        "contrast": contrast, "horizon": h, "learner": l, "pairs": r["pairs"], "test": r["test"],
                        "test_resistant": r["test_resistant"]}
                for arm in ("older", "fresh"):
                    arms.append({**base, "arm": arm,
                                 **{m: f"{r['learners'][l][arm][m]:.6f}" for m in MEASURES}})
                contrasts.append({**base, **{"d_" + m: f"{r['learners'][l]['d'][m]:.6f}" for m in MEASURES}})
    return arms, contrasts


def scan(files, ids):
    """Stop when a genome id, or the shape of an id, an accession or a sequence, stands in any text to write."""
    for name, text in files.items():
        tokens = set(re.findall(r"[A-Za-z0-9_.:-]+", text))
        need(not tokens & ids, f"{name}: a genome id stands in it; nothing is written")
        for shape in SHAPES:
            need(re.search(shape, text) is None, f"{name}: the shape of an id or a sequence stands in it")


def build(design_dir, results, listed_dir=None):
    """{file name: text} and the lines to print."""
    design = design_mod.read_json(Path(design_dir) / "design.json", "design")
    need(design.get("plan_sha256") == design_mod.PLAN_SHA256 and design.get("key") == design_mod.KEY,
         "design.json: not a design of the frozen plan")
    results = guard(results, "result dir")
    lists = read_lists(listed_dir) if listed_dir is not None else None
    readings_made = READINGS if lists is not None else READINGS[:1]
    values, lacking, ids = collect(design, results, lists)
    need(values, "no dataset with a pair to read")
    reading_rows, rule_rows, verdict_rows, verdicts, rules, reads = analyze(design, values, readings_made)
    stand_rows = stands(design, verdicts, readings_made)
    prediction_rows = predictions(verdicts, rules, reads, readings_made)
    arms, contrasts = arm_and_contrast_rows(values)
    files = {"arms.tsv": tsv(ARM_COLUMNS, arms), "contrasts.tsv": tsv(CONTRAST_COLUMNS, contrasts),
             "readings.tsv": tsv(READING_COLUMNS, reading_rows), "rule.tsv": tsv(RULE_COLUMNS, rule_rows),
             "verdicts.tsv": tsv(VERDICT_COLUMNS, verdict_rows), "stands.tsv": tsv(STAND_COLUMNS, stand_rows),
             "lacking.tsv": tsv(LACK_COLUMNS, lacking), "predictions.tsv": tsv(PREDICTION_COLUMNS, prediction_rows)}
    scan(files, ids)
    lines = [f"{r['clock']} {r['contrast']} horizon {r['horizon']}: first {r['first']}, second {r['second']}; "
             f"stands: {r['stands']}" for r in stand_rows]
    lines += [f"{r['prediction']} ({r['rule']}): {r['outcome']}" for r in prediction_rows]
    if lists is not None:
        named = {x["id"] for x in design["datasets"]}
        lines.append(f"lists: {len(lists)} files, {sum(1 for d in lists if d in named)} of them datasets of the design")
    return files, lines


def main(argv):
    args, listed = list(argv), None
    if "--listed" in args:
        i = args.index("--listed")
        listed = args[i + 1] if i + 1 < len(args) else None
        del args[i:i + 2]
    if len(args) != 3 or ("--listed" in argv and listed is None):
        print(__doc__.strip().split("\n\n")[1], file=sys.stderr)
        return 2
    try:
        out = design_mod.empty_out(args[2])
        files, lines = build(args[0], args[1], listed)
        out.mkdir(parents=True, exist_ok=True)
        for name, text in files.items():
            (out / name).write_text(text, encoding="utf-8")
        print("\n".join(lines))
        return 0
    except Stop as e:
        print(f"stagger_analyze: stopped: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
