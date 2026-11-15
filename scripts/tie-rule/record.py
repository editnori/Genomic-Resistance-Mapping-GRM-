"""Writes reference/tie-rule: the record of the tie rule experiment of 15 November 2026.

The experiment (lane `tierule`) fitted 169 datasets of release rev14-274-20261112
in five arms with a build of branch `exp/tie-rule-exact`:

    A  `main`'s binary
    B  the branch's binary, flag off: the present rule
    C  B, traced
    D  the branch's binary, `--exp-tie-rule exact`: the exact rule
    E  D, traced

This script reads those arms' reports and traces and the release's index, and
fits nothing. It writes plain text files, sorted, the same bytes on every run:

    python3 -I scripts/tie-rule/record.py RUNS RELEASE PLAN_DIR --out reference/tie-rule
    python3 -I scripts/tie-rule/record.py RUNS RELEASE PLAN_DIR --check reference/tie-rule

RUNS holds A/ B/ C/ D/ E/ (`<id>.json`, and `<id>.trace.jsonl` in C and E),
RELEASE the unpacked release (`index.json`, `datasets/`), and PLAN_DIR the
experiment's `PLAN.md` and `PLAN.sha256`. `--check` prints the files that
differ and exits 1 unless every file is the text this script writes.

The record's row files (datasets, measures, reproduction, seconds, traces,
inputs) come from the inputs; every other file is derived from the row files
alone, so `derive(read_rows(DIR))` gives the derived files again from a
committed record (tests/tie_rule). No string of a report other than ids,
species slugs and learner names is written.
"""

import argparse
import hashlib
import json
import math
import pathlib
import random
import re
import statistics
import sys

LEARNERS = ("scm", "cart", "rf", "gbt", "lasso")
TREES = ("cart", "rf", "gbt")
SPLITS = ("lineage", "random")
KINDS = ("cart", "rf", "gbt", "gbt-mic")
TALLY = ("searched", "exact_tied", "differ_tie", "differ_order", "differ_partition")
PAIRED_MEASURES = ("recall", "specificity", "balanced_accuracy", "log_loss")
THRESHOLD = 0.5
CLAMP = 1e-12
DRAWS = 4000
SEED = 20261115
SHARES = (0.01, 0.05)
SECONDS = re.compile(rb'"seconds":\s*-?[0-9][0-9.eE+-]*')

# The commits and binaries the arms ran with, as the experiment's PLAN.md and
# REPORT.md give them, and the settings this script uses.
SETTINGS = (
    ("arm_A", "main's binary, built at commit 852ebc91c3f649d3980b1d9a6844e1daec5df2d2"),
    ("arm_B", "the branch's binary, flag off (the present rule)"),
    ("arm_C", "arm B with --exp-tie-trace"),
    ("arm_D", "the branch's binary, --exp-tie-rule exact (the exact rule)"),
    ("arm_E", "arm D with --exp-tie-trace"),
    ("binary_branch_sha256", "f2423e4f2186ef7a38602e79c5e6d484ad1373c03c85bb5577621f3afafaf8b5"),
    ("binary_main_sha256", "8ab7cbff1145cc0f4c4cc51dfa659fdf00bf51d6bf127a9d947d3c29c71fb5e9"),
    ("branch", "exp/tie-rule-exact"),
    ("branch_commit_head", "35e1ae246a7f6806175ad66210f54670b967fd8f"),
    ("branch_commit_runs", "88af13052b8e66e34faea4172a8f63ddb81c91e9"),
    ("call_threshold", "0.5"),
    ("evaluation_revision", "14"),
    ("log_loss_clamp", "1e-12"),
    ("main_commit_runs", "852ebc91c3f649d3980b1d9a6844e1daec5df2d2"),
    ("paired_draws", str(DRAWS)),
    ("paired_interval", "the 101st and 3900th of the 4000 sorted draw means"),
    ("paired_resample", "species"),
    ("paired_seed", str(SEED)),
    ("release", "rev14-274-20261112"),
)

COLUMNS = {
    "datasets.tsv": ("dataset", "species", "genomes", "resistant", "features", "run"),
    "measures.tsv": (
        "dataset", "split", "learner", "species", "genomes",
        "recall_present", "recall_exact", "recall_change",
        "specificity_present", "specificity_exact", "specificity_change",
        "balanced_accuracy_present", "balanced_accuracy_exact", "balanced_accuracy_change",
        "log_loss_present", "log_loss_exact", "log_loss_change",
        "correct_present", "correct_exact", "calls_changed",
        "scores_changed", "metrics_changed", "fold_balanced_accuracy_changed",
    ),
    "seconds.tsv": ("dataset", "learner", "seconds_main", "seconds_present", "seconds_exact"),
    "traces.tsv": (
        "dataset", "kind", "fits",
        *(f"{k}_present" for k in TALLY), "fits_differing_present",
        *(f"{k}_exact" for k in TALLY), "fits_differing_exact",
        "trees_present", "trees_exact", "tree_positions", "trees_differ", "fit_mismatches",
    ),
    "reproduction.tsv": ("dataset", "main_equals_release", "present_equals_main", "traced_present_equals_present",
                         "traced_exact_equals_exact"),
    "inputs.sha256": None,
    "splits.tsv": (
        "dataset", "split", "best_present", "best_present_balanced_accuracy", "second_present_balanced_accuracy",
        "gap_present", "best_exact", "best_exact_balanced_accuracy", "best_changes", "tie_effect", "tie_effect_at_least_gap",
    ),
    "order.tsv": ("split", "first", "second", "datasets", "order_changes"),
    "paired.tsv": ("learner", "split", "measure", "datasets", "species", "mean", "low", "high", "higher", "lower", "equal"),
    "summary.tsv": ("key", "value"),
    "settings.tsv": ("key", "value"),
}
ROW_FILES = ("datasets.tsv", "measures.tsv", "reproduction.tsv", "seconds.tsv", "traces.tsv")
DERIVED_FILES = ("splits.tsv", "order.tsv", "paired.tsv", "summary.tsv")


def cell(v):
    """One value as the record writes it: an int, a float's shortest repr, a word, or empty for None."""
    if v is None:
        return ""
    if isinstance(v, bool):
        return "1" if v else "0"
    if isinstance(v, float):
        if math.isnan(v) or math.isinf(v):
            raise ValueError(f"not a finite value: {v}")
        return repr(v)
    return str(v)


def table(name, rows):
    cols = COLUMNS[name]
    lines = ["\t".join(cols)]
    for r in rows:
        lines.append("\t".join(cell(r[c]) for c in cols))
    return "\n".join(lines) + "\n"


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def masked(path):
    with open(path, "rb") as f:
        return SECONDS.sub(b'"seconds":0', f.read())


def log_loss(scores, labels):
    total, k = 0.0, 0
    for s, y in zip(scores, labels):
        if s is None:
            continue
        p = min(max(s, CLAMP), 1 - CLAMP)
        total += -(math.log(p) if y else math.log(1 - p))
        k += 1
    return total / k if k else None


def calls(scores):
    return [None if s is None else s >= THRESHOLD for s in scores]


def learner_rows(report):
    """{(split, learner): learner row} for the splits a report scored."""
    out = {}
    for s in report["splits"]:
        if s.get("unavailable_reason") is not None or not s.get("learners"):
            continue
        for row in s["learners"]:
            if row["learner"] in LEARNERS:
                out[(s["scheme"], row["learner"])] = row
    return out


def change(a, b):
    return None if a is None or b is None else b - a


# ---------------------------------------------------------------- row files


def run_ids(runs):
    """Datasets with a report in every arm and a trace in C and E."""
    ids = None
    for arm in "ABCDE":
        have = {p.name[: -len(".json")] for p in (runs / arm).glob("*.json") if not p.name.endswith(".trace.jsonl")}
        ids = have if ids is None else ids & have
    for arm in "CE":
        ids &= {p.name[: -len(".trace.jsonl")] for p in (runs / arm).glob("*.trace.jsonl")}
    return sorted(ids)


def collect(runs, release):
    runs, release = pathlib.Path(runs), pathlib.Path(release)
    index = json.loads((release / "index.json").read_text(encoding="utf-8"))
    ids = run_ids(runs)
    datasets = []
    for d in index["datasets"]:
        datasets.append(dict(dataset=d["id"], species=d["species_slug"], genomes=d["n"], resistant=d["n_resistant"],
                             features=d["features"], run=d["id"] in ids))
    datasets.sort(key=lambda r: r["dataset"])
    unknown = sorted(set(ids) - {r["dataset"] for r in datasets})
    if unknown:
        raise SystemExit(f"run datasets not in the release index: {unknown}")
    by_id = {r["dataset"]: r for r in datasets}

    measures, seconds, traces, inputs, reproduction = [], [], [], [], []
    for i in ids:
        paths = {arm: runs / arm / f"{i}.json" for arm in "ABCDE"}
        for arm in "ABD":
            inputs.append((f"{arm}/{i}.json", sha256(paths[arm])))
        for arm in "CE":
            inputs.append((f"{arm}/{i}.trace.jsonl", sha256(runs / arm / f"{i}.trace.jsonl")))
        a, b = masked(paths["A"]), masked(paths["B"])
        reproduction.append(dict(dataset=i, main_equals_release=a == masked(release / "datasets" / f"{i}.json"),
                                 present_equals_main=b == a, traced_present_equals_present=masked(paths["C"]) == b,
                                 traced_exact_equals_exact=masked(paths["E"]) == masked(paths["D"])))
        A, B, D = (json.loads(paths[arm].read_bytes()) for arm in "ABD")
        if B["counts"]["samples"] != by_id[i]["genomes"] or B["counts"]["features"] != by_id[i]["features"]:
            raise SystemExit(f"{i}: the report's counts differ from the release index")
        labels = [s["resistant"] for s in B["samples"]]
        ra, rb, rd = learner_rows(A), learner_rows(B), learner_rows(D)
        if set(rb) != set(rd) or set(ra) != set(rb):
            raise SystemExit(f"{i}: the arms scored different splits or learners")
        for (split, learner) in sorted(rb):
            p, e = rb[(split, learner)], rd[(split, learner)]
            mp, me = p["metrics"], e["metrics"]
            cp, ce = calls(p["scores"]), calls(e["scores"])
            measures.append(dict(
                dataset=i, split=split, learner=learner, species=B["species_slug"], genomes=mp["n"],
                recall_present=mp["sensitivity"], recall_exact=me["sensitivity"],
                recall_change=change(mp["sensitivity"], me["sensitivity"]),
                specificity_present=mp["specificity"], specificity_exact=me["specificity"],
                specificity_change=change(mp["specificity"], me["specificity"]),
                balanced_accuracy_present=mp["balanced_accuracy"], balanced_accuracy_exact=me["balanced_accuracy"],
                balanced_accuracy_change=change(mp["balanced_accuracy"], me["balanced_accuracy"]),
                log_loss_present=log_loss(p["scores"], labels), log_loss_exact=log_loss(e["scores"], labels),
                log_loss_change=None,
                correct_present=mp["tp"] + mp["tn"], correct_exact=me["tp"] + me["tn"],
                calls_changed=sum(x != y for x, y in zip(cp, ce)),
                scores_changed=p["scores"] != e["scores"], metrics_changed=mp != me,
                fold_balanced_accuracy_changed=p.get("fold_balanced_accuracy") != e.get("fold_balanced_accuracy"),
            ))
            measures[-1]["log_loss_change"] = change(measures[-1]["log_loss_present"], measures[-1]["log_loss_exact"])
        for learner in LEARNERS:
            sec = [sum(r[(s, learner)]["seconds"] for s in SPLITS if (s, learner) in r) for r in (ra, rb, rd)]
            seconds.append(dict(dataset=i, learner=learner, seconds_main=sec[0], seconds_present=sec[1], seconds_exact=sec[2]))
        traces.extend(trace_rows(i, runs / "C" / f"{i}.trace.jsonl", runs / "E" / f"{i}.trace.jsonl"))
    return dict(datasets=datasets, measures=measures, seconds=seconds, traces=traces, inputs=inputs, reproduction=reproduction)


def kind_of(fit):
    if fit["learner"] == "gbt":
        return "gbt" if fit.get("objective") == "logistic" else "gbt-mic"
    return fit["learner"]


def trace_rows(dataset, present_path, exact_path):
    fits = [[json.loads(line) for line in pathlib.Path(p).read_text(encoding="utf-8").splitlines()] for p in (present_path, exact_path)]
    rows = {k: dict(dataset=dataset, kind=k, fits=0, fits_differing_present=0, fits_differing_exact=0,
                    trees_present=0, trees_exact=0, tree_positions=0, trees_differ=0, fit_mismatches=0,
                    **{f"{t}_{arm}": 0 for t in TALLY for arm in ("present", "exact")}) for k in KINDS}
    if len(fits[0]) != len(fits[1]):
        raise SystemExit(f"{dataset}: the traces hold {len(fits[0])} and {len(fits[1])} fits")
    for c, e in zip(*fits):
        k = kind_of(c)
        r = rows[k]
        r["fits"] += 1
        if (c["learner"], c.get("objective"), c["seed"], c["fit_samples"]) != (e["learner"], e.get("objective"), e["seed"], e["fit_samples"]):
            r["fit_mismatches"] += 1
            continue
        for fit, arm in ((c, "present"), (e, "exact")):
            for t in TALLY:
                r[f"{t}_{arm}"] += fit["tally"][t]
            r[f"fits_differing_{arm}"] += any(fit["tally"][t] for t in ("differ_tie", "differ_order", "differ_partition"))
        r["trees_present"] += len(c["trees"])
        r["trees_exact"] += len(e["trees"])
        for j in range(max(len(c["trees"]), len(e["trees"]))):
            r["tree_positions"] += 1
            r["trees_differ"] += (c["trees"][j] if j < len(c["trees"]) else None) != (e["trees"][j] if j < len(e["trees"]) else None)
    return [rows[k] for k in KINDS]


# ---------------------------------------------------------------- reading a record back


def parse(v):
    if v == "":
        return None
    if re.fullmatch(r"-?[0-9]+", v):
        return int(v)
    try:
        return float(v)
    except ValueError:
        return v


def read_table(path):
    lines = pathlib.Path(path).read_text(encoding="utf-8").splitlines()
    cols = lines[0].split("\t")
    return [dict(zip(cols, (parse(v) for v in line.split("\t")))) for line in lines[1:]]


def read_rows(directory):
    d = pathlib.Path(directory)
    rows = {name[: -len(".tsv")]: read_table(d / name) for name in ROW_FILES}
    for r in rows["datasets"]:
        r["run"] = bool(r["run"])
    for r in rows["measures"]:
        for c in ("scores_changed", "metrics_changed", "fold_balanced_accuracy_changed"):
            r[c] = bool(r[c])
    for r in rows["seconds"]:
        for c in ("seconds_main", "seconds_present", "seconds_exact"):
            r[c] = float(r[c])
    for r in rows["reproduction"]:
        for c in COLUMNS["reproduction.tsv"][1:]:
            r[c] = bool(r[c])
    rows["inputs"] = [tuple(reversed(line.split("  ", 1))) for line in (d / "inputs.sha256").read_text(encoding="utf-8").splitlines()]
    return rows


# ---------------------------------------------------------------- derived files


def mean(v):
    return sum(v) / len(v)


def best(values):
    """The index's rule (src/report.rs `best_learner`): the highest balanced accuracy, a tie to the later learner."""
    top = None
    for learner in LEARNERS:
        v = values.get(learner)
        if v is not None and (top is None or v >= top[1]):
            top = (learner, v)
    return top


def sign(x):
    return (x > 0) - (x < 0)


def derive(rows):
    measures = rows["measures"]
    by = {}
    for r in measures:
        by.setdefault((r["learner"], r["split"]), []).append(r)
    for v in by.values():
        v.sort(key=lambda r: r["dataset"])
    per_dataset = {}
    for r in measures:
        per_dataset.setdefault((r["dataset"], r["split"]), {})[r["learner"]] = r

    splits = []
    for (dataset, split), m in sorted(per_dataset.items()):
        present = best({k: r["balanced_accuracy_present"] for k, r in m.items()})
        exact = best({k: r["balanced_accuracy_exact"] for k, r in m.items()})
        ranked = sorted((r["balanced_accuracy_present"] for r in m.values()), reverse=True)
        gap = ranked[0] - ranked[1]
        effect = max(abs(m[t]["balanced_accuracy_change"]) for t in TREES)
        splits.append(dict(dataset=dataset, split=split, best_present=present[0], best_present_balanced_accuracy=present[1],
                           second_present_balanced_accuracy=ranked[1], gap_present=gap, best_exact=exact[0],
                           best_exact_balanced_accuracy=exact[1], best_changes=present[0] != exact[0],
                           tie_effect=effect, tie_effect_at_least_gap=effect > 0 and effect >= gap))

    order = []
    for split in SPLITS:
        keys = sorted(k for k in per_dataset if k[1] == split)
        for a in range(len(LEARNERS)):
            for b in range(a + 1, len(LEARNERS)):
                x, y = LEARNERS[a], LEARNERS[b]
                flips = 0
                for k in keys:
                    m = per_dataset[k]
                    before = sign(m[x]["balanced_accuracy_present"] - m[y]["balanced_accuracy_present"])
                    after = sign(m[x]["balanced_accuracy_exact"] - m[y]["balanced_accuracy_exact"])
                    flips += before != after
                order.append(dict(split=split, first=x, second=y, datasets=len(keys), order_changes=flips))

    paired = []
    rng = random.Random(SEED)
    for learner in TREES:
        for split in SPLITS:
            rs = by[(learner, split)]
            groups = {}
            for r in rs:
                groups.setdefault(r["species"], []).append(r)
            species = sorted(groups)
            for measure in PAIRED_MEASURES:
                diffs = [r[f"{measure}_change"] for r in rs]
                draws = []
                for _ in range(DRAWS):
                    pick = [species[rng.randrange(len(species))] for _ in species]
                    v = [r[f"{measure}_change"] for s in pick for r in groups[s]]
                    draws.append(sum(v) / len(v))
                draws.sort()
                paired.append(dict(learner=learner, split=split, measure=measure, datasets=len(diffs), species=len(species),
                                   mean=mean(diffs), low=draws[int(0.025 * DRAWS)], high=draws[int(0.975 * DRAWS) - 1],
                                   higher=sum(x > 0 for x in diffs), lower=sum(x < 0 for x in diffs),
                                   equal=sum(x == 0 for x in diffs)))

    summary = summarize(rows, by, splits, order)
    return dict(splits=splits, order=order, paired=paired, summary=summary)


def summarize(rows, by, splits, order):
    s = {}
    datasets = rows["datasets"]
    run = [r for r in datasets if r["run"]]
    unrun = [r for r in datasets if not r["run"]]
    s["datasets.release"] = len(datasets)
    s["datasets.run"] = len(run)
    s["datasets.not_run"] = len(unrun)
    s["species.release"] = len({r["species"] for r in datasets})
    s["species.run"] = len({r["species"] for r in run})
    s["species.not_run"] = len({r["species"] for r in unrun} - {r["species"] for r in run})
    s["species.not_run.names"] = ",".join(sorted({r["species"] for r in unrun} - {r["species"] for r in run}))
    s["species.both"] = len({r["species"] for r in unrun} & {r["species"] for r in run})
    for name, group in (("run", run), ("not_run", unrun)):
        s[f"genomes.median.{name}"] = statistics.median(r["genomes"] for r in group)
        s[f"features.median.{name}"] = statistics.median(r["features"] for r in group)
    s["features.median.ratio"] = s["features.median.not_run"] / s["features.median.run"]

    for split in SPLITS:
        s[f"datasets.{split}"] = len(by[("cart", split)])
    for learner in LEARNERS:
        either_metric, either_ll, either_any = set(), set(), set()
        for split in SPLITS:
            rs = by[(learner, split)]
            p = f"{learner}.{split}"
            s[f"{p}.datasets"] = len(rs)
            s[f"{p}.genomes"] = sum(r["genomes"] for r in rs)
            s[f"{p}.calls_changed"] = sum(r["calls_changed"] for r in rs)
            s[f"{p}.metrics_changed"] = sum(r["metrics_changed"] for r in rs)
            s[f"{p}.scores_changed"] = sum(r["scores_changed"] for r in rs)
            s[f"{p}.net_correct"] = sum(r["correct_exact"] - r["correct_present"] for r in rs)
            either_metric |= {r["dataset"] for r in rs if r["metrics_changed"]}
            either_any |= {r["dataset"] for r in rs if r["metrics_changed"] or r["fold_balanced_accuracy_changed"]}
            for measure in PAIRED_MEASURES:
                d = [abs(r[f"{measure}_change"]) for r in rs]
                changed = [r for r in rs if r[f"{measure}_present"] != r[f"{measure}_exact"]]
                s[f"{p}.{measure}.changed"] = len(changed)
                s[f"{p}.{measure}.median_abs_change"] = statistics.median(d)
                top = max(range(len(d)), key=lambda j: (d[j], -j))
                s[f"{p}.{measure}.largest_abs_change"] = d[top]
                s[f"{p}.{measure}.largest_abs_change.dataset"] = rs[top]["dataset"]
                for t in SHARES:
                    s[f"{p}.{measure}.over_{t}"] = sum(x > t for x in d)
                if measure == "log_loss":
                    either_ll |= {r["dataset"] for r in changed}
        s[f"{learner}.either.metrics_changed"] = len(either_metric)
        s[f"{learner}.either.log_loss_changed"] = len(either_ll)
        s[f"{learner}.either.metrics_or_fold_changed"] = len(either_any)
        for measure in PAIRED_MEASURES:
            p = f"{learner}.all.{measure}"
            rs = [r for split in SPLITS for r in by[(learner, split)]]
            d = [abs(r[f"{measure}_change"]) for r in rs]
            s[f"{p}.largest_abs_change"] = max(d)
    trees_any = set()
    for t in TREES:
        trees_any |= {r["dataset"] for split in SPLITS for r in by[(t, split)] if r["metrics_changed"]}
    s["trees.either.metrics_changed"] = len(trees_any)
    ll_any = set()
    for t in TREES:
        ll_any |= {r["dataset"] for split in SPLITS for r in by[(t, split)] if r["log_loss_present"] != r["log_loss_exact"]}
    s["trees.either.log_loss_changed"] = len(ll_any)
    s["trees.calls_changed"] = sum(s[f"{t}.{split}.calls_changed"] for t in TREES for split in SPLITS)

    for split in SPLITS:
        rs = [r for r in splits if r["split"] == split]
        s[f"order.{split}.datasets"] = len(rs)
        s[f"order.{split}.best_changes"] = sum(r["best_changes"] for r in rs)
        s[f"order.{split}.tie_effect_at_least_gap"] = sum(r["tie_effect_at_least_gap"] for r in rs)
        s[f"order.{split}.median_gap"] = statistics.median(r["gap_present"] for r in rs)
        s[f"order.{split}.median_tie_effect"] = statistics.median(r["tie_effect"] for r in rs)
        s[f"order.{split}.tie_effect_positive"] = sum(r["tie_effect"] > 0 for r in rs)
        s[f"order.{split}.most_pair_changes"] = max(r["order_changes"] for r in order if r["split"] == split)

    for kind in KINDS:
        rs = [r for r in rows["traces"] if r["kind"] == kind]
        for c in COLUMNS["traces.tsv"][2:]:
            s[f"trace.{kind}.{c}"] = sum(r[c] for r in rs)
        s[f"trace.{kind}.datasets_differing_present"] = sum(r["fits_differing_present"] > 0 for r in rs)
        searched = s[f"trace.{kind}.searched_present"]
        if searched:
            s[f"trace.{kind}.exact_tied_share"] = s[f"trace.{kind}.exact_tied_present"] / searched
            differ = sum(s[f"trace.{kind}.{t}_present"] for t in ("differ_tie", "differ_order", "differ_partition"))
            s[f"trace.{kind}.differ_present"] = differ
            s[f"trace.{kind}.differ_share"] = differ / searched

    for learner in LEARNERS:
        rs = [r for r in rows["seconds"] if r["learner"] == learner]
        for arm in ("main", "present", "exact"):
            s[f"seconds.{learner}.{arm}"] = sum(r[f"seconds_{arm}"] for r in rs)
        s[f"seconds.{learner}.median_ratio_exact_present"] = statistics.median(
            r["seconds_exact"] / r["seconds_present"] for r in rs if r["seconds_present"] > 0)
        s[f"seconds.{learner}.median_ratio_present_main"] = statistics.median(
            r["seconds_present"] / r["seconds_main"] for r in rs if r["seconds_main"] > 0)

    for c in COLUMNS["reproduction.tsv"][1:]:
        s[f"reproduction.{c}"] = sum(r[c] for r in rows["reproduction"])
    s["inputs.files"] = len(rows["inputs"])
    return [dict(key=k, value=v) for k, v in sorted(s.items())]


# ---------------------------------------------------------------- writing


def render(rows, derived, plan_dir):
    files = {}
    files["datasets.tsv"] = table("datasets.tsv", rows["datasets"])
    files["measures.tsv"] = table("measures.tsv", sorted(rows["measures"], key=lambda r: (r["dataset"], r["split"], LEARNERS.index(r["learner"]))))
    files["seconds.tsv"] = table("seconds.tsv", sorted(rows["seconds"], key=lambda r: (r["dataset"], LEARNERS.index(r["learner"]))))
    files["traces.tsv"] = table("traces.tsv", sorted(rows["traces"], key=lambda r: (r["dataset"], KINDS.index(r["kind"]))))
    files["reproduction.tsv"] = table("reproduction.tsv", sorted(rows["reproduction"], key=lambda r: r["dataset"]))
    files["inputs.sha256"] = "".join(f"{h}  {p}\n" for p, h in sorted(rows["inputs"]))
    files["splits.tsv"] = table("splits.tsv", derived["splits"])
    files["order.tsv"] = table("order.tsv", derived["order"])
    files["paired.tsv"] = table("paired.tsv", derived["paired"])
    files["summary.tsv"] = table("summary.tsv", derived["summary"])
    files["settings.tsv"] = table("settings.tsv", [dict(key=k, value=v) for k, v in SETTINGS])
    plan = pathlib.Path(plan_dir)
    files["PLAN.md"] = (plan / "PLAN.md").read_bytes()
    files["PLAN.sha256"] = (plan / "PLAN.sha256").read_bytes()
    digest, name = files["PLAN.sha256"].decode("utf-8").split()[:2]
    if name != "PLAN.md" or hashlib.sha256(files["PLAN.md"]).hexdigest() != digest:
        raise SystemExit("PLAN.md does not have the sha256 PLAN.sha256 gives")
    return {k: v if isinstance(v, bytes) else v.encode("utf-8") for k, v in files.items()}


def sums(directory, files):
    """SHA256SUMS over every file of the directory but itself: the written ones and README.md."""
    d = pathlib.Path(directory)
    names = dict(files)
    readme = d / "README.md"
    if readme.exists():
        names["README.md"] = readme.read_bytes()
    return "".join(f"{hashlib.sha256(names[n]).hexdigest()}  {n}\n" for n in sorted(names)).encode("utf-8")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("runs")
    ap.add_argument("release")
    ap.add_argument("plan_dir")
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--out")
    g.add_argument("--check")
    a = ap.parse_args(argv)
    rows = collect(a.runs, a.release)
    files = render(rows, derive(rows), a.plan_dir)
    target = pathlib.Path(a.out or a.check)
    files["SHA256SUMS"] = sums(target, files)
    if a.check:
        bad = [n for n in sorted(files) if not (target / n).is_file() or (target / n).read_bytes() != files[n]]
        extra = sorted(p.name for p in target.iterdir() if p.name not in files and p.name != "README.md")
        for n in bad + extra:
            print(n)
        return 1 if bad or extra else 0
    target.mkdir(parents=True, exist_ok=True)
    for n, data in sorted(files.items()):
        (target / n).write_bytes(data)
    return 0


if __name__ == "__main__":
    sys.exit(main())
