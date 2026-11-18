"""Cells of a study-split record that must agree with each other, checked without the inputs.

problems(directory) gives one line per cell that does not follow from the
cells beside it or from another file of the record. It reads the row files,
plan copies, checksums and README; the derived files are held to derive() by
test_record.py.
"""

import hashlib
import math
import pathlib
import re
from collections import Counter, defaultdict
from fractions import Fraction

from .common import load_script

record = load_script()
STAMP = re.compile(r"^\d{4}-\d\d-\d\dT\d\d:\d\d:\d\d(\.\d\d)?Z$")


def when(stamp):
    """A time of order.tsv as (date, seconds), so that 00.58 sorts after 00; None if it is not one."""
    if not STAMP.match(stamp):
        return None
    day, clock = stamp[:-1].split("T")
    h, m, s = clock.split(":")
    return day, Fraction(int(h) * 3600 + int(m) * 60) + Fraction(s)


def load(directory):
    d = pathlib.Path(directory)
    files = {p.relative_to(d).as_posix(): p.read_bytes() for p in sorted(d.rglob("*")) if p.is_file()}
    tables = {}
    for name in record.COLUMNS:
        if name in files:
            tables[name] = record.parse(files[name].decode("utf-8"), name)
    return files, tables


def checksum_problems(files):
    out = []
    sums = {}
    for line in files["SHA256SUMS"].decode("utf-8").splitlines():
        digest, name = line.split("  ", 1)
        sums[name] = digest
    if set(sums) != set(files) - {"SHA256SUMS"}:
        out.append(f"SHA256SUMS names {sorted(set(sums) ^ (set(files) - {'SHA256SUMS'}))}")
    for name, digest in sums.items():
        if name in files and hashlib.sha256(files[name]).hexdigest() != digest:
            out.append(f"{name}: sha256 differs from SHA256SUMS")
    return out


def fit_problems(rows, where, ds):
    out = []
    for r in rows:
        tp, fn, tn, fp = (int(r[k]) for k in ("tp", "fn", "tn", "fp"))
        cells = record.fit_cells(tp, fn, tn, fp, record.num(r["log_loss"]))
        have = [r[k] for k in ("tp", "fn", "tn", "fp", "recall", "specificity", "balanced_accuracy", "log_loss",
                               "error_rate")]
        key = (r["dataset"], r.get("plan", "lineage"), r["learner"])
        if have != cells:
            out.append(f"{where} {key}: cells differ from the confusion counts")
        d = ds.get(r["dataset"])
        if d is None or (tp + fn, tn + fp) != (int(d["resistant"]), int(d["genomes"]) - int(d["resistant"])):
            out.append(f"{where} {key}: class counts differ from datasets.tsv")
        ll = float(r["log_loss"])
        lo, hi = -math.log(1 - record.CLIP), -math.log(record.CLIP)
        if not lo - 1e-12 <= ll <= hi + 1e-12:
            out.append(f"{where} {key}: log loss outside the clipped range")
        if r["learner"] == "scm":
            # A 0 or 1 score: each error loses -log(clip), each right call -log(1 - clip).
            e = Fraction(fn + fp, tp + fn + tn + fp)
            if abs(ll - (float(e) * hi + (1 - float(e)) * lo)) > 1e-9:
                out.append(f"{where} {key}: scm log loss differs from its error rate")
        if r["null_fit"] not in ("yes", "no") or (r["null_fit"] == "yes" and r["learner"] not in ("gbt", "lasso")):
            out.append(f"{where} {key}: null fit flag")
    return out


def problems(directory):
    files, t = load(directory)
    out = checksum_problems(files)
    S = {r["key"]: r["value"] for r in t["settings.tsv"]}
    ds = {r["dataset"]: r for r in t["datasets.tsv"]}
    n = {d: int(r["genomes"]) for d, r in ds.items()}
    res = {d: int(r["resistant"]) for d, r in ds.items()}
    # settings: the script's constants
    want = {"cutoff": record.fmt(record.CUTOFF), "clip": record.fmt(record.CLIP), "folds": str(record.FOLDS),
            "draws": str(len(record.DRAWS)), "threshold_SR": str(record.THRESHOLD["SR"]),
            "threshold_SM": str(record.THRESHOLD["SM"]), "threshold_MR": str(record.THRESHOLD["MR"]),
            "experiment_rule_learners": str(record.EXPERIMENT_LEARNERS), "no_study_min": str(record.NO_STUDY_MIN),
            "study_rule_line": record.RULE_LINE}
    for k, v in want.items():
        if S.get(k) != v:
            out.append(f"settings.tsv {k}: {S.get(k)} against the script's {v}")
    rule = dict(x.split(":") for x in S["study_rule"].split(","))
    # datasets.tsv
    for d, r in ds.items():
        if d != f"{r['species']}__" + re.sub(r"[^a-z0-9]+", "_", r["antibiotic"]):
            out.append(f"datasets.tsv {d}: species and antibiotic differ from the id")
        if r["largest_study_share"] != record.fmt(Fraction(int(r["largest_study"]), n[d])):
            out.append(f"datasets.tsv {d}: largest study share")
        if not 0 < res[d] < n[d] or int(r["no_study"]) > n[d] or int(r["largest_study"]) > n[d]:
            out.append(f"datasets.tsv {d}: counts out of range")
        if r["lineage_split"] not in ("scored", "unavailable"):
            out.append(f"datasets.tsv {d}: lineage split")
    # eligibility.tsv and variants.tsv
    ms, sh, gr = int(S["min_studies"]), Fraction(S["min_study_share"]), Fraction(S["max_group_share"])
    entered = set()
    for r in t["eligibility.tsv"]:
        gen, ws, st, lg = (int(r[k]) for k in ("genomes", "with_study", "studies", "largest_group"))
        flags = [st < ms, Fraction(ws, gen) < sh, Fraction(lg, gen) > gr]
        passes = not any(flags)
        have = [r["under_studies"], r["under_share"], r["over_group"], r["passes"]]
        if have != [record.fmt(x) for x in flags + [passes]]:
            out.append(f"eligibility.tsv {r['dataset']}: flags differ from the counts")
        if r["entered"] == "yes":
            entered.add(r["dataset"])
            if not passes or r["folds"] != "-":
                out.append(f"eligibility.tsv {r['dataset']}: entered against its counts or fold check")
        elif passes and r["folds"] == "-":
            out.append(f"eligibility.tsv {r['dataset']}: passes, not entered, and no fold check reason")
        if r["dataset"] in ds and (gen, st, ws) != (n[r["dataset"]], int(ds[r["dataset"]]["studies"]),
                                                    n[r["dataset"]] - int(ds[r["dataset"]]["no_study"])):
            out.append(f"eligibility.tsv {r['dataset']}: genomes or BioProjects differ from datasets.tsv")
    if entered != set(ds):
        out.append("eligibility.tsv: the datasets entered differ from datasets.tsv")
    el = t["eligibility.tsv"]
    for v in t["variants.tsv"]:
        a, b, c = int(v["min_studies"]), Fraction(v["min_study_share"]), Fraction(v["max_group_share"])
        keep = [r for r in el if int(r["studies"]) >= a and Fraction(int(r["with_study"]), int(r["genomes"])) >= b
                and Fraction(int(r["largest_group"]), int(r["genomes"])) <= c]
        got = [str(len(keep)), str(len({r["species"] for r in keep})), str(sum(int(r["genomes"]) for r in keep))]
        if got != [v["datasets"], v["species"], v["pairs"]]:
            out.append(f"variants.tsv {a} {b} {c}: {got} from eligibility.tsv")
    # plans.tsv
    plans = {(r["dataset"], r["plan"]): r for r in t["plans.tsv"]}
    names = {f"{k}-d{i}" for k in record.KINDS for i in record.DRAWS}
    for d in ds:
        if {p for (x, p) in plans if x == d} != names:
            out.append(f"plans.tsv {d}: plans")
    for (d, p), r in plans.items():
        kind, draw = p.split("-d")
        sizes = [int(r[f"size_{f}"]) for f in range(record.FOLDS)]
        rs = [int(r[f"resistant_{f}"]) for f in range(record.FOLDS)]
        if (r["kind"], r["draw"]) != (kind, draw):
            out.append(f"plans.tsv {d} {p}: kind or draw")
        if r["rule"] != ("-" if kind == "random" else rule[f"d{draw}"]):
            out.append(f"plans.tsv {d} {p}: rule differs from settings.tsv")
        if d in n and (sum(sizes) != n[d] or sum(rs) != res[d]):
            out.append(f"plans.tsv {d} {p}: fold counts do not add up to the dataset")
        if d in n and not all(0 < a < s and 0 < res[d] - a < n[d] - s for a, s in zip(rs, sizes)):
            out.append(f"plans.tsv {d} {p}: a fold lacks a class")
        if not re.fullmatch(r"[0-9a-f]{64}", r["fold_sha256"]):
            out.append(f"plans.tsv {d} {p}: fold sha256")
        if kind == "study" and r["study_mates"] != "0":
            out.append(f"plans.tsv {d} {p}: a study plan with study mates in training")
        if d in n:
            for k in ("study", "lineage"):
                if r[f"{k}_mate_share"] != record.fmt(Fraction(int(r[f"{k}_mates"]), n[d])):
                    out.append(f"plans.tsv {d} {p}: {k} mate share")
        if kind == "matched":
            s = plans.get((d, f"study-d{draw}"))
            if s is None or [s[f"size_{f}"] for f in range(5)] + [s[f"resistant_{f}"] for f in range(5)] != \
                    [r[f"size_{f}"] for f in range(5)] + [r[f"resistant_{f}"] for f in range(5)]:
                out.append(f"plans.tsv {d} {p}: fold counts differ from study draw {draw}")
    # fits.tsv and lineage.tsv
    out += fit_problems(t["fits.tsv"], "fits.tsv", ds)
    out += fit_problems(t["lineage.tsv"], "lineage.tsv", ds)
    fits = {(r["dataset"], r["plan"], r["learner"]): r for r in t["fits.tsv"]}
    want = {(d, p, l) for d in ds for p in names for l in record.LEARNERS}
    if set(fits) != want:
        out.append("fits.tsv: not one row per dataset, plan and learner")
    lin = {r["dataset"] for r in t["lineage.tsv"]}
    if lin != {d for d, r in ds.items() if r["lineage_split"] == "scored"}:
        out.append("lineage.tsv: datasets differ from datasets.tsv")
    if Counter(r["dataset"] for r in t["lineage.tsv"]) != Counter({d: len(record.LEARNERS) for d in lin}):
        out.append("lineage.tsv: not one row per learner")
    # subsets.tsv
    sub = defaultdict(dict)
    for r in t["subsets.tsv"]:
        sub[(r["dataset"], r["plan"], r["learner"])][r["subset"]] = r
    ns = {d for d, r in ds.items() if int(r["no_study"]) >= int(S["no_study_min"])}
    if {k[0] for k in sub} != ns or len(sub) != len(ns) * len(names) * len(record.LEARNERS):
        out.append("subsets.tsv: datasets differ from those with enough genomes without a BioProject")
    for key, parts in sub.items():
        f = fits.get(key)
        if set(parts) != {"study", "no-study"} or f is None:
            out.append(f"subsets.tsv {key}: subsets")
            continue
        for c in ("n", "resistant", "tp", "fn", "tn", "fp"):
            total = sum(int(x[c]) for x in parts.values())
            ref = n[key[0]] if c == "n" else res[key[0]] if c == "resistant" else int(f[c])
            if total != ref:
                out.append(f"subsets.tsv {key}: {c} does not add up to fits.tsv")
        for x in parts.values():
            if int(x["tp"]) + int(x["fn"]) != int(x["resistant"]) or \
                    sum(int(x[c]) for c in ("tp", "fn", "tn", "fp")) != int(x["n"]):
                out.append(f"subsets.tsv {key} {x['subset']}: counts")
        if int(parts["no-study"]["n"]) != int(ds[key[0]]["no_study"]):
            out.append(f"subsets.tsv {key}: genomes without a BioProject differ from datasets.tsv")
        pooled = sum(Fraction(float(x["log_loss"])) * int(x["n"]) for x in parts.values()) / n[key[0]]
        if abs(float(pooled) - float(f["log_loss"])) > 1e-12:
            out.append(f"subsets.tsv {key}: log loss does not pool to fits.tsv")
    # folds.tsv: the first experiment per study fold, pooled back to its two fits
    fold_rows = defaultdict(list)
    for r in t["folds.tsv"]:
        fold_rows[(r["dataset"], r["learner"])].append(r)
    if set(fold_rows) != {(d, l) for d in ds for l in record.LEARNERS}:
        out.append("folds.tsv: not every dataset and learner")
    for (d, l), rows in fold_rows.items():
        if d not in ds:
            continue
        p = plans[(d, "study-d0")]
        if [r["test"] for r in sorted(rows, key=lambda r: int(r["fold"]))] != [p[f"size_{f}"] for f in range(5)]:
            out.append(f"folds.tsv {d} {l}: test sizes differ from study draw 0")
        for arm, plan in (("study", "study-d0"), ("random", "random-d0")):
            f = fits[(d, plan, l)]
            for c in ("tp", "fn", "tn", "fp"):
                if sum(int(r[f"{arm}_{c}"]) for r in rows) != int(f[c]):
                    out.append(f"folds.tsv {d} {l}: {arm} {c} does not add up to fits.tsv")
            pooled = sum(Fraction(float(r[f"{arm}_log_loss"])) * int(r["test"]) for r in rows) / n[d]
            if abs(float(pooled) - float(f["log_loss"])) > 1e-12:
                out.append(f"folds.tsv {d} {l}: {arm} log loss does not pool to fits.tsv")
        for r in rows:
            if int(r["test_resistant"]) != int(r["study_tp"]) + int(r["study_fn"]) or \
                    int(r["test"]) != sum(int(r[f"random_{c}"]) for c in ("tp", "fn", "tn", "fp")) or \
                    int(r["test_resistant"]) != int(p[f"resistant_{r['fold']}"]):
                out.append(f"folds.tsv {d} {l} fold {r['fold']}: counts")
            rsz = [int(plans[(d, "random-d0")][f"size_{f}"]) for f in range(5)]
            if not int(r["test"]) * (n[d] - max(rsz)) <= int(r["train_random_sum"]) <= int(r["test"]) * (n[d] - min(rsz)):
                out.append(f"folds.tsv {d} {l} fold {r['fold']}: training sum out of range")
    # distances.tsv and variation.tsv
    dist = Counter((r["dataset"], r["kind"]) for r in t["distances.tsv"])
    if dist != Counter({(d, k): 10 for d in ds for k in record.KINDS}):
        out.append("distances.tsv: not 10 pairs of draws per dataset and kind")
    for r in t["distances.tsv"]:
        if r["dataset"] in n and not 0 <= int(r["changed"]) <= n[r["dataset"]] or int(r["draw_a"]) >= int(r["draw_b"]):
            out.append(f"distances.tsv {r['dataset']} {r['kind']} {r['draw_a']},{r['draw_b']}: out of range")
    if Counter(r["dataset"] for r in t["variation.tsv"]) != Counter({d: 8 for d in ds}):
        out.append("variation.tsv: not draws 1 to 4 under rules A and B per dataset")
    for r in t["variation.tsv"]:
        if r["dataset"] in n and not 0 <= int(r["changed"]) <= n[r["dataset"]]:
            out.append(f"variation.tsv {r['dataset']}: out of range")
    # order.tsv: the times the page puts in order
    o = {(r["experiment"], r["event"]): r["time"] for r in t["order.tsv"]}
    for (e, k), v in o.items():
        if not STAMP.match(v):
            out.append(f"order.tsv {e} {k}: time")
    chains = [
        ("first", ["fold_builder", "fold_files_first", "fold_files_last", "variant_table", "hash_time",
                   "first_fit_end", "last_fit_end"]),
        ("second", ["plan_builder", "stated_rule_choice", "hash_time", "fits_start", "first_fit_end",
                    "last_fit_end"]),
    ]
    for e, chain in chains:
        times = [when(o.get((e, k), "")) for k in chain]
        if None in times or times != sorted(times):
            out.append(f"order.tsv {e}: times out of order")
    if (when(o.get(("second", "rule_variation"), "")) or 0) < (when(o.get(("second", "plan_builder"), "")) or 1):
        out.append("order.tsv second: the movement table before the plan builder")
    # planfiles.tsv: each plan's sha256 is its line of inputs.sha256, and its hash time is order.tsv's
    listed = dict(reversed(x.split("  ", 1)) for x in files["inputs.sha256"].decode("utf-8").splitlines())
    for r in t["planfiles.tsv"]:
        if r["file"] != f"{r['experiment']}/PLAN.md" or listed.get(r["file"]) != r["sha256"]:
            out.append(f"planfiles.tsv {r['experiment']}: the plan's sha256 differs from inputs.sha256")
        if r["hash_time"] != o.get((r["experiment"], "hash_time")) or r["listed"] != "yes":
            out.append(f"planfiles.tsv {r['experiment']}: hash time differs from order.tsv")
    # checks.tsv: the counts that follow from the row files
    ck = {r["check"]: int(r["count"]) for r in t["checks.tsv"]}
    kinds = Counter(r["kind"] for r in t["plans.tsv"])
    want = {"fold_files_checked": len(plans), "copy_groups_whole_plans": len(plans),
            "study_plans_with_no_study_split": kinds["study"], "random_d0_equal_release_folds": len(ds),
            "matched_plans_equal_study_counts": kinds["matched"], "lineage_splits_clusters_whole": len(lin),
            "learner_fits_equal_evaluator_counts": len(t["fits.tsv"]) + len(t["lineage.tsv"]),
            "scm_fits_scores_0_or_1": len(plans), "genome_ids_in_record": 0}
    for k, v in want.items():
        if ck.get(k) != v:
            out.append(f"checks.tsv {k}: {ck.get(k)} against {v} from the row files")
    # totals.tsv: distinct genomes and BioProjects over the datasets
    tot = {r["key"]: int(r["value"]) for r in t["totals.tsv"]}
    if not max(n.values()) <= tot["genomes"] <= sum(n.values()):
        out.append("totals.tsv genomes: out of range")
    st = [int(r["studies"]) for r in ds.values()]
    if not max(st) <= tot["studies"] <= sum(st):
        out.append("totals.tsv studies: out of range")
    # inputs.sha256: sorted relative paths
    lines = files["inputs.sha256"].decode("utf-8").splitlines()
    if lines != sorted(lines, key=lambda x: x.split("  ", 1)[1]):
        out.append("inputs.sha256: not sorted by path")
    for line in lines:
        digest, path = line.split("  ", 1)
        if not re.fullmatch(r"[0-9a-f]{64}", digest) or path.startswith("/") or ".." in path.split("/"):
            out.append(f"inputs.sha256 {path}: line")
    # README: what the script writes
    if files.get("README.md") != record.readme(files, S):
        out.append("README.md: differs from what the script writes")
    return out
