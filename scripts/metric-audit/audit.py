#!/usr/bin/env python3
"""Read-only metric audit of GRM result reports, revisions 5 to 7.

From each report's saved per-sample labels and out-of-fold scores this
recomputes, over the called isolates, the pooled confusion matrices and
every derived metric, AUROC, AUPRC and ROC points, per-fold balanced
accuracy (outer folds rebuilt by foldgen with the pipeline's own seed and
rand crates, and compared with the saved test_fold), null-model folds, the
catalog baseline, catalog panels (from counted.json hits), probe-panel
no-call arithmetic, unsupported assays, MIC agreement with per-fold
breakpoints, eval_config, and the final forest's importance rows. It also
compares each file with the pre-run backup, the pipeline log, and the
published manifest. It writes only the output path it is given.

Schema-1 reports (revisions 5 and 6) score every isolate: a null-model fold
is one whose held-out scores are all exactly 0.5. Schema-2 reports
(revision 7) leave a null fit's isolates uncalled: a null-model fold is one
whose held-out scores are all null, and it has to be a gbt or lasso fold
with complexity 0. The best learner of a schema-2 report is the highest
random-split balanced accuracy among learners that call every isolate on
every split that ran, or none.

usage: python3 -I audit.py OUT.json [species-glob] --revision N [--root DIR]
           (--backup DIR --log FILE --launch UTC | --reshaped)
--revision is the evaluation revision the audited run writes. --backup is
the copy taken before launch (DIR/data/species/<slug>/results/<id>.json),
--log the run's pipeline log, and --launch its UTC start, such as
2026-10-24T12:43:12Z. --reshaped audits schema-2 reshapes of earlier
revisions, which have no run to compare with.
"""
import argparse
import base64
import datetime
import hashlib
import json
import math
import os
import re
import subprocess
import sys
import time
from collections import Counter, defaultdict
from fractions import Fraction

LANE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(LANE))
from result_scope import select_results

FOLDGEN = LANE + "/foldgen/target/release/foldgen"
ROOT = os.path.dirname(os.path.dirname(LANE))
REVISION = None
BACKUP = None
LOG = None
LAUNCH_EPOCH = None
RESHAPED = False
K, SEED = 5, 42
LEARNERS = ["scm", "cart", "rf", "gbt", "lasso"]
RANKED = {"cart", "rf", "gbt", "lasso"}
# src/probes.rs ASSAYS ids.
PROBE_IDS = {"xpert_rif", "xpert_xdr_inh", "xpert_xdr_eto", "xpert_xdr_fq", "xpert_xdr_amk",
             "gc_gyra", "ecoli_qrdr", "salmonella_qrdr"}
TOL = 1e-12
# EvalConfig as serialized: the forest and boosting seeds are derived per
# fit from `seed` and the outer fold, so the recorded config leaves them out.
EVAL_CONFIG_KEYS = {"folds", "inner_folds", "scm_max_rules", "cart_max_depth", "forest", "gbt",
                    "screened_features", "lasso_lambdas", "seed"}
FOREST_KEYS = {"n_trees", "subsample", "max_depth", "min_leaf"}
GBT_KEYS = {"rounds", "learning_rate", "max_depth", "lambda", "min_child_weight", "colsample", "patience"}
NO_WINNER = "no learner calls every isolate"
DONE = re.compile(r"\s+(\S+__\S+) done in (\S+) \((.*)\)\s*$")
FULL_COVERAGE = re.compile(r"best random-split learner among those that call every isolate: (\S+) (\S+)")
ANY_ROW = re.compile(r"best random-split balanced accuracy: (\S+) (\S+)")
RECOMPUTE = re.compile(r"\s+(\S+__\S+): (.+); recomputing report\s*$")


def options(description):
    """Parser for the run options the audit scripts share."""
    p = argparse.ArgumentParser(description=description)
    p.add_argument("--revision", type=int, required=True, help="evaluation revision the audited run writes")
    p.add_argument("--root", default=ROOT, help="directory holding data/ (default %(default)s)")
    p.add_argument("--backup", help="pre-launch copy holding data/species/<slug>/results/")
    p.add_argument("--log", help="the run's pipeline log")
    p.add_argument("--launch", help="UTC launch time, such as 2026-10-24T12:43:12Z")
    p.add_argument("--reshaped", action="store_true", help="schema-2 reshapes of earlier revisions; no run provenance")
    p.add_argument("--foldgen", default=FOLDGEN)
    return p


def configure(args, need=("backup", "log", "launch")):
    """Set the module's run settings from parsed options(); `need` names the
    provenance options a script cannot run without, unless --reshaped."""
    global ROOT, REVISION, BACKUP, LOG, LAUNCH_EPOCH, RESHAPED, FOLDGEN
    missing = [name for name in need if getattr(args, name) is None]
    if missing and not args.reshaped:
        sys.exit("%s required unless --reshaped" % ", ".join("--" + name for name in missing))
    ROOT, REVISION, BACKUP, LOG, RESHAPED, FOLDGEN = args.root, args.revision, args.backup, args.log, args.reshaped, args.foldgen
    LAUNCH_EPOCH = None
    if args.launch is not None:
        LAUNCH_EPOCH = datetime.datetime.fromisoformat(args.launch.replace("Z", "+00:00")).timestamp()


class FoldGen:
    def __init__(self):
        self.p = subprocess.Popen(
            [FOLDGEN], stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True
        )

    def _ask(self, line):
        self.p.stdin.write(line + "\n")
        self.p.stdin.flush()
        return [int(v) for v in self.p.stdout.readline().strip().split(",")]

    def stratified(self, y, k=K, seed=None):
        seed = SEED if seed is None else seed
        return self._ask("S %d %d %s" % (k, seed, "".join("1" if v else "0" for v in y)))

    def grouped(self, y, groups, k=K, seed=None):
        seed = SEED if seed is None else seed
        bits = "".join("1" if v else "0" for v in y)
        return self._ask("G %d %d %s %s" % (k, seed, bits, ",".join(str(g) for g in groups)))


def ratio(a, b):
    return 0.0 if b == 0 else a / b


def auroc(truth, scores):
    n_pos = sum(truth)
    n_neg = len(truth) - n_pos
    if n_pos == 0 or n_neg == 0:
        return None
    idx = sorted(range(len(scores)), key=lambda i: scores[i])
    rank_sum, i = 0.0, 0
    while i < len(idx):
        j = i
        while j + 1 < len(idx) and scores[idx[j + 1]] == scores[idx[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1.0
        rank_sum += sum(1 for k in idx[i : j + 1] if truth[k]) * avg
        i = j + 1
    return (rank_sum - (n_pos * (n_pos + 1)) / 2.0) / (float(n_pos) * float(n_neg))


def auprc(truth, scores):
    n_pos = sum(truth)
    if n_pos == 0:
        return None
    idx = sorted(range(len(scores)), key=lambda i: scores[i], reverse=True)
    tp = seen = 0
    ap, i = 0.0, 0
    while i < len(idx):
        j = i
        while j + 1 < len(idx) and scores[idx[j + 1]] == scores[idx[i]]:
            j += 1
        block_tp = sum(1 for k in idx[i : j + 1] if truth[k])
        seen += j - i + 1
        if block_tp > 0:
            tp += block_tp
            ap += block_tp / n_pos * (tp / seen)
        i = j + 1
    return ap


def roc_curve(truth, scores, max_points=60):
    if max_points < 2 or not any(truth) or all(truth):
        return []
    n_pos = float(max(sum(truth), 1))
    n_neg = float(max(len(truth) - sum(truth), 1))
    idx = sorted(range(len(scores)), key=lambda i: scores[i], reverse=True)
    pts = [(0.0, 0.0)]
    tp = fp = 0.0
    i = 0
    while i < len(idx):
        j = i
        while j + 1 < len(idx) and scores[idx[j + 1]] == scores[idx[i]]:
            j += 1
        for k in idx[i : j + 1]:
            if truth[k]:
                tp += 1.0
            else:
                fp += 1.0
        pts.append((fp / n_neg, tp / n_pos))
        i = j + 1
    if len(pts) > max_points:
        pts = [pts[i * (len(pts) - 1) // (max_points - 1)] for i in range(max_points)]
    return pts


def from_counts(tp, fp, tn, fn):
    n = tp + fp + tn + fn
    sens, spec, prec = ratio(tp, tp + fn), ratio(tn, tn + fp), ratio(tp, tp + fp)
    denom = math.sqrt(float(tp + fp) * float(tp + fn) * float(tn + fp) * float(tn + fn))
    return {
        "n": n, "tp": tp, "fp": fp, "tn": tn, "fn": fn,
        "accuracy": ratio(tp + tn, n),
        "balanced_accuracy": (sens + spec) / 2.0,
        "sensitivity": sens, "specificity": spec, "precision": prec,
        "f1": 0.0 if prec + sens == 0.0 else 2.0 * prec * sens / (prec + sens),
        "mcc": 0.0 if denom == 0.0 else (float(tp) * float(tn) - float(fp) * float(fn)) / denom,
        "vme": ratio(fn, tp + fn), "me": ratio(fp, fp + tn),
    }


def metrics(truth, scores, ranked):
    """Metrics of the called isolates; None when none is called. A null score is no call."""
    called = [(t, s) for t, s in zip(truth, scores) if s is not None]
    if not called:
        return None
    truth, scores = [t for t, _ in called], [s for _, s in called]
    tp = fp = tn = fn = 0
    for t, s in called:
        p = s >= 0.5
        if t and p:
            tp += 1
        elif p:
            fp += 1
        elif not t:
            tn += 1
        else:
            fn += 1
    m = from_counts(tp, fp, tn, fn)
    m["auroc"] = auroc(truth, scores) if ranked else None
    m["auprc"] = auprc(truth, scores) if ranked else None
    return m


def differs(a, b):
    """Whether two stored or computed values differ: null matches only null."""
    if a is None or b is None:
        return a is not b
    return abs(a - b) > TOL


def compare(stored, computed, where, out):
    """Append (where, field, stored, computed) for each disagreement; None means no metrics."""
    if stored is None or computed is None:
        if stored is not computed:
            out.append({"where": where, "field": "metrics", "stored": stored, "computed": computed})
        return stored is computed
    exact = True
    for key, want in computed.items():
        got = stored.get(key) if isinstance(stored, dict) else None
        if want is None or got is None or isinstance(want, int):
            ok = got == want
        else:
            ok = abs(got - want) <= TOL
            exact = exact and got == want
        if not ok:
            out.append({"where": where, "field": key, "stored": got, "computed": want})
    return exact


def rust_round(x):
    """f64::round: half away from zero, computed exactly."""
    fx = Fraction(x)
    f = math.floor(fx)
    d = fx - f
    if d > Fraction(1, 2):
        return f + 1
    if d < Fraction(1, 2):
        return f
    return f + 1 if x > 0 else f


def call_from_mic(mic, t):
    c, v = mic["censor"], mic["log2"]
    if c == "Exact":
        return "R" if v > t else "S"
    if c == "AtMost":
        return "S" if v <= t else None
    if c == "Above":
        return "R" if v >= t else None
    raise ValueError(c)


def infer_breakpoint(pairs):
    n = len(pairs)
    n_r = sum(1 for _, c in pairs if c == "R")
    if n < 50 or n_r < 10 or n - n_r < 10:
        return None
    lo = min(m["log2"] for m, _ in pairs)
    hi = max(m["log2"] for m, _ in pairs)
    best = None
    for t in range(lo - 1, hi + 1):
        correct = sum(1 for m, c in pairs if call_from_mic(m, t) == c)
        a = correct / n
        if best is None or a > best["agreement"]:
            best = {"r_above": t, "agreement": a, "n": n}
    return best if best is not None and best["agreement"] >= 0.95 else None


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fold_settings(rep):
    """Outer fold count and seed: eval_config's when the report records it."""
    cfg = rep.get("eval_config") or {}
    return cfg.get("folds", K), cfg.get("seed", SEED)


def fold_partition(rep, fg, truth):
    k, seed = fold_settings(rep)
    lineage = [s["lineage"] for s in rep["samples"]]
    return {"random": fg.stratified(truth, k, seed), "lineage": fg.grouped(truth, lineage, k, seed)}


def is_null_fold(schema, held):
    """A null-model fold: all held-out scores null (schema 2) or all exactly 0.5 (schema 1)."""
    if not held:
        return False
    if schema >= 2:
        return all(v is None for v in held)
    return all(v == 0.5 for v in held)


def null_folds(rep, truth, folds):
    """(scheme, learner, fold) of every null-model fold."""
    k = fold_settings(rep)[0]
    schema = rep.get("schema", 1)
    found = []
    for split in rep["splits"]:
        f = folds[split["scheme"]]
        for lr in split["learners"]:
            if lr["learner"] not in LEARNERS:
                continue
            for fold in range(k):
                if is_null_fold(schema, [lr["scores"][i] for i in range(len(truth)) if f[i] == fold]):
                    found.append((split["scheme"], lr["learner"], fold))
    return found


def learner_row(rep, scheme, name):
    split = next((s for s in rep["splits"] if s["scheme"] == scheme), None)
    return None if split is None else next((lr for lr in split["learners"] if lr["learner"] == name), None)


def calls_every_isolate(lr):
    nc = lr.get("no_call")
    return lr["metrics"] is not None and (nc is None or nc["resistant"] + nc["susceptible"] == 0)


def full_coverage_best(rep):
    """(learner, random BA) by the schema-2 rule: the highest random-split
    balanced accuracy among learners whose rows call every isolate on every
    split that ran, ties to the later learner; None when no learner does."""
    lineage_ran = any(split["scheme"] == "lineage" and split.get("unavailable_reason") is None for split in rep["splits"])
    best = None
    for name in LEARNERS:
        random, lineage = learner_row(rep, "random", name), learner_row(rep, "lineage", name)
        if random is None or not calls_every_isolate(random) or (lineage_ran and (lineage is None or not calls_every_isolate(lineage))):
            continue
        if best is None or random["metrics"]["balanced_accuracy"] >= best[1]:
            best = (name, random["metrics"]["balanced_accuracy"])
    return best


def max_random_best(rep):
    """(learner, random BA) by the schema-1 manifest rule: the highest
    random-split balanced accuracy among the five learners. Like revision-6
    grm report, it reads null metrics as balanced accuracy 0."""
    best = None
    for name in LEARNERS:
        random = learner_row(rep, "random", name)
        if random is None:
            continue
        ba = (random["metrics"] or {}).get("balanced_accuracy", 0.0)
        if best is None or ba >= best[1]:
            best = (name, ba)
    return best


def any_row_best(rep):
    """(learner, random BA) by the schema-1 log rule: Iterator::max_by over
    every random-split row, the catalog included; ties keep the last maximum."""
    best = None
    for lr in next(s for s in rep["splits"] if s["scheme"] == "random")["learners"]:
        if best is None or lr["metrics"]["balanced_accuracy"] >= best[1]:
            best = (lr["learner"], lr["metrics"]["balanced_accuracy"])
    return best


def lite(lr, schema):
    """The manifest's rates for one learner row, as grm report writes them for this schema."""
    if lr is None:
        return None
    m = lr["metrics"]
    keys = ("balanced_accuracy", "sensitivity", "specificity", "vme", "me")
    if m is None:
        out = dict.fromkeys(keys + ("auroc",))
    else:
        out = {k: m.get(k, 0.0) for k in keys}
        out["auroc"] = m.get("auroc")
    if schema >= 2 and lr.get("no_call") is not None:
        out["no_call"] = lr["no_call"]
    return out


def lite_differs(stored, want):
    if stored is None or want is None:
        return stored is not want
    if set(stored) != set(want):
        return True
    return any(stored[k] != want[k] if k == "no_call" else differs(stored[k], want[k]) for k in want)


def decode_bits(text, n):
    """Presence over samples from a feature card: base64 of little-endian u64 words."""
    raw = base64.b64decode(text)
    return [bool(raw[i >> 3] >> (i & 7) & 1) for i in range(n)]


species_finished = set()


def done_datasets():
    """Done lines and recompute reasons in LOG, by dataset; fills species_finished."""
    done, recomputed = {}, {}
    if LOG is None or not os.path.exists(LOG):
        return done, recomputed
    for line in open(LOG, errors="replace"):
        m = re.match(r"([a-z_]+): \d+ dataset reports written", line)
        if m:
            species_finished.add(m.group(1))
        m = DONE.match(line)
        if m:
            words = m.group(3)
            entry = {"form": "unrecognized", "learner": None, "ba": None, "words": words}
            if words == NO_WINNER:
                entry["form"] = "full coverage"
            elif w := FULL_COVERAGE.fullmatch(words):
                entry.update(form="full coverage", learner=w.group(1), ba=w.group(2))
            elif w := ANY_ROW.fullmatch(words):
                entry.update(form="any row", learner=w.group(1), ba=w.group(2))
            done[m.group(1)] = entry
        m = RECOMPUTE.match(line)
        if m:
            recomputed[m.group(1)] = m.group(2)
    return done, recomputed


def load_manifest(root):
    """The published manifest's schema, rows by dataset id, and panel rows by dataset id."""
    manifest = json.load(open(root + "/data/report/index.json"))
    panels = defaultdict(list)
    for p in manifest.get("panels", []):
        panels[p["dataset"]].append(p)
    return {"schema": manifest.get("schema", 1), "evaluation_revision": manifest["evaluation_revision"],
            "generated_at": manifest.get("generated_at"), "rows": {r["id"]: r for r in manifest["datasets"]},
            "panels": panels}


def audit_file(path, fg, counted_cache, log_done, log_recomputed, manifest):
    issues = []
    rec = {"file": path}
    st = os.stat(path)
    rec["mtime_epoch"] = st.st_mtime
    rec["sha256"] = sha256(path)
    rep = json.load(open(path))
    rid = rep["id"]
    slug = rep["species_slug"]
    schema = rep.get("schema", 1)
    rec.update(id=rid, species=slug, revision=rep["evaluation"]["revision"], schema=schema,
               valid=rep["evaluation"]["valid"], reasons=rep["evaluation"]["reasons"])
    samples = rep["samples"]
    n = len(samples)
    truth = [bool(s["resistant"]) for s in samples]
    n_r = sum(truth)
    c = rep["counts"]
    rec["counts"] = {"samples": c["samples"], "resistant": c["resistant"],
                     "susceptible": c["susceptible"], "labels_from_mic": c["labels_from_mic"]}

    def check(cond, what, where="counts/structure", **detail):
        if not cond:
            issues.append(dict(where=where, what=what, **detail))

    check(schema == (2 if RESHAPED or rec["revision"] >= 7 else 1), "schema does not match the revision", schema=schema)
    check(c["samples"] == n, "counts.samples != len(samples)", stored=c["samples"], computed=n)
    check(c["resistant"] == n_r, "counts.resistant", stored=c["resistant"], computed=n_r)
    check(c["susceptible"] == n - n_r, "counts.susceptible", stored=c["susceptible"], computed=n - n_r)
    check(c["resistant"] + c["susceptible"] == c["samples"], "R+S != samples")
    from_mic = sum(1 for s in samples if s["from_mic"])
    check(c["labels_from_mic"] == 0 and from_mic == 0, "inferred labels present",
          stored=c["labels_from_mic"], computed=from_mic)
    check(len({s["id"] for s in samples}) == n, "duplicate sample ids")
    lineages = len({s["lineage"] for s in samples})
    check(c["lineages"] == lineages, "counts.lineages", stored=c["lineages"], computed=lineages)
    # Reports that record the class minimum were gated on it; earlier ones on 25.
    min_class = c.get("min_class", 25)
    check(n_r >= min_class and n - n_r >= min_class, "class below the minimum", min_class=min_class)
    check(c["min_support"] == max(n // 200, 3), "min_support", stored=c["min_support"])

    cfg = rep.get("eval_config")
    rec["eval_config"] = cfg
    if cfg is None:
        check(rec["revision"] < 7, "report from revision 7 on without eval_config")
    else:
        check(set(cfg) == EVAL_CONFIG_KEYS, "eval_config fields", fields=sorted(cfg))
        check(set(cfg.get("forest", {})) == FOREST_KEYS, "eval_config.forest fields (a per-fit seed is not recorded)",
              fields=sorted(cfg.get("forest", {})))
        check(set(cfg.get("gbt", {})) == GBT_KEYS, "eval_config.gbt fields (a per-fit seed is not recorded)",
              fields=sorted(cfg.get("gbt", {})))
    k, seed = fold_settings(rep)

    folds = fold_partition(rep, fg, truth)
    # A lineage never spans folds, by construction; check both classes per
    # fold wherever the split was evaluated.
    available = {sp["scheme"] for sp in rep["splits"] if not sp.get("unavailable_reason")}
    check(len(rep["splits"]) == 2 and "random" in available, "split schemes", schemes=[sp["scheme"] for sp in rep["splits"]])
    for scheme in available:
        f = folds[scheme]
        for fold in range(k):
            r = sum(1 for i in range(n) if f[i] == fold and truth[i])
            s = sum(1 for i in range(n) if f[i] == fold and not truth[i])
            check(r > 0 and s > 0, "fold lacks a class", scheme=scheme, fold=fold, r=r, s=s)

    learner_rows = []
    exact_all = True
    tie_counts = Counter()
    nullf = []
    half_called = []
    constant_other = []
    small_folds = []
    for split in rep["splits"]:
        scheme = split["scheme"]
        f = folds[scheme]
        check(split["folds"] == k, "split.folds", scheme=scheme)
        check(split["n_groups"] == lineages, "split.n_groups", scheme=scheme, stored=split["n_groups"])
        if split.get("unavailable_reason"):
            check(split["learners"] == [], "unavailable split has learners", scheme=scheme)
            check("test_fold" not in split, "unavailable split saves test_fold", scheme=scheme)
            rec.setdefault("unavailable", []).append({"scheme": scheme, "reason": split["unavailable_reason"]})
            continue
        saved = split.get("test_fold")
        if schema >= 2:
            if saved != f:
                aligned = saved is not None and len(saved) == n
                issues.append({"where": scheme + "/test_fold", "what": "test_fold differs from the rebuilt folds",
                               "saved_len": None if saved is None else len(saved),
                               "positions": sum(a != b for a, b in zip(saved, f)) if aligned else None})
        else:
            check(saved is None, "schema-1 split saves test_fold", scheme=scheme)
        names = [lr["learner"] for lr in split["learners"]]
        want = LEARNERS + (["catalog"] if rep.get("catalog") else [])
        check(names == want, "learner list", scheme=scheme, stored=names, computed=want)
        for fold in range(k):
            r = sum(1 for i in range(n) if f[i] == fold and truth[i])
            s = sum(1 for i in range(n) if f[i] == fold and not truth[i])
            if min(r, s) < 5:
                small_folds.append({"scheme": scheme, "fold": fold, "r": r, "s": s})
        for lr in split["learners"]:
            name = lr["learner"]
            sc = lr["scores"]
            where = "%s/%s" % (scheme, name)
            check(len(sc) == n, "score length", where=where)
            check(all(v is None or (isinstance(v, (int, float)) and math.isfinite(v)) for v in sc), "non-finite score", where=where)
            nc = lr.get("no_call")
            uncalled = {"resistant": sum(1 for t, v in zip(truth, sc) if v is None and t),
                        "susceptible": sum(1 for t, v in zip(truth, sc) if v is None and not t)}
            if name == "catalog" or schema < 2:
                check(None not in sc, "null score in a row that calls every isolate", where=where)
                check(nc is None, "no_call on a row that calls every isolate", where=where)
            elif nc != uncalled:
                issues.append({"where": where + "/no_call", "what": "no_call differs from the null scores", "stored": nc, "computed": uncalled})
            ranked = name in RANKED
            called = [i for i in range(n) if sc[i] is not None]
            tc, scc = [truth[i] for i in called], [sc[i] for i in called]
            m = metrics(tc, scc, ranked)
            exact_all &= compare(lr["metrics"], m, where + "/pooled", issues)
            if ranked:
                roc = roc_curve(tc, scc)
                got = [tuple(p) for p in lr["roc"]]
                if len(got) != len(roc) or any(abs(a[0] - b[0]) > TOL or abs(a[1] - b[1]) > TOL for a, b in zip(got, roc)):
                    issues.append({"where": where + "/roc", "what": "roc points differ", "stored_len": len(got), "computed_len": len(roc)})
            else:
                check(lr["roc"] == [], "unranked learner has roc", where=where)
            check(lr["conformal"] == [], "conformal not empty", where=where)
            rates = ("tp", "fp", "tn", "fn", "balanced_accuracy")
            if name == "catalog":
                check(lr["fold_balanced_accuracy"] == [] and lr["fold_complexity"] == [], "catalog fold fields", where=where)
                flags = [1.0 if s["catalog_resistant"] else 0.0 for s in samples]
                check(sc == flags, "catalog scores != sample catalog_resistant", where=where)
                exact_all &= compare(rep["catalog"]["metrics"], m, "catalog/metrics", issues)
                learner_rows.append({"scheme": scheme, "learner": name, "called": len(called), **{k2: m and m[k2] for k2 in rates}})
                continue
            if name == "scm":
                check(all(v in (0.0, 1.0) for v in scc), "scm score not 0/1", where=where)
            complexity = lr["fold_complexity"] if len(lr["fold_complexity"]) == k else None
            fba = []
            for fold in range(k):
                idx = [i for i in range(n) if f[i] == fold]
                tk = [truth[i] for i in idx]
                sk = [sc[i] for i in idx]
                if is_null_fold(schema, sk):
                    fm = metrics(tk, sk, False)
                    fba.append(None if fm is None else fm["balanced_accuracy"])
                    nullf.append({"scheme": scheme, "learner": name, "fold": fold, "n": len(sk),
                                  "r": sum(tk), "s": len(tk) - sum(tk), "fold_ba": fba[-1],
                                  "complexity": None if complexity is None else complexity[fold]})
                    continue
                if None in sk:
                    issues.append({"where": "%s/fold %d" % (where, fold), "what": "fold called in part", "uncalled": sk.count(None)})
                fm = metrics(tk, sk, False)
                fba.append(None if fm is None else fm["balanced_accuracy"])
                held = [v for v in sk if v is not None]
                if schema >= 2 and held and all(v == 0.5 for v in held):
                    # A fitted model's 0.5 is a resistant call from revision 7.
                    half_called.append({"scheme": scheme, "learner": name, "fold": fold, "n": len(held),
                                        "complexity": None if complexity is None else complexity[fold]})
                elif held and len(set(held)) == 1:
                    constant_other.append({"scheme": scheme, "learner": name, "fold": fold, "value": held[0], "n": len(held)})
                else:
                    tie_counts[(scheme, name)] += sum(1 for v in held if v == 0.5)
            stored_fba = lr["fold_balanced_accuracy"]
            if len(stored_fba) != k or any(differs(a, b) for a, b in zip(stored_fba, fba)):
                issues.append({"where": where + "/fold_ba", "stored": stored_fba, "computed": fba})
            elif stored_fba != fba:
                exact_all = False
            check(len(lr["fold_complexity"]) == k, "fold_complexity length", where=where)
            if schema >= 2 and complexity is not None:
                for fold in range(k):
                    null = is_null_fold(schema, [sc[i] for i in range(n) if f[i] == fold])
                    zero = name in ("gbt", "lasso") and complexity[fold] == 0
                    if null != zero:
                        issues.append({"where": "%s/fold %d" % (where, fold), "complexity": complexity[fold],
                                       "what": "null fold without a zero gbt or lasso complexity" if null else "null fit scored"})
            learner_rows.append({"scheme": scheme, "learner": name, "called": len(called), "no_call": nc,
                                 **{k2: m and m[k2] for k2 in rates + ("vme", "me", "auroc")}})
    rec["learners"] = learner_rows
    rec["null_folds"] = nullf
    rec["half_called_folds"] = half_called
    rec["constant_non_half_folds"] = constant_other
    rec["half_ties_in_non_null_folds"] = {"%s/%s" % key: v for key, v in tie_counts.items() if v}
    rec["folds_with_class_below_5"] = small_folds
    rec["metrics_bit_exact"] = exact_all
    best = full_coverage_best(rep) if schema >= 2 else max_random_best(rep)
    rec["best"] = None if best is None else {"learner": best[0], "random": best[1]}

    # The final forest's importance names dataset rows, and the 30 the
    # report annotates have feature cards that match their presence bits.
    imp = rep["models"]["rf_importance"]
    ids = [w["feature"] for w in imp]
    check(len(imp) <= 50, "rf_importance longer than 50", n=len(imp))
    check(len(set(ids)) == len(ids), "rf_importance repeats a feature")
    check(all(isinstance(fid, int) and 0 <= fid < c["features"] for fid in ids), "rf_importance feature is not a dataset row",
          outside=[fid for fid in ids if not (isinstance(fid, int) and 0 <= fid < c["features"])][:5])
    positive = all(isinstance(w["weight"], (int, float)) and math.isfinite(w["weight"]) and w["weight"] > 0 for w in imp)
    check(positive, "rf_importance weight not positive")
    if positive and all(isinstance(fid, int) for fid in ids):
        check(all((a["weight"], -a["feature"]) > (b["weight"], -b["feature"]) for a, b in zip(imp, imp[1:])),
              "rf_importance not in weight then feature order")
    cards = {card["feature"]: card for card in rep["features"]}
    missing = [fid for fid in ids[:30] if fid not in cards]
    check(not missing, "rf_importance feature without a feature card", features=missing[:5])
    wrong = []
    for fid in ids[:30]:
        if fid in cards:
            bits = decode_bits(cards[fid]["bits"], n)
            counted = [sum(1 for b, t in zip(bits, truth) if b and t), sum(1 for b, t in zip(bits, truth) if b and not t)]
            if counted != [cards[fid]["n_resistant"], cards[fid]["n_susceptible"]]:
                wrong.append(fid)
    check(not wrong, "rf feature card counts differ from its bits", features=wrong[:5])
    rec["rf_cards_checked"] = sum(1 for fid in ids[:30] if fid in cards)

    # Catalog baseline and catalog panels from counted.json hits.
    hits = counted_cache.get(slug)
    if hits is None:
        cpath = "%s/data/species/%s/counted.json" % (ROOT, slug)
        hits = {g["genome_id"]: g["genes"] for g in json.load(open(cpath))["genomes"]}
        counted_cache[slug] = hits
    groups_names = rep["catalog_groups"]
    mapping_ok = all(
        s["id"] in hits and [groups_names[g] for g in s["groups"]] == [h["group"] for h in hits[s["id"]]]
        for s in samples
    )
    rec["counted_hits_match_samples"] = mapping_ok
    rec["samples_missing_from_counted"] = sum(1 for s in samples if s["id"] not in hits)
    if rep.get("catalog"):
        table = defaultdict(lambda: [0, 0])
        for s, t in zip(samples, truth):
            for g in s["groups"]:
                table[groups_names[g]][0 if t else 1] += 1
        stored = {row["group"]: [row["n_resistant"], row["n_susceptible"]] for row in rep["catalog"]["groups"]}
        if stored != dict(table):
            issues.append({"where": "catalog/groups", "what": "group carrier counts differ"})
    else:
        check(all(s["catalog_resistant"] is None for s in samples), "catalog flags without catalog")

    panel_rows = []
    probe_rows = []
    panel_ids = set()
    for p in rep["panels"]:
        panel_ids.add(p["id"])
        m = p["metrics"]
        nc = p.get("no_call")
        is_probe = p["id"] in PROBE_IDS
        if is_probe and nc is None:
            # Revision-5 probe rows predate no-call accounting.
            check(rec["revision"] < 6, "revision-%d probe panel without no_call" % rec["revision"], panel=p["id"])
            nc = {"resistant": 0, "susceptible": 0, "field_absent": True}
        if not is_probe:
            check(nc is None, "catalog panel carries no_call", panel=p["id"])
            row = {"id": p["id"], "kind": "catalog_panel", "tp": m["tp"], "fp": m["fp"], "tn": m["tn"], "fn": m["fn"]}
            check(m["n"] == n, "catalog panel n != samples", panel=p["id"])
            if mapping_ok:
                targets = p["targets"]
                det = defaultdict(lambda: [0, 0])
                scores = []
                for s, t in zip(samples, truth):
                    found = sorted({h["group"] for h in hits[s["id"]] if any(h["group"] == x or h["family"] == x for x in targets)})
                    scores.append(1.0 if found else 0.0)
                    for g in found:
                        det[g][0 if t else 1] += 1
                compare(m, metrics(truth, scores, False), "panel/" + p["id"], issues)
                want_det = [[g, v[0], v[1]] for g, v in sorted(det.items(), key=lambda kv: kv[0].encode())]
                if p["detected"] != want_det:
                    issues.append({"where": "panel/" + p["id"] + "/detected", "stored": p["detected"], "computed": want_det})
                row["recomputed_from_hits"] = True
            panel_rows.append(row)
        else:
            tp, fp, tn, fn = m["tp"], m["fp"], m["tn"], m["fn"]
            row = {"id": p["id"], "kind": "probe", "tp": tp, "fp": fp, "tn": tn, "fn": fn,
                   "no_call_r": nc["resistant"], "no_call_s": nc["susceptible"],
                   "no_call_field": not nc.get("field_absent", False)}
            check(tp + fn + nc["resistant"] == n_r, "probe R partition", panel=p["id"], row=row, cohort_r=n_r)
            check(fp + tn + nc["susceptible"] == n - n_r, "probe S partition", panel=p["id"], row=row, cohort_s=n - n_r)
            check(m["n"] == tp + fp + tn + fn, "probe n != called", panel=p["id"], n=m["n"])
            check(min(tp, fp, tn, fn, nc["resistant"], nc["susceptible"]) >= 0, "negative probe count", panel=p["id"])
            check(p["detected"] == [["mutation detected", tp, fp]], "probe detected row", panel=p["id"], detected=p["detected"])
            want = from_counts(tp, fp, tn, fn)
            want["auroc"] = None
            want["auprc"] = None
            compare(m, want, "probe/" + p["id"], issues)
            probe_rows.append(row)
    rec["panels"] = panel_rows
    rec["probe_panels"] = probe_rows
    # Revision-5 reports predate the field.
    unsupported = rep.get("unsupported_assays", [])
    rec["has_unsupported_field"] = "unsupported_assays" in rep
    rec["unsupported_assays"] = [{"id": u["id"], "reason": u["reason"]} for u in unsupported]
    for u in unsupported:
        check(bool(u["reason"]), "unsupported assay without reason", assay=u["id"])
        check(u["id"] not in panel_ids, "unsupported assay also scored", assay=u["id"])

    # MIC regression and categorical agreement.
    mic = rep.get("mic")
    has = [i for i, s in enumerate(samples) if s["mic"] is not None]
    n_exact = sum(1 for i in has if samples[i]["mic"]["censor"] == "Exact")
    rec["mic_samples"] = {"with_mic": len(has), "exact": n_exact}
    if mic is None:
        rec["mic"] = None
        # evaluate_mic needs >= 100 MICs and >= 50 exact; later early returns
        # (an empty fold or holdout) are possible but should be rare.
        if len(has) >= 100 and n_exact >= 50:
            rec["mic_absent_with_enough_mics"] = True
    else:
        check(mic["n"] == len(has), "mic.n", stored=mic["n"], computed=len(has))
        check(mic["n_exact"] == n_exact, "mic.n_exact", stored=mic["n_exact"], computed=n_exact)
        pred = mic["predicted"]
        check(len(pred) == n, "mic.predicted length")
        check(all((pred[i] is None) == (samples[i]["mic"] is None) for i in range(n)), "mic.predicted null pattern")
        local = fg.stratified([truth[i] for i in has], k, seed)
        mf = {i: local[j] for j, i in enumerate(has)}
        bps = []
        for fold in range(k):
            pairs = [(samples[i]["mic"], "R" if truth[i] else "S") for i in has if mf[i] != fold]
            bps.append(infer_breakpoint(pairs))
        if bps != mic["fold_breakpoints"]:
            issues.append({"where": "mic/fold_breakpoints", "stored": mic["fold_breakpoints"], "computed": bps})
        ea = ea_all = 0
        conf = Counter()
        for i in has:
            p_ = rust_round(pred[i])
            mm = samples[i]["mic"]
            if mm["censor"] == "Exact":
                conf[(mm["log2"], p_)] += 1
                ok = abs(p_ - mm["log2"]) <= 1
                ea += ok
            elif mm["censor"] == "AtMost":
                ok = p_ <= mm["log2"] + 1
            else:
                ok = p_ >= mm["log2"]
            ea_all += ok
        want_conf = [[o, p_, v] for (o, p_), v in sorted(conf.items())]
        if want_conf != mic["confusion"]:
            issues.append({"where": "mic/confusion", "what": "dilution confusion differs"})
        e1, e2 = ea / max(n_exact, 1), ea_all / len(has)
        if abs(mic["essential_agreement"] - e1) > TOL or abs(mic["essential_agreement_all"] - e2) > TOL:
            issues.append({"where": "mic/ea", "stored": [mic["essential_agreement"], mic["essential_agreement_all"]], "computed": [e1, e2]})
        cat = None
        if all(b is not None for b in mic["fold_breakpoints"]):
            calls = [1.0 if rust_round(pred[i]) > mic["fold_breakpoints"][mf[i]]["r_above"] else 0.0 for i in has]
            cat = metrics([truth[i] for i in has], calls, False)
        if (cat is None) != (mic["categorical"] is None):
            issues.append({"where": "mic/categorical", "what": "presence differs", "stored": mic["categorical"] is not None})
        elif cat is not None:
            compare(mic["categorical"], cat, "mic/categorical", issues)
        rec["mic"] = {"n": len(has), "n_exact": n_exact, "ea": e1, "ea_all": e2,
                      "categorical": None if cat is None else {k2: cat[k2] for k2 in ("tp", "fp", "tn", "fn")}}

    # Provenance: launch time, backup copy, pipeline log, manifest row.
    # grm writes results/<id>.json; the id also finds the backup of a copy
    # audited from elsewhere, as controls.py does.
    bpath = None if BACKUP is None else "%s/data/species/%s/results/%s.json" % (BACKUP, slug, rid)
    rec["backup_exists"] = bpath is not None and os.path.exists(bpath)
    rec["backup_sha256"] = sha256(bpath) if rec["backup_exists"] else None
    rec["identical_to_backup"] = rec["backup_sha256"] == rec["sha256"]
    rec["written_after_launch"] = None if LAUNCH_EPOCH is None else st.st_mtime >= LAUNCH_EPOCH
    rec["log_done"] = log_done.get(rid)
    rec["log_recomputed"] = log_recomputed.get(rid)
    if not RESHAPED and rec["revision"] == REVISION:
        check(rec["written_after_launch"], "revision-%d file older than launch" % REVISION)
        check(not rec["identical_to_backup"], "revision-%d file identical to its backup" % REVISION)
        check(rec["log_done"] is not None, "revision-%d file without a done log line" % REVISION)
        # grm logs a recompute only over a saved report; a dataset new to
        # this run has none.
        check((rec["log_recomputed"] is not None) == rec["backup_exists"],
              "revision-%d file without a recompute log line" % REVISION if rec["backup_exists"]
              else "recompute logged for a dataset without a backup report")
        logged = rec["log_done"]
        if logged is not None:
            if logged["form"] == "full coverage":
                check(schema >= 2, "full-coverage done line for a schema-1 file")
                want = (None, None) if best is None else (best[0], "%.3f" % best[1])
            else:
                check(logged["form"] == "any row" and schema < 2, "done line form does not match the schema", form=logged["form"])
                # Iterator::max_by over every random-split row, catalog included.
                b = any_row_best(rep) if schema < 2 else None
                want = (None, None) if b is None else (b[0], "%.3f" % b[1])
            check((logged["learner"], logged["ba"]) == want, "log best learner differs", stored=logged, computed=list(want))
    elif not RESHAPED:
        check(not rec["written_after_launch"] or rec["identical_to_backup"], "non-revision-%d file written after launch" % REVISION)
        check(rec["log_done"] is None, "log reports a recompute for a non-revision-%d file" % REVISION)
        check(slug not in species_finished, "species finished in this run but file is not revision %d" % REVISION)
    mrow = manifest["rows"].get(rid)
    rec["manifest_revision"] = None if mrow is None else mrow.get("evaluation", {}).get("revision")
    ppath = "%s/data/report/datasets/%s.json" % (ROOT, rid)
    if os.path.exists(ppath):
        pub = json.load(open(ppath))
        rec["published_revision"] = pub["evaluation"]["revision"]
        rec["published_identical"] = sha256(ppath) == rec["sha256"]
    else:
        rec["published_revision"] = None
        rec["published_identical"] = None
    # The manifest row was built from this file only when its published copy
    # is this file; otherwise the row is older or newer and says nothing here.
    rec["manifest_compared"] = mrow is not None and bool(rec["published_identical"])
    if rec["manifest_compared"]:
        check(manifest["schema"] >= schema, "schema-%d report in a schema-%d manifest" % (schema, manifest["schema"]))
        mbest = full_coverage_best(rep) if manifest["schema"] >= 2 else max_random_best(rep)
        if mbest is None:
            want_row = {"best_learner": None, "best_random": None, "best_lineage": None}
        else:
            lineage = learner_row(rep, "lineage", mbest[0])
            want_row = {"best_learner": mbest[0], "best_random": mbest[1],
                        "best_lineage": None if lineage is None or lineage["metrics"] is None else lineage["metrics"]["balanced_accuracy"]}
        got_row = {key: mrow.get(key) for key in want_row}
        if got_row["best_learner"] != want_row["best_learner"] or differs(got_row["best_random"], want_row["best_random"]) \
                or differs(got_row["best_lineage"], want_row["best_lineage"]):
            issues.append({"where": "manifest/best", "what": "manifest best learner differs", "stored": got_row, "computed": want_row})
        for entry in mrow["learners"]:
            for scheme in ("random", "lineage"):
                want_lite = lite(learner_row(rep, scheme, entry["learner"]), manifest["schema"])
                if lite_differs(entry[scheme], want_lite):
                    issues.append({"where": "manifest/%s/%s" % (scheme, entry["learner"]), "what": "manifest rates differ",
                                   "stored": entry[scheme], "computed": want_lite})
        check([entry["learner"] for entry in mrow["learners"]] == LEARNERS, "manifest learner rows")
        prows = manifest["panels"].get(rid, [])
        want_panels = [] if mbest is None else [p["id"] for p in rep["panels"]]
        best_model = None if mbest is None else lite(learner_row(rep, "random", mbest[0]), manifest["schema"])
        if [p["panel"] for p in prows] != want_panels or any(p["best_learner"] != mbest[0] or lite_differs(p["best_model"], best_model) for p in prows):
            issues.append({"where": "manifest/panels", "what": "panel rows do not follow the best learner",
                           "stored": [[p["panel"], p["best_learner"]] for p in prows], "computed": want_panels})

    # Backup comparison for the same dataset.
    if rec["backup_exists"] and not rec["identical_to_backup"]:
        old = json.load(open(bpath))
        old_ids = [s["id"] for s in old["samples"]]
        old_truth = [bool(s["resistant"]) for s in old["samples"]]
        same_ids = old_ids == [s["id"] for s in samples]
        rec["vs_backup"] = {
            "backup_revision": old["evaluation"]["revision"],
            "same_ordered_ids": same_ids,
            "same_labels": same_ids and old_truth == truth,
            "n_before": len(old_ids), "n_after": n,
            "r_before": sum(old_truth), "r_after": n_r,
            "ids_added": len(set(s["id"] for s in samples) - set(old_ids)),
            "ids_removed": len(set(old_ids) - set(s["id"] for s in samples)),
        }
        if same_ids:
            diff = []
            for sp_new, sp_old in zip(rep["splits"], old["splits"]):
                for a, b in zip(sp_new["learners"], sp_old["learners"]):
                    if a["scores"] != b["scores"]:
                        diff.append("%s/%s" % (sp_new["scheme"], a["learner"]))
                    elif a["learner"] == "rf" and rec["revision"] >= 7 and old["evaluation"]["revision"] < 7:
                        # Revision 7 screens the forest's features and splits
                        # its random streams, so every forest row changes.
                        issues.append({"where": "%s/rf" % sp_new["scheme"], "what": "forest scores unchanged from the backup"})
            rec["vs_backup"]["learners_with_changed_scores"] = diff
            rec["vs_backup"]["old_null_folds"] = len(null_folds(old, old_truth, fold_partition(old, fg, old_truth)))
            rec["vs_backup"]["old_panels"] = [[p["id"], p["metrics"]["tp"], p["metrics"]["fp"], p["metrics"]["tn"], p["metrics"]["fn"], p.get("no_call")] for p in old["panels"]]
            rec["vs_backup"]["old_unsupported"] = [u["id"] for u in old.get("unsupported_assays", [])]
    rec["issues"] = issues
    return rec


def main():
    p = options("Read-only metric audit of GRM result reports.")
    p.add_argument("out")
    p.add_argument("pattern", nargs="?", default="*", help="species glob")
    args = p.parse_args()
    configure(args)
    started = time.time()
    fg = FoldGen()
    log_done, log_recomputed = done_datasets()
    manifest = load_manifest(ROOT)
    counted_cache = {}
    files, scope = select_results(ROOT + "/data", args.pattern, LOG)
    records = []
    for path in files:
        try:
            rec = audit_file(str(path), fg, counted_cache, log_done, log_recomputed, manifest)
            if rec["id"] != path.stem or rec["species"] != path.parent.parent.name:
                rec["issues"].append({"where": "curation", "what": "report identity differs from its curated path"})
            records.append(rec)
        except Exception as e:  # report and keep going
            records.append({"file": str(path), "error": repr(e)})
        counted_cache = {k: v for k, v in counted_cache.items() if k == records[-1].get("species")}
    by_rev = Counter(r.get("revision") for r in records)
    summary = {
        "started_epoch": started,
        "finished_epoch": time.time(),
        "revision": REVISION,
        "reshaped": RESHAPED,
        "root": ROOT,
        "backup": BACKUP,
        "log": LOG,
        "launch_epoch": LAUNCH_EPOCH,
        "files": len(records),
        "by_revision": {str(k): v for k, v in by_rev.items()},
        "by_schema": {str(k): v for k, v in Counter(r.get("schema") for r in records).items()},
        "errors": [r for r in records if "error" in r] + scope["missing_curated_results"],
        "files_with_issues": [r["id"] for r in records if r.get("issues")],
        "issue_count": sum(len(r.get("issues", [])) for r in records),
        "manifest_schema": manifest["schema"],
        "manifest_evaluation_revision": manifest["evaluation_revision"],
        "manifest_generated_at": manifest["generated_at"],
        "manifest_compared": sum(1 for r in records if r.get("manifest_compared")),
        "log_done": len(log_done),
        "log_done_forms": dict(Counter(d["form"] for d in log_done.values())),
        "log_recompute_reasons": dict(Counter(log_recomputed.values())),
        "species_finished": sorted(species_finished),
        **scope,
    }
    json.dump({"summary": summary, "records": records}, open(args.out, "w"), indent=1, sort_keys=True)
    print(json.dumps(summary, indent=1, sort_keys=True))
    if scope["missing_curated_results"]:
        sys.exit(1)


if __name__ == "__main__":
    main()
