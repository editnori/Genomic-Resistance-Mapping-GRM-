"""Synthetic GRM dataset reports for the checker and audit tests.

One dataset of 80 isolates (40 resistant) in 16 lineages of five, scored by
the five learners and the catalog. Every stored metric comes from direct
formulas here, not from the code under test: confusion counts by
comparison, AUROC by counting resistant-susceptible pairs, AUPRC and ROC by
sweeping each distinct score as a cutoff. Scores are multiples of 1/32.

In the default revision-7 report gbt has the best random-split balanced
accuracy but its lineage fold 2 is a null fit, and lasso's random fold 4 is
a null fit, so the best learner that calls every isolate is rf.
"""
import base64
import copy
import json
import math
import os
import subprocess
import sys

# The audit directory sits beside the checker it tests.
LANE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FOLDGEN = os.path.join(LANE, "foldgen", "target", "release", "foldgen")
CHECKER = os.path.join(os.path.dirname(LANE), "check-results.py")
SLUG = "fixture_species"
ID = SLUG + "__fixturemycin"
N = 80
TRUTH = [i % 2 == 0 for i in range(N)]
LINEAGE = [i // 5 for i in range(N)]
LEARNERS = ["scm", "cart", "rf", "gbt", "lasso"]
RANKED = {"cart", "rf", "gbt", "lasso"}
EVAL_CONFIG = {
    "folds": 5, "inner_folds": 3, "scm_max_rules": 6, "cart_max_depth": 8,
    "forest": {"n_trees": 200, "subsample": 0.632, "max_depth": 30, "min_leaf": 1},
    "gbt": {"rounds": 400, "learning_rate": 0.1, "max_depth": 4, "lambda": 1.0,
            "min_child_weight": 1.0, "colsample": 0.5, "patience": 25},
    "screened_features": 4000, "lasso_lambdas": 30, "seed": 42,
}
WRONG = {"random": {"scm": 12, "cart": 10, "rf": 4, "gbt": 2, "lasso": 6},
         "lineage": {"scm": 14, "cart": 12, "rf": 8, "gbt": 6, "lasso": 10}}
COMPLEXITY = {"scm": 3.0, "cart": 7.0, "rf": 200.0, "gbt": 40.0, "lasso": 12.0}
NULL_FITS = (("lineage", "gbt", 2), ("random", "lasso", 4))
FEATURES = 20
RF_IMPORTANCE = [{"feature": 7, "weight": 0.25}, {"feature": 2, "weight": 0.125}, {"feature": 11, "weight": 0.125}]
SCOPE = ("Internal cross-validation on public BV-BRC genomes scored against reported R/S labels. "
         "This is not external or clinical validation.")


def folds():
    """Outer folds from foldgen, the pipeline's own fold code, at k 5 and seed 42."""
    bits = "".join("1" if t else "0" for t in TRUTH)
    text = "S 5 42 %s\nG 5 42 %s %s\n" % (bits, bits, ",".join(map(str, LINEAGE)))
    lines = subprocess.run([FOLDGEN], input=text, capture_output=True, text=True, check=True).stdout.splitlines()
    return {"random": [int(v) for v in lines[0].split(",")], "lineage": [int(v) for v in lines[1].split(",")]}


FOLDS = folds()


def wrong(scheme, learner, i):
    """Whether this learner calls isolate i wrongly; 37 is coprime to 80, so exactly WRONG of them are."""
    return (i * 37 + LEARNERS.index(learner) * 11) % N < WRONG[scheme][learner]


def score(scheme, learner, i, variant=0):
    resistant_call = TRUTH[i] != wrong(scheme, learner, i)
    if learner == "scm":
        return 1.0 if resistant_call else 0.0
    step = (i * 7 + LEARNERS.index(learner) * 3 + variant * 5) % 15
    return (17 + step) / 32 if resistant_call else (1 + step) / 32


def ratio(a, b):
    return a / b if b else 0.0


def rates(truth, scores, ranked):
    """Metrics of the called isolates by direct formulas; None when none is called."""
    pairs = [(t, s) for t, s in zip(truth, scores) if s is not None]
    if not pairs:
        return None
    tp = sum(1 for t, s in pairs if t and s >= 0.5)
    fp = sum(1 for t, s in pairs if not t and s >= 0.5)
    tn = sum(1 for t, s in pairs if not t and s < 0.5)
    fn = sum(1 for t, s in pairs if t and s < 0.5)
    sens, spec, prec = ratio(tp, tp + fn), ratio(tn, tn + fp), ratio(tp, tp + fp)
    denom = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    pos = [s for t, s in pairs if t]
    neg = [s for t, s in pairs if not t]
    auroc = auprc = None
    if ranked and pos and neg:
        auroc = sum(1.0 if p > q else 0.5 if p == q else 0.0 for p in pos for q in neg) / (len(pos) * len(neg))
    if ranked and pos:
        auprc, before = 0.0, 0.0
        for cut in sorted({s for _, s in pairs}, reverse=True):
            chosen = [t for t, s in pairs if s >= cut]
            recall = sum(chosen) / len(pos)
            auprc += (recall - before) * sum(chosen) / len(chosen)
            before = recall
    return {"n": len(pairs), "tp": tp, "fp": fp, "tn": tn, "fn": fn, "accuracy": ratio(tp + tn, len(pairs)),
            "balanced_accuracy": (sens + spec) / 2, "sensitivity": sens, "specificity": spec, "precision": prec,
            "f1": 0.0 if prec + sens == 0 else 2 * prec * sens / (prec + sens),
            "mcc": 0.0 if denom == 0 else (tp * tn - fp * fn) / denom,
            "vme": ratio(fn, tp + fn), "me": ratio(fp, fp + tn), "auroc": auroc, "auprc": auprc}


def roc(truth, scores):
    pairs = [(t, s) for t, s in zip(truth, scores) if s is not None]
    n_pos = sum(1 for t, _ in pairs if t)
    n_neg = len(pairs) - n_pos
    if not n_pos or not n_neg:
        return []
    points = [[0.0, 0.0]]
    for cut in sorted({s for _, s in pairs}, reverse=True):
        chosen = [t for t, s in pairs if s >= cut]
        points.append([(len(chosen) - sum(chosen)) / n_neg, sum(chosen) / n_pos])
    assert len(points) <= 60
    return points


def rescore(report):
    """Recompute every learner row's derived fields from its scores, as grm
    would: called-only metrics, ROC, fold balanced accuracy, and no_call."""
    truth = [s["resistant"] for s in report["samples"]]
    schema = report["schema"]
    for split in report["splits"]:
        test_fold = split.get("test_fold") or FOLDS[split["scheme"]]
        for row in split["learners"]:
            name, scores = row["learner"], row["scores"]
            row["metrics"] = rates(truth, scores, name in RANKED)
            row["roc"] = roc(truth, scores) if name in RANKED else []
            if name == "catalog":
                continue
            fold_ba = []
            for k in range(split["folds"]):
                held = [i for i in range(len(truth)) if test_fold[i] == k]
                m = rates([truth[i] for i in held], [scores[i] for i in held], False)
                fold_ba.append(None if m is None else m["balanced_accuracy"])
            row["fold_balanced_accuracy"] = fold_ba
            if schema >= 2:
                row["no_call"] = {"resistant": sum(1 for t, s in zip(truth, scores) if s is None and t),
                                  "susceptible": sum(1 for t, s in zip(truth, scores) if s is None and not t)}
    return report


def feature_card(f):
    present = [(i + f) % 4 == 0 for i in range(N)]
    words = [0, 0]
    for i, p in enumerate(present):
        if p:
            words[i >> 6] |= 1 << (i & 63)
    raw = b"".join(w.to_bytes(8, "little") for w in words)
    return {"feature": f, "n_kmers": 1, "kmer": "ACGT", "catalog_families": [], "catalog_products": [], "locus": None,
            "n_resistant": sum(1 for p, t in zip(present, TRUTH) if p and t),
            "n_susceptible": sum(1 for p, t in zip(present, TRUTH) if p and not t),
            "bits": base64.b64encode(raw).decode()}


def make_report(revision=7, null_fits=NULL_FITS, forest_variant=0, lineage_available=True):
    """A report as grm writes it for this revision: schema 2 from revision 7,
    where a null fit leaves its fold uncalled; schema 1 before, where it
    scored the fold 0.5."""
    schema = 2 if revision >= 7 else 1
    catalog_calls = [t != (i % 7 == 0) for i, t in enumerate(TRUTH)]
    samples = [{"id": "g%03d" % i, "resistant": TRUTH[i], "from_mic": False, "mic": None, "lineage": LINEAGE[i],
                "country": None, "year": None, "host": None, "mlst": None,
                "groups": [0] if catalog_calls[i] else [], "catalog_resistant": catalog_calls[i]} for i in range(N)]
    splits = []
    for scheme in ("random", "lineage"):
        split = {"scheme": scheme, "folds": 5, "n_groups": 16, "learners": []}
        if scheme == "lineage" and not lineage_available:
            split["unavailable_reason"] = "fixture: lineage split unavailable"
            splits.append(split)
            continue
        if schema >= 2:
            split["test_fold"] = list(FOLDS[scheme])
        for name in LEARNERS:
            variant = forest_variant if name == "rf" else 0
            scores = [score(scheme, name, i, variant) for i in range(N)]
            complexity = [COMPLEXITY[name]] * 5
            for null_scheme, null_name, fold in null_fits:
                if (null_scheme, null_name) == (scheme, name):
                    complexity[fold] = 0.0
                    for i in range(N):
                        if FOLDS[scheme][i] == fold:
                            scores[i] = None if schema >= 2 else 0.5
            split["learners"].append({"learner": name, "metrics": None, "fold_balanced_accuracy": [], "roc": [],
                                      "scores": scores, "conformal": [], "fold_complexity": complexity, "seconds": 1.0})
        split["learners"].append({"learner": "catalog", "metrics": None, "fold_balanced_accuracy": [], "roc": [],
                                  "scores": [1.0 if c else 0.0 for c in catalog_calls], "conformal": [],
                                  "fold_complexity": [], "seconds": 0.0})
        splits.append(split)
    report = {
        "schema": schema,
        "evaluation": {"revision": revision, "valid": True, "reasons": [], "scope": SCOPE},
        "id": ID, "species": "Fixture species", "species_slug": SLUG, "antibiotic": "fixturemycin",
        "drug_class": "fixture", "breakpoint": None,
        "counts": {"samples": N, "resistant": sum(TRUTH), "susceptible": N - sum(TRUTH), "labels_from_mic": 0,
                   "features": FEATURES, "patterns": FEATURES, "kmers": FEATURES, "lineages": 16, "min_support": 3},
        "lineage_threshold": 0.01,
        "catalog_groups": ["blaFIX"],
        "samples": samples,
        "splits": splits,
        "catalog": {"version": "fixture", "rule": [["blaFIX"]], "metrics": rates(TRUTH, [1.0 if c else 0.0 for c in catalog_calls], False),
                    "groups": [{"group": "blaFIX", "class": "BETA-LACTAM", "subclass": "FIXTURE", "counts_for_rule": True,
                                "n_resistant": sum(1 for c, t in zip(catalog_calls, TRUTH) if c and t),
                                "n_susceptible": sum(1 for c, t in zip(catalog_calls, TRUTH) if c and not t)}]},
        "panels": [],
        "unsupported_assays": [],
        "mic": None,
        "models": {"rf_importance": copy.deepcopy(RF_IMPORTANCE), "gbt_importance": []},
        "features": [feature_card(w["feature"]) for w in RF_IMPORTANCE],
        "seconds": 10.0,
    }
    if schema >= 2:
        report["eval_config"] = copy.deepcopy(EVAL_CONFIG)
    return rescore(report)


def isolated_report(revision=7):
    """make_report with one k-mer-identical pair isolated as grm does it: two
    same-lineage isolates in different random folds, neither in lasso's null
    fold 4, listed in kset_duplicates, the second moved to the first's fold
    (a tie goes to the lower-index member's fold), every row rescored.
    Returns the report and the counted.json kset digests the pair shares."""
    rep = make_report(revision)
    f = FOLDS["random"]
    a, b = next((a, b) for a in range(N) for b in range(a + 1, N)
                if LINEAGE[a] == LINEAGE[b] and f[a] != f[b] and 4 not in (f[a], f[b]))
    rep["kset_duplicates"] = [[a, b]]
    rep["splits"][0]["test_fold"][b] = f[a]
    shared = "kset of " + rep["samples"][a]["id"]
    return rescore(rep), {rep["samples"][i]["id"]: shared for i in (a, b)}


def row(report, scheme, learner):
    split = next(s for s in report["splits"] if s["scheme"] == scheme)
    return next(r for r in split["learners"] if r["learner"] == learner)


def lite(learner_row, schema=2):
    """The manifest rates grm report writes for one learner row."""
    if learner_row is None:
        return None
    m = learner_row["metrics"]
    keys = ("balanced_accuracy", "sensitivity", "specificity", "vme", "me", "auroc")
    out = dict.fromkeys(keys) if m is None else {k: m[k] for k in keys}
    if schema >= 2 and "no_call" in learner_row:
        out["no_call"] = learner_row["no_call"]
    return out


def manifest(report, best, schema=2, evaluation_revision=7):
    """index.json with one dataset row naming `best` (or None) as grm report would."""
    def find(scheme, name):
        split = next((s for s in report["splits"] if s["scheme"] == scheme), None)
        return None if split is None else next((r for r in split["learners"] if r["learner"] == name), None)

    best_lineage = None
    if best is not None and find("lineage", best) is not None and find("lineage", best)["metrics"] is not None:
        best_lineage = find("lineage", best)["metrics"]["balanced_accuracy"]
    dataset = {"id": report["id"], "species": report["species"], "species_slug": SLUG, "antibiotic": report["antibiotic"],
               "drug_class": report["drug_class"], "evaluation": report["evaluation"], "n": N, "n_resistant": sum(TRUTH),
               "labels_from_mic": 0, "features": FEATURES, "lineages": 16,
               "best_learner": best, "best_random": None if best is None else find("random", best)["metrics"]["balanced_accuracy"],
               "best_lineage": best_lineage,
               "learners": [{"learner": name, "random": lite(find("random", name), schema), "lineage": lite(find("lineage", name), schema)}
                            for name in LEARNERS],
               "catalog": None, "mic_essential_agreement": None, "scm_rules": [], "panels": [], "kover": None}
    return {"schema": schema, "generated_at": "2026-10-25T00:00:00Z", "evaluation_revision": evaluation_revision,
            "datasets": [dataset], "panels": [], "unsupported_panels": []}


def done_line(learner, ba, dataset=ID):
    if learner is None:
        return "    %s done in 3s (no learner calls every isolate)" % dataset
    return "    %s done in 3s (best random-split learner among those that call every isolate: %s %.3f)" % (dataset, learner, ba)


def old_done_line(learner, ba, dataset=ID):
    return "    %s done in 3s (best random-split balanced accuracy: %s %.3f)" % (dataset, learner, ba)


def recompute_line(reason="evaluation audit revision changed", dataset=ID):
    return "    %s: %s; recomputing report" % (dataset, reason)


def write_json(path, value):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        json.dump(value, f)


def checker_root(tmp, report):
    """A data root holding one report where check-results.py looks for it."""
    write_json(os.path.join(tmp, "curated", "datasets.json"), [{"id": ID, "species_slug": SLUG}])
    write_json(os.path.join(tmp, "species", SLUG, "results", ID + ".json"), report)
    return tmp


def run_checker(root, *args):
    return subprocess.run([sys.executable, "-I", CHECKER, root, *args], capture_output=True, text=True)


def audit_root(tmp, report, index=None, backup=None, log_lines=(), published=True, digests=None):
    """A /root/grm-shaped tree for audit.py: the report, counted.json for its
    catalog calls and kset digests (each genome its own unless `digests`
    names one), the manifest and its published copy, a backup, and a log."""
    digests = digests or {}
    root = os.path.join(tmp, "root")
    path = os.path.join(root, "data", "species", SLUG, "results", ID + ".json")
    write_json(path, report)
    write_json(os.path.join(root, "data", "curated", "datasets.json"), [{"id": ID, "species_slug": SLUG}])
    write_json(os.path.join(root, "data", "species", SLUG, "counted.json"),
               {"genomes": [{"genome_id": s["id"], "genes": [{"group": "blaFIX", "family": "blaFIX"}] if s["groups"] else [],
                             "kset_digest": digests.get(s["id"], "kset of " + s["id"])}
                            for s in report["samples"]]})
    write_json(os.path.join(root, "data", "report", "index.json"), index if index is not None else manifest(report, "rf"))
    if published:
        os.makedirs(os.path.join(root, "data", "report", "datasets"), exist_ok=True)
        with open(path, "rb") as src, open(os.path.join(root, "data", "report", "datasets", ID + ".json"), "wb") as dst:
            dst.write(src.read())
    backup_dir = os.path.join(tmp, "backup")
    if backup is not None:
        write_json(os.path.join(backup_dir, "data", "species", SLUG, "results", ID + ".json"), backup)
    log = os.path.join(tmp, "pipeline.log")
    with open(log, "w") as f:
        f.write("".join(line + "\n" for line in log_lines))
    return {"root": root, "report": path, "backup": backup_dir, "log": log}
