"""The nested search and the cutoff of plan block C5 on one outer fold (docs/settings-arms.md).

  python3 -I drive.py <grm> <data root> <species> <folds file> <default pick file> <fold> <work dir> <out dir>

<folds file> is what `grm learn --exp-settings-folds` wrote for the split:
its outer folds and, inside each outer fold's training genomes T_f, the inner
folds grouped by lineage. <default pick file> is the default's result for the
same dataset, split and seed (`grm learn --exp-pick --exp-pick-inner 0`).

On outer fold <fold>, for the five learners at once:

1. checks the folds file: every genome of T_f in one inner fold, no other
   genome in any, no lineage in two inner folds, both classes in each inner
   fold's test and training genomes;
2. writes the inner fold file (test: the three inner folds; unused: every
   genome outside T_f, the outer test genomes among them) and, per learner,
   the settings file of its 31 trials (trial 0 the release's setting, 1 to
   30 drawn by scripts/hyperopt/search.py at the protocol's sampler seed),
   and runs `grm learn --fold-file --exp-settings` once per learner, so each
   learner's wall time, CPU time and peak memory are its own;
3. per learner, picks a trial by the default-anchored rule on the mean of
   its three inner balanced accuracies (`pick`) and a cutoff on trial 0's
   inner calls (`choose_cutoff`), and writes pick.json;
4. only then writes the outer fold file (test: fold <fold>; train: T_f) and
   runs, per learner, the release's setting and the picked one once on it,
   at the seed the default fitted that fold with, and checks that the
   release's calls at 0.5 are the default's calls on that fold;
5. writes <out dir>/<dataset>.<scheme>-s<seed>.f<fold>.choices.json: per
   learner the trial picked and its setting, the rule's rank counts, the
   cutoff, and the outer calls of the default, the search and the cutoff
   arm, as counts. It names no genome. <work dir>/<stem>/timing.json holds
   each grm run's wall time, CPU time and peak memory and nothing else.

A restart skips a step whose files exist. Every file under <work dir> may
name genomes and stays on the machine. Exits 1 on any failure.
"""

import hashlib
import json
import os
import re
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from hyperopt import protocol as hyper_protocol  # noqa: E402
from hyperopt import search as hyper_search  # noqa: E402

LEARNERS = ["scm", "cart", "rf", "gbt", "lasso"]
TRIALS = 30
SAMPLER_SEED = 20261030
GAMMA = 0.04
GRID = [round(0.05 * k, 2) for k in range(1, 20)]
HALF = 0.5
INNER_FOLDS = 3
REVISION = 14
FOLD_FORMAT = "grm-fold-file-v1"
FOLDS_FORMAT = "grm-c5-folds-1"
SETTINGS_FORMAT = "grm-c5-settings-1"
TRIAL_FORMAT = "grm-c5-trial-1"
PICK_FORMAT = "grm-c9-pick-1"
CHOICE_FORMAT = "grm-c5-choices-1"
PROTOCOL = {"families": hyper_protocol.FAMILIES, "seeds": {"sampler": SAMPLER_SEED}}
GENOME_ID = re.compile(r"(?<![\w.])\d+\.\d+(?!\d)")


class DriveError(RuntimeError):
    pass


def r12(x):
    return None if x is None else round(x, 12)


def setting(learner, trial):
    """Trial 0 is the release's own fit; trials 1 to 30 are the protocol's draws."""
    return "release" if trial == 0 else hyper_search.sample_config(PROTOCOL, learner, trial)


def fold_seed(seed, fold):
    """The seed the release fits outer fold `fold` with: seed + 101 * fold, wrapping at 2^64."""
    return (seed + 101 * fold) % 2**64


def inner_seed(seed, fold):
    """The base seed of the inner run of outer fold `fold`."""
    return int.from_bytes(hashlib.sha256(f"c5-inner/{seed}/{fold}".encode()).digest()[:8], "big") >> 1


def balanced(c):
    pos, neg = c["tp"] + c["fn"], c["tn"] + c["fp"]
    if not pos or not neg:
        return None
    return (c["tp"] / pos + c["tn"] / neg) / 2


def add(a, b):
    return {k: a[k] + b[k] for k in ("n", "tp", "fp", "tn", "fn", "loss_sum")}


ZERO = {"n": 0, "tp": 0, "fp": 0, "tn": 0, "fn": 0, "loss_sum": 0.0}


# ---------- the folds ----------

def check_folds(doc, fold):
    """(T_f, test genomes) of outer fold `fold`, as sample indices, after the checks of step 1."""
    if doc.get("format") != FOLDS_FORMAT or doc.get("inner_folds") != INNER_FOLDS:
        raise DriveError("not a C5 folds file with three inner folds")
    if doc.get("unavailable_reason") is not None:
        raise DriveError("the split is unavailable")
    n = len(doc["genomes"])
    for key in ("resistant", "lineage", "outer", "used"):
        if len(doc[key]) != n:
            raise DriveError(f"the folds file's {key} has {len(doc[key])} entries for {n} genomes")
    if not 0 <= fold < doc["folds"] or len(doc["inner"]) != doc["folds"]:
        raise DriveError(f"no outer fold {fold}")
    train = [i for i in range(n) if doc["used"][i] and doc["outer"][i] != fold]
    test = [i for i in range(n) if doc["outer"][i] == fold]
    if not test or not train:
        raise DriveError(f"outer fold {fold} tests or trains no genome")
    inner = doc["inner"][fold]
    at = set(train)
    for i in range(n):
        j = inner[i]
        if (i in at) != (j is not None) or (j is not None and not 0 <= j < INNER_FOLDS):
            raise DriveError(f"the inner folds of outer fold {fold} do not cover exactly its training genomes")
    lineage_fold = {}
    for i in train:
        if lineage_fold.setdefault(doc["lineage"][i], inner[i]) != inner[i]:
            raise DriveError(f"an inner fold of outer fold {fold} splits a lineage")
    for j in range(INNER_FOLDS):
        for side in (lambda i: inner[i] == j, lambda i: inner[i] != j):
            rs = [doc["resistant"][i] for i in train if side(i)]
            if not rs or all(rs) or not any(rs):
                raise DriveError(f"inner fold {j} of outer fold {fold} lacks a class")
    return train, test


def scheme_tag(doc):
    return doc["scheme"]


def inner_file(doc, fold, train):
    ids = doc["genomes"]
    at = set(train)
    return {"format": FOLD_FORMAT, "dataset": doc["dataset"],
            "scheme": f"c5i-{scheme_tag(doc)}-s{doc['seed']}-f{fold}",
            "test": [[ids[i] for i in train if doc["inner"][fold][i] == j] for j in range(INNER_FOLDS)],
            "train": [],
            "unused": [ids[i] for i in range(len(ids)) if i not in at]}


def outer_file(doc, fold, train, test):
    ids = doc["genomes"]
    used = set(train) | set(test)
    return {"format": FOLD_FORMAT, "dataset": doc["dataset"],
            "scheme": f"c5o-{scheme_tag(doc)}-s{doc['seed']}-f{fold}",
            "test": [[ids[i] for i in test]],
            "train": [ids[i] for i in train],
            "unused": [ids[i] for i in range(len(ids)) if i not in used]}


def tuning(doc):
    """The lineage split's fits tune on lineages, a fold file's on copy groups, as the release's do."""
    return "lineage" if doc["scheme"] == "lineage" else "copies"


# ---------- the rule ----------

def history(docs, learner):
    """{trial: mean inner balanced accuracy, or None when the fit failed or an inner fold has no figure}."""
    out = {}
    for t, doc in docs.items():
        row = next(r for r in doc["learners"] if r["learner"] == learner)
        bas = row["fold_balanced_accuracy"]
        ok = row["status"] == "ok" and len(bas) == INNER_FOLDS and all(b is not None for b in bas)
        out[t] = r12(sum(bas) / INNER_FOLDS) if ok else None
    return out


def pick(means, gamma=GAMMA):
    """The default-anchored rule: keep trial 0 unless the best inner mean beats it by more than gamma.

    A tie at the best goes to the lowest trial id. Every mean and difference
    is compared at 12 decimals. The plain best is returned beside the pick.
    """
    ok = {t: m for t, m in means.items() if m is not None}
    counts = {"trials": len(means), "failed": len(means) - len(ok)}
    if 0 not in ok:
        return dict(counts, picked=0, plain_best=None, tied=[], status="release failed inside", release_mean=None,
                    best_mean=None, above_gamma=0, better=0, same=0, worse=0)
    best = max(ok.values())
    tied = sorted(t for t, m in ok.items() if m == best)
    lead = r12(best - ok[0])
    above = sum(1 for t, m in ok.items() if t and r12(m - ok[0]) > gamma)
    return dict(counts, picked=tied[0] if lead > gamma else 0, plain_best=tied[0], tied=tied, status="picked",
                release_mean=ok[0], best_mean=best, above_gamma=above,
                better=sum(1 for t, m in ok.items() if t and m > ok[0]),
                same=sum(1 for t, m in ok.items() if t and m == ok[0]),
                worse=sum(1 for t, m in ok.items() if t and m < ok[0]))


def grid_calls(doc, learner, folds=None):
    """{cutoff: calls summed over the result's folds} of one learner, or None when its fit failed."""
    row = next(r for r in doc["learners"] if r["learner"] == learner)
    if row["status"] != "ok":
        return None
    out = {}
    for k in row["calls"]:
        c = dict(ZERO)
        for f, calls in enumerate(k["folds"]):
            if folds is None or f in folds:
                c = add(c, calls)
        out[k["cutoff"]] = c
    return out


def choose_cutoff(by_cutoff):
    """(cutoff, its balanced accuracy, 0.5's) on the grid: the best pooled balanced accuracy at 12 decimals,
    a tie to the point nearest 0.5 and, of two equally near, the lower."""
    if by_cutoff is None:
        return HALF, None, None
    scored = {c: r12(balanced(by_cutoff[c])) for c in GRID if c in by_cutoff and balanced(by_cutoff[c]) is not None}
    if not scored:
        return HALF, None, None
    best = max(scored.values())
    cut = min((c for c, b in scored.items() if b == best), key=lambda c: (round(abs(c - HALF), 12), c))
    return cut, best, scored.get(HALF)


# ---------- grm ----------

def run_grm(grm, root, species, fold_file, settings, out, log, timing, step):
    """One grm learn --fold-file --exp-settings run; its wall, CPU and peak memory go to `timing`."""
    out.mkdir(parents=True, exist_ok=True)
    cmd = [str(grm), "--data", str(root), "learn", "--species", species, "--fold-file", str(fold_file),
           "--fold-out", str(out), "--exp-settings", str(settings), "--force"]
    t0 = time.monotonic()
    with open(log, "a") as f:
        f.write(f"{time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())} {step}: grm learn\n")
        f.flush()
        p = subprocess.Popen(cmd, stdout=f, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL)
        _, status, usage = os.wait4(p.pid, 0)
    rc = p.returncode = os.waitstatus_to_exitcode(status)
    timing.append({"step": step, "wall_s": round(time.monotonic() - t0, 3), "user_s": round(usage.ru_utime, 3),
                   "sys_s": round(usage.ru_stime, 3), "peak_kb": usage.ru_maxrss, "rc": rc})
    if rc != 0:
        raise DriveError(f"{step}: grm learn exited {rc}")


def read_trials(out, dataset, scheme, ids):
    docs = {}
    for t in ids:
        path = out / f"{dataset}.{scheme}.t{t}.json"
        with open(path) as f:
            doc = json.load(f)
        if (doc.get("format"), doc.get("evaluation_revision"), doc.get("dataset"), doc.get("scheme"), doc.get("trial")) \
                != (TRIAL_FORMAT, REVISION, dataset, scheme, t):
            raise DriveError(f"{path.name}: not the trial result asked for")
        docs[t] = doc
    return docs


def write_json(path, doc):
    tmp = path.with_name(path.name + ".tmp")
    with open(tmp, "w") as f:
        json.dump(doc, f, sort_keys=True)
    os.replace(tmp, path)


def load(path):
    with open(path) as f:
        return json.load(f)


# ---------- one outer fold ----------

def default_fold(pick_doc, folds_doc, fold, test, train):
    """The default's calls per learner on the fold, after checking the pick file is the same split and fold."""
    if (pick_doc.get("format"), pick_doc.get("dataset"), pick_doc.get("scheme"), pick_doc.get("seed")) \
            != (PICK_FORMAT, folds_doc["dataset"], folds_doc["scheme"], folds_doc["seed"]):
        raise DriveError("the default's pick file is not this dataset, split and seed")
    if pick_doc.get("unavailable_reason") is not None:
        raise DriveError("the default's split is unavailable")
    rec = pick_doc["folds"][fold]
    r = folds_doc["resistant"]
    shape = (len(train), sum(r[i] for i in train), len(test), sum(r[i] for i in test))
    if rec["fold"] != fold or (rec["train"], rec["train_resistant"], rec["test"], rec["test_resistant"]) != shape:
        raise DriveError(f"the default's fold {fold} is not the folds file's")
    out = {c["learner"]: {k: c[k] for k in ZERO} for c in rec["outer"]}
    if list(out) != LEARNERS:
        raise DriveError("the default's pick file does not hold the five learners in order")
    return out


def same_calls(a, b):
    return all(a[k] == b[k] for k in ("n", "tp", "fp", "tn", "fn")) and abs(a["loss_sum"] - b["loss_sum"]) <= 1e-9


def inner_settings(learner, seed, tune, ids):
    return {"format": SETTINGS_FORMAT, "seed": seed, "tuning": tune, "cutoffs": GRID,
            "trials": [{"id": t, "learners": {learner: setting(learner, t)}} for t in ids]}


def outer_settings(learner, c, seed, tune):
    """The release's setting (trial 0) with its cutoff, and the picked setting when it is another."""
    trials = [{"id": 0, "learners": {learner: "release"}, "cutoffs": {learner: c["cutoff"]}}]
    if c["picked"] != 0:
        trials.append({"id": c["picked"], "learners": {learner: setting(learner, c["picked"])}})
    # Calls at 0.5 and at the cutoff alone: the outer fold gives no figure at another grid point.
    return {"format": SETTINGS_FORMAT, "seed": seed, "tuning": tune, "cutoffs": [], "trials": trials}


def outer_calls(docs, trial, learner, cutoff):
    if trial not in docs:
        return None
    row = next((r for r in docs[trial]["learners"] if r["learner"] == learner), None)
    if row is None or row["status"] != "ok":
        return None
    for k in row["calls"]:
        if k["cutoff"] == cutoff:
            return k["folds"][0]
    raise DriveError(f"trial {trial} {learner}: no calls at {cutoff}")


def drive(grm, root, species, folds_path, pick_path, fold, work, out, trials=TRIALS):
    folds_doc, pick_doc = load(folds_path), load(pick_path)
    dataset, scheme, seed = folds_doc["dataset"], folds_doc["scheme"], folds_doc["seed"]
    stem = f"{dataset}.{scheme}-s{seed}.f{fold}"
    final = out / f"{stem}.choices.json"
    if final.is_file():
        return final
    train, test = check_folds(folds_doc, fold)
    default = default_fold(pick_doc, folds_doc, fold, test, train)
    tune = tuning(folds_doc)
    here = work / stem
    here.mkdir(parents=True, exist_ok=True)
    log, timing = here / "grm.log", []
    ids = list(range(trials + 1))

    timing_path = here / "timing.json"
    if timing_path.is_file():
        timing = load(timing_path)

    def grm_step(step, fold_file, settings_file, result_dir):
        if (here / f"{step}.done").is_file():
            return
        run_grm(grm, root, species, fold_file, settings_file, result_dir, log, timing, step)
        write_json(timing_path, timing)
        (here / f"{step}.done").write_text("")

    # Inner: T_f alone; the outer test genomes are unused.
    inner = inner_file(folds_doc, fold, train)
    write_json(here / "inner.folds.json", inner)
    inner_docs = {}
    for l in LEARNERS:
        write_json(here / f"inner-{l}.settings.json", inner_settings(l, inner_seed(seed, fold), tune, ids))
        grm_step(f"inner-{l}", here / "inner.folds.json", here / f"inner-{l}.settings.json", here / f"inner-{l}")
        inner_docs[l] = read_trials(here / f"inner-{l}", dataset, inner["scheme"], ids)

    # The pick and the cutoff, from inner scores alone, written before any outer fit.
    choice = {}
    for l in LEARNERS:
        rule = pick(history(inner_docs[l], l))
        cut, cut_ba, half_ba = choose_cutoff(grid_calls(inner_docs[l][0], l))
        choice[l] = dict(rule, cutoff=cut, cutoff_inner_balanced_accuracy=cut_ba,
                         half_inner_balanced_accuracy=half_ba)
    write_json(here / "pick.json", choice)

    # Outer: fold `fold` tested once, T_f trained, at the default's seed for that fold.
    outer = outer_file(folds_doc, fold, train, test)
    write_json(here / "outer.folds.json", outer)
    outer_docs = {}
    for l in LEARNERS:
        settings = outer_settings(l, choice[l], fold_seed(seed, fold), tune)
        write_json(here / f"outer-{l}.settings.json", settings)
        grm_step(f"outer-{l}", here / "outer.folds.json", here / f"outer-{l}.settings.json", here / f"outer-{l}")
        outer_docs[l] = read_trials(here / f"outer-{l}", dataset, outer["scheme"], [t["id"] for t in settings["trials"]])

    rows = []
    for l in LEARNERS:
        c = choice[l]
        release = outer_calls(outer_docs[l], 0, l, HALF)
        if release is None or not same_calls(release, default[l]):
            raise DriveError(f"{l}: the release's refit on fold {fold} is not the default's fit")
        search, why = release, None
        if c["picked"] != 0:
            search = outer_calls(outer_docs[l], c["picked"], l, HALF)
            if search is None:
                search, why = release, "the picked setting's outer fit failed; the release's is served"
        rows.append({
            "learner": l, "status": c["status"], "trials": c["trials"], "failed_trials": c["failed"],
            "picked_trial": c["picked"], "picked_setting": setting(l, c["picked"]),
            "plain_best_trial": c["plain_best"], "tied_at_best": len(c["tied"]),
            "release_inner_mean": c["release_mean"], "best_inner_mean": c["best_mean"],
            "trials_above_gamma": c["above_gamma"], "trials_better": c["better"], "trials_same": c["same"],
            "trials_worse": c["worse"], "cutoff": c["cutoff"],
            "cutoff_inner_balanced_accuracy": c["cutoff_inner_balanced_accuracy"],
            "half_inner_balanced_accuracy": c["half_inner_balanced_accuracy"], "search_note": why,
            "outer": {"default": default[l], "search": search, "cutoff": outer_calls(outer_docs[l], 0, l, c["cutoff"])},
            "seconds_inner": round(sum(d["learners"][0]["seconds"] for d in inner_docs[l].values()), 3),
            "seconds_outer": round(sum(d["learners"][0]["seconds"] for d in outer_docs[l].values()), 3),
        })
    r = folds_doc["resistant"]
    doc = {"format": CHOICE_FORMAT, "evaluation_revision": REVISION, "dataset": dataset, "scheme": scheme,
           "seed": seed, "fold": fold, "tuning": tune, "inner_seed": inner_seed(seed, fold),
           "outer_seed": fold_seed(seed, fold), "sampler_seed": SAMPLER_SEED, "gamma": GAMMA, "grid": GRID,
           "random_trials": trials,
           "counts": {"train": len(train), "train_resistant": sum(r[i] for i in train), "test": len(test),
                      "test_resistant": sum(r[i] for i in test)},
           "inner_folds": inner_docs["scm"][0]["folds"], "min_support_inner": inner_docs["scm"][0]["min_support"],
           "min_support_outer": outer_docs["scm"][0]["min_support"], "learners": rows, "timing": timing}
    text = json.dumps(doc, sort_keys=True)
    ids_here = set(folds_doc["genomes"])
    if any(m.group(0) in ids_here for m in GENOME_ID.finditer(text)):
        raise DriveError("the choices file would name a genome")
    out.mkdir(parents=True, exist_ok=True)
    write_json(final, doc)
    return final


def main(argv):
    if len(argv) != 8:
        print(__doc__.strip().splitlines()[2], file=sys.stderr)
        return 2
    grm, root, species, folds_path, pick_path, fold, work, out = argv
    try:
        path = drive(Path(grm), Path(root), species, Path(folds_path), Path(pick_path), int(fold), Path(work),
                     Path(out))
    except DriveError as e:
        print(f"drive: {e}", file=sys.stderr)
        return 1
    print(f"drive: {path.name}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
