"""combiner.py: can a combination of the five learners recover the misses of their mean? The reading frozen in
FREEZE-COMBINE (docs/combiner-result.md).

    python3 -I combiner.py plan
    python3 -I combiner.py read REPORTS_DIR --first FIRST --second SECOND --out READING_DIR [--expect 274]
    python3 -I combiner.py write READING_DIR --reports REPORTS_DIR --out RECORD_DIR
    python3 -I combiner.py print RECORD_DIR

FIRST is the learner overlap's record (reference/learner-overlap, its datasets.tsv) and SECOND its second pass's
(reference/learner-overlap/second-pass, its vote.tsv). The reports read are the ones the learner overlap pooled:
the five learners on the lineage split and 10 or more resistant genomes.

Three forms are read, and only these. The incumbent is the equal mean of the five stored held-out scores, called
when their math.fsum is 2.5 or more, as overlap_second.py calls it. The weighted mean has five non-negative weights
summing to one. The logistic combiner has an intercept and five coefficients on the five scores, with no
interaction term. Each arm is fitted per report and per stored outer fold on the genomes of the other folds, by
plain mean log loss, and applied to the held-out fold. Every form is called at 0.5. Nothing else enters a fit: no
catalogue call, lineage, species, year or country.

`read` runs the learner overlap's checks through learner_overlap.read_reports, which drops every key named `id`
as a file is parsed, then gates 1 to 6, and writes reading.json: per report, counts only. `write` runs gates 7 and
8 and writes datasets.tsv, pooled.tsv, rule.tsv, transitions.tsv, gates.txt and SHA256SUMS. `print` prints the
page's tables. A failed gate writes nothing and exits 1.

Python 3.12, standard library only; run with -I.
"""

import argparse
import collections
import fractions
import json
import math
import os
import shutil
import statistics
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common  # noqa: E402
import learner_overlap as lo  # noqa: E402

CUTOFF = 0.5                    # every form: a score of 0.5 or more is a resistant call
FIVE = lo.FIVE
VOTE_SUM = 2.5                  # the incumbent's call: math.fsum of the five scores, as overlap_second.py makes it
FORMS = ("mean", "weighted", "logistic")
ARMS = ("weighted", "logistic")
CLIP = 1e-4                     # log loss clips a score to [CLIP, 1 - CLIP] (src/cutoff.rs LOSS_CLIP)
MARGIN = common.MARGIN          # 0.02
W_CAP, W_TOL = 5000, 1e-9       # the weighted fit: at most W_CAP projected gradient steps, to a projected step of W_TOL
L_CAP, L_TOL, L_STEP = 100, 1e-10, 1e-8   # the logistic fit: at most L_CAP Newton steps, to a gradient and a step under these
PIVOT = 1e-13                   # a Newton system with a pivot under this has no solution: the fit has not converged
SUM_TOL = 1e-12                 # a fold's weights sum to 1 within this
MAX_LIST = 5                    # no list in a reading is longer than this, but the report list and the gate lines
MEASURES = ("balanced_accuracy", "recall", "specificity", "vme", "log_loss")
LOWER_BETTER = ("vme", "log_loss")
DECIMAL_LIMIT = 10              # a decimal cell under 10 in size has one digit before its point
READING = "reading.json"
RECORD_FILES = ("datasets.tsv", "pooled.tsv", "rule.tsv", "transitions.tsv", "gates.txt")
TRANSITIONS = ("recovered_1to4", "recovered_all5", "lost", "cleared", "new_wrong")
COUNTS = ("tp", "fn", "fp", "tn", "called")
DATASET_COLUMNS = (["species", "antibiotic", "genomes", "resistant", "susceptible", "mean_miss_by_1to4",
                    "mean_miss_by_all5"]
                   + [f"{f}_{c}" for f in FORMS for c in COUNTS + ("log_loss",)]
                   + [f"{a}_{t}" for a in ARMS for t in TRANSITIONS]
                   + [f"{a}_folds_converged" for a in ARMS])
POOLED_COLUMNS = ["scope", "species", "arm", "measure", "datasets", "median", "lo", "hi", "better", "worse", "same"]
RULE_COLUMNS = ["arm", "datasets", "species", "median_gain", "gain_lo", "gain_hi", "vme_median", "vme_lo", "vme_hi",
                "gain_at_least_margin", "reading"]
TRANSITION_COLUMNS = ["scope", "species", "arm", "datasets", "resistant", "mean_miss", "mean_miss_by_1to4",
                      "mean_miss_by_all5", "susceptible", "mean_wrong"] + list(TRANSITIONS)

Stop = lo.Stop
load_report = lo.load_report     # drops every key named `id` as the file is parsed


# ---- the forms -----------------------------------------------------------------------------------------------------

def is_call(score):
    return score >= CUTOFF


def mean_call(scores):
    """The incumbent: the five scores' math.fsum is 2.5 or more."""
    return math.fsum(scores) >= VOTE_SUM


def mean_score(scores):
    return math.fsum(scores) / len(scores)


def clipped(p):
    return min(max(p, CLIP), 1.0 - CLIP)


def log_loss(labels, probs):
    """Plain mean log loss, every score clipped to [CLIP, 1 - CLIP]."""
    return math.fsum(-math.log(clipped(p)) if t else -math.log(1.0 - clipped(p))
                     for t, p in zip(labels, probs)) / len(labels)


def groups_of(S, y, keep):
    """The genomes in `keep` as (scores, label, count): genomes with the same five scores and label are one term of
    the loss, so the fit is the same and faster."""
    return [(s, t, k) for (s, t), k in sorted(collections.Counter((S[i], y[i]) for i in keep).items())]


def project_simplex(v):
    """The point of the simplex (non-negative, summing to one) nearest v."""
    u = sorted(v, reverse=True)
    run, theta = 0.0, 0.0
    for j, x in enumerate(u, 1):
        run += x
        t = (run - 1.0) / j
        if x - t > 0:
            theta = t
    w = [max(x - theta, 0.0) for x in v]
    total = math.fsum(w)
    return [x / total for x in w]


def weighted_score(w, s):
    """The weighted mean, held inside the five scores' own range so that rounding cannot cross the bound."""
    return min(max(math.fsum(a * b for a, b in zip(w, s)), min(s)), max(s))


def weighted_loss_grad(w, groups):
    """Plain mean log loss of the weighted mean, clipped, and its gradient in the weights."""
    n = sum(k for _, _, k in groups)
    total, grad = 0.0, [0.0] * len(w)
    for s, t, k in groups:
        p = math.fsum(a * b for a, b in zip(w, s))
        c = clipped(p)
        total += k * (-math.log(c) if t else -math.log(1.0 - c))
        if CLIP < p < 1.0 - CLIP:
            d = k * (-1.0 / p if t else 1.0 / (1.0 - p))
            for j in range(len(w)):
                grad[j] += d * s[j]
    return total / n, [g / n for g in grad]


def projected_step(w, g):
    return max(abs(a - b) for a, b in zip(project_simplex([x - d for x, d in zip(w, g)]), w))


def fit_weighted(groups):
    """Five non-negative weights summing to one, by plain log loss: projected gradient with a Barzilai-Borwein step
    and backtracking, from the equal weights. Returns (weights, steps, converged)."""
    w = [1.0 / len(FIVE)] * len(FIVE)
    f, g = weighted_loss_grad(w, groups)
    alpha = 1.0
    for step in range(W_CAP):
        if projected_step(w, g) <= W_TOL:
            return w, step, True
        while True:
            trial = project_simplex([x - alpha * d for x, d in zip(w, g)])
            ft, gt = weighted_loss_grad(trial, groups)
            if ft <= f + 1e-4 * math.fsum(d * (a - b) for d, a, b in zip(g, trial, w)) or alpha < 1e-20:
                break
            alpha /= 2
        if trial == w:
            return w, step, False
        s_ = [a - b for a, b in zip(trial, w)]
        y_ = [a - b for a, b in zip(gt, g)]
        sy = math.fsum(a * b for a, b in zip(s_, y_))
        alpha = min(max(math.fsum(a * a for a in s_) / sy, 1e-10), 1e10) if sy > 0 else 1.0
        w, f, g = trial, ft, gt
    return w, W_CAP, projected_step(w, g) <= W_TOL


def sigmoid(z):
    if z >= 0:
        return 1.0 / (1.0 + math.exp(-z))
    e = math.exp(z)
    return e / (1.0 + e)


def logistic_score(beta, s):
    return sigmoid(beta[0] + math.fsum(b * x for b, x in zip(beta[1:], s)))


def logistic_nll(beta, groups):
    n = sum(k for _, _, k in groups)
    total = 0.0
    for s, t, k in groups:
        z = beta[0] + math.fsum(b * v for b, v in zip(beta[1:], s))
        total += k * ((max(z, 0.0) + math.log1p(math.exp(-abs(z)))) - (z if t else 0.0))
    return total / n


def solve(H, g):
    """H d = g by Gaussian elimination with partial pivoting; None when a pivot is under PIVOT."""
    n = len(g)
    a = [row[:] + [g[i]] for i, row in enumerate(H)]
    for col in range(n):
        piv = max(range(col, n), key=lambda r: abs(a[r][col]))
        if abs(a[piv][col]) < PIVOT:
            return None
        a[col], a[piv] = a[piv], a[col]
        for r in range(col + 1, n):
            m = a[r][col] / a[col][col]
            for c in range(col, n + 1):
                a[r][c] -= m * a[col][c]
    d = [0.0] * n
    for r in range(n - 1, -1, -1):
        d[r] = (a[r][n] - math.fsum(a[r][c] * d[c] for c in range(r + 1, n))) / a[r][r]
    return d


def fit_logistic(groups):
    """An intercept and five coefficients by maximum likelihood (plain log loss, no penalty): Newton's method with
    step halving, from zero. Converged when the mean gradient is under L_TOL and the last step under L_STEP times
    (1 + the largest coefficient). Returns (beta, steps, converged)."""
    n = sum(k for _, _, k in groups)
    size = 1 + len(FIVE)
    beta = [0.0] * size
    nll = logistic_nll(beta, groups)
    for step in range(L_CAP):
        grad, H = [0.0] * size, [[0.0] * size for _ in range(size)]
        for s, t, k in groups:
            x = (1.0,) + tuple(s)
            p = sigmoid(math.fsum(b * v for b, v in zip(beta, x)))
            r, wgt = k * (p - (1.0 if t else 0.0)), k * p * (1.0 - p)
            for i in range(size):
                grad[i] += r * x[i]
                for j in range(i, size):
                    H[i][j] += wgt * x[i] * x[j]
        grad = [v / n for v in grad]
        for i in range(size):
            for j in range(i, size):
                H[i][j] /= n
                H[j][i] = H[i][j]
        d = solve(H, grad)
        if d is None:
            return beta, step, False
        t_ = 1.0
        while True:
            trial = [b - t_ * v for b, v in zip(beta, d)]
            nt = logistic_nll(trial, groups)
            if nt <= nll or t_ < 1e-12:
                break
            t_ /= 2
        moved = max(abs(t_ * v) for v in d)
        beta, nll = trial, nt
        if max(abs(v) for v in grad) <= L_TOL and moved <= L_STEP * (1.0 + max(abs(b) for b in beta)):
            return beta, step + 1, True
    return beta, L_CAP, False


def fit_and_score(arm):
    """The arm's fit and its score, looked up when called."""
    return (fit_weighted, weighted_score) if arm == "weighted" else (fit_logistic, logistic_score)


# ---- cross-fitting over the stored folds ---------------------------------------------------------------------------

def training_genomes(folds, f):
    """The genomes an arm is fitted on for outer fold f: those of every other fold."""
    return [i for i, g in enumerate(folds) if g != f]


def crossfit(arm, S, y, folds):
    """Each fold's genomes scored by the arm fitted on the other folds' genomes. Returns the scores and, per fold, its
    record: the genomes trained on and held out, how many are both, the fit's steps and convergence, and for the
    weighted mean the smallest weight and the weights' sum."""
    fit, score = fit_and_score(arm)
    out = [None] * len(y)
    records = []
    for f in sorted(set(folds)):
        train = training_genomes(folds, f)
        test = [i for i, g in enumerate(folds) if g == f]
        model, steps, ok = fit(groups_of(S, y, train))
        for i in test:
            out[i] = score(model, S[i])
        rec = {"fold": f, "trained": len(train), "tested": len(test), "shared": len(set(train) & set(test)),
               "steps": steps, "converged": ok}
        if arm == "weighted":
            rec["weight_min"], rec["weight_sum"] = min(model), math.fsum(model)
        records.append(rec)
    return out, records


def fold_ok(rec, genomes):
    """Gate 2 on one fold's fit: no held-out genome among the training genomes, and the two add up to the report."""
    return rec["shared"] == 0 and rec["trained"] + rec["tested"] == genomes


# ---- counts --------------------------------------------------------------------------------------------------------

def confusion(labels, calls):
    """Counts at the cutoff over the genomes given a call; `called` is how many were given one."""
    c = {"tp": 0, "fn": 0, "fp": 0, "tn": 0, "called": 0}
    for t, k in zip(labels, calls):
        if k is None:
            continue
        c["called"] += 1
        c[("tp" if k else "fn") if t else ("fp" if k else "tn")] += 1
    return c


def coverage(c, genomes):
    """The share of a report's genomes given a call: never written as, or in place of, an accuracy."""
    return fractions.Fraction(c["called"], genomes)


def transitions(labels, mean_calls, arm_calls, called_by):
    """How an arm's calls part from the incumbent's. Resistant genomes the incumbent misses that the arm calls, where
    1 to 4 of the five call them and where none does; resistant genomes it calls that the arm misses; susceptible
    genomes it calls resistant that the arm does not; susceptible genomes the arm newly calls resistant."""
    out = {t: 0 for t in TRANSITIONS}
    for t, m, a, k in zip(labels, mean_calls, arm_calls, called_by):
        if t and not m and a:
            out["recovered_all5" if k == 0 else "recovered_1to4"] += 1
        elif t and m and not a:
            out["lost"] += 1
        elif not t and m and not a:
            out["cleared"] += 1
        elif not t and not m and a:
            out["new_wrong"] += 1
    return out


def report_counts(labels, S, folds):
    """Every count of one report, and its gate counts."""
    n = len(labels)
    called_by = [sum(1 for x in s if is_call(x)) for s in S]
    mean_calls = [mean_call(s) for s in S]
    probs = {"mean": [mean_score(s) for s in S]}
    calls = {"mean": mean_calls}
    gate = {"leaky": 0, "bound_called": 0, "bound_missed": 0, "bad_weights": 0, "uncovered": 0,
            "not_converged": {a: 0 for a in ARMS}, "fits": {a: 0 for a in ARMS}}
    row = {"genomes": n, "resistant": sum(labels), "susceptible": n - sum(labels),
           "mean_miss_by_1to4": sum(1 for t, m, k in zip(labels, mean_calls, called_by) if t and not m and k > 0),
           "mean_miss_by_all5": sum(1 for t, m, k in zip(labels, mean_calls, called_by) if t and not m and k == 0)}
    for arm in ARMS:
        p, recs = crossfit(arm, S, labels, folds)
        probs[arm], calls[arm] = p, [None if x is None else is_call(x) for x in p]
        for r in recs:
            gate["fits"][arm] += 1
            gate["leaky"] += not fold_ok(r, n)
            gate["not_converged"][arm] += not r["converged"]
            if arm == "weighted":
                gate["bad_weights"] += r["weight_min"] < 0 or abs(r["weight_sum"] - 1.0) > SUM_TOL
        row[f"{arm}_folds"] = recs
    for k, a in zip(called_by, calls["weighted"]):
        gate["bound_called"] += k == 0 and a is True
        gate["bound_missed"] += k == len(FIVE) and a is False
    for form in FORMS:
        c = confusion(labels, calls[form])
        if coverage(c, n) != 1:
            gate["uncovered"] += 1
        given = [(t, x) for t, x in zip(labels, probs[form]) if x is not None]
        c["log_loss"] = log_loss([t for t, _ in given], [x for _, x in given]) if given else None
        row[form] = c
    for arm in ARMS:
        row[f"{arm}_transitions"] = transitions(labels, mean_calls, calls[arm], called_by)
    return row, gate


# ---- measures from counts ------------------------------------------------------------------------------------------

def measure(c, name):
    """A measure of one form on one report, worked from its counts; log loss is stored as computed."""
    if name == "log_loss":
        return float(c["log_loss"])
    tp, fn, fp, tn = (int(c[k]) for k in ("tp", "fn", "fp", "tn"))
    recall, spec = fractions.Fraction(tp, tp + fn), fractions.Fraction(tn, tn + fp)
    return float({"balanced_accuracy": (recall + spec) / 2, "recall": recall, "specificity": spec,
                  "vme": fractions.Fraction(fn, tp + fn)}[name])


def form_counts(row, form):
    return {k: row[f"{form}_{k}"] for k in COUNTS + ("log_loss",)}


def difference(row, arm, name):
    """The arm's measure less the incumbent's, on one report, at 12 decimals."""
    return common.r12(measure(form_counts(row, arm), name) - measure(form_counts(row, "mean"), name))


def interval(by_species, name):
    """The species-clustered 95% interval of the median: resamples of species, seeded by the row's own name."""
    return common.species_interval(by_species, name)


def side(d, name):
    """'better', 'worse' or 'same' for the arm, at the margin."""
    if name in LOWER_BETTER:
        d = -d
    return "better" if d >= MARGIN else ("worse" if d <= -MARGIN else "same")


def pooled_row(rows, arm, name, scope, species=""):
    diffs = [difference(r, arm, name) for r in rows]
    out = {"scope": scope, "species": species, "arm": arm, "measure": name, "datasets": len(rows),
           "median": common.r12(statistics.median(diffs)), "lo": None, "hi": None}
    if scope == "all":
        by = collections.defaultdict(list)
        for r, d in zip(rows, diffs):
            by[r["species"]].append(d)
        out["lo"], out["hi"] = interval(dict(by), f"combiner|{arm}|{name}")
    for k in ("better", "worse", "same"):
        out[k] = sum(1 for d in diffs if side(d, name) == k)
    return out


def keep(median_gain, gain_lo, vme_lo):
    """FREEZE-COMBINE's rule: kept only if the median gain in balanced accuracy is at least the margin, its interval
    lies wholly above 0, and the very major error rate is not higher on its own interval."""
    return (median_gain >= MARGIN and gain_lo is not None and gain_lo > 0
            and not (vme_lo is not None and vme_lo > 0))


def rule(rows, arm):
    ba = pooled_row(rows, arm, "balanced_accuracy", "all")
    vme = pooled_row(rows, arm, "vme", "all")
    return {"arm": arm, "datasets": len(rows), "species": len({r["species"] for r in rows}),
            "median_gain": ba["median"], "gain_lo": ba["lo"], "gain_hi": ba["hi"], "vme_median": vme["median"],
            "vme_lo": vme["lo"], "vme_hi": vme["hi"], "gain_at_least_margin": ba["better"],
            "reading": "kept" if keep(ba["median"], ba["lo"], vme["lo"]) else "no difference shown"}


# ---- the read ------------------------------------------------------------------------------------------------------

def d1_reports(first_record):
    """The reports the learner overlap pooled, by its own rule."""
    out = set()
    for r in lo.read_tsv(os.path.join(first_record, "datasets.tsv")):
        if lo.pooled_in(dict(r, learner_count=int(r["learner_count"]), resistant=int(r["resistant"]))):
            out.add((r["species"], r["antibiotic"]))
    return out


def same_reports(wanted, read_keys):
    return wanted == read_keys


def incumbent_agrees(row, vote):
    """The incumbent counted again equals the second pass's vote.tsv row."""
    c = row["mean"]
    return (vote is not None and int(vote["mean_miss"]) == c["fn"] and int(vote["mean_wrong"]) == c["fp"]
            and int(vote["mean_miss_by_1to4"]) == row["mean_miss_by_1to4"]
            and int(vote["mean_miss_by_all5"]) == row["mean_miss_by_all5"])


def read(directory, first_record, second_record, expect=lo.EXPECT_FILES):
    """The learner overlap's checks, then gates 1 to 6; return the reading or raise Stop."""
    first = lo.read_reports(directory, expect)            # raises Stop on any of the learner overlap's checks
    wanted = d1_reports(first_record)
    used = {(d["species"], d["antibiotic"]) for d in first["datasets"]
            if lo.is_five(d) and d["resistant"] >= lo.MIN_RESISTANT}
    votes = {(r["species"], r["antibiotic"]): r for r in lo.read_tsv(os.path.join(second_record, "vote.tsv"))}
    total = collections.Counter()
    fits, unconverged = collections.Counter(), collections.Counter()
    rows, differ = [], []
    for name in lo.report_files(directory):
        report = load_report(os.path.join(directory, name))
        key = (report.get("species"), report.get("antibiotic"))
        if key not in wanted or key not in used:
            continue
        split = lo.lineage_splits(report)[0]
        entries = {lr["learner"]: lr["scores"] for lr in lo.learner_entries(split)}
        labels = [s.get("resistant") for s in report.get("samples") or []]
        S = [tuple(entries[n][i] for n in FIVE) for i in range(len(labels))]
        row, gate = report_counts(labels, S, split["test_fold"])
        del report, entries, S
        for k in ("leaky", "bound_called", "bound_missed", "bad_weights", "uncovered"):
            total[k] += gate[k]
        fits.update(gate["fits"])
        unconverged.update(gate["not_converged"])
        if not incumbent_agrees(row, votes.get(key)):
            differ.append(f"{key[0]} with {key[1]}")
        rows.append(dict(row, species=key[0], antibiotic=key[1]))
    keys = {(r["species"], r["antibiotic"]) for r in rows}
    lines = [
        f"gate 1 reports: {first['found']} files, the learner overlap's checks passed, {len(rows)} reports read, "
        f"{len(wanted)} pooled in the first record, {len(wanted ^ keys)} in one and not the other: "
        + ("pass" if same_reports(wanted, keys) and rows else "STOP"),
        f"gate 2 cross-fitting: {sum(fits.values())} fold fits, {total['leaky']} whose training and held-out genomes "
        f"overlap or do not add up to the report: "
        + ("pass" if total["leaky"] == 0 and sum(fits.values()) > 0 else "STOP"),
        f"gate 3 incumbent: {len(rows)} reports, {len(differ)} whose mean-of-five counts differ from vote.tsv: "
        + ("pass" if not differ else "STOP"),
        f"gate 4 convex bound: {total['bound_called']} genomes all five miss that the weighted mean calls, "
        f"{total['bound_missed']} genomes all five call that it misses, {total['bad_weights']} folds whose weights "
        f"are negative or do not sum to 1: "
        + ("pass" if total["bound_called"] == total["bound_missed"] == total["bad_weights"] == 0 else "STOP"),
        f"gate 5 convergence: {fits['weighted']} weighted fits, {unconverged['weighted']} not converged; "
        f"{fits['logistic']} logistic fits, {unconverged['logistic']} not converged: "
        + ("pass" if unconverged["weighted"] == unconverged["logistic"] == 0 else "STOP"),
        f"gate 6 cutoff: every form called at {CUTOFF}; {total['uncovered']} form and report pairs with a genome "
        f"given no call: " + ("pass" if total["uncovered"] == 0 and CUTOFF == 0.5 else "STOP"),
    ]
    lines += [f"  gate 3 differs: {d}" for d in sorted(differ)]
    if any("STOP" in line for line in lines):
        raise Stop("\n".join(first["checks"] + lines))
    return {"cutoff": CUTOFF, "forms": list(FORMS), "five": list(FIVE), "clip": CLIP, "margin": MARGIN,
            "first_checks": first["checks"], "checks": lines,
            "datasets": sorted(rows, key=lambda r: (r["species"], r["antibiotic"]))}


# ---- gates 7 and 8 -------------------------------------------------------------------------------------------------

def long_lists(reading):
    """Lists in the reading longer than MAX_LIST, but its report list and gate lines: a row per genome is one."""
    count = 0

    def walk(v):
        nonlocal count
        if isinstance(v, dict):
            for x in v.values():
                walk(x)
        elif isinstance(v, list):
            count += len(v) > MAX_LIST
            for x in v:
                walk(x)

    for k, v in reading.items():
        if k == "datasets":
            for r in v:
                walk(r)
        elif k not in ("checks", "first_checks"):
            walk(v)
    return count


def ids_of(directory):
    """The reports' genome ids, held in memory for gate 8 only; never printed or written."""
    ids = set()
    for name in lo.report_files(directory):
        with open(os.path.join(directory, name), encoding="utf-8") as f:
            for s in json.load(f).get("samples") or []:
                if isinstance(s.get("id"), str) and s["id"]:
                    ids.add(s["id"])
    return ids


def id_scan(paths, ids):
    """(files holding a genome id as a plain substring, id-shaped tokens and `id` or `samples` keys): counts only."""
    files, shapes = 0, 0
    for path in paths:
        with open(path, encoding="utf-8") as f:
            text = f.read()
        files += any(i in text for i in ids)
        shapes += len(lo.ID_SHAPE.findall(text)) + len(lo.BARRED_KEY.findall(text))
    return files, shapes


# ---- the record ----------------------------------------------------------------------------------------------------

def cell(x):
    """A decimal of the record at 12 places. Its size must stay under DECIMAL_LIMIT, so it has one digit before its
    point and cannot hold a genome id."""
    if x is None or x == "":
        return ""
    if abs(x) >= DECIMAL_LIMIT:
        raise Stop("gate 8 ids: a decimal of the record is 10 or more in size: STOP")
    return f"{common.r12(x):.12f}"


def dataset_row(r):
    out = {k: r[k] for k in ("species", "antibiotic", "genomes", "resistant", "susceptible", "mean_miss_by_1to4",
                             "mean_miss_by_all5")}
    for f in FORMS:
        for k in COUNTS:
            out[f"{f}_{k}"] = r[f][k]
        out[f"{f}_log_loss"] = cell(r[f]["log_loss"])
    for a in ARMS:
        for t in TRANSITIONS:
            out[f"{a}_{t}"] = r[f"{a}_transitions"][t]
        out[f"{a}_folds_converged"] = sum(1 for x in r[f"{a}_folds"] if x["converged"])
    return out


def scopes(rows):
    return [("all", "", rows)] + [("species", sp, [r for r in rows if r["species"] == sp])
                                  for sp in sorted({r["species"] for r in rows})]


def pooled_rows(rows):
    out = []
    for scope, sp, kept in scopes(rows):
        for arm in ARMS:
            for name in MEASURES:
                p = pooled_row(kept, arm, name, scope, sp)
                out.append(dict(p, median=cell(p["median"]), lo=cell(p["lo"]), hi=cell(p["hi"])))
    return out


def rule_rows(rows):
    out = []
    for arm in ARMS:
        r = rule(rows, arm)
        out.append(dict(r, **{k: cell(r[k]) for k in ("median_gain", "gain_lo", "gain_hi", "vme_median", "vme_lo",
                                                      "vme_hi")}))
    return out


def transition_rows(rows):
    out = []
    for scope, sp, kept in scopes(rows):
        for arm in ARMS:
            t = {"scope": scope, "species": sp, "arm": arm, "datasets": len(kept)}
            for col, src in (("resistant", "resistant"), ("mean_miss", "mean_fn"),
                             ("mean_miss_by_1to4", "mean_miss_by_1to4"), ("mean_miss_by_all5", "mean_miss_by_all5"),
                             ("susceptible", "susceptible"), ("mean_wrong", "mean_fp")):
                t[col] = sum(int(r[src]) for r in kept)
            for k in TRANSITIONS:
                t[k] = sum(int(r[f"{arm}_{k}"]) for r in kept)
            out.append(t)
    return out


def arithmetic_holds(r):
    """An arm's resistant calls are the incumbent's, plus those it recovers, less those it loses; its susceptible
    genomes called resistant are the incumbent's, less those it clears, plus those it newly calls."""
    for a in ARMS:
        if int(r[f"{a}_tp"]) != (int(r["mean_tp"]) + int(r[f"{a}_recovered_1to4"]) + int(r[f"{a}_recovered_all5"])
                                 - int(r[f"{a}_lost"])):
            return False
        if int(r[f"{a}_fp"]) != int(r["mean_fp"]) - int(r[f"{a}_cleared"]) + int(r[f"{a}_new_wrong"]):
            return False
    return True


def staged(out):
    lo.fresh_out(out)
    return tempfile.mkdtemp(prefix=".combiner-", dir=os.path.dirname(os.path.abspath(out)))


def write_reading(reading, out, ids):
    lists = long_lists(reading)
    if lists:
        raise Stop(f"gate 7 rows: {lists} lists in the reading longer than {MAX_LIST}: STOP")
    tmp = staged(out)
    try:
        with open(os.path.join(tmp, READING), "w", encoding="utf-8") as f:
            json.dump(reading, f, indent=1, sort_keys=True)
            f.write("\n")
        files, shapes = id_scan([os.path.join(tmp, READING)], ids)
        if files or shapes:
            raise Stop(f"gate 8 ids: {files} files holding a genome id, {shapes} id-shaped tokens or barred keys in "
                       "the reading: STOP")
        lo.publish(tmp, out)
    finally:
        if os.path.isdir(tmp):
            shutil.rmtree(tmp)


def write_record(reading_dir, reports_dir, out):
    with open(os.path.join(reading_dir, READING), encoding="utf-8") as f:
        reading = json.load(f)
    if (reading.get("cutoff"), reading.get("forms"), reading.get("five"), reading.get("clip")) != (
            CUTOFF, list(FORMS), list(FIVE), CLIP):
        raise Stop("the reading was made under other rules")
    rows = [dataset_row(r) for r in reading["datasets"]]
    bad = sum(1 for r in rows if not arithmetic_holds(r))
    species = len({r["species"] for r in rows})
    pooled, rules, trans = pooled_rows(rows), rule_rows(rows), transition_rows(rows)
    lists = long_lists(reading)
    line7 = (f"gate 7 rows: {len(rows)} rows in datasets.tsv for {len(reading['datasets'])} reports, {len(pooled)} in "
             f"pooled.tsv for {1 + species} scopes, {lists} lists in the reading longer than {MAX_LIST}, {bad} reports "
             f"whose arm counts do not follow from the incumbent's and the transitions: "
             + ("pass" if len(rows) == len(reading["datasets"]) and len(pooled) == (1 + species) * len(ARMS) * len(
                 MEASURES) and lists == 0 and bad == 0 else "STOP"))
    if "STOP" in line7:
        raise Stop(line7)
    ids = ids_of(reports_dir)
    tmp = staged(out)
    try:
        lo.write_tsv(os.path.join(tmp, "datasets.tsv"), DATASET_COLUMNS, rows)
        lo.write_tsv(os.path.join(tmp, "pooled.tsv"), POOLED_COLUMNS, pooled)
        lo.write_tsv(os.path.join(tmp, "rule.tsv"), RULE_COLUMNS, rules)
        lo.write_tsv(os.path.join(tmp, "transitions.tsv"), TRANSITION_COLUMNS, trans)
        scanned = lo.files_under(reading_dir) + [os.path.join(tmp, n) for n in RECORD_FILES[:4]]
        files, shapes = id_scan(scanned, ids)
        line8 = (f"gate 8 ids: {len(scanned)} files scanned against {len(ids):,} genome ids, {files} holding one as "
                 f"a substring, {shapes} id-shaped tokens or barred keys: " + ("pass" if files == shapes == 0 else "STOP"))
        if "STOP" in line8:
            raise Stop(line8)
        with open(os.path.join(tmp, "gates.txt"), "w", encoding="utf-8") as f:
            f.write("\n".join(reading["first_checks"] + reading["checks"] + [line7, line8]) + "\n")
            f.write("".join(f"reading {r['arm']}: {r['reading']}\n" for r in rules))
        with open(os.path.join(tmp, "SHA256SUMS"), "w", encoding="utf-8") as f:
            for name in RECORD_FILES:
                f.write(f"{lo.sha256(os.path.join(tmp, name))}  {name}\n")
        if any(id_scan([os.path.join(tmp, "gates.txt"), os.path.join(tmp, "SHA256SUMS")], ids)):
            raise Stop("gate 8 ids: gates.txt or SHA256SUMS holds a genome id or an id-shaped token: STOP")
        lo.publish(tmp, out)
    finally:
        if os.path.isdir(tmp):
            shutil.rmtree(tmp)
    return rows, rules, line7, line8


# ---- the page's tables ---------------------------------------------------------------------------------------------

def dec(x, places):
    """A page figure: a record decimal rounded once, half to even, to the places the table shows."""
    if x is None or x == "":
        return "none"
    q = round(fractions.Fraction(x) * 10 ** places)
    sign = "-" if q < 0 else ""
    whole, part = divmod(abs(q), 10 ** places)
    return f"{sign}{whole:,}.{part:0{places}d}"


def tables(record):
    rules = lo.read_tsv(os.path.join(record, "rule.tsv"))
    pooled = lo.read_tsv(os.path.join(record, "pooled.tsv"))
    trans = lo.read_tsv(os.path.join(record, "transitions.tsv"))
    c = lo.count
    return {
        "rule": ["| Arm | Reports | Species | Median gain in balanced accuracy | Its interval | Median change in the "
                 "very major error rate | Its interval | Reports gaining 0.02 or more | Reading |",
                 "|---|---:|---:|---:|---|---:|---|---:|---|"]
        + [f"| {r['arm']} | {c(r['datasets'])} | {c(r['species'])} | {dec(r['median_gain'], 4)} "
           f"| {dec(r['gain_lo'], 4)} to {dec(r['gain_hi'], 4)} | {dec(r['vme_median'], 4)} "
           f"| {dec(r['vme_lo'], 4)} to {dec(r['vme_hi'], 4)} | {c(r['gain_at_least_margin'])} | {r['reading']} |"
           for r in rules],
        "measures": ["| Arm | Measure | Reports | Median difference | Its interval | Better | Worse | Same |",
                     "|---|---|---:|---:|---|---:|---:|---:|"]
        + [f"| {p['arm']} | {p['measure']} | {c(p['datasets'])} | {dec(p['median'], 4)} | {dec(p['lo'], 4)} to "
           f"{dec(p['hi'], 4)} | {c(p['better'])} | {c(p['worse'])} | {c(p['same'])} |"
           for p in pooled if p["scope"] == "all"],
        "transitions": ["| Arm | Resistant the mean misses, 1 to 4 of the five call | Of those, the arm calls | "
                        "Resistant all five miss | Of those, the arm calls | Resistant the arm loses | Susceptible the "
                        "mean calls resistant | Of those, the arm clears | Susceptible the arm newly calls resistant |",
                        "|---|---:|---:|---:|---:|---:|---:|---:|---:|"]
        + [f"| {t['arm']} | {c(t['mean_miss_by_1to4'])} | {c(t['recovered_1to4'])} | {c(t['mean_miss_by_all5'])} "
           f"| {c(t['recovered_all5'])} | {c(t['lost'])} | {c(t['mean_wrong'])} | {c(t['cleared'])} "
           f"| {c(t['new_wrong'])} |" for t in trans if t["scope"] == "all"],
        "species": ["| Species | Arm | Reports | Median gain in balanced accuracy | Better | Worse | Same |",
                    "|---|---|---:|---:|---:|---:|---:|"]
        + [f"| *{p['species']}* | {p['arm']} | {c(p['datasets'])} | {dec(p['median'], 4)} | {c(p['better'])} "
           f"| {c(p['worse'])} | {c(p['same'])} |"
           for p in pooled if p["scope"] == "species" and p["measure"] == "balanced_accuracy"],
    }


# ---- the command line ----------------------------------------------------------------------------------------------

def plan():
    return {"cutoff": CUTOFF, "incumbent": "math.fsum(five scores) >= 2.5", "forms": list(FORMS), "arms": list(ARMS),
            "five": list(FIVE), "fit": "per report and stored outer fold, on the genomes of the other folds, by plain "
            "mean log loss", "clip": CLIP, "margin": MARGIN, "resamples": common.RESAMPLES,
            "interval_values": [common.LO_INDEX + 1, common.HI_INDEX + 1], "min_species": common.MIN_SPECIES,
            "weighted_fit": {"cap": W_CAP, "tolerance": W_TOL, "weights_sum_within": SUM_TOL},
            "logistic_fit": {"cap": L_CAP, "gradient": L_TOL, "step": L_STEP, "pivot": PIVOT},
            "measures": list(MEASURES), "lower_better": list(LOWER_BETTER), "datasets_columns": DATASET_COLUMNS,
            "pooled_columns": POOLED_COLUMNS, "rule_columns": RULE_COLUMNS, "transition_columns": TRANSITION_COLUMNS}


def main(argv=None):
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("plan")
    r = sub.add_parser("read")
    r.add_argument("reports")
    r.add_argument("--first", required=True)
    r.add_argument("--second", required=True)
    r.add_argument("--out", required=True)
    r.add_argument("--expect", type=int, default=lo.EXPECT_FILES)
    w = sub.add_parser("write")
    w.add_argument("reading")
    w.add_argument("--reports", required=True)
    w.add_argument("--out", required=True)
    p = sub.add_parser("print")
    p.add_argument("record")
    args = ap.parse_args(argv)
    try:
        if args.cmd == "plan":
            print(json.dumps(plan(), indent=1))
        elif args.cmd == "read":
            lo.fresh_out(args.out)
            reading = read(args.reports, args.first, args.second, args.expect)
            write_reading(reading, args.out, ids_of(args.reports))
            print("\n".join(reading["first_checks"] + reading["checks"]))
            print(f"read: {len(reading['datasets'])} reports; wrote {READING}")
        elif args.cmd == "write":
            rows, rules, line7, line8 = write_record(args.reading, args.reports, args.out)
            print(line7)
            print(line8)
            print(f"wrote datasets.tsv ({len(rows)} rows), pooled.tsv, rule.tsv, transitions.tsv, gates.txt, SHA256SUMS")
        else:
            for name, lines in tables(args.record).items():
                print(f"## {name}\n" + "\n".join(lines) + "\n")
    except Stop as e:
        print(str(e))
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
