"""Synthetic reports, exports and program calls for the D0 and C8 tests.

60 samples in 10 lineage clusters of 6 (sample i in cluster i // 6); cluster
c is in lineage fold c % 5, so each fold tests 12 samples and trains on 48.
The support floor is 3. The forest's screen in the export is worked out here
in plain Python (candidates, chi-square, ties to the lower row), apart from
the numpy code under test.
"""

import base64
import hashlib
import json
import os

N = 60
FOLDS = 5
LEARNERS = ("scm", "cart", "rf", "gbt", "lasso")


def cluster(i):
    return i // 6


def fold_of(i):
    return cluster(i) % FOLDS


TEST_FOLD = [fold_of(i) for i in range(N)]


def training(f):
    return [i for i in range(N) if TEST_FOLD[i] != f]


def flip(v, idx):
    out = list(v)
    for i in idx:
        out[i] = 1 - out[i]
    return out


def bits_text(pattern):
    value = sum(1 << i for i, v in enumerate(pattern) if v)
    words = (len(pattern) + 63) // 64
    return base64.b64encode(value.to_bytes(words * 8, "little")).decode()


def ba(truth, calls):
    tp = sum(t and c for t, c in zip(truth, calls))
    fn = sum(t and not c for t, c in zip(truth, calls))
    tn = sum(not t and not c for t, c in zip(truth, calls))
    fp = sum(not t and c for t, c in zip(truth, calls))
    sens = tp / (tp + fn) if tp + fn else 0.0
    spec = tn / (tn + fp) if tn + fp else 0.0
    return (sens + spec) / 2.0


def candidates(rows, train, floor):
    """`learn::training_candidates`: rows inside the floor, lowest row per
    training pattern."""
    seen, out = set(), []
    for r, p in enumerate(rows):
        key = tuple(p[i] for i in train)
        s = sum(key)
        if s < floor or s > len(train) - floor or key in seen:
            continue
        seen.add(key)
        out.append(r)
    return out


def chi(p, y, fit):
    a = float(sum(1 for i in fit if p[i] and y[i]))
    b = float(sum(1 for i in fit if p[i] and not y[i]))
    np_ = float(sum(1 for i in fit if y[i]))
    nn = float(len(fit)) - np_
    c, d = np_ - a, nn - b
    denom = (a + b) * (c + d) * np_ * nn
    return 0.0 if denom == 0.0 else (np_ + nn) * (a * d - b * c) ** 2 / denom


def screen(rows, y, train, fit, floor, k):
    cand = candidates(rows, train, floor)
    return sorted(cand, key=lambda r: (-chi(rows[r], y, fit), r))[:k]


def scores_for(truth, wrong_folds=(), wrong_samples=()):
    """Scores that call every sample right, except whole folds or samples."""
    out = []
    for i, t in enumerate(truth):
        bad = TEST_FOLD[i] in wrong_folds or i in wrong_samples
        out.append(float(t != bad))
    return out


def fold_model(learner, chosen):
    if chosen is None:
        if learner != "lasso":
            raise ValueError(f"only a lasso fold model can be null, not {learner}")
        return None
    rows = sorted(chosen)
    if learner == "scm":
        return {"kind": "conjunction", "p": 1.0, "rules": [{"feature": r, "presence": True} for r in rows]}
    if learner == "cart":
        nodes = [{"value": 0.5, "n": 48, "n_resistant": 24, "depth": 0, "split": {"feature": r, "gain": 0.1,
                  "absent": 1, "present": 2}} for r in rows]
        return {"nodes": nodes + [{"value": 0.0, "n": 1, "n_resistant": 0, "depth": 1}]}
    if learner in ("rf", "gbt"):
        return [{"feature": r, "weight": 1.0 / (j + 1)} for j, r in enumerate(rows)]
    if learner == "lasso":
        return {"intercept": 0.0, "coef": [[r, 1.0] for r in rows], "lambda": 0.1}
    raise ValueError(learner)


def holdouts(f):
    """Each fold's tuning holdouts: the lowest 2 (boosting) and 3 (lasso)
    training samples."""
    t = training(f)
    return {"gbt": t[:2], "lasso": t[:3]}


def make(directory, truth, program, rows, scores=None, chosen=None, cards=(), k=3,
         dataset="synthia_example__drugamycin", no_call=(), release=False):
    """Write report.json, export/{matrix.bin,dataset.json}, calls.tsv and
    reference.tsv under `directory`; return their paths.

    truth, program: 0/1 per sample. rows: matrix rows, 0/1 per sample.
    scores: learner to its 60 scores (default: right everywhere).
    chosen: learner to 5 sets of rows (default: empty sets).
    cards: (row, [families]) pairs for feature cards.
    """
    os.makedirs(directory, exist_ok=True)
    ids = [f"g{i:03d}" for i in range(N)]
    scores = dict(scores or {})
    chosen = dict(chosen or {})
    lineage = [cluster(i) for i in range(N)]
    learners = []
    for l in LEARNERS:
        sc = scores.get(l, scores_for(truth))
        fold_ba = []
        for f in range(FOLDS):
            idx = [i for i in range(N) if TEST_FOLD[i] == f]
            fold_ba.append(ba([truth[i] for i in idx], [sc[i] >= 0.5 for i in idx]))
        learners.append({"learner": l, "scores": sc, "fold_balanced_accuracy": fold_ba})
    fm = {"scheme": "lineage"}
    pattern_rows = set()
    for l in LEARNERS:
        sets = chosen.get(l, [set()] * FOLDS)
        fm[l] = [fold_model(l, s) for s in sets]
        if l in ("scm", "cart", "lasso"):
            for s in sets:
                pattern_rows |= set(s or ())
    report = {
        "schema": 2,
        "evaluation": {"revision": 14},
        "input_fingerprint": "00ff00ff00ff00ff",
        "eval_config": {"screened_features": k},
        "id": dataset,
        "species_slug": dataset.split("__")[0],
        "counts": {"min_support": 3, "samples": N},
        "samples": [{"id": ids[i], "resistant": bool(truth[i]), "lineage": lineage[i]} for i in range(N)],
        "splits": [{"scheme": "lineage", "folds": FOLDS, "n_groups": 10, "test_fold": TEST_FOLD,
                    "learners": learners}],
        "features": [{"feature": r, "bits": bits_text(rows[r]), "catalog_families": fams} for r, fams in cards],
        "fold_models": {"splits": [fm],
                        "features": [{"feature": r, "bits": bits_text(rows[r])} for r in sorted(pattern_rows)]},
    }
    if release:
        # A released report: no fold models, and its count of dataset rows.
        del report["fold_models"]
        report["counts"]["features"] = len(rows)
    rpath = os.path.join(directory, "report.json")
    with open(rpath, "w", encoding="utf-8") as f:
        json.dump(report, f)
    words = (N + 63) // 64
    raw = b"".join(sum(1 << i for i, v in enumerate(p) if v).to_bytes(words * 8, "little") for p in rows)
    epath = os.path.join(directory, "export")
    os.makedirs(epath, exist_ok=True)
    with open(os.path.join(epath, "matrix.bin"), "wb") as f:
        f.write(raw)
    screens = []
    plan_holdouts = []
    for fo in range(FOLDS):
        t = training(fo)
        plan_holdouts.append(holdouts(fo))
        screens.append({"fold": fo, "n_train": len(t), "n_candidates": len(candidates(rows, t, 3)),
                        "screened_rows": screen(rows, truth, t, t, 3, k)})
    meta = {
        "format": "grm-dataset-export-v1", "evaluation_revision": 14, "id": dataset,
        "species_slug": dataset.split("__")[0], "antibiotic": dataset.split("__")[1],
        "input_fingerprint": "00ff00ff00ff00ff", "screened_features": k, "folds": FOLDS, "inner_folds": 3,
        "seed": 42, "min_support": 3, "samples": ids, "resistant": [bool(t) for t in truth],
        "lineage": lineage, "duplicates": [], "rows": len(rows), "cols": N, "words": words,
        "matrix_sha256": hashlib.sha256(raw).hexdigest(),
        "splits": [{"plan": {"scheme": "lineage", "outer": TEST_FOLD, "holdouts": plan_holdouts},
                    "screens": screens}],
    }
    with open(os.path.join(epath, "dataset.json"), "w", encoding="utf-8") as f:
        json.dump(meta, f)
    cpath = os.path.join(directory, "calls.tsv")
    with open(cpath, "w", encoding="utf-8") as f:
        f.write("dataset_id\tgenome_id\tvariant\tcall\tpoint_screened\telements_used\n")
        for i in range(N):
            call = "no_call" if i in no_call else ("resistant" if program[i] else "susceptible")
            f.write(f"{dataset}\t{ids[i]}\ttool_default\t{call}\tyes\t\n")
            f.write(f"{dataset}\t{ids[i]}\tgenes_core\tsusceptible\tno\t\n")
    keep = [i for i in range(N) if i not in no_call]
    refpath = os.path.join(directory, "reference.tsv")
    with open(refpath, "w", encoding="utf-8") as f:
        f.write("dataset\tcaller\ttp\tfn\ttn\tfp\n")
        callers = [("program:tool_default", program)] + [
            (f"learner:{l['learner']}:lineage", [int(s >= 0.5) for s in l["scores"]]) for l in learners]
        for name, calls in callers:
            t = [truth[i] for i in keep]
            c = [calls[i] for i in keep]
            tp = sum(a and b for a, b in zip(t, c))
            fn = sum(a and not b for a, b in zip(t, c))
            tn = sum(not a and not b for a, b in zip(t, c))
            fp = sum(not a and b for a, b in zip(t, c))
            f.write(f"{dataset}\t{name}\t{tp}\t{fn}\t{tn}\t{fp}\n")
    return {"report": rpath, "export": epath, "calls": cpath, "reference": refpath}


def noise(count, seed):
    """Rows of a fixed pseudo-random draw, each present in about half the
    samples; none comes near any truth vector these tests use."""
    state = seed
    out = []
    for _ in range(count):
        row = []
        for _ in range(N):
            state = (state * 6364136223846793005 + 1442695040888963407) % (1 << 64)
            row.append((state >> 33) & 1)
        out.append(row)
    return out


def make_species(directory, rows, extra=5, seed=21):
    """A species matrix whose projection onto the 60 report genomes gives
    `rows` back, in order. Its 65 genomes put 5 others first and the report's
    in reverse; each row comes twice, the copy differing only on the other
    genomes, and a row present in 2 report genomes, under the floor, follows."""
    ids = [f"g{i:03d}" for i in range(N)]
    genomes = [f"x{j}" for j in range(extra)] + ids[::-1]
    col = {g: c for c, g in enumerate(genomes)}
    other = noise(2 * len(rows), seed)
    pats = []
    for r, row in enumerate(rows):
        for dup in range(2):
            p = [0] * len(genomes)
            for i in range(N):
                p[col[ids[i]]] = row[i]
            for j in range(extra):
                p[j] = other[2 * r + dup][j]
            pats.append(p)
        q = [0] * len(genomes)
        q[col[ids[r % N]]] = q[col[ids[(r + 1) % N]]] = q[0] = 1
        pats.append(q)
    words = (len(genomes) + 63) // 64
    os.makedirs(directory, exist_ok=True)
    with open(os.path.join(directory, "patterns.bin"), "wb") as f:
        for p in pats:
            f.write(sum(1 << c for c, v in enumerate(p) if v).to_bytes(words * 8, "little"))
    with open(os.path.join(directory, "matrix.json"), "w", encoding="utf-8") as f:
        json.dump({"k": 31, "genomes": genomes, "min_support": 2, "n_patterns": len(pats)}, f)
    return directory
