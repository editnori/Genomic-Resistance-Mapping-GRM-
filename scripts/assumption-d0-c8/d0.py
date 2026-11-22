"""D0 of the assumption plan: on a lineage fold where a learner misses and the
program does not, why?

    python3 -I d0.py dataset --report R.json --export DIR --calls calls.tsv \
        --reference callers.tsv --out DIR
    python3 -I d0.py combine --out DIR DATASET_DIR...

`dataset` reads one dataset's report (with `fold_models`), its `grm export`
directory (matrix.bin, dataset.json) and the program's per-genome calls, and
places each missed lineage fold of each learner in a cause, then the dataset.
`combine` adds the datasets up per learner. Nothing is fitted or drawn; the
screens are recomputed by the rule `learn::screen` states.

D0's sets are chosen by their result: its answer orders comparisons and tests
nothing. The output says so beside every count.

The causes, in the plan's order, which also settles a fold that meets two:

1. Absent. No matrix row agrees with the program's calls on the fold's
   training genomes at 0.95 or more, in either polarity.
2. Removed. Rows agree, but none is offered to the learner: its support on
   the training samples is outside [floor, n_train - floor], or (rf, gbt,
   lasso) the lowest row with its training pattern is past the 4,000 the
   learner's chi-square screen keeps.
3. Lineage. An offered matching row exists and the learner chose a feature
   that agrees with one at 0.95 or more on the training samples, and the
   genomes carrying that feature's resistant-indicating state lie in at most
   two lineage clusters.
4. Unused. An offered matching row exists and no chosen feature agrees with
   one at 0.95.

A fold whose chosen feature agrees with an offered match but spans more than
two clusters is outside the four: `chosen`. The report lists the choices
behind these readings (REPORT.md, "Open choices").
"""

import argparse
import csv
import hashlib
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import foldread as fr  # noqa: E402

CAUSES = ("absent", "removed", "lineage", "chosen", "unused")
FOUR = ("absent", "removed", "lineage", "unused")
# Read from a released report, which has no fold models.
OPEN = ("absent", "removed", "offered", "inside_floor")
OPEN_ONLY = ("offered", "inside_floor")
GAP = 0.10
EPS = 1e-9
SCREEN = 4000
CLUSTERS = 2
CHUNK = 1 << 16
READING = "sets chosen by their result; orders comparisons, tests nothing"
POP = np.array([bin(i).count("1") for i in range(256)], dtype=np.int64)
# Constants that mix each word into two independent 64-bit hashes of a
# row's training pattern, to find the lowest row of each pattern.
_rng = np.random.default_rng(20261121)
MIX = _rng.integers(1, 2**63, size=(2, 4096), dtype=np.uint64) | np.uint64(1)


class Gate(SystemExit):
    """An input disagrees with another; the dataset is not read."""


def gate(ok, message):
    if not ok:
        raise Gate(f"gate failed: {message}")


# ---------------------------------------------------------------- bits


def to_words(mask, n):
    """A list of booleans as little-endian 64-bit words."""
    words = np.zeros((n + 63) // 64, dtype=np.uint64)
    for i, v in enumerate(mask):
        if v:
            words[i >> 6] |= np.uint64(1) << np.uint64(i & 63)
    return words


def popcount(a):
    """Set bits per row of a 2-D uint64 array."""
    a = np.ascontiguousarray(a)
    return POP[a.view(np.uint8)].reshape(a.shape[0], -1).sum(axis=1)


def row_hash(masked):
    """Two 64-bit hashes per row; equal patterns hash alike."""
    w = masked.shape[1]
    out = np.zeros((masked.shape[0], 2), dtype=np.uint64)
    for k in range(2):
        x = masked ^ MIX[k, :w]
        x = (x ^ (x >> np.uint64(31))) * np.uint64(0x7FB5D329728EA185)
        x = (x ^ (x >> np.uint64(27))) * np.uint64(0x81DADEF4BC2DD44D)
        x = x ^ (x >> np.uint64(33))
        acc = np.zeros(masked.shape[0], dtype=np.uint64)
        for j in range(w):
            acc = (acc * np.uint64(0x9E3779B97F4A7C15)) ^ x[:, j]
        out[:, k] = acc
    return out


def chi_square(a, b, np_, nn):
    """`learn::screen`'s statistic, in its order of operations."""
    n = np_ + nn
    c = np_ - a
    d = nn - b
    denom = (a + b) * (c + d) * np_ * nn
    diff = a * d - b * c
    out = np.zeros(a.shape[0], dtype=np.float64)
    nz = denom != 0.0
    out[nz] = n * (diff[nz] * diff[nz]) / denom[nz]
    return out


# ---------------------------------------------------------------- inputs


def read_calls(path, dataset, variant):
    """Genome id to True (resistant), False (susceptible) or None (no call)."""
    out = {}
    with open(path, encoding="utf-8", newline="") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["dataset_id"] != dataset or r["variant"] != variant:
                continue
            gate(r["genome_id"] not in out, f"{dataset}: a genome has two {variant} calls")
            out[r["genome_id"]] = {"resistant": True, "susceptible": False, "no_call": None}[r["call"]]
    return out


def read_reference(path, dataset):
    """The repository's confusion counts per caller for one dataset."""
    out = {}
    with open(path, encoding="utf-8", newline="") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["dataset"] == dataset:
                out[r["caller"]] = {k: int(r[k]) for k in ("tp", "fn", "tn", "fp")}
    return out


def confusion(truth, calls):
    tp = sum(t and c for t, c in zip(truth, calls))
    fn = sum(t and not c for t, c in zip(truth, calls))
    tn = sum(not t and not c for t, c in zip(truth, calls))
    fp = sum(not t and c for t, c in zip(truth, calls))
    return {"tp": tp, "fn": fn, "tn": tn, "fp": fp}


def ratio(a, b):
    return a / b if b else 0.0


def ba(truth, calls):
    c = confusion(truth, calls)
    return (ratio(c["tp"], c["tp"] + c["fn"]) + ratio(c["tn"], c["tn"] + c["fp"])) / 2.0


class Export:
    """A `grm export` directory: the plan, the forest's screens, the matrix."""

    def __init__(self, path, check_sha=True):
        with open(os.path.join(path, "dataset.json"), encoding="utf-8") as f:
            self.meta = json.load(f)
        m = self.meta
        gate(m.get("format") == "grm-dataset-export-v1", f"{path}: not a grm-dataset-export-v1 file")
        size = m["rows"] * m["words"] * 8
        mpath = os.path.join(path, "matrix.bin")
        gate(os.path.getsize(mpath) == size, f"{path}: matrix.bin is not rows x words x 8 bytes")
        if check_sha:
            h = hashlib.sha256()
            with open(mpath, "rb") as f:
                for block in iter(lambda: f.read(1 << 24), b""):
                    h.update(block)
            gate(h.hexdigest() == m["matrix_sha256"], f"{path}: matrix.bin sha256 differs")
        self.x = np.memmap(mpath, dtype="<u8", mode="r", shape=(m["rows"], m["words"]))

    def split(self, scheme):
        for s in self.meta["splits"]:
            if s["plan"]["scheme"] == scheme:
                return s
        return None


# ---------------------------------------------------------------- one fold


def scan_fold(x, train, called, program, fits, y, floor):
    """One pass over the matrix for one outer fold.

    train: the fold's training samples (words). called: samples with a
    program call. program: samples the program calls resistant. fits: per
    learner with a screen, its fitting samples (words). y: resistant samples.
    Returns the candidates (lowest row per training pattern among rows inside
    the floor), the chi-square of each on each screening set, and the rows
    that match the program's calls, each with its candidate row (or -1).
    """
    tc = train & called
    n_tc = int(popcount(tc[None, :])[0])
    n_train = int(popcount(train[None, :])[0])
    lo, hi = floor, n_train - floor
    cand_rows, cand_hash = [], []
    match_rows, match_hash, match_inside = [], [], []
    for start in range(0, x.shape[0], CHUNK):
        block = np.asarray(x[start:start + CHUNK])
        masked = block & train
        support = popcount(masked)
        disagree = popcount((block ^ program) & tc)
        agree = n_tc - disagree
        best = np.maximum(agree, disagree)
        match = 20 * best >= 19 * n_tc if n_tc else np.zeros(block.shape[0], dtype=bool)
        inside = (support >= lo) & (support <= hi)
        h = row_hash(masked)
        idx = np.nonzero(inside)[0]
        cand_rows.append(idx + start)
        cand_hash.append(h[idx])
        midx = np.nonzero(match)[0]
        match_rows.append(midx + start)
        match_hash.append(h[midx])
        match_inside.append(inside[midx])
    cand_rows = np.concatenate(cand_rows)
    cand_hash = np.concatenate(cand_hash) if cand_rows.size else np.zeros((0, 2), dtype=np.uint64)
    # The lowest row of each training pattern: rows are in order, so the
    # first index np.unique returns is the lowest row.
    if cand_rows.size:
        _, first = np.unique(cand_hash, axis=0, return_index=True)
        first.sort()
        cand_rows = cand_rows[first]
        cand_hash = cand_hash[first]
    scores = {}
    for name, fit in fits.items():
        pos, neg = fit & y, fit & ~y
        np_ = float(popcount(pos[None, :])[0])
        nn = float(popcount(neg[None, :])[0])
        a = np.zeros(cand_rows.size, dtype=np.float64)
        b = np.zeros(cand_rows.size, dtype=np.float64)
        for s in range(0, cand_rows.size, CHUNK):
            block = np.asarray(x[cand_rows[s:s + CHUNK]])
            a[s:s + CHUNK] = popcount(block & pos)
            b[s:s + CHUNK] = popcount(block & neg)
        scores[name] = chi_square(a, b, np_, nn)
    match_rows = np.concatenate(match_rows)
    match_hash = np.concatenate(match_hash) if match_rows.size else np.zeros((0, 2), dtype=np.uint64)
    match_inside = np.concatenate(match_inside)
    # Each match inside the floor maps to its pattern's candidate row.
    lookup = {(int(h[0]), int(h[1])): int(r) for h, r in zip(cand_hash, cand_rows)}
    rep = np.array(
        [
            lookup[(int(h[0]), int(h[1]))] if ok else -1
            for h, ok in zip(match_hash, match_inside)
        ],
        dtype=np.int64,
    )
    return {
        "n_train": n_train,
        "n_called_train": n_tc,
        "candidates": cand_rows,
        "scores": scores,
        "match_rows": match_rows,
        "match_rep": rep,
    }


def top(scores, rows, k):
    """`learn::screen`'s top k: highest score first, equal scores by lower row."""
    order = np.lexsort((rows, -scores))
    return rows[order[:k]]


def screened_set(scan, name, k):
    return set(int(r) for r in top(scan["scores"][name], scan["candidates"], k))


def agrees(a, b, mask, n):
    """Whether pattern rows a (k x w) agree with pattern b on `mask` at 0.95
    or more, in either polarity."""
    d = popcount((a ^ b) & mask)
    return 20 * np.maximum(n - d, d) >= 19 * n


def clusters_of(x_row, program, train, called, lineage, n):
    """Lineage clusters among training genomes carrying the row's
    resistant-indicating state: presence when presence agrees with the
    program at least as often as absence does on called training genomes."""
    tc = train & called
    n_tc = int(popcount(tc[None, :])[0])
    d = int(popcount(((x_row ^ program) & tc)[None, :])[0])
    state = x_row if n_tc - d >= d else ~x_row
    carry = state & train
    seen = set()
    for i in range(n):
        if int(carry[i >> 6]) >> (i & 63) & 1:
            seen.add(("lineage", lineage[i]) if lineage[i] is not None else ("alone", i))
    return len(seen)


def fold_cause(x, scan, learner, chosen_rows, screen, program, train, called, lineage, n):
    """The cause of one learner's one fold, and the counts behind it."""
    rows = scan["match_rows"]
    rep = scan["match_rep"]
    out = {
        "matching_rows": int(rows.size),
        "inside_floor": int((rep >= 0).sum()),
        "offered": 0,
        "past_screen": 0,
        "chosen_agreeing": 0,
        "chosen_agreeing_in_2": 0,
    }
    if rows.size == 0:
        return "absent", out
    inside = rep[rep >= 0]
    if screen is None:
        offered = np.unique(inside)
    else:
        keep = np.array([int(r) in screen for r in inside], dtype=bool)
        offered = np.unique(inside[keep])
        out["past_screen"] = int(np.unique(inside[~keep]).size)
    out["offered"] = int(offered.size)
    if offered.size == 0:
        return "removed", out
    pats = np.asarray(x[offered]) & train
    n_train = scan["n_train"]
    agreeing = []
    for g in sorted(chosen_rows):
        row = np.asarray(x[g])
        if agrees(pats, row & train, train, n_train).any():
            agreeing.append(g)
    out["chosen_agreeing"] = len(agreeing)
    narrow = [g for g in agreeing if clusters_of(np.asarray(x[g]), program, train, called, lineage, n) <= CLUSTERS]
    out["chosen_agreeing_in_2"] = len(narrow)
    if narrow:
        return "lineage", out
    if agreeing:
        return "chosen", out
    return "unused", out


# ---------------------------------------------------------------- one dataset


def fold_open(scan, learner, screen):
    """Without fold models: how far the plan's causes can be read. A fold is
    `absent` or `removed` as with them; otherwise the matching row reaches
    the learner (`offered`: Lineage, chosen or Unused needs the fold model)
    or, for boosting and the lasso without the plan's holdouts, passes the
    floor with its screen not recomputed (`inside_floor`)."""
    rows, rep = scan["match_rows"], scan["match_rep"]
    out = {"matching_rows": int(rows.size), "inside_floor": int((rep >= 0).sum()), "offered": 0, "past_screen": 0}
    if rows.size == 0:
        return "absent", out
    inside = rep[rep >= 0]
    if inside.size == 0:
        return "removed", out
    if learner in ("gbt", "lasso") and screen is None:
        return "inside_floor", out
    if screen is None:
        offered = np.unique(inside)
    else:
        keep = np.array([int(r) in screen for r in inside], dtype=bool)
        offered = np.unique(inside[keep])
        out["past_screen"] = int(np.unique(inside[~keep]).size)
    out["offered"] = int(offered.size)
    return ("offered" if offered.size else "removed"), out


def common(report, calls, reference):
    """What every reading checks first: the calls cover the report's genomes
    and, with a reference, give its recorded counts."""
    rid = report["id"]
    ids = [s["id"] for s in report["samples"]]
    gate(sorted(calls) == sorted(ids), f"{rid}: the calls do not cover the report's genomes exactly")
    n = len(ids)
    floor = report["counts"]["min_support"]
    gate(floor == max(3, n // 200), f"{rid}: support floor differs from max(3, n / 200)")
    truth = [bool(s["resistant"]) for s in report["samples"]]
    keep = [i for i in range(n) if calls[ids[i]] is not None]
    truth_k = [truth[i] for i in keep]
    prog_all = [bool(calls[ids[i]]) for i in keep]
    if reference is not None:
        gate(reference.get("program:tool_default") == confusion(truth_k, prog_all),
             f"{rid}: program counts differ from the reference")
    return ids, n, floor, truth


def read_dataset(report, export, calls, reference=None, learners=None, gap=GAP, screen_size=None):
    """Every learner's fold rows and dataset row for one dataset, from its
    export and fold models. The screen keeps the report's
    `screened_features` rows; `screen_size`, when given, must equal it
    (4,000 in the release)."""
    rid = report["id"]
    ids, n, floor, truth = common(report, calls, reference)
    m = export.meta
    gate(m["samples"] == ids, f"{rid}: the export's samples are not the report's, in order")
    gate(m["cols"] == n and m["words"] == (n + 63) // 64, f"{rid}: matrix shape does not fit the samples")
    gate(m["input_fingerprint"] == report["input_fingerprint"], f"{rid}: input fingerprints differ")
    gate(m["evaluation_revision"] == report["evaluation"]["revision"], f"{rid}: evaluation revisions differ")
    gate(m["min_support"] == floor, f"{rid}: support floor differs from the export's")
    k = report["eval_config"]["screened_features"]
    gate(m["screened_features"] == k and (screen_size is None or k == screen_size), f"{rid}: screen size differs")
    gate(m["resistant"] == truth, f"{rid}: the export's labels differ from the report's")
    gate(m["lineage"] == [s["lineage"] for s in report["samples"]], f"{rid}: the export's lineages differ from the report's")
    x = export.x
    gate(x.shape == (m["rows"], m["words"]), f"{rid}: matrix shape")
    models = "fold_models" in report
    patterned = (report["fold_models"]["features"] if models else []) + report["features"]
    for f in patterned:
        row = int.from_bytes(np.asarray(x[int(f["feature"])]).astype("<u8").tobytes(), "little")
        gate(row == fr.decode_bits(f["bits"], n), f"{rid}: a fold-model or card pattern is not its matrix row")
    s = fr.split(report, "lineage")
    if s is None:
        return [], []
    es = export.split("lineage")
    gate(es is not None and not es["plan"].get("unavailable_reason"), f"{rid}: no lineage plan in the export")
    plan = es["plan"]
    gate(plan["outer"] == s["test_fold"], f"{rid}: the plan's outer folds are not the report's test_fold")
    inputs = {"x": x, "holdouts": plan["holdouts"], "screens": {sc["fold"]: sc for sc in es["screens"]},
              "floor": floor, "k": k, "models": models}
    return place(report, inputs, calls, learners, gap, reference)


def project(species_dir, ids, floor):
    """`Dataset::project` of the species matrix onto the report's genomes:
    species patterns read on those genomes, kept when their support lies in
    [floor, n - floor], one row per distinct projected pattern in the order
    of its first species pattern. Returns the rows and the species patterns
    kept."""
    with open(os.path.join(species_dir, "matrix.json"), encoding="utf-8") as f:
        meta = json.load(f)
    genomes = meta["genomes"]
    col = {g: i for i, g in enumerate(genomes)}
    gate(all(i in col for i in ids), "a report genome is not in the species matrix")
    ws = (len(genomes) + 63) // 64
    path = os.path.join(species_dir, "patterns.bin")
    gate(os.path.getsize(path) == meta["n_patterns"] * ws * 8, "patterns.bin size does not match matrix.json")
    pats = np.memmap(path, dtype="<u8", mode="r", shape=(meta["n_patterns"], ws))
    cols = np.array([col[i] for i in ids], dtype=np.int64)
    n = len(ids)
    w = (n + 63) // 64
    kept, chunks = 0, []
    for start in range(0, pats.shape[0], CHUNK):
        block = np.ascontiguousarray(pats[start:start + CHUNK])
        bits = np.unpackbits(block.view(np.uint8), axis=1, bitorder="little")[:, cols]
        support = bits.sum(axis=1, dtype=np.int64)
        inside = (support >= floor) & (support <= n - floor)
        kept += int(inside.sum())
        packed = np.packbits(bits[inside], axis=1, bitorder="little")
        pad = np.zeros((packed.shape[0], w * 8), dtype=np.uint8)
        pad[:, :packed.shape[1]] = packed
        chunks.append(pad.view("<u8"))
    rows = np.concatenate(chunks) if chunks else np.zeros((0, w), dtype="<u8")
    _, first = np.unique(row_hash(rows), axis=0, return_index=True)
    first.sort()
    return np.ascontiguousarray(rows[first]), kept


def read_release(report, species_dir, calls, reference=None, learners=None, gap=GAP, screen_size=None):
    """As `read_dataset` for a released report, which has no fold models and
    no export: the dataset matrix is projected from the species matrix and
    gated on the report's own counts and feature cards; the training samples
    are the report's `test_fold`. Boosting's and the lasso's screens need the
    plan's holdouts, so they are not recomputed."""
    rid = report["id"]
    ids, n, floor, truth = common(report, calls, reference)
    k = report["eval_config"]["screened_features"]
    gate(screen_size is None or k == screen_size, f"{rid}: screen size differs")
    x, kept = project(species_dir, ids, floor)
    gate(x.shape[0] == report["counts"]["features"], f"{rid}: {x.shape[0]} projected rows, the report has "
         f"{report['counts']['features']} features")
    for c in report["features"]:
        row = int.from_bytes(x[int(c["feature"])].astype("<u8").tobytes(), "little")
        gate(row == fr.decode_bits(c["bits"], n), f"{rid}: a feature card's pattern is not its projected row")
    s = fr.split(report, "lineage")
    if s is None:
        return [], []
    inputs = {"x": x, "holdouts": None, "screens": None, "floor": floor, "k": k, "models": False,
              "kept_patterns": kept}
    return place(report, inputs, calls, learners, gap, reference)


def place(report, inp, calls, learners, gap, reference):
    rid = report["id"]
    ids = [s["id"] for s in report["samples"]]
    n = len(ids)
    truth = [bool(s["resistant"]) for s in report["samples"]]
    lineage = [s["lineage"] for s in report["samples"]]
    keep = [i for i in range(n) if calls[ids[i]] is not None]
    truth_k = [truth[i] for i in keep]
    prog_ba = ba(truth_k, [bool(calls[ids[i]]) for i in keep])
    x, k, floor = inp["x"], inp["k"], inp["floor"]
    s = fr.split(report, "lineage")
    folds = s["folds"]
    called = to_words([calls[i] is not None for i in ids], n)
    program = to_words([bool(calls[i]) for i in ids], n)
    y = to_words(truth, n)
    rows_by_learner = {r["learner"]: r for r in s["learners"]}
    want = learners if learners is not None else fr.LEARNERS
    causes = CAUSES if inp["models"] else OPEN
    fold_out, data_out = [], []
    # The pooled gap of every learner, and its missed folds.
    plan_folds = {}
    for learner in want:
        lr = rows_by_learner[learner]
        lcalls = [lr["scores"][i] >= 0.5 for i in keep]
        if reference is not None:
            gate(reference.get(f"learner:{learner}:lineage") == confusion(truth_k, lcalls),
                 f"{rid}: {learner} lineage counts differ from the reference")
        lba = ba(truth_k, lcalls)
        dgap = prog_ba - lba
        folds_rows = []
        for f in range(folds):
            idx = [i for i in keep if s["test_fold"][i] == f]
            t = [truth[i] for i in idx]
            both = any(t) and not all(t)
            pf = ba(t, [bool(calls[ids[i]]) for i in idx])
            lf = ba(t, [lr["scores"][i] >= 0.5 for i in idx])
            if len(keep) == n:
                gate(abs(lf - lr["fold_balanced_accuracy"][f]) <= EPS, f"{rid}: {learner} fold {f} BA does not recompute")
            folds_rows.append((f, both, pf, lf, pf - lf))
        missed = [r[0] for r in folds_rows if r[1] and r[4] >= gap - EPS]
        fallback = False
        if dgap >= gap - EPS and not missed:
            scored = [r for r in folds_rows if r[1]]
            if scored:
                best = max(scored, key=lambda r: (r[4], -r[0]))
                missed = [best[0]]
                fallback = True
        plan_folds[learner] = (lba, dgap, folds_rows, missed, fallback)
    need = sorted({f for v in plan_folds.values() if v[1] >= gap - EPS for f in v[3]})
    scans = {}
    for f in need:
        train = to_words([t is not None and t != f for t in s["test_fold"]], n)
        fits = {"rf": train}
        if inp["holdouts"] is not None:
            hold = inp["holdouts"][f]
            for name in ("gbt", "lasso"):
                held = set(hold[name])
                ho = to_words([i in held for i in range(n)], n)
                gate(not (ho & ~train).any(), f"{rid}: a {name} holdout sample is not training in fold {f}")
                fits[name] = train & ~ho
        scan = scan_fold(x, train, called, program, fits, y, floor)
        if inp["screens"] is not None:
            sc = inp["screens"][f]
            gate(sc["n_train"] == scan["n_train"], f"{rid}: fold {f} training count differs from the export")
            gate(sc["n_candidates"] == scan["candidates"].size,
                 f"{rid}: fold {f} candidate count differs from the export")
            gate([int(r) for r in top(scan["scores"]["rf"], scan["candidates"], k)] == sc["screened_rows"],
                 f"{rid}: fold {f} forest screen differs from the export's")
        scan["train"] = train
        scans[f] = scan
    for learner in want:
        lba, dgap, folds_rows, missed, fallback = plan_folds[learner]
        in_miss = dgap >= gap - EPS
        sets = fr.fold_sets(report, "lineage", learner) if inp["models"] else None
        tally = dict.fromkeys(causes, 0)
        for f, both, pf, lf, fg in folds_rows:
            row = {"dataset": rid, "learner": learner, "fold": f, "in_miss": in_miss,
                   "both_classes": both, "program_fold_ba": pf, "learner_fold_ba": lf, "fold_gap": fg,
                   "missed": in_miss and f in missed, "fallback": fallback and in_miss and f in missed}
            if row["missed"]:
                scan = scans[f]
                screen = None
                if learner == "rf" or (learner in ("gbt", "lasso") and learner in scan["scores"]):
                    screen = screened_set(scan, learner, k)
                if inp["models"]:
                    cause, counts = fold_cause(x, scan, learner, sets[f], screen, program, scan["train"],
                                               called, lineage, n)
                else:
                    cause, counts = fold_open(scan, learner, screen)
                tc = scan["train"] & called
                n_tc = scan["n_called_train"]
                n_res = int(popcount((tc & program)[None, :])[0])
                row.update(counts)
                row["n_train"] = scan["n_train"]
                row["program_one_class"] = 20 * max(n_res, n_tc - n_res) >= 19 * n_tc
                row["cause"] = cause
                tally[cause] += 1
            fold_out.append(row)
        if in_miss:
            most = max(tally[c] for c in causes)
            top_causes = [c for c in causes if tally[c] == most and most > 0]
        else:
            top_causes = []
        data_out.append({"dataset": rid, "species": report["species_slug"], "learner": learner,
                         "program_ba": prog_ba, "learner_lineage_ba": lba, "gap": dgap, "in_miss": in_miss,
                         "screen": k, "source": "fold models" if inp["models"] else "release report, no fold models",
                         "missed_folds": len(missed) if in_miss else 0, "fallback": fallback if in_miss else False,
                         **{f"folds_{c}": tally.get(c) for c in CAUSES + OPEN_ONLY},
                         "cause": "+".join(top_causes) if top_causes else None})
    return fold_out, data_out


FOLD_HEADER = ["dataset", "learner", "fold", "in_miss", "both_classes", "program_fold_ba", "learner_fold_ba",
               "fold_gap", "missed", "fallback", "n_train", "program_one_class", "matching_rows", "inside_floor",
               "offered", "past_screen", "chosen_agreeing", "chosen_agreeing_in_2", "cause", "reading"]
DATA_HEADER = ["dataset", "species", "learner", "source", "program_ba", "learner_lineage_ba", "gap", "in_miss",
               "screen", "missed_folds", "fallback"] + [f"folds_{c}" for c in CAUSES + OPEN_ONLY] + \
    ["cause", "reading"]


def write_dataset(out, fold_rows, data_rows):
    os.makedirs(out, exist_ok=True)
    fr.write_tsv(os.path.join(out, "d0-folds.tsv"), FOLD_HEADER,
                 [[r.get(h) if h != "reading" else READING for h in FOLD_HEADER] for r in fold_rows])
    fr.write_tsv(os.path.join(out, "d0-datasets.tsv"), DATA_HEADER,
                 [[r.get(h) if h != "reading" else READING for h in DATA_HEADER] for r in data_rows])


# ---------------------------------------------------------------- combine


SUMMARY_HEADER = ["subset", "learner", "miss_datasets"] + [f"datasets_{c}" for c in CAUSES + OPEN_ONLY] + \
    [f"gap_share_{c}" for c in CAUSES + OPEN_ONLY] + ["ties", "reading"]


def combine(rows, subset="all"):
    """Per learner and over the five: datasets per cause (a tie counts in
    each tied cause) and each cause's share of the gap (a tie splits it).
    Also returns the five learners' counts and their miss pairs. `subset`
    names the rows in the table's first column."""
    out = []
    every, pairs = dict.fromkeys(CAUSES + OPEN_ONLY, 0), 0
    for learner in list(fr.LEARNERS) + ["all"]:
        mine = [r for r in rows if r["in_miss"] and (learner == "all" or r["learner"] == learner)]
        counts = dict.fromkeys(CAUSES + OPEN_ONLY, 0)
        share = dict.fromkeys(CAUSES + OPEN_ONLY, 0.0)
        ties = 0
        total = sum(r["gap"] for r in mine)
        for r in mine:
            causes = r["cause"].split("+") if r["cause"] else []
            ties += len(causes) > 1
            for c in causes:
                counts[c] += 1
                share[c] += r["gap"] / len(causes)
        if learner == "all":
            every, pairs = counts, len(mine)
        out.append([subset, learner, len(mine)] + [counts[c] for c in CAUSES + OPEN_ONLY] +
                   [share[c] / total if total else None for c in CAUSES + OPEN_ONLY] + [ties, READING])
    return out, every, pairs


def ordering(every, pairs):
    """The plan's ordering rule: the cause with the most datasets over the
    five learners sets the order when it holds 0.40 or more of the miss
    pairs (dataset and learner). A tie for most sets no order, and
    `chosen`, outside the four, never does."""
    if not pairs:
        return "no dataset in the miss set"
    still = sum(every.get(c, 0) for c in OPEN_ONLY)
    if still:
        known = ", ".join(f"{c} {every[c]}" for c in ("absent", "removed"))
        return (f"not settled: of {pairs} miss pairs, {known}; {every.get('offered', 0)} offered need the fold "
                f"models and {every.get('inside_floor', 0)} inside the floor need the plan's holdouts")
    most = max(every[c] for c in FOUR)
    tied = [c for c in FOUR if every[c] == most]
    share = most / pairs
    if most and len(tied) == 1 and share >= 0.40:
        best = tied[0]
        order = {"absent": "C3 first", "removed": "C4 floor arm first", "lineage": "C4 adjusted arms first",
                 "unused": "C5 first"}[best]
        return f"{best} holds {most} of {pairs} miss pairs ({share:.3f}): {order}"
    return f"no single cause holds 0.40 of {pairs} miss pairs: C4, C3, C5"


def complete(rows, expect_datasets=None, expect_pairs=None):
    """Gates on the rows combine reads: one row per dataset and learner, and
    the miss set's dataset and pair counts when given (the plan's 75 and
    204), so a species directory left out, or given twice, stops the run."""
    keys = [(r["dataset"], r["learner"]) for r in rows]
    gate(len(keys) == len(set(keys)), f"{len(keys) - len(set(keys))} dataset and learner rows come twice")
    miss = [r for r in rows if r["in_miss"]]
    if expect_datasets is not None:
        got = len({r["dataset"] for r in miss})
        gate(got == expect_datasets, f"{got} datasets in the miss set, expected {expect_datasets}")
    if expect_pairs is not None:
        gate(len(miss) == expect_pairs, f"{len(miss)} miss pairs, expected {expect_pairs}")


def read_rows(path):
    with open(path, encoding="utf-8", newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    for r in rows:
        r["in_miss"] = r["in_miss"] == "yes"
        r["gap"] = float(r["gap"])
        r["cause"] = None if r["cause"] == "NA" else r["cause"]
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="mode", required=True)
    d = sub.add_parser("dataset")
    d.add_argument("--report", required=True)
    d.add_argument("--export", required=True, help="the dataset's grm export directory")
    d.add_argument("--calls", required=True, help="the program adapter's calls.tsv")
    d.add_argument("--variant", default="tool_default")
    d.add_argument("--reference", help="reference/amrfinderplus-comparison/callers.tsv, for the count gates")
    d.add_argument("--screen", type=int, default=SCREEN, help="the screen size the report must hold (4000)")
    d.add_argument("--out", required=True)
    r = sub.add_parser("release", help="a released report with no fold models, read on the species matrix")
    r.add_argument("--report", required=True)
    r.add_argument("--species-matrix", required=True, help="the species' grm matrix directory")
    r.add_argument("--calls", required=True)
    r.add_argument("--variant", default="tool_default")
    r.add_argument("--reference")
    r.add_argument("--screen", type=int, default=SCREEN)
    r.add_argument("--out", required=True)
    c = sub.add_parser("combine")
    c.add_argument("--out", required=True)
    c.add_argument("--marks", help="intcount's cut table (datasets.tsv): adds the rows for M and not M")
    c.add_argument("--expect-datasets", type=int, help="datasets the miss set must hold (75)")
    c.add_argument("--expect-pairs", type=int, help="miss pairs, dataset and learner, it must hold (204)")
    c.add_argument("dirs", nargs="+", help="directories holding a d0-datasets.tsv")
    a = ap.parse_args(argv)
    if a.mode in ("dataset", "release"):
        report = fr.load_report(a.report, need_models=False)
        calls = read_calls(a.calls, report["id"], a.variant)
        ref = read_reference(a.reference, report["id"]) if a.reference else None
        if a.mode == "dataset":
            folds, data = read_dataset(report, Export(a.export), calls, ref, screen_size=a.screen)
        else:
            folds, data = read_release(report, a.species_matrix, calls, ref, screen_size=a.screen)
        write_dataset(a.out, folds, data)
        placed = [r for r in data if r["in_miss"]]
        print(f"d0: {report['id']}: {len(placed)} learners in the miss set; " +
              ", ".join(f"{r['learner']} {r['cause']}" for r in placed) + f" ({READING})")
    else:
        rows = []
        for p in a.dirs:
            rows += read_rows(os.path.join(p, "d0-datasets.tsv"))
        complete(rows, a.expect_datasets, a.expect_pairs)
        table, every, pairs = combine(rows)
        if a.marks:
            marks = fr.read_marks(a.marks)
            gone = sorted({r["dataset"] for r in rows} - set(marks))
            gate(not gone, f"{len(gone)} datasets are not in the marks table, the first {gone[0] if gone else ''}")
            for name, flag in (("not_M", False), ("M", True)):
                table += combine([r for r in rows if marks[r["dataset"]] is flag], name)[0]
        os.makedirs(a.out, exist_ok=True)
        fr.write_tsv(os.path.join(a.out, "d0-summary.tsv"), SUMMARY_HEADER, table)
        line = ordering(every, pairs)
        with open(os.path.join(a.out, "d0-ordering.txt"), "w", encoding="utf-8") as f:
            f.write(line + f" ({READING})\n")
        print(line)


if __name__ == "__main__":
    main()
