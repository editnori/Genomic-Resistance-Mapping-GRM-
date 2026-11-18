"""SYNTHETIC inputs for scripts/study-split/record.py: a reduced input set in the shape of the real one.

Nothing here is data. Three made-up species hold four datasets that enter:

  alphabacter_primus__onemycin     40 genomes, BioProjects PRJSYN01 to 08, 5 genomes each
  alphabacter_primus__twomycin     60 genomes: 20 shared with onemycin (PRJSYN05 to 08),
                                   30 in PRJSYN09 to 14 and 10 with no BioProject
  betabacter_secundus__onemycin    40 genomes, PRJSYN21 to 28; its first two genomes are one copy group
  gammabacter_tertius__onemycin    40 genomes, PRJSYN31 to 38; no lineage split

and four that do not: one under 5 BioProjects, one under the share and over
the group, one over the group, and deltabacter_quartus__onemycin, which
passes the counts and fails its fold check.

Each BioProject holds both classes. The first BioProject of every dataset
holds 2 resistant genomes, so the cells can be worked by hand:

  study plans      every learner is wrong on the genomes of the dataset's first BioProject
  random, matched  scm is wrong on the dataset's first genome (resistant); the others are right
  lineage split    every learner is wrong on the dataset's first two genomes

scm scores 0 or 1. cart, rf, gbt and lasso score 0.9, 0.8, 0.7 and 0.6 for
the class they call and one minus that for the other.

make(root, fault=None) writes the set under root and returns root. A fault
name breaks one input the way a test of a stop needs it (FAULTS).
"""

import hashlib
import itertools
import json
import os
import random
from fractions import Fraction

LEARNERS = ("scm", "cart", "rf", "gbt", "lasso")
SURE = {"cart": 0.9, "rf": 0.8, "gbt": 0.7, "lasso": 0.6}
DRAWS = range(5)
KINDS = ("study", "random", "matched")
DAY = "2026-01-02"
FINGERPRINT = "5y57e71c00000000"
META = "0" * 63 + "1"
RULE = {"min_studies": 5, "min_study_share": 0.8, "max_group_share": 0.5}
FIXED = "11:05:00Z"

# Resistant genomes in each BioProject of 5, by BioProject number; the first of each dataset holds 2.
RESISTANT = {1: 2, 2: 1, 3: 4, 4: 3, 5: 2, 6: 1, 7: 4, 8: 2, 9: 3, 10: 1, 11: 4, 12: 2, 13: 3, 14: 1,
             21: 2, 22: 4, 23: 1, 24: 3, 25: 2, 26: 4, 27: 1, 28: 3,
             31: 2, 32: 3, 33: 1, 34: 4, 35: 1, 36: 4, 37: 2, 38: 3}

FAULTS = {
    "fold_sha": "fold file sha256",
    "study_split": "BioProject split",
    "evaluator": "counts differ from the evaluator's",
    "copy_group": "copy group split",
    "matched": "matched draw 1 counts differ from study draw 1",
    "no_class": "lacks a class",
    "revision": "release revision",
    "plan_hash": "PLAN.md differs from its hash list",
    "run_log": "fit lines",
    "rule_time": "the stated time of the rule choice",
    "reasons": "reasons differ from the counts",
    "variant": "variant",
    "draw0": "draw 0 is not the first experiment's fold files",
    "lineage": "lineage cluster split",
    "labels": "labels differ from the release",
    "rule_line": "the rule's line reads ['0.20']",
    "genome_id": "genome ids of the inputs stand in the record: 1, 9001.003 first",
    # Inputs the script does not expect.
    "extra_input": "first/notes.txt: an input file the script does not expect",
    "extra_fit": "fits/alphabacter_primus__onemycin.study-d5.tsv: a fit file the script does not expect",
    "extra_column": "column 'weight', a column the script does not expect",
    "extra_dataset": "release/genomes/epsilonbacter_quintus__onemycin.tsv: a dataset table the script does not expect",
}


def sha(b):
    return hashlib.sha256(b).hexdigest()


def write(root, rel, data):
    p = os.path.join(root, rel)
    os.makedirs(os.path.dirname(p), exist_ok=True)
    with open(p, "wb") as f:
        f.write(data if isinstance(data, bytes) else data.encode())
    return p


def tsv(head, rows):
    return "\t".join(head) + "\n" + "".join("\t".join(str(c) for c in r) + "\n" for r in rows)


def genomes(d):
    """[(genome, label, study, lineage, copy_group)] of one dataset."""
    out = []

    def add(gid, k, pos):
        study = f"PRJSYN{k:02d}" if k else ""
        label = int(pos < RESISTANT[k]) if k else int(pos < 4)
        out.append([gid, label, study])
    if d == "alphabacter_primus__onemycin":
        for j in range(40):
            add(f"9001.{j:03d}", j // 5 + 1, j % 5)
    elif d == "alphabacter_primus__twomycin":
        for j in range(20, 80):
            add(f"9001.{j:03d}", j // 5 + 1 if j < 70 else 0, j % 5 if j < 70 else j - 70)
    elif d == "betabacter_secundus__onemycin":
        for j in range(40):
            add(f"9002.{j:03d}", j // 5 + 21, j % 5)
    else:
        for j in range(40):
            add(f"9003.{j:03d}", j // 5 + 31, j % 5)
    for i, g in enumerate(out):
        g.append(str(i % 6))
        g.append("C1" if d.startswith("betabacter") and i < 2 else "")
    return out


ENTERED = ["alphabacter_primus__onemycin", "alphabacter_primus__twomycin", "betabacter_secundus__onemycin",
           "gammabacter_tertius__onemycin"]
OTHERS = [
    # (dataset, genomes, with_study, studies, largest_group, reasons)
    ("alphabacter_primus__threemycin", 50, 50, 3, 20, ["3 studies, under 5"]),
    ("betabacter_secundus__twomycin", 50, 30, 6, 30,
     ["30 of 50 genomes with a study, under 80%", "largest group 30 of 50 genomes, over 50%"]),
    ("gammabacter_tertius__twomycin", 50, 50, 6, 30, ["largest group 30 of 50 genomes, over 50%"]),
    ("deltabacter_quartus__onemycin", 30, 30, 5, 8, ["study folds: fold 2 test genomes lack a class"]),
]
LINEAGE = ENTERED[:3]


def split_counts(units, counts):
    """Fill folds in turn with the given number of genomes from a list of genome indices."""
    fold, at = {}, 0
    for f, c in enumerate(counts):
        for i in units[at:at + c]:
            fold[i] = f
        at += c
    return fold


def keep_copy(fold, rows, members):
    """Move the second genome of the copy group to the first's fold, swapping it with a genome of its class there."""
    a, b = members
    if fold[a] == fold[b]:
        return
    c = next(i for i in sorted(fold) if fold[i] == fold[a] and i not in members and rows[i][1] == rows[b][1])
    fold[b], fold[c] = fold[a], fold[b]


def near_equal(total):
    return [total // 5 + (f < total % 5) for f in range(5)]


def plans_of(d, rows, fault):
    n = len(rows)
    res = [i for i in range(n) if rows[i][1]]
    sus = [i for i in range(n) if not rows[i][1]]
    copy = [i for i in range(n) if rows[i][4]]
    out = {}
    for i in DRAWS:
        r = random.Random(f"{d} random {i}")
        a, b = res[:], sus[:]
        r.shuffle(a)
        r.shuffle(b)
        fold = split_counts(a, near_equal(len(a)))
        fold.update(split_counts(b, near_equal(len(b))))
        if copy:
            keep_copy(fold, rows, copy)
            if fault == "copy_group" and i == 1:
                # The second genome of the copy group trades folds with a resistant genome of the next fold.
                a, b = copy
                c = next(j for j in res if fold[j] == (fold[a] + 1) % 5)
                fold[b], fold[c] = fold[c], fold[b]
        out[f"random-d{i}"] = fold
        r = random.Random(f"{d} study {i}")
        studies = sorted({g[2] for g in rows if g[2]})
        r.shuffle(studies)
        single = [j for j in range(n) if not rows[j][2]]
        r.shuffle(single)
        place = {s: k % 5 for k, s in enumerate(studies)}
        fold = {j: place[rows[j][2]] for j in range(n) if rows[j][2]}
        fold.update({j: (len(studies) + k) % 5 for k, j in enumerate(single)})
        if fault == "study_split" and i == 0 and d == ENTERED[0]:
            fold[0] = (fold[0] + 1) % 5
        if fault == "no_class" and i == 2 and d == ENTERED[0]:
            # Every resistant genome of one test fold moves to the next fold.
            f0 = fold[res[0]]
            for j in res:
                if fold[j] == f0:
                    fold[j] = (f0 + 1) % 5
        out[f"study-d{i}"] = fold
        sizes = [sum(1 for j in range(n) if fold[j] == f) for f in range(5)]
        rs = [sum(1 for j in res if fold[j] == f) for f in range(5)]
        r = random.Random(f"{d} matched {i}")
        a, b = res[:], sus[:]
        r.shuffle(a)
        r.shuffle(b)
        if fault == "matched" and i == 1 and d == ENTERED[0]:
            rs = [rs[0] - 1, rs[1] + 1] + rs[2:]
        mf = split_counts(a, rs)
        mf.update(split_counts(b, [s - x for s, x in zip(sizes, rs)]))
        if copy:
            keep_copy(mf, rows, copy)
        out[f"matched-d{i}"] = mf
    return out


def scores(rows, wrong):
    """{learner: [score]} with every learner wrong on the genome indices in wrong."""
    out = {}
    for l in LEARNERS:
        col = []
        for i, g in enumerate(rows):
            call = bool(g[1]) != (i in wrong)
            if l == "scm":
                col.append("1.0" if call else "0.0")
            else:
                col.append(repr(SURE[l] if call else 1 - SURE[l]))
        out[l] = col
    return out


def counts(rows, sc):
    cells = []
    for l in LEARNERS:
        tp = fn = tn = fp = 0
        for g, s in zip(rows, sc[l]):
            call = float(s) >= 0.5
            if g[1]:
                tp += call
                fn += not call
            else:
                fp += call
                tn += not call
        cells += [tp, fn, tn, fp, 0]
    return cells


COUNT_HEAD = [f"{l}_{k}" for l in LEARNERS for k in ("tp", "fn", "tn", "fp", "zero_folds")]


def plan_file(d, name):
    if name == "study-d0":
        return f"{d}.study-whole.json"
    if name == "random-d0":
        return f"{d}.release-random.json"
    return f"{d}.{name}.json"


def make(root, fault=None):
    assert fault is None or fault in FAULTS, fault
    os.makedirs(root, exist_ok=True)
    first_sets, second_sets, fit_index, release_index = [], [], [], []
    fit_lines, first_lines = [], []
    clock = 0
    for d in ENTERED:
        rows = genomes(d)
        n = len(rows)
        species, antibiotic = d.split("__")
        plans = plans_of(d, rows, fault)
        first_study = sorted({g[2] for g in rows if g[2]})[0]
        wrong_study = {i for i, g in enumerate(rows) if g[2] == first_study}
        ids = [g[0] for g in rows]
        entry = {"id": d, "species_slug": species, "genomes": n, "resistant": sum(g[1] for g in rows),
                 "groups": len({g[2] for g in rows if g[2]}) + sum(1 for g in rows if not g[2]),
                 "copy_groups": int(any(g[4] for g in rows)), "plans": {}}
        for name in sorted(plans):
            kind, draw = name.split("-d")
            fold = plans[name]
            doc = {"train": [], "unused": [], "test": [[ids[j] for j in range(n) if fold[j] == f] for f in range(5)]}
            raw = json.dumps(doc, indent=1).encode()
            digest = sha(raw)
            if fault == "fold_sha" and name == "random-d3" and d == ENTERED[1]:
                raw += b"\n"
            write(root, f"folds/{plan_file(d, name)}", raw)
            sizes = [sum(1 for j in range(n) if fold[j] == f) for f in range(5)]
            rs = [sum(1 for j in range(n) if fold[j] == f and rows[j][1]) for f in range(5)]
            entry["plans"][name] = {"kind": kind, "draw": int(draw), "scheme": plan_file(d, name)[len(d) + 1:-5],
                                    "file": plan_file(d, name), "sha256": digest, "sizes": sizes, "resistant": rs}
            wrong = wrong_study if kind == "study" else set()
            sc = scores(rows, wrong)
            if kind != "study":
                sc["scm"] = scores(rows, {0})["scm"]
            labels = [g[1] for g in rows]
            if fault == "labels" and name == "matched-d4" and d == ENTERED[2]:
                labels[5] = 1 - labels[5]
            write(root, f"fits/{d}.{name}.tsv",
                  tsv(["genome", "fold", "label"] + list(LEARNERS),
                      [[g[0], fold[j], labels[j]] + [sc[l][j] for l in LEARNERS] for j, g in enumerate(rows)]))
            cells = counts(rows, sc)
            if fault == "evaluator" and name == "study-d1" and d == ENTERED[3]:
                cells[0] += 1
            fit_index.append([d, name, digest, FINGERPRINT, "14"] + cells)
            clock += 1
            fit_lines.append(f"{DAY}T11:{40 + clock // 60:02d}:{clock % 60:02d}Z fit {d}.{name} done start "
                             f"{DAY}T11:40:00Z;  User time (seconds): 1.00")
        second_sets.append(entry)
        lin_ok = d in LINEAGE
        write(root, f"release/genomes/{d}.tsv",
              tsv(["genome", "label", "study", "lineage", "copy_group", "random_fold"],
                  [g + [plans["random-d0"][j]] for j, g in enumerate(rows)]))
        lin_cells = [0] * len(COUNT_HEAD)
        if lin_ok:
            lf = [int(g[3]) % 5 for g in rows]
            if fault == "lineage" and d == ENTERED[0]:
                lf[0] = (lf[0] + 1) % 5
            sc = scores(rows, {0, 1})
            write(root, f"release/lineage/{d}.tsv",
                  tsv(["genome", "fold", "label"] + list(LEARNERS),
                      [[g[0], lf[j], g[1]] + [sc[l][j] for l in LEARNERS] for j, g in enumerate(rows)]))
            lin_cells = counts(rows, sc)
        status = "scored" if lin_ok else "invalid lineage outer evaluation: fold 1 validation lacks one or both classes"
        release_index.append([d, species, antibiotic, FINGERPRINT, "13" if fault == "revision" else "14",
                              sha(d.encode()), status, len({g[3] for g in rows}), 5] + lin_cells)
        with_study = sum(1 for g in rows if g[2])
        studies = len({g[2] for g in rows if g[2]})
        largest = max(sum(1 for g in rows if g[2] == s) for s in {g[2] for g in rows if g[2]})
        reasons = []
        if fault == "reasons" and d == ENTERED[0]:
            reasons = ["8 studies, under 5"]
        first_sets.append({"id": d, "species_slug": species,
                           "counts": {"genomes": str(n), "with_study": str(with_study), "studies": str(studies),
                                      "largest_study": str(largest), "groups": str(entry["groups"]),
                                      "largest_group": str(largest), "dup_groups": "0"},
                           "eligible": True, "reasons": reasons,
                           "fold_files": {k: {"file": f"fold-files/{d}.{k}.json",
                                              "sha256": entry["plans"][p]["sha256"]}
                                          for k, p in (("release-random", "random-d0"), ("study-whole", "study-d0"))}})
        if fault == "draw0" and d == ENTERED[1]:
            first_sets[-1]["fold_files"]["study-whole"]["sha256"] = "f" * 64
        for k in ("release-random", "study-whole"):
            first_lines.append(f"{DAY}T10:{len(first_lines):02d}:00Z {d}.{k} done  User time (seconds): 1.00")
    for d, n, ws, st, lg, reasons in OTHERS:
        first_sets.append({"id": d, "species_slug": d.split("__")[0],
                           "counts": {"genomes": str(n), "with_study": str(ws), "studies": str(st),
                                      "largest_study": str(lg), "groups": str(st), "largest_group": str(lg),
                                      "dup_groups": "0"},
                           "eligible": False, "reasons": reasons})
    first_sets.sort(key=lambda e: e["id"])
    d1 = {"key": "9001.003" if fault == "genome_id" else "SYNTHETIC-first", "rule": RULE, "meta_sha256_gz": META, "run_order": ENTERED,
          "datasets": first_sets}
    d2 = {"key": "SYNTHETIC-second", "folds": 5, "draws": 5, "study_rule": {"0": "lane", "1": "B", "2": "B",
                                                                           "3": "B", "4": "B"},
          "meta_sha256_gz": META, "run_order": ENTERED, "datasets": second_sets}
    design1 = json.dumps(d1, indent=1).encode()
    design2 = json.dumps(d2, indent=1).encode()
    write(root, "first/design.json", design1)
    write(root, "second/design.json", design2)
    write(root, "fits/index.tsv", tsv(["dataset", "plan", "fold_sha256", "input_fingerprint", "revision"]
                                      + COUNT_HEAD, fit_index))
    write(root, "release/index.tsv", tsv(["dataset", "species", "antibiotic", "input_fingerprint", "revision",
                                          "report_sha256", "lineage_status", "lineage_groups", "lineage_folds"]
                                         + COUNT_HEAD, release_index))
    write(root, "release/release.tsv", tsv(["key", "value"], [["evaluation_revision", "14"],
                                                              ["generated_at", f"{DAY}T08:00:00Z"]]))
    write(root, "checks.tsv", tsv(["check", "count"], [["synthetic_inputs", 1]]))
    # Variants: the plan's rule and two others, counted from the eight datasets' counts.
    counted = [(e["species_slug"], *(int(e["counts"][k]) for k in ("genomes", "with_study", "studies",
                                                                     "largest_group"))) for e in first_sets]
    variants = []
    for ms, sh, gr in ((3, "0.5", "0.5"), (5, "0.8", "0.5"), (5, "0.8", "0.65")):
        keep = [c for c in counted
                if c[3] >= ms and Fraction(c[2], c[1]) >= Fraction(sh) and Fraction(c[4], c[1]) <= Fraction(gr)]
        variants.append([ms, sh, gr, len(keep), len({c[0] for c in keep}), sum(c[1] for c in keep)])
    if fault == "variant":
        variants[0][3] += 1
    write(root, "first/eligibility-variants.tsv",
          tsv(["min_studies", "min_study_share", "max_group_share", "datasets", "species", "cohort_genome_pairs"],
              variants))
    var = []
    for e in second_sets:
        for i in range(1, 5):
            for rule in ("A", "B"):
                changed = (0 if e["id"] == ENTERED[0] else i) if rule == "A" else 2 * i + 4
                var.append([e["id"], e["species_slug"], e["genomes"], e["groups"], i, rule,
                            f"{changed / e['genomes']:.6f}"])
    write(root, "second/study-rule-variation.tsv",
          tsv(["dataset", "species", "genomes", "groups", "draw", "rule", "changed_share_vs_d0"], var))
    write(root, "first/times.tsv", tsv(["role", "file", "mtime"], [
        ["fold_builder", "scripts/make_folds.py", f"{DAY}T09:00:00.00Z"],
        ["design", "results/design.json", f"{DAY}T09:10:00.00Z"],
        ["fold_files_first", "results/fold-files", f"{DAY}T09:09:00.00Z"],
        ["fold_files_last", "results/fold-files", f"{DAY}T09:10:00.00Z"],
        ["variant_counter", "scripts/eligibility_variants.py", f"{DAY}T09:20:00.00Z"],
        ["variant_table", "results/counts/eligibility-variants.tsv", f"{DAY}T09:20:00.00Z"],
        ["plan", "PLAN.md", f"{DAY}T09:30:00.00Z"],
        ["plan_hashes", "PLAN.sha256", f"{DAY}T09:40:00.00Z"]]))
    write(root, "second/times.tsv", tsv(["role", "file", "mtime"], [
        ["plan_builder", "scripts/make_plans.py", f"{DAY}T11:00:00.50Z"],
        ["design", "results/design.json", f"{DAY}T11:10:00.00Z"],
        ["plan_shapes", "results/plan-shapes.tsv", f"{DAY}T11:10:00.00Z"],
        ["rule_variation", "results/study-rule-variation.tsv", f"{DAY}T11:10:00.00Z"],
        ["plan", "PLAN.md", f"{DAY}T11:20:00.00Z"],
        ["plan_hashes", "PLAN.sha256", f"{DAY}T11:30:00.00Z"]]))
    plan1 = ("# SYNTHETIC first plan\n\nA SYNTHETIC plan for the tests; it fixes nothing. The species interval "
             "takes 10,000 resamples of species.\n").encode()
    plan2 = ("# SYNTHETIC second plan\n\nA SYNTHETIC plan for the tests. Rule A moves too little; I fixed at "
             f"{FIXED}, before measuring it, that a median under 0.10 over the 4 datasets is too little,\nand "
             "that rule B places draws 1 to 4.\n").encode()
    if fault == "rule_time":
        plan2 = plan2.replace(b"I fixed at", b"I chose at")
    if fault == "rule_line":
        plan2 = plan2.replace(b"under 0.10", b"under 0.20")
    for exp, plan, design, stamp in (("first", plan1, design1, "09:40:00"), ("second", plan2, design2, "11:30:00")):
        write(root, f"{exp}/PLAN.sha256", f"# SYNTHETIC hash list, taken at {DAY}T{stamp}Z\n"
                                          f"{sha(plan)}  PLAN.md\n{sha(design)}  results/design.json\n")
        if fault == "plan_hash" and exp == "second":
            plan = plan + b"An edit after the hash.\n"
        write(root, f"{exp}/PLAN.md", plan)
    write(root, "first/run.log", f"{DAY}T09:50:00Z start, SYNTHETIC\n" + "".join(x + "\n" for x in first_lines)
          + f"{DAY}T10:30:00Z END rc=0\n")
    if fault == "run_log":
        fit_lines = fit_lines[1:]
    write(root, "second/run.log", f"{DAY}T11:35:00Z start, SYNTHETIC\n{DAY}T11:40:00Z fits start, "
                                  f"{len(fit_lines)} in the queue\n" + "".join(x + "\n" for x in fit_lines)
          + f"{DAY}T12:00:00Z fits end\n")
    if fault == "extra_input":
        write(root, "first/notes.txt", "a file the set does not hold\n")
    if fault == "extra_fit":
        with open(os.path.join(root, f"fits/{ENTERED[0]}.study-d4.tsv"), "rb") as f:
            write(root, f"fits/{ENTERED[0]}.study-d5.tsv", f.read())
    if fault == "extra_column":
        path = os.path.join(root, f"fits/{ENTERED[1]}.random-d2.tsv")
        with open(path, encoding="utf-8") as f:
            lines = f.read().split("\n")
        write(root, f"fits/{ENTERED[1]}.random-d2.tsv",
              "\n".join(x + ("\tweight" if i == 0 else "\t1") if x else x for i, x in enumerate(lines)))
    if fault == "extra_dataset":
        with open(os.path.join(root, f"release/genomes/{ENTERED[3]}.tsv"), "rb") as f:
            write(root, "release/genomes/epsilonbacter_quintus__onemycin.tsv", f.read())
    return root


def brute_interval(groups):
    """The exact species interval by listing every one of the k**k ordered draws (small k only)."""
    k = len(groups)
    meds = []
    for draw in itertools.product(range(k), repeat=k):
        vals = sorted(v for s in draw for v in groups[s])
        m = len(vals)
        meds.append(Fraction(vals[(m - 1) // 2] + vals[m // 2]) / 2)
    meds.sort()
    tot = len(meds)
    lo = next(x for c, x in enumerate(meds, 1) if 40 * c > tot)
    hi = next(x for c, x in enumerate(meds, 1) if 40 * c >= 39 * tot)
    return lo, hi
