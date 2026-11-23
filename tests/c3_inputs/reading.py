"""A SYNTHETIC reading of plan block C3 for scripts/c3-inputs/page_tables.py and run_checks.py.

Every species, antibiotic, dataset, call, count, binary and build here is
made up in the test: three made-up species, seven made-up datasets, one of
them in M. The scorer's own functions write the reading's files, so they
have the form a real reading has. Each draw calls 40 resistant and 60
susceptible genomes; a step of 4 resistant genomes moves balanced accuracy
by 0.05. Most arm rows move every dataset by the same step, so the count's
species resamples run on a few rows only and the suite stays quick.
"""

import json
import math
import shutil
from pathlib import Path

from tests.c3_inputs import common

score = common.load("score")
page = common.load("page_tables")

L_COLUMNS = ["dataset", "species_slug", "species", "antibiotic", "samples", "resistant", "susceptible", "lineages",
             "in_M", "in_S", "in_A", "first_run", "first_run_A", "pairs_S+R", "labels_depend_on", "arms"]
# dataset, species, antibiotic, in M, in S, in A
L = [
    ("synthia_coli__drug_a", "Synthia coli", "drug a", 0, 1, 1),
    ("synthia_coli__drug_b", "Synthia coli", "drug b", 0, 0, 1),
    ("synthia_coli__drug_m", "Synthia coli", "drug m", 1, 1, 1),
    ("fakeria_maxima__drug_c", "Fakeria maxima", "drug c", 0, 1, 0),
    ("fakeria_maxima__drug_d", "Fakeria maxima", "drug d", 0, 0, 1),
    ("mockella_minor__drug_e", "Mockella minor", "drug e", 0, 0, 0),
    ("mockella_minor__drug_f", "Mockella minor", "drug f", 0, 0, 0),
]
NAMES = {d: (sp, ab) for d, sp, ab, *_ in L}
SPECIES = sorted({d.split("__")[0] for d, *_ in L})
FIRST = [d for d, _, _, m, _, _ in L if not m]
IN_A = [d for d, _, _, m, _, a in L if a and not m]
IN_S = [d for d, _, _, m, s, _ in L if s and not m]
NOT_A = [d for d in FIRST if d not in IN_A]
LEARNERS = common.LEARNERS
ARMS = page.ARMS
# The resistant genomes each arm calls resistant beyond the default's 28 of 40, on every dataset and learner.
STEP = {"k31": 0, "twin": 0, "unitigs": 4, "k31-counts": 0, "k31-aa": 0, "program": -4, "k31-program": 2}
# Rows that move by dataset: K31 plus counts with the forest, +4 on drug a and -4 on drug b on the lineage split,
# -4 on drug a's study run; the twin with the set covering machine, +1 at seed 46 on drug b.
VARY = {("k31-counts", "rf", "lineage", "synthia_coli__drug_a"): 4,
        ("k31-counts", "rf", "lineage", "synthia_coli__drug_b"): -4,
        ("k31-counts", "rf", "study", "synthia_coli__drug_a"): -4}
TWIN_VARY = {("scm", 46, "synthia_coli__drug_b"): 1}
# Drug e's lineage split did not run at seed 43 for the default and the arms; the program alone's study run on
# drug a was refused for no candidate.
UNSCORED = {("mockella_minor__drug_e", 43)}
REFUSED = {("synthia_coli__drug_a", "program")}
CODE = "c" * 64
GRM = "b" * 64


def calls(tp, tn=45):
    return {"n": 100, "tp": tp, "fp": 60 - tn, "tn": tn, "fn": 40 - tp, "loss_sum": 40.0 - tp / 2}


def constant_calls():
    """The constant model on folds whose training genomes are 40% resistant: every genome called S."""
    p = 0.4
    return {"n": 100, "tp": 0, "fp": 0, "tn": 60, "fn": 40, "loss_sum": -(40 * math.log(p) + 60 * math.log(1 - p))}


def draw(d, arm, split, seed, role, learner, withdrawn=()):
    sp, ab = NAMES[d]
    base = {"dataset": d, "species": sp, "antibiotic": ab, "arm": arm, "split": split, "seed": seed, "role": role,
            "learner": learner}
    if split == "lineage" and arm != "twin" and (d, seed) in UNSCORED:
        return dict(base, scored=0, constant=0, why="split unavailable")
    if split == "study" and (d, arm) in REFUSED:
        return dict(base, scored=1, constant=1, why=page.CONSTANT_WHY, **score.measures(constant_calls()))
    tp = 28 + STEP[arm] + VARY.get((arm, learner, split, d), 0)
    if arm == "twin":
        tp += TWIN_VARY.get((learner, seed, d), 0)
    return dict(base, scored=1, constant=0, why="", **score.measures(calls(tp)))


def draws_of(d, withdrawn=()):
    out = []
    in_a, in_s = d in IN_A, d in IN_S
    arms = ["k31"] + [a for a in ARMS if a not in withdrawn and (a not in page.PROGRAM_ARMS or in_a)]
    for arm in arms:
        for learner in LEARNERS:
            for seed in (42, 43, 44):
                out.append(draw(d, arm, "lineage", seed, "draw", learner))
            if in_s:
                out.append(draw(d, arm, "study", 42, "study", learner))
    for learner in LEARNERS:
        for seed in (45, 46, 47):
            out.append(draw(d, "twin", "lineage", seed, "twin", learner))
    return out


def write_l(path):
    with open(path, "w", encoding="utf-8") as f:
        f.write("\t".join(L_COLUMNS) + "\n")
        for d, sp, ab, m, s, a in L:
            arms = "k31,twin,k31-counts,unitigs,k31-aa" + (",program,k31-program" if a else "")
            cells = [d, d.split("__")[0], sp, ab, 150, 60, 90, 30, m, s, a, int(not m), int(a and not m), 0,
                     "none of the four", arms]
            f.write("\t".join(str(c) for c in cells) + "\n")


def build_record(slug, kind):
    rec = {"format": "grm-c3-build-1", "evaluation_revision": 14, "kind": kind, "species": slug, "genomes": 90,
           "space_k": {"counts": 31, "aa": 25, "unitigs": 16}[kind], "matrix_support": 2, "values_summed": 5000,
           "matrix_patterns": 300, "matrix_values_kept": 4000, "matrix_values_rare": 600, "matrix_values_core": 400}
    if kind == "aa":
        rec["proteins"] = 3000
    if kind == "unitigs":
        rec.update(graph="grm", graph_kmers=80000, unitigs=7000, pairs_held=5000, pairs_in_part=120,
                   kmers_without_unitig=0)
    return rec


def write_species(top, slug, withdrawn=(), code=CODE):
    """score/<slug>/ and run/<slug>/ of one SYNTHETIC species."""
    top = Path(top)
    folder = top / "score" / slug
    folder.mkdir(parents=True)
    draws = [r for d in sorted(FIRST) if d.startswith(slug + "__") for r in draws_of(d, withdrawn)]
    score.write_tsv(folder / "draws.tsv", draws, score.DRAW_COLUMNS)
    score.write_tsv(folder / "datasets.tsv", score.dataset_rows(draws), score.DATASET_COLUMNS)
    checks = [{"file": f"k31/{d}.lineage-s42.pick.json", "check": "exists", "passed": True, "detail": ""}
              for d in sorted({r["dataset"] for r in draws})]
    score.write_tsv(folder / "checks.tsv", checks, ["file", "check", "passed", "detail"])
    run = top / "run" / slug
    (run / "builds").mkdir(parents=True)
    (run / "binary.txt").write_text(f"commit {'a' * 40}\ngrm {GRM}\ncode {code}\nprodigal 2.6.3 {'d' * 64}\n")
    (run / "DONE").write_text(f"DONE {slug} 2026-11-23T10:00:00Z binary {GRM} code {code} results 12 "
                              f"SHA256SUMS {'e' * 64}\n")
    for arm, kind in page.BUILD_OF.items():
        if arm not in withdrawn:
            (run / "builds" / f"{kind}.build.json").write_text(json.dumps(build_record(slug, kind), indent=1))


def write_count(top, withdrawn=()):
    top = Path(top)
    rows = [r for slug in SPECIES for r in page.load(top / "score" / slug / "datasets.tsv")]
    counted = score.count(rows, withdrawn)
    (top / "count").mkdir(exist_ok=True)
    score.write_tsv(top / "count" / "count.tsv", counted, score.COUNT_COLUMNS)
    score.write_tsv(top / "count" / "count-no-constant.tsv", score.count(rows, withdrawn, keep_constant=False),
                    score.COUNT_COLUMNS)
    score.write_tsv(top / "count" / "conditions.tsv", score.conditions(counted), score.CONDITION_COLUMNS)


def build(top, withdrawn=(), codes=None):
    """A SYNTHETIC reading: top/reading and top/lists.tsv. codes maps a species to its binary's code sha256."""
    top = Path(top)
    top.mkdir(parents=True, exist_ok=True)
    write_l(top / "lists.tsv")
    reading = top / "reading"
    shutil.rmtree(reading, ignore_errors=True)
    reading.mkdir(parents=True)
    (reading / "reading.txt").write_text("A SYNTHETIC reading, made up in the test.\n")
    for slug in SPECIES:
        write_species(reading, slug, withdrawn, (codes or {}).get(slug, CODE))
    write_count(reading, withdrawn)
    return reading, top / "lists.tsv"


def rewrite_species(reading, slug, change):
    """Change one species' draws.tsv rows with change(rows) and write its datasets.tsv again from them."""
    folder = Path(reading) / "score" / slug
    rows = page.load(folder / "draws.tsv")
    rows = change(rows)
    score.write_tsv(folder / "draws.tsv", rows, score.DRAW_COLUMNS)
    score.write_tsv(folder / "datasets.tsv", score.dataset_rows([page.typed_draw(r) for r in rows]),
                    score.DATASET_COLUMNS)
