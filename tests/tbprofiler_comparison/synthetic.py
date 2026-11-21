"""A SYNTHETIC input set for scripts/tbprofiler-comparison/record.py: a call file, a compact file of reports, a
grade table, a plan with its addendum and hash list, a scoring log and the counts of one set-aside entry.

Every genome id, label, call, score, drug and count here is made up in the test. The ids have the shape of a real
genome id (digits, a dot, digits) so that the script's scan for them has something to find. The drugs are made-up
names, and the report ids carry the species prefix the script reads.
"""

import gzip
import hashlib
import json
import random
from pathlib import Path

from tests.tbprofiler_comparison.common import LEARNERS, SPLITS

SPECIES = "mycobacterium_tuberculosis"
GENOMES = 160
GROUP = 8  # genomes per lineage id; report 0 joins the first two ids, so 19 groups in all
# (antibiotic, kind, tool drugs, catalog row)
REPORTS = (("alphamycin", "primary", "alphamycin", True), ("betamycin", "primary", "betamycin", False),
           ("gammacillin", "primary", "gammacillin", False), ("deltazid", "primary", "deltazid", False),
           ("omegamycin", "proxy", "alphamycin,betamycin", True))
NOT_NAMED = "zetanide"
GRADES = (("alphamycin", 3, 2, 1, 0, 4), ("betamycin", 0, 4, 2, 1, 3), ("deltazid", 0, 0, 0, 0, 6),
          ("epsilonine", 1, 1, 0, 0, 0), ("gammacillin", 2, 0, 0, 0, 1))
SET_ASIDE = "betamycin"
WAYS = ("who_assoc", "who_assoc_only", "tool_default", "who_assoc_small")


def sha(b):
    return hashlib.sha256(b).hexdigest()


def genome_id(i):
    return f"{99000 + i % 7}.{1000 + i}"


def members(k):
    """The genomes of report k: each genome sits in four of the five reports with calls."""
    return [i for i in range(GENOMES) if (i + k) % 5 != 0]


def lineage(k, i):
    group = i // GROUP
    if k == 0 and group == 1:
        group = 0
    return 100 * k + group


def calls_for(rng, truth, drug_graded):
    """One row's call under each way of reading, its top grade and whether another grade was reported."""
    if rng.random() < 0.04:
        return {w: "no_call" for w in WAYS}, "-", "no"
    right = rng.random() < 0.8
    resistant = (truth if right else not truth) and drug_graded
    grade = ("Assoc w R" if rng.random() < 0.7 else "Assoc w R - Interim") if resistant else "-"
    other = "yes" if not resistant and rng.random() < 0.2 else "no"
    ways = {"who_assoc": resistant, "who_assoc_only": resistant and grade == "Assoc w R",
            "tool_default": resistant or other == "yes", "who_assoc_small": resistant and rng.random() > 0.02}
    return {w: "resistant" if v else "susceptible" for w, v in ways.items()}, grade, other


def score(rng, truth, learner):
    if learner == "scm":
        return 1.0 if (rng.random() < 0.75) == truth else 0.0
    return min(1.0, rng.random() * 0.6 + (0.35 if truth else 0.05))


def write(root, seed=7):
    """Write the whole input set under root; return the counts of the set-aside report as run."""
    root = Path(root)
    root.mkdir(parents=True, exist_ok=True)
    rng = random.Random(seed)
    graded = {g[0] for g in GRADES if g[1] + g[2] > 0}
    head = ["dataset_id", "genome_id", "kind", "tool_drugs"] + list(WAYS) + ["top_grade", "other_grade",
                                                                           "targets_covered", "large_only"]
    call_lines, compact, counts = ["\t".join(head)], [], None
    for k, (name, kind, drugs, has_catalog) in enumerate(REPORTS):
        did = f"{SPECIES}__{name}"
        ids = members(k)
        truth = [rng.random() < 0.5 for _ in ids]
        is_graded = any(d in graded for d in drugs.split(","))
        rows = [calls_for(rng, t, is_graded) for t in truth]
        for i, t, (ways, grade, other) in zip(ids, truth, rows):
            call_lines.append("\t".join([did, genome_id(i), kind, drugs] + [ways[w] for w in WAYS] + [
                grade, other, "no" if ways["who_assoc"] == "no_call" else "yes", "no"]))
        if name == SET_ASIDE:
            tp = sum(1 for t, r in zip(truth, rows) if t and r[0]["who_assoc"] == "resistant")
            fn = sum(1 for t, r in zip(truth, rows) if t and r[0]["who_assoc"] == "susceptible")
            tn = sum(1 for t, r in zip(truth, rows) if not t and r[0]["who_assoc"] == "susceptible")
            fp = sum(1 for t, r in zip(truth, rows) if not t and r[0]["who_assoc"] == "resistant")
            counts = [tp, fp, tn, fn]
        compact.append({
            "id": did, "species_slug": SPECIES, "antibiotic": name, "sid": [genome_id(i) for i in ids],
            "y": [int(t) for t in truth], "lineage": [lineage(k, i) for i in ids],
            "catalog_resistant": [rng.random() < 0.5 for _ in ids] if has_catalog else [None] * len(ids),
            "splits": {s: {"n_groups": 5, "scores": {l: [score(rng, t, l) for t in truth] for l in LEARNERS}}
                       for s in SPLITS}})
    ids = list(range(0, GENOMES, 4))
    compact.append({"id": f"{SPECIES}__{NOT_NAMED}", "species_slug": SPECIES, "antibiotic": NOT_NAMED,
                    "sid": [genome_id(i) for i in ids], "y": [i % 2 for i in ids], "lineage": [i // GROUP for i in ids],
                    "catalog_resistant": [None] * len(ids),
                    "splits": {s: {"scores": {l: [0.5] * len(ids) for l in LEARNERS}} for s in SPLITS}})
    compact.append({"id": "other_species__alphamycin", "species_slug": "other_species", "sid": ["1.1"]})
    (root / "calls.tsv").write_text("\n".join(call_lines) + "\n", encoding="utf-8")
    with gzip.open(root / "rev14_compact.jsonl.gz", "wt", encoding="utf-8") as f:
        for d in compact:
            f.write(json.dumps(d) + "\n")
    grade_head = ["tool_drug", "WHO catalogue v2: Assoc w R", "WHO catalogue v2: Assoc w R - Interim",
                  "tbdb: Uncertain significance", "tbdb: Not assoc w R - Interim", "tbdb: none",
                  "counted_grade_entries"]
    grade_lines = ["\t".join(grade_head)] + ["\t".join([g[0]] + [str(x) for x in g[1:]] + [str(g[1] + g[2])])
                                             for g in GRADES]
    (root / "database-grades-by-drug.tsv").write_text("\n".join(grade_lines) + "\n", encoding="utf-8")
    plan = b"# SYNTHETIC plan\n\nTool: TBProfiler 6.6.3 with its database, commit `3ed1cb99`.\n"
    addendum = b"# SYNTHETIC addendum\n"
    (root / "PLAN.md").write_bytes(plan)
    (root / "PLAN-ADDENDUM.md").write_bytes(addendum)
    (root / "PLAN.sha256").write_text(
        f"{sha(plan)}  PLAN.md\n# recorded 2026-01-01T10:00:00Z, before any phenotype was read\n"
        f"{sha(addendum)}  PLAN-ADDENDUM.md\n# addendum recorded 2026-01-01T10:05:00Z, before any phenotype was read\n",
        encoding="utf-8")
    (root / "finish.log").write_text(
        "2026-01-01T11:00:00Z plan and addendum match PLAN.sha256\n2026-01-01T11:00:04Z adapter done\n"
        "2026-01-01T11:00:05Z input limit done\na free line of the log\n2026-01-01T11:01:00Z evaluator done\n",
        encoding="utf-8")
    tp, fp, tn, fn = counts
    (root / "exploratory.json").write_text(json.dumps({
        "label": "EXPLORATORY: SYNTHETIC",
        "rule_that_sets_an_entry_aside": {
            "result_files_with_a_note": 2, "entries_set_aside_by_drug_and_grade": {f"{SET_ASIDE} | Assoc w R": 1},
            "primary_endpoints_touched": {f"{SPECIES}__{SET_ASIDE}": {
                "counts_as_run": [tp, fp, tn, fn], "counts_with_the_entry_counted": [tp, fp + 1, tn - 1, fn]}}}},
        indent=1), encoding="utf-8")
    return counts
