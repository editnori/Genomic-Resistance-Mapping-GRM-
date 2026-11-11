"""Hand-built inputs with known answers for the outside-study kit's tests.

A made-up species, `testus_bacterium`, with three dataset reports:

  drugA  eligible; lineage and random splits; SCM, CART and a complete lasso;
         a catalog rule [["C:BETA-LACTAM"]]
  drugB  eligible; the lineage split unavailable; the lasso incomplete (it
         reads a feature the report does not list); no catalog rule
  drugC  pinned but not eligible

Each report's cohort is the eight genomes 999.1 to 999.8. The frozen list
is 120 labelled genomes (1001.1 to 1001.120) in three BioProjects, plus
genomes that test the exclusions and the label rules:

  1001.121  R on drugA, no sequence
  1001.122  S on drugA, species check failed
  1001.123  R on drugA, fetched and passed, not counted
  1001.124  call I on drugA          (BioProject P4)
  1001.125  call conflict on drugA   (BioProject P4)
  1001.126  call none on drugA       (BioProject P4)
  1001.1    borderline in the species check: kept, scored, flagged

Every model here calls resistant exactly when feature 1 is present, except
drugB's CART, which is one leaf at 0.3 (always susceptible). So the counts
follow from which genomes carry feature 1:

  project  R   FN  S   FP
  P1       30  3   30  3
  P2       20  4   20  4
  P3       10  5   10  0

drugA SCM, CART and lasso: TP 48, FN 12, FP 7, TN 53.
The known-gene calls are resistant for the R genomes of P1 and P2 only:
TP 50, FN 10, FP 0, TN 60.
drugB labels the first 25 R and first 25 S genomes of P1 with the same
labels: SCM TP 22, FN 3, FP 3, TN 22; CART TP 0, FN 25, FP 0, TN 25.
"""

import base64
import hashlib
import importlib.util
import json
import math
import shutil
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
KIT = ROOT / "scripts" / "outside-study"


def load(name):
    spec = importlib.util.spec_from_file_location(f"outside_study_{name}_under_test", KIT / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


S = load("score")

SPECIES = "testus_bacterium"
NAME = "Testus bacterium"
A, B, C = (f"{SPECIES}__drug{x}" for x in "ABC")
K = 5
KMER = {1: "AAAAC", 2: "AAAAG", 3: "AAAAT"}
PRIMARY = "primary (other-study)"
COHORT = [f"999.{i}" for i in range(1, 9)]
COHORT_R = {"999.1", "999.2", "999.3", "999.4"}
# Feature presence in the cohort (published bits).
COHORT_F = {1: {"999.1", "999.2", "999.5"}, 2: {"999.1", "999.3"}, 3: {"999.4"}}
COHORT_ST = {g: ("MLST.testus.1" if i < 4 else "MLST.testus.2") for i, g in enumerate(COHORT)}
PROJECTS = (("P1", 30, 3, 30, 3), ("P2", 20, 4, 20, 4), ("P3", 10, 5, 10, 0))
MLST_CYCLE = ("MLST.testus.1", "MLST.testus.7", "")

# The reports' cross-validation rows: (sensitivity, specificity) per split
# and report learner, chosen so each reading has a known answer against the
# new intervals (scm 48/60 and 53/60: about [0.677, 0.892] and [0.774, 0.952];
# known gene 50/60 and 60/60: about [0.715, 0.917] and [0.940, 1]).
REFS = {
    A: {
        "lineage": {"scm": (0.95, 0.90), "cart": (0.80, 0.50), "lasso": (0.60, 0.99), "catalog": (0.50, 1.0)},
        "random": {"scm": (0.85, 0.80), "cart": (0.81, 0.82), "lasso": (0.83, 0.84), "catalog": (0.86, 0.87)},
    },
    B: {
        "lineage": None,
        "random": {"scm": (0.70, 0.71), "cart": (0.72, 0.73), "lasso": (0.74, 0.75), "catalog": (0.76, 0.77)},
    },
    C: {
        "lineage": {"scm": (0.5, 0.5), "cart": (0.5, 0.5), "lasso": (0.5, 0.5), "catalog": (0.5, 0.5)},
        "random": {"scm": (0.5, 0.5), "cart": (0.5, 0.5), "lasso": (0.5, 0.5), "catalog": (0.5, 0.5)},
    },
}


def bits(present):
    """Base64 of little-endian u64 words, bit i for cohort sample i."""
    words = [0] * ((len(COHORT) + 63) // 64)
    for i, g in enumerate(COHORT):
        if g in present:
            words[i // 64] |= 1 << (i % 64)
    return base64.b64encode(b"".join(w.to_bytes(8, "little") for w in words)).decode()


def ref_metrics(sens, spec):
    # Counts out of 100 R and 100 S, so the counts carried beside the
    # estimate are known too.
    tp, tn = round(100 * sens), round(100 * spec)
    return {"n": 200, "tp": tp, "fn": 100 - tp, "tn": tn, "fp": 100 - tn,
            "sensitivity": sens, "specificity": spec, "vme": 1 - sens, "me": 1 - spec}


# Gene calls of the cohort, as grm count writes them, sorted by group.
def hit(group, cls, subclass):
    return {"family": group, "group": group, "allele": group, "containment": 1.0, "class": cls, "subclass": subclass}


BLAZ = hit("blaZ", "BETA-LACTAM", "PENICILLIN")
MECA = hit("mecA", "BETA-LACTAM", "METHICILLIN")
TET = hit("tet(38)", "TETRACYCLINE", "TETRACYCLINE")
COHORT_HITS = {
    "999.1": [BLAZ, MECA], "999.2": [MECA], "999.3": [TET], "999.4": [],
    "999.5": [BLAZ], "999.6": [TET], "999.7": [], "999.8": [TET],
}
GROUP_NAMES = ["blaZ", "mecA", "tet(38)"]
# The catalog table and counts, worked by hand from COHORT_HITS.
CATALOG_GROUPS = [
    {"group": "tet(38)", "class": "TETRACYCLINE", "subclass": "TETRACYCLINE", "counts_for_rule": False,
     "n_resistant": 1, "n_susceptible": 2},
    {"group": "blaZ", "class": "BETA-LACTAM", "subclass": "PENICILLIN", "counts_for_rule": True,
     "n_resistant": 1, "n_susceptible": 1},
    {"group": "mecA", "class": "BETA-LACTAM", "subclass": "METHICILLIN", "counts_for_rule": True,
     "n_resistant": 2, "n_susceptible": 0},
]
CATALOG_RESISTANT = {"999.1": True, "999.2": True, "999.5": True}
RULE = [["C:BETA-LACTAM"]]


def make_report(dataset):
    features = [1, 2] if dataset != B else [1]
    rows = []
    for f in features:
        rows.append({"feature": f, "n_kmers": 10 * f, "kmer": KMER[f],
                     "n_resistant": len(COHORT_F[f] & COHORT_R),
                     "n_susceptible": len(COHORT_F[f] - COHORT_R), "bits": bits(COHORT_F[f])})
    if dataset == B:
        cart = {"nodes": [{"value": 0.3, "n": 8, "n_resistant": 4, "depth": 0, "split": None}]}
        lasso = {"intercept": -1.0, "coef": [[1, 2.5], [99, 0.3]], "lambda": 0.1}
    else:
        cart = {"nodes": [
            {"value": 0.5, "n": 8, "n_resistant": 4, "depth": 0,
             "split": {"feature": 1, "gain": 0.3, "absent": 1, "present": 2}},
            {"value": 0.1, "n": 5, "n_resistant": 2, "depth": 1, "split": None},
            {"value": 0.9, "n": 3, "n_resistant": 2, "depth": 1, "split": None},
        ]}
        lasso = {"intercept": -1.0, "coef": [[1, 2.5], [2, 0.3]], "lambda": 0.1}
    splits = []
    for scheme in ("random", "lineage"):
        refs = REFS[dataset][scheme]
        if refs is None:
            splits.append({"scheme": scheme, "unavailable_reason": "fold 4 validation lacks one or both classes",
                           "learners": []})
            continue
        learners = [{"learner": lr, "metrics": ref_metrics(*refs[lr])} for lr in ("scm", "cart", "lasso", "catalog")]
        learners.insert(2, {"learner": "rf", "metrics": ref_metrics(0.5, 0.5)})
        splits.append({"scheme": scheme, "unavailable_reason": None, "learners": learners})
    samples = []
    for g in COHORT:
        samples.append({"id": g, "resistant": g in COHORT_R, "mlst": COHORT_ST[g],
                        "country": "Testland", "groups": [GROUP_NAMES.index(h["group"]) for h in COHORT_HITS[g]],
                        "catalog_resistant": g in CATALOG_RESISTANT})
    catalog = None
    if dataset != B:
        catalog = {"version": "2026-08-07.1", "rule": RULE, "groups": CATALOG_GROUPS,
                   "metrics": {"n": 8, "tp": 2, "fn": 2, "fp": 1, "tn": 3}}
    return {
        "schema": 2, "id": dataset, "species": NAME, "species_slug": SPECIES,
        "antibiotic": dataset.split("__")[1],
        "input_fingerprint": "f00d" + dataset[-1],
        "evaluation": {"revision": 12, "valid": True, "reasons": []},
        "counts": {"samples": 8, "resistant": 4, "susceptible": 4, "k": K},
        "features": rows,
        "models": {"scm": {"kind": "conjunction", "rules": [{"feature": 1, "presence": True}]},
                   "cart": cart, "cart_depth": 2, "lasso": lasso},
        "splits": splits,
        "catalog": catalog,
        "catalog_groups": GROUP_NAMES,
        "samples": samples,
    }


def report_bytes(report):
    return (json.dumps(report, sort_keys=True, indent=1) + "\n").encode()


def sha(data):
    return hashlib.sha256(data).hexdigest()


# What `grm predict` would write, worked out here without score.py.
def model_call(dataset, learner, found):
    f1 = found.get(1, False)
    if learner == "scm":
        return ("resistant" if f1 else "susceptible"), None
    if learner == "cart":
        if dataset == B:
            return "susceptible", 0.3
        return ("resistant", 0.9) if f1 else ("susceptible", 0.1)
    if learner == "lasso":
        eta = -1.0 + (2.5 if f1 else 0.0) + (0.3 if found.get(2, False) else 0.0)
        p = 1.0 / (1.0 + math.exp(-eta))
        return ("resistant" if p >= 0.5 else "susceptible"), p
    raise ValueError(learner)


MODEL_FEATURES = {
    (A, "scm"): [1], (A, "cart"): [1], (A, "lasso"): [1, 2],
    (B, "scm"): [1], (B, "cart"): [],
    (C, "scm"): [1], (C, "cart"): [1], (C, "lasso"): [1, 2],
}
LEARNERS = {A: ("scm", "cart", "lasso"), B: ("scm", "cart"), C: ("scm", "cart", "lasso")}


def make_row(dataset, report_sha, genome, learner, found):
    call, score = model_call(dataset, learner, found)
    feats = [{"feature": f, "kmer": KMER[f], "found": bool(found.get(f, False)), "n_kmers": 10 * f,
              "kmers_found": 10 * f if found.get(f, False) else 0}
             for f in MODEL_FEATURES[(dataset, learner)]]
    return {"contract": "grm-predict-1", "dataset": dataset, "report_sha256": report_sha,
            "input_fingerprint": "f00d" + dataset[-1], "evaluation_revision": 12, "k": K,
            "genome": genome, "source": "kset", "assembly": None, "learner": learner,
            "call": call, "score": score, "features": feats, "partial_features": 0}


def new_genomes():
    """The 120 labelled frozen genomes, with their known properties."""
    out, n = [], 0
    for proj, nr, fn, ns, fp in PROJECTS:
        for i in range(nr):
            n += 1
            out.append({"id": f"1001.{n}", "project": proj, "resistant": True, "f1": i >= fn,
                        "f2": proj == "P1", "known": proj in ("P1", "P2")})
        for i in range(ns):
            n += 1
            out.append({"id": f"1001.{n}", "project": proj, "resistant": False, "f1": i < fp,
                        "f2": proj == "P1", "known": False})
    for g in out:
        k = int(g["id"].split(".")[1])
        g["country"] = "Atlantis" if g["project"] == "P1" else "Lemuria"
        g["mlst"] = MLST_CYCLE[k % 3]
    return out


EXTRA = {
    "1001.121": {"label": "R", "fasta": False},
    "1001.122": {"label": "S", "verdict": "fail"},
    "1001.123": {"label": "R", "counted": False},
    "1001.124": {"label": "I"},
    "1001.125": {"label": "conflict"},
    "1001.126": {"label": "none"},
}
BORDERLINE = "1001.1"
SAME_STUDY = ["2001.1", "2001.2", "2001.3", "2001.4"]
DROPPED = "3001.1"


def pair(g, d, call, pop, eligible, frozen, same):
    return {"genome_id": g, "dataset_id": d, "call": call,
            "records_R": int(call in ("R", "conflict")), "records_S": int(call in ("S", "conflict")),
            "records_I": int(call == "I"), "records_no_call": int(call == "none"),
            "genome_population": pop, "genome_drop_reason": "" if pop != "dropped" else "older study",
            "dataset_eligible_in_population": "yes" if eligible else "no",
            "in_frozen_test": "yes" if frozen else "no", "in_same_study_set": "yes" if same else "no"}


def candidate(g, pop, project, country, mlst, on):
    return {"genome_id": g, "species": NAME, "species_slug": SPECIES, "bioproject": project,
            "isolation_country": country, "mlst": mlst, "population": pop,
            "on_eligible_dataset": "yes" if on else "no"}


def world(drug_b_resistant=25):
    """Every input the scoring needs, as Python objects. `drug_b_resistant`
    sets how many R genomes drugB labels (25 makes it just reportable)."""
    reports = {d: make_report(d) for d in (A, B, C)}
    shas = {d: sha(report_bytes(r)) for d, r in reports.items()}
    gens = new_genomes()
    cands, pairs = [], []
    for g in gens:
        cands.append(candidate(g["id"], PRIMARY, g["project"], g["country"], g["mlst"], True))
        pairs.append(pair(g["id"], A, "R" if g["resistant"] else "S", PRIMARY, True, True, False))
    p1r = [g for g in gens if g["project"] == "P1" and g["resistant"]][:drug_b_resistant]
    p1s = [g for g in gens if g["project"] == "P1" and not g["resistant"]][:25]
    for g in p1r + p1s:
        pairs.append(pair(g["id"], B, "R" if g["resistant"] else "S", PRIMARY, True, True, False))
    # A non-eligible dataset's labels never count.
    pairs.append(pair(gens[0]["id"], C, "S", PRIMARY, False, False, False))
    for g, e in EXTRA.items():
        cands.append(candidate(g, PRIMARY, "P4" if e["label"] in ("I", "conflict", "none") else "P1",
                               "Atlantis", "", True))
        lab = e["label"]
        pairs.append(pair(g, A, lab, PRIMARY, True, lab in ("R", "S"), False))
    for i, g in enumerate(SAME_STUDY):
        cands.append(candidate(g, "same-study", "P9", "Atlantis", "", True))
        pairs.append(pair(g, A, "R" if i % 2 == 0 else "S", "same-study", True, False, True))
    cands.append(candidate(DROPPED, "dropped", "P8", "Atlantis", "", False))
    pairs.append(pair(DROPPED, A, "R", "dropped", False, False, False))
    datasets = [
        {"dataset_id": A, "species_slug": SPECIES, "eligible_primary": "yes", "eligible_same_study": "yes"},
        {"dataset_id": B, "species_slug": SPECIES, "eligible_primary": "yes", "eligible_same_study": "no"},
        {"dataset_id": C, "species_slug": SPECIES, "eligible_primary": "no", "eligible_same_study": "no"},
    ]
    frozen = [g["id"] for g in gens] + list(EXTRA)
    availability, verdicts, counted = {}, {}, {}
    for g in frozen:
        e = EXTRA.get(g, {})
        fasta = e.get("fasta", True)
        availability[g] = {"genome_id": g, "fasta": "yes" if fasta else "no",
                           "detail": "" if fasta else "no contigs in BV-BRC"}
        if fasta:
            v = e.get("verdict", "borderline" if g == BORDERLINE else "pass")
            verdicts[g] = {"genome_id": g, "verdict": v,
                           "reason": "nearer another species" if v == "fail" else "", "distance": "0.01", "gc": "0.33"}
            ok = e.get("counted", True) and v != "fail"
            counted[g] = {"genome_id": g, "counted": "yes" if ok else "no",
                          "detail": "" if ok else "count failed: truncated gzip"}
        else:
            counted[g] = {"genome_id": g, "counted": "no", "detail": "no FASTA"}
    found = {g["id"]: {1: g["f1"], 2: g["f2"]} for g in gens}
    for g in ("1001.124", "1001.125", "1001.126"):
        found[g] = {1: True, 2: False}
    survivors = [g for g in frozen if g in found]
    rows = {A: [], B: []}
    for d in (A, B):
        for g in survivors:
            for lr in LEARNERS[d]:
                rows[d].append(make_row(d, shas[d], g, lr, found[g]))
    known = {g["id"]: g["known"] for g in gens}
    kg = {"unavailable": None,
          "datasets": {A: {"dataset": A, "rule": json.dumps(RULE, separators=(",", ":")),
                           "catalog_version": "2026-08-07.1", "replay": "pass", "samples": "8"},
                       B: {"dataset": B, "rule": "", "catalog_version": "", "replay": "no rule published",
                           "samples": "0"}},
          "calls": {(A, g): {"dataset": A, "genome_id": g, "groups": "", "rule_groups": "",
                             "catalog_resistant": "yes" if known.get(g) else "no"} for g in survivors}}
    return {
        "reports": reports, "shas": shas, "genomes": gens, "candidates": cands, "pairs": pairs,
        "datasets": datasets, "availability": availability, "verdicts": verdicts, "counted": counted,
        "rows": rows, "known_gene": kg, "survivors": survivors, "found": found,
    }


def inputs(w, rehearsal=False, population="primary"):
    """The dict `score.score` takes."""
    return {
        "population": population, "species": SPECIES, "rehearsal": rehearsal,
        "candidates": w["candidates"], "pairs": w["pairs"], "datasets": w["datasets"],
        "reports": {d: (r, w["shas"][d]) for d, r in w["reports"].items()},
        "species_reports": [w["reports"][d] for d in sorted(w["reports"])],
        "availability": w["availability"], "verdicts": w["verdicts"], "counted": w["counted"],
        "rows": w["rows"], "known_gene": w["known_gene"], "hashes": {"fixture": "1"},
    }


CAND_COLS = ["genome_id", "species", "species_slug", "bioproject", "isolation_country", "mlst",
             "population", "on_eligible_dataset"]
PAIR_COLS = ["genome_id", "dataset_id", "call", "records_R", "records_S", "records_I", "records_no_call",
             "genome_population", "genome_drop_reason", "dataset_eligible_in_population", "in_frozen_test",
             "in_same_study_set"]
DS_COLS = ["dataset_id", "species_slug", "eligible_primary", "eligible_same_study"]


def write_world(w, root):
    """Writes the world as files in the layout score.py's `load_inputs`
    reads. Returns the paths and the overrides `score.main` takes."""
    root = Path(root)
    reports, tables, work = root / "reports", root / "tables", root / "work"
    for p in (reports, tables, work / "fetch", work / "species-check", work / "count", work / "predict",
              work / "known-gene"):
        p.mkdir(parents=True, exist_ok=True)
    pins = []
    for d, r in sorted(w["reports"].items()):
        data = report_bytes(r)
        (reports / f"{d}.json").write_bytes(data)
        pins.append(f"{sha(data)}  {d}.json")
    (root / "reports.sha256").write_text("\n".join(pins) + "\n")
    S.write_tsv(tables / "candidates.tsv", CAND_COLS, w["candidates"])
    S.write_tsv(tables / "pairs.tsv", PAIR_COLS, w["pairs"])
    S.write_tsv(tables / "datasets.tsv", DS_COLS, w["datasets"])
    (root / "tables.sha256").write_text("".join(
        f"{S.sha256_file(tables / n)}  {n}\n" for n in ("candidates.tsv", "pairs.tsv", "datasets.tsv")))
    S.write_tsv(work / "fetch/availability.tsv", ["genome_id", "fasta", "detail"], w["availability"].values())
    S.write_tsv(work / "species-check/verdicts.tsv", ["genome_id", "verdict", "reason", "distance", "gc"],
                w["verdicts"].values())
    S.write_tsv(work / "count/counted.tsv", ["genome_id", "counted", "detail"], w["counted"].values())
    for d, rows in w["rows"].items():
        (work / "predict" / f"{d}.jsonl").write_text("".join(json.dumps(r, sort_keys=True) + "\n" for r in rows))
    kg = w["known_gene"]
    if kg["unavailable"] is not None:
        (work / "known-gene/unavailable.txt").write_text(kg["unavailable"] + "\n")
    else:
        S.write_tsv(work / "known-gene/datasets.tsv", ["dataset", "rule", "catalog_version", "replay", "samples"],
                    kg["datasets"].values())
        S.write_tsv(work / "known-gene/calls.tsv", ["dataset", "genome_id", "groups", "rule_groups",
                                                    "catalog_resistant"], kg["calls"].values())
    return {"reports": reports, "tables": tables, "work": work, "pins": root / "reports.sha256",
            "table_pins": root / "tables.sha256"}


def git(root, *args):
    return subprocess.run(["git", "-C", str(root), "-c", "user.name=t", "-c", "user.email=t@example.invalid",
                           "-c", "commit.gpgsign=false", *args], capture_output=True, text=True, check=True)


def frozen_kit(root):
    """A git repository holding stand-in frozen files and a FROZEN.sha256
    that pins them, committed: a kit root the gate passes."""
    root = Path(root)
    (root / "docs").mkdir(parents=True)
    (root / "scripts/outside-study").mkdir(parents=True)
    (root / "docs/outside-study-test.md").write_text("plan\n")
    (root / "scripts/outside-study/score.py").write_text("score\n")
    (root / "scripts/outside-study/run.py").write_text("run\n")
    (root / "scripts/outside-study/FROZEN.sha256").write_text("".join(
        f"{S.sha256_file(root / rel)}  {rel}\n" for rel in S.FROZEN_FILES))
    git(root, "init", "-q")
    git(root, "add", "-A")
    git(root, "commit", "-q", "-m", "freeze")
    return root


def copy_kit(dst):
    """The kit's own scripts and pins, copied, for tests that edit them."""
    shutil.copytree(KIT, dst, ignore=shutil.ignore_patterns("__pycache__"))
    return Path(dst)


def rehearsal_world():
    """A rehearsal: the frozen list is the cohort's own eight genomes,
    labelled as the reports label them, with rows built from their bits."""
    w = world()
    ids = COHORT
    w["candidates"] = [candidate(g, PRIMARY, "REHEARSAL-" + "ABC"[i % 3], "Testland", COHORT_ST[g], True)
                       for i, g in enumerate(ids)]
    w["pairs"] = [pair(g, d, "R" if g in COHORT_R else "S", PRIMARY, True, True, False)
                  for d in (A, B) for g in ids]
    w["availability"] = {g: {"genome_id": g, "fasta": "yes", "detail": ""} for g in ids}
    w["verdicts"] = {g: {"genome_id": g, "verdict": "pass", "reason": ""} for g in ids}
    w["counted"] = {g: {"genome_id": g, "counted": "yes", "detail": ""} for g in ids}
    found = {g: {f: g in COHORT_F[f] for f in (1, 2)} for g in ids}
    w["rows"] = {d: [make_row(d, w["shas"][d], g, lr, found[g]) for g in ids for lr in LEARNERS[d]] for d in (A, B)}
    w["known_gene"]["calls"] = {(A, g): {"catalog_resistant": "yes" if g in CATALOG_RESISTANT else "no"}
                                for g in ids}
    w["survivors"], w["found"] = list(ids), found
    return w
