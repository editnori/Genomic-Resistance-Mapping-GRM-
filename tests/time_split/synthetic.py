"""A SYNTHETIC input set for scripts/time-split/record.py: a design, result files, a study count, plans and logs.

Every genome id, score, day and count here is made up in the test. The ids have the shape of a real genome id
(digits, a dot, digits) so that the script's scan for them has something to find.
"""

import hashlib
import json
import random
from pathlib import Path

from tests.time_split.common import ARMS, LEARNERS

SPECIES = ("alpha_species", "beta_species", "gamma_species")
DRUGS = ("drug_one", "drug_two", "drug_three", "drug_four")


def sha(b):
    return hashlib.sha256(b).hexdigest()


def design(n_per_species=4, test_n=8, test_resistant=4, train_r=30, train_s=40):
    datasets = []
    for i, sp in enumerate(SPECIES):
        for j in range(n_per_species):
            did = f"{sp}__{DRUGS[j]}"
            train = train_r + train_s
            datasets.append({
                "id": did, "species_slug": sp, "antibiotic": DRUGS[j], "cohort": 200 + 10 * j, "cohort_resistant": 90,
                "cut_day": f"2020-0{1 + j}-1{i}", "share_on_or_before_cut": 0.6123, "test_studies": 3,
                "test": {"n": test_n, "resistant": test_resistant}, "train_resistant": train_r,
                "train_susceptible": train_s,
                "train_from": {"past": {"past": train, "later_other_studies": 0, "later_test_studies": 0},
                               "random": {"past": train - 15, "later_other_studies": 10, "later_test_studies": 5},
                               "study-out": {"past": train - 12, "later_other_studies": 12, "later_test_studies": 0}},
                "eligible": True,
                "fold_files": {a: {"file": f"{did}.{a}.json", "sha256": sha(f"{did}|{a}".encode())} for a in ARMS},
            })
    datasets.append({"id": "gamma_species__drug_none", "species_slug": "gamma_species", "eligible": False})
    return {"cut_share": 0.6, "min_test_class": 20, "min_train_class": 25, "datasets": datasets}


def genome_ids(d, k, n):
    return [f"{990000 + 7 * k}.{1000 + i:04d}" for i in range(n)]


def score(rng, truth, learner, arm):
    """A made-up score: scm and part of cart at 0 or 1, the others in between."""
    shift = {"past": -0.15, "random": 0.0, "study-out": -0.08}[arm]
    if learner == "scm" or (learner == "cart" and rng.random() < 0.5):
        right = rng.random() < 0.8 + shift
        return 1.0 if right == truth else 0.0
    p = rng.random() * 0.5 + (0.45 if truth else 0.05) + (shift if truth else -shift)
    return min(max(p, 0.0), 1.0)


def result(d, arm, ids, truth, seed, run):
    rng = random.Random(f"{seed}|{run}|{d['id']}|{arm}")
    learners = []
    for learner in LEARNERS:
        scores = [score(rng, t, learner, arm) for t in truth]
        tp = sum(1 for t, s in zip(truth, scores) if t and s >= 0.5)
        fn = sum(1 for t, s in zip(truth, scores) if t and s < 0.5)
        tn = sum(1 for t, s in zip(truth, scores) if not t and s < 0.5)
        fp = sum(1 for t, s in zip(truth, scores) if not t and s >= 0.5)
        learners.append({"learner": learner,
                         "metrics": {"n": len(truth), "tp": tp, "fn": fn, "tn": tn, "fp": fp,
                                     "sensitivity": tp / (tp + fn), "specificity": tn / (tn + fp)},
                         "scores": scores, "fold_complexity": [3]})
    train = d["train_resistant"] + d["train_susceptible"]
    return {"format": "grm-fold-file-result-v1", "evaluation_revision": 14,
            "fold_file_sha256": d["fold_files"][arm]["sha256"], "id": d["id"], "species_slug": d["species_slug"],
            "scheme": f"first-public-{arm}", "eval_config": {"seed": 42, "screened_features": 4000, "inner_folds": 3},
            "samples": [{"genome_id": g, "resistant": t, "fold": 0} for g, t in zip(ids, truth)],
            "folds": [{"train": train, "train_resistant": d["train_resistant"], "test": len(truth),
                       "test_resistant": sum(truth), "candidates": 1000}],
            "split": {"scheme": f"first-public-{arm}", "folds": 1, "learners": learners}}


def write(root, des=None, seed=1):
    """Write the whole input set under root and return the design."""
    root = Path(root)
    des = des or design()
    eligible = [d for d in des["datasets"] if d.get("eligible")]
    for sub in ("results", "before", "plan/scripts", "rerun", "logs"):
        (root / sub).mkdir(parents=True, exist_ok=True)
    (root / "design.json").write_text(json.dumps(des, indent=1))
    per_dataset, species = {}, {}
    for k, d in enumerate(eligible):
        n, r = d["test"]["n"], d["test"]["resistant"]
        ids = genome_ids(d, k, n)
        truth = [i < r for i in range(n)]
        for arm in ARMS:
            for folder, run in (("results", "after"), ("before", "before")):
                body = result(d, arm, ids, truth, seed, run)
                (root / folder / f"{d['id']}.first-public-{arm}.json").write_text(json.dumps(body))
        per_dataset[d["id"]] = {"drawn_test_studies": d["test_studies"], "studies_with_a_tested_genome": 2,
                                "tested_genomes": n}
        s = species.setdefault(d["species_slug"], {"datasets": 0, "genomes": 0})
        s["datasets"] += 1
        s["genomes"] += n
    studies = {"datasets": len(eligible), "tested_pairs": sum(d["test"]["n"] for d in eligible), "species": {},
               "per_dataset": per_dataset}
    for slug, s in sorted(species.items()):
        studies["species"][slug] = {"datasets": s["datasets"], "distinct_test_genomes": s["genomes"],
                                    "distinct_test_studies": 3,
                                    "test_genomes_per_study": [s["genomes"] - 3, 2, 1]}
    studies["distinct_test_genomes"] = sum(s["genomes"] for s in species.values())
    studies["sum_of_per_species_study_counts"] = 3 * len(species)
    studies["studies_tested_in_two_species"] = {"SYNTHETIC-STUDY-A": ["alpha_species", "beta_species"]}
    studies["distinct_test_studies"] = 3 * len(species) - 1
    (root / "studies.json").write_text(json.dumps(studies, indent=1))
    plan = b"# Synthetic plan\n"
    analyze = b"# synthetic analysis script\n"
    tables = b"# synthetic layout script\n"
    rerun = b"# Synthetic rerun plan\n"
    (root / "plan" / "PLAN.md").write_bytes(plan)
    (root / "plan" / "scripts" / "analyze.py").write_bytes(analyze)
    (root / "plan" / "scripts" / "tables.py").write_bytes(tables)
    (root / "plan" / "PLAN.sha256").write_text(
        "# sha256 of the plan and of the scripts that apply it, taken at 2026-01-01T10:00:00Z\n"
        f"{sha(plan)}  PLAN.md\n{sha(analyze)}  scripts/analyze.py\n")
    (root / "plan" / "tables.sha256").write_text(
        f"# sha256 of the table layout script, taken at 2026-01-01T10:22:30Z, before any score was read.\n"
        f"{sha(tables)}  scripts/tables.py\n")
    (root / "rerun" / "RERUN-PLAN.md").write_bytes(rerun)
    (root / "rerun" / "RERUN-PLAN.sha256").write_text(
        f"# sha256 of the rerun plan, taken at 2026-01-02T09:00:00Z\n{sha(rerun)}  RERUN-PLAN.md\n")
    first = ["2026-01-01T10:20:00Z binary sha256 0"]
    rerun_log = ["2026-01-02T09:05:00Z START"]
    t = 0
    for d in eligible:
        for arm in ARMS:
            t += 1
            first.append(f"2026-01-01T10:{20 + t // 4:02d}:{(t * 7) % 60:02d}Z {d['id']}.{arm} done in 5s")
            a, b = f"2026-01-02T09:{5 + t // 4:02d}:{(t * 3) % 60:02d}Z", f"2026-01-02T09:{6 + t // 4:02d}:{(t * 3) % 60:02d}Z"
            rerun_log.append(f"{a} {b} rerun-{d['id']}.{arm} rc=0 load=1.0  Elapsed 0:05")
    first.append("2026-01-01T11:00:00Z END rc=0")
    (root / "logs" / "first-run.log").write_text("\n".join(first) + "\n")
    (root / "logs" / "rerun.log").write_text("\n".join(rerun_log) + "\n")
    return des
