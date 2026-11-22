"""Shared loaders and SYNTHETIC fixtures for the small-cohort scripts (scripts/c1-cutoff, docs/cohort-draw.md).

Every sample, report and result here is made up in the test. No genome or
call comes from a real dataset; the ids are synthetic numbers in the BV-BRC
id shape and the dataset names are made up.
"""

import importlib.util
import json
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[2]
KIT = ROOT / "scripts" / "c1-cutoff"
DATASET = "synthia_coli__drug_a"
LEARNERS = ["scm", "cart", "rf", "gbt", "lasso"]
CONFIG = {"folds": 5, "inner_folds": 3, "seed": 42, "screened_features": 4000}


def load(name):
    spec = importlib.util.spec_from_file_location(f"c1_cutoff_{name}", KIT / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def samples(resistant=60, susceptible=80, lineages=30, seed=1):
    """SYNTHETIC samples: ids 9000.<n>, calls, lineage numbers spread unevenly."""
    out = []
    state = seed
    for i in range(resistant + susceptible):
        state = (state * 6364136223846793005 + 1442695040888963407) % (1 << 64)
        out.append({"id": f"9000.{i + 1}", "resistant": i < resistant, "lineage": (state >> 33) % lineages})
    return out


def report(smp, dataset=DATASET, fingerprint="fp-synthetic", lineage_ran=True, n_groups=30):
    r = sum(1 for s in smp if s["resistant"])
    return {
        "id": dataset, "species_slug": dataset.split("__")[0], "species": "Synthia coli",
        "antibiotic": dataset.split("__")[1].replace("_", " "), "input_fingerprint": fingerprint,
        "eval_config": dict(CONFIG), "counts": {"samples": len(smp), "resistant": r, "susceptible": len(smp) - r},
        "splits": [{"scheme": "random", "n_groups": len(smp)},
                   {"scheme": "lineage", "n_groups": n_groups, "unavailable_reason": None if lineage_ran else "too few lineages"}],
        "samples": smp,
    }


def write_release(root, reports):
    (root / "datasets").mkdir(parents=True, exist_ok=True)
    for rep in reports:
        (root / "datasets" / f"{rep['id']}.json").write_text(json.dumps(rep))


def calls(n_pos, n_neg, tp, tn, loss=1.0):
    return {"n": n_pos + n_neg, "tp": tp, "fn": n_pos - tp, "tn": tn, "fp": n_neg - tn, "loss_sum": loss}


def result(r=25, draw=0, seed=42, dataset=DATASET, samples_=140, resistant=60, copies=(2, 1), new=(30, 10),
           study=False, cohort_sha="c" * 64, fingerprint="fp-synthetic", truth_from_folds=False,
           lineage_unavailable=None, est_tp=20, est_tn=21, folds=(0.7, 0.8, 0.75, 0.9, 0.85), ba=0.82):
    """A SYNTHETIC result of grm learn --exp-cohort with consistent counts."""
    scored = (samples_ - 2 * r - copies[0], resistant - r - copies[1])
    est = {"learner": None, "calls": calls(r, r, est_tp, est_tn), "balanced_accuracy": ba,
           "fold_balanced_accuracy": list(folds), "no_call": None}

    def split(scheme, unavailable=None):
        return {"scheme": scheme, "n_groups": 20, "unavailable": unavailable,
                "learners": [] if unavailable else [dict(est, learner=name) for name in LEARNERS]}

    estimates = [split("random"), split("lineage", lineage_unavailable)]
    if study:
        estimates.append(split("study-whole"))
    if truth_from_folds:
        all_ = calls(r, r, est_tp, est_tn)
    else:
        all_ = calls(scored[1], scored[0] - scored[1], scored[1] - 5, scored[0] - scored[1] - 10, loss=12.5)
    nl = calls(new[1], new[0] - new[1], new[1] - 2, new[0] - new[1] - 4, loss=6.0)
    return {
        "format": "grm-cohort-draw-result-v1", "evaluation_revision": 14, "release_fingerprint": fingerprint,
        "eval_config": dict(CONFIG, seed=seed), "id": dataset, "species_slug": dataset.split("__")[0],
        "antibiotic": "drug a", "r": r, "draw": draw, "seed": seed, "cohort_file_sha256": cohort_sha,
        "dataset": {"samples": samples_, "resistant": resistant}, "refused": None,
        "cohort": {"samples": 2 * r, "resistant": r, "lineages": 20, "copy_groups": 0, "min_support": 3, "features": 100},
        "scored": {"samples": scored[0], "resistant": scored[1]},
        "new_lineage": {"samples": new[0], "resistant": new[1]},
        "copies_left_out": {"samples": copies[0], "resistant": copies[1]},
        "estimates": estimates,
        "truth": {"unavailable": None, "learners": [
            {"learner": name, "null_fit": False, "all": dict(all_), "new_lineage": dict(nl)} for name in LEARNERS]},
        "seconds": 1.0,
    }


def expect_row(r=25, draw=0, seed=42, dataset=DATASET, study="0", cohort_sha="c" * 64):
    return {"dataset": dataset, "species_slug": dataset.split("__")[0], "r": str(r), "draw": str(draw),
            "seed": str(seed), "study": study, "cohort": f"{dataset}.r{r}.d{draw}.cohort.json",
            "cohort_sha256": cohort_sha, "file": f"{dataset}.r{r}.d{draw}.s{seed}.cohort.json"}
