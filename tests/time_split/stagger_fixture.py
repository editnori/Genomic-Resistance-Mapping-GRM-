"""SYNTHETIC inputs for the staggered split's scripts: a data root, a fetch directory, a built root, result files
and lists, all made up in the test.

Every genome id here starts with `syn-` and every study with `syn-study-`, so no real BV-BRC genome id or BioProject
accession can equal one. Days, years, labels, digests and scores are made up too.
"""

import contextlib
import datetime
import hashlib
import importlib.util
import io
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT / "scripts" / "time-split"
DESIGN = SCRIPTS / "stagger_design.py"
ANALYZE = SCRIPTS / "stagger_analyze.py"
FAULTS = SCRIPTS / "stagger_faults.py"
PAGE = ROOT / "docs" / "time-split-stagger.md"
LEARNERS = ("scm", "cart", "rf", "gbt", "lasso")


def load(path, name):
    """A script as a module, loaded once."""
    module = sys.modules.get(name)
    if module is None:
        spec = importlib.util.spec_from_file_location(name, path)
        module = importlib.util.module_from_spec(spec)
        sys.modules[name] = module
        spec.loader.exec_module(module)
    return module


def design_script():
    return load(DESIGN, "stagger_design_under_test")


def analyze_script():
    return load(ANALYZE, "stagger_analyze_under_test")


def run(module, *argv):
    """main() with its output caught: (exit code, stdout, stderr)."""
    out, err = io.StringIO(), io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        code = module.main([str(a) for a in argv])
    return code, out.getvalue(), err.getvalue()


def sha(text):
    return hashlib.sha256(text.encode()).hexdigest()


# ----------------------------------------------------------------- a made-up world on disk

SPECIES = ("alpha_species", "beta_species", "gamma_species")
DRUGS = ("drug_one", "drug_two", "drug_three", "drug_four", "drug_five")


def gid(slug, i):
    return f"syn-{slug[0]}-{i:04d}"


def genome_rows(slug, n=400):
    """One species' genomes: day, year, study, digest and whether the built matrix keeps it."""
    start = datetime.date(2015, 1, 1)
    rows = []
    for i in range(n):
        g = gid(slug, i)
        rows.append({
            "id": g,
            "day": (start + datetime.timedelta(days=i // 2)).isoformat() + "T08:00:00.000Z",
            "year": None if i % 13 == 7 else 2010 + i * 12 // n,
            "study": None if i % 17 == 5 else f"syn-study-{slug}-{i // 20:03d}",
            "digest": f"syn-digest-{slug}-{i - 1 if i % 40 == 1 else i}",
            "in_matrix": i % 31 != 0,
        })
    return rows


def label(slug, drug, i):
    """R, S or None: made-up calls. drug_five is mostly susceptible; drug_two leaves some genomes unlabelled."""
    if drug == "drug_two" and i % 11 == 3:
        return None
    h = sha(f"{slug}|{drug}|{i}")
    if drug == "drug_five":
        return "R" if h[0] == "0" else "S"
    return "R" if h[0] in "01234567" else "S"


def write_custom(base, genomes, datasets):
    """A data root, a fetch directory and a built root under base from explicit rows.

    genomes: dicts with id, slug, day (text or None), year, study, digest and in_matrix (False keeps a genome out
    of its built matrix; None for every genome of a species writes no matrix for it). datasets: dicts with id,
    species_slug, antibiotic and labels ({genome id: R or S}). Returns (root, fetch, built root)."""
    base = Path(base)
    root, fetch, built_root = base / "root", base / "fetch", base / "built"
    for p in (root / "curated", root / "bvbrc", fetch, built_root):
        p.mkdir(parents=True, exist_ok=True)
    slugs = list(dict.fromkeys(g["slug"] for g in genomes))
    for slug in slugs:
        rows = [g for g in genomes if g["slug"] == slug]
        (root / "species" / slug).mkdir(parents=True, exist_ok=True)
        (root / "species" / slug / "genomes.txt").write_text("".join(r["id"] + "\n" for r in rows))
        if any(r["in_matrix"] is not None for r in rows):
            kept = [r for r in rows if r["in_matrix"]]
            matrix = built_root / "species" / slug / "matrix"
            matrix.mkdir(parents=True, exist_ok=True)
            (matrix / "matrix.json").write_text(json.dumps(
                {"genomes": [r["id"] for r in kept], "kset_digests": [r["digest"] for r in kept]}))
    specs, phenotypes = [], []
    for d in datasets:
        name = d["species_slug"].replace("_", " ").capitalize()
        specs.append({"id": d["id"], "species": name, "species_slug": d["species_slug"],
                      "antibiotic": d["antibiotic"]})
        for g, call in d["labels"].items():
            phenotypes.append({"genome_id": g, "species": name, "antibiotic": d["antibiotic"], "reported": call,
                               "conflict": False, "label_from_mic": False})
    (root / "curated" / "datasets.json").write_text(json.dumps(specs))
    (root / "curated" / "phenotypes.jsonl").write_text("".join(json.dumps(p) + "\n" for p in phenotypes))
    (root / "bvbrc" / "genomes.jsonl").write_text(
        "".join(json.dumps({"genome_id": g["id"], "collection_year": g["year"]}) + "\n" for g in genomes))
    (fetch / "batch-0001.json").write_text(json.dumps(
        [{"genome_id": g["id"], "date_inserted": g["day"], "bioproject_accession": g["study"]} for g in genomes]))
    return root, fetch, built_root


def write_world(base, n=400, built=True, species=SPECIES, drugs=DRUGS):
    """The made-up world of 3 species and 5 drugs. Returns (root, fetch, built root, labels)."""
    genomes, datasets, labels = [], [], {}
    for slug in species:
        rows = genome_rows(slug, n)
        for r in rows:
            genomes.append({"id": r["id"], "slug": slug, "day": r["day"], "year": r["year"], "study": r["study"],
                            "digest": r["digest"], "in_matrix": r["in_matrix"] if built else None})
        for drug in drugs:
            calls = {r["id"]: label(slug, drug, i) for i, r in enumerate(rows) if label(slug, drug, i)}
            datasets.append({"id": f"{slug}__{drug}", "species_slug": slug, "antibiotic": drug, "labels": calls})
            labels.update({(g, drug): c for g, c in calls.items()})
    root, fetch, built_root = write_custom(base, genomes, datasets)
    return root, fetch, built_root, labels


# ----------------------------------------------------------------- result files

def metrics(truth, scores):
    tp = sum(1 for t, s in zip(truth, scores) if t and s >= 0.5)
    fn = sum(1 for t, s in zip(truth, scores) if t and s < 0.5)
    tn = sum(1 for t, s in zip(truth, scores) if not t and s < 0.5)
    fp = sum(1 for t, s in zip(truth, scores) if not t and s >= 0.5)
    return {"n": len(truth), "tp": tp, "fn": fn, "tn": tn, "fp": fp, "sensitivity": tp / (tp + fn),
            "specificity": tn / (tn + fp)}


def result(dataset, scheme, sha256, samples, scores_by_learner, train, train_resistant):
    """A result file in the shape `grm learn --fold-file` writes (docs/fold-file.md)."""
    truth = [t for _, t in samples]
    learners = [{"learner": l, "metrics": metrics(truth, scores_by_learner[l]), "scores": scores_by_learner[l],
                 "fold_complexity": [3]} for l in LEARNERS]
    return {"format": "grm-fold-file-result-v1", "evaluation_revision": 14, "fold_file_sha256": sha256,
            "id": dataset, "species_slug": dataset.split("__")[0], "scheme": scheme,
            "eval_config": {"seed": 42, "screened_features": 4000, "inner_folds": 3},
            "samples": [{"genome_id": g, "resistant": t, "fold": 0} for g, t in samples],
            "folds": [{"train": train, "train_resistant": train_resistant, "test": len(truth),
                       "test_resistant": sum(truth), "candidates": 1000}],
            "split": {"scheme": scheme, "folds": 1, "learners": learners}}


def made_up_score(scheme, genome, learner, truth, origin, window):
    """A made-up score: an older origin is a little further from the truth."""
    u = int(sha(f"{scheme}|{genome}|{learner}")[:8], 16) / 2 ** 32
    lean = 0.08 * (window - origin)
    p = 0.55 * u + (0.4 - lean if truth else 0.05 + lean)
    return 1.0 if learner == "scm" and p >= 0.5 else 0.0 if learner == "scm" else min(max(p, 0.0), 1.0)


def write_results(design_dir, result_dir, labels):
    """A result file for every fold file of a design, from the fold files themselves."""
    design = json.loads((Path(design_dir) / "design.json").read_text())
    result_dir = Path(result_dir)
    result_dir.mkdir(parents=True, exist_ok=True)
    for d in design["datasets"]:
        if d["why_not"]:
            continue
        drug = d["antibiotic"]
        for w in d["windows"]:
            for f in w["fits"]:
                if f.get("same_file_as"):
                    continue
                body = json.loads((Path(design_dir) / "fold-files" / f["file"]).read_text())
                test = body["test"][0]
                samples = [(g, labels[g, drug] == "R") for g in test]
                scores = {l: [made_up_score(f["scheme"], g, l, t, f["origin"], w["window"]) for g, t in samples]
                          for l in LEARNERS}
                train_r = sum(labels[g, drug] == "R" for g in body["train"])
                (result_dir / f"{d['id']}.{f['scheme']}.json").write_text(json.dumps(
                    result(d["id"], f["scheme"], f["sha256"], samples, scores, len(body["train"]), train_r)))
    return design
