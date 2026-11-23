"""A wider arm's study fold files (docs/cap-arms.md, The held-out genomes, study split).

  python3 -I study.py <reviewed folds dir> <study.txt> <cohorts dir> <arm> <curated dir>
                      <genome meta .jsonl.gz[,...]> <structure.json> <out dir> <out tsv>

For each dataset of <study.txt>, the reviewed plan
<reviewed folds dir>/<dataset>.study-whole.json is the default's fold file
(grm-fold-file-v1). The arm's file keeps its test folds, train and unused
lists and adds the arm's added genomes of the dataset's cohort (in the
arm's cohort list and not the default's, with a reported R or S call on the
dataset's drug, not a conflict and not from a MIC, as cohorts.py reads
them). Each added genome goes to train when its BioProject (from the genome
row files, a comma list: the master genome rows and bioprojects.py's
answers for the added genomes, a genome two files disagree on unknown) is
known and no tested genome has it, and no tested genome has
its k-mer digest (from the species grouping); else to unused. An added
genome is never tested: a fold file cannot hold a genome out of one fold
without testing it there, so it trains in every fold or in none. When no
added genome trains, no file is written and the dataset is not fitted for
the arm.

Writes <out dir>/<dataset>.study-whole.json per dataset with a training
genome added, and <out tsv>: dataset, arm, added, train_added,
unused_bioproject, unused_unknown, unused_duplicate, written. Prints the
totals.
"""

import importlib.util
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
COLUMNS = ["dataset", "arm", "added", "train_added", "unused_bioproject", "unused_unknown", "unused_duplicate",
           "written"]


def _load(name):
    spec = importlib.util.spec_from_file_location(f"c2_{name}", HERE / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


COHORTS = _load("cohorts")
FOLDS = _load("check_folds")


def labelled(phenotypes, species, antibiotic):
    """The genomes a dataset's cohort can hold: a reported R or S call on its drug, not a conflict, not from a MIC."""
    return {p["genome_id"] for p in phenotypes
            if p["species"] == species and p["antibiotic"] == antibiotic and COHORTS.reported(p) is not None}


def place(f0, added, bp, digest):
    """(train, unused by reason) for the added genomes against the default's file f0."""
    tested = set().union(*map(set, f0["test"])) if f0["test"] else set()
    test_bp = {bp[g] for g in tested if g in bp}
    test_digest = {digest[g] for g in tested if g in digest}
    train, unused = [], {"bioproject": [], "unknown": [], "duplicate": []}
    for g in sorted(added):
        if g not in bp:
            unused["unknown"].append(g)
        elif bp[g] in test_bp:
            unused["bioproject"].append(g)
        elif digest.get(g) in test_digest:
            unused["duplicate"].append(g)
        else:
            train.append(g)
    return train, unused


def fold_file(f0, train, unused):
    return {"format": f0["format"], "dataset": f0["dataset"], "scheme": f0["scheme"], "test": f0["test"],
            "train": list(f0["train"]) + train,
            "unused": list(f0["unused"]) + sorted(unused["bioproject"] + unused["unknown"] + unused["duplicate"])}


def main(argv):
    if len(argv) != 10:
        raise SystemExit(__doc__)
    folds, study_txt, cohorts, arm, curated, meta, structure, out, out_tsv = argv[1:]
    study = FOLDS.read_ids(study_txt)
    base = set(FOLDS.read_ids(Path(cohorts) / "cap-2000.txt"))
    wide = set(FOLDS.read_ids(Path(cohorts) / f"{arm}.txt"))
    specs = {d["id"]: d for d in FOLDS.load_json(Path(curated) / "datasets.json")}
    phenotypes = COHORTS.read_jsonl(Path(curated) / "phenotypes.jsonl") if study else []
    bp = FOLDS.bioprojects(meta)
    _, digest = FOLDS.grouping(structure)
    Path(out).mkdir(parents=True, exist_ok=True)
    rows = []
    for d in study:
        f0 = FOLDS.load_json(Path(folds) / f"{d}.study-whole.json")
        spec = specs[d]
        added = (wide - base) & labelled(phenotypes, spec["species"], spec["antibiotic"])
        train, unused = place(f0, added, bp, digest)
        if train:
            (Path(out) / f"{d}.study-whole.json").write_text(json.dumps(fold_file(f0, train, unused)) + "\n",
                                                             encoding="utf-8")
        rows.append([d, arm, len(added), len(train), len(unused["bioproject"]), len(unused["unknown"]),
                     len(unused["duplicate"]), int(bool(train))])
    FOLDS.write(out_tsv, COLUMNS, rows)
    print(f"{arm}: study datasets {len(rows)}, fold files written {sum(r[-1] for r in rows)}, "
          f"added genomes {sum(r[2] for r in rows)}, in train {sum(r[3] for r in rows)}")


if __name__ == "__main__":
    main(sys.argv)
