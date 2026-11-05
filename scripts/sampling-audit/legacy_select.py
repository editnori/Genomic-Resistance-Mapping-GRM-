#!/usr/bin/env python3
"""A Python port of `select_genomes` in src/curate.rs (unchanged from
7711459 through c90ebaf), for read-only audits of how curation picks the
genomes each species downloads.

Input contract. The port reads what `grm curate` writes and what
`select_genomes` is handed: `curated/phenotypes.jsonl` (one consolidated
genome x antibiotic row after QC) and `curated/datasets.json` (the species x
antibiotic pairs that passed `--min-class`). Of each phenotype row it uses
genome_id, species, antibiotic, reported, conflict, label_from_mic and label;
of each dataset, species and antibiotic. Everything that decides which
genomes reach those files happens earlier, in `curate`, and is not redone
here (see README.md, "What precedes selection").

The Rust function, step by step, as ported:
  1. targets: the antibiotics with a dataset for the species.
  2. For each phenotype row of the species on a target antibiotic, in file
     order, take `reported_label()`: the reported call when it is R or S,
     the pair is not conflicted and the label was not inferred from MIC.
     Rows without one (I, conflicted, MIC-inferred, unlabelled) are skipped.
  3. class_sizes[(antibiotic, call)] counts those rows; each genome keeps
     its (antibiotic, call) list in file order.
  4. A genome's score is the left-to-right f64 sum of 1 / class_size over
     its list. Python's built-in sum() compensates float error from 3.12 on,
     so the port adds in a loop to reproduce Rust's `Iterator::sum`.
  5. Sort by score descending (`f64::total_cmp`), then genome id ascending
     by bytes (`str::cmp`), and keep the first `cap`.
A genome with no reported R/S label on a target antibiotic is never scored,
so it is not picked at any cap, however large.

usage: python3 -I legacy_select.py CURATED_DIR SPECIES CAP
prints the selected ids, one per line, as genomes.txt holds them.
"""
import json
import os
import sys

R, S = "R", "S"


def reported_label(row):
    """`Phenotype::reported_label`: the unconflicted, not MIC-inferred R or S call."""
    call = row.get("reported")
    if call in (R, S) and not row.get("conflict", False) and not row.get("label_from_mic", False):
        return call
    return None


def slug(species):
    """`curate::slug`: lower case, runs of non-alphanumerics to one underscore."""
    s = "".join(c if c.isascii() and c.isalnum() else "_" for c in species.lower())
    return "_".join(part for part in s.split("_") if part)


def targets(datasets, species):
    return {d["antibiotic"] for d in datasets if d["species"] == species}


def labels_by_genome(phenotypes, datasets, species):
    """Steps 1-3: (class_sizes, {genome_id: [(antibiotic, call), ...]})."""
    drugs = targets(datasets, species)
    class_sizes = {}
    labels = {}
    for p in phenotypes:
        if p["species"] != species or p["antibiotic"] not in drugs:
            continue
        call = reported_label(p)
        if call is None:
            continue
        key = (p["antibiotic"], call)
        class_sizes[key] = class_sizes.get(key, 0) + 1
        labels.setdefault(p["genome_id"], []).append(key)
    return class_sizes, labels


def score(label_list, class_sizes):
    """Step 4: Rust's sequential f64 sum, not Python 3.12's compensated sum()."""
    total = 0.0
    for key in label_list:
        total += 1.0 / class_sizes[key]
    return total


def ranking(phenotypes, datasets, species):
    """Every scored genome in selection order, as [(score, genome_id)]."""
    class_sizes, labels = labels_by_genome(phenotypes, datasets, species)
    scored = [(score(ls, class_sizes), gid) for gid, ls in labels.items()]
    # Scores are finite and positive, so negation orders them as
    # b.total_cmp(a) does; ids compare as Rust's str::cmp, by UTF-8 bytes.
    scored.sort(key=lambda t: (-t[0], t[1].encode("utf-8")))
    return scored


def select_genomes(phenotypes, datasets, species, cap):
    return [gid for _, gid in ranking(phenotypes, datasets, species)[:cap]]


PHENOTYPE_FIELDS = ("genome_id", "species", "antibiotic", "reported", "conflict", "label_from_mic", "label")


def read_phenotypes(path, species=None):
    """phenotypes.jsonl rows, trimmed to the fields the audit uses. With
    `species` (a set), rows of other species are dropped while reading."""
    rows = []
    intern = sys.intern
    with open(path, encoding="utf-8") as f:
        for line in f:
            if not line.strip():
                continue
            p = json.loads(line)
            if species is not None and p["species"] not in species:
                continue
            rows.append({
                "genome_id": intern(p["genome_id"]),
                "species": intern(p["species"]),
                "antibiotic": intern(p["antibiotic"]),
                "reported": p.get("reported"),
                "conflict": bool(p.get("conflict", False)),
                "label_from_mic": bool(p.get("label_from_mic", False)),
                "label": p.get("label"),
            })
    return rows


def read_datasets(path):
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def main(argv):
    if len(argv) != 4:
        sys.exit(__doc__.rsplit("usage:", 1)[1].strip())
    curated, species, cap = argv[1], argv[2], int(argv[3])
    phenotypes = read_phenotypes(os.path.join(curated, "phenotypes.jsonl"), {species})
    datasets = read_datasets(os.path.join(curated, "datasets.json"))
    for gid in select_genomes(phenotypes, datasets, species, cap):
        print(gid)


if __name__ == "__main__":
    main(sys.argv)
