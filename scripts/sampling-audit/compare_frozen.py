#!/usr/bin/env python3
"""Compares the ported selection with a frozen `grm curate` run: for each
species with a dataset, the port's selection from RUN/curated/phenotypes.jsonl
and RUN/curated/datasets.json at CAP against RUN/species/<slug>/genomes.txt,
id for id and in order. Reads RUN only and writes only OUT.json.

usage: python3 -I compare_frozen.py OUT.json RUN CAP [--sums SHA256SUMS]
--sums checks each genomes.txt against a sha256sum listing whose paths are
relative to the listing's directory (./RUN-NAME/species/<slug>/genomes.txt).
"""
import argparse
import hashlib
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import legacy_select as ls  # noqa: E402


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_sums(path):
    sums = {}
    base = os.path.dirname(os.path.abspath(path))
    with open(path, encoding="utf-8") as f:
        for line in f:
            digest, _, name = line.rstrip("\n").partition("  ")
            sums[os.path.normpath(os.path.join(base, name))] = digest
    return sums


def compare(run, cap, sums=None):
    curated = os.path.join(run, "curated")
    phen_path = os.path.join(curated, "phenotypes.jsonl")
    ds_path = os.path.join(curated, "datasets.json")
    datasets = ls.read_datasets(ds_path)
    phenotypes = ls.read_phenotypes(phen_path)
    out = {
        "run": os.path.abspath(run),
        "cap": cap,
        "inputs": {"phenotypes.jsonl": sha256(phen_path), "datasets.json": sha256(ds_path)},
        "species": [],
    }
    for species in sorted({d["species"] for d in datasets}):
        path = os.path.join(run, "species", ls.slug(species), "genomes.txt")
        with open(path, encoding="utf-8") as f:
            frozen = f.read().split("\n")
        if frozen and frozen[-1] == "":
            frozen.pop()
        ported = ls.select_genomes(phenotypes, datasets, species, cap)
        first_diff = next((i for i, (a, b) in enumerate(zip(ported, frozen)) if a != b), None)
        if first_diff is None and len(ported) != len(frozen):
            first_diff = min(len(ported), len(frozen))
        row = {
            "species": species,
            "frozen": len(frozen),
            "ported": len(ported),
            "identical": ported == frozen,
            "same_set": set(ported) == set(frozen),
            "first_difference": first_diff,
            "genomes_txt_sha256": sha256(path),
        }
        if sums is not None:
            row["sha256_matches_listing"] = sums.get(os.path.normpath(os.path.abspath(path))) == row["genomes_txt_sha256"]
        out["species"].append(row)
    out["identical_species"] = sum(r["identical"] for r in out["species"])
    out["total_species"] = len(out["species"])
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("out")
    ap.add_argument("run")
    ap.add_argument("cap", type=int)
    ap.add_argument("--sums")
    a = ap.parse_args(argv)
    result = compare(a.run, a.cap, read_sums(a.sums) if a.sums else None)
    with open(a.out, "w", encoding="utf-8") as f:
        json.dump(result, f, indent=1)
        f.write("\n")
    print(f"{a.run}: {result['identical_species']}/{result['total_species']} species identical")


if __name__ == "__main__":
    main()
