"""scan.py: fail when a file names a genome of the species or a BioProject
accession, or when a table could stand for 1 to 4 genomes.

    python3 -I scan.py --ids SPECIES.json PATH [PATH ...]

SPECIES.json is the release's species file (its `genomes[].id`). Every file
under each PATH is read. A token of digits, a dot and digits that is a genome
id of the list, or any text of the form PRJ<letters><digits>, is a hit.

Each table the mechanism scripts write is checked by its name (TABLES): a
genome count is 0, "under 5", 5 and more or blank; a recall stands only over 10 or
more resistant genomes; a share is blank where its count is "under 5"; a
voted M2 row that does not pass writes no share; a family or gene in under 5
genomes is not named. A .tsv of another name, other than the run's timing
log, is a hit, and so is a known table without a column the checks read.
The output gives counts per file and the rule's name only, never the hit
itself. Exit 1 on any hit.

Python 3.12, standard library only; run with -I.
"""

import argparse
import json
import os
import re
import sys

TOKEN = re.compile(rb"(?<![\w.])\d+\.\d+(?![\w])")
BIOPROJECT = re.compile(rb"PRJ[A-Z]{1,3}\d+")
UNDER = "under 5"
FOLDED = {"m2-families.tsv": ("family", "families seen in under 5 genomes"),
          "m4-genes.tsv": ("gene", "genes seen in under 5 genomes")}
_M2_SHARES = {"earlier_point": ("earlier_genomes",), "earlier_gene": ("earlier_genomes",),
              "later_point": ("later_genomes",), "later_gene": ("later_genomes",),
              "d_point": ("earlier_genomes", "later_genomes"), "d_gene": ("earlier_genomes", "later_genomes"),
              "tvd": ("earlier_genomes", "later_genomes"), "p_point": ("earlier_genomes", "later_genomes"),
              "p_tvd": ("earlier_genomes", "later_genomes")}
# name -> genome count columns, {recall: its count}, {share: its counts}
TABLES = {
    "strata.tsv": (("resistant", "resistant_called", "susceptible", "susceptible_called"),
                   {"recall": "resistant"}, {}),
    "pooled.tsv": (("resistant_pairs", "resistant_called", "susceptible_pairs", "susceptible_called"),
                   {"recall": "resistant_pairs"}, {}),
    "datasets.tsv": (("genomes", "with_table"), {}, {}),
    "nocall.tsv": (("no_call_resistant", "no_call_susceptible"), {}, {}),
    "m2-cells.tsv": (("genomes", "earlier_genomes", "later_genomes", "earlier_units", "later_units",
                      "earlier_point_disrupt", "later_point_disrupt"), {}, _M2_SHARES),
    "m2-families.tsv": (("earlier_genomes", "later_genomes"), {},
                        {"earlier_mix": ("earlier_genomes",), "later_mix": ("later_genomes",)}),
    "m2-years.tsv": (("genomes",), {}, {}),
    "m2-lineage.tsv": ((), {}, {}),
    "timesplit.tsv": (("resistant_later_test", "with_determinant", "with_unseen_determinant", "no_determinant"),
                      {}, {"share_unseen": ("with_determinant",)}),
    "m4-curve.tsv": (("resistant",), {}, {}),
    "m4-coverage.tsv": (("resistant", "without_determinant"), {}, {"share_without": ("resistant",)}),
    "m4-convergence.tsv": (("resistant_covered",), {}, {}),
    "m4-genes.tsv": (("resistant_with_gene",), {}, {}),
}
NOT_TABLES = ("timing.tsv",)
RECALL_FLOOR = 10


def load_ids(path):
    with open(path, encoding="utf-8") as f:
        return {g["id"].encode() for g in json.load(f)["genomes"]}


def scan_file(path, ids):
    with open(path, "rb") as f:
        data = f.read()
    genomes = sum(1 for m in TOKEN.finditer(data) if m.group(0) in ids)
    projects = len(BIOPROJECT.findall(data))
    return genomes, projects


def is_count(text):
    """A genome count as the tables write it, or blank where none applies."""
    return text in ("", "0", UNDER) or (text.isascii() and text.isdigit() and int(text) >= 5 and text[0] != "0")


def check_table(path):
    """Rule name -> hits for one table, by its file name."""
    name = os.path.basename(path)
    hits = {}

    def hit(rule):
        hits[rule] = hits.get(rule, 0) + 1

    if not name.endswith(".tsv") or name in NOT_TABLES:
        return hits
    if name not in TABLES:
        hit("a table the scan does not know")
        return hits
    counts, recalls, shares = TABLES[name]
    with open(path, encoding="utf-8", newline="") as f:
        lines = f.read().split("\n")
    if lines and lines[-1] == "":
        lines.pop()
    if not lines:
        hit("no header")
        return hits
    header = lines[0].split("\t")
    needed = set(counts) | set(recalls) | set(recalls.values()) | set(shares)
    needed |= {c for cs in shares.values() for c in cs}
    if name == "m2-cells.tsv":
        needed |= {"guard", "passes"}
    if name in FOLDED:
        needed.add(FOLDED[name][0])
    if needed - set(header):
        hit("a column the scan reads is missing")
        return hits
    for line in lines[1:]:
        cells = line.split("\t")
        if len(cells) != len(header):
            hit("a row of another width")
            continue
        r = dict(zip(header, cells))
        for col in counts:
            if not is_count(r[col]):
                hit(f"{col}: not 0, under 5 or 5 and more")
        for col, of in recalls.items():
            n = r[of]
            if r[col] != "" and not (n.isascii() and n.isdigit() and int(n) >= RECALL_FLOOR):
                hit(f"{col}: written over under 10")
        for col, ofs in shares.items():
            if r[col] != "" and any(r[o] == UNDER for o in ofs):
                hit(f"{col}: written beside a count under 5")
        if name == "m2-cells.tsv" and r["guard"] != "pooled" and r["passes"] != "1":
            if any(r[col] != "" for col in _M2_SHARES):
                hit("a voted row that does not pass writes a share")
        if name in FOLDED:
            col, label = FOLDED[name]
            if r[col] != label and any(r[c] == UNDER for c in counts):
                hit(f"a {col} under 5 genomes is named")
    return hits


def files_under(paths):
    for p in paths:
        if os.path.isfile(p):
            yield p
            continue
        for root, _, names in os.walk(p):
            for n in sorted(names):
                yield os.path.join(root, n)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--ids", required=True)
    ap.add_argument("paths", nargs="+")
    args = ap.parse_args(argv)
    ids = load_ids(args.ids)
    files = hits = 0
    tables = 0
    for path in files_under(args.paths):
        files += 1
        g, b = scan_file(path, ids)
        small = check_table(path)
        tables += os.path.basename(path) in TABLES
        if g or b or small:
            hits += 1
            print(f"scan: {path}: {g} genome ids, {b} BioProject accessions")
            for rule in sorted(small):
                print(f"scan: {path}: {small[rule]} rows: {rule}")
    print(f"scan: {files} files, {tables} tables, {hits} with a hit, against {len(ids)} genome ids")
    return 1 if hits else 0


if __name__ == "__main__":
    sys.exit(main())
