"""scan.py: fail when a file names a genome of the species or a BioProject
accession, holds a sequence or a position, writes a rate under its floor or
names a catalogue family or gene carried by under 5 genomes.

    python3 -I scan.py --ids SPECIES.json PATH [PATH ...]

SPECIES.json is the release's species file (its `genomes[].id`). Every file
under each PATH is read. A token of digits, a dot and digits that is a genome
id of the list, or any text of the form PRJ<letters><digits>, is a hit, and so
is a word of 20 or more bases or 30 or more amino acids, a point symbol of the
program (gene_S83L) or a change as TBProfiler writes it (c.-15C>T, p.Ser450Leu).
A bare number is not read as a position.

Each table the mechanism scripts write is checked by its name (TABLES), under
the lead's ruling on question 12 (23 Nov): a genome count is written as it is,
0 to 4 included, or blank; a recall stands only over 10 or more resistant
genomes; an M2 share stands only in a cell whose BioProject vote passes the
plan's entry rule (a voted row by its own `passes`, the pooled row and the
family mix by the BioProject row of their cell); a time-split or M4 share
stands only over 5 or more genomes; a family or gene is named only where 5 or
more genomes carry it, in each period for a family. A .tsv of another name,
other than the run's timing log, is a hit, and so is a known table without a
column the checks read. The output gives counts per file and the rule's name
only, never the hit itself. Exit 1 on any hit.

Python 3.12, standard library only; run with -I.
"""

import argparse
import json
import os
import re
import sys

TOKEN = re.compile(rb"(?<![\w.])\d+\.\d+(?![\w])")
BIOPROJECT = re.compile(rb"PRJ[A-Z]{1,3}\d+")
# 20 or more bases, or 30 or more amino acids, as one word
SEQUENCE = re.compile(rb"(?<![A-Za-z])(?:[ACGTN]{20,}|[acgtn]{20,}|[ACDEFGHIKLMNPQRSTVWY]{30,})(?![A-Za-z])")
# a position: the program's point symbol (gene_S83L, 23S_A2058G, gene_C-42T, gene_Q11Ter) or a change as
# TBProfiler writes it (c.-15C>T, p.Ser450Leu, n.1401A>G); a bare number is not caught
POSITION = re.compile(rb"(?<![A-Za-z0-9])(?:[A-Za-z0-9][A-Za-z0-9()']{1,15}_[A-Z]-?\d{1,6}(?:[A-Z*]|Ter|del|ins|dup|fs)"
                      rb"|[cnpgmr]\.\(?[-*]?(?:[A-Z][a-z]{2})?\d+)")
FLOOR = 5
FOLDED = {"m2-families.tsv": ("family", "other families"),
          "m4-genes.tsv": ("gene", "other genes")}
_M2_SHARES = ("earlier_point", "earlier_gene", "later_point", "later_gene", "d_point", "d_gene", "tvd", "p_point",
              "p_tvd")
_M2_CELL = ("species_slug", "cell", "class", "dataset", "periods")
# name -> genome count columns, {recall: its count}, {share: its counts, each 5 or more}
TABLES = {
    "strata.tsv": (("resistant", "resistant_called", "susceptible", "susceptible_called"),
                   {"recall": "resistant"}, {}),
    "pooled.tsv": (("resistant_pairs", "resistant_called", "susceptible_pairs", "susceptible_called"),
                   {"recall": "resistant_pairs"}, {}),
    "datasets.tsv": (("genomes", "with_table"), {}, {}),
    "nocall.tsv": (("no_call_resistant", "no_call_susceptible"), {}, {}),
    "m2-cells.tsv": (("genomes", "earlier_genomes", "later_genomes", "earlier_units", "later_units",
                      "earlier_point_disrupt", "later_point_disrupt"), {}, {}),
    "m2-families.tsv": (("earlier_genomes", "later_genomes"), {}, {}),
    "m2-years.tsv": (("genomes",), {}, {}),
    "m2-lineage.tsv": ((), {}, {}),
    "timesplit.tsv": (("resistant_later_test", "with_determinant", "with_unseen_determinant", "no_determinant"),
                      {}, {"share_unseen": ("with_determinant",)}),
    "m4-curve.tsv": (("resistant",), {}, {}),
    "m4-coverage.tsv": (("resistant", "without_determinant"), {}, {"share_without": ("resistant",)}),
    "m4-convergence.tsv": (("resistant_covered",), {}, {}),
    "m4-genes.tsv": (("resistant_with_gene",), {}, {}),
    "page-endpoints.tsv": (("resistant_with_result", "resistant_no_result", "without_determinant", "no_main", "no_sub",
                            "no_cluster"), {}, {"share_without": ("resistant_with_result",)}),
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


def scan_text(path):
    """Sequences and positions in one file, counted, never printed."""
    with open(path, "rb") as f:
        data = f.read()
    return len(SEQUENCE.findall(data)), len(POSITION.findall(data))


def is_count(text):
    """A genome count as the tables write it, or blank where none applies."""
    return text in ("", "0") or (text.isascii() and text.isdigit() and text[0] != "0")


def at_least(text, floor):
    return text.isascii() and text.isdigit() and int(text) >= floor


def read_rows(path):
    with open(path, encoding="utf-8", newline="") as f:
        lines = f.read().split("\n")
    if lines and lines[-1] == "":
        lines.pop()
    return lines


def passing_cells(path):
    """The cells of an m2-cells.tsv whose BioProject vote passes, or None."""
    if not os.path.isfile(path):
        return None
    lines = read_rows(path)
    if not lines:
        return None
    header = lines[0].split("\t")
    if not set(_M2_CELL) | {"guard", "passes"} <= set(header):
        return None
    out = set()
    for line in lines[1:]:
        r = dict(zip(header, line.split("\t")))
        if r.get("guard") == "bioproject" and r.get("passes") == "1":
            out.add(tuple(r.get(c, "") for c in _M2_CELL))
    return out


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
    lines = read_rows(path)
    if not lines:
        hit("no header")
        return hits
    header = lines[0].split("\t")
    needed = set(counts) | set(recalls) | set(recalls.values()) | set(shares)
    needed |= {c for cs in shares.values() for c in cs}
    if name == "m2-cells.tsv":
        needed |= {"guard", "passes"} | set(_M2_SHARES) | set(_M2_CELL)
    if name == "m2-families.tsv":
        needed |= {"species_slug", "class", "periods", "earlier_mix", "later_mix"}
    if name in FOLDED:
        needed.add(FOLDED[name][0])
    if needed - set(header):
        hit("a column the scan reads is missing")
        return hits
    rows = []
    for line in lines[1:]:
        cells = line.split("\t")
        if len(cells) != len(header):
            hit("a row of another width")
            continue
        rows.append(dict(zip(header, cells)))
    if name == "m2-cells.tsv":
        passing = {tuple(r[c] for c in _M2_CELL) for r in rows if r["guard"] == "bioproject" and r["passes"] == "1"}
    elif name == "m2-families.tsv":
        passing = passing_cells(os.path.join(os.path.dirname(path), "m2-cells.tsv"))
        if passing is None:
            hit("a family mix with no cells table beside it")
            return hits
    for r in rows:
        for col in counts:
            if not is_count(r[col]):
                hit(f"{col}: not a count")
        for col, of in recalls.items():
            if r[col] != "" and not at_least(r[of], RECALL_FLOOR):
                hit(f"{col}: written over under 10")
        for col, ofs in shares.items():
            if r[col] != "" and not all(at_least(r[o], FLOOR) for o in ofs):
                hit(f"{col}: written over under 5 genomes")
        if name == "m2-cells.tsv" and any(r[col] != "" for col in _M2_SHARES):
            if r["guard"] == "pooled":
                if tuple(r[c] for c in _M2_CELL) not in passing:
                    hit("a pooled share in a cell whose BioProject vote does not pass")
            elif r["passes"] != "1":
                hit("a voted row that does not pass writes a share")
        if name == "m2-families.tsv" and (r["earlier_mix"] != "" or r["later_mix"] != ""):
            if (r["species_slug"], "species_class", r["class"], "", r["periods"]) not in passing:
                hit("a family mix in a cell whose BioProject vote does not pass")
        if name in FOLDED:
            col, label = FOLDED[name]
            ns = [r[c] for c in counts]
            if r[col] != label and (any(n not in ("", "0") and not at_least(n, FLOOR) for n in ns)
                                    or sum(int(n) for n in ns if n.isascii() and n.isdigit()) < FLOOR):
                hit(f"a {col} carried by under 5 genomes is named")
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
        seq, pos = scan_text(path)
        small = check_table(path)
        tables += os.path.basename(path) in TABLES
        if g or b or seq or pos or small:
            hits += 1
            print(f"scan: {path}: {g} genome ids, {b} BioProject accessions, {seq} sequences, {pos} positions")
            for rule in sorted(small):
                print(f"scan: {path}: {small[rule]} rows: {rule}")
    print(f"scan: {files} files, {tables} tables, {hits} with a hit, against {len(ids)} genome ids")
    return 1 if hits else 0


if __name__ == "__main__":
    sys.exit(main())
