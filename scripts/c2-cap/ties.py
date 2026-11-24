"""Count, per species, the genomes that share the selection rule's score at the 4,096th place (docs/cap-arms-result.md, Limits).

  python3 -I ties.py <curated dir> <out tsv>

<curated dir> holds the rebuilt revision 14 table's species.json,
datasets.json and phenotypes.jsonl, as cohorts.py reads them. For each
species of species.tsv, it scores every genome as cohorts.py's rank() does
(`curate::select_genomes`: the sum of 1 / (size of the call's class) over a
genome's reported R or S calls, summed in the phenotypes' order), holds its
order to rank()'s, and takes the score of the 4,096th genome. A genome ties
at the cut when its score equals that score exactly, as the rule's sort
compares them; genomes of one score sort by id, so the cut keeps the tied
genomes of lower id.

Writes <out tsv>, one row per species, counts only: species, slug, ranked,
tied (the genomes sharing the score at the 4,096th place), kept (those of
them inside the top 4,096) and left (those past it). A species whose rule
ranks 4,096 genomes or fewer has no cut: tied, kept and left are empty and
cut is "no cut". No genome id is written or printed.
"""

import importlib.util
import json
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
CUT = 4096
COLUMNS = ["species", "species_slug", "ranked", "cut", "tied", "kept", "left"]


def _cohorts():
    spec = importlib.util.spec_from_file_location("c2_cohorts_for_ties", HERE / "cohorts.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


cohorts = _cohorts()


def scores(phenotypes, datasets, species):
    """{genome: score} as cohorts.rank() sums it."""
    targets = {d["antibiotic"] for d in datasets if d["species"] == species}
    sizes = defaultdict(int)
    labels = defaultdict(list)
    for p in phenotypes:
        if p["species"] != species or p["antibiotic"] not in targets:
            continue
        call = cohorts.reported(p)
        if call is None:
            continue
        sizes[(p["antibiotic"], call)] += 1
        labels[p["genome_id"]].append((p["antibiotic"], call))
    out = {}
    for g, ls in labels.items():
        t = 0.0
        for a, c in ls:
            t += 1.0 / sizes[(a, c)]
        out[g] = t
    return out


def tie_row(phenotypes, datasets, species, slug):
    score = scores(phenotypes, datasets, species)
    order = cohorts.rank(phenotypes, datasets, species)
    if order != sorted(score, key=lambda g: (-score[g], g)):
        raise SystemExit(f"{slug}: the scores do not give cohorts.py's ranking")
    row = {"species": species, "species_slug": slug, "ranked": len(order)}
    if len(order) <= CUT:
        return dict(row, cut="no cut", tied="", kept="", left="")
    at = score[order[CUT - 1]]
    kept = sum(1 for g in order[:CUT] if score[g] == at)
    left = sum(1 for g in order[CUT:] if score[g] == at)
    return dict(row, cut=str(CUT), tied=kept + left, kept=kept, left=left)


def main(argv):
    if len(argv) != 2:
        sys.stderr.write(__doc__)
        return 2
    curated, out = Path(argv[0]), Path(argv[1])
    slugs = [line.split("\t")[0] for line in (HERE / "species.tsv").read_text(encoding="utf-8").splitlines()[1:]
             if line.strip()]
    names = {s["slug"]: s["species"] for s in json.loads((curated / "species.json").read_text(encoding="utf-8"))}
    datasets = json.loads((curated / "datasets.json").read_text(encoding="utf-8"))
    phenotypes = cohorts.read_jsonl(curated / "phenotypes.jsonl")
    rows = [tie_row(phenotypes, datasets, names[slug], slug) for slug in slugs]
    text = "\t".join(COLUMNS) + "\n" + "".join("\t".join(str(r[c]) for c in COLUMNS) + "\n" for r in rows)
    out.write_text(text, encoding="utf-8")
    sys.stdout.write(text)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
