"""The cohorts of plan block C2 (docs/cap-arms.md): the release's 2,000 and the wider caps.

  python3 -I cohorts.py <curated dir> <release species json> <species slug> <out dir> [<cap> ...]

<curated dir> holds the rebuilt revision 14 table's datasets.json and
phenotypes.jsonl; <release species json> is the release's
species/<slug>.json, whose genomes are the default's cohort. The caps
(default 4000 4096) are the wider arms, cap-<n> for each, the largest also
named cap-all. Each wider cohort is the top <cap> genomes by the release's
selection rule, `curate::select_genomes` (src/curate.rs), on the rebuilt
table: every genome with a reported R or S call on a drug of the species'
datasets, not a conflict and not inferred from a MIC, scores the sum of
1 / (size of that call's class) over those calls, summed in the
phenotypes' order as the Rust sums them; genomes sort by score, highest
first, then by id. The ranking does not depend on the cap, so the cohorts
nest.

Refuses, writing nothing: a cap above 4,096, the reviewed matrix gate
(docs/matrix-scaling.md); a wider cohort that leaves out a genome of the
release's 2,000; a cohort genome the rebuilt table does not rank; a cohort
with a genome listed twice.

Writes <out dir>/<arm>.txt per arm, one genome id per line sorted by id
(the canonical order of scripts/input-comparison/run.py), which stays on
the box; <out dir>/cohorts.json, counts only: per arm its genomes, the
genomes it adds to the default, how many of the release's 2,000 it holds,
and its list's line count and sha256; and <out dir>/added-labelled.tsv,
counts only, per dataset of the species in the rebuilt table and wider
arm: the species and antibiotic in full, the arm's genomes with a call on
the dataset's drug (labelled, as the rule reads a call), and its added
genomes with one, resistant and susceptible. An arm with no added
labelled genome on a dataset is none-added there (none_added 1): it is not
fitted and stays out of the count, on both splits. A dataset with few is
fitted and counted; no threshold is set. Prints cohorts.json.
"""

import hashlib
import json
import sys
from collections import defaultdict
from pathlib import Path

GATE = 4096
DEFAULT = "cap-2000"
ALL = "cap-all"


def read_jsonl(path):
    with open(path, encoding="utf-8") as f:
        return [json.loads(line) for line in f if line.strip()]


def reported(p):
    """The call the rule counts, as `Phenotype::reported_label` gives it: R or S, not a conflict, not from a MIC."""
    if p.get("conflict") or p.get("label_from_mic"):
        return None
    return p.get("reported") if p.get("reported") in ("R", "S") else None


def rank(phenotypes, datasets, species):
    """Every genome the rule ranks for `species`, best first."""
    targets = {d["antibiotic"] for d in datasets if d["species"] == species}
    sizes = defaultdict(int)
    labels = defaultdict(list)
    for p in phenotypes:
        if p["species"] != species or p["antibiotic"] not in targets:
            continue
        call = reported(p)
        if call is None:
            continue
        sizes[(p["antibiotic"], call)] += 1
        labels[p["genome_id"]].append((p["antibiotic"], call))
    scored = []
    for g, ls in labels.items():
        t = 0.0
        for a, c in ls:
            t += 1.0 / sizes[(a, c)]
        scored.append((-t, g))
    scored.sort()
    return [g for _, g in scored]


def cohorts(order, release, caps):
    """{arm: sorted ids}: the default's release list and the top <cap> of `order` for each cap."""
    if len(set(release)) != len(release):
        raise ValueError("the release's list holds a genome twice")
    if len(set(order)) != len(order):
        raise ValueError("the ranking holds a genome twice")
    ranked = set(order)
    gone = [g for g in release if g not in ranked]
    if gone:
        raise ValueError(f"{len(gone)} genomes of the release's list are not ranked on the rebuilt table")
    out = {DEFAULT: sorted(release)}
    caps = sorted(set(caps))
    for cap in caps:
        if cap > GATE:
            raise ValueError(f"cap {cap} is past {GATE:,}, the reviewed matrix gate")
        if cap <= len(release):
            raise ValueError(f"cap {cap} is not wider than the default's {len(release):,}")
    for cap in caps:
        cohort = order[:cap]
        out[f"cap-{cap}"] = sorted(cohort)
    out[ALL] = out.pop(f"cap-{caps[-1]}")
    return out


def check(arms, release, order):
    """Raises unless each arm holds no genome twice, holds only ranked genomes, stays inside the gate, and each
    wider arm holds the release's whole list."""
    ranked, base = set(order), set(release)
    for arm, ids in arms.items():
        if len(ids) > GATE:
            raise ValueError(f"{arm} holds {len(ids):,} genomes, past {GATE:,}, the reviewed matrix gate")
        if len(set(ids)) != len(ids):
            raise ValueError(f"{arm} holds a genome twice")
        stray = [g for g in ids if g not in ranked]
        if stray:
            raise ValueError(f"{arm} holds {len(stray)} genomes the rebuilt table does not rank")
        if arm != DEFAULT:
            dropped = base - set(ids)
            if dropped:
                raise ValueError(f"{arm} drops {len(dropped)} genomes of the default's 2,000")


def list_text(ids):
    return "".join(g + "\n" for g in ids)


def counts(arms, order):
    base = set(arms[DEFAULT])
    return {"ranked": len(order),
            "arms": {arm: {"genomes": len(ids), "added": len(set(ids) - base), "holds_default": len(base & set(ids)),
                           "lines": len(ids), "sha256": hashlib.sha256(list_text(ids).encode("utf-8")).hexdigest()}
                     for arm, ids in arms.items()}}


LABELLED_COLUMNS = ["dataset", "species", "antibiotic", "arm", "labelled", "added", "added_labelled", "added_resistant",
                    "added_susceptible", "none_added"]


def added_labelled(arms, phenotypes, datasets, species):
    """added-labelled.tsv's rows: per dataset of `species` and wider arm, the arm's labelled genomes and its added
    genomes with a call on the dataset's drug, by class."""
    calls = defaultdict(lambda: defaultdict(set))
    for p in phenotypes:
        call = reported(p) if p["species"] == species else None
        if call is not None:
            calls[p["antibiotic"]][call].add(p["genome_id"])
    base = set(arms[DEFAULT])
    out = []
    for d in sorted((d for d in datasets if d["species"] == species), key=lambda d: d["id"]):
        r, s = calls[d["antibiotic"]]["R"], calls[d["antibiotic"]]["S"]
        for arm in sorted(a for a in arms if a != DEFAULT):
            ids = set(arms[arm])
            added = ids - base
            n = len(added & (r | s))
            out.append({"dataset": d["id"], "species": species, "antibiotic": d["antibiotic"], "arm": arm,
                        "labelled": len(ids & (r | s)), "added": len(added), "added_labelled": n,
                        "added_resistant": len(added & r), "added_susceptible": len(added & s),
                        "none_added": int(n == 0)})
    return out


def main(argv):
    if len(argv) < 5:
        raise SystemExit(__doc__)
    curated, release_json, slug, out = Path(argv[1]), Path(argv[2]), argv[3], Path(argv[4])
    caps = [int(c) for c in argv[5:]] or [4000, 4096]
    with open(release_json, encoding="utf-8") as f:
        rel = json.load(f)
    release = [g["id"] for g in rel["genomes"]]
    with open(curated / "species.json", encoding="utf-8") as f:
        names = {s["slug"]: s["species"] for s in json.load(f)}
    if slug not in names:
        raise SystemExit(f"the rebuilt table has no species {slug}")
    with open(curated / "datasets.json", encoding="utf-8") as f:
        datasets = json.load(f)
    phenotypes = read_jsonl(curated / "phenotypes.jsonl")
    order = rank(phenotypes, datasets, names[slug])
    try:
        arms = cohorts(order, release, caps)
        check(arms, release, order)
    except ValueError as e:
        raise SystemExit(f"{slug}: {e}") from None
    out.mkdir(parents=True, exist_ok=True)
    for arm, ids in arms.items():
        (out / f"{arm}.txt").write_text(list_text(ids), encoding="utf-8")
    doc = {"species": names[slug], "slug": slug, **counts(arms, order)}
    (out / "cohorts.json").write_text(json.dumps(doc, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    rows = added_labelled(arms, phenotypes, datasets, names[slug])
    with open(out / "added-labelled.tsv", "w", encoding="utf-8") as f:
        f.write("\t".join(LABELLED_COLUMNS) + "\n")
        for r in rows:
            f.write("\t".join(str(r[c]) for c in LABELLED_COLUMNS) + "\n")
    print(json.dumps(doc, sort_keys=True))


if __name__ == "__main__":
    main(sys.argv)
