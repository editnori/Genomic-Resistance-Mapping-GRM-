"""The list and the floor of plan block C3's run on its datasets of M (docs/input-arms.md, The run on the datasets of M).

  python3 -I m_list.py <data root> <records> <lists.tsv> <species slug> <datasets> <list dir> <floor.tsv>

Run on the run's machine, before any score is read: it reads the pinned
label table, the records that table was built from and the species' matrix,
and no fit, report or score.

The pinned table (curated/phenotypes.jsonl, from curated.tar) holds each
genome and antibiotic pair's consolidated call, `reported` and `conflict`,
and no record's call. So the records file <records> (the kit's
genome_amr.jsonl: the 1,285,111 lab records the pinned table was rebuilt
from on 21 Nov 2026, which the run script pins by its own sha256) is named
with its sha256 before the list is built: that line is printed and
flushed first, before any check, so a stop still shows it. Each record's raw
antibiotic name goes through curated/drug_names.json, as `grm curate`
resolves it; a record whose name is not a drug, or whose genome has no row
in the table, is passed over. Its call is read with parse_call's words
(src/curate.rs). Every pair of the listed datasets' antibiotics must then
consolidate under r-or-s-wins to the table's own `reported` and `conflict`,
and each name's record count must equal drug_names.json's: a records file
that is not the one the table was built from stops the list.

A dataset's cohort is grm's (src/analysis.rs, cohort): the table's genomes
for its antibiotic among the matrix genomes, with an unconflicted reported R
or S and no label from an MIC. Its size and classes must equal lists.tsv's.
The list holds the cohort genomes whose records for that antibiotic hold an
intermediate call; each is R or S under r-or-s-wins. The floor: a dataset
passes it (`enters` 1) when each class keeps at least 25 genomes, each
counted once, among the cohort genomes (all of which the lineage split
scores) less the list.

Writes <list dir>/<dataset>.txt (the listed genome ids, one a line; it
stays on the run's machine) and <floor.tsv> (counts only), and prints the
records file's sha256, then the counts per dataset. Exits 1, writing
nothing, on any difference.
"""

import csv
import hashlib
import json
import sys
from pathlib import Path

FLOOR = 25
R_WORDS = {"resistant", "nonsusceptible", "non-susceptible"}
S_WORDS = {"susceptible"}
I_WORDS = {"intermediate", "susceptible-dose dependent", "reduced susceptibility"}
FLOOR_COLUMNS = ["dataset", "species_slug", "species", "antibiotic", "cohort", "resistant", "susceptible",
                 "listed", "listed_resistant", "listed_susceptible", "second_resistant", "second_susceptible",
                 "floor", "enters"]


class Stop(Exception):
    pass


def value_text(v):
    """A record field as grm reads it (src/bvbrc.rs, value_text): trimmed text, a list joined with ';', None if blank."""
    if v is None:
        return None
    if isinstance(v, str):
        return v.strip() or None
    if isinstance(v, list):
        parts = [p for p in (value_text(x) for x in v) if p is not None]
        return ";".join(parts) if parts else None
    if isinstance(v, bool):
        return "true" if v else "false"
    return repr(v) if isinstance(v, float) else str(v)


def parse_call(text):
    """src/curate.rs, parse_call: R, S, I or None, on the trimmed lowercase text."""
    if text is None:
        return None
    t = text.strip().lower()
    return "R" if t in R_WORDS else "S" if t in S_WORDS else "I" if t in I_WORDS else None


def consolidate(calls):
    """r-or-s-wins (src/curate.rs, LabelRule::ROrSWins): (call, conflict). S with I is S, I with R is R, R with S none."""
    if "R" in calls and "S" in calls:
        return None, True
    for c in ("R", "S", "I"):
        if c in calls:
            return c, False
    return None, False


def read_jsonl(path):
    """Each line's JSON, one at a time: the records file is too large to hold whole."""
    with open(path, encoding="utf-8") as f:
        for line in f:
            if line.strip():
                yield json.loads(line)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def records_line(records_path):
    """The line that names the records file with its sha256: main prints it before any check."""
    return (f"the pinned table curated/phenotypes.jsonl holds no record's call; the calls are read from "
            f"{Path(records_path).name}, sha256 {sha256(records_path)}")


def build(root, records_path, lists_rows, slug, datasets):
    """(floor rows, {dataset: listed genome ids}, printout lines). Raises Stop on any difference."""
    root = Path(root)
    lines = []
    listed_rows = {r["dataset"]: r for r in lists_rows}
    specs = {d["id"]: d for d in json.loads((root / "curated" / "datasets.json").read_text())}
    for d in datasets:
        if d not in listed_rows or listed_rows[d]["species_slug"] != slug or listed_rows[d]["in_M"] != "1":
            raise Stop(f"{d} is not a dataset of M of {slug} in lists.tsv")
        if d not in specs:
            raise Stop(f"{d} is not in curated/datasets.json")
    drugs = {specs[d]["antibiotic"] for d in datasets}
    names = json.loads((root / "curated" / "drug_names.json").read_text())["names"]
    resolve = {n["raw"]: (n["kind"], n["canonical"]) for n in names}
    want_records = {n["raw"]: int(n["records"]) for n in names}
    genomes = json.loads((root / "species" / slug / "matrix" / "matrix.json").read_text())["genomes"]
    column = {g: i for i, g in enumerate(genomes)}
    table = list(read_jsonl(root / "curated" / "phenotypes.jsonl"))
    passing = {p["genome_id"] for p in table}
    pairs = {(p["genome_id"], p["antibiotic"]): p for p in table if p["antibiotic"] in drugs}

    calls, seen = {}, {}
    for rec in read_jsonl(records_path):
        raw = value_text(rec.get("antibiotic"))
        if raw is None:
            raise Stop("a record has no antibiotic name")
        seen[raw] = seen.get(raw, 0) + 1
        kind, canonical = resolve.get(raw, (None, None))
        if kind is None:
            raise Stop("a record's antibiotic name is not in curated/drug_names.json")
        if kind != "drug" or canonical not in drugs:
            continue
        genome = value_text(rec.get("genome_id"))
        if genome not in passing:
            continue
        call = parse_call(value_text(rec.get("resistant_phenotype")))
        if call is not None:
            calls.setdefault((genome, canonical), set()).add(call)
    if seen != want_records:
        diff = sorted(set(seen) ^ set(want_records)) or sorted(k for k in seen if seen[k] != want_records.get(k))
        raise Stop(f"the records' name counts differ from curated/drug_names.json for {len(diff)} names")
    bad = 0
    for key in set(calls) | {k for k, p in pairs.items() if p.get("reported") is not None}:
        p = pairs.get(key)
        got = consolidate(calls.get(key, set()))
        if p is None or got != (p.get("reported"), bool(p.get("conflict"))):
            bad += 1
    if bad:
        raise Stop(f"{bad} genome and antibiotic pairs of the records do not consolidate to the table's call")
    lines.append(f"records: every pair of {len(drugs)} antibiotics consolidates to the table's call under r-or-s-wins")

    rows, lists = [], {}
    for d in datasets:
        spec, lrow = specs[d], listed_rows[d]
        cohort = sorted((p for (g, a), p in pairs.items() if a == spec["antibiotic"] and g in column
                         and p.get("reported") in ("R", "S") and not p.get("conflict") and not p.get("label_from_mic")),
                        key=lambda p: column[p["genome_id"]])
        ids = [p["genome_id"] for p in cohort]
        if len(set(ids)) != len(ids):
            raise Stop(f"{d}: a cohort genome is listed twice")
        r = sum(1 for p in cohort if p["reported"] == "R")
        s = len(cohort) - r
        if (len(cohort), r, s) != (int(lrow["samples"]), int(lrow["resistant"]), int(lrow["susceptible"])):
            raise Stop(f"{d}: the cohort holds {len(cohort)} genomes, {r} R and {s} S; "
                       f"lists.tsv gives {lrow['samples']}, {lrow['resistant']} and {lrow['susceptible']}")
        held = [p for p in cohort if "I" in calls.get((p["genome_id"], spec["antibiotic"]), set())]
        lr = sum(1 for p in held if p["reported"] == "R")
        ls = len(held) - lr
        row = {"dataset": d, "species_slug": slug, "species": lrow["species"], "antibiotic": lrow["antibiotic"],
               "cohort": len(cohort), "resistant": r, "susceptible": s, "listed": len(held),
               "listed_resistant": lr, "listed_susceptible": ls, "second_resistant": r - lr,
               "second_susceptible": s - ls, "floor": FLOOR}
        row["enters"] = int(row["second_resistant"] >= FLOOR and row["second_susceptible"] >= FLOOR)
        rows.append(row)
        lists[d] = [p["genome_id"] for p in held]
        lines.append(f"{d}: cohort {len(cohort)} ({r} R, {s} S); listed {len(held)} ({lr} R, {ls} S); "
                     f"second reading {r - lr} R and {s - ls} S, {'enters' if row['enters'] else 'under the floor of 25'}")
    return rows, lists, lines


def write_floor(path, rows):
    with open(path, "w", encoding="utf-8", newline="") as f:
        f.write("\t".join(FLOOR_COLUMNS) + "\n")
        for r in rows:
            f.write("\t".join(str(r[c]) for c in FLOOR_COLUMNS) + "\n")


def read_tsv(path):
    with open(path, encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def main(argv):
    if len(argv) != 7:
        print(__doc__.strip().splitlines()[2], file=sys.stderr)
        return 2
    root, records, lists_path, slug, datasets, list_dir, floor_path = argv
    print(records_line(records), flush=True)
    try:
        rows, lists, lines = build(root, records, read_tsv(lists_path), slug, datasets.split(","))
    except Stop as e:
        print(f"m_list: stopped: {e}")
        return 1
    list_dir = Path(list_dir)
    list_dir.mkdir(parents=True, exist_ok=True)
    for d, ids in lists.items():
        (list_dir / f"{d}.txt").write_text("".join(f"{g}\n" for g in ids))
    Path(floor_path).parent.mkdir(parents=True, exist_ok=True)
    write_floor(floor_path, rows)
    print("\n".join(lines))
    print(f"floor: {sum(r['enters'] for r in rows)} of {len(rows)} datasets enter the second reading")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
