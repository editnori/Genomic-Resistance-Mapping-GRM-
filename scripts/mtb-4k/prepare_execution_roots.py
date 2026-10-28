#!/usr/bin/env python3
"""Prepare the matched MTB roots. NOT EXECUTED by the wiring or fix lanes.

usage: python3 -I prepare_execution_roots.py RUN_ROOT
       python3 -I prepare_execution_roots.py --verify-frozen RUN_ROOT

Cohorts, each its own data root under RUN_ROOT/cohorts/<name>/data:

  pilot256        rows 1-256 of the frozen order (the small biological pilot)
  pilot256-stock  the same 256, for the stock count + matrix comparator
  prefix          rows 1-2000 (the frozen 2k candidate, byte for byte)
  added           rows 2001-4000
  full            all 4,000
  matched2k       the size-matched control: the 2,000 of the 4,000 with the
                  smallest SHA-256 of MATCHED_SEED + ":" + genome_id, in
                  frozen order (fixed here, before any result exists)

Every root gets real `species/<slug>` and `fasta` directories (FASTA are
links to the manifest-verified frozen files, never copies) and links to the
frozen candidate's read-only `catalog` and `curated`. `bvbrc` is a real,
lane-local directory: `genome_amr.jsonl` and `genomes.jsonl` link to the
frozen records, and `features/` is a lane-local copy of the candidate's
annotation cache, each file hash-checked against the manifest written here.
`grm learn` writes newly fetched annotations only there (and refuses a
`bvbrc` or `features` link; see src/annotate.rs), so training never writes
into the frozen candidate root.

RUN_ROOT/frozen-inputs.json records the SHA-256 of every frozen input this
lane reads or links (catalog, curated, BV-BRC records, annotation cache).
--verify-frozen recomputes them after any stage; a difference is a stop.
"""

import argparse
import csv
import hashlib
import json
import os
import shutil
import sys
from pathlib import Path

EVIDENCE = Path("/root/grm-research/grm-research-mtb-4k-cohort-r2")
CANDIDATE = Path("/root/grm-extra-species-analysis-20261027/data")
SLUG = "mycobacterium_tuberculosis"
MATCHED_SEED = "grm-mtb-4k-matched2k-20261028"
PILOT = 256
EXPECTED = {
    "full": ("80faef89571ecd2fd658204abd39fd4b55e5a86bfe8e4f379274162868f44eaf", 4000),
    "prefix": ("57af1c86c35fb0a4c62d3c0784a09804a919ba75941159d0c7b3cf56870aa2e3", 2000),
    "added": ("a980d9e014fb90972a46b3fc01e65a9654deb58935d3ddf0f54018c19e7def85", 2000),
}
ALLOWED_FASTA_PARENTS = {
    Path("/root/grm/data/fasta"),
    Path("/root/grm-research/grm-research-mtb-scaleup-r1/fasta"),
}
BVBRC_RECORDS = ["genome_amr.jsonl", "genomes.jsonl"]


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def matched_selection(full_ids):
    """The size-matched control: half of the full cohort, chosen by a hash
    fixed before any result, kept in frozen order."""
    rank = sorted(full_ids, key=lambda g: hashlib.sha256(f"{MATCHED_SEED}:{g}".encode()).hexdigest())
    chosen = set(rank[: len(full_ids) // 2])
    return [g for g in full_ids if g in chosen]


def frozen_manifest(candidate: Path):
    """SHA-256 of every frozen input file the roots read or link."""
    out = {}
    for sub in ["catalog", "curated"]:
        d = candidate / sub
        assert d.is_dir() and not d.is_symlink(), d
        for p in sorted(d.rglob("*")):
            if p.is_file():
                out[str(p)] = sha256_file(p)
    for name in BVBRC_RECORDS:
        p = candidate / "bvbrc" / name
        assert p.is_file() and not p.is_symlink(), p
        out[str(p)] = sha256_file(p)
    features = candidate / "bvbrc" / "features"
    for p in sorted(features.iterdir()) if features.is_dir() else []:
        assert p.is_file() and not p.is_symlink(), p
        out[str(p)] = sha256_file(p)
    return out


def prepare_bvbrc(data: Path, candidate: Path, manifest: dict):
    """A lane-local bvbrc: frozen records linked, the annotation cache copied
    and checked file by file against `manifest`."""
    bvbrc = data / "bvbrc"
    bvbrc.mkdir()
    for name in BVBRC_RECORDS:
        os.symlink(candidate / "bvbrc" / name, bvbrc / name)
    features = bvbrc / "features"
    features.mkdir()
    source = candidate / "bvbrc" / "features"
    copied = 0
    for p in sorted(source.iterdir()) if source.is_dir() else []:
        target = features / p.name
        shutil.copyfile(p, target)
        assert sha256_file(target) == manifest[str(p)], f"{p} changed while it was copied"
        copied += 1
    return copied


def verify_frozen(run_root: Path) -> int:
    recorded = json.loads((run_root / "frozen-inputs.json").read_text())
    candidate = Path(recorded["candidate"])
    now = frozen_manifest(candidate)
    changed = sorted(k for k in set(now) | set(recorded["sha256"]) if now.get(k) != recorded["sha256"].get(k))
    print(json.dumps({"frozen_files": len(now), "changed": changed}, indent=2))
    return 1 if changed else 0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_root", type=Path)
    parser.add_argument("--verify-frozen", action="store_true")
    args = parser.parse_args()
    if args.verify_frozen:
        return verify_frozen(args.run_root)
    if os.path.lexists(args.run_root):
        raise SystemExit(f"refusing existing run root: {args.run_root}")

    full_bytes = (EVIDENCE / "results/genomes-4000.txt").read_bytes()
    prefix_bytes = (EVIDENCE / "proofs/frozen-candidate-2000.txt").read_bytes()
    added_bytes = (EVIDENCE / "results/added-2000.txt").read_bytes()
    lines = full_bytes.splitlines(keepends=True)
    assert len(lines) == 4000
    assert b"".join(lines[:2000]) == prefix_bytes
    assert b"".join(lines[2000:]) == added_bytes
    cohort_bytes = {"full": full_bytes, "prefix": prefix_bytes, "added": added_bytes}
    cohort_ids = {}
    for name, data in cohort_bytes.items():
        digest, count = EXPECTED[name]
        ids = data.decode().splitlines()
        assert len(ids) == count and len(set(ids)) == count
        assert sha256_bytes(data) == digest
        assert all(Path(x).name == x and x not in {"", ".", ".."} for x in ids)
        cohort_ids[name] = ids
    cohort_ids["matched2k"] = matched_selection(cohort_ids["full"])
    cohort_ids["pilot256"] = cohort_ids["prefix"][:PILOT]
    cohort_ids["pilot256-stock"] = cohort_ids["pilot256"]
    for name in ("matched2k", "pilot256", "pilot256-stock"):
        cohort_bytes[name] = "".join(f"{g}\n" for g in cohort_ids[name]).encode()
    assert cohort_bytes["pilot256"] == b"".join(lines[:PILOT])
    assert len(cohort_ids["matched2k"]) == 2000

    with (EVIDENCE / "results/fasta-manifest.tsv").open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 4000
    assert [row["genome_id"] for row in rows] == cohort_ids["full"]
    sources = {}
    for row in rows:
        assert row["read_hash_reproducible"] == "True"
        assert row["unchanged_since_resolution"] == "True"
        source = Path(row["canonical_source_path"])
        assert source.parent in ALLOWED_FASTA_PARENTS
        assert source.name == row["genome_id"] + ".fna.gz"
        assert source.is_file() and not source.is_symlink()
        assert source.resolve(strict=True) == source
        assert source.stat().st_size == int(row["size_bytes"])
        sources[row["genome_id"]] = source

    manifest = frozen_manifest(CANDIDATE)
    args.run_root.mkdir()
    (args.run_root / "frozen-inputs.json").write_text(
        json.dumps({"candidate": str(CANDIDATE), "sha256": manifest}, indent=1) + "\n"
    )
    seeded = {}
    for name, ids in cohort_ids.items():
        data = args.run_root / "cohorts" / name / "data"
        species = data / "species" / SLUG
        fasta = data / "fasta"
        species.mkdir(parents=True)
        fasta.mkdir()
        for directory in ["catalog", "curated"]:
            target = CANDIDATE / directory
            assert target.is_dir()
            os.symlink(target, data / directory, target_is_directory=True)
        seeded[name] = prepare_bvbrc(data, CANDIDATE, manifest)
        with (species / "genomes.txt").open("xb") as handle:
            handle.write(cohort_bytes[name])
        for genome_id in ids:
            os.symlink(sources[genome_id], fasta / f"{genome_id}.fna.gz")

    for directory in ["logs", "backups", "audit", "target"]:
        (args.run_root / directory).mkdir()
    prefix = set(cohort_ids["prefix"])
    receipt = {
        "run_root": str(args.run_root),
        "cohorts": {name: len(ids) for name, ids in cohort_ids.items()},
        "genome_sha256": {name: sha256_bytes(data) for name, data in cohort_bytes.items()},
        "matched2k": {"seed": MATCHED_SEED, "from_prefix": sum(g in prefix for g in cohort_ids["matched2k"]),
                      "from_added": sum(g not in prefix for g in cohort_ids["matched2k"])},
        "fasta_links": sum(len(ids) for ids in cohort_ids.values()),
        "fasta_bytes_copied": 0,
        "input_root": str(CANDIDATE),
        "frozen_inputs": len(manifest),
        "bvbrc": "lane-local directory per cohort; records linked read-only; features copied and hash-checked",
        "features_seeded_per_cohort": seeded,
    }
    (args.run_root / "preparation.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(receipt, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
