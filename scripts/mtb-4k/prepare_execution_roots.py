#!/usr/bin/env python3
"""Prepare the matched MTB roots. NOT EXECUTED by the wiring or fix lanes.

usage: python3 -I prepare_execution_roots.py RUN_ROOT
       python3 -I prepare_execution_roots.py --verify-frozen RUN_ROOT
       python3 -I prepare_execution_roots.py --check-lane COHORT_DATA_ROOT

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
frozen records, and `features/` holds lane-local copies of the candidate's
cached annotations for that cohort's own genomes only, each hash-checked
against the manifest written here. `grm learn` reads and fetches annotations
only for genomes in its cohort, writes newly fetched ones only there, and
refuses a `bvbrc` or `features` link (see src/annotate.rs), so training
never writes into the frozen candidate root. The two pilot roots never
train and get an empty `features/`.

Before it creates anything, the script sums the annotation bytes it will
copy and refuses unless at least MIN_FREE_AFTER bytes stay free on the run
root's filesystem afterwards (the plan's 45 GiB gate).

RUN_ROOT/frozen-inputs.json records the SHA-256 of every frozen input this
lane reads or links: catalog, curated, BV-BRC records, the candidate's whole
annotation cache, every FASTA in the linked store (checked here against the
cohort's FASTA manifest), each cohort's genomes.txt, and where each FASTA
link points. --verify-frozen recomputes all of them after any stage; a
difference is a stop.

--check-lane refuses a cohort data root with a link anywhere but where this
script made one (catalog, curated, the two BV-BRC record files, and each
FASTA), or with any entry that is neither a real directory nor a regular
file. GRM confines its own matrix, count and learn writes; this check also
covers `grm report`, which the plan runs before the run manifest.
"""

import argparse
import csv
import hashlib
import json
import os
import shutil
import stat
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
    # matched_selection of the frozen full list, measured read-only on Vast on 2026-10-28 and
    # frozen in report_contract.RULES, so a changed seed or selection is a stop here too.
    "matched2k": ("212a21a5d2574e0d62d379fccddb0839c4f44f5a7cd36e43a558dc2c0926c3d1", 2000),
}
ALLOWED_FASTA_PARENTS = {
    Path("/root/grm/data/fasta"),
    Path("/root/grm-research/grm-research-mtb-scaleup-r1/fasta"),
}
BVBRC_RECORDS = ["genome_amr.jsonl", "genomes.jsonl"]
# Cohorts that train and so may read or fetch annotations; the pilots never do.
TRAINING = ["prefix", "added", "full", "matched2k"]
MIN_FREE_AFTER = 45 * 1024**3


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def free_bytes(path: Path) -> int:
    s = os.statvfs(path)
    return s.f_bavail * s.f_frsize


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


def cached_features(candidate: Path, ids):
    """The candidate's cached annotation files for `ids`, in `ids` order:
    `annotate::load_features` reads `bvbrc/features/<genome_id>.json` only
    for genomes of the cohort it trains on."""
    source = candidate / "bvbrc" / "features"
    out = []
    for gid in ids:
        p = source / f"{gid}.json"
        if p.is_file() and not p.is_symlink():
            out.append(p)
    return out


def prepare_bvbrc(data: Path, candidate: Path, manifest: dict, ids):
    """A lane-local bvbrc: frozen records linked, and the cached annotations
    of `ids` (none for a pilot) copied and checked file by file against
    `manifest`. Returns the files and bytes copied."""
    bvbrc = data / "bvbrc"
    bvbrc.mkdir()
    for name in BVBRC_RECORDS:
        os.symlink(candidate / "bvbrc" / name, bvbrc / name)
    features = bvbrc / "features"
    features.mkdir()
    copied = {"files": 0, "bytes": 0}
    for p in cached_features(candidate, ids):
        target = features / p.name
        shutil.copyfile(p, target)
        assert sha256_file(target) == manifest[str(p)], f"{p} changed while it was copied"
        copied["files"] += 1
        copied["bytes"] += target.stat().st_size
    return copied


def verify_frozen(run_root: Path) -> int:
    recorded = json.loads((run_root / "frozen-inputs.json").read_text())
    candidate = Path(recorded["candidate"])
    now = frozen_manifest(candidate)
    changed = sorted(k for k in set(now) | set(recorded["sha256"]) if now.get(k) != recorded["sha256"].get(k))
    fasta_changed = []
    for gid, rec in sorted(recorded["fasta"].items()):
        source = Path(rec["source"])
        if not source.is_file() or source.is_symlink() or sha256_file(source) != rec["sha256"]:
            fasta_changed.append(str(source))
    links_changed = []
    for name, digest in sorted(recorded["cohort_genomes_sha256"].items()):
        data = run_root / "cohorts" / name / "data"
        genomes = data / "species" / SLUG / "genomes.txt"
        if not genomes.is_file() or genomes.is_symlink() or sha256_file(genomes) != digest:
            links_changed.append(f"{genomes}: cohort list changed")
            continue
        ids = genomes.read_text().splitlines()
        fasta = data / "fasta"
        names = sorted(os.listdir(fasta)) if fasta.is_dir() and not fasta.is_symlink() else []
        if names != sorted(f"{g}.fna.gz" for g in ids):
            links_changed.append(f"{fasta}: entries differ from the cohort")
        for gid in ids:
            link = fasta / f"{gid}.fna.gz"
            if not link.is_symlink() or os.readlink(link) != recorded["fasta"][gid]["source"]:
                links_changed.append(str(link))
    print(json.dumps({"frozen_files": len(now), "fasta_files": len(recorded["fasta"]), "changed": changed,
                      "fasta_changed": fasta_changed, "links_changed": links_changed}, indent=2))
    return 1 if changed or fasta_changed or links_changed else 0


def check_lane(data: Path):
    """Problems with a cohort data root: a link anywhere but where this
    script made one, or an entry that is not a real directory or file."""
    problems = []
    if data.is_symlink() or not data.is_dir():
        return [f"{data} is not a real directory"]
    allowed = {data / "catalog": CANDIDATE / "catalog", data / "curated": CANDIDATE / "curated"}
    for name in BVBRC_RECORDS:
        allowed[data / "bvbrc" / name] = CANDIDATE / "bvbrc" / name
    for d in [data / "species", data / "species" / SLUG, data / "fasta", data / "bvbrc", data / "bvbrc" / "features"]:
        if d.is_symlink() or not d.is_dir():
            problems.append(f"{d} must be a real directory")
    for top, dirs, files in os.walk(data, followlinks=False):
        for name in dirs + files:
            p = Path(top) / name
            mode = os.lstat(p).st_mode
            if stat.S_ISLNK(mode):
                target = Path(os.readlink(p))
                if p in allowed:
                    if target != allowed[p]:
                        problems.append(f"{p} links to {target}, not {allowed[p]}")
                elif p.parent == data / "fasta" and name.endswith(".fna.gz"):
                    if target.parent not in ALLOWED_FASTA_PARENTS or target.name != name:
                        problems.append(f"{p} links outside the frozen FASTA stores: {target}")
                else:
                    problems.append(f"{p} is a link where GRM reads or writes lane files")
            elif not (stat.S_ISDIR(mode) or stat.S_ISREG(mode)):
                problems.append(f"{p} is neither a directory nor a regular file")
    return problems


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("path", type=Path, help="RUN_ROOT, or a cohort data root with --check-lane")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--verify-frozen", action="store_true")
    mode.add_argument("--check-lane", action="store_true")
    args = parser.parse_args()
    if args.verify_frozen:
        try:
            return verify_frozen(args.path)
        except (OSError, ValueError, KeyError) as exc:  # a missing or unreadable record is a stop, not a traceback
            print(f"prepare_execution_roots: STOP: --verify-frozen cannot read {args.path}: {type(exc).__name__}: {exc}",
                  file=sys.stderr)
            return 1
    if args.check_lane:
        problems = check_lane(args.path)
        print(json.dumps({"data_root": str(args.path), "problems": problems}, indent=2))
        return 1 if problems else 0
    run_root = args.path
    if os.path.lexists(run_root):
        raise SystemExit(f"refusing existing run root: {run_root}")

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
    assert (sha256_bytes(cohort_bytes["matched2k"]), len(cohort_ids["matched2k"])) == EXPECTED["matched2k"], "matched2k differs from its frozen list"

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
        assert sha256_file(source) == row["sha256"], f"{source} differs from the FASTA manifest"
        sources[row["genome_id"]] = {"source": str(source), "sha256": row["sha256"]}

    # The annotation copies are the only bytes this script writes in bulk.
    copies = {name: cached_features(CANDIDATE, cohort_ids[name]) for name in TRAINING}
    copy_bytes = sum(p.stat().st_size for files in copies.values() for p in files)
    free_before = free_bytes(run_root.parent)
    if free_before - copy_bytes < MIN_FREE_AFTER:
        raise SystemExit(f"refusing: copying {copy_bytes} annotation bytes would leave {free_before - copy_bytes} bytes free, "
                         f"below the {MIN_FREE_AFTER}-byte floor")

    manifest = frozen_manifest(CANDIDATE)
    run_root.mkdir()
    (run_root / "frozen-inputs.json").write_text(json.dumps({
        "candidate": str(CANDIDATE), "sha256": manifest, "fasta": sources,
        "cohort_genomes_sha256": {name: sha256_bytes(data) for name, data in cohort_bytes.items()},
    }, indent=1) + "\n")
    seeded = {}
    for name, ids in cohort_ids.items():
        data = run_root / "cohorts" / name / "data"
        species = data / "species" / SLUG
        fasta = data / "fasta"
        species.mkdir(parents=True)
        fasta.mkdir()
        for directory in ["catalog", "curated"]:
            target = CANDIDATE / directory
            assert target.is_dir()
            os.symlink(target, data / directory, target_is_directory=True)
        seeded[name] = prepare_bvbrc(data, CANDIDATE, manifest, ids if name in TRAINING else [])
        with (species / "genomes.txt").open("xb") as handle:
            handle.write(cohort_bytes[name])
        for genome_id in ids:
            os.symlink(sources[genome_id]["source"], fasta / f"{genome_id}.fna.gz")

    for directory in ["logs", "backups", "audit", "target"]:
        (run_root / directory).mkdir()
    prefix = set(cohort_ids["prefix"])
    receipt = {
        "run_root": str(run_root),
        "cohorts": {name: len(ids) for name, ids in cohort_ids.items()},
        "genome_sha256": {name: sha256_bytes(data) for name, data in cohort_bytes.items()},
        "matched2k": {"seed": MATCHED_SEED, "from_prefix": sum(g in prefix for g in cohort_ids["matched2k"]),
                      "from_added": sum(g not in prefix for g in cohort_ids["matched2k"]),
                      "sha256": sha256_bytes(cohort_bytes["matched2k"]), "frozen_sha256": EXPECTED["matched2k"][0]},
        "fasta_links": sum(len(ids) for ids in cohort_ids.values()),
        "fasta_bytes_copied": 0,
        "fasta_sha256_checked": len(sources),
        "input_root": str(CANDIDATE),
        "frozen_inputs": len(manifest),
        "bvbrc": "lane-local directory per cohort; records linked read-only; each training cohort's own cached annotations copied and hash-checked; pilots get none",
        "features_seeded_per_cohort": seeded,
        "disk": {"free_bytes_before": free_before, "annotation_bytes_copied": sum(s["bytes"] for s in seeded.values()),
                 "free_bytes_after_estimate": free_before - copy_bytes, "min_free_after": MIN_FREE_AFTER},
    }
    (run_root / "preparation.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(receipt, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
