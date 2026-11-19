#!/usr/bin/env python3
"""Verify the cohort-r2 evidence files a cleaned run reads, as push-data carries them: flat under the data root's
bindings/ (plan revision 7). Run from the reviewed checkout, never from the evidence.

The evidence root (`grm-research-mtb-4k-cohort-r2`, 88 files) holds scripts, logs and proofs that push-data does
not carry, so its own scripts/verify_evidence.py runs on the laptop against the surviving copy, and its receipt is
kept with the run. On Boat this script checks what preparation reads:

- the evidence's SHA256SUMS, pushed as mtb-cohort-r2-sha256sums.txt, hashes to EVIDENCE_SUMS_SHA256;
- each file prepare_execution_roots.EVIDENCE_FILES names hashes as SHA256SUMS lists its evidence path, and as
  prepare_execution_roots.EXPECTED pins it;
- the cohort checks the evidence script makes on these files: 4,000 unique genomes, the frozen 2k is their exact
  first rows, added is the rest, and the FASTA manifest lists the 4,000 in order, each second read reproducing.

usage: python3 -I verify_evidence.py --root BINDINGS_DIR
Prints a JSON receipt; exits 0 when verified, 2 otherwise."""
import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path

LANE = Path(__file__).resolve().parent
if sys.path[:1] != [str(LANE)]:
    sys.path.insert(0, str(LANE))
import prepare_execution_roots as prep  # noqa: E402

EVIDENCE_SUMS = "mtb-cohort-r2-sha256sums.txt"
EVIDENCE_SUMS_SHA256 = "d58008c92ccdc07fde911c6b0cf055d5e473c7416dc4887185c8500b82fef833"
# each EVIDENCE_FILES key's path in the evidence root
EVIDENCE_PATHS = {"full": "results/genomes-4000.txt", "prefix": "proofs/frozen-candidate-2000.txt",
                  "added": "results/added-2000.txt", "fasta-manifest": "results/fasta-manifest.tsv"}


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def verify(root: Path):
    """{check: bool} and the problems found."""
    problems = []
    sums_bytes = (root / EVIDENCE_SUMS).read_bytes()
    checks = {"sha256sums_pinned": sha256(sums_bytes) == EVIDENCE_SUMS_SHA256}
    listed = {}
    for line in sums_bytes.decode().splitlines():
        digest, sep, name = line.partition("  ")
        listed[name[2:] if name.startswith("./") else name] = digest
    data = {}
    for key, flat in prep.EVIDENCE_FILES.items():
        data[key] = (root / flat).read_bytes()
        digest = sha256(data[key])
        checks[f"{key}_as_sha256sums_lists"] = listed.get(EVIDENCE_PATHS[key]) == digest
        checks[f"{key}_as_preparation_pins"] = prep.EXPECTED[key][0] == digest
    ids = data["full"].decode().splitlines()
    frozen = data["prefix"].decode().splitlines()
    added = data["added"].decode().splitlines()
    rows = list(csv.DictReader(data["fasta-manifest"].decode().splitlines(), delimiter="\t"))
    checks.update({
        "cohort_4000_unique": len(ids) == 4000 and len(set(ids)) == 4000,
        "frozen_2k_exact_prefix": ids[:2000] == frozen and data["full"].startswith(data["prefix"]),
        "added_2000_exact": added == ids[2000:] and len(added) == 2000,
        "fasta_manifest_4000_ordered": len(rows) == 4000 and [r["genome_id"] for r in rows] == ids,
        "all_fasta_second_reads_reproduce": all(
            r["sha256"] == r["second_pass_sha256"] and r["size_bytes"] == r["second_pass_size_bytes"]
            and r["read_hash_reproducible"] == "True" for r in rows),
    })
    problems += [name for name, ok in checks.items() if not ok]
    return checks, problems


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0], allow_abbrev=False)
    ap.add_argument("--root", required=True, type=Path, help="the folder holding the flat evidence files (bindings/)")
    args = ap.parse_args(argv)
    try:
        checks, problems = verify(args.root)
    except (OSError, ValueError, KeyError) as exc:
        checks, problems = {}, [f"cannot read the evidence under {args.root}: {type(exc).__name__}: {exc}"]
    print(json.dumps({"verified": not problems, "root": str(args.root), "checks": checks, "problems": problems},
                     indent=2, sort_keys=True))
    if problems:
        print(f"verify_evidence: STOP: {problems[0]}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
