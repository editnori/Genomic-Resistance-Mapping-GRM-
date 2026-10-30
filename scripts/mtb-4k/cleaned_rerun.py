#!/usr/bin/env python3
"""The cleaned-cohort MTB rerun's cohort bindings (plan revision 6): every cohort input that is
unknown until revision N's species check, selections and FASTA fetch exist. Each starts as an
@@NAME@@ placeholder and fails closed until a reviewed commit fills it. Nothing here is a guess.
Revision N's validator values (the revision, checker pins, fields and convergence rule) are
revision_n.py's; this file adds only what depends on the cleaned cohorts.

usage: python3 -I cleaned_rerun.py status [--group prepare|reference]
       python3 -I cleaned_rerun.py get NAME

status lists the placeholders (all, or one group) and whether each is filled; it exits 0 only
when every one listed is filled and well formed. get prints one filled value and exits 2 while it
is a placeholder or malformed.

Groups. prepare: everything prepare_execution_roots.py --cohorts cleaned and report_contract.py's
profile rN-cleaned read; the plan's sections 3-8 need it. reference: revision N's MTB 2k
candidate (its trained matrix and reports), which the plan's pending-reference checks read. It
may be filled after training; those checks stay required for the final result.

The admin's rulings of 2026-10-29 that this file and prepare_execution_roots.py enforce: the
cleaned 2k is revision N's own fresh MTB 2k selection (not the frozen prefix minus two) and must
be exactly the first rows of revision N's cap-4000 selection, which is the cleaned full cohort;
added is the rows after it; matched2k is matched_selection of the cleaned full cohort at the
prefix's size. The species check's MTB "fail" genomes are the exclusions, and they include
1773.41701 and 1773.40949; "borderline" genomes are recorded and never excluded; every cohort
genome must have a verdict in the bound post-fetch species check and none may be "fail", and each
selection's curated/species.json must report selected_unchecked 0 for MTB.

To fill: replace placeholder strings in VALUES with the reviewed values, in a commit that changes
VALUES and nothing else. The tests hold for a wholly unfilled or a wholly filled file, so no test
changes with a fill. SPECIES_CHECK_FORMAT is a decision this file checks, not free text: it must be
SPECIES_CHECK_FORMAT_PARSED, meaning that the committed species-check code writes the schemas
prepare_execution_roots.py parses (species_check.json and curated/species.json as relayed on
2026-10-29). Any other value means the parser needs a reviewed change first.
"""
import sys
from pathlib import Path, PurePosixPath

LANE = Path(__file__).resolve().parent
if sys.path[:1] != [str(LANE)]:
    sys.path.insert(0, str(LANE))
import revision_n as rev_n  # noqa: E402

REQUIRED_EXCLUSIONS = ("1773.41701", "1773.40949")  # the 2026-10-29 ruling; both fail the species check
SPECIES_CHECK_FORMAT_PARSED = "species-check-relayed-20261029"

# name: (kind, group, what it is). The plan's placeholder table lists the same names.
TABLE = {
    "SPECIES_CHECK_PATH": ("path", "prepare", "the post-fetch (second) species check's species_check.json"),
    "SPECIES_CHECK_SHA256": ("sha256", "prepare", "its SHA-256"),
    "SPECIES_CHECK_FORMAT": ("format", "prepare", f"{SPECIES_CHECK_FORMAT_PARSED} once the committed schemas match the parser"),
    "EXCLUSIONS_COUNT": ("exclusions", "prepare", "how many MTB genomes it fails (at least the two of the ruling)"),
    "RN_PREFIX_GENOMES_PATH": ("path", "prepare", "revision N's MTB 2k selection, species/mycobacterium_tuberculosis/genomes.txt"),
    "RN_PREFIX_GENOMES_SHA256": ("sha256", "prepare", "its SHA-256"),
    "RN_PREFIX_SPECIES_JSON_PATH": ("path", "prepare", "the curated/species.json of that selection"),
    "RN_PREFIX_SPECIES_JSON_SHA256": ("sha256", "prepare", "its SHA-256"),
    "RN_FULL_GENOMES_PATH": ("path", "prepare", "revision N's MTB cap-4000 selection, in genomes.txt form"),
    "RN_FULL_GENOMES_SHA256": ("sha256", "prepare", "its SHA-256"),
    "RN_FULL_SPECIES_JSON_PATH": ("path", "prepare", "the curated/species.json of that selection"),
    "RN_FULL_SPECIES_JSON_SHA256": ("sha256", "prepare", "its SHA-256"),
    "RN_INPUT_ROOT": ("path", "prepare", "the data root whose catalog, curated, BV-BRC records and annotations the cohorts link"),
    "RN_INPUT_DATASETS_SHA256": ("sha256", "prepare", "SHA-256 of RN_INPUT_ROOT/curated/datasets.json"),
    "RN_INPUT_PHENOTYPES_SHA256": ("sha256", "prepare", "SHA-256 of RN_INPUT_ROOT/curated/phenotypes.jsonl (the labels)"),
    "RN_INPUT_CATALOG_VERSION_SHA256": ("sha256", "prepare", "SHA-256 of RN_INPUT_ROOT/catalog/version.txt"),
    "RN_INPUT_CATALOG_CDS_SHA256": ("sha256", "prepare", "SHA-256 of RN_INPUT_ROOT/catalog/AMR_CDS.fa"),
    "RN_INPUT_CATALOG_REFERENCE_SHA256": ("sha256", "prepare", "SHA-256 of RN_INPUT_ROOT/catalog/ReferenceGeneCatalog.txt"),
    "RN_INPUT_MANIFEST_SHA256": ("sha256", "prepare", "SHA-256 of the sorted path<TAB>sha256 manifest of every RN_INPUT_ROOT file the "
                                 "cohort roots can read (prepare --derive-cleaned prints it)"),
    "RN_FASTA_MANIFEST_PATH": ("path", "prepare", "a fasta-manifest.tsv (the cohort evidence's columns) for cohort genomes outside the frozen 4k"),
    "RN_FASTA_MANIFEST_SHA256": ("sha256", "prepare", "its SHA-256"),
    "RN_FASTA_STORE": ("path", "prepare", "the directory holding those FASTA"),
    "CLEANED_PREFIX_SHA256": ("sha256", "prepare", "SHA-256 of the cleaned prefix list (prepare --derive-cleaned)"),
    "CLEANED_PREFIX_COUNT": ("count", "prepare", "its length"),
    "CLEANED_ADDED_SHA256": ("sha256", "prepare", "SHA-256 of the cleaned added list"),
    "CLEANED_ADDED_COUNT": ("count", "prepare", "its length"),
    "CLEANED_FULL_SHA256": ("sha256", "prepare", "SHA-256 of the cleaned full list"),
    "CLEANED_FULL_COUNT": ("count", "prepare", "its length"),
    "CLEANED_MATCHED2K_SHA256": ("sha256", "prepare", "SHA-256 of the cleaned matched2k list"),
    "CLEANED_MATCHED2K_COUNT": ("count", "prepare", "its length (the cleaned prefix's)"),
    "RN_REFERENCE_ROOT": ("path", "reference", "revision N's MTB 2k candidate run root, holding data/species/mycobacterium_tuberculosis"),
}

VALUES = {
    "SPECIES_CHECK_PATH": "/root/grm-research/grm-mtb-4k-cleaned-r10-20261030/bindings/species_check.json",
    "SPECIES_CHECK_SHA256": "d80eef2949036ff5495f5ab097569f8b8bdfbf3dcaa92dee95493866b1645b42",
    "SPECIES_CHECK_FORMAT": "species-check-relayed-20261029",
    "EXCLUSIONS_COUNT": 2,
    "RN_PREFIX_GENOMES_PATH": "/root/grm-research/grm-mtb-4k-cleaned-r10-20261030/bindings/prefix-genomes.txt",
    "RN_PREFIX_GENOMES_SHA256": "b104d8bd7ffaa2d0dc3032fe30ca0c0b32670a7099a02bcd49860e5c77f59a64",
    "RN_PREFIX_SPECIES_JSON_PATH": "/root/grm-research/grm-mtb-4k-cleaned-r10-20261030/bindings/prefix-species.json",
    "RN_PREFIX_SPECIES_JSON_SHA256": "1ab2855747ca1c079cfcdbb38f04022c3142454f443dae9f9e7684bf5a4037f4",
    "RN_FULL_GENOMES_PATH": "/root/grm-research/grm-mtb-4k-cleaned-r10-20261030/bindings/full-genomes.txt",
    "RN_FULL_GENOMES_SHA256": "f2a18a43d1743b46fbbb7f1b136a4921ba9e884f19f7942c131716a90845fd75",
    "RN_FULL_SPECIES_JSON_PATH": "/root/grm-research/grm-mtb-4k-cleaned-r10-20261030/bindings/full-species.json",
    "RN_FULL_SPECIES_JSON_SHA256": "97e80c085b5bb425d34b41c8141dddab1c2ec547b9a4f9a592639c240af0252c",
    "RN_INPUT_ROOT": "/root/grm-research/grm-mtb-4k-cleaned-r10-20261030/inputs",
    "RN_INPUT_DATASETS_SHA256": "5519c61d281a685b58b33e1870f7b4186911394773224295ec81c26c8d0096a8",
    "RN_INPUT_PHENOTYPES_SHA256": "c56481f99d037ebcc44a2cafe585a20460c872317350167af49ae3678c18415d",
    "RN_INPUT_CATALOG_VERSION_SHA256": "43ad0274317a2cbf80e01626e3ae6084a02f471be2f870d57fdbe2d35761911a",
    "RN_INPUT_CATALOG_CDS_SHA256": "bcc57417c5a170fb4a664de6898455990fcc296dbfe214b6874d761dad63962c",
    "RN_INPUT_CATALOG_REFERENCE_SHA256": "fa41ade01712841e639c91163243bace1c34755b927e69e0c8678540268a7fd4",
    "RN_INPUT_MANIFEST_SHA256": "1a9b9f19c57faf23af7116af62936440b5c14bd0473a3ec2e927a0e53822db44",
    "RN_FASTA_MANIFEST_PATH": "/root/grm-research/grm-mtb-4k-cleaned-r10-20261030/bindings/replacement-fasta-manifest.tsv",
    "RN_FASTA_MANIFEST_SHA256": "c54d3033de31b8db670f7294b4268c0e8460192779e93a0605864fba614928fa",
    "RN_FASTA_STORE": "/root/grm-research/grm-research-mtb-scaleup-r1/fasta",
    "CLEANED_PREFIX_SHA256": "b104d8bd7ffaa2d0dc3032fe30ca0c0b32670a7099a02bcd49860e5c77f59a64",
    "CLEANED_PREFIX_COUNT": 2000,
    "CLEANED_ADDED_SHA256": "51199cd34b123bca55a0dee2827b398a650774bfbef3df09d2264b79f617b765",
    "CLEANED_ADDED_COUNT": 2000,
    "CLEANED_FULL_SHA256": "f2a18a43d1743b46fbbb7f1b136a4921ba9e884f19f7942c131716a90845fd75",
    "CLEANED_FULL_COUNT": 4000,
    "CLEANED_MATCHED2K_SHA256": "efe0e90769b460e8d5949ab38881f24f8b137d29912d3ec3dfb7eaeb870be9ad",
    "CLEANED_MATCHED2K_COUNT": 2000,
    "RN_REFERENCE_ROOT": "/root/grm-r10-275-run-20261030",
}

# The input root's files the cohorts' labels, gates, gene calls and fingerprints come from, each bound by hash
# (review P1; the re-check's N1 added the two catalog files grm matrix loads, which GRM identifies only by version.txt).
INPUT_ROOT_FILES = {"curated/datasets.json": "RN_INPUT_DATASETS_SHA256", "curated/phenotypes.jsonl": "RN_INPUT_PHENOTYPES_SHA256",
                    "catalog/version.txt": "RN_INPUT_CATALOG_VERSION_SHA256", "catalog/AMR_CDS.fa": "RN_INPUT_CATALOG_CDS_SHA256",
                    "catalog/ReferenceGeneCatalog.txt": "RN_INPUT_CATALOG_REFERENCE_SHA256"}
COHORT_PINS = {name: (f"CLEANED_{name.upper()}_SHA256", f"CLEANED_{name.upper()}_COUNT")
               for name in ("prefix", "added", "full", "matched2k")}
GROUPS = {g: [n for n, (_, group, _) in TABLE.items() if group == g] for g in ("prepare", "reference")}
CONTRACT_NAMES = [n for pair in COHORT_PINS.values() for n in pair]
# Values derived from the bound inputs, which preparation recomputes and requires (--derive-cleaned prints them).
DERIVED_NAMES = CONTRACT_NAMES + ["RN_INPUT_MANIFEST_SHA256"]
INPUT_NAMES = [n for n in GROUPS["prepare"] if n not in DERIVED_NAMES]
Unfilled, Malformed = rev_n.Unfilled, rev_n.Malformed


def checked(name, value):
    """The value, or Malformed when it is not of its kind."""
    kind, _, what = TABLE[name]
    integer = isinstance(value, int) and not isinstance(value, bool)
    ok = {
        "sha256": lambda v: isinstance(v, str) and bool(rev_n.SHA256.fullmatch(v)),
        "count": lambda v: integer and v >= 1,
        "exclusions": lambda v: integer and v >= len(REQUIRED_EXCLUSIONS),
        "path": lambda v: isinstance(v, str) and v.startswith("/") and ".." not in PurePosixPath(v).parts
        and str(PurePosixPath(v)) == v,
        "format": lambda v: isinstance(v, str) and bool(v),
    }[kind](value)
    if not ok:
        raise Malformed(f"{name} = {value!r} is not a valid {kind} ({what})")
    if kind == "format" and value != SPECIES_CHECK_FORMAT_PARSED:
        raise Malformed(f"{name} is {value!r}, not {SPECIES_CHECK_FORMAT_PARSED!r}: the species-check schemas need a reviewed parser change first")
    return value


def require(names, values=None):
    """{name: value} for `names`; Unfilled naming every placeholder among them, else Malformed."""
    table = {n: (kind, what) for n, (kind, _, what) in TABLE.items()}
    return rev_n.require(names, VALUES if values is None else values, table, checked)


def main(argv=None):
    import argparse
    import json
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0], allow_abbrev=False)
    sub = ap.add_subparsers(dest="command", required=True)
    sub.add_parser("status").add_argument("--group", choices=sorted(GROUPS))
    sub.add_parser("get").add_argument("name", choices=sorted(TABLE))
    args = ap.parse_args(argv)
    try:
        if args.command == "status":
            names = GROUPS[args.group] if args.group else list(TABLE)
            missing = [n for n in names if rev_n.is_placeholder(VALUES.get(n, "@@X@@"))]
            print(json.dumps({n: n not in missing for n in names}, indent=1))
            require(names)
            return 0
        value = require([args.name])[args.name]
        print(value)
        return 0
    except (Unfilled, Malformed) as exc:
        print(f"cleaned_rerun: STOP: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
