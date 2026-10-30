#!/usr/bin/env python3
"""Prepare the matched MTB roots. NOT EXECUTED by the wiring or fix lanes.

usage: python3 -I prepare_execution_roots.py RUN_ROOT --cohorts frozen|cleaned
       python3 -I prepare_execution_roots.py --derive-cleaned
       python3 -I prepare_execution_roots.py --verify-frozen RUN_ROOT
       python3 -I prepare_execution_roots.py --check-lane COHORT_DATA_ROOT

Cohorts, each its own data root under RUN_ROOT/cohorts/<name>/data. With --cohorts frozen
(plan revision 5, the current run):

  pilot256        rows 1-256 of the frozen order (the small biological pilot)
  pilot256-stock  the same 256, for the stock count + matrix comparator
  prefix          rows 1-2000 (the frozen 2k candidate, byte for byte)
  added           rows 2001-4000
  full            all 4,000
  matched2k       the size-matched control: the 2,000 of the 4,000 with the
                  smallest SHA-256 of MATCHED_SEED + ":" + genome_id, in
                  frozen order (fixed here, before any result exists)

With --cohorts cleaned (plan revision 6, the cleaned-cohort rerun on revision N), every value
comes from cleaned_rerun.py, and the script refuses, before it creates anything, while one it
reads is an unfilled placeholder. Each input is read only at its bound SHA-256: the post-fetch
species check (species_check.json), revision N's MTB 2k and cap-4000 selections (genomes.txt
form) with each selection's curated/species.json, a FASTA manifest for the cohort genomes the
frozen 4k list does not hold, and RN_INPUT_ROOT's curated/datasets.json, curated/phenotypes.jsonl,
catalog/version.txt, catalog/AMR_CDS.fa and catalog/ReferenceGeneCatalog.txt. RN_INPUT_MANIFEST_SHA256
binds the sorted path<TAB>sha256 manifest of every input-root file a cohort root can read
(input_root_files: catalog/ and curated/ whole, the two record files, the cohort genomes' cached
annotations); --derive-cleaned prints it with the CLEANED_* pins. Then, by the admin's rulings of 2026-10-29:

  prefix          revision N's own 2k selection (a fresh selection, not the frozen rows)
  full            revision N's cap-4000 selection, whose first rows must be exactly the prefix:
                  otherwise a stop, before anything is created
  added           the full selection's rows after the prefix
  matched2k       matched_selection of the full selection, at the prefix's size

The exclusions are the MTB genomes the species check fails; they must include 1773.41701 and
1773.40949 and number EXCLUSIONS_COUNT. Borderline genomes are recorded and never excluded. Every
full genome (so every cohort genome) must have a species-check verdict, none may fail or be in the
check's unchecked list, and each selection's species.json must report selected_unchecked 0 and
its length. Each list must equal its CLEANED_* pin (hash and count); --derive-cleaned prints those
pins from the bound inputs and writes nothing. FASTA come from the frozen manifest, or for a
replacement from the bound manifest under RN_FASTA_STORE, each hash-checked. catalog, curated,
the BV-BRC records and cached annotations come from RN_INPUT_ROOT. No pilot roots are made: the
section 4 pilot is not repeated. preparation.json records the exclusions, the borderline genomes
per cohort, the replacements, and how the prefix differs from the frozen rows 1-2000;
frozen-inputs.json binds every input above, and --verify-frozen re-hashes them. The run root also
gets fasta-manifest.tsv, the manifest of every cohort genome, for the duplicate screen.

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
import importlib
import json
import os
import shutil
import stat
import sys
from pathlib import Path

# M1: python3 -I (3.11+ implies -P) leaves this script's directory off sys.path. Put back that one
# directory, first, and require cleaned_rerun to load from it.
LANE = Path(__file__).resolve().parent
if sys.path[:1] != [str(LANE)]:
    sys.path.insert(0, str(LANE))


def _sibling(name):
    module = importlib.import_module(name)
    if Path(module.__file__).resolve().parent != LANE:
        raise ImportError(f"{name} loaded from {module.__file__}, not from {LANE}")
    return module


try:
    rerun = _sibling("cleaned_rerun")
except ImportError as _exc:
    print(f"prepare_execution_roots: STOP: cannot load cleaned_rerun: {_exc}", file=sys.stderr)
    sys.exit(2)

EVIDENCE = Path("/root/grm-research/grm-research-mtb-4k-cohort-r2")
CANDIDATE = Path("/root/grm-extra-species-analysis-20261027/data")
SLUG = "mycobacterium_tuberculosis"
MATCHED_SEED = "grm-mtb-4k-matched2k-20261028"
PILOT = 256
PREFIX_ROWS = 2000  # the frozen 2k candidate's rows
VERDICTS = ("pass", "borderline", "fail")
MANIFEST_COLUMNS = ["genome_id", "size_bytes", "sha256", "canonical_source_path"]
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


def matched_selection(full_ids, size=None):
    """The size-matched control: `size` genomes of the full cohort (half by
    default; the prefix's size), chosen by a hash fixed before any result,
    kept in frozen order."""
    size = len(full_ids) // 2 if size is None else size
    rank = sorted(full_ids, key=lambda g: hashlib.sha256(f"{MATCHED_SEED}:{g}".encode()).hexdigest())
    chosen = set(rank[:size])
    return [g for g in full_ids if g in chosen]


def refuse(message):
    raise SystemExit(f"refusing: {message}")


def id_list(data: bytes, what: str):
    """Genome IDs, one per LF-terminated line: no blank, padded, duplicate or path-like entry."""
    try:
        text = data.decode()
    except UnicodeDecodeError:
        refuse(f"{what} is not UTF-8")
    if not text.endswith("\n") or "\r" in text:
        refuse(f"{what} must be LF-terminated lines")
    ids = text[:-1].split("\n")
    if len(set(ids)) != len(ids) or any(not g or g != g.strip() or Path(g).name != g or g in {".", ".."} for g in ids):
        refuse(f"{what} holds a blank, padded, duplicate or path-like genome ID")
    return ids


def read_bound(path, digest: str, what: str) -> bytes:
    """The file's bytes, only when they hash to the bound SHA-256."""
    path = Path(path)
    if not path.is_file():
        refuse(f"{what} {path} is missing")
    data = path.read_bytes()
    if sha256_bytes(data) != digest:
        refuse(f"{what} {path} has SHA-256 {sha256_bytes(data)}, not the bound {digest}")
    return data


def species_verdicts(data: bytes):
    """{genome_id: verdict} for MTB from the species check's species_check.json, the check's MTB
    unchecked list, and its MTB summary. Any shape the parser does not know is a refusal."""
    try:
        doc = json.loads(data)
    except ValueError:
        refuse("the species check's output is not JSON")
    if not isinstance(doc, dict) or not isinstance(doc.get("species"), list) or not isinstance(doc.get("genomes"), list):
        refuse("the species check's output lacks its species and genomes lists")
    entries = [s for s in doc["species"] if isinstance(s, dict) and s.get("slug") == SLUG]
    if len(entries) != 1 or not isinstance(entries[0].get("species"), str) or not isinstance(entries[0].get("unchecked"), list):
        refuse("the species check's output has no single MTB entry with a species name and an unchecked list")
    entry = entries[0]
    verdicts = {}
    for g in doc["genomes"]:
        if not isinstance(g, dict):
            refuse("the species check's output holds a genome row that is not an object")
        if g.get("species") != entry["species"]:
            continue
        gid, verdict = g.get("genome_id"), g.get("verdict")
        if not isinstance(gid, str) or not gid or gid in verdicts or verdict not in VERDICTS:
            refuse(f"the species check's MTB row {gid!r} has a missing, repeated or unknown verdict {verdict!r}")
        verdicts[gid] = verdict
    for key, verdict in (("failed", "fail"), ("borderline", "borderline")):
        if entry.get(key) != sum(v == verdict for v in verdicts.values()):
            refuse(f"the species check's MTB summary {key}={entry.get(key)!r} disagrees with its {verdict} rows")
    return verdicts, list(entry["unchecked"]), {k: v for k, v in entry.items() if not isinstance(v, list)}


def species_selection(data: bytes, what: str, selected: int):
    """The MTB entry of a selection's curated/species.json: selected_unchecked must be 0 and
    selected_genomes the selection's length."""
    try:
        doc = json.loads(data)
    except ValueError:
        refuse(f"{what} is not JSON")
    entries = [s for s in doc if isinstance(s, dict) and s.get("slug") == SLUG] if isinstance(doc, list) else []
    if len(entries) != 1:
        refuse(f"{what} is not a list with one MTB selection")
    unchecked = entries[0].get("selected_unchecked")
    if not isinstance(unchecked, int) or isinstance(unchecked, bool) or unchecked != 0:
        refuse(f"{what}: selected_unchecked is {unchecked!r}, not 0, so a selected replacement has not been species-checked")
    if entries[0].get("selected_genomes") != selected:
        refuse(f"{what}: selected_genomes is {entries[0].get('selected_genomes')!r}, not the selection's {selected}")
    return entries[0]


def cleaned_cohorts(frozen_full, verdicts, unchecked, prefix, full):
    """The cleaned prefix, added, full and matched2k lists, and what preparation.json records."""
    if full[:len(prefix)] != prefix:
        refuse("the cleaned 2k is not the first rows of the cleaned 4k: revision N's 2k selection is not a prefix of its "
               "cap-4000 selection, so prefix, added and full cannot be built; stop and bring it back to the admin")
    fails = [g for g, v in verdicts.items() if v == "fail"]
    absent = [g for g in rerun.REQUIRED_EXCLUSIONS if verdicts.get(g) != "fail"]
    if absent:
        refuse(f"the species check does not fail {absent}, which the 2026-10-29 ruling excludes")
    missing = [g for g in full if g not in verdicts]
    if missing:
        refuse(f"{len(missing)} cohort genomes have no verdict in the bound species check: {missing[:10]}")
    excluded = [g for g in full if verdicts[g] == "fail"]
    if excluded:
        refuse(f"excluded genomes are in the cohorts: {excluded[:10]}")
    pending = sorted(set(full) & set(unchecked))
    if pending:
        refuse(f"cohort genomes are in the species check's unchecked list: {pending[:10]}")
    cohorts = {"prefix": list(prefix), "added": full[len(prefix):], "full": list(full),
               "matched2k": matched_selection(full, len(prefix))}
    frozen = set(frozen_full)
    frozen_prefix = frozen_full[:PREFIX_ROWS]
    info = {"exclusions": fails,
            "borderline": {name: [g for g in ids if verdicts[g] == "borderline"] for name, ids in cohorts.items()},
            "replacements": [g for g in full if g not in frozen],
            "prefix_vs_frozen_rows_1_2000": {"left": [g for g in frozen_prefix if g not in set(prefix)],
                                             "joined": [g for g in prefix if g not in set(frozen_prefix)],
                                             "same_order_on_common": [g for g in prefix if g in set(frozen_prefix)]
                                             == [g for g in frozen_prefix if g in set(prefix)]},
            "full_vs_frozen_4k": {"left": [g for g in frozen_full if g not in set(full)], "joined": [g for g in full if g not in frozen]}}
    return cohorts, info


def replacement_fasta(data: bytes, ids, store: Path):
    """{genome_id: {source, sha256}} from the bound manifest for the cohort genomes outside the frozen
    4k list: exactly those genomes, each a regular file directly under `store` with its manifest
    size and SHA-256."""
    rows = list(csv.DictReader(data.decode().splitlines(), delimiter="\t"))
    if not rows and ids:
        refuse("the replacement FASTA manifest is empty")
    if rows and not set(MANIFEST_COLUMNS) <= set(rows[0]):
        refuse(f"the replacement FASTA manifest lacks the columns {MANIFEST_COLUMNS}")
    got = [row["genome_id"] for row in rows]
    if len(set(got)) != len(got) or set(got) != set(ids):
        refuse(f"the replacement FASTA manifest covers {sorted(set(got))[:10]}, not exactly the replacements {sorted(ids)[:10]}")
    sources = {}
    for row in rows:
        source = Path(row["canonical_source_path"])
        if (source.parent != store or source.name != row["genome_id"] + ".fna.gz" or not source.is_file()
                or source.is_symlink() or source.resolve(strict=True) != source):
            refuse(f"replacement FASTA {source} is not a regular {row['genome_id']}.fna.gz directly under {store}")
        if source.stat().st_size != int(row["size_bytes"]) or sha256_file(source) != row["sha256"]:
            refuse(f"{source} differs from the replacement FASTA manifest")
        sources[row["genome_id"]] = {"source": str(source), "sha256": row["sha256"], "size_bytes": int(row["size_bytes"])}
    return sources


def input_root_files(input_root: Path, full_ids):
    """Every RN_INPUT_ROOT file a cleaned cohort root can read, by the rules that make the roots: all of
    catalog/ and curated/ (each linked whole), the two BV-BRC record files (linked) and the cached annotations
    of the cohort genomes (copied by prepare_bvbrc; every cohort is part of full)."""
    files = []
    for sub in ("catalog", "curated"):
        files += [p for p in (input_root / sub).rglob("*") if p.is_file()]
    files += [input_root / "bvbrc" / name for name in BVBRC_RECORDS]
    files += cached_features(input_root, full_ids)
    return sorted(files)


def input_root_manifest(input_root: Path, full_ids):
    """The sorted `path<TAB>sha256` manifest of input_root_files (paths relative to the root), and {path: sha256}."""
    hashes = {str(p): sha256_file(p) for p in input_root_files(input_root, full_ids)}
    lines = sorted(f"{Path(p).relative_to(input_root)}\t{h}\n" for p, h in hashes.items())
    return "".join(lines).encode(), hashes


def derived_pins(cohorts, input_root: Path):
    """{name: value} of every DERIVED_NAMES value: the four cleaned lists and the input-root manifest."""
    manifest, _ = input_root_manifest(input_root, cohorts["full"])
    return {**cohort_pins(cohorts), "RN_INPUT_MANIFEST_SHA256": sha256_bytes(manifest)}


def load_cleaned(frozen_full):
    """The cleaned cohorts from the bound inputs, the values that bind them, and what preparation.json records."""
    try:
        v = rerun.require(rerun.INPUT_NAMES)
    except (rerun.Unfilled, rerun.Malformed) as exc:
        refuse(str(exc))
    bound = {}
    input_root = Path(v["RN_INPUT_ROOT"])
    for rel, sha_name in rerun.INPUT_ROOT_FILES.items():  # the labels and catalog the cohorts are learned from (review P1)
        bound[str(input_root / rel)] = v[sha_name]
        read_bound(input_root / rel, v[sha_name], f"RN_INPUT_ROOT's {rel}")

    def read(path_name, sha_name, what):
        bound[v[path_name]] = v[sha_name]
        return read_bound(v[path_name], v[sha_name], what)

    verdicts, unchecked, summary = species_verdicts(read("SPECIES_CHECK_PATH", "SPECIES_CHECK_SHA256", "the species check's output"))
    prefix = id_list(read("RN_PREFIX_GENOMES_PATH", "RN_PREFIX_GENOMES_SHA256", "revision N's 2k selection"), "revision N's 2k selection")
    full = id_list(read("RN_FULL_GENOMES_PATH", "RN_FULL_GENOMES_SHA256", "revision N's cap-4000 selection"), "revision N's cap-4000 selection")
    species_selection(read("RN_PREFIX_SPECIES_JSON_PATH", "RN_PREFIX_SPECIES_JSON_SHA256", "the 2k selection's species.json"),
                      "the 2k selection's species.json", len(prefix))
    species_selection(read("RN_FULL_SPECIES_JSON_PATH", "RN_FULL_SPECIES_JSON_SHA256", "the cap-4000 selection's species.json"),
                      "the cap-4000 selection's species.json", len(full))
    cohorts, info = cleaned_cohorts(frozen_full, verdicts, unchecked, prefix, full)
    if len(info["exclusions"]) != v["EXCLUSIONS_COUNT"]:
        refuse(f"the species check fails {len(info['exclusions'])} MTB genomes, not EXCLUSIONS_COUNT {v['EXCLUSIONS_COUNT']}")
    store = Path(v["RN_FASTA_STORE"])
    fasta = replacement_fasta(read("RN_FASTA_MANIFEST_PATH", "RN_FASTA_MANIFEST_SHA256", "the replacement FASTA manifest"),
                              info["replacements"], store)
    # The re-check's N1: every input-root file the cohort roots can read, recorded for --verify-frozen.
    manifest, hashes = input_root_manifest(input_root, cohorts["full"])
    bound.update(hashes)
    record = {"bound": bound, "species_check_format": v["SPECIES_CHECK_FORMAT"], "species_check_mtb": summary,
              "input_root": v["RN_INPUT_ROOT"], "fasta_store": str(store),
              "input_manifest": {"sha256": sha256_bytes(manifest), "files": len(hashes)}}
    return cohorts, {**record, **info}, fasta


def cohort_pins(cohorts):
    """{CLEANED_* name: value} for the cleaned lists."""
    out = {}
    for name, (sha_name, count_name) in rerun.COHORT_PINS.items():
        out[sha_name] = sha256_bytes("".join(f"{g}\n" for g in cohorts[name]).encode())
        out[count_name] = len(cohorts[name])
    return out


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
    cleaned = recorded.get("cleaned")
    if cleaned is not None:
        # The cleaned rerun's bound inputs: the species check, both selections, their species.json and the replacement manifest.
        for path, digest in sorted(cleaned["bound"].items()):
            p = Path(path)
            if not p.is_file() or sha256_file(p) != digest:
                changed.append(path)
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
    # The input root and FASTA stores this run's preparation recorded (a cleaned run's are revision N's), else the frozen ones.
    candidate, parents = CANDIDATE, set(ALLOWED_FASTA_PARENTS)
    record = data.parent.parent.parent / "frozen-inputs.json"
    if data.parent.parent.name == "cohorts" and record.is_file():
        recorded = json.loads(record.read_text())
        candidate = Path(recorded["candidate"])
        parents |= {Path(p) for p in recorded.get("fasta_parents", [])}
    allowed = {data / "catalog": candidate / "catalog", data / "curated": candidate / "curated"}
    for name in BVBRC_RECORDS:
        allowed[data / "bvbrc" / name] = candidate / "bvbrc" / name
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
                    if target.parent not in parents or target.name != name:
                        problems.append(f"{p} links outside the frozen FASTA stores: {target}")
                else:
                    problems.append(f"{p} is a link where GRM reads or writes lane files")
            elif not (stat.S_ISDIR(mode) or stat.S_ISREG(mode)):
                problems.append(f"{p} is neither a directory nor a regular file")
    return problems


def frozen_cohorts():
    """The frozen lists from the cohort evidence, each checked against EXPECTED: ({name: ids}, {name: bytes})."""
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
    return cohort_ids, cohort_bytes


def main() -> int:
    parser = argparse.ArgumentParser(allow_abbrev=False)
    parser.add_argument("path", type=Path, nargs="?", help="RUN_ROOT, or a cohort data root with --check-lane")
    parser.add_argument("--cohorts", choices=["frozen", "cleaned"], help="which cohorts to prepare (required to prepare)")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--verify-frozen", action="store_true")
    mode.add_argument("--check-lane", action="store_true")
    mode.add_argument("--derive-cleaned", action="store_true", help="print the CLEANED_* pins from the bound inputs; write nothing")
    args = parser.parse_args()
    if args.derive_cleaned:
        if args.path is not None or args.cohorts is not None:
            parser.error("--derive-cleaned takes no RUN_ROOT or --cohorts")
        try:
            rerun.require(rerun.INPUT_NAMES)
        except (rerun.Unfilled, rerun.Malformed) as exc:
            refuse(str(exc))
        cleaned_ids, info, _ = load_cleaned(frozen_cohorts()[0]["full"])
        pins = derived_pins(cleaned_ids, Path(info["input_root"]))
        print(json.dumps({"pins": pins, **{k: v for k, v in info.items() if k != "bound"}}, indent=1))
        return 0
    if args.path is None or (args.cohorts is not None) == (args.verify_frozen or args.check_lane):
        parser.error("give RUN_ROOT with --cohorts to prepare, or a path with --verify-frozen or --check-lane")
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
    if args.cohorts == "cleaned":  # every cleaned value, before any input is read
        try:
            rerun.require(rerun.INPUT_NAMES + rerun.DERIVED_NAMES)
        except (rerun.Unfilled, rerun.Malformed) as exc:
            refuse(str(exc))

    cohort_ids, cohort_bytes = frozen_cohorts()
    frozen_full = cohort_ids["full"]
    pinned_matched = EXPECTED["matched2k"][0]
    candidate, fasta_parents, cleaned, replacement_sources = CANDIDATE, set(ALLOWED_FASTA_PARENTS), None, {}
    if args.cohorts == "cleaned":
        cleaned_ids, cleaned, replacement_sources = load_cleaned(frozen_full)
        try:
            pins = rerun.require(rerun.DERIVED_NAMES)
        except (rerun.Unfilled, rerun.Malformed) as exc:
            refuse(str(exc))
        derived = derived_pins(cleaned_ids, Path(cleaned["input_root"]))
        differ = sorted(n for n in derived if derived[n] != pins[n])
        if differ:
            refuse(f"the cleaned cohorts or the input-root manifest differ from their pins {differ}: derived {[derived[n] for n in differ]}")
        # The section 4 pilot is not repeated: a cleaned run root holds the four training cohorts only.
        cohort_ids = dict(cleaned_ids)
        cohort_bytes = {name: "".join(f"{g}\n" for g in ids).encode() for name, ids in cohort_ids.items()}
        pinned_matched = pins["CLEANED_MATCHED2K_SHA256"]
        candidate = Path(cleaned["input_root"])
        fasta_parents.add(Path(cleaned["fasta_store"]))

    with (EVIDENCE / "results/fasta-manifest.tsv").open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 4000
    assert [row["genome_id"] for row in rows] == frozen_full
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
        sources[row["genome_id"]] = {"source": str(source), "sha256": row["sha256"], "size_bytes": int(row["size_bytes"])}
    sources.update(replacement_sources)
    linked = {g for ids in cohort_ids.values() for g in ids}
    if args.cohorts == "cleaned":
        sources = {g: rec for g, rec in sources.items() if g in linked}  # every cohort genome, and only those
    assert linked <= set(sources), sorted(linked - set(sources))[:10]

    # The annotation copies are the only bytes this script writes in bulk.
    copies = {name: cached_features(candidate, cohort_ids[name]) for name in TRAINING}
    copy_bytes = sum(p.stat().st_size for files in copies.values() for p in files)
    free_before = free_bytes(run_root.parent)
    if free_before - copy_bytes < MIN_FREE_AFTER:
        raise SystemExit(f"refusing: copying {copy_bytes} annotation bytes would leave {free_before - copy_bytes} bytes free, "
                         f"below the {MIN_FREE_AFTER}-byte floor")

    manifest = frozen_manifest(candidate)
    run_root.mkdir()
    (run_root / "frozen-inputs.json").write_text(json.dumps({
        "candidate": str(candidate), "sha256": manifest, "fasta": sources,
        "cohort_genomes_sha256": {name: sha256_bytes(data) for name, data in cohort_bytes.items()},
        "cohort_selection": args.cohorts, "fasta_parents": sorted(str(p) for p in fasta_parents), "cleaned": cleaned,
    }, indent=1) + "\n")
    with (run_root / "fasta-manifest.tsv").open("x", newline="") as handle:  # every cohort genome, for the duplicate screen
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(MANIFEST_COLUMNS)
        for gid in cohort_ids["full"]:
            writer.writerow([gid, sources[gid]["size_bytes"], sources[gid]["sha256"], sources[gid]["source"]])
    seeded = {}
    for name, ids in cohort_ids.items():
        data = run_root / "cohorts" / name / "data"
        species = data / "species" / SLUG
        fasta = data / "fasta"
        species.mkdir(parents=True)
        fasta.mkdir()
        for directory in ["catalog", "curated"]:
            target = candidate / directory
            assert target.is_dir()
            os.symlink(target, data / directory, target_is_directory=True)
        seeded[name] = prepare_bvbrc(data, candidate, manifest, ids if name in TRAINING else [])
        with (species / "genomes.txt").open("xb") as handle:
            handle.write(cohort_bytes[name])
        for genome_id in ids:
            os.symlink(sources[genome_id]["source"], fasta / f"{genome_id}.fna.gz")

    for directory in ["logs", "backups", "audit", "target"]:
        (run_root / directory).mkdir()
    prefix = set(cohort_ids["prefix"])
    receipt = {
        "run_root": str(run_root),
        "cohort_selection": args.cohorts,
        "cohorts": {name: len(ids) for name, ids in cohort_ids.items()},
        "genome_sha256": {name: sha256_bytes(data) for name, data in cohort_bytes.items()},
        "matched2k": {"seed": MATCHED_SEED, "from_prefix": sum(g in prefix for g in cohort_ids["matched2k"]),
                      "from_added": sum(g not in prefix for g in cohort_ids["matched2k"]),
                      "sha256": sha256_bytes(cohort_bytes["matched2k"]), "frozen_sha256": pinned_matched},
        "cleaned": None if cleaned is None else {k: v for k, v in cleaned.items() if k != "bound"},
        "fasta_links": sum(len(ids) for ids in cohort_ids.values()),
        "fasta_bytes_copied": 0,
        "fasta_sha256_checked": len(sources),
        "input_root": str(candidate),
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
