#!/usr/bin/env python3
"""Full-pipeline input comparison: every arm through `grm`, then J and K.

One command runs a comparison spec (spec.json beside this file holds the
decided arms, DECISIONS V-1): each arm changes one input of the default arm,
goes through `grm curate`, `grm count --k`, `grm structure`, `grm matrix` and
`grm arm` (blocks, encoders, the learners' screen and every learner, once per
outer fold of the reference plans), and its prediction rows are paired with
the default arm's by J (`scripts/ablation/stats.py`) under the arm's pairing
rule, then reported by K (`scripts/ablation/report.py`).

An arm is a data root, so arms share a root when they share a cohort and a k:
blocks and budget arms score on the default root; a k arm counts the default
cohort again at its k; a cap or selection arm curates its own cohort. Each
root links the raw inputs of --raw (bvbrc/, fasta/, catalog/ and species-check/) and
curates with --sigdrug-vocabulary, as revision 12's S0 did.

Pairing follows from the factor an arm changes:
    cap, selection   common_outer_test_set_balanced (another cohort)
    k                k31_groups (the default cohort, the default's k = 31 plans)
    blocks, budget   same_plan
An arm that changes no factor, or more than one, is refused.

Decisions come from the lineage split only. Random-split contrasts are planned in family D:
J gives their estimate, interval and AUROC delta but no p and no direction, and K prints them
as descriptive (main's ruling after the calibration study, 5 Nov 2026).

Every data root puts its genomes in canonical order (sorted by genome id) after `grm curate` and
before count, structure, matrix, fold plans and evaluation (REVIEW-INPUT M1), so two arms with the
same genome set give bit-identical predictions whatever order their selection rule drew them in.
A cohort arm whose genome set equals the default's is not run: its runs are `identical_cohort`
with the reason, and K reports it with no contrasts (`--run-identical-cohorts` runs it anyway,
for the order-twin proof only). The spec's `panel_floor` sets the eligible genomes a species
needs (more than 2,000, so cap-4000 and grr-v1 change the set); the runner refuses a species
below it unless `--below-panel-floor` is given (synthetic studies only).

Every spec run gets a sentinel arm, `sentinel-seed-twin`: the default's learner-seed twin
(learner seed 1001), which is exchangeable with the default. The manifest names it in
`sentinel_arms`. When its lineage-split contrast for a learner is anything but
inconclusive, or isn't produced, K marks every decision of that learner on the species
unreliable. summary.json then shows those decisions as "unreliable" (the decision itself is in
`raw_decisions`) and leaves them out of `any_improves` and `any_worse`.

Calibration arms (calibration.json) change no input: they are the default made
to look like an arm, so a decision on one is a false positive of the decision
procedure. An arm with "calibration": "cohort_path" sends the default's inputs
through the cohort-arm path (common_outer_test_set_balanced, with its own
"balance_seed"); "calibration": "learner_seed" is a same_plan twin of the
default whose learners draw from another "learner_seed".

Kset sharing (H-1 gap G1): every genome of a root counted at the pool's k is a pool genome, so its kset is the
pool's, byte for byte. After the pool is counted, each such root gets a hard link to the pool's finished kset and
sidecar for every genome it doesn't already hold, `grm count` reuses them, and the disk check counts them once.
After the count every linked kset must still be the pool's file (a recount would write a new one); the root's
shared.json records how many of its ksets are the pool's and how many it counted itself. --no-share-ksets counts
every root's own, as before, for the identity proof only.

Resume (G6): a reused WORK skips an arm on a dataset when its resume.json, written after the arm finished, still
matches: the same inputs.json, spec, species grouping and reference plans (the default arm's plans for that
dataset), the same --trace-training and --run-identical-cohorts, and every output file with its recorded sha256.
Anything else reruns the arm.

`run.py combine` (G2) makes one J record set and one K report from runs of one spec split by datasets: each part is
a WORK whose spec equals the full spec but for its datasets. The parts' datasets must be disjoint and cover the
full spec's, the parts must share the species pool, the grm binary, the vocabulary and the source commit, and every
RAW file each part hashed must be in the pinned master manifest. J's family C is then one BH family per species
over every part's lineage contrasts. Each copied arm output must match its part's resume.json, and
the spec's data_root.manifest_sha256 pins the master manifest.

usage: python3 -B run.py SPEC --raw RAW_ROOT --work WORK --grm GRM_BINARY --source-commit SHA
           [--jobs N] [--estimate-only] [--free-disk-bytes N] [--below-panel-floor] [--trace-training]
           [--no-share-ksets]
       python3 -B run.py combine SPEC --work-from WORK [--work-from WORK ...] --master-manifest FILE --out DIR

Writes WORK/roots/<root>/ (data roots), WORK/arms/<arm>/<dataset>/ (`grm arm`
output), WORK/j/{contrasts,missing,planned}.json, WORK/k/manifest.json and
WORK/report/ablation-report.{json,md}, plus WORK/summary.json: each arm's
decision per learner and split. WORK/spec.json and WORK/raw-manifest.txt keep
the spec and RAW's file list for combine. Needs numpy (J). Reads no network:
every genome's FASTA must already be in RAW_ROOT/fasta.
"""
import argparse
import concurrent.futures
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "ablation"))
import stats  # noqa: E402  J
sys.path.insert(0, str(ROOT))
from scripts.ablation import report as k_report  # noqa: E402  K
from scripts.ablation import schema as k_schema  # noqa: E402

FACTORS = {"cohort": ("cap", "selection", "selection_seed", "stratum"), "k": ("k",),
           "blocks": ("blocks",), "budget": ("screened_features",)}
PAIRING = {"cohort": "common_outer_test_set_balanced", "k": "k31_groups", "blocks": "same_plan",
           "budget": "same_plan"}
DEFAULT_ARM = "default"
SENTINEL = {"name": "sentinel-seed-twin", "calibration": "learner_seed", "learner_seed": 1001}
# Calibration arm kind -> (its seed key, its pairing rule).
CALIBRATION = {"cohort_path": ("balance_seed", "common_outer_test_set_balanced"),
               "learner_seed": ("learner_seed", "same_plan")}
LEARNERS = ("scm", "cart", "rf", "gbt", "lasso")
SPLITS = ("random", "lineage")


def child_env():
    """The environment every grm child gets: this one without GRM_SIGDRUG_VOCABULARY, so a vocabulary reaches curate
    only through --sigdrug-vocabulary, which inputs.json hashes (REVIEW-INPUT8 L3)."""
    return {k: v for k, v in os.environ.items() if k != "GRM_SIGDRUG_VOCABULARY"}


def datasets_view(path):
    """The sha256 of a datasets.json with every entry's `selection` record removed, serialized as grm writes it
    (two-space indent, no trailing newline): revision 12's bytes from this tree's curate."""
    entries = json.loads(Path(path).read_text(encoding="utf-8"))
    for e in entries:
        e.pop("selection", None)
    return sha256_bytes(json.dumps(entries, indent=2, ensure_ascii=False).encode())


# A dataset's eligible-pool counts, which the species check moves without changing what the default learns (main's
# ruling after H1DATA): left out of the gate's dataset spec, as `selection` is.
COUNT_FIELDS = ("n_resistant", "n_susceptible", "n_intermediate")


def substance(root, slugs):
    """What the revision-12 gate compares of a curated root (main's ruling after H1DATA, the rule of H1DATA's
    tools/substantive_check.py), for each species in `slugs`:
    - `genome_set`: the sha256 of its species/<slug>/genomes.txt sorted, one id per line (a duplicate id changes it);
    - `datasets`: for each of its datasets in curated/datasets.json, the sha256 of the canonical JSON of
      {"spec": the entry without COUNT_FIELDS and `selection`, "rows": the curated/phenotypes.jsonl row of every
      selected genome for the dataset's antibiotic, by genome id}. The whole row fixes the learned label
      (`Phenotype::reported_label`: reported call, conflict, label_from_mic) and the rest of what curate recorded.
    Revision 12's values come from its master root through this same function."""
    root = Path(root)
    entries = [d for d in json.loads((root / "curated" / "datasets.json").read_text(encoding="utf-8"))
               if d["species_slug"] in slugs]
    ids = {slug: sorted(line.strip() for line in (root / "species" / slug / "genomes.txt")
                        .read_text(encoding="utf-8").splitlines() if line.strip()) for slug in slugs}
    names = {d["species_slug"]: d["species"] for d in entries}
    want = {(names[slug], g) for slug in slugs if slug in names for g in ids[slug]}
    rows = {}
    with open(root / "curated" / "phenotypes.jsonl", encoding="utf-8") as f:
        for line in f:
            p = json.loads(line)
            if (p["species"], p["genome_id"]) in want:
                key = (p["genome_id"], p["antibiotic"])
                if key in rows:
                    raise SystemExit(f"{root}: two phenotype rows for genome {key[0]} and {key[1]}")
                rows[key] = p
    out = {"genome_set": {slug: sha256_bytes("".join(g + "\n" for g in ids[slug]).encode()) for slug in slugs},
           "datasets": {slug: {} for slug in slugs}}
    for d in entries:
        spec = {k: v for k, v in d.items() if k not in COUNT_FIELDS + ("selection",)}
        selected = [rows[(g, d["antibiotic"])] for g in ids[d["species_slug"]] if (g, d["antibiotic"]) in rows]
        out["datasets"][d["species_slug"]][d["id"]] = sha256_bytes(canonical({"spec": spec, "rows": selected}).encode())
    return out


def gate_line(gate, spec):
    """The revision-12 gate's decision as the run prints it on stderr at the start, before anything is counted
    (main's ruling on REVIEW-INTEG7 L2): passed, with what it compared, or not required on a synthetic RAW."""
    if gate is None:
        return ("revision-12 gate: NOT REQUIRED: the spec says \"synthetic\": true and RAW is a synthetic root "
                "(truth.json); this run is ungated")
    species = ", ".join(sorted(gate["genome_set"]))
    n = sum(len(v) for v in gate["datasets"].values())
    recorded = ("as the spec records them" if gate["recorded"]["as_the_spec_records"]
                else "NOT as the spec records them (see the gate record)")
    return (f"revision-12 gate: PASSED: the default's genome sets for {species} and the spec and selected genomes' "
            f"phenotype rows of their {n} datasets equal revision 12's; recorded, not gated: the drawn order and the "
            f"datasets.json view, {recorded}"
            + (" (the spec also says synthetic; the gate ran)" if spec.get("synthetic") is True else ""))


def disk_bytes(path, skip=()):
    """The disk a folder's files take, in allocated blocks as df counts them; symlinks are not followed, files named
    in `skip` are left out, and a file with several hard links in the folder (a shared kset, G1) counts once."""
    blocks = {}
    for p in Path(path).rglob("*"):
        if p.is_file() and not p.is_symlink() and p.name not in skip:
            st = p.lstat()
            blocks[(st.st_dev, st.st_ino)] = st.st_blocks * 512
    return sum(blocks.values())


def file_sha256(path):
    """A file's sha256, read in 1 MiB blocks."""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(2 ** 20), b""):
            h.update(block)
    return h.hexdigest()


def sha256_bytes(data):
    return hashlib.sha256(data).hexdigest()


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True, allow_nan=False)


def load_spec(path):
    spec = json.loads(Path(path).read_text(encoding="utf-8"))
    if spec.get("format") != "grm-input-comparison-spec" or spec.get("version") != 1:
        raise SystemExit(f"{path}: not a grm-input-comparison-spec version 1")
    if "data_root" in spec:
        pinned = (spec["data_root"] or {}).get("manifest_sha256")
        if not isinstance(pinned, str) or not re.fullmatch(r"[0-9a-f]{64}", pinned):
            raise SystemExit(f"{path}: data_root.manifest_sha256 is {pinned!r}, not a sha256: pin the data root's "
                             "manifest (main's ruling after REVIEW-H1DATA)")
    default = spec["default"]
    datasets = [d["id"] for d in spec["datasets"]]
    if len(set(datasets)) != len(datasets):
        raise SystemExit(f"{path}: a dataset is listed twice: {datasets}")
    for ds in datasets:
        if not NAME.fullmatch(ds):
            raise SystemExit(f"{path}: dataset id {ds!r} is not a plain name ({NAME.pattern})")
    arms = []
    names = {DEFAULT_ARM}
    for arm in spec["arms"]:
        name = arm["name"]
        if not NAME.fullmatch(name):
            raise SystemExit(f"{path}: arm name {name!r} is not a plain name ({NAME.pattern})")
        if name in names:
            raise SystemExit(f"{path}: arm {name} appears twice")
        names.add(name)
        if "calibration" in arm:
            kind = arm["calibration"]
            seed_key = CALIBRATION.get(kind, (None,))[0]
            if seed_key is None or set(arm) != {"name", "calibration", seed_key} \
                    or type(arm[seed_key]) is not int:
                raise SystemExit(f"{path}: calibration arm {name} must be {{name, calibration, <seed>}} with "
                                 f"calibration one of {sorted(CALIBRATION)} and an integer seed")
            arms.append({**default, **arm, "factor": "calibration", "pairing": CALIBRATION[kind][1]})
            continue
        changes = {k for k in arm if k != "name"}
        unknown = changes - {k for keys in FACTORS.values() for k in keys}
        if unknown:
            raise SystemExit(f"{path}: arm {name} sets unknown inputs {sorted(unknown)}")
        changed = {f for f, keys in FACTORS.items() if any(k in arm and arm[k] != default.get(k) for k in keys)}
        if len(changed) != 1:
            raise SystemExit(f"{path}: arm {name} changes {sorted(changed) or 'nothing'}; each arm changes one input "
                             "of the default")
        factor = changed.pop()
        arms.append({**default, **arm, "factor": factor, "pairing": PAIRING[factor]})
        if factor == "cohort":
            arms[-1]["kind"] = cohort_kind(default, arms[-1])
    cohort_paired = [a["name"] for a in arms if a["pairing"] == "common_outer_test_set_balanced"]
    if spec.get("repeats", 1) != 1 and cohort_paired:
        raise SystemExit(f"{path}: repeats {spec['repeats']} with a cohort-paired arm ({', '.join(cohort_paired)}): a "
                         "cohort arm's balanced subsample can differ by repeat, and J pairs one sample set per arm "
                         "across repeats; every decided spec uses repeats 1 (PR #69's review)")
    unknown = set(spec.get("base_rates", {})) - set(COHORT_KINDS)
    if unknown:
        raise SystemExit(f"{path}: base_rates for {sorted(unknown)}; the kinds are {list(COHORT_KINDS)}")
    if SENTINEL["name"] in names:
        raise SystemExit(f"{path}: the arm name {SENTINEL['name']} is the runner's sentinel")
    arms.append({**default, **SENTINEL, "factor": "calibration", "pairing": "same_plan", "sentinel": True})
    return spec, {**default, "name": DEFAULT_ARM, "factor": None, "pairing": "same_plan"}, arms


# The comparison's disk check (main's rulings on REVIEW-INPUT4 to REVIEW-INPUT6 M1) when the spec states none: a
# k31 kset's size, the share of free disk the run may take, the matrix and sidecar margin on each data root, and
# the arm outputs (WORK/arms) as a share of the roots' ksets. All are stated figures until H-1 measures real
# per-species ksets and matrices.
POOL_BUDGET = {"default_kset_bytes": 27.5 * 2 ** 20, "max_disk_share": 0.8, "root_overhead_share": 0.10,
               "arm_output_share": 0.05, "work_files_bytes": 64 * 2 ** 20}
# What a WORK is tied to: its RAW input manifest, its grm binary and its SigDrug vocabulary (main's rulings on
# REVIEW-INPUT6 and REVIEW-INPUT7 L3).
INPUTS = "inputs.json"
# The spec a WORK ran (canonical JSON, so its sha256 is K's protocol_sha256) and the RAW file list inputs.json hashes:
# what `combine` checks a part against.
SPEC_COPY = "spec.json"
RAW_MANIFEST = "raw-manifest.txt"
# Per arm and dataset: what its outputs were made from (G6); per data root: its ksets shared with the pool (G1).
RESUME = "resume.json"
SHARED = "shared.json"
ARM_OUTPUTS = {"runs.json", "predictions.jsonl", "arm.json"}
# The RAW folders a data root links; species-check/ only when RAW has one.
RAW_PARTS = ("bvbrc", "fasta", "catalog", "species-check")
# Provenance side files (docs/provenance-fields.md) that no arm reads: `grm fetch genome-dates`' in RAW, left out of
# its input manifest, and `grm curate`'s, removed from each root after it curates and never copied or sized. Curate's
# is about 310 MB on the real root, once per data root.
RAW_SIDE_FILES = ("bvbrc/genome_dates.jsonl", "bvbrc/genome_dates.manifest.json")
CURATED_SIDE_FILES = ("provenance.jsonl",)
# synthetic_root.py's truth.json format: the marker a "synthetic": true spec needs in RAW (REVIEW-INTEG7 L2).
SYNTHETIC_FORMAT = "grm-input-comparison-synthetic"
# Arm and dataset names become path components under WORK (REVIEW-INPUT7 L2).
NAME = re.compile(r"[a-z0-9][a-z0-9._-]{0,63}")
# A cohort arm's kind, for its measured no-signal base rate (main's ruling on REVIEW-INPUT3 I-c).
COHORT_KINDS = ("cap-down", "cap-up", "selection-swap")


def cohort_kind(default, arm):
    """cap-down or cap-up for a cohort arm that changes the cap, selection-swap for one that changes the rule."""
    if arm["selection"] != default["selection"] or arm.get("selection_seed") != default.get("selection_seed") \
            or arm.get("stratum") != default.get("stratum"):
        return "selection-swap"
    return "cap-down" if arm["cap"] < default["cap"] else "cap-up"


def base_rate(spec, arm):
    """The manifest's base_rates record of a cohort arm: its kind's measured rate, or not measured."""
    rate = spec.get("base_rates", {}).get(arm["kind"])
    if rate is None:
        return {"arm": arm["name"], "kind": arm["kind"], "rate": None, "ci95": None, "seeds": None,
                "source": "not measured"}
    return {"arm": arm["name"], "kind": arm["kind"], **rate}


def root_key(arm):
    cohort = [arm["cap"], arm["min_class"], arm["selection"], arm.get("selection_seed"), arm.get("stratum")]
    return sha256_bytes(canonical([cohort, arm["k"]]).encode())[:12], cohort


class Runner:
    def __init__(self, args, spec):
        # The binary subprocess runs: a bare name through PATH, as subprocess finds it; that file is hashed.
        grm = shutil.which(args.grm)
        if grm is None:
            raise SystemExit(f"--grm {args.grm}: no such executable (a bare name is looked up on PATH)")
        self.grm, self.raw, self.work, self.spec = str(Path(grm).resolve()), Path(args.raw).resolve(), Path(args.work), spec
        vocabulary = getattr(args, "sigdrug_vocabulary", None)
        self.vocabulary = Path(vocabulary).resolve() if vocabulary else None
        self.curate_inputs = spec.get("curate_inputs", {})
        self.slug = spec["species"]
        self.log = self.work / "logs"
        self.log.mkdir(parents=True, exist_ok=True)
        self.seconds = {}
        self.trace = getattr(args, "trace_training", False)
        self.budget = {**POOL_BUDGET, **spec.get("pool_budget", {})}
        self.run_identical = getattr(args, "run_identical_cohorts", False)
        self.share = not getattr(args, "no_share_ksets", False)
        # The pool's kset folder, once `pool` has it counted; arms skipped by resume.
        self.pool_ksets = None
        self.resumed = []

    def call(self, name, root, *argv):
        started = time.monotonic()
        cmd = [self.grm, "--data", str(root), *argv]
        with open(self.log / f"{name}.log", "ab") as log:
            log.write((" ".join(cmd) + "\n").encode())
            log.flush()
            rc = subprocess.call(cmd, stdout=log, stderr=subprocess.STDOUT, env=child_env())
        self.seconds[name] = self.seconds.get(name, 0.0) + time.monotonic() - started
        if rc != 0:
            raise RuntimeError(f"{' '.join(cmd)} exited {rc}; see {self.log / (name + '.log')}")

    def link_raw(self, root):
        """Links RAW's folders into a data root as revision 12's S0 laid its root out: bvbrc/, fasta/, catalog/,
        and species-check/ when RAW has one, so grm curate drops the genomes the species check fails."""
        for part in RAW_PARTS:
            link = root / part
            if (self.raw / part).exists() and not link.exists():
                link.symlink_to(self.raw / part, target_is_directory=True)

    def check_curate_inputs(self):
        """What the spec's curate_inputs require: RAW's species check and a SigDrug vocabulary, as S0 curated. A spec
        must state them (main's ruling on REVIEW-INPUT8 L3)."""
        if "curate_inputs" not in self.spec:
            raise SystemExit("the spec states no curate_inputs (species_check, sigdrug_vocabulary): say how its roots "
                             "curate")
        if self.curate_inputs.get("species_check") and not (self.raw / "species-check" / "species_check.json").is_file():
            raise SystemExit(f"the spec curates with the species check, and RAW {self.raw} has no "
                             "species-check/species_check.json")
        if self.curate_inputs.get("sigdrug_vocabulary") and self.vocabulary is None:
            raise SystemExit("the spec curates with GRM's SigDrug vocabulary copy: pass --sigdrug-vocabulary DIR (the "
                             "copy revision 12's S0 used)")

    def check_synthetic(self):
        """A spec that says "synthetic": true runs only on a synthetic RAW (main's ruling on REVIEW-INTEG7 L2): RAW's
        truth.json must be synthetic_root.py's (its format) and list exactly the genomes of RAW's bvbrc/genomes.jsonl.
        The format alone could sit in a copied file; the genome list ties the marker to this RAW's own data, which
        a real root's BV-BRC metadata can't match by accident."""
        if self.spec.get("synthetic") is not True:
            return

        def refuse(why):
            raise SystemExit(f'the spec says "synthetic": true, but RAW {self.raw} is not a synthetic root: {why}. A '
                             'real run is never ungated: drop "synthetic" and give the spec its revision12_gate')
        try:
            truth = json.loads((self.raw / "truth.json").read_text())
        except (OSError, ValueError) as e:
            refuse(f"no readable truth.json ({e.__class__.__name__})")
        if not isinstance(truth, dict) or truth.get("format") != SYNTHETIC_FORMAT:
            refuse(f"truth.json's format is not {SYNTHETIC_FORMAT}")
        try:
            listed = {g["genome_id"] for g in truth["genomes"]}
            with open(self.raw / "bvbrc" / "genomes.jsonl", encoding="utf-8") as f:
                held = {json.loads(line)["genome_id"] for line in f if line.strip()}
        except (OSError, ValueError, KeyError, TypeError) as e:
            refuse(f"truth.json's genomes can't be checked against bvbrc/genomes.jsonl ({e.__class__.__name__})")
        if not listed or listed != held:
            refuse(f"truth.json lists {len(listed)} genomes and bvbrc/genomes.jsonl holds {len(held)}, "
                   f"{len(listed ^ held)} of them in only one")

    def check_revision12(self):
        """The H-1 gate: before anything is counted or learned, and in --estimate-only, the default arm must be
        revision 12's cohort in substance (main's ruling after H1DATA, replacing the order gate of REVIEW-INPUT7 to
        REVIEW-INPUT9). A spec without the gate refuses unless it says "synthetic": true. With it, for every species
        `genome_set` names (the run's own species must be one), `substance` of the default root must give:
        - its `genome_set` value: revision 12's genomes, as a sorted set;
        - for each dataset `datasets[<slug>]` pins, that dataset, with the same `datasets` value: the same spec
          outside the eligible-pool counts and the same phenotype row for every selected genome;
        - no dataset of those species that the gate doesn't pin.
        The drawn order (the run's own species from curated.json `drawn_order_sha256`, the others from their
        untouched species/<slug>/genomes.txt) and the datasets.json view (`datasets_view`) are recorded beside the
        spec's `recorded` values, and decide nothing. Returns the record."""
        gate = self.spec.get("revision12_gate")
        if not gate:
            if self.spec.get("synthetic") is True:
                return None
            raise SystemExit('the spec has no revision12_gate and doesn\'t say "synthetic": true; a real run is never '
                             "ungated")
        sets, pinned = gate.get("genome_set") or {}, gate.get("datasets") or {}
        if self.slug not in sets:
            raise SystemExit(f"the revision12_gate names no genome set for {self.slug} (it names "
                             f"{sorted(sets) or 'none'}): add the sha256 of revision 12's species/{self.slug}/"
                             "genomes.txt sorted, from run.substance")
        root = self.curated(self.default)
        for slug in sorted(sets):
            if not (root / "species" / slug / "genomes.txt").is_file():
                raise SystemExit(f"the default root has no species/{slug}/genomes.txt for the gate to compare")
            if not pinned.get(slug):
                raise SystemExit(f"the revision12_gate pins no dataset of {slug}")
        have = substance(root, sorted(sets))
        for slug, want in sorted(sets.items()):
            if have["genome_set"][slug] != want:
                raise SystemExit(f"the default arm's genome set for {slug} is {have['genome_set'][slug]}, not "
                                 f"revision 12's {want}: this RAW or curate isn't revision 12's, so no arm runs")
        for slug in sorted(sets):
            for ds, want in sorted(pinned[slug].items()):
                got = have["datasets"][slug].get(ds)
                if got is None:
                    raise SystemExit(f"the default arm has no dataset {ds}, which revision 12 has, so no arm runs")
                if got != want:
                    raise SystemExit(f"the default arm's dataset {ds} is {got}, not revision 12's {want}: its spec "
                                     "(outside the eligible-pool counts) or a selected genome's phenotype row differs, "
                                     "so no arm runs")
            # An extra dataset would be learned and reported though revision 12 has none (main's ruling after
            # REVIEW-H1DATA, L1).
            extra = sorted(set(have["datasets"][slug]) - set(pinned[slug]))
            if extra:
                raise SystemExit(f"the default arm has {len(extra)} dataset(s) of {slug} that revision 12 hasn't: "
                                 f"{', '.join(extra)}, so no arm runs")
        order = {slug: (json.loads((root / "curated.json").read_text())["drawn_order_sha256"] if slug == self.slug
                        else file_sha256(root / "species" / slug / "genomes.txt")) for slug in sorted(sets)}
        recorded = {"genome_order": order, "datasets_json_sha256": datasets_view(root / "curated" / "datasets.json")}
        spec_recorded = gate.get("recorded") or {}
        recorded["as_the_spec_records"] = all(spec_recorded.get(k) == v for k, v in recorded.items())
        return {"genome_set": {slug: have["genome_set"][slug] for slug in sorted(sets)},
                "datasets": {slug: {ds: have["datasets"][slug][ds] for ds in sorted(pinned[slug])} for slug in sorted(sets)},
                "recorded": recorded}

    def genomes_txt(self, root):
        return root / "species" / self.slug / "genomes.txt"

    def curated(self, arm):
        """The arm's data root, curated, with genomes.txt in canonical order."""
        key, _ = root_key(arm)
        root = self.work / "roots" / f"{key}-k{arm['k']}"
        if (root / "curated.json").exists():
            return root
        root.mkdir(parents=True, exist_ok=True)
        self.link_raw(root)
        curate = ["curate", "--cap", str(arm["cap"]), "--min-class", str(arm["min_class"]),
                  "--selection", arm["selection"]]
        if self.vocabulary:
            curate += ["--sigdrug-vocabulary", str(self.vocabulary)]
        if arm["selection"] != "legacy":
            curate += ["--selection-seed", str(arm["selection_seed"])]
            if arm.get("stratum"):
                curate += ["--stratum", str(self.raw / arm["stratum"])]
        self.call(f"root-{root.name}", root, *curate)
        for name in CURATED_SIDE_FILES:
            (root / "curated" / name).unlink(missing_ok=True)
        drawn = self.canonicalize(root)
        (root / "curated.json").write_text(canonical({"arm": arm["name"], "drawn_order_sha256": drawn}))
        return root

    def canonicalize(self, root):
        """Rewrites genomes.txt sorted by genome id; returns the drawn order's sha256."""
        path = self.genomes_txt(root)
        text = path.read_text(encoding="utf-8")
        ids = [line.strip() for line in text.splitlines() if line.strip()]
        path.write_text("\n".join(sorted(ids)) + "\n", encoding="utf-8")
        return sha256_bytes(text.encode())

    def build(self, root, k):
        """Count, structure and matrix of a canonicalized root."""
        name = f"root-{root.name}"
        self.call(name, root, "count", "--species", self.slug, "--k", str(k))
        self.call(name, root, "structure", "--species", self.slug)
        self.call(name, root, "matrix", "--species", self.slug)

    def union(self, arms):
        """The species pool's genomes: the union of every arm's curated cohort (main's ruling on REVIEW-INPUT4
        M1), sorted by id. Curates every arm's root first."""
        return sorted(set().union(*(self.genome_set(arm) for arm in arms)))

    def kset_estimate(self, n):
        """(bytes, bytes per genome, where that came from) for n genomes' ksets: the mean size of the species'
        own counted ksets in RAW when there are any, else the spec's stated default."""
        sizes = [p.stat().st_size for p in (self.raw / "species" / self.slug / "ksets").glob("*.kset")]
        if sizes:
            per, source = sum(sizes) / len(sizes), f"mean of the species' {len(sizes)} counted ksets in RAW"
        else:
            per = self.budget["default_kset_bytes"]
            source = f"the spec's stated default ({per / 2 ** 20:.1f} MiB per genome)"
        return n * per, per, source

    def kset_dir(self, root):
        return root / "species" / self.slug / "ksets"

    def finished_ksets(self, root):
        """The genomes whose kset `grm count` finished in this root: the kset and its sidecar both written."""
        ksets = self.kset_dir(root)
        if not ksets.is_dir():
            return set()
        return {p.stem for p in ksets.glob("*.kset") if (ksets / f"{p.stem}.json").exists()}

    def shares(self, arm):
        """Whether the arm's root takes its ksets from the pool (G1): it is counted at the pool's k, the default's."""
        return self.share and arm["k"] == self.default["k"]

    def pool_kset_ids(self, arms):
        """The genomes whose pool kset a root can link: every genome of the union while the pool is still to be
        counted (the pool's pending ksets count it), else the pool's finished ksets."""
        union = self.union(arms)
        root = self.pool_root(union)
        if (root / "species" / self.slug / "structure.json").exists():
            return self.finished_ksets(root)
        return set(union)

    def counted_roots(self, arms):
        """Every distinct data root the run counts, by name, with its genome count, k and the ksets it still has
        to write (`pending`; 0 once built.json exists): the default's and each arm's that has its own (an
        identical cohort is not built; blocks, budget and calibration arms reuse the default's root). A root that
        shares the pool's ksets (G1) writes none of those (`shared`). Curates them; counts nothing."""
        out = {}
        for arm in arms:
            if (arm["factor"] == "cohort" and not self.run_identical
                    and self.genome_set(arm) == self.genome_set(self.default)):
                continue
            key, _ = root_key(arm)
            name = f"{key}-k{arm['k']}"
            if name not in out:
                root, ids = self.work / "roots" / name, self.genome_set(arm)
                built = (root / "built.json").exists()
                pending = set() if built else ids - self.finished_ksets(root)
                shared = pending & self.pool_kset_ids(arms) if self.shares(arm) else set()
                out[name] = {"genomes": len(ids), "k": arm["k"], "built": built, "arms": [],
                             "shared": len(shared), "pending": len(pending - shared)}
            out[name]["arms"].append(arm["name"])
        return out

    def share_ksets(self, root):
        """Gives every genome of a root counted at the pool's k the pool's finished kset as a hard link and its
        sidecar as a copy, in the root's kset folder, where `grm count` reuses them (G1). `grm count` writes a kset
        only to a .part file it renames, which replaces the link, but it rewrites a reused sidecar in place (when its
        catalog version differs or it has no kset digest), so a linked sidecar would carry that write into the pool
        (REVIEW-H1CODE L1): the sidecar, a few kilobytes, is copied. A genome the root already has a kset or sidecar
        for keeps its own. Returns the genomes shared."""
        mine, pool = self.kset_dir(root), self.pool_ksets
        mine.mkdir(parents=True, exist_ok=True)
        have = self.finished_ksets(pool.parents[2])
        linked = []
        for g in sorted(set(self.genomes_txt(root).read_text(encoding="utf-8").split()) & have):
            if any(os.path.lexists(mine / f"{g}.{ext}") for ext in ("kset", "json")):
                continue
            try:
                os.link(pool / f"{g}.kset", mine / f"{g}.kset")
                shutil.copyfile(pool / f"{g}.json", mine / f"{g}.json")
            except OSError as e:
                raise SystemExit(f"can't share the pool's kset of {g} with root {root.name} ({e}); the disk check "
                                 "counted it once") from None
            linked.append(g)
        return linked

    def record_shared(self, root, linked):
        """After `grm count`: every shared kset must still be the pool's file, and its sidecar copy must still hold
        the pool's bytes. A recount writes a new kset and sidecar, and a sidecar rewritten for another catalog means
        the root counted otherwise than the pool, so either refuses. Writes the root's shared.json: its genomes, how
        many ksets are the pool's, and how many it counted itself."""
        mine, pool = self.kset_dir(root), self.pool_ksets

        def same_kset(g):
            a, b = mine / f"{g}.kset", pool / f"{g}.kset"
            return a.exists() and b.exists() and os.path.samefile(a, b)

        def same_sidecar(g):
            a, b = mine / f"{g}.json", pool / f"{g}.json"
            return a.is_file() and b.is_file() and a.read_bytes() == b.read_bytes()

        rewritten = [g for g in linked if not (same_kset(g) and same_sidecar(g))]
        if rewritten:
            raise SystemExit(f"root {root.name}: grm count rewrote the kset or sidecar of {len(rewritten)} of the "
                             f"{len(linked)} genomes shared with the species pool (e.g. {rewritten[0]}), which the "
                             "disk check counted once")
        ids = self.genomes_txt(root).read_text(encoding="utf-8").split()
        shared = sum(same_kset(g) for g in ids)
        (root / SHARED).write_text(canonical({"genomes": len(ids), "shared": shared, "own": len(ids) - shared,
                                              "linked_now": len(linked), "pool": pool.parents[2].name}))

    def missing_fastas(self, union):
        """The union's genomes with no FASTA in RAW, which the pool refuses."""
        return [g for g in union if not (self.raw / "fasta" / f"{g}.fna.gz").exists()]

    def pool_root(self, union):
        return self.work / "roots" / f"pool-{sha256_bytes(canonical(union).encode())[:12]}"

    def pool_pending(self, union):
        """The pool ksets still to write: none once its structure.json exists."""
        root = self.pool_root(union)
        if (root / "species" / self.slug / "structure.json").exists():
            return 0
        return len(set(union) - self.finished_ksets(root))

    def disk_check(self, union, free=None, roots=None, measure=True, pool_pending=None, pool_curated=0):
        """The disk the run still has to write, against max_disk_share of the free disk. On a fresh WORK that is
        the whole comparison; on a reused one, only what isn't finished:
        - the pool's pending ksets (`pool_pending`, default all of the union);
        - per data root (`counted_roots`), its pending ksets and, unless built.json exists, its matrix and
          sidecars at root_overhead_share of all its genomes' ksets;
        - the arm outputs (WORK/arms), rewritten on every run, at arm_output_share of every root's ksets;
        - the pool's copy of the default's curated/ (`pool_curated`, its size in disk blocks, until the pool has it);
        - WORK-level files (logs, inputs.json, disk-estimate.json, species-groups.json, j/, k/, the report) at the
          stated work_files_bytes.
        `grm count` reuses a kset only from its own root's kset folder, so a root counted at the pool's k has the
        pool's ksets linked in (G1) and its pending ksets leave those out; a root at another k counts its own. A
        k21 or k25 root is sized as a k31 one, an upper bound. `free` is the box's free disk (--free-disk-bytes);
        `measure` also reads WORK's volume now and takes the smaller. A dict with ok False when it doesn't fit."""
        roots = roots or {}
        _, per, source = self.kset_estimate(1)
        pool = (len(union) if pool_pending is None else pool_pending) * per
        sized = {}
        for name, r in sorted(roots.items()):
            pending = r.get("pending", r["genomes"])
            matrix = 0 if r.get("built") else r["genomes"] * per * self.budget["root_overhead_share"]
            sized[name] = {**r, "bytes": round(pending * per + matrix)}
        arm_outputs = sum(r["genomes"] for r in roots.values()) * per * self.budget["arm_output_share"]
        work_files = self.budget["work_files_bytes"]
        need = pool + sum(r["bytes"] for r in sized.values()) + arm_outputs + pool_curated + work_files
        measured = None
        if measure:
            self.work.mkdir(parents=True, exist_ok=True)
            measured = shutil.disk_usage(self.work).free
        free = min(f for f in (free, measured) if f is not None)
        share = self.budget["max_disk_share"]
        return {"genomes": len(union), "pool_bytes": round(pool), "roots": sized, "arm_output_bytes": round(arm_outputs),
                "pool_curated_bytes": round(pool_curated), "work_files_bytes": round(work_files),
                "kset_bytes": round(need), "per_genome_bytes": round(per), "source": source, "free_bytes": free,
                "max_disk_share": share, "ok": need <= share * free}

    def plan_disk(self, arms, free=None, measure=True):
        """`disk_check` over this run: the arms' union, the pool's pending ksets and every data root."""
        union = self.union(arms)
        return self.disk_check(union, free, self.counted_roots(arms), measure, self.pool_pending(union),
                               self.pool_curated(union))

    def pool_curated(self, union):
        """The pool's copy of the default's curated/ still to make, in disk blocks; 0 once the pool has it."""
        if (self.pool_root(union) / "curated").exists():
            return 0
        return disk_bytes(self.curated(self.default) / "curated", CURATED_SIDE_FILES)

    def check_disk(self, arms, free=None):
        """Before anything is counted, on every run: refuses when what is still to be written passes the disk
        share; writes WORK/disk-estimate.json."""
        disk = self.plan_disk(arms, free)
        (self.work / "disk-estimate.json").write_text(json.dumps(disk, indent=1, sort_keys=True))
        if not disk["ok"]:
            raise SystemExit(f"the comparison still has to write about {disk['kset_bytes'] / 2 ** 30:.1f} GiB (the "
                             f"pool's pending ksets of {disk['genomes']} genomes, {len(disk['roots'])} data roots and "
                             f"the arm outputs, {disk['source']}), over {disk['max_disk_share']:.0%} of "
                             f"{disk['free_bytes'] / 2 ** 30:.1f} GiB free")
        return disk

    def check_inputs(self):
        """Ties WORK to its inputs: the sha256 of RAW's input manifest (every file under bvbrc/, fasta/, catalog/
        and species-check/ but RAW_SIDE_FILES, and each arm's stratum file, with its sha256 and size; RAW is walked
        without following
        links and a symlink in it is refused), of the grm binary subprocess runs, and of the SigDrug vocabulary.
        A fresh WORK records them in inputs.json; a reused WORK whose inputs changed, or that has roots but no
        record, is refused. RAW's manifest lines go to raw-manifest.txt, for `combine`."""
        def manifest(top, parts, label):
            lines = []
            for part in parts:
                base = top / part
                if not base.exists():
                    continue
                if base.is_symlink():
                    raise SystemExit(f"{label} {base} is a symlink; the input manifest follows no links")
                if base.is_file():
                    lines.append(f"{file_sha256(base)}  {base.stat().st_size}  {part}")
                    continue
                for here, dirs, files in os.walk(base, followlinks=False):
                    for name in sorted(dirs + files):
                        p = Path(here) / name
                        if p.is_symlink():
                            raise SystemExit(f"{label} holds a symlink, {p}; the input manifest follows no links")
                    dirs.sort()
                    for name in sorted(files):
                        p = Path(here) / name
                        if top == self.raw and p.relative_to(top).as_posix() in RAW_SIDE_FILES:
                            continue
                        lines.append(f"{file_sha256(p)}  {p.stat().st_size}  {p.relative_to(top)}")
            return lines

        # No symlink anywhere in RAW, not only under the hashed parts (REVIEW-INPUT8 L5).
        for here, dirs, files in os.walk(self.raw, followlinks=False):
            for name in sorted(dirs + files):
                if (Path(here) / name).is_symlink():
                    raise SystemExit(f"RAW holds a symlink, {Path(here) / name}; the input manifest follows no links")
        strata = sorted({a["stratum"] for a in self.all_arms if a.get("stratum")})
        for s in strata:
            if not (self.raw / s).resolve().is_relative_to(self.raw):
                raise SystemExit(f"stratum {s} is outside RAW {self.raw}")
        lines = manifest(self.raw, [*RAW_PARTS, *strata], "RAW")
        now = {"raw_manifest_sha256": sha256_bytes("\n".join(lines).encode()), "raw_files": len(lines),
               "grm_sha256": file_sha256(Path(self.grm)), "grm": self.grm,
               "vocabulary_sha256": (sha256_bytes("\n".join(manifest(self.vocabulary, ["."], "the vocabulary"))
                                                  .encode()) if self.vocabulary else None)}
        record = self.work / INPUTS
        if record.exists():
            was = json.loads(record.read_text())
            changed = [k for k in ("raw_manifest_sha256", "grm_sha256", "vocabulary_sha256")
                       if was.get(k) != now[k]]
            if changed:
                raise SystemExit(f"WORK {self.work} was made from another {' and '.join(changed)} "
                                 f"({', '.join(f'{k} {was.get(k)} now {now[k]}' for k in changed)}); use a fresh WORK")
        elif (self.work / "roots").exists():
            raise SystemExit(f"WORK {self.work} has data roots but no {INPUTS}, so its RAW and grm are unknown; "
                             "use a fresh WORK")
        else:
            self.work.mkdir(parents=True, exist_ok=True)
            record.write_text(canonical(now))
        (self.work / RAW_MANIFEST).write_text("\n".join(lines), encoding="utf-8")
        return now

    def pool(self, arms):
        """The species grouping every arm splits on: one root over the union of every arm's curated cohort,
        counted at the default's k and structured once, before any grm arm (after `check_disk`). Refuses a
        genome of the union without a FASTA (with its id), and a pool whose structure.json doesn't hold exactly
        the union. Returns its structure.json (grm arm --groups)."""
        union = self.union(arms)
        root = self.pool_root(union)
        structure = root / "species" / self.slug / "structure.json"
        if not structure.exists():
            missing = self.missing_fastas(union)
            if missing:
                raise SystemExit(f"the species pool needs a FASTA for every genome of the arms' union; "
                                 f"{len(missing)} of {len(union)} have none: {', '.join(missing[:20])}")
            (root / "species" / self.slug).mkdir(parents=True, exist_ok=True)
            self.link_raw(root)
            if not (root / "curated").exists():
                shutil.copytree(self.curated(self.default) / "curated", root / "curated",
                                ignore=shutil.ignore_patterns(*CURATED_SIDE_FILES))
            self.genomes_txt(root).write_text("\n".join(union) + "\n", encoding="utf-8")
            (root / "curated.json").write_text(canonical({"arm": "species-pool", "genomes": len(union)}))
            name = f"root-{root.name}"
            self.call(name, root, "count", "--species", self.slug, "--k", str(self.default["k"]))
            self.call(name, root, "structure", "--species", self.slug)
        s = json.loads(structure.read_text())
        if s["genomes"] != union:
            lost = sorted(set(union) - set(s["genomes"]))
            raise SystemExit(f"the species pool holds {len(s['genomes'])} genomes, not the arms' union of "
                             f"{len(union)}; missing {', '.join(lost[:20])} (see {self.log}/{root.name}.log)")
        if s["k"] != self.default["k"]:
            raise SystemExit(f"the species pool was counted at k = {s['k']}, the default at k = {self.default['k']}")
        self.pool_ksets = self.kset_dir(root)
        return structure

    def root(self, arm):
        """The arm's data root, built; a root counted at the pool's k first links the pool's ksets (G1), once the
        pool is counted."""
        root = self.curated(arm)
        done = root / "built.json"
        if not done.exists():
            linked = self.share_ksets(root) if self.shares(arm) and self.pool_ksets else []
            self.build(root, arm["k"])
            if self.pool_ksets:
                self.record_shared(root, linked)
            done.write_text(canonical({"arm": arm["name"], "k": arm["k"], "cohort": root_key(arm)[1]}))
        return root

    def genome_set(self, arm):
        return set(self.genomes_txt(self.curated(arm)).read_text(encoding="utf-8").split())

    def eligible(self, arm):
        """The species' eligible genomes: an unconflicted reported R or S call on one of its dataset drugs."""
        root = self.curated(arm)
        drugs = {d["antibiotic"] for d in json.loads((root / "curated" / "datasets.json").read_text())
                 if d["species_slug"] == self.slug}
        out = set()
        with open(root / "curated" / "phenotypes.jsonl", encoding="utf-8") as f:
            for line in f:
                p = json.loads(line)
                if (p["antibiotic"] in drugs and p.get("reported") in ("R", "S") and not p.get("conflict")
                        and not p.get("label_from_mic")):
                    out.add(p["genome_id"])
        return out

    def identical(self, arm, dataset, out):
        """Records a cohort arm with the default's genome set as identical_cohort, without running it."""
        n = len(self.genome_set(arm))
        reason = (f"identical_cohort: arm {arm['name']} selects the default's genome set ({n} genomes); "
                  "with canonical genome order it would give the default's predictions, so it is not run")
        runs = [{"dataset": dataset, "arm": arm["name"], "blocks": arm["blocks"], "learner": learner,
                 "split": split, "repeats": self.spec.get("repeats", 1), "fold_plan_sha256": None,
                 "status": "identical_cohort", "reason": reason, "metrics": {}, "fit_count": 0, "seconds": 0.0,
                 "predictions": None, "pairing": arm["pairing"], "arm_plan_sha256": None, "cohort_drops": None}
                for split in SPLITS for learner in LEARNERS]
        (out / "runs.json").write_text(json.dumps(runs, indent=1))
        (out / "predictions.jsonl").write_text("")
        (out / "arm.json").write_text(json.dumps({"arm": arm["name"], "dataset": dataset, "blocks": [],
                                                  "unavailable": [], "identical_cohort": reason}, indent=1))

    def arm(self, arm, dataset, run_identical=False):
        out = self.work / "arms" / arm["name"] / dataset
        arms_dir = (self.work / "arms").resolve()
        if not out.resolve().is_relative_to(arms_dir) or out.resolve() == arms_dir:
            raise SystemExit(f"arm {arm['name']!r}, dataset {dataset!r}: {out} is outside {arms_dir}")
        # Plans from an earlier run of this arm never reach this one (REVIEW-INPUT6 L3); a plans folder that is a
        # link, or that can't be removed, stops the run (REVIEW-INPUT7 L2).
        plans = out / "plans"
        if plans.is_symlink():
            raise SystemExit(f"{plans} is a symlink; refusing to clear it")
        record = self.resume_record(arm, dataset, run_identical)
        if self.complete(out, record):
            self.resumed.append(f"{arm['name']}/{dataset}")
            return out
        (out / RESUME).unlink(missing_ok=True)
        if plans.exists():
            shutil.rmtree(plans)
        out.mkdir(parents=True, exist_ok=True)
        if (arm["factor"] == "cohort" and not run_identical
                and self.genome_set(arm) == self.genome_set(self.default)):
            self.identical(arm, dataset, out)
            self.finished(out, record)
            return out
        root = self.root(arm)
        argv = ["arm", "--species", self.slug, "--dataset", dataset, "--arm", arm["name"],
                "--blocks", ",".join(arm["blocks"]), "--screened-features", str(arm["screened_features"]),
                "--repeats", str(self.spec.get("repeats", 1)), "--pairing", arm["pairing"], "--out", str(out)]
        if arm["name"] != DEFAULT_ARM:
            argv += ["--reference", str(self.work / "arms" / DEFAULT_ARM / dataset)]
        for key in ("balance_seed", "learner_seed"):
            if key in arm:
                argv += [f"--{key.replace('_', '-')}", str(arm[key])]
        if self.trace:
            argv.append("--trace-training")
        if arm["name"] == DEFAULT_ARM or arm["pairing"] == "common_outer_test_set_balanced":
            argv += ["--groups", str(self.groups)]
        self.call(f"arm-{arm['name']}-{dataset}", root, *argv)
        self.finished(out, record)
        return out

    def resume_record(self, arm, dataset, run_identical):
        """What an arm's output on a dataset was made from (G6): WORK's inputs.json (RAW's manifest, the grm binary,
        the vocabulary), the spec, the species grouping, the reference plans (the default arm's plans for this
        dataset; none for the default itself) and the flags that change `grm arm`'s output."""
        def sha(path):
            return file_sha256(path) if path and Path(path).is_file() else None

        ref = self.work / "arms" / DEFAULT_ARM / dataset / "plans"
        return {"arm": arm["name"], "dataset": dataset, "inputs_sha256": sha(self.work / INPUTS),
                "spec_sha256": sha256_bytes(canonical(self.spec).encode()),
                "groups_sha256": sha(getattr(self, "groups", None)),
                "reference_plans": {} if arm["name"] == DEFAULT_ARM
                else {p.name: file_sha256(p) for p in sorted(ref.glob("*.json"))},
                "trace_training": bool(self.trace), "run_identical_cohorts": bool(run_identical)}

    @staticmethod
    def output_files(out):
        return {str(p.relative_to(out)) for p in out.rglob("*") if p.is_file() and p.relative_to(out) != Path(RESUME)}

    def finished(self, out, record):
        """Writes resume.json after an arm finished: its record and the sha256 of every file it wrote."""
        outputs = {name: file_sha256(out / name) for name in sorted(self.output_files(out))}
        (out / RESUME).write_text(canonical({**record, "outputs": outputs}))

    def complete(self, out, record):
        """Whether a reused WORK may skip this arm on this dataset: its resume.json holds the same record, its
        outputs include runs.json, predictions.jsonl and arm.json, the folder holds exactly the files it lists, none
        a symlink, and each still has its sha256."""
        try:
            was = json.loads((out / RESUME).read_text())
            outputs = was.pop("outputs")
        except (OSError, ValueError, KeyError, AttributeError):
            return False
        if was != record or not isinstance(outputs, dict) or not ARM_OUTPUTS <= set(outputs):
            return False
        if any(p.is_symlink() for p in out.rglob("*")) or self.output_files(out) != set(outputs):
            return False
        return all(file_sha256(out / name) == h for name, h in outputs.items())


def read_rows(path):
    with open(path, encoding="utf-8") as f:
        return [json.loads(line) for line in f if line.strip()]


def contrasts(spec, arms, work):
    """J: one family C log_loss contrast per (dataset, arm, learner, split) against the default."""
    margins = {"log_loss": spec["margin"]["log_loss"]}
    development = [d["id"] for d in spec["datasets"] if d["partition"] == "development"]
    records, planned, missing = [], [], {}
    for d in spec["datasets"]:
        ds = d["id"]
        runs = {}
        rows = []
        for arm in [DEFAULT_ARM] + [a["name"] for a in arms]:
            out = work / "arms" / arm / ds
            for u in json.loads((out / "runs.json").read_text()):
                runs[(arm, u["learner"], u["split"])] = u
            rows += read_rows(out / "predictions.jsonl")
        for arm in arms:
            for learner in spec.get("learners", LEARNERS):
                for split in spec.get("splits", SPLITS):
                    decides = d["partition"] == "confirmatory" and split not in k_schema.DESCRIPTIVE_SPLITS
                    key = stats.contrast_key("C" if decides else "D", learner, split,
                                             DEFAULT_ARM, arm["name"], "log_loss")
                    pid = f"{ds}/{key}"
                    planned.append(pid)
                    down = [runs[(a, learner, split)] for a in (DEFAULT_ARM, arm["name"])
                            if runs[(a, learner, split)]["status"] != "ok"]
                    if down:
                        missing[pid] = f"arm {down[0]['arm']} has no {learner} {split} result: {down[0]['reason']}"
                        continue
                    try:
                        paired = stats.pair_arms(rows, dataset=ds, learner=learner, split=split,
                                                 baseline=DEFAULT_ARM, treatment=arm["name"],
                                                 pairing=arm["pairing"])
                    except stats.PairingError as e:
                        missing[pid] = f"pairing refused: {e}"
                        continue
                    records += stats.compare(paired, family=key.split(":")[0], endpoints=["log_loss"],
                                             margins=margins, development=development,
                                             seed=spec.get("seed", 0))
    return stats.finalize(records, planned, missing), planned, missing


def manifest(spec, default, arms, work, binary_sha256, source_commit, seconds):
    """K's run manifest: every arm's runs, with each paired arm's rule and reference."""
    blocks, runs, datasets = {}, [], []
    for d in spec["datasets"]:
        ds = d["id"]
        record = json.loads((work / "arms" / DEFAULT_ARM / ds / "arm.json").read_text())
        lineages = len({c for c in json.loads((work / "arms" / DEFAULT_ARM / ds / "plans" / "lineage-0.json")
                                              .read_text())["cluster"]})
        datasets.append({"id": ds, "species": record["species"], "partition": d["partition"],
                         "n_samples": record["samples"], "group_counts": {"lineage": lineages}})
        for arm in [default] + arms:
            out = work / "arms" / arm["name"] / ds
            a = json.loads((out / "arm.json").read_text())
            # One spec per block id: the default root's. A k or cohort arm's own K spec (another k
            # or cohort) is in its arm.json.
            for b in a["blocks"]:
                if arm["factor"] in (None, "blocks", "budget"):
                    blocks.setdefault(b["id"], b)
            for bid, reason in a["unavailable"]:
                blocks.setdefault(bid, {"id": bid, "kind": {"C": "CuratedExact", "H": "Hierarchical",
                                                            "P": "Presence", "N": "Count", "X": "Context",
                                                            "K": "Anonymous"}[bid],
                                        "version": "", "source_sha256": "",
                                        "priority": "CHPNXK".index(bid), "screen_cap": None, "catalog_date": None,
                                        "circularity_risk": bid != "K", "unavailable_reason": reason})
            for u in json.loads((out / "runs.json").read_text()):
                run = {k: u[k] for k in ("dataset", "arm", "blocks", "learner", "split", "repeats",
                                         "fold_plan_sha256", "status", "reason", "metrics", "fit_count", "seconds")}
                run["id"] = f"{u['arm']}.{ds}.{u['learner']}.{u['split']}"
                run["predictions"] = (f"arms/{arm['name']}/{ds}/predictions.jsonl" if u["status"] == "ok"
                                      else None)
                if arm["name"] != DEFAULT_ARM:
                    run["pairing"] = {"rule": u["pairing"], "reference_arm": DEFAULT_ARM,
                                      "arm_plan_sha256": u["arm_plan_sha256"]}
                runs.append(run)
    return {
        "schema_version": k_schema.SCHEMA_VERSION,
        "protocol_sha256": sha256_bytes(canonical(spec).encode()),
        "source_commit": source_commit,
        "binary_sha256": binary_sha256,
        "blocks": [blocks[b] for b in sorted(blocks, key="CHPNXK".index)],
        "datasets": datasets,
        "missingness": [],
        "group_proxy_audit": [],
        "runs": runs,
        "compute": {"planned_fit_units": sum(u["fit_count"] for u in runs),
                    "actual_seconds": round(sum(seconds.values()), 3), "tiers_completed": ["full"]},
        "sentinel_arms": [a["name"] for a in arms if a.get("sentinel")],
        "base_rates": [base_rate(spec, a) for a in arms if a["factor"] == "cohort"],
    }


def refused_plans(work, datasets, arms):
    """Every refused (unavailable) fold or cohort plan of this run's `arms`, with its reason (REVIEW-INPUT3 I-a);
    a folder a reused WORK keeps for an arm no longer in the spec is not read."""
    out = []
    for ds in datasets:
        for plans in sorted(work / "arms" / arm / ds / "plans" for arm in arms):
            if not plans.is_dir():
                continue
            for path in sorted(plans.glob("*.json")):
                p = json.loads(path.read_text())
                reason = p["plan"]["unavailable_reason"] if "plan" in p else p["unavailable_reason"]
                if reason is not None:
                    out.append({"dataset": ds, "arm": plans.parent.parent.name, "plan": path.name, "reason": reason})
    return out


def summarize(spec, arms, contrasts, missing, unreliable=()):
    """Each arm's decision per dataset, learner and split. A random-split contrast reads "descriptive"
    and a contrast K marked unreliable reads "unreliable"; raw_decisions keeps J's value."""
    flagged = {i for u in unreliable for i in u["contrasts"]}
    out = {}
    for arm in arms:
        rows, raw = {}, {}
        for x in contrasts:
            if x["treatment"] == arm["name"]:
                key = f"{x['dataset']}:{x['learner']}:{x['split']}"
                raw[key] = x["decision"]
                rows[key] = ("descriptive" if x["family"] == "D" and x["split"] in k_schema.DESCRIPTIVE_SPLITS
                             else "unreliable" if x["id"] in flagged else x["decision"])
        for pid, reason in missing.items():
            parts = k_schema.parse_key(k_schema.contrast_key(pid))
            if parts["treatment"] == arm["name"]:
                rows[f"{k_schema.contrast_dataset(pid)}:{parts['learner']}:{parts['split']}"] = \
                    f"not produced: {reason}"
        decisions = {v for v in rows.values()}
        out[arm["name"]] = {"factor": arm["factor"], "pairing": arm["pairing"], "decisions": rows,
                            "raw_decisions": raw, "sentinel": bool(arm.get("sentinel")),
                            "any_improves": "improves" in decisions, "any_worse": "worse" in decisions}
    return out


def main(argv):
    if argv[:1] == ["combine"]:
        return combine(argv[1:])
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("spec")
    p.add_argument("--raw", required=True)
    p.add_argument("--work", required=True)
    p.add_argument("--grm", required=True)
    p.add_argument("--source-commit", required=True)
    p.add_argument("--jobs", type=int, default=1)
    p.add_argument("--run-identical-cohorts", action="store_true",
                   help="run cohort arms with the default's genome set anyway (the order-twin proof only)")
    p.add_argument("--trace-training", action="store_true",
                   help="have every grm arm write training.jsonl (each fold's training genome ids)")
    p.add_argument("--below-panel-floor", action="store_true",
                   help="allow a species below the spec's panel floor (synthetic studies only)")
    p.add_argument("--estimate-only", action="store_true",
                   help="curate every arm, print the species pool's kset estimate against --free-disk-bytes "
                        "and exit (1 when it doesn't fit): the check to make before any box")
    p.add_argument("--sigdrug-vocabulary",
                   help="GRM's SigDrug vocabulary copy for grm curate (the one revision 12's S0 used)")
    p.add_argument("--free-disk-bytes", type=int,
                   help="the box's free disk; with --estimate-only it replaces this machine's, otherwise the "
                        "smaller of it and WORK's volume is used")
    p.add_argument("--no-share-ksets", action="store_true",
                   help="count every root's own ksets instead of linking the pool's (the G1 identity proof only)")
    a = p.parse_args(argv)
    spec, default, arms = load_spec(a.spec)
    work = Path(a.work)
    runner = Runner(a, spec)
    runner.default = default
    runner.all_arms = [default, *arms]
    runner.check_curate_inputs()
    runner.check_synthetic()
    runner.check_inputs()
    floor = spec.get("panel_floor", {}).get("min_eligible_genomes", 0)
    eligible = len(runner.eligible(default))
    if eligible < floor and not a.below_panel_floor:
        raise SystemExit(f"{spec['species']}: {eligible} eligible genomes, below the panel floor of {floor} "
                         "(spec panel_floor): cap-4000 and grr-v1 would not change the genome set")
    gate = runner.check_revision12()
    print(gate_line(gate, spec), file=sys.stderr)
    # The gate's record, drawn order and datasets.json view included (K's manifest, a closed record, keeps "passed").
    (work / "revision12-gate.json").write_text(json.dumps(gate, indent=1, sort_keys=True))
    if a.estimate_only:
        # Before any box: the gate, then the box's free disk (--free-disk-bytes, not this machine's), and the union's
        # FASTAs, which the pool refuses without (it reports them; the exit code is the disk's).
        disk = runner.plan_disk([default, *arms], a.free_disk_bytes, measure=a.free_disk_bytes is None)
        union = runner.union([default, *arms])
        missing = runner.missing_fastas(union)
        disk["union_without_fasta"], disk["union_without_fasta_examples"] = len(missing), missing[:20]
        disk["revision12_gate"] = gate
        print(json.dumps(disk, indent=1, sort_keys=True))
        if missing:
            print(f"{len(missing)} of the union's {len(union)} genomes have no FASTA in RAW: the run would refuse at "
                  f"the species pool (e.g. {', '.join(missing[:5])})", file=sys.stderr)
        raise SystemExit(0 if disk["ok"] else 1)
    names = {DEFAULT_ARM, *(arm["name"] for arm in arms)}
    stale = sorted(p.name for p in (work / "arms").glob("*") if p.is_dir() and p.name not in names)
    if stale:
        print(f"WORK keeps folders of arms not in this spec, not read: {', '.join(stale)}", file=sys.stderr)
    started = time.monotonic()
    (work / SPEC_COPY).write_text(canonical(spec))
    runner.check_disk([default, *arms], a.free_disk_bytes)
    runner.groups = runner.pool([default, *arms])
    grouping = json.loads(runner.groups.read_text())
    for arm in [default, *arms]:
        outside = runner.genome_set(arm) - set(grouping["genomes"])
        if outside:
            raise SystemExit(f"arm {arm['name']} holds {len(outside)} genomes outside the species pool, "
                             f"e.g. {sorted(outside)[0]}")
    species_groups = {"structure": str(runner.groups), "sha256": sha256_bytes(runner.groups.read_bytes()),
                      "genomes": len(grouping["genomes"]), "lineages": grouping["lineages"]["n_clusters"],
                      "arms": [arm["name"] for arm in [default, *arms]]}
    (work / "species-groups.json").write_text(canonical(species_groups))
    datasets = [d["id"] for d in spec["datasets"]]
    for ds in datasets:
        runner.arm(default, ds)
    # Roots are built one at a time (count uses every core); arms then run in parallel.
    for arm in arms:
        if arm["factor"] != "cohort" or a.run_identical_cohorts \
                or runner.genome_set(arm) != runner.genome_set(default):
            runner.root(arm)
    with concurrent.futures.ThreadPoolExecutor(max_workers=a.jobs) as pool:
        for f in [pool.submit(runner.arm, arm, ds, a.run_identical_cohorts) for arm in arms for ds in datasets]:
            f.result()
    if runner.resumed:
        print(f"resumed {len(runner.resumed)} arm outputs from this WORK: {', '.join(sorted(runner.resumed))}",
              file=sys.stderr)
    meta = {"synthetic": spec.get("synthetic") is True, "revision12_gate": "not_required" if gate is None else "passed",
            "species_pool": {"genomes": species_groups["genomes"], "arms": species_groups["arms"],
                             "structure_sha256": species_groups["sha256"]}}
    summary = report_run(spec, default, arms, work, file_sha256(Path(runner.grm)), a.source_commit, runner.seconds,
                         meta)
    summary.update({"seconds": round(time.monotonic() - started, 1), "stage_seconds": runner.seconds,
                    "eligible_genomes": eligible, "panel_floor": floor, "species_groups": species_groups,
                    "stale_arm_folders": stale, "resumed_arms": sorted(runner.resumed), "revision12_gate": gate})
    write_summary(work, summary)


def report_run(spec, default, arms, work, binary_sha256, source_commit, seconds, meta):
    """J and K over WORK's arm outputs (WORK/j, WORK/k/manifest.json with `meta`'s synthetic, revision12_gate and
    species_pool, WORK/report), as a run and `combine` both write them; returns summary.json's decision fields."""
    records, planned, missing = contrasts(spec, arms, work)
    (work / "j").mkdir(exist_ok=True)
    for name, value in (("contrasts", records), ("planned", planned), ("missing", missing)):
        (work / "j" / f"{name}.json").write_text(json.dumps(value, indent=1, sort_keys=True))
    m = manifest(spec, default, arms, work, binary_sha256, source_commit, seconds)
    m.update(meta)
    (work / "k").mkdir(exist_ok=True)
    (work / "k" / "manifest.json").write_text(json.dumps(m, indent=1, sort_keys=True))
    report = k_report.build_report(m, records, [], [], planned, [], not_produced=missing,
                                   out_dir=str(work / "report"), prediction_root=str(work))
    datasets = [d["id"] for d in spec["datasets"]]
    return {"contrast_counts": report["contrast_counts"], "unreliable": report["unreliable"],
            "identical_cohorts": sorted({u["arm"] for u in report["identical_cohorts"]}),
            "refused_plans": refused_plans(work, datasets, [DEFAULT_ARM, *(a["name"] for a in arms)]),
            "arms": summarize(spec, arms, records, missing, report["unreliable"])}


def write_summary(work, summary):
    (work / "summary.json").write_text(json.dumps(summary, indent=1, sort_keys=True))
    for name, s in summary["arms"].items():
        counts = {}
        for v in s["decisions"].values():
            counts[v.split(":")[0]] = counts.get(v.split(":")[0], 0) + 1
        print(f"{name:<14} {s['factor']:<7} {s['pairing']:<22} {counts}")
    print(f"refused plans: {len(summary['refused_plans'])}; report: {work / 'report'}; {summary['seconds']} s")


def read_part(work):
    """A finished run's records that `combine` checks: its spec.json, inputs.json, raw-manifest.txt (which must be
    the manifest inputs.json hashed), K's manifest, species-groups.json and summary.json."""
    def load(name):
        try:
            return json.loads((work / name).read_text(encoding="utf-8"))
        except (OSError, ValueError) as e:
            raise SystemExit(f"part {work}: no readable {name} ({e.__class__.__name__}); combine takes finished runs "
                             "of this runner") from None
    part = {"work": work, "spec": load(SPEC_COPY), "inputs": load(INPUTS), "manifest": load("k/manifest.json"),
            "groups": load("species-groups.json"), "summary": load("summary.json")}
    try:
        text = (work / RAW_MANIFEST).read_bytes()
    except OSError as e:
        raise SystemExit(f"part {work}: no readable {RAW_MANIFEST} ({e.__class__.__name__})") from None
    if sha256_bytes(text) != part["inputs"].get("raw_manifest_sha256"):
        raise SystemExit(f"part {work}: {RAW_MANIFEST} is not the RAW manifest its {INPUTS} hashed")
    part["raw"] = [line.split("  ", 2) for line in text.decode("utf-8").split("\n") if line]
    return part


def read_master(path):
    """A sha256sum-format manifest ("<sha256>  <path>" or "<sha256> *<path>") as {path: sha256}."""
    out = {}
    for n, line in enumerate(Path(path).read_text(encoding="utf-8").splitlines(), 1):
        if not line.strip():
            continue
        digest, sep, name = line[:64], line[64:66], line[66:]
        if not re.fullmatch(r"[0-9a-f]{64}", digest) or sep not in ("  ", " *") or not name:
            raise SystemExit(f"master manifest {path}, line {n}: not a sha256sum line")
        out[name.removeprefix("./")] = digest
    return out


def check_parts(spec_path, spec, parts, master_path):
    """The refusals of `combine` (G2): see the module docstring."""
    def without_datasets(s):
        return {k: v for k, v in s.items() if k != "datasets"}

    entries = {d["id"]: d for d in spec["datasets"]}
    owner = {}
    for part in parts:
        s = part["spec"]
        if not isinstance(s, dict) or canonical(without_datasets(s)) != canonical(without_datasets(spec)):
            differ = sorted(k for k in set(s) | set(spec) if k != "datasets" and s.get(k) != spec.get(k)) \
                if isinstance(s, dict) else ["everything"]
            raise SystemExit(f"part {part['work']}: its spec differs from {spec_path} in {', '.join(differ)}, not only "
                             "in its datasets")
        if part["manifest"].get("protocol_sha256") != sha256_bytes(canonical(s).encode()):
            raise SystemExit(f"part {part['work']}: K's manifest was not made from its {SPEC_COPY}")
        for d in s.get("datasets", []):
            if entries.get(d.get("id")) != d:
                raise SystemExit(f"part {part['work']}: dataset {d.get('id')} with partition {d.get('partition')} is "
                                 f"not in {spec_path}")
            if d["id"] in owner:
                raise SystemExit(f"dataset {d['id']} is in two parts: {owner[d['id']]} and {part['work']}")
            owner[d["id"]] = part["work"]
    left = [ds for ds in entries if ds not in owner]
    if left:
        raise SystemExit(f"no part holds {', '.join(left)} of {spec_path}")

    def same(label, get):
        if len({canonical(get(p)) for p in parts}) != 1:
            raise SystemExit(f"the parts differ in {label}: "
                             + "; ".join(f"{p['work']}: {get(p)}" for p in parts))
    same("the species pool (species-groups.json but its structure path)",
         lambda p: {k: v for k, v in p["groups"].items() if k != "structure"})
    same("the grm binary (inputs.json grm_sha256)", lambda p: p["inputs"].get("grm_sha256"))
    same("the vocabulary (inputs.json vocabulary_sha256)", lambda p: p["inputs"].get("vocabulary_sha256"))
    same("the source commit (K's manifest)", lambda p: p["manifest"].get("source_commit"))
    same("synthetic and the revision-12 gate (K's manifest)",
         lambda p: [p["manifest"].get("synthetic"), p["manifest"].get("revision12_gate")])
    same("eligible genomes and the panel floor (summary.json)",
         lambda p: [p["summary"].get("eligible_genomes"), p["summary"].get("panel_floor")])
    for p in parts:
        if p["manifest"].get("binary_sha256") != p["inputs"].get("grm_sha256"):
            raise SystemExit(f"part {p['work']}: K's binary_sha256 is not the grm its {INPUTS} hashed")
    # The one pinned master manifest (REVIEW-H1CODE L6): the spec names its sha256 in data_root.
    pinned = (spec.get("data_root") or {}).get("manifest_sha256")
    if pinned is None:
        raise SystemExit(f"{spec_path} pins no data_root.manifest_sha256, so combine can't tell which master manifest "
                         "the parts must come from")
    if file_sha256(master_path) != pinned:
        raise SystemExit(f"master manifest {master_path} is {file_sha256(master_path)}, not the {pinned} "
                         f"{spec_path} pins")
    master = read_master(master_path)
    for p in parts:
        outside = [name for digest, _, name in p["raw"] if master.get(name) != digest]
        if outside:
            raise SystemExit(f"part {p['work']}: {len(outside)} of its {len(p['raw'])} RAW files are not in the master "
                             f"manifest {master_path} with that sha256, e.g. {', '.join(outside[:5])}")


def verify_arm_output(part, src, arm, dataset):
    """A part's arm output as `combine` takes it (REVIEW-H1CODE L2): its resume.json names this arm and dataset, was
    made from the part's inputs.json and spec, lists runs.json, predictions.jsonl and arm.json, and the folder holds
    exactly the files it lists, none a symlink, each with its recorded sha256. Returns the record."""
    def refuse(why):
        raise SystemExit(f"part {part['work']}: arm {arm} on {dataset}: {why}")
    try:
        rec = json.loads((src / RESUME).read_text())
        outputs = rec["outputs"]
    except (OSError, ValueError, KeyError, TypeError):
        refuse(f"no readable {RESUME} with its outputs")
    if not isinstance(rec, dict) or not isinstance(outputs, dict):
        refuse(f"{RESUME} is not a resume record")
    if (rec.get("arm"), rec.get("dataset")) != (arm, dataset):
        refuse(f"{RESUME} is the record of {rec.get('arm')}/{rec.get('dataset')}")
    if rec.get("inputs_sha256") != file_sha256(part["work"] / INPUTS):
        refuse(f"{RESUME} was not made from the part's {INPUTS}")
    if rec.get("spec_sha256") != part["manifest"].get("protocol_sha256"):
        refuse(f"{RESUME} was not made from the part's spec")
    if not ARM_OUTPUTS <= set(outputs):
        refuse(f"{RESUME} lists no {', '.join(sorted(ARM_OUTPUTS - set(outputs)))}")
    if any(p.is_symlink() for p in src.rglob("*")):
        refuse("its folder holds a symlink")
    files = Runner.output_files(src)
    if files != set(outputs):
        refuse(f"its files differ from {RESUME}'s list: {', '.join(sorted(files ^ set(outputs)))}")
    changed = sorted(n for n, h in outputs.items() if file_sha256(src / n) != h)
    if changed:
        refuse(f"{', '.join(changed)} changed after the arm finished")
    return rec


def combine(argv):
    """One J record set and one K report over the dataset parts of one spec (G2)."""
    p = argparse.ArgumentParser(prog="run.py combine", description=combine.__doc__)
    p.add_argument("spec", help="the full spec; each part ran it with only its own datasets")
    p.add_argument("--work-from", action="append", required=True, metavar="WORK",
                   help="a finished part's WORK; once per part")
    p.add_argument("--master-manifest", required=True,
                   help="the pinned master manifest (sha256sum lines) that holds every part's RAW files")
    p.add_argument("--out", required=True, help="a new or empty folder for the combined outputs")
    a = p.parse_args(argv)
    spec, default, arms = load_spec(a.spec)
    out = Path(a.out)
    if out.exists() and (not out.is_dir() or any(out.iterdir())):
        raise SystemExit(f"--out {out} is not an empty folder")
    parts = [read_part(Path(w)) for w in a.work_from]
    check_parts(a.spec, spec, parts, a.master_manifest)
    names = [DEFAULT_ARM, *(arm["name"] for arm in arms)]
    flags = {}
    for part in parts:
        for d in part["spec"]["datasets"]:
            for name in names:
                src = part["work"] / "arms" / name / d["id"]
                if not (src / "runs.json").is_file():
                    raise SystemExit(f"part {part['work']} has no output of arm {name} on {d['id']}")
                rec = verify_arm_output(part, src, name, d["id"])
                flags.setdefault((rec["trace_training"], rec["run_identical_cohorts"]), f"{part['work']}: {name}/{d['id']}")
    if len(flags) > 1:
        raise SystemExit("the parts' arms ran with different --trace-training / --run-identical-cohorts: "
                         + "; ".join(f"{where} ran with {list(f)}" for f, where in sorted(flags.items())))
    for part in parts:
        for d in part["spec"]["datasets"]:
            for name in names:
                shutil.copytree(part["work"] / "arms" / name / d["id"], out / "arms" / name / d["id"], symlinks=True)
    first = parts[0]
    seconds = {"parts": sum(p["manifest"]["compute"]["actual_seconds"] for p in parts)}
    meta = {k: first["manifest"][k] for k in ("synthetic", "revision12_gate", "species_pool") if k in first["manifest"]}
    record = {"spec_sha256": sha256_bytes(canonical(spec).encode()),
              "master_manifest": str(Path(a.master_manifest).resolve()),
              "master_manifest_sha256": file_sha256(a.master_manifest),
              "grm_sha256": first["inputs"]["grm_sha256"], "source_commit": first["manifest"]["source_commit"],
              "parts": [{"work": str(p["work"].resolve()), "spec_sha256": p["manifest"]["protocol_sha256"],
                         "datasets": [d["id"] for d in p["spec"]["datasets"]],
                         "raw_manifest_sha256": p["inputs"]["raw_manifest_sha256"],
                         "species_groups_sha256": p["groups"]["sha256"]} for p in parts]}
    summary = report_run(spec, default, arms, out, first["manifest"]["binary_sha256"], record["source_commit"],
                         seconds, meta)
    (out / "combine.json").write_text(json.dumps(record, indent=1, sort_keys=True))
    summary.update({"seconds": round(sum(p["summary"].get("seconds", 0) for p in parts), 1), "stage_seconds": {},
                    "eligible_genomes": first["summary"]["eligible_genomes"],
                    "panel_floor": first["summary"]["panel_floor"], "species_groups": first["groups"],
                    "stale_arm_folders": [], "resumed_arms": [], "combine": record,
                    "revision12_gate": first["summary"].get("revision12_gate")})
    write_summary(out, summary)


if __name__ == "__main__":
    main(sys.argv[1:])
