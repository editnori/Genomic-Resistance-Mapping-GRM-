#!/usr/bin/env python3
"""Scores the outside-study test: a retrospective test on newly deposited
studies from the same database (docs/outside-study-test.md).

Reads the run script's work directory (its receipts and the fetch,
species-check, count, known-gene and predict outputs), the frozen population
tables in scripts/outside-study/frozen/ (candidates.tsv.gz, pairs.tsv.gz and
datasets.tsv), and the served revision-12 dataset reports. Writes TSV tables,
results.json with every number, and reading.txt. Standard library only; the
same inputs give byte-identical outputs.

The run script (run.py) loads this file by path for the shared pieces: the
report pins, the final models replayed from published bits, and the
known-gene rule. Nothing here fits, tunes or recalibrates a model.
"""

import argparse
import base64
import collections
import gzip
import hashlib
import json
import math
import os
import re
import subprocess
import sys
from pathlib import Path

NAME = "a retrospective test on newly deposited studies from the same database"
CONTRACT = "grm-predict-1"
REVISION = 12
REPORT_SCHEMA = 2
CUTOFF = 0.5
MIN_CLASS = 25
ALPHA = 0.05
MODEL_LEARNERS = ("scm", "cart", "lasso")
KNOWN = "known_gene"
LEARNERS = MODEL_LEARNERS + (KNOWN,)
# The report's learner row each learner here is compared with.
REPORT_LEARNER = {"scm": "scm", "cart": "cart", "lasso": "lasso", KNOWN: "catalog"}
SPLITS = ("lineage", "random")
METRICS_READ = ("sensitivity", "specificity")
GROUPINGS = ("bioproject", "country", "mlst_in_cohort")
CALLS_VERSION_SUFFIX = "|group-calls-v2"
POPULATIONS = {
    "primary": {
        "candidates": "primary (other-study)",
        "pair_flag": "in_frozen_test",
        "dataset_flag": "eligible_primary",
    },
    "same-study": {
        "candidates": "same-study",
        "pair_flag": "in_same_study_set",
        "dataset_flag": "eligible_same_study",
    },
}
# Exclusions after sequence, label-blind, in this order (docs, item 3).
EXCLUSION_ORDER = ("no sequence", "species check failed", "not counted")
FORBIDDEN_WORDS = ("validated", "generalizes", "clinical", "proves", "proven")
STAGES = ("setup", "fetch", "species-check", "count", "known-gene", "predict")
# The kit root is where this file is; no option or argument points it
# elsewhere.
KIT_ROOT = Path(__file__).resolve().parents[2]
KIT = KIT_ROOT / "scripts" / "outside-study"
SCORE_REL = "scripts/outside-study/score.py"
RUN_REL = "scripts/outside-study/run.py"
FROZEN_REL = "scripts/outside-study/FROZEN.sha256"
FROZEN_FILES = (
    "docs/outside-study-test.md",
    SCORE_REL,
    RUN_REL,
    "scripts/outside-study/reports.sha256",
    "scripts/outside-study/tables.sha256",
    "scripts/outside-study/frozen/RULES.md",
)
# git is never looked up on PATH: the first of these that exists.
GIT_PATHS = ("/usr/bin/git", "/bin/git", "/usr/local/bin/git")
COMMIT = re.compile(r"^[0-9a-f]{40}$")
TABLE_FILES = ("candidates.tsv", "pairs.tsv", "datasets.tsv")
# How each table is stored in the kit's frozen/ folder.
TABLE_STORED = {"candidates.tsv": "candidates.tsv.gz", "pairs.tsv": "pairs.tsv.gz", "datasets.tsv": "datasets.tsv"}
GENOME_ID = re.compile(r"^[0-9]+\.[0-9]+$")
NO_BIOPROJECT = "no BioProject recorded"
NO_COUNTRY = "no country recorded"
REFERENCE_KIND = {KNOWN: "fixed rule, whole cohort"}
CV_REFERENCE = "cross-validation, lineage split"


class Refusal(Exception):
    """An input the plan does not allow; nothing is written."""


def refuse(message):
    raise Refusal(message)


# ------------------------------------------------------------------ files


def sha256_bytes(data):
    return hashlib.sha256(data).hexdigest()


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_tsv(path):
    """Rows of a tab-separated file with a header, as the frozen tables and run.py have
    them: no quoting, tabs never inside a value."""
    with open(path, newline="") as f:
        return tsv_rows(f.read(), path)


def tsv_rows(text, path):
    lines = text.split("\n")
    if lines and lines[-1] == "":
        lines.pop()
    if not lines:
        refuse(f"{path}: empty")
    header = lines[0].split("\t")
    rows = []
    for n, line in enumerate(lines[1:], start=2):
        values = line.split("\t")
        if len(values) != len(header):
            refuse(f"{path}:{n}: {len(values)} fields, header has {len(header)}")
        rows.append(dict(zip(header, values)))
    return rows


def read_pins(path):
    """`<sha256>  <name>` lines, as sha256sum writes them."""
    return parse_pins(Path(path).read_text(), path)


def parse_pins(text, path):
    pins = {}
    for line in text.splitlines():
        if not line.strip():
            continue
        sha, name = line.split(None, 1)
        name = name.lstrip("*").strip()
        if name in pins:
            refuse(f"{path}: {name} pinned twice")
        pins[name] = sha
    return pins


def genome_key(g):
    taxon, _, suffix = g.partition(".")
    return (int(taxon), int(suffix)) if taxon.isdigit() and suffix.isdigit() else (1 << 62, 0, g)


def fmt(v):
    if v is None:
        return ""
    if isinstance(v, bool):
        return "yes" if v else "no"
    if isinstance(v, float):
        return f"{v:.6f}"
    return str(v)


def write_tsv(path, header, rows):
    out = ["\t".join(header)]
    for r in rows:
        cells = [fmt(r.get(h)) for h in header]
        for c in cells:
            if "\t" in c or "\n" in c:
                refuse(f"{path}: a value holds a tab or newline: {c!r}")
        out.append("\t".join(cells))
    Path(path).write_text("\n".join(out) + "\n")


def parse_json(text, where):
    """JSON text, refused (not a traceback) when it does not parse."""
    try:
        return json.loads(text)
    except json.JSONDecodeError as e:
        refuse(f"{where}: not JSON ({e})")


def dumps(obj):
    return json.dumps(obj, sort_keys=True, indent=1, allow_nan=False) + "\n"


# ------------------------------------------------------------------ gate


def git_binary():
    for path in GIT_PATHS:
        if os.path.isfile(path) and os.access(path, os.X_OK):
            return path
    refuse(f"no git at {', '.join(GIT_PATHS)}; the gate never looks git up on PATH")


def git_env():
    """The whole environment git runs with: no inherited GIT_* variable,
    no system or user configuration."""
    return {"PATH": "/usr/bin:/bin", "HOME": os.environ.get("HOME") or "/", "LC_ALL": "C",
            "GIT_CONFIG_NOSYSTEM": "1", "GIT_CONFIG_GLOBAL": "/dev/null", "GIT_OPTIONAL_LOCKS": "0"}


def git(*args):
    return subprocess.run([git_binary(), "-C", str(KIT_ROOT), "--no-replace-objects", *args],
                          capture_output=True, env=git_env())


def check_running(path, rel):
    """Refuses unless `path`, the file that is running, is the kit's own
    `rel`."""
    if Path(path).resolve() != KIT_ROOT / rel:
        refuse(f"the running script {Path(path).resolve()} is not the kit's own {KIT_ROOT / rel}")


def work_tree_bytes(rel):
    """A frozen file as it is on disk, refused when it or a folder above it
    is a symbolic link."""
    path = KIT_ROOT
    for part in Path(rel).parts:
        path = path / part
        if path.is_symlink():
            refuse(f"{path} is a symbolic link; a frozen path never is")
    if not path.is_file():
        refuse(f"{rel} is missing from the kit")
    return path.read_bytes()


def committed_bytes(rel):
    """A frozen file as committed at HEAD, read from the commit itself, so
    neither the index nor its flags can stand in for it."""
    tree = git("ls-tree", "HEAD", "--", rel)
    if tree.returncode != 0 or not tree.stdout.strip():
        refuse(f"{rel} is not committed at HEAD")
    if tree.stdout.split()[0] != b"100644" and tree.stdout.split()[0] != b"100755":
        refuse(f"{rel} is committed at HEAD as mode {tree.stdout.split()[0].decode()}, not as a file")
    blob = git("cat-file", "blob", f"HEAD:{rel}")
    if blob.returncode != 0:
        refuse(f"{rel} cannot be read from HEAD")
    return blob.stdout


def frozen_gate(expect_commit=None, rehearsal=False):
    """Refuses unless the code that is running is the kit's committed,
    pinned code, at the expected commit. In order: this file is the kit's
    score.py; git (an absolute path, a clean environment) puts the kit root
    at the top of its work tree; FROZEN.sha256 is the same at HEAD and on
    disk; every frozen file has its pinned SHA-256 both at HEAD and on disk;
    HEAD is `expect_commit`, which only a rehearsal may leave out. Returns
    the commit, the expected commit, the SHA-256 of FROZEN.sha256 and the
    pins. It stops an accident or a casual edit; it cannot stop someone who
    rewrites this function and commits that (docs/outside-study-test.md)."""
    if expect_commit is None and not rehearsal:
        refuse("a real run needs --expect-commit, the published commit of the kit")
    if expect_commit is not None and not COMMIT.match(expect_commit):
        refuse(f"--expect-commit {expect_commit!r} is not 40 lower-case hex characters")
    check_running(__file__, SCORE_REL)
    top = git("rev-parse", "--show-toplevel")
    if top.returncode != 0:
        refuse(f"{KIT_ROOT} is not in a git work tree: the plan and scoring script are not frozen")
    if Path(top.stdout.decode().strip()).resolve() != KIT_ROOT:
        refuse(f"git puts the top of the work tree at {top.stdout.decode().strip()}, not at the kit root {KIT_ROOT}")
    head = git("rev-parse", "--verify", "HEAD^{commit}")
    if head.returncode != 0:
        refuse("no git commit at the kit root")
    commit = head.stdout.decode().strip()
    frozen_disk = work_tree_bytes(FROZEN_REL) if (KIT_ROOT / FROZEN_REL).exists() else None
    if frozen_disk is None:
        refuse(f"{KIT_ROOT / FROZEN_REL} is missing: the plan and scoring script are not frozen")
    if committed_bytes(FROZEN_REL) != frozen_disk:
        refuse("FROZEN.sha256 differs from the one committed at HEAD")
    pins = parse_pins(frozen_disk.decode(), FROZEN_REL)
    for rel in pins:
        if rel not in FROZEN_FILES:
            refuse(f"FROZEN.sha256 pins {rel}, which is not a frozen file")
    for rel in FROZEN_FILES:
        if rel not in pins:
            refuse(f"FROZEN.sha256 does not pin {rel}")
        got = sha256_bytes(committed_bytes(rel))
        if got != pins[rel]:
            refuse(f"{rel} at HEAD has sha256 {got}, FROZEN.sha256 pins {pins[rel]}")
        got = sha256_bytes(work_tree_bytes(rel))
        if got != pins[rel]:
            refuse(f"{rel} has sha256 {got}, FROZEN.sha256 pins {pins[rel]}")
    if expect_commit is not None and commit != expect_commit:
        refuse(f"HEAD is {commit}, not the expected commit {expect_commit}")
    return {"commit": commit, "expect_commit": expect_commit, "frozen_sha256": sha256_bytes(frozen_disk),
            "files": {rel: pins[rel] for rel in FROZEN_FILES}}


def load_tables(rehearsal, tables_dir=None):
    """The frozen population. A real run reads the kit's own frozen/ folder
    and checks each table's uncompressed bytes against tables.sha256; only
    a rehearsal may name another folder of plain TSV files, whose pins are
    not checked."""
    raw = {}
    if tables_dir is None:
        tpins = read_pins(KIT / "tables.sha256")
        for name in TABLE_FILES:
            data = (KIT / "frozen" / TABLE_STORED[name]).read_bytes()
            if TABLE_STORED[name].endswith(".gz"):
                data = gzip.decompress(data)
            got = sha256_bytes(data)
            if got != tpins[name]:
                refuse(f"the frozen table {name} has sha256 {got}, the kit pins {tpins[name]}")
            raw[name] = data
        checked = True
    else:
        if not rehearsal:
            refuse("--tables is for a rehearsal only; a real run reads the kit's own frozen tables")
        for name in TABLE_FILES:
            raw[name] = (Path(tables_dir) / name).read_bytes()
        checked = False
    out = {"hashes": {f"tables/{n}": sha256_bytes(raw[n]) for n in TABLE_FILES}, "pins_checked": checked}
    for name in TABLE_FILES:
        out[name[: -len(".tsv")]] = tsv_rows(raw[name].decode(), name)
    return out


# ------------------------------------------------------------------ reports


def report_path(reports_dir, dataset):
    return Path(reports_dir) / f"{dataset}.json"


def load_report(reports_dir, dataset, pins):
    """A served revision-12 report, refused unless its bytes match the pin."""
    path = report_path(reports_dir, dataset)
    data = path.read_bytes()
    sha = sha256_bytes(data)
    want = pins.get(f"{dataset}.json")
    if want is None:
        refuse(f"{dataset}: no pinned report")
    if sha != want:
        refuse(f"{dataset}: report sha256 {sha} is not the pinned {want}")
    report = json.loads(data)
    if report.get("schema") != REPORT_SCHEMA:
        refuse(f"{dataset}: report schema {report.get('schema')}")
    if report["id"] != dataset:
        refuse(f"{path}: id {report['id']}")
    ev = report["evaluation"]
    if ev["revision"] != REVISION or not ev["valid"]:
        refuse(f"{dataset}: evaluation revision {ev['revision']} valid {ev['valid']}")
    return report, sha


def listed_features(report):
    return {f["feature"]: f for f in report["features"]}


def readable_kmer(kmer, k):
    return len(kmer) == k and set(kmer.upper()) <= set("ACGT")


def learner_features(report, learner):
    """The features a final model reads, in grm predict's order, or None
    when the report has no such model."""
    m = report["models"]
    if learner == "scm":
        return [r["feature"] for r in m["scm"]["rules"]]
    if learner == "cart":
        out = []
        for node in m["cart"]["nodes"]:
            s = node.get("split")
            if s is not None and node["depth"] < m["cart_depth"] and s["feature"] not in out:
                out.append(s["feature"])
        return out
    if learner == "lasso":
        lasso = m.get("lasso")
        return None if lasso is None else [c[0] for c in lasso["coef"]]
    raise ValueError(learner)


def complete_learners(report):
    """The learners whose published form is a complete model
    (docs/prediction-contract.md): every feature listed with a readable
    k-mer. Forests and boosted trees never are."""
    listed = listed_features(report)
    k = report["counts"]["k"]
    out = []
    for learner in MODEL_LEARNERS:
        ids = learner_features(report, learner)
        if ids is None:
            continue
        if all(f in listed and readable_kmer(listed[f]["kmer"], k) for f in ids):
            out.append(learner)
    return out


def final_call(report, learner, found):
    """The frozen final model applied to feature values (`found`: feature
    id -> bool). Returns (call, score) exactly as `grm predict` does."""
    m = report["models"]
    if learner == "scm":
        holds = [found[r["feature"]] == r["presence"] for r in m["scm"]["rules"]]
        if m["scm"]["kind"] == "conjunction":
            resistant = all(holds)
        elif m["scm"]["kind"] == "disjunction":
            resistant = any(holds)
        else:
            refuse(f"{report['id']}: SCM kind {m['scm']['kind']}")
        return ("resistant" if resistant else "susceptible"), None
    if learner == "cart":
        nodes, depth, i = m["cart"]["nodes"], m["cart_depth"], 0
        while True:
            node = nodes[i]
            s = node.get("split")
            if s is not None and node["depth"] < depth:
                i = s["present"] if found[s["feature"]] else s["absent"]
                continue
            score = node["value"]
            return ("resistant" if score >= CUTOFF else "susceptible"), score
    if learner == "lasso":
        lasso = m["lasso"]
        if not lasso["coef"]:
            return None, None
        eta = lasso["intercept"]
        for f, w in lasso["coef"]:
            if found[f]:
                eta += w
        score = 1.0 / (1.0 + math.exp(-eta))
        return ("resistant" if score >= CUTOFF else "susceptible"), score
    raise ValueError(learner)


def sample_bit(bits_b64, i):
    """Bit i of a feature's published presence bits: base64 of the
    little-endian u64 words, bit i for sample i."""
    data = base64.b64decode(bits_b64)
    return bool((data[i // 8] >> (i % 8)) & 1)


def check_row(report, sha, row, learners):
    """Refuses a prediction row that is not this report's frozen model under
    grm-predict-1, or whose call does not follow from its feature values."""
    d = report["id"]
    if row.get("contract") != CONTRACT:
        refuse(f"{d}: a row under contract {row.get('contract')!r}, not {CONTRACT}")
    if row.get("evaluation_revision") != REVISION:
        refuse(f"{d}: a row from evaluation revision {row.get('evaluation_revision')}")
    if row.get("report_sha256") != sha:
        refuse(f"{d}: a row names report sha256 {row.get('report_sha256')}, the pinned report is {sha}")
    if row.get("dataset") != d:
        refuse(f"{d}: a row for dataset {row.get('dataset')}")
    if row.get("k") != report["counts"]["k"] or row.get("input_fingerprint") != report["input_fingerprint"]:
        refuse(f"{d}: a row with another k or input fingerprint")
    learner = row.get("learner")
    if learner not in learners:
        refuse(f"{d}: a row for learner {learner!r}, scored learners are {list(learners)}")
    listed = listed_features(report)
    ids = learner_features(report, learner)
    got = [f["feature"] for f in row["features"]]
    if sorted(got) != sorted(set(ids)) or len(got) != len(set(got)):
        refuse(f"{d}: {row['genome']} {learner}: features {got}, the model reads {ids}")
    found = {}
    for f in row["features"]:
        if f["kmer"] != listed[f["feature"]]["kmer"]:
            refuse(f"{d}: feature {f['feature']} k-mer {f['kmer']} is not the published one")
        if not isinstance(f["found"], bool):
            refuse(f"{d}: feature {f['feature']} found {f['found']!r}")
        found[f["feature"]] = f["found"]
    call, score = final_call(report, learner, found)
    if row.get("call") != call:
        refuse(f"{d}: {row['genome']} {learner}: call {row.get('call')}, the model gives {call}")
    if (score is None) != (row.get("score") is None) or (
        score is not None and abs(score - row["score"]) > 1e-12
    ):
        refuse(f"{d}: {row['genome']} {learner}: score {row.get('score')}, the model gives {score}")
    return found


def replay_problems(report, sha, rows):
    """On cohort genomes: every feature's value must be its published bit,
    and every call the final model's call on the published bits. Returns the
    mismatches (empty when the replay holds) and the genomes checked."""
    index = {s["id"]: i for i, s in enumerate(report["samples"])}
    listed = listed_features(report)
    learners = complete_learners(report)
    problems, genomes = [], set()
    for row in rows:
        g = row["genome"]
        if g not in index:
            problems.append(f"{g}: not a sample of {report['id']}")
            continue
        genomes.add(g)
        try:
            check_row(report, sha, row, learners)
        except Refusal as e:
            problems.append(str(e))
            continue
        i = index[g]
        bits = {f: sample_bit(listed[f]["bits"], i) for f in learner_features(report, row["learner"])}
        for f in row["features"]:
            if f["found"] != bits[f["feature"]]:
                problems.append(f"{g} feature {f['feature']}: found {f['found']}, published bit {bits[f['feature']]}")
        call, _ = final_call(report, row["learner"], bits)
        if row["call"] != call:
            problems.append(f"{g} {row['learner']}: call {row['call']}, final model on the bits {call}")
    for g in sorted(genomes, key=genome_key):
        have = sorted(r["learner"] for r in rows if r["genome"] == g)
        if have != sorted(learners):
            problems.append(f"{g}: rows for {have}, complete learners {learners}")
    return problems, len(genomes)


def read_rows(path):
    rows = []
    with open(path) as f:
        for n, line in enumerate(f, start=1):
            if line.strip():
                try:
                    rows.append(json.loads(line))
                except json.JSONDecodeError:
                    refuse(f"{path}:{n}: not JSON")
    return rows


# ------------------------------------------------------------------ known-gene rule


def token_match(token, cls, subclass):
    """`C:<class token>` or `S:<subclass token>` against a hit's class or
    subclass, split on '/' (`catalog::token_match`)."""
    kind, want = token[:2], token[2:]
    field = cls if kind == "C:" else subclass
    return any(t.strip() == want for t in field.split("/"))


def explain(rule, hits):
    """(all rule groups met, groups used): `catalog::explain`."""
    used, met = set(), True
    for group in rule:
        matching = [h for h in hits if any(token_match(t, h["class"], h["subclass"]) for t in group)]
        if not matching:
            met = False
        used.update(h["group"] for h in matching)
    return met, sorted(used)


def known_gene_replay(report, hits_by_genome, every_sample=True):
    """Recomputes, from `grm count`'s gene calls of the cohort's own
    genomes, every sample's `groups` and `catalog_resistant`, the
    `catalog.groups` table and the catalog's confusion counts. Returns the
    mismatches; empty means the replay holds. Only a rehearsal passes
    `every_sample=False`: it then checks the samples it has calls for and
    skips the cohort-wide table and counts."""
    cat = report.get("catalog")
    if cat is None:
        return ["the report publishes no catalog rule"]
    rule, names = cat["rule"], report["catalog_groups"]
    problems, table, conf = [], {}, {"tp": 0, "fp": 0, "tn": 0, "fn": 0}
    for s in report["samples"]:
        hits = hits_by_genome.get(s["id"])
        if hits is None:
            if every_sample:
                problems.append(f"{s['id']}: no gene calls counted")
            continue
        got = [h["group"] for h in hits]
        want = [names[i] for i in s["groups"]]
        if got != want:
            problems.append(f"{s['id']}: groups {got}, published {want}")
        call, _ = explain(rule, hits)
        if call != s["catalog_resistant"]:
            problems.append(f"{s['id']}: catalog_resistant {call}, published {s['catalog_resistant']}")
        key = ("tp" if call else "fn") if s["resistant"] else ("fp" if call else "tn")
        conf[key] += 1
        for h in hits:
            row = table.get(h["group"])
            if row is None:
                row = table[h["group"]] = {
                    "group": h["group"],
                    "class": h["class"],
                    "subclass": h["subclass"],
                    "counts_for_rule": h["group"] in explain(rule, [h])[1],
                    "n_resistant": 0,
                    "n_susceptible": 0,
                }
            row["n_resistant" if s["resistant"] else "n_susceptible"] += 1
    if not problems and every_sample:
        rows = [table[g] for g in sorted(table)]
        rows.sort(key=lambda r: -(r["n_resistant"] + r["n_susceptible"]))
        if rows != cat["groups"]:
            problems.append("the catalog.groups table does not replay")
        for key, v in conf.items():
            if cat["metrics"][key] != v:
                problems.append(f"catalog {key} {v}, published {cat['metrics'][key]}")
    return problems


# ------------------------------------------------------------------ intervals


_LANCZOS = (
    0.9999999999998099,
    676.5203681218851,
    -1259.1392167224028,
    771.3234287776531,
    -176.6150291621406,
    12.507343278686905,
    -0.13857109526572012,
    9.984369578019572e-6,
    1.5056327351493116e-7,
)


def ln_gamma(x):
    """ln Γ(x), the Lanczos series of src/intervals.rs."""
    if x < 0.5:
        return math.log(math.pi / math.sin(math.pi * x)) - ln_gamma(1.0 - x)
    x = x - 1.0
    a = _LANCZOS[0]
    t = x + 7.0 + 0.5
    for i in range(1, 9):
        a += _LANCZOS[i] / (x + i)
    return 0.5 * math.log(2.0 * math.pi) + (x + 0.5) * math.log(t) - t + math.log(a)


def _beta_cf(a, b, x):
    tiny = 1e-300
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < tiny:
        d = tiny
    d = 1.0 / d
    h = d
    for m in range(1, 10001):
        m = float(m)
        m2 = 2.0 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < tiny:
            d = tiny
        c = 1.0 + aa / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < tiny:
            d = tiny
        c = 1.0 + aa / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        dl = d * c
        h *= dl
        if abs(dl - 1.0) < 1e-16:
            break
    return h


def beta_cdf(x, a, b):
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    front = math.exp(
        ln_gamma(a + b) - ln_gamma(a) - ln_gamma(b) + a * math.log(x) + b * math.log(1.0 - x)
    )
    if x < (a + 1.0) / (a + b + 2.0):
        return front * _beta_cf(a, b, x) / a
    return 1.0 - front * _beta_cf(b, a, 1.0 - x) / b


def beta_quantile(p, a, b):
    lo, hi = 0.0, 1.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if mid <= lo or mid >= hi:
            break
        if beta_cdf(mid, a, b) < p:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def clopper_pearson(x, n, alpha=ALPHA):
    """Two-sided exact interval, as `intervals::clopper_pearson`; None when
    n is 0."""
    if n == 0 or x > n:
        return None
    xf, nf = float(x), float(n)
    lo = 0.0 if x == 0 else beta_quantile(alpha / 2.0, xf, nf - xf + 1.0)
    hi = 1.0 if x == n else beta_quantile(1.0 - alpha / 2.0, xf + 1.0, nf - xf)
    return lo, hi


# ------------------------------------------------------------------ metrics


def ratio(a, b):
    return a / b if b else None


def confusion(pairs):
    """pairs: (resistant label, call) with call 'resistant', 'susceptible'
    or None (no call)."""
    c = {"tp": 0, "fn": 0, "fp": 0, "tn": 0, "no_call_resistant": 0, "no_call_susceptible": 0}
    for resistant, call in pairs:
        if call is None:
            c["no_call_resistant" if resistant else "no_call_susceptible"] += 1
        elif resistant:
            c["tp" if call == "resistant" else "fn"] += 1
        else:
            c["fp" if call == "resistant" else "tn"] += 1
    return c


def estimates(c, with_intervals):
    sens = ratio(c["tp"], c["tp"] + c["fn"])
    spec = ratio(c["tn"], c["tn"] + c["fp"])
    out = {
        "sensitivity": sens,
        "specificity": spec,
        "balanced_accuracy": None if sens is None or spec is None else (sens + spec) / 2.0,
        "precision": ratio(c["tp"], c["tp"] + c["fp"]),
        "vme": ratio(c["fn"], c["tp"] + c["fn"]),
        "me": ratio(c["fp"], c["fp"] + c["tn"]),
    }
    if with_intervals:
        for name, x, n in (
            ("sensitivity", c["tp"], c["tp"] + c["fn"]),
            ("specificity", c["tn"], c["tn"] + c["fp"]),
        ):
            ci = clopper_pearson(x, n)
            out[f"{name}_lo"] = None if ci is None else ci[0]
            out[f"{name}_hi"] = None if ci is None else ci[1]
    return out


def reading(reference, lo, hi):
    """The fixed descriptive reading: 'below' when the report's estimate is
    above the new interval, 'above' when it is under it."""
    if reference > hi:
        return "below"
    if reference < lo:
        return "above"
    return "consistent"


def reference_row(report, split, learner):
    """The report's cross-validation row for a learner in one split, or the
    reason there is none."""
    s = next((s for s in report["splits"] if s["scheme"] == split), None)
    if s is None:
        return None, f"the report has no {split} split"
    if s.get("unavailable_reason"):
        return None, f"{split} split unavailable: {s['unavailable_reason']}"
    row = next((r for r in s["learners"] if r["learner"] == REPORT_LEARNER[learner]), None)
    if row is None or row.get("metrics") is None:
        return None, f"no {REPORT_LEARNER[learner]} row in the {split} split"
    m = row["metrics"]
    return {
        "n": m["n"], "tp": m["tp"], "fn": m["fn"], "fp": m["fp"], "tn": m["tn"],
        "sensitivity": m["sensitivity"], "specificity": m["specificity"],
    }, None


# ------------------------------------------------------------------ the frozen population


def population(candidates, species, which):
    want = POPULATIONS[which]["candidates"]
    out = {}
    for c in candidates:
        if c["species_slug"] == species and c["population"] == want and c["on_eligible_dataset"] == "yes":
            g = c["genome_id"]
            if g in out:
                refuse(f"candidates: {g} listed twice")
            out[g] = c
    if not out:
        refuse(f"no {which} genome of {species} on an eligible dataset")
    return out


def eligible(dataset_rows, species, which):
    flag = POPULATIONS[which]["dataset_flag"]
    return sorted(d["dataset_id"] for d in dataset_rows if d["species_slug"] == species and d[flag] == "yes")


def labels(pairs, pop, which, datasets):
    """(dataset, genome) -> True for R, False for S. Only reported R and S
    calls are labels; intermediate, conflict and no-call pairs carry none.
    Refuses pairs whose frozen-set flag disagrees with the call or the
    population, which would mean the two populations were mixed."""
    flag = POPULATIONS[which]["pair_flag"]
    ds = set(datasets)
    out, seen = {}, set()
    for p in pairs:
        key = (p["dataset_id"], p["genome_id"])
        if key in seen:
            refuse(f"pairs: {key} twice")
        seen.add(key)
        if p["dataset_id"] not in ds:
            continue
        in_pop = p["genome_id"] in pop
        if p[flag] == "yes" and not in_pop:
            refuse(f"pairs: {key} is in the {which} frozen set but its genome is not in the {which} population")
        if not in_pop:
            continue
        if p["call"] in ("R", "S"):
            if p[flag] != "yes":
                refuse(f"pairs: {key} has call {p['call']} but is not flagged {flag}")
            out[key] = p["call"] == "R"
        elif p[flag] == "yes":
            refuse(f"pairs: {key} has call {p['call']!r} yet is flagged {flag}")
    return out


def exclusions(ids, availability, verdicts, counted):
    """Label-blind, in EXCLUSION_ORDER. Returns genome -> (reason, detail)
    and the borderline genomes, which stay in."""
    out, borderline = {}, set()
    for g in ids:
        a = availability.get(g)
        if a is None:
            refuse(f"the fetch record has no row for {g}")
        if a["fasta"] != "yes":
            out[g] = ("no sequence", a["detail"])
            continue
        v = verdicts.get(g)
        if v is None:
            refuse(f"the species check has no verdict for {g}, whose sequence was fetched")
        if v["verdict"] == "fail":
            out[g] = ("species check failed", v["reason"])
            continue
        if v["verdict"] == "borderline":
            borderline.add(g)
        elif v["verdict"] != "pass":
            refuse(f"species check verdict {v['verdict']!r} for {g}")
        c = counted.get(g)
        if c is None:
            refuse(f"the count record has no row for {g}")
        if c["counted"] != "yes":
            out[g] = ("not counted", c["detail"])
    return out, borderline


def cohort_sts(reports):
    """MLST types of the species' cohort: every sample of every report."""
    return {s["mlst"].strip() for r in reports for s in r["samples"] if s.get("mlst") and s["mlst"].strip()}


def no_st(value):
    """An MLST value that names no sequence type: empty, or whose ST field
    (after the last dot) is '-' or empty."""
    v = value.strip()
    return not v or v.rsplit(".", 1)[-1].strip() in ("", "-")


def largest(counter):
    """(name, count) with the largest count, ties to the first name; None
    when empty."""
    if not counter:
        return None
    return min(counter.items(), key=lambda kv: (-kv[1], kv[0]))


# ------------------------------------------------------------------ scoring


def score(inputs):
    """All numbers, from loaded inputs (see `load_inputs`). Pure: no file
    is read or written here."""
    which, species = inputs["population"], inputs["species"]
    pop = population(inputs["candidates"], species, which)
    ids = sorted(pop, key=genome_key)
    datasets = eligible(inputs["datasets"], species, which)
    if not datasets:
        refuse(f"no dataset of {species} is eligible in the {which} population")
    reports = inputs["reports"]
    for d in datasets:
        if d not in reports:
            refuse(f"no report for eligible dataset {d}")
    lab = labels(inputs["pairs"], pop, which, datasets)
    requested = set(inputs["availability"])
    if requested != set(pop):
        extra = sorted(requested - set(pop), key=genome_key)
        missing = sorted(set(pop) - requested, key=genome_key)
        if extra:
            refuse(f"genome {extra[0]} is outside the frozen list ({len(extra)} such)")
        refuse(f"frozen genome {missing[0]} was never requested ({len(missing)} such)")
    excluded, borderline = exclusions(ids, inputs["availability"], inputs["verdicts"], inputs["counted"])
    survivors = [g for g in ids if g not in excluded]
    if inputs["predict_survivors"] != survivors:
        refuse(f"predict/survivors.txt lists {len(inputs['predict_survivors'])} genomes, not the "
               f"{len(survivors)} survivors of the availability, species-check and count records")
    sts = cohort_sts(inputs["species_reports"])

    def group_of(g, grouping):
        c = pop[g]
        if grouping == "bioproject":
            return c["bioproject"].strip() or NO_BIOPROJECT
        if grouping == "country":
            return c["isolation_country"].strip() or NO_COUNTRY
        st = c["mlst"].strip()
        return "no ST" if no_st(st) else ("yes" if st in sts else "no")

    genomes = []
    for g in ids:
        c = pop[g]
        genomes.append({
            "genome_id": g,
            "bioproject": group_of(g, "bioproject"),
            "country": group_of(g, "country"),
            "mlst": c["mlst"].strip(),
            "mlst_in_cohort": group_of(g, "mlst_in_cohort"),
            "fasta": inputs["availability"][g]["fasta"] == "yes",
            "species_check": inputs["verdicts"].get(g, {}).get("verdict"),
            "borderline": g in borderline,
            "scored": g not in excluded,
            "exclusion": excluded[g][0] if g in excluded else None,
            "exclusion_detail": excluded[g][1] if g in excluded else None,
        })

    # Prediction rows: checked, then indexed by (dataset, learner, genome).
    surv = set(survivors)
    calls, found_by = {}, {}
    learners_by = {}
    for d in datasets:
        report, sha = reports[d]
        learners = complete_learners(report)
        learners_by[d] = list(learners)
        for row in inputs["rows"].get(d, []):
            g = row.get("genome")
            if g not in pop:
                refuse(f"{d}: genome {g} is outside the frozen list")
            if g in excluded:
                refuse(f"{d}: excluded genome {g} ({excluded[g][0]}) has a prediction row")
            found = check_row(report, sha, row, learners)
            key = (d, row["learner"], g)
            if key in calls:
                refuse(f"{d}: two rows for {g} {row['learner']}")
            calls[key] = row["call"]
            found_by[key] = found
        if inputs["rehearsal"]:
            # Stand-in genomes are cohort samples: on each report they are
            # samples of, their rows must be the report's own bits and
            # final-model calls, and every labelled one must be among them.
            members = {s["id"] for s in report["samples"]}
            mine = [r for r in inputs["rows"].get(d, []) if r["genome"] in members]
            problems, _ = replay_problems(report, sha, mine)
            if problems:
                refuse(f"rehearsal replay on {d}: {problems[0]} ({len(problems)} mismatches)")
            for (dd, g) in lab:
                if dd == d and g in surv and g not in members:
                    refuse(f"rehearsal: {g} is labelled on {d} but is not one of its samples")
        for learner in learners:
            for g in survivors:
                if (d, learner, g) not in calls:
                    refuse(f"{d}: no {learner} row for scored genome {g}")
    for d in inputs["rows"]:
        if d not in datasets:
            refuse(f"prediction rows for {d}, which is not eligible in the {which} population")

    # Known-gene comparator.
    kg = inputs["known_gene"]
    known_status = {}
    for d in datasets:
        if kg["unavailable"] is not None:
            known_status[d] = f"not available: {kg['unavailable']}"
            continue
        status = kg["datasets"].get(d)
        if status is None:
            refuse(f"the known-gene record has no row for {d}")
        if status["replay"] != "pass":
            known_status[d] = status["replay"]
            continue
        rule = json.loads(status["rule"])
        report = reports[d][0]
        if report.get("catalog") is None or rule != report["catalog"]["rule"]:
            refuse(f"{d}: the known-gene rule is not the report's catalog.rule")
        known_status[d] = "scored"
        learners_by[d].append(KNOWN)
        for g in survivors:
            row = kg["calls"].get((d, g))
            if row is None:
                refuse(f"{d}: no known-gene call for scored genome {g}")
            calls[(d, KNOWN, g)] = "resistant" if row["catalog_resistant"] == "yes" else "susceptible"
    if kg["unavailable"] is None:
        for (d, g) in kg["calls"]:
            if g not in pop:
                refuse(f"{d}: a known-gene call for {g}, outside the frozen list")
            if g in excluded:
                refuse(f"{d}: a known-gene call for excluded genome {g}")

    dataset_rows, metric_rows, group_rows, loso_rows, feature_rows = [], [], [], [], []
    for d in datasets:
        report = reports[d][0]
        labelled = [g for g in ids if (d, g) in lab]
        scored = [g for g in labelled if g in surv]
        n_r = sum(lab[(d, g)] for g in scored)
        n_s = len(scored) - n_r
        reportable = n_r >= MIN_CLASS and n_s >= MIN_CLASS
        drow = {
            "dataset": d,
            "labelled_resistant": sum(lab[(d, g)] for g in labelled),
            "labelled_susceptible": sum(not lab[(d, g)] for g in labelled),
            "scored_resistant": n_r,
            "scored_susceptible": n_s,
            "reportable": reportable,
            "known_gene": known_status[d],
            "lasso_complete": "lasso" in learners_by[d],
        }
        for reason in EXCLUSION_ORDER:
            drow[f"excluded_{reason.replace(' ', '_')}"] = sum(
                1 for g in labelled if excluded.get(g, ("",))[0] == reason
            )
        drow["scored_borderline"] = sum(1 for g in scored if g in borderline)
        # Description only: how much of each class one study holds, and how
        # many leave-one-study-out folds keep 25 of each class.
        studies_d = sorted({group_of(g, "bioproject") for g in scored} - {NO_BIOPROJECT})
        for cls, name, n in ((True, "resistant", n_r), (False, "susceptible", n_s)):
            top = largest(collections.Counter(group_of(g, "bioproject") for g in scored
                                              if lab[(d, g)] == cls and group_of(g, "bioproject") != NO_BIOPROJECT))
            drow[f"top_study_{name}"] = None if top is None else top[0]
            drow[f"top_study_{name}_count"] = 0 if top is None else top[1]
            drow[f"top_study_{name}_share"] = None if top is None else ratio(top[1], n)
        drow["studies"] = len(studies_d)
        folds = 0
        for b in studies_d:
            rest_r = sum(1 for g in scored if lab[(d, g)] and group_of(g, "bioproject") != b)
            rest_s = sum(1 for g in scored if not lab[(d, g)] and group_of(g, "bioproject") != b)
            folds += rest_r >= MIN_CLASS and rest_s >= MIN_CLASS
        drow["loso_folds_reportable"] = folds
        dataset_rows.append(drow)
        for learner in learners_by[d]:
            c = confusion((lab[(d, g)], calls[(d, learner, g)]) for g in scored)
            row = {"dataset": d, "learner": learner, "reportable": reportable,
                   "scored_resistant": n_r, "scored_susceptible": n_s, **c,
                   "reference_kind": REFERENCE_KIND.get(learner, CV_REFERENCE)}
            if reportable:
                row.update(estimates(c, True))
            for split in SPLITS:
                ref, why = reference_row(report, split, learner)
                for k in ("n", "tp", "fn", "fp", "tn", "sensitivity", "specificity"):
                    row[f"{split}_{k}"] = None if ref is None else ref[k]
                row[f"{split}_note"] = why
            row["reference_at_0_or_1"] = any(
                row[f"{split}_{m}"] in (0.0, 1.0) for split in SPLITS for m in METRICS_READ
                if row[f"{split}_{m}"] is not None)
            for m in METRICS_READ:
                lo, hi = row.get(f"{m}_lo"), row.get(f"{m}_hi")
                ref = row["lineage_" + m]
                if not reportable:
                    row[f"reading_{m}"] = "not read: counts only"
                elif ref is None:
                    row[f"reading_{m}"] = "not read: no lineage-split estimate"
                elif lo is None:
                    row[f"reading_{m}"] = "not read: no interval"
                else:
                    row[f"reading_{m}"] = reading(ref, lo, hi)
            # Each BioProject left out in turn.
            # "no BioProject recorded" is never a fold and never a study.
            projects = sorted({group_of(g, "bioproject") for g in scored} - {NO_BIOPROJECT})
            loso = {m: [] for m in METRICS_READ}
            for b in projects:
                rest = [g for g in scored if group_of(g, "bioproject") != b]
                rc = confusion((lab[(d, g)], calls[(d, learner, g)]) for g in rest)
                rr = sum(lab[(d, g)] for g in rest)
                rs = len(rest) - rr
                lrow = {"dataset": d, "learner": learner, "left_out": b,
                        "scored_resistant": rr, "scored_susceptible": rs,
                        "reportable": rr >= MIN_CLASS and rs >= MIN_CLASS, **rc}
                if lrow["reportable"]:
                    e = estimates(rc, False)
                    for m in ("sensitivity", "specificity", "balanced_accuracy"):
                        lrow[m] = e[m]
                    for m in METRICS_READ:
                        if e[m] is not None:
                            loso[m].append(e[m])
                loso_rows.append(lrow)
            row["studies"] = len(projects)
            for m in METRICS_READ:
                row[f"loso_{m}_n"] = len(loso[m])
                row[f"loso_{m}_min"] = min(loso[m]) if loso[m] else None
                row[f"loso_{m}_max"] = max(loso[m]) if loso[m] else None
            metric_rows.append(row)
            for grouping in GROUPINGS:
                values = sorted({group_of(g, grouping) for g in scored})
                for v in values:
                    part = [g for g in scored if group_of(g, grouping) == v]
                    gc = confusion((lab[(d, g)], calls[(d, learner, g)]) for g in part)
                    pr = sum(lab[(d, g)] for g in part)
                    group_rows.append({"dataset": d, "learner": learner, "grouping": grouping,
                                       "group": v, "scored_resistant": pr,
                                       "scored_susceptible": len(part) - pr, **gc})
            if learner == KNOWN:
                continue
            listed = listed_features(report)
            seen = []
            for f in learner_features(report, learner):
                if f in seen:
                    continue
                seen.append(f)
                carry_r = sum(1 for g in scored if lab[(d, g)] and found_by[(d, learner, g)][f])
                carry_s = sum(1 for g in scored if not lab[(d, g)] and found_by[(d, learner, g)][f])
                lf = listed[f]
                feature_rows.append({
                    "dataset": d, "learner": learner, "feature": f, "kmer": lf["kmer"],
                    "n_kmers": lf["n_kmers"],
                    "cohort_resistant_carriers": lf["n_resistant"],
                    "cohort_resistant": report["counts"]["resistant"],
                    "cohort_resistant_share": ratio(lf["n_resistant"], report["counts"]["resistant"]),
                    "cohort_susceptible_carriers": lf["n_susceptible"],
                    "cohort_susceptible": report["counts"]["susceptible"],
                    "cohort_susceptible_share": ratio(lf["n_susceptible"], report["counts"]["susceptible"]),
                    "new_resistant_carriers": carry_r, "new_resistant": n_r,
                    "new_resistant_share": ratio(carry_r, n_r),
                    "new_susceptible_carriers": carry_s, "new_susceptible": n_s,
                    "new_susceptible_share": ratio(carry_s, n_s),
                })

    # Each learner beside the known-gene rule, genome by genome: counts
    # only, over the scored labelled genomes that have both calls.
    paired_rows = []
    for d in datasets:
        if known_status[d] != "scored":
            continue
        labelled = [g for g in ids if (d, g) in lab and g in surv]
        for learner in learners_by[d]:
            if learner == KNOWN:
                continue
            prow = {"dataset": d, "learner": learner}
            for cls, name in ((True, "resistant"), (False, "susceptible")):
                both = learner_only = known_only = neither = 0
                for g in labelled:
                    if lab[(d, g)] != cls or calls[(d, learner, g)] is None:
                        continue
                    a = calls[(d, learner, g)] == "resistant"
                    k = calls[(d, KNOWN, g)] == "resistant"
                    both += a and k
                    learner_only += a and not k
                    known_only += k and not a
                    neither += not a and not k
                prow.update({f"{name}_genomes": both + learner_only + known_only + neither,
                             f"{name}_both": both, f"{name}_learner_only": learner_only,
                             f"{name}_known_gene_only": known_only, f"{name}_neither": neither})
            paired_rows.append(prow)

    # Scored genomes with the same k-mer set as a cohort genome or as another
    # scored genome (grm count's kset digest). A list only: nobody is
    # excluded by it.
    dup_rows = []
    new_digest = inputs["digests"]["new"]
    by_digest = collections.defaultdict(list)
    for g, dg in inputs["digests"]["cohort"].items():
        by_digest[dg].append((g, "cohort"))
    for g in survivors:
        if g not in new_digest:
            refuse(f"the count record has no k-mer set digest for scored genome {g}")
        by_digest[new_digest[g]].append((g, "scored new genome"))
    for g in survivors:
        others = [(o, where) for o, where in by_digest[new_digest[g]] if o != g]
        for o, where in sorted(others, key=lambda x: (x[1], genome_key(x[0]))):
            dup_rows.append({"genome_id": g, "other_genome_id": o, "other_is": where, "kset_digest": new_digest[g]})
    duplicate_counts = {
        "scored genomes with the k-mer set of a cohort genome":
            len({r["genome_id"] for r in dup_rows if r["other_is"] == "cohort"}),
        "scored genomes with the k-mer set of another scored genome":
            len({r["genome_id"] for r in dup_rows if r["other_is"] == "scored new genome"}),
    }

    tallies = []
    for learner in LEARNERS:
        rows = [r for r in metric_rows if r["learner"] == learner]
        if not rows:
            continue
        t = {"learner": learner, "below": 0, "consistent": 0, "above": 0, "not_read": 0}
        for r in rows:
            for m in METRICS_READ:
                v = r[f"reading_{m}"]
                t[v if v in ("below", "consistent", "above") else "not_read"] += 1
        tallies.append(t)

    exclusion_counts = {reason: sum(1 for v in excluded.values() if v[0] == reason) for reason in EXCLUSION_ORDER}
    return {
        "name": NAME,
        "rehearsal": inputs["rehearsal"],
        "species": species,
        "population": which,
        "contract": CONTRACT,
        "evaluation_revision": REVISION,
        "cutoff": CUTOFF,
        "min_class": MIN_CLASS,
        "confidence": 1.0 - ALPHA,
        "datasets": datasets,
        "reports": {d: reports[d][1] for d in datasets},
        "genomes_frozen": len(ids),
        "genomes_scored": len(survivors),
        "genomes_borderline_scored": sum(1 for g in survivors if g in borderline),
        "exclusions": exclusion_counts,
        "genomes": genomes,
        "dataset_counts": dataset_rows,
        "metrics": metric_rows,
        "reading_tallies": tallies,
        "groups": group_rows,
        "leave_one_study_out": loso_rows,
        "features": feature_rows,
        "paired": paired_rows,
        "duplicates": dup_rows,
        "duplicate_counts": duplicate_counts,
        "attempts": inputs["attempts"],
        "table_pins_checked": inputs["table_pins_checked"],
        "grm_sha256": inputs["grm_sha256"],
        "inputs": inputs["hashes"],
    }


# ------------------------------------------------------------------ outputs

GENOME_COLS = ["genome_id", "bioproject", "country", "mlst", "mlst_in_cohort", "fasta",
               "species_check", "borderline", "scored", "exclusion", "exclusion_detail"]
DATASET_COLS = ["dataset", "labelled_resistant", "labelled_susceptible",
                "excluded_no_sequence", "excluded_species_check_failed", "excluded_not_counted",
                "scored_resistant", "scored_susceptible", "scored_borderline", "reportable",
                "lasso_complete", "known_gene",
                "top_study_resistant", "top_study_resistant_count", "top_study_resistant_share",
                "top_study_susceptible", "top_study_susceptible_count", "top_study_susceptible_share",
                "studies", "loso_folds_reportable"]
METRIC_COLS = (
    ["dataset", "learner", "reportable", "scored_resistant", "scored_susceptible",
     "tp", "fn", "fp", "tn", "no_call_resistant", "no_call_susceptible",
     "sensitivity", "sensitivity_lo", "sensitivity_hi",
     "specificity", "specificity_lo", "specificity_hi",
     "balanced_accuracy", "precision", "vme", "me", "reference_kind", "reference_at_0_or_1"]
    + [f"{s}_{k}" for s in SPLITS for k in ("n", "tp", "fn", "fp", "tn", "sensitivity", "specificity", "note")]
    + ["reading_sensitivity", "reading_specificity", "studies"]
    + [f"loso_{m}_{k}" for m in METRICS_READ for k in ("n", "min", "max")]
)
TALLY_COLS = ["learner", "below", "consistent", "above", "not_read"]
GROUP_COLS = ["dataset", "learner", "grouping", "group", "scored_resistant", "scored_susceptible",
              "tp", "fn", "fp", "tn", "no_call_resistant", "no_call_susceptible"]
LOSO_COLS = ["dataset", "learner", "left_out", "scored_resistant", "scored_susceptible", "reportable",
             "tp", "fn", "fp", "tn", "no_call_resistant", "no_call_susceptible",
             "sensitivity", "specificity", "balanced_accuracy"]
FEATURE_COLS = ["dataset", "learner", "feature", "kmer", "n_kmers",
                "cohort_resistant_carriers", "cohort_resistant", "cohort_resistant_share",
                "cohort_susceptible_carriers", "cohort_susceptible", "cohort_susceptible_share",
                "new_resistant_carriers", "new_resistant", "new_resistant_share",
                "new_susceptible_carriers", "new_susceptible", "new_susceptible_share"]
PAIRED_COLS = ["dataset", "learner"] + [
    f"{c}_{k}" for c in ("resistant", "susceptible")
    for k in ("genomes", "both", "learner_only", "known_gene_only", "neither")]
DUPLICATE_COLS = ["genome_id", "other_genome_id", "other_is", "kset_digest"]


def pct(v):
    return "n/a" if v is None else f"{100.0 * v:.1f}%"


def reading_text(res):
    lines = []
    head = f"{'REHEARSAL ON COHORT GENOMES: MEASURES NOTHING. ' if res['rehearsal'] else ''}"
    lines.append(f"{head}Outside-study test, {res['species']}, {res['population']} population")
    lines.append(f"This is {NAME}.")
    lines.append("It is not a lockbox: the batch's metadata and label totals were read before the rules were frozen.")
    lines.append("The genomes come from the same database as the cohort.")
    lines.append("")
    lines.append(f"Frozen genomes: {res['genomes_frozen']}. Scored: {res['genomes_scored']} "
                 f"({res['genomes_borderline_scored']} of them borderline in the species check, kept and flagged).")
    for reason in EXCLUSION_ORDER:
        lines.append(f"Genomes excluded, {reason}: {res['exclusions'][reason]}.")
    lines.append("")
    lines.append("Each model is the report's final model, fitted once on the whole cohort and applied unchanged.")
    lines.append("The report's figures are for models refitted in each cross-validation fold and scored on held-out cohort genomes.")
    lines.append(f"Calls use the {CUTOFF} cutoff. Intervals are exact Clopper-Pearson {int(round(100 * (1 - ALPHA)))}% intervals.")
    lines.append(f"A dataset with fewer than {MIN_CLASS} resistant or {MIN_CLASS} susceptible scored genomes is reported as counts only.")
    lines.append("Sensitivity and specificity are over the genomes a learner called; genomes with no call are counted in their own column.")
    lines.append("The known-gene rule's reference figure is a fixed rule applied to the whole cohort: nothing was fitted, "
                 "so there is no cross-validation figure for it.")
    lines.append("The datasets share genomes and studies, so they are not independent of each other.")
    lines.append("Precision depends on the share of resistant genomes, which differs between the cohort and the new genomes.")
    dc = res["duplicate_counts"]
    lines.append(f"Scored genomes with the same k-mer set as a cohort genome: "
                 f"{dc['scored genomes with the k-mer set of a cohort genome']}; as another scored genome: "
                 f"{dc['scored genomes with the k-mer set of another scored genome']} (listed in duplicates.tsv; none is excluded for it).")
    lines.append("Stage attempts: " + ", ".join(f"{s} {res['attempts'][s]}" for s in STAGES) + ".")
    lines.append("")
    for d in res["dataset_counts"]:
        lines.append(f"{d['dataset']}: {d['scored_resistant']} resistant and {d['scored_susceptible']} susceptible genomes scored"
                     f" ({'estimates' if d['reportable'] else 'counts only'}); known-gene rule: {d['known_gene']}.")
        shares = []
        for name in ("resistant", "susceptible"):
            if d[f"top_study_{name}"] is None:
                shares.append(f"no {name} genome is in a study")
            else:
                shares.append(f"{d[f'top_study_{name}']} holds {d[f'top_study_{name}_count']} of the "
                              f"{d[f'scored_{name}']} {name} genomes ({pct(d[f'top_study_{name}_share'])})")
        lines.append(f"  By study: {shares[0]}; {shares[1]}. Of {d['studies']} studies, "
                     f"{d['loso_folds_reportable']} can be left out with {MIN_CLASS} of each class remaining.")
        if d["loso_folds_reportable"] == 0:
            lines.append("  No leave-one-study-out fold is reportable, so no study-to-study range can be shown, "
                         "and the estimate is in effect one study's.")
        for r in res["metrics"]:
            if r["dataset"] != d["dataset"]:
                continue
            counts = f"TP {r['tp']}, FN {r['fn']}, FP {r['fp']}, TN {r['tn']}"
            if r["no_call_resistant"] or r["no_call_susceptible"]:
                counts += f", no call {r['no_call_resistant'] + r['no_call_susceptible']}"
            if not r["reportable"]:
                lines.append(f"  {r['learner']}: {counts}.")
                continue
            lines.append(
                f"  {r['learner']}: {counts}. Sensitivity {pct(r['sensitivity'])} "
                f"[{pct(r['sensitivity_lo'])}, {pct(r['sensitivity_hi'])}], lineage split {pct(r['lineage_sensitivity'])}, "
                f"random split {pct(r['random_sensitivity'])}: {r['reading_sensitivity']}. "
                f"Specificity {pct(r['specificity'])} [{pct(r['specificity_lo'])}, {pct(r['specificity_hi'])}], "
                f"lineage split {pct(r['lineage_specificity'])}, random split {pct(r['random_specificity'])}: "
                f"{r['reading_specificity']}."
            )
            if r["lineage_note"]:
                lines.append(f"    No lineage-split figure: {r['lineage_note']}.")
            lines.append(f"    Reference figure: {r['reference_kind']}"
                         f"{'; a reference sensitivity or specificity is exactly 0 or 1' if r['reference_at_0_or_1'] else ''}.")
            if r["loso_sensitivity_n"]:
                lines.append(
                    f"    With each of {r['studies']} BioProjects left out: sensitivity {pct(r['loso_sensitivity_min'])} to "
                    f"{pct(r['loso_sensitivity_max'])}, specificity {pct(r['loso_specificity_min'])} to "
                    f"{pct(r['loso_specificity_max'])} (over {r['loso_sensitivity_n']} estimates)."
                )
    lines.append("")
    lines.append("Readings per learner (each dataset contributes sensitivity and specificity):")
    for t in res["reading_tallies"]:
        lines.append(f"  {t['learner']}: below {t['below']}, consistent {t['consistent']}, above {t['above']}, not read {t['not_read']}.")
    lines.append("The reading is descriptive. It has no multiplicity correction and is not a hypothesis test.")
    lines.append("paired.tsv sets each learner beside the known-gene rule on the same genomes (counts only, no test).")
    lines.append("")
    lines.append("Not run here: sequence-distance novelty, partial-feature counts (--matrix), calibration,")
    lines.append("any hypothesis test, forests and boosted trees (not published in full).")
    text = "\n".join(lines) + "\n"
    for word in FORBIDDEN_WORDS:
        if word in text.lower():
            refuse(f"the reading would contain {word!r}")
    return text


OUTPUT_NAMES = ("genomes.tsv", "datasets.tsv", "metrics.tsv", "reading_counts.tsv", "groups.tsv",
                "leave_one_study_out.tsv", "features.tsv", "paired.tsv", "duplicates.tsv", "reading.txt",
                "results.json")


def output_name(res, name):
    """A rehearsal's every output starts with 'rehearsal-'."""
    return ("rehearsal-" if res["rehearsal"] else "") + name


def write_outputs(res, out):
    out = Path(out)
    text = reading_text(res)
    out.mkdir(parents=True, exist_ok=True)
    n = lambda name: out / output_name(res, name)  # noqa: E731
    write_tsv(n("genomes.tsv"), GENOME_COLS, res["genomes"])
    write_tsv(n("datasets.tsv"), DATASET_COLS, res["dataset_counts"])
    write_tsv(n("metrics.tsv"), METRIC_COLS, res["metrics"])
    write_tsv(n("reading_counts.tsv"), TALLY_COLS, res["reading_tallies"])
    write_tsv(n("groups.tsv"), GROUP_COLS, res["groups"])
    write_tsv(n("leave_one_study_out.tsv"), LOSO_COLS, res["leave_one_study_out"])
    write_tsv(n("features.tsv"), FEATURE_COLS, res["features"])
    write_tsv(n("paired.tsv"), PAIRED_COLS, res["paired"])
    write_tsv(n("duplicates.tsv"), DUPLICATE_COLS, res["duplicates"])
    n("reading.txt").write_text(text)
    n("results.json").write_text(dumps(res))


# ------------------------------------------------------------------ loading


def load_receipts(work, args, gate):
    """The six stage receipts, refused unless they come from one run of
    this kit at HEAD, with one grm binary, for these arguments. Returns
    the SHA-256 each work file had when its stage wrote it, the grm
    binary's SHA-256 and the number of attempts per stage."""
    receipts = {}
    for stage in STAGES:
        p = work / "receipts" / f"{stage}.json"
        if not p.is_file():
            refuse(f"the {stage} stage has no receipt; the run is not complete")
        receipts[stage] = parse_json(p.read_text(), f"receipts/{stage}.json")
    scripts = {"run.py": gate["files"][RUN_REL], "score.py": gate["files"][SCORE_REL]}
    grm = {r["grm"]["sha256"] for r in receipts.values()}
    if len(grm) != 1:
        refuse("the receipts name more than one grm binary")
    written = {}
    for stage, r in receipts.items():
        if r["kit"]["commit"] != gate["commit"]:
            refuse(f"the {stage} receipt names commit {r['kit']['commit']}, HEAD is {gate['commit']}")
        for name, sha in scripts.items():
            if r["kit"][name] != sha:
                refuse(f"the {stage} stage ran a {name} with sha256 {r['kit'][name]}, the kit pins {sha}")
        for key, want in (("species", args.species), ("population", args.population), ("rehearsal", args.rehearsal)):
            if r["args"][key] != want:
                refuse(f"the {stage} stage ran with {key} {r['args'][key]!r}, this scoring has {want!r}")
        for key, sha in r["outputs"].items():
            if key.startswith("work/"):
                written[key[len("work/"):]] = (stage, sha)
    attempts = {stage: 0 for stage in STAGES}
    log = work / "receipts" / "attempts.jsonl"
    if not log.is_file():
        refuse("the run left no receipts/attempts.jsonl")
    for line in log.read_text().splitlines():
        if line.strip():
            a = parse_json(line, "receipts/attempts.jsonl")
            if a["event"] == "start" and a["stage"] in attempts:
                attempts[a["stage"]] += 1
    return written, grm.pop(), attempts


def load_inputs(args, gate):
    work, reports_dir = Path(args.work).resolve(), Path(args.reports)
    tables = load_tables(args.rehearsal, args.tables)
    hashes = dict(tables["hashes"])
    pins = read_pins(KIT / "reports.sha256")
    species = args.species
    species_ids = sorted(n[: -len(".json")] for n in pins if n.startswith(species + "__"))
    species_reports, reports = [], {}
    for d in species_ids:
        r, sha = load_report(reports_dir, d, pins)
        species_reports.append(r)
        reports[d] = (r, sha)
        hashes[f"reports/{d}.json"] = sha
    candidates = tables["candidates"]
    if args.rehearsal:
        cohort = {s["id"] for r in species_reports for s in r["samples"]}
        for c in candidates:
            if c["genome_id"] not in cohort:
                refuse(f"rehearsal: {c['genome_id']} is not a cohort genome; a rehearsal scores cohort genomes only")
    else:
        cohort = {s["id"] for r in species_reports for s in r["samples"]}
        for c in candidates:
            if c["species_slug"] == species and c["genome_id"] in cohort:
                refuse(f"{c['genome_id']} is a cohort genome")

    written, grm_sha, attempts = load_receipts(work, args, gate)
    for stage in STAGES:
        hashes[f"work/receipts/{stage}.json"] = sha256_file(work / "receipts" / f"{stage}.json")

    def read(rel):
        """A work file's text, refused unless it has the SHA-256 the stage
        that wrote it recorded."""
        data = (work / rel).read_bytes()
        got = sha256_bytes(data)
        if rel not in written:
            refuse(f"work/{rel} is in no stage's receipt")
        stage, want = written[rel]
        if got != want:
            refuse(f"work/{rel} has sha256 {got}, the {stage} stage recorded {want}")
        hashes[f"work/{rel}"] = got
        return data.decode()

    def tsv(rel, key):
        return {r[key]: r for r in tsv_rows(read(rel), rel)}

    availability = tsv("fetch/availability.tsv", "genome_id")
    verdicts = tsv("species-check/verdicts.tsv", "genome_id")
    counted = tsv("count/counted.tsv", "genome_id")
    digests = {"new": {g: r["kset_digest"] for g, r in counted.items() if r["counted"] == "yes"},
               "cohort": {g: r["kset_digest"] for g, r in tsv("count/cohort-digests.tsv", "genome_id").items()}}
    rows = {}
    listed = sorted(rel for rel in written if rel.startswith("predict/") and rel.endswith(".jsonl"))
    present = sorted(f"predict/{p.name}" for p in (work / "predict").glob("*.jsonl"))
    if listed != present:
        refuse(f"the predict folder holds {present}, the predict stage recorded {listed}")
    for rel in present:
        d = rel[len("predict/"): -len(".jsonl")]
        rows[d] = [parse_json(line, f"work/{rel}") for line in read(rel).splitlines() if line.strip()]
    predict_survivors = [g for g in read("predict/survivors.txt").splitlines() if g.strip()]
    kg = {"unavailable": None, "datasets": {}, "calls": {}}
    if (work / "known-gene" / "unavailable.txt").exists():
        if not args.rehearsal:
            refuse("the known-gene stage recorded no calls; only a rehearsal may go on without them")
        kg["unavailable"] = read("known-gene/unavailable.txt").strip()
    else:
        kg["datasets"] = tsv("known-gene/datasets.tsv", "dataset")
        for r in tsv_rows(read("known-gene/calls.tsv"), "known-gene/calls.tsv"):
            key = (r["dataset"], r["genome_id"])
            if key in kg["calls"]:
                refuse(f"known-gene calls: {key} twice")
            kg["calls"][key] = r
    return {
        "population": args.population,
        "species": species,
        "rehearsal": args.rehearsal,
        "candidates": candidates,
        "pairs": tables["pairs"],
        "datasets": tables["datasets"],
        "reports": reports,
        "species_reports": species_reports,
        "availability": availability,
        "verdicts": verdicts,
        "counted": counted,
        "digests": digests,
        "rows": rows,
        "predict_survivors": predict_survivors,
        "known_gene": kg,
        "attempts": attempts,
        "table_pins_checked": tables["pins_checked"],
        "grm_sha256": grm_sha,
        "hashes": hashes,
    }


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--work", required=True, help="the run script's work directory")
    ap.add_argument("--reports", required=True, help="the served revision-12 dataset reports")
    ap.add_argument("--species", required=True)
    ap.add_argument("--population", required=True, choices=sorted(POPULATIONS))
    ap.add_argument("--expect-commit", help="the published commit of the kit (40 hex characters); a real run needs it")
    ap.add_argument("--out", required=True)
    ap.add_argument("--rehearsal", action="store_true",
                    help="stand-in inputs of cohort genomes only, to show the kit runs; measures nothing")
    ap.add_argument("--tables", help="a rehearsal's stand-in candidates.tsv, pairs.tsv and datasets.tsv")
    args = ap.parse_args(argv)
    # The kit root, the pins and the frozen tables are the kit's own: no
    # option points them elsewhere, and --tables needs --rehearsal.
    try:
        if args.tables is not None and not args.rehearsal:
            refuse("--tables is for a rehearsal only; a real run reads the kit's own frozen tables")
        gate = frozen_gate(args.expect_commit, args.rehearsal)
        res = score(load_inputs(args, gate))
        res["kit"] = gate
        write_outputs(res, args.out)
    except Refusal as e:
        print(f"score.py: refused: {e}", file=sys.stderr)
        return 2
    print(f"score.py: {res['genomes_scored']} of {res['genomes_frozen']} genomes scored on "
          f"{len(res['datasets'])} datasets; outputs in {args.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
