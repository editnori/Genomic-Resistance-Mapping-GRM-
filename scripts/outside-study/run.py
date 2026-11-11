#!/usr/bin/env python3
"""The run script of the outside-study test, a retrospective test on
newly deposited studies from the same database (docs/outside-study-test.md).

Stages, each callable on its own and in this order:

  setup          derive the frozen ids and the eligible datasets from the
                 kit's frozen tables; build the scratch data root (links and
                 copies; the master root is only read)
  fetch          grm fetch fasta for the frozen ids; a database answer of no
                 sequence is the availability result, any other failure
                 stops the stage
  species-check  grm species-check against the master root's recorded
                 references; the recorded verdicts must come out unchanged
  count          grm count for the frozen genomes and the cohort genomes the
                 replay checks need
  known-gene     the known-gene comparator, after its replay on the cohort
  predict        grm predict per dataset report over the surviving genomes,
                 after the final models replay on the cohort

Each stage writes work/receipts/<stage>.json (counts, SHA-256 of inputs and
outputs, the kit commit, the grm binary's SHA-256, seconds) and refuses to
run when an earlier receipt is missing or a file it lists has changed, or
when its own receipt exists. Every start and end of a stage is appended to
work/receipts/attempts.jsonl. Every stage passes the frozen gate
(score.frozen_gate) first. Standard library only.
"""

import argparse
import concurrent.futures
import datetime
import importlib.util
import json
import shutil
import subprocess
import sys
import time
from pathlib import Path

# The kit root is where this file is; no option points it elsewhere.
KIT_ROOT = Path(__file__).resolve().parents[2]
KIT = KIT_ROOT / "scripts" / "outside-study"
# What `grm species-check` decides a fail and a boundary flag from. When a
# recorded genome keeps all of them, only the GC rule can move its verdict.
SAME_DISTANCES = ("distance", "nearest_species", "nearest_distance", "inside_species")
# The last part of `grm fetch fasta`'s failure message when BV-BRC answered
# and the answer held no sequence (src/bvbrc.rs, `fetch_fasta`). Any other
# failure (a timeout, a refused connection, an HTTP error status) is not an
# answer about the genome.
NO_SEQUENCE_ANSWER = "empty FASTA response"


def _load_score():
    spec = importlib.util.spec_from_file_location("outside_study_score", KIT / "score.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


S = _load_score()
Refusal, refuse = S.Refusal, S.refuse
STAGES = S.STAGES


def check_running():
    """This file and the score.py it loaded are the kit's own."""
    if Path(__file__).resolve() != KIT / "run.py":
        refuse(f"the running script {Path(__file__).resolve()} is not the kit's own {KIT / 'run.py'}")
    if Path(S.__file__).resolve() != KIT / "score.py":
        refuse(f"the loaded score.py {Path(S.__file__).resolve()} is not {KIT / 'score.py'}")
    S.check_running(__file__, S.RUN_REL)


def no_sequence_answer(message):
    m = message.strip()
    return m == NO_SEQUENCE_ANSWER or m.endswith(": " + NO_SEQUENCE_ANSWER)


def log(msg):
    print(f"run.py: {msg}", file=sys.stderr, flush=True)


def read_ids(path):
    ids = [line.strip() for line in Path(path).read_text().splitlines() if line.strip()]
    if len(ids) != len(set(ids)):
        refuse(f"{path}: an id is listed twice")
    for g in ids:
        if not S.GENOME_ID.match(g):
            refuse(f"{path}: {g!r} is not a BV-BRC genome id")
    return ids


def listing_digest(directory, pattern):
    """SHA-256 of `name<TAB>size` lines for the files a directory holds:
    names the master root's thousands of FASTA files without reading them."""
    d = Path(directory)
    if not d.is_dir():
        return "absent"
    lines = sorted(f"{p.name}\t{p.stat().st_size}" for p in d.glob(pattern))
    return S.sha256_bytes(("\n".join(lines) + "\n").encode())


class Run:
    def __init__(self, args):
        self.args = args
        self.grm = Path(args.grm).resolve()
        self.master = Path(args.master).resolve()
        self.work = Path(args.work).resolve()
        self.root = self.work / "root"
        self.species = args.species
        self.population = args.population
        self.out_slug = f"outside__{args.species}"
        self.replay_slug = f"replay__{args.species}"
        self.rehearsal = args.rehearsal
        self.jobs = max(1, args.jobs)
        if self.work == self.master or self.master in self.work.parents:
            refuse("the work directory is inside the master root, which this script never writes")

    # -------------------------------------------------------------- receipts

    def receipt_path(self, stage):
        return self.work / "receipts" / f"{stage}.json"

    def attempt(self, stage, event, message=""):
        """Appends one line to receipts/attempts.jsonl: every start and every
        end (done, refused or failed) of a stage."""
        p = self.work / "receipts" / "attempts.jsonl"
        p.parent.mkdir(parents=True, exist_ok=True)
        line = {"stage": stage, "event": event, "message": message,
                "utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")}
        with open(p, "a") as f:
            f.write(json.dumps(line, sort_keys=True) + "\n")

    def key(self, path):
        """How a receipt names a file: relative to the work directory when it
        is inside it, so a saved work directory can be scored anywhere."""
        p = Path(path)
        if p.is_absolute() and self.work in p.parents:
            return "work/" + str(p.relative_to(self.work))
        return str(p)

    def digest(self, key):
        if key.startswith("listing:"):
            _, directory, pattern = key.split(":", 2)
            return listing_digest(directory, pattern)
        p = self.work / key[len("work/"):] if key.startswith("work/") else Path(key)
        return S.sha256_file(p) if p.is_file() else "absent"

    def gate(self):
        """The frozen gate, before every stage, real or rehearsal."""
        check_running()
        return S.frozen_gate(self.args.expect_commit, self.rehearsal)

    def kit_state(self, gate):
        return dict(gate, **{
            "run.py": S.sha256_file(Path(__file__).resolve()),
            "score.py": S.sha256_file(Path(S.__file__).resolve()),
        })

    def check_earlier(self, stage, gate):
        """This stage has no receipt yet; every earlier receipt exists and
        every file it lists is unchanged; the same grm binary, kit commit
        and scripts; the same arguments."""
        if self.receipt_path(stage).exists():
            refuse(f"{stage}: this stage already has a receipt in {self.work}; a stage never runs twice "
                   f"in one work directory, so start a new work directory")
        grm_sha = S.sha256_file(self.grm)
        kit = self.kit_state(gate)
        for earlier in STAGES[: STAGES.index(stage)]:
            p = self.receipt_path(earlier)
            if not p.is_file():
                refuse(f"{stage}: the {earlier} stage has no receipt; run it first")
            r = json.loads(p.read_text())
            if r["grm"]["sha256"] != grm_sha:
                refuse(f"{stage}: the grm binary differs from the one the {earlier} stage ran")
            if r["kit"]["run.py"] != kit["run.py"] or r["kit"]["score.py"] != kit["score.py"]:
                refuse(f"{stage}: the kit scripts changed since the {earlier} stage")
            if r["kit"]["commit"] != kit["commit"]:
                refuse(f"{stage}: the kit commit changed since the {earlier} stage")
            if r["args"] != self.arg_record():
                refuse(f"{stage}: the arguments differ from the {earlier} stage's")
            for section in ("inputs", "outputs"):
                for key, sha in r[section].items():
                    if self.digest(key) != sha:
                        refuse(f"{stage}: {key} changed since the {earlier} stage")

    def arg_record(self):
        a = self.args
        opt = lambda v: None if v is None else str(Path(v).resolve())  # noqa: E731
        return {
            "grm": str(self.grm), "master": str(self.master), "work": str(self.work),
            "species": self.species, "population": self.population,
            "reports": str(Path(a.reports).resolve()), "expect_commit": a.expect_commit,
            "tables": opt(a.tables), "ids": opt(a.ids), "datasets": opt(a.datasets),
            "rehearsal": self.rehearsal,
        }

    def write_receipt(self, stage, started, inputs, outputs, counts, gate, extra=None):
        r = {
            "stage": stage,
            "rehearsal": self.rehearsal,
            "args": self.arg_record(),
            "jobs": self.jobs,
            "grm": {"path": str(self.grm), "sha256": S.sha256_file(self.grm)},
            "kit": self.kit_state(gate),
            "inputs": {self.key(k): self.digest(self.key(k)) for k in sorted(inputs)},
            "outputs": {self.key(k): self.digest(self.key(k)) for k in sorted(outputs)},
            "counts": counts,
            "seconds": round(time.monotonic() - started, 3),
        }
        r.update(extra or {})
        p = self.receipt_path(stage)
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(json.dumps(r, sort_keys=True, indent=1) + "\n")
        log(f"{stage}: done in {r['seconds']:.0f} s; receipt {p}")
        for k, v in counts.items():
            log(f"  {k}: {v}")

    # -------------------------------------------------------------- inputs

    def frozen_lists(self):
        """The frozen ids and the eligible datasets. A real run derives both
        from the kit's frozen tables, with score.py's definitions; only a
        rehearsal may give them as files."""
        a = self.args
        if a.ids is not None:
            ids, ds = read_ids(a.ids), read_ids_plain(a.datasets)
        else:
            tables = S.load_tables(self.rehearsal, a.tables)
            ids = list(S.population(tables["candidates"], self.species, self.population))
            ds = S.eligible(tables["datasets"], self.species, self.population)
            if not ds:
                refuse(f"no dataset of {self.species} is eligible in the {self.population} population")
        for d in ds:
            if not d.startswith(self.species + "__"):
                refuse(f"{d} is not a {self.species} dataset")
        return sorted(ids, key=S.genome_key), ds

    def ids(self):
        return read_ids(self.work / "ids.txt")

    def datasets(self):
        return read_ids_plain(self.work / "datasets.txt")

    def reports(self, datasets=None):
        pins = S.read_pins(KIT / "reports.sha256")
        out = {}
        for d in datasets or self.datasets():
            out[d] = S.load_report(self.args.reports, d, pins)
        ks = {r["counts"]["k"] for r, _ in out.values()}
        if len(ks) != 1:
            refuse(f"the reports use several k: {sorted(ks)}")
        return out

    def report_inputs(self, datasets=None):
        return [str(S.report_path(self.args.reports, d).resolve()) for d in datasets or self.datasets()]

    def cohort_ids(self, reports):
        ids = set()
        for r, _ in reports.values():
            ids.update(s["id"] for s in r["samples"])
        return ids

    def grm_run(self, *argv, log_name):
        cmd = [str(self.grm), *argv]
        log("$ " + " ".join(cmd if len(cmd) < 20 else cmd[:12] + [f"... ({len(cmd)} arguments)"]))
        logp = self.work / "logs" / log_name
        logp.parent.mkdir(parents=True, exist_ok=True)
        with open(logp, "w") as f:
            p = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT)
        tail = logp.read_text().splitlines()[-15:]
        for line in tail:
            print(f"  | {line}", file=sys.stderr)
        if p.returncode != 0:
            refuse(f"{' '.join(argv[:3])} exited {p.returncode}; see {logp}")
        return logp

    def fasta(self, g):
        return self.root / "fasta" / f"{g}.fna.gz"

    # -------------------------------------------------------------- stages

    def setup(self, gate):
        started = time.monotonic()
        if self.root.exists() and any(self.root.iterdir()):
            refuse(f"{self.root} exists and is not empty; start a new work directory")
        ids, datasets = self.frozen_lists()
        reports = self.reports(datasets)
        cohort = self.cohort_ids(reports)
        if self.rehearsal:
            stray = [g for g in ids if g not in cohort]
            if stray:
                refuse(f"rehearsal: {stray[0]} is not a cohort genome; a rehearsal runs on cohort genomes only")
        else:
            inside = [g for g in ids if g in cohort]
            if inside:
                refuse(f"{inside[0]} is a cohort genome of these reports; the frozen list holds none")
        mfasta = self.master / "fasta"
        held = [g for g in ids if (mfasta / f"{g}.fna.gz").exists()]
        if held:
            refuse(f"{held[0]} already has a FASTA file in the master root; fetch must find the frozen genomes' availability itself")
        for sub in ("fasta", "catalog", "curated", "species-check/sketches", "species"):
            (self.root / sub).mkdir(parents=True, exist_ok=True)
        linked = 0
        if mfasta.is_dir():
            for p in mfasta.glob("*.fna.gz"):
                (self.root / "fasta" / p.name).symlink_to(p)
                linked += 1
        sketches = 0
        msk = self.master / "species-check" / "sketches"
        if msk.is_dir():
            for p in msk.glob("*.json"):
                (self.root / "species-check" / "sketches" / p.name).symlink_to(p)
                sketches += 1
        ids_file, ds_file = self.work / "ids.txt", self.work / "datasets.txt"
        ids_file.write_text("".join(g + "\n" for g in ids))
        ds_file.write_text("".join(d + "\n" for d in datasets))
        inputs = [*self.report_inputs(datasets), f"listing:{mfasta}:*.fna.gz", f"listing:{msk}:*.json"]
        if self.args.ids is not None:
            inputs += [str(Path(self.args.ids).resolve()), str(Path(self.args.datasets).resolve())]
        catalog_files = 0
        for name in ("AMR_CDS.fa", "ReferenceGeneCatalog.txt", "version.txt"):
            src = self.master / "catalog" / name
            if src.is_file():
                shutil.copyfile(src, self.root / "catalog" / name)
                inputs.append(str(src))
                catalog_files += 1
        if catalog_files not in (0, 3):
            refuse("the master catalog folder is incomplete")
        if catalog_files == 0 and not self.rehearsal:
            refuse("the master root has no catalog; counting needs it for the known-gene comparator")
        for rel in ("curated/species.json", "species-check/species_check.json"):
            src = self.master / rel
            if not src.is_file():
                refuse(f"the master root has no {rel}")
            inputs.append(str(src))
        for p in sorted((self.master / "species").glob("*/genomes.txt")):
            inputs.append(str(p))
        out_list = self.root / "species" / self.out_slug / "genomes.txt"
        out_list.parent.mkdir(parents=True, exist_ok=True)
        out_list.write_text("\n".join(ids) + "\n")
        replay = sorted(cohort, key=S.genome_key)
        missing = [g for g in replay if not (mfasta / f"{g}.fna.gz").exists()]
        if missing and not self.rehearsal:
            refuse(f"{len(missing)} cohort genomes of these reports have no FASTA in the master root, e.g. {missing[0]}; the replay checks need every one")
        replay = [g for g in replay if g not in missing]
        rp = self.root / "species" / self.replay_slug / "genomes.txt"
        rp.parent.mkdir(parents=True, exist_ok=True)
        rp.write_text("\n".join(replay) + ("\n" if replay else ""))
        self.write_receipt("setup", started, inputs, [str(ids_file), str(ds_file), str(out_list), str(rp)], {
            "frozen ids requested": len(ids),
            "datasets": len(reports),
            "cohort genomes for the replay checks": len(replay),
            "cohort genomes without FASTA (rehearsal only)": len(missing),
            "master FASTA files linked": linked,
            "master sketches linked": sketches,
            "catalog files copied": catalog_files,
        }, gate, {"frozen_ids": ids, "datasets": datasets})

    def fetch(self, gate):
        started = time.monotonic()
        ids = self.ids()
        failures_log = self.root / "species" / self.out_slug / "fasta_failures.tsv"
        if failures_log.exists():
            failures_log.unlink()
        self.grm_run("--data", str(self.root), "fetch", "fasta", "--species", self.out_slug,
                     "--workers", str(self.jobs), log_name="fetch.log")
        failed = {}
        if failures_log.exists():
            for line in failures_log.read_text().splitlines():
                g, _, why = line.partition("\t")
                failed[g] = why
        rows, outputs, unanswered = [], [], []
        for g in ids:
            f = self.fasta(g)
            if f.is_file() and not f.is_symlink():
                rows.append({"genome_id": g, "fasta": "yes", "detail": ""})
                outputs.append(str(f))
            elif g in failed and no_sequence_answer(failed[g]):
                rows.append({"genome_id": g, "fasta": "no", "detail": " ".join(failed[g].split())})
            else:
                unanswered.append((g, " ".join(failed.get(g, "no FASTA file after fetch, and no failure recorded").split())))
        if unanswered:
            g, why = unanswered[0]
            refuse(f"fetch: {len(unanswered)} genomes have no FASTA and no answer of no sequence from the database "
                   f"(e.g. {g}: {why}); the stage stops without a receipt and may be run again in this work directory")
        out = self.work / "fetch" / "availability.tsv"
        out.parent.mkdir(parents=True, exist_ok=True)
        S.write_tsv(out, ["genome_id", "fasta", "detail"], rows)
        n_yes = sum(r["fasta"] == "yes" for r in rows)
        self.write_receipt("fetch", started, [], [str(out), *outputs], {
            "requested": len(ids),
            "fetched": n_yes,
            "no sequence (the database answered with none)": len(ids) - n_yes,
            "bytes fetched (gzip)": sum(Path(p).stat().st_size for p in outputs),
        }, gate)

    def species_check(self, gate):
        started = time.monotonic()
        ids = self.ids()
        avail = {r["genome_id"]: r for r in S.read_tsv(self.work / "fetch" / "availability.tsv")}
        fetched = [g for g in ids if avail[g]["fasta"] == "yes"]
        # Fresh copies every run: the recorded file and the cohort lists are
        # the master's, never an earlier run's.
        recorded_path = self.master / "species-check" / "species_check.json"
        dst = self.root / "species-check" / "species_check.json"
        shutil.copyfile(recorded_path, dst)
        shutil.copyfile(self.master / "curated" / "species.json", self.root / "curated" / "species.json")
        species_rows = json.loads((self.master / "curated" / "species.json").read_text())
        name = None
        for row in species_rows:
            src = self.master / "species" / row["slug"] / "genomes.txt"
            text = src.read_text() if src.is_file() else ""
            if row["slug"] == self.species:
                name = row["species"]
                text = text.rstrip("\n") + ("\n" if text.strip() else "") + "".join(g + "\n" for g in fetched)
            d = self.root / "species" / row["slug"]
            d.mkdir(parents=True, exist_ok=True)
            (d / "genomes.txt").write_text(text)
        if name is None:
            refuse(f"{self.species} is not in the master root's curated/species.json")
        recorded = json.loads(recorded_path.read_text())
        rec_refs = {s["species"]: s["reference"] for s in recorded["species"]}
        if name not in rec_refs:
            refuse(f"the recorded species check has no reference for {name}")
        self.grm_run("--data", str(self.root), "species-check", log_name="species-check.log")
        now = json.loads(dst.read_text())
        changes = []
        for s in now["species"]:
            if s["species"] in rec_refs and s["reference"] != rec_refs[s["species"]]:
                changes.append({"genome_id": s["reference"], "species": s["species"],
                                "recorded": f"reference {rec_refs[s['species']]}", "now": f"reference {s['reference']}"})
        now_by = {(g["species"], g["genome_id"]): g for g in now["genomes"]}
        # A recorded genome's verdict may move only between pass and
        # borderline, and only when its distances are the recorded ones: then
        # the GC rule alone moved it, because the species' GC median and
        # spread are now taken over the cohort and the new genomes together.
        # Such a move is listed and the run goes on. Any other change stops it.
        moves = []
        for g in recorded["genomes"]:
            n = now_by.get((g["species"], g["genome_id"]))
            if n is not None and n["verdict"] == g["verdict"]:
                continue
            row = {"genome_id": g["genome_id"], "species": g["species"],
                   "recorded": g["verdict"], "now": None if n is None else n["verdict"]}
            gc_only = (
                n is not None
                and {g["verdict"], n["verdict"]} == {"pass", "borderline"}
                and all(k in g and k in n and g[k] == n[k] for k in SAME_DISTANCES)
            )
            (moves if gc_only else changes).append(row)
        cdir = self.work / "species-check"
        cdir.mkdir(parents=True, exist_ok=True)
        for stale in ("cohort-changes.tsv", "cohort-gc-moves.tsv"):
            if (cdir / stale).exists():
                (cdir / stale).unlink()
        if changes:
            p = cdir / "cohort-changes.tsv"
            S.write_tsv(p, ["genome_id", "species", "recorded", "now"], changes)
            refuse(f"{len(changes)} recorded species-check results changed when the frozen genomes joined; see {p}")
        moves_path = cdir / "cohort-gc-moves.tsv"
        S.write_tsv(moves_path, ["genome_id", "species", "recorded", "now"], moves)
        rows = []
        for g in fetched:
            n = now_by.get((name, g))
            if n is None:
                refuse(f"the species check gave no verdict for fetched genome {g}")
            rows.append({"genome_id": g, "verdict": n["verdict"],
                         "reason": " ".join((n.get("reason") or "").split()),
                         "distance": n["distance"], "gc": n["gc"]})
        out = cdir / "verdicts.tsv"
        S.write_tsv(out, ["genome_id", "verdict", "reason", "distance", "gc"], rows)
        count = {v: sum(r["verdict"] == v for r in rows) for v in ("pass", "borderline", "fail")}
        self.write_receipt("species-check", started, [str(recorded_path)], [str(out), str(moves_path), str(dst)], {
            "fetched genomes judged": len(rows),
            "pass": count["pass"], "borderline (kept, flagged)": count["borderline"], "fail (excluded)": count["fail"],
            "recorded genomes with unchanged verdicts": len(recorded["genomes"]) - len(moves),
            "recorded genomes whose GC flag alone moved (listed)": len(moves),
            "reference": rec_refs[name],
        }, gate)

    def count(self, gate):
        started = time.monotonic()
        ids = self.ids()
        reports = self.reports()
        k = next(iter(reports.values()))[0]["counts"]["k"]
        argv = ["--data", str(self.root), "count", "--species", f"{self.out_slug},{self.replay_slug}", "--k", str(k)]
        if not (self.root / "catalog" / "version.txt").is_file():
            if not self.rehearsal:
                refuse("no catalog in the scratch root")
            argv.append("--allow-missing-catalog")
        self.grm_run(*argv, log_name="count.log")
        outputs, counts = [], {"k": k}
        for slug in (self.out_slug, self.replay_slug):
            cj = self.root / "species" / slug / "counted.json"
            c = json.loads(cj.read_text())
            outputs.append(str(cj))
            for g in c["genomes"]:
                outputs.append(str(self.root / "species" / slug / "ksets" / f"{g['genome_id']}.kset"))
                if not g.get("kset_digest"):
                    refuse(f"{slug}: counted.json records no k-mer set digest for {g['genome_id']}")
            counts[f"{slug}: counted"] = len(c["genomes"])
            counts[f"{slug}: missing FASTA"] = len(c["missing_fasta"])
            counts[f"{slug}: failed"] = len(c["failed"])
            if slug == self.replay_slug:
                if c["missing_fasta"] or c["failed"]:
                    refuse(f"{len(c['missing_fasta']) + len(c['failed'])} cohort genomes for the replay checks were not counted")
                cohort_rows = [{"genome_id": g["genome_id"], "kset_digest": g["kset_digest"]}
                               for g in sorted(c["genomes"], key=lambda g: S.genome_key(g["genome_id"]))]
            if slug == self.out_slug:
                got = {g["genome_id"]: g for g in c["genomes"]}
                why = dict((g, w) for g, w in c["failed"])
                rows = []
                for g in ids:
                    if g in got:
                        rows.append({"genome_id": g, "counted": "yes", "detail": "",
                                     "kset_digest": got[g]["kset_digest"]})
                    else:
                        rows.append({"genome_id": g, "counted": "no",
                                     "detail": " ".join(why.get(g, "no FASTA").split()), "kset_digest": ""})
        out = self.work / "count" / "counted.tsv"
        out.parent.mkdir(parents=True, exist_ok=True)
        S.write_tsv(out, ["genome_id", "counted", "detail", "kset_digest"], rows)
        cd = self.work / "count" / "cohort-digests.tsv"
        S.write_tsv(cd, ["genome_id", "kset_digest"], cohort_rows)
        outputs += [str(out), str(cd)]
        counts["kset bytes"] = sum(Path(p).stat().st_size for p in outputs if p.endswith(".kset"))
        self.write_receipt("count", started, [], outputs, counts, gate)

    def survivors(self):
        ids = self.ids()
        avail = {r["genome_id"]: r for r in S.read_tsv(self.work / "fetch" / "availability.tsv")}
        verdicts = {r["genome_id"]: r for r in S.read_tsv(self.work / "species-check" / "verdicts.tsv")}
        counted = {r["genome_id"]: r for r in S.read_tsv(self.work / "count" / "counted.tsv")}
        excluded, _ = S.exclusions(ids, avail, verdicts, counted)
        return [g for g in ids if g not in excluded], excluded

    def counted_genes(self, slug):
        c = json.loads((self.root / "species" / slug / "counted.json").read_text())
        return {g["genome_id"]: g for g in c["genomes"]}

    def known_gene(self, gate):
        started = time.monotonic()
        reports = self.reports()
        survivors, _ = self.survivors()
        kdir = self.work / "known-gene"
        if kdir.exists():
            shutil.rmtree(kdir)
        kdir.mkdir(parents=True)
        version_file = self.root / "catalog" / "version.txt"
        if not version_file.is_file():
            if not self.rehearsal:
                refuse("no catalog in the scratch root")
            out = kdir / "unavailable.txt"
            out.write_text("no AMRFinderPlus catalog in the stand-in master root, so grm count recorded no gene calls\n")
            self.write_receipt("known-gene", started, [], [str(out)], {"known-gene comparator": "not available (rehearsal without a catalog)"}, gate)
            return
        version = version_file.read_text().strip()
        cohort = self.counted_genes(self.replay_slug)
        new = self.counted_genes(self.out_slug)
        want_version = version + S.CALLS_VERSION_SUFFIX
        for g in list(cohort.values()) + list(new.values()):
            if g["catalog_version"] != want_version:
                refuse(f"{g['genome_id']}: gene calls from catalog {g['catalog_version']!r}, not {want_version!r}")
        hits = {g: c["genes"] for g, c in cohort.items()}
        ds_rows, call_rows, problems, counts = [], [], [], {}
        for d, (report, _) in sorted(reports.items()):
            cat = report.get("catalog")
            if cat is None:
                ds_rows.append({"dataset": d, "rule": "", "catalog_version": "", "replay": "no rule published", "samples": 0})
                counts[f"{d}"] = "no rule published"
                continue
            if cat["version"] != version:
                refuse(f"{d}: the report's catalog is {cat['version']}, the scratch root's {version}")
            p = S.known_gene_replay(report, hits, every_sample=not self.rehearsal)
            problems.extend(f"{d}: {x}" for x in p)
            ds_rows.append({"dataset": d, "rule": json.dumps(cat["rule"], separators=(",", ":")),
                            "catalog_version": version, "replay": "pass" if not p else "fail",
                            "samples": len(report["samples"])})
            counts[d] = f"replay {'pass' if not p else 'FAIL'} on {len(report['samples'])} samples"
            for g in survivors:
                met, used = S.explain(cat["rule"], new[g]["genes"])
                call_rows.append({"dataset": d, "genome_id": g,
                                  "groups": ";".join(h["group"] for h in new[g]["genes"]),
                                  "rule_groups": ";".join(used),
                                  "catalog_resistant": "yes" if met else "no"})
        if problems:
            p = kdir / "replay-problems.txt"
            p.write_text("\n".join(problems) + "\n")
            refuse(f"the known-gene replay failed on the cohort ({len(problems)} mismatches); see {p}")
        S.write_tsv(kdir / "datasets.tsv", ["dataset", "rule", "catalog_version", "replay", "samples"], ds_rows)
        S.write_tsv(kdir / "calls.tsv", ["dataset", "genome_id", "groups", "rule_groups", "catalog_resistant"], call_rows)
        counts["scored genomes with calls"] = len(survivors)
        self.write_receipt("known-gene", started, [str(version_file)],
                           [str(kdir / "datasets.tsv"), str(kdir / "calls.tsv")], counts, gate)

    def predict_jobs(self, jobs):
        def one(job):
            d, genomes, kdir, out = job
            argv = [str(self.grm), "predict", "--report", str(S.report_path(self.args.reports, d)), "--out", str(out)]
            for g in genomes:
                argv += ["--kset", str(kdir / f"{g}.kset")]
            with open(out.with_suffix(".stderr"), "w") as f:
                p = subprocess.run(argv, stdout=subprocess.DEVNULL, stderr=f)
            return d, out, p.returncode

        with concurrent.futures.ThreadPoolExecutor(max_workers=self.jobs) as pool:
            done = list(pool.map(one, jobs))
        for d, out, code in done:
            if code != 0:
                refuse(f"grm predict on {d} exited {code}; see {out.with_suffix('.stderr')}")

    def predict(self, gate):
        started = time.monotonic()
        reports = self.reports()
        survivors, _ = self.survivors()
        if not survivors:
            refuse("no genome survived the exclusions")
        pdir, rdir = self.work / "predict", self.work / "predict-replay"
        for d in (pdir, rdir):
            if d.exists():
                shutil.rmtree(d)
        rdir.mkdir(parents=True)
        okset = self.root / "species" / self.out_slug / "ksets"
        rkset = self.root / "species" / self.replay_slug / "ksets"
        cohort_counted = set(self.counted_genes(self.replay_slug))

        # First the cohort's own genomes, judged before any new genome is
        # predicted.
        replay_jobs = []
        for d, (report, _) in sorted(reports.items()):
            samples = [s["id"] for s in report["samples"] if s["id"] in cohort_counted]
            if samples:
                replay_jobs.append((d, samples, rkset, rdir / f"{d}.jsonl"))
        self.predict_jobs(replay_jobs)
        counts, problems, checked = {}, [], 0
        for d, (report, sha) in sorted(reports.items()):
            rp = rdir / f"{d}.jsonl"
            if rp.exists():
                p, n = S.replay_problems(report, sha, S.read_rows(rp))
                problems.extend(f"{d}: {x}" for x in p)
                checked += n
                counts[f"{d}: cohort genomes replayed"] = n
        if problems:
            p = rdir / "replay-problems.txt"
            p.write_text("\n".join(problems) + "\n")
            refuse(f"the final models do not replay on the cohort ({len(problems)} mismatches); see {p}")

        # Then the new genomes.
        pdir.mkdir(parents=True)
        (pdir / "survivors.txt").write_text("\n".join(survivors) + "\n")
        self.predict_jobs([(d, survivors, okset, pdir / f"{d}.jsonl") for d in sorted(reports)])
        for d in sorted(reports):
            rows = S.read_rows(pdir / f"{d}.jsonl")
            learners = sorted({r["learner"] for r in rows})
            counts[f"{d}: rows"] = f"{len(rows)} ({', '.join(learners)})"
        counts["genomes scored"] = len(survivors)
        counts["cohort genome-report replays"] = checked
        outputs = [str(pdir / "survivors.txt")] + [str(p) for p in sorted(pdir.glob("*.jsonl"))] \
            + [str(p) for p in sorted(rdir.glob("*.jsonl"))]
        self.write_receipt("predict", started, [], outputs, counts, gate)


def read_ids_plain(path):
    ds = [line.strip() for line in Path(path).read_text().splitlines() if line.strip()]
    if not ds or len(ds) != len(set(ds)):
        refuse(f"{path}: empty or a dataset listed twice")
    return ds


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("stage", choices=STAGES + ("all",))
    ap.add_argument("--grm", required=True, help="the grm binary, built at the kit's commit")
    ap.add_argument("--master", required=True, help="the master data root; only read")
    ap.add_argument("--work", required=True, help="a new work directory; the scratch root is <work>/root")
    ap.add_argument("--species", required=True, help="species slug, e.g. staphylococcus_aureus")
    ap.add_argument("--population", required=True, choices=sorted(S.POPULATIONS))
    ap.add_argument("--reports", required=True, help="the served revision-12 dataset reports")
    ap.add_argument("--expect-commit", help="the published commit of the kit (40 hex characters); a real run needs it")
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--rehearsal", action="store_true",
                    help="cohort genomes as a stand-in list, to show the kit runs; measures nothing")
    ap.add_argument("--tables", help="rehearsal only: stand-in candidates.tsv, pairs.tsv and datasets.tsv")
    ap.add_argument("--ids", help="rehearsal only: the stand-in genome ids, one per line")
    ap.add_argument("--datasets", help="rehearsal only: the dataset ids, one per line")
    args = ap.parse_args(argv)
    stages = STAGES if args.stage == "all" else (args.stage,)
    try:
        if not args.rehearsal:
            for opt in ("tables", "ids", "datasets"):
                if getattr(args, opt) is not None:
                    refuse(f"--{opt} is for a rehearsal only; a real run derives the lists from the kit's frozen tables")
        if (args.ids is None) != (args.datasets is None):
            refuse("--ids and --datasets go together")
        run = Run(args)
    except Refusal as e:
        log(f"refused: {e}")
        return 2
    for stage in stages:
        log(f"{stage}: start")
        run.attempt(stage, "start")
        try:
            gate = run.gate()
            run.check_earlier(stage, gate)
            getattr(run, stage.replace("-", "_"))(gate)
        except Refusal as e:
            run.attempt(stage, "refused", str(e))
            log(f"refused: {e}")
            return 2
        except Exception as e:
            run.attempt(stage, "failed", f"{type(e).__name__}: {e}")
            raise
        run.attempt(stage, "done")
    return 0


if __name__ == "__main__":
    sys.exit(main())
