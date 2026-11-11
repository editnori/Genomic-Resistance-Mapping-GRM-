#!/usr/bin/env python3
"""Box-side run script for the outside-study test, a retrospective test on
newly deposited studies from the same database (docs/outside-study-test.md).

Stages, each callable on its own and in this order:

  setup          build the scratch data root (links and copies; the master
                 root is only read)
  fetch          grm fetch fasta for the frozen ids; the failures are the
                 availability result
  species-check  grm species-check against the master root's recorded
                 references; the recorded verdicts must come out unchanged
  count          grm count for the frozen genomes and the cohort genomes the
                 replay checks need
  known-gene     the known-gene comparator, after its replay on the cohort
  predict        grm predict per dataset report over the surviving genomes,
                 after the final models replay on the cohort

Each stage writes work/receipts/<stage>.json (counts, SHA-256 of inputs and
outputs, the kit commit, the grm binary's SHA-256, seconds) and refuses to
run when an earlier receipt is missing or a file it lists has changed.
known-gene and predict refuse unless the plan and the scripts match the
committed FROZEN.sha256. Standard library only.
"""

import argparse
import concurrent.futures
import importlib.util
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

KIT = Path(__file__).resolve().parent
STAGES = ("setup", "fetch", "species-check", "count", "known-gene", "predict")
GATED = ("known-gene", "predict")


def _load_score():
    spec = importlib.util.spec_from_file_location("outside_study_score", KIT / "score.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


S = _load_score()
Refusal, refuse = S.Refusal, S.refuse


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
        self.out_slug = f"outside__{args.species}"
        self.replay_slug = f"replay__{args.species}"
        self.kit_root = Path(args.kit_root).resolve()
        self.rehearsal = args.rehearsal
        self.jobs = max(1, args.jobs)
        if self.work == self.master or self.master in self.work.parents:
            refuse("the work directory is inside the master root, which this script never writes")

    # -------------------------------------------------------------- receipts

    def receipt_path(self, stage):
        return self.work / "receipts" / f"{stage}.json"

    def digest(self, key):
        if key.startswith("listing:"):
            _, directory, pattern = key.split(":", 2)
            return listing_digest(directory, pattern)
        p = Path(key)
        return S.sha256_file(p) if p.is_file() else "absent"

    def kit_state(self):
        head = subprocess.run(["git", "-C", str(self.kit_root), "rev-parse", "HEAD"],
                              capture_output=True, text=True)
        return {
            "commit": head.stdout.strip() if head.returncode == 0 else None,
            "run.py": S.sha256_file(KIT / "run.py"),
            "score.py": S.sha256_file(KIT / "score.py"),
        }

    def check_earlier(self, stage):
        """Every earlier receipt exists and every file it lists is unchanged;
        the same grm binary and kit scripts; the same arguments."""
        grm_sha = S.sha256_file(self.grm)
        kit = self.kit_state()
        for earlier in STAGES[: STAGES.index(stage)]:
            p = self.receipt_path(earlier)
            if not p.is_file():
                refuse(f"{stage}: the {earlier} stage has no receipt; run it first")
            r = json.loads(p.read_text())
            if r["grm"]["sha256"] != grm_sha:
                refuse(f"{stage}: the grm binary differs from the one the {earlier} stage ran")
            if r["kit"]["run.py"] != kit["run.py"] or r["kit"]["score.py"] != kit["score.py"]:
                refuse(f"{stage}: the kit scripts changed since the {earlier} stage")
            if r["args"] != self.arg_record():
                refuse(f"{stage}: the arguments differ from the {earlier} stage's")
            for section in ("inputs", "outputs"):
                for key, sha in r[section].items():
                    if self.digest(key) != sha:
                        refuse(f"{stage}: {key} changed since the {earlier} stage")

    def arg_record(self):
        a = self.args
        return {
            "grm": str(self.grm), "master": str(self.master), "work": str(self.work),
            "species": self.species, "ids": str(Path(a.ids).resolve()),
            "reports": str(Path(a.reports).resolve()), "datasets": str(Path(a.datasets).resolve()),
            "rehearsal": self.rehearsal,
        }

    def write_receipt(self, stage, started, inputs, outputs, counts, frozen=None):
        r = {
            "stage": stage,
            "rehearsal": self.rehearsal,
            "args": self.arg_record(),
            "jobs": self.jobs,
            "grm": {"path": str(self.grm), "sha256": S.sha256_file(self.grm)},
            "kit": self.kit_state(),
            "frozen": frozen,
            "inputs": {k: self.digest(k) for k in sorted(inputs)},
            "outputs": {k: self.digest(k) for k in sorted(outputs)},
            "counts": counts,
            "seconds": round(time.monotonic() - started, 3),
        }
        p = self.receipt_path(stage)
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(json.dumps(r, sort_keys=True, indent=1) + "\n")
        log(f"{stage}: done in {r['seconds']:.0f} s; receipt {p}")
        for k, v in counts.items():
            log(f"  {k}: {v}")

    # -------------------------------------------------------------- inputs

    def ids(self):
        return read_ids(self.args.ids)

    def datasets(self):
        ds = read_ids_plain(self.args.datasets)
        for d in ds:
            if not d.startswith(self.species + "__"):
                refuse(f"{self.args.datasets}: {d} is not a {self.species} dataset")
        return ds

    def reports(self):
        pins = S.read_pins(self.args.pins)
        out = {}
        for d in self.datasets():
            out[d] = S.load_report(self.args.reports, d, pins)
        ks = {r["counts"]["k"] for r, _ in out.values()}
        if len(ks) != 1:
            refuse(f"the reports use several k: {sorted(ks)}")
        return out

    def report_inputs(self):
        return [str(S.report_path(self.args.reports, d).resolve()) for d in self.datasets()]

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

    def setup(self):
        started = time.monotonic()
        ids = self.ids()
        reports = self.reports()
        cohort = self.cohort_ids(reports)
        if self.rehearsal:
            stray = [g for g in ids if g not in cohort]
            if stray:
                refuse(f"rehearsal: {stray[0]} is not a cohort genome; a rehearsal runs on cohort genomes only")
        else:
            inside = [g for g in ids if g in cohort]
            if inside:
                refuse(f"{inside[0]} is a cohort genome of these reports; the ids file must hold the frozen list")
        mfasta = self.master / "fasta"
        held = [g for g in ids if (mfasta / f"{g}.fna.gz").exists()]
        if held:
            refuse(f"{held[0]} already has a FASTA file in the master root; fetch must find the frozen genomes' availability itself")
        if self.root.exists() and any(self.root.iterdir()):
            log(f"setup: reusing {self.root}; links and copies are refreshed")
        for sub in ("fasta", "catalog", "curated", "species-check/sketches", "species"):
            (self.root / sub).mkdir(parents=True, exist_ok=True)
        linked = 0
        if mfasta.is_dir():
            for p in mfasta.glob("*.fna.gz"):
                dst = self.root / "fasta" / p.name
                if not dst.is_symlink() and not dst.exists():
                    dst.symlink_to(p)
                    linked += 1
        sketches = 0
        msk = self.master / "species-check" / "sketches"
        if msk.is_dir():
            for p in msk.glob("*.json"):
                dst = self.root / "species-check" / "sketches" / p.name
                if not dst.is_symlink() and not dst.exists():
                    dst.symlink_to(p)
                    sketches += 1
        inputs = [str(Path(self.args.ids).resolve()), str(Path(self.args.datasets).resolve()),
                  *self.report_inputs(),
                  f"listing:{mfasta}:*.fna.gz", f"listing:{msk}:*.json"]
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
        self.write_receipt("setup", started, inputs, [str(out_list), str(rp)], {
            "frozen ids requested": len(ids),
            "datasets": len(reports),
            "cohort genomes for the replay checks": len(replay),
            "cohort genomes without FASTA (rehearsal only)": len(missing),
            "master FASTA files linked now": linked,
            "master sketches linked now": sketches,
            "catalog files copied": catalog_files,
        })

    def fetch(self):
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
        rows, outputs = [], []
        for g in ids:
            f = self.fasta(g)
            if f.is_file() and not f.is_symlink():
                rows.append({"genome_id": g, "fasta": "yes", "detail": ""})
                outputs.append(str(f))
            else:
                why = failed.get(g, "no FASTA file after fetch, and no failure recorded")
                rows.append({"genome_id": g, "fasta": "no", "detail": " ".join(why.split())})
        out = self.work / "fetch" / "availability.tsv"
        out.parent.mkdir(parents=True, exist_ok=True)
        S.write_tsv(out, ["genome_id", "fasta", "detail"], rows)
        n_yes = sum(r["fasta"] == "yes" for r in rows)
        self.write_receipt("fetch", started, [], [str(out), *outputs], {
            "requested": len(ids),
            "fetched": n_yes,
            "not fetched": len(ids) - n_yes,
            "bytes fetched (gzip)": sum(Path(p).stat().st_size for p in outputs),
        })

    def species_check(self):
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
        for g in recorded["genomes"]:
            n = now_by.get((g["species"], g["genome_id"]))
            if n is None or n["verdict"] != g["verdict"]:
                changes.append({"genome_id": g["genome_id"], "species": g["species"],
                                "recorded": g["verdict"], "now": None if n is None else n["verdict"]})
        if changes:
            p = self.work / "species-check" / "cohort-changes.tsv"
            p.parent.mkdir(parents=True, exist_ok=True)
            S.write_tsv(p, ["genome_id", "species", "recorded", "now"], changes)
            refuse(f"{len(changes)} recorded species-check results changed when the frozen genomes joined; see {p}")
        rows = []
        for g in fetched:
            n = now_by.get((name, g))
            if n is None:
                refuse(f"the species check gave no verdict for fetched genome {g}")
            rows.append({"genome_id": g, "verdict": n["verdict"],
                         "reason": " ".join((n.get("reason") or "").split()),
                         "distance": n["distance"], "gc": n["gc"]})
        out = self.work / "species-check" / "verdicts.tsv"
        out.parent.mkdir(parents=True, exist_ok=True)
        S.write_tsv(out, ["genome_id", "verdict", "reason", "distance", "gc"], rows)
        count = {v: sum(r["verdict"] == v for r in rows) for v in ("pass", "borderline", "fail")}
        self.write_receipt("species-check", started, [str(recorded_path)], [str(out), str(dst)], {
            "fetched genomes judged": len(rows),
            "pass": count["pass"], "borderline (kept, flagged)": count["borderline"], "fail (excluded)": count["fail"],
            "recorded genomes with unchanged verdicts": len(recorded["genomes"]),
            "reference": rec_refs[name],
        })

    def count(self):
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
            counts[f"{slug}: counted"] = len(c["genomes"])
            counts[f"{slug}: missing FASTA"] = len(c["missing_fasta"])
            counts[f"{slug}: failed"] = len(c["failed"])
            if slug == self.replay_slug and (c["missing_fasta"] or c["failed"]):
                refuse(f"{len(c['missing_fasta']) + len(c['failed'])} cohort genomes for the replay checks were not counted")
            if slug == self.out_slug:
                got = {g["genome_id"] for g in c["genomes"]}
                why = dict((g, w) for g, w in c["failed"])
                rows = []
                for g in ids:
                    if g in got:
                        rows.append({"genome_id": g, "counted": "yes", "detail": ""})
                    else:
                        rows.append({"genome_id": g, "counted": "no",
                                     "detail": " ".join(why.get(g, "no FASTA").split())})
        out = self.work / "count" / "counted.tsv"
        out.parent.mkdir(parents=True, exist_ok=True)
        S.write_tsv(out, ["genome_id", "counted", "detail"], rows)
        outputs.append(str(out))
        counts["kset bytes"] = sum(Path(p).stat().st_size for p in outputs if p.endswith(".kset"))
        self.write_receipt("count", started, [], outputs, counts)

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

    def known_gene(self):
        started = time.monotonic()
        frozen = S.frozen_gate(self.kit_root)
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
            self.write_receipt("known-gene", started, [], [str(out)], {"known-gene comparator": "not available (rehearsal without a catalog)"}, frozen)
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
                           [str(kdir / "datasets.tsv"), str(kdir / "calls.tsv")], counts, frozen)

    def predict(self):
        started = time.monotonic()
        frozen = S.frozen_gate(self.kit_root)
        reports = self.reports()
        survivors, _ = self.survivors()
        if not survivors:
            refuse("no genome survived the exclusions")
        pdir, rdir = self.work / "predict", self.work / "predict-replay"
        for d in (pdir, rdir):
            if d.exists():
                shutil.rmtree(d)
            d.mkdir(parents=True)
        (pdir / "survivors.txt").write_text("\n".join(survivors) + "\n")
        okset = self.root / "species" / self.out_slug / "ksets"
        rkset = self.root / "species" / self.replay_slug / "ksets"
        cohort_counted = set(self.counted_genes(self.replay_slug))

        def one(job):
            d, genomes, kdir, out = job
            argv = [str(self.grm), "predict", "--report", str(S.report_path(self.args.reports, d)), "--out", str(out)]
            for g in genomes:
                argv += ["--kset", str(kdir / f"{g}.kset")]
            with open(out.with_suffix(".stderr"), "w") as f:
                p = subprocess.run(argv, stdout=subprocess.DEVNULL, stderr=f)
            return d, out, p.returncode

        jobs = []
        for d, (report, _) in sorted(reports.items()):
            jobs.append((d, survivors, okset, pdir / f"{d}.jsonl"))
            samples = [s["id"] for s in report["samples"] if s["id"] in cohort_counted]
            if samples:
                jobs.append((d, samples, rkset, rdir / f"{d}.jsonl"))
        with concurrent.futures.ThreadPoolExecutor(max_workers=self.jobs) as pool:
            done = list(pool.map(one, jobs))
        for d, out, code in done:
            if code != 0:
                refuse(f"grm predict on {d} exited {code}; see {out.with_suffix('.stderr')}")
        counts, problems, checked = {}, [], 0
        for d, (report, sha) in sorted(reports.items()):
            rp = rdir / f"{d}.jsonl"
            if rp.exists():
                p, n = S.replay_problems(report, sha, S.read_rows(rp))
                problems.extend(f"{d}: {x}" for x in p)
                checked += n
                counts[f"{d}: cohort genomes replayed"] = n
            rows = S.read_rows(pdir / f"{d}.jsonl")
            learners = sorted({r["learner"] for r in rows})
            counts[f"{d}: rows"] = f"{len(rows)} ({', '.join(learners)})"
        if problems:
            p = rdir / "replay-problems.txt"
            p.write_text("\n".join(problems) + "\n")
            refuse(f"the final models do not replay on the cohort ({len(problems)} mismatches); see {p}")
        counts["genomes scored"] = len(survivors)
        counts["cohort genome-report replays"] = checked
        outputs = [str(pdir / "survivors.txt")] + [str(p) for p in sorted(pdir.glob("*.jsonl"))] \
            + [str(p) for p in sorted(rdir.glob("*.jsonl"))]
        self.write_receipt("predict", started, [], outputs, counts, frozen)


def read_ids_plain(path):
    ds = [line.strip() for line in Path(path).read_text().splitlines() if line.strip()]
    if not ds or len(ds) != len(set(ds)):
        refuse(f"{path}: empty or a dataset listed twice")
    return ds


def main(argv=None, overrides=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("stage", choices=STAGES + ("all",))
    ap.add_argument("--grm", required=True, help="the grm binary, built at the kit's commit")
    ap.add_argument("--master", required=True, help="the master data root; only read")
    ap.add_argument("--work", required=True, help="the work directory; the scratch root is <work>/root")
    ap.add_argument("--species", required=True, help="species slug, e.g. staphylococcus_aureus")
    ap.add_argument("--ids", required=True, help="the frozen genome ids, one per line")
    ap.add_argument("--reports", required=True, help="the served revision-12 dataset reports")
    ap.add_argument("--datasets", required=True, help="the dataset ids to score, one per line")
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--rehearsal", action="store_true",
                    help="cohort genomes as a stand-in list, to show the kit runs; measures nothing")
    args = ap.parse_args(argv)
    # The report pins and the kit root are fixed; only the tests set other
    # values, through `overrides`.
    args.pins, args.kit_root = KIT / "reports.sha256", S.KIT_ROOT
    for k, v in (overrides or {}).items():
        setattr(args, k, v)
    stages = STAGES if args.stage == "all" else (args.stage,)
    try:
        run = Run(args)
        for stage in stages:
            log(f"{stage}: start")
            run.check_earlier(stage)
            if stage in GATED:
                S.frozen_gate(run.kit_root)
            getattr(run, stage.replace("-", "_"))()
    except Refusal as e:
        log(f"refused: {e}")
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
