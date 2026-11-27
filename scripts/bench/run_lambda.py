"""Lambda execution runner skeleton for Track A AMR tool benchmarks.

Features:
- Writes pinned `run-manifest.json` (software/DB/adapter/env-lock SHA-256s,
  split scheme, time_provenance, selection_role, offline flag).
- Enforces `--offline` execution (blocks proxy vars, sets offline env flags,
  requires local DB directory/file and optional expected DB SHA-256 verification).
- Atomic per-unit and run-level `DONE` / `FAILED` sentinel checkpointing with
  idempotent resume and config-drift refusal.
- Writes wall time, CPU time, and peak RSS strictly to separate `logs/timing.tsv`
  so `grm-bench-result-v1` JSON bytes stay 100% deterministic.
- Runs `score.check_zero_id_hygiene` on every result artifact before writing
  final run-level `DONE`.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable

import numpy as np

from scripts.bench.protocol import (
    ADAPTERS_DIR,
    ALL_TOOLS,
    AdapterSpec,
    BenchError,
    MODE_PRETRAINED,
    SELECTION_ROLE_INTERNAL_CV,
    canonical_json_bytes,
    load_adapter,
    require,
    sha256_bytes,
    sha256_file,
    validate_time_provenance,
)
from scripts.bench.score import check_zero_id_hygiene

try:
    import resource
except ImportError:  # pragma: no cover - non-POSIX fallback
    resource = None  # type: ignore[assignment]

RUN_MANIFEST_FORMAT = "grm-bench-run-manifest-v1"
TIMING_COLUMNS = ["unit_id", "tool", "mode", "status", "wall_s", "cpu_s", "peak_rss_kb"]
_GENOME_ID_SCRUB_RE = re.compile(r"(?<![\w.])\d+\.\d+(?!\d)")


def _peak_rss_kb() -> int:
    if resource is None:
        return 0
    return int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)


def _cpu_seconds() -> float:
    t = os.times()
    return float(t.user + t.system)


def scrub_genome_ids(text: str, known_genome_ids: set[str] | None = None) -> str:
    """Scrub known genome IDs from failure logs/sentinels before writing to disk."""
    if not known_genome_ids:
        return text
    return _GENOME_ID_SCRUB_RE.sub(
        lambda m: "[genome_id]" if m.group(0) in known_genome_ids else m.group(0),
        text,
    )


def hash_path_deterministic(path: Path | str) -> str:
    """Compute a deterministic SHA-256 over a file or directory tree (sorted relative paths)."""
    p = Path(path)
    require(p.exists(), f"path to hash does not exist: {p}")
    if p.is_file():
        return sha256_file(p)
    require(p.is_dir(), f"expected file or directory: {p}")
    h = hashlib.sha256()
    files = sorted(f for f in p.rglob("*") if f.is_file())
    require(len(files) > 0, f"database/asset directory is empty: {p}")
    for f in files:
        rel = f.relative_to(p).as_posix().encode("utf-8")
        digest = sha256_file(f).encode("ascii")
        h.update(rel + b"\0" + digest + b"\n")
    return h.hexdigest()


def build_offline_env(base_env: dict[str, str] | None = None, *, offline: bool = True) -> dict[str, str]:
    """Build subprocess environment with strict offline guards when `offline=True`."""
    env = dict(os.environ if base_env is None else base_env)
    if offline:
        for proxy_var in (
            "HTTP_PROXY",
            "HTTPS_PROXY",
            "ALL_PROXY",
            "http_proxy",
            "https_proxy",
            "all_proxy",
            "FTP_PROXY",
            "ftp_proxy",
        ):
            env.pop(proxy_var, None)
        env["GRM_BENCH_OFFLINE"] = "1"
        env["HF_HUB_OFFLINE"] = "1"
        env["CONDA_OFFLINE"] = "1"
        env["PIP_NO_INDEX"] = "1"
        env["NO_NETWORK"] = "1"
    return env


def verify_offline_db(
    adapter: AdapterSpec,
    db_path: Path | str | None,
    *,
    offline: bool,
    expected_db_sha256: str | None = None,
) -> str | None:
    """Verify local DB existence and optional SHA-256 pin before running any unit."""
    if adapter.mode == MODE_PRETRAINED:
        if offline:
            require(
                db_path is not None,
                f"offline pretrained run for {adapter.tool} requires --db-path pointing to local database",
            )
        if db_path is not None:
            actual = hash_path_deterministic(db_path)
            if expected_db_sha256 is not None:
                require(
                    actual == expected_db_sha256,
                    f"{adapter.tool}: database SHA-256 mismatch (got {actual[:12]}, expected {expected_db_sha256[:12]})",
                )
            return actual
        return None
    if db_path is not None:
        actual = hash_path_deterministic(db_path)
        if expected_db_sha256 is not None:
            require(actual == expected_db_sha256, f"{adapter.tool}: asset SHA-256 mismatch")
        return actual
    return None


def build_run_manifest(
    adapter: AdapterSpec,
    *,
    split_scheme: str,
    selection_role: str,
    time_provenance: dict[str, Any],
    offline: bool,
    db_sha256: str | None = None,
    env_lock_path: Path | str | None = None,
    unit_ids: list[str] | None = None,
) -> dict[str, Any]:
    """Create a deterministic, zero-path, pinned `run-manifest.json`."""
    prov = validate_time_provenance(
        split_scheme=split_scheme,
        declared_label=time_provenance["label"],
        fold_cutoff_date=time_provenance.get("fold_cutoff_date"),
        software_release_date=time_provenance.get("software_release_date", adapter.manifest["software_pin"].get("release_date")),
        database_release_date=time_provenance.get("database_release_date", (adapter.manifest["database_pin"] or {}).get("release_date")),
        weights_release_date=time_provenance.get("weights_release_date"),
        training_labels_max_date=time_provenance.get("training_labels_max_date"),
        mode=adapter.mode,
        is_bvbrc_insertion_time_split=bool(time_provenance.get("is_bvbrc_insertion_time_split", False)),
    )
    env_lock_sha = sha256_file(Path(env_lock_path)) if env_lock_path is not None else None
    return {
        "format": RUN_MANIFEST_FORMAT,
        "tool": adapter.tool,
        "mode": adapter.mode,
        "comparison_kind": adapter.comparison_kind,
        "primary_reading": adapter.primary_reading,
        "split_scheme": split_scheme,
        "selection_role": selection_role,
        "time_provenance": prov,
        "offline": bool(offline),
        "runtime_pins": {
            "python_version": ".".join(str(x) for x in sys.version_info[:3]),
            "numpy_version": str(np.__version__),
            "env_lock_sha256": env_lock_sha,
        },
        "software_pin": adapter.manifest["software_pin"],
        "database_pin": adapter.manifest["database_pin"],
        "database_sha256": db_sha256,
        "adapter_hashes": {
            "manifest_sha256": adapter.manifest_sha256,
            "drug_map_sha256": adapter.drug_map_sha256,
            "species_support_sha256": adapter.species_support_sha256,
        },
        "unit_count": len(unit_ids) if unit_ids is not None else 0,
    }


def _write_atomic(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_bytes(data)
    os.replace(tmp, path)


def _append_timing_row(timing_path: Path, row: dict[str, Any]) -> None:
    timing_path.parent.mkdir(parents=True, exist_ok=True)
    write_header = not timing_path.exists()
    with timing_path.open("a", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=TIMING_COLUMNS, delimiter="\t", lineterminator="\n")
        if write_header:
            writer.writeheader()
        writer.writerow({k: row[k] for k in TIMING_COLUMNS})


@dataclass(frozen=True)
class WorkUnit:
    unit_id: str
    payload: Any


class LambdaRunner:
    """Checkpoint-safe, offline-capable benchmark runner for Lambda boxes."""

    def __init__(
        self,
        adapter: AdapterSpec,
        out_dir: Path | str,
        *,
        split_scheme: str,
        time_provenance: dict[str, Any],
        selection_role: str = SELECTION_ROLE_INTERNAL_CV,
        offline: bool = True,
        db_path: Path | str | None = None,
        expected_db_sha256: str | None = None,
        env_lock_path: Path | str | None = None,
        known_genome_ids: set[str] | None = None,
    ) -> None:
        self.adapter = adapter
        self.out_dir = Path(out_dir)
        self.checkpoints_dir = self.out_dir / "checkpoints"
        self.results_dir = self.out_dir / "results"
        self.logs_dir = self.out_dir / "logs"
        self.timing_path = self.logs_dir / "timing.tsv"
        self.manifest_path = self.out_dir / "run-manifest.json"
        self.done_path = self.out_dir / "DONE"
        self.failed_path = self.out_dir / "FAILED"

        self.split_scheme = split_scheme
        self.time_provenance = time_provenance
        self.selection_role = selection_role
        self.offline = offline
        self.db_path = Path(db_path) if db_path is not None else None
        self.expected_db_sha256 = expected_db_sha256
        self.env_lock_path = Path(env_lock_path) if env_lock_path is not None else None
        self.known_genome_ids = known_genome_ids
        self.env = build_offline_env(offline=offline)

    def _unit_done(self, unit_id: str) -> Path:
        return self.checkpoints_dir / f"{unit_id}.DONE"

    def _unit_failed(self, unit_id: str) -> Path:
        return self.checkpoints_dir / f"{unit_id}.FAILED"

    def _write_failed(self, step: str, detail: str) -> None:
        safe = scrub_genome_ids(detail, self.known_genome_ids)
        body = f"step={step}\nerror={safe}\n"
        _write_atomic(self.failed_path, body.encode("utf-8"))

    def run(
        self,
        units: list[WorkUnit],
        executor: Callable[[WorkUnit, dict[str, str], Path], Path | None],
        *,
        result_paths_to_verify: list[Path] | None = None,
    ) -> dict[str, Any]:
        """Execute all work units with checkpointing, resume, timing, and zero-ID gate.

        `executor(unit, env, results_dir)` runs one unit and optionally returns a
        result JSON path (which will be checked by the zero-ID hygiene gate).
        """
        sorted_units = sorted(units, key=lambda u: u.unit_id)
        require(len({u.unit_id for u in sorted_units}) == len(sorted_units), "duplicate unit_id")
        for u in sorted_units:
            require(
                bool(re.fullmatch(r"[A-Za-z0-9._-]+", u.unit_id)),
                f"invalid unit_id for checkpoint filename: {u.unit_id!r}",
            )

        self.out_dir.mkdir(parents=True, exist_ok=True)
        self.checkpoints_dir.mkdir(parents=True, exist_ok=True)
        self.results_dir.mkdir(parents=True, exist_ok=True)
        self.logs_dir.mkdir(parents=True, exist_ok=True)

        # Verify offline DB / SHA-256 before writing manifest or running anything
        try:
            db_sha = verify_offline_db(
                self.adapter,
                self.db_path,
                offline=self.offline,
                expected_db_sha256=self.expected_db_sha256,
            )
            manifest = build_run_manifest(
                self.adapter,
                split_scheme=self.split_scheme,
                selection_role=self.selection_role,
                time_provenance=self.time_provenance,
                offline=self.offline,
                db_sha256=db_sha,
                env_lock_path=self.env_lock_path,
                unit_ids=[u.unit_id for u in sorted_units],
            )
        except Exception as exc:  # noqa: BLE001
            self._write_failed("preflight", str(exc))
            raise

        manifest_bytes = canonical_json_bytes(manifest)
        manifest_sha = sha256_bytes(manifest_bytes)
        if self.manifest_path.exists():
            prev_bytes = self.manifest_path.read_bytes()
            if prev_bytes != manifest_bytes:
                msg = "refused: existing run-manifest.json differs; cannot resume with changed config or pins"
                self._write_failed("manifest_resume_check", msg)
                raise BenchError(msg)
        else:
            _write_atomic(self.manifest_path, manifest_bytes)

        # Clear stale run-level FAILED sentinel before executing
        if self.failed_path.exists():
            self.failed_path.unlink()

        executed = 0
        skipped = 0
        produced_results: list[Path] = list(result_paths_to_verify or [])

        for unit in sorted_units:
            u_done = self._unit_done(unit.unit_id)
            u_fail = self._unit_failed(unit.unit_id)
            if u_done.is_file():
                sentinel_text = u_done.read_text(encoding="utf-8")
                if f"manifest_sha256={manifest_sha}" in sentinel_text:
                    skipped += 1
                    for line in sentinel_text.splitlines():
                        if line.startswith("result_rel="):
                            rel = line.split("=", 1)[1].strip()
                            if rel:
                                produced_results.append(self.out_dir / rel)
                    continue

            if u_fail.exists():
                u_fail.unlink()

            t0 = time.monotonic()
            c0 = _cpu_seconds()
            try:
                res_path = executor(unit, dict(self.env), self.results_dir)
            except Exception as exc:  # noqa: BLE001
                wall = round(time.monotonic() - t0, 6)
                cpu = round(_cpu_seconds() - c0, 6)
                _append_timing_row(
                    self.timing_path,
                    {
                        "unit_id": unit.unit_id,
                        "tool": self.adapter.tool,
                        "mode": self.adapter.mode,
                        "status": "FAILED",
                        "wall_s": wall,
                        "cpu_s": cpu,
                        "peak_rss_kb": _peak_rss_kb(),
                    },
                )
                safe_err = scrub_genome_ids(f"{type(exc).__name__}: {exc}", self.known_genome_ids)
                _write_atomic(u_fail, f"unit={unit.unit_id}\nerror={safe_err}\n".encode("utf-8"))
                self._write_failed(f"unit:{unit.unit_id}", safe_err)
                raise BenchError(f"unit {unit.unit_id} failed: {safe_err}") from exc

            wall = round(time.monotonic() - t0, 6)
            cpu = round(_cpu_seconds() - c0, 6)
            _append_timing_row(
                self.timing_path,
                {
                    "unit_id": unit.unit_id,
                    "tool": self.adapter.tool,
                    "mode": self.adapter.mode,
                    "status": "DONE",
                    "wall_s": wall,
                    "cpu_s": cpu,
                    "peak_rss_kb": _peak_rss_kb(),
                },
            )
            rel_line = ""
            if res_path is not None:
                rp = Path(res_path).resolve()
                produced_results.append(rp)
                try:
                    rel_line = f"result_rel={rp.relative_to(self.out_dir.resolve()).as_posix()}\n"
                except ValueError:
                    rel_line = ""
            _write_atomic(
                u_done,
                f"unit={unit.unit_id}\nmanifest_sha256={manifest_sha}\n{rel_line}".encode("utf-8"),
            )
            executed += 1

        # Mandatory zero-ID hygiene scan over all produced result JSONs before writing DONE
        unique_results = sorted({Path(p).resolve() for p in produced_results if Path(p).is_file()})
        for rp in unique_results:
            try:
                check_zero_id_hygiene(rp, known_genome_ids=self.known_genome_ids)
            except Exception as exc:  # noqa: BLE001
                self._write_failed(f"hygiene:{rp.name}", str(exc))
                raise

        # Write run-level DONE sentinel with SHA-256 sums of manifest + result files
        sum_lines = [f"{manifest_sha}  run-manifest.json"]
        for rp in unique_results:
            try:
                rel = rp.relative_to(self.out_dir.resolve()).as_posix()
            except ValueError:
                rel = rp.name
            sum_lines.append(f"{sha256_file(rp)}  {rel}")
        _write_atomic(self.done_path, ("\n".join(sum_lines) + "\n").encode("utf-8"))

        return {
            "status": "DONE",
            "tool": self.adapter.tool,
            "manifest_sha256": manifest_sha,
            "units_total": len(sorted_units),
            "units_executed": executed,
            "units_skipped": skipped,
            "results_verified": len(unique_results),
            "done_path": str(self.done_path),
            "timing_path": str(self.timing_path),
        }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="run_lambda.py", description="Track A Lambda benchmark runner (preflight/manifest CLI).")
    ap.add_argument("--tool", required=True, choices=ALL_TOOLS)
    ap.add_argument("--out", required=True, help="Output directory for run-manifest.json, checkpoints, results, logs.")
    ap.add_argument("--split-scheme", required=True, help="Split scheme (lineage, random, study-whole, time, external:<name>).")
    ap.add_argument("--time-provenance", required=True, help="Provenance label (retrospective_modern_db, historical_prospective, etc.).")
    ap.add_argument("--fold-cutoff-date", default=None)
    ap.add_argument("--software-release-date", default=None)
    ap.add_argument("--database-release-date", default=None)
    ap.add_argument("--training-labels-max-date", default=None)
    ap.add_argument("--bvbrc-time-split", action="store_true")
    ap.add_argument("--selection-role", default=SELECTION_ROLE_INTERNAL_CV)
    ap.add_argument("--offline", action="store_true", help="Enforce offline execution and local DB verification.")
    ap.add_argument("--db-path", default=None)
    ap.add_argument("--expected-db-sha256", default=None)
    ap.add_argument("--env-lock", default=None)
    ap.add_argument("--dry-run", action="store_true", help="Validate adapter, offline DB, and write run-manifest.json without running units.")
    args = ap.parse_args(argv)

    try:
        adapter = load_adapter(args.tool, adapters_dir=ADAPTERS_DIR)
        prov = {
            "label": args.time_provenance,
            "fold_cutoff_date": args.fold_cutoff_date,
            "software_release_date": args.software_release_date,
            "database_release_date": args.database_release_date,
            "training_labels_max_date": args.training_labels_max_date,
            "is_bvbrc_insertion_time_split": args.bvbrc_time_split,
        }
        runner = LambdaRunner(
            adapter,
            args.out,
            split_scheme=args.split_scheme,
            time_provenance=prov,
            selection_role=args.selection_role,
            offline=args.offline,
            db_path=args.db_path,
            expected_db_sha256=args.expected_db_sha256,
            env_lock_path=args.env_lock,
        )
        if args.dry_run:
            summary = runner.run([], executor=lambda u, env, rdir: None)
            print(json.dumps(summary, sort_keys=True))
            return 0
        print("run_lambda.py: pass --dry-run for preflight or invoke LambdaRunner from a task driver", file=sys.stderr)
        return 2
    except BenchError as exc:
        print(f"run_lambda.py: refused: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
