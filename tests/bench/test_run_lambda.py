"""Synthetic unit tests for Lambda execution runner (`scripts/bench/run_lambda.py`):
- Pinned run-manifest.json (DB/env/adapter SHA-256s)
- --offline mode guard & DB SHA-256 verification
- Atomic DONE/FAILED sentinel checkpointing, idempotent resume, config-drift refusal
- Separate logs/timing.tsv
- Mandatory zero-ID hygiene gate before run-level DONE
"""

from __future__ import annotations

import csv
import json
import tempfile
import unittest
from pathlib import Path

from scripts.bench.pretrained import ProjectedEndpoint
from scripts.bench.protocol import (
    BenchError,
    TIER_SUPPORTED_PRIMARY,
    TIME_PROVENANCE_RETRO_MODERN_DB,
    TIME_PROVENANCE_RETRO_SAME_DB,
    load_adapter,
)
from scripts.bench.run_lambda import (
    LambdaRunner,
    WorkUnit,
    build_offline_env,
    hash_path_deterministic,
    main as lambda_main,
)
from scripts.bench.score import score_pretrained_cohort, write_result_atomic


def _clean_result_doc(adapter, known_ids: set[str]):
    gids = sorted(known_ids)
    n = len(gids)
    truth = tuple(i < n // 2 for i in range(n))
    calls = tuple("R" if y else "S" for y in truth)
    ep = ProjectedEndpoint(
        dataset="mycobacterium_tuberculosis__rifampin",
        species_slug="mycobacterium_tuberculosis",
        canonical_drug="rifampin",
        split_scheme="lineage",
        selection_role="internal_cv",
        support_tier=TIER_SUPPORTED_PRIMARY,
        unsupported_reason="-",
        reading=adapter.primary_reading,
        genome_ids=tuple(gids),
        truth=truth,
        lineages=tuple(i % 10 for i in range(n)),
        calls=calls,
        fold_of=tuple(0 for _ in range(n)),
    )
    return score_pretrained_cohort(
        adapter,
        [ep],
        split_scheme="lineage",
        time_provenance={"label": TIME_PROVENANCE_RETRO_SAME_DB},
        known_genome_ids=known_ids,
    )


class LambdaRunnerTests(unittest.TestCase):
    def test_offline_env_strips_proxies_and_sets_offline_flags(self):
        env = build_offline_env(
            {"PATH": "/usr/bin", "HTTP_PROXY": "http://127.0.0.1:8080", "https_proxy": "http://127.0.0.1:8080"},
            offline=True,
        )
        self.assertNotIn("HTTP_PROXY", env)
        self.assertNotIn("https_proxy", env)
        self.assertEqual(env["GRM_BENCH_OFFLINE"], "1")
        self.assertEqual(env["PIP_NO_INDEX"], "1")
        self.assertEqual(env["NO_NETWORK"], "1")

    def test_checkpoint_resume_timing_tsv_and_manifest_pins(self):
        adapter = load_adapter("tbprofiler")
        known_ids = {f"777000.{100 + i}" for i in range(60)}
        result_doc = _clean_result_doc(adapter, known_ids)

        with tempfile.TemporaryDirectory() as tmp:
            tdir = Path(tmp)
            db_dir = tdir / "tbdb"
            db_dir.mkdir()
            (db_dir / "catalogue.tsv").write_text("gene\tdrug\nrpoB\trifampicin\n", encoding="utf-8")
            (db_dir / "version.txt").write_text("3ed1cb99\n", encoding="utf-8")
            db_sha = hash_path_deterministic(db_dir)

            env_lock = tdir / "requirements.lock"
            env_lock.write_text("numpy==1.26.4\n", encoding="utf-8")

            out_dir = tdir / "run_out"
            units = [WorkUnit("u1", "a"), WorkUnit("u2", "b")]
            ran: list[str] = []

            def exec_unit(unit: WorkUnit, env: dict[str, str], rdir: Path) -> Path:
                self.assertEqual(env["GRM_BENCH_OFFLINE"], "1")
                ran.append(unit.unit_id)
                out_res = rdir / f"{unit.unit_id}.result.json"
                write_result_atomic(out_res, result_doc, known_genome_ids=known_ids)
                return out_res

            runner = LambdaRunner(
                adapter,
                out_dir,
                split_scheme="lineage",
                time_provenance={"label": TIME_PROVENANCE_RETRO_SAME_DB},
                offline=True,
                db_path=db_dir,
                expected_db_sha256=db_sha,
                env_lock_path=env_lock,
                known_genome_ids=known_ids,
            )
            s1 = runner.run(units, exec_unit)
            self.assertEqual(s1["status"], "DONE")
            self.assertEqual(s1["units_executed"], 2)
            self.assertEqual(s1["units_skipped"], 0)
            self.assertEqual(ran, ["u1", "u2"])
            self.assertTrue((out_dir / "DONE").is_file())
            self.assertFalse((out_dir / "FAILED").exists())

            # Inspect run-manifest.json pins
            manifest = json.loads((out_dir / "run-manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(manifest["database_sha256"], db_sha)
            self.assertIsNotNone(manifest["runtime_pins"]["env_lock_sha256"])
            self.assertTrue(manifest["offline"])

            # Inspect separate logs/timing.tsv
            with (out_dir / "logs" / "timing.tsv").open(encoding="utf-8", newline="") as fh:
                trows = list(csv.DictReader(fh, delimiter="\t"))
            self.assertEqual([r["unit_id"] for r in trows], ["u1", "u2"])
            self.assertTrue(all(r["status"] == "DONE" for r in trows))

            # Resume skips already-completed units
            ran.clear()
            s2 = runner.run(units, exec_unit)
            self.assertEqual(s2["units_executed"], 0)
            self.assertEqual(s2["units_skipped"], 2)
            self.assertEqual(ran, [])

            # Mutating DB file and resuming MUST fail (config/pin drift refusal!)
            (db_dir / "version.txt").write_text("tampered\n", encoding="utf-8")
            with self.assertRaises(BenchError):
                runner.run(units, exec_unit)
            self.assertTrue((out_dir / "FAILED").is_file())

    def test_offline_missing_db_and_sha_mismatch_refused(self):
        adapter = load_adapter("tbprofiler")
        with tempfile.TemporaryDirectory() as tmp:
            tdir = Path(tmp)
            # Missing db_path in offline mode refused
            r1 = LambdaRunner(
                adapter,
                tdir / "o1",
                split_scheme="lineage",
                time_provenance={"label": TIME_PROVENANCE_RETRO_SAME_DB},
                offline=True,
                db_path=None,
            )
            with self.assertRaisesRegex(BenchError, "--db-path"):
                r1.run([], lambda u, e, r: None)
            self.assertTrue((tdir / "o1" / "FAILED").is_file())

            # Mismatched expected_db_sha256 refused
            db_file = tdir / "db.bin"
            db_file.write_bytes(b"synthetic-db")
            r2 = LambdaRunner(
                adapter,
                tdir / "o2",
                split_scheme="lineage",
                time_provenance={"label": TIME_PROVENANCE_RETRO_SAME_DB},
                offline=True,
                db_path=db_file,
                expected_db_sha256="0" * 64,
            )
            with self.assertRaisesRegex(BenchError, "SHA-256 mismatch"):
                r2.run([], lambda u, e, r: None)

    def test_unit_failure_writes_scrubbed_failed_sentinel_and_hygiene_gate_blocks_done(self):
        adapter = load_adapter("tbprofiler")
        known_ids = {f"777000.{100 + i}" for i in range(60)}
        with tempfile.TemporaryDirectory() as tmp:
            tdir = Path(tmp)
            db_file = tdir / "db.bin"
            db_file.write_bytes(b"synthetic-db")

            # 1. Failing unit scrubs known genome ID from FAILED sentinel
            out1 = tdir / "fail_run"
            runner1 = LambdaRunner(
                adapter,
                out1,
                split_scheme="lineage",
                time_provenance={"label": TIME_PROVENANCE_RETRO_SAME_DB},
                offline=True,
                db_path=db_file,
                known_genome_ids=known_ids,
            )

            def boom(unit, env, rdir):
                raise RuntimeError("crash on genome 777000.100")

            with self.assertRaises(BenchError):
                runner1.run([WorkUnit("u1", None)], boom)
            self.assertFalse((out1 / "DONE").exists())
            failed_text = (out1 / "FAILED").read_text(encoding="utf-8")
            self.assertNotIn("777000.100", failed_text)
            self.assertIn("[genome_id]", failed_text)

            # 2. Leaky result JSON blocked by zero-ID gate before writing run-level DONE
            out2 = tdir / "leak_run"
            runner2 = LambdaRunner(
                adapter,
                out2,
                split_scheme="lineage",
                time_provenance={"label": TIME_PROVENANCE_RETRO_SAME_DB},
                offline=True,
                db_path=db_file,
                known_genome_ids=known_ids,
            )
            leaky_doc = _clean_result_doc(adapter, known_ids)
            leaky_doc["genotype_to_drug_rule"] += " 777000.101"

            def write_leaky(unit, env, rdir):
                p = rdir / "leaky.json"
                p.write_text(json.dumps(leaky_doc), encoding="utf-8")
                return p

            with self.assertRaisesRegex(BenchError, "genome id hit"):
                runner2.run([WorkUnit("u1", None)], write_leaky)
            self.assertFalse((out2 / "DONE").exists())
            self.assertTrue((out2 / "FAILED").is_file())

    def test_cli_dry_run_preflight(self):
        import contextlib
        import io

        with tempfile.TemporaryDirectory() as tmp:
            tdir = Path(tmp)
            db_file = tdir / "db.bin"
            db_file.write_bytes(b"synthetic-db")
            out_dir = tdir / "cli_out"
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                rc = lambda_main([
                    "--tool", "amrfinderplus",
                    "--out", str(out_dir),
                    "--split-scheme", "time",
                    "--time-provenance", TIME_PROVENANCE_RETRO_MODERN_DB,
                    "--fold-cutoff-date", "2023-01-01",
                    "--bvbrc-time-split",
                    "--offline",
                    "--db-path", str(db_file),
                    "--dry-run",
                ])
            self.assertEqual(rc, 0)
            payload = json.loads(buf.getvalue())
            self.assertEqual(payload["status"], "DONE")
            self.assertTrue((out_dir / "run-manifest.json").is_file())
            self.assertTrue((out_dir / "DONE").is_file())


if __name__ == "__main__":
    unittest.main()
