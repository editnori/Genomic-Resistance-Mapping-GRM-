"""Unit, planted-fault, and end-to-end mock tests for runners.py and cli.py."""

import contextlib
import io
import json
import os
import pathlib
import subprocess
import sys
import tempfile
import unittest

from tests.kmer_bench import common

parity = common.load_module("parity")
synthetic = common.load_module("synthetic")
runners = common.load_module("runners")
cli = common.load_module("cli")


class RunnersAndCliTests(unittest.TestCase):
    def _make_cohort(self, root: pathlib.Path):
        spec = synthetic.CohortSpec(
            seed=17,
            n_genomes=4,
            k=21,
            min_support=2,
            core_len=900,
            n_core_contigs=2,
            n_accessory_genes=3,
            accessory_len=240,
            snp_rate=0.008,
            n_plasmids=1,
            plasmid_len=320,
            plasmid_copy_numbers=(5,),
            read_coverage=12.0,
            read_length=75,
            insert_size=130,
            error_rate=0.005,
        )
        return synthetic.generate_cohort(root, spec)

    def test_parse_gnu_time_v_mmss_and_hhmmss_and_fields(self):
        sample_mmss = """\
\tCommand being timed: "grm count --species escherichia_coli"
\tUser time (seconds): 512.34
\tSystem time (seconds): 14.56
\tPercent of CPU this job got: 1075%
\tElapsed (wall clock) time (h:mm:ss or m:ss): 0:49.02
\tMinor (reclaiming a frame) page faults: 123456
\tMajor (requiring I/O) page faults: 7
\tMaximum resident set size (kbytes): 2930764
\tFile system inputs: 8192
\tFile system outputs: 16384
\tExit status: 0
"""
        tm = runners.parse_gnu_time_v(sample_mmss)
        self.assertAlmostEqual(tm.wall_seconds, 49.02, places=6)
        self.assertAlmostEqual(tm.user_seconds, 512.34, places=6)
        self.assertAlmostEqual(tm.sys_seconds, 14.56, places=6)
        self.assertEqual(tm.max_rss_kb, 2_930_764)
        self.assertAlmostEqual(tm.cpu_percent, 1075.0, places=6)
        self.assertEqual(tm.minor_page_faults, 123456)
        self.assertEqual(tm.major_page_faults, 7)
        self.assertEqual(tm.fs_inputs, 8192)
        self.assertEqual(tm.fs_outputs, 16384)
        self.assertEqual(tm.exit_status, 0)

        sample_hhmmss = """\
\tUser time (seconds): 9001.50
\tSystem time (seconds): 120.25
\tElapsed (wall clock) time (h:mm:ss or m:ss): 1:02:03
\tMaximum resident set size (kbytes): 12582912
\tExit status: 137
"""
        tm2 = runners.parse_gnu_time_v(sample_hhmmss)
        self.assertAlmostEqual(tm2.wall_seconds, 3723.0, places=6)
        self.assertEqual(tm2.max_rss_kb, 12_582_912)
        self.assertEqual(tm2.exit_status, 137)

    def test_parse_gnu_time_v_rejects_missing_required_fields(self):
        with self.assertRaisesRegex(ValueError, "missing required fields"):
            runners.parse_gnu_time_v("User time (seconds): 1.0\nExit status: 0\n")
        with self.assertRaisesRegex(ValueError, "unrecognized elapsed"):
            runners.parse_gnu_time_v(
                "User time (seconds): 1.0\nSystem time (seconds): 0.1\n"
                "Elapsed (wall clock) time (h:mm:ss or m:ss): 99\n"
                "Maximum resident set size (kbytes): 1024\n"
            )

    def test_compute_disk_metrics_and_bytes_per_kmer(self):
        dm = runners.compute_disk_metrics(
            peak_temp_bytes=16_000_000,
            final_artifact_bytes=6_000_000,
            kmers_counted=1_000_000,
        )
        self.assertAlmostEqual(dm.bytes_per_kmer, 6.0, places=6)
        self.assertEqual(dm.peak_total_bytes, 22_000_000)
        dm_zero = runners.compute_disk_metrics(0, 128, 0)
        self.assertIsNone(dm_zero.bytes_per_kmer)
        with self.assertRaises(ValueError):
            runners.compute_disk_metrics(-1, 0, 10)

    def test_disk_sampler_tracks_temp_and_artifact_bytes(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            tdir = root / "tmp"
            adir = root / "art"
            tdir.mkdir()
            adir.mkdir()
            sampler = runners.DiskSampler([tdir], [adir], poll_seconds=0.01)
            sampler.sample_once()
            self.assertEqual(sampler.peak_temp_bytes, 0)
            (tdir / "spill.bin").write_bytes(b"x" * 4096)
            (adir / "out.kset").write_bytes(b"y" * 1024)
            sampler.sample_once()
            self.assertEqual(sampler.peak_temp_bytes, 4096)
            self.assertEqual(sampler.peak_total_bytes, 5120)
            (tdir / "spill.bin").unlink()
            sampler.sample_once()
            # Peak retains the high-water mark after temp cleanup.
            self.assertEqual(sampler.peak_temp_bytes, 4096)
            self.assertEqual(runners.directory_bytes([adir]), 1024)

    def test_command_plans_cover_all_nine_tools_across_workloads(self):
        expected_tools = {
            "grm", "fastk", "kmc3", "jellyfish", "ntcard",
            "kmtricks", "ggcat", "cuttlefish2", "bcalm2",
        }
        self.assertEqual(set(runners.TOOL_REGISTRY.keys()), expected_tools)
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            cdir = root / "cohort"
            rdir = root / "run"
            count = 0
            for tool, meta in runners.TOOL_REGISTRY.items():
                for wl in meta["workloads"]:
                    plan = runners.build_command_plan(tool, wl, cdir, rdir, k=31, threads=4, mem_gb=8)
                    self.assertEqual(plan.tool, tool)
                    self.assertEqual(plan.workload, wl)
                    self.assertGreaterEqual(len(plan.steps), 1)
                    self.assertTrue(all(isinstance(arg, str) and arg for step in plan.steps for arg in step))
                    count += 1
            self.assertEqual(count, 20)
            # Check GRM T3 variants emit distinct flags.
            p_cb = runners.build_command_plan("grm", "T3", cdir, rdir, variant="column-blocks-64")
            p_ds = runners.build_command_plan("grm", "T3", cdir, rdir, variant="disk-sparse-256")
            self.assertIn("--column-blocks", p_cb.steps[0])
            self.assertIn("--disk-sparse", p_ds.steps[0])

    def test_command_plan_rejects_unknown_tool_or_unsupported_workload(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = pathlib.Path(tmp)
            with self.assertRaisesRegex(ValueError, "unknown tool"):
                runners.build_command_plan("nonexistent", "T1", p, p)
            with self.assertRaisesRegex(ValueError, "unknown workload"):
                runners.build_command_plan("grm", "T99", p, p)
            with self.assertRaisesRegex(ValueError, "does not support"):
                runners.build_command_plan("bcalm2", "T2", p, p)

    def test_real_mode_gracefully_skips_missing_external_binary(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            cdir = root / "cohort"
            manifest = self._make_cohort(cdir)
            plan = runners.build_command_plan(
                "fastk", "T1", cdir, root / "work",
                k=21, binary_override="definitely_missing_kmer_binary_xyz_9999",
            )
            # Force provenance check against the missing override via capture_tool_provenance.
            prov = runners.capture_tool_provenance("fastk", binary_override="definitely_missing_kmer_binary_xyz_9999")
            self.assertFalse(prov.installed)
            self.assertIsNone(prov.binary_sha256)
            # In real mode, if tool is absent from PATH, skip_missing=True returns skipped_missing_binary.
            old_path = os.environ.get("PATH", "")
            os.environ["PATH"] = ""
            try:
                rec = runners.execute_case(plan, manifest, cdir, mode="real", skip_missing=True)
                self.assertEqual(rec["status"], "skipped_missing_binary")
                with self.assertRaises(FileNotFoundError):
                    runners.execute_case(plan, manifest, cdir, mode="real", skip_missing=False)
            finally:
                os.environ["PATH"] = old_path

    def test_mock_end_to_end_t1_through_t5_and_parity_pass(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            cdir = root / "cohort"
            odir = root / "out"
            manifest = self._make_cohort(cdir)
            results = []
            for wl in ("T1", "T2", "T3", "T4", "T5"):
                for tool, meta in runners.TOOL_REGISTRY.items():
                    if wl not in meta["workloads"]:
                        continue
                    plan = runners.build_command_plan(tool, wl, cdir, odir / "work", k=21, threads=4)
                    rec = runners.execute_case(plan, manifest, cdir, mode="mock")
                    self.assertEqual(rec["status"], "ok", (wl, tool, rec))
                    self.assertTrue(rec["parity"]["ok"], (wl, tool))
                    results.append(rec)
            report = runners.build_bench_report(results, manifest, strict_hygiene=True)
            self.assertEqual(report["format"], "grm-kmer-bench-v1")
            jpath, tpath = runners.write_bench_outputs(odir, report, strict_hygiene=True)
            self.assertTrue(jpath.is_file() and tpath.is_file())
            tsv_lines = tpath.read_text(encoding="utf-8").splitlines()
            self.assertEqual(len(tsv_lines), len(results) + 1)
            self.assertEqual(tsv_lines[0], "\t".join(runners.TSV_COLUMNS))

    def test_mock_planted_faults_fail_parity_for_t1_t2_t3_t4_t5(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            cdir = root / "cohort"
            manifest = self._make_cohort(cdir)
            cases = (
                ("T1", "kmc3", "drop_kmer"),
                ("T1", "kmc3", "alter_count"),
                ("T2", "fastk", "drop_kmer"),
                ("T2", "fastk", "alter_count"),
                ("T3", "grm", "flip_bit"),
                ("T4", "bcalm2", "merge_unitig"),
                ("T5", "kmtricks", "flip_bit"),
            )
            for wl, tool, fault in cases:
                with self.subTest(workload=wl, tool=tool, fault=fault):
                    plan = runners.build_command_plan(tool, wl, cdir, root / "work", k=21)
                    rec = runners.execute_case(plan, manifest, cdir, mode="mock", inject_fault=fault)
                    self.assertEqual(rec["status"], "parity_failed")
                    self.assertFalse(rec["parity"]["ok"])

    def test_cli_subcommands_generate_verify_plan_run_summarize(self):
        with tempfile.TemporaryDirectory() as tmp, contextlib.redirect_stdout(io.StringIO()):
            root = pathlib.Path(tmp)
            cdir = root / "cohort"
            self.assertEqual(
                cli.main([
                    "generate", "--out", str(cdir), "--seed", "19", "--n-genomes", "4",
                    "--k", "21", "--core-len", "600", "--n-accessory-genes", "3",
                    "--accessory-len", "160", "--n-plasmids", "1", "--plasmid-len", "240",
                    "--plasmid-copies", "4", "--coverage", "10", "--read-length", "75",
                    "--insert-size", "130",
                ]),
                0,
            )
            # verify-parity set via .kset
            fna = cdir / "assemblies" / "syn_g0000.fna"
            ref = parity.count_file(fna, k=21, min_abundance=1)
            kset_path = root / "syn_g0000.kset"
            parity.write_kset(kset_path, "syn_g0000", parity.KmerSpace(21), ref.mixed_sorted)
            par_json = root / "par_set.json"
            self.assertEqual(
                cli.main([
                    "verify-parity", "--mode", "set", "--k", "21", "--min-abundance", "1",
                    "--reference", str(fna), "--candidate-kset", str(kset_path),
                    "--json-out", str(par_json),
                ]),
                0,
            )
            self.assertTrue(json.loads(par_json.read_text(encoding="utf-8"))["ok"])

            # verify-parity abundance via truth dump
            ab_json = root / "par_ab.json"
            self.assertEqual(
                cli.main([
                    "verify-parity", "--mode", "abundance", "--k", "21", "--min-abundance", "2",
                    "--reference", str(cdir / "reads" / "syn_g0000_R1.fastq"),
                    str(cdir / "reads" / "syn_g0000_R2.fastq"),
                    "--candidate-dump", str(cdir / "truth" / "syn_g0000.reads.counts.tsv"),
                    "--json-out", str(ab_json),
                ]),
                0,
            )
            self.assertTrue(json.loads(ab_json.read_text(encoding="utf-8"))["ok"])

            # plan
            plan_json = root / "plan.json"
            self.assertEqual(
                cli.main([
                    "plan", "--cohort-dir", str(cdir), "--out-dir", str(root / "pout"),
                    "--workloads", "all", "--tools", "all", "--k", "21",
                    "--include-grm-variants", "--json-out", str(plan_json),
                ]),
                0,
            )
            plan_doc = json.loads(plan_json.read_text(encoding="utf-8"))
            self.assertEqual(plan_doc["format"], "grm-kmer-bench-plan-v1")
            self.assertEqual(len(plan_doc["cases"]), 22)

            # run --dry-run
            self.assertEqual(
                cli.main([
                    "run", "--cohort-dir", str(cdir), "--out-dir", str(root / "dry"),
                    "--workloads", "T1,T3", "--tools", "grm,kmc3", "--k", "21", "--dry-run",
                ]),
                0,
            )
            dry_doc = json.loads((root / "dry" / "bench.json").read_text(encoding="utf-8"))
            self.assertTrue(all(r["status"] == "dry_run" for r in dry_doc["results"]))

            # run --mock
            self.assertEqual(
                cli.main([
                    "run", "--cohort-dir", str(cdir), "--out-dir", str(root / "mock"),
                    "--workloads", "all", "--tools", "all", "--k", "21", "--mock",
                ]),
                0,
            )

            # summarize
            sum_tsv = root / "sum.tsv"
            sum_md = root / "sum.md"
            self.assertEqual(
                cli.main([
                    "summarize", "--inputs", str(root / "mock" / "bench.json"),
                    "--tsv-out", str(sum_tsv), "--md-out", str(sum_md),
                ]),
                0,
            )
            self.assertIn("| Workload | Tool |", sum_md.read_text(encoding="utf-8"))
            self.assertEqual(sum_tsv.read_text(encoding="utf-8"), (root / "mock" / "summary.tsv").read_text(encoding="utf-8"))

    def test_cli_verify_parity_exits_nonzero_on_planted_fault(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            cdir = root / "cohort"
            self._make_cohort(cdir)
            truth = cdir / "truth" / "syn_g0000.assembly.counts.tsv"
            lines = truth.read_text(encoding="utf-8").splitlines()
            # Drop first k-mer -> exit 1.
            bad_dump = root / "dropped.tsv"
            bad_dump.write_text("\n".join(lines[1:]) + "\n", encoding="utf-8")
            out_json = root / "fail.json"
            rc = cli.main([
                "verify-parity", "--mode", "set", "--k", "21", "--min-abundance", "1",
                "--reference", str(cdir / "assemblies" / "syn_g0000.fna"),
                "--candidate-dump", str(bad_dump), "--json-out", str(out_json),
            ])
            self.assertEqual(rc, 1)
            self.assertFalse(json.loads(out_json.read_text(encoding="utf-8"))["ok"])

    def test_mock_bench_output_deterministic_across_pythonhashseed(self):
        cli_path = common.SCRIPT_DIR / "cli.py"
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            cdir = root / "cohort"
            self._make_cohort(cdir)
            json_blobs = {}
            tsv_blobs = {}
            for seed in ("0", "1", "4242"):
                odir = root / f"run_{seed}"
                env = {
                    "PATH": "/usr/bin:/bin",
                    "PYTHONHASHSEED": seed,
                }
                proc = subprocess.run(
                    [
                        sys.executable, "-s", "-P", "-B", "-W", "error", str(cli_path),
                        "run", "--cohort-dir", str(cdir), "--out-dir", str(odir),
                        "--workloads", "all", "--tools", "all", "--k", "21", "--mock",
                    ],
                    env=env,
                    capture_output=True,
                    text=True,
                    check=False,
                )
                self.assertEqual(proc.returncode, 0, proc.stderr)
                json_blobs[seed] = (odir / "bench.json").read_text(encoding="utf-8")
                tsv_blobs[seed] = (odir / "summary.tsv").read_text(encoding="utf-8")
            self.assertEqual(json_blobs["0"], json_blobs["1"])
            self.assertEqual(json_blobs["0"], json_blobs["4242"])
            self.assertEqual(tsv_blobs["0"], tsv_blobs["1"])
            self.assertEqual(tsv_blobs["0"], tsv_blobs["4242"])

    def test_hygiene_scan_catches_genome_ids_em_dashes_and_private_paths(self):
        self.assertEqual(runners.scan_for_leaks("clean syn_g0001 report 0.0245"), [])
        self.assertIn("bvbrc_id_shape", runners.scan_for_leaks("genome 562.12345 leaked"))
        self.assertIn("known_genome_id", runners.scan_for_leaks("id 12.34 here", forbidden_ids={"12.34"}))
        self.assertIn("em_dash", runners.scan_for_leaks("text \u2014 dash"))
        self.assertIn("private_home_path", runners.scan_for_leaks("/" + "home/alice/work"))
        self.assertEqual(runners.scan_for_leaks("/" + "home/user/work"), [])
        self.assertIn("windows_user_path", runners.scan_for_leaks("C:\\" + "Users\\alice\\work"))
        self.assertIn("thread_link", runners.scan_for_leaks("t3" + "-thread://" + "abc"))
        self.assertIn("box_id", runners.scan_for_leaks("bx" + "_a1b2c3d4"))


if __name__ == "__main__":
    unittest.main()
