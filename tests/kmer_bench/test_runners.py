"""Unit, planted-fault, stub real-mode, and end-to-end tests for runners.py and cli.py."""

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


def _write_exec_script(path: pathlib.Path, body: str) -> pathlib.Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(f"#!{sys.executable}\n{body}", encoding="utf-8")
    os.chmod(path, 0o755)
    return path


class RunnersAndCliTests(unittest.TestCase):
    def _make_cohort(self, root: pathlib.Path, **overrides):
        params = dict(
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
        params.update(overrides)
        return synthetic.generate_cohort(root, synthetic.CohortSpec(**params))

    def test_parse_gnu_time_v_mmss_past_one_minute_and_hhmmss_and_fields(self):
        # Covers <1 min, >1 min m:ss (catches F30), and h:mm:ss (F31).
        sample_mmss = """\
\tCommand being timed: "grm count --species synthetic_sp"
\tUser time (seconds): 512.34
\tSystem time (seconds): 14.56
\tPercent of CPU this job got: 1075%
\tElapsed (wall clock) time (h:mm:ss or m:ss): 2:05.50
\tMinor (reclaiming a frame) page faults: 123456
\tMajor (requiring I/O) page faults: 7
\tMaximum resident set size (kbytes): 2930764
\tFile system inputs: 8192
\tFile system outputs: 16384
\tExit status: 0
"""
        tm = runners.parse_gnu_time_v(sample_mmss)
        self.assertAlmostEqual(tm.wall_seconds, 125.50, places=6)
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
            self.assertEqual(sampler.peak_temp_bytes, 4096)
            self.assertEqual(runners.directory_bytes([adir]), 1024)

    def test_command_plans_cover_all_nine_tools_flags_and_comparability(self):
        # Covers all 9 tools, T5 floor=2 (F34), grm T1 --k (F35), FastK -t2 (F36),
        # Jellyfish -C (F37), shared grm_data (C2), and grm T2/T5 not_comparable (C12).
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
                    plan = runners.build_command_plan(tool, wl, cdir, rdir, k=21, threads=4, mem_gb=8)
                    self.assertEqual(plan.tool, tool)
                    self.assertEqual(plan.workload, wl)
                    self.assertGreaterEqual(len(plan.steps), 1)
                    self.assertTrue(all(isinstance(arg, str) and arg for step in plan.steps for arg in step))
                    want_min = 1 if wl in ("T1", "T3", "T4") else 2
                    self.assertEqual(plan.min_abundance, want_min, (tool, wl))
                    count += 1
            self.assertEqual(count, 20)

            # F35: grm T1 passes requested --k (21 and 31).
            grm_t1_21 = runners.build_command_plan("grm", "T1", cdir, rdir, k=21)
            grm_t1_31 = runners.build_command_plan("grm", "T1", cdir, rdir, k=31)
            self.assertEqual(grm_t1_21.steps[0][grm_t1_21.steps[0].index("--k") + 1], "21")
            self.assertEqual(grm_t1_31.steps[0][grm_t1_31.steps[0].index("--k") + 1], "31")
            self.assertEqual(grm_t1_31.output_scope, "cohort")

            # C2: T1, T3 (default & column-blocks-64), and T4 share one grm_data directory.
            grm_t3 = runners.build_command_plan("grm", "T3", cdir, rdir, k=31)
            grm_t4 = runners.build_command_plan("grm", "T4", cdir, rdir, k=31)
            self.assertEqual(grm_t1_31.grm_data_dir, grm_t3.grm_data_dir)
            self.assertEqual(grm_t1_31.grm_data_dir, grm_t4.grm_data_dir)

            # C12: grm T2 and T5 are marked not_comparable.
            grm_t2 = runners.build_command_plan("grm", "T2", cdir, rdir, k=31)
            grm_t5 = runners.build_command_plan("grm", "T5", cdir, rdir, k=31)
            self.assertFalse(grm_t2.comparable)
            self.assertFalse(grm_t5.comparable)
            self.assertIn("min_abundance=1", grm_t2.not_comparable_reason or "")

            # F36: FastK floor is -t1 for T1 and -t2 for T2.
            fk1 = runners.build_command_plan("fastk", "T1", cdir, rdir, k=21)
            fk2 = runners.build_command_plan("fastk", "T2", cdir, rdir, k=21)
            self.assertIn("-t1", fk1.steps[0])
            self.assertIn("-t2", fk2.steps[0])

            # F37: Jellyfish uses -C (canonical) in both T1 and T2.
            jf1 = runners.build_command_plan("jellyfish", "T1", cdir, rdir, k=21)
            jf2 = runners.build_command_plan("jellyfish", "T2", cdir, rdir, k=21)
            self.assertIn("-C", jf1.steps[0])
            self.assertIn("-C", jf2.steps[0])

            # F34: T5 floor is 2 in kmtricks (--hard-min 2) and ggcat (-s 2).
            kt5 = runners.build_command_plan("kmtricks", "T5", cdir, rdir, k=21)
            gg5 = runners.build_command_plan("ggcat", "T5", cdir, rdir, k=21)
            self.assertEqual(kt5.steps[0][kt5.steps[0].index("--hard-min") + 1], "2")
            self.assertEqual(gg5.steps[0][gg5.steps[0].index("-s") + 1], "2")

            # GRM T3 variants emit distinct flags.
            p_cb = runners.build_command_plan("grm", "T3", cdir, rdir, variant="column-blocks-64")
            p_ds = runners.build_command_plan("grm", "T3", cdir, rdir, variant="disk-sparse-256")
            self.assertIn("--column-blocks", p_cb.steps[0])
            self.assertIn("--disk-sparse", p_ds.steps[0])

    def test_gzip_reads_manifest_filenames_and_dash_leading_cohort_path(self):
        # Covers C9 (--gzip-reads filenames from manifest) and C10 (dash-leading cohort directory name).
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            dash_cdir = root / "-dash-cohort"
            manifest = self._make_cohort(dash_cdir, n_genomes=2, gzip_reads=True)
            rdir = root / "run"

            # Call with a relative dash-leading path to confirm build_command_plan resolves it (C10).
            old_cwd = os.getcwd()
            os.chdir(root)
            try:
                rel_cohort = pathlib.Path("-dash-cohort")
                for tool, meta in runners.TOOL_REGISTRY.items():
                    for wl in meta["workloads"]:
                        plan = runners.build_command_plan(tool, wl, rel_cohort, rdir, k=21)
                        for step in plan.steps:
                            for arg in step:
                                self.assertFalse(arg.startswith("-dash-cohort"), (tool, wl, arg))
                                self.assertFalse(arg.startswith("@-dash-cohort"), (tool, wl, arg))
                # C9: T2 plans for fastk, jellyfish, ntcard reference existing .fastq.gz files.
                for tool in ("fastk", "jellyfish", "ntcard"):
                    plan = runners.build_command_plan(tool, "T2", rel_cohort, rdir, k=21, cohort_manifest=manifest)
                    step = plan.steps[0]
                    gz_args = [a for a in step if a.endswith(".fastq.gz")]
                    self.assertEqual(len(gz_args), 2, (tool, step))
                    for p_str in gz_args:
                        self.assertTrue(pathlib.Path(p_str).is_file(), p_str)
                # C4 + C9: materialized FOF lists hold resolved absolute .fastq.gz paths.
                kmc_plan = runners.build_command_plan("kmc3", "T2", rel_cohort, rdir, k=21, cohort_manifest=manifest)
                lists = runners.materialize_run_lists(kmc_plan, manifest, rel_cohort)
                kmc_lines = lists["kmc_fof.txt"].read_text(encoding="utf-8").splitlines()
                self.assertEqual(len(kmc_lines), 2)
                for ln in kmc_lines:
                    self.assertTrue(pathlib.Path(ln).is_absolute() and ln.endswith(".fastq.gz") and pathlib.Path(ln).is_file())
            finally:
                os.chdir(old_cwd)

    def test_grm_provenance_uses_sha256_without_version_flag(self):
        # Covers C11: grm binary fingerprint records SHA-256 even when binary rejects --version.
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            stub = _write_exec_script(root / "grm", "import sys\nsys.exit(2)\n")
            prov = runners.capture_tool_provenance("grm", binary_override=str(stub))
            self.assertTrue(prov.installed)
            self.assertIsNotNone(prov.binary_sha256)
            self.assertEqual(len(prov.binary_sha256), 64)
            self.assertIsNotNone(prov.version_line)
            self.assertIn(f"sha256:{prov.binary_sha256[:16]}", prov.version_line)

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
            prov = runners.capture_tool_provenance("fastk", binary_override="definitely_missing_kmer_binary_xyz_9999")
            self.assertFalse(prov.installed)
            self.assertIsNone(prov.binary_sha256)
            rec = runners.execute_case(plan, manifest, cdir, mode="real", skip_missing=True)
            self.assertEqual(rec["status"], "skipped_missing_binary")
            with self.assertRaises(FileNotFoundError):
                runners.execute_case(plan, manifest, cdir, mode="real", skip_missing=False)

    def test_real_mode_stubs_empty_fail_stop_parity_fail_pass_and_unchecked(self):
        # Covers C1, C2, C4, C5, C12, F40, F41, F61 in real execution mode using deterministic Python stubs.
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            cdir = root / "cohort"
            manifest = self._make_cohort(cdir, n_genomes=4, k=21, generate_reads=False)

            # (a) Zero-byte exit-0 stub MUST NOT return "ok" (C1, F41).
            noop_stub = _write_exec_script(root / "bins" / "noop", "import sys\nsys.exit(0)\n")
            for tool, wl in (("bcalm2", "T4"), ("jellyfish", "T1"), ("kmc3", "T1"), ("grm", "T1")):
                plan = runners.build_command_plan(
                    tool, wl, cdir, root / f"run_noop_{tool}", k=21, binary_override=str(noop_stub), cohort_manifest=manifest
                )
                rec = runners.execute_case(plan, manifest, cdir, mode="real", skip_missing=False)
                self.assertEqual(rec["mode"], "real")
                self.assertEqual(rec["status"], "failed_empty_output", (tool, wl, rec))
                self.assertFalse(rec["parity"]["ok"])

            # (b) Multi-step failing stub stops immediately after step 1 (F40, F41).
            marker = root / "steps_called.log"
            fail_stub = _write_exec_script(
                root / "bins" / "failstep",
                (
                    f"import sys\n"
                    f"if sys.argv[1:] != ['-version']:\n"
                    f"    with open({str(marker)!r}, 'a') as f: f.write(' '.join(sys.argv[1:]) + '\\n')\n"
                    f"sys.exit(1)\n"
                ),
            )
            base_plan = runners.build_command_plan(
                "bcalm2", "T4", cdir, root / "run_fail", k=21, binary_override=str(fail_stub), cohort_manifest=manifest
            )
            two_step = runners.CommandPlan(
                **{
                    **base_plan.__dict__,
                    "steps": ((str(fail_stub), "step1"), (str(fail_stub), "step2")),
                }
            )
            rec_fail = runners.execute_case(two_step, manifest, cdir, mode="real", skip_missing=False)
            self.assertEqual(rec_fail["status"], "failed_exit_1")
            self.assertEqual(marker.read_text(encoding="utf-8").splitlines(), ["step1"])

            # (c) Parity-mismatching stub for grm T1 -> status='parity_failed', cli.py run exits 1 (C1, F41, F61).
            bad_grm = _write_exec_script(
                root / "bins" / "bad_grm",
                f"""import pathlib, sys
sys.path.insert(0, {str(common.SCRIPT_DIR)!r})
import parity
args = sys.argv[1:]
data_dir = pathlib.Path(args[args.index('--data') + 1])
sp = args[args.index('--species') + 1]
k = int(args[args.index('--k') + 1])
space = parity.KmerSpace(k)
kdir = data_dir / 'species' / sp / 'ksets'
kdir.mkdir(parents=True, exist_ok=True)
gids = [ln.strip() for ln in (data_dir / 'species' / sp / 'genomes.txt').read_text().splitlines() if ln.strip()]
for idx, gid in enumerate(gids):
    res = parity.count_file(data_dir / 'fasta' / f'{{gid}}.fna.gz', k=k, min_abundance=1)
    vals = list(res.mixed_sorted)
    if idx == 0 and vals:
        vals = vals[1:]  # Drop 1 k-mer -> parity failure!
    parity.write_kset(kdir / f'{{gid}}.kset', gid, space, vals)
""",
            )
            with contextlib.redirect_stdout(io.StringIO()):
                rc_bad = cli.main([
                    "run", "--cohort-dir", str(cdir), "--out-dir", str(root / "cli_bad_grm"),
                    "--workloads", "T1", "--tools", "grm", "--k", "21",
                    "--grm-bin", str(bad_grm), "--fail-missing",
                ])
            self.assertEqual(rc_bad, 1)
            bad_doc = json.loads((root / "cli_bad_grm" / "bench.json").read_text(encoding="utf-8"))
            self.assertEqual(bad_doc["results"][0]["status"], "parity_failed")
            self.assertFalse(bad_doc["results"][0]["parity"]["ok"])

            # (d) Passing grm T1 stub: verifies C1 (real parity), C2 (grm_data staging), C5 (all-genome bytes/k-mer).
            good_grm = _write_exec_script(
                root / "bins" / "good_grm",
                f"""import json, pathlib, sys
sys.path.insert(0, {str(common.SCRIPT_DIR)!r})
import parity
args = sys.argv[1:]
data_dir = pathlib.Path(args[args.index('--data') + 1])
sp = args[args.index('--species') + 1]
k = int(args[args.index('--k') + 1])
space = parity.KmerSpace(k)
kdir = data_dir / 'species' / sp / 'ksets'
kdir.mkdir(parents=True, exist_ok=True)
gids = [ln.strip() for ln in (data_dir / 'species' / sp / 'genomes.txt').read_text().splitlines() if ln.strip()]
for gid in gids:
    res = parity.count_file(data_dir / 'fasta' / f'{{gid}}.fna.gz', k=k, min_abundance=1)
    parity.write_kset(kdir / f'{{gid}}.kset', gid, space, res.mixed_sorted)
""",
            )
            with contextlib.redirect_stdout(io.StringIO()):
                rc_good = cli.main([
                    "run", "--cohort-dir", str(cdir), "--out-dir", str(root / "cli_good_grm"),
                    "--workloads", "T1,T2", "--tools", "grm", "--k", "21",
                    "--grm-bin", str(good_grm), "--fail-missing",
                ])
            self.assertEqual(rc_good, 0)
            good_doc = json.loads((root / "cli_good_grm" / "bench.json").read_text(encoding="utf-8"))
            by_wl = {r["workload"]: r for r in good_doc["results"]}
            t1_rec = by_wl["T1"]
            self.assertEqual(t1_rec["status"], "ok")
            self.assertTrue(t1_rec["parity"]["ok"])
            self.assertEqual(t1_rec["parity"]["genomes_passed"], 4)
            # C5: denominator is sum of distinct k-mers across ALL 4 genomes, not genome 0 alone.
            total_k = sum(g["assembly_truth"]["distinct_kmers_min1"] for g in manifest["genomes"])
            g0_k = manifest["genomes"][0]["assembly_truth"]["distinct_kmers_min1"]
            self.assertEqual(t1_rec["disk"]["kmers_counted"], total_k)
            self.assertGreater(total_k, g0_k * 3)
            # C12: T2 is not_comparable.
            self.assertEqual(by_wl["T2"]["status"], "not_comparable")
            self.assertFalse(by_wl["T2"]["comparable"])

            # (e) Non-empty binary unparsed output (FastK .ktab binary bytes) -> status='unchecked', never 'ok' (C1).
            bin_fastk = _write_exec_script(
                root / "bins" / "FastK",
                """import pathlib, sys
for arg in sys.argv[1:]:
    if arg.startswith('-N'):
        p = pathlib.Path(arg[2:] + '.ktab')
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_bytes(b'\\x00FASTKBIN\\x00\\x01\\x02' * 32)
""",
            )
            fk_plan = runners.build_command_plan(
                "fastk", "T1", cdir, root / "run_fk_bin", k=21, binary_override=str(bin_fastk), cohort_manifest=manifest
            )
            rec_fk = runners.execute_case(fk_plan, manifest, cdir, mode="real", skip_missing=False)
            self.assertEqual(rec_fk["status"], "unchecked")
            self.assertIsNone(rec_fk["parity"]["ok"])
            self.assertFalse(rec_fk["parity"]["checked"])

            # (f) Passing bcalm2 T4 stub run from an unrelated cwd verifies C4 (absolute paths in assemblies.fof)
            # and real unitig FASTA parity parsing (C1).
            good_bcalm = _write_exec_script(
                root / "bins" / "bcalm",
                f"""import pathlib, sys
sys.path.insert(0, {str(common.SCRIPT_DIR)!r})
import parity
args = sys.argv[1:]
fof = pathlib.Path(args[args.index('-in') + 1])
out_prefix = args[args.index('-out') + 1]
paths = [pathlib.Path(ln.strip()) for ln in fof.read_text().splitlines() if ln.strip()]
assert all(p.is_absolute() and p.is_file() for p in paths), paths
# Write each contig from the single-contig test cohort as a unitig.
out_fa = pathlib.Path(out_prefix + '.unitigs.fa')
out_fa.parent.mkdir(parents=True, exist_ok=True)
seqs = []
for p in paths:
    for name, seq in parity.iter_sequence_records(p):
        seqs.append((name, seq.decode('ascii')))
out_fa.write_text(''.join(f'>{{n}}\\n{{s}}\\n' for n, s in seqs), encoding='ascii')
""",
            )
            # Use a tiny 2-genome cohort with 1 non-overlapping contig per genome so contigs == unitigs.
            u_cdir = root / "u_cohort"
            (u_cdir / "assemblies").mkdir(parents=True)
            rng_u = synthetic.DeterministicRng(505)
            s_a = rng_u.random_dna(120)
            s_b = rng_u.random_dna(120)
            (u_cdir / "assemblies" / "g0.fna").write_text(f">u0\n{s_a}\n", encoding="ascii")
            (u_cdir / "assemblies" / "g1.fna").write_text(f">u1\n{s_b}\n", encoding="ascii")
            u_manifest = {
                "format": synthetic.COHORT_FORMAT,
                "spec": {"seed": 1, "n_genomes": 2, "k": 21, "min_support": 2, "core_len": 120, "n_plasmids": 0, "generate_reads": False, "read_coverage": 0, "error_rate": 0},
                "cohort_truth": {"unitigs": {"graph_kmers": 200, "n_unitigs": 2}},
                "genomes": [
                    {"genome": "g0", "fasta": "assemblies/g0.fna", "assembly_truth": {"distinct_kmers_min1": 100}},
                    {"genome": "g1", "fasta": "assemblies/g1.fna", "assembly_truth": {"distinct_kmers_min1": 100}},
                ],
            }
            unrelated = root / "unrelated_cwd"
            unrelated.mkdir()
            old_cwd = os.getcwd()
            os.chdir(unrelated)
            try:
                bc_plan = runners.build_command_plan(
                    "bcalm2", "T4", u_cdir, root / "run_bcalm", k=21, binary_override=str(good_bcalm), cohort_manifest=u_manifest
                )
                rec_bc = runners.execute_case(bc_plan, u_manifest, u_cdir, mode="real", skip_missing=False)
            finally:
                os.chdir(old_cwd)
            self.assertEqual(rec_bc["status"], "ok", rec_bc)
            self.assertTrue(rec_bc["parity"]["ok"])

    def test_mock_end_to_end_self_check_honesty_and_exclusion_from_real_summary(self):
        # Covers C3 (mode='mock_self_check', self_check labels, excluded from real comparison summary),
        # C5 (grm T1 denominator = all genomes), F39 (grm mock 6 B/k-mer vs 5 B/k-mer), and F42 (support=2).
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
                    plan = runners.build_command_plan(tool, wl, cdir, odir / "work", k=21, threads=4, cohort_manifest=manifest)
                    rec = runners.execute_case(plan, manifest, cdir, mode="mock")
                    self.assertEqual(rec["mode"], "mock_self_check")
                    self.assertTrue(rec["synthetic_placeholder"])
                    if plan.comparable:
                        self.assertEqual(rec["status"], "ok", (wl, tool, rec))
                        self.assertTrue(rec["parity"]["ok"], (wl, tool))
                        self.assertTrue(rec["parity"]["self_check"])
                        self.assertEqual(rec["parity"]["check_type"], "mock_self_check")
                    else:
                        self.assertEqual(rec["status"], "not_comparable", (wl, tool, rec))
                    results.append(rec)

            # F39 & C5: grm T1 mock uses 6 B/k-mer over all genomes; kmc3 T1 uses 5 B/k-mer over genome 0.
            by_key = {(r["workload"], r["tool"]): r for r in results}
            total_k = sum(g["assembly_truth"]["distinct_kmers_min1"] for g in manifest["genomes"])
            g0_k = manifest["genomes"][0]["assembly_truth"]["distinct_kmers_min1"]
            self.assertEqual(by_key[("T1", "grm")]["disk"]["kmers_counted"], total_k)
            self.assertAlmostEqual(by_key[("T1", "grm")]["disk"]["bytes_per_kmer"], 6.0, places=6)
            self.assertEqual(by_key[("T1", "kmc3")]["disk"]["kmers_counted"], g0_k)
            self.assertAlmostEqual(by_key[("T1", "kmc3")]["disk"]["bytes_per_kmer"], 5.0, places=6)

            # F42: mock T3 and T5 use min_support=2 (kept < distinct).
            asm_mat = manifest["cohort_truth"]["assembly_matrix"]
            self.assertLess(asm_mat["kmers_kept"], asm_mat["kmers_distinct"])
            self.assertEqual(by_key[("T3", "grm")]["parity"]["matrix"]["expected_kept"], asm_mat["kmers_kept"])
            rd_mat = manifest["cohort_truth"]["reads_matrix_min2"]
            self.assertEqual(by_key[("T5", "kmtricks")]["parity"]["matrix"]["expected_kept"], rd_mat["kmers_kept"])

            report = runners.build_bench_report(results, manifest, strict_hygiene=True)
            self.assertEqual(report["format"], "grm-kmer-bench-v1")
            jpath, tpath = runners.write_bench_outputs(odir, report, strict_hygiene=True)
            self.assertTrue(jpath.is_file() and tpath.is_file())
            tsv_lines = tpath.read_text(encoding="utf-8").splitlines()
            self.assertEqual(len(tsv_lines), len(results) + 1)
            self.assertEqual(tsv_lines[0], "\t".join(runners.TSV_COLUMNS))
            self.assertIn("mode", runners.TSV_COLUMNS)
            mode_idx = runners.TSV_COLUMNS.index("mode")
            for ln in tsv_lines[1:]:
                self.assertEqual(ln.split("\t")[mode_idx], "mock_self_check")

            # C3: real comparison summary MUST exclude mock_self_check rows by default.
            self.assertEqual(runners.filter_real_comparable_results(results), [])
            md_real_only = runners.render_summary_markdown(results)
            self.assertNotIn("mock_self_check`", md_real_only)
            self.assertIn("Excluded", md_real_only)
            md_with_mock = runners.render_summary_markdown(results, include_mock=True)
            self.assertIn("`mock_self_check`", md_with_mock)

    def test_mock_planted_faults_fail_parity_and_tsv_parity_column_tracks_false(self):
        # Covers mock fault injection AND F44 (TSV parity column reflects true/false/empty accurately).
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
            failed_recs = []
            for wl, tool, fault in cases:
                with self.subTest(workload=wl, tool=tool, fault=fault):
                    plan = runners.build_command_plan(tool, wl, cdir, root / "work", k=21, cohort_manifest=manifest)
                    rec = runners.execute_case(plan, manifest, cdir, mode="mock", inject_fault=fault)
                    self.assertEqual(rec["status"], "parity_failed")
                    self.assertFalse(rec["parity"]["ok"])
                    failed_recs.append(rec)

            # F44: TSV parity_ok column must be 'false' for failed rows, 'true' for passing, '' for None.
            ok_plan = runners.build_command_plan("kmc3", "T1", cdir, root / "work", k=21, cohort_manifest=manifest)
            ok_rec = runners.execute_case(ok_plan, manifest, cdir, mode="mock")
            nc_plan = runners.build_command_plan("grm", "T2", cdir, root / "work", k=21, cohort_manifest=manifest)
            nc_rec = runners.execute_case(nc_plan, manifest, cdir, mode="mock")
            tsv = runners.render_summary_tsv([failed_recs[0], ok_rec, nc_rec])
            par_idx = runners.TSV_COLUMNS.index("parity_ok")
            stat_idx = runners.TSV_COLUMNS.index("status")
            by_status = {ln.split("\t")[stat_idx]: ln.split("\t")[par_idx] for ln in tsv.splitlines()[1:]}
            self.assertEqual(by_status["parity_failed"], "false")
            self.assertEqual(by_status["ok"], "true")
            self.assertEqual(by_status["not_comparable"], "")

    def test_cli_subcommands_generate_verify_plan_run_summarize_and_flags(self):
        # Covers F55 (generate collapses plasmids by default), F56 (verify-parity rejects mismatched kset k),
        # F58 (--max-cap respected), F60 (plan & run propagate --k), and C3 (summarize excludes mock by default).
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
            # F55: default generate collapses plasmid copies (1 plasmid_0 contig in FASTA).
            m_default = json.loads((cdir / "manifest.json").read_text(encoding="utf-8"))
            self.assertTrue(m_default["spec"]["collapse_plasmids_in_assembly"])
            recs_def = list(parity.iter_sequence_records(cdir / "assemblies" / "syn_g0000.fna"))
            self.assertEqual(sum(1 for name, _ in recs_def if "plasmid_0" in name), 1)

            cdir_keep = root / "cohort_keep"
            self.assertEqual(
                cli.main([
                    "generate", "--out", str(cdir_keep), "--seed", "19", "--n-genomes", "2",
                    "--k", "21", "--core-len", "600", "--n-accessory-genes", "1",
                    "--accessory-len", "160", "--n-plasmids", "1", "--plasmid-len", "240",
                    "--plasmid-copies", "4", "--keep-plasmid-copies-in-fasta", "--no-reads",
                ]),
                0,
            )
            m_keep = json.loads((cdir_keep / "manifest.json").read_text(encoding="utf-8"))
            self.assertFalse(m_keep["spec"]["collapse_plasmids_in_assembly"])
            recs_keep = list(parity.iter_sequence_records(cdir_keep / "assemblies" / "syn_g0000.fna"))
            self.assertEqual(sum(1 for name, _ in recs_keep if "plasmid_0" in name), 4)

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

            # F56: verify-parity rejects a .kset built at a different k.
            ref15 = parity.count_file(fna, k=15, min_abundance=1)
            kset_15 = root / "syn_g0000_k15.kset"
            parity.write_kset(kset_15, "syn_g0000", parity.KmerSpace(15), ref15.mixed_sorted)
            with self.assertRaisesRegex(ValueError, "does not match expected k"):
                cli.main([
                    "verify-parity", "--mode", "set", "--k", "21", "--min-abundance", "1",
                    "--reference", str(fna), "--candidate-kset", str(kset_15),
                ])

            # F58: verify-parity --max-cap respected.
            rep_fa = root / "repeat15.fna"
            unit_seq = synthetic.DeterministicRng(33).random_dna(30)
            rep_fa.write_text(f">c0\n{unit_seq}\n" * 10, encoding="ascii")
            ref_rep = parity.count_file(rep_fa, k=15, min_abundance=1)
            space15 = parity.KmerSpace(15)
            capped_dump = root / "capped.tsv"
            capped_dump.write_text(
                "".join(f"{space15.decode(km)}\t{min(cnt, 5)}\n" for km, cnt in sorted(ref_rep.counts.items())),
                encoding="ascii",
            )
            self.assertEqual(
                cli.main([
                    "verify-parity", "--mode", "abundance", "--k", "15", "--min-abundance", "1",
                    "--reference", str(rep_fa), "--candidate-dump", str(capped_dump),
                ]),
                1,
            )
            self.assertEqual(
                cli.main([
                    "verify-parity", "--mode", "abundance", "--k", "15", "--min-abundance", "1",
                    "--max-cap", "5", "--reference", str(rep_fa), "--candidate-dump", str(capped_dump),
                ]),
                0,
            )

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

            # F60: plan propagates --k=21 into every case.
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
            self.assertTrue(all(c["k"] == 21 for c in plan_doc["cases"]))

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
            self.assertTrue(all(r["k"] == 21 for r in dry_doc["results"]))

            # F60: run --mock propagates --k=21.
            self.assertEqual(
                cli.main([
                    "run", "--cohort-dir", str(cdir), "--out-dir", str(root / "mock"),
                    "--workloads", "all", "--tools", "all", "--k", "21", "--mock",
                ]),
                0,
            )
            mock_doc = json.loads((root / "mock" / "bench.json").read_text(encoding="utf-8"))
            self.assertTrue(all(r["k"] == 21 for r in mock_doc["results"]))

            # C3: summarize excludes mock_self_check rows by default, includes with --include-mock.
            sum_real_tsv = root / "sum_real.tsv"
            sum_real_md = root / "sum_real.md"
            self.assertEqual(
                cli.main([
                    "summarize", "--inputs", str(root / "mock" / "bench.json"),
                    "--tsv-out", str(sum_real_tsv), "--md-out", str(sum_real_md),
                ]),
                0,
            )
            self.assertEqual(len(sum_real_tsv.read_text(encoding="utf-8").splitlines()), 1)
            self.assertIn("Excluded", sum_real_md.read_text(encoding="utf-8"))

            sum_all_tsv = root / "sum_all.tsv"
            sum_all_md = root / "sum_all.md"
            self.assertEqual(
                cli.main([
                    "summarize", "--inputs", str(root / "mock" / "bench.json"),
                    "--include-mock", "--include-not-comparable",
                    "--tsv-out", str(sum_all_tsv), "--md-out", str(sum_all_md),
                ]),
                0,
            )
            self.assertIn("| Workload | Tool |", sum_all_md.read_text(encoding="utf-8"))
            self.assertEqual(
                sum_all_tsv.read_text(encoding="utf-8"),
                (root / "mock" / "summary.tsv").read_text(encoding="utf-8"),
            )

    def test_cli_verify_parity_exits_nonzero_on_planted_fault(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            cdir = root / "cohort"
            self._make_cohort(cdir)
            truth = cdir / "truth" / "syn_g0000.assembly.counts.tsv"
            lines = truth.read_text(encoding="utf-8").splitlines()
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

    def test_hygiene_scan_all_id_forms_near_misses_dashes_and_report_gate(self):
        # Covers C6 (all 6+ ID forms, 20+ near misses, em/en dashes in decoded JSON),
        # C8 (runtime-generated non-species synthetic IDs), F43, and F62.
        rng = synthetic.DeterministicRng(909)

        # Generate 10 synthetic IDs for each of 7 forms at runtime (no real taxon or accession literals).
        planted_cases = []
        for i in range(10):
            # 3-7 digit synthetic non-species taxon prefix (9999900..9999999)
            t_long = 9999900 + i
            g_num = 1000 + rng.randrange(89000)
            planted_cases.append(("bvbrc_id_shape", f"entry {t_long}.{g_num} found"))
            # 1-2 digit synthetic taxon prefix
            t_short = 11 + i
            g_num2 = 1000 + rng.randrange(89000)
            planted_cases.append(("bvbrc_id_shape", f"entry {t_short}.{g_num2} found"))
            # GCA_ assembly accession
            acc_a = 900000000 + rng.randrange(99999999)
            planted_cases.append(("assembly_accession", f"asm GCA_{acc_a}.1 here"))
            # GCF_ assembly accession
            acc_f = 900000000 + rng.randrange(99999999)
            planted_cases.append(("assembly_accession", f"asm GCF_{acc_f}.2 here"))
            # Run accession ([SED]RR)
            pfx = ("SRR", "ERR", "DRR")[i % 3]
            planted_cases.append(("run_accession", f"run {pfx}{9000000 + i * 37} here"))
            # BioSample (SAM[NED])
            spfx = ("SAMN", "SAMEA", "SAMD")[i % 3]
            planted_cases.append(("biosample_accession", f"biosample {spfx}{90000000 + i * 53} here"))
            # BioProject (PRJ[NED])
            ppfx = ("PRJNA", "PRJEB", "PRJDB")[i % 3]
            planted_cases.append(("bioproject_accession", f"project {ppfx}{900000 + i * 19} here"))

        for want_tag, sample in planted_cases:
            with self.subTest(tag=want_tag):
                self.assertIn(want_tag, runners.scan_for_leaks(sample))

        # Explicit forbidden_ids catches short numeric tokens and accessions alike.
        short_id = f"{12}.{34}"
        self.assertIn("known_genome_id", runners.scan_for_leaks(f"tag {short_id} here", forbidden_ids={short_id}))

        # 25 near-miss strings across decimals, dates, versions, times, and metric labels -> 0 false alarms.
        near_misses = [
            "clean syn_g0001 report 0.0245",
            "elapsed 0:49.02",
            "elapsed 1:23.45",
            "elapsed 12:34.567",
            "elapsed 1:02:03.456",
            "version 2.2.3",
            "rustc 1.98.1",
            "Python 3.12.3",
            "Ubuntu 24.04",
            "bcalm2 2.2.3",
            "ver 12.345",
            "v100.200",
            "date 2026-11-27",
            "date 2026.11.27",
            "date 15.01.2026",
            "date 2026.11",
            "throughput 208.832 Mbp/s",
            "runtime 123.456 s",
            "size 512.340 MB",
            "size 5.97 B/k-mer",
            "accuracy 99.95%",
            "delta -12.345",
            "wall_seconds: 512.340000",
            "score=123.456",
            "row 1.7 part 6.2 step-100.200 T1.100",
        ]
        for text in near_misses:
            with self.subTest(near_miss=text):
                self.assertEqual(runners.scan_for_leaks(text), [], f"false alarm on {text!r}")

        # Em dash (\u2014) and en dash (\u2013) caught in plain text AND JSON-escaped text (C6).
        self.assertIn("em_dash", runners.scan_for_leaks("text \u2014 dash"))
        self.assertNotIn("en_dash", runners.scan_for_leaks("text \u2014 dash"))
        self.assertIn("en_dash", runners.scan_for_leaks("text \u2013 dash"))
        self.assertIn("em_dash", runners.scan_for_leaks(json.dumps({"note": "a \u2014 b"}, ensure_ascii=True)))
        self.assertIn("en_dash", runners.scan_for_leaks(json.dumps({"note": "a \u2013 b"}, ensure_ascii=True)))

        # Private paths, thread links, box IDs.
        self.assertIn("private_home_path", runners.scan_for_leaks("/" + "home/alice/work"))
        self.assertEqual(runners.scan_for_leaks("/" + "home/user/work"), [])
        self.assertIn("windows_user_path", runners.scan_for_leaks("C:\\" + "Users\\alice\\work"))
        self.assertIn("thread_link", runners.scan_for_leaks("t3" + "-thread://" + "abc"))
        self.assertIn("box_id", runners.scan_for_leaks("bx" + "_a1b2c3d4"))

        # F43 & C6: build_bench_report AND write_bench_outputs reject em dash, en dash, and leaked IDs
        # inside JSON strings, while allowing 100+ second float timings in JSON and summary.tsv.
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            cdir = root / "cohort"
            manifest = self._make_cohort(cdir, n_genomes=2, generate_reads=False)
            plan = runners.build_command_plan("grm", "T1", cdir, root / "work", k=21, cohort_manifest=manifest)
            clean_rec = runners.execute_case(plan, manifest, cdir, mode="mock")
            clean_rec["timing"]["user_seconds"] = 512.340123
            clean_rec["timing"]["wall_seconds"] = 123.456789
            rep_ok = runners.build_bench_report([clean_rec], manifest, strict_hygiene=True)
            runners.write_bench_outputs(root / "clean_out", rep_ok, strict_hygiene=True)

            syn_id = f"{9999991}.{12345}"
            for label, bad_note in (
                ("em_dash", "note with \u2014 dash"),
                ("en_dash", "note with \u2013 dash"),
                ("bvbrc_id_shape", f"leaked {syn_id}"),
                ("assembly_accession", f"leaked GCA_{901234567}.1"),
            ):
                with self.subTest(poison=label):
                    bad_rec = dict(clean_rec)
                    bad_rec["notes"] = bad_note
                    with self.assertRaisesRegex(ValueError, label):
                        runners.build_bench_report([bad_rec], manifest, strict_hygiene=True)
                    unverified = runners.build_bench_report([bad_rec], manifest, strict_hygiene=False)
                    with self.assertRaisesRegex(ValueError, label):
                        runners.write_bench_outputs(root / f"bad_{label}", unverified, strict_hygiene=True)


if __name__ == "__main__":
    unittest.main()
