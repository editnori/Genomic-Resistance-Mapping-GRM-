"""Unit and determinism tests for scripts/kmer-bench/synthetic.py."""

import hashlib
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


def _tree_hashes(root: pathlib.Path):
    out = {}
    for p in sorted(root.rglob("*")):
        if p.is_file():
            rel = p.relative_to(root).as_posix()
            out[rel] = hashlib.sha256(p.read_bytes()).hexdigest()
    return out


class SyntheticCohortTests(unittest.TestCase):
    def _small_spec(self, **overrides):
        params = dict(
            seed=42,
            n_genomes=6,
            k=21,
            min_support=2,
            core_len=1500,
            n_core_contigs=3,
            n_accessory_genes=4,
            accessory_len=300,
            snp_rate=0.006,
            n_plasmids=2,
            plasmid_len=450,
            plasmid_copy_numbers=(4, 8),
            read_coverage=14.0,
            read_length=80,
            insert_size=140,
            error_rate=0.005,
        )
        params.update(overrides)
        return synthetic.CohortSpec(**params)

    def test_determinism_in_process_across_repeated_calls(self):
        spec = self._small_spec()
        with tempfile.TemporaryDirectory() as t1, tempfile.TemporaryDirectory() as t2:
            m1 = synthetic.generate_cohort(pathlib.Path(t1), spec)
            m2 = synthetic.generate_cohort(pathlib.Path(t2), spec)
            self.assertEqual(m1, m2)
            self.assertEqual(_tree_hashes(pathlib.Path(t1)), _tree_hashes(pathlib.Path(t2)))

    def test_determinism_across_pythonhashseed_0_1_4242(self):
        cli_path = common.SCRIPT_DIR / "cli.py"
        snapshots = {}
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            for seed in ("0", "1", "4242"):
                out_dir = root / f"h{seed}"
                env = {
                    "PATH": os.environ.get("PATH", "/usr/bin:/bin"),
                    "PYTHONHASHSEED": seed,
                }
                proc = subprocess.run(
                    [
                        sys.executable, "-s", "-P", "-B", "-W", "error", str(cli_path),
                        "generate", "--out", str(out_dir),
                        "--seed", "101", "--n-genomes", "4", "--k", "21",
                        "--core-len", "900", "--n-accessory-genes", "3", "--accessory-len", "200",
                        "--n-plasmids", "1", "--plasmid-len", "300", "--plasmid-copies", "5",
                        "--coverage", "10", "--read-length", "75", "--insert-size", "130",
                    ],
                    env=env,
                    capture_output=True,
                    text=True,
                    check=False,
                )
                self.assertEqual(proc.returncode, 0, proc.stderr)
                snapshots[seed] = _tree_hashes(out_dir)
            self.assertEqual(snapshots["0"], snapshots["1"])
            self.assertEqual(snapshots["0"], snapshots["4242"])

    def test_assembly_features_contigs_n_runs_lowercase_and_collapsed_plasmids(self):
        spec = self._small_spec(collapse_plasmids_in_assembly=True)
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            manifest = synthetic.generate_cohort(root, spec)
            g0 = manifest["genomes"][0]
            self.assertEqual(g0["plasmids_true_copy_number"], {"plasmid_0": 4, "plasmid_1": 8})
            fna_text = (root / g0["fasta"]).read_text(encoding="utf-8")
            self.assertIn("NNNNN", fna_text)
            self.assertTrue(any(ch in "acgt" for ch in fna_text))
            records = list(parity.iter_sequence_records(root / g0["fasta"]))
            self.assertEqual(len(records), g0["contigs"])
            # Collapsed assembly stores each plasmid as 1 contig, so assembly [>=2] k-mers are near 0.
            self.assertEqual(
                sum(1 for name, _ in records if "plasmid_0" in name),
                1,
            )
            self.assertLess(g0["assembly_truth"]["distinct_kmers_min2"], 25)

    def test_noisy_paired_fastq_recovers_plasmid_copy_number_ratio(self):
        spec = self._small_spec(
            n_genomes=4,
            core_len=2400,
            plasmid_len=600,
            plasmid_copy_numbers=(4, 9),
            read_coverage=16.0,
            error_rate=0.008,
        )
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            manifest = synthetic.generate_cohort(root, spec)
            g0 = manifest["genomes"][0]
            rt = g0["reads_truth"]
            self.assertGreater(rt["injected_base_errors"], 0)
            # Sequencing errors create many singleton k-mers that min_abundance=2 removes.
            self.assertGreater(rt["distinct_kmers_min1"], rt["distinct_kmers_min2"])
            # Observed plasmid-to-chromosome read depth ratios track true copy numbers (4x and 9x).
            r0 = rt["plasmid_observed_depth_ratio"]["plasmid_0"]
            r1 = rt["plasmid_observed_depth_ratio"]["plasmid_1"]
            self.assertGreater(r0, 2.5)
            self.assertLess(r0, 5.5)
            self.assertGreater(r1, 6.5)
            self.assertLess(r1, 11.5)
            self.assertGreater(r1, r0)

    def test_manifest_ground_truth_matches_independent_recount(self):
        spec = self._small_spec(n_genomes=4, core_len=1200)
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            manifest = synthetic.generate_cohort(root, spec)
            for rec in manifest["genomes"]:
                asm = parity.count_file(root / rec["fasta"], k=spec.k, min_abundance=1)
                self.assertEqual(asm.fnv1a64_digest, rec["assembly_truth"]["fnv1a64_min1"])
                self.assertEqual(asm.sha256_canonical, rec["assembly_truth"]["sha256_canonical_min1"])
                self.assertEqual(asm.sha256_spectrum, rec["assembly_truth"]["sha256_spectrum_min1"])
                # Truth TSV table matches exact counts.
                tsv_counts = parity.parse_kmer_text_dump(
                    root / "truth" / f"{rec['genome']}.assembly.counts.tsv", k=spec.k, min_abundance=1
                )
                self.assertTrue(parity.verify_abundance_parity(asm.counts, tsv_counts, k=spec.k).ok)
                rd = parity.count_files(
                    [root / rec["reads_r1"], root / rec["reads_r2"]], k=spec.k, min_abundance=2
                )
                self.assertEqual(rd.fnv1a64_digest, rec["reads_truth"]["fnv1a64_min2"])
                self.assertEqual(rd.sha256_spectrum, rec["reads_truth"]["sha256_spectrum_min2"])

    def test_gzipped_fastq_reads_match_uncompressed_truth(self):
        spec_plain = self._small_spec(n_genomes=2, core_len=900, gzip_reads=False)
        spec_gz = self._small_spec(n_genomes=2, core_len=900, gzip_reads=True)
        with tempfile.TemporaryDirectory() as t1, tempfile.TemporaryDirectory() as t2:
            m1 = synthetic.generate_cohort(pathlib.Path(t1), spec_plain)
            m2 = synthetic.generate_cohort(pathlib.Path(t2), spec_gz)
            self.assertEqual(m1["genomes"][0]["reads_truth"], m2["genomes"][0]["reads_truth"])
            self.assertTrue((pathlib.Path(t2) / m2["genomes"][0]["reads_r1"]).name.endswith(".fastq.gz"))

    def test_generated_cohort_contains_zero_genome_ids_or_private_leaks(self):
        spec = self._small_spec(n_genomes=3, core_len=900)
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            synthetic.generate_cohort(root, spec)
            for p in sorted(root.rglob("*")):
                if p.is_file():
                    text = p.read_text(encoding="utf-8")
                    leaks = runners.scan_for_leaks(text)
                    self.assertEqual(leaks, [], f"{p.name}: {leaks}")


if __name__ == "__main__":
    unittest.main()
