"""Unit, determinism, and planted-fault tests for scripts/kmer-bench/synthetic.py."""

import hashlib
import os
import pathlib
import statistics
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

    def test_mutate_dna_always_mutates_to_different_nucleotide(self):
        # Directly catches F46: at snp_rate=1.0, EVERY base must change to a different nucleotide.
        for base in "ACGT":
            rng = synthetic.DeterministicRng(10 + "ACGT".index(base))
            mutated = synthetic.mutate_dna(base * 1000, snp_rate=1.0, rng=rng)
            self.assertEqual(len(mutated), 1000)
            self.assertNotIn(base, mutated)
            self.assertEqual(set(mutated), set("ACGT") - {base})

    def test_read_errors_mutate_base_and_r2_is_reverse_complement(self):
        # Directly catches F47 (every injected Q2 '#' read error changes the base) and
        # F49 (R2 is reverse complement of fragment end, not same strand as R1).
        rng_err = synthetic.DeterministicRng(107)
        r1_err, r2_err, n_err = synthetic._simulate_paired_reads(
            [("homopol", "A" * 300, 1)],
            coverage=20.0,
            read_length=60,
            insert_size=120,
            error_rate=0.20,
            rng=rng_err,
        )
        self.assertGreater(n_err, 50)
        checked_errors = 0
        for (_h1, s1, q1), (_h2, s2, q2) in zip(r1_err, r2_err):
            # On an all-A template, every clean ('I') base in R1 is 'A' (forward) or 'T' (reverse),
            # and R2 is on the opposite strand.
            clean_r1 = [ch for ch, q in zip(s1, q1) if q == "I"]
            self.assertTrue(len(clean_r1) > 0)
            r1_expected_base = clean_r1[0]
            self.assertIn(r1_expected_base, ("A", "T"))
            self.assertTrue(all(ch == r1_expected_base for ch in clean_r1))
            r2_expected_base = "T" if r1_expected_base == "A" else "A"
            for ch, q in zip(s1, q1):
                if q == "#":
                    self.assertNotEqual(ch, r1_expected_base)
                    checked_errors += 1
            for ch, q in zip(s2, q2):
                if q == "#":
                    self.assertNotEqual(ch, r2_expected_base)
                    checked_errors += 1
                else:
                    self.assertEqual(ch, r2_expected_base)
        self.assertEqual(checked_errors, n_err)

        # Clean reads on a random non-palindromic template: R2 must be revcomp of fragment end (F49).
        rng_rc = synthetic.DeterministicRng(211)
        tmpl = rng_rc.random_dna(500)
        tmpl_rc = synthetic.revcomp(tmpl)
        r1_clean, r2_clean, zero_err = synthetic._simulate_paired_reads(
            [("c0", tmpl, 1)],
            coverage=12.0,
            read_length=60,
            insert_size=120,
            error_rate=0.0,
            rng=rng_rc,
        )
        self.assertEqual(zero_err, 0)
        self.assertGreater(len(r1_clean), 20)
        for (_h1, s1, _q1), (_h2, s2, _q2) in zip(r1_clean, r2_clean):
            s2_rc = synthetic.revcomp(s2)
            if s1 in tmpl:
                pos = tmpl.index(s1)
                self.assertEqual(s2_rc, tmpl[pos + 60 : pos + 120])
                self.assertNotIn(s2, tmpl)
            else:
                pos = tmpl_rc.index(s1)
                self.assertEqual(s2_rc, tmpl_rc[pos + 60 : pos + 120])
                self.assertNotIn(s2, tmpl_rc)

    def test_assembly_features_contigs_n_runs_lowercase_and_collapsed_plasmids(self):
        # Explicitly inspects sequence bytes (not FASTA headers!) to catch F50.
        spec = self._small_spec(collapse_plasmids_in_assembly=True, include_lowercase=True)
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            manifest = synthetic.generate_cohort(root, spec)
            g0 = manifest["genomes"][0]
            self.assertEqual(g0["plasmids_true_copy_number"], {"plasmid_0": 4, "plasmid_1": 8})
            records = list(parity.iter_sequence_records(root / g0["fasta"]))
            self.assertEqual(len(records), g0["contigs"])
            seq_by_name = {name: seq.decode("ascii") for name, seq in records}
            self.assertIn("NNNNN", seq_by_name["syn_g0000_core_0"])
            # Last core contig MUST have a 50-bp lowercase stretch in its sequence bytes (F50).
            last_core = seq_by_name["syn_g0000_core_2"]
            self.assertTrue(all(ch in "acgt" for ch in last_core[60:110]))
            self.assertEqual(sum(1 for name, _ in records if "plasmid_0" in name), 1)
            self.assertLess(g0["assembly_truth"]["distinct_kmers_min2"], 25)

        # With include_lowercase=False, ZERO sequence bases may be lowercase.
        spec_no_lc = self._small_spec(n_genomes=2, include_lowercase=False, generate_reads=False)
        with tempfile.TemporaryDirectory() as tmp2:
            root2 = pathlib.Path(tmp2)
            m2 = synthetic.generate_cohort(root2, spec_no_lc)
            for _name, seq in parity.iter_sequence_records(root2 / m2["genomes"][0]["fasta"]):
                self.assertFalse(any(ch in b"acgt" for ch in seq))

    def test_plasmid_carriers_chrom_median_depth_and_kmtricks_fof(self):
        # Directly catches F51 (exact plasmid carrier indices), F52 (median vs mean chrom depth),
        # and F53 (reads.kmtricks.fof includes both R1 and R2).
        spec = self._small_spec(
            seed=42,
            n_genomes=6,
            core_len=1800,
            plasmid_len=450,
            plasmid_copy_numbers=(4, 8),
            read_coverage=14.0,
            error_rate=0.006,
        )
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            manifest = synthetic.generate_cohort(root, spec)
            genomes = manifest["genomes"]

            # F51: plasmid_0 carried by g % 2 == 0 -> [0, 2, 4]; plasmid_1 by g % 3 == 0 -> [0, 3].
            p0_carriers = [i for i, g in enumerate(genomes) if "plasmid_0" in g["plasmids_true_copy_number"]]
            p1_carriers = [i for i, g in enumerate(genomes) if "plasmid_1" in g["plasmids_true_copy_number"]]
            self.assertEqual(p0_carriers, [0, 2, 4])
            self.assertEqual(p1_carriers, [0, 3])
            self.assertEqual(genomes[1]["plasmids_true_copy_number"], {})
            self.assertEqual(genomes[5]["plasmids_true_copy_number"], {})

            # F52: chrom_median_depth matches exact median (and differs from mean).
            g0 = genomes[0]
            space = parity.KmerSpace(spec.k)
            rd_counts = parity.count_files(
                [root / g0["reads_r1"], root / g0["reads_r2"]], k=spec.k, min_abundance=1
            ).counts
            # Recompute chromosome k-mers from core contigs (note N-run in core_0 was stripped from read template,
            # so using core contigs gives >99% of chromosome k-mers; also verify via exact generator logic).
            core_kmers = set()
            for cname, seq in parity.iter_sequence_records(root / g0["fasta"]):
                if "_core_" in cname:
                    core_kmers.update(space.iter_canonical(seq))
            depths = [rd_counts[km] for km in core_kmers if km in rd_counts]
            exact_med = round(float(statistics.median(depths)), 4)
            exact_mean = round(float(statistics.mean(depths)), 4)
            self.assertNotEqual(exact_med, exact_mean)
            self.assertEqual(g0["reads_truth"]["chrom_median_depth"], exact_med)
            self.assertNotEqual(g0["reads_truth"]["chrom_median_depth"], exact_mean)

            # F53: reads.kmtricks.fof MUST include both R1 and R2 for every genome.
            fof_lines = (root / "reads.kmtricks.fof").read_text(encoding="ascii").splitlines()
            self.assertEqual(len(fof_lines), len(genomes))
            for line, rec in zip(fof_lines, genomes):
                self.assertEqual(line, f"{rec['genome']} : {rec['reads_r1']} ; {rec['reads_r2']}")
                self.assertIn("_R1.fastq", line)
                self.assertIn("_R2.fastq", line)

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
            self.assertGreater(rt["distinct_kmers_min1"], rt["distinct_kmers_min2"])
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
