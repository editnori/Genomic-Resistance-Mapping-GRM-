"""scripts/c2-cap/fetch_check.py, scrub.py, timing.py and run_species.sh's text (docs/cap-arms.md, The fetch and Files).

Every genome, FASTA, list and log here is SYNTHETIC.
"""

import gzip
import hashlib
import io
import json
import re
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

from tests.c2_cap import common

fetch = common.load("fetch_check")
scrub = common.load("scrub")
timing = common.load("timing")

RUNNER = common.KIT / "run_species.sh"


def fasta(path, contigs, header="contig"):
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="ascii") as f:
        for i, n in enumerate(contigs):
            f.write(f">{header}{i}\n{'ACGT' * (n // 4)}\n")


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


class FetchCheck(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        t = self.t = Path(self.tmp.name)
        self.release = [f"5{i:03d}.1" for i in range(6)]
        self.added = [f"5{i:03d}.2" for i in range(4)]
        sizes = {g: [400, 200, 80] for g in self.release + self.added}
        for g in self.release + self.added:
            fasta(t / "fasta" / f"{g}.fna.gz", sizes[g])
        common.write_json(t / "release.json", {"genomes": [{"id": g, "length": 680, "contigs": 3, "n50": 400}
                                                           for g in self.release]})
        common.write_ids(t / "cohorts" / "cap-2000.txt", self.release)
        common.write_ids(t / "cohorts" / "cap-4000.txt", self.release + self.added)
        common.write_ids(t / "cohorts" / "cap-all.txt", self.release + self.added)
        common.write_json(t / "cohorts" / "cohorts.json", {"arms": {
            "cap-2000": {"genomes": 6, "added": 0}, "cap-4000": {"genomes": 10, "added": 4},
            "cap-all": {"genomes": 10, "added": 4}}})
        # The master list holds the release's six and one added genome.
        listed = self.release + self.added[:1]
        (t / "master.sha256").write_text("".join(f"{sha(t / 'fasta' / f'{g}.fna.gz')}  fasta/{g}.fna.gz\n"
                                                 for g in listed), encoding="utf-8")
        common.meta(t / "meta.jsonl.gz", {g: "PRJ1" for g in listed})

    def tearDown(self):
        self.tmp.cleanup()

    def run_check(self):
        t = self.t
        rc, text = common.quiet(fetch.main, ["fetch_check.py", str(t / "release.json"), str(t / "cohorts"),
                                             str(t / "fasta"), str(t / "master.sha256"), str(t / "meta.jsonl.gz"),
                                             str(t / "fasta_failures.tsv"), str(t / "fetch.json"),
                                             str(t / "fetched.sha256")])
        return rc, json.loads((t / "fetch.json").read_text(encoding="utf-8"))

    def test_a_clean_fetch_passes_and_records_the_added_genomes_sha256(self):
        rc, doc = self.run_check()
        self.assertEqual(rc, 0)
        self.assertEqual((doc["genomes"], doc["release"], doc["added"], doc["failed_release"], doc["failed_added"]),
                         (10, 6, 4, 0, 0))
        self.assertEqual((doc["added_in_master"], doc["added_sha256_recorded"]), (1, 3))
        self.assertEqual(len((self.t / "fetched.sha256").read_text(encoding="utf-8").splitlines()), 3)

    def test_release_bytes_that_differ_with_the_same_sequence_are_counted_not_failed(self):
        fasta(self.t / "fasta" / f"{self.release[2]}.fna.gz", [400, 200, 80], header="renamed")
        rc, doc = self.run_check()
        self.assertEqual(rc, 0)
        self.assertEqual((doc["release_bytes_differ_from_master"], doc["release_stats_differ_from_release"]), (1, 0))

    def test_a_release_genome_whose_length_contigs_or_n50_differ_fails(self):
        fasta(self.t / "fasta" / f"{self.release[3]}.fna.gz", [400, 200, 40, 40])
        rc, doc = self.run_check()
        self.assertEqual(rc, 1)
        self.assertEqual(doc["release_stats_differ_from_release"], 1)

    def test_a_failed_release_genome_fails(self):
        (self.t / "fasta" / f"{self.release[0]}.fna.gz").unlink()
        rc, doc = self.run_check()
        self.assertEqual((rc, doc["failed_release"]), (1, 1))

    def test_a_failed_added_genome_leaves_every_wider_cohort(self):
        (self.t / "fasta" / f"{self.added[2]}.fna.gz").unlink()
        (self.t / "fasta_failures.tsv").write_text(f"{self.added[3]}\tHTTP 500\n", encoding="utf-8")
        rc, doc = self.run_check()
        self.assertEqual(rc, 0)
        self.assertEqual(doc["failed_added"], 2)
        self.assertEqual(doc["left_wider_cohorts"], {"cap-4000": 2, "cap-all": 2})
        for arm in ("cap-4000", "cap-all"):
            ids = (self.t / "cohorts" / f"{arm}.txt").read_text(encoding="utf-8").split()
            self.assertEqual(ids, self.release + self.added[:2])
        self.assertEqual((self.t / "cohorts" / "cap-2000.txt").read_text(encoding="utf-8").split(), self.release)
        arms = json.loads((self.t / "cohorts" / "cohorts.json").read_text(encoding="utf-8"))["arms"]
        self.assertEqual((arms["cap-all"]["genomes"], arms["cap-all"]["added"], arms["cap-all"]["fetch_failed"]),
                         (8, 2, 2))

    def test_n50_is_the_contig_at_half_the_length(self):
        p = self.t / "x.fna.gz"
        fasta(p, [100, 300, 200, 400])
        self.assertEqual(fetch.stats(p), (1000, 4, 300))


class Scrub(unittest.TestCase):
    def test_known_genome_ids_are_written_as_genome_id_and_counted_by_the_scan(self):
        with tempfile.TemporaryDirectory() as tmp:
            t = Path(tmp)
            common.write_json(t / "release.json", {"genomes": [{"id": "1352.1001"}]})
            common.write_ids(t / "cohorts" / "cap-all.txt", ["1352.1001", "1352.2002"])
            known = scrub.ids(t / "release.json", t / "cohorts")
            out = io.BytesIO()
            scrub.scrub([b"genome 1352.2002 and 1352.1001, k 31, took 12.5 s, 41352.2002x\n"], out, known)
            self.assertEqual(out.getvalue(), b"genome [genome id] and [genome id], k 31, took 12.5 s, 41352.2002x\n")
            (t / "logs").mkdir()
            (t / "logs" / "a.log").write_text("clean 12.5 s\n", encoding="utf-8")
            self.assertEqual(scrub.scan([t / "logs"], known), 0)
            (t / "logs" / "b.tsv").write_text("x\t1352.2002\n", encoding="utf-8")
            self.assertEqual(scrub.scan([t / "logs"], known), 1)


class Timing(unittest.TestCase):
    def test_each_time_file_gives_wall_cpu_and_peak(self):
        with tempfile.TemporaryDirectory() as tmp:
            t = Path(tmp)
            (t / "fit-a.time").write_text("Command exited with non-zero status 1\n"
                                          "fit-a 2.89 s wall, 1466544 KB peak, 4.33 s user, 1.07 s sys\n",
                                          encoding="utf-8")
            (t / "count.time").write_text("count 129.53 s wall, 1668636 KB peak, 1485.74 s user, 227.27 s sys\n",
                                          encoding="utf-8")
            (t / "empty.time").write_text("", encoding="utf-8")
            self.assertEqual(timing.rows(t), [("count", "129.53", "1713.01", "1668636"), ("empty", "-", "-", "-"),
                                              ("fit-a", "2.89", "5.40", "1466544")])


class Runner(unittest.TestCase):
    def setUp(self):
        self.text = RUNNER.read_text(encoding="utf-8")

    @unittest.skipUnless(shutil.which("bash"), "no bash")
    def test_the_runner_parses(self):
        subprocess.run(["bash", "-n", str(RUNNER)], check=True)

    def test_no_bash_special_variable_is_set_or_read(self):
        # bash's GROUPS holds the user's group ids; $GROUPS reads 1000, not a path.
        for name in ("GROUPS", "RANDOM", "SECONDS", "HOSTNAME", "UID", "PPID", "LINENO", "BASHPID", "PWD"):
            self.assertIsNone(re.search(rf"(?<![\w$]){name}=|\${{?{name}\b", self.text), name)

    def test_scores_go_under_score_and_the_end_scan_covers_results_score_and_logs(self):
        self.assertIn('"$OUT/score" > "$OUT/score/score.log"', self.text)
        self.assertNotRegex(self.text, r"score\.py[^\n]*\$OUT/logs")
        self.assertIn("scrub.py\" --scan \"$R/species/$SP.json\" \"$COH\" results score logs", self.text)

    def test_with_score_no_the_score_files_move_where_nothing_reads_them(self):
        self.assertIn("for f in predictions.jsonl runs.json repeats.json; do", self.text)
        self.assertIn('mv -f "$FITS/$1/$2/$f" "$UNREAD/$1/$2/$f"', self.text)

    def test_every_fit_checks_its_folds_before_anything_is_scored(self):
        self.assertLess(self.text.index('check_folds.py" lineage'), self.text.index('"$HERE/score.py"'))
        self.assertLess(self.text.index('check_folds.py" study'), self.text.index('"$HERE/score.py"'))
        self.assertIn("--pairing common_outer_test_set_balanced --reference \"$FITS/cap-2000/$d\"", self.text)

    def test_cap_4000_is_withdrawn_from_the_runs_default_arms(self):
        self.assertIn("ARMS=${8:-cap-all}", self.text)
        self.assertIn("grep -Eq '^(cap-4000|cap-all)(,(cap-4000|cap-all))?$'", self.text)

    def test_lineage_chains_run_at_once_under_the_memory_guard_and_at_half_after_a_memory_kill(self):
        t = self.text
        self.assertIn('chain "$1"; shift', t)
        self.assertLess(t.index('chain "$1"; shift'), t.index('chain "$d" &'))
        self.assertIn('while [ "$(jobs -rp | wc -l)" -ge "$(at_once)" ]; do wait -n || failed=1; done', t)
        self.assertIn("n=$(( (total - 4194304) / peak ))", t)
        self.assertIn('[ "$n" -le "$LIMIT" ] || n=$LIMIT', t)
        self.assertIn('touch "$OUT/flags/memory-killed"', t)
        self.assertIn("LIMIT=$(( LIMIT / 2 )); [ \"$LIMIT\" -ge 1 ] || LIMIT=1", t)
        self.assertIn("echo \"$AT_ONCE\" | grep -Eq '^[1-4]$'", t)
        study = t[t.index('fit "study-cap-2000-$d"'):t.index("STEP=timing")]
        self.assertNotRegex(study, r"&\s*$|& *\n")

    def test_a_none_added_arm_is_not_fitted_and_its_check_reads_the_cohort_count(self):
        t = self.text
        self.assertLess(t.index('if none_added "$arm" "$d"; then'), t.index('--arm "$arm" --blocks K'))
        self.assertIn('"$OUT/results/folds-lineage-$d.tsv" "$OUT/results/added-labelled.tsv"', t)
        self.assertIn('"$OUT/results/roots.json" \\\n  "$OUT/results/added-labelled.tsv" > "$OUT/expect.tsv.part"', t)

    def test_the_bioproject_answers_stay_in_work_and_only_counts_reach_results(self):
        t = self.text
        self.assertIn('"$HERE/bioprojects.py" "$COH" "$KIT/master-genome-meta.jsonl.gz" "$W/bioprojects.part.jsonl.gz"', t)
        self.assertIn('"$OUT/results/bioprojects.json"', t)
        self.assertIn("META=$META,$W/bioprojects.jsonl.gz", t)
        self.assertLess(t.index("step bioprojects bioprojects"), t.index("step plans plans"))


if __name__ == "__main__":
    unittest.main()
