"""pilot_interrupt.py against the real GRM binary (Linux).

usage: GRM=target/release/grm python3 -I scripts/mtb-4k/tests/test_pilot_interrupt.py

A disk-sparse run at width 64 over 128 genomes is held after block 0 folds:
genome 64's FASTA is a named pipe nobody writes, so the first read of block 1
blocks. A wrapper standing in for the resource guard runs GRM, so its command
line also contains `GRM --data ROOT matrix`, as the guard's and
/usr/bin/time's do. The helper must stop and kill GRM only, at
folded_columns = 64, and the run must then resume to the matrix an
uninterrupted run gives. It must refuse when nothing runs, and when the run
ends before the kill point.
"""
import json
import os
import random
import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
HELPER = HERE.parent / "pilot_interrupt.py"
sys.path.insert(0, str(HERE.parent))
import pilot_interrupt as pi  # noqa: E402

SLUG = "mycobacterium_tuberculosis"
FILES = ["matrix.json", "patterns.bin", "support.bin", "kmer_offsets.bin", "kmers.bin"]
WRAPPER = "import subprocess, sys; r = subprocess.run(sys.argv[2:]); open(sys.argv[1], 'w').write(str(r.returncode))"


def grm():
    path = os.environ.get("GRM")
    if not path or not Path(path).is_file():
        raise SystemExit("set GRM to a built grm binary (cargo build --release)")
    return str(Path(path).resolve())


def build_root(base: Path, n=128, seed=3):
    rng = random.Random(seed)
    root = base / "data"
    sp = root / "species" / SLUG
    sp.mkdir(parents=True)
    (root / "fasta").mkdir()
    core = "".join(rng.choice("ACGT") for _ in range(600))
    ids = [f"G{i:03d}" for i in range(n)]
    texts = {}
    for i, g in enumerate(ids):
        seq = list(core)
        for pos in rng.sample(range(len(core)), 8):
            seq[pos] = rng.choice("ACGT")
        texts[g] = (">c1\n" + "".join(seq) + "\n>c2\n" + "".join(rng.choice("ACGT") for _ in range(40 + i % 7)) + "\n").encode()
        (root / "fasta" / f"{g}.fna.gz").write_bytes(texts[g])
    (sp / "genomes.txt").write_text("\n".join(ids) + "\n")
    return root, ids, texts


def sparse(root):
    return [grm(), "--data", str(root), "matrix", "--species", SLUG, "--min-support", "2", "--disk-sparse",
            "--block-width", "64", "--allow-missing-catalog"]


def helper(root, receipt, folded=64, timeout=120):
    return subprocess.run([sys.executable, "-I", str(HELPER), "--data", str(root), "--grm", grm(), "--folded", str(folded),
                           "--receipt", str(receipt), "--timeout", str(timeout), "--poll", "0.05"],
                          capture_output=True, text=True)


def matrix_bytes(root):
    return {f: (root / "species" / SLUG / "matrix" / f).read_bytes() for f in FILES}


class Interrupt(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.base = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def gated(self, name):
        """A root whose genome 64 is a pipe, and its text."""
        root, ids, texts = build_root(self.base / name)
        fifo = root / "fasta" / f"{ids[64]}.fna.gz"
        fifo.unlink()
        os.mkfifo(fifo)
        return root, fifo, texts[ids[64]]

    def test_find_run_matches_only_the_exact_command_line(self):
        g, d = "/x/grm", "/run/cohorts/pilot256/data"
        table = {
            "11": ["python3", "/guard.py", "--log", "l", "--", "/usr/bin/time", "-v", g, "--data", d, "matrix"],
            "12": ["/usr/bin/time", "-v", g, "--data", d, "matrix", "--species", SLUG],
            "13": [g, "--data", d, "matrix", "--species", SLUG, "--disk-sparse"],
            "14": [g, "--data", d + "-stock", "matrix"],
            "15": [g, "--data", d, "structure"],
            "16": None,
        }
        self.assertEqual(pi.find_run(g, d, pids=list(table), read=table.get), [13])
        table["17"] = list(table["13"])
        self.assertEqual(pi.find_run(g, d, pids=list(table), read=table.get), [13, 17])

    def test_the_run_is_stopped_after_block_zero_folds_and_only_grm_is_killed(self):
        reference, _, _ = build_root(self.base / "reference")
        subprocess.run(sparse(reference), check=True, capture_output=True)
        root, fifo, text = self.gated("pilot")
        rc_file = self.base / "wrapper.rc"
        wrapper = subprocess.Popen([sys.executable, "-c", WRAPPER, str(rc_file), *sparse(root)],
                                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        receipt = self.base / "interrupt.json"
        p = helper(root, receipt)
        self.assertEqual(p.returncode, 0, p.stderr)
        self.assertEqual(wrapper.wait(timeout=60), 0)  # the guard stand-in was never signalled
        self.assertEqual(rc_file.read_text(), "-9")      # and saw its child killed
        got = json.loads(receipt.read_text())
        self.assertTrue(got["interrupted"] and got["gone"])
        state = got["stopped_state"]
        self.assertEqual((state["folded_columns"], state["saved"], state["block_width"]), (64, False, 64))
        self.assertIn("block-00000.partial", state["run_files"])
        self.assertIn("block-00000.partial.receipt.json", state["run_files"])
        self.assertEqual(got["cmdline"][:4], [grm(), "--data", str(root), "matrix"])
        # A receipt is never overwritten.
        before = receipt.read_bytes()
        self.assertNotEqual(helper(root, receipt, timeout=1).returncode, 0)
        self.assertEqual(receipt.read_bytes(), before)
        # The resume, with the real FASTA in place, gives the uninterrupted matrix.
        real = fifo.with_name(fifo.name + ".real")
        real.write_bytes(text)
        os.replace(real, fifo)
        subprocess.run(sparse(root), check=True, capture_output=True)
        self.assertEqual(matrix_bytes(root), matrix_bytes(reference))
        self.assertFalse((root / "species" / SLUG / "matrix-partials").exists())

    def test_nothing_running_is_a_stop(self):
        root, _, _ = build_root(self.base / "idle")
        p = helper(root, self.base / "idle.json", timeout=1)
        self.assertEqual(p.returncode, 1)
        self.assertIn("found []", json.loads((self.base / "idle.json").read_text())["error"])

    def test_a_run_that_ends_before_the_kill_point_is_a_stop(self):
        root, fifo, text = self.gated("early")
        run = subprocess.Popen(sparse(root), stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        receipt = self.base / "early.json"
        watcher = subprocess.Popen([sys.executable, "-I", str(HELPER), "--data", str(root), "--grm", grm(),
                                    "--folded", "4096", "--receipt", str(receipt), "--timeout", "120", "--poll", "0.05"],
                                   stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        # Let block 0 fold, then release the pipe: hand it its bytes once and
        # put the real file at its name for every later read.
        deadline = time.time() + 60
        while time.time() < deadline:
            _, state = pi.run_state(root / "species" / SLUG)
            if state and state.get("folded_columns") == 64:
                break
            time.sleep(0.05)
        fd = os.open(fifo, os.O_WRONLY)
        real = fifo.with_name(fifo.name + ".real")
        real.write_bytes(text)
        os.replace(real, fifo)
        os.write(fd, text)
        os.close(fd)
        self.assertEqual(run.wait(timeout=120), 0)
        _, err = watcher.communicate(timeout=120)
        self.assertEqual(watcher.returncode, 1, err)
        self.assertIn("ended before the kill point", json.loads(receipt.read_text())["error"])


if __name__ == "__main__":
    unittest.main()
