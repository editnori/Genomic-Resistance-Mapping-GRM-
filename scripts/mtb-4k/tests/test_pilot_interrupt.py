"""pilot_interrupt.py against the real GRM binary (Linux).

usage: GRM=target/release/grm python3 -I scripts/mtb-4k/tests/test_pilot_interrupt.py

A disk-sparse run at width 64 over 128 genomes is held after block 0 folds:
genome 64's FASTA is a named pipe nobody writes, so the first read of block 1
blocks. A stand-in for the resource guard is the launcher, a file named as the
guard and run as the plan runs it (`python3 .../resource_guard.py --log L
--resources R -- COMMAND`): it runs COMMAND as its child (so its command line
also contains `GRM --data ROOT matrix`, as the guard's and /usr/bin/time's
do), and when the child ends it writes a record in the guard's form
({"command", "rc", "guard_reason"}). The helper must stop and kill GRM only,
at folded_columns = 64, verify the kill in that record, and the run must then
resume to the matrix an uninterrupted run gives.

L1 (plan-validators review): the helper binds only to a process below the
launcher. A stale exact-argv run outside it (a real GRM blocked on another
species of the same root) is never bound or signalled and is a stop. On every
stop (nothing found, timeout, the run ending first, a foreign match, a bad
guard record, an existing receipt) the launched run is stopped too, and a
--launcher that is not a sibling job is refused with nothing signalled.

The focused re-check: a sibling that is not this run's guard is refused with
nothing signalled (N1); every stop lists the matching runs still going and
signals none (N2); a signal during the stop after a taken receipt path does
not abort it (N5); the stop kills deepest first, and a matrix.json in the
stopped state is not an interrupted run (N10).
"""
import json
import os
import random
import signal
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
# The guard stand-in: `resource_guard.py --reason R --log L --resources PATH -- COMMAND...` (R: a guard reason to
# record, "" for none).
GUARD_STAND_IN = """import json, subprocess, sys
cut = sys.argv.index("--")
opts = dict(zip(sys.argv[1:cut:2], sys.argv[2:cut:2]))
r = subprocess.run(sys.argv[cut + 1:])
open(opts["--resources"], "w").write(json.dumps({"command": sys.argv[cut + 1:], "rc": r.returncode,
                                                  "guard_reason": opts["--reason"] or None}))
"""


def guard_command(base: Path, resources, command, reason=""):
    """The plan's `python3 $GUARD --log ... --resources ... -- COMMAND`, with the stand-in as $GUARD."""
    script = base / "guard-stand-in" / pi.GUARD
    if not script.exists():
        script.parent.mkdir(parents=True, exist_ok=True)
        script.write_text(GUARD_STAND_IN)
    return [sys.executable, str(script), "--reason", reason, "--log", str(base / "guard.log"), "--resources", str(resources),
            "--", *command]


def idle(root):
    """A guarded command that holds `GRM --data ROOT matrix` in its argv but runs only a sleep below the guard."""
    return ["sh", "-c", "exec sleep 120", "--", grm(), "--data", str(root), "matrix"]


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


def sparse(root, slug=SLUG):
    return [grm(), "--data", str(root), "matrix", "--species", slug, "--min-support", "2", "--disk-sparse",
            "--block-width", "64", "--allow-missing-catalog"]


def helper(root, receipt, launcher, resources, folded=64, timeout=120):
    return subprocess.run([sys.executable, "-I", "-B", str(HELPER), "--launcher", str(launcher), "--guard-resources", str(resources),
                           "--data", str(root), "--grm", grm(), "--folded", str(folded), "--receipt", str(receipt),
                           "--timeout", str(timeout), "--poll", "0.05"], capture_output=True, text=True)


def matrix_bytes(root):
    return {f: (root / "species" / SLUG / "matrix" / f).read_bytes() for f in FILES}


def wait_for(predicate, seconds=60):
    deadline = time.time() + seconds
    while time.time() < deadline:
        if predicate():
            return True
        time.sleep(0.05)
    return False


class Interrupt(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.base = Path(self.tmp.name)
        self.started = []

    def tearDown(self):
        # Only processes this test started, and their children.
        for p in self.started:
            if p.poll() is None:
                for pid, (_, start) in sorted(pi.descendants((p.pid, 0)).items(), key=lambda kv: -kv[1][0]):
                    pi.signal_ident((pid, start), signal.SIGKILL)
                p.kill()
            p.wait()
        self.tmp.cleanup()

    def launch(self, resources, command, reason=""):
        return self.start(guard_command(self.base, resources, command, reason))

    def start(self, argv):
        p = subprocess.Popen(argv, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        self.started.append(p)
        return p

    def kill_orphan(self, pid):
        """A run whose guard was killed: it is no longer below anything this test started."""
        st = pi.proc_stat(pid)
        if st is not None:
            pi.signal_ident((pid, st[2]), signal.SIGKILL)

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
        # L1: only a match below the launcher can be bound; one outside its tree is foreign.
        self.assertEqual(pi.split_matches([13, 17], {12: (1, 5), 13: (2, 7)}), ([(13, 7)], [17]))

    def test_the_run_is_stopped_after_block_zero_folds_and_only_grm_is_killed(self):
        reference, _, _ = build_root(self.base / "reference")
        subprocess.run(sparse(reference), check=True, capture_output=True)
        root, fifo, text = self.gated("pilot")
        resources = self.base / "guard.json"
        wrapper = self.launch(resources, sparse(root))
        receipt = self.base / "interrupt.json"
        p = helper(root, receipt, wrapper.pid, resources)
        self.assertEqual(p.returncode, 0, p.stderr)
        self.assertEqual(wrapper.wait(timeout=60), 0)  # the guard stand-in was never signalled
        self.assertEqual(json.loads(resources.read_text())["rc"], -9)  # and saw its child killed
        got = json.loads(receipt.read_text())
        self.assertTrue(got["interrupted"] and got["gone"] and got["launcher_exited"])
        self.assertEqual(got["guard_record"]["rc"], -9)
        self.assertEqual(got["launcher"]["pid"], wrapper.pid)
        state = got["stopped_state"]
        self.assertEqual((state["folded_columns"], state["saved"], state["block_width"], state["matrix_json"]), (64, False, 64, False))
        self.assertIn("block-00000.partial", state["run_files"])
        self.assertIn("block-00000.partial.receipt.json", state["run_files"])
        self.assertEqual(got["after_kill_state"]["run_dir"], state["run_dir"])
        self.assertEqual(got["cmdline"][:4], [grm(), "--data", str(root), "matrix"])
        self.assertEqual(got["matches_at_stop"], {"below_launcher": [got["pid"]], "outside": []})
        # A receipt is never overwritten.
        before = receipt.read_bytes()
        self.assertNotEqual(helper(root, receipt, wrapper.pid, self.base / "guard-2.json", timeout=1).returncode, 0)
        self.assertEqual(receipt.read_bytes(), before)
        # The resume, with the real FASTA in place, gives the uninterrupted matrix.
        real = fifo.with_name(fifo.name + ".real")
        real.write_bytes(text)
        os.replace(real, fifo)
        subprocess.run(sparse(root), check=True, capture_output=True)
        self.assertEqual(matrix_bytes(root), matrix_bytes(reference))
        self.assertFalse((root / "species" / SLUG / "matrix-partials").exists())

    def test_nothing_running_below_the_launcher_is_a_stop_that_stops_the_launch(self):
        root, _, _ = build_root(self.base / "idle")
        resources = self.base / "idle-guard.json"
        wrapper = self.launch(resources, idle(root))
        self.assertTrue(wait_for(lambda: pi.descendants((wrapper.pid, 0))))
        sleeper = list(pi.descendants((wrapper.pid, 0)))[0]
        p = helper(root, self.base / "idle.json", wrapper.pid, resources, timeout=1)
        self.assertEqual(p.returncode, 1)
        got = json.loads((self.base / "idle.json").read_text())
        self.assertIn("appeared below the launcher", got["error"])
        self.assertEqual(got["launched_run_stopped"]["killed"], [sleeper])
        self.assertTrue(got["launched_run_stopped"]["launcher_exited"])
        self.assertEqual(wrapper.wait(timeout=60), 0)
        self.assertEqual(json.loads(resources.read_text())["rc"], -9)

    def test_a_run_that_ends_before_the_kill_point_is_a_stop(self):
        root, fifo, text = self.gated("early")
        resources = self.base / "early-guard.json"
        wrapper = self.launch(resources, sparse(root))
        receipt = self.base / "early.json"
        watcher = subprocess.Popen([sys.executable, "-I", "-B", str(HELPER), "--launcher", str(wrapper.pid), "--guard-resources", str(resources),
                                    "--data", str(root), "--grm", grm(), "--folded", "4096", "--receipt", str(receipt),
                                    "--timeout", "120", "--poll", "0.05"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        self.started.append(watcher)
        # Let block 0 fold, then release the pipe: hand it its bytes once and
        # put the real file at its name for every later read.
        self.assertTrue(wait_for(lambda: (pi.run_state(root / "species" / SLUG)[1] or {}).get("folded_columns") == 64))
        fd = os.open(fifo, os.O_WRONLY)
        real = fifo.with_name(fifo.name + ".real")
        real.write_bytes(text)
        os.replace(real, fifo)
        os.write(fd, text)
        os.close(fd)
        self.assertEqual(wrapper.wait(timeout=120), 0)
        self.assertEqual(json.loads(resources.read_text())["rc"], 0)
        _, err = watcher.communicate(timeout=120)
        self.assertEqual(watcher.returncode, 1, err)
        self.assertIn("ended before the kill point", json.loads(receipt.read_text())["error"])

    def test_a_stale_run_outside_the_launcher_is_never_bound_or_signalled(self):
        # The review's decoy: an exact-argv match already running against the same root when the
        # helper starts. Here it is a real GRM on another species of the root, blocked on a pipe.
        root, fifo, text = self.gated("decoy")
        decoy_sp = root / "species" / "decoy"
        decoy_sp.mkdir(parents=True)
        for i in range(4):
            os.mkfifo(root / "fasta" / f"D{i:03d}.fna.gz")
        (decoy_sp / "genomes.txt").write_text("".join(f"D{i:03d}\n" for i in range(4)))
        decoy = subprocess.Popen(sparse(root, "decoy"), stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        self.started.append(decoy)
        self.assertTrue(wait_for(lambda: decoy.pid in pi.find_run(grm(), str(root))))
        resources = self.base / "decoy-guard.json"
        wrapper = self.launch(resources, sparse(root))
        p = helper(root, self.base / "decoy.json", wrapper.pid, resources)
        self.assertEqual(p.returncode, 1, p.stderr)
        got = json.loads((self.base / "decoy.json").read_text())
        self.assertFalse(got["interrupted"])
        self.assertIn("not launched by this run", got["error"])
        self.assertIn(str(decoy.pid), got["error"])
        stopped = got["launched_run_stopped"]
        self.assertNotIn(decoy.pid, stopped["frozen"] + stopped["killed"])
        self.assertTrue(stopped["launcher_exited"])
        self.assertIsNone(decoy.poll())  # untouched, and not stopped:
        self.assertNotEqual(pi.proc_stat(decoy.pid)[0], "T")
        self.assertEqual(wrapper.wait(timeout=60), 0)  # the launched run is stopped
        self.assertEqual(json.loads(resources.read_text())["rc"], -9)
        self.assertEqual(pi.find_run(grm(), str(root)), [decoy.pid])

    def test_a_match_that_appears_after_binding_fails_the_interrupt_at_the_stop(self):
        # Block 0 itself is held (genome 0 is a pipe), so the helper binds long before the kill
        # point. A stale-looking run then starts outside the launch; when block 0 is released and
        # folds, the scan at the stop must see it and refuse, leaving it untouched.
        root, ids, texts = build_root(self.base / "late")
        gate = root / "fasta" / f"{ids[0]}.fna.gz"
        gate.unlink()
        os.mkfifo(gate)
        resources = self.base / "late-guard.json"
        wrapper = self.launch(resources, sparse(root))
        watcher = subprocess.Popen([sys.executable, "-I", "-B", str(HELPER), "--launcher", str(wrapper.pid), "--guard-resources", str(resources),
                                    "--data", str(root), "--grm", grm(), "--folded", "64", "--receipt", str(self.base / "late.json"),
                                    "--timeout", "120", "--poll", "0.05"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        self.started.append(watcher)
        self.assertTrue(wait_for(lambda: pi.find_run(grm(), str(root))))
        time.sleep(1.0)  # the helper polls every 0.05 s; it has bound the run below the wrapper by now
        late_sp = root / "species" / "late"
        late_sp.mkdir(parents=True)
        for name in ("L000", "L001"):
            os.mkfifo(root / "fasta" / f"{name}.fna.gz")
        (late_sp / "genomes.txt").write_text("L000\nL001\n")
        late = subprocess.Popen(sparse(root, "late"), stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        self.started.append(late)
        self.assertTrue(wait_for(lambda: late.pid in pi.find_run(grm(), str(root))))
        fd = os.open(gate, os.O_WRONLY)  # release block 0: its first read gets the bytes,
        real = gate.with_name(gate.name + ".real")  # and the real file stands at its name
        real.write_bytes(texts[ids[0]])
        os.replace(real, gate)
        os.write(fd, texts[ids[0]])
        os.close(fd)
        _, err = watcher.communicate(timeout=120)
        self.assertEqual(watcher.returncode, 1, err)
        got = json.loads((self.base / "late.json").read_text())
        self.assertFalse(got["interrupted"])
        self.assertEqual(got["matches_at_stop"], {"below_launcher": [got["pid"]], "outside": [late.pid]})
        self.assertIn("at the stop the matching processes", got["error"])
        self.assertIsNone(late.poll())
        self.assertNotEqual(pi.proc_stat(late.pid)[0], "T")
        self.assertEqual(wrapper.wait(timeout=60), 0)

    def test_a_timeout_stops_the_launched_run(self):
        root, _, _ = self.gated("timeout")
        resources = self.base / "timeout-guard.json"
        wrapper = self.launch(resources, sparse(root))
        p = helper(root, self.base / "timeout.json", wrapper.pid, resources, folded=128, timeout=3)
        self.assertEqual(p.returncode, 1)
        got = json.loads((self.base / "timeout.json").read_text())
        self.assertIn("did not reach 128", got["error"])
        self.assertEqual(got["launched_run_stopped"]["killed"], [got["pid"]])
        self.assertEqual(wrapper.wait(timeout=60), 0)
        self.assertEqual(json.loads(resources.read_text())["rc"], -9)
        self.assertEqual(pi.find_run(grm(), str(root)), [])

    def test_a_launcher_that_is_not_a_sibling_job_is_refused_and_nothing_is_signalled(self):
        root, _, _ = build_root(self.base / "wrong")
        resources = self.base / "wrong-guard.json"
        wrapper = self.launch(resources, idle(root))
        self.assertTrue(wait_for(lambda: pi.descendants((wrapper.pid, 0))))
        child = list(pi.descendants((wrapper.pid, 0)))[0]
        for i, launcher in enumerate([1, os.getpid(), child, 0, -1, 2 ** 22 + 7]):
            with self.subTest(launcher=launcher):
                p = helper(root, self.base / f"wrong-{i}.json", launcher, resources, timeout=1)
                self.assertEqual(p.returncode, 1)
                got = json.loads((self.base / f"wrong-{i}.json").read_text())
                self.assertIn("not a running sibling job", got["error"])
                self.assertEqual(got["launched_run_stopped"], {"frozen": [], "killed": [], "launcher_killed": False, "launcher_exited": False})
        self.assertIsNone(wrapper.poll())
        self.assertEqual(pi.proc_stat(child)[0] in ("S", "R"), True)

    def test_an_existing_receipt_is_kept_and_the_launch_is_stopped_before_any_interrupt(self):
        root, _, _ = self.gated("receipt")
        resources = self.base / "receipt-guard.json"
        receipt = self.base / "receipt.json"
        receipt.write_text('{"earlier": true}\n')
        wrapper = self.launch(resources, sparse(root))
        self.assertTrue(wait_for(lambda: pi.find_run(grm(), str(root))))
        p = helper(root, receipt, wrapper.pid, resources)
        self.assertEqual(p.returncode, 1)
        self.assertIn("cannot reserve the receipt", p.stderr)
        self.assertEqual(receipt.read_text(), '{"earlier": true}\n')
        self.assertEqual(wrapper.wait(timeout=60), 0)
        self.assertEqual(json.loads(resources.read_text())["rc"], -9)

    def test_a_guard_record_that_is_not_a_clean_kill_fails_the_interrupt(self):
        root, _, _ = self.gated("reason")
        resources = self.base / "reason-guard.json"
        wrapper = self.launch(resources, sparse(root), reason="disk floor crossed")
        p = helper(root, self.base / "reason.json", wrapper.pid, resources)
        self.assertEqual(p.returncode, 1)
        got = json.loads((self.base / "reason.json").read_text())
        self.assertFalse(got["interrupted"])
        self.assertIn("not a SIGKILLed run", got["error"])
        self.assertEqual(wrapper.wait(timeout=60), 0)

    def test_an_earlier_guard_record_is_a_stop(self):
        root, _, _ = self.gated("stale-record")
        resources = self.base / "stale-guard.json"
        resources.write_text('{"rc": 137, "guard_reason": null, "command": []}\n')
        wrapper = self.launch(resources, sparse(root))  # the launch's own record will land on the earlier one
        self.assertTrue(wait_for(lambda: pi.find_run(grm(), str(root))))
        p = helper(root, self.base / "stale.json", wrapper.pid, resources)
        self.assertEqual(p.returncode, 1)
        self.assertIn("a record from another launch", json.loads((self.base / "stale.json").read_text())["error"])
        self.assertEqual(wrapper.wait(timeout=60), 0)
        self.assertEqual(pi.find_run(grm(), str(root)), [])

    def test_guard_problem_accepts_only_this_runs_guard(self):
        g, d, r = "/x/grm", "/run/cohorts/pilot256/data", "/run/logs/pilot-sparse-1.json"
        run = ["/usr/bin/time", "-v", g, "--data", d, "matrix", "--species", SLUG]
        good = ["python3", "/tools/resource_guard.py", "--log", "l", "--resources", r, "--", *run]
        self.assertIsNone(pi.guard_problem(good, r, g, d))
        self.assertIsNone(pi.guard_problem(["/usr/bin/python3.12", *good[1:]], r, g, d))
        bad = {
            "not python": ["sleep", "300"],
            "not the guard": ["python3", "/tools/other.py", "--resources", r, "--", *run],
            "interpreter flag first": ["python3", "-I", "/tools/resource_guard.py", "--resources", r, "--", *run],
            "other resources": ["python3", "/tools/resource_guard.py", "--resources", r + ".x", "--", *run],
            "resources twice": ["python3", "/tools/resource_guard.py", "--resources", r, "--resources", r, "--", *run],
            "resources after --": ["python3", "/tools/resource_guard.py", "--", "--resources", r, *run],
            "no --": ["python3", "/tools/resource_guard.py", "--resources", r, *run],
            "another root": ["python3", "/tools/resource_guard.py", "--resources", r, "--", g, "--data", d + "-stock", "matrix"],
            "another subcommand": ["python3", "/tools/resource_guard.py", "--resources", r, "--", g, "--data", d, "structure"],
            "empty": [],
            "unreadable": None,
        }
        for name, argv in bad.items():
            with self.subTest(name):
                self.assertIsNotNone(pi.guard_problem(argv, r, g, d))

    def test_a_sibling_that_is_not_this_runs_guard_is_refused_and_nothing_is_signalled(self):
        # N1: the run, a guarded run on another root, and a plain sleep are all sibling jobs of the helper.
        root, _, _ = self.gated("n1")
        other_root, _, _ = self.gated("n1-other")
        resources, other_resources = self.base / "n1-guard.json", self.base / "n1-other-guard.json"
        real = self.launch(resources, sparse(root))
        other = self.launch(other_resources, sparse(other_root))
        plain = self.start(["sleep", "120"])
        self.assertTrue(wait_for(lambda: pi.find_run(grm(), str(root)) and pi.find_run(grm(), str(other_root))))
        runs = {str(root): pi.find_run(grm(), str(root)), str(other_root): pi.find_run(grm(), str(other_root))}
        cases = {"a plain sleep": (plain.pid, resources), "another run's guard": (other.pid, resources),
                 "another run's guard and its own resources": (other.pid, other_resources),
                 "this run's guard, other resources": (real.pid, self.base / "n1-elsewhere.json")}
        for i, (name, (launcher, res)) in enumerate(cases.items()):
            with self.subTest(name):
                p = helper(root, self.base / f"n1-{i}.json", launcher, res, timeout=1)
                self.assertEqual(p.returncode, 1, p.stderr)
                got = json.loads((self.base / f"n1-{i}.json").read_text())
                self.assertIn("is not the resource guard of this run", got["error"])
                self.assertEqual(got["launched_run_stopped"], {"frozen": [], "killed": [], "launcher_killed": False, "launcher_exited": False})
                self.assertEqual([r["pid"] for r in got["still_running"]], runs[str(root)])
                self.assertIn("still running, not signalled", got["error"])
        for p in (real, other, plain):
            self.assertIsNone(p.poll())
        for pids in runs.values():
            self.assertNotEqual(pi.proc_stat(pids[0])[0], "T")
        self.assertEqual({k: pi.find_run(grm(), k) for k in runs}, runs)

    def test_a_guard_that_dies_before_the_bind_leaves_its_run_listed_and_unsignalled(self):
        # N2: the guard is killed before the helper starts; its GRM keeps running outside any launch.
        root, _, _ = self.gated("n2")
        resources = self.base / "n2-guard.json"
        wrapper = self.launch(resources, sparse(root))
        self.assertTrue(wait_for(lambda: pi.find_run(grm(), str(root))))
        [run] = pi.find_run(grm(), str(root))
        try:
            wrapper.kill()
            wrapper.wait()
            p = helper(root, self.base / "n2.json", wrapper.pid, resources, timeout=1)
            self.assertEqual(p.returncode, 1)
            got = json.loads((self.base / "n2.json").read_text())
            self.assertIn("not a running sibling job", got["error"])
            self.assertEqual([r["pid"] for r in got["still_running"]], [run])
            self.assertIn(f"still running, not signalled: [{run}]", p.stderr)
            self.assertNotEqual(pi.proc_stat(run)[0], "T")
            listed = subprocess.run([sys.executable, "-I", "-B", str(HELPER), "--list-runs", "--data", str(root), "--grm", grm()],
                                    capture_output=True, text=True)
            self.assertEqual(listed.returncode, 1, listed.stderr)
            self.assertEqual([r["pid"] for r in json.loads(listed.stdout)], [run])
            self.assertNotEqual(pi.proc_stat(run)[0], "T")
        finally:
            self.kill_orphan(run)
        self.assertTrue(wait_for(lambda: not pi.find_run(grm(), str(root))))
        listed = subprocess.run([sys.executable, "-I", "-B", str(HELPER), "--list-runs", "--data", str(root), "--grm", grm()],
                                capture_output=True, text=True)
        self.assertEqual((listed.returncode, json.loads(listed.stdout)), (0, []))

    def test_a_signal_during_the_stop_after_a_taken_receipt_path_does_not_abort_it(self):
        # N5: the launcher becomes the guard only after 2 s; the helper, stopping for a taken receipt path,
        # gets SIGTERM while it waits for that. It must still stop the launch, without a traceback.
        root, _, _ = self.gated("n5")
        resources, receipt = self.base / "n5-guard.json", self.base / "n5.json"
        receipt.write_text('{"earlier": true}\n')
        guard = guard_command(self.base, resources, sparse(root))
        delayed = self.start([sys.executable, "-c", "import os, sys, time; time.sleep(2); os.execv(sys.executable, [sys.executable, *sys.argv[1:]])",
                              *guard[1:]])
        watcher = subprocess.Popen([sys.executable, "-I", "-B", str(HELPER), "--launcher", str(delayed.pid), "--guard-resources", str(resources),
                                    "--data", str(root), "--grm", grm(), "--folded", "64", "--receipt", str(receipt)],
                                   stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        self.started.append(watcher)
        time.sleep(1.0)
        watcher.send_signal(signal.SIGTERM)
        _, err = watcher.communicate(timeout=120)
        self.assertEqual(watcher.returncode, 1, err)
        self.assertNotIn("Traceback", err)
        self.assertIn("cannot reserve the receipt", err)
        self.assertEqual(receipt.read_text(), '{"earlier": true}\n')
        self.assertEqual(delayed.wait(timeout=60), 0)
        self.assertEqual(pi.find_run(grm(), str(root)), [])

    def test_the_stop_kills_deepest_first(self):
        # N10: a process between the guard and GRM, as /usr/bin/time is; GRM is killed before it.
        root, _, _ = self.gated("deep")
        resources = self.base / "deep-guard.json"
        middle = [sys.executable, "-c", "import subprocess, sys; sys.exit(subprocess.run(sys.argv[1:]).returncode)"]
        wrapper = self.launch(resources, [*middle, *sparse(root)])
        p = helper(root, self.base / "deep.json", wrapper.pid, resources, folded=128, timeout=3)
        self.assertEqual(p.returncode, 1)
        got = json.loads((self.base / "deep.json").read_text())
        self.assertIn("did not reach 128", got["error"])
        killed = got["launched_run_stopped"]["killed"]
        self.assertEqual(len(killed), 2)
        self.assertEqual(killed[0], got["pid"])  # GRM, the deepest, first
        self.assertEqual(wrapper.wait(timeout=60), 0)

    def test_a_matrix_json_in_the_stopped_state_is_not_an_interrupted_run(self):
        # N10 (I5): every other field is the interrupted state.
        state = {"folded_columns": 64, "saved": False, "matrix_json": False, "genomes": 128, "block_width": 64,
                 "run_files": ["block-00000.partial", "block-00000.partial.receipt.json", "run.json"]}
        self.assertEqual(pi.interrupted_problems(state, 64), [])
        self.assertEqual(pi.interrupted_problems({**state, "matrix_json": True}, 64), ["a matrix.json exists"])


if __name__ == "__main__":
    unittest.main()
