"""scan.py under a scripted fake binary, and, when smoke.sh provides them,
the real lasso-scan binary on a synthetic data root.

usage: python3 -I -B scripts/lasso-scan/test_scan.py
Set GRM_LASSO_SCAN_BIN and GRM_LASSO_SCAN_FIXTURE (as smoke.sh does) to also
run the RealBinary tests.
"""
import hashlib
import json
import os
import shutil
import signal
import subprocess
import sys
import tempfile
import time
import unittest
from dataclasses import replace
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import scan  # noqa: E402

# Prints a dataset's events as the real binary would, following a per-dataset
# plan in $FAKE_PLAN: sleeps, claimed fit seconds, errors, unconverged points,
# capped points, the number of fits, and the exit status.
FAKE = r'''
import json, os, sys, time
args = sys.argv[1:]
inputs_only = args[0] == "--inputs-only"
root, species, dataset = args[1:] if inputs_only else args
plan = json.load(open(os.environ["FAKE_PLAN"]))[dataset]
with open(plan["pids"], "a") as f:
    f.write(f"{os.getpid()}\n")

def emit(event):
    print(json.dumps(dict(event, dataset=dataset)), flush=True)

time.sleep(plan.get("load_sleep", 0))
emit({"event": "inputs", "samples": 100, "resistant": 40, "features": 50, "lineages": 20,
      "min_support": 3, "input_fingerprint": plan["fingerprint"], "seconds": 0.01})
if inputs_only:
    sys.exit(plan.get("exit", 0))
sites = [("random", f) for f in range(5)] + [("lineage", f) for f in range(5)] + [("final", None)]
for k, (scheme, fold) in enumerate(sites[:plan.get("fits", 11)]):
    step = plan.get("steps", {}).get(str(k), {})
    emit({"event": "start", "scheme": scheme, "outer_fold": fold})
    time.sleep(step.get("sleep", plan.get("fit_sleep", 0)))
    ok = {"irls_steps": 3, "sweep_cap_hits": 0, "converged": True}
    capped = step.get("capped", plan.get("capped", 0))
    path = [dict(ok, sweep_cap_hits=int(i < capped)) for i in range(10)]
    audit = {"path_points": 10, "path_unconverged": step.get("path_unconverged", 0),
             "selected_point": 2, "selected_convergence": ok, "refit_convergence": ok}
    error = step.get("error")
    emit({"event": "fit", "scheme": scheme, "outer_fold": fold,
          "seconds": step.get("seconds", plan.get("fit_seconds", 0.001)),
          "audit": None if error else audit, "error": error,
          "trace": {"path": path, "refit": [ok] * 3},
          "model": None if error else {"intercept": 0.1, "coef": [[1, 0.5]], "lambda": 0.01}})
    if error:
        sys.exit(1)
sys.exit(plan.get("exit", 0))
'''

TEST_RULE = replace(scan.REGISTERED, fit_floor_seconds=0.5)
REVISION7 = {"reports": 275, "lasso_seconds": 14989.0, "report_seconds": 25088.0}


def entry(dataset, lasso_seconds=1.0, rank=1):
    """A list entry; r7_fit is a tenth of `lasso_seconds`."""
    return {"rank": rank, "species": dataset.split("__")[0], "id": dataset, "samples": 100, "features": 50,
            "lasso_seconds": lasso_seconds, "outer_fits": 10, "report_seconds": 3.0,
            "input_fingerprint": hashlib.sha256(dataset.encode()).hexdigest()[:16]}


def write_report(root, e, fingerprint=None):
    path = Path(root, "species", e["species"], "results", e["id"] + ".json")
    path.parent.mkdir(parents=True, exist_ok=True)

    def split(scheme):
        return {"scheme": scheme, "folds": 5, "learners": [{"learner": "lasso", "seconds": e["lasso_seconds"] / 2}]}

    path.write_text(json.dumps({"evaluation": {"revision": 7, "valid": True}, "id": e["id"], "species_slug": e["species"],
                                "input_fingerprint": fingerprint or e["input_fingerprint"], "seconds": 3.0,
                                "splits": [split("random"), split("lineage")]}))


class Harness(unittest.TestCase):
    """A temporary input root, list, fake binary, and output path."""

    def setUp(self):
        tmp = tempfile.TemporaryDirectory(prefix="grm-lasso-scan-test.")
        self.addCleanup(tmp.cleanup)
        self.tmp = Path(tmp.name)
        self.root = self.tmp / "root"
        self.out = self.tmp / "out"
        self.pids = self.tmp / "pids"
        self.fake = self.tmp / "fake-lasso-scan"
        self.fake.write_text(f"#!{sys.executable} -I\n{FAKE}")
        self.fake.chmod(0o755)
        self.plan_file = self.tmp / "fake-plan.json"
        self.addCleanup(setattr, scan, "free_bytes", scan.free_bytes)

    def setup_scan(self, entries, fake=None, revision7=REVISION7):
        for e in entries:
            write_report(self.root, e)
        self.list = self.tmp / "scan-list.json"
        self.list.write_text(json.dumps({"revision7": revision7, "datasets": entries}))
        plans = {e["id"]: {"fingerprint": e["input_fingerprint"], "pids": str(self.pids)} for e in entries}
        for dataset, extra in (fake or {}).items():
            plans[dataset].update(extra)
        self.plan_file.write_text(json.dumps(plans))
        os.environ["FAKE_PLAN"] = str(self.plan_file)

    def args(self, *extra):
        return scan.parse(["--input-root", str(self.root), "--out", str(self.out), "--binary", str(self.fake),
                           "--list", str(self.list), *extra])

    def run_scan(self, *extra, rule=TEST_RULE):
        started = time.monotonic()
        rc = scan.run(self.args(*extra), rule)
        return rc, time.monotonic() - started

    def summary(self):
        return json.loads((self.out / "summary.json").read_text())

    def stop(self):
        return json.loads((self.out / "stop.json").read_text())

    def assert_no_child_alive(self):
        for pid in map(int, self.pids.read_text().split()):
            for _ in range(50):
                try:
                    os.kill(pid, 0)
                except ProcessLookupError:
                    break
                time.sleep(0.1)
            else:
                self.fail(f"binary {pid} outlived the stop")

    def assert_stopped(self, rc, kind):
        self.assertEqual(rc, 2)
        self.assertEqual(self.stop()["kind"], kind, self.stop())
        self.assertEqual(self.summary()["status"], "stopped")
        self.assert_no_child_alive()

    def need_floor_space(self):
        """Runs under the registered rule need its 30 GiB disk floor free."""
        if scan.free_bytes(self.tmp) < scan.REGISTERED.disk_floor_bytes + 2**30:
            self.skipTest("less than the registered disk floor is free here")


class RegisteredRule(unittest.TestCase):
    def test_the_registered_numbers(self):
        r = scan.REGISTERED
        self.assertEqual((r.fit_floor_seconds, r.fit_factor, r.dataset_factor), (60.0, 20.0, 5.0))
        self.assertEqual((r.max_ratio, r.max_median_ratio, r.max_projected_hours, r.cap_review_share), (5.0, 3.0, 16.0, 0.05))
        self.assertEqual((r.final_fit_share, r.disk_floor_bytes, r.rayon_threads), (1.1, 30 * 2**30, 64))

    def test_the_committed_list_and_its_budgets(self):
        scan_list = json.loads((HERE / "scan-list.json").read_text())
        datasets = scan_list["datasets"]
        self.assertEqual([d["rank"] for d in datasets], [1, 2, 3, 4, 9, 11, 13, 14, 16, 17, 20, 21, 25, 26, 29, 37, 38, 65])
        self.assertEqual(scan_list["revision7"]["species_by_lasso_seconds"],
                         ["mycobacterium_tuberculosis", "klebsiella_pneumoniae", "salmonella_enterica", "escherichia_coli", "staphylococcus_aureus"])
        self.assertEqual(round(scan_list["revision7"]["lasso_seconds"]), 14989)
        self.assertEqual(round(scan_list["revision7"]["report_seconds"]), 25088)
        b = [scan.budgets(d, scan.REGISTERED, 1.0) for d in datasets]
        self.assertEqual(datasets[0]["id"], "mycobacterium_tuberculosis__delamanid")
        self.assertEqual((round(b[0]["fit_seconds"]), round(b[0]["dataset_seconds"])), (842, 2315))
        self.assertEqual(round(sum(x["dataset_seconds"] for x in b)), 20026)
        need = scan.FIXED_BYTES + sum(scan.disk_need(d) for d in datasets)
        self.assertEqual(need, 75_249_664)
        # The projections the go/no-go limit was set from.
        hours = lambda m: (25088 + (m - 1) * 1.1 * 14989) / 3600  # noqa: E731
        self.assertEqual([round(hours(m), 1) for m in (1.92, 3, 5)], [11.2, 16.1, 25.3])


class Scans(Harness):
    def test_a_clean_scan_is_a_go_and_records_every_fit(self):
        entries = [entry("alpha__one", 1.0, 1), entry("beta__two", 0.8, 2)]
        self.setup_scan(entries, {"beta__two": {"capped": 1}})
        rc, _ = self.run_scan()
        self.assertEqual(rc, 0)
        s = self.summary()
        self.assertEqual((s["status"], s["decision"], s["reasons"]), ("complete", "go", []))
        self.assertEqual([r["id"] for r in s["datasets"]], ["alpha__one", "beta__two"])
        self.assertAlmostEqual(s["datasets"][0]["ratio"], 0.01)
        self.assertEqual((s["points"], s["capped_points"]), (2 * 11 * 13, 11))
        events = [json.loads(line) for line in (self.out / "events.jsonl").read_text().splitlines()]
        self.assertEqual([e["event"] for e in events].count("fit"), 22)
        plan = json.loads((self.out / "plan.json").read_text())
        self.assertEqual(plan["binary"]["sha256"], scan.sha256(self.fake))
        self.assertEqual(plan["rule"]["fit_floor_seconds"], 0.5)

    def test_a_lasso_error_stops_at_once_and_runs_nothing_after(self):
        error = "lasso refit failure at dataset alpha__one split lineage outer fold 1"
        self.setup_scan([entry("alpha__one"), entry("beta__two", rank=2)], {"alpha__one": {"steps": {"6": {"error": error}}}})
        rc, _ = self.run_scan()
        self.assert_stopped(rc, "lasso error")
        self.assertIn("lineage outer fold 1", self.stop()["detail"])
        self.assertEqual(self.stop()["completed"], [])
        self.assertEqual(len(self.pids.read_text().split()), 1)

    def test_an_unconverged_point_stops_even_when_another_was_selected(self):
        self.setup_scan([entry("alpha__one")], {"alpha__one": {"steps": {"3": {"path_unconverged": 1}}, "fit_sleep": 0.2}})
        rc, elapsed = self.run_scan()
        self.assert_stopped(rc, "unconverged point")
        # It stopped on the fourth fit, not after the remaining seven.
        self.assertLess(elapsed, 2.0)

    def test_a_fit_over_its_budget_is_killed_at_the_budget(self):
        # r7_fit 0.1 s: the fit budget is max(0.5, 20 * 0.1) = 2 s.
        self.setup_scan([entry("alpha__one")], {"alpha__one": {"steps": {"2": {"sleep": 30}}}})
        rc, elapsed = self.run_scan()
        self.assert_stopped(rc, "fit budget overrun")
        self.assertLess(elapsed, 6)

    def test_a_claimed_fit_time_over_budget_stops(self):
        self.setup_scan([entry("alpha__one")], {"alpha__one": {"steps": {"0": {"seconds": 2.5}}}})
        rc, _ = self.run_scan()
        self.assert_stopped(rc, "fit budget overrun")

    def test_a_dataset_over_its_budget_stops(self):
        # The second dataset's r7_fit is 0.02 s: its fits get the 0.5 s floor
        # and the dataset 1.1 s. The first leaves the whole scan 55 s of slack.
        # Three 0.3 s fits leave 0.2 s, so the fourth is stopped mid-fit by the
        # dataset budget, before its own; claimed times alone also count.
        cases = [
            ({"fit_sleep": 0.3, "fit_seconds": 0.3, "steps": {"3": {"sleep": 30}}}, "while running"),
            ({"fit_seconds": 0.3}, "from claimed times"),
        ]
        for fake, case in cases:
            with self.subTest(case=case):
                self.setUp()
                self.setup_scan([entry("alpha__one", 10.0), entry("beta__two", 0.2, 2)], {"beta__two": fake})
                rc, elapsed = self.run_scan()
                self.assert_stopped(rc, "dataset budget overrun")
                self.assertEqual(self.stop()["completed"], ["alpha__one"])
                self.assertLess(elapsed, 4)

    def test_the_whole_scan_budget_counts_loading_too(self):
        # One dataset with a 0.55 s budget that takes 30 s to load.
        self.setup_scan([entry("alpha__one", 0.1)], {"alpha__one": {"load_sleep": 30}})
        rc, elapsed = self.run_scan()
        self.assert_stopped(rc, "scan budget overrun")
        self.assertLess(elapsed, 5)

    def test_inputs_that_are_not_the_revision7_inputs_stop(self):
        self.setup_scan([entry("alpha__one")], {"alpha__one": {"fingerprint": "ffffffffffffffff"}})
        rc, _ = self.run_scan()
        self.assert_stopped(rc, "inputs differ from revision 7")

    def test_a_failed_binary_or_a_short_run_stops(self):
        for fake, kind in (({"exit": 101}, "scan binary failed"), ({"fits": 10}, "unexpected output")):
            with self.subTest(fake=fake):
                self.setUp()
                self.setup_scan([entry("alpha__one")], {"alpha__one": fake})
                rc, _ = self.run_scan()
                self.assert_stopped(rc, kind)

    def test_no_go_limits(self):
        # With r7_fit 1 s, a claimed outer-fit time of r seconds sets the
        # dataset's ratio to r, inside every budget.
        cases = [
            ([5.2, 1.0, 1.0], ["worst ratio 5.20 is above 5.0"]),
            ([3.5, 3.5, 1.0], ["median ratio 3.50 is above 3.0", "the 275-report run projects to 18.4 h, above 16.0 h"]),
            ([3.0, 3.0, 3.0], ["the 275-report run projects to 16.1 h, above 16.0 h"]),
            ([2.9, 2.9, 5.0], []),
        ]
        for ratios, reasons in cases:
            with self.subTest(ratios=ratios):
                self.setUp()
                entries = [entry(f"s{i}__d", 10.0, i + 1) for i in range(3)]
                self.setup_scan(entries, {e["id"]: {"fit_seconds": r, "steps": {"10": {"seconds": 0.001}}} for e, r in zip(entries, ratios)})
                rc, _ = self.run_scan()
                s = self.summary()
                self.assertEqual((rc, s["status"], s["reasons"]), (1 if reasons else 0, "complete", reasons))
                self.assertEqual(s["decision"], "no-go" if reasons else "go")

    def test_sweep_cap_hits_on_over_five_percent_of_points_call_for_review(self):
        # 143 solved points: 7 capped is 4.9%, 8 is 5.6%.
        for fits, rc_want in ((7, 0), (8, 1)):
            with self.subTest(fits=fits):
                self.setUp()
                self.setup_scan([entry("alpha__one")], {"alpha__one": {"steps": {str(k): {"capped": 1} for k in range(fits)}}})
                rc, _ = self.run_scan()
                self.assertEqual(rc, rc_want)
        self.assertEqual(self.summary()["reasons"], ["8 of 143 solved points (5.6%) hit the sweep cap, above 5%: review first"])

    def test_inputs_only_loads_and_checks_without_fitting(self):
        self.setup_scan([entry("alpha__one"), entry("beta__two", rank=2)])
        rc, _ = self.run_scan("--inputs-only")
        self.assertEqual(rc, 0)
        s = self.summary()
        self.assertEqual([r["input_fingerprint"] for r in s["datasets"]],
                         [entry("alpha__one")["input_fingerprint"], entry("beta__two")["input_fingerprint"]])
        events = [json.loads(line)["event"] for line in (self.out / "events.jsonl").read_text().splitlines()]
        self.assertEqual(events, ["inputs", "inputs"])


class Refusals(Harness):
    def refused(self, *extra):
        rc, _ = self.run_scan(*extra)
        self.assertEqual(rc, 3)
        self.assertFalse(self.pids.exists(), "nothing may run after a refusal")

    def test_a_report_that_is_not_the_listed_revision7_report_is_refused(self):
        self.setup_scan([entry("alpha__one")])
        write_report(self.root, entry("alpha__one"), fingerprint="0000000000000000")
        self.refused()
        self.assertFalse(self.out.exists())

    def test_a_used_or_nested_output_directory_is_refused(self):
        self.setup_scan([entry("alpha__one")])
        self.out.mkdir()
        (self.out / "old").write_text("x")
        self.refused()
        self.out = self.root / "scan"
        self.refused()

    def test_free_space_minus_the_estimate_must_stay_above_the_floor(self):
        self.setup_scan([entry("alpha__one")])
        need = scan.FIXED_BYTES + scan.disk_need(entry("alpha__one"))
        scan.free_bytes = lambda path: TEST_RULE.disk_floor_bytes + need - 1
        self.refused()
        scan.free_bytes = lambda path: TEST_RULE.disk_floor_bytes + need
        rc, _ = self.run_scan()
        self.assertEqual(rc, 0)

    def test_the_plan_runs_nothing(self):
        self.setup_scan([entry("alpha__one")])
        rc, _ = self.run_scan("--plan")
        self.assertEqual(rc, 0)
        self.assertFalse(self.out.exists())
        self.assertFalse(self.pids.exists())


class DiskFloor(Harness):
    def test_free_space_is_checked_before_each_dataset(self):
        entries = [entry("alpha__one"), entry("beta__two", rank=2)]
        self.setup_scan(entries)
        need = scan.FIXED_BYTES + sum(scan.disk_need(e) for e in entries)
        calls = []

        def free(path):
            calls.append(path)
            # Enough at launch and before the first dataset; then another job
            # takes the headroom.
            return TEST_RULE.disk_floor_bytes + need if len(calls) <= 2 else TEST_RULE.disk_floor_bytes

        scan.free_bytes = free
        rc, _ = self.run_scan()
        self.assert_stopped(rc, "disk floor")
        self.assertEqual(self.stop()["completed"], ["alpha__one"])
        self.assertIn("before beta__two", self.stop()["detail"])


class Signals(Harness):
    """The driver as a process, as the launch runs it, with the registered rule."""

    def start_driver(self):
        self.need_floor_space()
        self.setup_scan([entry("alpha__one", 100.0)], {"alpha__one": {"steps": {"1": {"sleep": 60}}}})
        cmd = [sys.executable, "-I", "-B", str(HERE / "scan.py"), "--input-root", str(self.root), "--out", str(self.out),
               "--binary", str(self.fake), "--list", str(self.list)]
        driver = subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, env=os.environ)
        events = self.out / "events.jsonl"
        for _ in range(200):
            if events.exists() and events.read_text().count('"event": "start"') >= 2:
                return driver
            time.sleep(0.1)
        driver.kill()
        self.fail("the driver did not reach the second fit")

    def test_sigterm_stops_the_scan_on_the_record(self):
        driver = self.start_driver()
        driver.send_signal(signal.SIGTERM)
        self.assertEqual(driver.wait(timeout=30), 2)
        self.assert_stopped(2, "driver stopped")

    def test_a_killed_driver_takes_the_binary_with_it(self):
        driver = self.start_driver()
        driver.kill()
        driver.wait(timeout=30)
        self.assert_no_child_alive()


@unittest.skipUnless(os.environ.get("GRM_LASSO_SCAN_BIN") and os.environ.get("GRM_LASSO_SCAN_FIXTURE"), "set by smoke.sh")
class RealBinary(Harness):
    """The lasso-scan binary, registered rule, on the synthetic root that
    analysis::tests::write_lasso_scan_fixture_root writes."""

    def setUp(self):
        super().setUp()
        self.need_floor_space()
        self.fake = Path(os.environ["GRM_LASSO_SCAN_BIN"])
        shutil.copytree(os.environ["GRM_LASSO_SCAN_FIXTURE"], self.root)
        self.report = self.root / "species/fixture/results/fixture__meropenem.json"
        self.pids.write_text("")

    def use_list(self, lasso_split_seconds=None):
        report = json.loads(self.report.read_text())
        if lasso_split_seconds is not None:
            for split in report["splits"]:
                split["learners"][0]["seconds"] = lasso_split_seconds
            self.report.write_text(json.dumps(report))
        e = {"rank": 1, "species": "fixture", "id": "fixture__meropenem", "samples": 150, "features": 5,
             "lasso_seconds": sum(s["learners"][0]["seconds"] for s in report["splits"]),
             "outer_fits": sum(s["folds"] for s in report["splits"]), "report_seconds": report["seconds"],
             "input_fingerprint": report["input_fingerprint"]}
        self.list = self.tmp / "scan-list.json"
        self.list.write_text(json.dumps({"revision7": REVISION7, "datasets": [e]}))

    def test_the_real_binary_scans_the_fixture_to_a_go(self):
        self.use_list()
        rc, elapsed = self.run_scan(rule=scan.REGISTERED)
        self.assertEqual(rc, 0, (self.out / "fixture__meropenem.stderr.log").read_text())
        s = self.summary()
        self.assertEqual((s["decision"], s["datasets"][0]["valid"], s["datasets"][0]["samples"]), ("go", True, 150))
        events = [json.loads(line) for line in (self.out / "events.jsonl").read_text().splitlines()]
        self.assertEqual([e["event"] for e in events], ["inputs"] + ["start", "fit"] * 11)
        fits = [e for e in events if e["event"] == "fit"]
        self.assertTrue(all(f["error"] is None and f["audit"]["path_unconverged"] == 0 for f in fits))
        self.assertEqual([f["scheme"] for f in fits], ["random"] * 5 + ["lineage"] * 5 + ["final"])
        written = sum(p.stat().st_size for p in self.out.iterdir())
        print(f"\nreal binary: {elapsed:.2f} s, {written} bytes written", file=sys.stderr)

    def test_the_real_binary_checks_inputs_only(self):
        self.use_list()
        rc, _ = self.run_scan("--inputs-only", rule=scan.REGISTERED)
        self.assertEqual(rc, 0)

    def test_a_budget_overrun_kills_the_real_binary(self):
        self.use_list(lasso_split_seconds=1e-6)
        rc, _ = self.run_scan(rule=scan.REGISTERED)
        self.assertEqual(rc, 2)
        self.assertTrue(self.stop()["kind"].endswith("budget overrun"), self.stop())

    def test_other_inputs_stop_the_real_binary(self):
        report = json.loads(self.report.read_text())
        report["input_fingerprint"] = "0123456789abcdef"
        self.report.write_text(json.dumps(report))
        self.use_list()
        rc, _ = self.run_scan(rule=scan.REGISTERED)
        self.assertEqual(rc, 2)
        self.assertEqual(self.stop()["kind"], "inputs differ from revision 7")


if __name__ == "__main__":
    unittest.main(verbosity=2)
