#!/usr/bin/env python3
"""Mutation check for the revised checker, audit, and audit scripts.

Each mutant breaks one rule in a scratch copy of scripts/check-results.py
and this lane (scripts and tests, with foldgen linked) and runs the test
file that owns the rule. A mutant is killed only when named tests in that
file fail; a file that cannot run counts as not killed. The unmutated copy
must pass all three files first. Every `before` text must occur exactly
once.

usage: python3 -I tests/mutate.py OUT.json
"""
import json
import os
import shutil
import subprocess
import sys
import tempfile

LANE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
NAME = os.path.basename(LANE)
CHECK, AUDIT, SCRIPTS = "tests/test_check_results.py", "tests/test_audit.py", "tests/test_scripts.py"
# Mutated files, relative to the scratch data/logs.
C, A = "check-results.py", NAME + "/audit.py"
SUMMARIZE, TRANSITIONS, NULLIMPACT = NAME + "/summarize.py", NAME + "/transitions.py", NAME + "/nullimpact.py"

MUTANTS = [
    # check-results.py
    ("schema-1 rows may hold null scores", C, [('if learner["learner"] == "catalog" or schema == 1:', 'if learner["learner"] == "catalog":')], CHECK),
    ("no_call need not match the null scores", C, [('assert no_call == uncalled, ("no_call differs from the null scores",) + identity', "pass")], CHECK),
    ("scm, cart, and rf may make no call", C, [('assert learner["learner"] in ("gbt", "lasso") or no_call == {"resistant": 0, "susceptible": 0}, ("only gbt and lasso fits can be null",) + identity', "pass")], CHECK),
    ("metrics n counts uncalled isolates", C, [('assert metrics["n"] == len(called) and', 'assert metrics["n"] == n and')], CHECK),
    ("balanced accuracy scores uncalled isolates 0.5", C, [("balanced_accuracy([truth[i] for i in called], [scores[i] for i in called])",
                                                         "balanced_accuracy(truth, [0.5 if v is None else v for v in scores])")], CHECK),
    ("a row without calls may keep metrics", C, [('assert metrics is None and learner["roc"] == []', 'assert learner["roc"] == []')], CHECK),
    ("a fold may be called in part", C, [('assert None not in fold and fold_ba[k] is not None, ("a fold is called in part",) + identity + (k,)', "pass")], CHECK),
    ("fold balanced accuracy is not recomputed", C, [('assert math.isclose(fold_ba[k], balanced_accuracy([truth[i] for i in held], fold), abs_tol=1e-12), ("fold balanced accuracy",) + identity + (k,)', "pass")], CHECK),
    ("test_fold shape is not checked", C, [('assert folds is not None and len(folds) == n and all(0 <= f < split["folds"] for f in folds), ("test_fold", path, split["scheme"])',
                                            'assert folds is not None, ("test_fold", path, split["scheme"])')], CHECK),
    ("a test fold may lack a class", C, [('assert any(held) and not all(held), ("test fold lacks a class", path, split["scheme"], k)', "pass")], CHECK),
    ("a lineage may span test folds", C, [('assert all(len(v) == 1 for v in lineage_folds.values()), ("a lineage spans test folds", path)', "pass")], CHECK),
    ("a schema-1 split may save test_fold", C, [('assert folds is None, ("test_fold in a schema-1 report", path, split["scheme"])', "pass")], CHECK),
    ("an unavailable split may save test_fold", C, [('and not split["learners"] and "test_fold" not in split, ("unavailable split", path)',
                                                     'and not split["learners"], ("unavailable split", path)')], CHECK),
    ("any schema passes with any revision", C, [('assert schema == (2 if args.reshaped or audit["revision"] >= 7 else 1)', "assert schema in (1, 2)")], CHECK),
    ("revision-6 learner-row count", C, [('assert names == LEARNERS + (["catalog"] if report.get("catalog") else []), ("learner rows"',
                                          'assert len(names) >= 5, ("learner rows"')], CHECK),
    ("winner ignores lineage coverage", C, [("if random is None or not calls_every_isolate(random) or (lineage_ran and (lineage is None or not calls_every_isolate(lineage))):",
                                             "if random is None or not calls_every_isolate(random):")], CHECK),
    ("winner tie goes to the earlier learner", C, [('if best is None or random["metrics"]["balanced_accuracy"] >= best[1]:',
                                                    'if best is None or random["metrics"]["balanced_accuracy"] > best[1]:')], CHECK),
    ("schema-2 winner counts the catalog", C, [("    for name in LEARNERS:\n", '    for name in LEARNERS + ["catalog"]:\n')], CHECK),
    ("no-winner done line is not parsed", C, [("if words == NO_WINNER:", "if False:")], CHECK),
    ("done-line form is not checked", C, [('assert line is None or line[0] == ("full coverage" if schema == 2 else "any row")', "assert True")], CHECK),
    ("revision-7 report needs no done line", C, [('assert line is not None or audit["revision"] < 7, ("revision-7 report without a done line", path)', "pass")], CHECK),
    # audit/audit.py
    ("audit scores uncalled isolates 0.5", A, [("called = [(t, s) for t, s in zip(truth, scores) if s is not None]",
                                                "called = [(t, 0.5 if s is None else s) for t, s in zip(truth, scores)]")], AUDIT),
    ("schema-2 null folds found by 0.5 scores", A, [("return all(v is None for v in held)", "return all(v == 0.5 for v in held)")], AUDIT),
    ("null fold and zero complexity not compared", A, [("if null != zero:", "if False:")], AUDIT),
    ("test_fold not compared with the rebuilt folds", A, [("if saved != f:", "if False:")], AUDIT),
    ("no_call not compared with the null scores", A, [("elif nc != uncalled:", "elif False:")], AUDIT),
    ("null matches a number", A, [("return a is not b", "return False")], AUDIT),
    ("folds rebuilt from the default seed", A, [('return cfg.get("folds", K), cfg.get("seed", SEED)', "return K, SEED")], AUDIT),
    ("schema and revision not paired", A, [('check(schema == (2 if RESHAPED or rec["revision"] >= 7 else 1), "schema does not match the revision"',
                                            'check(True, "schema does not match the revision"')], AUDIT),
    ("eval_config optional from revision 7", A, [('check(rec["revision"] < 7, "report from revision 7 on without eval_config")', 'check(True, "")')], AUDIT),
    ("a recorded forest seed is accepted", A, [('check(set(cfg.get("forest", {})) == FOREST_KEYS,', 'check(set(cfg.get("forest", {})) <= FOREST_KEYS | {"seed"},')], AUDIT),
    ("no-winner done line is not parsed", A, [("if words == NO_WINNER:", "if False:")], AUDIT),
    ("any-row done line accepted for schema 2", A, [('check(logged["form"] == "any row" and schema < 2,', 'check(logged["form"] == "any row",')], AUDIT),
    ("winner ignores lineage coverage", A, [("if random is None or not calls_every_isolate(random) or (lineage_ran and (lineage is None or not calls_every_isolate(lineage))):",
                                             "if random is None or not calls_every_isolate(random):")], AUDIT),
    ("schema-2 winner by the schema-1 manifest rule", A, [("best = full_coverage_best(rep) if schema >= 2 else max_random_best(rep)", "best = max_random_best(rep)")], AUDIT),
    ("manifest winner by the schema-1 rule", A, [('mbest = full_coverage_best(rep) if manifest["schema"] >= 2 else max_random_best(rep)', "mbest = max_random_best(rep)")], AUDIT),
    ("schema-1 manifest rule fails on null metrics", A, [('ba = (random["metrics"] or {}).get("balanced_accuracy", 0.0)', 'ba = random["metrics"]["balanced_accuracy"]')], AUDIT),
    ("manifest not compared", A, [('if rec["manifest_compared"]:', "if False:")], AUDIT),
    ("null manifest rates read as 0", A, [('out = dict.fromkeys(keys + ("auroc",))', 'out = dict.fromkeys(keys + ("auroc",), 0.0)')], AUDIT),
    ("panel rows not checked", A, [('issues.append({"where": "manifest/panels"', '[].append({"where": "manifest/panels"')], AUDIT),
    ("revision-6 probe no_call rule", A, [('check(rec["revision"] < 6, "revision-%d probe panel without no_call"', 'check(rec["revision"] != 6, "revision-%d probe panel without no_call"')], AUDIT),
    ("launch time not checked", A, [('check(rec["written_after_launch"], "revision-%d file older than launch" % REVISION)', 'check(True, "")')], AUDIT),
    ("recompute line not required", A, [('check((rec["log_recomputed"] is not None) == rec["backup_exists"],', "check(True,")], AUDIT),
    ("unchanged forest accepted", A, [('elif a["learner"] == "rf" and rec["revision"] >= 7 and old["evaluation"]["revision"] < 7:', "elif False:")], AUDIT),
    ("rf feature rows not range-checked", A, [('check(all(isinstance(fid, int) and 0 <= fid < c["features"] for fid in ids), "rf_importance feature is not a dataset row",',
                                               'check(True, "rf_importance feature is not a dataset row",')], AUDIT),
    ("rf card counts not checked", A, [('if counted != [cards[fid]["n_resistant"], cards[fid]["n_susceptible"]]:', "if False:")], AUDIT),
    ("revision-6 check closure", A, [('def check(cond, what, where="counts/structure", **detail):', "def check(cond, what, **detail):"),
                                     ("issues.append(dict(where=where, what=what, **detail))", 'issues.append(dict(where="counts/structure", what=what, **detail))')], AUDIT),
    # audit scripts
    ("summarize sorts null beside numbers", SUMMARIZE, [('"fold_ba_values": null_first(', '"fold_ba_values": sorted(')], SCRIPTS),
    ("summarize keeps one row's no-calls per key", SUMMARIZE,
     [('uncalled["%s/%s" % (x["scheme"], x["learner"])] += ', 'uncalled["%s/%s" % (x["scheme"], x["learner"])] = ')], SCRIPTS),
    ("transitions misses null folds", TRANSITIONS, [("if schema >= 2 and audit.is_null_fold(schema, held):", "if False:")], SCRIPTS),
    ("nullimpact counts no metrics as 0", NULLIMPACT, [("return None if kept is None else kept[key]", "return 0.0 if kept is None else kept[key]")], SCRIPTS),
]


def run_tests(lane, test):
    return subprocess.run([sys.executable, "-I", "-W", "ignore::ResourceWarning", test], cwd=lane, capture_output=True, text=True)


def main():
    out = sys.argv[1]
    with tempfile.TemporaryDirectory() as tmp:
        logs = os.path.join(tmp, "logs")
        lane = os.path.join(logs, NAME)
        os.makedirs(lane)
        shutil.copy2(os.path.join(os.path.dirname(LANE), "check-results.py"), logs)
        for name in ("audit.py", "nullimpact.py", "summarize.py", "transitions.py"):
            shutil.copy2(os.path.join(LANE, name), lane)
        os.symlink(os.path.join(LANE, "foldgen"), os.path.join(lane, "foldgen"))
        shutil.copytree(os.path.join(LANE, "tests"), os.path.join(lane, "tests"), ignore=shutil.ignore_patterns("__pycache__"))
        baseline = {test: run_tests(lane, test).returncode for test in (CHECK, AUDIT, SCRIPTS)}
        assert all(code == 0 for code in baseline.values()), baseline
        results = []
        for name, relative, edits, test in MUTANTS:
            path = os.path.join(logs, relative)
            original = open(path).read()
            text = original
            for before, after in edits:
                assert text.count(before) == 1, (name, before, text.count(before))
                text = text.replace(before, after)
            open(path, "w").write(text)
            try:
                result = run_tests(lane, test)
            finally:
                open(path, "w").write(original)
            failed = [line.split(" (")[0] for line in result.stderr.splitlines() if line.startswith(("FAIL: ", "ERROR: "))]
            killed = result.returncode != 0 and bool(failed)
            results.append({"mutant": name, "file": relative, "test_file": test, "exit": result.returncode,
                            "killed": killed, "failing_tests": failed})
            print("%-48s %-38s %s  %d failing: %s" % (name, relative, "killed" if killed else "NOT KILLED", len(failed),
                                                       ", ".join(f.split(": ")[1] for f in failed[:3])), flush=True)
    summary = {"baseline": baseline, "mutants": len(results), "killed": sum(r["killed"] for r in results),
               "survived": [r["mutant"] for r in results if not r["killed"]], "results": results}
    json.dump(summary, open(out, "w"), indent=1)
    print(json.dumps({k: v for k, v in summary.items() if k != "results"}))
    sys.exit(1 if summary["survived"] else 0)


if __name__ == "__main__":
    main()
