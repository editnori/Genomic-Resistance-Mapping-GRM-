"""Plant each fault of docs/outside-learners.md's table and run the test meant to catch it.

  python3 -I planted_faults.py <out dir>

Run it with the python of the venv that holds requirements.txt's wheels: each
test runs under that same interpreter (sys.executable). Copies this
checkout's last commit to <out dir>/tree (git archive) and, one fault at a
time, changes one piece of text there, runs the one test named for it, and
puts the text back, as scripts/c9-pick/planted_faults.py does. A fault is
caught when its test fails; each test must also pass with no fault planted.
Writes <out dir>/faults.tsv (fault, what it plants, file, test, result) and
<out dir>/logs/, prints one line per fault, and exits 1 when a fault is not
caught. Every test it runs is SYNTHETIC; it reads no data.
"""

import importlib.util
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
RUNNER = "scripts/c7-outside/outside.py"
SCORE = "scripts/c7-outside/score.py"
T = "tests.c7_outside.test_runner"
S = "tests.c7_outside.test_score"

FAULTS = [
    ("1", "a library arm that fits on the outer test genomes too", RUNNER,
     '    fit = np.flatnonzero(fo["fit"])\n',
     '    fit = np.flatnonzero(fo["fit"] | fo["test"])\n',
     "py", f"{T}.Arms.test_no_outer_test_label_moves_an_xgboost_score"),
    ("2", "a study fold's boosting holdout drawn over its outer test genomes too", RUNNER,
     "        fit, hold = tuning_holdout(train, ex.y, groups, GBT_HOLDOUT,",
     "        fit, hold = tuning_holdout(train | test, ex.y, groups, GBT_HOLDOUT,",
     "py", f"{T}.Plans.test_study_holdouts_and_inner_folds_lie_inside_each_folds_training_genomes"),
    ("3", "early stopping read on the outer test fold", RUNNER,
     '    hold = np.flatnonzero(fo["hold"])\n',
     '    hold = np.flatnonzero(fo["test"])\n',
     "py", f"{T}.Arms.test_no_outer_test_label_moves_a_lightgbm_score"),
    ("4", "the elastic net's mix and penalty picked on the outer test fold", RUNNER,
     "    pick_on = np.flatnonzero(train)",
     "    pick_on = np.flatnonzero(test)",
     "py", f"{T}.Arms.test_no_outer_test_label_moves_an_elastic_net_score"),
    ("5", "a study report read where a lineage one belongs", SCORE,
     ' and doc.get("scheme") == exp["scheme"]',
     "",
     "py", f"{S}.Checks.test_a_study_report_read_as_a_lineage_one_fails"),
    ("6", "the elastic net scored against GRM's boosted trees, not its lasso", RUNNER,
     '"elasticnet": "lasso"}',
     '"elasticnet": "gbt"}',
     "py", f"{S}.Checks.test_each_arm_pairs_with_its_own_grm_learner"),
    ("7", "the twin scored as an arm", SCORE,
     '    for arm in ("xgboost", "lightgbm", "elasticnet"):\n        learner = ARMS[arm]',
     '    for arm in ("xgboost", "lightgbm", "elasticnet", TWIN):\n        learner = ARMS.get(arm, "gbt")',
     "py", f"{S}.Count.test_conditions_are_read_for_the_three_arms_and_never_for_the_twin"),
    ("8", "a missing draw read as a zero", SCORE,
     'rows.append(dict(base, learner=learner, scored=0, constant=0, why="draw missing"))',
     'rows.append(dict(base, learner=learner, scored=1, constant=0, why="", **{m: 0.0 for m in MEASURES}))',
     "py", f"{S}.Checks.test_a_missing_draw_fails_and_is_never_read_as_a_zero"),
    ("9", "a report that carries feature names", RUNNER,
     '{"rounds_run": len(history["validation"]["logloss"]), "rounds_kept": kept}',
     '{"rounds_run": len(history["validation"]["logloss"]), "rounds_kept": kept, "features": ["row 1", "row 2"]}',
     "py", f"{T}.Arms.test_the_report_holds_only_its_listed_keys_and_no_genome"),
]


def _c9():
    spec = importlib.util.spec_from_file_location("c9_faults_for_c7", HERE.parent / "c9-pick" / "planted_faults.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def run_test(tree, kind, name, log, env):
    """C9's python test run, as CI runs the suite, with this interpreter, whose venv holds the pinned wheels."""
    suite = "/".join(name.split(".")[:2])
    cmd = [sys.executable, "-s", "-P", "-B", "-W", "error", "-m", "unittest", "discover", "-s", suite, "-t", ".",
           "-k", name]
    env = {"PATH": "/usr/bin:/bin", "HOME": env.get("HOME", "/tmp"), "PYTHONHASHSEED": "0"}
    with open(log, "w") as f:
        p = subprocess.run(cmd, cwd=tree, env=env, stdout=f, stderr=subprocess.STDOUT)
    return p.returncode == 0, "Ran 1 test" in Path(log).read_text(errors="replace")


def main(argv):
    c9 = _c9()
    c9.FAULTS = FAULTS
    c9.run_test = run_test
    c9.REPO = HERE.parents[1]
    return c9.main(argv)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
