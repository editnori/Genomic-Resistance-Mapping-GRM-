"""Plant each fault of docs/input-arms.md's table and run the test meant to catch it.

  python3 -I planted_faults.py <out dir>

Copies this checkout's last commit to <out dir>/tree (git archive) and, one
fault at a time, changes one piece of text there, runs the one test named
for it, and puts the text back. A fault is caught when its test fails and
the tree still builds; each test must also pass with no fault planted. Rust
tests build in <out dir>/target, whatever CARGO_TARGET_DIR says, as
scripts/c9-pick/planted_faults.py explains. Writes <out dir>/faults.tsv
(fault, what it plants, file, test, result) and <out dir>/logs/, prints one
line per fault, and exits 1 when a fault is not caught. Every test it runs
is SYNTHETIC; it reads no data.
"""

import importlib.util
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ARMS = "src/input_arms.rs"
SCORE = "scripts/c3-inputs/score.py"
T = "tests.c3_inputs.test_score"

FAULTS = [
    ("1", "a cohort filter that reads the outer test genomes' calls", ARMS,
     "        if s < min_support || s > max {",
     "        if s < 2 * min_support || s > max {",
     "rust", "input_arms::tests::no_outer_test_genome_changes_a_fits_candidates"),
    ("2a", "a joined arm that drops the K31 block", ARMS,
     "    if arm == Arm::K31Program {\n        for r in 0..ds.x.rows {",
     "    if arm == Arm::K31Program && ds.x.rows == 0 {\n        for r in 0..ds.x.rows {",
     "rust", "input_arms::tests::the_joined_arm_holds_each_block_once_the_program_first"),
    ("2b", "a joined arm that holds the program block twice", ARMS,
     "    for r in 0..program.x.rows {\n        x.push_row(program.x.row(r));",
     "    for r in (0..program.x.rows).chain(0..program.x.rows) {\n        x.push_row(program.x.row(r));",
     "rust", "input_arms::tests::the_joined_arm_holds_each_block_once_the_program_first"),
    ("2d", "a joined arm that keeps a K31 copy of a program column", ARMS,
     "            if held.contains(ds.x.row(r)) {",
     "            if ds.x.rows == 0 && held.contains(ds.x.row(r)) {",
     "rust", "input_arms::tests::a_k31_row_equal_to_a_program_column_enters_once_as_the_programs"),
    ("2c", "the scorer passing an inputs file whose rows double a block", SCORE,
     'inputs["arm_rows"] == inputs["program_rows"] + (inputs["k31_rows"] - dropped if joined else 0)',
     'inputs["arm_rows"] >= inputs["program_rows"]',
     "py", f"{T}.Checks.test_a_joined_arm_that_drops_or_doubles_a_block_fails"),
    ("3", "a program column built from the label", ARMS,
     "        x.push_row(program.x.row(r));\n        held.insert",
     "        let leak: Vec<u64> = program.x.row(r).iter().zip(ds.y.words()).map(|(a, b)| a | b).collect();\n"
     "        x.push_row(&leak);\n        held.insert",
     "rust", "input_arms::tests::an_arm_matrix_does_not_change_with_the_labels"),
    ("4a", "an arm scored on folds other than the default's", SCORE,
     "base is not None and fold_shape(base) == fold_shape(doc)",
     "base is not None",
     "py", f"{T}.Checks.test_an_arm_on_other_folds_fails"),
    ("4b", "an arm paired with the default at another seed", SCORE,
     'base_file = base_of.get((e["dataset"], e["split"], e["seed"]))',
     'base_file = base_of.get((e["dataset"], e["split"], str(RELEASE_SEED)))',
     "py", f"{T}.Checks.test_an_arm_is_held_to_the_default_at_its_own_seed"),
    ("5a", "a twin at a draw's seed", SCORE,
     "seed in twin_seeds and seed not in draw_seeds",
     "seed in twin_seeds",
     "py", f"{T}.Checks.test_a_twin_at_a_draw_seed_fails"),
    ("5b", "a twin listed at the draws' seeds", "scripts/c3-inputs/expect.py",
     "    if set(seeds) & set(twin):",
     "    if False:",
     "py", "tests.c3_inputs.test_lists.Expect.test_a_twin_at_a_draw_seed_is_refused"),
    ("5c", "a twin fitted with an input arm", SCORE,
     "inputs is None and stray_inputs == 0",
     "inputs is None",
     "py", f"{T}.Checks.test_a_twin_with_an_inputs_file_fails"),
    ("6", "a default not held to the release", SCORE,
     "out = C9.check(e9, doc, release if is_default else None, served_learner if is_default else None)",
     "out = C9.check(e9, doc, None, None)",
     "py", f"{T}.Checks.test_a_default_that_differs_from_the_release_fails"),
    ("7", "a difference taken as the default's value less the arm's", SCORE,
     'diffs.append(r12(float(a[m]) - float(b[m])))',
     'diffs.append(r12(float(b[m]) - float(a[m])))',
     "py", f"{T}.Count.test_every_difference_is_the_arms_value_less_the_defaults"),
    ("8", "Holm over one arm's five learners", SCORE,
     '    for r, p in zip(family, C9.holm([r["sign_p"] for r in family])):\n        r["holm_p"] = p',
     '    for a in FAMILY:\n        part = [r for r in family if r["pair"] == f"{a}-{DEFAULT}"]\n'
     '        for r, p in zip(part, C9.holm([r["sign_p"] for r in part])):\n            r["holm_p"] = p',
     "py", f"{T}.Count.test_holm_runs_over_the_thirty_lineage_pairs_of_the_primary_measure"),
    ("9", "a count taken on part of the family", SCORE,
     "    if missing:\n        raise SystemExit",
     "    if False:\n        raise SystemExit",
     "py", f"{T}.Count.test_the_count_refuses_part_of_the_family"),
    ("10", "an input arm that runs without --force, so a kept report could skip it", "src/main.rs",
     'requires_all = ["exp_pick", "exp_program_tables", "force"]',
     'requires_all = ["exp_pick", "exp_program_tables"]',
     "rust-bin", "tests::learn_takes_the_c3_input_arm_only_when_asked"),
]


def _c9():
    spec = importlib.util.spec_from_file_location("c9_faults_for_c3", HERE.parent / "c9-pick" / "planted_faults.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main(argv):
    if len(argv) != 1:
        print(__doc__.strip().splitlines()[2], file=sys.stderr)
        return 2
    c9 = _c9()
    run_test = c9.run_test

    def run(tree, kind, name, log, env):
        if kind != "rust-bin":
            return run_test(tree, kind, name, log, env)
        # A test of the binary's argument parser lives in src/main.rs.
        import subprocess
        cmd = ["cargo", "test", "--release", "--locked", "--bin", "grm", "--", "--exact", name]
        with open(log, "w") as f:
            p = subprocess.run(cmd, cwd=tree, env=env, stdout=f, stderr=subprocess.STDOUT)
        text = Path(log).read_text(errors="replace")
        built = not any(line.startswith("error[E") or line.startswith("error: could not compile")
                        for line in text.splitlines())
        return p.returncode == 0, built and (" 1 passed" in text or " 0 passed" in text)

    c9.FAULTS = FAULTS
    c9.run_test = run
    c9.REPO = HERE.parents[1]
    return c9.main(argv)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
