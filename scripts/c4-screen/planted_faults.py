"""Plant each fault of docs/screen-arms.md's table and run the test meant to catch it.

  python3 -I planted_faults.py <out dir>

scripts/c3-inputs/planted_faults.py's runner, imported, with this block's
faults: it copies this checkout's last commit to <out dir>/tree (git
archive) and, one fault at a time, changes one piece of text there, runs
the one test named for it, and puts the text back. A fault is caught when
its test fails and the tree still builds; each test must also pass with no
fault planted. Writes <out dir>/faults.tsv (fault, what it plants, file,
test, result) and <out dir>/logs/, prints one line per fault, and exits 1
when a fault is not caught. Every test it runs is SYNTHETIC; it reads no
data.
"""

import importlib.util
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ARM = "src/learn/screen_arm.rs"
SCORE = "scripts/c4-screen/score.py"
EXPECT = "scripts/c4-screen/expect.py"
R = "learn::screen_arm::tests"
T = "tests.c4_screen.test_score"
L = "tests.c4_screen.test_lists"

FAULTS = [
    ("1", "a screen that reads an outer test genome: the strata count every genome, not the fitting ones", ARM,
     "    let mut by: BTreeMap<u32, (u64, u64)> = BTreeMap::new();\n    for i in fit.iter_ones() {",
     "    let mut by: BTreeMap<u32, (u64, u64)> = BTreeMap::new();\n    for i in 0..fit.len() {",
     "rust", f"{R}::no_sample_outside_the_fit_changes_the_screen"),
    ("2", "strata taken from labels: a resistant genome's stratum moved by its label", ARM,
     "    let mut by: BTreeMap<u32, (u64, u64)> = BTreeMap::new();\n    for i in fit.iter_ones() {",
     "    let lineage: Vec<u32> = (0..fit.len())\n"
     "        .map(|i| lineage[i] + 1_000 * (y.contains(i) && i % 2 == 0) as u32)\n        .collect();\n"
     "    let mut by: BTreeMap<u32, (u64, u64)> = BTreeMap::new();\n    for i in fit.iter_ones() {",
     "rust", f"{R}::the_strata_do_not_change_with_the_labels"),
    ("3", "a stratum of one class counted as informative", ARM,
     "        if n >= 2 && r > 0 && r < n {",
     "        if n >= 2 {",
     "rust", f"{R}::a_stratum_of_one_class_or_one_sample_changes_no_statistic"),
    ("4a", "the floor arm run at another budget", "src/main.rs",
     "        floor.is_none() || screened_features == DEFAULT_SCREENED_FEATURES,",
     "        floor.is_none() || screened_features > 0,",
     "rust-bin", "tests::learn_takes_the_c4_screen_and_floor_arms_only_when_asked"),
    ("4b", "the scorer passing a floor arm's inputs file at another budget", SCORE,
     '(inputs.get("rule"), inputs.get("budget")) == (rule, budget)',
     'inputs.get("rule") == rule',
     "py", f"{T}.Checks.test_a_floor_arm_at_another_budget_fails"),
    ("5a", "an arm fitted with the default's screen", ARM,
     "                cmh_screen(x, y, fit, &ctx.lineage, k)\n            });",
     "                let _ = &ctx.lineage;\n                Screened {\n                    feats: super::screen(x, y, fit, k),\n"
     "                    kept_at_zero: 0,\n                    strata: 0,\n                    informative: 0,\n"
     "                }\n            });",
     "rust", f"{R}::a_fit_under_the_adjusted_arm_takes_its_screen_and_with_no_arm_the_chi_square"),
    ("5b", "the scorer passing an arm whose screens are the default's", SCORE,
     'all(s["rule"] == rule for s in screens)',
     'all(s["rule"] in ("cmh", "chi-square") for s in screens)',
     "py", f"{T}.Checks.test_an_arm_screened_by_the_defaults_rule_fails"),
    ("6a", "the twin checked as an arm, not as the default at seeds of its own", SCORE,
     "    if arm in (DEFAULT, TWIN):\n        out = C3.check(",
     "    if arm == DEFAULT:\n        out = C3.check(",
     "py", f"{T}.Checks.test_the_twin_passes_as_the_default_at_its_own_seeds_and_fails_with_an_inputs_file"),
    ("6b", "the twin listed with an arm's inputs file", EXPECT,
     'f"{TWIN}/{d}.lineage-s{s}.pick.json", "-"])',
     'f"{TWIN}/{d}.lineage-s{s}.pick.json", inputs_name(TWIN, d, "drawn", s)])',
     "py", f"{L}.Expect.test_only_an_arm_names_an_inputs_file"),
    ("7", "a missing draw read as a zero", SCORE,
     '"scored": 0, "constant": 0, "why": "missing"})',
     '"scored": 1, "constant": 0, "why": "", **{m: 0.0 for m in MEASURES}})',
     "py", f"{T}.Scores.test_a_dataset_with_a_missing_draw_is_not_scored"),
    ("8a", "a study result listed as a lineage draw", EXPECT,
     '            out.append([d, species, arm, "study", STUDY_SCHEME, RELEASE_SEED, "study",',
     '            out.append([d, species, arm, "lineage", STUDY_SCHEME, RELEASE_SEED, "study",',
     "py", f"{L}.Expect.test_each_result_is_on_its_own_split"),
    ("8b", "the scorer passing a study result read as a lineage result", SCORE,
     '(exp["split"], exp["scheme"]) in (("lineage", "lineage"), ("study", STUDY_SCHEME))\n       and ',
     "",
     "py", f"{T}.Checks.test_a_study_result_read_as_a_lineage_result_fails"),
    ("9", "Holm over each arm's own learners, not the 16 pairs", SCORE,
     '    for r, p in zip(fam, C9.holm([r["sign_p"] for r in fam])):\n        r["holm_p"] = p',
     '    for a in FAMILY:\n        part = [r for r in fam if arm_of(r["pair"]) == a]\n'
     '        for r, p in zip(part, C9.holm([r["sign_p"] for r in part])):\n            r["holm_p"] = p',
     "py", f"{T}.Count.test_holm_over_sixteen_is_held_to_the_plain_step_down"),
    ("10", "a screen arm's scm and cart counted", SCORE,
     "    return SCREENED if arm in SCREEN_ARMS else LEARNERS",
     "    return LEARNERS",
     "py", f"{T}.Count.test_a_screen_arms_scm_and_cart_are_not_counted"),
    ("11", "a floor arm unlike the default below 800 genomes passed", SCORE,
     "untimed(doc) == untimed(base))",
     "C3.fold_shape(doc) == C3.fold_shape(base))",
     "py", f"{T}.Checks.test_a_floor_arm_below_800_genomes_unlike_the_default_fails_and_at_800_may_differ"),
    ("12", "a screen arm whose scm or cart differs from the default's passed", SCORE,
     "same_calls(base, doc, [l for l in LEARNERS if l not in SCREENED]))",
     "same_calls(base, doc, []))",
     "py", f"{T}.Checks.test_a_screen_arm_whose_scm_or_cart_differs_from_the_default_fails"),
    ("13", "the floor arm run without --force, so a kept report could skip it", "src/main.rs",
     '            value_parser = RangedU64ValueParser::<u32>::new().range(1..),\n'
     '            requires_all = ["exp_pick", "force"],',
     '            value_parser = RangedU64ValueParser::<u32>::new().range(1..),\n'
     '            requires_all = ["exp_pick"],',
     "rust-bin", "tests::learn_takes_the_c4_screen_and_floor_arms_only_when_asked"),
]


def _c3():
    spec = importlib.util.spec_from_file_location("c3_faults_for_c4", HERE.parent / "c3-inputs" / "planted_faults.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main(argv):
    c3 = _c3()
    # C3's runner plants its module's FAULTS: this block's, for the call.
    c3.FAULTS = FAULTS
    return c3.main(argv)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
