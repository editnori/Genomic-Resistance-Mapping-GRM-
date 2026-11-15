"""What the learner cross-check suite shares: the kit's modules, the committed result and the tables its README holds.

The kit's folder, scripts/learner-crosscheck, is no package name, so its modules are loaded by path. Only the modules
that need no library are loaded: panels, scm_reference, exact, table and check_env."""

import ast
import collections
import importlib.util
import json
import os
import pathlib
import subprocess

ROOT = pathlib.Path(__file__).resolve().parents[2]
KIT = ROOT / "scripts" / "learner-crosscheck"
REFERENCE = ROOT / "reference" / "learner-crosscheck"
FIXTURES = ROOT / "tests" / "data" / "lasso-stopping-fixtures.txt"

NAMES = {
    ("scm", "fits"): "Set covering machine",
    ("cart", "nodes"): "Decision tree",
    ("rf", "trees"): "Random forest",
    ("gbt", "tree by tree"): "Boosted trees, tree by tree",
    ("gbt", "end to end"): "Boosted trees, end to end",
}


def module(name):
    spec = importlib.util.spec_from_file_location(f"learner_crosscheck_{name}", KIT / f"{name}.py")
    loaded = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(loaded)
    return loaded


table = module("table")
panels = module("panels")


def git(*args):
    """A git command's output in the repository, stripped. It fails the test that calls it if git does."""
    env = {"PATH": os.environ.get("PATH", "/usr/bin:/bin"), "HOME": os.environ.get("HOME", "/"), "LC_ALL": "C"}
    done = subprocess.run(["git", "-C", str(ROOT), *args], capture_output=True, text=True, env=env, check=False)
    if done.returncode:
        raise AssertionError(f"git {' '.join(args)} exited {done.returncode}: {done.stderr.strip()}")
    return done.stdout.strip()


def load(name):
    return json.loads((REFERENCE / name).read_text(encoding="utf-8"))


def summary():
    return table.read_summary(REFERENCE / "summary.tsv")


def differences():
    return table.read_differences(REFERENCE / "differences.tsv")


def trees():
    """end-to-end-trees.tsv, the record added on 15 November 2026: each end-to-end boosted tree that differs."""
    return table.read_trees(REFERENCE / "end-to-end-trees.tsv")


def constants(name):
    """The top-level names a kit file sets to a literal, or to another such name, read without running the file."""
    found = {}
    for node in ast.parse((KIT / f"{name}.py").read_text(encoding="utf-8")).body:
        if isinstance(node, ast.Assign) and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name):
            if isinstance(node.value, ast.Name) and node.value.id in found:
                found[node.targets[0].id] = found[node.value.id]
                continue
            try:
                found[node.targets[0].id] = ast.literal_eval(node.value)
            except ValueError:
                pass
    return found


def detail_sum(check, key):
    """A count of detail.json summed over the panels that have it."""
    return sum(panel[check][key] for panel in load("detail.json").values() if key in panel.get(check, {}))


def totals_table():
    """The README's table of totals, from totals.json."""
    lines = [
        "| Learner and check | Unit | Compared | Matched | Shown cause | No shown cause | Panels matched | Panels with a shown cause | Panels with none shown | Panels not run |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in load("totals.json"):
        by = row["panels"]
        lines.append(
            f"| {NAMES[(row['learner'], row['check'])]} | {row['units']} | {row['compared']:,} | {row['units_matched']:,} | {row['units_explained']:,} "
            f"| {row['units_unexplained']:,} | {by['matched']} | {by['explained']} | {by['unexplained']} | {by['not applicable']} |"
        )
    return lines


def match_table():
    """The README's table of what a match holds although the two sides differ, from the result files."""
    totals = {(row["learner"], row["check"]): row for row in load("totals.json")}
    matched = collections.defaultdict(list)
    for row in summary():
        if row["outcome"] == "matched":
            matched[(row["learner"], row["check"])].append(row)
    end = [row["panel"] for row in matched[("gbt", "end to end")]]
    record = [row for row in trees() if row["panel"] in end]
    kept = sum(load("exploratory.json")["end_to_end_trees"][panel]["trees_both_keep"] for panel in end)
    tie = "scikit-learn split on another feature of the same exact tie for the largest gain"
    rows = (
        ("cart", "nodes", f"{totals[('cart', 'nodes')]['units_matched']:,} nodes", detail_sum("cart nodes", "nodes_other_member_of_a_tie"), tie),
        ("rf", "trees", f"{detail_sum('rf trees', 'nodes_agree'):,} nodes that agreed", detail_sum("rf trees", "nodes_other_member_of_a_tie"), tie),
        ("gbt", "tree by tree", f"{len(matched[('gbt', 'tree by tree')])} panels", sum(row["compared"] == 0 for row in matched[("gbt", "tree by tree")]),
         "neither side kept a tree, so no tree was compared"),
        ("gbt", "end to end", f"{len(end)} panels", len({row["panel"] for row in record}),
         f"{len(record)} of the {kept:,} trees that both sides keep there are not identical (`end-to-end-trees.tsv`)"),
    )  # fmt: skip
    lines = ["| Learner and check | Matched | Of those | In which |", "|---|---|---:|---|"]
    for learner, check, units, count, what in rows:
        lines.append(f"| {NAMES[(learner, check)]} | {units} | {count:,} | {what} |")
    return lines


def panel_table():
    """The README's table by panel, from summary.tsv: the units that matched of those compared."""
    lines = ["| Panel | Samples | Features | " + " | ".join(NAMES[check[:2]] for check in table.CHECKS) + " |", "|---|---:|---:|" + "---:|" * len(table.CHECKS)]
    rows = collections.OrderedDict()
    for row in summary():
        rows.setdefault((row["panel"], row["samples"], row["features"]), []).append(row)
    for (panel, samples, features), checks in rows.items():
        cells = ["not run" if row["outcome"] == "not applicable" else f"{row['matched']:,} of {row['compared']:,}" for row in checks]
        lines.append(f"| `{panel}` | {samples} | {features} | " + " | ".join(cells) + " |")
    return lines


def cause_table():
    """The README's table of causes, from differences.tsv: the rows and the panels of each cause."""
    rows, where = collections.Counter(), collections.defaultdict(set)
    for row in differences():
        key = (row["learner"], row["check"], row["cause"])
        rows[key] += 1
        where[key].add(row["panel"])
    lines = ["| Learner and check | Cause | Difference rows | Panels |", "|---|---|---:|---:|"]
    for learner, check, _ in table.CHECKS:
        for key in sorted(k for k in rows if k[:2] == (learner, check)):
            lines.append(f"| {NAMES[key[:2]]} | {key[2] or 'none shown'} | {rows[key]:,} | {len(where[key])} |")
    return lines


# (label, check in detail.json, key): the further counts the README gives, each a sum over the panels.
FURTHER = (
    ("Set covering machine: fits the independent implementation gave rule for rule", "scm fits", "independent_same"),
    ("Set covering machine: rules in GRM's fits", "scm fits", "rules"),
    ("Set covering machine: steps read against pyscm-ml", "scm fits", "steps_read_against_library"),
    ("Set covering machine: steps where GRM's rule was among the library's best", "scm fits", "steps_certified"),
    ("Set covering machine: of those, by the same effect only", "scm fits", "steps_certified_same_effect_only"),
    ("Set covering machine: steps where the library preferred a rule GRM does not have", "scm fits", "steps_dominated"),
    ("Set covering machine: steps of any other kind", "scm fits", "steps_other"),
    ("Set covering machine: fits pyscm-ml refused (one class)", "scm fits", "fits_library_refuses"),
    ("Set covering machine: fits pyscm-ml ran freely", "scm fits", "library_free_fits"),
    ("Set covering machine: free fits with GRM's rules", "scm fits", "library_free_same_rules"),
    ("Set covering machine: free fits with GRM's calls on every sample", "scm fits", "library_free_same_calls"),
    ("Decision tree: nodes", "cart nodes", "nodes"),
    ("Decision tree: nodes searched for a split", "cart nodes", "nodes_searched"),
    ("Decision tree: both made a leaf", "cart nodes", "nodes_both_leaf"),
    ("Decision tree: both split on the same feature", "cart nodes", "nodes_same_feature"),
    ("Decision tree: they split on two features of one exact tie, GRM on the lowest", "cart nodes", "nodes_other_member_of_a_tie"),
    ("Decision tree: nodes with an exact tie for the largest gain", "cart nodes", "nodes_with_an_exact_tie"),
    ("Decision tree: stop rules broken", "cart nodes", "stop_rule_failures"),
    ("Decision tree: node records wrong", "cart nodes", "node_record_failures"),
    ("Decision tree: panels where scikit-learn's own tree was the same tree", "cart nodes", "free_tree_same"),
    ("Decision tree: places where its own tree parted from GRM's at an exact tie", "cart nodes", "free_tree_parts_at_exact_ties"),
    ("Decision tree: places where it parted anywhere else", "cart nodes", "free_tree_parts_elsewhere"),
    ("Decision tree: panels with the same tree at scikit-learn's default floor of 0", "cart nodes", "default_floor_same"),
    ("Random forest: nodes searched for a split", "rf trees", "nodes_searched"),
    ("Random forest: nodes that agreed", "rf trees", "nodes_agree"),
    ("Random forest: both made a leaf", "rf trees", "nodes_both_leaf"),
    ("Random forest: both split on the same feature", "rf trees", "nodes_same_feature"),
    ("Random forest: they split on two features of one exact tie, GRM on the lowest", "rf trees", "nodes_other_member_of_a_tie"),
    ("Random forest: nodes with an exact tie for the largest gain", "rf trees", "nodes_with_an_exact_tie"),
    ("Boosted trees, tree by tree: rounds grown", "gbt tree by tree", "rounds_grown"),
    ("Boosted trees, tree by tree: trees kept", "gbt tree by tree", "trees_kept"),
    ("Boosted trees, end to end: trees GRM kept", "gbt end to end", "grm_trees"),
    ("Boosted trees, end to end: trees XGBoost kept", "gbt end to end", "xgboost_trees"),
    ("Boosted trees, end to end: trees GRM kept at the production column share", "gbt end to end", "grm_trees_at_production_colsample"),
)


def further_table():
    lines = ["| Count, over the 16 panels | Value |", "|---|---:|"]
    for label, check, key in FURTHER:
        lines.append(f"| {label} | {detail_sum(check, key):,} |")
    return lines


def exploratory_table():
    e = load("exploratory.json")
    ties, splits, trees = e["rounding_ties"], e["equal_splits_total"], e["end_to_end_trees_total"]
    gains = e["stored_gain"]
    rows = (
        ("Decision tree: split nodes", gains["decision_tree"]["split_nodes"]),
        ("Decision tree: split nodes whose stored gain is the gain the kit recomputes", gains["decision_tree"]["stored_gain_is_the_recomputed_one"]),
        ("Random forest: split nodes", gains["forest"]["split_nodes"]),
        ("Random forest: split nodes whose stored gain is the gain the kit recomputes", gains["forest"]["stored_gain_is_the_recomputed_one"]),
        ("Decision tree: nodes where rounding broke an exact tie", ties["cart"]["nodes"]),
        ("Decision tree: of those, the two features are present on opposite samples of the node", ties["cart"]["features_present_where_the_other_is_absent"]),
        ("Random forest: nodes where rounding broke an exact tie", ties["rf"]["nodes"]),
        ("Random forest: of those, the two features are present on opposite samples of the node", ties["rf"]["features_present_where_the_other_is_absent"]),
        ("Boosted trees: equal splits where the two features are present on the same samples", splits["same_way_round"]),
        ("Boosted trees: equal splits where they are present on opposite samples", splits["other_way_round"]),
        ("Boosted trees: equal splits where GRM names the higher-numbered feature", splits["grm_names_the_higher_feature"]),
        ("Boosted trees, end to end: trees both sides keep", trees["trees_both_keep"]),
        ("Boosted trees, end to end: of those, identical to XGBoost's tree of the same round", trees["identical"]),
    )
    return ["| Exploratory count, over the 16 panels | Value |", "|---|---:|"] + [f"| {label} | {value:,} |" for label, value in rows]


def tree_record_table():
    """The README's table of end-to-end-trees.tsv: per panel and cause, the trees and the margin gaps before them."""
    gaps = collections.defaultdict(list)
    for row in trees():
        gaps[(row["panel"], row["cause"] or "none shown")].append(row["fitting_margin_gap"])
    order = [panel["name"] for panel in panels.read(FIXTURES)]
    lines = ["| Panel | Cause of the tree | Trees | Smallest margin gap | Largest margin gap |", "|---|---|---:|---:|---:|"]
    for panel, cause in sorted(gaps, key=lambda key: (order.index(key[0]), key[1])):
        found = gaps[(panel, cause)]
        lines.append(f"| `{panel}` | {cause} | {len(found):,} | {min(found):.3g} | {max(found):.3g} |")
    return lines


TABLES = (totals_table, match_table, panel_table, cause_table, tree_record_table, further_table, exploratory_table)

if __name__ == "__main__":
    for make in TABLES:
        print("\n".join(make()), end="\n\n")
