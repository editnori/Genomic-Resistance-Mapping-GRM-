"""What the learner cross-check suite shares: the kit's modules, the committed result and the tables its README holds.

The kit's folder, scripts/learner-crosscheck, is no package name, so its modules are loaded by path. Only the modules
that need no library are loaded: panels, scm_reference, exact, table and check_env."""

import collections
import importlib.util
import json
import pathlib

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


def load(name):
    return json.loads((REFERENCE / name).read_text(encoding="utf-8"))


def summary():
    return table.read_summary(REFERENCE / "summary.tsv")


def differences():
    return table.read_differences(REFERENCE / "differences.tsv")


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


TABLES = (totals_table, panel_table, cause_table, further_table, exploratory_table)

if __name__ == "__main__":
    for make in TABLES:
        print("\n".join(make()), end="\n\n")
