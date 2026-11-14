"""The decision tree and the random forest against scikit-learn.

scikit-learn's split search is the reference, node by node: a one-split tree
fitted on the training samples a node holds, with the tree's class weights
and, for the forest, the features GRM drew at that node. exact.py settles a
tied node, where scikit-learn keeps the first best feature of a seeded random
order and GRM the lowest one.
"""

import math

import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.tree import DecisionTreeClassifier

import exact
import panels as P

GAIN_FLOOR = 1e-12  # GRM splits a node only on a gain above this (src/learn/tree.rs)
VALUE_TOLERANCE = 1e-9
SCORE_TOLERANCE = 1e-12


def dense(rows, n):
    x = np.zeros((n, len(rows)), dtype=np.float32)
    for feature, row in enumerate(rows):
        x[P.members(row), feature] = 1.0
    return x


def node_samples(tree, rows, root):
    """The samples of `root` that reach each node of a GRM tree."""
    sets = [0] * len(tree["n"])
    sets[0] = root
    stack = [0]
    while stack:
        node = stack.pop()
        if tree["feature"][node] is None:
            continue
        row = rows[tree["feature"][node]]
        sets[tree["absent"][node]] = sets[node] & ~row
        sets[tree["present"][node]] = sets[node] & row
        stack += [tree["absent"][node], tree["present"][node]]
    return sets


def float_gain(weights, p, n, pp, pn):
    """A split's gain in float64, operation for operation as `fit_tree` computes it."""

    def gini(a, b):
        t = a + b
        if t == 0.0:
            return 0.0
        q = a / t
        return 2.0 * q * (1.0 - q)

    total = float(weights.n_resistant + weights.n_susceptible)
    w_pos = total / (2.0 * max(float(weights.n_resistant), 1.0))
    w_neg = total / (2.0 * max(float(weights.n_susceptible), 1.0))
    np_, nn = p * w_pos, n * w_neg
    parent = (np_ + nn) * gini(np_, nn)
    lp, ln = pp * w_pos, pn * w_neg
    rp, rn = np_ - lp, nn - ln
    return parent - (lp + ln) * gini(lp, ln) - (rp + rn) * gini(rp, rn)


def stump(x, y, weights, samples, candidates):
    """scikit-learn's one-split tree on a node: None for a leaf, else
    (feature, resistant share where absent, resistant share where present)."""
    members = P.members(samples)
    labels = np.array([(y >> i) & 1 for i in members], dtype=np.int64)
    weight = np.where(labels == 1, float(weights.resistant), float(weights.susceptible))
    model = DecisionTreeClassifier(
        criterion="gini",
        splitter="best",
        max_depth=1,
        min_samples_split=2,
        min_samples_leaf=1,
        max_features=None,
        random_state=0,
        min_impurity_decrease=GAIN_FLOOR / float(weight.sum()),
    )
    columns = sorted(candidates)
    model.fit(x[np.ix_(members, columns)], labels, sample_weight=weight)
    t = model.tree_
    if t.node_count == 1:
        return None
    if t.node_count != 3 or t.threshold[0] != 0.5:
        raise AssertionError("a one-split tree on 0/1 features splits once, at 0.5")
    classes = model.classes_.tolist()
    share = lambda node: float(t.value[node, 0, classes.index(1)]) if 1 in classes else 0.0
    return columns[int(t.feature[0])], share(t.children_left[0]), share(t.children_right[0])


def check_tree(tree, rows, x, y, train, max_depth, draws=None):
    """One GRM tree, node by node. `draws` maps a node to the features offered
    to it; without it every feature is offered to every node."""
    resistant, susceptible = train & y, train & ~y
    weights = exact.Weights(P.count(resistant), P.count(susceptible))
    sets = node_samples(tree, rows, train)
    out = {
        "nodes": len(sets),
        "searched": 0,
        "agree": 0,
        "both_leaf": 0,
        "same_feature": 0,
        "tie_other_member": 0,
        "tied_nodes": 0,
        "stop_rule_failures": [],
        "value_failures": [],
        "disagreements": [],
    }
    all_features = range(len(rows))
    for node, samples in enumerate(sets):
        pos, neg = samples & y, samples & ~y
        p, n = P.count(pos), P.count(neg)
        depth = tree["depth"][node]
        feature = tree["feature"][node]
        if (tree["n"][node], tree["n_resistant"][node]) != (p + n, p):
            out["value_failures"].append(f"node {node}: counts")
        if not math.isclose(tree["value"][node], float(weights.value(p, n)), rel_tol=0, abs_tol=SCORE_TOLERANCE):
            out["value_failures"].append(f"node {node}: value")
        if feature is not None and tree["depth"][tree["absent"][node]] != depth + 1:
            out["value_failures"].append(f"node {node}: child depth")
        if p + n == 0:
            out["value_failures"].append(f"node {node}: no training sample reaches it")
            continue
        must_stop = depth >= max_depth or p == 0 or n == 0 or p + n < 2
        if must_stop and feature is not None:
            out["stop_rule_failures"].append(f"node {node}: split at depth {depth} with {p} resistant and {n} susceptible samples")
        if depth >= max_depth:
            continue  # the depth limit is no split search; the free tree covers it
        if draws is None:
            candidates = all_features
        elif node in draws:
            candidates = draws[node]
        elif must_stop:
            candidates = all_features  # a pure node draws nothing; any feature set gives a leaf
        else:
            out["stop_rule_failures"].append(f"node {node}: searched but has no draw")
            continue
        out["searched"] += 1
        top, at = exact.best(weights, pos, neg, rows, candidates)
        if len(at) > 1 and top is not None and top > 0:
            out["tied_nodes"] += 1
        if feature is None and top is not None and top > 0:
            out["stop_rule_failures"].append(f"node {node}: a leaf, and feature {at[0]} has the positive gain {float(top):.6g}")
        reference = stump(x, y, weights, samples, candidates)
        if feature is None and reference is None:
            out["agree"] += 1
            out["both_leaf"] += 1
            continue
        if feature is None or reference is None:
            out["disagreements"].append(
                {
                    "node": node,
                    "what": f"GRM {'splits on feature ' + str(feature) if feature is not None else 'makes a leaf'}, scikit-learn "
                    f"{'splits on feature ' + str(reference[0]) if reference is not None else 'makes a leaf'}",
                    "cause": "",
                    "evidence": f"largest exact gain {float(top) if top is not None else 'none'}",
                }
            )
            continue
        theirs, absent_share, present_share = reference
        in_set = theirs in at
        lowest = bool(at) and feature == at[0]
        mine_present, theirs_present = samples & rows[feature], samples & rows[theirs]
        values = True
        if theirs_present == mine_present or theirs_present == samples & ~mine_present:
            if theirs_present != mine_present:
                absent_share, present_share = present_share, absent_share
            values = (
                abs(tree["value"][tree["absent"][node]] - absent_share) <= VALUE_TOLERANCE
                and abs(tree["value"][tree["present"][node]] - present_share) <= VALUE_TOLERANCE
            )
        if in_set and lowest and values:
            out["agree"] += 1
            out["same_feature" if theirs == feature else "tie_other_member"] += 1
            continue
        what, cause, evidence = [], "", ""
        if not in_set:
            what.append(f"scikit-learn's feature {theirs} has not the largest exact gain")
        if not lowest:
            what.append(f"GRM's feature {feature} is not the lowest feature with the largest exact gain")
            if feature in at:
                def as_grm(f):
                    return float_gain(weights, p, n, P.count(pos & rows[f]), P.count(neg & rows[f]))

                if as_grm(feature) != as_grm(at[0]):
                    cause = "float64 rounding broke an exact tie: GRM's own gains for the two features differ"
                    evidence = f"features {at[0]} and {feature}: exact gains equal, float64 gains {as_grm(at[0])!r} and {as_grm(feature)!r}"
                    if samples & rows[at[0]] == samples & ~rows[feature]:
                        evidence += "; each is present on the node's samples exactly where the other is absent"
        if not values:
            what.append("child values differ")
        out["disagreements"].append({"node": node, "what": "; ".join(what), "cause": cause, "evidence": evidence})
    return out


def walk_scores(tree, rows, n, max_depth=None):
    """Each sample's score: the value of the leaf it reaches, nodes at `max_depth` being leaves."""
    scores = [0.0] * n
    stack = [(0, (1 << n) - 1)]
    while stack:
        node, samples = stack.pop()
        feature = tree["feature"][node]
        if feature is None or (max_depth is not None and tree["depth"][node] >= max_depth):
            for i in P.members(samples):
                scores[i] = tree["value"][node]
            continue
        row = rows[feature]
        stack.append((tree["absent"][node], samples & ~row))
        stack.append((tree["present"][node], samples & row))
    return scores


def largest_gap(a, b):
    return max((abs(u - v) for u, v in zip(a, b)), default=0.0)


def free_tree(tree, rows, x, y, train, max_depth, floor):
    """scikit-learn's own tree, grown freely, walked beside GRM's by the
    training samples each node holds."""
    members = P.members(train)
    labels = np.array([(y >> i) & 1 for i in members], dtype=np.int64)
    model = DecisionTreeClassifier(
        criterion="gini",
        splitter="best",
        max_depth=max_depth,
        min_samples_split=2,
        min_samples_leaf=1,
        max_features=None,
        random_state=0,
        class_weight="balanced",
        min_impurity_decrease=floor / len(members) if floor else 0.0,
    )
    model.fit(x[members], labels)
    t = model.tree_
    classes = model.classes_.tolist()
    weights = exact.Weights(P.count(train & y), P.count(train & ~y))
    out = {"reference_nodes": int(t.node_count), "same_feature": 0, "equal_split": 0, "leaves": 0, "parts": [], "value_gap": 0.0}
    stack = [(0, 0, train)]
    while stack:
        mine, theirs, samples = stack.pop()
        feature = tree["feature"][mine]
        their_leaf = t.children_left[theirs] == -1
        if feature is None and their_leaf:
            out["leaves"] += 1
            share = float(t.value[theirs, 0, classes.index(1)]) if 1 in classes else 0.0
            out["value_gap"] = max(out["value_gap"], abs(share - tree["value"][mine]))
            continue
        pos, neg = samples & y, samples & ~y
        if feature is None or their_leaf:
            top, at = exact.best(weights, pos, neg, rows, range(len(rows)))
            out["parts"].append(
                {
                    "node": mine,
                    "what": f"GRM {'makes a leaf' if feature is None else 'splits on feature ' + str(feature)}, scikit-learn "
                    f"{'makes a leaf' if their_leaf else 'splits on feature ' + str(int(t.feature[theirs]))}",
                    "equal_gain": False,
                    "largest_exact_gain": float(top) if top is not None else None,
                }
            )
            continue
        other = int(t.feature[theirs])
        mine_present, their_present = samples & rows[feature], samples & rows[other]
        left, right = int(t.children_left[theirs]), int(t.children_right[theirs])
        if their_present == mine_present:
            out["same_feature" if other == feature else "equal_split"] += 1
            stack += [(tree["absent"][mine], left, samples & ~mine_present), (tree["present"][mine], right, mine_present)]
        elif their_present == samples & ~mine_present:
            out["equal_split"] += 1
            stack += [(tree["absent"][mine], right, samples & ~mine_present), (tree["present"][mine], left, mine_present)]
        else:
            equal = exact.gain(weights, pos, neg, rows[feature]) == exact.gain(weights, pos, neg, rows[other])
            out["parts"].append({"node": mine, "what": f"GRM splits on feature {feature}, scikit-learn on feature {other}", "equal_gain": equal})
    out["same_tree"] = not out["parts"]
    return out


def check_decision_tree(panel, grm):
    n, rows, y, train = panel["n"], panel["rows"], panel["y"], panel["fit"]
    x = dense(rows, n)
    tree, max_depth = grm["tree"], grm["max_depth"]
    nodes = check_tree(tree, rows, x, y, train, max_depth)
    score_gap = max(largest_gap(grm["scores_by_depth"][d - 1], walk_scores(tree, rows, n, d)) for d in range(1, max_depth + 1))
    free = free_tree(tree, rows, x, y, train, max_depth, GAIN_FLOOR)
    default = free_tree(tree, rows, x, y, train, max_depth, 0.0)
    differences = [{"unit": f"node {d['node']}", "what": d["what"], "cause": d["cause"], "evidence": d["evidence"]} for d in nodes["disagreements"]]
    for text in nodes["stop_rule_failures"] + nodes["value_failures"]:
        differences.append({"unit": text.split(":")[0], "what": "stop rule or node record: " + text, "cause": "", "evidence": ""})
    if score_gap > SCORE_TOLERANCE:
        differences.append({"unit": "scores", "what": f"scores differ from the tree's leaf values by up to {score_gap:.3g}", "cause": "", "evidence": ""})
    for part in free["parts"]:
        if not part["equal_gain"]:
            differences.append({"unit": f"free tree node {part['node']}", "what": part["what"], "cause": "", "evidence": "the two splits have unequal exact gains"})
    if free["value_gap"] > VALUE_TOLERANCE:
        differences.append({"unit": "free tree", "what": f"leaf values differ by up to {free['value_gap']:.3g}", "cause": "", "evidence": ""})
    unexplained = sum(not d["cause"] for d in differences)
    return {
        "outcome": "matched" if not differences else "unexplained" if unexplained else "explained",
        "units": "nodes",
        "compared": nodes["searched"],
        "matched": nodes["agree"],
        "explained": sum(bool(d["cause"]) for d in nodes["disagreements"]),
        "unexplained": sum(not d["cause"] for d in nodes["disagreements"]),
        "detail": {
            "nodes": nodes["nodes"],
            "nodes_searched": nodes["searched"],
            "nodes_both_leaf": nodes["both_leaf"],
            "nodes_same_feature": nodes["same_feature"],
            "nodes_other_member_of_a_tie": nodes["tie_other_member"],
            "nodes_with_an_exact_tie": nodes["tied_nodes"],
            "stop_rule_failures": len(nodes["stop_rule_failures"]),
            "node_record_failures": len(nodes["value_failures"]),
            "largest_score_gap": score_gap,
            "free_tree_same": free["same_tree"],
            "free_tree_nodes_same_feature": free["same_feature"],
            "free_tree_nodes_equal_split": free["equal_split"],
            "free_tree_parts_at_exact_ties": sum(part["equal_gain"] for part in free["parts"]),
            "free_tree_parts_elsewhere": sum(not part["equal_gain"] for part in free["parts"]),
            "free_tree_reference_nodes": free["reference_nodes"],
            "default_floor_same": default["same_tree"],
            "default_floor_reference_nodes": default["reference_nodes"],
        },
        "differences": differences,
        "free_tree": free,
        "default_floor_tree": default,
    }


def library_forest(x, y, train, take, max_features, max_depth, seed):
    """scikit-learn's own forest: reported, not judged."""
    members = P.members(train)
    labels = np.array([(y >> i) & 1 for i in members], dtype=np.int64)
    model = RandomForestClassifier(
        n_estimators=200,
        criterion="gini",
        max_depth=max_depth,
        min_samples_split=2,
        min_samples_leaf=1,
        max_features=max_features,
        bootstrap=True,
        max_samples=take,
        class_weight="balanced_subsample",
        min_impurity_decrease=GAIN_FLOOR / take,
        random_state=seed,
        n_jobs=1,
    )
    model.fit(x[members], labels)
    classes = model.classes_.tolist()
    return model.predict_proba(x)[:, classes.index(1)].tolist() if 1 in classes else [0.0] * x.shape[0]


def check_forest(panel, grm):
    n, y, train = panel["n"], panel["y"], panel["fit"]
    order = grm["screen"]
    rows = [panel["rows"][f] for f in order]  # the forest names positions in the screen's order
    x = dense(rows, n)
    params = grm["params"]
    take = math.floor(P.count(train) * params["subsample"] + 0.5)  # Rust's round, half away from zero
    max_features = max(math.ceil(math.sqrt(len(rows))), 1)
    problems = []
    if sorted(order) != list(range(len(panel["rows"]))):
        problems.append("the screen's order is no ordering of the panel's features")
    if (grm["take"], grm["max_features"], len(grm["trees"])) != (take, max_features, params["n_trees"]):
        problems.append("subsample size, draw size or tree count")
    totals = {key: 0 for key in ("nodes", "searched", "agree", "both_leaf", "same_feature", "tie_other_member", "tied_nodes")}
    differences, trees_ok = [], 0
    sums = [0.0] * n
    for index, entry in enumerate(grm["trees"]):
        tree = entry["tree"]
        sample = P.decode_bits(n, entry["subsample"])
        bad = []
        if sample & ~train or P.count(sample) != take:
            bad.append({"what": f"the subsample is not {take} distinct training samples", "cause": "", "evidence": ""})
        draws = {node: features for node, features in entry["draws"]}
        for node, features in draws.items():
            if len(features) != max_features or len(set(features)) != len(features) or not all(0 <= f < len(rows) for f in features):
                bad.append({"what": f"node {node}: the draw is not {max_features} distinct features", "cause": "", "evidence": ""})
        result = check_tree(tree, rows, x, y, sample, params["max_depth"], draws)
        for key in totals:
            totals[key] += result[key]
        bad += [{"what": f"node {d['node']}: {d['what']}", "cause": d["cause"], "evidence": d["evidence"]} for d in result["disagreements"]]
        bad += [{"what": text, "cause": "", "evidence": ""} for text in result["stop_rule_failures"] + result["value_failures"]]
        trees_ok += not bad
        differences += [{"unit": f"tree {index}", **d} for d in bad]
        for i, value in enumerate(walk_scores(tree, rows, n)):
            sums[i] += value
    score_gap = largest_gap(grm["scores"], [s / len(grm["trees"]) for s in sums])
    if score_gap > SCORE_TOLERANCE:
        differences.append({"unit": "scores", "what": f"scores differ from the mean of the trees' leaf values by up to {score_gap:.3g}", "cause": "", "evidence": ""})
    differences += [{"unit": "forest", "what": text, "cause": "", "evidence": ""} for text in problems]

    library = [library_forest(x, y, train, take, max_features, params["max_depth"], seed) for seed in (0, 1)]
    calls = lambda scores: [s >= 0.5 for s in scores]
    gap = lambda a, b: sum(abs(u - v) for u, v in zip(a, b)) / n
    flips = lambda a, b: sum(u != v for u, v in zip(calls(a), calls(b)))
    bad_trees = {d["unit"] for d in differences if d["unit"].startswith("tree ")}
    explained_trees = {unit for unit in bad_trees if all(d["cause"] for d in differences if d["unit"] == unit)}
    unexplained = sum(not d["cause"] for d in differences)
    return {
        "outcome": "matched" if not differences else "unexplained" if unexplained else "explained",
        "units": "trees",
        "compared": len(grm["trees"]),
        "matched": trees_ok,
        "explained": len(explained_trees),
        "unexplained": len(bad_trees) - len(explained_trees),
        "detail": {
            "nodes": totals["nodes"],
            "nodes_searched": totals["searched"],
            "nodes_agree": totals["agree"],
            "nodes_both_leaf": totals["both_leaf"],
            "nodes_same_feature": totals["same_feature"],
            "nodes_other_member_of_a_tie": totals["tie_other_member"],
            "nodes_with_an_exact_tie": totals["tied_nodes"],
            "subsample_size": take,
            "features_drawn_per_node": max_features,
            "largest_score_gap": score_gap,
            "library_forest_mean_abs_score_gap": gap(grm["scores"], library[0]),
            "library_forest_calls_differing": flips(grm["scores"], library[0]),
            "library_forest_seed_to_seed_mean_abs_score_gap": gap(library[0], library[1]),
            "library_forest_seed_to_seed_calls_differing": flips(library[0], library[1]),
        },
        "differences": differences,
    }
