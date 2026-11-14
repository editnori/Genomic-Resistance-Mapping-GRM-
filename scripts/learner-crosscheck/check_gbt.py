"""The boosted trees against XGBoost's exact grower.

Tree by tree: for each tree GRM grew, XGBoost grows one tree from the same
state, GRM's margin before that round as `base_margin`, on the features that
round drew. End to end: XGBoost's own loop with early stopping, against GRM's
fit with every feature offered to every tree.

A tree that is not identical gets a cause from quantities recomputed here in
float64 from GRM's margin, and again as XGBoost holds them, in float32.
"""

import json
import math

import numpy as np
import xgboost

import panels as P

LEAF_TOLERANCE = 1e-5
PROBABILITY_TOLERANCE = 1e-4
XGB_MIN_SPLIT = 1e-6  # XGBoost's kRtEps: a split needs a loss change above it
GRM_STOP_MARGIN = 1e-9  # src/learn/gbt.rs: an improvement must exceed it

# GRM's settings (GbtParams) in XGBoost's names. Each tree's features are
# handed over as that round's columns, so XGBoost itself samples nothing.
PARAMS = {
    "objective": "binary:logistic",
    "tree_method": "exact",
    "max_depth": 4,
    "eta": 0.1,
    "lambda": 1.0,
    "alpha": 0.0,
    "gamma": 0.0,
    "min_child_weight": 1.0,
    "max_delta_step": 0.0,
    "subsample": 1.0,
    "colsample_bytree": 1.0,
    "colsample_bylevel": 1.0,
    "colsample_bynode": 1.0,
    "base_score": 0.5,
    "nthread": 1,
    "seed": 0,
    "verbosity": 0,
}


def settings(params):
    """PARAMS, after checking that GRM ran with the values they map."""
    mapped = {"eta": params["learning_rate"], "max_depth": params["max_depth"], "lambda": params["lambda"], "min_child_weight": params["min_child_weight"]}
    for name, value in mapped.items():
        if PARAMS[name] != value:
            raise AssertionError(f"GRM ran with {name} = {value}, and the kit maps {PARAMS[name]}")
    return dict(PARAMS)


def dense(panel):
    x = np.zeros((panel["n"], len(panel["rows"])), dtype=np.float32)
    for feature, row in enumerate(panel["rows"]):
        x[P.members(row), feature] = 1.0
    return x


def class_weights(y, fit, n):
    """GRM's sample weights: balanced on the fitting samples, given to every sample."""
    pos = P.count(fit & y)
    neg = P.count(fit) - pos
    total = pos + neg
    w_pos, w_neg = total / (2.0 * max(pos, 1.0)), total / (2.0 * max(neg, 1.0))
    return np.array([w_pos if (y >> i) & 1 else w_neg for i in range(n)])


def nested(tree, node=0):
    """A GRM tree as nested dicts: {"leaf": v} or {"feature", "absent", "present"}."""
    if tree["feature"][node] is None:
        return {"leaf": tree["value"][node]}
    return {
        "feature": tree["feature"][node],
        "absent": nested(tree, tree["absent"][node]),
        "present": nested(tree, tree["present"][node]),
    }


def from_dump(node, columns):
    """One tree of XGBoost's JSON dump, in the same form. On 0/1 features
    its split `f < 0.5` sends the samples lacking the feature to `yes`."""
    if "leaf" in node:
        return {"leaf": node["leaf"]}
    if node["split_condition"] != 0.5:
        raise AssertionError("a split of a 0/1 feature is at 0.5")
    children = {child["nodeid"]: child for child in node["children"]}
    return {
        "feature": columns[int(node["split"][1:])],
        "absent": from_dump(children[node["yes"]], columns),
        "present": from_dump(children[node["no"]], columns),
    }


def identical(mine, theirs):
    if "leaf" in mine or "leaf" in theirs:
        return "leaf" in mine and "leaf" in theirs and abs(mine["leaf"] - theirs["leaf"]) <= LEAF_TOLERANCE
    return mine["feature"] == theirs["feature"] and identical(mine["absent"], theirs["absent"]) and identical(mine["present"], theirs["present"])


def add_tree(margin, tree, x):
    """Adds a tree's leaf values to the margin, as `gbt::fit` does."""

    def walk(node, index):
        if "leaf" in node:
            margin[index] += node["leaf"]
            return
        present = x[index, node["feature"]] == 1.0
        walk(node["absent"], index[~present])
        walk(node["present"], index[present])

    walk(tree, np.arange(x.shape[0]))


class State:
    """The gradients of one round, in float64 as GRM holds them and in
    float32 as XGBoost does. Sums over a node are in float64 in both."""

    def __init__(self, margin, labels, weight):
        p = 1.0 / (1.0 + np.exp(-margin))
        self.g = (p - labels) * weight
        self.h = np.maximum(p * (1.0 - p), 1e-6) * weight
        # GRM's gradients with XGBoost's Hessian floor, still in float64.
        self.h_floor = np.maximum(p * (1.0 - p), 1e-16) * weight
        m32, y32, w32 = margin.astype(np.float32), labels.astype(np.float32), weight.astype(np.float32)
        p32 = np.float32(1.0) / (np.float32(1.0) + np.exp(-m32))
        self.g32 = ((p32 - y32) * w32).astype(np.float64)
        self.h32 = (np.maximum(p32 * (np.float32(1.0) - p32), np.float32(1e-16)) * w32).astype(np.float64)

    def loss_change(self, x, samples, feature, lam, single=False):
        """A split's loss change, XGBoost's scale: twice GRM's gain."""
        g, h = (self.g32, self.h32) if single else (self.g, self.h)
        present = samples[x[samples, feature] == 1.0]
        gl, hl = g[present].sum(), h[present].sum()
        gt, ht = g[samples].sum(), h[samples].sum()
        gr, hr = gt - gl, ht - hl
        return gl * gl / (hl + lam) + gr * gr / (hr + lam) - gt * gt / (ht + lam), min(hl, hr)


def cause_of(mine, theirs, state, x, fitting, candidates, params):
    """Why two trees are not identical: (cause, what, evidence). The cause is "" when none is shown."""
    lam, rate, floor = params["lambda"], params["learning_rate"], params["min_child_weight"]
    queue = [(mine, theirs, fitting, "root")]
    leaves = []
    while queue:
        a, b, samples, path = queue.pop(0)
        a_leaf, b_leaf = "leaf" in a, "leaf" in b
        if a_leaf and b_leaf:
            if abs(a["leaf"] - b["leaf"]) > LEAF_TOLERANCE:
                leaves.append((a, b, samples, path))
            continue
        if not a_leaf and not b_leaf and a["feature"] == b["feature"]:
            present = x[samples, a["feature"]] == 1.0
            queue.append((a["absent"], b["absent"], samples[~present], path + ".absent"))
            queue.append((a["present"], b["present"], samples[present], path + ".present"))
            continue
        # The first node where the two trees part.
        if a_leaf:
            return "", f"at {path} GRM makes a leaf and XGBoost splits on feature {b['feature']}", ""
        mine_change, _ = state.loss_change(x, samples, a["feature"], lam)
        mine_single, _ = state.loss_change(x, samples, a["feature"], lam, single=True)
        if b_leaf:
            nu = abs(mine_change - mine_single)
            what = f"at {path} GRM splits on feature {a['feature']} and XGBoost makes a leaf"
            evidence = f"loss change {mine_change:.6g} in float64 and {mine_single:.6g} in float32; XGBoost needs more than {XGB_MIN_SPLIT:g}"
            if mine_change <= XGB_MIN_SPLIT + 4 * nu:
                return "minimum split gain", what, evidence
            return "", what, evidence
        what = f"at {path} GRM splits on feature {a['feature']} and XGBoost on feature {b['feature']}"
        mine_present = x[samples, a["feature"]] == 1.0
        their_present = x[samples, b["feature"]] == 1.0
        if (mine_present == their_present).all() or (mine_present != their_present).all():
            return "equal split", what, "the two features split the node's fitting samples alike"
        their_change, _ = state.loss_change(x, samples, b["feature"], lam)
        their_single, _ = state.loss_change(x, samples, b["feature"], lam, single=True)
        nu = max(abs(mine_change - mine_single), abs(their_change - their_single))
        valid = [state.loss_change(x, samples, f, lam) for f in candidates]
        largest = max((change for change, smaller in valid if smaller >= floor), default=-math.inf)
        evidence = (
            f"float64 loss changes {mine_change!r} and {their_change!r}, float32 {mine_single!r} and {their_single!r}; "
            f"largest float64 loss change among the round's features {largest!r}"
        )
        # The slack covers this file's own summation order, which is not GRM's.
        if mine_change >= largest - 1e-9 * max(1.0, abs(largest)) and 0 <= mine_change - their_change <= 4 * nu:
            return "float32 rounding", what, evidence
        return "", what, evidence
    # The shapes and the features agree; every leaf that differs needs a cause.
    causes, notes = [], []
    for a, b, samples, path in leaves:
        g = state.g[samples].sum()
        h_single = state.h32[samples].sum()
        refit = -g / (state.h_floor[samples].sum() + lam) * rate
        notes.append(f"{path}: GRM {a['leaf']!r}, XGBoost {b['leaf']!r}, Hessian sum {h_single:.6g}")
        if b["leaf"] == 0.0 and h_single < floor:
            causes.append("leaf set to 0")
        elif abs(refit - b["leaf"]) <= LEAF_TOLERANCE:
            causes.append("Hessian floor")
        else:
            causes.append("")
    what = f"{len(leaves)} leaf value{'s' if len(leaves) != 1 else ''} differ, the first at {leaves[0][3]}"
    if all(causes):
        return causes[0], what, "; ".join(notes[:3])
    return "", what, "; ".join(notes[:3])


def one_tree(x, labels, weight, margin, fitting, candidates, settings_):
    """The tree XGBoost grows from GRM's state on the round's features."""
    columns = sorted(candidates)
    data = xgboost.DMatrix(
        x[np.ix_(fitting, columns)],
        label=labels[fitting],
        weight=weight[fitting],
        base_margin=margin[fitting],
        missing=np.nan,
        nthread=1,
    )
    model = xgboost.train(settings_, data, num_boost_round=1)
    (dump,) = model.get_dump(dump_format="json")
    return from_dump(json.loads(dump), columns)


def tree_by_tree(panel, grm):
    n, y = panel["n"], panel["y"]
    x = dense(panel)
    fit = panel["fit"]
    fitting = np.array(P.members(fit), dtype=np.int64)
    labels = np.array([(y >> i) & 1 for i in range(n)], dtype=np.float64)
    weight = class_weights(y, fit, n)
    params = grm["params"]
    settings_ = settings(params)
    margin = np.zeros(n)
    differences, same, causes = [], 0, {}
    for round_, (tree, draw) in enumerate(zip(grm["trees"], grm["draws"])):
        mine = nested(tree)
        theirs = one_tree(x, labels, weight, margin, fitting, draw, settings_)
        if identical(mine, theirs):
            same += 1
        else:
            cause, what, evidence = cause_of(mine, theirs, State(margin, labels, weight), x, fitting, sorted(draw), params)
            causes[cause or "none"] = causes.get(cause or "none", 0) + 1
            differences.append({"unit": f"round {round_}", "what": what, "cause": cause, "evidence": evidence})
        add_tree(margin, mine, x)
    margin_gap = max((abs(a - b) for a, b in zip(margin.tolist(), grm["margin"])), default=0.0)
    if margin_gap > 1e-9:
        raise AssertionError(f"the kit's margin differs from GRM's by up to {margin_gap:.3g}: it does not read GRM's trees as GRM does")
    unexplained = sum(not d["cause"] for d in differences)
    return {
        "outcome": "matched" if not differences else "unexplained" if unexplained else "explained",
        "units": "trees",
        "compared": len(grm["trees"]),
        "matched": same,
        "explained": sum(bool(d["cause"]) for d in differences),
        "unexplained": unexplained,
        "detail": {"rounds_grown": params["rounds"], "trees_kept": grm["kept"], "features_per_tree": len(grm["draws"][0]) if grm["draws"] else None, "causes": causes, "largest_margin_gap": margin_gap},
        "differences": differences,
    }


def replay_stop(losses, patience):
    """GRM's stopping rule on a history of validation losses: the trees it would keep,
    or None if the history ends before the rule does."""
    best, kept, since = math.log(2.0), 0, 0
    for round_, loss in enumerate(losses):
        if loss < best - GRM_STOP_MARGIN:
            best, kept, since = loss, round_ + 1, 0
        else:
            since += 1
            if since >= patience:
                return kept
    return None


def end_to_end(panel, grm):
    if "error" in grm:
        return {
            "outcome": "not applicable",
            "units": "panels",
            "compared": 0,
            "matched": 0,
            "explained": 0,
            "unexplained": 0,
            "detail": {"reason": "production's two-class gate on the tuning holdout: " + grm["error"]},
            "differences": [],
        }
    n, y = panel["n"], panel["y"]
    x = dense(panel)
    fit, validation = P.decode_bits(n, grm["fit"]), P.decode_bits(n, grm["validation"])
    fitting, held = np.array(P.members(fit), dtype=np.int64), np.array(P.members(validation), dtype=np.int64)
    labels = np.array([(y >> i) & 1 for i in range(n)], dtype=np.float64)
    weight = class_weights(y, fit, n)
    mine = grm["colsample_1"]
    params = mine["params"]
    settings_ = settings(params)
    settings_["eval_metric"] = "logloss"
    train = xgboost.DMatrix(x[fitting], label=labels[fitting], weight=weight[fitting], missing=np.nan, nthread=1)
    valid = xgboost.DMatrix(x[held], label=labels[held], weight=weight[held], missing=np.nan, nthread=1)
    history = {}
    model = xgboost.train(
        settings_,
        train,
        num_boost_round=params["rounds"],
        evals=[(valid, "validation")],
        early_stopping_rounds=params["patience"],
        evals_result=history,
        verbose_eval=False,
    )
    losses = history["validation"]["logloss"]
    kept = model.best_iteration + 1
    probability = model.predict(xgboost.DMatrix(x, missing=np.nan, nthread=1), iteration_range=(0, kept))
    gap = float(np.max(np.abs(probability.astype(np.float64) - np.array(mine["probability"]))))
    columns = list(range(x.shape[1]))
    theirs = [from_dump(json.loads(dump), columns) for dump in model.get_dump(dump_format="json")]
    ours = [nested(tree) for tree in mine["trees"]]
    detail = {
        "grm_trees": mine["kept"],
        "xgboost_trees": kept,
        "xgboost_rounds_run": len(losses),
        "largest_probability_gap": gap,
        "fitting_samples": len(fitting),
        "validation_samples": len(held),
        "grm_trees_at_production_colsample": grm["colsample_production"]["kept"],
    }
    out = {"units": "panels", "compared": 1, "matched": 0, "explained": 0, "unexplained": 0, "detail": detail, "differences": []}
    if mine["kept"] == kept and gap <= PROBABILITY_TOLERANCE:
        out.update(outcome="matched", matched=1)
        return out
    what = f"GRM keeps {mine['kept']} trees and XGBoost {kept}; probabilities differ by up to {gap:.3g}"
    cause, evidence = "", ""
    first = next((r for r, (a, b) in enumerate(zip(ours, theirs)) if not identical(a, b)), None)
    detail["first_tree_not_identical"] = first
    if mine["kept"] == 0 and min(losses[: params["patience"]]) >= math.log(2.0) - 1e-6:
        cause = "no tree kept"
        evidence = f"XGBoost's lowest validation loss in its first {params['patience']} rounds is {min(losses[: params['patience']])!r}; with no tree it is ln 2 = {math.log(2.0)!r}"
    elif first is not None:
        margin = np.zeros(n)
        for tree in ours[:first]:
            add_tree(margin, tree, x)
        tree_cause, tree_what, tree_evidence = cause_of(ours[first], theirs[first], State(margin, labels, weight), x, fitting, columns, params)
        what += f"; the first tree that is not identical is round {first}: {tree_what}"
        cause = f"a tree differs ({tree_cause})" if tree_cause else ""
        evidence = tree_evidence
    elif mine["kept"] != kept:
        replayed = replay_stop(losses, params["patience"])
        evidence = f"every tree both keep is identical; GRM's stopping rule on XGBoost's validation losses keeps {replayed} trees"
        if replayed == mine["kept"]:
            cause = "stopping rule"
    else:
        evidence = "every tree is identical and both keep the same number"
    out["differences"].append({"unit": "fit", "what": what, "cause": cause, "evidence": evidence})
    out.update(outcome="explained" if cause else "unexplained", explained=int(bool(cause)), unexplained=int(not cause))
    return out
