"""exploratory.py FIXTURES GRM_DIR OUT_DIR

Exploratory. Written after the first comparison was read, so nothing here is
part of the plan and nothing here decides an outcome. It reads the same
inputs as compare.py, and compare.py's differences.tsv in OUT_DIR, and writes
OUT_DIR/exploratory.json. It answers four questions the first reading raised.

1. Is the float64 gain this kit recomputes for a tree split the gain GRM
   stored? It counts the split nodes of every decision tree and forest tree
   where `exact.float_gain` divided by the tree's sample count is
   exactly GRM's stored gain.
2. Where rounding broke an exact tie in a tree, are the two features present
   on the node's samples exactly where the other is absent?
3. Where a boosted tree and XGBoost's differ by an "equal split", do the two
   features split the node's fitting samples the same way round or the other
   way round, and which of the two names the lower feature?
4. End to end, how many of the trees GRM kept are identical to XGBoost's tree
   of the same round, on every panel, matched or not?

It needs the pinned libraries, as compare.py does.
"""

import os

for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_name] = "1"

import json
import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import numpy as np
import xgboost

import check_gbt
import check_trees
import exact
import panels as P
import table


def stored_gains(tree, rows, y, train):
    """(split nodes, those whose stored gain is the recomputed float64 gain)."""
    weights = exact.Weights(P.count(train & y), P.count(train & ~y))
    total = float(P.count(train))
    sets = check_trees.node_samples(tree, rows, train)
    splits = same = 0
    for node, samples in enumerate(sets):
        feature = tree["feature"][node]
        if feature is None:
            continue
        pos, neg = samples & y, samples & ~y
        gain = exact.float_gain(weights, P.count(pos), P.count(neg), P.count(pos & rows[feature]), P.count(neg & rows[feature]))
        splits += 1
        same += gain / total == tree["gain"][node]
    return splits, same


def equal_splits(panel, grm, rounds):
    """For the named rounds: how the two features at the first differing node relate."""
    n, y, fit = panel["n"], panel["y"], panel["fit"]
    x = check_gbt.dense(panel)
    fitting = np.array(P.members(fit), dtype=np.int64)
    labels = np.array([(y >> i) & 1 for i in range(n)], dtype=np.float64)
    weight = check_gbt.class_weights(y, fit, n)
    settings = check_gbt.settings(grm["params"])
    margin = np.zeros(n)
    out = {"same_way_round": 0, "other_way_round": 0, "grm_names_the_higher_feature": 0}
    for round_, (tree, draw) in enumerate(zip(grm["trees"], grm["draws"])):
        mine = check_gbt.nested(tree)
        if round_ in rounds:
            theirs = check_gbt.one_tree(x, labels, weight, margin, fitting, draw, settings)
            queue = [(mine, theirs, fitting)]
            while queue:
                a, b, samples = queue.pop(0)
                if "leaf" in a or "leaf" in b:
                    continue
                present = x[samples, a["feature"]] == 1.0
                if a["feature"] == b["feature"]:
                    queue += [(a["absent"], b["absent"], samples[~present]), (a["present"], b["present"], samples[present])]
                    continue
                other = x[samples, b["feature"]] == 1.0
                out["same_way_round" if (present == other).all() else "other_way_round"] += 1
                out["grm_names_the_higher_feature"] += a["feature"] > b["feature"]
                break
        check_gbt.add_tree(margin, mine, x)
    return out


def end_to_end_trees(panel, grm):
    """GRM's kept trees against XGBoost's own run, round by round."""
    if "error" in grm:
        return None
    n, y = panel["n"], panel["y"]
    x = check_gbt.dense(panel)
    fit, validation = P.decode_bits(n, grm["fit"]), P.decode_bits(n, grm["validation"])
    fitting, held = np.array(P.members(fit), dtype=np.int64), np.array(P.members(validation), dtype=np.int64)
    labels = np.array([(y >> i) & 1 for i in range(n)], dtype=np.float64)
    weight = check_gbt.class_weights(y, fit, n)
    mine = grm["colsample_1"]
    settings = check_gbt.settings(mine["params"])
    settings["eval_metric"] = "logloss"
    train = xgboost.DMatrix(x[fitting], label=labels[fitting], weight=weight[fitting], missing=np.nan, nthread=1)
    valid = xgboost.DMatrix(x[held], label=labels[held], weight=weight[held], missing=np.nan, nthread=1)
    model = xgboost.train(
        settings,
        train,
        num_boost_round=mine["params"]["rounds"],
        evals=[(valid, "validation")],
        early_stopping_rounds=mine["params"]["patience"],
        verbose_eval=False,
    )
    columns = list(range(x.shape[1]))
    theirs = [check_gbt.from_dump(json.loads(dump), columns) for dump in model.get_dump(dump_format="json")]
    ours = [check_gbt.nested(tree) for tree in mine["trees"]]
    both = min(len(ours), model.best_iteration + 1)
    return {"trees_both_keep": both, "identical": sum(check_gbt.identical(a, b) for a, b in zip(ours[:both], theirs[:both]))}


def main(argv):
    if len(argv) != 4:
        sys.exit(__doc__.splitlines()[0])
    fixtures, grm_dir, out = argv[1], pathlib.Path(argv[2]), pathlib.Path(argv[3])
    differences = table.read_differences(out / "differences.tsv")
    result = {
        "label": "exploratory: written after the first comparison was read; no outcome depends on it",
        "stored_gain": {"decision_tree": [0, 0], "forest": [0, 0]},
        "rounding_ties": {},
        "equal_splits": {},
        "end_to_end_trees": {},
    }
    for learner in ("cart", "rf"):
        rows = [d for d in differences if d["learner"] == learner and d["cause"].startswith("float64 rounding")]
        opposite = sum("exactly where the other is absent" in d["evidence"] for d in rows)
        result["rounding_ties"][learner] = {"nodes": len(rows), "features_present_where_the_other_is_absent": opposite, "other_exact_ties": len(rows) - opposite}
    for panel in P.read(fixtures):
        name = panel["name"]
        grm = json.loads((grm_dir / f"{name}.json").read_text(encoding="utf-8"))
        n, y, train = panel["n"], panel["y"], panel["fit"]
        splits, same = stored_gains(grm["decision_tree"]["tree"], panel["rows"], y, train)
        result["stored_gain"]["decision_tree"][0] += splits
        result["stored_gain"]["decision_tree"][1] += same
        screened = [panel["rows"][f] for f in grm["forest"]["screen"]]
        for entry in grm["forest"]["trees"]:
            splits, same = stored_gains(entry["tree"], screened, y, P.decode_bits(n, entry["subsample"]))
            result["stored_gain"]["forest"][0] += splits
            result["stored_gain"]["forest"][1] += same
        rounds = {int(d["unit"].split()[1]) for d in differences if d["panel"] == name and d["check"] == "tree by tree" and d["cause"] == "equal split"}
        if rounds:
            result["equal_splits"][name] = equal_splits(panel, grm["gbt_tree_by_tree"], rounds)
        trees = end_to_end_trees(panel, grm["gbt_end_to_end"])
        if trees is not None:
            result["end_to_end_trees"][name] = trees
    for learner in ("decision_tree", "forest"):
        splits, same = result["stored_gain"][learner]
        result["stored_gain"][learner] = {"split_nodes": splits, "stored_gain_is_the_recomputed_one": same}
    result["equal_splits_total"] = {key: sum(v[key] for v in result["equal_splits"].values()) for key in ("same_way_round", "other_way_round", "grm_names_the_higher_feature")}
    result["end_to_end_trees_total"] = {key: sum(v[key] for v in result["end_to_end_trees"].values()) for key in ("trees_both_keep", "identical")}
    (out / "exploratory.json").write_text(json.dumps(result, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({k: v for k, v in result.items() if k not in ("equal_splits", "end_to_end_trees")}, indent=1))


if __name__ == "__main__":
    main(sys.argv)
