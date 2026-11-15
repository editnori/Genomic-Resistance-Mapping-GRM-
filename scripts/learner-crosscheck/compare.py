"""compare.py FIXTURES GRM_DIR OUT_DIR [--jobs N]

The reference side of the learner cross-check. FIXTURES is the panel file
(tests/data/lasso-stopping-fixtures.txt), GRM_DIR what the Rust harness wrote
(`cargo run --release --example learner-crosscheck FIXTURES GRM_DIR`).

For each panel it sets GRM's fits against the references, by the rules of
reference/learner-crosscheck/PLAN.md:

- the set covering machine against scm_reference.py and pyscm-ml (check_scm.py);
- the decision tree and the random forest against scikit-learn (check_trees.py);
- the boosted trees against XGBoost (check_gbt.py).

It writes OUT_DIR/summary.tsv, differences.tsv, totals.json (table.py),
detail.json with each check's further counts, versions.json, and
OUT_DIR/panels/<panel>.json with everything a check recorded. Since 15
November 2026 it also writes end-to-end-trees.tsv (table.py), a record added
after the plan.

It needs the pinned libraries of requirements.txt and requirements-scm.txt.
It runs N panels at a time (default 4, at most 8), one thread each. It reads
its two inputs and writes OUT_DIR only.
"""

import os

# One thread per process; set before numpy and XGBoost load.
for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_name] = "1"

import concurrent.futures
import hashlib
import json
import pathlib
import platform
import sys

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))  # `python3 -I` leaves the script's folder out

import panels as P
import table

KEEP = ("outcome", "units", "compared", "matched", "explained", "unexplained")


def sha256(path):
    return hashlib.sha256(pathlib.Path(path).read_bytes()).hexdigest()


def compare_panel(fixtures, grm_dir, name):
    import check_gbt
    import check_scm
    import check_trees

    (panel,) = [p for p in P.read(fixtures) if p["name"] == name]
    grm = json.loads((pathlib.Path(grm_dir) / f"{name}.json").read_text(encoding="utf-8"))
    if grm["panel"] != name or grm["counts"] != P.counts(panel):
        raise SystemExit(f"{name}: the Rust harness and panels.py decode the panel differently")
    return {
        "panel": name,
        "counts": grm["counts"],
        "selected": grm["selected"],
        "checks": {
            "scm fits": check_scm.check(panel, grm["scm"]),
            "cart nodes": check_trees.check_decision_tree(panel, grm["decision_tree"]),
            "rf trees": check_trees.check_forest(panel, grm["forest"]),
            "gbt tree by tree": check_gbt.tree_by_tree(panel, grm["gbt_tree_by_tree"]),
            "gbt end to end": check_gbt.end_to_end(panel, grm["gbt_end_to_end"]),
        },
    }


def versions(fixtures, grm_dir, names):
    import joblib
    import numpy
    import scipy
    import sklearn
    import xgboost
    from importlib import metadata

    index = json.loads((pathlib.Path(grm_dir) / "index.json").read_text(encoding="utf-8"))
    return {
        "python": platform.python_version(),
        "platform": platform.platform(),
        "libraries": {
            "numpy": numpy.__version__,
            "scipy": scipy.__version__,
            "scikit-learn": sklearn.__version__,
            "joblib": joblib.__version__,
            "threadpoolctl": metadata.version("threadpoolctl"),
            "six": metadata.version("six"),
            "xgboost": xgboost.__version__,
            "pyscm-ml": metadata.version("pyscm-ml"),
        },
        "xgboost_build": {key: xgboost.build_info()[key] for key in ("USE_CUDA", "USE_OPENMP", "GCC_VERSION", "GLIBC_VERSION") if key in xgboost.build_info()},
        "evaluation_revision": index["evaluation_revision"],
        "fixtures_sha256": sha256(fixtures),
        "grm_fits_sha256": {name: sha256(pathlib.Path(grm_dir) / f"{name}.json") for name in names},
    }


def main(argv):
    args = argv[1:]
    jobs = 4
    if "--jobs" in args:
        at = args.index("--jobs")
        jobs = int(args[at + 1])
        del args[at : at + 2]
    if len(args) != 3 or not 1 <= jobs <= 8:
        sys.exit(__doc__.splitlines()[0])
    fixtures, grm_dir, out = args[0], args[1], pathlib.Path(args[2])
    (out / "panels").mkdir(parents=True, exist_ok=True)
    panels = P.read(fixtures)
    names = [panel["name"] for panel in panels]

    with concurrent.futures.ProcessPoolExecutor(max_workers=jobs) as pool:
        results = list(pool.map(compare_panel, [fixtures] * len(names), [grm_dir] * len(names), names))

    summary, differences, trees, detail = [], [], [], {}
    for panel, result in zip(panels, results):
        name = panel["name"]
        (out / "panels" / f"{name}.json").write_text(json.dumps(result, indent=1, sort_keys=True) + "\n", encoding="utf-8")
        detail[name] = {"selected": result["selected"]}
        for learner, check, _ in table.CHECKS:
            entry = result["checks"][f"{learner} {check}"]
            summary.append(
                {
                    "panel": name,
                    "samples": panel["n"],
                    "features": len(panel["rows"]),
                    "fit_samples": P.count(panel["fit"]),
                    "learner": learner,
                    "check": check,
                    **{key: entry[key] for key in KEEP},
                }
            )
            differences += [{"panel": name, "learner": learner, "check": check, **row} for row in entry["differences"]]
            detail[name][f"{learner} {check}"] = entry["detail"]
        trees += [{"panel": name, **row} for row in result["checks"]["gbt end to end"]["trees"]]
    table.write(out / "summary.tsv", table.SUMMARY_COLUMNS, summary)
    table.write(out / "differences.tsv", table.DIFFERENCE_COLUMNS, differences)
    table.write(out / "end-to-end-trees.tsv", table.TREE_COLUMNS, trees)
    summary = table.read_summary(out / "summary.tsv")
    found = table.problems(summary, table.read_differences(out / "differences.tsv"), names)
    if found:
        sys.exit("the tables are not consistent:\n" + "\n".join(found))
    dump = lambda value: json.dumps(value, indent=1, sort_keys=True) + "\n"
    (out / "totals.json").write_text(dump(table.totals(summary)), encoding="utf-8")
    (out / "detail.json").write_text(dump(detail), encoding="utf-8")
    (out / "versions.json").write_text(dump(versions(fixtures, grm_dir, names)), encoding="utf-8")
    for row in table.totals(summary):
        print(row["learner"], row["check"], row["panels"], f"{row['units_matched']} of {row['compared']} {row['units']}", sep="\t")


if __name__ == "__main__":
    main(sys.argv)
