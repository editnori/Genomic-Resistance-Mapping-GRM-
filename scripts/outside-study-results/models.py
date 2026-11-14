"""models.py ROOT --digests REPORTS
models.py OLDER LATER
models.py DIGESTS --reports REPORTS
models.py --fields OLDER_REPORT LATER_REPORT

ROOT is reference/outside-study. Each run of record scored the final models of
the reports it pinned. A later release can keep those models or change them,
and only in the first case is a run's estimate an estimate for a model the
site serves. This script says which, from the reports' own fields.

With --digests it prints DIGESTS. REPORTS is the `datasets/` directory of an
unpacked release. For each dataset of ROOT's `metrics.tsv` files it writes six
rows: the five learners (scm, cart, rf, gbt, lasso) and the known-gene rule.
A row holds the report's sha256, evaluation revision and input fingerprint,
and two digests.

`model_sha256` is the sha256 of one JSON text, keys sorted and no spaces, of
three things: the learner's name, everything the report holds for it under
`models`, and each feature id those fields name with the `kmer`, `n_kmers`
and `bits` of its entry in the report's `features` list (three nulls where the
report does not list it). For the forest and the boosted trees the report
publishes their largest importances and no tree, so the digest covers the
importances only. For the known-gene rule it covers `catalog.rule` and
`catalog.version`.

`read_features_sha256` is the sha256 of the features the kit reads for a
model it can score (`score.learner_features`), each with its `kmer` and
`n_kmers`, in id order. The runs wrote the same three columns to
`features.tsv`, so the suite rebuilds this digest from the committed runs.
The cell is empty for the forest, the boosted trees and the known-gene rule.

With two DIGESTS files it prints, learner by learner, on how many datasets
the two releases have the same `model_sha256`, then on how many they have the
same input fingerprint, then each row that differs. It exits 1 if a digest
or a fingerprint differs.

With --reports it reads every row of DIGESTS again from REPORTS and exits 1
unless each row is the one the reports give.

With --fields it prints every field that differs between two reports of one
dataset, in the parts the digest covers, as `learner<TAB>path`. It compares
the fields themselves and uses no digest.

A digest that is equal means every field it covers is equal, to the last
printed digit. It says nothing of a part the report does not publish. A
feature the report does not list is in the digest as its id alone. The id is
a row of the dataset matrix, and the report's `input_fingerprint` identifies
that matrix (`input_fingerprint`, src/analysis.rs), so an id stands for the
same row in two reports with one fingerprint.
"""
import csv
import hashlib
import importlib.util
import io
import json
import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent


def _module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


revision = _module("outside_study_revision", HERE / "revision.py")
tables, score = revision.tables, revision.score

LEARNERS = ("scm", "cart", "rf", "gbt", "lasso")
KNOWN = score.KNOWN
ROWS = LEARNERS + (KNOWN,)
COLUMNS = ["dataset", "learner", "report_sha256", "evaluation_revision", "input_fingerprint", "model_sha256", "read_features_sha256"]
# The columns that are the same on every row of a dataset.
OF_THE_REPORT = ("report_sha256", "input_fingerprint")


def model(report, learner):
    """Everything the report holds for one learner's final model, or for the known-gene rule."""
    m = report["models"]
    if learner == "scm":
        return {"scm": m["scm"]}
    if learner == "cart":
        return {"cart": m["cart"], "cart_depth": m["cart_depth"]}
    if learner == "lasso":
        return {"lasso": m.get("lasso"), "lasso_fit": m.get("lasso_fit")}
    if learner in ("rf", "gbt"):
        return {f"{learner}_importance": m[f"{learner}_importance"]}
    if learner == KNOWN:
        catalog = report.get("catalog")
        return None if catalog is None else {"rule": catalog["rule"], "version": catalog["version"]}
    raise ValueError(learner)


def named(report, learner):
    """Every feature id the learner's fields name, in id order."""
    m = report["models"]
    if learner == "scm":
        ids = [r["feature"] for r in m["scm"]["rules"]] + [e["feature"] for s in m["scm"]["steps"] for e in s["equivalent"]]
    elif learner == "cart":
        ids = [n["split"]["feature"] for n in m["cart"]["nodes"] if n.get("split") is not None]
    elif learner == "lasso":
        ids = [] if m.get("lasso") is None else [c[0] for c in m["lasso"]["coef"]]
    elif learner in ("rf", "gbt"):
        ids = [w["feature"] for w in m[f"{learner}_importance"]]
    else:
        ids = []
    return sorted(set(ids))


def record(report, learner):
    """What `model_sha256` covers: the learner, its fields, and what each feature id they name stands for."""
    listed = score.listed_features(report)
    features = []
    for f in named(report, learner):
        entry = listed.get(f)
        features.append([f, None, None, None] if entry is None else [f, entry["kmer"], entry["n_kmers"], entry["bits"]])
    return {"learner": learner, "model": model(report, learner), "features": features}


def read_features(report, learner):
    """[[feature, kmer, n_kmers]] of the features the kit reads for a model, in id order.

    None for a row the kit scores no model of: the forest, the boosted trees,
    the known-gene rule, and a report with no lasso."""
    if learner not in score.MODEL_LEARNERS:
        return None
    ids = score.learner_features(report, learner)
    if ids is None:
        return None
    listed = score.listed_features(report)
    return [[f, None, None] if f not in listed else [f, listed[f]["kmer"], listed[f]["n_kmers"]] for f in sorted(set(ids))]


def sha(value):
    """The sha256 of a value's JSON text: keys sorted, no spaces, ASCII, a number as Python prints it."""
    text = json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True, allow_nan=False)
    return hashlib.sha256(text.encode("ascii")).hexdigest()


def row(report, digest, learner):
    """One DIGESTS row as a dict."""
    read = read_features(report, learner)
    return {
        "dataset": report["id"],
        "learner": learner,
        "report_sha256": digest,
        "evaluation_revision": str(report["evaluation"]["revision"]),
        "input_fingerprint": report["input_fingerprint"],
        "model_sha256": sha(record(report, learner)),
        "read_features_sha256": "" if read is None else sha(read),
    }


def datasets(root):
    """The datasets of ROOT's metrics files, in name order."""
    data = tables.load(pathlib.Path(root))
    return sorted({r["dataset"] for runs in data.values() for _, _, metrics, _ in runs.values() for r in metrics})


def run_features(root):
    """{(dataset, learner): [[feature, kmer, n_kmers]]} as the runs under ROOT recorded them.

    One entry for each model a run scored, from its `metrics.tsv` row and its
    `features.tsv` rows. A lasso that keeps no weight has a row and no feature."""
    out = {}
    root = pathlib.Path(root)
    for run in sorted(p.parent for p in root.glob("*/*/metrics.tsv")):
        got = {(r["dataset"], r["learner"]): [] for r in tables.rows(run / "metrics.tsv") if r["learner"] in score.MODEL_LEARNERS}
        for r in tables.rows(run / "features.tsv"):
            got[r["dataset"], r["learner"]].append([int(r["feature"]), r["kmer"], int(r["n_kmers"])])
        for key, features in got.items():
            features.sort()
            if out.setdefault(key, features) != features:
                raise SystemExit(f"{key[0]} {key[1]}: two runs recorded different features")
    return out


def digests(root, reports):
    """DIGESTS as text: six rows for each dataset of ROOT's metrics files."""
    out = io.StringIO()
    w = csv.DictWriter(out, COLUMNS, delimiter="\t", lineterminator="\n")
    w.writeheader()
    for dataset in datasets(root):
        report, digest = revision.admitted(reports, dataset)
        for learner in ROWS:
            w.writerow(row(report, digest, learner))
    return out.getvalue()


def load(path):
    """(evaluation revision, {(dataset, learner): row}) of a DIGESTS file."""
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if reader.fieldnames != COLUMNS:
            raise SystemExit(f"{path}: the columns are not the ones this script writes")
        rows = list(reader)
    by_key = {(r["dataset"], r["learner"]): r for r in rows}
    revisions = {r["evaluation_revision"] for r in rows}
    if len(by_key) != len(rows) or len(revisions) != 1:
        raise SystemExit(f"{path}: a row is listed twice, or the reports are not of one evaluation revision")
    for dataset in {d for d, _ in by_key}:
        if any((dataset, learner) not in by_key for learner in ROWS):
            raise SystemExit(f"{path}: {dataset} lacks a row")
        if any(len({by_key[dataset, learner][c] for learner in ROWS}) != 1 for c in OF_THE_REPORT):
            raise SystemExit(f"{path}: the rows of {dataset} name more than one report or input fingerprint")
    if any(learner not in ROWS for _, learner in by_key):
        raise SystemExit(f"{path}: a row of a learner this script does not write")
    return int(revisions.pop()), by_key


def compare(older, later):
    """([row key with the same `model_sha256`], [row key where it differs]) of two loaded DIGESTS files."""
    if set(older) != set(later):
        raise SystemExit("the two files do not hold the same datasets and learners")
    same = [key for key in sorted(older) if older[key]["model_sha256"] == later[key]["model_sha256"]]
    return same, [key for key in sorted(older) if key not in set(same)]


def same_input(older, later):
    """The datasets whose reports have one input fingerprint in two loaded DIGESTS files."""
    return sorted({d for (d, learner), r in older.items() if r["input_fingerprint"] == later[d, learner]["input_fingerprint"]})


def summary(older_path, later_path):
    """(text, whether a digest or a fingerprint differs) of the comparison of two DIGESTS files."""
    (old_rev, older), (new_rev, later) = load(older_path), load(later_path)
    if old_rev == new_rev:
        raise SystemExit(f"both files are of evaluation revision {old_rev}")
    same, differ = compare(older, later)
    n = len({d for d, _ in older})
    out = [f"revision {old_rev} and revision {new_rev}, {n} datasets"]
    for learner in ROWS:
        kept = sum(1 for _, name in same if name == learner)
        out.append(f"{learner}: the same {'rule' if learner == KNOWN else 'final model'} on {kept} of {n} datasets")
    kept = sum(1 for _, name in same if name != KNOWN)
    out.append(f"the five learners: {kept} of {len(LEARNERS) * n} dataset and learner pairs")
    inputs = same_input(older, later)
    out.append(f"input fingerprint: the same on {len(inputs)} of {n} datasets")
    out += [f"differs: {dataset} {learner}" for dataset, learner in differ]
    other_input = sorted({d for d, _ in older} - set(inputs))
    out += [f"input fingerprint differs: {dataset}" for dataset in other_input]
    return "\n".join(out) + "\n", bool(differ or other_input)


def stale(path, reports):
    """The rows of a DIGESTS file that are not the rows REPORTS give, as (dataset, learner, column)."""
    _, by_key = load(path)
    out = []
    for dataset in sorted({d for d, _ in by_key}):
        report, digest = revision.admitted(reports, dataset)
        for learner in ROWS:
            want = row(report, digest, learner)
            out += [(dataset, learner, c) for c in COLUMNS if by_key[dataset, learner][c] != want[c]]
    return out


def walk(a, b, path=""):
    """The paths at which two JSON values differ. A list of another length differs as a whole."""
    if isinstance(a, dict) and isinstance(b, dict):
        out = []
        for k in sorted(set(a) | set(b)):
            out += walk(a[k], b[k], f"{path}.{k}") if k in a and k in b else [f"{path}.{k}"]
        return out
    if isinstance(a, list) and isinstance(b, list) and len(a) == len(b):
        return [p for i, (x, y) in enumerate(zip(a, b)) for p in walk(x, y, f"{path}[{i}]")]
    # 1 and 1.0 are equal numbers and other JSON texts, and the digest tells them apart.
    return [] if type(a) is type(b) and a == b else [path]


def fields(older, later):
    """[(learner, path)] of every field that differs between two reports, in the parts the digest covers."""
    out = []
    for learner in ROWS:
        out += [(learner, p.lstrip(".")) for p in walk(record(older, learner), record(later, learner))]
    return out


def main(argv):
    args = argv[1:]
    if len(args) == 3 and args[1] == "--digests":
        sys.stdout.write(digests(args[0], args[2]))
        return 0
    if len(args) == 3 and args[1] == "--reports":
        wrong = stale(args[0], args[2])
        for dataset, learner, column in wrong:
            print(f"differs from its report: {dataset} {learner} {column}")
        return 1 if wrong else 0
    if len(args) == 3 and args[0] == "--fields":
        older, later = (json.loads(pathlib.Path(p).read_bytes()) for p in args[1:])
        if older.get("id") != later.get("id"):
            raise SystemExit(f"the reports are of two datasets: {older.get('id')!r} and {later.get('id')!r}")
        for learner, path in fields(older, later):
            print(f"{learner}\t{path}")
        return 0
    if len(args) != 2 or any(a.startswith("--") for a in args):
        raise SystemExit(__doc__)
    text, failed = summary(args[0], args[1])
    sys.stdout.write(text)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
