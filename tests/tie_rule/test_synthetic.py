"""scripts/tie-rule/record.py on SYNTHETIC reports, traces, release and plan written here.

No real report, trace or release is read. Every dataset id starts with
SYNTHETIC, every species slug with synthetic_, and the plan says SYNTHETIC.
"""

import contextlib
import hashlib
import importlib.util
import io
import json
import math
import pathlib
import subprocess
import sys
import tempfile
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "tie-rule" / "record.py"


def load():
    spec = importlib.util.spec_from_file_location("tie_rule_record_synthetic", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


R = load()
LEARNERS = ("scm", "cart", "rf", "gbt", "lasso")
SPLITS = ("lineage", "random")
MARK = "SYNTHETIC-NOT-COPIED"

LABELS = {
    "SYNTHETIC-a1": [1, 1, 1, 0, 0, 0],
    "SYNTHETIC-a2": [1, 0, 1, 0, 1, 0],
    "SYNTHETIC-b1": [1, 1, 0, 0, 0, 0],
    "SYNTHETIC-b2": [1, 0, 0, 1],
    "SYNTHETIC-c1": [1, 0, 0],
}
SPECIES = {
    "SYNTHETIC-a1": "synthetic_alpha",
    "SYNTHETIC-a2": "synthetic_alpha",
    "SYNTHETIC-b1": "synthetic_beta",
    "SYNTHETIC-b2": "synthetic_beta",
    "SYNTHETIC-c1": "synthetic_gamma",
}
FEATURES = {"SYNTHETIC-a1": 40, "SYNTHETIC-a2": 60, "SYNTHETIC-b1": 80, "SYNTHETIC-b2": 100, "SYNTHETIC-c1": 500}
ONLY_RANDOM = {"SYNTHETIC-a2"}
RUN = ("SYNTHETIC-a1", "SYNTHETIC-a2", "SYNTHETIC-b1")
# SYNTHETIC-b2 has every report and C's trace but not E's; SYNTHETIC-c1 has no report.
HAND = {
    ("SYNTHETIC-a1", "lineage", "cart"): ([0.9, 0.4, 0.7, 0.2, 0.6, 0.1], [0.9, 0.6, 0.7, 0.2, 0.6, 0.1]),
}
SECONDS = {"A": 1.0, "B": 1.25, "C": 2.0, "D": 1.5, "E": 3.0, "release": 9.0}


def present_scores(dataset, split, learner):
    d = sorted(LABELS).index(dataset)
    k = LEARNERS.index(learner)
    s = SPLITS.index(split)
    out = [round(((i * 3 + k * 5 + s * 7 + d * 11) % 19) / 20 + 0.025, 3) for i in range(len(LABELS[dataset]))]
    if dataset == "SYNTHETIC-a2" and learner == "scm":
        out[0] = None
    return out


def scores(dataset, split, learner, rule):
    if (dataset, split, learner) in HAND:
        return list(HAND[(dataset, split, learner)][rule == "exact"])
    out = present_scores(dataset, split, learner)
    if rule == "exact":
        if learner == "rf":
            out[0] = round(out[0] + 0.01, 3)
        if learner == "gbt" and dataset == "SYNTHETIC-b1":
            out[0] = round(1 - out[0], 3)
    return out


def metrics(sc, labels):
    tp = sum(1 for s, y in zip(sc, labels) if s is not None and s >= 0.5 and y)
    fn = sum(1 for s, y in zip(sc, labels) if s is not None and s < 0.5 and y)
    tn = sum(1 for s, y in zip(sc, labels) if s is not None and s < 0.5 and not y)
    fp = sum(1 for s, y in zip(sc, labels) if s is not None and s >= 0.5 and not y)
    sens, spec = tp / (tp + fn), tn / (tn + fp)
    return dict(n=tp + fn + tn + fp, tp=tp, fn=fn, tn=tn, fp=fp, sensitivity=sens, specificity=spec,
                balanced_accuracy=(sens + spec) / 2)


def report(dataset, arm):
    labels = LABELS[dataset]
    rule = "exact" if arm in "DE" else "present"
    splits = []
    for split in SPLITS:
        if split == "lineage" and dataset in ONLY_RANDOM:
            splits.append({"scheme": split, "unavailable_reason": "SYNTHETIC: one lineage group", "learners": []})
            continue
        rows = []
        for k, learner in enumerate(LEARNERS):
            sc = scores(dataset, split, learner, rule)
            m = metrics(sc, labels)
            rows.append({"learner": learner, "metrics": m, "scores": sc, "seconds": SECONDS[arm] * (k + 1),
                         "fold_balanced_accuracy": [m["balanced_accuracy"]], "note": MARK})
        splits.append({"scheme": split, "unavailable_reason": None, "learners": rows})
    return {"label": "SYNTHETIC", "id": dataset, "species_slug": SPECIES[dataset], "note": MARK,
            "counts": {"samples": len(labels), "features": FEATURES[dataset]},
            "samples": [{"id": f"{dataset}-g{i}", "resistant": bool(y)} for i, y in enumerate(labels)],
            "splits": splits}


def fit_lines(rule):
    lines = []
    for learner, objective in (("cart", None), ("rf", None), ("gbt", "logistic"), ("gbt", "censored")):
        for seed in (1, 2):
            if rule == "present":
                tally = dict(searched=10, exact_tied=4, differ_tie=0, differ_order=0, differ_partition=0)
                tally["differ_partition" if learner == "gbt" else "differ_tie"] = int(seed == 1)
            else:
                tally = dict(searched=11, exact_tied=5, differ_tie=0, differ_order=0, differ_partition=0)
            trees = {
                "cart": ["c-p"] if rule == "present" or seed == 2 else ["c-e"],
                "rf": ["r1", "r2"] if rule == "present" else ["r1", "r3", "r4"],
                "gbt": ["g1", "g2"],
            }[learner]
            fit = {"learner": learner, "seed": seed, "fit_samples": 5, "tally": tally, "trees": trees, "note": MARK}
            if objective:
                fit["objective"] = objective
            lines.append(json.dumps(fit))
    return "\n".join(lines) + "\n"


def write_json(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=1) + "\n", encoding="utf-8")


def world(base):
    """RUNS, RELEASE and PLAN folders under base, all SYNTHETIC."""
    base = pathlib.Path(base)
    runs, release, plan = base / "runs", base / "release", base / "plan"
    for d in ("SYNTHETIC-a1", "SYNTHETIC-a2", "SYNTHETIC-b1", "SYNTHETIC-b2"):
        for arm in "ABCDE":
            write_json(runs / arm / f"{d}.json", report(d, arm))
        (runs / "C" / f"{d}.trace.jsonl").write_text(fit_lines("present"), encoding="utf-8")
        if d != "SYNTHETIC-b2":
            (runs / "E" / f"{d}.trace.jsonl").write_text(fit_lines("exact"), encoding="utf-8")
        rel = report(d, "A")
        for s in rel["splits"]:
            for row in s["learners"]:
                row["seconds"] = SECONDS["release"]
        write_json(release / "datasets" / f"{d}.json", rel)
    write_json(release / "index.json", {"label": "SYNTHETIC", "datasets": [
        {"id": d, "species_slug": SPECIES[d], "n": len(LABELS[d]), "n_resistant": sum(LABELS[d]), "features": FEATURES[d],
         "note": MARK} for d in sorted(LABELS, reverse=True)]})
    plan.mkdir(parents=True)
    (plan / "PLAN.md").write_text("# SYNTHETIC plan\n\nWritten by tests/tie_rule.\n", encoding="utf-8")
    digest = hashlib.sha256((plan / "PLAN.md").read_bytes()).hexdigest()
    (plan / "PLAN.sha256").write_text(f"{digest}  PLAN.md\n", encoding="utf-8")
    return runs, release, plan


def run_main(*argv):
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        rc = R.main([str(a) for a in argv])
    return rc, out.getvalue()


def files_of(directory):
    return {p.name: p.read_bytes() for p in sorted(pathlib.Path(directory).iterdir())}


class Synthetic(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory(prefix="SYNTHETIC-tie-rule-")
        self.addCleanup(tmp.cleanup)
        self.base = pathlib.Path(tmp.name)
        self.runs, self.release, self.plan = world(self.base)
        self.out = self.base / "out"
        rc, printed = run_main(self.runs, self.release, self.plan, "--out", self.out)
        self.assertEqual((rc, printed), (0, ""))

    def rows(self, name):
        return R.read_table(self.out / name)

    def summary(self):
        return {r["key"]: r["value"] for r in self.rows("summary.tsv")}

    def measure(self, dataset, split, learner):
        found = [r for r in self.rows("measures.tsv") if (r["dataset"], r["split"], r["learner"]) == (dataset, split, learner)]
        self.assertEqual(len(found), 1)
        return found[0]

    # ---------------------------------------------------------------- the same bytes

    def test_a_second_run_writes_the_same_bytes(self):
        again = self.base / "again"
        self.assertEqual(run_main(self.runs, self.release, self.plan, "--out", again), (0, ""))
        self.assertEqual(files_of(self.out), files_of(again))

    def test_two_hash_seeds_write_the_same_bytes(self):
        written = []
        for seed in ("0", "4242"):
            target = self.base / f"seed-{seed}"
            env = {"PATH": "/usr/bin:/bin", "PYTHONHASHSEED": seed}
            done = subprocess.run([sys.executable, "-s", "-P", "-B", str(SCRIPT), str(self.runs), str(self.release),
                                   str(self.plan), "--out", str(target)], env=env, capture_output=True, text=True, check=False)
            self.assertEqual((done.returncode, done.stdout, done.stderr), (0, "", ""))
            written.append(files_of(target))
        self.assertEqual(written[0], written[1])
        self.assertEqual(written[0], files_of(self.out))

    def test_the_files_are_the_record_s_files_and_plain_text(self):
        want = set(R.COLUMNS) | {"PLAN.md", "PLAN.sha256", "SHA256SUMS"}
        self.assertEqual(set(files_of(self.out)), want)
        for name, data in files_of(self.out).items():
            text = data.decode("utf-8")
            self.assertTrue(text.endswith("\n"), name)
            self.assertNotIn("\r", text, name)
            if name.endswith(".tsv"):
                lines = text.splitlines()
                self.assertEqual(tuple(lines[0].split("\t")), R.COLUMNS[name], name)
                self.assertTrue(all(len(line.split("\t")) == len(R.COLUMNS[name]) for line in lines), name)

    def test_no_string_of_a_report_but_ids_species_and_learners_is_written(self):
        for name, data in files_of(self.out).items():
            self.assertNotIn(MARK, data.decode("utf-8"), name)
            self.assertNotIn("SYNTHETIC: one lineage group", data.decode("utf-8"), name)
            self.assertNotIn("-g0", data.decode("utf-8"), name)

    # ---------------------------------------------------------------- row files

    def test_datasets_rows_follow_the_index_and_the_arms(self):
        rows = self.rows("datasets.tsv")
        self.assertEqual([r["dataset"] for r in rows], sorted(LABELS))
        for r in rows:
            d = r["dataset"]
            self.assertEqual((r["species"], r["genomes"], r["resistant"], r["features"]),
                             (SPECIES[d], len(LABELS[d]), sum(LABELS[d]), FEATURES[d]))
            self.assertEqual(r["run"], int(d in RUN), d)

    def test_a_dataset_missing_from_an_arm_is_not_run(self):
        for name in ("measures.tsv", "seconds.tsv", "traces.tsv", "reproduction.tsv"):
            self.assertEqual(sorted({r["dataset"] for r in self.rows(name)}), list(RUN), name)
        text = (self.out / "inputs.sha256").read_text(encoding="utf-8")
        self.assertNotIn("SYNTHETIC-b2", text)
        self.assertNotIn("SYNTHETIC-c1", text)

    def test_the_cells_of_one_learner_by_hand(self):
        r = self.measure("SYNTHETIC-a1", "lineage", "cart")
        # present calls: 1 0 1 | 0 1 0 against labels 1 1 1 | 0 0 0; exact calls: 1 1 1 | 0 1 0
        self.assertEqual((r["recall_present"], r["recall_exact"]), (2 / 3, 1.0))
        self.assertEqual((r["specificity_present"], r["specificity_exact"]), (2 / 3, 2 / 3))
        self.assertEqual((r["balanced_accuracy_present"], r["balanced_accuracy_exact"]), (2 / 3, (1.0 + 2 / 3) / 2))
        self.assertEqual(r["recall_change"], 1.0 - 2 / 3)
        self.assertEqual(r["specificity_change"], 0.0)
        self.assertEqual((r["correct_present"], r["correct_exact"], r["calls_changed"]), (4, 5, 1))
        self.assertEqual((r["scores_changed"], r["metrics_changed"], r["fold_balanced_accuracy_changed"]), (1, 1, 1))
        self.assertEqual((r["species"], r["genomes"]), ("synthetic_alpha", 6))
        present = -(math.log(0.9) + math.log(0.4) + math.log(0.7) + math.log(0.8) + math.log(0.4) + math.log(0.9)) / 6
        exact = -(math.log(0.9) + math.log(0.6) + math.log(0.7) + math.log(0.8) + math.log(0.4) + math.log(0.9)) / 6
        self.assertAlmostEqual(r["log_loss_present"], present, delta=1e-15)
        self.assertAlmostEqual(r["log_loss_exact"], exact, delta=1e-15)
        self.assertAlmostEqual(r["log_loss_change"], exact - present, delta=1e-15)

    def test_an_unchanged_learner_has_no_change(self):
        for d, split in (("SYNTHETIC-a1", "random"), ("SYNTHETIC-b1", "lineage")):
            for learner in ("scm", "lasso", "cart"):
                r = self.measure(d, split, learner)
                self.assertEqual((r["scores_changed"], r["metrics_changed"], r["calls_changed"]), (0, 0, 0))
                self.assertEqual((r["recall_change"], r["log_loss_change"]), (0.0, 0.0))

    def test_a_score_that_moves_without_crossing_changes_log_loss_only(self):
        r = self.measure("SYNTHETIC-b1", "random", "rf")
        self.assertEqual((r["scores_changed"], r["metrics_changed"], r["calls_changed"]), (1, 0, 0))
        self.assertNotEqual(r["log_loss_change"], 0.0)

    def test_a_missing_score_is_left_out(self):
        r = self.measure("SYNTHETIC-a2", "random", "scm")
        self.assertEqual(r["genomes"], 5)
        sc = present_scores("SYNTHETIC-a2", "random", "scm")
        labels = LABELS["SYNTHETIC-a2"]
        want = -sum(math.log(s) if y else math.log(1 - s) for s, y in zip(sc[1:], labels[1:])) / 5
        self.assertAlmostEqual(r["log_loss_present"], want, delta=1e-15)

    def test_log_loss_clamps_scores_of_0_and_1(self):
        self.assertEqual(R.log_loss([0.0], [True]), -math.log(1e-12))
        self.assertEqual(R.log_loss([1.0], [False]), -math.log(1 - (1 - 1e-12)))
        self.assertEqual(R.log_loss([1.0, None], [True, False]), -math.log(1 - 1e-12))
        self.assertIsNone(R.log_loss([None], [True]))
        self.assertEqual(R.calls([0.5, 0.4999, None]), [True, False, None])

    def test_an_unavailable_split_is_left_out(self):
        self.assertEqual({r["split"] for r in self.rows("measures.tsv") if r["dataset"] == "SYNTHETIC-a2"}, {"random"})
        s = self.summary()
        self.assertEqual((s["datasets.lineage"], s["datasets.random"]), (2, 3))

    def test_seconds_are_summed_over_the_splits(self):
        rows = {(r["dataset"], r["learner"]): r for r in self.rows("seconds.tsv")}
        for d in RUN:
            nsplits = 1 if d in ONLY_RANDOM else 2
            for k, learner in enumerate(LEARNERS):
                r = rows[(d, learner)]
                self.assertEqual((r["seconds_main"], r["seconds_present"], r["seconds_exact"]),
                                 tuple(float(SECONDS[a] * (k + 1) * nsplits) for a in "ABD"))

    def test_reproduction_masks_the_seconds_only(self):
        for r in self.rows("reproduction.tsv"):
            self.assertEqual(tuple(r[c] for c in R.COLUMNS["reproduction.tsv"][1:]), (1, 1, 1, 1), r["dataset"])
        path = self.runs / "C" / "SYNTHETIC-a1.json"
        r = json.loads(path.read_text(encoding="utf-8"))
        r["splits"][0]["learners"][1]["scores"][0] = 0.91
        write_json(path, r)
        rep = {x["dataset"]: x for x in R.collect(self.runs, self.release)["reproduction"]}
        self.assertFalse(rep["SYNTHETIC-a1"]["traced_present_equals_present"])
        self.assertTrue(rep["SYNTHETIC-a1"]["present_equals_main"])
        self.assertTrue(rep["SYNTHETIC-a2"]["traced_present_equals_present"])

    def test_the_release_report_is_compared_with_arm_a(self):
        path = self.release / "datasets" / "SYNTHETIC-b1.json"
        r = json.loads(path.read_text(encoding="utf-8"))
        r["counts"]["features"] += 1
        write_json(path, r)
        rep = {x["dataset"]: x for x in R.collect(self.runs, self.release)["reproduction"]}
        self.assertEqual((rep["SYNTHETIC-b1"]["main_equals_release"], rep["SYNTHETIC-a1"]["main_equals_release"]), (False, True))

    def test_inputs_give_the_sha256_of_each_report_and_trace(self):
        lines = (self.out / "inputs.sha256").read_text(encoding="utf-8").splitlines()
        self.assertEqual(len(lines), 5 * len(RUN))
        names = []
        for line in lines:
            digest, name = line.split("  ")
            suffix = ".trace.jsonl" if name[0] in "CE" else ".json"
            self.assertTrue(name.endswith(suffix), name)
            self.assertEqual(hashlib.sha256((self.runs / name).read_bytes()).hexdigest(), digest, name)
            names.append(name)
        self.assertEqual(names, sorted(names))
        self.assertEqual({n.split("/")[0] for n in names}, set("ABCDE"))

    def test_trace_rows_by_hand(self):
        rows = {(r["dataset"], r["kind"]): r for r in self.rows("traces.tsv")}
        self.assertEqual(len(rows), 4 * len(RUN))
        for d in RUN:
            cart, rf, gbt, mic = (rows[(d, k)] for k in ("cart", "rf", "gbt", "gbt-mic"))
            for r in (cart, rf, gbt, mic):
                self.assertEqual((r["fits"], r["searched_present"], r["exact_tied_present"], r["searched_exact"],
                                  r["exact_tied_exact"], r["fits_differing_present"], r["fits_differing_exact"],
                                  r["fit_mismatches"]), (2, 20, 8, 22, 10, 1, 0, 0))
            self.assertEqual((cart["differ_tie_present"], cart["differ_partition_present"]), (1, 0))
            self.assertEqual((gbt["differ_tie_present"], gbt["differ_partition_present"]), (0, 1))
            self.assertEqual((mic["differ_partition_present"], mic["differ_order_present"]), (1, 0))
            self.assertEqual((cart["trees_present"], cart["trees_exact"], cart["tree_positions"], cart["trees_differ"]), (2, 2, 2, 1))
            self.assertEqual((rf["trees_present"], rf["trees_exact"], rf["tree_positions"], rf["trees_differ"]), (4, 6, 6, 4))
            self.assertEqual((gbt["tree_positions"], gbt["trees_differ"]), (4, 0))

    def test_a_fit_that_does_not_line_up_is_counted_and_not_tallied(self):
        present, exact = self.base / "p.jsonl", self.base / "e.jsonl"
        present.write_text(fit_lines("present"), encoding="utf-8")
        lines = fit_lines("exact").splitlines()
        fit = json.loads(lines[0])
        fit["seed"] = 99
        lines[0] = json.dumps(fit)
        exact.write_text("\n".join(lines) + "\n", encoding="utf-8")
        cart = R.trace_rows("SYNTHETIC-x", present, exact)[0]
        self.assertEqual((cart["kind"], cart["fits"], cart["fit_mismatches"], cart["searched_present"], cart["trees_present"]),
                         ("cart", 2, 1, 10, 1))
        exact.write_text("\n".join(lines[1:]) + "\n", encoding="utf-8")
        with self.assertRaises(SystemExit):
            R.trace_rows("SYNTHETIC-x", present, exact)

    # ---------------------------------------------------------------- derived files

    def test_the_derived_files_come_from_the_row_files_alone(self):
        derived = R.derive(R.read_rows(self.out))
        for name in R.DERIVED_FILES:
            self.assertEqual(R.table(name, derived[name[: -len(".tsv")]]).encode("utf-8"), (self.out / name).read_bytes(), name)

    def test_the_best_learner_ties_to_the_later_learner(self):
        self.assertEqual(R.best({"scm": 0.8, "cart": 0.7, "rf": 0.8, "gbt": 0.6, "lasso": 0.5}), ("rf", 0.8))
        self.assertEqual(R.best({"scm": 0.8, "cart": 0.7, "rf": 0.8, "gbt": 0.6, "lasso": 0.8}), ("lasso", 0.8))
        self.assertEqual(R.best({"scm": 0.9, "cart": None}), ("scm", 0.9))

    def test_splits_rows_follow_the_measures(self):
        measures = self.rows("measures.tsv")
        for r in self.rows("splits.tsv"):
            m = {x["learner"]: x for x in measures if (x["dataset"], x["split"]) == (r["dataset"], r["split"])}
            for rule in ("present", "exact"):
                ba = {k: m[k][f"balanced_accuracy_{rule}"] for k in LEARNERS}
                top = max(reversed(LEARNERS), key=lambda k: ba[k])
                self.assertEqual((r[f"best_{rule}"], r[f"best_{rule}_balanced_accuracy"]), (top, ba[top]))
            ranked = sorted((m[k]["balanced_accuracy_present"] for k in LEARNERS), reverse=True)
            self.assertEqual((r["second_present_balanced_accuracy"], r["gap_present"]), (ranked[1], ranked[0] - ranked[1]))
            effect = max(abs(m[k]["balanced_accuracy_change"]) for k in ("cart", "rf", "gbt"))
            self.assertEqual(r["tie_effect"], effect)
            self.assertEqual(r["tie_effect_at_least_gap"], int(effect > 0 and effect >= r["gap_present"]))
            self.assertEqual(r["best_changes"], int(r["best_present"] != r["best_exact"]))
        hand = [r for r in self.rows("splits.tsv") if (r["dataset"], r["split"]) == ("SYNTHETIC-a1", "lineage")][0]
        self.assertEqual(hand["tie_effect"], (1.0 + 2 / 3) / 2 - 2 / 3)

    def test_order_counts_a_pair_whose_sign_changes(self):
        measures = self.rows("measures.tsv")
        rows = self.rows("order.tsv")
        self.assertEqual(len(rows), 2 * 10)
        for r in rows:
            n = 0
            keys = sorted({x["dataset"] for x in measures if x["split"] == r["split"]})
            for d in keys:
                m = {x["learner"]: x for x in measures if (x["dataset"], x["split"]) == (d, r["split"])}
                before = m[r["first"]]["balanced_accuracy_present"] - m[r["second"]]["balanced_accuracy_present"]
                after = m[r["first"]]["balanced_accuracy_exact"] - m[r["second"]]["balanced_accuracy_exact"]
                n += (before > 0, before < 0) != (after > 0, after < 0)
            self.assertEqual((r["datasets"], r["order_changes"]), (len(keys), n), r)
        self.assertGreater(sum(r["order_changes"] for r in rows), 0)

    def test_paired_rows_by_hand(self):
        measures = self.rows("measures.tsv")
        rows = self.rows("paired.tsv")
        self.assertEqual([(r["learner"], r["split"], r["measure"]) for r in rows],
                         [(t, s, m) for t in ("cart", "rf", "gbt") for s in SPLITS for m in R.PAIRED_MEASURES])
        for r in rows:
            rs = [x for x in measures if (x["learner"], x["split"]) == (r["learner"], r["split"])]
            d = [x[f"{r['measure']}_change"] for x in sorted(rs, key=lambda x: x["dataset"])]
            groups = {}
            for x in rs:
                groups.setdefault(x["species"], []).append(x[f"{r['measure']}_change"])
            means = [sum(v) / len(v) for v in groups.values()]
            self.assertEqual((r["datasets"], r["species"]), (len(d), len(groups)))
            self.assertEqual(r["mean"], sum(d) / len(d))
            self.assertEqual((r["higher"], r["lower"], r["equal"]), (sum(x > 0 for x in d), sum(x < 0 for x in d), sum(x == 0 for x in d)))
            self.assertLessEqual(min(means) - 1e-12, r["low"])
            self.assertLessEqual(r["low"], r["high"])
            self.assertLessEqual(r["high"], max(means) + 1e-12)

    def test_one_species_gives_an_interval_of_one_point(self):
        rows = R.read_rows(self.out)
        for r in rows["measures"]:
            r["species"] = "synthetic_alpha"
        for r in R.derive(rows)["paired"]:
            self.assertEqual((r["species"], r["low"], r["high"]), (1, r["mean"], r["mean"]))

    def test_summary_counts(self):
        s = self.summary()
        self.assertEqual((s["datasets.release"], s["datasets.run"], s["datasets.not_run"]), (5, 3, 2))
        self.assertEqual((s["species.release"], s["species.run"], s["species.not_run"], s["species.both"]), (3, 2, 1, 1))
        self.assertEqual(s["species.not_run.names"], "synthetic_gamma")
        self.assertEqual((s["genomes.median.run"], s["genomes.median.not_run"]), (6, 3.5))
        self.assertEqual((s["features.median.run"], s["features.median.not_run"], s["features.median.ratio"]), (60, 300.0, 5.0))
        self.assertEqual([s[f"{k}.either.metrics_changed"] for k in LEARNERS], [0, 1, 0, 1, 0])
        self.assertEqual([s[f"{k}.either.log_loss_changed"] for k in LEARNERS], [0, 1, 3, 1, 0])
        self.assertEqual((s["trees.either.metrics_changed"], s["trees.either.log_loss_changed"]), (2, 3))
        self.assertEqual(s["cart.lineage.calls_changed"], 1)
        self.assertEqual(s["cart.lineage.net_correct"], 1)
        self.assertEqual(s["cart.lineage.recall.over_0.05"], 1)
        self.assertEqual(s["cart.lineage.recall.largest_abs_change.dataset"], "SYNTHETIC-a1")
        self.assertEqual(s["rf.random.scores_changed"], 3)
        self.assertEqual(s["trees.calls_changed"], sum(s[f"{t}.{x}.calls_changed"] for t in ("cart", "rf", "gbt") for x in SPLITS))
        self.assertEqual(s["trace.cart.searched_present"], 20 * 3)
        self.assertEqual(s["trace.gbt-mic.differ_present"], 3)
        self.assertEqual(s["trace.rf.exact_tied_share"], 0.4)
        self.assertEqual((s["seconds.cart.median_ratio_exact_present"], s["seconds.lasso.median_ratio_present_main"]), (1.2, 1.25))
        self.assertEqual((s["reproduction.main_equals_release"], s["inputs.files"]), (3, 15))

    def test_settings_are_the_script_s(self):
        self.assertEqual([(r["key"], str(r["value"])) for r in self.rows("settings.tsv")], list(R.SETTINGS))

    def test_a_value_that_is_not_finite_is_refused(self):
        for v in (math.nan, math.inf, -math.inf):
            with self.assertRaises(ValueError):
                R.cell(v)
        self.assertEqual([R.cell(v) for v in (None, True, False, 0.1, 3, "x")], ["", "1", "0", "0.1", "3", "x"])

    # ---------------------------------------------------------------- what stops the script

    def test_the_plan_must_have_its_sha256(self):
        (self.plan / "PLAN.md").write_text("# SYNTHETIC plan, changed\n", encoding="utf-8")
        with self.assertRaises(SystemExit):
            run_main(self.runs, self.release, self.plan, "--out", self.base / "other")
        self.assertFalse((self.base / "other").exists())

    def test_counts_that_differ_from_the_index_stop_it(self):
        path = self.runs / "B" / "SYNTHETIC-a1.json"
        r = json.loads(path.read_text(encoding="utf-8"))
        r["counts"]["samples"] += 1
        write_json(path, r)
        with self.assertRaises(SystemExit):
            R.collect(self.runs, self.release)

    def test_arms_that_scored_other_learners_stop_it(self):
        path = self.runs / "D" / "SYNTHETIC-b1.json"
        r = json.loads(path.read_text(encoding="utf-8"))
        r["splits"][1]["learners"].pop()
        write_json(path, r)
        with self.assertRaises(SystemExit):
            R.collect(self.runs, self.release)

    def test_a_run_dataset_missing_from_the_index_stops_it(self):
        index = json.loads((self.release / "index.json").read_text(encoding="utf-8"))
        index["datasets"] = [d for d in index["datasets"] if d["id"] != "SYNTHETIC-a2"]
        write_json(self.release / "index.json", index)
        with self.assertRaises(SystemExit):
            R.collect(self.runs, self.release)

    # ---------------------------------------------------------------- --check

    def test_check_passes_the_written_record(self):
        self.assertEqual(run_main(self.runs, self.release, self.plan, "--check", self.out), (0, ""))

    def test_check_names_a_changed_byte(self):
        path = self.out / "measures.tsv"
        data = bytearray(path.read_bytes())
        at = data.index(b"0.6666666666666666")
        data[at + 2] = ord("7")
        path.write_bytes(bytes(data))
        self.assertEqual(run_main(self.runs, self.release, self.plan, "--check", self.out), (1, "measures.tsv\n"))

    def test_check_names_a_missing_and_an_extra_file(self):
        (self.out / "order.tsv").unlink()
        (self.out / "notes.txt").write_text("SYNTHETIC\n", encoding="utf-8")
        self.assertEqual(run_main(self.runs, self.release, self.plan, "--check", self.out), (1, "order.tsv\nnotes.txt\n"))

    def test_the_checksums_cover_the_readme(self):
        (self.out / "README.md").write_text("# SYNTHETIC record\n", encoding="utf-8")
        self.assertEqual(run_main(self.runs, self.release, self.plan, "--check", self.out), (1, "SHA256SUMS\n"))
        self.assertEqual(run_main(self.runs, self.release, self.plan, "--out", self.out), (0, ""))
        self.assertEqual(run_main(self.runs, self.release, self.plan, "--check", self.out), (0, ""))
        sums = (self.out / "SHA256SUMS").read_text(encoding="utf-8").splitlines()
        names = [line.split("  ")[1] for line in sums]
        self.assertEqual(names, sorted(set(files_of(self.out)) - {"SHA256SUMS"}))
        for line in sums:
            digest, name = line.split("  ")
            self.assertEqual(hashlib.sha256((self.out / name).read_bytes()).hexdigest(), digest, name)

    def test_the_command_line_fails_on_a_changed_record(self):
        (self.out / "paired.tsv").write_text("SYNTHETIC\n", encoding="utf-8")
        env = {"PATH": "/usr/bin:/bin"}
        done = subprocess.run([sys.executable, "-I", "-B", str(SCRIPT), str(self.runs), str(self.release), str(self.plan),
                               "--check", str(self.out)], env=env, capture_output=True, text=True, check=False)
        self.assertEqual((done.returncode, done.stdout, done.stderr), (1, "paired.tsv\n", ""))


if __name__ == "__main__":
    unittest.main()
