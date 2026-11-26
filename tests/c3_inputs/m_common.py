"""SYNTHETIC fits, score files, kept reports and fold results for the run on the datasets of M (scripts/c3-inputs).

Every genome id, class, fold and score here is made up in the test; the ids
are numbers of the form a genome id takes, so the id checks have something
to find. Each result's per-fold counts are taken from its genomes' scores by
the rule grm's cutoff uses, so a fixture's counts come back from its scores
unless a test changes one.
"""

import json
import pathlib

from tests.c3_inputs import common

checks = common.load("m_checks")
LEARNERS = common.LEARNERS
SPECIES = common.SPECIES
M_SET = "synthia_coli__drug_m"
FOLDS = 5


def ids(n, base=90007):
    return [f"{base}.{i + 1}" for i in range(n)]


def genomes(n=80, resistant_every=2, folds=FOLDS, base=90007):
    """[(genome id, resistant, fold, [five scores])], made up: scores spread over [0, 1], some exactly 0.5 and 0."""
    out = []
    for i, g in enumerate(ids(n, base)):
        r = i % resistant_every == 0
        scores = []
        for j in range(len(LEARNERS)):
            s = ((i * 37 + j * 11) % 101) / 100
            s = min(s, 1.0)
            if (i + j) % 17 == 0:
                s = 0.5
            scores.append(s if not r else min(1.0, s + 0.25))
        out.append((g, r, i % folds, scores))
    return out


def pick_from(gs, dataset=M_SET, scheme="lineage", seed=42, folds=FOLDS):
    """A result as grm learn --exp-pick --exp-pick-inner 0 writes it, its counts taken from the genomes' scores."""
    doc = common.pick(dataset=dataset, scheme=scheme, seed=seed)
    n, n_r = len(gs), sum(1 for g in gs if g[1])
    recs = []
    for f, (got, test, test_r) in enumerate(checks.fold_calls(gs, folds)):
        rec = {"fold": f, "train": n - test, "train_resistant": n_r - test_r, "test": test, "test_resistant": test_r,
               "outer": [dict(learner=l, **c) for l, c in zip(LEARNERS, got)],
               "inner": {"status": "not run", "reason": None, "inner_folds": [], "learners": [], "tied": [],
                         "pick": None}}
        recs.append(rec)
    doc["folds"] = recs
    doc["tested"], doc["tested_resistant"] = n, n_r
    doc["unavailable_reason"] = None
    return doc


def scores_file(gs, dataset=M_SET, seed=42, arm="k31", folds=FOLDS):
    """What grm learn --exp-scores writes for one drawn fit."""
    return {"format": "grm-c3-scores-1", "evaluation_revision": 14, "dataset": dataset, "scheme": "lineage",
            "seed": seed, "arm": arm, "learners": list(LEARNERS), "folds": folds, "unavailable_reason": None,
            "genomes": [{"genome_id": g, "resistant": r, "fold": f, "scores": s} for g, r, f, s in gs]}


def report(gs, dataset=M_SET):
    """A report's lineage split as grm learn writes it, holding each genome's score in sample order."""
    return {"id": dataset, "samples": [{"id": g, "resistant": r} for g, r, _, _ in gs],
            "splits": [{"scheme": "random", "learners": []},
                       {"scheme": "lineage", "unavailable_reason": None, "test_fold": [f for _, _, f, _ in gs],
                        "learners": [{"learner": l, "scores": [s[j] for _, _, _, s in gs]}
                                     for j, l in enumerate(LEARNERS)]}]}


def fold_result(gs, dataset=M_SET, scheme="study-whole"):
    """A study fit's fold result (grm learn --fold-file --fold-out), its samples by tested position."""
    return {"format": "grm-fold-file-result-v1", "evaluation_revision": 14, "id": dataset, "scheme": scheme,
            "samples": [{"genome_id": g, "resistant": r, "fold": f} for g, r, f, _ in gs],
            "split": {"scheme": scheme, "learners": [{"learner": l, "scores": [s[j] for _, _, _, s in gs]}
                                                    for j, l in enumerate(LEARNERS)]}}


def write_json(path, doc):
    path = pathlib.Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(doc))


def write_tsv(path, columns, rows):
    path = pathlib.Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        f.write("\t".join(columns) + "\n")
        for r in rows:
            f.write("\t".join(str(r[c]) for c in columns) + "\n")


def lay_out(top, rows, gs_of, refused=(), study_gs_of=None):
    """Write each expect row's result under top/picks, its scores under top/scores or top/study, and the
    default's and the twin's kept reports under top/reports. gs_of(row) gives a row's genomes; refused names
    files written as refusals for no candidate. Returns (expect, picks, scores, reports, study)."""
    top = pathlib.Path(top)
    expect = []
    for r in rows:
        e = dict(zip(common.COLUMNS, [str(x) for x in r]))
        expect.append(e)
        gs = gs_of(e)
        if e["file"] in refused:
            write_json(top / "picks" / (e["file"][: -len(".pick.json")] + ".refused.json"),
                       {"format": "grm-c3-refused-1", "dataset": e["dataset"], "arm": e["arm"], "scheme": e["scheme"],
                        "seed": int(e["seed"]), "reason": "no candidate"})
            continue
        write_json(top / "picks" / e["file"], pick_from(gs, e["dataset"], e["scheme"], int(e["seed"])))
        if e["split"] == "lineage":
            arm = "k31" if e["arm"] == "twin" else e["arm"]
            write_json(checks.score_path(top / "scores", e), scores_file(gs, e["dataset"], int(e["seed"]), arm))
            if e["arm"] in ("k31", "twin"):
                write_json(checks.report_path(top / "reports", e), report(gs, e["dataset"]))
        else:
            write_json(checks.study_path(top / "study", e), fold_result(gs, e["dataset"], e["scheme"]))
    return expect, top / "picks", top / "scores", top / "reports", top / "study"
