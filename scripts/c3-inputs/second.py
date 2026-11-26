"""The second reading of plan block C3's run on its datasets of M (docs/input-arms.md, The run on the datasets of M).

  python3 -I second.py <expect.tsv> <picks dir> <scores dir> <study dir> <release dir> <floor.tsv> <list dir> <out dir>

Run on the run's machine after score.py's first reading. It refits
nothing. For each dataset that <floor.tsv> enters (m_list.py's table, made
from the label table before any score was read), every result expect.tsv
lists for it is read again with the dataset's listed genomes
(<list dir>/<dataset>.txt, those whose records for its antibiotic hold an
intermediate call) left out of every outer test fold: of every arm, the
default and the twin, at every seed and on both splits. A fold that holds a
listed genome gets each learner's calls taken again from the genomes'
scores (m_checks.py's calls, as grm's cutoff makes them); a fold that holds
none keeps the counts the fit wrote. The training counts stay as fitted.
An arm's study result that grm learn refused for no candidate is the
constant model (score.py, ruling 16) on the default's study folds less the
listed genomes, with the training counts as fitted. A split that could not
be drawn stays unscored, as in the first reading.

The scores are read as m_checks.py's counts form reads them: the lineage
fits' score files (<scores dir>/<arm>/<dataset>.lineage-s<seed>.scores.json)
and the study fits' fold results (<study dir>/<arm>/<dataset>.study-whole.json).
Before any genome is left out, every result's counts must come back from
its scores, as Same counts asks; each listed genome must be a cohort genome
of every lineage score file; and the genomes each lineage score file keeps
must give floor.tsv's counts by class. Any difference stops it, writing
nothing.

A dataset that floor.tsv does not enter is not scored in this reading;
floor.tsv names it with its counts by class. Writes <out dir>/draws.tsv and
<out dir>/datasets.tsv in score.py's columns, and prints counts only.
"""

import importlib.util
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent


def _load(name, file):
    spec = importlib.util.spec_from_file_location(name, HERE / file)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


SCORE = _load("c3_score_for_second", "score.py")
CHECKS = _load("c3_checks_for_second", "m_checks.py")
LEARNERS = SCORE.LEARNERS
DEFAULT = SCORE.DEFAULT
COUNTED = ("n", "tp", "fp", "tn", "fn", "loss_sum")


class Stop(Exception):
    pass


def read_list(path):
    """The listed genome ids of one dataset, one a line."""
    ids = [line.strip() for line in Path(path).read_text().splitlines() if line.strip()]
    if len(set(ids)) != len(ids):
        raise Stop(f"{Path(path).name}: a genome is listed twice")
    return set(ids)


def again(doc, genomes, left_out):
    """A result read again with left_out out of every outer test fold; each fold's training counts as fitted."""
    if doc.get("unavailable_reason") is not None:
        return doc
    per = CHECKS.fold_calls(genomes, len(doc["folds"]), left_out)
    folds = []
    for f, (rec, (got, test, test_r)) in enumerate(zip(doc["folds"], per)):
        if any(g[2] == f and g[0] in left_out for g in genomes):
            rec = dict(rec, test=test, test_resistant=test_r,
                       outer=[dict(w, **{k: c[k] for k in COUNTED}) for w, c in zip(rec["outer"], got)])
        folds.append(rec)
    gone = [g for g in genomes if g[0] in left_out]
    return dict(doc, folds=folds, tested=doc["tested"] - len(gone),
                tested_resistant=doc["tested_resistant"] - sum(1 for g in gone if g[1]))


def floor_of(rows, lists, expect):
    """The datasets that enter, after holding floor.tsv to the list files and to expect.tsv."""
    datasets = {e["dataset"] for e in expect}
    named = [r["dataset"] for r in rows]
    if len(set(named)) != len(named) or set(named) != datasets:
        raise Stop("floor.tsv's datasets are not expect.tsv's, each once")
    for r in rows:
        d = r["dataset"]
        if d not in lists:
            raise Stop(f"{d}: no list file")
        listed = int(r["listed"])
        if len(lists[d]) != listed or int(r["listed_resistant"]) + int(r["listed_susceptible"]) != listed:
            raise Stop(f"{d}: floor.tsv lists {listed} genomes, the list file {len(lists[d])}")
        second = (int(r["resistant"]) - int(r["listed_resistant"]), int(r["susceptible"]) - int(r["listed_susceptible"]))
        if second != (int(r["second_resistant"]), int(r["second_susceptible"])):
            raise Stop(f"{d}: floor.tsv's second-reading counts are not the cohort's less the list")
        enters = int(min(second) >= int(r["floor"]))
        if int(r["floor"]) != 25 or enters != int(r["enters"]):
            raise Stop(f"{d}: floor.tsv's entry is not the floor of 25 a class")
    return {r["dataset"] for r in rows if r["enters"] == "1"}, {r["dataset"]: r for r in rows}


def second(expect, picks, scores, study, names, floor_rows, lists):
    """(draws rows, entering datasets, datasets under the floor). Raises Stop on any difference."""
    picks = Path(picks)
    enters, by = floor_of(floor_rows, lists, expect)
    rows_in = [e for e in expect if e["dataset"] in enters]
    docs, read, refused = {}, {}, {}
    for e in rows_in:
        path = picks / e["file"]
        if not path.is_file():
            if e["arm"] != DEFAULT and (picks / SCORE.refused_name(e["file"])).is_file():
                continue
            raise Stop(f"{e['file']}: no result")
        doc = CHECKS.load(path)
        try:
            genomes = CHECKS.genomes_of(e, scores, study)
        except FileNotFoundError:
            raise Stop(f"{e['file']}: its genomes' scores were not kept")
        diffs = CHECKS.compare(doc, genomes, e["file"])
        if diffs:
            raise Stop(f"{len(diffs)} counts do not come back from the scores; " + "; ".join(diffs[:3]))
        left_out = lists[e["dataset"]]
        if e["split"] == "lineage" and doc.get("unavailable_reason") is None:
            ids = {g[0] for g in genomes}
            if not left_out <= ids:
                raise Stop(f"{e['file']}: {len(left_out - ids)} listed genomes are not in the fit's cohort")
            kept = [g for g in genomes if g[0] not in left_out]
            r = by[e["dataset"]]
            got = (sum(1 for g in kept if g[1]), sum(1 for g in kept if not g[1]))
            if got != (int(r["second_resistant"]), int(r["second_susceptible"])):
                raise Stop(f"{e['file']}: the fit keeps {got[0]} R and {got[1]} S, floor.tsv "
                           f"{r['second_resistant']} and {r['second_susceptible']}")
        docs[e["file"]] = again(doc, genomes, left_out)
        read[e["file"]] = genomes
    base_of = {(e["dataset"], e["split"], e["seed"]): e["file"] for e in rows_in if e["arm"] == DEFAULT}
    for e in rows_in:
        if e["file"] in docs:
            continue
        base = base_of.get((e["dataset"], e["split"], e["seed"]))
        if base not in docs:
            raise Stop(f"{e['file']}: refused, and the default's result at its dataset, split and seed is not read")
        refused[e["file"]] = docs[base]
    draws = SCORE.score_rows(rows_in, docs, names, refused)
    under = sorted(set(by) - enters)
    return draws, sorted(enters), under


def main(argv):
    if len(argv) != 8:
        print(__doc__.strip().splitlines()[2], file=sys.stderr)
        return 2
    expect_path, picks, scores, study, release, floor_path, list_dir, out = argv
    expect = SCORE.read_tsv(expect_path)
    floor_rows = SCORE.read_tsv(floor_path)
    names = {}
    for d in sorted({e["dataset"] for e in expect}):
        with open(Path(release) / "datasets" / f"{d}.json") as f:
            rel = json.load(f)
        names[d] = (rel["species"], rel["antibiotic"])
    try:
        lists = {r["dataset"]: read_list(Path(list_dir) / f"{r['dataset']}.txt") for r in floor_rows
                 if (Path(list_dir) / f"{r['dataset']}.txt").is_file()}
        draws, enters, under = second(expect, picks, scores, study, names, floor_rows, lists)
    except Stop as e:
        print(f"second: stopped: {e}")
        return 1
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    SCORE.write_tsv(out / "draws.tsv", draws, SCORE.DRAW_COLUMNS)
    SCORE.write_tsv(out / "datasets.tsv", SCORE.dataset_rows(draws), SCORE.DATASET_COLUMNS)
    print(f"second reading: {len(enters)} of {len(enters) + len(under)} datasets enter, {len(draws)} draw rows; "
          f"under the floor of 25 a class: {', '.join(under) if under else 'none'}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
