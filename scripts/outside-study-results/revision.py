"""revision.py ROOT FIGURES [--check FILE]
revision.py ROOT --figures REPORTS

ROOT is reference/outside-study. Each run of record read its new estimates
against the lineage-split figure of the report it pinned. A later release can
change that figure and keep the final model, and then the estimates are still
the served model's while the readings are against a figure the site no longer
shows. This script reads the same estimates against a later release's figures.

With --figures it prints FIGURES. REPORTS is the `datasets/` directory of an
unpacked release. For each dataset and learner of ROOT's `metrics.tsv` files
it writes one row: the report's sha256 and evaluation revision, and the
report's lineage-split row as the kit reads it (scripts/outside-study/score.py,
`reference_row`), or the kit's note where the report has no such row. Where
that row has `called_metrics`, the figure over the fitted folds only, its five
counts follow. It stops on a report of another dataset, of another schema than
the kit's, or with an evaluation marked invalid. It checks no hash pin: the
kit pins the reports of one revision only.

Without --figures it prints the tables of the section "Revision 14" of
docs/outside-study-results.md: for each population and learner, the datasets
where the learner holds against the run's figure and against the later one,
every row read whose figure differs, measure by measure, and for each such
dataset the genomes of null-fit folds that the later figure counts.

It estimates nothing. Each new estimate and interval is the run's. A row whose
five lineage counts are the run's keeps the run's readings. A row whose counts
differ is read with the kit's own rule (`reading`) on the later figure and the
interval as `metrics.tsv` prints it, to six decimals. A later figure within
one millionth of a printed interval end cannot be read that way, and the
script stops there.

With --check FILE it prints nothing and exits 1 unless FILE holds, between the
lines `<!-- revision:begin -->` and `<!-- revision:end -->`, exactly the text
this script prints.
"""
import csv
import hashlib
import importlib.util
import io
import json
import pathlib
import re
import sys

HERE = pathlib.Path(__file__).resolve().parent


def _module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


tables = _module("outside_study_tables", HERE / "tables.py")
score = _module("outside_study_score", HERE.parent / "outside-study" / "score.py")

BEGIN, END = "<!-- revision:begin -->", "<!-- revision:end -->"
COUNTS = ("n", "tp", "fn", "fp", "tn")
MEASURES = ("sensitivity", "specificity")
COLUMNS = (
    ["dataset", "learner", "report_sha256", "evaluation_revision"]
    + [f"lineage_{k}" for k in COUNTS + MEASURES]
    + ["lineage_note"]
    + [f"called_{k}" for k in COUNTS]
)
# metrics.tsv prints an interval end to six decimals.
MARGIN = 1e-6


def called_row(report, learner):
    """The `called_metrics` of the report's lineage-split row for a learner, or None."""
    split = next((s for s in report["splits"] if s["scheme"] == "lineage"), None)
    rows = [] if split is None else split.get("learners") or []
    row = next((r for r in rows if r["learner"] == score.REPORT_LEARNER[learner]), None)
    return None if row is None else row.get("called_metrics")


def figures(root, reports):
    """FIGURES as text: one row for each dataset and learner of ROOT's metrics files."""
    data = tables.load(pathlib.Path(root))
    wanted = sorted({(r["dataset"], r["learner"]) for runs in data.values() for _, _, metrics, _ in runs.values() for r in metrics})
    out = io.StringIO()
    w = csv.writer(out, delimiter="\t", lineterminator="\n")
    w.writerow(COLUMNS)
    loaded = report = digest = None
    for dataset, learner in wanted:
        if dataset != loaded:
            raw = (pathlib.Path(reports) / f"{dataset}.json").read_bytes()
            loaded, report, digest = dataset, json.loads(raw), hashlib.sha256(raw).hexdigest()
            if report.get("id") != dataset:
                raise SystemExit(f"{dataset}.json is the report of {report.get('id')!r}")
            if report.get("schema") != score.REPORT_SCHEMA or report["evaluation"].get("valid") is not True:
                raise SystemExit(f"{dataset}.json: schema {report.get('schema')!r}, evaluation valid {report['evaluation'].get('valid')!r}")
        ref, note = score.reference_row(report, "lineage", learner)
        called = called_row(report, learner)
        w.writerow(
            [dataset, learner, digest, report["evaluation"]["revision"]]
            + ["" if ref is None else repr(ref[k]) for k in COUNTS + MEASURES]
            + ["" if note is None else note]
            + ["" if called is None else repr(called[k]) for k in COUNTS]
        )
    return out.getvalue()


def load_figures(path):
    """(evaluation revision, {(dataset, learner): row}) of a FIGURES file."""
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if reader.fieldnames != COLUMNS:
            raise SystemExit(f"{path}: the columns are not the ones this script writes")
        rows = list(reader)
    by_key = {(r["dataset"], r["learner"]): r for r in rows}
    revisions = {r["evaluation_revision"] for r in rows}
    if len(by_key) != len(rows) or len(revisions) != 1:
        raise SystemExit(f"{path}: a row is listed twice, or the reports are not of one evaluation revision")
    return int(revisions.pop()), by_key


def counts(row):
    """A row's five lineage counts as text, the same in `metrics.tsv` and in FIGURES."""
    return tuple(row[f"lineage_{k}"] for k in COUNTS)


def reread(r, measure, figure):
    """The kit's reading of a run's interval against another figure of the report."""
    lo, hi = float(r[measure + "_lo"]), float(r[measure + "_hi"])
    if min(abs(figure - lo), abs(figure - hi)) < MARGIN:
        raise SystemExit(f"{r['dataset']} {r['learner']} {measure}: the figure is within one millionth of a printed interval end")
    return score.reading(figure, lo, hi)


def later(data, figs):
    """{population: [row]} of the rows the pages read, each with its later figure and readings.

    A row here is the run's `metrics.tsv` row with four keys added: `species`,
    `figure` (its FIGURES row), `differs` and `later` ({measure: reading})."""
    out = {}
    for pop, runs in data.items():
        out[pop] = []
        for species, (_, _, metrics, _) in runs.items():
            for r in metrics:
                if r["learner"] not in tables.LEARNERS or tables.holds(r) is None:
                    continue
                f = figs.get((r["dataset"], r["learner"]))
                if f is None or f["lineage_sensitivity"] == "":
                    raise SystemExit(f"{r['dataset']} {r['learner']}: read by the run, and the later reports give no lineage-split figure")
                differs = counts(f) != counts(r)
                readings = {m: reread(r, m, float(f["lineage_" + m])) if differs else r["reading_" + m] for m in MEASURES}
                out[pop].append({**r, "species": species, "figure": f, "differs": differs, "later": readings})
    return out


def holds_later(row):
    return all(row["later"][m] != "below" for m in MEASURES)


def null_folds(rows):
    """[(row, genomes, in null-fit folds, called resistant there, called susceptible there)]

    One entry for each dataset and learner read whose FIGURES row has `called_`
    counts, the smallest share of genomes first. The genomes of null-fit folds
    are the ones the figure over every genome counts and the figure over the
    fitted folds does not."""
    out = {}
    for pop, _ in tables.POPULATIONS:
        for r in rows.get(pop, []):
            f = r["figure"]
            if f["called_n"] == "" or (r["dataset"], r["learner"]) in out:
                continue
            added = {k: int(f[f"lineage_{k}"]) - int(f[f"called_{k}"]) for k in COUNTS}
            if min(added.values()) < 0 or added["n"] != added["tp"] + added["fn"] + added["fp"] + added["tn"]:
                raise SystemExit(f"{r['dataset']} {r['learner']}: the figure over the fitted folds is no part of the figure over every genome")
            out[r["dataset"], r["learner"]] = (r, int(f["lineage_n"]), added["n"], added["tp"] + added["fp"], added["fn"] + added["tn"])
    return sorted(out.values(), key=lambda e: (e[2] / e[1], e[0]["dataset"], e[0]["learner"]))


def figure_cell(row, measure):
    """'40.8%, 29 of 71': a report figure with its counts, from a metrics.tsv row or a FIGURES row."""
    hit, miss = ("tp", "fn") if measure == "sensitivity" else ("tn", "fp")
    hit, miss = int(row["lineage_" + hit]), int(row["lineage_" + miss])
    return f"{tables.pct(row['lineage_' + measure])}, {tables.num(hit)} of {tables.num(hit + miss)}"


def build(root, figures_path):
    data = tables.load(pathlib.Path(root))
    revision, figs = load_figures(figures_path)
    rows = later(data, figs)
    out = ["**Datasets where a learner holds, against each report figure**", ""]
    out.append(f"| Population | Learner | Datasets read | Report figure differs | Holds, runs of record | Holds, revision {revision} |")
    out.append("|---|---|---:|---:|---:|---:|")
    for pop, _ in tables.POPULATIONS:
        for learner in tables.LEARNERS:
            read = [r for r in rows.get(pop, []) if r["learner"] == learner]
            if read:
                out.append(
                    f"| {pop} | {tables.NAMES[learner]} | {len(read)} | {sum(r['differs'] for r in read)} "
                    f"| {sum(tables.holds(r) for r in read)} | {sum(holds_later(r) for r in read)} |"
                )
    out += ["", "**Rows read whose report figure differs**", ""]
    out.append("Each new estimate is followed by its 95% interval and its counts, and each report figure by its counts.")
    out.append("")
    out.append(f"| Population | Species | Drug | Learner | Measure | New estimate | Report, runs of record | Reading | Report, revision {revision} | Reading |")
    out.append("|---|---|---|---|---|---|---|---|---|---|")
    for pop, _ in tables.POPULATIONS:
        for r in rows.get(pop, []):
            if not r["differs"]:
                continue
            for m in MEASURES:
                out.append(
                    f"| {pop} | {tables.title(r['species'])} | {tables.drug(r['dataset'])} | {tables.NAMES[r['learner']]} | {m} | {tables.estimate(r, m)} "
                    f"| {figure_cell(r, m)} | {r['reading_' + m]} | {figure_cell(r['figure'], m)} | {r['later'][m]} |"
                )
    out += ["", f"**Genomes of null-fit folds in the revision {revision} figure**", ""]
    out.append("Each dataset read whose report keeps a second figure, over the fitted folds only. The calls are the ones the genomes of null-fit folds get.")
    out.append("")
    out.append("| Species | Drug | Learner | Genomes | In null-fit folds | Share | Called resistant | Called susceptible |")
    out.append("|---|---|---|---:|---:|---:|---:|---:|")
    for r, n, null, resistant, susceptible in null_folds(rows):
        out.append(
            f"| {tables.title(r['species'])} | {tables.drug(r['dataset'])} | {tables.NAMES[r['learner']]} | {tables.num(n)} | {tables.num(null)} "
            f"| {tables.pct(null / n)} | {tables.num(resistant)} | {tables.num(susceptible)} |"
        )
    return "\n".join(out).rstrip("\n") + "\n"


def check(root, figures_path, document):
    """True when the document's revision block is exactly what `build` prints."""
    doc = pathlib.Path(document).read_text(encoding="utf-8")
    m = re.search(re.escape(BEGIN) + r"\n(.*?)" + re.escape(END), doc, re.S)
    return bool(m) and m.group(1).strip("\n") == build(root, figures_path).strip("\n")


def main(argv):
    args = argv[1:]
    if len(args) == 3 and args[1] == "--figures":
        sys.stdout.write(figures(args[0], args[2]))
        return 0
    if len(args) == 4 and args[2] == "--check":
        if check(args[0], args[1], args[3]):
            return 0
        print("the revision tables in the document differ from what the files give", file=sys.stderr)
        return 1
    if len(args) != 2 or args[1].startswith("--"):
        raise SystemExit(__doc__)
    sys.stdout.write(build(args[0], args[1]))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
