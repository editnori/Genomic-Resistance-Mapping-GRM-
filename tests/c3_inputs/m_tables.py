"""SYNTHETIC tables, records and lists for scripts/c3-inputs/m_page.py's tests (test_m_page.py, test_m_record.py).

The 7 tables are m_reading.py's own, written by its floor and count forms
from test_m_reading.py's made-up species directories: 6 made-up datasets of M
in two made-up species and 8 made-up first-run datasets, no genome. The
record holds a made-up README.md and first-run-times.tsv. No dataset, figure
or time here comes from a run.
"""

from pathlib import Path

from tests.c3_inputs import common, m_common
from tests.c3_inputs import test_m_reading as mr

m_page = common.load("m_page")
M_SET = {r["dataset"] for r in mr.LISTS if r["in_M"] == "1"}
# Each class lists no genome: every dataset of M is in the second reading, the one under the floor too, since 30
# resistant and 40 susceptible genomes stay.
ALL_IN = {d: (0, 0) for d in mr.LISTED}
TIMES = "species_slug\tstep\tseconds\tsource_sha256\nsynth_c\tbuild\t72\t" + "a" * 64 + "\n"


def seven(top, listed=mr.LISTED):
    """m_reading.py's 7 tables from test_m_reading.py's species directories: <top>/run/tables. The second reading's
    file holds every dataset that passes the floor, as second.py writes it."""
    top = Path(top)
    fs = mr.folders(top / "species", listed)
    if listed[mr.UNDER] == (0, 0):
        start = len(mr.RECORD) + len(mr.M_SETS["synth_c"])
        m_common.write_tsv(fs[1] / "score" / "second" / "datasets.tsv", mr.score.DATASET_COLUMNS,
                           mr.rows_of(mr.M_SETS["synth_d"], mr.MEASURES, start, shift=4))
    mr.tables(top / "run", fs, mr.rows_of(mr.RECORD, mr.TWO, 0))
    return top / "run" / "tables"


def floor_only(top, listed=mr.LISTED):
    """The floor form's two tables alone, from floor.tsv files alone: <top>/run/tables."""
    top = Path(top)
    fs = mr.folders(top / "species", listed, scores=False)
    out = top / "run" / "tables"
    mr.reading.floor_form(out, fs, mr.LISTS)
    return out


def record(top):
    """A record directory as reference/input-arms-m stands before the run: README.md and first-run-times.tsv."""
    rec = Path(top) / "record"
    rec.mkdir(parents=True)
    (rec / "README.md").write_text("# SYNTHETIC record\n", encoding="utf-8")
    (rec / "first-run-times.tsv").write_text(TIMES, encoding="utf-8")
    return rec


def lists_file(path):
    """A lists.tsv of the made-up datasets with the two columns m_page.py reads."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("dataset\tin_M\n" + "".join(f"{r['dataset']}\t{r['in_M']}\n" for r in mr.LISTS),
                    encoding="utf-8")
    return path


def edit(record_dir, name, change, sums=True):
    """Rewrite one table of a record through change(text) -> text, and its SHA256SUMS with it unless sums is False,
    so a gate other than the sums' is the one that stops."""
    path = Path(record_dir) / name
    path.write_text(change(path.read_text(encoding="utf-8")), encoding="utf-8")
    if sums:
        (Path(record_dir) / "SHA256SUMS").write_text(m_page.sums_text(record_dir), encoding="utf-8")
