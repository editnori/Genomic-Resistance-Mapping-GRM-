"""The committed export of stage 2's one real reading and the run's facts beside it: each file's sha256 against
SHA256SUMS, the sums against the values the run printed, and run.tsv against the export it describes.

The export holds counts and estimates in whole parts per million, and no genome id. This module reads it and
reads no score."""

import fractions
import hashlib
import unittest

from tests.matched_accuracy import stage2_common as c

# The 6 sums the run printed, equal to the 6 lines of the export's SHA256SUMS.
EXPORT_SHA256 = {
    "cells.tsv": "bd7c1ed20c55ce19b95bf27299e1ffa0bd9965ebc054227b292b223d41dbe7ff",
    "inputs.tsv": "dbeb7577ca538c81dd99bf8accc452e76ca10274801318063106e362855cf57e",
    "interval.tsv": "18531cf1a2f7ae625227c3eea8a1dae8822b2319130e24f1e0238273f8088f61",
    "measure.tsv": "fd440424b9beeb7c5857c8871788623d1df2338d7f402422f2e01d2ecdb27bf1",
    "null_fits.tsv": "e1ae5271a07bc5a4f99b9794d5355f51585d6804f978ce333bc50119dcfdc252",
    "typing.tsv": "e974d6cbb38119866e114c048b19e78749008f7a147f0f259cfcd86d503c5590",
}
# The merge commit ruling R66 names, at which the reading ran.
MERGE = "338695a85738375546e1904dddd8672b794501d3"
RUN_FACTS = 19
SCAN_FACTS = 11


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run_facts():
    """{(part, fact): value} of run.tsv, split on tabs: three fields a line, a header, no fact twice."""
    lines = (c.STAGE2_COPY / "run.tsv").read_text(encoding="utf-8").split("\n")
    assert lines[0] == "part\tfact\tvalue" and lines[-1] == "", lines[:1]
    out = {}
    for line in lines[1:-1]:
        part, fact, value = line.split("\t")
        assert (part, fact) not in out, (part, fact)
        out[(part, fact)] = value
    return out


def whole(value):
    """A count as run.tsv writes it, with commas between thousands."""
    assert value == f"{int(value.replace(',', '')):,}", value
    return int(value.replace(",", ""))


def inputs():
    """{input: count} of the export's inputs table."""
    text = (c.STAGE2_COPY / "inputs.tsv").read_text(encoding="utf-8")
    return {line.split("\t")[0]: line.split("\t")[2] for line in text.split("\n")[1:] if line.count("\t") == 3}


class Export(unittest.TestCase):
    def test_the_directory_holds_the_six_files_their_sums_and_the_runs_facts(self):
        self.assertEqual(sorted(p.name for p in c.STAGE2_COPY.iterdir()),
                         sorted(list(EXPORT_SHA256) + ["SHA256SUMS", "run.tsv"]))
        self.assertEqual(sorted(EXPORT_SHA256), sorted(c.EXPORT))

    def test_each_files_sha256_is_its_line_of_sha256sums(self):
        lines = (c.STAGE2_COPY / "SHA256SUMS").read_text(encoding="utf-8").split("\n")
        self.assertEqual(lines[-1], "")
        listed = dict(reversed(line.split("  ")) for line in lines[:-1])
        self.assertEqual(listed, EXPORT_SHA256)
        for name, digest in listed.items():
            self.assertEqual(sha(c.STAGE2_COPY / name), digest, name)

    def test_inputs_ends_with_the_line(self):
        text = (c.STAGE2_COPY / "inputs.tsv").read_text(encoding="utf-8")
        self.assertTrue(text.endswith("\n" + c.stage2.OPENED + "\n"))


class RunFacts(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.facts = run_facts()

    def test_one_row_a_fact_of_the_run_and_of_the_scan(self):
        parts = [p for p, _ in self.facts]
        self.assertEqual((parts.count("run"), parts.count("scan"), len(parts)), (RUN_FACTS, SCAN_FACTS,
                                                                                 RUN_FACTS + SCAN_FACTS))

    def test_the_run_ran_once_at_the_merge_and_ended_with_exit_0_and_no_stop(self):
        f = self.facts
        self.assertEqual(f[("run", "Code")], f"main {MERGE}")
        self.assertEqual([f[("run", k)] for k in ("Runs", "Exit", "Stopped lines", "Traceback lines",
                                                  "Changed files in the clone")], ["1", "0", "0", "0", "0"])
        self.assertEqual(f[("run", 'The line "scores opened: yes"')], "printed 1 time")
        self.assertEqual(f[("run", "Pinned inputs with an equal sha256")], "6 of 6")

    def test_the_export_facts_are_the_committed_files(self):
        f = self.facts
        files = [(c.STAGE2_COPY / name).read_bytes() for name in EXPORT_SHA256]
        self.assertEqual(whole(f[("run", "Export files")]), len(files))
        self.assertEqual(whole(f[("run", "Export lines")]), sum(data.count(b"\n") for data in files))
        self.assertEqual(whole(f[("run", "Export bytes")]), sum(len(data) for data in files))
        self.assertEqual(whole(f[("run", "Full stops in the export")]), sum(data.count(b".") for data in files))

    def test_the_reports_and_numpy_are_the_inputs_table(self):
        f, i = self.facts, inputs()
        self.assertEqual(whole(f[("run", "Reports in the directory")]), int(i["reports"]))
        self.assertEqual(f[("run", "numpy")], ".".join(i[f"numpy version, {p}"] for p in ("major", "minor", "micro")))
        self.assertEqual(f[("run", "numpy")], c.stage2.NUMPY_VERSION)

    def test_the_scan_found_nothing_and_its_control_was_found_in_every_text(self):
        f = self.facts
        for fact in ("Whole ids found", "Ids found inside a longer run of digits", "Accessions found",
                     "Tokens of an id's shape", "Tokens of an accession's shape"):
            self.assertEqual(whole(f[("scan", fact)]), 0, fact)
        texts = whole(f[("scan", "Texts")])
        self.assertEqual(f[("scan", "Control: a planted id and a planted accession")],
                         f"found in {texts} of {texts} texts")
        # The 6 export files and the run's printout of 3 lines.
        self.assertEqual(texts, len(EXPORT_SHA256) + 1)

    def test_the_scans_texts_are_the_export_and_the_three_lines_the_run_printed(self):
        # The three printed lines rebuilt from the export: the line, then said()'s line for each reading.
        measure = {r["reading"]: r for r in c.read_tsv(c.STAGE2_COPY / "measure.tsv") if r["level"] == "overall"}
        bounds = {r["reading"]: r for r in c.read_tsv(c.STAGE2_COPY / "interval.tsv") if r["read"] == "primary"}
        result = {}
        for name in c.stage2.READINGS:
            o, b = measure[name], bounds[name]
            result[name] = {"primary": {
                "levels": {"overall": {"all": dict({k: int(o[k]) for k in ("groups", "pairs", "datasets", "species")},
                                                   value=fractions.Fraction(int(o["estimate_ppm"]), c.stage2.PPM))}},
                "interval": tuple(fractions.Fraction(int(b[k]), c.stage2.PPM) for k in ("lower_ppm", "upper_ppm")),
                "species": int(b["species"])}}
        printed = [c.stage2.OPENED] + c.stage2.said(result)
        self.assertEqual(whole(self.facts[("scan", "Bytes of the 7 texts")]),
                         whole(self.facts[("run", "Export bytes")]) + sum(len(x.encode("utf-8")) + 1 for x in printed))

    def test_the_times_are_utc_and_in_order(self):
        f = self.facts
        times = [f[("run", "Start")], f[("run", "End")], f[("scan", "Start")]]
        for t in times:
            self.assertRegex(t, r"^2026-11-27T\d\d:\d\d:\d\dZ$")
        self.assertEqual(times, sorted(times))


if __name__ == "__main__":
    unittest.main()
