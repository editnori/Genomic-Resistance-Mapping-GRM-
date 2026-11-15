"""The reference files the clindamycin writer reads, and what the bundle may hold.

Each file of the folder the bundle names as a source must equal its line in
the folder's SHA256SUMS, and its bytes must give that sha256; the one source
beside the folder, the release's checksum lines in reference/lineage-control,
must equal its line in that directory's README. Each of the six released
reports' sha256 the bundle carries must equal its line in the release's own
checksum lines, and none is null. The bundle holds no row per sample and no
genome id: its fields are a fixed list of scalars.
"""

import hashlib
import json
import re
import unittest

from . import common

TOP = ["page", "reference", "release", "evaluation_revision", "min_class", "threshold", "sources", "datasets"]
FIELDS = [
    "kind",
    "pmid",
    "samples",
    "study_samples",
    "disagreements",
    "study_disagreements",
    "no_pmid_samples",
    "no_pmid_disagreements",
    "other_pmid_samples",
    "other_pmid_disagreements",
    "left_resistant",
    "left_susceptible",
    "scored_without",
    "rows_passing",
    "rows_tested",
    "distinct_passing",
    "distinct_tested",
    "direction",
    "report_sha256",
]
# A BV-BRC genome id: taxon and serial joined by a dot.
GENOME_ID = re.compile(r"\b\d+\.\d+\b")


def sums():
    """SHA256SUMS of the reference folder: path to sha256."""
    lines = {}
    for line in (common.REFERENCE / "SHA256SUMS").read_text(encoding="utf-8").splitlines():
        digest, rel = line.split("  ")
        lines[rel] = digest
    return lines


def release_lines():
    """The release's own checksum line of each report: dataset id to sha256."""
    lines = {}
    for line in common.RELEASE_SUMS.read_text(encoding="utf-8").splitlines():
        digest, name = line.split("  ")
        lines[name[len("./datasets/") : -len(".json")]] = digest
    return lines


class ReferenceFiles(unittest.TestCase):
    def test_each_file_of_the_folder_gives_its_sha256sums_line(self):
        listed = sums()
        self.assertGreater(len(listed), 0)
        for rel, digest in listed.items():
            self.assertEqual(hashlib.sha256((common.REFERENCE / rel).read_bytes()).hexdigest(), digest, rel)

    def test_each_source_of_the_bundle_is_a_listed_file_with_its_listed_sha256(self):
        listed = sums()
        sources = common.bundled()["sources"]
        beside = f"../../lineage-control/{common.RELEASE_SUMS.name}"
        self.assertIn(beside, sources)
        self.assertGreater(len(sources), 1)
        for rel, digest in sources.items():
            self.assertEqual(hashlib.sha256((common.REFERENCE / rel).read_bytes()).hexdigest(), digest, rel)
            if rel == beside:
                self.assertEqual((common.REFERENCE / rel).resolve(), common.RELEASE_SUMS.resolve())
                readme = (common.RELEASE_SUMS.parent / "README.md").read_text(encoding="utf-8")
                self.assertIn(f"{digest}  {common.RELEASE_SUMS.name}", readme)
                continue
            self.assertIn(rel, listed, f"{rel} is not in SHA256SUMS")
            self.assertEqual(digest, listed[rel], rel)

    def test_each_report_sha256_is_the_releases_own_line(self):
        data = common.bundled()
        self.assertEqual(common.RELEASE_SUMS.name, f"{data['release']}.reports.sha256")
        lines = release_lines()
        self.assertEqual(len(data["datasets"]), 6)
        for d, entry in data["datasets"].items():
            self.assertIsNotNone(entry["report_sha256"], d)
            self.assertEqual(entry["report_sha256"], lines[d], d)

    def test_the_checks_own_report_lines_are_the_releases(self):
        lines = release_lines()
        kept = (common.REFERENCE / "results" / "logs" / "released-reports.sha256").read_text(encoding="utf-8").splitlines()
        self.assertGreater(len(kept), 0)
        for line in kept:
            digest, name = line.split("  ")
            self.assertEqual(digest, lines[name[: -len(".json")]], name)


class NoSampleField(unittest.TestCase):
    def test_the_bundle_holds_the_fixed_fields_and_scalars_only(self):
        data = common.bundled()
        self.assertEqual(list(data), TOP)
        self.assertTrue(all(isinstance(v, str) for v in data["sources"].values()))
        for d, entry in data["datasets"].items():
            self.assertEqual(list(entry), FIELDS, d)
            for field, value in entry.items():
                self.assertIsInstance(value, (str, int, bool, type(None)), f"{d} {field}")

    def test_no_string_of_the_bundle_is_a_genome_id(self):
        strings = []

        def walk(v):
            if isinstance(v, dict):
                for k, x in v.items():
                    strings.append(k)
                    walk(x)
            elif isinstance(v, list):
                for x in v:
                    walk(x)
            elif isinstance(v, str):
                strings.append(v)

        walk(common.bundled())
        self.assertEqual([s for s in strings if GENOME_ID.search(s)], [])

    def test_the_bundle_is_small(self):
        # Six datasets of 19 fields and ten source hashes; a row per sample would not fit.
        self.assertLess(common.BUNDLED.stat().st_size, 8192)
        self.assertEqual(len(json.loads(common.BUNDLED.read_text(encoding="utf-8"))["datasets"]), 6)


if __name__ == "__main__":
    unittest.main()
