"""fingerprint.py: FNV-1a known answers, and that every input `grm learn` hashes changes the
fingerprint while what it does not read leaves it alone.

usage: FOLDGEN=<foldgen binary> python3 -I scripts/mtb-4k/tests/test_fingerprint.py

The byte-exact oracle is not here: it is the 22 saved revision-7 MTB reports on Vast, each of
which fingerprint.py must reproduce from the frozen 2k root (the lane receipt records that run).
"""
import json
import os
import shutil
import struct
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
import fingerprint as fp  # noqa: E402
import fixtures as fx  # noqa: E402

DID = "mycobacterium_tuberculosis__amikacin"


class Fnv(unittest.TestCase):
    def test_published_fnv1a_values_as_in_src_fnv_rs(self):
        for text, digest in [("", 0xCBF29CE484222325), ("a", 0xAF63DC4C8601EC8C), ("foobar", 0x85944171F73967E8)]:
            h = fp.Fnv()
            h.bytes(text.encode())
            self.assertEqual(h.h, digest, text)
        self.assertEqual(fp.Fnv().hex(), "cbf29ce484222325")


class Inputs(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.rules = fx.make_rules({"prefix": ("amikacin",)})
        cls.root = fx.build_cohort(Path(cls.tmp.name) / "pristine", "prefix", cls.rules)

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def key(self, root, config=None):
        return fp.Inputs(root / "data", fx.vc.SPECIES).fingerprint(DID, config or fx.CONFIG)

    def copy(self, label):
        root = Path(self.tmp.name) / label
        shutil.copytree(self.root, root)
        return root

    def edit_json(self, path, fn):
        value = json.loads(path.read_text())
        fn(value)
        path.write_text(json.dumps(value))

    def edit_phenotype(self, root, fn, genome="g000002"):
        path = root / "data/curated/phenotypes.jsonl"
        lines = []
        for line in path.read_text().splitlines():
            p = json.loads(line)
            if p["genome_id"] == genome and p["antibiotic"] == "amikacin":
                fn(p)
            lines.append(json.dumps(p))
        path.write_text("\n".join(lines) + "\n")

    def test_the_fixture_report_carries_the_recomputed_fingerprint(self):
        report = json.loads((self.root / "data/species" / fx.vc.SPECIES / "results" / f"{DID}.json").read_text())
        self.assertEqual(report["input_fingerprint"], self.key(self.root))
        inputs = fp.Inputs(self.root / "data", fx.vc.SPECIES)
        self.assertEqual(inputs.class_counts(inputs.datasets[DID]), (30, 30))

    def test_every_input_the_evaluation_reads_changes_it(self):
        sp = lambda root: root / "data/species" / fx.vc.SPECIES  # noqa: E731
        base = self.key(self.root)
        edits = {
            "a matrix word": lambda r: (sp(r) / "matrix/kmers.bin").write_bytes(b"\7" + (sp(r) / "matrix/kmers.bin").read_bytes()[1:]),
            "a support count": lambda r: (sp(r) / "matrix/support.bin").write_bytes(struct.pack("<I", 9) + (sp(r) / "matrix/support.bin").read_bytes()[4:]),
            "matrix summary min_support": lambda r: self.edit_json(sp(r) / "matrix/matrix.json", lambda m: m.update(min_support=3)),
            "a kset digest": lambda r: self.edit_json(sp(r) / "matrix/matrix.json", lambda m: m["kset_digests"].__setitem__(0, "f" * 16)),
            "a reported call": lambda r: self.edit_phenotype(r, lambda p: p.update(reported="S")),
            "an MIC": lambda r: self.edit_phenotype(r, lambda p: p.update(mic={"log2": 3, "censor": "Exact"})),
            "MIC censoring": lambda r: (self.edit_phenotype(r, lambda p: p.update(mic={"log2": 3, "censor": "Above"}))),
            "an MIC-inferred label on a matrix genome": lambda r: self.edit_phenotype(r, lambda p: p.update(label_from_mic=True)),
            "a lineage label": lambda r: self.edit_json(sp(r) / "structure.json", lambda s: s["lineages"]["labels"].__setitem__(2, 99)),
            "the lineage threshold": lambda r: self.edit_json(sp(r) / "structure.json", lambda s: s["lineages"].update(threshold=0.001)),
            "a gene call's containment": lambda r: self.edit_json(sp(r) / "counted.json", lambda c: c["genomes"][2]["genes"][0].update(containment=0.95)),
            "a gene call's allele": lambda r: self.edit_json(sp(r) / "counted.json", lambda c: c["genomes"][2]["genes"][0].update(allele="other")),
            "the class minimum": lambda r: self.edit_json(r / "data/curated/datasets.json", lambda ds: next(d for d in ds if d["id"] == DID).update(min_class=24)),
            "the breakpoint": lambda r: self.edit_json(r / "data/curated/datasets.json", lambda ds: next(d for d in ds if d["id"] == DID).update(breakpoint={"r_above": 2, "agreement": 0.9, "n": 40})),
        }
        seen = {base}
        for label, edit in edits.items():
            with self.subTest(label):
                root = self.copy("in-" + label.replace(" ", "-").replace("'", ""))
                edit(root)
                got = self.key(root)
                self.assertNotEqual(got, base, label)
                seen.add(got)
        for label, change in [("seed", lambda c: c.update(seed=7)), ("lasso_lambdas", lambda c: c.update(lasso_lambdas=31)),
                              ("forest subsample", lambda c: c["forest"].update(subsample=0.5)), ("gbt max_depth", lambda c: c["gbt"].update(max_depth=5))]:
            config = json.loads(json.dumps(fx.CONFIG))
            change(config)
            self.assertNotEqual(self.key(self.root, config), base, label)

    def test_a_catalog_version_that_moved_from_the_gene_calls_is_a_stop(self):
        root = self.copy("catalog-moved")
        (root / "data/catalog/version.txt").write_text("other\n")
        with self.assertRaises(ValueError):
            fp.Inputs(root / "data", fx.vc.SPECIES)

    def test_what_the_evaluation_does_not_read_leaves_it_alone(self):
        base = self.key(self.root)
        root = self.copy("unread")
        sp = root / "data/species" / fx.vc.SPECIES
        self.edit_json(sp / "counted.json", lambda c: c["genomes"][2].update(n_kmers=12345, stats={"contigs": 9}))
        with (root / "data/curated/phenotypes.jsonl").open("a") as f:
            f.write(json.dumps({"genome_id": "g000002", "species": "Another species", "antibiotic": "amikacin", "reported": "R",
                                "conflict": False, "mic": None, "label_from_mic": False}) + "\n")
            f.write(json.dumps({"genome_id": "not-in-matrix", "species": fx.SPECIES_NAME, "antibiotic": "amikacin", "reported": "R",
                                "conflict": False, "mic": None, "label_from_mic": False}) + "\n")
        self.edit_json(root / "data/curated/datasets.json", lambda ds: next(d for d in ds if d["id"] == DID).update(drug_class="other", n_resistant=999))
        self.assertEqual(self.key(root), base)


class RevisionSevenView(unittest.TestCase):
    """L5: the plan's revision-7 pre-check builds its view with --make-view, then --check
    --expect-reports. The fixture splits one cohort root as the frozen roots are split: the matrix,
    counted.json and structure.json in one (the preflight root), curated, catalog and the reports in
    the other (the candidate)."""

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        base = Path(cls.tmp.name)
        rules = fx.make_rules({"prefix": ("amikacin", "isoniazid")})
        whole = fx.build_cohort(base / "whole", "prefix", rules) / "data"
        sp = f"species/{fx.vc.SPECIES}"
        cls.pre, cls.cand = base / "pre" / "data", base / "cand" / "data"
        shutil.copytree(whole / sp / "matrix", cls.pre / sp / "matrix")
        for name in ("counted.json", "structure.json"):
            shutil.copyfile(whole / sp / name, cls.pre / sp / name)
        shutil.copytree(whole / sp / "results", cls.cand / sp / "results")
        for name in ("curated", "catalog"):
            shutil.copytree(whole / name, cls.cand / name)
        (base / "empty" / "data" / sp / "results").mkdir(parents=True)
        for name in ("curated", "catalog"):
            shutil.copytree(whole / name, base / "empty" / "data" / name)
        cls.empty = base / "empty" / "data"

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def run_cli(self, *args):
        return subprocess.run([sys.executable, "-I", "-B", str(HERE.parent / "fingerprint.py"), *map(str, args)],
                              cwd=self.tmp.name, capture_output=True, text=True)

    def view(self, label):
        return Path(self.tmp.name) / label / "rev7-view" / "data"

    def test_the_view_holds_links_only_and_every_report_matches(self):
        view, out = self.view("good"), Path(self.tmp.name) / "good.json"
        p = self.run_cli(view, "--make-view", "--matrix-from", self.pre, "--inputs-from", self.cand, "--check", "--expect-reports", 2, "--out", out)
        self.assertEqual(p.returncode, 0, p.stderr)
        got = json.loads(out.read_text())
        self.assertEqual((len(got["reports"]), got["matched"]), (2, 2))
        self.assertEqual(len(got["view"]), 6)
        links = sorted(str(q.relative_to(view)) for q in view.rglob("*") if q.is_symlink())
        dirs = sorted(str(q.relative_to(view)) for q in view.rglob("*") if q.is_dir() and not q.is_symlink())
        self.assertEqual(links, sorted(got["view"]))
        self.assertEqual(dirs, ["species", f"species/{fx.vc.SPECIES}"])
        for rel, target in got["view"].items():
            self.assertEqual(os.readlink(view / rel), target)
        # Built once: an existing view is refused, and nothing in it changes.
        p = self.run_cli(view, "--make-view", "--matrix-from", self.pre, "--inputs-from", self.cand, "--check")
        self.assertEqual(p.returncode, 2)
        self.assertIn("exists", p.stderr)

    def test_a_count_other_than_expected_or_an_empty_results_dir_fails(self):
        p = self.run_cli(self.view("count"), "--make-view", "--matrix-from", self.pre, "--inputs-from", self.cand, "--check", "--expect-reports", 22)
        self.assertEqual(p.returncode, 1)
        self.assertIn("not the expected 22", p.stderr)
        p = self.run_cli(self.view("empty"), "--make-view", "--matrix-from", self.pre, "--inputs-from", self.empty, "--check")
        self.assertEqual(p.returncode, 0)  # --check alone passes an empty results directory ...
        p = self.run_cli(self.view("empty"), "--check", "--expect-reports", 22)
        self.assertEqual(p.returncode, 1)  # ... --expect-reports does not

    def test_any_existing_entry_at_the_view_path_is_refused_and_left_alone(self):
        # N4: an empty directory (which mkdir(parents=True) alone would build into) and a dangling link.
        empty = self.view("existing-dir")
        empty.mkdir(parents=True)
        dangling = self.view("existing-link")
        dangling.parent.mkdir(parents=True)
        os.symlink(Path(self.tmp.name) / "nowhere", dangling)
        for view in (empty, dangling):
            with self.subTest(view.parent.parent.name):
                p = self.run_cli(view, "--make-view", "--matrix-from", self.pre, "--inputs-from", self.cand, "--check")
                self.assertEqual(p.returncode, 2, p.stderr)
                self.assertIn("exists; the view is built once", p.stderr)
        self.assertEqual(os.listdir(empty), [])
        self.assertEqual(os.readlink(dangling), str(Path(self.tmp.name) / "nowhere"))

    def test_reports_from_takes_the_results_of_a_third_root(self):
        # Plan revision 6, section 7b: revision N's reports against revision N's matrix and the rerun's inputs.
        view, out = self.view("third"), Path(self.tmp.name) / "third.json"
        p = self.run_cli(view, "--make-view", "--matrix-from", self.pre, "--inputs-from", self.empty, "--reports-from", self.cand,
                         "--check", "--expect-reports", 2, "--out", out)
        self.assertEqual(p.returncode, 0, p.stderr)
        got = json.loads(out.read_text())
        self.assertEqual((len(got["reports"]), got["matched"]), (2, 2))
        sp = f"species/{fx.vc.SPECIES}"
        self.assertEqual(got["view"][f"{sp}/results"], str(self.cand / sp / "results"))
        self.assertEqual(got["view"]["curated"], str(self.empty / "curated"))
        p = self.run_cli(self.view("third-alone"), "--check", "--reports-from", self.cand)
        self.assertEqual(p.returncode, 2)
        self.assertIn("--reports-from needs --make-view", p.stderr)
        self.assertFalse(os.path.lexists(self.view("third-alone")))

    def test_a_missing_source_or_a_partial_flag_set_stops_before_anything_is_made(self):
        view = self.view("missing")
        p = self.run_cli(view, "--make-view", "--matrix-from", self.cand, "--inputs-from", self.cand, "--check")
        self.assertEqual(p.returncode, 2)
        self.assertIn("view sources missing", p.stderr)
        self.assertFalse(os.path.lexists(view))
        p = self.run_cli(view, "--make-view", "--matrix-from", self.pre, "--check")
        self.assertEqual(p.returncode, 2)
        self.assertIn("--make-view needs", p.stderr)
        self.assertFalse(os.path.lexists(view))


if __name__ == "__main__":
    unittest.main()
