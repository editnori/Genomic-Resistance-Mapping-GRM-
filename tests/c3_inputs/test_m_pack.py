"""scripts/c3-inputs/pack.py: what leaves the run's machine holds counts per dataset and no genome.

Every OUT here is SYNTHETIC: made-up results, tables and logs, and ids of
the form a genome id takes, written by the test (m_common.py).
"""

import contextlib
import hashlib
import io
import json
import os
import tarfile
import tempfile
import unittest
from pathlib import Path

from tests.c3_inputs import common, m_common

pack = common.load("pack")
D = m_common.M_SET


def out_dir(top):
    """An OUT as run_m.sh leaves it before the pack: work/ and flags/ hold what stays on the machine."""
    out = Path(top) / "out"
    gs = m_common.genomes(30)
    m_common.write_json(out / "results" / "picks" / "k31" / f"{D}.lineage-s42.pick.json", m_common.pick_from(gs))
    m_common.write_json(out / "results" / "inputs" / "program" / f"{D}.drawn.inputs.json", common.inputs(D))
    m_common.write_tsv(out / "results" / "m" / "floor.tsv", ["dataset", "resistant", "susceptible"],
                       [{"dataset": D, "resistant": 15, "susceptible": 15}])
    m_common.write_tsv(out / "score" / "datasets.tsv", ["dataset", "balanced_accuracy"],
                       [{"dataset": D, "balanced_accuracy": 0.75}])
    (out / "logs").mkdir(parents=True)
    (out / "logs" / "fit-k31-s42-m.log").write_text("read 30 genomes; [genome id] has no call; grm 1.98.1\n")
    (out / "binary.txt").write_text("grm " + "a" * 64 + "\n")
    (out / "SHA256SUMS").write_text("b" * 64 + "  binary.txt\n")
    m_common.write_json(out / "work" / "scores" / "k31" / f"{D}.lineage-s42.scores.json", m_common.scores_file(gs))
    (out / "work" / "m" / "list").mkdir(parents=True)
    (out / "work" / "m" / "list" / f"{D}.txt").write_text(gs[0][0] + "\n")
    (out / "flags").mkdir()
    (out / "flags" / "release.done").write_text(gs[1][0] + "\n")
    (out / "run.lock").write_text("")
    species = Path(top) / "species.json"
    species.write_text(json.dumps({"genomes": [{"id": g[0]} for g in gs]}))
    return out, species, gs


def run(argv):
    with contextlib.redirect_stdout(io.StringIO()) as o:
        rc = pack.main([str(a) for a in argv])
    return rc, o.getvalue()


class Check(unittest.TestCase):
    def test_a_pack_of_counts_and_tables_passes_and_leaves_work_and_flags(self):
        with tempfile.TemporaryDirectory() as d:
            out, species, _ = out_dir(d)
            files = pack.check(out, pack.species_ids(species))
            self.assertEqual(files, sorted(["SHA256SUMS", "binary.txt", "logs/fit-k31-s42-m.log",
                                            f"results/inputs/program/{D}.drawn.inputs.json", "results/m/floor.tsv",
                                            f"results/picks/k31/{D}.lineage-s42.pick.json", "score/datasets.tsv"]))
            self.assertEqual(run(["check", out, species])[0], 0)

    def refused(self, plant, why):
        with tempfile.TemporaryDirectory() as d:
            out, species, gs = out_dir(d)
            plant(out, gs)
            with self.assertRaisesRegex(pack.Stop, why):
                pack.check(out, pack.species_ids(species))
            rc, said = run(["check", out, species])
            self.assertEqual(rc, 1)
            self.assertFalse(any(g[0] in said for g in gs), said)

    def test_a_score_file_is_refused_by_its_name_and_by_its_format(self):
        # By its name alone, even holding nothing; and by its format under another name.
        self.refused(lambda out, gs: m_common.write_json(out / "score" / f"{D}.lineage-s43.scores.json", {}),
                     "a score file")
        self.refused(lambda out, gs: m_common.write_json(out / "check" / "copy.json", m_common.scores_file(gs)),
                     "a score file")

    def test_a_study_fold_result_is_refused(self):
        self.refused(lambda out, gs: m_common.write_json(out / "results" / "study.json", m_common.fold_result(gs)),
                     "a study fit's fold result")

    def test_a_kept_report_is_refused(self):
        self.refused(lambda out, gs: m_common.write_json(out / "results" / "reports" / "k31.json", m_common.report(gs)),
                     "a report")

    def test_a_genome_id_in_a_table_or_a_log_is_refused(self):
        self.refused(lambda out, gs: m_common.write_tsv(out / "score" / "rows.tsv", ["dataset", "g"],
                                                        [{"dataset": D, "g": gs[4][0]}]), "names 1 genome ids")
        self.refused(lambda out, gs: (out / "logs" / "m_list.log").write_text(f"{gs[2][0]},{gs[3][0]}\n"),
                     "names 2 genome ids")

    def test_a_link_is_refused_to_a_file_or_a_folder(self):
        def to_file(out, gs):
            (out / "check").mkdir()
            os.symlink(out / "binary.txt", out / "check" / "b.txt")
        self.refused(lambda out, gs: os.symlink(out / "work" / "scores", out / "results" / "scores"), "a link")
        self.refused(to_file, "a link")

    def test_a_list_file_of_genomes_is_refused_even_empty(self):
        def plant(out, gs):
            (out / "results" / "m" / "list").mkdir(parents=True)
            (out / "results" / "m" / "list" / f"{D}.txt").write_text("")
        self.refused(plant, "a list file")

    def test_an_id_followed_by_an_underscore_or_a_point_is_refused(self):
        self.refused(lambda out, gs: (out / "logs" / "a.log").write_text(f"wrote {gs[5][0]}_scores\n"),
                     "names 1 genome ids")
        self.refused(lambda out, gs: (out / "logs" / "b.log").write_text(f"wrote {gs[5][0]}.fna\n"),
                     "names 1 genome ids")

    def test_an_id_after_an_underscore_or_a_point_is_refused(self):
        self.refused(lambda out, gs: (out / "logs" / "c.log").write_text(f"set_{gs[6][0]} and x.{gs[7][0]}\n"),
                     "names 2 genome ids")
        # Inside a longer run of points and digits, where a shorter match would hide it.
        self.refused(lambda out, gs: (out / "logs" / "d.log").write_text(f"0.{gs[8][0]}\n"), "names 1 genome ids")

    def test_an_id_in_a_file_or_folder_name_is_refused_and_not_printed(self):
        self.refused(lambda out, gs: (out / "score" / f"{gs[9][0]}.tsv").write_text("dataset\n"),
                     r"score/\[genome id\]\.tsv: its path names 1 genome ids")
        def folder(out, gs):
            (out / "check" / f"run_{gs[10][0]}").mkdir(parents=True)
            (out / "check" / f"run_{gs[10][0]}" / "x.txt").write_text("0.75\n")
        self.refused(folder, r"check/run_\[genome id\]/x\.txt: its path names 1 genome ids")

    def test_a_compressed_file_is_refused_by_its_bytes_or_its_name(self):
        import gzip
        self.refused(lambda out, gs: (out / "logs" / "fit.log").write_bytes(gzip.compress(gs[11][0].encode())),
                     r"a compressed file \(gzip\)")
        self.refused(lambda out, gs: (out / "score" / "rows.tsv.xz").write_bytes(b"dataset\n"),
                     r"a compressed file \(by its name\)")

    def test_a_file_that_is_not_utf_8_or_holds_a_nul_is_refused(self):
        self.refused(lambda out, gs: (out / "logs" / "w.log").write_bytes(f"\ufeff{gs[12][0]}".encode("utf-16")),
                     "not UTF-8 text")
        self.refused(lambda out, gs: (out / "logs" / "n.log").write_bytes(gs[12][0].encode("utf-16-le")),
                     "holds a NUL byte")
        self.refused(lambda out, gs: (out / "logs" / "l.log").write_bytes(b"read \xe9t\xe9\n"), "not UTF-8 text")

    def test_figures_and_versions_of_the_form_of_an_id_that_name_no_genome_pass(self):
        with tempfile.TemporaryDirectory() as d:
            out, species, gs = out_dir(d)
            self.assertNotIn("90007.31", {g[0] for g in gs})
            (out / "logs" / "other.log").write_text(f"0.75 1.98 90007.31 9{gs[0][0]} 90007.1000 v{gs[1][0]}\n")
            self.assertIn("logs/other.log", pack.check(out, pack.species_ids(species)))


class Tar(unittest.TestCase):
    def test_the_tar_holds_the_listed_files_alone_with_no_owner_or_time(self):
        with tempfile.TemporaryDirectory() as d:
            out, species, _ = out_dir(d)
            files = pack.check(out, pack.species_ids(species))
            rc, said = run(["tar", out, species, out / "pack.tar"])
            self.assertEqual(rc, 0)
            data = (out / "pack.tar").read_bytes()
            self.assertIn(f"sha256 {hashlib.sha256(data).hexdigest()}", said)
            with tarfile.open(out / "pack.tar") as tar:
                members = tar.getmembers()
                self.assertEqual([m.name for m in members], files)
                for m in members:
                    self.assertEqual((m.mtime, m.uid, m.gid, m.uname, m.gname, m.mode), (0, 0, 0, "", "", 0o644))
                    self.assertEqual(tar.extractfile(m).read(), (out / m.name).read_bytes())
            self.assertEqual(pack.check(out, pack.species_ids(species)), files)
            rc, again = run(["tar", out, species, Path(d) / "again.tar"])
            self.assertEqual((rc, (Path(d) / "again.tar").read_bytes()), (0, data))

    def test_a_tar_path_inside_out_other_than_pack_tar_is_refused(self):
        with tempfile.TemporaryDirectory() as d:
            out, species, _ = out_dir(d)
            rc, _ = run(["tar", out, species, out / "results" / "pack.tar"])
            self.assertEqual(rc, 2)
            self.assertFalse((out / "results" / "pack.tar").exists())

    def test_a_refused_pack_writes_no_tar(self):
        with tempfile.TemporaryDirectory() as d:
            out, species, gs = out_dir(d)
            (out / "score" / "ids.tsv").write_text(gs[7][0] + "\n")
            self.assertEqual(run(["tar", out, species, out / "pack.tar"])[0], 1)
            self.assertFalse((out / "pack.tar").exists())


if __name__ == "__main__":
    unittest.main()
