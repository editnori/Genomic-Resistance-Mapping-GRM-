"""Synthetic unit tests for external AMR portal metadata normalizers."""

from __future__ import annotations

import gzip
import importlib.util
import json
import pathlib
import subprocess
import sys
import tempfile
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "dedup" / "portals.py"


def _load_portals():
    spec = importlib.util.spec_from_file_location("dedup_portals_mod", SCRIPT)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


portals = _load_portals()


class PortalNormalizersTest(unittest.TestCase):
    def test_ebi_amr_merges_multi_row_ast_and_emits_insdc_synonyms(self):
        rows = [
            {
                "BioSample_ID": "SAMEA100001",
                "secondary_sample_accession": "ERS500001",
                "run_accession": "ERR700001",
                "assembly_accession": "GCA_000111222.1",
                "bioproject_accession": "PRJEB12345",
                "secondary_study_accession": "ERP012345",
                "scientific_name": "Escherichia coli str. K-12",
                "strain": "  EC-01  ",
                "pmid": "PMID:31234567",
                "antibiotic": "ampicillin",
            },
            {
                "BioSample_ID": "SAMEA100001",
                "secondary_sample_accession": "ERS500001",
                "run_accession": "ERR700002",
                "assembly_accession": "NA",
                "bioproject_accession": "PRJEB12345",
                "secondary_study_accession": "ERP012345",
                "scientific_name": "Escherichia coli",
                "strain": "EC-01",
                "pmid": "31234567;31234568",
                "antibiotic": "ciprofloxacin",
            },
        ]
        recs, syns = portals.normalize_rows("cabbage", rows)
        self.assertEqual(len(recs), 1)
        r = recs[0]
        self.assertEqual(r["genome_id"], "SAMEA100001")
        self.assertEqual(r["portal"], "ebi_amr")
        self.assertEqual(r["species"], "Escherichia coli")
        self.assertEqual(r["biosample_accession"], "ERS500001")
        self.assertEqual(r["biosample_aliases"], ["SAMEA100001"])
        self.assertEqual(r["assembly_accession"], "GCA_000111222.1")
        self.assertEqual(r["sra_accession"], "ERR700001;ERR700002")
        self.assertEqual(r["bioproject_accession"], "ERP012345")
        self.assertEqual(r["bioproject_aliases"], ["PRJEB12345"])
        self.assertEqual(r["pmid"], "31234567;31234568")
        self.assertEqual(r["strain"], "EC-01")
        self.assertIn(("sample", "ERS500001", "SAMEA100001"), syns)
        self.assertIn(("run_to_sample", "ERR700001", "SAMEA100001"), syns)
        self.assertIn(("run_to_sample", "ERR700002", "SAMEA100001"), syns)
        self.assertIn(("study", "ERP012345", "PRJEB12345"), syns)

    def test_cryptic_sets_default_mtb_species_and_longitudinal_subject_id(self):
        rows = [
            {
                "UNIQUEID": "site.02.subj.0010.lab.1",
                "SITEID": "02",
                "SUBJID": "0010",
                "BIOSAMPLE": "SAMEA90001",
                "ENA_SAMPLE": "ERS90001",
                "ENA_RUN": "ERR90001",
                "BIOPROJECT": "PRJEB41116",
                "ENA_STUDY": "ERP124852",
            },
            {
                "UNIQUEID": "site.02.subj.0010.lab.2",
                "SITEID": "02",
                "SUBJID": "0010",
                "BIOSAMPLE": "SAMEA90002",
                "ENA_SAMPLE": "ERS90002",
                "ENA_RUN": "ERR90002",
                "BIOPROJECT": "PRJEB41116",
                "ENA_STUDY": "ERP124852",
            },
        ]
        recs, syns = portals.normalize_rows("cryptic", rows)
        self.assertEqual(len(recs), 2)
        for rec in recs:
            self.assertEqual(rec["species"], "Mycobacterium tuberculosis")
            self.assertEqual(rec["subject_id"], "02:0010")
        self.assertIn(("sample", "ERS90001", "SAMEA90001"), syns)
        self.assertIn(("run_to_sample", "ERR90002", "SAMEA90002"), syns)
        self.assertIn(("study", "ERP124852", "PRJEB41116"), syns)

    def test_ncbi_pd_narms_links_genbank_refseq_and_run_to_biosample(self):
        rows = [
            {
                "target_acc": "PDT000123456.1",
                "biosample_acc": "SAMN12345678",
                "sra_sample": "SRS1234567",
                "asm_acc": "GCA_001234567.1",
                "refseq_asm_acc": "GCF_001234567.1",
                "Run": "SRR9990001",
                "bioproject_acc": "PRJNA292661",
                "scientific_name": "Salmonella enterica subsp. enterica",
                "strain": "SAL-42",
                "host_subject_id": "Not applicable",
            }
        ]
        recs, syns = portals.normalize_rows("narms", rows)
        self.assertEqual(len(recs), 1)
        r = recs[0]
        self.assertEqual(r["genome_id"], "PDT000123456.1")
        self.assertEqual(r["portal"], "ncbi_pd")
        self.assertEqual(r["species"], "Salmonella enterica")
        self.assertEqual(r["assembly_accession"], "GCA_001234567.1;GCF_001234567.1")
        self.assertNotIn("subject_id", r)
        self.assertIn(("sample", "SAMN12345678", "SRS1234567"), syns)
        self.assertIn(("run_to_sample", "SRR9990001", "SAMN12345678"), syns)
        self.assertIn(("assembly", "GCA_001234567.1", "GCF_001234567.1"), syns)

    def test_enterobase_and_pubmlst_normalizers(self):
        eb_rows = [
            {
                "Uberstrain": "ESC_AA0001AA",
                "Name": "EC_2025_01",
                "Sample ID": "SAMN555001",
                "Secondary Sample Accession": "SRS555001",
                "Accession": "SRR555001",
                "Assembly Accession": "GCA_005550001.1",
                "Project ID": "PRJNA555001",
                "Study Accession": "SRP555001",
                "Species": "Escherichia coli",
            }
        ]
        eb_recs, eb_syns = portals.normalize_rows("enterobase", eb_rows)
        self.assertEqual(eb_recs[0]["genome_id"], "ESC_AA0001AA")
        self.assertEqual(eb_recs[0]["strain"], "EC_2025_01")
        self.assertIn(("sample", "SAMN555001", "SRS555001"), eb_syns)
        self.assertIn(("study", "PRJNA555001", "SRP555001"), eb_syns)

        pm_rows = [
            {
                "id": "4201",
                "isolate": "NG-4201",
                "species": "Neisseria gonorrhoeae",
                "biosample_accession": "SAMEA777001",
                "ena_sample": "ERS777001",
                "ena_run_accession": "ERR777001",
                "bioproject_accession": "PRJEB77701",
                "ena_study": "ERP077701",
                "patient_id": "PT-9",
                "pubmed_id": "29000001",
            }
        ]
        pm_recs, pm_syns = portals.normalize_rows("pubmlst", pm_rows)
        self.assertEqual(pm_recs[0]["genome_id"], "pubmlst:4201")
        self.assertEqual(pm_recs[0]["subject_id"], "PT-9")
        self.assertEqual(pm_recs[0]["pmid"], "29000001")
        self.assertIn(("sample", "ERS777001", "SAMEA777001"), pm_syns)
        self.assertIn(("run_to_sample", "ERR777001", "SAMEA777001"), pm_syns)
        self.assertIn(("study", "ERP077701", "PRJEB77701"), pm_syns)

    def test_determinism_under_row_reordering_and_gzip_cli_round_trip(self):
        with tempfile.TemporaryDirectory() as tmp:
            tdir = pathlib.Path(tmp)
            tsv_gz = tdir / "ebi.tsv.gz"
            header = "BioSample_ID\tsecondary_sample_accession\trun_accession\tspecies\tbioproject_accession\tsecondary_study_accession\n"
            row_a = "SAMEA2\tERS2\tERR2\tKlebsiella pneumoniae\tPRJEB2\tERP2\n"
            row_b = "SAMEA1\tERS1\tERR1\tKlebsiella pneumoniae\tPRJEB1\tERP1\n"
            with gzip.open(tsv_gz, "wt", encoding="utf-8", newline="\n") as fh:
                fh.write(header + row_a + row_b)

            tsv_rev = tdir / "ebi_rev.tsv"
            tsv_rev.write_text(header + row_b + row_a, encoding="utf-8")

            out1 = tdir / "q1.jsonl"
            syn1 = tdir / "s1.tsv"
            out2 = tdir / "q2.jsonl"
            syn2 = tdir / "s2.tsv"

            for inp, q_out, s_out in ((tsv_gz, out1, syn1), (tsv_rev, out2, syn2)):
                proc = subprocess.run(
                    [
                        sys.executable,
                        "-I",
                        "-B",
                        str(SCRIPT),
                        "--portal",
                        "ebi_amr",
                        "--input",
                        str(inp),
                        "--out",
                        str(q_out),
                        "--synonym-out",
                        str(s_out),
                    ],
                    capture_output=True,
                    text=True,
                    check=True,
                )
                self.assertEqual(proc.returncode, 0)

            self.assertEqual(out1.read_bytes(), out2.read_bytes())
            self.assertEqual(syn1.read_bytes(), syn2.read_bytes())
            lines = [json.loads(line) for line in out1.read_text(encoding="utf-8").splitlines()]
            self.assertEqual([x["genome_id"] for x in lines], ["SAMEA1", "SAMEA2"])

            merged_syn = tdir / "merged.tsv"
            proc = subprocess.run(
                [
                    sys.executable,
                    "-I",
                    "-B",
                    str(SCRIPT),
                    "--merge-synonyms",
                    str(syn1),
                    str(syn2),
                    "--synonym-out",
                    str(merged_syn),
                ],
                capture_output=True,
                text=True,
                check=True,
            )
            self.assertEqual(proc.returncode, 0)
            self.assertEqual(merged_syn.read_bytes(), syn1.read_bytes())

    def test_conflicting_subject_id_and_missing_identifiers_raise(self):
        with self.assertRaises(ValueError):
            portals.normalize_rows(
                "cryptic",
                [
                    {"UNIQUEID": "u1", "SUBJID": "A"},
                    {"UNIQUEID": "u1", "SUBJID": "B"},
                ],
            )
        with self.assertRaises(ValueError):
            portals.normalize_rows("ebi_amr", [{"species": "Escherichia coli"}])
        with self.assertRaises(ValueError):
            portals.canonical_portal("not_a_portal")


if __name__ == "__main__":
    unittest.main()
