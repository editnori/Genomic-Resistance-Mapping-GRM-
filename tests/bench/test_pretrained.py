"""Synthetic unit tests for pretrained mode: deduplicated FASTA execution,
all 7 tool output parsers, genome_calls.jsonl, and projection onto rev14,
grm-fold-file-v1 (study-whole/time), and external deduplicated cohorts."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from scripts.bench.pretrained import (
    GenomeTask,
    combine_proxy_calls,
    parse_amrfinderplus_tsv,
    parse_card_rgi_tsv,
    parse_kleborate_tsv,
    parse_mykrobe_json,
    parse_resfinder_json,
    parse_staramr_summary_tsv,
    parse_tbprofiler_json,
    project_onto_external_cohort,
    project_onto_fold_file,
    project_onto_rev14_dataset,
    read_genome_calls_jsonl,
    run_pretrained_unique_genomes,
    write_genome_calls_jsonl,
)
from scripts.bench.protocol import (
    BenchError,
    EXTERNAL_COHORT_FORMAT,
    TIER_PROXY_CLASS,
    TIER_SUPPORTED_PRIMARY,
    TIER_SUPPORTED_UNGRADED,
    TIER_UNSUPPORTED,
    load_adapter,
    validate_external_cohort,
    validate_genome_call_record,
)
from tests.bench.common import build_synthetic_export, build_synthetic_fold_file


class PretrainedParsersTests(unittest.TestCase):
    def test_amrfinderplus_parser_readings(self):
        tsv = (
            "Element symbol\tElement type\tElement subtype\tScope\tMethod\tClass\tSubclass\n"
            "blaKPC-2\tAMR\tAMR\tcore\tEXACTX\tBETA-LACTAM\tCARBAPENEM\n"
            "gyrA_S83L\tAMR\tPOINT\tcore\tPOINTX\tQUINOLONE\tQUINOLONE\n"
            "tetA\tAMR\tAMR\tcore\tPARTIALX\tTETRACYCLINE\tTETRACYCLINE\n"
            "effluxX\tAMR\tAMR\tplus\tBLASTX\tMACROLIDE\tMACROLIDE\n"
        )
        drugs = ["CARBAPENEM", "QUINOLONE", "TETRACYCLINE", "MACROLIDE", "COLISTIN"]
        by_def = {r["tool_drug"]: r["call"] for r in parse_amrfinderplus_tsv(tsv, "g1", "escherichia_coli", drugs, "tool_default")}
        self.assertEqual(by_def, {"CARBAPENEM": "R", "QUINOLONE": "R", "TETRACYCLINE": "R", "MACROLIDE": "S", "COLISTIN": "S"})

        by_intact = {r["tool_drug"]: r["call"] for r in parse_amrfinderplus_tsv(tsv, "g1", "escherichia_coli", drugs, "tool_default_intact")}
        self.assertEqual(by_intact["TETRACYCLINE"], "S")
        self.assertEqual(by_intact["CARBAPENEM"], "R")

        by_gcore = {r["tool_drug"]: r["call"] for r in parse_amrfinderplus_tsv(tsv, "g1", "escherichia_coli", drugs, "genes_core")}
        self.assertEqual(by_gcore["QUINOLONE"], "S")
        self.assertEqual(by_gcore["CARBAPENEM"], "R")

        by_all = {r["tool_drug"]: r["call"] for r in parse_amrfinderplus_tsv(tsv, "g1", "escherichia_coli", drugs, "genes_all_scopes")}
        self.assertEqual(by_all["MACROLIDE"], "R")
        self.assertEqual(by_all["QUINOLONE"], "S")

    def test_resfinder_parser_and_missing_drug_no_call(self):
        payload = json.dumps({
            "phenotypes": {
                "ciprofloxacin": {"amr_resistant": True, "seq_regions": [], "point_mutations": ["gyrA_S83L"]},
                "ampicillin": {"amr_resistant": True, "seq_regions": ["blaTEM-1B"]},
                "colistin": {"amr_resistant": False, "seq_regions": []},
            }
        })
        recs = {r["tool_drug"]: r for r in parse_resfinder_json(payload, "g1", "escherichia_coli", ["ciprofloxacin", "ampicillin", "colistin", "meropenem"], "pheno_default")}
        self.assertEqual(recs["ciprofloxacin"]["call"], "R")
        self.assertEqual(recs["ampicillin"]["call"], "R")
        self.assertEqual(recs["colistin"]["call"], "S")
        self.assertEqual(recs["meropenem"]["call"], "no_call")

        acq = {r["tool_drug"]: r["call"] for r in parse_resfinder_json(payload, "g1", "escherichia_coli", ["ciprofloxacin", "ampicillin"], "acquired_only")}
        self.assertEqual(acq["ciprofloxacin"], "S")
        self.assertEqual(acq["ampicillin"], "R")

    def test_card_rgi_parser_readings(self):
        tsv = (
            "ORF_ID\tCut_Off\tBest_Hit_ARO\tDrug Class\tResistance Mechanism\tAntibiotic\n"
            "orf1\tPerfect\tNDM-1\tcarbapenem\tantibiotic inactivation\tmeropenem\n"
            "orf2\tStrict\tAcrB\tfluoroquinolone antibiotic\tantibiotic efflux\tciprofloxacin\n"
            "orf3\tLoose\tmcr-1\tpeptide antibiotic\ttarget alteration\tcolistin\n"
        )
        drugs = ["meropenem", "ciprofloxacin", "colistin"]
        sp = {r["tool_drug"]: r["call"] for r in parse_card_rgi_tsv(tsv, "g1", "escherichia_coli", drugs, "strict_perfect")}
        self.assertEqual(sp, {"meropenem": "R", "ciprofloxacin": "R", "colistin": "S"})
        ne = {r["tool_drug"]: r["call"] for r in parse_card_rgi_tsv(tsv, "g1", "escherichia_coli", drugs, "strict_perfect_no_efflux")}
        self.assertEqual(ne, {"meropenem": "R", "ciprofloxacin": "S", "colistin": "S"})
        po = {r["tool_drug"]: r["call"] for r in parse_card_rgi_tsv(tsv, "g1", "escherichia_coli", drugs, "perfect_only")}
        self.assertEqual(po, {"meropenem": "R", "ciprofloxacin": "S", "colistin": "S"})

    def test_staramr_parser_and_qc_fail_no_call(self):
        tsv = "Isolate ID\tQuality Module\tPredicted Phenotype\ng1\tFailed\tampicillin, tetracycline\n"
        recs = {r["tool_drug"]: r for r in parse_staramr_summary_tsv(tsv, "g1", "salmonella_enterica", ["ampicillin", "ciprofloxacin"])}
        self.assertEqual(recs["ampicillin"]["call"], "R")
        self.assertEqual(recs["ciprofloxacin"]["call"], "no_call")

    def test_kleborate_parser(self):
        tsv = "strain\tspecies\tBla_Carb_acquired\tFlq_mutations\tCol_mutations\ng1\tKlebsiella pneumoniae\tKPC-2\tGyrA-83I\t-\n"
        recs = {r["tool_drug"]: r["call"] for r in parse_kleborate_tsv(tsv, "g1", "klebsiella_pneumoniae", ["Bla_Carb_acquired", "Flq_mutations", "Col_mutations", "Fcyn_acquired"], "tool_default")}
        self.assertEqual(recs, {"Bla_Carb_acquired": "R", "Flq_mutations": "R", "Col_mutations": "S", "Fcyn_acquired": "no_call"})

    def test_tbprofiler_parser_grades_and_coverage_precedence(self):
        payload = json.dumps({
            "dr_variants": [
                {"type": "missense_variant", "indel_len": 0, "drugs": [{"drug": "rifampicin", "confidence": "Assoc w R"}]},
                {"type": "large_deletion", "indel_len": 120, "drugs": [{"drug": "isoniazid", "confidence": "Assoc w R - Interim"}]},
                {"type": "missense_variant", "indel_len": 0, "drugs": [{"drug": "ethambutol", "confidence": "Uncertain significance"}]},
            ],
            "qc": {
                "target_qc": [
                    {"percent_depth_pass": 85.0, "drugs": ["rifampicin", "pyrazinamide"]},
                ]
            },
        })
        drugs = ["rifampicin", "isoniazid", "ethambutol", "pyrazinamide", "amikacin"]
        wa = {r["tool_drug"]: r["call"] for r in parse_tbprofiler_json(payload, "g1", "mycobacterium_tuberculosis", drugs, "who_assoc")}
        # Resistant rifampicin stands despite incomplete coverage; pyrazinamide becomes no_call!
        self.assertEqual(wa, {"rifampicin": "R", "isoniazid": "R", "ethambutol": "S", "pyrazinamide": "no_call", "amikacin": "S"})

        wao = {r["tool_drug"]: r["call"] for r in parse_tbprofiler_json(payload, "g1", "mycobacterium_tuberculosis", drugs, "who_assoc_only")}
        self.assertEqual(wao["isoniazid"], "S")

        td = {r["tool_drug"]: r["call"] for r in parse_tbprofiler_json(payload, "g1", "mycobacterium_tuberculosis", drugs, "tool_default")}
        self.assertEqual(td["ethambutol"], "R")

        sm = {r["tool_drug"]: r["call"] for r in parse_tbprofiler_json(payload, "g1", "mycobacterium_tuberculosis", drugs, "who_assoc_small")}
        self.assertEqual(sm["isoniazid"], "S")

    def test_mykrobe_parser_never_maps_minor_r_or_uncalled_to_s(self):
        payload = json.dumps({
            "g1": {
                "susceptibility": {
                    "Rifampicin": {"predict": "R"},
                    "Isoniazid": {"predict": "r"},
                    "Ethambutol": {"predict": "N"},
                    "Pyrazinamide": {"predict": "S"},
                }
            }
        })
        recs = {r["tool_drug"]: r for r in parse_mykrobe_json(payload, "g1", "mycobacterium_tuberculosis", ["Rifampicin", "Isoniazid", "Ethambutol", "Pyrazinamide", "Amikacin"])}
        self.assertEqual(recs["Rifampicin"]["call"], "R")
        self.assertEqual(recs["Isoniazid"]["call"], "I")
        self.assertEqual(recs["Ethambutol"]["call"], "no_call")
        self.assertEqual(recs["Pyrazinamide"]["call"], "S")
        self.assertEqual(recs["Amikacin"]["call"], "no_call")


class PretrainedExecutionAndProjectionTests(unittest.TestCase):
    def test_proxy_call_combination_precedence(self):
        self.assertEqual(combine_proxy_calls([("R", None), ("no_call", "cov")])[0], "R")
        self.assertEqual(combine_proxy_calls([("R", None), ("failure", "err")])[0], "R")
        self.assertEqual(combine_proxy_calls([("failure", "err"), ("S", None)])[0], "failure")
        self.assertEqual(combine_proxy_calls([("no_call", "cov"), ("S", None)])[0], "no_call")
        self.assertEqual(combine_proxy_calls([("I", "minor"), ("S", None)])[0], "I")
        self.assertEqual(combine_proxy_calls([("S", None), ("S", None)])[0], "S")

    def test_deduplicated_unique_fasta_execution_and_failure_accounting(self):
        adapter = load_adapter("tbprofiler")
        with tempfile.TemporaryDirectory() as tmp:
            tdir = Path(tmp)
            f1 = tdir / "g1.fa"
            f2 = tdir / "g2.fa"
            f1.write_text(">c1\nATCGATCG\n", encoding="utf-8")
            f2.write_text(">c1\nGCTAGCTA\n", encoding="utf-8")

            # Same genome listed 3 times (as if across 3 drugs) -> runs ONCE!
            tasks = [
                GenomeTask("g1", "mycobacterium_tuberculosis", f1),
                GenomeTask("g1", "mycobacterium_tuberculosis", f1),
                GenomeTask("g1", "mycobacterium_tuberculosis", f1),
                GenomeTask("g2", "mycobacterium_tuberculosis", f2),
            ]
            seen_runs: list[str] = []

            def fake_runner(spec, task):
                seen_runs.append(task.genome_id)
                if task.genome_id == "g2":
                    raise RuntimeError("simulated tool crash")
                return [
                    {"genome_id": "g1", "species_slug": "mycobacterium_tuberculosis", "reading": "who_assoc", "tool_drug": "rifampicin", "call": "R"},
                    {"genome_id": "g1", "species_slug": "mycobacterium_tuberculosis", "reading": "who_assoc", "tool_drug": "isoniazid", "call": "S"},
                ]

            out_jsonl = tdir / "genome_calls.jsonl"
            stats = run_pretrained_unique_genomes(adapter, tasks, fake_runner, out_jsonl)
            self.assertEqual(seen_runs, ["g1", "g2"])
            self.assertEqual(stats["input_tasks"], 4)
            self.assertEqual(stats["unique_genomes_run"], 2)
            self.assertEqual(stats["failed_genomes"], 1)

            recs = read_genome_calls_jsonl(out_jsonl)
            g2_calls = [r for r in recs if r["genome_id"] == "g2"]
            self.assertGreater(len(g2_calls), 0)
            self.assertTrue(all(r["call"] == "failure" and "RuntimeError" in r["reason"] for r in g2_calls))

    def test_genome_calls_validation_requires_reason_on_non_binary_calls(self):
        with self.assertRaises(BenchError):
            validate_genome_call_record({
                "genome_id": "g1",
                "species_slug": "escherichia_coli",
                "reading": "tool_default",
                "tool_drug": "CIPROFLOXACIN",
                "call": "no_call",
            })
        with self.assertRaises(BenchError):
            validate_genome_call_record({
                "genome_id": "g1",
                "species_slug": "escherichia_coli",
                "reading": "tool_default",
                "tool_drug": "CIPROFLOXACIN",
                "call": "I",
            })

    def test_projection_onto_rev14_fold_file_and_external_cohort(self):
        adapter = load_adapter("tbprofiler")
        with tempfile.TemporaryDirectory() as tmp:
            tdir = Path(tmp)
            # 1. rev14 export on supported_primary (rifampin), supported_ungraded (cycloserine),
            #    proxy_class (ciprofloxacin), and unsupported (amoxicillin)
            _, meta_rif, _ = build_synthetic_export(tdir, species_slug="mycobacterium_tuberculosis", antibiotic="rifampin", n_samples=60)
            _, meta_cyc, _ = build_synthetic_export(tdir, species_slug="mycobacterium_tuberculosis", antibiotic="cycloserine", n_samples=60)
            _, meta_cip, _ = build_synthetic_export(tdir, species_slug="mycobacterium_tuberculosis", antibiotic="ciprofloxacin", n_samples=60)
            _, meta_amx, _ = build_synthetic_export(tdir, species_slug="mycobacterium_tuberculosis", antibiotic="amoxicillin", n_samples=60)

            records = []
            for idx, gid in enumerate(meta_rif["samples"]):
                # Leave last genome missing to verify missing supported call -> failure (never S!)
                if idx == 59:
                    continue
                rif_call = "no_call" if idx == 0 else ("R" if meta_rif["resistant"][idx] else "S")
                rec_rif = {"genome_id": gid, "species_slug": "mycobacterium_tuberculosis", "reading": "who_assoc", "tool_drug": "rifampicin", "call": rif_call}
                if rif_call == "no_call":
                    rec_rif["reason"] = "incomplete_target_coverage"
                records.append(rec_rif)
                records.append({"genome_id": gid, "species_slug": "mycobacterium_tuberculosis", "reading": "who_assoc", "tool_drug": "cycloserine", "call": "S"})
                records.append({"genome_id": gid, "species_slug": "mycobacterium_tuberculosis", "reading": "who_assoc", "tool_drug": "levofloxacin", "call": "R" if idx % 3 == 0 else "S"})
                records.append({"genome_id": gid, "species_slug": "mycobacterium_tuberculosis", "reading": "who_assoc", "tool_drug": "moxifloxacin", "call": "no_call" if idx == 1 else "S", **({"reason": "cov"} if idx == 1 else {})})

            jsonl_path = tdir / "calls.jsonl"
            write_genome_calls_jsonl(jsonl_path, records)
            loaded_recs = read_genome_calls_jsonl(jsonl_path)

            ep_rif = project_onto_rev14_dataset(adapter, loaded_recs, meta_rif, "lineage")
            self.assertEqual(ep_rif.support_tier, TIER_SUPPORTED_PRIMARY)
            self.assertEqual(ep_rif.calls[0], "no_call")
            self.assertEqual(ep_rif.calls[59], "failure")

            ep_cyc = project_onto_rev14_dataset(adapter, loaded_recs, meta_cyc, "lineage")
            self.assertEqual(ep_cyc.support_tier, TIER_SUPPORTED_UNGRADED)

            ep_cip = project_onto_rev14_dataset(adapter, loaded_recs, meta_cip, "lineage")
            self.assertEqual(ep_cip.support_tier, TIER_PROXY_CLASS)
            # idx=1 has levofloxacin=S and moxifloxacin=no_call -> proxy union yields no_call (NEVER S!)
            self.assertEqual(ep_cip.calls[1], "no_call")
            # idx=0 has levofloxacin=R and moxifloxacin=S -> proxy union yields R
            self.assertEqual(ep_cip.calls[0], "R")

            ep_amx = project_onto_rev14_dataset(adapter, loaded_recs, meta_amx, "lineage")
            self.assertEqual(ep_amx.support_tier, TIER_UNSUPPORTED)
            self.assertEqual(ep_amx.unsupported_reason, "drug_not_in_tbdb")
            self.assertTrue(all(c == "unsupported" for c in ep_amx.calls))

            # 2. grm-fold-file-v1 (study-whole and time)
            ff_study = build_synthetic_fold_file(tdir / "study.json", meta_rif, scheme="study-whole", n_test_folds=2, n_unused=4)
            ep_study = project_onto_fold_file(adapter, loaded_recs, meta_rif, ff_study)
            self.assertEqual(ep_study.split_scheme, "study-whole")
            self.assertEqual(len(ep_study.genome_ids), 56)
            self.assertIsNotNone(ep_study.fold_file_sha256)

            ff_time = build_synthetic_fold_file(tdir / "time.json", meta_rif, scheme="time", n_test_folds=1, n_unused=4)
            ep_time = project_onto_fold_file(adapter, loaded_recs, meta_rif, ff_time)
            self.assertEqual(ep_time.split_scheme, "time")
            self.assertLess(len(ep_time.genome_ids), 56)

            # 3. External deduplicated cohort (evaluation_only)
            ext_cohort = {
                "format": EXTERNAL_COHORT_FORMAT,
                "cohort_name": "synth_ext_mtb",
                "dataset": "mycobacterium_tuberculosis__rifampin",
                "species_slug": "mycobacterium_tuberculosis",
                "antibiotic": "rifampin",
                "selection_role": "evaluation_only",
                "deduplicated": True,
                "samples": meta_rif["samples"][:20],
                "resistant": meta_rif["resistant"][:20],
                "lineage": meta_rif["lineage"][:20],
            }
            ep_ext = project_onto_external_cohort(adapter, loaded_recs, ext_cohort)
            self.assertEqual(ep_ext.selection_role, "evaluation_only")
            self.assertEqual(ep_ext.split_scheme, "external:synth_ext_mtb")

            # External cohort guards: reject selection, non-deduplicated, or duplicate genome IDs
            with self.assertRaises(BenchError):
                validate_external_cohort(ext_cohort, allow_selection=True)
            bad_role = dict(ext_cohort, selection_role="internal_cv")
            with self.assertRaises(BenchError):
                project_onto_external_cohort(adapter, loaded_recs, bad_role)
            bad_dup = dict(ext_cohort, samples=ext_cohort["samples"][:19] + [ext_cohort["samples"][0]])
            with self.assertRaises(BenchError):
                project_onto_external_cohort(adapter, loaded_recs, bad_dup)


if __name__ == "__main__":
    unittest.main()
