"""Unit tests for src/drugs.rs vocabulary parsing and all 11 adapter registries."""

from __future__ import annotations

import csv
import json
import shutil
import tempfile
import unittest
from pathlib import Path

from scripts.bench.protocol import (
    ADAPTERS_DIR,
    ALL_TOOLS,
    BenchError,
    COMPARISON_END_TO_END,
    COMPARISON_FIXED_LEARNER,
    MODE_PRETRAINED,
    MODE_TRAINABLE,
    PRETRAINED_TOOLS,
    TIER_PROXY_CLASS,
    TIER_SUPPORTED_PRIMARY,
    TIER_SUPPORTED_UNGRADED,
    TIER_UNSUPPORTED,
    TRAINABLE_TOOLS,
    load_adapter,
    load_all_adapters,
)
from scripts.bench.vocab import DRUGS_RS, REPO_ROOT, dataset_id, drug_slug, load_vocabulary


class VocabularyTests(unittest.TestCase):
    def test_parses_src_drugs_rs_vocabulary_and_inhibitors(self):
        v = load_vocabulary(DRUGS_RS)
        self.assertGreaterEqual(len(v.drugs), 130)
        self.assertIn("meropenem", v.drugs)
        self.assertIn("rifampin", v.drugs)
        self.assertIn("para-aminosalicylic acid", v.drugs)
        self.assertIn("clavulanate", v.inhibitors)
        self.assertIn("avibactam", v.inhibitors)
        self.assertEqual(v.overrides["rifampicin"], "rifampin")

    def test_canonical_rules_match_src_drugs_rs(self):
        v = load_vocabulary(DRUGS_RS)
        # Valid single and combination names
        for good in (
            "meropenem",
            "amoxicillin/clavulanate",
            "ceftazidime/avibactam",
            "sulfamethoxazole/trimethoprim",
            "dalfopristin/quinupristin",
            "para-aminosalicylic acid",
            "cefpodoxime/ceftriaxone",
        ):
            self.assertTrue(v.is_canonical(good), good)

        # Invalid orderings, non-canonical overrides, uppercase, duplicates, unknown
        for bad in (
            "rifampicin",
            "trimethoprim/sulfamethoxazole",
            "clavulanate/amoxicillin",
            "quinupristin/dalfopristin",
            "Meropenem",
            "meropenem/meropenem",
            "clavulanate",
            "not_a_drug",
            "",
        ):
            self.assertFalse(v.is_canonical(bad), bad)

        self.assertEqual(v.canonicalize("Rifampicin"), "rifampin")
        self.assertEqual(v.canonicalize("trimethoprim/sulfamethoxazole"), "sulfamethoxazole/trimethoprim")
        self.assertEqual(v.canonicalize("amoxicillin/clavulanic acid"), "amoxicillin/clavulanate")
        self.assertEqual(v.class_of("meropenem"), "carbapenem")
        self.assertEqual(v.class_of("amoxicillin/clavulanate"), "penicillin + inhibitor")
        self.assertEqual(v.class_of("ceftazidime/avibactam"), "cephalosporin + inhibitor")
        self.assertEqual(v.class_of("sulfamethoxazole/trimethoprim"), "folate pathway combination")

    def test_rev14_lists_match_canonical_vocab_and_slugs(self):
        v = load_vocabulary(DRUGS_RS)
        lists_tsv = REPO_ROOT / "scripts" / "c3-inputs" / "lists.tsv"
        with lists_tsv.open(encoding="utf-8", newline="") as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        self.assertEqual(len(rows), 268)
        for r in rows:
            ab = r["antibiotic"]
            self.assertTrue(v.is_canonical(ab), f"rev14 antibiotic not canonical: {ab}")
            self.assertEqual(dataset_id(r["species_slug"], ab), r["dataset"])
            self.assertEqual(r["dataset"].split("__", 1)[1], drug_slug(ab))


class AdapterRegistryTests(unittest.TestCase):
    def test_all_11_adapters_load_and_validate_against_src_drugs_rs(self):
        adapters = load_all_adapters(ADAPTERS_DIR)
        self.assertEqual(set(adapters), set(ALL_TOOLS))
        self.assertEqual(len(adapters), 11)

        for tool in PRETRAINED_TOOLS:
            spec = adapters[tool]
            self.assertEqual(spec.mode, MODE_PRETRAINED)
            self.assertEqual(spec.comparison_kind, COMPARISON_END_TO_END)
            self.assertIsNotNone(spec.manifest["database_pin"])
            self.assertEqual(len(spec.species_support), 274, tool)

        for tool in TRAINABLE_TOOLS:
            spec = adapters[tool]
            self.assertEqual(spec.mode, MODE_TRAINABLE)
            self.assertEqual(len(spec.species_support), 274, tool)
            if tool == "kover":
                self.assertEqual(spec.comparison_kind, COMPARISON_END_TO_END)
            else:
                self.assertEqual(spec.comparison_kind, COMPARISON_FIXED_LEARNER)

    def test_tbprofiler_matches_committed_rev14_tier_structure(self):
        spec = load_adapter("tbprofiler")
        self.assertEqual(spec.bootstrap_cluster_unit, "lineage")
        mtb = {ab: e for (sp, ab), e in spec.species_support.items() if sp == "mycobacterium_tuberculosis"}
        self.assertEqual(len(mtb), 22)
        prim = sorted(ab for ab, e in mtb.items() if e.support_tier == TIER_SUPPORTED_PRIMARY)
        ungr = sorted(ab for ab, e in mtb.items() if e.support_tier == TIER_SUPPORTED_UNGRADED)
        prox = sorted(ab for ab, e in mtb.items() if e.support_tier == TIER_PROXY_CLASS)
        unsup = sorted(ab for ab, e in mtb.items() if e.support_tier == TIER_UNSUPPORTED)
        self.assertEqual(len(prim), 14)
        self.assertEqual(ungr, ["cycloserine", "para-aminosalicylic acid"])
        self.assertEqual(prox, ["ciprofloxacin", "ofloxacin", "prothionamide", "rifabutin"])
        self.assertEqual(unsup, ["amoxicillin", "nicotinamide"])
        # All non-Mtb species must be unsupported with explicit reason
        non_mtb = [e for (sp, _), e in spec.species_support.items() if sp != "mycobacterium_tuberculosis"]
        self.assertEqual(len(non_mtb), 252)
        self.assertTrue(all(e.support_tier == TIER_UNSUPPORTED and e.unsupported_reason == "species_not_mtb" for e in non_mtb))

    def test_amrfinderplus_matches_253_published_rule_endpoints(self):
        spec = load_adapter("amrfinderplus")
        self.assertEqual(spec.bootstrap_cluster_unit, "species")
        prim = [e for e in spec.species_support.values() if e.support_tier == TIER_SUPPORTED_PRIMARY]
        prox = [e for e in spec.species_support.values() if e.support_tier == TIER_PROXY_CLASS]
        unsup = [e for e in spec.species_support.values() if e.support_tier == TIER_UNSUPPORTED]
        self.assertEqual(len(prim), 241)
        self.assertEqual(len(prox), 12)
        self.assertEqual(len(unsup), 21)

    def test_mykrobe_declares_explicit_intermediate_handling_and_lineage_bootstrap(self):
        spec = load_adapter("mykrobe")
        self.assertEqual(spec.bootstrap_cluster_unit, "lineage")
        self.assertTrue(spec.manifest["i_handling"]["tool_emits_intermediate"])
        self.assertEqual(spec.manifest["i_handling"]["policy"], "separate_count_never_map_to_s")
        # S. enterica non-typhi must be explicitly unsupported
        salm = [e for (sp, _), e in spec.species_support.items() if sp == "salmonella_enterica"]
        self.assertTrue(all(e.support_tier == TIER_UNSUPPORTED and e.unsupported_reason == "mykrobe_panel_typhi_only" for e in salm))

    def test_unknown_species_drug_pair_defaults_to_unsupported_never_susceptible(self):
        spec = load_adapter("kleborate")
        entry = spec.support_for("vibrio_cholerae", "meropenem")
        self.assertEqual(entry.support_tier, TIER_UNSUPPORTED)
        self.assertEqual(entry.unsupported_reason, "not_in_precommitted_support_registry")

    def test_corrupted_adapter_registries_are_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "adapters"
            shutil.copytree(ADAPTERS_DIR / "tbprofiler", root / "tbprofiler")

            # 1. Invalid canonical drug in drug_map.tsv
            dm_path = root / "tbprofiler" / "drug_map.tsv"
            orig_dm = dm_path.read_text(encoding="utf-8")
            dm_path.write_text(orig_dm + "fake_drug\trifampicin\texact\tbad spelling\n", encoding="utf-8")
            with self.assertRaises(BenchError):
                load_adapter("tbprofiler", adapters_dir=root)
            dm_path.write_text(orig_dm, encoding="utf-8")

            # 2. Missing unsupported_reason on unsupported row
            sp_path = root / "tbprofiler" / "species_support.tsv"
            orig_sp = sp_path.read_text(encoding="utf-8")
            sp_path.write_text(orig_sp.replace("drug_not_in_tbdb", "-"), encoding="utf-8")
            with self.assertRaises(BenchError):
                load_adapter("tbprofiler", adapters_dir=root)
            sp_path.write_text(orig_sp, encoding="utf-8")

            # 3. Bad i_handling policy in manifest.json
            m_path = root / "tbprofiler" / "manifest.json"
            m_doc = json.loads(m_path.read_text(encoding="utf-8"))
            m_doc["i_handling"]["policy"] = "map_i_to_s"
            m_path.write_text(json.dumps(m_doc), encoding="utf-8")
            with self.assertRaises(BenchError):
                load_adapter("tbprofiler", adapters_dir=root)


if __name__ == "__main__":
    unittest.main()
