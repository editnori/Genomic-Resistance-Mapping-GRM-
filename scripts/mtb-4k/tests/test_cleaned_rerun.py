#!/usr/bin/env python3
"""cleaned_rerun.py: the cleaned cohorts' placeholders refuse until filled, are checked by kind, and
split into the prepare group (sections 3-8) and the reference group (section 7b, pending).

usage: python3 -I scripts/mtb-4k/tests/test_cleaned_rerun.py
"""
import json,re,subprocess,sys,tempfile,unittest
from pathlib import Path
LANE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(LANE))
import cleaned_rerun as cr

H="a"*64
FILLED={"SPECIES_CHECK_PATH":"/data/species-check/species_check.json","SPECIES_CHECK_SHA256":H,
 "SPECIES_CHECK_FORMAT":cr.SPECIES_CHECK_FORMAT_PARSED,"EXCLUSIONS_COUNT":2,
 "RN_PREFIX_GENOMES_PATH":"/r9/data/species/mycobacterium_tuberculosis/genomes.txt","RN_PREFIX_GENOMES_SHA256":H,
 "RN_PREFIX_SPECIES_JSON_PATH":"/r9/data/curated/species.json","RN_PREFIX_SPECIES_JSON_SHA256":H,
 "RN_FULL_GENOMES_PATH":"/r9-4k/data/species/mycobacterium_tuberculosis/genomes.txt","RN_FULL_GENOMES_SHA256":H,
 "RN_FULL_SPECIES_JSON_PATH":"/r9-4k/data/curated/species.json","RN_FULL_SPECIES_JSON_SHA256":H,
 "RN_INPUT_ROOT":"/r9/data","RN_FASTA_MANIFEST_PATH":"/r9/fasta-manifest.tsv","RN_FASTA_MANIFEST_SHA256":H,
 "RN_FASTA_STORE":"/r9/data/fasta","CLEANED_PREFIX_SHA256":H,"CLEANED_PREFIX_COUNT":2000,"CLEANED_ADDED_SHA256":H,
 "CLEANED_ADDED_COUNT":2000,"CLEANED_FULL_SHA256":H,"CLEANED_FULL_COUNT":4000,"CLEANED_MATCHED2K_SHA256":H,
 "CLEANED_MATCHED2K_COUNT":2000,"RN_REFERENCE_ROOT":"/r9"}

def run(script,*args,cwd):
 return subprocess.run([sys.executable,"-I","-B",str(script),*map(str,args)],cwd=cwd,capture_output=True,text=True)

def with_values(directory,values):
 """A copy of the lane's cleaned_rerun.py and revision_n.py whose VALUES block is rewritten."""
 for name in ("revision_n.py","cleaned_rerun.py"):
  (Path(directory)/name).write_text((LANE/name).read_text())
 path=Path(directory)/"cleaned_rerun.py"; text=path.read_text()
 start=text.index("\nVALUES = {\n")+1; end=text.index("\n}\n",start)+3
 path.write_text(text[:start]+"VALUES = {\n"+"".join(f"    {n!r}: {values[n]!r},\n" for n in cr.TABLE)+"}\n"+text[end:])
 return path

class Placeholders(unittest.TestCase):
 def test_the_table_values_and_groups_agree(self):
  self.assertEqual(set(cr.VALUES),set(cr.TABLE)); self.assertEqual(set(FILLED),set(cr.TABLE))
  self.assertEqual(sorted(cr.GROUPS["prepare"]+cr.GROUPS["reference"]),sorted(cr.TABLE))
  self.assertEqual(cr.GROUPS["reference"],["RN_REFERENCE_ROOT"])
  self.assertEqual(sorted(cr.INPUT_NAMES+cr.CONTRACT_NAMES),sorted(cr.GROUPS["prepare"]))
  self.assertEqual(cr.REQUIRED_EXCLUSIONS,("1773.41701","1773.40949"))
 def test_the_committed_file_is_wholly_unfilled_and_refuses_or_wholly_filled(self):
  missing=[n for n in cr.TABLE if cr.rev_n.is_placeholder(cr.VALUES[n])]
  self.assertIn(missing,[list(cr.TABLE),[]],"cleaned_rerun.py is partly filled")
  if missing:
   with self.assertRaises(cr.Unfilled) as stop: cr.require(list(cr.TABLE))
   self.assertEqual(stop.exception.names,list(cr.TABLE))
  else: cr.require(list(cr.TABLE))
 def test_unfilled_values_refuse_naming_every_one_and_the_cli_exits_2(self):
  unfilled={n:f"@@{n}@@" for n in cr.TABLE}
  with self.assertRaises(cr.Unfilled) as stop: cr.require(cr.INPUT_NAMES,unfilled)
  self.assertEqual(stop.exception.names,cr.INPUT_NAMES)
  with self.assertRaises(cr.Unfilled) as stop: cr.require(list(cr.TABLE),{**FILLED,"RN_REFERENCE_ROOT":"@@RN_REFERENCE_ROOT@@"})
  self.assertEqual(stop.exception.names,["RN_REFERENCE_ROOT"])
  with tempfile.TemporaryDirectory() as tmp:
   script=with_values(tmp,unfilled)
   p=run(script,"status",cwd=tmp)
   self.assertEqual(p.returncode,2); self.assertIn("cleaned_rerun: STOP: unfilled placeholder(s)",p.stderr)
   self.assertEqual(json.loads(p.stdout),{n:False for n in cr.TABLE})
   p=run(script,"get","RN_INPUT_ROOT",cwd=tmp); self.assertEqual((p.returncode,p.stdout),(2,"")); self.assertIn("@@RN_INPUT_ROOT@@",p.stderr)
 def test_the_groups_let_training_proceed_while_the_reference_is_pending(self):
  with tempfile.TemporaryDirectory() as tmp:
   script=with_values(tmp,{**FILLED,"RN_REFERENCE_ROOT":"@@RN_REFERENCE_ROOT@@"})
   self.assertEqual(run(script,"status","--group","prepare",cwd=tmp).returncode,0)
   p=run(script,"status","--group","reference",cwd=tmp); self.assertEqual(p.returncode,2); self.assertIn("@@RN_REFERENCE_ROOT@@",p.stderr)
   self.assertEqual(run(script,"status",cwd=tmp).returncode,2)
   p=run(script,"get","CLEANED_FULL_COUNT",cwd=tmp); self.assertEqual((p.returncode,p.stdout.strip()),(0,"4000"))
   script=with_values(tmp,FILLED)
   self.assertEqual(run(script,"status",cwd=tmp).returncode,0)

class Kinds(unittest.TestCase):
 def test_filled_values_are_well_formed(self):
  self.assertEqual(cr.require(list(cr.TABLE),FILLED),FILLED)
 def test_values_are_checked_by_kind(self):
  for name,bad in [("SPECIES_CHECK_PATH","relative/path.json"),("SPECIES_CHECK_PATH","/a/../b"),("SPECIES_CHECK_PATH","/a/"),
                   ("SPECIES_CHECK_SHA256","A"*64),("EXCLUSIONS_COUNT",1),("EXCLUSIONS_COUNT",True),("CLEANED_FULL_COUNT",0),
                   ("CLEANED_PREFIX_COUNT","2000"),("RN_REFERENCE_ROOT",3)]:
   with self.subTest(name=name,bad=bad):
    with self.assertRaises(cr.Malformed): cr.require([name],{**FILLED,name:bad})
 def test_another_species_check_format_needs_a_parser_change(self):
  with self.assertRaises(cr.Malformed) as stop: cr.require(["SPECIES_CHECK_FORMAT"],{**FILLED,"SPECIES_CHECK_FORMAT":"species-check-v2"})
  self.assertIn("reviewed parser change",str(stop.exception))

if __name__=="__main__":
 result=unittest.main(exit=False,verbosity=2).result
 sys.exit(0 if result.wasSuccessful() else 1)
