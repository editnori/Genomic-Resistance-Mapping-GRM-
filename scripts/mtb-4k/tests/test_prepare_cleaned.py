#!/usr/bin/env python3
"""prepare_execution_roots.py --cohorts cleaned on a synthetic frozen 4k, a fresh selection and a
post-fetch species check, in the shapes relayed on 2026-10-29.

usage: python3 -I scripts/mtb-4k/tests/test_prepare_cleaned.py
"""
import contextlib,csv,gzip,hashlib,io,json,os,shutil,sys,tempfile,unittest
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent)); sys.path.insert(0,str(HERE))
import prepare_execution_roots as prep  # noqa: E402
import test_prepare_execution_roots as tp  # noqa: E402  (the frozen-cohort fixture)

MTB="Mycobacterium tuberculosis"
RULED=("1773.41701","1773.40949")
BORDERLINE="1773.100114"  # a frozen prefix genome the check calls borderline: kept, recorded
REPLACEMENTS=["1773.900001","1773.900002"]

def sha(b): return hashlib.sha256(b).hexdigest()
def lines(ids): return "".join(f"{g}\n" for g in ids).encode()

class Cleaned(unittest.TestCase):
 def setUp(self):
  self.saved=(prep.EVIDENCE,prep.CANDIDATE,prep.ALLOWED_FASTA_PARENTS,prep.EXPECTED,prep.free_bytes,dict(prep.rerun.VALUES))
  self.tmp=tempfile.TemporaryDirectory(); self.base=Path(self.tmp.name)
  ids=[f"1773.{100000+i}" for i in range(4000)]; ids[352],ids[526]=RULED  # rows 353 and 527, as in the frozen list
  self.frozen=ids
  tp.Prepare.redirect(self,self.base,ids=ids)
  prep.free_bytes=lambda _:prep.MIN_FREE_AFTER+10**9
  # The fresh selection: the frozen rows 1-2000 without the two, cap-filled by two replacements, two neighbours swapped.
  prefix=[g for g in ids[:2000] if g not in RULED]; prefix[10],prefix[11]=prefix[11],prefix[10]; prefix+=REPLACEMENTS
  self.prefix=prefix; self.full=prefix+ids[2000:]
  self.verdicts={g:"pass" for g in ids+REPLACEMENTS}; self.verdicts.update({g:"fail" for g in RULED}); self.verdicts[BORDERLINE]="borderline"
  self.unchecked=[]
  self.store=self.base/"r9-fasta"; self.store.mkdir()
  rows=[]
  for g in REPLACEMENTS:
   p=self.store/f"{g}.fna.gz"; p.write_bytes(gzip.compress(f">c\nTTGA{g}\n".encode(),mtime=0)); b=p.read_bytes()
   rows.append({"genome_id":g,"size_bytes":len(b),"sha256":sha(b),"canonical_source_path":str(p)})
  self.fasta_rows=rows
  (prep.CANDIDATE/"curated"/"phenotypes.jsonl").write_text('{"genome_id":"1773.100000","antibiotic":"isoniazid","label":"R"}\n')
  (prep.CANDIDATE/"catalog"/"AMR_CDS.fa").write_text(">blaX\nATGAAA\n")
  (prep.CANDIDATE/"catalog"/"ReferenceGeneCatalog.txt").write_text("allele\tgene_family\nblaX-1\tblaX\n")
  (prep.CANDIDATE/"curated"/"funnel.json").write_text("{}\n")
  # Every cohort genome has a usable cached annotation, as revision N's staging must provide.
  for g in ids+REPLACEMENTS: (prep.CANDIDATE/"bvbrc"/"features"/f"{g}.json").write_text(f'[{{"gene": "{g}"}}]')
  self.write_inputs()

 def tearDown(self):
  prep.EVIDENCE,prep.CANDIDATE,prep.ALLOWED_FASTA_PARENTS,prep.EXPECTED,prep.free_bytes,values=self.saved
  prep.rerun.VALUES.clear(); prep.rerun.VALUES.update(values); self.tmp.cleanup()

 def write_inputs(self,prefix=None,full=None,prefix_unchecked=0,full_unchecked=0,count=2,fmt=None):
  prefix=self.prefix if prefix is None else prefix; full=self.full if full is None else full
  d=self.base/"r9"; d.mkdir(exist_ok=True)
  check={"k":21,"sketch_size":1000,"config":{"fail_distance":0.13,"species_distance":0.05,"gc_z":6.0,"gc_points":1.0},"pairs":[],
   "species":[{"species":"Other","slug":"other","unchecked":[],"failed":0,"borderline":0},
              {"species":MTB,"slug":prep.SLUG,"reference":"1267357.3","checked":len(self.verdicts),"unchecked":self.unchecked,
               "failed":sum(v=="fail" for v in self.verdicts.values()),"borderline":sum(v=="borderline" for v in self.verdicts.values())}],
   "genomes":[{"genome_id":g,"species":MTB,"verdict":v,"gc":0.656,"distance":0.0003} for g,v in sorted(self.verdicts.items())]
             +[{"genome_id":"1773.41701","species":"Other","verdict":"pass"}]}
  files={"species_check.json":json.dumps(check).encode(),"prefix-genomes.txt":lines(prefix),"full-genomes.txt":lines(full),
   "prefix-species.json":json.dumps([{"species":"Other","slug":"other","selected_genomes":5,"selected_unchecked":3},
                                      {"species":MTB,"slug":prep.SLUG,"datasets":22,"labeled_genomes":5000,"selected_genomes":len(prefix),"selected_unchecked":prefix_unchecked}]).encode(),
   "full-species.json":json.dumps([{"species":MTB,"slug":prep.SLUG,"datasets":22,"labeled_genomes":5000,"selected_genomes":len(full),"selected_unchecked":full_unchecked}]).encode()}
  buf=io.StringIO(); w=csv.DictWriter(buf,fieldnames=prep.MANIFEST_COLUMNS,delimiter="\t",lineterminator="\n"); w.writeheader(); w.writerows(self.fasta_rows)
  files["fasta-manifest.tsv"]=buf.getvalue().encode()
  for name,data in files.items(): (d/name).write_bytes(data)
  h={name:sha(data) for name,data in files.items()}
  values={"SPECIES_CHECK_PATH":str(d/"species_check.json"),"SPECIES_CHECK_SHA256":h["species_check.json"],
   "SPECIES_CHECK_FORMAT":fmt or prep.rerun.SPECIES_CHECK_FORMAT_PARSED,"EXCLUSIONS_COUNT":count,
   "RN_PREFIX_GENOMES_PATH":str(d/"prefix-genomes.txt"),"RN_PREFIX_GENOMES_SHA256":h["prefix-genomes.txt"],
   "RN_PREFIX_SPECIES_JSON_PATH":str(d/"prefix-species.json"),"RN_PREFIX_SPECIES_JSON_SHA256":h["prefix-species.json"],
   "RN_FULL_GENOMES_PATH":str(d/"full-genomes.txt"),"RN_FULL_GENOMES_SHA256":h["full-genomes.txt"],
   "RN_FULL_SPECIES_JSON_PATH":str(d/"full-species.json"),"RN_FULL_SPECIES_JSON_SHA256":h["full-species.json"],
   "RN_INPUT_ROOT":str(prep.CANDIDATE),"RN_FASTA_MANIFEST_PATH":str(d/"fasta-manifest.tsv"),"RN_FASTA_MANIFEST_SHA256":h["fasta-manifest.tsv"],
   "RN_FASTA_STORE":str(self.store),"RN_REFERENCE_ROOT":"@@RN_REFERENCE_ROOT@@",
   **{name:sha((prep.CANDIDATE/rel).read_bytes()) for rel,name in prep.rerun.INPUT_ROOT_FILES.items()}}
  prep.rerun.VALUES.update(values)
  full_list=list(full); n=len(prefix)
  derived={"prefix":list(prefix),"added":full_list[n:],"full":full_list,"matched2k":prep.matched_selection(full_list,n)}
  prep.rerun.VALUES.update(prep.derived_pins(derived,prep.CANDIDATE))
  return derived

 def main(self,*argv):
  saved=sys.argv; sys.argv=["prepare_execution_roots.py",*map(str,argv)]; out=io.StringIO()
  try:
   with contextlib.redirect_stdout(out): return prep.main(),out.getvalue()
  finally: sys.argv=saved

 def refused(self,needle):
  run=self.base/"run"
  with self.assertRaises(SystemExit) as stop: self.main(run,"--cohorts","cleaned")
  self.assertIn(needle,str(stop.exception)); self.assertFalse(os.path.lexists(run),"a refused preparation created the run root")
  return str(stop.exception)

 def test_the_cleaned_roots_follow_the_ruled_design(self):
  rc,out=self.main("--derive-cleaned"); self.assertEqual(rc,0)
  derived=json.loads(out)
  self.assertEqual(derived["pins"],{n:prep.rerun.VALUES[n] for n in derived["pins"]})
  run=self.base/"run"; rc,_=self.main(run,"--cohorts","cleaned"); self.assertEqual(rc,0)
  receipt=json.loads((run/"preparation.json").read_text())
  self.assertEqual(receipt["cohorts"],{"prefix":2000,"added":2000,"full":4000,"matched2k":2000})
  self.assertEqual(sorted(os.listdir(run/"cohorts")),["added","full","matched2k","prefix"])  # no pilot roots
  c=receipt["cleaned"]
  self.assertEqual(sorted(c["exclusions"]),sorted(RULED)); self.assertEqual(c["replacements"],REPLACEMENTS)
  self.assertEqual(c["borderline"]["prefix"],[BORDERLINE]); self.assertEqual(c["borderline"]["full"],[BORDERLINE])
  self.assertEqual(c["prefix_vs_frozen_rows_1_2000"]["left"],list(RULED)); self.assertEqual(c["prefix_vs_frozen_rows_1_2000"]["joined"],REPLACEMENTS)
  self.assertFalse(c["prefix_vs_frozen_rows_1_2000"]["same_order_on_common"])
  for name,ids in {"prefix":self.prefix,"added":self.full[2000:],"full":self.full}.items():
   self.assertEqual((run/"cohorts"/name/"data/species"/prep.SLUG/"genomes.txt").read_bytes(),lines(ids))
   self.assertEqual(sorted(os.listdir(run/"cohorts"/name/"data/fasta")),sorted(f"{g}.fna.gz" for g in ids))
   rc,out=self.main("--check-lane",run/"cohorts"/name/"data"); self.assertEqual(rc,0,out)
  self.assertEqual(os.readlink(run/"cohorts/prefix/data/fasta"/f"{REPLACEMENTS[0]}.fna.gz"),str(self.store/f"{REPLACEMENTS[0]}.fna.gz"))
  self.assertNotIn(RULED[0],(run/"cohorts/full/data/species"/prep.SLUG/"genomes.txt").read_text().split())
  with open(run/"fasta-manifest.tsv",newline="") as f: rows=list(csv.DictReader(f,delimiter="\t"))
  self.assertEqual([r["genome_id"] for r in rows],self.full)
  self.assertEqual(self.main("--verify-frozen",run)[0],0)
  # A bound input changed afterwards is a stop.
  sc=Path(prep.rerun.VALUES["SPECIES_CHECK_PATH"]); sc.write_bytes(sc.read_bytes()+b" ")
  rc,out=self.main("--verify-frozen",run); self.assertEqual(rc,1); self.assertIn(str(sc),json.loads(out)["changed"])

 def test_a_borderline_genome_is_never_excluded(self):
  self.verdicts.update({g:"borderline" for g in self.prefix[:5]}); self.write_inputs()
  run=self.base/"run"; self.assertEqual(self.main(run,"--cohorts","cleaned")[0],0)
  kept=(run/"cohorts/prefix/data/species"/prep.SLUG/"genomes.txt").read_text().split()
  for g in self.prefix[:5]+[BORDERLINE]: self.assertIn(g,kept)
  self.assertEqual(sorted(json.loads((run/"preparation.json").read_text())["cleaned"]["exclusions"]),sorted(RULED))

 def test_the_prefix_must_be_the_first_rows_of_full(self):
  self.write_inputs(full=self.full[1:]+self.full[:1])
  self.refused("the cleaned 2k is not the first rows of the cleaned 4k")

 def test_each_species_check_condition_stops_before_anything_is_created(self):
  self.verdicts[self.full[3000]]="fail"; self.write_inputs(count=3); self.refused("excluded genomes are in the cohorts")
  self.verdicts[self.full[3000]]="pass"; del self.verdicts[self.full[3001]]; self.write_inputs(); self.refused("have no verdict")
  self.verdicts[self.full[3001]]="pass"; self.unchecked=[REPLACEMENTS[1]]; self.write_inputs(); self.refused("unchecked list")
  self.unchecked=[]; self.verdicts[RULED[1]]="borderline"; self.write_inputs()  # EXCLUSIONS_COUNT stays 2: at least the two
  self.refused("does not fail ['1773.40949']")
  self.verdicts[RULED[1]]="fail"; self.write_inputs(count=3); self.refused("not EXCLUSIONS_COUNT 3")
  self.write_inputs(prefix_unchecked=1); self.refused("selected_unchecked is 1, not 0")
  self.write_inputs(full_unchecked=None); self.refused("selected_unchecked is None, not 0")

 def test_bindings_pins_fasta_and_placeholders_stop(self):
  self.write_inputs()
  p=Path(prep.rerun.VALUES["RN_FULL_GENOMES_PATH"]); p.write_bytes(p.read_bytes()+b"x\n"); self.refused("not the bound")
  self.write_inputs(); prep.rerun.VALUES["CLEANED_ADDED_SHA256"]="0"*64; self.refused("differ from their pins ['CLEANED_ADDED_SHA256']")
  self.write_inputs(); (self.store/f"{REPLACEMENTS[0]}.fna.gz").write_bytes(b"edited"); self.refused("differs from the replacement FASTA manifest")
  self.setUp_store_again(); self.write_inputs(); prep.rerun.VALUES["RN_FASTA_STORE"]=str(self.base); self.refused("directly under")
  self.write_inputs(); prep.rerun.VALUES["RN_INPUT_ROOT"]="@@RN_INPUT_ROOT@@"; self.refused("@@RN_INPUT_ROOT@@")
  self.write_inputs(); prep.rerun.VALUES["CLEANED_MATCHED2K_COUNT"]="@@CLEANED_MATCHED2K_COUNT@@"; self.refused("@@CLEANED_MATCHED2K_COUNT@@")
  self.write_inputs(fmt="species-check-v2"); self.refused("reviewed parser change")

 def test_the_input_roots_labels_and_catalog_are_bound_before_anything_is_read(self):
  # Review P1: RN_INPUT_ROOT's curated/datasets.json, curated/phenotypes.jsonl and catalog/version.txt each have a bound hash.
  for rel,name in prep.rerun.INPUT_ROOT_FILES.items():
   with self.subTest(rel=rel):
    self.write_inputs(); path=prep.CANDIDATE/rel; original=path.read_bytes()
    path.write_bytes(original+b"\n"); self.refused(f"RN_INPUT_ROOT's {rel}")
    path.write_bytes(original)
    self.write_inputs(); prep.rerun.VALUES[name]=f"@@{name}@@"; self.refused(f"@@{name}@@")
  self.write_inputs(); run=self.base/"run"; self.assertEqual(self.main(run,"--cohorts","cleaned")[0],0)
  bound=json.loads((run/"frozen-inputs.json").read_text())["cleaned"]["bound"]
  for rel in prep.rerun.INPUT_ROOT_FILES: self.assertIn(str(prep.CANDIDATE/rel),bound)
  ph=prep.CANDIDATE/"curated"/"phenotypes.jsonl"; ph.write_bytes(ph.read_bytes()+b"\n")
  rc,out=self.main("--verify-frozen",run); self.assertEqual(rc,1); self.assertIn(str(ph),json.loads(out)["changed"])

 def setUp_store_again(self):
  for row in self.fasta_rows:
   p=Path(row["canonical_source_path"]); p.write_bytes(gzip.compress(f">c\nTTGA{row['genome_id']}\n".encode(),mtime=0))

 def test_the_input_manifest_is_exactly_what_the_cohort_roots_can_read_and_any_change_stops(self):
  # The re-check's N1. The set: walk every prepared cohort root, follow its links into RN_INPUT_ROOT, and map each
  # copied annotation back to its source. It must equal input_root_files, which the manifest hashes.
  run=self.base/"run"; self.assertEqual(self.main(run,"--cohorts","cleaned")[0],0)
  root=prep.CANDIDATE.resolve(); seen=set()
  for name in ("prefix","added","full","matched2k"):
   data=run/"cohorts"/name/"data"
   stores=[s.resolve() for s in prep.ALLOWED_FASTA_PARENTS|{self.store}]
   for top,dirs,files in os.walk(data,followlinks=True):
    for f in files:
     p=Path(top)/f; real=p.resolve()
     if real.is_relative_to(root): seen.add(real); continue  # a link into RN_INPUT_ROOT
     rel=p.relative_to(data)
     if rel==Path("species")/prep.SLUG/"genomes.txt" or any(real.is_relative_to(s) for s in stores): continue  # written by prepare; FASTA
     # Review A3: anything else must be a copy of the input-root file at the same path, and it counts as read.
     source=root/rel
     self.assertTrue(source.is_file() and source.read_bytes()==p.read_bytes(),
                     f"{p} is no link into RN_INPUT_ROOT, no copy of one of its files and not written by preparation")
     seen.add(source)
  expected={p.resolve() for p in prep.input_root_files(prep.CANDIDATE,self.full)}
  self.assertEqual(seen,expected)
  self.assertIn(root/"curated"/"funnel.json",expected); self.assertIn(root/"bvbrc"/"genomes.jsonl",expected)
  self.assertTrue(any(p.parent.name=="features" for p in expected))
  frozen=json.loads((run/"frozen-inputs.json").read_text())["cleaned"]
  self.assertEqual(frozen["input_manifest"]["files"],len(expected))
  for p in expected: self.assertIn(str(prep.CANDIDATE/p.relative_to(root)),frozen["bound"])
  # A change to any file in the set stops preparation; one outside it does not.
  feature=next(p for p in prep.input_root_files(prep.CANDIDATE,self.full) if p.parent.name=="features")
  shutil.rmtree(run)  # a refused preparation makes no run root, so this is removed once
  for path in (prep.CANDIDATE/"curated"/"funnel.json",prep.CANDIDATE/"bvbrc"/"genomes.jsonl",feature,prep.CANDIDATE/"catalog"/"AMR_CDS.fa"):
   with self.subTest(path=path.name):
    original=path.read_bytes(); path.write_bytes(original+b"\n")
    self.refused("AMR_CDS.fa" if path.name=="AMR_CDS.fa" else "RN_INPUT_MANIFEST_SHA256")
    path.write_bytes(original)
  (prep.CANDIDATE/"curated"/"new.json").write_text("{}")
  self.refused("RN_INPUT_MANIFEST_SHA256"); (prep.CANDIDATE/"curated"/"new.json").unlink()
  outside=prep.CANDIDATE/"bvbrc"/"features"/"other.0.json"; outside.write_text(outside.read_text()+"x")
  self.assertEqual(self.main(run,"--cohorts","cleaned")[0],0)
  # And after preparation, --verify-frozen catches a change to a manifest file.
  path=prep.CANDIDATE/"curated"/"funnel.json"; path.write_text("{\"x\":1}\n")
  rc,out=self.main("--verify-frozen",run); self.assertEqual(rc,1); self.assertIn(str(path),json.loads(out)["changed"])

 def test_a_missing_or_unusable_cached_annotation_stops_before_anything_is_created(self):
  # grm learn fetches an annotation it cannot read or parse (src/annotate.rs:101-117); the cleaned run must fetch none.
  added=self.full[2500]; replacement=REPLACEMENTS[0]
  for gid,edit in ((added,"delete"),(replacement,"delete"),(added,"garble"),(added,"object")):
   with self.subTest(gid=gid,edit=edit):
    path=prep.CANDIDATE/"bvbrc"/"features"/f"{gid}.json"; original=path.read_bytes()
    if edit=="delete": path.unlink()
    else: path.write_text("not json" if edit=="garble" else '{"gene": "x"}')
    self.write_inputs()
    message=self.refused("1 cohort genomes have no usable cached annotation")
    self.assertIn(gid,message)
    path.write_bytes(original)
  for g in self.full[3000:3030]: (prep.CANDIDATE/"bvbrc"/"features"/f"{g}.json").unlink()
  self.write_inputs(); message=self.refused("30 cohort genomes have no usable cached annotation")
  self.assertIn("first 20",message); self.assertEqual(sum(g in message for g in self.full[3000:3030]),20)
  for g in self.full[3000:3030]: (prep.CANDIDATE/"bvbrc"/"features"/f"{g}.json").write_text(f'[{{"gene": "{g}"}}]')
  self.write_inputs(); self.assertEqual(self.main(self.base/"run","--cohorts","cleaned")[0],0)

 def test_a_missing_record_file_is_a_clean_refusal_and_verify_frozen_reports_it(self):
  # Review A2: no traceback, nothing created.
  self.write_inputs(); rec=prep.CANDIDATE/"bvbrc"/"genomes.jsonl"; original=rec.read_bytes(); rec.unlink()
  self.refused("has no regular bvbrc/genomes.jsonl")
  rec.write_bytes(original); self.write_inputs()
  run=self.base/"run"; self.assertEqual(self.main(run,"--cohorts","cleaned")[0],0)
  rec.unlink()
  rc,out=self.main("--verify-frozen",run); self.assertEqual(rc,1); self.assertIn(str(rec),json.loads(out)["changed"])
  rec.write_bytes(original)

 def test_matched_selection_takes_the_prefix_size(self):
  ids=[f"g{i}" for i in range(4000)]
  self.assertEqual(prep.matched_selection(ids),prep.matched_selection(ids,2000))
  m=prep.matched_selection(ids[:3998],1999); self.assertEqual(len(m),1999); self.assertEqual(m,[g for g in ids if g in set(m)])

if __name__=="__main__":
 result=unittest.main(exit=False,verbosity=2).result
 sys.exit(0 if result.wasSuccessful() else 1)
