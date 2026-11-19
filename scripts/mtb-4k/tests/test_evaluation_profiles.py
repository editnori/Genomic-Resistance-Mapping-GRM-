#!/usr/bin/env python3
"""Evaluation profiles of report_contract.py: r8 and rN (revision N on the frozen cohort rules).

usage: FOLDGEN=<foldgen binary> [MTB4K_TEST_PROFILE=r8|rN] python3 -I scripts/mtb-4k/tests/test_evaluation_profiles.py

The suite runs under MTB4K_TEST_PROFILE (profile_boot.py; r8 by default). These tests hold under
either: "the tree" below is this checkout's own checker blobs, and OTHER is a revision other than
the suite's (9 under r8, 10 under a revision-9 rN), so the revision a profile binds is always
seen to move.

- --profile is read from the process's arguments at import, as argparse reads it;
- rN refuses, and leaves the active profile as it was, while revision_n.py's values hold a
  placeholder (also with MTB4K_TEST_PROFILE=rN set, in a subprocess), when a value is malformed,
  when revision N drops or renames a field, when its checker judges convergence statuses unlike
  690815a, and when a pin differs from the tree's blob;
- rN takes exactly the revision-9 facts of 2026-10-29 (b2fd433): they activate on a tree that
  holds b2fd433's blobs and stop at the pin check on any other;
- with the tree's pins, rN activates at OTHER, and that revision reaches every check:
  validate_manifest, the runner's checker and audit calls, the postflight, the summary coverage
  and the comparator's receipt check. A root restamped to OTHER stops at the tree's checker, which
  accepts revisions up to the suite's only.
"""
import contextlib,copy,hashlib,io,json,os,re,shutil,subprocess,sys,tempfile,unittest
from pathlib import Path
LANE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(LANE)); sys.path.insert(0,str(LANE/"tests"))
import profile_boot as pb  # noqa: E402  (MTB4K_TEST_PROFILE, before report_contract loads)
import report_contract as vc
import revision_n
import cleaned_rerun
import run_revision8_checks as wrapper
import compare_cohorts as cc
import fixtures as fx

BASE=pb.PROFILE
# The revision-9 facts relayed on 2026-10-29 (r9/lasso-exchange-fix, b2fd433), in revision_n.py's form.
R9={"REVISION_N":9,"RN_CHECKER_COMMIT":"b2fd4337ce88be64767a3a552e87df21e64803c4",
 "RN_PIN_CHECK_RESULTS":"b7e5a8a067179bc4578d6487aebe2baad51f4fe419a44aa2e3a942ae3d10a347",
 "RN_PIN_METRIC_AUDIT":"c1af2ead6fbedb869a067509efd2a431a218a9ddc520fdfb3625d2054e4af98b",
 "RN_PIN_SUMMARIZE":"fd5fdfab1bedb8118c95eaadef19452fc89786f371a5d1757bd784844eb3295f",
 "RN_PIN_RESULT_SCOPE":"3bcb2db34b23accd9186a595e6e0002cbd3bc292ab548442d8a53ba94b65ae2a",
 "RN_FIELD_CHANGES":{},"RN_CONVERGENCE_STATUS_RULE":"as-690815a"}
TREE_PINS={name:vc.file_sha256(vc.SCRIPTS/rel) for rel,name in revision_n.PIN_NAMES.items()}
# A revision the tree's checker does not accept: 9 on revision 8's blobs, 10 on revision 9's, 11 on revision 10's.
OTHER=max(vc.REVISION,8)+1
TEN=vc.DUPLICATE_ISOLATION_REVISION
DECLARED={"OPTIONAL_TOP_KEYS":["kset_duplicates"]}  # revision 10's RN_FIELD_CHANGES (c90ebaf)
IN_TREE={**R9,**TREE_PINS,"REVISION_N":OTHER,"RN_FIELD_CHANGES":DECLARED if OTHER>=TEN else {}}
AT_TEN={**R9,**TREE_PINS,"REVISION_N":TEN,"RN_FIELD_CHANGES":DECLARED}
UNFILLED={n:f"@@{n}@@" for n in revision_n.TABLE}

def quiet(fn):
 with contextlib.redirect_stdout(io.StringIO()),contextlib.redirect_stderr(io.StringIO()) as err:
  return fn(),err

@contextlib.contextmanager
def profile(name,values=None):
 vc.activate(name,values)
 try: yield
 finally: vc.activate(BASE)

def state(): return (vc.PROFILE,vc.REVISION,dict(vc.PINS),vc.RULES,vc.upstream_checker,set(vc.FIT_KEYS))

def restamp(root,revision):
 """A fixture root whose reports, index and run manifest carry another evaluation revision."""
 for p in (root/"data/species"/vc.SPECIES/"results").glob("*.json"):
  r=json.loads(p.read_text()); r["evaluation"]["revision"]=revision; fx.write_json(p,r)
 fx.refresh_root(root,"prefix")
 index=root/"data/report/index.json"; i=json.loads(index.read_text()); i["evaluation_revision"]=revision; fx.write_json(index,i)
 m=json.loads((root/"validation-manifest.json").read_text()); m["evaluation_revision"]=revision
 for e in m["inputs"]: e["sha256"]=fx.sha(root/e["path"])
 fx.write_json(root/"validation-manifest.json",m)

def unfilled_lane(directory):
 """A copy of scripts/mtb-4k (no tests but profile_boot.py) and the four pinned checker files, with
 every revision_n.py value reset to its placeholder, whatever this tree holds."""
 scripts=Path(directory)/"scripts"
 shutil.copytree(LANE,scripts/"mtb-4k",ignore=shutil.ignore_patterns("__pycache__"))
 for rel in revision_n.PIN_NAMES:
  (scripts/rel).parent.mkdir(parents=True,exist_ok=True); shutil.copyfile(vc.SCRIPTS/rel,scripts/rel)
 path=scripts/"mtb-4k/revision_n.py"; text=path.read_text()
 start=text.index("\nVALUES = {\n")+1; end=text.index("\n}\n",start)+3
 path.write_text(text[:start]+"VALUES = {\n"+"".join(f'    "{n}": "@@{n}@@",\n' for n in revision_n.TABLE)+"}\n"+text[end:])
 return scripts/"mtb-4k"

class Selection(unittest.TestCase):
 def test_profile_is_read_from_argv_as_argparse_reads_it(self):
  cases=[([],"r8"),(["out.json","--root","x"],"r8"),(["--profile","rN"],"rN"),(["--profile=rN"],"rN"),
   (["--profile","rN","--profile","r8"],"r8"),(["--","--profile","rN"],"r8"),(["--profile"],"r8")]
  for argv,want in cases:
   with self.subTest(argv=argv): self.assertEqual(vc.profile_from_argv(argv),want)
 def test_this_process_runs_the_suites_profile(self):
  self.assertEqual(vc.PROFILE,BASE)
  self.assertEqual(vc.PINS,{rel:TREE_PINS[name] for rel,name in revision_n.PIN_NAMES.items()})
  self.assertIs(vc.RULES,vc.R8_RULES)
  if BASE=="r8": self.assertEqual((vc.REVISION,vc.PINS),(8,vc.R8_PINS))
 def test_an_unknown_profile_is_refused(self):
  before=state()
  with self.assertRaises(vc.ValidationError) as stop: vc.activate("r9")
  self.assertIn("unknown evaluation profile",str(stop.exception)); self.assertEqual(state(),before)

class Refusals(unittest.TestCase):
 def refused(self,values,*needles):
  before=state()
  with self.assertRaises(vc.ValidationError) as stop: vc.activate("rN",values)
  for needle in needles: self.assertIn(needle,str(stop.exception))
  self.assertEqual(state(),before,"a refused activation changed the active profile")
  return str(stop.exception)
 def test_unfilled_values_refuse_and_name_every_one(self):
  message=self.refused(UNFILLED,"unfilled placeholder")
  for name in revision_n.CONTRACT_NAMES: self.assertIn(f"@@{name}@@",message)
 def test_an_unfilled_rn_refuses_with_the_switch_set(self):
  with tempfile.TemporaryDirectory() as tmp:
   lane=unfilled_lane(tmp)
   p=subprocess.run([sys.executable,"-I","-B","-c","import sys; sys.path.insert(0,sys.argv[1]); import profile_boot",str(lane/"tests")],
                    cwd=tmp,env={**os.environ,"MTB4K_TEST_PROFILE":"rN"},capture_output=True,text=True)
   self.assertNotEqual(p.returncode,0,p.stdout)
   self.assertIn("profile rN: unfilled placeholder(s): @@REVISION_N@@",p.stderr)
 def test_one_placeholder_left_is_enough_to_refuse(self):
  for name in revision_n.CONTRACT_NAMES:
   with self.subTest(name=name): self.refused({**IN_TREE,name:f"@@{name}@@"},f"@@{name}@@")
 def test_malformed_values_refuse(self):
  for name,bad in [("REVISION_N",8),("REVISION_N","9"),("REVISION_N",True),("RN_PIN_CHECK_RESULTS","B7E5"+"0"*60),
                   ("RN_PIN_METRIC_AUDIT","0"*63),("RN_FIELD_CHANGES",{"NEW_KEYS":["x"]}),("RN_FIELD_CHANGES",[]),
                   ("RN_CONVERGENCE_STATUS_RULE","")]:
   with self.subTest(name=name,bad=bad): self.refused({**IN_TREE,name:bad},name)
 def test_a_dropped_or_renamed_field_needs_a_contract_change(self):
  renamed=sorted(vc.BASE_KEYS["CONVERGENCE_KEYS"]-{"sweep_cap_hits"}|{"cap_hits"})
  self.refused({**IN_TREE,"RN_FIELD_CHANGES":{"CONVERGENCE_KEYS":renamed}},"drops or renames ['sweep_cap_hits']")
 def test_other_convergence_semantics_need_a_contract_change(self):
  self.refused({**IN_TREE,"RN_CONVERGENCE_STATUS_RULE":"every-step-capped-is-unconverged"},"validate_convergence needs a reviewed change")
 def test_a_pin_that_is_not_the_trees_blob_refuses(self):
  self.refused({**IN_TREE,"RN_PIN_CHECK_RESULTS":"0"*64},"check-results.py","differs from the reviewed pin 0000")
 def test_the_revision_9_facts_are_taken_and_activate_only_on_b2fd433s_blobs(self):
  self.assertEqual(revision_n.require(list(revision_n.TABLE),R9),R9)
  if all(TREE_PINS[n]==R9[n] for n in TREE_PINS):
   with profile("rN",R9): self.assertEqual((vc.PROFILE,vc.REVISION),("rN",9))
  else:
   self.refused(R9,"differs from the reviewed pin")

class Cleaned(unittest.TestCase):
 """Profile rN-cleaned: revision N's values with the cleaned cohort pins of plan revision 6."""
 def pins(self,**edit):
  v={n:f"@@{n}@@" for n in cleaned_rerun.TABLE}
  for name in fx.COHORTS:
   sha,count=cleaned_rerun.COHORT_PINS[name]
   v[sha]=hashlib.sha256(fx.genome_bytes(name)).hexdigest(); v[count]=len(fx.cohort_genomes(name))
  v.update(edit); return v
 def test_unfilled_values_of_either_file_refuse_naming_both(self):
  before=state()
  with self.assertRaises(vc.ValidationError) as stop: vc.activate("rN-cleaned",UNFILLED,{n:f"@@{n}@@" for n in cleaned_rerun.TABLE})
  self.assertIn("@@REVISION_N@@",str(stop.exception)); self.assertIn("@@CLEANED_MATCHED2K_COUNT@@",str(stop.exception))
  with self.assertRaises(vc.ValidationError) as stop: vc.activate("rN-cleaned",IN_TREE,self.pins(CLEANED_ADDED_COUNT="@@CLEANED_ADDED_COUNT@@"))
  self.assertIn("@@CLEANED_ADDED_COUNT@@",str(stop.exception)); self.assertNotIn("@@REVISION_N@@",str(stop.exception))
  self.assertEqual(state(),before)
 def test_the_cleaned_rules_bind_the_pins_and_keep_the_drug_sets(self):
  vc.activate("rN-cleaned",IN_TREE,self.pins())
  try:
   self.assertEqual((vc.PROFILE,vc.REVISION),("rN-cleaned",OTHER))
   for name in fx.COHORTS:
    self.assertEqual((vc.RULES[name].size,vc.RULES[name].genomes_sha256),(len(fx.cohort_genomes(name)),__import__("hashlib").sha256(fx.genome_bytes(name)).hexdigest()))
   self.assertEqual(vc.RULES["prefix"].drugs,vc.ALL_DRUGS); self.assertEqual(vc.RULES["full"].drugs,vc.ALL_DRUGS)
   self.assertEqual((vc.RULES["added"].drugs,vc.RULES["added"].unavailable),(vc.ADDED_DRUGS,None))
   self.assertIsNone(vc.RULES["matched2k"].drugs)
   with tempfile.TemporaryDirectory() as tmp:
    root=fx.build_cohort(Path(tmp)/"added","added",fx.make_rules())
    restamp(root,OTHER)
    cohort=vc.validate_manifest(root/"validation-manifest.json",vc.RULES)  # the eight receipts derived from the phenotypes
    self.assertEqual(len(cohort.reports),14)
    vc.activate("rN-cleaned",IN_TREE,self.pins(CLEANED_ADDED_SHA256="0"*64))
    with self.assertRaises(vc.ValidationError) as stop: vc.validate_manifest(root/"validation-manifest.json",vc.RULES)
    self.assertIn("cohort genome-list hash differs from frozen rule",str(stop.exception))
  finally: vc.activate(BASE)

class KsetDuplicates(unittest.TestCase):
 """Revision 10's kset_duplicates (duplicate isolation; report_contract.py cites c90ebaf's audit.py and check-results.py)."""
 DID=vc.DRUG_IDS["amikacin"]
 def setUp(self):
  self.tmp=tempfile.TemporaryDirectory(); vc.activate("rN",AT_TEN)
  self.rules=fx.make_rules({"prefix":("amikacin",)}); self.root=fx.build_cohort(Path(self.tmp.name)/"prefix","prefix",self.rules)
  restamp(self.root,TEN); self.report_path=self.root/"data/species"/vc.SPECIES/"results"/(self.DID+".json")
 def tearDown(self):
  vc.activate(BASE); self.tmp.cleanup()
 def report(self): return json.loads(self.report_path.read_text())
 def pair(self,same_random=True,avoid=()):
  """Two samples of one lineage, in the same random test fold or not."""
  r=self.report(); rf=next(s for s in r["splits"] if s["scheme"]=="random")["test_fold"]; s=r["samples"]
  for i in range(len(s)):
   for j in range(i+1,len(s)):
    if s[i]["lineage"]==s[j]["lineage"] and (rf[i]==rf[j])==same_random and not {i,j}&set(avoid): return i,j
  raise AssertionError("no such pair in the fixture")
 def share(self,i,j):
  """Genome j takes genome i's kset digest in matrix.json, structure.json and counted.json, as identical ksets would."""
  r=self.report(); a,b=r["samples"][i]["id"],r["samples"][j]["id"]; sp=self.root/"data/species"/vc.SPECIES
  for f in (sp/"matrix/matrix.json",sp/"structure.json"):
   m=json.loads(f.read_text()); m["kset_digests"][m["genomes"].index(b)]=m["kset_digests"][m["genomes"].index(a)]; fx.write_json(f,m)
  c=json.loads((sp/"counted.json").read_text()); by={g["genome_id"]:g for g in c["genomes"]}
  by[b]["kset_digest"]=by[a]["kset_digest"]; fx.write_json(sp/"counted.json",c)
 def store(self,edit):
  """The report edited, its fingerprint recomputed from the root, and the run manifest rebound."""
  r=self.report(); edit(r); r["input_fingerprint"]=fx.fpmod.Inputs(self.root/"data",vc.SPECIES).fingerprint(self.DID,r["eval_config"])
  fx.write_json(self.report_path,r); fx.refresh_root(self.root,"prefix")
  m=json.loads((self.root/"validation-manifest.json").read_text())
  m["reports"][self.DID].update(input_fingerprint=r["input_fingerprint"],sha256=fx.sha(self.report_path))
  for e in m["inputs"]: e["sha256"]=fx.sha(self.root/e["path"])
  for key in ("matrix","genomes"): m[key]["sha256"]=fx.sha(self.root/m[key]["path"])
  fx.write_json(self.root/"validation-manifest.json",m)
 def validate(self): return vc.validate_manifest(self.root/"validation-manifest.json",self.rules)
 def refused(self,needle):
  with self.assertRaises(vc.ValidationError) as stop: self.validate()
  self.assertIn(needle,str(stop.exception)); return str(stop.exception)
 def test_without_a_group_the_field_is_absent_and_the_report_passes(self):
  self.store(lambda r:None); self.validate()
  self.store(lambda r:r.update(kset_duplicates=[])); self.refused("must be a non-empty list when present")
 def test_a_group_is_listed_exactly_as_the_digests_give(self):
  i,j=self.pair(); self.share(i,j)
  self.store(lambda r:None); self.refused("kset_duplicates is missing: the root's counted.json digests give 1 groups")
  self.store(lambda r:r.update(kset_duplicates=[[i,j]])); cohort=self.validate(); self.assertIn(self.DID,cohort.reports)
  k,l=self.pair(avoid=(i,j))  # one lineage, one random fold, but distinct digests
  cases=[("wrong group",[[k,l]],"differs from the root's counted.json kset digests"),
         ("unsorted group",[[j,i]],"must be an ascending list"),
         ("overlapping groups",[[i,j],[i,j]],"groups out of order, overlapping or out of range"),
         ("a one-member group",[[i]],"must be an ascending list"),
         ("an index out of range",[[i,10**6]],"groups out of order, overlapping or out of range"),
         ("not a list",{"0":[i,j]},"must be a non-empty list when present"),
         ("empty list",[],"must be a non-empty list when present")]
  for label,value,needle in cases:
   with self.subTest(label):
    self.store(lambda r,value=value:r.update(kset_duplicates=value)); self.refused(needle)
 def test_members_split_across_random_folds_refuse(self):
  i,j=self.pair(same_random=False); self.share(i,j)
  self.store(lambda r:r.update(kset_duplicates=[[i,j]])); self.refused("a group spans test folds in the random split")
 def test_the_field_below_the_threshold_refuses(self):
  i,j=self.pair(); self.store(lambda r:r.update(kset_duplicates=[[i,j]]))
  vc.activate("rN",{**R9,**TREE_PINS,"REVISION_N":TEN-1,"RN_FIELD_CHANGES":{}}); restamp(self.root,TEN-1)
  self.store(lambda r:r.update(kset_duplicates=[[i,j]])); self.refused("unknown=['kset_duplicates']")
 def test_the_declaration_is_required_exactly_from_the_threshold(self):
  before=state()
  for values,needle in (({**AT_TEN,"RN_FIELD_CHANGES":{}},"must declare it under OPTIONAL_TOP_KEYS"),
                        ({**AT_TEN,"REVISION_N":TEN-1},"must not declare it under OPTIONAL_TOP_KEYS"),
                        ({**AT_TEN,"RN_FIELD_CHANGES":{"OPTIONAL_TOP_KEYS":["kset_duplicates","other"]}},"cannot validate: ['other']")):
   with self.subTest(needle):
    with self.assertRaises(vc.ValidationError) as stop: vc.activate("rN",values)
    self.assertIn(needle,str(stop.exception)); self.assertEqual(state(),before)
 def test_an_l1_failure_is_the_invalid_report_stop(self):
  # Known limit L1: a move can empty a class from a random test fold, and the dataset's evaluation fails. The
  # comparator loads strictly (L3), so such a report is the existing invalid-report stop, never a crash.
  self.store(lambda r:r["evaluation"].update(valid=False,reasons=["invalid random outer evaluation"]))
  with self.assertRaises(vc.ValidationError): self.validate()
  c=vc.validate_manifest(self.root/"validation-manifest.json",self.rules,collect_invalid=True); self.assertIn(self.DID,c.invalid)
 def test_the_grouping_matches_audit_pys_own_when_the_tree_has_it(self):
  samples=[vc.Sample(f"g{i}",True,0) for i in range(8)]; digests={"g0":"a","g1":"b","g2":"a","g3":"c","g4":"b","g5":"d","g6":"a","g7":"e"}
  want=[[0,2,6],[1,4]]; self.assertEqual(vc.derived_kset_duplicates(samples,digests,"synthetic"),want)
  if hasattr(vc.upstream_audit,"kset_duplicates"):  # audit.py from c90ebaf on
   self.assertEqual(vc.upstream_audit.kset_duplicates([{"id":s.id} for s in samples],digests),want)

class Active(unittest.TestCase):
 def test_rn_binds_revision_pins_rules_and_modules(self):
  with profile("rN",IN_TREE):
   self.assertEqual((vc.PROFILE,vc.REVISION),("rN",OTHER))
   self.assertEqual(vc.PINS,{rel:TREE_PINS[name] for rel,name in revision_n.PIN_NAMES.items()}); self.assertIs(vc.RULES,vc.R8_RULES)
   self.assertIsNotNone(vc.upstream_checker.check_lasso_fit)
  self.assertEqual(vc.PROFILE,BASE)
 def test_an_added_field_extends_its_key_set_and_the_base_profile_restores_it(self):
  fit_keys=sorted(vc.BASE_KEYS["FIT_KEYS"]|{"exchange_steps"})
  with profile("rN",{**IN_TREE,"RN_FIELD_CHANGES":{**IN_TREE["RN_FIELD_CHANGES"],"FIT_KEYS":fit_keys}}):
   self.assertIn("exchange_steps",vc.FIT_KEYS)
   fit=fx.fit(); fit["exchange_steps"]=0
   vc.expect_keys(fit,vc.FIT_KEYS,"synthetic")
  self.assertNotIn("exchange_steps",vc.FIT_KEYS)
  with self.assertRaises(vc.ValidationError): vc.expect_keys(fit,vc.FIT_KEYS,"synthetic")
 def test_the_revision_reaches_the_contract(self):
  with tempfile.TemporaryDirectory() as tmp:
   rules=fx.make_rules(); root=fx.build_cohort(Path(tmp)/"prefix","prefix",rules)
   vc.validate_manifest(root/"validation-manifest.json",rules)
   with profile("rN",IN_TREE):
    with self.assertRaises(vc.ValidationError) as stop: vc.validate_manifest(root/"validation-manifest.json",rules)
    self.assertIn("wrong species/revision/schema",str(stop.exception))
    restamp(root,OTHER)
    cohort=vc.validate_manifest(root/"validation-manifest.json",rules)
    self.assertEqual({r.receipt["evaluation_revision"] for r in cohort.reports.values()},{OTHER})
   with self.assertRaises(vc.ValidationError): vc.validate_manifest(root/"validation-manifest.json",rules)
 def test_the_revision_reaches_the_runner_and_the_trees_checker_refuses_it(self):
  with tempfile.TemporaryDirectory() as tmp:
   rules=fx.make_rules(); root=fx.build_cohort(Path(tmp)/"prefix","prefix",rules)
   args=[str(root/"audit.json"),"--root",str(root),"--run-manifest",str(root/"validation-manifest.json"),
    "--expected-source-commit",fx.SOURCE,"--binary",str(root/"synthetic-binary"),"--backup",str(root/"backup"),
    "--log",str(root/"pipeline.log"),"--launch","2000-01-01T00:00:00Z","--foldgen",fx.foldgen_path(),
    "--expected-foldgen-sha256",fx.foldgen_sha256(),"--checker-out",str(root/"checker.json"),"--summary-out",str(root/"summary.json"),
    "--receipt",str(root/"receipt.json")]
   with profile("rN",IN_TREE):
    restamp(root,OTHER)
    code,err=quiet(lambda:wrapper.main(args+["--profile","r8"],rules=rules))
    self.assertEqual(code,2); self.assertIn("--profile r8 is not the active profile rN",err.getvalue())
    code,err=quiet(lambda:wrapper.main(args+["--profile","rN"],rules=rules))
    self.assertEqual(code,2); self.assertIn(f"revision-{OTHER} checker exited",err.getvalue())
    self.assertIn("('revision',",err.getvalue())
    self.assertFalse((root/"audit.json").exists() or (root/"receipt.json").exists())
 def test_the_revision_reaches_the_postflight_summary_and_comparator(self):
  with tempfile.TemporaryDirectory() as tmp:
   rules=fx.make_rules(); root=fx.build_cohort(Path(tmp)/"prefix","prefix",rules)
   cohort=vc.validate_manifest(root/"validation-manifest.json",rules)
   bodies=fx.passing_outputs(root,cohort.manifest)
   paths={}
   for key,body in bodies.items(): paths[key]=root/f"{key}.json"; fx.write_json(paths[key],body)
   self.assertEqual(wrapper.postflight(paths["audit"],cohort,paths["checker"]),[])
   self.assertEqual(wrapper.summary_errors(paths["summary"],cohort),[])
   receipt=fx.write_receipt(root)
   with profile("rN",IN_TREE):
    errors=wrapper.postflight(paths["audit"],cohort,paths["checker"])
    for needle in (f"auditor revision is not {OTHER}","historical/wrong revisions in audited set","manifest schema/revision differs",
                   "checker accepted mixed revisions"): self.assertIn(needle,errors)
    self.assertTrue(any("wrong revision/schema/validity" in e for e in errors))
    self.assertEqual(wrapper.summary_errors(paths["summary"],cohort),["summarizer coverage/issues differ"])
    with self.assertRaises(vc.ValidationError) as stop: cc.load_receipt(receipt,cohort,fx.SOURCE,fx.foldgen_sha256())
    self.assertIn(f"receipt is not for evaluation revision {OTHER}",str(stop.exception))
    # The same outputs restamped to OTHER pass that revision's postflight.
    other=copy.deepcopy(bodies); s=other["audit"]["summary"]; key=str(OTHER)
    s.update(revision=OTHER,by_revision={key:s["files"]},manifest_evaluation_revision=OTHER)
    for rec in other["audit"]["records"]: rec.update(revision=OTHER,published_revision=OTHER)
    for rec in other["checker"]["reports"]: rec["revision"]=OTHER
    other["summary"]["coverage"]={key:other["summary"]["coverage"][str(cohort.manifest["evaluation_revision"])]}
    for key_,body in other.items(): fx.write_json(paths[key_],body)
    self.assertEqual(wrapper.postflight(paths["audit"],cohort,paths["checker"]),[])
    self.assertEqual(wrapper.summary_errors(paths["summary"],cohort),[])

if __name__=="__main__":
 result=unittest.main(exit=False,verbosity=2).result
 sys.exit(0 if result.wasSuccessful() else 1)
