#!/usr/bin/env python3
"""Evaluation profiles of report_contract.py: r8 and rN (revision N on the frozen cohort rules).

usage: FOLDGEN=<foldgen binary> python3 -I scripts/mtb-4k/tests/test_evaluation_profiles.py

- --profile is read from the process's arguments at import, as argparse reads it;
- rN refuses, and leaves the active profile as it was, while revision_n.py holds a placeholder,
  when a value is malformed, when revision N drops or renames a report or status field, when its
  checker judges convergence statuses unlike 690815a, and when a pin differs from the in-tree blob;
- rN takes exactly the revision-9 facts of 2026-10-29 (b2fd433), and against this tree, whose
  checker blobs are revision 8's, those facts stop at the pin check;
- with values whose pins are this tree's blobs, rN activates, and its revision reaches every
  check: validate_manifest, the runner's checker and audit calls, the postflight, the summary
  coverage and the comparator's receipt check. This tree's check-results.py accepts revisions 3-8
  only, so a revision-9 root stops at the checker; a passing rN run needs revision N's checker.
"""
import contextlib,copy,io,json,sys,tempfile,unittest
from pathlib import Path
LANE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(LANE)); sys.path.insert(0,str(LANE/"tests"))
import report_contract as vc
import revision_n
import run_revision8_checks as wrapper
import compare_cohorts as cc
import fixtures as fx

# The revision-9 facts relayed on 2026-10-29 (r9/lasso-exchange-fix, b2fd433), in revision_n.py's form.
R9={"REVISION_N":9,"RN_CHECKER_COMMIT":"b2fd4337ce88be64767a3a552e87df21e64803c4",
 "RN_PIN_CHECK_RESULTS":"b7e5a8a067179bc4578d6487aebe2baad51f4fe419a44aa2e3a942ae3d10a347",
 "RN_PIN_METRIC_AUDIT":"c1af2ead6fbedb869a067509efd2a431a218a9ddc520fdfb3625d2054e4af98b",
 "RN_PIN_SUMMARIZE":"fd5fdfab1bedb8118c95eaadef19452fc89786f371a5d1757bd784844eb3295f",
 "RN_PIN_RESULT_SCOPE":"3bcb2db34b23accd9186a595e6e0002cbd3bc292ab548442d8a53ba94b65ae2a",
 "RN_FIELD_CHANGES":{},"RN_CONVERGENCE_STATUS_RULE":"as-690815a"}
# The same with this tree's (revision-8) checker blobs as pins, so rN can activate here.
IN_TREE={**R9,**{name:vc.R8_PINS[rel] for rel,name in revision_n.PIN_NAMES.items()}}

def quiet(fn):
 with contextlib.redirect_stdout(io.StringIO()),contextlib.redirect_stderr(io.StringIO()) as err:
  return fn(),err

@contextlib.contextmanager
def profile(name,values=None):
 vc.activate(name,values)
 try: yield
 finally: vc.activate("r8")

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

class Selection(unittest.TestCase):
 def test_profile_is_read_from_argv_as_argparse_reads_it(self):
  cases=[([],"r8"),(["out.json","--root","x"],"r8"),(["--profile","rN"],"rN"),(["--profile=rN"],"rN"),
   (["--profile","rN","--profile","r8"],"r8"),(["--","--profile","rN"],"r8"),(["--profile"],"r8")]
  for argv,want in cases:
   with self.subTest(argv=argv): self.assertEqual(vc.profile_from_argv(argv),want)
 def test_this_process_runs_r8(self):
  self.assertEqual((vc.PROFILE,vc.REVISION,vc.PINS,vc.RULES),("r8",8,vc.R8_PINS,vc.R8_RULES))
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
 def test_the_committed_placeholders_refuse_and_name_every_one(self):
  message=self.refused(None,"unfilled placeholder")
  for name in revision_n.CONTRACT_NAMES: self.assertIn(f"@@{name}@@",message)
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
 def test_the_revision_9_facts_are_taken_and_stop_at_this_trees_revision_8_checker(self):
  self.assertEqual(revision_n.require(list(revision_n.TABLE),R9),R9)
  self.refused(R9,"check-results.py","differs from the reviewed pin b7e5a8a0")

class Active(unittest.TestCase):
 def test_rn_binds_revision_pins_rules_and_modules(self):
  with profile("rN",IN_TREE):
   self.assertEqual((vc.PROFILE,vc.REVISION),("rN",9))
   self.assertEqual(vc.PINS,vc.R8_PINS); self.assertIs(vc.RULES,vc.R8_RULES)
   self.assertIsNotNone(vc.upstream_checker.check_lasso_fit)
  self.assertEqual((vc.PROFILE,vc.REVISION),("r8",8))
 def test_an_added_field_extends_its_key_set_and_r8_restores_it(self):
  fit_keys=sorted(vc.BASE_KEYS["FIT_KEYS"]|{"exchange_steps"})
  with profile("rN",{**IN_TREE,"RN_FIELD_CHANGES":{"FIT_KEYS":fit_keys}}):
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
    restamp(root,9)
    cohort=vc.validate_manifest(root/"validation-manifest.json",rules)
    self.assertEqual({r.receipt["evaluation_revision"] for r in cohort.reports.values()},{9})
   with self.assertRaises(vc.ValidationError): vc.validate_manifest(root/"validation-manifest.json",rules)
 def test_the_revision_reaches_the_runner_and_this_trees_checker_refuses_it(self):
  with tempfile.TemporaryDirectory() as tmp:
   rules=fx.make_rules(); root=fx.build_cohort(Path(tmp)/"prefix","prefix",rules)
   args=[str(root/"audit.json"),"--root",str(root),"--run-manifest",str(root/"validation-manifest.json"),
    "--expected-source-commit",fx.SOURCE,"--binary",str(root/"synthetic-binary"),"--backup",str(root/"backup"),
    "--log",str(root/"pipeline.log"),"--launch","2000-01-01T00:00:00Z","--foldgen",fx.foldgen_path(),
    "--expected-foldgen-sha256",fx.foldgen_sha256(),"--checker-out",str(root/"checker.json"),"--summary-out",str(root/"summary.json"),
    "--receipt",str(root/"receipt.json")]
   with profile("rN",IN_TREE):
    restamp(root,9)
    code,err=quiet(lambda:wrapper.main(args+["--profile","r8"],rules=rules))
    self.assertEqual(code,2); self.assertIn("--profile r8 is not the active profile rN",err.getvalue())
    code,err=quiet(lambda:wrapper.main(args+["--profile","rN"],rules=rules))
    self.assertEqual(code,2); self.assertIn("revision-9 checker exited",err.getvalue())
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
    for needle in ("auditor revision is not 9","historical/wrong revisions in audited set","manifest schema/revision differs",
                   "checker accepted mixed revisions"): self.assertIn(needle,errors)
    self.assertTrue(any("wrong revision/schema/validity" in e for e in errors))
    self.assertEqual(wrapper.summary_errors(paths["summary"],cohort),["summarizer coverage/issues differ"])
    with self.assertRaises(vc.ValidationError) as stop: cc.load_receipt(receipt,cohort,fx.SOURCE,fx.foldgen_sha256())
    self.assertIn("receipt is not for evaluation revision 9",str(stop.exception))
    # The same outputs restamped to revision 9 pass the revision-9 postflight.
    nine=copy.deepcopy(bodies); s=nine["audit"]["summary"]
    s.update(revision=9,by_revision={"9":s["files"]},manifest_evaluation_revision=9)
    for rec in nine["audit"]["records"]: rec.update(revision=9,published_revision=9)
    for rec in nine["checker"]["reports"]: rec["revision"]=9
    nine["summary"]["coverage"]={"9":nine["summary"]["coverage"]["8"]}
    for key,body in nine.items(): fx.write_json(paths[key],body)
    self.assertEqual(wrapper.postflight(paths["audit"],cohort,paths["checker"]),[])
    self.assertEqual(wrapper.summary_errors(paths["summary"],cohort),[])

if __name__=="__main__":
 result=unittest.main(exit=False,verbosity=2).result
 sys.exit(0 if result.wasSuccessful() else 1)
