#!/usr/bin/env python3
"""run_revision8_checks.py: package K's WrapperTests with the real foldgen, plus the receipt.

usage: FOLDGEN=<foldgen binary> python3 -I scripts/mtb-4k/tests/test_run_revision8_checks.py

Ported from K's tests/test_validators.py WrapperTests. Path-adapted: the in-tree checker and audit
run instead of vendored copies, every run also writes --receipt, and --foldgen is the committed
foldgen, which these fixtures' folds come from (K used an in-process double that emitted its index
folds). The double stays for one injected-failure test. Every K assertion is kept.
IRLS0 (marked): "wrapper final refit zero-step convergence" changed, not dropped. Its refit
status is now 0 steps with 1 cap hit (still a stop before the checker), and it is named "wrapper
final refit zero-step convergence with a cap hit". 690815a writes irls_steps 0 when the warm
start already meets the KKT tolerance (src/learn/lasso.rs:70), the pinned check-results.py:114
accepts any non-negative step count, and the K review treated that status as impossible because
it was at 855bdba.
R1 (marked): two statuses changed, not dropped. "wrapper final selected impossible status and
excess path" now carries {1,2,true} (was {1,1,true}; still a stop), and "wrapper fold refit
all-steps-capped convergence" ({1,1,true}) became "wrapper fold refit cap hits exceed steps"
({2,3,true}; still a stop before the checker). A new positive control passes the complete
wrapper with {1,1,true} and {3,3,true} on every fit: 690815a can write them, and the pinned
check-results.py accepts them.

Follow-up tests (plan-validators review): I4, the runner refuses a foldgen that is not the
recorded build (--expected-foldgen-sha256, required); L2, verify_outputs, which the runner runs
before it writes a receipt and compare_cohorts.py runs again on the outputs a receipt names,
accepts a real passing run's outputs and rejects each edit that makes them record a failure.
"""
import contextlib,io,json,sys,tempfile,unittest
from pathlib import Path
LANE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(LANE)); sys.path.insert(0,str(LANE/"tests"))
import profile_boot as pb  # noqa: E402  (MTB4K_TEST_PROFILE, before report_contract loads)
import report_contract as vc
import run_revision8_checks as wrapper
import fixtures as fx

NEGATIVE=[]
def quiet(fn):
 with contextlib.redirect_stdout(io.StringIO()): return fn()

class WrapperTests(unittest.TestCase):
 def setUp(self):
  self.tmp=tempfile.TemporaryDirectory(); self.rules=fx.make_rules(); self.root=fx.build_cohort(Path(self.tmp.name)/"prefix","prefix",self.rules)
  self.did=vc.DRUG_IDS["amikacin"]
 def tearDown(self): self.tmp.cleanup()
 def args(self): return [str(self.root/"audit.json"),"--root",str(self.root),"--run-manifest",str(self.root/"validation-manifest.json"),
  "--expected-source-commit",fx.SOURCE,"--binary",str(self.root/"synthetic-binary"),"--backup",str(self.root/"backup"),
  "--log",str(self.root/"pipeline.log"),"--launch","2000-01-01T00:00:00Z","--foldgen",fx.foldgen_path(),
  "--expected-foldgen-sha256",fx.foldgen_sha256(),"--checker-out",str(self.root/"checker.json"),"--summary-out",str(self.root/"summary.json"),"--receipt",str(self.root/"receipt.json")]
 def run_wrapper(self,args=None,factory=None): return quiet(lambda:wrapper.main(pb.argv(args or self.args()),rules=self.rules,foldgen_factory=factory))
 def reject_report_mutation(self,label,mutate):
  root=fx.build_cohort(Path(self.tmp.name)/(label.replace(" ","-").replace("/","-")),"prefix",self.rules)
  fx.mutate_report(root,self.did,mutate)
  old=self.root; self.root=root
  try: code=self.run_wrapper()
  finally: self.root=old
  self.assertEqual(code,2); self.assertFalse((root/"checker.json").exists()); self.assertFalse((root/"receipt.json").exists())
  NEGATIVE.append({"control":label,"result":"wrapper preflight rejected before checker"})
 def test_complete_wrapper_and_summarizer_pass(self):
  code=self.run_wrapper(); self.assertEqual(code,0)
  audit=json.loads((self.root/"audit.json").read_text()); self.assertTrue(audit["summary"]["revision8_wrapper"]["valid"]); self.assertEqual(audit["summary"]["by_revision"],{str(vc.REVISION):22})
  summary=json.loads((self.root/"summary.json").read_text()); self.assertEqual(summary["coverage"][str(vc.REVISION)]["files"],22); self.assertEqual(summary["coverage"]["8"]["issues"],0)
  receipt=json.loads((self.root/"receipt.json").read_text())
  self.assertEqual((receipt["schema"],receipt["cohort"],receipt["valid"],receipt["errors"],receipt["reports"]),("grm-mtb-4k-revision8-check-receipt-v1","prefix",True,[],22))
  self.assertEqual(receipt["run_manifest_sha256"],vc.file_sha256(self.root/"validation-manifest.json")); self.assertEqual(receipt["pins"],vc.PINS)
  self.assertEqual(receipt["foldgen_sha256"],vc.file_sha256(fx.foldgen_path()))
  for key,entry in receipt["outputs"].items(): self.assertEqual(vc.file_sha256(entry["path"]),entry["sha256"],key)
 def test_foldgen_rebuild_trap_fails_postflight(self):
  did=vc.DRUG_IDS["amikacin"]; path=self.root/"data/species"/vc.SPECIES/"results"/(did+".json"); report=json.loads(path.read_text()); split=report["splits"][0]
  i=next(i for i,f in enumerate(split["test_fold"]) if f==0 and report["samples"][i]["resistant"]); j=next(j for j,f in enumerate(split["test_fold"]) if f==1 and report["samples"][j]["resistant"])
  split["test_fold"][i],split["test_fold"][j]=split["test_fold"][j],split["test_fold"][i]; fx.rescore(report); fx.write_json(path,report); fx.refresh_root(self.root,"prefix")
  code=self.run_wrapper(); self.assertNotEqual(code,0)
  audit=json.loads((self.root/"audit.json").read_text()); self.assertFalse(audit["summary"]["revision8_wrapper"]["valid"]); self.assertGreater(audit["summary"]["issue_count"],0)
  self.assertTrue(any("test_fold" in json.dumps(rec.get("issues")) for rec in audit["records"])); self.assertFalse((self.root/"receipt.json").exists())
  NEGATIVE.append({"control":"foldgen rebuilt-fold mismatch (real foldgen)","result":"audit postflight rejected"})
 def test_unknown_status_stops_before_checker(self):
  did=vc.DRUG_IDS["amikacin"]; path=self.root/"data/species"/vc.SPECIES/"results"/(did+".json"); report=json.loads(path.read_text())
  report["models"]["lasso_fit"]["selected_convergence"]["unknown"]="done"; fx.write_json(path,report); fx.refresh_root(self.root,"prefix")
  code=self.run_wrapper(); self.assertEqual(code,2); self.assertFalse((self.root/"checker.json").exists())
  NEGATIVE.append({"control":"wrapper unknown final status","result":"stopped before checker"})
 def test_final_impossible_statuses_stop_before_checker(self):
  self.reject_report_mutation("wrapper final selected impossible status and excess path",lambda r:r["models"]["lasso_fit"].update(
   path_points=31,selected_convergence={"irls_steps":1,"sweep_cap_hits":2,"converged":True}))  # CHANGED ASSERTION (R1): was {1,1,true}
  # CHANGED ASSERTION (IRLS0): the refit status was {0,0,true}, which 690815a writes on a valid warm start.
  self.reject_report_mutation("wrapper final refit zero-step convergence with a cap hit",lambda r:r["models"]["lasso_fit"].update(
   refit_convergence={"irls_steps":0,"sweep_cap_hits":1,"converged":True}))
 def test_per_fold_impossible_statuses_stop_before_checker(self):
  self.reject_report_mutation("wrapper fold selected impossible status and excess path",lambda r:r["splits"][0]["learners"][4]["lasso_fits"][0].update(
   path_points=31,selected_convergence={"irls_steps":1,"sweep_cap_hits":2,"converged":True}))
  # CHANGED ASSERTION (R1): was "wrapper fold refit all-steps-capped convergence", {1,1,true}, which 690815a can write.
  self.reject_report_mutation("wrapper fold refit cap hits exceed steps",lambda r:r["splits"][0]["learners"][4]["lasso_fits"][0].update(
   refit_convergence={"irls_steps":2,"sweep_cap_hits":3,"converged":True}))
 def test_valid_intermediate_cap_passes_complete_wrapper(self):
  def mutate(r):
   fits=[r["models"]["lasso_fit"]]
   for split in r["splits"]: fits.extend(split["learners"][4]["lasso_fits"])
   for fit in fits:
    for key in ("selected_convergence","refit_convergence"):
     if fit[key] is not None: fit[key]={"irls_steps":2,"sweep_cap_hits":1,"converged":True}
  fx.mutate_report(self.root,self.did,mutate)
  code=self.run_wrapper(); self.assertEqual(code,0)
 def test_all_steps_capped_converged_passes_complete_wrapper(self):
  # R1: {1,1,true} and {3,3,true} on every fit pass the contract, the pinned checker and the audit.
  def mutate(r):
   fits=[r["models"]["lasso_fit"]]
   for split in r["splits"]: fits.extend(split["learners"][4]["lasso_fits"])
   for i,fit in enumerate(fits):
    for key in ("selected_convergence","refit_convergence"):
     if fit[key] is not None: fit[key]={"irls_steps":1+2*(i%2),"sweep_cap_hits":1+2*(i%2),"converged":True}
  fx.mutate_report(self.root,self.did,mutate)
  code=self.run_wrapper(); self.assertEqual(code,0)
  self.assertTrue(json.loads((self.root/"receipt.json").read_text())["valid"])
 def test_added_14_reports_and_eight_skips_pass(self):
  root=fx.build_cohort(Path(self.tmp.name)/"added","added",self.rules); old=self.root; self.root=root
  try: code=self.run_wrapper()
  finally: self.root=old
  self.assertEqual(code,0); audit=json.loads((root/"audit.json").read_text())
  self.assertEqual(audit["summary"]["files"],14); self.assertEqual(len(audit["summary"]["class_minimum_skipped"]),8)
  self.assertEqual(audit["summary"]["by_revision"],{str(vc.REVISION):14})
  self.assertEqual(json.loads((root/"receipt.json").read_text())["class_gate_skips"],8)
 def test_future_launch_stale_report_rejected(self):
  args=self.args(); args[args.index("--launch")+1]="2999-01-01T00:00:00Z"; code=self.run_wrapper(args)
  self.assertEqual(code,2); audit=json.loads((self.root/"audit.json").read_text()); self.assertFalse(audit["summary"]["revision8_wrapper"]["valid"])
  self.assertFalse((self.root/"receipt.json").exists())
  NEGATIVE.append({"control":"stale pre-launch report timestamp","result":"audit postflight rejected"})

class FollowupTests(unittest.TestCase):
 setUp=WrapperTests.setUp; tearDown=WrapperTests.tearDown; args=WrapperTests.args; run_wrapper=WrapperTests.run_wrapper
 def test_i4_a_foldgen_that_is_not_the_recorded_build_stops_before_the_checker(self):
  args=self.args(); args[args.index("--expected-foldgen-sha256")+1]="0"*64
  self.assertEqual(self.run_wrapper(args),2); self.assertFalse((self.root/"checker.json").exists()); self.assertFalse((self.root/"receipt.json").exists())
  args=self.args(); args[args.index("--expected-foldgen-sha256")+1]="not-a-digest"
  self.assertEqual(self.run_wrapper(args),2); self.assertFalse((self.root/"checker.json").exists())
  args=self.args(); i=args.index("--expected-foldgen-sha256"); del args[i:i+2]
  with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit) as stop: wrapper.main(pb.argv(args),rules=self.rules)
  self.assertEqual(stop.exception.code,2)
  NEGATIVE.append({"control":"I4 foldgen is not the recorded build","result":"stopped before checker"})
 def test_l2_verify_outputs_accepts_a_real_pass_and_rejects_each_failure(self):
  self.assertEqual(self.run_wrapper(),0)
  cohort=vc.validate_manifest(self.root/"validation-manifest.json",self.rules)
  real={"checker":self.root/"checker.json","audit":self.root/"audit.json","summary":self.root/"summary.json"}
  binary=vc.file_sha256(self.root/"synthetic-binary"); fg=fx.foldgen_sha256()
  self.assertEqual(wrapper.verify_outputs(real,cohort,fx.SOURCE,binary,fg),[])
  self.assertTrue(wrapper.verify_outputs(real,cohort,"2"*40,binary,fg)); self.assertTrue(wrapper.verify_outputs(real,cohort,fx.SOURCE,"0"*64,fg))
  self.assertTrue(wrapper.verify_outputs(real,cohort,fx.SOURCE,binary,"0"*64))
  def audit_edit(fn):
   def apply(path):
    value=json.loads(path.read_text()); fn(value); path.write_text(json.dumps(value))
   return apply
  cases=[
   ("audit",audit_edit(lambda a:a["summary"]["revision8_wrapper"].update(valid=False)),"valid"),
   ("audit",audit_edit(lambda a:a["summary"].update(issue_count=1)),"issues"),
   ("audit",audit_edit(lambda a:a["records"][0].update(sha256="0"*64)),"SHA-256 differs"),
   ("audit",audit_edit(lambda a:a["records"][0].update(id="other")),"ids differ"),
   ("checker",audit_edit(lambda c:c["reports"][0].update(id="other")),"checker report ids"),
   ("summary",audit_edit(lambda s:s["coverage"][str(vc.REVISION)].update(issues=1)),"summarizer"),
  ]
  for key,fn,match in cases:
   with self.subTest(match):
    copy=self.root/f"edited-{key}.json"; copy.write_bytes(real[key].read_bytes()); fn(copy)
    errors=wrapper.verify_outputs({**real,key:copy},cohort,fx.SOURCE,binary,fg)
    self.assertTrue(any(match in e for e in errors),errors)
    NEGATIVE.append({"control":f"L2 verify_outputs: {match}","result":"rejected"})

class RealFoldgenTests(unittest.TestCase):
 """Experiment E2: with folds from the real foldgen the wrapper passes, a swap of two same-class
 random folds fails the rebuild, and K's index folds (the double) fail it everywhere."""
 def test_matched2k_and_full_pass_and_the_index_fold_double_is_caught(self):
  with tempfile.TemporaryDirectory() as tmp:
   rules=fx.make_rules()
   for name in ("full","matched2k"):
    root=fx.build_cohort(Path(tmp)/name,name,rules)
    args=[str(root/"audit.json"),"--root",str(root),"--run-manifest",str(root/"validation-manifest.json"),"--expected-source-commit",fx.SOURCE,
          "--binary",str(root/"synthetic-binary"),"--backup",str(root/"backup"),"--log",str(root/"pipeline.log"),"--launch","2000-01-01T00:00:00Z",
          "--foldgen",fx.foldgen_path(),"--expected-foldgen-sha256",fx.foldgen_sha256(),
          "--checker-out",str(root/"checker.json"),"--summary-out",str(root/"summary.json"),"--receipt",str(root/"receipt.json")]
    self.assertEqual(quiet(lambda:wrapper.main(pb.argv(args),rules=rules)),0,name)
   root=fx.build_cohort(Path(tmp)/"double","prefix",rules)
   args=[str(root/"audit.json"),"--root",str(root),"--run-manifest",str(root/"validation-manifest.json"),"--expected-source-commit",fx.SOURCE,
         "--binary",str(root/"synthetic-binary"),"--backup",str(root/"backup"),"--log",str(root/"pipeline.log"),"--launch","2000-01-01T00:00:00Z",
         "--foldgen",fx.foldgen_path(),"--expected-foldgen-sha256",fx.foldgen_sha256(),
         "--checker-out",str(root/"checker.json"),"--summary-out",str(root/"summary.json"),"--receipt",str(root/"receipt.json")]
   self.assertEqual(quiet(lambda:wrapper.main(pb.argv(args),rules=rules,foldgen_factory=fx.FakeFoldGen)),2)
   audit=json.loads((root/"audit.json").read_text()); self.assertGreater(audit["summary"]["issue_count"],0)
   NEGATIVE.append({"control":"index-fold double against real-foldgen folds","result":"audit postflight rejected"})

if __name__=="__main__":
 result=unittest.main(exit=False,verbosity=2).result
 print(json.dumps({"negative_controls":len(NEGATIVE)}))
 sys.exit(0 if result.wasSuccessful() else 1)
