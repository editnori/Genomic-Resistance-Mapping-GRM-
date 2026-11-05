#!/usr/bin/env python3
"""report_contract.py: package K's guard, solver-status and known-answer tests, plus the
reconciliation's F1 rule, its 38-case negative-control matrix, the fingerprint binding and the
checker pins.

usage: FOLDGEN=<foldgen binary> python3 -I scripts/mtb-4k/tests/test_report_contract.py

Ported from K's tests/test_validators.py (GuardTests, SolverStatusTests, KnownAnswerTests) with
paths adapted; every K assertion is kept. The fixtures draw folds from the real foldgen.

IRLS0 (marked): two solver-status controls changed, not dropped. SolverStatusTests' zero-step
case is now 0 steps with 1 cap hit (still rejected), and 0 steps with none, converged, is
accepted; the E4 entry is now "convergence irls_steps=0 with sweep_cap_hits 1" (still rejected,
38 in all). 690815a writes irls_steps 0 when the warm start already meets the KKT tolerance
(src/learn/lasso.rs:70), the pinned check-results.py:114 accepts any non-negative step count,
and the K review treated that status as impossible because it was at 855bdba.

R1 (marked): the contract no longer rejects a converged status whose every IRLS step hit the cap.
SolverStatusTests' "every solver step capped" ({3,3,true}) became "cap hits exceed steps after
two steps" ({2,3,true}, still rejected), and a new test accepts {1,1,true}, {3,3,true} and
{3,1,true} for selected and refit, in the contract and in the pinned checker. The E4 entry
"every IRLS step capped yet converged" ({1,1,true} on the refit) became "refit sweep_cap_hits >
irls_steps" ({3,4,true}, still rejected, 38 in all). Under 690815a a step whose inner solve hits
its iteration guard counts as a cap hit and the point can still converge
(docs/lasso-stopping.md:39), which check-results.py accepts.
"""
import copy,hashlib,json,math,os,shutil,subprocess,sys,tempfile,unittest
from pathlib import Path
LANE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(LANE)); sys.path.insert(0,str(LANE/"tests"))
import profile_boot as pb  # noqa: E402  (MTB4K_TEST_PROFILE, before report_contract loads)
import report_contract as vc
import fixtures as fx

NEGATIVE=[]
def reject(case,label,match,fn):
 try: fn()
 except vc.ValidationError as exc:
  case.assertIn(match,str(exc)); NEGATIVE.append({"control":label,"matched":match,"result":"rejected"}); return
 case.fail("guard accepted negative control: "+label)

class GuardTests(unittest.TestCase):
 def setUp(self):
  self.tmp=tempfile.TemporaryDirectory(); self.rules=fx.make_rules({"prefix":("amikacin",)})
  self.root=fx.build_cohort(Path(self.tmp.name)/"prefix","prefix",self.rules); self.did=vc.DRUG_IDS["amikacin"]
 def tearDown(self): self.tmp.cleanup()
 def validate(self): return vc.validate_manifest(self.root/"validation-manifest.json",self.rules)
 def rehash(self):
  path=self.root/"validation-manifest.json"; manifest=json.loads(path.read_text()); entry=manifest["reports"][self.did]
  entry["sha256"]=fx.sha(self.root/entry["path"]); fx.write_json(path,manifest)
 def change(self,fn,rescore=False):
  path=self.root/"data/species"/vc.SPECIES/"results"/(self.did+".json"); report=json.loads(path.read_text()); fn(report)
  if rescore: fx.rescore(report)
  fx.write_json(path,report); self.rehash(); return report
 def raw(self,fn,rehash=True):
  path=self.root/"data/species"/vc.SPECIES/"results"/(self.did+".json"); path.write_text(fn(path.read_text()))
  if rehash: self.rehash()
 def test_historical_revision7_rejected(self):
  self.change(lambda r:r["evaluation"].update(revision=7)); reject(self,"historical revision 7",f"revision {vc.REVISION}",self.validate)
 def test_wrong_schema_rejected(self):
  self.change(lambda r:r.update(schema=3)); reject(self,"wrong report schema","schema 2",self.validate)
 def test_invalid_evaluation_rejected(self):
  self.change(lambda r:r["evaluation"].update(valid=False,reasons=["synthetic invalid"])); reject(self,"invalid evaluation","valid=true",self.validate)
 def test_unknown_critical_field_rejected(self):
  self.change(lambda r:r.update(unknown_revision8_field=1)); reject(self,"unknown report field","unknown",self.validate)
 def test_unknown_lasso_status_field_rejected(self):
  def mutate(r): r["splits"][0]["learners"][4]["lasso_fits"][0]["selected_convergence"]["vendor_status"]="finished"
  self.change(mutate); reject(self,"unknown convergence status field","field mismatch",self.validate)
 def test_missing_lasso_status_rejected(self):
  self.change(lambda r:r["splits"][0]["learners"][4].pop("lasso_fits")); reject(self,"missing fold lasso status","lasso_fits",self.validate)
 def test_unconverged_path_rejected(self):
  self.change(lambda r:r["splits"][0]["learners"][4]["lasso_fits"][0].update(path_unconverged=1)); reject(self,"unconverged path","unconverged",self.validate)
 def test_unconverged_selected_refit_rejected(self):
  self.change(lambda r:r["splits"][0]["learners"][4]["lasso_fits"][0]["refit_convergence"].update(converged=False)); reject(self,"unconverged refit","unconverged",self.validate)
 def test_missing_final_status_rejected(self):
  self.change(lambda r:r["models"].pop("lasso_fit")); reject(self,"missing final lasso status","models.lasso_fit",self.validate)
 def test_empty_path_score_mismatch_rejected(self):
  self.change(lambda r:r["splits"][0]["learners"][4]["lasso_fits"][0].update(fx.fit(empty=True))); reject(self,"empty-path/no-call mismatch","empty-path",self.validate)
 def test_corrupted_score_rejected(self):
  def mutate(r):
   row=r["splits"][0]["learners"][0]; row["scores"][0]=1.0-row["scores"][0]
  self.change(mutate); reject(self,"corrupted held-out score","recomputed",self.validate)
 def test_non_probability_score_rejected(self):
  self.change(lambda r:r["splits"][0]["learners"][1]["scores"].__setitem__(0,1.25)); reject(self,"score outside probability range","outside [0,1]",self.validate)
 def test_corrupted_metric_rejected(self):
  self.change(lambda r:r["splits"][0]["learners"][1]["metrics"].update(balanced_accuracy=0.123)); reject(self,"corrupted metric","balanced_accuracy",self.validate)
 def test_corrupted_fold_ba_rejected(self):
  self.change(lambda r:r["splits"][0]["learners"][1]["fold_balanced_accuracy"].__setitem__(0,0.123)); reject(self,"corrupted fold metric","fold_balanced_accuracy",self.validate)
 def test_no_call_count_rejected(self):
  self.change(lambda r:r["splits"][1]["learners"][3]["no_call"].update(resistant=0)); reject(self,"wrong no-call denominator","no_call",self.validate)
 def test_partial_fold_no_call_rejected(self):
  def mutate(r):
   split=r["splits"][1]; row=split["learners"][3]; i=next(i for i,f in enumerate(split["test_fold"]) if f==0); row["scores"][i]=None
  self.change(mutate,rescore=True); reject(self,"partial-fold no-call","called in part",self.validate)
 def test_fingerprint_mismatch_rejected(self):
  path=self.root/"validation-manifest.json"; manifest=json.loads(path.read_text()); manifest["reports"][self.did]["input_fingerprint"]="0"*16; fx.write_json(path,manifest)
  reject(self,"stale input fingerprint","fingerprint",self.validate)
 def test_stale_input_hash_rejected(self):
  path=self.root/"data/species"/vc.SPECIES/"counted.json"; counted=json.loads(path.read_text()); counted["genomes"][0]["n_kmers"]=1; fx.write_json(path,counted)
  reject(self,"stale counted input","stale input SHA-256",self.validate)
 def test_report_hash_mismatch_rejected(self):
  self.raw(lambda s:s.replace('"seconds":10.0','"seconds":11.0',1),rehash=False); reject(self,"report hash binding","SHA-256",self.validate)
 def test_unknown_split_status_rejected(self):
  self.change(lambda r:r["splits"][1].update(scheme="temporal")); reject(self,"unknown split scheme","unknown or duplicate",self.validate)
 def test_random_unavailable_rejected(self):
  def mutate(r):
   split=r["splits"][0]; split["learners"]=[]; split.pop("test_fold"); split["unavailable_reason"]="synthetic"
  self.change(mutate); reject(self,"unavailable random split","random split cannot be unavailable",self.validate)
 def test_malformed_json_rejected(self):
  self.raw(lambda s:s[:100]); reject(self,"malformed JSON","malformed JSON",self.validate)
 def test_nan_rejected(self):
  self.raw(lambda s:s.replace('"seconds":10.0','"seconds":NaN',1)); reject(self,"NaN JSON","nonfinite",self.validate)
 def test_overflow_nonfinite_rejected(self):
  self.raw(lambda s:s.replace('"seconds":10.0','"seconds":1e999',1)); reject(self,"overflow nonfinite JSON","nonfinite",self.validate)
 def test_duplicate_json_field_rejected(self):
  self.raw(lambda s:s.replace('"schema":2','"schema":2,"schema":2',1)); reject(self,"duplicate JSON field","duplicate JSON field",self.validate)
 def test_duplicate_genome_rejected(self):
  def mutate(r):
   r["samples"][1]=copy.deepcopy(r["samples"][0])
   for split in r["splits"]:
    split["test_fold"][1]=split["test_fold"][0]
    for row in split["learners"]: row["scores"][1]=row["scores"][0]
  self.change(mutate,rescore=True); reject(self,"duplicate genome ID","duplicate sample",self.validate)
 def test_fold_out_of_range_rejected(self):
  self.change(lambda r:r["splits"][0].__setitem__("test_fold",[99]+r["splits"][0]["test_fold"][1:])); reject(self,"fold outside range","outside configured range",self.validate)
 def test_lineage_spanning_fold_rejected(self):
  def mutate(r):
   split=r["splits"][1]; split["test_fold"][0]=(split["test_fold"][0]+1)%5
  self.change(mutate,rescore=True); reject(self,"lineage spans folds","lineage spans",self.validate)
 def test_catalog_arithmetic_rejected(self):
  self.change(lambda r:r["catalog"]["metrics"].update(tp=r["catalog"]["metrics"]["tp"]+1)); reject(self,"catalog metric arithmetic","catalog.metrics",self.validate)
 def test_probe_arithmetic_rejected(self):
  def mutate(r):
   m=vc.upstream_audit.from_counts(1,2,3,4); m["auroc"]=None; m["auprc"]=None
   r["panels"].append({"id":"xpert_rif","name":"synthetic","mirrors":"synthetic","targets":[],"metrics":m,"detected":[],
                       "no_call":{"resistant":0,"susceptible":0}})
  self.change(mutate); reject(self,"probe no-call partition","partition",self.validate)
 def test_wrong_frozen_cohort_hash_rejected(self):
  reject(self,"frozen cohort hash binding","frozen rule",lambda:vc.validate_manifest(self.root/"validation-manifest.json",vc.RULES))
 def test_float_confusion_count_rejected(self):
  self.change(lambda r:r["splits"][0]["learners"][1]["metrics"].update(tp=float(r["splits"][0]["learners"][1]["metrics"]["tp"])))
  reject(self,"non-integer confusion count","recomputed integer",self.validate)
 def test_duplicate_panel_rejected(self):
  def mutate(r):
   m=vc.upstream_audit.from_counts(10,5,20,25); m["auroc"]=None; m["auprc"]=None
   panel={"id":"synthetic_catalog_panel","name":"synthetic","mirrors":"synthetic","targets":[],"metrics":m,"detected":[]}
   r["panels"]=[panel,copy.deepcopy(panel)]
  self.change(mutate); reject(self,"duplicate panel ID","duplicate panel",self.validate)
 def test_negative_complexity_rejected(self):
  self.change(lambda r:r["splits"][0]["learners"][1]["fold_complexity"].__setitem__(0,-1.0)); reject(self,"negative fold complexity","negative model complexity",self.validate)
 def test_float_schema_rejected(self):
  self.change(lambda r:r.update(schema=2.0)); reject(self,"non-integer schema","schema 2",self.validate)

class SolverStatusTests(unittest.TestCase):
 def test_valid_intermediate_cap_and_short_or_empty_paths(self):
  fit=fx.fit(); fit["path_points"]=7
  for key in ("selected_convergence","refit_convergence"):
   fit[key]={"irls_steps":2,"sweep_cap_hits":1,"converged":True}
  vc.validate_lasso_fit(fit,"valid intermediate cap",30)
  vc.validate_lasso_fit(fx.fit(empty=True),"valid empty path",30)
 def test_impossible_convergence_rejected_for_selected_and_refit(self):
  cases=[("zero solver steps with a cap hit",0,1,"exceeds irls_steps"),  # CHANGED ASSERTION (IRLS0): was 0 steps, 0 hits
         ("cap hits exceed steps",1,2,"exceeds irls_steps"),
         ("cap hits exceed steps after two steps",2,3,"exceeds irls_steps")]  # CHANGED ASSERTION (R1): was 3 steps, 3 hits
  for label,steps,hits,match in cases:
   for key in ("selected_convergence","refit_convergence"):
    with self.subTest(label=label,status=key):
     fit=fx.fit(); fit[key]={"irls_steps":steps,"sweep_cap_hits":hits,"converged":True}
     reject(self,label+"/"+key,match,lambda fit=fit:vc.validate_lasso_fit(fit,"synthetic",30))
 def test_capped_steps_that_converge_are_accepted_for_selected_and_refit(self):
  # R1: under 690815a a capped step can still end in a converged point, and check-results.py accepts it.
  for steps,hits in ((1,1),(3,3),(3,1)):
   for key in ("selected_convergence","refit_convergence"):
    with self.subTest(steps=steps,hits=hits,status=key):
     fit=fx.fit(); fit[key]={"irls_steps":steps,"sweep_cap_hits":hits,"converged":True}
     self.assertTrue(vc.validate_lasso_fit(fit,"synthetic",30))
     vc.upstream_checker.check_lasso_fit(copy.deepcopy(fit),("synthetic",))
 def test_zero_step_warm_start_accepted_for_selected_and_refit(self):
  # IRLS0: 690815a takes 0 IRLS steps where the warm start already meets the KKT tolerance.
  for steps in (0,3):
   for key in ("selected_convergence","refit_convergence"):
    with self.subTest(steps=steps,status=key):
     fit=fx.fit(); fit[key]={"irls_steps":steps,"sweep_cap_hits":0,"converged":True}
     self.assertTrue(vc.validate_lasso_fit(fit,"synthetic",30))
 def test_unconverged_zero_step_status_is_left_to_the_checker(self):
  # IRLS0: the contract checks a status's keys and bounds; the pinned checker judges convergence.
  status={"irls_steps":0,"sweep_cap_hits":0,"converged":False}
  vc.validate_convergence(status,"synthetic")
  fit=fx.fit(); fit["refit_convergence"]=status
  with self.assertRaises(AssertionError) as raised: vc.upstream_checker.check_lasso_fit(copy.deepcopy(fit),("synthetic",))
  self.assertIn("unconverged two-class lasso fit",raised.exception.args[0])
  reject(self,"unconverged zero-step refit","invalid revision-8 lasso status",lambda:vc.validate_lasso_fit(fit,"synthetic",30))
 def test_excess_path_points_rejected(self):
  fit=fx.fit(); fit["path_points"]=31
  reject(self,"path exceeds configured lambdas","path_points exceeds eval_config.lasso_lambdas",
         lambda:vc.validate_lasso_fit(fit,"synthetic",30))

class KnownAnswerTests(unittest.TestCase):
 def test_clopper_pearson_closed_forms(self):
  x=vc.clopper_pearson(0,10); self.assertEqual(x["lower"],0.0); self.assertAlmostEqual(x["upper"],1-0.025**0.1,places=14)
  x=vc.clopper_pearson(10,10); self.assertEqual(x["upper"],1.0); self.assertAlmostEqual(x["lower"],0.025**0.1,places=14)
  x=vc.clopper_pearson(1,2); self.assertAlmostEqual(x["lower"],1-math.sqrt(0.975),places=14); self.assertAlmostEqual(x["upper"],math.sqrt(0.975),places=14)
  x=vc.clopper_pearson(0,0); self.assertEqual((x["lower"],x["upper"],x["status"]),(None,None,"undefined_no_called_denominator"))
  for n in (1,2,10,25,4000):
   self.assertAlmostEqual(vc.clopper_pearson(0,n)["upper"],1-0.025**(1/n),places=14)
   self.assertAlmostEqual(vc.clopper_pearson(n,n)["lower"],0.025**(1/n),places=14)
 def test_clopper_pearson_exact_binomial_tail_definitions(self):
  def cumulative(n,p,lo,hi): return sum(math.comb(n,k)*p**k*(1-p)**(n-k) for k in range(lo,hi+1))
  for n,x in ((2,1),(10,3),(10,0),(10,10),(20,7)):
   interval=vc.clopper_pearson(x,n)
   if x>0: self.assertAlmostEqual(cumulative(n,interval["lower"],x,n),0.025,places=12)
   if x<n: self.assertAlmostEqual(cumulative(n,interval["upper"],0,x),0.025,places=12)
 def test_extreme_beta_quantile_complements(self):
  p=1-1e-12; q=1-p
  self.assertAlmostEqual(vc.beta_ppf(p,1,5),1-q**0.2,delta=2e-12)
  self.assertAlmostEqual(vc.clopper_pearson(0,5,alpha=2e-12)["upper"],1-(1e-12)**0.2,delta=2e-12)
 def test_known_confusion_no_call_and_fold_denominators(self):
  with tempfile.TemporaryDirectory() as tmp:
   rules=fx.make_rules({"prefix":("amikacin",)}); root=fx.build_cohort(Path(tmp)/"prefix","prefix",rules); c=vc.validate_manifest(root/"validation-manifest.json",rules)
   report=c.reports[vc.DRUG_IDS["amikacin"]]; row=report.splits["lineage"].learners["gbt"].receipt
   self.assertEqual(row["metrics"]["n"],48); self.assertEqual(row["uncalled"],12)
   self.assertEqual(row["no_call"],{"resistant":6,"susceptible":6})
   self.assertEqual(row["metrics"]["tp"]+row["metrics"]["fn"]+row["no_call"]["resistant"],30)
   self.assertEqual(row["metrics"]["tn"]+row["metrics"]["fp"]+row["no_call"]["susceptible"],30)
   null_fold=row["folds"][2]; self.assertEqual(null_fold["called"],0); self.assertEqual(null_fold["held_out"],12); self.assertIsNone(null_fold["metrics"])
   self.assertEqual(row["clopper_pearson_95"]["vme"]["denominator"],row["metrics"]["tp"]+row["metrics"]["fn"])
   # Both orientations from one computation: sensitivity's interval mirrors VME's.
   vme,sens=row["clopper_pearson_95"]["vme"],row["clopper_pearson_95"]["sensitivity"]
   self.assertEqual((sens["lower"],sens["upper"]),(1-vme["upper"],1-vme["lower"]))
 def test_null_metrics_are_not_zero_or_dropped_from_five_learner_set(self):
  with tempfile.TemporaryDirectory() as tmp:
   rules=fx.make_rules({"prefix":("amikacin",)}); root=fx.build_cohort(Path(tmp)/"prefix","prefix",rules); did=vc.DRUG_IDS["amikacin"]
   path=root/"data/species"/vc.SPECIES/"results"/(did+".json"); report=json.loads(path.read_text()); split=report["splits"][0]; row=split["learners"][4]
   row["scores"]=[None]*60; row["fold_complexity"]=[0.0]*5; row["lasso_fits"]=[fx.fit(empty=True) for _ in range(5)]; fx.rescore(report); fx.write_json(path,report)
   manifest=json.loads((root/"validation-manifest.json").read_text()); manifest["reports"][did]["sha256"]=fx.sha(path); fx.write_json(root/"validation-manifest.json",manifest)
   c=vc.validate_manifest(root/"validation-manifest.json",rules); summary=c.reports[did].splits["random"].receipt
   lasso=summary["learners"][4]; self.assertIsNone(lasso["metrics"]); self.assertEqual(lasso["uncalled"],60)
   median=summary["five_learner_medians"]["balanced_accuracy"]; self.assertIsNone(median["median"]); self.assertEqual(median["nonnull_learners"],4); self.assertEqual(median["null_learners"],["lasso"])
 def test_unavailable_lineage_stays_unavailable(self):
  with tempfile.TemporaryDirectory() as tmp:
   rules=fx.make_rules({"prefix":("amikacin",)}); root=fx.build_cohort(Path(tmp)/"prefix","prefix",rules); did=vc.DRUG_IDS["amikacin"]
   path=root/"data/species"/vc.SPECIES/"results"/(did+".json"); report=json.loads(path.read_text()); split=report["splits"][1]
   split["learners"]=[]; split.pop("test_fold"); split["unavailable_reason"]="synthetic unavailable"; fx.write_json(path,report)
   manifest=json.loads((root/"validation-manifest.json").read_text()); manifest["reports"][did]["sha256"]=fx.sha(path); fx.write_json(root/"validation-manifest.json",manifest)
   c=vc.validate_manifest(root/"validation-manifest.json",rules); lineage=c.reports[did].splits["lineage"]
   self.assertFalse(lineage.available); self.assertEqual(lineage.receipt["learners"],[]); self.assertIsNone(lineage.receipt["five_learner_medians"])

def null_fold(r,complexity=0.0,fold=0):
 """E3: random-split lasso fold `fold` uncalled, with `complexity`, its path record kept."""
 split=r["splits"][0]; row=split["learners"][4]
 for i,f in enumerate(split["test_fold"]):
  if f==fold: row["scores"][i]=None
 row["fold_complexity"][fold]=complexity

class F1LassoFoldRuleTests(unittest.TestCase):
 """RECONCILE F1 and experiment E3: K required path_points == 0 exactly when a lasso fold has no
 calls, which rejected a valid report; the checker's rule (855bdba) is the contract now."""
 def setUp(self):
  self.tmp=tempfile.TemporaryDirectory(); self.rules=fx.make_rules({"prefix":("amikacin",)}); self.did=vc.DRUG_IDS["amikacin"]
 def tearDown(self): self.tmp.cleanup()
 def variant(self,label,fn):
  root=fx.build_cohort(Path(self.tmp.name)/label,"prefix",self.rules); fx.mutate_report(root,self.did,fn,rescore_after=True)
  # The checker runs on a copy whose curated list holds this one dataset (as the reconciliation's
  # scripts/06c_checker_isolated.py), since the single-drug fixture has no other results.
  data=Path(self.tmp.name)/(label+"-checker-data"); shutil.copytree(root/"data",data)
  curated=json.loads((data/"curated/datasets.json").read_text()); fx.write_json(data/"curated/datasets.json",[r for r in curated if r["id"]==self.did])
  checker=subprocess.run([sys.executable,"-I",str(vc.pinned("check-results.py")),str(data),"--log",str(root/"pipeline.log"),
                          "--out",str(Path(self.tmp.name)/(label+"-checker.json"))],capture_output=True,text=True)
  return root,checker.returncode
 def test_nonempty_path_whose_refit_keeps_no_weight_is_a_valid_null_fold(self):
  root,checker=self.variant("null-fold",lambda r:null_fold(r))
  self.assertEqual(checker,0)
  c=vc.validate_manifest(root/"validation-manifest.json",self.rules)
  self.assertEqual(c.reports[self.did].splits["random"].learners["lasso"].receipt["folds"][0]["called"],0)
 def test_the_three_invalid_variants_are_rejected_as_by_the_checker(self):
  for label,fn,match in [
    ("called fold with zero complexity",lambda r:r["splits"][0]["learners"][4]["fold_complexity"].__setitem__(0,0.0),"zero complexity"),
    ("empty path with called scores",lambda r:r["splits"][0]["learners"][4]["lasso_fits"].__setitem__(0,fx.fit(empty=True)),"empty-path"),
    ("uncalled fold with positive complexity",lambda r:null_fold(r,complexity=12.0),"null fit is inconsistent")]:
   with self.subTest(label):
    root,checker=self.variant(label.replace(" ","-"),fn)
    self.assertNotEqual(checker,0,label)
    reject(self,"F1 "+label,match,lambda root=root:vc.validate_manifest(root/"validation-manifest.json",self.rules))

L=lambda r,s,i:r["splits"][s]["learners"][i]  # noqa: E731
LF=lambda r:L(r,0,4)["lasso_fits"][0]  # noqa: E731
def _partial(r):
 split=r["splits"][1]; row=split["learners"][3]; i=next(i for i,f in enumerate(split["test_fold"]) if f==0); row["scores"][i]=None
def _dup(r):
 r["samples"][1]=copy.deepcopy(r["samples"][0])
 for split in r["splits"]:
  split["test_fold"][1]=split["test_fold"][0]
  for row in split["learners"]: row["scores"][1]=row["scores"][0]
def _nan_score(s):
 rep=json.loads(s); rep["splits"][0]["learners"][2]["scores"][0]=float("nan")
 return json.dumps(rep,separators=(",",":"),allow_nan=True)+"\n"
# The reconciliation's E4 matrix (scripts/06_cross.py e4_mutation_matrix): K's guard mutations and
# five raw-JSON cases. K rejected all 38, W 12, the checker at 2d30bba 21.
# IRLS0: one entry's status changed, not dropped (the module docstring says why).
MUTATIONS={
 "historical revision 7":(lambda r:r["evaluation"].update(revision=7),False),
 "wrong report schema 3":(lambda r:r.update(schema=3),False),
 "evaluation.valid=false":(lambda r:r["evaluation"].update(valid=False,reasons=["synthetic invalid"]),False),
 "unknown top-level field":(lambda r:r.update(unknown_revision8_field=1),False),
 "unknown convergence status field":(lambda r:LF(r)["selected_convergence"].update(vendor_status="finished"),False),
 "missing per-fold lasso_fits":(lambda r:L(r,0,4).pop("lasso_fits"),False),
 "path_unconverged=1":(lambda r:LF(r).update(path_unconverged=1),False),
 "refit converged=false":(lambda r:LF(r)["refit_convergence"].update(converged=False),False),
 "missing models.lasso_fit":(lambda r:r["models"].pop("lasso_fit"),False),
 "empty path with called fold":(lambda r:L(r,0,4)["lasso_fits"].__setitem__(0,fx.fit(empty=True)),False),
 "corrupted held-out score (metrics stale)":(lambda r:L(r,0,0)["scores"].__setitem__(0,1.0-L(r,0,0)["scores"][0]),False),
 "score 1.25 outside [0,1] (rescored)":(lambda r:L(r,0,1)["scores"].__setitem__(0,1.25),True),
 "stored balanced_accuracy wrong (confusion right)":(lambda r:L(r,0,1)["metrics"].update(balanced_accuracy=0.123),False),
 "stored fold_balanced_accuracy wrong":(lambda r:L(r,0,1)["fold_balanced_accuracy"].__setitem__(0,0.123),False),
 "no_call count wrong":(lambda r:L(r,1,3)["no_call"].update(resistant=0),False),
 "partial-fold no-call (rescored)":(_partial,True),
 "duplicate genome ID (rescored)":(_dup,True),
 "test_fold out of range":(lambda r:r["splits"][0].__setitem__("test_fold",[99]+r["splits"][0]["test_fold"][1:]),False),
 "lineage spans folds (rescored)":(lambda r:r["splits"][1]["test_fold"].__setitem__(0,(r["splits"][1]["test_fold"][0]+1)%5),True),
 "catalog metric arithmetic":(lambda r:r["catalog"]["metrics"].update(tp=r["catalog"]["metrics"]["tp"]+1),False),
 "float confusion count":(lambda r:L(r,0,1)["metrics"].update(tp=float(L(r,0,1)["metrics"]["tp"])),False),
 "negative fold complexity":(lambda r:L(r,0,1)["fold_complexity"].__setitem__(0,-1.0),False),
 "float schema 2.0":(lambda r:r.update(schema=2.0),False),
 "convergence irls_steps=0 with sweep_cap_hits 1":(lambda r:LF(r).update(selected_convergence={"irls_steps":0,"sweep_cap_hits":1,"converged":True}),False),  # CHANGED ASSERTION (IRLS0): was {0,0,true}
 "sweep_cap_hits > irls_steps":(lambda r:LF(r).update(selected_convergence={"irls_steps":1,"sweep_cap_hits":2,"converged":True}),False),
 "refit sweep_cap_hits > irls_steps":(lambda r:LF(r).update(refit_convergence={"irls_steps":3,"sweep_cap_hits":4,"converged":True}),False),  # CHANGED ASSERTION (R1): was "every IRLS step capped yet converged", {1,1,true}
 "path_points 31 > lasso_lambdas 30":(lambda r:LF(r).update(path_points=31),False),
 "selected_point == path_points":(lambda r:LF(r).update(selected_point=30),False),
 "convergence object without irls_steps/sweep_cap_hits (W fixture shape)":(lambda r:LF(r).update(selected_convergence={"converged":True}),False),
 "evaluation.scope lacks non-clinical limitation":(lambda r:r["evaluation"].update(scope="Internal cross-validation."),False),
 "lasso_fits on a non-lasso row":(lambda r:L(r,0,1).update(lasso_fits=[fx.fit()]*5),False),
 "input_fingerprint malformed":(lambda r:r.update(input_fingerprint="not-hex"),False),
 "missing no_call on a row with null fold":(lambda r:L(r,1,3).pop("no_call"),False),
}
RAW={
 "NaN seconds (raw JSON)":lambda s:s.replace('"seconds":10.0','"seconds":NaN',1),
 "1e999 seconds (raw JSON)":lambda s:s.replace('"seconds":10.0','"seconds":1e999',1),
 "duplicate JSON key (raw)":lambda s:s.replace('"schema":2','"schema":2,"schema":2',1),
 "truncated JSON (raw)":lambda s:s[:100],
 "NaN held-out score (raw JSON)":_nan_score,
}

class MutationMatrixTests(unittest.TestCase):
 def test_all_38_reconciliation_corruptions_are_rejected(self):
  self.assertEqual(len(MUTATIONS)+len(RAW),38)
  rules=fx.make_rules({"prefix":("amikacin",)}); did=vc.DRUG_IDS["amikacin"]; accepted=[]
  with tempfile.TemporaryDirectory() as tmp:
   pristine=fx.build_cohort(Path(tmp)/"pristine","prefix",rules)
   vc.validate_manifest(pristine/"validation-manifest.json",rules)
   for i,label in enumerate(list(MUTATIONS)+list(RAW)):
    with self.subTest(label):
     root=Path(tmp)/f"m{i:02d}"; shutil.copytree(pristine,root)
     path=root/"data/species"/vc.SPECIES/"results"/(did+".json")
     if label in RAW: path.write_text(RAW[label](path.read_text()))
     else:
      fn,rs=MUTATIONS[label]; rep=json.loads(path.read_text()); fn(rep)
      if rs: fx.rescore(rep)
      fx.write_json(path,rep)
     manifest=json.loads((root/"validation-manifest.json").read_text()); manifest["reports"][did]["sha256"]=fx.sha(path); fx.write_json(root/"validation-manifest.json",manifest)
     try: vc.validate_manifest(root/"validation-manifest.json",rules); accepted.append(label)
     except vc.ValidationError: NEGATIVE.append({"control":"E4 "+label,"result":"rejected"})
  self.assertEqual(accepted,[])

class FingerprintBindingTests(unittest.TestCase):
 """RECONCILE section 4: recomputed = report = manifest, from the root's own files."""
 def setUp(self):
  self.tmp=tempfile.TemporaryDirectory(); self.rules=fx.make_rules({"prefix":("amikacin",)})
  self.root=fx.build_cohort(Path(self.tmp.name)/"prefix","prefix",self.rules); self.did=vc.DRUG_IDS["amikacin"]
 def tearDown(self): self.tmp.cleanup()
 def rehash_inputs(self):
  path=self.root/"validation-manifest.json"; manifest=json.loads(path.read_text())
  for entry in manifest["inputs"]: entry["sha256"]=fx.sha(self.root/entry["path"])
  fx.write_json(path,manifest)
 def test_the_fixture_reports_carry_the_recomputed_fingerprint(self):
  c=vc.validate_manifest(self.root/"validation-manifest.json",self.rules)
  self.assertTrue(c.reports[self.did].receipt["input_fingerprint_recomputed"])
 def test_inputs_changed_after_the_reports_are_a_binding_stop_even_when_rehashed(self):
  for label,edit in [
    ("a phenotype's call flipped",lambda root:(root/"data/curated/phenotypes.jsonl").write_text((root/"data/curated/phenotypes.jsonl").read_text().replace('"reported": "R"','"reported": "S"',1))),
    ("a matrix word changed",lambda root:(root/"data/species"/vc.SPECIES/"matrix/kmers.bin").write_bytes(b"\1"+(root/"data/species"/vc.SPECIES/"matrix/kmers.bin").read_bytes()[1:])),
    ("the catalog version changed",lambda root:(root/"data/catalog/version.txt").write_text("other\n"))]:
   with self.subTest(label):
    saved={p:p.read_bytes() for p in [self.root/"data/curated/phenotypes.jsonl",self.root/"data/species"/vc.SPECIES/"matrix/kmers.bin",self.root/"data/catalog/version.txt"]}
    edit(self.root)
    with self.assertRaises(vc.ValidationError): vc.validate_manifest(self.root/"validation-manifest.json",self.rules)
    self.rehash_inputs()
    with self.assertRaises(vc.BindingError) as stop: vc.validate_manifest(self.root/"validation-manifest.json",self.rules,collect_invalid=True)
    NEGATIVE.append({"control":"fingerprint binding: "+label,"result":"rejected","message":str(stop.exception)[-120:]})
    for p,b in saved.items(): p.write_bytes(b)
    self.rehash_inputs()
  vc.validate_manifest(self.root/"validation-manifest.json",self.rules)
 def test_a_class_gate_receipt_must_match_the_roots_phenotypes(self):
  rules=fx.make_rules()
  root=fx.build_cohort(Path(self.tmp.name)/"added","added",rules)
  path=root/"validation-manifest.json"; manifest=json.loads(path.read_text())
  row=next(r for r in manifest["class_gate_unavailable"] if r["antibiotic"]=="para-aminosalicylic acid")
  self.assertEqual((row["resistant"],row["susceptible"]),(5,55))
  bad=dict(rules); bad["added"]=vc.CohortRule(2000,rules["added"].genomes_sha256,vc.ADDED_DRUGS,unavailable={**fx.FIXTURE_ADDED_UNAVAILABLE,"para-aminosalicylic acid":(6,54)})
  row.update(resistant=6,susceptible=54); fx.write_json(path,manifest)
  reject(self,"class-gate receipt disagrees with the root's phenotypes","differs from the root's phenotypes",lambda:vc.validate_manifest(path,bad))

class CollectInvalidTests(unittest.TestCase):
 def test_a_contract_failure_is_collected_and_a_binding_failure_still_stops(self):
  with tempfile.TemporaryDirectory() as tmp:
   rules=fx.make_rules({"prefix":("amikacin",)}); did=vc.DRUG_IDS["amikacin"]
   root=fx.build_cohort(Path(tmp)/"prefix","prefix",rules)
   fx.mutate_report(root,did,lambda r:r["splits"][0]["learners"][1]["metrics"].update(balanced_accuracy=0.123))
   with self.assertRaises(vc.ValidationError): vc.validate_manifest(root/"validation-manifest.json",rules)
   c=vc.validate_manifest(root/"validation-manifest.json",rules,collect_invalid=True)
   self.assertEqual(list(c.invalid),[did]); self.assertEqual(c.reports,{})
   path=root/"data/species"/vc.SPECIES/"results"/(did+".json"); path.write_text(path.read_text().replace('"seconds":10.0','"seconds":12.0',1))
   with self.assertRaises(vc.BindingError): vc.validate_manifest(root/"validation-manifest.json",rules,collect_invalid=True)
 def test_i1_a_bound_report_corrupted_into_invalid_json_is_a_binding_stop(self):
  # I1: the report's SHA-256 is checked before it is parsed. A bound report cut short is a
  # binding failure even when collecting, where it used to be recorded as an invalid report.
  with tempfile.TemporaryDirectory() as tmp:
   rules=fx.make_rules({"prefix":("amikacin",)}); did=vc.DRUG_IDS["amikacin"]
   root=fx.build_cohort(Path(tmp)/"prefix","prefix",rules)
   path=root/"data/species"/vc.SPECIES/"results"/(did+".json"); path.write_bytes(path.read_bytes()[:-40])
   with self.assertRaises(vc.BindingError) as stop: vc.validate_manifest(root/"validation-manifest.json",rules,collect_invalid=True)
   self.assertIn("report SHA-256 differs",str(stop.exception))
   # Bound again to the corrupted bytes, it is a contract failure (malformed JSON), not a binding one.
   manifest=json.loads((root/"validation-manifest.json").read_text()); manifest["reports"][did]["sha256"]=fx.sha(path); fx.write_json(root/"validation-manifest.json",manifest)
   c=vc.validate_manifest(root/"validation-manifest.json",rules,collect_invalid=True)
   self.assertIn("malformed JSON",c.invalid[did])
   NEGATIVE.append({"control":"I1 bound report corrupted into invalid JSON","result":"binding stop"})

class MatchedAndPinTests(unittest.TestCase):
 def test_the_frozen_matched2k_list_and_rules_account_for_all_22(self):
  self.assertEqual(vc.RULES["matched2k"].genomes_sha256,"212a21a5d2574e0d62d379fccddb0839c4f44f5a7cd36e43a558dc2c0926c3d1")
  for name in ("prefix","added","full"):
   rule=vc.RULES[name]; self.assertEqual(set(rule.drugs)|set(rule.unavailable),set(vc.ALL_DRUGS))
  self.assertEqual(vc.RULES["added"].unavailable,vc.ADDED_UNAVAILABLE)
 def test_a_matched2k_list_that_is_not_the_frozen_selection_is_rejected(self):
  with tempfile.TemporaryDirectory() as tmp:
   rules=fx.make_rules(); root=fx.build_cohort(Path(tmp)/"matched2k","matched2k",rules)
   vc.validate_manifest(root/"validation-manifest.json",rules)
   reject(self,"matched2k list differs from the frozen selection","frozen rule",lambda:vc.validate_manifest(root/"validation-manifest.json",vc.RULES))
 def test_a_checker_blob_that_moved_from_its_pin_stops_the_contract(self):
  with tempfile.TemporaryDirectory() as tmp:
   scripts=Path(tmp)/"scripts"
   for rel in list(vc.PINS)+["mtb-4k/report_contract.py","mtb-4k/fingerprint.py","mtb-4k/revision_n.py","mtb-4k/cleaned_rerun.py"]:
    (scripts/rel).parent.mkdir(parents=True,exist_ok=True); shutil.copyfile(vc.SCRIPTS/rel,scripts/rel)
   probe=[sys.executable,"-I","-B","-c",f"import sys; sys.argv=['probe','--profile',{pb.PROFILE!r}]; sys.path.insert(0,{str(scripts/'mtb-4k')!r}); import report_contract; print('loaded')"]
   self.assertEqual(subprocess.run(probe,capture_output=True,text=True).returncode,0)
   with (scripts/"check-results.py").open("a") as f: f.write("\n# edited\n")
   p=subprocess.run(probe,capture_output=True,text=True)
   self.assertNotEqual(p.returncode,0); self.assertIn("differs from the reviewed pin",p.stderr)
   NEGATIVE.append({"control":"checker blob moved from its pin","result":"rejected at import"})

if __name__=="__main__":
 result=unittest.main(exit=False,verbosity=2).result
 print(json.dumps({"negative_controls":len(NEGATIVE)}))
 sys.exit(0 if result.wasSuccessful() else 1)
