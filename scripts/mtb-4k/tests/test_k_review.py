#!/usr/bin/env python3
"""Independent synthetic/numerical review; no biological execution.

The K review's prior harness (Vast grm-mtb-4k-validators-review-kimi-20261028-r3,
independent/test_review.py, SHA-256 77f559fe...), adapted at call signatures only to run against
the canonical scripts/mtb-4k: the modules it imports, the temporary directories, the mutant
copy (a pinned scripts tree instead of K's upstream/), the wrapper's --receipt, a skip when the
Vast-only evidence is absent, and one assertion changed by decision D1 (marked).
IRLS0 (marked): {0 steps, 0 cap hits, converged} moved from the impossible statuses to an
accepted one. 690815a writes irls_steps 0 when the warm start already meets the KKT tolerance
(src/learn/lasso.rs:70), the pinned check-results.py:114 accepts any non-negative step count,
and the K review treated that status as impossible because it was at 855bdba.

usage: FOLDGEN=<foldgen binary> python3 -I scripts/mtb-4k/tests/test_k_review.py --receipt OUT.json
"""
import argparse,copy,csv,decimal,importlib.util,json,math,os,shutil,sys,tempfile,time,unittest
from pathlib import Path
from types import SimpleNamespace
LANE=Path(__file__).resolve().parents[1]; REVIEW=None
sys.path.insert(0,str(LANE)); sys.path.insert(0,str(LANE/"tests"))
import report_contract as vc, compare_cohorts as cc, run_revision8_checks as wrapper, fixtures as fx
EVIDENCE=Path("/root/grm-research/grm-research-mtb-scaleup-rev7cand-r1/results/audit.json")
def readable(path):
 try: return path.is_file() and os.access(path,os.R_OK)
 except OSError: return False
OBS={"cp_cases":[],"beta_cases":[],"mutations":[],"invalid_solver_status":[]}
D=decimal.Decimal; decimal.getcontext().prec=80

def cdf(n,k,p):
 if k<0:return D(0)
 if k>=n:return D(1)
 if k>n//2:return 1-survival(n,k+1,p)
 p=D(p); q=1-p; term=q**n; total=term
 for i in range(k):term=term*(n-i)*p/((i+1)*q); total+=term
 return total

def survival(n,k,p):
 if k<=0:return D(1)
 if k>n:return D(0)
 if k<=n//2:return 1-cdf(n,k-1,p)
 p=D(p); q=1-p; term=p**n; total=term
 for i in range(n,k,-1):term=term*i*q/((n-i+1)*p); total+=term
 return total

def survival_quantile(target,n,k):
 lo,hi=D(0),D(1); target=D(target)
 for _ in range(110):
  mid=(lo+hi)/2
  if survival(n,k,mid)<target:lo=mid
  else:hi=mid
 return (lo+hi)/2

def reference_cp(e,n):
 lower=D(0) if e==0 else survival_quantile("0.025",n,e); upper=D(1)
 if e<n:
  lo,hi=D(0),D(1)
  for _ in range(110):
   mid=(lo+hi)/2
   if cdf(n,e,mid)>D("0.025"):lo=mid
   else:hi=mid
  upper=(lo+hi)/2
 return lower,upper

def ratio(a,b):return a/b if b else 0.0

def metrics(truth,scores):
 pairs=[(t,s) for t,s in zip(truth,scores) if s is not None]
 tp=sum(t and s>=.5 for t,s in pairs); fp=sum(not t and s>=.5 for t,s in pairs)
 tn=sum(not t and s<.5 for t,s in pairs); fn=sum(t and s<.5 for t,s in pairs)
 sens,spec,prec=ratio(tp,tp+fn),ratio(tn,tn+fp),ratio(tp,tp+fp)
 den=math.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)); pos=[s for t,s in pairs if t]; neg=[s for t,s in pairs if not t]
 auroc=sum((p>q)+.5*(p==q) for p in pos for q in neg)/(len(pos)*len(neg)) if pos and neg else None
 auprc=None
 if pos:
  auprc=before=0.0
  for cut in sorted({s for _,s in pairs},reverse=True):
   chosen=[t for t,s in pairs if s>=cut]; recall=sum(chosen)/len(pos)
   auprc+=(recall-before)*sum(chosen)/len(chosen); before=recall
 return {"n":len(pairs),"tp":tp,"fp":fp,"tn":tn,"fn":fn,"accuracy":ratio(tp+tn,len(pairs)),
  "balanced_accuracy":(sens+spec)/2,"sensitivity":sens,"specificity":spec,"precision":prec,
  "f1":0.0 if prec+sens==0 else 2*prec*sens/(prec+sens),"mcc":0.0 if den==0 else (tp*tn-fp*fn)/den,
  "vme":ratio(fn,tp+fn),"me":ratio(fp,fp+tn),"auroc":auroc,"auprc":auprc}

def roc(truth,scores):
 pairs=[(t,s) for t,s in zip(truth,scores) if s is not None]; np=sum(t for t,_ in pairs); nn=len(pairs)-np
 points=[[0.0,0.0]]
 for cut in sorted({s for _,s in pairs},reverse=True):
  chosen=[t for t,s in pairs if s>=cut]; points.append([(len(chosen)-sum(chosen))/nn,sum(chosen)/np])
 return points

def fit(status=None,empty=False):
 if empty:return {"path_points":0,"path_unconverged":0,"selected_point":None,"selected_convergence":None,"refit_convergence":None}
 status=status or {"irls_steps":2,"sweep_cap_hits":0,"converged":True}
 return {"path_points":1,"path_unconverged":0,"selected_point":0,"selected_convergence":copy.deepcopy(status),"refit_convergence":copy.deepcopy(status)}

class NumericalTests(unittest.TestCase):
 def test_cp_decimal_binomial_tails(self):
  cases=[]
  for n in (1,2,3,10,25,100,1000,4000):
   for e in sorted({0,1,n//4,n//2,max(0,n-1),n}):cases.append((e,n))
  worst=D(0)
  for e,n in cases:
   got=vc.clopper_pearson(e,n); lower,upper=reference_cp(e,n)
   errors=[abs(D(got["lower"])-lower),abs(D(got["upper"])-upper)]; worst=max(worst,*errors)
   OBS["cp_cases"].append({"errors":e,"denominator":n,"lower":got["lower"],"upper":got["upper"],
    "reference_lower":str(lower),"reference_upper":str(upper),"absolute_errors":[str(x) for x in errors]})
   self.assertLessEqual(float(errors[0]),2e-12,(e,n,"lower")); self.assertLessEqual(float(errors[1]),2e-12,(e,n,"upper"))
  OBS["cp_max_absolute_error"]=str(worst); OBS["cp_case_count"]=len(cases)

 def test_beta_cdf_quantile_tails(self):
  failures=[]; shapes=((1,1),(1,5),(5,1),(2,3),(3,10),(10,3),(500,500),(1,4000),(4000,1))
  for a,b in shapes:
   n=a+b-1
   for x in (1e-12,1e-6,.01,.1,.5,.9,.999999):
    got=vc.regularized_beta(x,a,b); want=survival(n,a,D(x)); error=abs(D(got)-want)
    OBS["beta_cases"].append({"kind":"cdf","a":a,"b":b,"x":x,"error":str(error)})
    if float(error)>2e-12:failures.append(("cdf",a,b,x,float(error)))
  for a,b in ((1,1),(1,5),(5,1),(2,3),(3,10),(1,4000),(4000,1)):
   n=a+b-1
   for p in (1e-12,.001,.025,.5,.975,.999,1-1e-12):
    got=vc.beta_ppf(p,a,b); want=survival_quantile(p,n,a); error=abs(D(got)-want)
    OBS["beta_cases"].append({"kind":"ppf","a":a,"b":b,"p":p,"error":str(error)})
    if float(error)>2e-12:failures.append(("ppf",a,b,p,float(error)))
  self.assertEqual(failures,[])
  self.assertEqual(vc.regularized_beta(0.0,2,3),0.0); self.assertEqual(vc.regularized_beta(1.0,2,3),1.0)

 def test_cp_edges_and_denominator_guards(self):
  for n in (1,2,10,25,4000):
   self.assertAlmostEqual(vc.clopper_pearson(0,n)["upper"],1-.025**(1/n),places=14)
   self.assertAlmostEqual(vc.clopper_pearson(n,n)["lower"],.025**(1/n),places=14)
  undefined=vc.clopper_pearson(0,0)
  self.assertEqual((undefined["lower"],undefined["upper"]),(None,None)); self.assertEqual(undefined["status"],"undefined_no_called_denominator")
  for e,n in ((-1,10),(11,10),(True,10),(1.5,10),(1,-1),(1,10.0)):
   with self.subTest(errors=e,denominator=n):
    with self.assertRaises(vc.ValidationError):vc.clopper_pearson(e,n)

 def test_beta_bisection_mutation_detected(self):
  with tempfile.TemporaryDirectory(dir=REVIEW) as tmp:
   needle="for _ in range(160):"
   self.assertEqual((vc.LANE/"report_contract.py").read_text().count(needle),1)
   mutant=fx.contract_copy(Path(tmp)/"mutant",needle,"for _ in range(4):")
   spec=importlib.util.spec_from_file_location("review_cp_mutant",mutant); module=importlib.util.module_from_spec(spec)
   sys.modules[spec.name]=module
   try:spec.loader.exec_module(module); got=module.clopper_pearson(3,10)
   finally:sys.modules.pop(spec.name,None)
  lower,upper=reference_cp(3,10); discrepancy=max(abs(D(got["lower"])-lower),abs(D(got["upper"])-upper))
  self.assertGreater(float(discrepancy),1e-3)
  OBS["mutations"].append({"mutation":"beta_ppf 160 -> 4 iterations","detected":True,"discrepancy":str(discrepancy)})

class MetricPairingTests(unittest.TestCase):
 def test_known_metrics_no_calls_fold_denominators(self):
  truth=[True,True,True,True,False,False,False,False]; scores=[.9,.4,None,.5,.6,.2,None,.5]; folds=[0,0,1,0,0,0,1,0]
  samples=[vc.Sample(f"g{i}",truth[i],i//2) for i in range(8)]; expected=metrics(truth,scores)
  self.assertEqual((expected["tp"],expected["fp"],expected["tn"],expected["fn"]),(2,2,1,1))
  row={"learner":"lasso","metrics":expected,"fold_balanced_accuracy":[.5,None],"roc":roc(truth,scores),
   "scores":scores,"conformal":[],"fold_complexity":[2.0,0.0],"seconds":0.0,"no_call":{"resistant":1,"susceptible":1},
   "lasso_fits":[fit(),fit(empty=True)]}
  got=vc.validate_learner(row,"random",folds,samples,2,False,30,"independent-known-answer").receipt
  self.assertEqual(got["metrics"]["n"],6); self.assertEqual(got["metrics"]["sensitivity"],2/3)
  self.assertEqual(got["metrics"]["specificity"],1/3); self.assertEqual(got["metrics"]["balanced_accuracy"],.5)
  self.assertEqual(got["metrics"]["vme"],1/3); self.assertEqual(got["metrics"]["me"],2/3)
  self.assertEqual(got["no_call"],{"resistant":1,"susceptible":1})
  self.assertEqual((got["clopper_pearson_95"]["vme"]["errors"],got["clopper_pearson_95"]["vme"]["denominator"]),(1,3))
  self.assertEqual((got["clopper_pearson_95"]["me"]["errors"],got["clopper_pearson_95"]["me"]["denominator"]),(2,3))
  self.assertEqual((got["folds"][0]["held_out"],got["folds"][0]["called"]),(6,6))
  self.assertEqual((got["folds"][1]["held_out"],got["folds"][1]["called"]),(2,0))
  self.assertEqual(got["folds"][1]["no_call"],{"resistant":1,"susceptible":1}); self.assertIsNone(got["folds"][1]["metrics"])

 @staticmethod
 def pairing_report(order,truth,folds,scores,scheme="random"):
  samples={g:vc.Sample(g,truth[g],0) for g in order}; learner=vc.LearnerData("cart",scheme,scores,folds,{})
  split=vc.SplitData(scheme,True,None,folds,{"cart":learner},{})
  return SimpleNamespace(samples=samples,sample_order=order,splits={scheme:split})

 def test_genome_fold_pairing_not_row_order(self):
  order=["g1","g2","g3","g4"]; truth={"g1":True,"g2":False,"g3":True,"g4":False}
  pf={"g1":0,"g2":1,"g3":2,"g4":3}; ff={"g1":0,"g2":2,"g3":2,"g4":3}
  ps={"g1":.4,"g2":.2,"g3":None,"g4":.8}; fs={"g1":.6,"g2":.8,"g3":.7,"g4":None}
  prefix=self.pairing_report(order,truth,pf,ps); full=self.pairing_report(list(reversed(order)),truth,ff,fs)
  row,issues=cc.pair_learner(prefix,full,"random","cart",set(order),{g:True for g in order})
  self.assertEqual(row["pairing"]["fold_mismatch"],["g2"]); self.assertEqual(row["pairing"]["eligible_matched_pairs"],3)
  self.assertEqual((row["matched_fold"]["fixed"],row["matched_fold"]["broken"]),(1,0))
  self.assertEqual(row["matched_fold"]["prefix_no_call_full_call"],1); self.assertEqual(row["matched_fold"]["prefix_call_full_no_call"],1)
  self.assertEqual(row["matched_fold"]["score_delta_full_minus_prefix"]["n"],1)
  self.assertAlmostEqual(row["matched_fold"]["score_delta_full_minus_prefix"]["mean"],.2)
  # CHANGED ASSERTION (decision D1): K's comparator returned "random folds differ" as a blocking
  # issue; the canonical one invalidates this row's matched random control instead.
  self.assertEqual(row["genome_aligned_all_folds"]["broken"],1); self.assertFalse(row["matched_fold"]["control_valid"]); self.assertEqual(issues,[])

 def test_changed_lineage_group_excluded(self):
  order=["g1","g2","g3","g4"]; truth={"g1":True,"g2":False,"g3":True,"g4":False}; folds={"g1":0,"g2":1,"g3":2,"g4":3}
  ps={"g1":.4,"g2":.2,"g3":None,"g4":.8}; fs={"g1":.6,"g2":.8,"g3":.7,"g4":None}
  prefix=self.pairing_report(order,truth,folds,ps,"lineage"); full=self.pairing_report(list(reversed(order)),truth,folds,fs,"lineage")
  stable={"g1":True,"g2":True,"g3":False,"g4":True}; row,issues=cc.pair_learner(prefix,full,"lineage","cart",set(order),stable)
  self.assertEqual(row["pairing"]["fold_mismatch"],[]); self.assertEqual(row["pairing"]["changed_lineage_group"],["g3"])
  self.assertEqual(row["matched_fold"]["pairs"],3); self.assertEqual((row["matched_fold"]["fixed"],row["matched_fold"]["broken"]),(1,1)); self.assertEqual(issues,[])

 def test_five_learner_null_median(self):
  rows=[]
  for i,name in enumerate(vc.LEARNERS):
   m=None if name=="lasso" else {"sensitivity":.5+i/10,"specificity":.8,"balanced_accuracy":.5+i/10,"vme":.2,"me":.1}
   rows.append({"learner":name,"metrics":m,"no_call_rate":i/10})
  med=vc.learner_medians(rows)["balanced_accuracy"]
  self.assertIsNone(med["median"]); self.assertEqual(med["null_learners"],["lasso"]); self.assertEqual(med["prespecified_learners"],5)
  self.assertAlmostEqual(med["available_only_median"],.65)

class FrozenEvidenceTests(unittest.TestCase):
 @unittest.skipUnless(readable(EVIDENCE),"frozen MTB evidence exists only on Vast; run this file there")
 def test_added_gate_counts_derive_from_frozen_evidence(self):
  with open("/root/grm-research/grm-research-mtb-scaleup-rev7cand-r1/results/audit.json") as f:rows=json.load(f)["class_counts"]
  with open("/root/grm-research/grm-research-mtb-4k-cohort-r2/results/class-gates.tsv") as f:gates={r["dataset_id"]:r for r in csv.DictReader(f,delimiter="\t")}
  self.assertEqual(len(rows),22); derived={}
  for row in rows:
   drug=row["drug"]; gate=gates[vc.DRUG_IDS[drug]]
   self.assertEqual((row["hyp4k_R"],row["hyp4k_S"]),(int(gate["resistant"]),int(gate["susceptible"])))
   self.assertTrue(row["cand2k_gate"] and row["hyp4k_gate"]); derived[drug]=(row["hyp4k_R"]-row["cand2k_R"],row["hyp4k_S"]-row["cand2k_S"])
  unavailable={d:rs for d,rs in derived.items() if min(rs)<25}; eligible=tuple(d for d in vc.ALL_DRUGS if d not in unavailable)
  self.assertEqual(unavailable,vc.ADDED_UNAVAILABLE); self.assertEqual(eligible,vc.ADDED_DRUGS)
  OBS["derived_added_class_counts"]=derived; OBS["derived_added_unavailable"]=unavailable

class SolverStatusTests(unittest.TestCase):
 def test_valid_intermediate_cap_hit_accepted(self):
  self.assertTrue(vc.validate_lasso_fit(fit({"irls_steps":2,"sweep_cap_hits":1,"converged":True}),"valid-intermediate-cap",30))

 def test_impossible_iteration_cap_counts_fail_closed(self):
  # CHANGED ASSERTION (IRLS0): {0,0,true} left the impossible set; 690815a writes it on a valid warm start.
  accepted=[]; statuses=({"irls_steps":1,"sweep_cap_hits":1,"converged":True},{"irls_steps":1,"sweep_cap_hits":2,"converged":True})
  for status in statuses:
   try:vc.validate_lasso_fit(fit(status),"impossible-status",30)
   except vc.ValidationError:continue
   accepted.append(status)
  OBS["invalid_solver_status"].append({"direct_accepted":accepted}); self.assertEqual(accepted,[])
  self.assertTrue(vc.validate_lasso_fit(fit({"irls_steps":0,"sweep_cap_hits":0,"converged":True}),"zero-step-warm-start",30))

 def test_full_wrapper_rejects_impossible_status_excess_path(self):
  with tempfile.TemporaryDirectory(dir=REVIEW) as tmp:
   rules=fx.make_rules(); root=fx.build_cohort(Path(tmp)/"prefix","prefix",rules); did=vc.DRUG_IDS["amikacin"]
   path=root/"data/species"/vc.SPECIES/"results"/(did+".json"); report=json.loads(path.read_text())
   report["models"]["lasso_fit"]["path_points"]=report["eval_config"]["lasso_lambdas"]+1
   report["models"]["lasso_fit"]["selected_convergence"]={"irls_steps":1,"sweep_cap_hits":1,"converged":True}
   fx.write_json(path,report); fx.refresh_root(root,"prefix")
   args=[str(root/"audit.json"),"--root",str(root),"--run-manifest",str(root/"validation-manifest.json"),
    "--expected-source-commit",fx.SOURCE,"--binary",str(root/"synthetic-binary"),"--backup",str(root/"backup"),
    "--log",str(root/"pipeline.log"),"--launch","2000-01-01T00:00:00Z","--foldgen",str(root/"synthetic-foldgen"),
    "--expected-foldgen-sha256",fx.sha(root/"synthetic-foldgen"),  # path-adapted (I4)
    "--checker-out",str(root/"checker.json"),"--summary-out",str(root/"summary.json"),"--receipt",str(root/"receipt.json")]
   code=wrapper.main(args,rules=rules,foldgen_factory=fx.FakeFoldGen); valid=None
   if (root/"audit.json").exists():valid=json.loads((root/"audit.json").read_text())["summary"]["revision8_wrapper"]["valid"]
   OBS["invalid_solver_status"].append({"wrapper_exit":code,"wrapper_valid":valid,"path_points":31,
    "configured_lasso_lambdas":30,"selected_convergence":{"irls_steps":1,"sweep_cap_hits":1,"converged":True}})
   self.assertNotEqual(code,0,"full checker/auditor/wrapper accepted impossible revision-8 solver evidence")

def flatten(suite):
 for item in suite:
  if isinstance(item,unittest.TestSuite):yield from flatten(item)
  else:yield item

def main():
 parser=argparse.ArgumentParser(); parser.add_argument("--receipt",required=True,type=Path); args=parser.parse_args(); started=time.time()
 suite=unittest.defaultTestLoader.loadTestsFromModule(sys.modules[__name__]); ids=[t.id() for t in flatten(suite)]
 result=unittest.TextTestRunner(verbosity=2).run(suite); failed={t.id() for t,_ in result.failures+result.errors}
 receipt={"schema":"grm-mtb-4k-independent-review-tests-v1","BIOLOGICAL_EXECUTION":"NO","RUN_4K":"NO","synthetic_tests_only":True,
  "tests_run":result.testsRun,"failures":len(result.failures),"errors":len(result.errors),"skipped":len(result.skipped),"passed":result.wasSuccessful(),
  "tests":[{"id":i,"status":"failed" if i in failed else "passed"} for i in ids],
  "failure_details":[{"test":t.id(),"traceback":tb} for t,tb in result.failures+result.errors],
  "wall_seconds":time.time()-started,"observed":OBS}
 args.receipt.write_text(json.dumps(receipt,indent=1,sort_keys=True,allow_nan=False)+"\n"); return 0 if result.wasSuccessful() else 1
if __name__=="__main__":sys.exit(main())
