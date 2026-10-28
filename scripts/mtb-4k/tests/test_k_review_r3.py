#!/usr/bin/env python3
"""Fresh R3 synthetic and numerical review; no biological execution.

The K review's R3 harness (Vast grm-mtb-4k-validators-review-kimi-20261028-r3,
independent/test_r3.py, SHA-256 bbc50baa...), adapted at call signatures only to run against the
canonical scripts/mtb-4k: the modules it imports, the temporary directories, the wrapper's
--receipt and the real foldgen (these fixtures' folds come from it; K's runs used an index-fold
double), and the mutant copies (a pinned scripts tree instead of K's upstream/). No assertion
changed.

usage: FOLDGEN=<foldgen binary> python3 -I scripts/mtb-4k/tests/test_k_review_r3.py --receipt OUT.json
"""
import argparse,contextlib,copy,decimal,functools,importlib.util,io,json,math,shutil,sys,tempfile,time,unittest
from pathlib import Path
LANE=Path(__file__).resolve().parents[1]; REVIEW=None
sys.path.insert(0,str(LANE)); sys.path.insert(0,str(LANE/"tests"))
import report_contract as vc, run_revision8_checks as wrapper, fixtures as fx
D=decimal.Decimal; decimal.getcontext().prec=100
TOL=D("2e-12"); EXTREME_P=1-1e-12
OBS={"beta_cases":[],"cp95_cases":[],"former_extreme_beta":[],"closed_form_cases":[],
     "direct_negative_controls":[],"wrapper_negative_controls":[],"wrapper_positive_controls":[],"mutations":[]}

@functools.lru_cache(maxsize=None)
def beta_norm(a,b):
    return D(math.factorial(a+b-1))/(D(math.factorial(a-1))*D(math.factorial(b-1)))

def decimal_betacf(a,b,x):
    tiny=D("1e-900"); eps=D("1e-80"); qab=a+b; qap=a+1; qam=a-1
    c=D(1); d=1-qab*x/qap
    if abs(d)<tiny: d=tiny
    d=1/d; h=d
    for m in range(1,2001):
        m2=2*m
        aa=D(m)*(b-m)*x/((qam+m2)*(a+m2)); d=1+aa*d
        if abs(d)<tiny: d=tiny
        c=1+aa/c
        if abs(c)<tiny: c=tiny
        d=1/d; h*=d*c
        aa=-(a+m)*(qab+m)*x/((a+m2)*(qap+m2)); d=1+aa*d
        if abs(d)<tiny: d=tiny
        c=1+aa/c
        if abs(c)<tiny: c=tiny
        d=1/d; delta=d*c; h*=delta
        if abs(delta-1)<=eps: return h
    raise ArithmeticError(("decimal beta continued fraction",a,b,x))

def oracle_beta(x,a,b):
    """100-digit continued fraction with exact integer-beta normalization."""
    x=D(x); a=int(a); b=int(b)
    if x==0 or x==1: return x
    front=x**a*(1-x)**b*beta_norm(a,b)
    da=D(a); db=D(b)
    if x<(da+1)/(da+db+2): return front*decimal_betacf(da,db,x)/da
    return 1-front*decimal_betacf(db,da,1-x)/db

@functools.lru_cache(maxsize=None)
def oracle_ppf(p,a,b):
    """Direct 100-digit bisection; deliberately no candidate-style complement."""
    p=D(p); lo=D(0); hi=D(1)
    for _ in range(180):
        mid=(lo+hi)/2
        if mid==lo or mid==hi: break
        if oracle_beta(mid,a,b)<p: lo=mid
        else: hi=mid
    return (lo+hi)/2

def reference_cp(e,n,alpha=0.05):
    q=D(alpha)/2
    lower=D(0) if e==0 else oracle_ppf(q,e,n-e+1)
    upper=D(1) if e==n else oracle_ppf(1-q,e+1,n-e)
    return lower,upper

def error(got,want): return abs(D(got)-want)

def wrapper_args(root):
    return [str(root/"audit.json"),"--root",str(root),"--run-manifest",str(root/"validation-manifest.json"),
            "--expected-source-commit",fx.SOURCE,"--binary",str(root/"synthetic-binary"),"--backup",str(root/"backup"),
            "--log",str(root/"pipeline.log"),"--launch","2000-01-01T00:00:00Z","--foldgen",fx.foldgen_path(),
            "--checker-out",str(root/"checker.json"),"--summary-out",str(root/"summary.json"),"--receipt",str(root/"receipt.json")]

class NumericalTests(unittest.TestCase):
    def test_beta_112_prior_probes(self):
        failures=[]
        shapes=((1,1),(1,5),(5,1),(2,3),(3,10),(10,3),(500,500),(1,4000),(4000,1))
        for a,b in shapes:
            for x in (1e-12,1e-6,.01,.1,.5,.9,.999999):
                got=vc.regularized_beta(x,a,b); want=oracle_beta(D(x),a,b); err=error(got,want)
                OBS["beta_cases"].append({"kind":"cdf","a":a,"b":b,"x":x,"got":got,"reference":str(want),"absolute_error":str(err)})
                if err>TOL: failures.append(("cdf",a,b,x,str(err)))
        former={(1,5),(2,3),(3,10),(1,4000)}
        for a,b in ((1,1),(1,5),(5,1),(2,3),(3,10),(1,4000),(4000,1)):
            for p in (1e-12,.001,.025,.5,.975,.999,EXTREME_P):
                got=vc.beta_ppf(p,a,b); want=oracle_ppf(D(p),a,b); err=error(got,want)
                row={"kind":"ppf","a":a,"b":b,"p":p,"got":got,"reference":str(want),"absolute_error":str(err)}
                OBS["beta_cases"].append(row)
                if p==EXTREME_P and (a,b) in former: OBS["former_extreme_beta"].append(row)
                if err>TOL: failures.append(("ppf",a,b,p,str(err)))
        self.assertEqual(vc.regularized_beta(0.0,2,3),0.0)
        self.assertEqual(vc.regularized_beta(1.0,2,3),1.0)
        OBS["beta_probe_count"]=len(OBS["beta_cases"])
        OBS["beta_max_absolute_error"]=str(max(D(r["absolute_error"]) for r in OBS["beta_cases"]))
        self.assertEqual(OBS["beta_probe_count"],112)
        self.assertEqual(len(OBS["former_extreme_beta"]),4)
        self.assertEqual(failures,[])

    def test_cp95_39_cases_and_closed_forms(self):
        failures=[]; cases=[]
        for n in (1,2,3,10,25,100,1000,4000):
            for e in sorted({0,1,n//4,n//2,max(0,n-1),n}): cases.append((e,n))
        worst=D(0)
        for e,n in cases:
            got=vc.clopper_pearson(e,n); lower,upper=reference_cp(e,n)
            errors=[error(got["lower"],lower),error(got["upper"],upper)]; worst=max(worst,*errors)
            OBS["cp95_cases"].append({"errors":e,"denominator":n,"lower":got["lower"],"upper":got["upper"],
                "reference_lower":str(lower),"reference_upper":str(upper),"absolute_errors":[str(v) for v in errors]})
            if errors[0]>TOL: failures.append((e,n,"lower",str(errors[0])))
            if errors[1]>TOL: failures.append((e,n,"upper",str(errors[1])))
        OBS["cp95_case_count"]=len(cases); OBS["cp95_endpoint_count"]=2*len(cases)
        OBS["cp95_max_absolute_error"]=str(worst)
        self.assertEqual(len(cases),39); self.assertEqual(failures,[])
        q=D("0.025")
        for n in (1,2,10,25,4000):
            for label,got,want in (
                ("zero_error_upper",vc.clopper_pearson(0,n)["upper"],1-q**(D(1)/D(n))),
                ("all_error_lower",vc.clopper_pearson(n,n)["lower"],q**(D(1)/D(n)))):
                err=error(got,want); OBS["closed_form_cases"].append({"case":label,"n":n,"absolute_error":str(err)})
                self.assertLessEqual(err,TOL,(label,n))
        alpha=2e-12; qext=D(alpha)/2
        got=vc.beta_ppf(EXTREME_P,1,5); want=1-(1-D(EXTREME_P))**(D(1)/D(5)); err=error(got,want)
        OBS["closed_form_cases"].append({"case":"beta_ppf_1_minus_1e-12_1_5","got":got,"reference":str(want),"absolute_error":str(err)})
        self.assertLessEqual(err,TOL)
        got=vc.clopper_pearson(0,5,alpha=alpha)["upper"]; want=1-qext**(D(1)/D(5)); err=error(got,want)
        OBS["closed_form_cases"].append({"case":"cp_zero_of_five_alpha_2e-12","got":got,"reference":str(want),"absolute_error":str(err)})
        self.assertLessEqual(err,TOL)
        OBS["closed_form_case_count"]=len(OBS["closed_form_cases"])

class SolverTests(unittest.TestCase):
    def test_direct_status_and_path_matrix(self):
        for points in (1,7,29,30):
            fit=fx.fit(); fit["path_points"]=points; fit["selected_point"]=0 if points==1 else 2
            for key in ("selected_convergence","refit_convergence"):
                fit[key]={"irls_steps":2,"sweep_cap_hits":1,"converged":True}
            self.assertTrue(vc.validate_lasso_fit(fit,"direct-valid-short",30))
        self.assertFalse(vc.validate_lasso_fit(fx.fit(empty=True),"direct-valid-empty",30))
        statuses=(("zero_step",{"irls_steps":0,"sweep_cap_hits":0,"converged":True}),
                  ("excess_cap",{"irls_steps":1,"sweep_cap_hits":2,"converged":True}),
                  ("all_capped",{"irls_steps":1,"sweep_cap_hits":1,"converged":True}))
        for label,status in statuses:
            for key in ("selected_convergence","refit_convergence"):
                fit=fx.fit(); fit[key]=copy.deepcopy(status)
                with self.subTest(label=label,key=key):
                    with self.assertRaises(vc.ValidationError): vc.validate_lasso_fit(fit,"direct-invalid",30)
                OBS["direct_negative_controls"].append({"kind":label,"status":key,"rejected":True})
        fit=fx.fit(); fit["path_points"]=31
        with self.assertRaises(vc.ValidationError): vc.validate_lasso_fit(fit,"direct-excess-path",30)
        OBS["direct_negative_controls"].append({"kind":"excess_path","rejected":True})
        self.assertEqual(len(OBS["direct_negative_controls"]),7)

class WrapperTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory(dir=REVIEW); self.rules=fx.make_rules(); self.did=vc.DRUG_IDS["amikacin"]
    def tearDown(self): self.tmp.cleanup()

    def run_case(self,label,cohort="prefix",mutate=None):
        with tempfile.TemporaryDirectory(dir=self.tmp.name) as case:
            root=fx.build_cohort(Path(case)/cohort,cohort,self.rules)
            if mutate is not None: fx.mutate_report(root,self.did,mutate)
            stdout=io.StringIO(); stderr=io.StringIO()
            with contextlib.redirect_stdout(stdout),contextlib.redirect_stderr(stderr):
                code=wrapper.main(wrapper_args(root),rules=self.rules)
            result={"label":label,"cohort":cohort,"exit":code,"checker_exists":(root/"checker.json").exists(),
                    "audit_exists":(root/"audit.json").exists(),"summary_exists":(root/"summary.json").exists(),
                    "stdout":stdout.getvalue(),"stderr":stderr.getvalue()}
            if result["audit_exists"]:
                audit=json.loads((root/"audit.json").read_text()); summary=audit["summary"]
                result["valid"]=summary["revision8_wrapper"]["valid"]; result["by_revision"]=summary["by_revision"]
                result["class_minimum_skipped"]=len(summary["class_minimum_skipped"])
            if result["checker_exists"]:
                checker=json.loads((root/"checker.json").read_text())
                result["checker_reports"]=len(checker["reports"]); result["checker_skips"]=len(checker["class_minimum_skipped"])
            return result

    def negative(self,case):
        self.assertEqual(case["exit"],2,case)
        self.assertFalse(case["checker_exists"],case)
        self.assertFalse(case["audit_exists"],case)
        self.assertFalse(case["summary_exists"],case)
        self.assertIn("VALIDATION FAILED:",case["stderr"])
        OBS["wrapper_negative_controls"].append(case)

    def positive(self,case):
        self.assertEqual(case["exit"],0,case)
        self.assertTrue(case["checker_exists"] and case["audit_exists"] and case["summary_exists"],case)
        self.assertTrue(case["valid"],case)
        OBS["wrapper_positive_controls"].append(case)

    def test_final_and_per_fold_fail_before_checker(self):
        statuses=(("zero_step",{"irls_steps":0,"sweep_cap_hits":0,"converged":True}),
                  ("excess_cap",{"irls_steps":1,"sweep_cap_hits":2,"converged":True}),
                  ("all_capped",{"irls_steps":1,"sweep_cap_hits":1,"converged":True}))
        cases=[]
        for location in ("final","per_fold"):
            cases.append((location,"selected_convergence",statuses[2][0],statuses[2][1],31,"prior_path31_selected_all_capped"))
            for key in ("selected_convergence","refit_convergence"):
                for label,status in statuses: cases.append((location,key,label,status,None,label))
            cases.append((location,None,"excess_path_only",None,31,"excess_path_only"))
        self.assertEqual(len(cases),16)
        for location,key,label,status,points,kind in cases:
            with self.subTest(location=location,key=key,kind=kind):
                def mutate(r,location=location,key=key,status=status,points=points):
                    if location=="final": fit=r["models"]["lasso_fit"]
                    else:
                        split=next(s for s in r["splits"] if s["scheme"]=="random")
                        row=next(row for row in split["learners"] if row["learner"]=="lasso")
                        fit=row["lasso_fits"][0]
                    if points is not None: fit["path_points"]=points
                    if key is not None: fit[key]=copy.deepcopy(status)
                case=self.run_case(f"{location}/{key}/{kind}",mutate=mutate)
                case.update({"location":location,"status_field":key,"kind":kind,"path_points":points,"status":status})
                self.negative(case)
        self.assertEqual(len(OBS["wrapper_negative_controls"]),16)

    def test_valid_intermediate_short_and_empty_pass(self):
        empty=[]
        def valid_short(r):
            fits=[r["models"]["lasso_fit"]]
            for split in r["splits"]:
                if "unavailable_reason" in split: continue
                row=next(row for row in split["learners"] if row["learner"]=="lasso")
                fits.extend(row["lasso_fits"])
            for fit in fits:
                if fit["path_points"]==0: empty.append(True); continue
                fit["path_points"]=7
                for key in ("selected_convergence","refit_convergence"):
                    fit[key]={"irls_steps":2,"sweep_cap_hits":1,"converged":True}
        case=self.run_case("valid_intermediate_cap_short_and_per_fold_empty",mutate=valid_short)
        self.positive(case); self.assertTrue(empty); self.assertEqual(case["by_revision"],{"8":22})
        def final_empty(r):
            r["models"]["lasso_fit"]=fx.fit(empty=True); r["models"]["lasso"]=None
        case=self.run_case("valid_final_empty_path",mutate=final_empty)
        self.positive(case); self.assertEqual(case["by_revision"],{"8":22})

    def test_added_14_and_eight_skips_wrapper(self):
        case=self.run_case("added_14_and_eight_skips",cohort="added")
        self.positive(case)
        self.assertEqual(case["by_revision"],{"8":14})
        self.assertEqual(case["class_minimum_skipped"],8)
        self.assertEqual(case["checker_reports"],14)
        self.assertEqual(case["checker_skips"],8)

class MutationTests(unittest.TestCase):
    @contextlib.contextmanager
    def mutated(self,name,old,new):
        with tempfile.TemporaryDirectory(dir=REVIEW) as tmp:
            mutant=fx.contract_copy(Path(tmp)/"mutant",old,new)
            spec=importlib.util.spec_from_file_location("r3_"+name,mutant)
            module=importlib.util.module_from_spec(spec); sys.modules[spec.name]=module
            try:
                spec.loader.exec_module(module); yield module
            finally: sys.modules.pop(spec.name,None)

    def test_beta_guard_source_mutations_detected(self):
        with self.mutated("without_beta_symmetry","if p>0.5: return 1.0-beta_ppf(1.0-p,b,a)","if False: return 1.0-beta_ppf(1.0-p,b,a)") as mutant:
            detected=[]
            for a,b in ((1,5),(2,3),(3,10),(1,4000)):
                got=mutant.beta_ppf(EXTREME_P,a,b); err=error(got,oracle_ppf(D(EXTREME_P),a,b))
                self.assertGreater(err,TOL,(a,b,err)); detected.append(str(err))
            OBS["mutations"].append({"mutation":"remove beta upper-tail symmetry","detected":True,"former_extreme_errors":detected})
        with self.mutated("without_survival_complement","return 1.0-beta_ppf(q,b,a)","return beta_ppf(q,b,a)") as mutant:
            got=mutant.clopper_pearson(0,5,alpha=2e-12)["upper"]
            want=1-(D(2e-12)/2)**(D(1)/D(5)); err=error(got,want)
            self.assertGreater(err,TOL)
            OBS["mutations"].append({"mutation":"remove survival outer complement","detected":True,"absolute_error":str(err)})

    def test_solver_guard_source_mutations_detected(self):
        with self.mutated("without_positive_steps",'steps=integer(status["irls_steps"],where+".irls_steps",1)','steps=integer(status["irls_steps"],where+".irls_steps")') as mutant:
            status={"irls_steps":0,"sweep_cap_hits":0,"converged":False}
            with self.assertRaises(vc.ValidationError): vc.validate_convergence(status,"baseline-positive-steps")
            mutant.validate_convergence(status,"mutant-positive-steps")
            OBS["mutations"].append({"mutation":"remove positive IRLS requirement","detected":True})
        with self.mutated("without_cap_bound","if hits>steps: fail","if False: fail") as mutant:
            fit=fx.fit(); fit["selected_convergence"]={"irls_steps":1,"sweep_cap_hits":2,"converged":True}
            self.assertTrue(mutant.validate_lasso_fit(fit,"mutant-cap-bound",30))
            OBS["mutations"].append({"mutation":"remove sweep_cap_hits upper bound","detected":True})
        with self.mutated("without_final_cap_guard","if converged and hits==steps: fail","if False: fail") as mutant:
            fit=fx.fit(); fit["selected_convergence"]={"irls_steps":1,"sweep_cap_hits":1,"converged":True}
            self.assertTrue(mutant.validate_lasso_fit(fit,"mutant-final-cap",30))
            OBS["mutations"].append({"mutation":"remove final-step cap convergence guard","detected":True})
        with self.mutated("without_path_bound",'if fit["path_points"]>lambda_count: fail','if False: fail') as mutant:
            fit=fx.fit(); fit["path_points"]=31
            self.assertTrue(mutant.validate_lasso_fit(fit,"mutant-path-bound",30))
            OBS["mutations"].append({"mutation":"remove configured lambda path bound","detected":True})

def flatten(suite):
    for item in suite:
        if isinstance(item,unittest.TestSuite): yield from flatten(item)
        else: yield item

def main():
    parser=argparse.ArgumentParser(); parser.add_argument("--receipt",required=True,type=Path); args=parser.parse_args()
    started=time.time(); suite=unittest.defaultTestLoader.loadTestsFromModule(sys.modules[__name__]); ids=[t.id() for t in flatten(suite)]
    result=unittest.TextTestRunner(verbosity=2).run(suite); failed={t.id() for t,_ in result.failures+result.errors}
    receipt={"schema":"grm-mtb-4k-independent-review-r3-tests-v1","BIOLOGICAL_EXECUTION":"NO","RUN_4K":"NO","synthetic_tests_only":True,
             "tests_run":result.testsRun,"failures":len(result.failures),"errors":len(result.errors),"skipped":len(result.skipped),"passed":result.wasSuccessful(),
             "direct_negative_control_count":len(OBS["direct_negative_controls"]),"wrapper_negative_control_count":len(OBS["wrapper_negative_controls"]),
             "wrapper_positive_control_count":len(OBS["wrapper_positive_controls"]),"source_mutation_count":len(OBS["mutations"]),
             "tests":[{"id":i,"status":"failed" if i in failed else "passed"} for i in ids],
             "failure_details":[{"test":t.id(),"traceback":tb} for t,tb in result.failures+result.errors],
             "wall_seconds":time.time()-started,"observed":OBS}
    args.receipt.write_text(json.dumps(receipt,indent=1,sort_keys=True,allow_nan=False)+"\n")
    return 0 if result.wasSuccessful() else 1
if __name__=="__main__": sys.exit(main())
