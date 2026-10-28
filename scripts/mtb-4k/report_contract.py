#!/usr/bin/env python3
"""Strict revision-8 report contract over saved metadata and held-out scores only.

Port of package K's validator_core.py (Vast grm-mtb-4k-validators-kimi-20261028-r2), reconciled
with package W (RECONCILE.md section 5). It vendors nothing: the repo's own checker
(scripts/check-results.py) and metric audit (scripts/metric-audit/audit.py) are loaded in place,
after their SHA-256 match PINS, the blobs reviewed at 855bdba (unchanged through 3483065). A pin
that moves is a stop until the pins are reviewed again.

Changes from K: the per-fold lasso rule is the checker's (F1: a fold is called exactly when its
complexity is positive, and a fold whose record has no selected point has no calls); a fourth
cohort, matched2k, whose genome list is frozen and whose drug set follows the class gate (F6);
every class-gate receipt must equal the counts recomputed from the root's own phenotypes; every
report's input_fingerprint must equal fingerprint.py's recomputation from the root's files, as
well as the run manifest; and validate_manifest(collect_invalid=True) records a report that fails
the contract instead of stopping, so the comparator can exclude it (W's rule). Binding failures
(hashes, fingerprint, identity) always stop.
"""
from __future__ import annotations
import hashlib, importlib.util, json, math, re, statistics, sys
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path, PurePosixPath

LANE=Path(__file__).resolve().parent
SCRIPTS=LANE.parent
sys.path.insert(0,str(LANE))
import fingerprint as fpmod  # noqa: E402
SCHEMA,REVISION=2,8
SPECIES="mycobacterium_tuberculosis"
LEARNERS=("scm","cart","rf","gbt","lasso")
RANKED={"cart","rf","gbt","lasso"}
SCHEMES=("random","lineage")
CLASS_SKIP="too few labeled genomes after download, skipped"
FINGERPRINT=re.compile(r"[0-9a-f]{16}\Z")
SHA256=re.compile(r"[0-9a-f]{64}\Z")
COMMIT=re.compile(r"[0-9a-f]{40}\Z")
TOL=1e-12
# In-tree checker blobs this contract runs, as reviewed (RECONCILE.md section 3).
PINS={"check-results.py":"e42392a9a6473e10a59854a312a2891a783e697cf8ebd2304fcc4204631a3774",
 "metric-audit/audit.py":"1f5a0468949616d6fa6609c9ee48c8507224827726cad0fa9ca4d7e09936f42f",
 "metric-audit/summarize.py":"fd5fdfab1bedb8118c95eaadef19452fc89786f371a5d1757bd784844eb3295f",
 "result_scope.py":"3bcb2db34b23accd9186a595e6e0002cbd3bc292ab548442d8a53ba94b65ae2a"}
VALIDATOR_FILES=("report_contract.py","compare_cohorts.py","run_revision8_checks.py","fingerprint.py",
 "prepare_execution_roots.py","write_run_manifest.py")
DRUG_IDS={
 "amikacin":"mycobacterium_tuberculosis__amikacin","amoxicillin":"mycobacterium_tuberculosis__amoxicillin",
 "capreomycin":"mycobacterium_tuberculosis__capreomycin","ciprofloxacin":"mycobacterium_tuberculosis__ciprofloxacin",
 "clofazimine":"mycobacterium_tuberculosis__clofazimine","cycloserine":"mycobacterium_tuberculosis__cycloserine",
 "delamanid":"mycobacterium_tuberculosis__delamanid","ethambutol":"mycobacterium_tuberculosis__ethambutol",
 "ethionamide":"mycobacterium_tuberculosis__ethionamide","isoniazid":"mycobacterium_tuberculosis__isoniazid",
 "kanamycin":"mycobacterium_tuberculosis__kanamycin","levofloxacin":"mycobacterium_tuberculosis__levofloxacin",
 "linezolid":"mycobacterium_tuberculosis__linezolid","moxifloxacin":"mycobacterium_tuberculosis__moxifloxacin",
 "nicotinamide":"mycobacterium_tuberculosis__nicotinamide","ofloxacin":"mycobacterium_tuberculosis__ofloxacin",
 "para-aminosalicylic acid":"mycobacterium_tuberculosis__para_aminosalicylic_acid",
 "prothionamide":"mycobacterium_tuberculosis__prothionamide","pyrazinamide":"mycobacterium_tuberculosis__pyrazinamide",
 "rifabutin":"mycobacterium_tuberculosis__rifabutin","rifampin":"mycobacterium_tuberculosis__rifampin",
 "streptomycin":"mycobacterium_tuberculosis__streptomycin"}
ID_TO_DRUG={v:k for k,v in DRUG_IDS.items()}
ALL_DRUGS=tuple(DRUG_IDS)
ADDED_DRUGS=("amikacin","capreomycin","delamanid","ethambutol","ethionamide","isoniazid","kanamycin",
 "moxifloxacin","ofloxacin","prothionamide","pyrazinamide","rifabutin","rifampin","streptomycin")
ADDED_UNAVAILABLE={"amoxicillin":(0,113),"ciprofloxacin":(0,87),"clofazimine":(0,774),
 "cycloserine":(0,383),"levofloxacin":(774,1),"linezolid":(0,774),"nicotinamide":(0,0),
 "para-aminosalicylic acid":(5,1152)}
@dataclass(frozen=True)
class CohortRule:
 size:int
 genomes_sha256:str
 drugs:object  # tuple of drugs with reports; None: whatever the class gate gives on the root's genomes
 min_class:int=25
 unavailable:object=None  # frozen {drug:(R,S)} receipts; None: derived like drugs
RULES={
 "prefix":CohortRule(2000,"57af1c86c35fb0a4c62d3c0784a09804a919ba75941159d0c7b3cf56870aa2e3",ALL_DRUGS,unavailable={}),
 "added":CohortRule(2000,"a980d9e014fb90972a46b3fc01e65a9654deb58935d3ddf0f54018c19e7def85",ADDED_DRUGS,unavailable=ADDED_UNAVAILABLE),
 "full":CohortRule(4000,"80faef89571ecd2fd658204abd39fd4b55e5a86bfe8e4f379274162868f44eaf",ALL_DRUGS,unavailable={}),
 # prepare_execution_roots.matched_selection of the frozen 4,000 (measured read-only on Vast, 2026-10-28).
 "matched2k":CohortRule(2000,"212a21a5d2574e0d62d379fccddb0839c4f44f5a7cd36e43a558dc2c0926c3d1",None)}
for _rule in RULES.values():  # every production cohort accounts for all 22 datasets
 assert _rule.drugs is None or set(_rule.drugs)|set(_rule.unavailable)==set(ALL_DRUGS)
TOP_KEYS={"schema","evaluation","input_fingerprint","eval_config","id","species","species_slug","antibiotic",
 "drug_class","breakpoint","counts","lineage_threshold","catalog_groups","samples","splits","catalog","panels",
 "unsupported_assays","mic","models","features","seconds"}
METRIC_KEYS={"n","tp","fp","tn","fn","accuracy","balanced_accuracy","sensitivity","specificity","precision",
 "f1","mcc","vme","me","auroc","auprc"}
CONFIG_KEYS={"folds","inner_folds","scm_max_rules","cart_max_depth","forest","gbt","screened_features","lasso_lambdas","seed"}
FOREST_KEYS={"n_trees","subsample","max_depth","min_leaf"}
GBT_KEYS={"rounds","learning_rate","max_depth","lambda","min_child_weight","colsample","patience"}
FIT_KEYS={"path_points","path_unconverged","selected_point","selected_convergence","refit_convergence"}
CONVERGENCE_KEYS={"irls_steps","sweep_cap_hits","converged"}
MATRIX_BINS=("patterns.bin","support.bin","kmer_offsets.bin","kmers.bin")

class ValidationError(ValueError): pass
class BindingError(ValidationError):
 """A report or input that differs from what the run manifest binds: never only an exclusion."""
def fail(where,message): raise ValidationError(f"{where}: {message}")
def bind_fail(where,message): raise BindingError(f"{where}: {message}")
def file_sha256(path):
 h=hashlib.sha256()
 with Path(path).open("rb") as stream:
  for chunk in iter(lambda:stream.read(1<<20),b""): h.update(chunk)
 return h.hexdigest()
def pinned(rel):
 path=SCRIPTS/rel
 actual=file_sha256(path)
 if actual!=PINS[rel]: raise ValidationError(f"{path}: SHA-256 {actual} differs from the reviewed pin {PINS[rel]}")
 return path
def check_pins(): return {rel:str(pinned(rel)) for rel in PINS}
def validator_sha256(): return {name:file_sha256(LANE/name) for name in VALIDATOR_FILES if (LANE/name).is_file()}
def _load_module(name,path):
 spec=importlib.util.spec_from_file_location(name,path)
 if spec is None or spec.loader is None: raise RuntimeError(f"cannot load in-tree source {path}")
 module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module); return module
check_pins()
upstream_audit=_load_module("grm_in_tree_metric_audit",pinned("metric-audit/audit.py"))
upstream_checker=_load_module("grm_in_tree_result_checker",pinned("check-results.py"))
def expect_keys(value,required,where,optional=None):
 if not isinstance(value,dict): fail(where,"expected object")
 optional=optional or set(); missing=sorted(required-set(value)); unknown=sorted(set(value)-required-optional)
 if missing or unknown: fail(where,f"field mismatch; missing={missing}, unknown={unknown}")
def integer(value,where,minimum=0):
 if not isinstance(value,int) or isinstance(value,bool) or value<minimum: fail(where,f"expected integer >= {minimum}, got {value!r}")
 return value
def number(value,where):
 if not isinstance(value,(int,float)) or isinstance(value,bool) or not math.isfinite(value): fail(where,f"expected finite number, got {value!r}")
 return float(value)
def boolean(value,where):
 if not isinstance(value,bool): fail(where,f"expected boolean, got {value!r}")
 return value
def text(value,where,nonempty=True):
 if not isinstance(value,str) or (nonempty and not value): fail(where,"expected string")
 return value
def digest(value,where,pattern=SHA256):
 value=text(value,where)
 if not pattern.fullmatch(value): fail(where,f"malformed digest {value!r}")
 return value
def no_duplicate_objects(pairs):
 out={}
 for key,value in pairs:
  if key in out: raise ValidationError(f"duplicate JSON field {key!r}")
  out[key]=value
 return out
def reject_constant(value): raise ValidationError(f"nonfinite JSON constant {value}")
def check_finite_tree(value,where="json"):
 if isinstance(value,float) and not math.isfinite(value): fail(where,"nonfinite number")
 if isinstance(value,list):
  for i,item in enumerate(value): check_finite_tree(item,f"{where}[{i}]")
 elif isinstance(value,dict):
  for key,item in value.items(): check_finite_tree(item,f"{where}.{key}")
def read_json(path):
 raw=Path(path).read_bytes(); sha=hashlib.sha256(raw).hexdigest()
 try: value=json.loads(raw,object_pairs_hook=no_duplicate_objects,parse_constant=reject_constant)
 except (json.JSONDecodeError,UnicodeDecodeError) as exc: fail(str(path),f"malformed JSON: {exc}")
 check_finite_tree(value,str(path)); return value,sha
def safe_relative(root,value,where):
 value=text(value,where); p=PurePosixPath(value)
 if p.is_absolute() or ".." in p.parts or "." in p.parts or not p.parts: fail(where,f"unsafe relative path {value!r}")
 return root.joinpath(*p.parts)
def close_enough(a,b):
 if a is None or b is None: return a is b
 if not isinstance(a,(int,float)) or isinstance(a,bool) or not isinstance(b,(int,float)) or isinstance(b,bool): return False
 return abs(a-b)<=TOL
def compare_metrics(stored,computed,where):
 if computed is None:
  if stored is not None: fail(where,"metrics must be null when no isolate is called")
  return
 if stored is None: fail(where,"metrics are null despite called isolates")
 expect_keys(stored,METRIC_KEYS,where)
 for key,want in computed.items():
  got=stored.get(key)
  if isinstance(want,int):
   if not isinstance(got,int) or isinstance(got,bool) or got!=want: fail(f"{where}.{key}",f"stored {got!r}, recomputed integer {want!r}")
  elif want is None or got is None:
   if got!=want: fail(f"{where}.{key}",f"stored {got!r}, recomputed {want!r}")
  elif not isinstance(got,(int,float)) or isinstance(got,bool) or not close_enough(got,want):
   fail(f"{where}.{key}",f"stored {got!r}, recomputed {want!r}")

def betacf(a,b,x):
 tiny,eps=1e-300,3e-16; qab,qap,qam=a+b,a+1.0,a-1.0; c=1.0; d=1.0-qab*x/qap
 if abs(d)<tiny: d=tiny
 d=1.0/d; h=d
 for m in range(1,201):
  m2=2*m; aa=m*(b-m)*x/((qam+m2)*(a+m2)); d=1.0+aa*d
  if abs(d)<tiny: d=tiny
  c=1.0+aa/c
  if abs(c)<tiny: c=tiny
  d=1.0/d; h*=d*c; aa=-(a+m)*(qab+m)*x/((a+m2)*(qap+m2)); d=1.0+aa*d
  if abs(d)<tiny: d=tiny
  c=1.0+aa/c
  if abs(c)<tiny: c=tiny
  d=1.0/d; delta=d*c; h*=delta
  if abs(delta-1.0)<=eps: return h
 raise ArithmeticError(f"incomplete beta did not converge for a={a},b={b},x={x}")
def regularized_beta(x,a,b):
 if not 0<=x<=1 or a<=0 or b<=0: raise ValueError("invalid beta arguments")
 if x in (0.0,1.0): return x
 front=math.exp(math.lgamma(a+b)-math.lgamma(a)-math.lgamma(b)+a*math.log(x)+b*math.log1p(-x))
 if x<(a+1)/(a+b+2): return front*betacf(a,b,x)/a
 return 1.0-front*betacf(b,a,1.0-x)/b
def beta_ppf(p,a,b):
 if not 0<p<1 or a<=0 or b<=0: raise ValueError("invalid beta quantile")
 if p>0.5: return 1.0-beta_ppf(1.0-p,b,a)
 lo,hi=0.0,1.0
 for _ in range(160):
  mid=(lo+hi)/2
  if mid in (lo,hi): break
  if regularized_beta(mid,a,b)<p: lo=mid
  else: hi=mid
 return (lo+hi)/2
def beta_ppf_survival(q,a,b):
 if not 0<q<1 or a<=0 or b<=0: raise ValueError("invalid beta survival quantile")
 return 1.0-beta_ppf(q,b,a)
def clopper_pearson(errors,denominator,alpha=0.05):
 integer(errors,"CP errors"); integer(denominator,"CP denominator")
 if errors>denominator: fail("CP","errors exceed denominator")
 if denominator==0: return {"errors":0,"denominator":0,"confidence":1-alpha,"lower":None,"upper":None,"status":"undefined_no_called_denominator"}
 lower=0.0 if errors==0 else beta_ppf(alpha/2,errors,denominator-errors+1)
 upper=1.0 if errors==denominator else beta_ppf_survival(alpha/2,errors+1,denominator-errors)
 return {"errors":errors,"denominator":denominator,"confidence":1-alpha,"lower":lower,"upper":upper,"status":"defined"}
def complement(interval):
 """The success-rate interval (sensitivity from VME, specificity from ME) from the same computation."""
 if interval["status"]!="defined": return {"successes":0,"denominator":0,"confidence":interval["confidence"],"lower":None,"upper":None,"status":interval["status"]}
 return {"successes":interval["denominator"]-interval["errors"],"denominator":interval["denominator"],"confidence":interval["confidence"],
  "lower":1.0-interval["upper"],"upper":1.0-interval["lower"],"status":"defined"}

def validate_config(config,where="eval_config"):
 expect_keys(config,CONFIG_KEYS,where)
 for key in ("folds","inner_folds","scm_max_rules","cart_max_depth","screened_features","lasso_lambdas"): integer(config[key],f"{where}.{key}",1)
 integer(config["seed"],f"{where}.seed")
 expect_keys(config["forest"],FOREST_KEYS,f"{where}.forest")
 for key in ("n_trees","max_depth","min_leaf"): integer(config["forest"][key],f"{where}.forest.{key}",1)
 number(config["forest"]["subsample"],f"{where}.forest.subsample")
 expect_keys(config["gbt"],GBT_KEYS,f"{where}.gbt")
 for key in ("rounds","max_depth"): integer(config["gbt"][key],f"{where}.gbt.{key}",1)
 integer(config["gbt"]["patience"],f"{where}.gbt.patience")
 for key in ("learning_rate","lambda","min_child_weight","colsample"): number(config["gbt"][key],f"{where}.gbt.{key}")
 if config["folds"]<2 or config["inner_folds"]<2: fail(where,"outer and inner folds must be at least two")
 return config

def validate_convergence(status,where):
 expect_keys(status,CONVERGENCE_KEYS,where)
 steps=integer(status["irls_steps"],where+".irls_steps",1)
 hits=integer(status["sweep_cap_hits"],where+".sweep_cap_hits")
 converged=boolean(status["converged"],where+".converged")
 if hits>steps: fail(where,"sweep_cap_hits exceeds irls_steps")
 if converged and hits==steps: fail(where,"converged status cannot have every IRLS step hit the sweep cap")
def validate_lasso_fit(fit,where,lambda_count):
 expect_keys(fit,FIT_KEYS,where)
 for key in ("path_points","path_unconverged"): integer(fit[key],f"{where}.{key}")
 if fit["path_points"]>lambda_count: fail(where,"path_points exceeds eval_config.lasso_lambdas")
 if fit["selected_point"] is not None: integer(fit["selected_point"],where+".selected_point")
 for key in ("selected_convergence","refit_convergence"):
  if fit[key] is not None: validate_convergence(fit[key],f"{where}.{key}")
 try: return upstream_checker.check_lasso_fit(fit,(where,))
 except AssertionError as exc: fail(where,f"invalid revision-8 lasso status: {exc}")

@dataclass
class Sample:
 id:str
 resistant:bool
 lineage:int
@dataclass
class LearnerData:
 name:str
 scheme:str
 scores:dict
 folds:dict
 receipt:dict
@dataclass
class SplitData:
 scheme:str
 available:bool
 reason:object
 folds:dict
 learners:dict
 receipt:dict
@dataclass
class ReportData:
 id:str
 drug:str
 path:Path
 sha256:str
 input_fingerprint:str
 eval_config:dict
 lineage_threshold:float
 samples:dict
 sample_order:list
 splits:dict
 receipt:dict
@dataclass
class CohortData:
 name:str
 root:Path
 manifest_path:Path
 manifest:dict
 genomes:list
 lineage:dict
 reports:dict
 rule:CohortRule
 inputs:object=None
 invalid:dict=field(default_factory=dict)
 manifest_sha256:str=""

def error_intervals(metrics):
 if metrics is None: return {"vme":None,"me":None,"sensitivity":None,"specificity":None}
 vme=clopper_pearson(metrics["fn"],metrics["tp"]+metrics["fn"]); me=clopper_pearson(metrics["fp"],metrics["tn"]+metrics["fp"])
 return {"vme":vme,"me":me,"sensitivity":complement(vme),"specificity":complement(me)}

def validate_learner(row,scheme,folds,samples,k,catalog_present,lambda_count,where):
 expect_keys(row,{"learner","metrics","fold_balanced_accuracy","roc","scores","conformal","fold_complexity","seconds"},where,{"no_call","lasso_fits"})
 name=text(row["learner"],where+".learner")
 if name not in set(LEARNERS)|{"catalog"}: fail(where,f"unknown learner {name!r}")
 if name=="catalog" and not catalog_present: fail(where,"catalog row without catalog baseline")
 if name!="lasso" and "lasso_fits" in row: fail(where,"lasso status on a non-lasso row")
 number(row["seconds"],where+".seconds")
 if row["conformal"]!=[]: fail(where+".conformal","exploratory conformal field must be empty")
 scores_raw=row["scores"]
 if not isinstance(scores_raw,list) or len(scores_raw)!=len(samples): fail(where+".scores","score vector does not cover samples exactly")
 scores=[]
 for i,value in enumerate(scores_raw):
  if value is None: scores.append(None)
  else:
   value=number(value,f"{where}.scores[{i}]")
   if not 0<=value<=1: fail(f"{where}.scores[{i}]","probability score outside [0,1]")
   scores.append(value)
 truth=[s.resistant for s in samples]
 uncalled={"resistant":sum(v is None and t for v,t in zip(scores,truth)),"susceptible":sum(v is None and not t for v,t in zip(scores,truth))}
 if name=="catalog":
  if "no_call" in row or None in scores: fail(where,"catalog must call every isolate and omit no_call")
 else:
  if "no_call" not in row: fail(where+".no_call","schema-2 learner rows require explicit no_call counts")
  expect_keys(row["no_call"],{"resistant","susceptible"},where+".no_call")
  integer(row["no_call"]["resistant"],where+".no_call.resistant"); integer(row["no_call"]["susceptible"],where+".no_call.susceptible")
  if row["no_call"]!=uncalled: fail(where+".no_call",f"stored {row['no_call']}, recomputed {uncalled}")
  if name not in ("gbt","lasso") and any(v is None for v in scores): fail(where,f"{name} cannot have null-fit no-calls")
 complexity=row["fold_complexity"]
 if name=="catalog":
  if complexity!=[] or row["fold_balanced_accuracy"]!=[]: fail(where,"catalog fold fields must be empty")
 else:
  if not isinstance(complexity,list) or len(complexity)!=k: fail(where+".fold_complexity","expected one value per outer fold")
  for i,value in enumerate(complexity):
   number(value,f"{where}.fold_complexity[{i}]")
   if value<0: fail(f"{where}.fold_complexity[{i}]","negative model complexity")
 fold_receipts=[]; computed_fba=[]
 for fold in range(k):
  idx=[i for i,f in enumerate(folds) if f==fold]; held_scores=[scores[i] for i in idx]; held_truth=[truth[i] for i in idx]
  nulls=[v is None for v in held_scores]
  if any(nulls) and not all(nulls): fail(where,f"fold {fold} is called in part")
  if all(nulls):
   if name not in ("gbt","lasso") or complexity[fold]!=0: fail(where,f"fold {fold} null fit is inconsistent with learner/complexity")
  elif name in ("gbt","lasso") and complexity[fold]==0: fail(where,f"fold {fold} has zero complexity but called scores")
  m=upstream_audit.metrics(held_truth,held_scores,False); computed_fba.append(None if m is None else m["balanced_accuracy"])
  called=[i for i in idx if scores[i] is not None]
  cm=upstream_audit.metrics([truth[i] for i in called],[scores[i] for i in called],False)
  fold_receipts.append({"fold":fold,"held_out":len(idx),"resistant":sum(held_truth),"susceptible":len(idx)-sum(held_truth),
   "called":len(called),"called_resistant":sum(truth[i] for i in called),"called_susceptible":len(called)-sum(truth[i] for i in called),
   "no_call":{"resistant":sum(scores[i] is None and truth[i] for i in idx),"susceptible":sum(scores[i] is None and not truth[i] for i in idx)},
   "complexity":complexity[fold] if name!="catalog" else None,"metrics":cm})
 if name!="catalog":
  stored_fba=row["fold_balanced_accuracy"]
  if not isinstance(stored_fba,list) or len(stored_fba)!=k: fail(where+".fold_balanced_accuracy","expected one value per outer fold")
  for fold,(got,want) in enumerate(zip(stored_fba,computed_fba)):
   if not close_enough(got,want): fail(where+".fold_balanced_accuracy",f"fold {fold}: stored {got!r}, recomputed {want!r}")
 computed=upstream_audit.metrics(truth,scores,name in RANKED); compare_metrics(row["metrics"],computed,where+".metrics")
 ranked=name in RANKED
 expected_roc=upstream_audit.roc_curve([t for t,v in zip(truth,scores) if v is not None],[v for v in scores if v is not None]) if ranked and computed else []
 got_roc=row["roc"]
 if not isinstance(got_roc,list) or len(got_roc)!=len(expected_roc): fail(where+".roc","ROC point count differs")
 for i,(got,want) in enumerate(zip(got_roc,expected_roc)):
  if not isinstance(got,list) or len(got)!=2 or any(not close_enough(a,b) for a,b in zip(got,want)): fail(where+".roc",f"point {i} differs")
 if name=="scm" and any(v not in (0.0,1.0) for v in scores): fail(where,"SCM scores must be 0/1")
 fits=None
 if name=="lasso":
  fits=row.get("lasso_fits")
  if not isinstance(fits,list) or len(fits)!=k: fail(where+".lasso_fits","revision 8 requires one lasso status per outer fold")
  for fold,fit in enumerate(fits):
   has_refit=validate_lasso_fit(fit,f"{where}.lasso_fits[{fold}]",lambda_count)
   idx=[i for i,f in enumerate(folds) if f==fold]
   # F1, the checker's rule (check-results.py from 855bdba): a lasso fold is called exactly when its
   # complexity, the refit's weight count, is positive (the null-fit checks above), and a fold whose
   # record has no selected point (an empty path) has no calls. A non-empty path whose refit keeps
   # no weight is a valid null fold.
   if not has_refit and any(scores[i] is not None for i in idx): fail(where,f"fold {fold} empty-path lasso record has called scores")
 ids=[s.id for s in samples]
 receipt={"learner":name,"scheme":scheme,"metrics":computed,"no_call":uncalled if name!="catalog" else None,
  "called":sum(v is not None for v in scores),"uncalled":sum(v is None for v in scores),
  "no_call_rate":upstream_audit.ratio(sum(v is None for v in scores),len(scores)),
  "fold_balanced_accuracy":computed_fba if name!="catalog" else [],"folds":fold_receipts,
  "clopper_pearson_95":error_intervals(computed),"lasso_fits":fits}
 return LearnerData(name,scheme,dict(zip(ids,scores)),dict(zip(ids,folds)),receipt)

def validate_report(path,expected_id,expected,cohort):
 report,sha=read_json(path)
 if sha!=expected["sha256"]: bind_fail(str(path),f"report SHA-256 differs from run manifest: {sha} != {expected['sha256']}")
 if not isinstance(report,dict) or not isinstance(report.get("evaluation"),dict): fail(str(path),"missing evaluation object")
 if not isinstance(report["evaluation"].get("revision"),int) or isinstance(report["evaluation"].get("revision"),bool) or report["evaluation"].get("revision")!=REVISION: fail(str(path),f"matched controls require evaluation revision 8, got {report['evaluation'].get('revision')!r}")
 expect_keys(report,TOP_KEYS,str(path))
 if not isinstance(report["schema"],int) or isinstance(report["schema"],bool) or report["schema"]!=SCHEMA: fail(str(path),f"matched controls require schema 2, got {report['schema']!r}")
 expect_keys(report["evaluation"],{"revision","valid","reasons","scope"},str(path)+".evaluation")
 if report["evaluation"]["valid"] is not True or report["evaluation"]["reasons"]!=[]: fail(str(path),"only evaluation.valid=true with empty reasons is accepted")
 scope=text(report["evaluation"]["scope"],str(path)+".evaluation.scope")
 if "not external or clinical validation" not in scope: fail(str(path),"evaluation scope lacks internal-CV/non-clinical limitation")
 if report["id"]!=expected_id or report["species_slug"]!=SPECIES: bind_fail(str(path),"report identity/species differs from manifest")
 drug=ID_TO_DRUG.get(expected_id)
 if drug is None or report["antibiotic"]!=drug: bind_fail(str(path),"report antibiotic differs from curated identifier")
 fingerprint=digest(report["input_fingerprint"],str(path)+".input_fingerprint",FINGERPRINT)
 if fingerprint!=expected["input_fingerprint"]: bind_fail(str(path),"input fingerprint differs from run manifest (stale or mismatched input)")
 config=validate_config(report["eval_config"],str(path)+".eval_config")
 if config!=cohort.manifest["eval_config"]: bind_fail(str(path),"eval_config differs from run manifest")
 recomputed=cohort.inputs.fingerprint(expected_id,config)
 if fingerprint!=recomputed: bind_fail(str(path),f"input fingerprint {fingerprint} differs from {recomputed} recomputed from the root's files (stale report or changed inputs)")
 threshold=number(report["lineage_threshold"],str(path)+".lineage_threshold")
 if threshold!=cohort.manifest["lineage_threshold"]: bind_fail(str(path),"lineage threshold differs from manifest/structure")
 expect_keys(report["counts"],{"samples","resistant","susceptible","labels_from_mic","labels_from_mic_excluded","features","patterns","kmers","k","lineages","min_support","min_class"},str(path)+".counts")
 counts=report["counts"]
 for key in counts: integer(counts[key],str(path)+f".counts.{key}")
 if counts["k"]!=31 or counts["min_class"]!=cohort.rule.min_class: fail(str(path),"counts k or class minimum differs from cohort rule")
 if not isinstance(report["samples"],list): fail(str(path),"samples must be a list")
 if not isinstance(report["catalog_groups"],list) or any(not isinstance(x,str) for x in report["catalog_groups"]): fail(str(path),"catalog_groups must be strings")
 samples=[]
 for i,raw in enumerate(report["samples"]):
  where=f"{path}.samples[{i}]"
  expect_keys(raw,{"id","resistant","from_mic","mic","lineage","country","year","host","mlst","groups","catalog_resistant"},where)
  sid=text(raw["id"],where+".id"); resistant=boolean(raw["resistant"],where+".resistant")
  if raw["from_mic"] is not False: fail(where,"MIC-inferred labels are excluded")
  lineage=integer(raw["lineage"],where+".lineage")
  if not isinstance(raw["groups"],list) or any(integer(g,where+".groups")>=len(report["catalog_groups"]) for g in raw["groups"]): fail(where+".groups","invalid catalog group index")
  if raw["mic"] is not None:
   expect_keys(raw["mic"],{"log2","censor"},where+".mic"); integer(raw["mic"]["log2"],where+".mic.log2",-1000000)
   if raw["mic"]["censor"] not in ("Exact","AtMost","Above"): fail(where+".mic.censor","unknown censor status")
  if raw["catalog_resistant"] is not None: boolean(raw["catalog_resistant"],where+".catalog_resistant")
  if sid not in cohort.genomes: fail(where,f"genome {sid!r} is outside frozen cohort")
  if cohort.lineage[sid]!=lineage: fail(where,f"genome {sid!r} lineage differs from structure.json")
  samples.append(Sample(sid,resistant,lineage))
 ids=[s.id for s in samples]
 if len(set(ids))!=len(ids): fail(str(path),"duplicate sample genome IDs")
 n=len(samples); nr=sum(s.resistant for s in samples); ns=n-nr
 if (counts["samples"],counts["resistant"],counts["susceptible"],counts["labels_from_mic"])!=(n,nr,ns,0): fail(str(path),f"label/count mismatch: recomputed n={n},R={nr},S={ns}")
 if nr<cohort.rule.min_class or ns<cohort.rule.min_class: fail(str(path),"evaluated class count below required gate")
 if counts["lineages"]!=len({s.lineage for s in samples}): fail(str(path),"counts.lineages differs from samples")
 if counts["min_support"]!=max(n//200,3): fail(str(path),"dataset projection min_support differs from committed policy")
 if not isinstance(report["splits"],list) or len(report["splits"])!=2: fail(str(path),"exactly random and lineage splits are required")
 splits={}
 for si,split in enumerate(report["splits"]):
  where=f"{path}.splits[{si}]"; expect_keys(split,{"scheme","folds","n_groups","learners"},where,{"test_fold","unavailable_reason"})
  scheme=split["scheme"]
  if scheme not in SCHEMES or scheme in splits: fail(where,f"unknown or duplicate split scheme {scheme!r}")
  k=integer(split["folds"],where+".folds",2)
  integer(split["n_groups"],where+".n_groups",1)
  if k!=config["folds"] or split["n_groups"]!=counts["lineages"]: fail(where,"fold or lineage group count differs")
  if "unavailable_reason" in split:
   if scheme!="lineage": fail(where,"random split cannot be unavailable")
   reason=text(split["unavailable_reason"],where+".unavailable_reason")
   if split["learners"]!=[] or "test_fold" in split: fail(where,"unavailable split must have no learners or folds")
   splits[scheme]=SplitData(scheme,False,reason,{},{},{"scheme":scheme,"status":"unavailable","unavailable_reason":reason,"learners":[],"five_learner_medians":None})
   continue
  folds=split.get("test_fold")
  if not isinstance(folds,list) or len(folds)!=n: fail(where,"available split must save one test fold per sample")
  for i,fold in enumerate(folds):
   integer(fold,f"{where}.test_fold[{i}]")
   if fold>=k: fail(where,"test fold outside configured range")
  for fold in range(k):
   held=[s.resistant for s,f in zip(samples,folds) if f==fold]
   if not held or all(held) or not any(held): fail(where,f"fold {fold} lacks a class")
  if scheme=="lineage":
   groups=defaultdict(set)
   for s,fold in zip(samples,folds): groups[s.lineage].add(fold)
   if any(len(v)!=1 for v in groups.values()): fail(where,"a lineage spans held-out folds")
  expected_names=list(LEARNERS)+(["catalog"] if report["catalog"] is not None else [])
  learners_raw=split["learners"]
  if not isinstance(learners_raw,list) or [r.get("learner") for r in learners_raw]!=expected_names: fail(where,f"learner rows must be exactly {expected_names}")
  learners={}
  for li,row in enumerate(learners_raw):
   learner=validate_learner(row,scheme,folds,samples,k,report["catalog"] is not None,config["lasso_lambdas"],f"{where}.learners[{li}]"); learners[learner.name]=learner
  receipt={"scheme":scheme,"status":"available","folds":k,"n_groups":split["n_groups"],"learners":[learners[name].receipt for name in expected_names]}
  receipt["five_learner_medians"]=learner_medians([learners[name].receipt for name in LEARNERS])
  splits[scheme]=SplitData(scheme,True,None,dict(zip(ids,folds)),learners,receipt)
 if set(splits)!=set(SCHEMES): fail(str(path),"random and lineage analyses must be separate and complete")
 models=report["models"]
 if not isinstance(models,dict) or "lasso_fit" not in models: fail(str(path),"revision 8 requires final models.lasso_fit")
 has_refit=validate_lasso_fit(models["lasso_fit"],str(path)+".models.lasso_fit",config["lasso_lambdas"])
 if (models.get("lasso") is not None)!=has_refit: fail(str(path),"final lasso model inconsistent with solver evidence")
 validate_catalog(report,samples,str(path))
 validate_panels(report,nr,ns,str(path))
 receipt={"id":expected_id,"antibiotic":drug,"sha256":sha,"input_fingerprint":fingerprint,"input_fingerprint_recomputed":True,
  "schema":SCHEMA,"evaluation_revision":REVISION,"valid":True,"counts":{"samples":n,"resistant":nr,"susceptible":ns,"lineages":counts["lineages"],
  "labels_from_mic":0,"min_class":counts["min_class"],"min_support":counts["min_support"]},"lineage_threshold":threshold,
  "eval_config":config,"splits":{scheme:splits[scheme].receipt for scheme in SCHEMES}}
 return ReportData(expected_id,drug,path,sha,fingerprint,config,threshold,{s.id:s for s in samples},ids,splits,receipt)

def validate_catalog(report,samples,where):
 catalog=report["catalog"]
 if catalog is None:
  if any(s_raw["catalog_resistant"] is not None for s_raw in report["samples"]): fail(where,"catalog flags without catalog")
  return
 expect_keys(catalog,{"version","rule","metrics","groups"},where+".catalog")
 flags=[s["catalog_resistant"] for s in report["samples"]]
 if any(v is None for v in flags): fail(where+".catalog","catalog baseline requires a flag per sample")
 truth=[s.resistant for s in samples]; scores=[1.0 if v else 0.0 for v in flags]
 computed=upstream_audit.metrics(truth,scores,False); compare_metrics(catalog["metrics"],computed,where+".catalog.metrics")
 for split in report["splits"]:
  if "unavailable_reason" in split: continue
  row=next(r for r in split["learners"] if r["learner"]=="catalog")
  if row["scores"]!=scores: fail(where+".catalog","catalog learner scores differ from sample flags")
 table=defaultdict(lambda:[0,0])
 for raw,t in zip(report["samples"],truth):
  for group_index in raw["groups"]: table[report["catalog_groups"][group_index]][0 if t else 1]+=1
 stored={row["group"]:[row["n_resistant"],row["n_susceptible"]] for row in catalog["groups"]}
 if len(stored)!=len(catalog["groups"]): fail(where+".catalog.groups","duplicate catalog group rows")
 if stored!=dict(table): fail(where+".catalog.groups","catalog carrier counts differ")

def validate_panels(report,nr,ns,where):
 panel_ids=set()
 for i,panel in enumerate(report["panels"]):
  w=f"{where}.panels[{i}]"; expect_keys(panel,{"id","name","mirrors","targets","metrics","detected"},w,{"no_call"})
  pid=text(panel["id"],w+".id")
  if pid in panel_ids: fail(w,"duplicate panel ID")
  panel_ids.add(pid); is_probe=pid in upstream_audit.PROBE_IDS
  if not isinstance(panel["metrics"],dict): fail(w+".metrics","expected object")
  for key in ("tp","fp","tn","fn"): integer(panel["metrics"][key],w+".metrics."+key)
  m=panel["metrics"]
  if is_probe:
   expect_keys(panel["no_call"],{"resistant","susceptible"},w+".no_call"); nc=panel["no_call"]
   integer(nc["resistant"],w+".no_call.resistant"); integer(nc["susceptible"],w+".no_call.susceptible")
   if m["tp"]+m["fn"]+nc["resistant"]!=nr or m["fp"]+m["tn"]+nc["susceptible"]!=ns: fail(w,"probe calls/no-calls do not partition class denominators")
  elif "no_call" in panel or m["n"]!=nr+ns: fail(w,"catalog panel no_call/n arithmetic is inconsistent")
  derived=upstream_audit.from_counts(m["tp"],m["fp"],m["tn"],m["fn"]); derived["auroc"]=None; derived["auprc"]=None
  compare_metrics(m,derived,w+".metrics")
 for i,assay in enumerate(report["unsupported_assays"]):
  w=f"{where}.unsupported_assays[{i}]"; expect_keys(assay,{"id","name","mirrors","reason"},w); text(assay["reason"],w+".reason")
  if assay["id"] in panel_ids: fail(w,"unsupported assay also has a scored panel")

def learner_medians(rows):
 if len(rows)!=5 or {r["learner"] for r in rows}!=set(LEARNERS): fail("five-learner median","prespecified learner set incomplete")
 out={}
 for field_name in ("sensitivity","specificity","balanced_accuracy","vme","me"):
  values=[None if r["metrics"] is None else r["metrics"][field_name] for r in rows]; present=[v for v in values if v is not None]
  out[field_name]={"median":statistics.median(values) if len(present)==5 else None,"available_only_median":statistics.median(present) if present else None,
   "prespecified_learners":5,"nonnull_learners":len(present),"null_learners":[r["learner"] for r,v in zip(rows,values) if v is None]}
 out["no_call_rate"]={"median":statistics.median([r["no_call_rate"] for r in rows]),"available_only_median":statistics.median([r["no_call_rate"] for r in rows]),
  "prespecified_learners":5,"nonnull_learners":5,"null_learners":[]}
 return out

def required_inputs(gm,mm):
 base=f"data/species/{SPECIES}"
 return {gm["path"],mm["path"],f"{base}/structure.json",f"{base}/counted.json","data/curated/datasets.json","data/report/index.json",
  "data/curated/phenotypes.jsonl","data/catalog/version.txt"}|{f"{base}/matrix/{name}" for name in MATRIX_BINS}

def validate_manifest(path,rules,collect_invalid=False):
 path=Path(path)
 manifest,manifest_sha=read_json(path)
 expect_keys(manifest,{"schema","cohort","species_slug","source_commit","binary_sha256","evaluation_revision","report_schema","eval_config","lineage_threshold","genomes","matrix","inputs","reports","class_gate_unavailable"},str(path))
 if manifest["schema"]!="grm-mtb-4k-run-manifest-v1": fail(str(path),"unknown run-manifest schema")
 name=manifest["cohort"]
 if name not in rules: fail(str(path),f"unknown cohort {name!r}")
 rule=rules[name]
 if manifest["species_slug"]!=SPECIES or manifest["evaluation_revision"]!=REVISION or manifest["report_schema"]!=SCHEMA: fail(str(path),"wrong species/revision/schema")
 digest(manifest["source_commit"],str(path)+".source_commit",COMMIT); digest(manifest["binary_sha256"],str(path)+".binary_sha256")
 validate_config(manifest["eval_config"],str(path)+".eval_config"); number(manifest["lineage_threshold"],str(path)+".lineage_threshold")
 root=path.parent; cohort=CohortData(name,root,path,manifest,[],{},{},rule,manifest_sha256=manifest_sha)
 gm,mm=manifest["genomes"],manifest["matrix"]
 expect_keys(gm,{"path","sha256"},str(path)+".genomes"); expect_keys(mm,{"path","sha256","k","min_support"},str(path)+".matrix")
 if mm["k"]!=31 or mm["min_support"]!=2: fail(str(path),"matched matrix must use k=31 and min_support=2")
 if not isinstance(manifest["inputs"],list): fail(str(path),"inputs must be a list")
 input_paths=set()
 for i,entry in enumerate(manifest["inputs"]):
  expect_keys(entry,{"path","sha256"},f"{path}.inputs[{i}]"); p=safe_relative(root,entry["path"],f"{path}.inputs[{i}]")
  if entry["path"] in input_paths: fail(str(path),f"duplicate input path {entry['path']}")
  input_paths.add(entry["path"]); actual=file_sha256(p)
  if actual!=digest(entry["sha256"],f"{path}.inputs[{i}].sha256"): bind_fail(str(p),f"stale input SHA-256: {actual} != {entry['sha256']}")
 required=required_inputs(gm,mm)
 if not required<=input_paths: fail(str(path),f"manifest inputs omit {sorted(required-input_paths)}")
 genomes_path=safe_relative(root,gm["path"],str(path)+".genomes.path")
 if file_sha256(genomes_path)!=gm["sha256"] or gm["sha256"]!=rule.genomes_sha256: fail(str(genomes_path),"cohort genome-list hash differs from frozen rule")
 text_data=genomes_path.read_text(encoding="utf-8")
 if not text_data.endswith("\n"): fail(str(genomes_path),"genome list lacks final newline")
 genomes=text_data.splitlines()
 if len(genomes)!=rule.size or len(set(genomes))!=len(genomes) or any(not g or g!=g.strip() for g in genomes): fail(str(genomes_path),"wrong/duplicate/blank genome IDs")
 cohort.genomes=genomes
 matrix_path=safe_relative(root,mm["path"],str(path)+".matrix.path"); matrix,matrix_sha=read_json(matrix_path)
 if matrix_sha!=mm["sha256"]: fail(str(matrix_path),"matrix metadata digest differs from manifest")
 if matrix.get("k")!=31 or matrix.get("min_support")!=2 or matrix.get("genomes")!=genomes: fail(str(matrix_path),"matrix k/min_support/genome order differs")
 digests=matrix.get("kset_digests")
 if not isinstance(digests,list) or len(digests)!=len(genomes) or any(not isinstance(x,str) or not FINGERPRINT.fullmatch(x) for x in digests): fail(str(matrix_path),"incomplete matrix kset digests")
 structure_path=root/"data/species"/SPECIES/"structure.json"; structure,_=read_json(structure_path)
 if structure.get("genomes")!=genomes or structure.get("k")!=31 or structure.get("kset_digests")!=digests: fail(str(structure_path),"structure is stale relative to matrix")
 lineages=structure.get("lineages")
 if not isinstance(lineages,dict): fail(str(structure_path),"missing lineages")
 labels=lineages.get("labels")
 if not isinstance(labels,list) or len(labels)!=len(genomes): fail(str(structure_path),"lineage labels do not cover cohort")
 for i,value in enumerate(labels): integer(value,f"{structure_path}.lineages.labels[{i}]")
 n_clusters=len(set(labels))
 if lineages.get("n_clusters")!=n_clusters or lineages.get("threshold")!=manifest["lineage_threshold"]: fail(str(structure_path),"lineage threshold/cluster count differs")
 largest=max(Counter(labels).values())/len(labels)
 if not close_enough(lineages.get("largest_share"),largest): fail(str(structure_path),"lineage largest_share differs")
 cohort.lineage=dict(zip(genomes,labels))
 counted_path=root/"data/species"/SPECIES/"counted.json"; counted,_=read_json(counted_path)
 if counted.get("k")!=31 or [g.get("genome_id") for g in counted.get("genomes",[])]!=genomes: fail(str(counted_path),"counted genomes/k differ")
 if [g.get("kset_digest") for g in counted["genomes"]]!=digests or counted.get("missing_fasta")!=[] or counted.get("failed")!=[] or counted.get("pending",[])!=[]: fail(str(counted_path),"counted digests differ or pending/failed genomes remain")
 curated_path=root/"data/curated/datasets.json"; curated,_=read_json(curated_path)
 if not isinstance(curated,list): fail(str(curated_path),"expected dataset list")
 curated_by_id={row.get("id"):row for row in curated if isinstance(row,dict) and row.get("species_slug")==SPECIES}
 if set(curated_by_id)!=set(DRUG_IDS.values()) or len(curated_by_id)!=sum(1 for row in curated if isinstance(row,dict) and row.get("species_slug")==SPECIES): fail(str(curated_path),"curated set must contain exact 22 IDs")
 for drug,dataset_id in DRUG_IDS.items():
  row=curated_by_id[dataset_id]
  if row.get("species_slug")!=SPECIES or row.get("antibiotic")!=drug or row.get("min_class")!=rule.min_class: fail(str(curated_path),f"curated dataset mismatch for {drug}")
 try: cohort.inputs=fpmod.Inputs(root/"data",SPECIES)
 except (OSError,KeyError,ValueError) as exc: bind_fail(str(root),f"fingerprint inputs unreadable: {exc}")
 # Class-gate counts recomputed from the root's own phenotypes on its matrix genomes.
 gate={drug:cohort.inputs.class_counts(curated_by_id[did]) for drug,did in DRUG_IDS.items()}
 if not isinstance(manifest["reports"],dict): fail(str(path),"reports must be an object")
 rows=manifest["class_gate_unavailable"]
 if not isinstance(rows,list): fail(str(path),"class_gate_unavailable must be a list")
 eligible={d for d,(r,s) in gate.items() if r>=rule.min_class and s>=rule.min_class}
 expected_drugs=set(rule.drugs) if rule.drugs is not None else eligible
 if expected_drugs!=eligible: fail(str(path),f"the root's phenotypes give eligible drugs {sorted(eligible)}, not the frozen {sorted(expected_drugs)}")
 expected_ids={DRUG_IDS[d] for d in expected_drugs}
 if set(manifest["reports"])!=expected_ids: fail(str(path),f"report set mismatch: expected {len(expected_ids)}, got {len(manifest['reports'])}")
 expected_unavailable=set(rule.unavailable) if rule.unavailable is not None else set(ALL_DRUGS)-expected_drugs
 seen=set()
 for i,row in enumerate(rows):
  where=f"{path}.class_gate_unavailable[{i}]"; expect_keys(row,{"id","antibiotic","min_class","resistant","susceptible","reason"},where)
  drug=row["antibiotic"]
  if drug not in expected_unavailable or drug in seen: fail(where,f"unexpected or duplicate unavailable drug {drug!r}")
  seen.add(drug)
  for key in ("min_class","resistant","susceptible"): integer(row[key],where+"."+key)
  r,s=rule.unavailable[drug] if rule.unavailable is not None else gate[drug]
  if (row["id"],row["min_class"],row["resistant"],row["susceptible"],row["reason"])!=(DRUG_IDS[drug],25,r,s,CLASS_SKIP): fail(where,"class-gate receipt differs from frozen evidence")
  if (r,s)!=gate[drug]: fail(where,f"class-gate receipt {(r,s)} differs from the root's phenotypes {gate[drug]}")
  if row["resistant"]>=25 and row["susceptible"]>=25: fail(where,"unavailable receipt unexpectedly passes gate")
 if seen!=expected_unavailable: fail(str(path),f"missing unavailable receipts: {sorted(expected_unavailable-seen)}")
 index_path=root/"data/report/index.json"; index,_=read_json(index_path)
 if index.get("schema")!=SCHEMA or index.get("evaluation_revision")!=REVISION: fail(str(index_path),"index has wrong schema/revision")
 index_rows=index.get("datasets")
 if not isinstance(index_rows,list) or len(index_rows)!=len(expected_ids) or {r.get("id") for r in index_rows}!=expected_ids: fail(str(index_path),"index dataset set differs or contains duplicates")
 for dataset_id in sorted(expected_ids):
  entry=manifest["reports"][dataset_id]; expect_keys(entry,{"path","sha256","input_fingerprint"},f"{path}.reports.{dataset_id}")
  conventional=f"data/species/{SPECIES}/results/{dataset_id}.json"
  if entry["path"]!=conventional: fail(str(path),f"noncanonical report path for {dataset_id}")
  digest(entry["sha256"],f"{path}.reports.{dataset_id}.sha256"); digest(entry["input_fingerprint"],f"{path}.reports.{dataset_id}.input_fingerprint",FINGERPRINT)
  try:
   report=validate_report(safe_relative(root,entry["path"],str(path)),dataset_id,entry,cohort)
   if (report.receipt["counts"]["resistant"],report.receipt["counts"]["susceptible"])!=gate[report.drug]: fail(str(report.path),f"report class counts differ from the root's phenotypes {gate[report.drug]}")
   cohort.reports[dataset_id]=report
  except BindingError: raise
  except ValidationError as exc:
   if not collect_invalid: raise
   cohort.invalid[dataset_id]=str(exc)
 return cohort

def compare_lineage_partitions(prefix,full):
 full_set=set(full.genomes); common=[g for g in prefix.genomes if g in full_set]
 if len(common)!=len(prefix.genomes): fail("lineage comparison","prefix genomes are not all present in full cohort")
 pairs=Counter((prefix.lineage[g],full.lineage[g]) for g in common); p_to_f=defaultdict(set); f_to_p=defaultdict(set)
 for p,f in pairs: p_to_f[p].add(f); f_to_p[f].add(p)
 same=all(len(v)==1 for v in p_to_f.values()) and all(len(v)==1 for v in f_to_p.values())
 rows=[{"prefix_group":p,"full_groups":sorted(p_to_f[p]),"common_members":sum(n for (pp,_),n in pairs.items() if pp==p)} for p in sorted(p_to_f)]
 full_sizes=Counter(full.lineage[g] for g in full.genomes)
 return {"common_genomes":len(common),"changed_raw_group_ids":sum(prefix.lineage[g]!=full.lineage[g] for g in common),
  "same_common_genome_partition":same,"split_prefix_groups":{str(p):sorted(v) for p,v in p_to_f.items() if len(v)>1},
  "merged_full_groups_on_common":{str(f):sorted(v) for f,v in f_to_p.items() if len(v)>1},"prefix_group_mapping":rows,
  "full_group_sizes_including_added":{str(f):full_sizes[f] for f in sorted(f_to_p)},
  "note":"Partition relations use common genomes; full-group sizes also show added members. Raw renumbering alone is not a partition change."}

def cohort_public_receipt(cohort):
 return {"cohort":cohort.name,"root":str(cohort.root),"run_manifest":str(cohort.manifest_path),"run_manifest_sha256":cohort.manifest_sha256,
  "source_commit":cohort.manifest["source_commit"],"binary_sha256":cohort.manifest["binary_sha256"],
  "genomes":len(cohort.genomes),"genomes_sha256":cohort.manifest["genomes"]["sha256"],"matrix":{"k":31,"min_support":2},
  "eval_config":cohort.manifest["eval_config"],"lineage_threshold":cohort.manifest["lineage_threshold"],
  "eligible_reports":len(cohort.reports)+len(cohort.invalid),"valid_reports":len(cohort.reports),
  "eligible_drugs":sorted(ID_TO_DRUG[i] for i in set(cohort.reports)|set(cohort.invalid)),"invalid_reports":dict(sorted(cohort.invalid.items())),
  "class_gate_unavailable":cohort.manifest["class_gate_unavailable"],"reports":[cohort.reports[i].receipt for i in sorted(cohort.reports)]}
