"""Small synthetic revision-8 roots; no biological input or learned result.

Port of package K's tests/fixtures.py. Differences: every report's random and lineage test_fold
comes from the committed foldgen (FOLDGEN, built from scripts/metric-audit/foldgen), as `grm learn`
assigns them and as the metric audit rebuilds them; a fourth cohort, matched2k; roots carry what
fingerprint.py reads (matrix files, curated phenotypes, the catalog version, full gene calls), and
every report's input_fingerprint is recomputed from them; and the added cohort's eight class-gate
receipts follow the fixture's own labels (FIXTURE_ADDED_UNAVAILABLE), which the full cohort's labels
on the added genomes reproduce. FakeFoldGen, K's index-fold double, stays only for injected-failure
tests: on these fixtures it disagrees with every saved fold.
"""
import atexit,copy,hashlib,json,math,os,struct,subprocess,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import report_contract as vc
import fingerprint as fpmod
import prepare_execution_roots as prep

CONFIG={"folds":5,"inner_folds":3,"scm_max_rules":6,"cart_max_depth":8,
 "forest":{"n_trees":200,"subsample":0.632,"max_depth":30,"min_leaf":1},
 "gbt":{"rounds":400,"learning_rate":0.1,"max_depth":4,"lambda":1.0,"min_child_weight":1.0,"colsample":0.5,"patience":25},
 "screened_features":4000,"lasso_lambdas":30,"seed":42}
SCOPE="Internal cross-validation on synthetic labels. This is not external or clinical validation."
SOURCE="1"*40
SPECIES_NAME="Synthetic tuberculosis"
CATALOG_VERSION="synthetic"
COHORTS=("prefix","added","full","matched2k")
# The added genomes' labels for the eight drugs the added cohort cannot evaluate, in the frozen shape
# (no resistant isolate; 5 resistant; nearly all resistant; no label at all).
FIXTURE_ADDED_UNAVAILABLE={"amoxicillin":(0,60),"ciprofloxacin":(0,60),"clofazimine":(0,60),"cycloserine":(0,60),
 "levofloxacin":(59,1),"linezolid":(0,60),"nicotinamide":(0,0),"para-aminosalicylic acid":(5,55)}

def foldgen_path():
 path=os.environ.get("FOLDGEN")
 if not path or not Path(path).is_file() or not os.access(path,os.X_OK):
  raise SystemExit("set FOLDGEN to the foldgen binary built from scripts/metric-audit/foldgen (cargo build --release --locked)")
 return path

def foldgen_sha256():
 """The --expected-foldgen-sha256 every runner and comparator call in these tests passes (I4)."""
 return hashlib.sha256(Path(foldgen_path()).read_bytes()).hexdigest()

class FoldGen:
 """Client for the committed foldgen, with metric-audit/audit.py's protocol."""
 def __init__(self,path=None):
  self.p=subprocess.Popen([path or foldgen_path()],stdin=subprocess.PIPE,stdout=subprocess.PIPE,text=True)
 def ask(self,line):
  self.p.stdin.write(line+"\n"); self.p.stdin.flush()
  return [int(v) for v in self.p.stdout.readline().strip().split(",")]
 def stratified(self,y,k=5,seed=42): return self.ask("S %d %d %s"%(k,seed,"".join("1" if v else "0" for v in y)))
 def grouped(self,y,groups,k=5,seed=42): return self.ask("G %d %d %s %s"%(k,seed,"".join("1" if v else "0" for v in y),",".join(map(str,groups))))
 def close(self):
  if self.p.poll() is None:
   self.p.stdin.close(); self.p.wait(timeout=10)
_FG=None
def foldgen():
 global _FG
 if _FG is None:
  _FG=FoldGen(); atexit.register(_FG.close)
 return _FG

def ratio(a,b): return a/b if b else 0.0

def metrics_direct(truth,scores,ranked):
 pairs=[(t,s) for t,s in zip(truth,scores) if s is not None]
 if not pairs: return None
 tp=sum(t and s>=0.5 for t,s in pairs); fp=sum(not t and s>=0.5 for t,s in pairs)
 tn=sum(not t and s<0.5 for t,s in pairs); fn=sum(t and s<0.5 for t,s in pairs)
 sens,spec,prec=ratio(tp,tp+fn),ratio(tn,tn+fp),ratio(tp,tp+fp)
 denom=math.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)); pos=[s for t,s in pairs if t]; neg=[s for t,s in pairs if not t]
 auroc=auprc=None
 if ranked and pos and neg: auroc=sum(1.0 if p>q else 0.5 if p==q else 0.0 for p in pos for q in neg)/(len(pos)*len(neg))
 if ranked and pos:
  auprc,before=0.0,0.0
  for cut in sorted({s for _,s in pairs},reverse=True):
   chosen=[t for t,s in pairs if s>=cut]; recall=sum(chosen)/len(pos)
   auprc+=(recall-before)*sum(chosen)/len(chosen); before=recall
 return {"n":len(pairs),"tp":tp,"fp":fp,"tn":tn,"fn":fn,"accuracy":ratio(tp+tn,len(pairs)),
  "balanced_accuracy":(sens+spec)/2,"sensitivity":sens,"specificity":spec,"precision":prec,
  "f1":0.0 if prec+sens==0 else 2*prec*sens/(prec+sens),"mcc":0.0 if denom==0 else (tp*tn-fp*fn)/denom,
  "vme":ratio(fn,tp+fn),"me":ratio(fp,fp+tn),"auroc":auroc,"auprc":auprc}

def roc_direct(truth,scores):
 pairs=[(t,s) for t,s in zip(truth,scores) if s is not None]; np=sum(t for t,_ in pairs); nn=len(pairs)-np
 if not np or not nn: return []
 points=[[0.0,0.0]]
 for cut in sorted({s for _,s in pairs},reverse=True):
  chosen=[t for t,s in pairs if s>=cut]; points.append([(len(chosen)-sum(chosen))/nn,sum(chosen)/np])
 return points

def status(): return {"irls_steps":3,"sweep_cap_hits":0,"converged":True}
def fit(empty=False):
 if empty: return {"path_points":0,"path_unconverged":0,"selected_point":None,"selected_convergence":None,"refit_convergence":None}
 return {"path_points":30,"path_unconverged":0,"selected_point":2,"selected_convergence":status(),"refit_convergence":status()}
def null_fit(scheme,learner,fold): return (scheme=="lineage" and learner=="gbt" and fold==2) or (scheme=="random" and learner=="lasso" and fold==4)
def random_fold(i): return (i//2)%5

def full_genomes(): return [f"g{i:06d}" for i in range(4000)]
def cohort_genomes(name):
 if name=="full": return full_genomes()
 if name=="prefix": return full_genomes()[:2000]
 if name=="added": return full_genomes()[2000:]
 return prep.matched_selection(full_genomes())
def lineage_of(g): return (int(g[1:])%2000%600)//6
def lineage_labels(name): return [lineage_of(g) for g in cohort_genomes(name)]
def genome_bytes(name): return ("\n".join(cohort_genomes(name))+"\n").encode()
def make_rules(drugs=None):
 drugs=drugs or {}
 sha=lambda name:hashlib.sha256(genome_bytes(name)).hexdigest()  # noqa: E731
 return {"prefix":vc.CohortRule(2000,sha("prefix"),drugs.get("prefix",vc.ALL_DRUGS),unavailable={}),
  "added":vc.CohortRule(2000,sha("added"),drugs.get("added",vc.ADDED_DRUGS),unavailable=FIXTURE_ADDED_UNAVAILABLE),
  "full":vc.CohortRule(4000,sha("full"),drugs.get("full",vc.ALL_DRUGS),unavailable={}),
  "matched2k":vc.CohortRule(2000,sha("matched2k"),None)}
LABELLED=full_genomes()[:60]+full_genomes()[2000:2060]
def label(drug,g):
 """A genome's reported call for a drug, the same in every cohort; None when it has none."""
 if g not in LABELLED: return None
 i=int(g[1:])
 if i<2000 or drug not in FIXTURE_ADDED_UNAVAILABLE: return i%2==0
 j=i-2000
 return {"nicotinamide":None,"levofloxacin":j!=0,"para-aminosalicylic acid":j<5}.get(drug,False)
def sample_ids(name,drug="amikacin"):
 order={"prefix":full_genomes()[:60],"added":full_genomes()[2000:2060],"full":LABELLED}.get(name)
 if order is None: order=[g for g in cohort_genomes(name) if g in set(LABELLED)]
 return [g for g in order if label(drug,g) is not None]
def drugs_of(name,rules):
 rule=rules[name]
 if rule.drugs is not None: return tuple(rule.drugs)
 return tuple(d for d in vc.ALL_DRUGS if min(class_counts(name,d))>=rule.min_class)
def class_counts(name,drug):
 labels=[label(drug,g) for g in sample_ids(name,drug)]
 return sum(labels),len(labels)-sum(labels)

def score(scheme,learner,i,n):
 li=vc.LEARNERS.index(learner); wrong=(i*13+li*5+(7 if scheme=="lineage" else 0))%n<li+1
 truth=i%2==0; call=(not truth) if wrong else truth
 if learner=="scm": return 1.0 if call else 0.0
 step=(i*7+li*3)%5
 return (27+step)/32 if call else (1+step)/32

def rescore(report):
 truth=[s["resistant"] for s in report["samples"]]; n=len(truth)
 for split in report["splits"]:
  if "unavailable_reason" in split: continue
  folds=split["test_fold"]
  for row in split["learners"]:
   name,scores=row["learner"],row["scores"]; ranked=name in vc.RANKED
   row["metrics"]=metrics_direct(truth,scores,ranked); row["roc"]=roc_direct(truth,scores) if ranked else []
   if name=="catalog": continue
   row["no_call"]={"resistant":sum(v is None and t for v,t in zip(scores,truth)),"susceptible":sum(v is None and not t for v,t in zip(scores,truth))}
   fba=[]
   for fold in range(split["folds"]):
    idx=[i for i,f in enumerate(folds) if f==fold]; m=metrics_direct([truth[i] for i in idx],[scores[i] for i in idx],False)
    fba.append(None if m is None else m["balanced_accuracy"])
   row["fold_balanced_accuracy"]=fba
 if report.get("catalog"):
  row=next(r for s in report["splits"] if s["scheme"]=="random" for r in s["learners"] if r["learner"]=="catalog")
  report["catalog"]["metrics"]=copy.deepcopy(row["metrics"])
  flags=[bool(v) for v in row["scores"]]
  report["catalog"]["groups"][0]["n_resistant"]=sum(f and t for f,t in zip(flags,truth))
  report["catalog"]["groups"][0]["n_susceptible"]=sum(f and not t for f,t in zip(flags,truth))
 return report

def catalog_call(g): return LABELLED.index(g)%3!=0 if g in LABELLED else False

def make_report(cohort,drug,lineage_by_id,lineage_available=True,fg=None,fingerprint="0"*16):
 fg=fg or foldgen()
 ids=sample_ids(cohort,drug); n=len(ids); did=vc.DRUG_IDS[drug]; truth=[label(drug,g) for g in ids]
 catalog_calls=[catalog_call(g) for g in ids]
 samples=[{"id":g,"resistant":truth[i],"from_mic":False,"mic":None,"lineage":lineage_by_id[g],"country":None,"year":None,
  "host":None,"mlst":None,"groups":[0] if catalog_calls[i] else [],"catalog_resistant":catalog_calls[i]} for i,g in enumerate(ids)]
 splits=[]
 for scheme in vc.SCHEMES:
  split={"scheme":scheme,"folds":5,"n_groups":len({s["lineage"] for s in samples}),"learners":[]}
  if scheme=="lineage" and not lineage_available:
   split["unavailable_reason"]="synthetic: insufficient independent lineage groups"; splits.append(split); continue
  folds=fg.stratified(truth,5,42) if scheme=="random" else fg.grouped(truth,[s["lineage"] for s in samples],5,42); split["test_fold"]=folds
  for name in vc.LEARNERS:
   scores=[score(scheme,name,i,n) for i in range(n)]; complexity=[{"scm":3.0,"cart":7.0,"rf":200.0,"gbt":40.0,"lasso":12.0}[name]]*5
   for fold in range(5):
    if null_fit(scheme,name,fold):
     complexity[fold]=0.0
     for i in range(n):
      if folds[i]==fold: scores[i]=None
   row={"learner":name,"metrics":None,"fold_balanced_accuracy":[],"roc":[],"scores":scores,"conformal":[],
        "fold_complexity":complexity,"seconds":1.0}
   if name=="lasso": row["lasso_fits"]=[fit(null_fit(scheme,name,fold)) for fold in range(5)]
   split["learners"].append(row)
  split["learners"].append({"learner":"catalog","metrics":None,"fold_balanced_accuracy":[],"roc":[],
   "scores":[1.0 if c else 0.0 for c in catalog_calls],"conformal":[],"fold_complexity":[],"seconds":0.0})
  splits.append(split)
 report={"schema":2,"evaluation":{"revision":8,"valid":True,"reasons":[],"scope":SCOPE},
  "input_fingerprint":fingerprint,"eval_config":copy.deepcopy(CONFIG),
  "id":did,"species":SPECIES_NAME,"species_slug":vc.SPECIES,"antibiotic":drug,"drug_class":"synthetic","breakpoint":None,
  "counts":{"samples":n,"resistant":sum(truth),"susceptible":n-sum(truth),"labels_from_mic":0,"labels_from_mic_excluded":0,
   "features":20,"patterns":20,"kmers":20,"k":31,"lineages":len({s["lineage"] for s in samples}),"min_support":max(n//200,3),"min_class":25},
  "lineage_threshold":0.0001,"catalog_groups":["fixtureGroup"],"samples":samples,"splits":splits,
  "catalog":{"version":CATALOG_VERSION,"rule":[["fixtureGroup"]],"metrics":None,"groups":[{"group":"fixtureGroup","class":"SYNTHETIC",
   "subclass":"FIXTURE","counts_for_rule":True,"n_resistant":0,"n_susceptible":0}]},"panels":[],"unsupported_assays":[],"mic":None,
  "models":{"scm":{},"cart":{},"cart_depth":0,"lasso":{"intercept":-1.0,"coef":[[4,0.5]],"lambda":0.01},"lasso_fit":fit(),
   "rf_importance":[],"gbt_importance":[]},"features":[],"seconds":10.0}
 return rescore(report)

def write_json(path,value):
 path.parent.mkdir(parents=True,exist_ok=True); path.write_text(json.dumps(value,separators=(",",":"),allow_nan=False)+"\n")
def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def done_line(report):
 best=vc.upstream_audit.full_coverage_best(report)
 if best is None: return f"    {report['id']} done in 3s (no learner calls every isolate)"
 return f"    {report['id']} done in 3s (best random-split learner among those that call every isolate: {best[0]} {best[1]:.3f})"
def lite(report,scheme,name):
 split=next(s for s in report["splits"] if s["scheme"]==scheme)
 if "unavailable_reason" in split: return None
 return vc.upstream_audit.lite(next(r for r in split["learners"] if r["learner"]==name),2)

def input_names(slug=vc.SPECIES):
 base=f"data/species/{slug}"
 return [f"{base}/genomes.txt",f"{base}/matrix/matrix.json",f"{base}/structure.json",f"{base}/counted.json",
  "data/curated/datasets.json","data/report/index.json","data/curated/phenotypes.jsonl","data/catalog/version.txt"]+\
  [f"{base}/matrix/{name}" for name in vc.MATRIX_BINS]

def refresh_root(root,name):
 slug=vc.SPECIES; results=root/"data/species"/slug/"results"; reports={}
 for path in sorted(results.glob("*.json")):
  report=json.loads(path.read_text()); reports[report["id"]]=report
  write_json(root/"data/report/datasets"/(report["id"]+".json"),report)
 best_rows=[]
 for did,report in sorted(reports.items()):
  best=vc.upstream_audit.full_coverage_best(report)
  lineage_best=None
  if best and next(s for s in report["splits"] if s["scheme"]=="lineage").get("unavailable_reason") is None:
   row=lite(report,"lineage",best[0]); lineage_best=None if row is None else row["balanced_accuracy"]
  best_rows.append({"id":did,"evaluation":report["evaluation"],"best_learner":None if best is None else best[0],
   "best_random":None if best is None else best[1],"best_lineage":lineage_best,
   "learners":[{"learner":n,"random":lite(report,"random",n),"lineage":lite(report,"lineage",n)} for n in vc.LEARNERS]})
 write_json(root/"data/report/index.json",{"schema":2,"evaluation_revision":8,"generated_at":"synthetic","datasets":best_rows,"panels":[],"unsupported_panels":[]})
 log=[done_line(report) for _,report in sorted(reports.items())]
 manifest_path=root/"validation-manifest.json"
 if manifest_path.exists():
  manifest=json.loads(manifest_path.read_text())
  for row in manifest["class_gate_unavailable"]: log.append(f"    {row['id']}: {vc.CLASS_SKIP}")
 else:
  for drug in sorted(set(vc.ALL_DRUGS)-{r["antibiotic"] for r in reports.values()}):
   if name=="added" and drug in FIXTURE_ADDED_UNAVAILABLE or name=="matched2k": log.append(f"    {vc.DRUG_IDS[drug]}: {vc.CLASS_SKIP}")
 log.append(f"{slug}: {len(reports)} dataset reports written"); (root/"pipeline.log").write_text("\n".join(log)+"\n")
 if manifest_path.exists():
  manifest=json.loads(manifest_path.read_text())
  for did,entry in manifest["reports"].items(): entry["sha256"]=sha(root/entry["path"])
  for entry in manifest["inputs"]: entry["sha256"]=sha(root/entry["path"])
  write_json(manifest_path,manifest)

def write_inputs(root,name,rules):
 """Everything fingerprint.py and the run manifest read, before any report exists."""
 slug=vc.SPECIES; sp=root/"data/species"/slug; genomes=cohort_genomes(name); labels=lineage_labels(name)
 sp.mkdir(parents=True,exist_ok=True)
 (sp/"genomes.txt").write_bytes(genome_bytes(name))
 digests=[hashlib.sha256(g.encode()).hexdigest()[:16] for g in genomes]
 n=len(genomes); words=(n+63)//64
 write_json(sp/"matrix/matrix.json",{"k":31,"genomes":genomes,"min_support":2,"n_patterns":20,"kmers_distinct":20,"kmers_kept":20,"kmers_rare":0,"kmers_core":0,"kset_digests":digests})
 (sp/"matrix/patterns.bin").write_bytes(b"".join(struct.pack("<Q",(p*0x9E3779B97F4A7C15+w)&(2**64-1)) for p in range(20) for w in range(words)))
 (sp/"matrix/support.bin").write_bytes(struct.pack("<20I",*range(2,22)))
 (sp/"matrix/kmer_offsets.bin").write_bytes(struct.pack("<21Q",*range(21)))
 (sp/"matrix/kmers.bin").write_bytes(struct.pack("<20Q",*range(1000,1020)))
 write_json(sp/"structure.json",{"genomes":genomes,"k":31,"kset_digests":digests,
  "lineages":{"threshold":0.0001,"labels":labels,"n_clusters":len(set(labels)),"largest_share":max(labels.count(x) for x in set(labels))/len(labels)},
  "merges":[],"leaf_order":[],"coords":[],"nearest":[],"sketch_size":1000})
 hit={"family":"fixtureGroup","group":"fixtureGroup","allele":"fixtureAllele","containment":1.0,"class":"SYNTHETIC","subclass":"FIXTURE"}
 counted={"k":31,"genomes":[{"genome_id":g,"n_kmers":0,"stats":{},"catalog_version":CATALOG_VERSION+fpmod.CALLS_SHAPE,"kset_digest":digests[i],
  "genes":[dict(hit)] if catalog_call(g) else []} for i,g in enumerate(genomes)],"missing_fasta":[],"failed":[]}
 write_json(sp/"counted.json",counted)
 curated=[{"id":did,"species":SPECIES_NAME,"species_slug":vc.SPECIES,"antibiotic":drug,"drug_class":"synthetic","n_resistant":0,"n_susceptible":0,
  "n_intermediate":0,"n_with_mic":0,"n_labels_from_mic":0,"breakpoint":None,"min_class":25} for drug,did in vc.DRUG_IDS.items()]
 write_json(root/"data/curated/datasets.json",curated)
 rule=rules[name]
 drugs=set(vc.ALL_DRUGS) if rule.drugs is None else set(rule.drugs)|set(rule.unavailable or {})
 lines=[]
 for drug in vc.ALL_DRUGS:
  if drug not in drugs: continue
  for g in genomes:
   call=label(drug,g)
   if call is None: continue
   lines.append(json.dumps({"genome_id":g,"species":SPECIES_NAME,"antibiotic":drug,"reported":"R" if call else "S","conflict":False,"mic":None,
    "mic_spread":0,"label":"R" if call else "S","label_from_mic":False,"methods":[],"standards":[],"n_records":1}))
 (root/"data/curated/phenotypes.jsonl").write_text("\n".join(lines)+"\n")
 (root/"data/catalog").mkdir(parents=True,exist_ok=True); (root/"data/catalog/version.txt").write_text(CATALOG_VERSION+"\n")
 return dict(zip(genomes,labels))

def build_cohort(root,name,rules,lineage_available=True,fg=None):
 root=Path(root); root.mkdir(parents=True,exist_ok=True)
 lineage_by_id=write_inputs(root,name,rules)
 inputs=fpmod.Inputs(root/"data",vc.SPECIES)
 drugs=drugs_of(name,rules)
 reports={drug:make_report(name,drug,lineage_by_id,lineage_available,fg,inputs.fingerprint(vc.DRUG_IDS[drug],CONFIG)) for drug in drugs}
 for drug,report in reports.items(): write_json(root/"data/species"/vc.SPECIES/"results"/(report["id"]+".json"),report)
 (root/"backup").mkdir(exist_ok=True); (root/"synthetic-binary").write_text("synthetic only\n"); (root/"synthetic-foldgen").write_text("#!/bin/sh\nexit 99\n"); (root/"synthetic-foldgen").chmod(0o755)
 refresh_root(root,name)
 report_entries={}
 for drug,report in reports.items():
  rel=f"data/species/{vc.SPECIES}/results/{report['id']}.json"
  report_entries[report["id"]]={"path":rel,"sha256":sha(root/rel),"input_fingerprint":report["input_fingerprint"]}
 unavailable=[]
 for drug in vc.ALL_DRUGS:
  if drug in drugs: continue
  rule=rules[name]
  if rule.unavailable is not None and drug not in rule.unavailable: continue
  r,s=class_counts(name,drug)
  unavailable.append({"id":vc.DRUG_IDS[drug],"antibiotic":drug,"min_class":25,"resistant":r,"susceptible":s,"reason":vc.CLASS_SKIP})
 names=input_names()
 manifest={"schema":"grm-mtb-4k-run-manifest-v1","cohort":name,"species_slug":vc.SPECIES,"source_commit":SOURCE,
  "binary_sha256":sha(root/"synthetic-binary"),"evaluation_revision":8,"report_schema":2,"eval_config":copy.deepcopy(CONFIG),
  "lineage_threshold":0.0001,"genomes":{"path":names[0],"sha256":sha(root/names[0])},
  "matrix":{"path":names[1],"sha256":sha(root/names[1]),"k":31,"min_support":2},
  "inputs":[{"path":rel,"sha256":sha(root/rel)} for rel in names],"reports":report_entries,"class_gate_unavailable":unavailable}
 write_json(root/"validation-manifest.json",manifest)
 refresh_root(root,name)
 return root

def passing_outputs(root,manifest,source=SOURCE):
 """The checker, audit and summary outputs of a passing run_revision8_checks.py run on this root,
 in the fields its postflight and verify_outputs read, for comparator tests that do not run the
 checks (the runner's own outputs are tested end to end). Since L2 the comparator re-reads them."""
 root=Path(root); n=len(manifest["reports"]); base=f"{root}/data/species/{vc.SPECIES}/results"
 skips=[{"file":f"{base}/{row['id']}.json"} for row in manifest["class_gate_unavailable"]]
 records=[{"id":did,"file":f"{base}/{did}.json","sha256":entry["sha256"],"issues":[],"revision":8,"schema":2,"valid":True,
  "written_after_launch":True,"identical_to_backup":False,"log_done":3.0,"published_revision":8,"published_identical":True,
  "manifest_compared":True,"counted_hits_match_samples":True,"samples_missing_from_counted":0,"counts":{"labels_from_mic":0}}
  for did,entry in sorted(manifest["reports"].items())]
 summary={"revision":8,"files":n,"by_revision":{"8":n},"by_schema":{"2":n},"errors":[],"files_with_issues":[],"not_in_current_curation":[],
  "not_run":[],"missing_curated_results":[],"issue_count":0,"manifest_schema":2,"manifest_evaluation_revision":8,"manifest_compared":n,
  "log_done":n,"species_finished":[vc.SPECIES],"class_minimum_skipped":skips,
  "revision8_wrapper":{"schema":"grm-revision8-audit-wrapper-v2","valid":True,"errors":[],"checker_pins":dict(vc.PINS),
   "run_source_commit":source,"run_manifest_sha256":sha(root/"validation-manifest.json"),"binary_sha256":manifest["binary_sha256"],
   "foldgen_sha256":foldgen_sha256(),"biological_execution":False}}
 checker={"reports":[{"id":did,"revision":8} for did in sorted(manifest["reports"])],"class_minimum_skipped":skips,
  "missing_curated_results":[],"not_in_current_curation":[]}
 return {"checker":checker,"audit":{"summary":summary,"records":records},"summary":{"coverage":{"8":{"files":n,"issues":0}}}}

def write_receipt(root,name=None,source=SOURCE,edit=None,outputs_edit=None):
 """A receipt in run_revision8_checks.py's form, bound to the root's current run manifest and to
 passing_outputs(); outputs_edit(outputs) may change those outputs before they are hashed, as a
 forged receipt over a failed run would bind them (L2)."""
 root=Path(root); manifest_path=root/"validation-manifest.json"; manifest=json.loads(manifest_path.read_text())
 bodies=passing_outputs(root,manifest,source)
 if outputs_edit: outputs_edit(bodies)
 outputs={}
 for key in ("checker","audit","summary"):
  p=root/f"synthetic-{key}.json"; p.write_text(json.dumps(bodies[key])+"\n"); outputs[key]={"path":str(p),"sha256":sha(p)}
 receipt={"schema":"grm-mtb-4k-revision8-check-receipt-v1","cohort":name or manifest["cohort"],"valid":True,"errors":[],"evaluation_revision":8,
  "reports":len(manifest["reports"]),"class_gate_skips":len(manifest["class_gate_unavailable"]),"run_manifest":str(manifest_path),
  "run_manifest_sha256":sha(manifest_path),"source_commit":source,"binary_sha256":manifest["binary_sha256"],"foldgen":foldgen_path(),
  "foldgen_sha256":sha(Path(foldgen_path())),"pins":dict(vc.PINS),"validator_sha256":vc.validator_sha256(),"outputs":outputs,
  "biological_execution":False}
 if edit: receipt=edit(receipt)
 path=root/"receipt.json"; path.write_text(receipt if isinstance(receipt,str) else json.dumps(receipt,indent=1)+"\n")
 return path

def contract_copy(directory,old=None,new=None):
 """report_contract.py, optionally with `old` replaced by `new` (exactly once), in a copy of the
 scripts tree it pins, so a source mutant loads against byte-identical checker blobs. The
 review harnesses' mutation tests copied K's vendored upstream/ the same way."""
 import shutil
 scripts=Path(directory)/"scripts"
 for rel in list(vc.PINS)+["mtb-4k/fingerprint.py","mtb-4k/revision_n.py"]:
  (scripts/rel).parent.mkdir(parents=True,exist_ok=True); shutil.copyfile(vc.SCRIPTS/rel,scripts/rel)
 source=(vc.LANE/"report_contract.py").read_text()
 if old is not None:
  if source.count(old)!=1: raise RuntimeError(("mutation needle count",source.count(old),old))
  source=source.replace(old,new)
 path=scripts/"mtb-4k"/"report_contract.py"; path.write_text(source)
 return path

def mutate_report(root,did,fn,rescore_after=False):
 path=Path(root)/"data/species"/vc.SPECIES/"results"/(did+".json"); report=json.loads(path.read_text()); fn(report)
 if rescore_after: rescore(report)
 write_json(path,report); refresh_root(Path(root),json.loads((Path(root)/"validation-manifest.json").read_text())["cohort"])
 return report

class FakeFoldGen:
 """In-process foldgen test double that emits K's index folds; for injected-failure tests only."""
 p=None
 def stratified(self,y,k=5,seed=42): return [random_fold(i) for i in range(len(y))]
 def grouped(self,y,groups,k=5,seed=42): return [g%k for g in groups]
