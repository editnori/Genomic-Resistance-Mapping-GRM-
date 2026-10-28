#!/usr/bin/env python3
"""Fail-closed revision-8 gate for one cohort around the repo's own checker and metric audit.

usage: python3 -I run_revision8_checks.py AUDIT_OUT.json --root ROOT --run-manifest ROOT/validation-manifest.json
           --expected-source-commit SHA --binary GRM --backup DIR --log PIPELINE_LOG --launch ISO8601
           --foldgen FOLDGEN --expected-foldgen-sha256 SHA256
           --checker-out CHECKER.json --summary-out SUMMARY.json --receipt RECEIPT.json

Port of package K's revision8_audit.py. It is not a parallel checker: it runs the in-tree
scripts/check-results.py, scripts/metric-audit/audit.py --revision 8 and summarize.py, after
report_contract.py has checked their SHA-256 against the reviewed pins.

1. Preflight: the run manifest and every report pass report_contract.validate_manifest strictly
   (any failure stops before the checker runs), the manifest's source commit is the expected one,
   the binary's SHA-256 is the manifest's, and the foldgen's SHA-256 is the one the plan recorded
   when it built foldgen (--expected-foldgen-sha256; I4).
2. The checker, then the metric audit with the real foldgen, which rebuilds every saved test_fold
   of every report; a fold that differs is an audit issue.
3. Postflight (K's): revision 8 only, the expected report and class-gate skip counts, every
   audited record bound to the run manifest's report by id and SHA-256, no issues, log, backup,
   launch and published-copy evidence, and summarizer coverage.
4. Only when all of that passes, RECEIPT.json: the receipt compare_cohorts.py requires and
   parses, binding the run manifest, source, binary, foldgen and checker pins to the outputs'
   SHA-256. verify_outputs() is the same pass check over the outputs a receipt names, run here
   before the receipt is written and again by compare_cohorts.py (L2), which re-reads them.

Exit 0 when the receipt is written; 2 otherwise (the audit output keeps the reasons), including
when the validators cannot load (a moved pin).
"""
import argparse,datetime,importlib,json,os,subprocess,sys,time
from pathlib import Path
# M1: python3 -I (3.11+ implies -P) leaves this script's directory off sys.path. Put back that one
# directory, first, and nothing else, and require report_contract to load from it.
LANE=Path(__file__).resolve().parent
if sys.path[:1]!=[str(LANE)]: sys.path.insert(0,str(LANE))
def _sibling(name):
 module=importlib.import_module(name)
 if Path(module.__file__).resolve().parent!=LANE: raise ImportError(f"{name} loaded from {module.__file__}, not from {LANE}")
 return module
try:
 vc=_sibling("report_contract")
 CHECKER=vc.pinned("check-results.py")
 AUDIT=vc.pinned("metric-audit/audit.py")
 SUMMARIZE=vc.pinned("metric-audit/summarize.py")
except (ImportError,ValueError,OSError) as _exc:  # ValueError, OSError: report_contract's pin check
 print(f"run_revision8_checks: STOP: cannot load the validators: {_exc}",file=sys.stderr); sys.exit(2)
RECEIPT_SCHEMA="grm-mtb-4k-revision8-check-receipt-v1"
WRAPPER_SCHEMA="grm-revision8-audit-wrapper-v2"
def run_command(command,where):
 result=subprocess.run(command,capture_output=True,text=True,env=dict(os.environ,PYTHONDONTWRITEBYTECODE="1"))
 if result.returncode: raise vc.ValidationError(f"{where} exited {result.returncode}: {result.stderr[-4000:]}")
 return result

def postflight(audit_out,cohort,checker_out):
 # The run manifest's report set, not cohort.reports: a cohort loaded by another caller is
 # strict too, but the count must never follow what happened to load.
 run,_=vc.read_json(audit_out); summary=run.get("summary",{}); records=run.get("records",[]); expected=len(cohort.manifest["reports"]); errors=[]
 def need(condition,message):
  if not condition: errors.append(message)
 need(summary.get("revision")==8,"auditor revision is not 8")
 need(summary.get("files")==expected,"audited report count differs")
 need(summary.get("by_revision")=={"8":expected},"historical/wrong revisions in audited set")
 need(summary.get("by_schema")=={"2":expected},"wrong schemas in audited set")
 for key in ("errors","files_with_issues","not_in_current_curation","not_run","missing_curated_results"): need(summary.get(key)==[],f"audit field {key} is nonempty")
 need(summary.get("issue_count")==0,"auditor found metric/fold/catalog/probe issues")
 need(summary.get("manifest_schema")==2 and summary.get("manifest_evaluation_revision")==8,"manifest schema/revision differs")
 need(summary.get("manifest_compared")==expected,"manifest/published comparison incomplete")
 need(summary.get("log_done")==expected,"pipeline done-line coverage differs")
 need(summary.get("species_finished")==[vc.SPECIES],"species completion log missing")
 expected_skips={row["id"] for row in cohort.manifest["class_gate_unavailable"]}
 actual_skips={Path(x["file"]).stem for x in summary.get("class_minimum_skipped",[])}
 need(actual_skips==expected_skips,"class-gate skips differ")
 need(len(records)==expected,"record count differs")
 # Each audited record is one of the run manifest's reports, byte for byte.
 bound={did:entry["sha256"] for did,entry in cohort.manifest["reports"].items()}
 need(sorted(rec.get("id") for rec in records if isinstance(rec,dict))==sorted(bound),"audited report ids differ from the run manifest's")
 for rec in records:
  rid=rec.get("id",rec.get("file","unknown"))
  need("error" not in rec,f"{rid}: auditor exception")
  need(rec.get("sha256")==bound.get(rid),f"{rid}: audited report SHA-256 differs from the run manifest")
  need(rec.get("issues")==[],f"{rid}: report issues present")
  need(rec.get("revision")==8 and rec.get("schema")==2 and rec.get("valid") is True,f"{rid}: wrong revision/schema/validity")
  need(rec.get("written_after_launch") is True,f"{rid}: stale launch timestamp")
  need(rec.get("identical_to_backup") is False,f"{rid}: stale identical backup")
  need(rec.get("log_done") is not None,f"{rid}: no done line")
  need(rec.get("published_revision")==8 and rec.get("published_identical") is True,f"{rid}: published copy differs")
  need(rec.get("manifest_compared") is True,f"{rid}: manifest not compared")
  need(rec.get("counted_hits_match_samples") is True and rec.get("samples_missing_from_counted")==0,f"{rid}: counted calls differ")
  need(rec.get("counts",{}).get("labels_from_mic")==0,f"{rid}: MIC-inferred labels present")
 checker,_=vc.read_json(checker_out)
 need(len(checker.get("reports",[]))==expected,"checker report count differs")
 need(sorted(r.get("id") for r in checker.get("reports",[]))==sorted(bound),"checker report ids differ from the run manifest's")
 need({r.get("revision") for r in checker.get("reports",[])}=={8},"checker accepted mixed revisions")
 need({Path(x["file"]).stem for x in checker.get("class_minimum_skipped",[])}==expected_skips,"checker skip accounting differs")
 need(checker.get("missing_curated_results")==[] and checker.get("not_in_current_curation")==[],"checker scope incomplete")
 return errors

def summary_errors(summary_out,cohort):
 summary,_=vc.read_json(summary_out); coverage=summary.get("coverage",{}).get("8",{})
 if coverage.get("files")!=len(cohort.manifest["reports"]) or coverage.get("issues")!=0: return ["summarizer coverage/issues differ"]
 return []

def verify_outputs(outputs,cohort,source_commit,binary_sha256,foldgen_sha256):
 """L2: re-read the checker, audit and summary outputs a receipt names and require that they
 record a pass for this cohort's run manifest: the postflight again, the wrapper block the runner
 wrote (valid, no errors, the same run manifest, source, binary, foldgen and pins), and the
 summarizer's coverage. Returns the reasons it does not; [] when it does."""
 errors=postflight(outputs["audit"],cohort,outputs["checker"])
 run,_=vc.read_json(outputs["audit"]); wrapper=run.get("summary",{}).get("revision8_wrapper")
 if not isinstance(wrapper,dict): return errors+["the audit output has no revision8_wrapper block"]
 expected={"schema":WRAPPER_SCHEMA,"valid":True,"errors":[],"checker_pins":vc.PINS,"run_source_commit":source_commit,
  "run_manifest_sha256":cohort.manifest_sha256,"binary_sha256":binary_sha256,"foldgen_sha256":foldgen_sha256,"biological_execution":False}
 for key,want in expected.items():
  if wrapper.get(key)!=want: errors.append(f"audit wrapper {key} is {wrapper.get(key)!r}, not {want!r}")
 return errors+summary_errors(outputs["summary"],cohort)

def close_foldgens(instances):
 for fg in instances:
  process=getattr(fg,"p",None)
  if process is None: continue
  try:
   if process.stdin: process.stdin.close()
  except OSError: pass
  try: process.wait(timeout=5)
  except subprocess.TimeoutExpired:
   process.terminate()
   try: process.wait(timeout=5)
   except subprocess.TimeoutExpired: process.kill(); process.wait()
  if process.poll() is None: raise vc.ValidationError("foldgen remained running")

def main(argv=None,rules=None,foldgen_factory=None):
 parser=argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
 parser.add_argument("out",type=Path); parser.add_argument("pattern",nargs="?",default=vc.SPECIES)
 parser.add_argument("--root",required=True,type=Path); parser.add_argument("--run-manifest",required=True,type=Path)
 parser.add_argument("--expected-source-commit",required=True); parser.add_argument("--binary",required=True,type=Path)
 parser.add_argument("--backup",required=True,type=Path); parser.add_argument("--log",required=True,type=Path)
 parser.add_argument("--launch",required=True); parser.add_argument("--foldgen",required=True,type=Path)
 parser.add_argument("--expected-foldgen-sha256",required=True)
 parser.add_argument("--checker-out",required=True,type=Path); parser.add_argument("--summary-out",required=True,type=Path)
 parser.add_argument("--receipt",required=True,type=Path)
 args=parser.parse_args(argv); started=time.time(); outputs=[args.out,args.checker_out,args.summary_out,args.receipt]
 try:
  if len({str(p) for p in outputs})!=4 or any(p.exists() for p in outputs): raise vc.ValidationError("output paths must be distinct and absent")
  if args.pattern!=vc.SPECIES: raise vc.ValidationError("this lane audits only MTB")
  if args.run_manifest.parent.resolve()!=args.root.resolve(): raise vc.ValidationError("run manifest must be inside audited root")
  datetime.datetime.fromisoformat(args.launch.replace("Z","+00:00")).astimezone(datetime.timezone.utc)
  vc.digest(args.expected_foldgen_sha256,"--expected-foldgen-sha256")
  cohort=vc.validate_manifest(args.run_manifest,rules or vc.RULES)
  if cohort.manifest["source_commit"]!=args.expected_source_commit: raise vc.ValidationError("source commit differs from expected source")
  binary_hash=vc.file_sha256(args.binary)
  if binary_hash!=cohort.manifest["binary_sha256"]: raise vc.ValidationError("binary SHA-256 differs from manifest")
  if not args.foldgen.is_file() or not os.access(args.foldgen,os.X_OK): raise vc.ValidationError("foldgen must be executable")
  foldgen_hash=vc.file_sha256(args.foldgen)
  if foldgen_hash!=args.expected_foldgen_sha256: raise vc.ValidationError(f"foldgen SHA-256 {foldgen_hash} differs from the recorded build {args.expected_foldgen_sha256}")
  run_command([sys.executable,"-I",str(CHECKER),str(args.root/"data"),"--log",str(args.log),"--skip-log",str(args.log),"--out",str(args.checker_out)],"revision-8 checker")
  audit=vc.upstream_audit; old_argv=sys.argv[:]; old_factory=audit.FoldGen; instances=[]
  factory=foldgen_factory or old_factory
  def tracked_factory():
   fg=factory(); instances.append(fg); return fg
  audit.FoldGen=tracked_factory; audit.species_finished.clear(); audit_code=0
  sys.argv=[str(AUDIT),str(args.out),args.pattern,"--revision","8","--root",str(args.root),"--backup",str(args.backup),
            "--log",str(args.log),"--launch",args.launch,"--foldgen",str(args.foldgen)]
  try:
   try: audit.main()
   except SystemExit as exc: audit_code=int(exc.code or 0)
  finally:
   sys.argv=old_argv; audit.FoldGen=old_factory; close_foldgens(instances)
  if audit_code or not args.out.exists(): raise vc.ValidationError(f"in-tree auditor exited {audit_code} without complete receipt")
  errors=postflight(args.out,cohort,args.checker_out)
  run,_=vc.read_json(args.out)
  run["summary"]["revision8_wrapper"]={"schema":WRAPPER_SCHEMA,"valid":not errors,"errors":errors,
   "checker_pins":vc.PINS,"validator_sha256":vc.validator_sha256(),"run_source_commit":args.expected_source_commit,
   "run_manifest":str(args.run_manifest),"run_manifest_sha256":cohort.manifest_sha256,"binary":str(args.binary),"binary_sha256":binary_hash,
   "foldgen":str(args.foldgen),"foldgen_sha256":foldgen_hash,"elapsed_seconds":time.time()-started,
   "biological_execution":False,"independent_review":"REQUIRED"}
  args.out.write_text(json.dumps(run,indent=1,sort_keys=True,allow_nan=False)+"\n")
  if errors: raise vc.ValidationError("postflight rejected audit: "+"; ".join(errors))
  run_command([sys.executable,"-I",str(SUMMARIZE),str(args.out),str(args.summary_out)],"revision-8 summarizer")
  bound={"checker":args.checker_out,"audit":args.out,"summary":args.summary_out}
  errors=verify_outputs(bound,cohort,args.expected_source_commit,binary_hash,foldgen_hash)
  if errors: raise vc.ValidationError("the outputs do not record a pass: "+"; ".join(errors))
  receipt={"schema":RECEIPT_SCHEMA,"cohort":cohort.name,"valid":True,"errors":[],"evaluation_revision":vc.REVISION,
   "reports":len(cohort.reports),"class_gate_skips":len(cohort.manifest["class_gate_unavailable"]),
   "run_manifest":str(args.run_manifest),"run_manifest_sha256":cohort.manifest_sha256,"source_commit":args.expected_source_commit,
   "binary_sha256":binary_hash,"foldgen":str(args.foldgen),"foldgen_sha256":foldgen_hash,"pins":vc.PINS,"validator_sha256":vc.validator_sha256(),
   "outputs":{key:{"path":str(p),"sha256":vc.file_sha256(p)} for key,p in bound.items()},
   "biological_execution":False}
  with args.receipt.open("x") as f: f.write(json.dumps(receipt,indent=1,sort_keys=True,allow_nan=False)+"\n")
  print(json.dumps({"valid":True,"reports":len(cohort.reports),"audit":str(args.out),"checker":str(args.checker_out),
                    "summary":str(args.summary_out),"receipt":str(args.receipt)},sort_keys=True))
  return 0
 except (vc.ValidationError,OSError,KeyError,TypeError,ValueError,ArithmeticError) as exc:
  print(f"VALIDATION FAILED: {exc}",file=sys.stderr); return 2
if __name__=="__main__": sys.exit(main())
