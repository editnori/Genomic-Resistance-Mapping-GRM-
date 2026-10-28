#!/usr/bin/env python3
"""Fail-closed revision-8 gate for one cohort around the repo's own checker and metric audit.

usage: python3 -I run_revision8_checks.py AUDIT_OUT.json --root ROOT --run-manifest ROOT/validation-manifest.json
           --expected-source-commit SHA --binary GRM --backup DIR --log PIPELINE_LOG --launch ISO8601
           --foldgen FOLDGEN --checker-out CHECKER.json --summary-out SUMMARY.json --receipt RECEIPT.json

Port of package K's revision8_audit.py. It is not a parallel checker: it runs the in-tree
scripts/check-results.py, scripts/metric-audit/audit.py --revision 8 and summarize.py, after
report_contract.py has checked their SHA-256 against the reviewed pins.

1. Preflight: the run manifest and every report pass report_contract.validate_manifest strictly
   (any failure stops before the checker runs), the manifest's source commit is the expected one,
   and the binary's SHA-256 is the manifest's.
2. The checker, then the metric audit with the real foldgen, which rebuilds every saved test_fold
   of every report; a fold that differs is an audit issue.
3. Postflight (K's): revision 8 only, the expected report and class-gate skip counts, no issues,
   log, backup, launch and published-copy evidence, and summarizer coverage.
4. Only when all of that passes, RECEIPT.json: the receipt compare_cohorts.py requires and
   parses, binding the run manifest, source, binary, foldgen and checker pins to the outputs'
   SHA-256.

Exit 0 when the receipt is written; 2 otherwise (the audit output keeps the reasons).
"""
import argparse,datetime,json,os,subprocess,sys,time
from pathlib import Path
import report_contract as vc
CHECKER=vc.pinned("check-results.py")
AUDIT=vc.pinned("metric-audit/audit.py")
SUMMARIZE=vc.pinned("metric-audit/summarize.py")
RECEIPT_SCHEMA="grm-mtb-4k-revision8-check-receipt-v1"
def run_command(command,where):
 result=subprocess.run(command,capture_output=True,text=True,env=dict(os.environ,PYTHONDONTWRITEBYTECODE="1"))
 if result.returncode: raise vc.ValidationError(f"{where} exited {result.returncode}: {result.stderr[-4000:]}")
 return result

def postflight(audit_out,cohort,checker_out):
 run,_=vc.read_json(audit_out); summary=run.get("summary",{}); records=run.get("records",[]); expected=len(cohort.reports); errors=[]
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
 for rec in records:
  rid=rec.get("id",rec.get("file","unknown"))
  need("error" not in rec,f"{rid}: auditor exception")
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
 need({r.get("revision") for r in checker.get("reports",[])}=={8},"checker accepted mixed revisions")
 need({Path(x["file"]).stem for x in checker.get("class_minimum_skipped",[])}==expected_skips,"checker skip accounting differs")
 need(checker.get("missing_curated_results")==[] and checker.get("not_in_current_curation")==[],"checker scope incomplete")
 return errors

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
 parser.add_argument("--checker-out",required=True,type=Path); parser.add_argument("--summary-out",required=True,type=Path)
 parser.add_argument("--receipt",required=True,type=Path)
 args=parser.parse_args(argv); started=time.time(); outputs=[args.out,args.checker_out,args.summary_out,args.receipt]
 try:
  if len({str(p) for p in outputs})!=4 or any(p.exists() for p in outputs): raise vc.ValidationError("output paths must be distinct and absent")
  if args.pattern!=vc.SPECIES: raise vc.ValidationError("this lane audits only MTB")
  if args.run_manifest.parent.resolve()!=args.root.resolve(): raise vc.ValidationError("run manifest must be inside audited root")
  datetime.datetime.fromisoformat(args.launch.replace("Z","+00:00")).astimezone(datetime.timezone.utc)
  cohort=vc.validate_manifest(args.run_manifest,rules or vc.RULES)
  if cohort.manifest["source_commit"]!=args.expected_source_commit: raise vc.ValidationError("source commit differs from expected source")
  binary_hash=vc.file_sha256(args.binary)
  if binary_hash!=cohort.manifest["binary_sha256"]: raise vc.ValidationError("binary SHA-256 differs from manifest")
  if not args.foldgen.is_file() or not os.access(args.foldgen,os.X_OK): raise vc.ValidationError("foldgen must be executable")
  foldgen_hash=vc.file_sha256(args.foldgen)
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
  run["summary"]["revision8_wrapper"]={"schema":"grm-revision8-audit-wrapper-v2","valid":not errors,"errors":errors,
   "checker_pins":vc.PINS,"validator_sha256":vc.validator_sha256(),"run_source_commit":args.expected_source_commit,
   "run_manifest":str(args.run_manifest),"run_manifest_sha256":cohort.manifest_sha256,"binary":str(args.binary),"binary_sha256":binary_hash,
   "foldgen":str(args.foldgen),"foldgen_sha256":foldgen_hash,"elapsed_seconds":time.time()-started,
   "biological_execution":False,"independent_review":"REQUIRED"}
  args.out.write_text(json.dumps(run,indent=1,sort_keys=True,allow_nan=False)+"\n")
  if errors: raise vc.ValidationError("postflight rejected audit: "+"; ".join(errors))
  run_command([sys.executable,"-I",str(SUMMARIZE),str(args.out),str(args.summary_out)],"revision-8 summarizer")
  summary,_=vc.read_json(args.summary_out); coverage=summary.get("coverage",{}).get("8",{})
  if coverage.get("files")!=len(cohort.reports) or coverage.get("issues")!=0: raise vc.ValidationError("summarizer coverage/issues differ")
  receipt={"schema":RECEIPT_SCHEMA,"cohort":cohort.name,"valid":True,"errors":[],"evaluation_revision":vc.REVISION,
   "reports":len(cohort.reports),"class_gate_skips":len(cohort.manifest["class_gate_unavailable"]),
   "run_manifest":str(args.run_manifest),"run_manifest_sha256":cohort.manifest_sha256,"source_commit":args.expected_source_commit,
   "binary_sha256":binary_hash,"foldgen":str(args.foldgen),"foldgen_sha256":foldgen_hash,"pins":vc.PINS,"validator_sha256":vc.validator_sha256(),
   "outputs":{key:{"path":str(p),"sha256":vc.file_sha256(p)} for key,p in (("checker",args.checker_out),("audit",args.out),("summary",args.summary_out))},
   "biological_execution":False}
  with args.receipt.open("x") as f: f.write(json.dumps(receipt,indent=1,sort_keys=True,allow_nan=False)+"\n")
  print(json.dumps({"valid":True,"reports":len(cohort.reports),"audit":str(args.out),"checker":str(args.checker_out),
                    "summary":str(args.summary_out),"receipt":str(args.receipt)},sort_keys=True))
  return 0
 except (vc.ValidationError,OSError,KeyError,TypeError,ValueError,ArithmeticError) as exc:
  print(f"VALIDATION FAILED: {exc}",file=sys.stderr); return 2
if __name__=="__main__": sys.exit(main())
