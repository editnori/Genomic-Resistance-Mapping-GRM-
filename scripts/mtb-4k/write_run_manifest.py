#!/usr/bin/env python3
"""Write a cohort's run manifest (grm-mtb-4k-run-manifest-v1) after `grm learn` and `grm report`.

usage: python3 -I write_run_manifest.py ROOT --cohort NAME --source-commit SHA --binary GRM [--log PIPELINE_LOG]

ROOT is the cohort's execution root (RUN/cohorts/<name>, holding data/). The manifest binds what
compare_cohorts.py and run_revision8_checks.py read, each by SHA-256: the genome list, the matrix
summary and its four binaries, structure.json, counted.json, curated datasets and phenotypes, the
catalog version and the dashboard index; the source commit and binary; the EvalConfig and lineage
threshold; and each report with its input_fingerprint. Before writing anything it:

- recomputes every report's input_fingerprint from the root's files (fingerprint.py) and stops on
  any difference (a stale or mismatched report);
- recomputes each curated dataset's class counts from the root's phenotypes on its matrix genomes,
  as `grm learn` gates them: a dataset that passes must have a report, one that does not gets a
  class-gate receipt with those counts and must have no report, and with --log its skip line must
  be in the pipeline log;
- checks the finished manifest with report_contract.validate_manifest against the frozen cohort
  rules, and writes it (validation-manifest.json, never over an existing one) only if it passes.

Exit 0 when written; 2 otherwise.
"""
import argparse,copy,json,sys
from pathlib import Path
import report_contract as vc
import fingerprint as fpmod

def build(root,cohort,source_commit,binary,log=None):
 data=root/"data"; slug=vc.SPECIES; base=f"data/species/{slug}"
 inputs=fpmod.Inputs(data,slug)
 reports={}; configs=[]
 for path in sorted((data/"species"/slug/"results").glob("*.json")):
  report,sha=vc.read_json(path)
  did=report["id"]
  if path.stem!=did: raise vc.ValidationError(f"{path}: report id {did} differs from its file name")
  want=inputs.fingerprint(did,report["eval_config"])
  if report["input_fingerprint"]!=want: raise vc.ValidationError(f"{path}: input_fingerprint {report['input_fingerprint']} differs from {want} recomputed from the root (stale report)")
  reports[did]={"path":f"{base}/results/{did}.json","sha256":sha,"input_fingerprint":want}
  configs.append(report["eval_config"])
 if not configs: raise vc.ValidationError(f"{root}: no reports")
 if any(c!=configs[0] for c in configs[1:]): raise vc.ValidationError(f"{root}: reports used different eval_config")
 skipped=None
 if log is not None:
  skipped={line.strip().split(":",1)[0] for line in Path(log).read_text(errors="replace").splitlines() if line.rstrip().endswith(vc.CLASS_SKIP)}
 unavailable=[]
 for did,spec in sorted(inputs.datasets.items()):
  r,s=inputs.class_counts(spec); minimum=spec.get("min_class",25)
  if r>=minimum and s>=minimum:
   if did not in reports: raise vc.ValidationError(f"{did}: {r} resistant and {s} susceptible pass the class gate, but there is no report")
   continue
  if did in reports: raise vc.ValidationError(f"{did}: a report exists although {r} resistant and {s} susceptible fail the class gate")
  if skipped is not None and did not in skipped: raise vc.ValidationError(f"{did}: no class-gate skip line in {log}")
  unavailable.append({"id":did,"antibiotic":spec["antibiotic"],"min_class":minimum,"resistant":r,"susceptible":s,"reason":vc.CLASS_SKIP})
 structure=json.loads((data/"species"/slug/"structure.json").read_text())
 names=[f"{base}/genomes.txt",f"{base}/matrix/matrix.json",f"{base}/structure.json",f"{base}/counted.json","data/curated/datasets.json",
  "data/report/index.json","data/curated/phenotypes.jsonl","data/catalog/version.txt"]+[f"{base}/matrix/{n}" for n in vc.MATRIX_BINS]
 summary=json.loads((data/"species"/slug/"matrix/matrix.json").read_text())
 return {"schema":"grm-mtb-4k-run-manifest-v1","cohort":cohort,"species_slug":slug,"source_commit":source_commit,
  "binary_sha256":vc.file_sha256(binary),"evaluation_revision":vc.REVISION,"report_schema":vc.SCHEMA,"eval_config":copy.deepcopy(configs[0]),
  "lineage_threshold":structure["lineages"]["threshold"],"genomes":{"path":names[0],"sha256":vc.file_sha256(root/names[0])},
  "matrix":{"path":names[1],"sha256":vc.file_sha256(root/names[1]),"k":summary["k"],"min_support":summary["min_support"]},
  "inputs":[{"path":rel,"sha256":vc.file_sha256(root/rel)} for rel in names],"reports":reports,"class_gate_unavailable":unavailable}

def main(argv=None,rules=None):
 ap=argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
 ap.add_argument("root",type=Path); ap.add_argument("--cohort",required=True)
 ap.add_argument("--source-commit",required=True); ap.add_argument("--binary",required=True,type=Path); ap.add_argument("--log",type=Path)
 args=ap.parse_args(argv)
 final=args.root/"validation-manifest.json"; part=args.root/"validation-manifest.json.part"
 try:
  if final.exists() or part.exists(): raise vc.ValidationError(f"{final} (or its .part) exists; a run manifest is written once")
  vc.digest(args.source_commit,"source commit",vc.COMMIT)
  manifest=build(args.root,args.cohort,args.source_commit,args.binary,args.log)
  with part.open("x") as f: f.write(json.dumps(manifest,indent=1,sort_keys=True,allow_nan=False)+"\n")
  try: cohort=vc.validate_manifest(part,rules or vc.RULES)
  except BaseException:
   part.unlink(); raise
  part.rename(final)
 except (vc.ValidationError,OSError,KeyError,TypeError,ValueError) as exc:
  print(f"write_run_manifest: STOP: {exc}",file=sys.stderr); return 2
 print(json.dumps({"manifest":str(final),"sha256":vc.file_sha256(final),"cohort":cohort.name,"reports":len(cohort.reports),
                   "class_gate_unavailable":len(manifest["class_gate_unavailable"])},sort_keys=True))
 return 0
if __name__=="__main__": sys.exit(main())
