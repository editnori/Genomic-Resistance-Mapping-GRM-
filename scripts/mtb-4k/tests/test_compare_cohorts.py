#!/usr/bin/env python3
"""compare_cohorts.py on the canonical strict fixtures (real foldgen folds).

usage: FOLDGEN=<foldgen binary> python3 -I scripts/mtb-4k/tests/test_compare_cohorts.py

- ComparisonTests: package K's ComparisonTests, path-adapted to four cohorts, roots and receipts.
  Changed assertions, each named in its test: the D1 decision (a random-fold mismatch is recorded
  and invalidates that row's matched random control, not the comparison), counts that now follow
  the real foldgen folds instead of K's index folds, and a missing genome that the merged contract
  now stops as an invalid report (K's coverage check is kept at function level).
- WScenarioTests: package W's six comparator scenarios on these fixtures.
- FixTests: controls for F2-F7 and the re-review's N6.
- FollowupTests: the plan-validators review's L2 (a receipt whose bound outputs record a failure
  is refused, however it is forged), L3 (any invalid report stops the comparison), I4 (every
  receipt's foldgen is the recorded build) and I9 (a matched control with no pair is invalid).
- EndToEndRealFoldsTests: experiment E2 as a passing run: run_revision8_checks.py with the real
  foldgen on all four cohorts, then the comparator on their receipts.

Since L2 the comparator re-reads the outputs a receipt binds, so fx.write_receipt binds outputs in
the runner's form that record a pass (fx.passing_outputs); the runner's real outputs pass the same
check in EndToEndRealFoldsTests.
"""
import contextlib,copy,hashlib,io,json,shutil,sys,tempfile,unittest
from pathlib import Path
LANE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(LANE)); sys.path.insert(0,str(LANE/"tests"))
import report_contract as vc
import compare_cohorts as cc
import run_revision8_checks as runner
import fixtures as fx

NEGATIVE=[]
def reject(case,label,match,fn):
 try: fn()
 except vc.ValidationError as exc:
  case.assertIn(match,str(exc)); NEGATIVE.append({"control":label,"matched":match,"result":"rejected"}); return
 case.fail("accepted negative control: "+label)

def screen(base,exact=0):
 path=Path(base)/f"screen-{exact}.json"; path.write_text(json.dumps({"cross_cohort_exact":exact,"cross_cohort_near":0})+"\n"); return path

def run_compare(roots,rules,source=fx.SOURCE,receipt_source=None,duplicates=None,accept=None,edit=None,outputs_edit=None,foldgen=None):
 receipts={n:str(fx.write_receipt(r,source=receipt_source or source,edit=(edit.get(n) if edit else None),
                                  outputs_edit=(outputs_edit.get(n) if outputs_edit else None))) for n,r in roots.items()}
 return cc.compare({n:str(r) for n,r in roots.items()},receipts,str(duplicates or screen(Path(roots["prefix"]).parent)),source,rules,accept,
                   expected_foldgen=foldgen or fx.foldgen_sha256())

def run_main(roots,base,rules,extra=(),outputs_edit=None):
 """compare_cohorts.main on these roots: (exit code, the output or None, stderr)."""
 receipts=[f"--checker-receipt={n}={fx.write_receipt(r,outputs_edit=(outputs_edit.get(n) if outputs_edit else None))}" for n,r in roots.items()]
 out=Path(base)/f"out-{len(list(Path(base).glob('out-*')))}.json"
 argv=[f"--cohort={n}={r}" for n,r in roots.items()]+receipts+["--duplicates",str(screen(base)),"--expected-source-commit",fx.SOURCE,
       "--expected-foldgen-sha256",fx.foldgen_sha256(),"--out",str(out),*extra]
 saved=cc.vc.RULES
 cc.vc.RULES=rules
 err=io.StringIO()
 try:
  with contextlib.redirect_stdout(io.StringIO()),contextlib.redirect_stderr(err): code=cc.main(argv)
 finally: cc.vc.RULES=saved
 return code,(json.loads(out.read_text()) if out.exists() else None),err.getvalue()

def report_path(root,did): return Path(root)/"data/species"/vc.SPECIES/"results"/(did+".json")
def rehash(root,did):
 path=Path(root)/"validation-manifest.json"; manifest=json.loads(path.read_text()); entry=manifest["reports"][did]; entry["sha256"]=fx.sha(Path(root)/entry["path"]); fx.write_json(path,manifest)

def relabel(root,drug,genome,resistant):
 """A cohort self-consistent with another label for one genome: its phenotype, report sample,
 metrics, fingerprint and manifest all change together."""
 root=Path(root); pheno=root/"data/curated/phenotypes.jsonl"; lines=[]
 for line in pheno.read_text().splitlines():
  p=json.loads(line)
  if p["genome_id"]==genome and p["antibiotic"]==drug: p["reported"]=p["label"]="R" if resistant else "S"
  lines.append(json.dumps(p))
 pheno.write_text("\n".join(lines)+"\n")
 did=vc.DRUG_IDS[drug]
 def edit(r):
  for s in r["samples"]:
   if s["id"]==genome: s["resistant"]=resistant
  r["counts"]["resistant"]=sum(s["resistant"] for s in r["samples"]); r["counts"]["susceptible"]=len(r["samples"])-r["counts"]["resistant"]
 fx.mutate_report(root,did,edit,rescore_after=True)
 refingerprint(root)

def refingerprint(root):
 root=Path(root); manifest=json.loads((root/"validation-manifest.json").read_text()); inputs=fx.fpmod.Inputs(root/"data",vc.SPECIES)
 for did,entry in manifest["reports"].items():
  path=root/entry["path"]; r=json.loads(path.read_text()); r["input_fingerprint"]=entry["input_fingerprint"]=inputs.fingerprint(did,r["eval_config"]); fx.write_json(path,r)
 for row in manifest["class_gate_unavailable"]:
  row["resistant"],row["susceptible"]=inputs.class_counts(inputs.datasets[row["id"]])
 fx.write_json(root/"validation-manifest.json",manifest); fx.refresh_root(root,manifest["cohort"])

class Base(unittest.TestCase):
 @classmethod
 def setUpClass(cls):
  cls.tmp=tempfile.TemporaryDirectory(); cls.rules=fx.make_rules(); cls.roots={}
  for name in fx.COHORTS: cls.roots[name]=fx.build_cohort(Path(cls.tmp.name)/name,name,cls.rules)
  cls.baseline=run_compare(cls.roots,cls.rules)
 @classmethod
 def tearDownClass(cls): cls.tmp.cleanup()
 def fresh(self,label,name):
  root=Path(self.tmp.name)/label; shutil.copytree(self.roots[name],root); return root
 def with_root(self,name,root): return {**self.roots,name:root}
 def row(self,out,analysis,drug,scheme,learner=0): return next(r for r in out["analyses"][analysis]["rows"] if r["drug"]==drug and r["scheme"]==scheme)["comparisons"][learner]

class ComparisonTests(Base):
 def test_22_14_22_and_eight_unavailable_accounting(self):
  out=self.baseline; self.assertTrue(out["validation"]["valid"],out["validation"]["blocking_issues"])
  # Path-adapted: four cohorts; matched2k's 14 follow the fixture's class gate on its genomes.
  self.assertEqual({n:c["eligible_reports"] for n,c in out["cohorts"].items()},{"prefix":22,"added":14,"full":22,"matched2k":14})
  unavailable=out["cohorts"]["added"]["class_gate_unavailable"]; self.assertEqual(len(unavailable),8); self.assertEqual({r["antibiotic"] for r in unavailable},set(vc.ADDED_UNAVAILABLE))
  self.assertEqual(out["cohorts"]["prefix"]["class_gate_unavailable"],[]); self.assertEqual(out["cohorts"]["full"]["class_gate_unavailable"],[])
  self.assertEqual(len(out["analyses"]["full_vs_prefix"]["rows"]),44); self.assertEqual(len(out["analyses"]["added_vs_prefix"]["rows"]),28)
  for row in out["analyses"]["full_vs_prefix"]["rows"]:
   self.assertEqual(len(row["comparisons"]),5)
   for cohort in ("prefix","full"):
    if row[cohort]["status"]=="available": self.assertEqual(row[cohort]["five_learner_medians"]["balanced_accuracy"]["prespecified_learners"],5)
 def test_disjoint_and_full_added_semantics(self):
  out=self.baseline; pa=out["analyses"]["added_vs_prefix"]; self.assertEqual(pa["analysis_type"],"added_vs_prefix_side_by_side")
  self.assertFalse(pa["paired_isolates"]); self.assertFalse(pa["frozen_model_external_validation"])
  fa=out["analyses"]["full_cohort_cv_added_subset"]; self.assertFalse(fa["frozen_model_external_validation"])
  row=fa["rows"][0]; self.assertEqual(row["added_only_cv"]["learners"][0]["metrics"]["n"],60)
  self.assertEqual(row["full_cohort_cv_added_subset"]["learners"][0]["metrics"]["n"],60)
  self.assertIn("not a frozen",fa["note"]); self.assertNotIn("matched_fold",pa["rows"][0])
 def test_row_order_does_not_change_pairing(self):
  full=self.fresh("row-order-full","full"); did=vc.DRUG_IDS["amikacin"]; path=report_path(full,did); report=json.loads(path.read_text())
  report["samples"].reverse()
  for split in report["splits"]:
   split["test_fold"].reverse()
   for row in split["learners"]: row["scores"].reverse()
  fx.write_json(path,report); rehash(full,did); out=run_compare(self.with_root("full",full),self.rules)
  self.assertTrue(out["validation"]["valid"],out["validation"]["blocking_issues"])
  self.assertEqual(self.row(out,"full_vs_prefix","amikacin","random")["matched_fold"],self.row(self.baseline,"full_vs_prefix","amikacin","random")["matched_fold"])
 def test_fold_mismatch_excluded_and_blocks_matched_random(self):
  full=self.fresh("fold-mismatch-full","full"); did=vc.DRUG_IDS["amikacin"]; path=report_path(full,did); report=json.loads(path.read_text()); split=report["splits"][0]
  i=next(i for i,f in enumerate(split["test_fold"]) if f==0 and report["samples"][i]["resistant"])
  j=next(j for j,f in enumerate(split["test_fold"]) if f==1 and report["samples"][j]["resistant"])
  split["test_fold"][i],split["test_fold"][j]=split["test_fold"][j],split["test_fold"][i]; fx.rescore(report); fx.write_json(path,report); rehash(full,did)
  out=run_compare(self.with_root("full",full),self.rules)
  # CHANGED ASSERTION (decision D1): K asserted valid=False with "random folds differ" as a
  # blocking issue. The mismatch is now recorded per row and invalidates that row's matched
  # random control only; fold integrity is enforced by run_revision8_checks.py's foldgen rebuild
  # (test_run_revision8_checks.py: test_foldgen_rebuild_trap_fails_postflight, real foldgen).
  self.assertTrue(out["validation"]["valid"],out["validation"]["blocking_issues"])
  self.assertFalse(any("random folds differ" in x["issue"] for x in out["validation"]["blocking_issues"]))
  self.assertTrue(any("random folds differ" in w["issue"] for w in out["validation"]["warnings"]))
  row=self.row(out,"full_vs_prefix","amikacin","random")
  # Counts from the two reports' saved folds (K's 2 and 58 assumed identical index folds).
  prefix=json.loads(report_path(self.roots["prefix"],did).read_text()); pf=dict(zip([s["id"] for s in prefix["samples"]],prefix["splits"][0]["test_fold"]))
  ff=dict(zip([s["id"] for s in report["samples"]],split["test_fold"])); want=[g for g in pf if pf[g]!=ff[g]]
  self.assertEqual(row["pairing"]["fold_mismatch"],want); self.assertEqual(row["pairing"]["fold_mismatch_count"],len(want))
  self.assertEqual(row["pairing"]["eligible_matched_pairs"],60-len(want)); self.assertEqual(row["genome_aligned_all_folds"]["pairs"],60)
  self.assertFalse(row["matched_fold"]["control_valid"]); self.assertIn("random folds differ",row["matched_fold"]["control_invalid_reason"])
  NEGATIVE.append({"control":"held-out fold mismatch","result":"excluded from matched_fold; matched random control invalid (D1)"})
 def test_missing_full_genome_is_coverage_failure(self):
  full=self.fresh("genome-mismatch-full","full"); did=vc.DRUG_IDS["amikacin"]; path=report_path(full,did); report=json.loads(path.read_text())
  report["samples"].pop(0); report["counts"]["samples"]=119; report["counts"]["resistant"]=59
  for split in report["splits"]:
   split["test_fold"].pop(0)
   for row in split["learners"]: row["scores"].pop(0)
  fx.rescore(report); fx.write_json(path,report); rehash(full,did)
  # CHANGED ASSERTION: the merged contract checks report counts against the root's own
  # phenotypes, so this report is invalid at load where K reached a coverage issue. Since L3 an
  # invalid report stops the comparison (it was excluded before). K's coverage assertions are kept
  # on pair_learner itself below.
  reject(self,"report counts differ from the root's phenotypes (L3 stop)","root's phenotypes",lambda:run_compare(self.with_root("full",full),self.rules))
  p=vc.validate_manifest(self.roots["prefix"]/"validation-manifest.json",self.rules)
  row,issues=cc.pair_learner(p.reports[did],self.coverage_report(did),"random","scm",set(p.genomes),{g:True for g in p.genomes})
  self.assertTrue(any("genome coverage" in x for x in issues))
  self.assertEqual(row["pairing"]["common_samples"],59); self.assertEqual(row["pairing"]["prefix_only"],["g000000"]); NEGATIVE.append({"control":"genome coverage mismatch","result":"blocked"})
 def coverage_report(self,did):
  """The full report with g000000 removed, as K's test had it reach pair_learner."""
  c=vc.validate_manifest(self.roots["full"]/"validation-manifest.json",self.rules); r=c.reports[did]
  r.samples.pop("g000000"); r.sample_order.remove("g000000")
  for split in r.splits.values():
   for learner in split.learners.values(): learner.scores.pop("g000000"); learner.folds.pop("g000000")
  return r
 def test_changed_lineage_groups_recorded_not_assumed_stable(self):
  full=self.fresh("lineage-change-full","full"); sp=full/"data/species"/vc.SPECIES/"structure.json"; structure=json.loads(sp.read_text()); labels=structure["lineages"]["labels"]; labels[0]=100
  structure["lineages"]["n_clusters"]=len(set(labels)); structure["lineages"]["largest_share"]=max(labels.count(x) for x in set(labels))/len(labels); fx.write_json(sp,structure)
  manifest=json.loads((full/"validation-manifest.json").read_text())
  for entry in manifest["inputs"]:
   if entry["path"].endswith("structure.json"): entry["sha256"]=fx.sha(sp)
  fx.write_json(full/"validation-manifest.json",manifest)
  for drug in vc.ALL_DRUGS:
   did=vc.DRUG_IDS[drug]; path=report_path(full,did); report=json.loads(path.read_text())
   sample=next((s for s in report["samples"] if s["id"]=="g000000"),None)
   if sample is None: continue
   sample["lineage"]=100
   report["counts"]["lineages"]=len({s["lineage"] for s in report["samples"]})
   for split in report["splits"]: split["n_groups"]=report["counts"]["lineages"]
   fx.write_json(path,report); rehash(full,did)
  refingerprint(full)  # path-adapted: the lineage labels are fingerprint inputs
  out=run_compare(self.with_root("full",full),self.rules); self.assertTrue(out["validation"]["valid"],out["validation"]["blocking_issues"])
  lineage=out["lineage_comparison"]["prefix_in_full"]; self.assertFalse(lineage["same_common_genome_partition"]); self.assertEqual(lineage["changed_raw_group_ids"],1)
  self.assertEqual(lineage["split_prefix_groups"],{"0":[0,100]}); self.assertTrue(out["validation"]["warnings"])
  row=self.row(out,"full_vs_prefix","amikacin","lineage")
  self.assertGreater(row["pairing"]["changed_lineage_group_count"],0); self.assertLess(row["pairing"]["eligible_matched_pairs"],60)
 def test_disjoint_overlap_negative_control(self):
  p=vc.validate_manifest(self.roots["prefix"]/"validation-manifest.json",self.rules); a=vc.validate_manifest(self.roots["added"]/"validation-manifest.json",self.rules); a.genomes[0]=p.genomes[0]
  # Path-adapted (L3): no exclusion argument any more.
  _,issues=cc.compare_prefix_added(p,a,set()); self.assertTrue(any("not disjoint" in x["issue"] for x in issues)); NEGATIVE.append({"control":"non-disjoint prefix/added","result":"blocked"})
 def test_missing_or_wrong_unavailable_receipt_rejected(self):
  original=json.loads((self.roots["added"]/"validation-manifest.json").read_text())
  for label,mutate in [("missing eight-drug receipt",lambda m:m["class_gate_unavailable"].pop()),("wrong gate counts",lambda m:m["class_gate_unavailable"][0].update(resistant=25))]:
   manifest=copy.deepcopy(original); mutate(manifest); path=self.roots["added"]/(label.replace(" ","-")+".json"); fx.write_json(path,manifest)
   reject(self,label,"receipt" if label.startswith("wrong") else "missing unavailable",lambda path=path:vc.validate_manifest(path,self.rules))
 def test_source_and_binary_mismatch_blocked(self):
  # Path-adapted: receipts now bind the source, so they are written for the expected source.
  out=run_compare(self.roots,self.rules,source="2"*40); self.assertFalse(out["validation"]["valid"]); self.assertTrue(any("source commit" in x["issue"] for x in out["validation"]["blocking_issues"]))
  full=self.fresh("bad-binary-full","full"); manifest=json.loads((full/"validation-manifest.json").read_text()); manifest["binary_sha256"]="0"*64; fx.write_json(full/"validation-manifest.json",manifest)
  out=run_compare(self.with_root("full",full),self.rules); self.assertFalse(out["validation"]["valid"]); self.assertTrue(any("binaries differ" in x["issue"] for x in out["validation"]["blocking_issues"]))
  NEGATIVE.append({"control":"source/binary mismatch","result":"blocked"})

class WScenarioTests(Base):
 """Package W's six comparator scenarios (its tests/test_compare_cohorts.py), on strict fixtures."""
 def test_statistics(self):
  self.assertAlmostEqual(vc.complement(vc.clopper_pearson(10,10))["upper"],0.30850,places=4)
  ci=vc.complement(vc.clopper_pearson(5,10)); self.assertAlmostEqual(ci["lower"],0.18709,places=4); self.assertAlmostEqual(ci["upper"],0.81291,places=4)
  self.assertAlmostEqual(vc.complement(vc.clopper_pearson(0,10))["lower"],0.69150,places=4)
  self.assertIsNone(vc.complement(vc.clopper_pearson(0,0))["lower"])
  self.assertAlmostEqual(cc.mcnemar_exact(0,5),0.0625); self.assertEqual(cc.mcnemar_exact(0,0),1.0)
  self.assertEqual(cc.holm({"a":0.01,"b":0.04,"c":0.03}),{"a":0.03,"b":0.06,"c":0.06})
 def test_invalid_in_one_cohort_stops_the_comparison(self):
  matched=self.fresh("w-invalid-matched","matched2k"); did=vc.DRUG_IDS["amikacin"]
  fx.mutate_report(matched,did,lambda r:r["evaluation"].update(valid=False,reasons=["final model has an unconverged two-class lasso fit"]))
  code,out,err=run_main(self.with_root("matched2k",matched),self.tmp.name,self.rules)
  # CHANGED ASSERTION (L3): W's rule excluded the drug in matched2k and exited 3. Nothing upstream
  # can produce that state (grm learn and the runner both stop on an invalid report), so the
  # comparator now stops too: exit 2, no output, the report named.
  self.assertEqual(code,2); self.assertIsNone(out); self.assertIn("evaluation.valid=true",err); self.assertIn(did,err)
  NEGATIVE.append({"control":"invalid in one cohort","result":"stop, exit 2, no output (L3)"})
 def test_not_interpretable_drugs_are_derived_flagged_and_kept_out_of_holm(self):
  out=self.baseline; flagged=set(out["not_interpretable_for_before_after"]); self.assertEqual(flagged,set(vc.ADDED_UNAVAILABLE))
  for analysis in ("full_vs_prefix","full_vs_matched2k"):
   for row in out["analyses"][analysis]["rows"]:
    self.assertEqual(row["before_after_interpretable"],row["drug"] not in flagged)
    for comp in row["comparisons"]:
     if comp["genome_aligned_all_folds"]: self.assertEqual("mcnemar_p_holm" in comp["genome_aligned_all_folds"],row["drug"] not in flagged)
  self.assertEqual(out["added_class_counts_from_full"]["para-aminosalicylic acid"],[5,55])
 def test_a_report_whose_metrics_do_not_match_its_scores_is_invalid(self):
  added=self.fresh("w-metrics-added","added"); did=vc.DRUG_IDS["amikacin"]
  fx.mutate_report(added,did,lambda r:r["splits"][0]["learners"][2]["metrics"].update(tp=r["splits"][0]["learners"][2]["metrics"]["tp"]+1))
  code,out,err=run_main(self.with_root("added",added),self.tmp.name,self.rules)
  # CHANGED ASSERTION (L3): a stop with no output, where W excluded the drug and exited 3.
  self.assertEqual(code,2); self.assertIsNone(out); self.assertIn("recomputed",err)
 def test_exact_cross_cohort_duplicates_stop_unless_accepted(self):
  dup=screen(self.tmp.name,exact=1)
  reject(self,"exact cross-cohort duplicate","exact cross-cohort duplicates",lambda:run_compare(self.roots,self.rules,duplicates=dup))
  out=run_compare(self.roots,self.rules,duplicates=dup,accept=hashlib.sha256(dup.read_bytes()).hexdigest())
  self.assertTrue(out["validation"]["valid"]); self.assertTrue(out["duplicate_screen"]["accepted_by_sha256"])
 def test_missing_checker_receipt_stops(self):
  receipts={n:str(fx.write_receipt(r)) for n,r in self.roots.items() if n!="matched2k"}
  reject(self,"missing checker receipt","one root and one checker receipt",lambda:cc.compare({n:str(r) for n,r in self.roots.items()},receipts,str(screen(self.tmp.name)),fx.SOURCE,self.rules,
         expected_foldgen=fx.foldgen_sha256()))

class FixTests(Base):
 def test_f3_receipts_must_match_the_full_cohorts_labels_on_the_added_genomes(self):
  # E1b: K accepted a receipt contradicting the full cohort's labels; W checked names only.
  full=self.fresh("f3-full","full"); relabel(full,"amoxicillin","g002003",True)
  reject(self,"F3 receipt vs full labels","class-gate receipt (0, 60) differs from the full cohort's labels on its genomes (1, 59)",
         lambda:run_compare(self.with_root("full",full),self.rules))
 def test_f5_receipts_must_prove_their_cohort_passed(self):
  for label,edit,match in [
    ("placeholder {}",lambda r:"{}","field mismatch"),
    ("failed checks",lambda r:{**r,"valid":False,"errors":["auditor found issues"]},"did not pass"),
    ("another cohort's receipt",lambda r:{**r,"cohort":"prefix"},"not 'added'"),
    ("an older run manifest",lambda r:{**r,"run_manifest_sha256":"0"*64},"another run manifest"),
    ("other checker blobs",lambda r:{**r,"pins":{**r["pins"],"check-results.py":"0"*64}},"reviewed pins"),
    ("revision 7",lambda r:{**r,"evaluation_revision":7},"revision 8"),
    ("a report count that differs",lambda r:{**r,"reports":13},"count differs"),
    ("biological execution claimed",lambda r:{**r,"biological_execution":True},"biological")]:
   with self.subTest(label):
    reject(self,"F5 "+label,match,lambda edit=edit:run_compare(self.roots,self.rules,edit={"added":edit}))
  receipt=fx.write_receipt(self.roots["added"]); out=json.loads(receipt.read_text())["outputs"]["audit"]["path"]
  Path(out).write_text("{}\n")
  receipts={n:str(receipt) if n=="added" else str(fx.write_receipt(r)) for n,r in self.roots.items()}
  reject(self,"F5 output changed after the checks","changed after the checks",lambda:cc.compare({n:str(r) for n,r in self.roots.items()},receipts,str(screen(self.tmp.name)),fx.SOURCE,self.rules,
         expected_foldgen=fx.foldgen_sha256()))
 def test_f6_matched2k_must_be_the_frozen_selection(self):
  # E1: W accepted the prefix root as matched2k.
  saved=fx.cohort_genomes
  fx.cohort_genomes=lambda name:saved("prefix") if name=="matched2k" else saved(name)
  try:
   rules=fx.make_rules(); root=fx.build_cohort(Path(self.tmp.name)/"f6-matched","matched2k",rules)
  finally: fx.cohort_genomes=saved
  reject(self,"F6 matched2k is the prefix","matched_selection",lambda:run_compare(self.with_root("matched2k",root),rules))
 def test_f4_and_f7_flags_and_no_constant_claims(self):
  out=self.baseline; text=json.dumps(out)
  # F7: W's constant `fold_changes_recorded: True` is gone; the accounting is computed.
  self.assertNotIn("fold_changes_recorded",text)
  for r in out["analyses"]["full_vs_prefix"]["rows"]:
   if r["drug"]=="amoxicillin": self.assertFalse(r["before_after_interpretable"]); self.assertIn("added genomes",r["not_interpretable_reason"])
  for analysis in ("full_vs_prefix","full_vs_matched2k"):
   for r in out["analyses"][analysis]["rows"]:
    for comp in r["comparisons"]:
     if comp["status"]=="available":
      pairing=comp["pairing"]; self.assertEqual(pairing["fold_mismatch_count"],len(pairing["fold_mismatch"]))
      self.assertEqual(pairing["eligible_matched_pairs"],comp["matched_fold"]["pairs"])
      fold_receipts=r["full"]["learners"][vc.LEARNERS.index(comp["learner"])]["folds"]
      self.assertEqual(sum(f["held_out"] for f in fold_receipts),r["full"]["learners"][vc.LEARNERS.index(comp["learner"])]["called"]+r["full"]["learners"][vc.LEARNERS.index(comp["learner"])]["uncalled"])
 def test_n6_an_unavailable_primary_split_is_listed_and_stops_for_review(self):
  full=self.fresh("n6-full","full"); did=vc.DRUG_IDS["amikacin"]
  def mutate(r):
   split=r["splits"][1]; split["learners"]=[]; split.pop("test_fold"); split["unavailable_reason"]="too few lineage groups"
  fx.mutate_report(full,did,mutate)
  code,out,_=run_main(self.with_root("full",full),self.tmp.name,self.rules)
  self.assertEqual(code,3); self.assertEqual(out["primary_split_unavailable"],{"amikacin":[{"cohort":"full","reason":"too few lineage groups"}]})
  self.assertNotIn("amikacin",[r["drug"] for r in out["size_matched_decomposition"]])
  self.assertEqual(self.row(out,"full_vs_prefix","amikacin","lineage")["status"],"unavailable")

def failed_audit(outputs):
 """What the runner's outputs hold after a failed audit: the wrapper invalid, one issue."""
 summary=outputs["audit"]["summary"]; summary["issue_count"]=1
 summary["revision8_wrapper"].update(valid=False,errors=["auditor found metric/fold/catalog/probe issues"])
 outputs["audit"]["records"][0]["issues"]=[{"where":"splits[0].test_fold","what":"saved fold differs from foldgen"}]

class FollowupTests(Base):
 def test_l2_a_forged_receipt_over_a_failed_audit_stops(self):
  # The review's repro: receipt valid=true and every hash right, but the audit output it binds
  # records the failure. compare exited 0; it must stop, with no output.
  code,out,err=run_main(self.roots,self.tmp.name,self.rules,outputs_edit={"full":failed_audit})
  self.assertEqual(code,2); self.assertIsNone(out); self.assertIn("do not record a pass",err); self.assertIn("valid is False",err)
  NEGATIVE.append({"control":"L2 forged receipt over a failed audit","result":"stop, exit 2, no output"})
 def test_l2_each_way_the_bound_outputs_can_disagree_with_a_pass_stops(self):
  def wrapper(**kw): return lambda o:o["audit"]["summary"]["revision8_wrapper"].update(**kw)
  def summary(**kw): return lambda o:o["audit"]["summary"].update(**kw)
  def record(**kw): return lambda o:o["audit"]["records"][0].update(**kw)
  cases=[
   ("wrapper errors listed",wrapper(errors=["x"]),"audit wrapper errors"),
   ("wrapper bound to another run manifest",wrapper(run_manifest_sha256="0"*64),"run_manifest_sha256"),
   ("wrapper for another source",wrapper(run_source_commit="2"*40),"run_source_commit"),
   ("wrapper for another binary",wrapper(binary_sha256="0"*64),"binary_sha256"),
   ("wrapper for another foldgen",wrapper(foldgen_sha256="0"*64),"foldgen_sha256"),
   ("wrapper ran other checker pins",wrapper(checker_pins={**vc.PINS,"check-results.py":"0"*64}),"checker_pins"),
   ("no wrapper block",lambda o:o["audit"]["summary"].pop("revision8_wrapper"),"no revision8_wrapper"),
   ("an audit issue counted",summary(issue_count=2),"issue"),
   ("an audited report missing",lambda o:o["audit"]["records"].pop(),"record count differs"),
   ("an audited record for another report",record(id="mycobacterium_tuberculosis__other"),"ids differ"),
   ("an audited record of other bytes",record(sha256="0"*64),"SHA-256 differs from the run manifest"),
   ("an audited record with issues",record(issues=[{"where":"x","what":"y"}]),"issues present"),
   ("revision 7 records",summary(by_revision={"7":22}),"wrong revisions"),
   ("checker reports missing",lambda o:o["checker"]["reports"].pop(),"checker report count"),
   ("summarizer issues",lambda o:o["summary"]["coverage"]["8"].update(issues=1),"summarizer coverage"),
  ]
  for label,edit,match in cases:
   with self.subTest(label):
    reject(self,"L2 "+label,match,lambda edit=edit:run_compare(self.roots,self.rules,outputs_edit={"prefix":edit}))
 def test_l2_the_runners_own_outputs_after_a_failed_audit_stop_compare(self):
  # A real failed run: the foldgen rebuild trap (two same-class random folds swapped). The runner
  # writes no receipt; a receipt forged over its real outputs is refused.
  root=self.fresh("l2-real-trap","prefix"); did=vc.DRUG_IDS["amikacin"]; path=report_path(root,did); report=json.loads(path.read_text()); split=report["splits"][0]
  i=next(i for i,f in enumerate(split["test_fold"]) if f==0 and report["samples"][i]["resistant"]); j=next(j for j,f in enumerate(split["test_fold"]) if f==1 and report["samples"][j]["resistant"])
  split["test_fold"][i],split["test_fold"][j]=split["test_fold"][j],split["test_fold"][i]; fx.rescore(report); fx.write_json(path,report); rehash(root,did); fx.refresh_root(root,"prefix")
  args=runner_args(root,"real")
  with contextlib.redirect_stdout(io.StringIO()),contextlib.redirect_stderr(io.StringIO()): self.assertEqual(runner.main(args,rules=self.rules),2)
  self.assertFalse((root/"real-receipt.json").exists())
  outputs={"checker":root/"real-checker.json","audit":root/"real-audit.json","summary":root/"real-summary.json"}
  if not outputs["summary"].exists(): outputs["summary"].write_text(json.dumps({"coverage":{"8":{"files":22,"issues":0}}})+"\n")
  forged=fx.write_receipt(root,edit=lambda r:{**r,"outputs":{k:{"path":str(v),"sha256":fx.sha(v)} for k,v in outputs.items()}})
  receipts={n:str(forged) if n=="prefix" else str(fx.write_receipt(r)) for n,r in self.roots.items()}
  reject(self,"L2 receipt forged over the runner's real failed outputs","do not record a pass",
         lambda:cc.compare({**{n:str(r) for n,r in self.roots.items()},"prefix":str(root)},receipts,str(screen(self.tmp.name)),fx.SOURCE,self.rules,
                           expected_foldgen=fx.foldgen_sha256()))
 def test_l3_any_invalid_report_in_any_cohort_stops(self):
  for name in fx.COHORTS:
   with self.subTest(name):
    root=self.fresh(f"l3-{name}",name); did=vc.DRUG_IDS["isoniazid"]
    fx.mutate_report(root,did,lambda r:r["evaluation"].update(valid=False,reasons=["unconverged"]))
    reject(self,f"L3 invalid report in {name}","evaluation.valid=true",lambda root=root,name=name:run_compare(self.with_root(name,root),self.rules))
 def test_i4_a_receipt_for_another_foldgen_build_stops(self):
  reject(self,"I4 receipt foldgen is not the recorded build","not the recorded build",lambda:run_compare(self.roots,self.rules,foldgen="0"*64))
  reject(self,"I4 malformed expected foldgen","malformed digest",lambda:run_compare(self.roots,self.rules,foldgen="not-a-digest"))
  with contextlib.redirect_stderr(io.StringIO()) as err, self.assertRaises(SystemExit) as stop:
   cc.main([f"--cohort={n}={r}" for n,r in self.roots.items()]+["--duplicates","x","--expected-source-commit",fx.SOURCE,"--out","x"])
  self.assertEqual(stop.exception.code,2); self.assertIn("--expected-foldgen-sha256",err.getvalue())
 def test_i9_a_matched_control_with_no_pair_is_invalid(self):
  order=["g1","g2","g3","g4"]; truth={"g1":True,"g2":False,"g3":True,"g4":False}
  prefix=pairing_report(order,truth,{"g1":0,"g2":1,"g3":2,"g4":3},"lineage")
  full=pairing_report(order,truth,{"g1":1,"g2":0,"g3":3,"g4":2},"lineage")
  row,issues=cc.pair_learner(prefix,full,"lineage","cart",set(order),{g:True for g in order})
  self.assertEqual(row["pairing"]["eligible_matched_pairs"],0); self.assertEqual(issues,[])
  self.assertFalse(row["matched_fold"]["control_valid"]); self.assertIn("no common genome",row["matched_fold"]["control_invalid_reason"])
  row,_=cc.pair_learner(prefix,prefix,"lineage","cart",set(order),{g:True for g in order})
  self.assertTrue(row["matched_fold"]["control_valid"]); self.assertNotIn("control_invalid_reason",row["matched_fold"])
  NEGATIVE.append({"control":"I9 lineage matched control with no pair","result":"control_valid false"})

def pairing_report(order,truth,folds,scheme):
 """A minimal ReportData for pair_learner: one learner, given folds, fixed scores."""
 learner=vc.LearnerData("cart",scheme,{g:.7 if truth[g] else .3 for g in order},dict(folds),{})
 split=vc.SplitData(scheme,True,None,dict(folds),{"cart":learner},{})
 return vc.ReportData("x","amikacin",Path("x"),"0"*64,"0"*16,{},0.0,{g:vc.Sample(g,truth[g],0) for g in order},list(order),{scheme:split},{})

def runner_args(root,tag):
 return [str(root/f"{tag}-audit.json"),"--root",str(root),"--run-manifest",str(root/"validation-manifest.json"),"--expected-source-commit",fx.SOURCE,
         "--binary",str(root/"synthetic-binary"),"--backup",str(root/"backup"),"--log",str(root/"pipeline.log"),"--launch","2000-01-01T00:00:00Z",
         "--foldgen",fx.foldgen_path(),"--expected-foldgen-sha256",fx.foldgen_sha256(),"--checker-out",str(root/f"{tag}-checker.json"),
         "--summary-out",str(root/f"{tag}-summary.json"),"--receipt",str(root/f"{tag}-receipt.json")]

class EndToEndRealFoldsTests(unittest.TestCase):
 def test_e2_real_foldgen_receipts_and_comparison_pass(self):
  with tempfile.TemporaryDirectory() as tmp:
   rules=fx.make_rules(); roots={n:fx.build_cohort(Path(tmp)/n,n,rules) for n in fx.COHORTS}; receipts={}
   for name,root in roots.items():
    # Path-adapted (I4): the runner requires the recorded foldgen build.
    with contextlib.redirect_stdout(io.StringIO()): self.assertEqual(runner.main(runner_args(root,"real"),rules=rules),0,name)
    receipts[name]=str(root/"real-receipt.json")
   # L2: the comparator re-reads the runner's real outputs and accepts them.
   out=cc.compare({n:str(r) for n,r in roots.items()},receipts,str(screen(tmp)),fx.SOURCE,rules,expected_foldgen=fx.foldgen_sha256())
   self.assertTrue(out["validation"]["valid"],out["validation"]["blocking_issues"])
   rates=[v["fold_mismatch"]/v["common"] for k,v in out["random_fold_accounting"].items() if k.startswith("full_vs_prefix/")]
   self.assertTrue(rates and max(rates)>0.2,rates)  # real folds differ, as grm's do
   self.assertEqual(len({r["foldgen_sha256"] for r in out["checker_receipts"].values()}),1)

if __name__=="__main__":
 result=unittest.main(exit=False,verbosity=2).result
 print(json.dumps({"negative_controls":len(NEGATIVE)}))
 sys.exit(0 if result.wasSuccessful() else 1)
