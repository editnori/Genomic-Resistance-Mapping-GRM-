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
- EndToEndRealFoldsTests: experiment E2 as a passing run: run_revision8_checks.py with the real
  foldgen on all four cohorts, then the comparator on their receipts.
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

def run_compare(roots,rules,source=fx.SOURCE,receipt_source=None,duplicates=None,accept=None,edit=None):
 receipts={n:str(fx.write_receipt(r,source=receipt_source or source,edit=(edit.get(n) if edit else None))) for n,r in roots.items()}
 return cc.compare({n:str(r) for n,r in roots.items()},receipts,str(duplicates or screen(Path(roots["prefix"]).parent)),source,rules,accept)

def run_main(roots,base,rules,extra=()):
 receipts=[f"--checker-receipt={n}={fx.write_receipt(r)}" for n,r in roots.items()]
 out=Path(base)/f"out-{len(list(Path(base).glob('out-*')))}.json"
 argv=[f"--cohort={n}={r}" for n,r in roots.items()]+receipts+["--duplicates",str(screen(base)),"--expected-source-commit",fx.SOURCE,"--out",str(out),*extra]
 saved=cc.vc.RULES
 cc.vc.RULES=rules
 try:
  with contextlib.redirect_stdout(io.StringIO()): code=cc.main(argv)
 finally: cc.vc.RULES=saved
 return code,(json.loads(out.read_text()) if out.exists() else None)

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
  fx.rescore(report); fx.write_json(path,report); rehash(full,did); out=run_compare(self.with_root("full",full),self.rules)
  # CHANGED ASSERTION: the merged contract checks report counts against the root's own
  # phenotypes, so this report is now invalid at load and the drug is excluded where K reached a
  # coverage issue. K's coverage assertions are kept on pair_learner itself below.
  self.assertEqual(out["excluded"]["amikacin"][0]["cohort"],"full"); self.assertIn("root's phenotypes",out["excluded"]["amikacin"][0]["reasons"][0])
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
  _,issues=cc.compare_prefix_added(p,a,set(),{}); self.assertTrue(any("not disjoint" in x["issue"] for x in issues)); NEGATIVE.append({"control":"non-disjoint prefix/added","result":"blocked"})
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
 def test_invalid_in_one_cohort_is_excluded_there_and_stops_for_review(self):
  matched=self.fresh("w-invalid-matched","matched2k"); did=vc.DRUG_IDS["amikacin"]
  fx.mutate_report(matched,did,lambda r:r["evaluation"].update(valid=False,reasons=["final model has an unconverged two-class lasso fit"]))
  code,out=run_main(self.with_root("matched2k",matched),self.tmp.name,self.rules)
  self.assertEqual(code,3); self.assertEqual(list(out["excluded"]),["amikacin"]); self.assertEqual(out["excluded"]["amikacin"][0]["cohort"],"matched2k")
  drugs=lambda name:{r["drug"] for r in out["analyses"][name]["rows"]}  # noqa: E731
  self.assertIn("amikacin",drugs("full_vs_prefix")); self.assertIn("amikacin",drugs("added_vs_prefix"))
  self.assertNotIn("amikacin",drugs("matched2k_vs_prefix")); self.assertNotIn("amikacin",drugs("full_vs_matched2k"))
  self.assertNotIn("amikacin",[r["drug"] for r in out["size_matched_decomposition"]])
  NEGATIVE.append({"control":"invalid in one cohort","result":"excluded there, exit 3"})
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
  code,out=run_main(self.with_root("added",added),self.tmp.name,self.rules)
  self.assertEqual(code,3); self.assertIn("recomputed",out["excluded"]["amikacin"][0]["reasons"][0])
 def test_exact_cross_cohort_duplicates_stop_unless_accepted(self):
  dup=screen(self.tmp.name,exact=1)
  reject(self,"exact cross-cohort duplicate","exact cross-cohort duplicates",lambda:run_compare(self.roots,self.rules,duplicates=dup))
  out=run_compare(self.roots,self.rules,duplicates=dup,accept=hashlib.sha256(dup.read_bytes()).hexdigest())
  self.assertTrue(out["validation"]["valid"]); self.assertTrue(out["duplicate_screen"]["accepted_by_sha256"])
 def test_missing_checker_receipt_stops(self):
  receipts={n:str(fx.write_receipt(r)) for n,r in self.roots.items() if n!="matched2k"}
  reject(self,"missing checker receipt","one root and one checker receipt",lambda:cc.compare({n:str(r) for n,r in self.roots.items()},receipts,str(screen(self.tmp.name)),fx.SOURCE,self.rules))

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
  reject(self,"F5 output changed after the checks","changed after the checks",lambda:cc.compare({n:str(r) for n,r in self.roots.items()},receipts,str(screen(self.tmp.name)),fx.SOURCE,self.rules))
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
  code,out=run_main(self.with_root("full",full),self.tmp.name,self.rules)
  self.assertEqual(code,3); self.assertEqual(out["primary_split_unavailable"],{"amikacin":[{"cohort":"full","reason":"too few lineage groups"}]})
  self.assertNotIn("amikacin",[r["drug"] for r in out["size_matched_decomposition"]])
  self.assertEqual(self.row(out,"full_vs_prefix","amikacin","lineage")["status"],"unavailable")

class EndToEndRealFoldsTests(unittest.TestCase):
 def test_e2_real_foldgen_receipts_and_comparison_pass(self):
  with tempfile.TemporaryDirectory() as tmp:
   rules=fx.make_rules(); roots={n:fx.build_cohort(Path(tmp)/n,n,rules) for n in fx.COHORTS}; receipts={}
   for name,root in roots.items():
    args=[str(root/"audit.json"),"--root",str(root),"--run-manifest",str(root/"validation-manifest.json"),"--expected-source-commit",fx.SOURCE,
          "--binary",str(root/"synthetic-binary"),"--backup",str(root/"backup"),"--log",str(root/"pipeline.log"),"--launch","2000-01-01T00:00:00Z",
          "--foldgen",fx.foldgen_path(),"--checker-out",str(root/"checker.json"),"--summary-out",str(root/"summary.json"),"--receipt",str(root/"real-receipt.json")]
    with contextlib.redirect_stdout(io.StringIO()): self.assertEqual(runner.main(args,rules=rules),0,name)
    receipts[name]=str(root/"real-receipt.json")
   out=cc.compare({n:str(r) for n,r in roots.items()},receipts,str(screen(tmp)),fx.SOURCE,rules)
   self.assertTrue(out["validation"]["valid"],out["validation"]["blocking_issues"])
   rates=[v["fold_mismatch"]/v["common"] for k,v in out["random_fold_accounting"].items() if k.startswith("full_vs_prefix/")]
   self.assertTrue(rates and max(rates)>0.2,rates)  # real folds differ, as grm's do
   self.assertEqual(len({r["foldgen_sha256"] for r in out["checker_receipts"].values()}),1)

if __name__=="__main__":
 result=unittest.main(exit=False,verbosity=2).result
 print(json.dumps({"negative_controls":len(NEGATIVE)}))
 sys.exit(0 if result.wasSuccessful() else 1)
