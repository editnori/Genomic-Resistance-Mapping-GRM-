#!/usr/bin/env python3
"""Validate and compare the matched MTB revision-8 cohorts (stdlib only; saved scores and metadata only).

usage: python3 -I compare_cohorts.py
           --cohort prefix=ROOT --cohort added=ROOT --cohort full=ROOT --cohort matched2k=ROOT
           --checker-receipt prefix=PATH ... (one per cohort, from run_revision8_checks.py)
           --duplicates SCREEN.json [--accept-reviewed-duplicates SHA256]
           --expected-source-commit SHA --expected-foldgen-sha256 SHA256 --out OUT.json

ROOT is a cohort's execution root (RUN/cohorts/<name>): validation-manifest.json from
write_run_manifest.py, beside data/. The canonical merge of packages K and W (RECONCILE.md): K's
strict loader and accounting (report_contract.py), W's four-cohort design.

1. Loading. Every manifest, input hash, report hash and fingerprint, and every report, is
   checked strictly by report_contract.validate_manifest; any of those failing, a receipt that
   does not prove its cohort passed (F5), a duplicate screen with an exact cross-cohort duplicate
   not accepted by its SHA-256, or cohorts that are not the frozen prefix, added, full and
   matched_selection(full) lists (F6) stop the run with no output (exit 2). A receipt proves its
   cohort passed only when the checker, audit and summary outputs it names still hash as it
   says and, read again, record that pass for this run manifest, source, binary and foldgen
   (L2, run_revision8_checks.verify_outputs), and its foldgen is the build the plan recorded
   (--expected-foldgen-sha256; I4).
2. Class gates and labels. The full cohort's report labels on each other cohort's genomes must
   give that cohort's report counts and its class-gate receipts exactly (F3): the added
   cohort's eight receipts and matched2k's. From them, the drugs not interpretable for
   before/after claims (fewer than 25 of a class among the added genomes) must be the frozen
   eight, and the added cohort's eligible drugs the frozen 14.
3. Any invalid report is a stop (L3). W's rule excluded a drug invalid in one cohort and exited
   3, but nothing upstream can produce that state: revision-8 `grm learn` stops a species at its
   first invalid report, and run_revision8_checks.py, around the reviewed in-tree checker, gives
   a cohort with any invalid report no receipt. A per-drug exclusion would need a per-drug
   checker, which the plan forbids. So an invalid or missing eligible report stops here too
   (exit 2, no output), and the parent decides.
   A drug whose primary (lineage) split is unavailable in a cohort has no primary result
   there: listed under `primary_split_unavailable`, left out of the decomposition; exit 3.
4. Comparisons, on the primary split (lineage) and the secondary (random), each row with
   recomputed metrics, Clopper-Pearson 95% intervals for VME and ME and, from the same
   computation, sensitivity and specificity, and the null-aware five-learner medians:
   - full vs prefix and full vs matched2k, paired by genome ID, never row position:
     `genome_aligned_all_folds` transitions (fixed, broken, no-call moves, score deltas) with an
     exact McNemar p-value on fixed vs broken, Holm-adjusted within learner x split across the
     interpretable drugs (descriptive); and `matched_fold`, the same restricted to genomes held
     out in the same fold and, for lineage, an unchanged group relation. grm's random folds are
     a seeded shuffle over each cohort's own samples, so they differ between cohorts: the
     mismatch is counted per row and marks that row's matched random control invalid, but does
     not invalidate the comparison (decision D1, F2). Fold integrity itself is enforced per
     cohort by run_revision8_checks.py, which rebuilds every saved fold with the real foldgen.
   - added vs prefix and matched2k vs prefix, side by side (disjoint or independent cohorts);
   - the full cohort's CV scores on the added isolates beside added-only CV (not frozen-model
     scoring);
   - the size-matched decomposition of the primary split's median balanced accuracy:
     total (full - prefix) = composition (matched2k - prefix) + size (full - matched2k).
   The eight not-interpretable drugs are flagged `before_after_interpretable: false` in every
   comparison with full and kept out of Holm (F4). Per-fold receipts are recomputed from the
   saved scores (F7: no constant claims).
5. Blocking issues (output written, exit 1): source commit, binary, foldgen or EvalConfig
   differing across cohorts or from the expected source; genome coverage or label mismatch in a
   paired comparison; prefix and added not disjoint; full not the exact prefix + added.

Exit status: 0 complete; 1 output written but the comparison is invalid; 2 input or integrity
stop, no output, including when the validators cannot load (a moved pin); 3 output written and
valid, but a primary split is unavailable (stop for review before any before/after claim).
"""
import argparse, importlib, json, math, statistics, sys
from collections import defaultdict
from pathlib import Path
# M1: python3 -I (3.11+ implies -P) leaves this script's directory off sys.path. Put back that one
# directory, first, and nothing else, and require each sibling module to load from it.
LANE=Path(__file__).resolve().parent
if sys.path[:1]!=[str(LANE)]: sys.path.insert(0,str(LANE))
def _sibling(name):
 module=importlib.import_module(name)
 if Path(module.__file__).resolve().parent!=LANE: raise ImportError(f"{name} loaded from {module.__file__}, not from {LANE}")
 return module
try: vc=_sibling("report_contract"); prep=_sibling("prepare_execution_roots"); runner=_sibling("run_revision8_checks")
except (ImportError,ValueError,OSError) as _exc:  # ValueError, OSError: report_contract's pin check
 print(f"compare_cohorts: STOP: cannot load the validators: {_exc}",file=sys.stderr); sys.exit(2)

COHORTS=("prefix","added","full","matched2k")
PRIMARY="lineage"
MIN_CLASS=25
RECEIPT_SCHEMA="grm-mtb-4k-revision8-check-receipt-v1"
EXPECTED_ADDED_ELIGIBLE=set(vc.ADDED_DRUGS)
EXPECTED_SPARSE_ADDED=set(vc.ADDED_UNAVAILABLE)

def mcnemar_exact(b,c):
 """Two-sided exact McNemar p-value from discordant counts b and c (W)."""
 n=b+c
 if n==0: return 1.0
 k=min(b,c)
 return min(1.0,2*sum(math.comb(n,i) for i in range(k+1))/2**n)

def holm(pvalues):
 """Holm step-down adjustment of {key: p} (W)."""
 order=sorted(pvalues,key=lambda key:pvalues[key]); m,running,out=len(order),0.0,{}
 for i,key in enumerate(order):
  running=max(running,min(1.0,(m-i)*pvalues[key])); out[key]=running
 return out

def subset_row(report,scheme,name,ids=None):
 split=report.splits[scheme]
 if not split.available: return None
 learner=split.learners[name]
 if ids is None: return learner.receipt
 selected=[g for g in report.sample_order if g in ids]
 truth=[report.samples[g].resistant for g in selected]
 scores=[learner.scores[g] for g in selected]
 metrics=vc.upstream_audit.metrics(truth,scores,name in vc.RANKED)
 no_call={"resistant":sum(v is None and t for v,t in zip(scores,truth)),"susceptible":sum(v is None and not t for v,t in zip(scores,truth))}
 folds=[]; k=report.eval_config["folds"]
 for fold in range(k):
  idx=[i for i,g in enumerate(selected) if learner.folds[g]==fold]
  called=[i for i in idx if scores[i] is not None]
  cm=vc.upstream_audit.metrics([truth[i] for i in called],[scores[i] for i in called],False)
  folds.append({"fold":fold,"status":"not_in_subset" if not idx else "present","held_out":len(idx),
   "resistant":sum(truth[i] for i in idx),"susceptible":len(idx)-sum(truth[i] for i in idx),"called":len(called),
   "called_resistant":sum(truth[i] for i in called),"called_susceptible":len(called)-sum(truth[i] for i in called),
   "no_call":{"resistant":sum(scores[i] is None and truth[i] for i in idx),"susceptible":sum(scores[i] is None and not truth[i] for i in idx)},
   "metrics":cm})
 return {"learner":name,"scheme":scheme,"metrics":metrics,"no_call":no_call,"called":sum(v is not None for v in scores),
  "uncalled":sum(v is None for v in scores),"no_call_rate":vc.upstream_audit.ratio(sum(v is None for v in scores),len(scores)),
  "folds":folds,"clopper_pearson_95":vc.error_intervals(metrics),"subset_of_report":True}

def report_summary(report,scheme,ids=None):
 split=report.splits[scheme]
 if not split.available:
  return {"status":"unavailable","unavailable_reason":split.reason,"learners":[],"five_learner_medians":None}
 rows=[subset_row(report,scheme,name,ids) for name in vc.LEARNERS]
 return {"status":"available","learners":rows,"five_learner_medians":vc.learner_medians(rows)}

def median_ba(report,scheme):
 split=report.splits[scheme]
 return split.receipt["five_learner_medians"]["balanced_accuracy"]["median"] if split.available else None

def transitions(pairs):
 out={"pairs":len(pairs),"both_called":0,"fixed":0,"broken":0,"unchanged_correct":0,"unchanged_wrong":0,
      "prefix_no_call_full_call":0,"prefix_call_full_no_call":0,"both_no_call":0}
 deltas=[]
 for truth,p_score,f_score in pairs:
  if p_score is None and f_score is None: out["both_no_call"]+=1; continue
  if p_score is None: out["prefix_no_call_full_call"]+=1; continue
  if f_score is None: out["prefix_call_full_no_call"]+=1; continue
  out["both_called"]+=1; deltas.append(f_score-p_score)
  pc=(p_score>=0.5)==truth; fc=(f_score>=0.5)==truth
  if not pc and fc: out["fixed"]+=1
  elif pc and not fc: out["broken"]+=1
  elif pc and fc: out["unchanged_correct"]+=1
  else: out["unchanged_wrong"]+=1
 out["score_delta_full_minus_prefix"]={"n":len(deltas),"mean":statistics.fmean(deltas) if deltas else None,
  "median":statistics.median(deltas) if deltas else None,"min":min(deltas) if deltas else None,
  "max":max(deltas) if deltas else None,"mean_absolute":statistics.fmean(map(abs,deltas)) if deltas else None}
 out["mcnemar_p"]=mcnemar_exact(out["fixed"],out["broken"])
 out["definition"]="fixed=base incorrect/full correct; broken=base correct/full incorrect at score>=0.5; no-call transitions are separate; 'prefix' names the base cohort"
 return out

def stable_lineage_genomes(base,full):
 common=set(base.genomes); p_to_f=defaultdict(set); f_to_p=defaultdict(set)
 for g in common: p_to_f[base.lineage[g]].add(full.lineage[g]); f_to_p[full.lineage[g]].add(base.lineage[g])
 return {g:(len(p_to_f[base.lineage[g]])==1 and len(f_to_p[full.lineage[g]])==1) for g in common}

def pair_learner(prefix,full,scheme,name,prefix_genomes,stable_groups):
 ps,fs=prefix.splits[scheme],full.splits[scheme]
 base={"learner":name,"scheme":scheme,"status":"available","genome_aligned_all_folds":None,"matched_fold":None}
 if not ps.available or not fs.available:
  base.update(status="unavailable",prefix_reason=ps.reason,full_reason=fs.reason,
              note="Unavailable split remains unavailable; no zero scores or error counts are imputed")
  return base,[]
 pids=set(prefix.samples); fids=set(full.samples); common=[g for g in prefix.sample_order if g in fids]
 prefix_only=[g for g in prefix.sample_order if g not in fids]
 full_only_prefix=sorted(g for g in fids-pids if g in prefix_genomes)
 label_mismatch=[g for g in common if prefix.samples[g].resistant!=full.samples[g].resistant]
 pl,fl=ps.learners[name],fs.learners[name]
 fold_mismatch=[g for g in common if pl.folds[g]!=fl.folds[g]]
 group_mismatch=[g for g in common if scheme=="lineage" and not stable_groups[g]] if scheme=="lineage" else []
 eligible=[g for g in common if g not in set(fold_mismatch) and g not in set(group_mismatch)]
 all_pairs=[(prefix.samples[g].resistant,pl.scores[g],fl.scores[g]) for g in common]
 matched_pairs=[(prefix.samples[g].resistant,pl.scores[g],fl.scores[g]) for g in eligible]
 base["pairing"]={"prefix_samples":len(pids),"full_samples":len(fids),"common_samples":len(common),
  "coverage_of_prefix":vc.upstream_audit.ratio(len(common),len(pids)),"prefix_only":prefix_only,
  "full_only_from_prefix_cohort":full_only_prefix,"label_mismatch":label_mismatch,"fold_mismatch":fold_mismatch,
  "fold_mismatch_count":len(fold_mismatch),"changed_lineage_group":group_mismatch,
  "changed_lineage_group_count":len(group_mismatch),"eligible_matched_pairs":len(eligible),
  "eligible_coverage_of_common":vc.upstream_audit.ratio(len(eligible),len(common)),
  "alignment":"genome ID; matched_fold additionally requires equal held-out fold and, for lineage, an unchanged common-genome group relation"}
 base["genome_aligned_all_folds"]=transitions(all_pairs)
 base["genome_aligned_all_folds"]["interpretation"]="descriptive across independently trained CV runs; different folds are included and counted above"
 base["matched_fold"]=transitions(matched_pairs)
 base["matched_fold"]["interpretation"]="same genome and held-out fold only; lineage pairs also require an unchanged group relation"
 # D1 (F2): a random-fold mismatch is expected (each cohort shuffles its own samples). It is
 # counted above and invalidates this row's matched random control, not the comparison.
 # I9: a matched control with no eligible pair holds nothing to compare, for either split.
 control_valid=not (scheme=="random" and fold_mismatch) and bool(eligible)
 base["matched_fold"]["control_valid"]=control_valid
 if scheme=="random" and fold_mismatch: base["matched_fold"]["control_invalid_reason"]=f"random folds differ for {len(fold_mismatch)} of {len(common)} common genomes"
 elif not eligible:
  base["matched_fold"]["control_invalid_reason"]="no common genome is held out in the same fold"+(" with an unchanged lineage group" if scheme=="lineage" else "")
 issues=[]
 if prefix_only or full_only_prefix or label_mismatch: issues.append("genome coverage/label mismatch")
 return base,issues

def paired_rows(base_cohort,full,base_name,eligible,flagged):
 rows=[]; issues=[]; stable=stable_lineage_genomes(base_cohort,full)
 for drug in sorted(eligible):
  did=vc.DRUG_IDS[drug]
  b,f=base_cohort.reports[did],full.reports[did]
  for scheme in vc.SCHEMES:
   row={"drug":drug,"dataset_id":did,"scheme":scheme,"primary":scheme==PRIMARY,"before_after_interpretable":drug not in flagged,
        base_name:report_summary(b,scheme),"full":report_summary(f,scheme),"comparisons":[]}
   if drug in flagged: row["not_interpretable_reason"]=flagged[drug]
   for name in vc.LEARNERS:
    comparison,problems=pair_learner(b,f,scheme,name,set(base_cohort.genomes),stable)
    row["comparisons"].append(comparison)
    for problem in problems: issues.append({"analysis":f"full_vs_{base_name}","drug":drug,"scheme":scheme,"learner":name,"issue":problem})
   rows.append(row)
 # Holm within learner x split, across the interpretable drugs, on the genome-aligned McNemar p-values.
 for scheme in vc.SCHEMES:
  for li,name in enumerate(vc.LEARNERS):
   ps={r["drug"]:r["comparisons"][li]["genome_aligned_all_folds"]["mcnemar_p"] for r in rows
       if r["scheme"]==scheme and r["before_after_interpretable"] and r["comparisons"][li]["genome_aligned_all_folds"]}
   for drug,p in holm(ps).items():
    next(r for r in rows if r["drug"]==drug and r["scheme"]==scheme)["comparisons"][li]["genome_aligned_all_folds"]["mcnemar_p_holm"]=p
 return rows,issues

def compare_paired(base_cohort,full,base_name,eligible,flagged):
 rows,issues=paired_rows(base_cohort,full,base_name,eligible,flagged)
 return {"analysis_type":f"full_vs_{base_name}_internal_cv","paired_by":"genome ID and held-out-fold metadata, never row position",
  "multiplicity":"exact McNemar on genome-aligned fixed vs broken; Holm within learner x split across interpretable drugs; descriptive",
  "external_validation":False,"rows":rows},issues

def compare_side_by_side(a,b,a_name,b_name,eligible,flagged):
 rows=[]
 for drug in sorted(eligible):
  did=vc.DRUG_IDS[drug]; ra,rb=a.reports[did],b.reports[did]
  for scheme in vc.SCHEMES:
   ma,mb=median_ba(ra,scheme),median_ba(rb,scheme)
   row={"drug":drug,"dataset_id":did,"scheme":scheme,"primary":scheme==PRIMARY,a_name:report_summary(ra,scheme),b_name:report_summary(rb,scheme),
        "median_balanced_accuracy":{a_name:ma,b_name:mb,"delta":None if ma is None or mb is None else mb-ma}}
   if flagged is not None:
    row["before_after_interpretable"]=drug not in flagged
    if drug in flagged: row["not_interpretable_reason"]=flagged[drug]
   rows.append(row)
 return {"analysis_type":f"{b_name}_vs_{a_name}_side_by_side","paired_isolates":False,"frozen_model_external_validation":False,
  "note":"Independent cohort CV estimates; no isolate pairing and no external/prospective/clinical claim","rows":rows}

def compare_prefix_added(prefix,added,eligible):
 """Disjoint-cohort CV side by side (K), with K's disjointness check."""
 issues=[]
 overlap=sorted(set(prefix.genomes)&set(added.genomes))
 if overlap: issues.append({"analysis":"added_vs_prefix","issue":"cohorts are not disjoint","overlap":overlap[:20]})
 return compare_side_by_side(prefix,added,"prefix","added",eligible,None),issues

def compare_full_added(added,full,eligible):
 rows=[]; issues=[]; added_set=set(added.genomes)
 for drug in sorted(eligible):
  did=vc.DRUG_IDS[drug]; a,f=added.reports[did],full.reports[did]
  for scheme in vc.SCHEMES:
   rows.append({"drug":drug,"dataset_id":did,"scheme":scheme,"primary":scheme==PRIMARY,"added_only_cv":report_summary(a,scheme),
                "full_cohort_cv_added_subset":report_summary(f,scheme,added_set)})
 return {"analysis_type":"full_cohort_cv_scores_on_added_isolates_side_by_side_with_added_only_cv",
  "paired_isolates":False,"frozen_model_external_validation":False,
  "note":"Full rows are subsets of full-cohort held-out CV scores, not a frozen prefix/added model externally applied to added isolates. Fold denominators are subset denominators; folds with no added isolates are not_in_subset, not zero-score folds","rows":rows},issues

def load_receipt(path,cohort,expected_source,expected_foldgen):
 """F5: the run_revision8_checks.py receipt for this cohort, bound to its run manifest, source,
 binary, foldgen (I4), pins and outputs. A placeholder or a failed check is a stop, and so is a
 receipt whose outputs, read again, do not record the pass it claims (L2)."""
 path=Path(path); rec,sha=vc.read_json(path); where=str(path)
 vc.expect_keys(rec,{"schema","cohort","valid","errors","evaluation_revision","reports","class_gate_skips","run_manifest","run_manifest_sha256",
  "source_commit","binary_sha256","foldgen","foldgen_sha256","pins","validator_sha256","outputs","biological_execution"},where)
 if rec["schema"]!=RECEIPT_SCHEMA: vc.fail(where,"not a run_revision8_checks.py receipt")
 if rec["cohort"]!=cohort.name: vc.fail(where,f"receipt is for cohort {rec['cohort']!r}, not {cohort.name!r}")
 if rec["valid"] is not True or rec["errors"]!=[]: vc.fail(where,"the revision-8 checks did not pass")
 if rec["evaluation_revision"]!=vc.REVISION: vc.fail(where,"receipt is not for evaluation revision 8")
 if rec["reports"]!=len(cohort.manifest["reports"]) or rec["class_gate_skips"]!=len(cohort.manifest["class_gate_unavailable"]): vc.fail(where,"receipt report or skip count differs from the run manifest")
 if rec["run_manifest_sha256"]!=cohort.manifest_sha256: vc.fail(where,"receipt was written for another run manifest")
 if rec["source_commit"]!=expected_source or rec["binary_sha256"]!=cohort.manifest["binary_sha256"]: vc.fail(where,"receipt source or binary differs")
 vc.digest(rec["foldgen_sha256"],where+".foldgen_sha256")
 if rec["foldgen_sha256"]!=expected_foldgen: vc.fail(where,f"receipt foldgen {rec['foldgen_sha256']} is not the recorded build {expected_foldgen} (I4)")
 if rec["pins"]!=vc.PINS: vc.fail(where,"receipt ran other checker blobs than the reviewed pins")
 if not isinstance(rec["outputs"],dict) or set(rec["outputs"])!={"checker","audit","summary"}: vc.fail(where,"receipt outputs incomplete")
 for key,entry in rec["outputs"].items():
  vc.expect_keys(entry,{"path","sha256"},f"{where}.outputs.{key}")
  if vc.file_sha256(entry["path"])!=entry["sha256"]: vc.fail(where,f"{key} output changed after the checks: {entry['path']}")
 if rec["biological_execution"] is not False: vc.fail(where,"receipt claims biological execution")
 # L2: the receipt's `valid` is a claim; the outputs it binds must record the pass themselves.
 errors=runner.verify_outputs({key:Path(entry["path"]) for key,entry in rec["outputs"].items()},cohort,expected_source,
                              cohort.manifest["binary_sha256"],expected_foldgen)
 if errors: vc.fail(where,"the outputs it binds do not record a pass: "+"; ".join(errors[:10]))
 return {"path":str(path),"sha256":sha,"foldgen_sha256":rec["foldgen_sha256"],"validator_sha256":rec["validator_sha256"]}

def label_counts(full_report,genomes):
 labels=[s.resistant for g,s in full_report.samples.items() if g in genomes]
 return sum(labels),len(labels)-sum(labels)

def cross_check_labels(cohorts):
 """F3: the full cohort's labels on each cohort's genomes reproduce its report counts and its
 class-gate receipts. Returns problems (stops) and per-drug counts on the added genomes."""
 full=cohorts["full"]; problems=[]; added_counts={}
 for name in ("prefix","added","matched2k"):
  c=cohorts[name]; genomes=set(c.genomes)
  receipts={row["antibiotic"]:(row["resistant"],row["susceptible"]) for row in c.manifest["class_gate_unavailable"]}
  for drug in vc.ALL_DRUGS:
   fr=full.reports.get(vc.DRUG_IDS[drug])
   if fr is None: continue  # not eligible in full (never, under the frozen rule): nothing to cross-check
   derived=label_counts(fr,genomes)
   if name=="added": added_counts[drug]=derived
   report=c.reports.get(vc.DRUG_IDS[drug])
   if drug in receipts and receipts[drug]!=derived: problems.append(f"{name} {drug}: class-gate receipt {receipts[drug]} differs from the full cohort's labels on its genomes {derived}")
   if report is not None:
    got=(report.receipt["counts"]["resistant"],report.receipt["counts"]["susceptible"])
    if got!=derived or {g for g in report.samples}!={g for g in fr.samples if g in genomes}: problems.append(f"{name} {drug}: report labels {got} differ from the full cohort's labels on its genomes {derived}")
 return problems,added_counts

def compare(paths,receipts,duplicates,expected_source,rules=None,accept_duplicates=None,*,expected_foldgen):
 rules=rules or vc.RULES
 if set(paths)!=set(COHORTS) or set(receipts)!=set(COHORTS): raise vc.ValidationError(f"need one root and one checker receipt for each of {list(COHORTS)}")
 vc.digest(expected_source,"expected source commit",vc.COMMIT); vc.digest(expected_foldgen,"expected foldgen SHA-256")
 # L3: strict. Any invalid report stops here (validate_manifest raises), with no output.
 cohorts={name:vc.validate_manifest(Path(paths[name])/"validation-manifest.json",rules) for name in COHORTS}
 receipt_rows={name:load_receipt(receipts[name],cohorts[name],expected_source,expected_foldgen) for name in COHORTS}
 screen_path=Path(duplicates); screen,screen_sha=vc.read_json(screen_path)
 if not isinstance(screen,dict) or not isinstance(screen.get("cross_cohort_exact"),int): vc.fail(str(screen_path),"not a screen_duplicates.py output")
 if screen["cross_cohort_exact"] and accept_duplicates!=screen_sha: vc.fail(str(screen_path),f"the duplicate screen reports {screen['cross_cohort_exact']} exact cross-cohort duplicates; stop for a parent decision")
 p,a,f,m=(cohorts[x] for x in COHORTS)
 if m.genomes!=prep.matched_selection(f.genomes): vc.fail("matched2k","cohort is not matched_selection of the full cohort (F6)")
 problems,added_counts=cross_check_labels(cohorts)
 if problems: vc.fail("class gates (F3)","; ".join(problems[:10]))
 flagged={d:f"{r} resistant and {s} susceptible isolates among the added genomes (fewer than {MIN_CLASS} of a class)"
          for d,(r,s) in sorted(added_counts.items()) if r<MIN_CLASS or s<MIN_CLASS}
 if set(flagged)!=EXPECTED_SPARSE_ADDED: vc.fail("not interpretable",f"derived {sorted(flagged)}, expected the frozen {sorted(EXPECTED_SPARSE_ADDED)}")
 eligible={name:{vc.ID_TO_DRUG[i] for i in set(c.manifest["reports"])} for name,c in cohorts.items()}
 if eligible["added"]!=EXPECTED_ADDED_ELIGIBLE: vc.fail("added eligibility",f"{sorted(eligible['added'])} is not the frozen 14")
 primary_unavailable={}
 for name,c in cohorts.items():
  for did,report in sorted(c.reports.items()):
   if not report.splits[PRIMARY].available: primary_unavailable.setdefault(report.drug,[]).append({"cohort":name,"reason":report.splits[PRIMARY].reason})
 issues=[]
 for name,c in cohorts.items():
  if c.manifest["source_commit"]!=expected_source: issues.append({"issue":f"{name} source commit differs from expected reviewed source"})
 if len({c.manifest["binary_sha256"] for c in cohorts.values()})!=1: issues.append({"issue":"cohort binaries differ"})
 if len({r["foldgen_sha256"] for r in receipt_rows.values()})!=1: issues.append({"issue":"cohort foldgen binaries differ"})
 configs=[c.manifest["eval_config"] for c in cohorts.values()]
 if any(cf!=configs[0] for cf in configs[1:]): issues.append({"issue":"EvalConfig/seeds differ across cohorts"})
 if f.genomes[:len(p.genomes)]!=p.genomes or f.genomes[len(p.genomes):]!=a.genomes: issues.append({"issue":"full cohort is not exact ordered prefix+added"})
 both=lambda x,y: eligible[x]&eligible[y]  # noqa: E731
 full_prefix,i1=compare_paired(p,f,"prefix",both("prefix","full"),flagged)
 full_matched,i2=compare_paired(m,f,"matched2k",both("matched2k","full"),flagged)
 full_added,i3=compare_full_added(a,f,both("added","full"))
 added_prefix,i4=compare_prefix_added(p,a,both("prefix","added"))
 issues+=i1+i2+i3+i4
 analyses={"full_vs_prefix":full_prefix,"full_vs_matched2k":full_matched,"added_vs_prefix":added_prefix,
  "matched2k_vs_prefix":compare_side_by_side(p,m,"prefix","matched2k",both("prefix","matched2k"),flagged),
  "full_cohort_cv_added_subset":full_added}
 decomposition=[]
 for drug in sorted(both("prefix","full")&eligible["matched2k"]-set(primary_unavailable)):
  did=vc.DRUG_IDS[drug]; meds={name:median_ba(cohorts[name].reports[did],PRIMARY) for name in ("prefix","matched2k","full")}
  decomposition.append({"drug":drug,"before_after_interpretable":drug not in flagged,"median_balanced_accuracy":meds,
   "total":None if None in (meds["full"],meds["prefix"]) else meds["full"]-meds["prefix"],
   "composition":None if None in (meds["matched2k"],meds["prefix"]) else meds["matched2k"]-meds["prefix"],
   "size":None if None in (meds["full"],meds["matched2k"]) else meds["full"]-meds["matched2k"]})
 warnings=[]
 thresholds={name:c.manifest["lineage_threshold"] for name,c in cohorts.items()}
 if len(set(thresholds.values()))!=1: warnings.append({"issue":"lineage thresholds differ across independently rebuilt structures","thresholds":thresholds})
 lineage={"prefix_in_full":vc.compare_lineage_partitions(p,f),"matched2k_in_full":vc.compare_lineage_partitions(m,f)}
 for key,value in lineage.items():
  if not value["same_common_genome_partition"]: warnings.append({"issue":f"lineage partition changed on common genomes ({key}); affected lineage pairs are excluded from matched_fold"})
 random_mismatch={}
 for key in ("full_vs_prefix","full_vs_matched2k"):
  for row in analyses[key]["rows"]:
   if row["scheme"]=="random" and row["comparisons"][0]["status"]=="available":
    pr=row["comparisons"][0]["pairing"]; random_mismatch[f"{key}/{row['drug']}"]={"fold_mismatch":pr["fold_mismatch_count"],"common":pr["common_samples"]}
 if any(v["fold_mismatch"] for v in random_mismatch.values()):
  warnings.append({"issue":"random folds differ between cohorts, as expected for per-cohort stratified shuffles; matched random controls are marked invalid row by row (decision D1), and fold integrity is enforced per cohort by the foldgen rebuild"})
 out={"schema":"grm-mtb-4k-comparison-v3","run_source_commit":expected_source,"foldgen_sha256":expected_foldgen,"validator_sha256":vc.validator_sha256(),"checker_pins":vc.PINS,
  "validation":{"valid":not issues,"blocking_issues":issues,"warnings":warnings},
  "scope":"Internal cross-validation only; no external, prospective, clinical, or causal validation. McNemar p-values are descriptive, Holm-adjusted within learner and split across interpretable drugs.",
  "validator_execution":{"biological_execution":False,"saved_scores_and_metadata_only":True},
  "VALIDATORS_INDEPENDENT_REVIEW":"REQUIRED","WIRING_APPROVAL":"NOT_ASSERTED",
  "cohorts":{name:vc.cohort_public_receipt(c) for name,c in cohorts.items()},"checker_receipts":receipt_rows,
  "duplicate_screen":{"path":str(screen_path),"sha256":screen_sha,"cross_cohort_exact":screen["cross_cohort_exact"],
   "cross_cohort_near":screen.get("cross_cohort_near"),"accepted_by_sha256":bool(screen["cross_cohort_exact"])},
  "primary_split":PRIMARY,"threshold":0.5,"not_interpretable_for_before_after":flagged,"added_class_counts_from_full":{d:list(v) for d,v in sorted(added_counts.items())},
  "eligible":{name:sorted(v) for name,v in eligible.items()},"not_eligible":{name:sorted(set(vc.ALL_DRUGS)-v) for name,v in eligible.items()},
  "primary_split_unavailable":primary_unavailable,"random_fold_accounting":random_mismatch,
  "lineage_comparison":lineage,"analyses":analyses,"size_matched_decomposition":decomposition}
 return out

def parse_pairs(specs,what):
 out={}
 for spec in specs:
  name,_,value=spec.partition("=")
  if name not in COHORTS or name in out or not value: raise vc.ValidationError(f"--{what} {spec}: expected one each of {list(COHORTS)} as NAME=PATH")
  out[name]=value
 return out

def main(argv=None):
 parser=argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
 parser.add_argument("--cohort",action="append",default=[],required=True)
 parser.add_argument("--checker-receipt",action="append",default=[])
 parser.add_argument("--duplicates",required=True)
 parser.add_argument("--accept-reviewed-duplicates")
 parser.add_argument("--expected-source-commit",required=True)
 parser.add_argument("--expected-foldgen-sha256",required=True)
 parser.add_argument("--out",required=True,type=Path)
 args=parser.parse_args(argv)
 try:
  result=compare(parse_pairs(args.cohort,"cohort"),parse_pairs(args.checker_receipt,"checker-receipt"),args.duplicates,
                 args.expected_source_commit,accept_duplicates=args.accept_reviewed_duplicates,expected_foldgen=args.expected_foldgen_sha256)
 except (vc.ValidationError,OSError,KeyError,TypeError,ValueError,ArithmeticError) as exc:
  print(f"compare_cohorts: STOP: {exc}",file=sys.stderr); return 2
 with args.out.open("x") as f:
  f.write(json.dumps(result,indent=1,sort_keys=True,allow_nan=False)+"\n")
 print(json.dumps({"out":str(args.out),"valid":result["validation"]["valid"],"blocking_issues":len(result["validation"]["blocking_issues"]),
                   "warnings":len(result["validation"]["warnings"]),
                   "primary_split_unavailable":sorted(result["primary_split_unavailable"]),
                   "reports":{n:c["eligible_reports"] for n,c in result["cohorts"].items()}},sort_keys=True))
 if not result["validation"]["valid"]: return 1
 if result["primary_split_unavailable"]:
  print("compare_cohorts: stop for review before any before/after claim",file=sys.stderr); return 3
 return 0
if __name__=="__main__": sys.exit(main())
