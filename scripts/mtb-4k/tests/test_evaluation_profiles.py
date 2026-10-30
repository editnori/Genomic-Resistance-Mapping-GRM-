#!/usr/bin/env python3
"""Evaluation profiles of report_contract.py: r8 and rN (revision N on the frozen cohort rules).

usage: FOLDGEN=<foldgen binary> [MTB4K_TEST_PROFILE=r8|rN] python3 -I scripts/mtb-4k/tests/test_evaluation_profiles.py

The suite runs under MTB4K_TEST_PROFILE (profile_boot.py; r8 by default). These tests hold under
either: "the tree" below is this checkout's own checker blobs, and OTHER is a revision other than
the suite's (9 under r8, 10 under a revision-9 rN), so the revision a profile binds is always
seen to move.

- --profile is read from the process's arguments at import, as argparse reads it;
- rN refuses, and leaves the active profile as it was, while revision_n.py's values hold a
  placeholder (also with MTB4K_TEST_PROFILE=rN set, in a subprocess), when a value is malformed,
  when revision N drops or renames a field, when its checker judges convergence statuses unlike
  690815a, and when a pin differs from the tree's blob;
- rN takes exactly the revision-9 facts of 2026-10-29 (b2fd433): they activate on a tree that
  holds b2fd433's blobs and stop at the pin check on any other;
- with the tree's pins, rN activates at OTHER, and that revision reaches every check:
  validate_manifest, the runner's checker and audit calls, the postflight, the summary coverage
  and the comparator's receipt check. A root restamped to OTHER stops at the tree's checker, which
  accepts revisions up to the suite's only.
"""
import contextlib,copy,io,json,os,re,shutil,subprocess,sys,tempfile,unittest
from pathlib import Path
LANE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(LANE)); sys.path.insert(0,str(LANE/"tests"))
import profile_boot as pb  # noqa: E402  (MTB4K_TEST_PROFILE, before report_contract loads)
import report_contract as vc
import revision_n
import run_revision8_checks as wrapper
import compare_cohorts as cc
import fixtures as fx

BASE=pb.PROFILE
# The revision-9 facts relayed on 2026-10-29 (r9/lasso-exchange-fix, b2fd433), in revision_n.py's form.
R9={"REVISION_N":9,"RN_CHECKER_COMMIT":"b2fd4337ce88be64767a3a552e87df21e64803c4",
 "RN_PIN_CHECK_RESULTS":"b7e5a8a067179bc4578d6487aebe2baad51f4fe419a44aa2e3a942ae3d10a347",
 "RN_PIN_METRIC_AUDIT":"c1af2ead6fbedb869a067509efd2a431a218a9ddc520fdfb3625d2054e4af98b",
 "RN_PIN_SUMMARIZE":"fd5fdfab1bedb8118c95eaadef19452fc89786f371a5d1757bd784844eb3295f",
 "RN_PIN_RESULT_SCOPE":"3bcb2db34b23accd9186a595e6e0002cbd3bc292ab548442d8a53ba94b65ae2a",
 "RN_FIELD_CHANGES":{},"RN_CONVERGENCE_STATUS_RULE":"as-690815a"}
TREE_PINS={name:vc.file_sha256(vc.SCRIPTS/rel) for rel,name in revision_n.PIN_NAMES.items()}
OTHER=10 if vc.REVISION==9 else 9
IN_TREE={**R9,**TREE_PINS,"REVISION_N":OTHER}
UNFILLED={n:f"@@{n}@@" for n in revision_n.TABLE}

def quiet(fn):
 with contextlib.redirect_stdout(io.StringIO()),contextlib.redirect_stderr(io.StringIO()) as err:
  return fn(),err

@contextlib.contextmanager
def profile(name,values=None):
 vc.activate(name,values)
 try: yield
 finally: vc.activate(BASE)

def state(): return (vc.PROFILE,vc.REVISION,dict(vc.PINS),vc.RULES,vc.upstream_checker,set(vc.FIT_KEYS))

def restamp(root,revision):
 """A fixture root whose reports, index and run manifest carry another evaluation revision."""
 for p in (root/"data/species"/vc.SPECIES/"results").glob("*.json"):
  r=json.loads(p.read_text()); r["evaluation"]["revision"]=revision; fx.write_json(p,r)
 fx.refresh_root(root,"prefix")
 index=root/"data/report/index.json"; i=json.loads(index.read_text()); i["evaluation_revision"]=revision; fx.write_json(index,i)
 m=json.loads((root/"validation-manifest.json").read_text()); m["evaluation_revision"]=revision
 for e in m["inputs"]: e["sha256"]=fx.sha(root/e["path"])
 fx.write_json(root/"validation-manifest.json",m)

def unfilled_lane(directory):
 """A copy of scripts/mtb-4k (no tests but profile_boot.py) and the four pinned checker files, with
 every revision_n.py value reset to its placeholder, whatever this tree holds."""
 scripts=Path(directory)/"scripts"
 shutil.copytree(LANE,scripts/"mtb-4k",ignore=shutil.ignore_patterns("__pycache__"))
 for rel in revision_n.PIN_NAMES:
  (scripts/rel).parent.mkdir(parents=True,exist_ok=True); shutil.copyfile(vc.SCRIPTS/rel,scripts/rel)
 path=scripts/"mtb-4k/revision_n.py"; text=path.read_text()
 start=text.index("\nVALUES = {\n")+1; end=text.index("\n}\n",start)+3
 path.write_text(text[:start]+"VALUES = {\n"+"".join(f'    "{n}": "@@{n}@@",\n' for n in revision_n.TABLE)+"}\n"+text[end:])
 return scripts/"mtb-4k"

class Selection(unittest.TestCase):
 def test_profile_is_read_from_argv_as_argparse_reads_it(self):
  cases=[([],"r8"),(["out.json","--root","x"],"r8"),(["--profile","rN"],"rN"),(["--profile=rN"],"rN"),
   (["--profile","rN","--profile","r8"],"r8"),(["--","--profile","rN"],"r8"),(["--profile"],"r8")]
  for argv,want in cases:
   with self.subTest(argv=argv): self.assertEqual(vc.profile_from_argv(argv),want)
 def test_this_process_runs_the_suites_profile(self):
  self.assertEqual(vc.PROFILE,BASE)
  self.assertEqual(vc.PINS,{rel:TREE_PINS[name] for rel,name in revision_n.PIN_NAMES.items()})
  self.assertIs(vc.RULES,vc.R8_RULES)
  if BASE=="r8": self.assertEqual((vc.REVISION,vc.PINS),(8,vc.R8_PINS))
 def test_an_unknown_profile_is_refused(self):
  before=state()
  with self.assertRaises(vc.ValidationError) as stop: vc.activate("r9")
  self.assertIn("unknown evaluation profile",str(stop.exception)); self.assertEqual(state(),before)

class Refusals(unittest.TestCase):
 def refused(self,values,*needles):
  before=state()
  with self.assertRaises(vc.ValidationError) as stop: vc.activate("rN",values)
  for needle in needles: self.assertIn(needle,str(stop.exception))
  self.assertEqual(state(),before,"a refused activation changed the active profile")
  return str(stop.exception)
 def test_unfilled_values_refuse_and_name_every_one(self):
  message=self.refused(UNFILLED,"unfilled placeholder")
  for name in revision_n.CONTRACT_NAMES: self.assertIn(f"@@{name}@@",message)
 def test_an_unfilled_rn_refuses_with_the_switch_set(self):
  with tempfile.TemporaryDirectory() as tmp:
   lane=unfilled_lane(tmp)
   p=subprocess.run([sys.executable,"-I","-B","-c","import sys; sys.path.insert(0,sys.argv[1]); import profile_boot",str(lane/"tests")],
                    cwd=tmp,env={**os.environ,"MTB4K_TEST_PROFILE":"rN"},capture_output=True,text=True)
   self.assertNotEqual(p.returncode,0,p.stdout)
   self.assertIn("profile rN: unfilled placeholder(s): @@REVISION_N@@",p.stderr)
 def test_one_placeholder_left_is_enough_to_refuse(self):
  for name in revision_n.CONTRACT_NAMES:
   with self.subTest(name=name): self.refused({**IN_TREE,name:f"@@{name}@@"},f"@@{name}@@")
 def test_malformed_values_refuse(self):
  for name,bad in [("REVISION_N",8),("REVISION_N","9"),("REVISION_N",True),("RN_PIN_CHECK_RESULTS","B7E5"+"0"*60),
                   ("RN_PIN_METRIC_AUDIT","0"*63),("RN_FIELD_CHANGES",{"NEW_KEYS":["x"]}),("RN_FIELD_CHANGES",[]),
                   ("RN_CONVERGENCE_STATUS_RULE","")]:
   with self.subTest(name=name,bad=bad): self.refused({**IN_TREE,name:bad},name)
 def test_a_dropped_or_renamed_field_needs_a_contract_change(self):
  renamed=sorted(vc.BASE_KEYS["CONVERGENCE_KEYS"]-{"sweep_cap_hits"}|{"cap_hits"})
  self.refused({**IN_TREE,"RN_FIELD_CHANGES":{"CONVERGENCE_KEYS":renamed}},"drops or renames ['sweep_cap_hits']")
 def test_other_convergence_semantics_need_a_contract_change(self):
  self.refused({**IN_TREE,"RN_CONVERGENCE_STATUS_RULE":"every-step-capped-is-unconverged"},"validate_convergence needs a reviewed change")
 def test_a_pin_that_is_not_the_trees_blob_refuses(self):
  self.refused({**IN_TREE,"RN_PIN_CHECK_RESULTS":"0"*64},"check-results.py","differs from the reviewed pin 0000")
 def test_the_revision_9_facts_are_taken_and_activate_only_on_b2fd433s_blobs(self):
  self.assertEqual(revision_n.require(list(revision_n.TABLE),R9),R9)
  if all(TREE_PINS[n]==R9[n] for n in TREE_PINS):
   with profile("rN",R9): self.assertEqual((vc.PROFILE,vc.REVISION),("rN",9))
  else:
   self.refused(R9,"differs from the reviewed pin")

class Active(unittest.TestCase):
 def test_rn_binds_revision_pins_rules_and_modules(self):
  with profile("rN",IN_TREE):
   self.assertEqual((vc.PROFILE,vc.REVISION),("rN",OTHER))
   self.assertEqual(vc.PINS,{rel:TREE_PINS[name] for rel,name in revision_n.PIN_NAMES.items()}); self.assertIs(vc.RULES,vc.R8_RULES)
   self.assertIsNotNone(vc.upstream_checker.check_lasso_fit)
  self.assertEqual(vc.PROFILE,BASE)
 def test_an_added_field_extends_its_key_set_and_the_base_profile_restores_it(self):
  fit_keys=sorted(vc.BASE_KEYS["FIT_KEYS"]|{"exchange_steps"})
  with profile("rN",{**IN_TREE,"RN_FIELD_CHANGES":{"FIT_KEYS":fit_keys}}):
   self.assertIn("exchange_steps",vc.FIT_KEYS)
   fit=fx.fit(); fit["exchange_steps"]=0
   vc.expect_keys(fit,vc.FIT_KEYS,"synthetic")
  self.assertNotIn("exchange_steps",vc.FIT_KEYS)
  with self.assertRaises(vc.ValidationError): vc.expect_keys(fit,vc.FIT_KEYS,"synthetic")
 def test_the_revision_reaches_the_contract(self):
  with tempfile.TemporaryDirectory() as tmp:
   rules=fx.make_rules(); root=fx.build_cohort(Path(tmp)/"prefix","prefix",rules)
   vc.validate_manifest(root/"validation-manifest.json",rules)
   with profile("rN",IN_TREE):
    with self.assertRaises(vc.ValidationError) as stop: vc.validate_manifest(root/"validation-manifest.json",rules)
    self.assertIn("wrong species/revision/schema",str(stop.exception))
    restamp(root,OTHER)
    cohort=vc.validate_manifest(root/"validation-manifest.json",rules)
    self.assertEqual({r.receipt["evaluation_revision"] for r in cohort.reports.values()},{OTHER})
   with self.assertRaises(vc.ValidationError): vc.validate_manifest(root/"validation-manifest.json",rules)
 def test_the_revision_reaches_the_runner_and_the_trees_checker_refuses_it(self):
  with tempfile.TemporaryDirectory() as tmp:
   rules=fx.make_rules(); root=fx.build_cohort(Path(tmp)/"prefix","prefix",rules)
   args=[str(root/"audit.json"),"--root",str(root),"--run-manifest",str(root/"validation-manifest.json"),
    "--expected-source-commit",fx.SOURCE,"--binary",str(root/"synthetic-binary"),"--backup",str(root/"backup"),
    "--log",str(root/"pipeline.log"),"--launch","2000-01-01T00:00:00Z","--foldgen",fx.foldgen_path(),
    "--expected-foldgen-sha256",fx.foldgen_sha256(),"--checker-out",str(root/"checker.json"),"--summary-out",str(root/"summary.json"),
    "--receipt",str(root/"receipt.json")]
   with profile("rN",IN_TREE):
    restamp(root,OTHER)
    code,err=quiet(lambda:wrapper.main(args+["--profile","r8"],rules=rules))
    self.assertEqual(code,2); self.assertIn("--profile r8 is not the active profile rN",err.getvalue())
    code,err=quiet(lambda:wrapper.main(args+["--profile","rN"],rules=rules))
    self.assertEqual(code,2); self.assertIn(f"revision-{OTHER} checker exited",err.getvalue())
    self.assertIn("('revision',",err.getvalue())
    self.assertFalse((root/"audit.json").exists() or (root/"receipt.json").exists())
 def test_the_revision_reaches_the_postflight_summary_and_comparator(self):
  with tempfile.TemporaryDirectory() as tmp:
   rules=fx.make_rules(); root=fx.build_cohort(Path(tmp)/"prefix","prefix",rules)
   cohort=vc.validate_manifest(root/"validation-manifest.json",rules)
   bodies=fx.passing_outputs(root,cohort.manifest)
   paths={}
   for key,body in bodies.items(): paths[key]=root/f"{key}.json"; fx.write_json(paths[key],body)
   self.assertEqual(wrapper.postflight(paths["audit"],cohort,paths["checker"]),[])
   self.assertEqual(wrapper.summary_errors(paths["summary"],cohort),[])
   receipt=fx.write_receipt(root)
   with profile("rN",IN_TREE):
    errors=wrapper.postflight(paths["audit"],cohort,paths["checker"])
    for needle in (f"auditor revision is not {OTHER}","historical/wrong revisions in audited set","manifest schema/revision differs",
                   "checker accepted mixed revisions"): self.assertIn(needle,errors)
    self.assertTrue(any("wrong revision/schema/validity" in e for e in errors))
    self.assertEqual(wrapper.summary_errors(paths["summary"],cohort),["summarizer coverage/issues differ"])
    with self.assertRaises(vc.ValidationError) as stop: cc.load_receipt(receipt,cohort,fx.SOURCE,fx.foldgen_sha256())
    self.assertIn(f"receipt is not for evaluation revision {OTHER}",str(stop.exception))
    # The same outputs restamped to OTHER pass that revision's postflight.
    other=copy.deepcopy(bodies); s=other["audit"]["summary"]; key=str(OTHER)
    s.update(revision=OTHER,by_revision={key:s["files"]},manifest_evaluation_revision=OTHER)
    for rec in other["audit"]["records"]: rec.update(revision=OTHER,published_revision=OTHER)
    for rec in other["checker"]["reports"]: rec["revision"]=OTHER
    other["summary"]["coverage"]={key:other["summary"]["coverage"][str(cohort.manifest["evaluation_revision"])]}
    for key_,body in other.items(): fx.write_json(paths[key_],body)
    self.assertEqual(wrapper.postflight(paths["audit"],cohort,paths["checker"]),[])
    self.assertEqual(wrapper.summary_errors(paths["summary"],cohort),[])

if __name__=="__main__":
 result=unittest.main(exit=False,verbosity=2).result
 sys.exit(0 if result.wasSuccessful() else 1)
