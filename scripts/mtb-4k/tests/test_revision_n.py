#!/usr/bin/env python3
"""revision_n.py, whether the committed file is unfilled (this branch) or filled (a fill commit, such
as the revision-9 fallback merge's, which changes VALUES and no test).

- The committed file is wholly unfilled or wholly filled. Unfilled, every CLI call refuses; filled,
  `status` and `verify-checkout` of this checkout exit 0, and every value has its kind and format.
- Unfilled values refuse, in require and in the CLI, on a copy whose VALUES this file resets.
- Values are checked by kind; a copy whose VALUES this file sets to the revision-9 facts runs.
- verify-checkout on a scratch git repository.

usage: python3 -I scripts/mtb-4k/tests/test_revision_n.py
"""
import hashlib,json,os,re,shutil,subprocess,sys,tempfile,unittest
from pathlib import Path
LANE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(LANE))
import revision_n as rn
REPO=LANE.parents[1]  # the checkout this file belongs to

R9={"REVISION_N":9,"RN_CHECKER_COMMIT":"b2fd4337ce88be64767a3a552e87df21e64803c4",
 "RN_PIN_CHECK_RESULTS":"b7e5a8a067179bc4578d6487aebe2baad51f4fe419a44aa2e3a942ae3d10a347",
 "RN_PIN_METRIC_AUDIT":"c1af2ead6fbedb869a067509efd2a431a218a9ddc520fdfb3625d2054e4af98b",
 "RN_PIN_SUMMARIZE":"fd5fdfab1bedb8118c95eaadef19452fc89786f371a5d1757bd784844eb3295f",
 "RN_PIN_RESULT_SCOPE":"3bcb2db34b23accd9186a595e6e0002cbd3bc292ab548442d8a53ba94b65ae2a",
 "RN_FIELD_CHANGES":{},"RN_CONVERGENCE_STATUS_RULE":"as-690815a"}

# Revision 10 as relayed on 2026-10-30 (r10/integration, c90ebaf).
R10={**R9,"REVISION_N":10,"RN_CHECKER_COMMIT":"c90ebaf48998dc244806ac435c3706bfd87c9a9b",
 "RN_PIN_CHECK_RESULTS":"28b65bdbffb40ce6d4c010d6b7509855e9fa1b7c78d337889334492a2bcfa9df",
 "RN_PIN_METRIC_AUDIT":"7b76746b896a7b4dab96ac360b4250e2bcdefa3bcb2dc5615d0414d8edfd1f30",
 "RN_FIELD_CHANGES":{"OPTIONAL_TOP_KEYS":["kset_duplicates"]}}

def run(script,*args,cwd):
 return subprocess.run([sys.executable,"-I","-B",str(script),*map(str,args)],cwd=cwd,capture_output=True,text=True)

def with_values(directory,values):
 """A copy of revision_n.py whose whole VALUES block is rewritten to `values`, whatever the committed
 block holds (placeholders, filled values, any layout)."""
 text=(LANE/"revision_n.py").read_text()
 start=text.index("\nVALUES = {\n")+1; end=text.index("\n}\n",start)+3
 block="VALUES = {\n"+"".join(f"    {name!r}: {values[name]!r},\n" for name in rn.TABLE)+"}\n"
 path=Path(directory)/"revision_n.py"; path.write_text(text[:start]+block+text[end:]); return path

def unfilled_copy(directory): return with_values(directory,{n:f"@@{n}@@" for n in rn.TABLE})

def filled_copy(directory,values): return with_values(directory,values)

class Placeholders(unittest.TestCase):
 def test_the_committed_file_is_wholly_unfilled_and_refuses_or_wholly_filled_and_verifies(self):
  self.assertEqual(set(rn.VALUES),set(rn.TABLE))
  missing=rn.unfilled()
  self.assertIn(missing,[list(rn.TABLE),[]],"revision_n.py is partly filled")
  with tempfile.TemporaryDirectory() as cwd:
   status=run(LANE/"revision_n.py","status",cwd=cwd)
   verify=run(LANE/"revision_n.py","verify-checkout",REPO,cwd=cwd)
   if missing:  # this branch: everything refuses
    self.assertEqual(rn.VALUES,{n:f"@@{n}@@" for n in rn.TABLE})
    self.assertEqual((status.returncode,verify.returncode),(2,2),status.stderr+verify.stderr)
    for name in rn.TABLE:
     with self.subTest(name=name): self.assertEqual(run(LANE/"revision_n.py","get",name,cwd=cwd).returncode,2)
    with self.assertRaises(rn.Unfilled): rn.require(rn.CONTRACT_NAMES)
   else:  # a fill commit: everything is well formed and this checkout verifies
    self.assertEqual(status.returncode,0,status.stderr)
    values=rn.require(list(rn.TABLE))
    self.assertTrue(isinstance(values["REVISION_N"],int) and not isinstance(values["REVISION_N"],bool) and values["REVISION_N"]>=9)
    self.assertRegex(values["RN_CHECKER_COMMIT"],r"\A[0-9a-f]{40}\Z")
    for name in rn.PIN_NAMES.values(): self.assertRegex(values[name],r"\A[0-9a-f]{64}\Z")
    self.assertIsInstance(values["RN_FIELD_CHANGES"],dict)
    self.assertEqual(verify.returncode,0,verify.stdout+verify.stderr)
  # The contract reads every value but the checker commit, which verify-checkout reads.
  self.assertEqual(set(rn.CONTRACT_NAMES)|{"RN_CHECKER_COMMIT"},set(rn.TABLE))
 def test_require_names_every_unfilled_placeholder(self):
  with self.assertRaises(rn.Unfilled) as stop: rn.require(rn.CONTRACT_NAMES,{n:f"@@{n}@@" for n in rn.TABLE})
  self.assertEqual(stop.exception.names,rn.CONTRACT_NAMES)
  with self.assertRaises(rn.Unfilled) as stop: rn.require(rn.CONTRACT_NAMES,{**R9,"RN_PIN_SUMMARIZE":"@@RN_PIN_SUMMARIZE@@"})
  self.assertEqual(stop.exception.names,["RN_PIN_SUMMARIZE"])
  with self.assertRaises(rn.Unfilled): rn.require(["REVISION_N"],{})
  self.assertEqual(rn.unfilled({n:f"@@{n}@@" for n in rn.TABLE}),list(rn.TABLE))
 def test_the_cli_refuses_while_unfilled(self):
  with tempfile.TemporaryDirectory() as cwd:
   script=unfilled_copy(cwd)
   p=run(script,"status",cwd=cwd)
   self.assertEqual(p.returncode,2,p.stderr); self.assertIn("revision_n: STOP: unfilled placeholder(s)",p.stderr)
   self.assertEqual(json.loads(p.stdout),{n:False for n in rn.TABLE})
   for name in rn.TABLE:
    with self.subTest(name=name):
     p=run(script,"get",name,cwd=cwd)
     self.assertEqual((p.returncode,p.stdout),(2,"")); self.assertIn(f"@@{name}@@",p.stderr)
   p=run(script,"verify-checkout",cwd,cwd=cwd)
   self.assertEqual(p.returncode,2); self.assertIn("@@RN_CHECKER_COMMIT@@",p.stderr)

class Filled(unittest.TestCase):
 def test_the_revision_9_facts_are_well_formed(self):
  self.assertEqual(rn.require(list(rn.TABLE),R9),R9)
 def test_the_revision_10_facts_are_well_formed(self):
  self.assertEqual(rn.require(list(rn.TABLE),R10),R10)
  with self.assertRaises(rn.Malformed): rn.require(["RN_FIELD_CHANGES"],{**R10,"RN_FIELD_CHANGES":{"OPTIONAL_TOP_KEYS":[]}})
 def test_values_are_checked_by_kind(self):
  for name,bad in [("REVISION_N",8),("REVISION_N",9.0),("REVISION_N",True),("REVISION_N","9"),("RN_CHECKER_COMMIT","b2fd433"),
                   ("RN_CHECKER_COMMIT","B2FD4337CE88BE64767A3A552E87DF21E64803C4"),("RN_PIN_CHECK_RESULTS","0"*63),
                   ("RN_PIN_RESULT_SCOPE","g"*64),("RN_FIELD_CHANGES",None),("RN_FIELD_CHANGES",{"FIT_KEYS":[]}),
                   ("RN_FIELD_CHANGES",{"OTHER_KEYS":["x"]}),("RN_CONVERGENCE_STATUS_RULE",1)]:
   with self.subTest(name=name,bad=bad):
    with self.assertRaises(rn.Malformed): rn.require([name],{**R9,name:bad})
 def test_a_copy_filled_with_the_revision_9_facts_runs(self):
  with tempfile.TemporaryDirectory() as tmp:
   script=filled_copy(tmp,R9)
   p=run(script,"status",cwd=tmp)
   self.assertEqual(p.returncode,0,p.stderr); self.assertEqual(json.loads(p.stdout),{n:True for n in rn.TABLE})
   for name,want in [("REVISION_N","9"),("RN_CHECKER_COMMIT",R9["RN_CHECKER_COMMIT"]),("RN_FIELD_CHANGES","{}"),
                     ("RN_CONVERGENCE_STATUS_RULE","as-690815a"),("RN_PIN_SUMMARIZE",R9["RN_PIN_SUMMARIZE"])]:
    p=run(script,"get",name,cwd=tmp); self.assertEqual((p.returncode,p.stdout.strip()),(0,want),p.stderr)

class VerifyCheckout(unittest.TestCase):
 def git(self,repo,*args):
  return subprocess.run(["git","-C",str(repo),"-c","user.name=t","-c","user.email=t@example.invalid","-c","commit.gpgsign=false",*args],
                        capture_output=True,text=True,check=True).stdout.strip()
 def test_a_clean_descendant_with_the_pinned_blobs_passes_and_each_fault_is_named(self):
  with tempfile.TemporaryDirectory() as tmp:
   repo=Path(tmp)/"ws"; repo.mkdir(); self.git(repo,"init","-q")
   for rel in rn.PIN_NAMES:
    (repo/"scripts"/rel).parent.mkdir(parents=True,exist_ok=True); (repo/"scripts"/rel).write_text(f"# {rel}\n")
   self.git(repo,"add","-A"); self.git(repo,"commit","-q","-m","checker")
   checker=self.git(repo,"rev-parse","HEAD")
   (repo/"other.txt").write_text("merge\n"); self.git(repo,"add","-A"); self.git(repo,"commit","-q","-m","merge")
   values={**R9,"RN_CHECKER_COMMIT":checker,**{n:hashlib.sha256((repo/"scripts"/rel).read_bytes()).hexdigest() for rel,n in rn.PIN_NAMES.items()}}
   self.assertEqual(rn.verify_checkout(repo,values),[])
   (repo/"untracked.txt").write_text("x\n")
   self.assertEqual(rn.verify_checkout(repo,values),[f"{repo} is not a clean git checkout"])
   (repo/"untracked.txt").unlink()
   (repo/"scripts/result_scope.py").write_text("# edited\n")
   problems=rn.verify_checkout(repo,values)
   self.assertEqual(len(problems),2); self.assertTrue(any("result_scope.py does not hold the pinned SHA-256" in p for p in problems))
   self.git(repo,"checkout","-q","--","scripts/result_scope.py")
   branch=self.git(repo,"rev-parse","--abbrev-ref","HEAD")
   self.git(repo,"checkout","-q","--orphan","other"); self.git(repo,"commit","-q","-m","unrelated")
   problems=rn.verify_checkout(repo,values)
   self.assertIn(f"HEAD does not descend from RN_CHECKER_COMMIT {checker}",problems)
   self.git(repo,"checkout","-q",branch)
   self.assertEqual(rn.verify_checkout(repo,values),[])
   problems=rn.verify_checkout(repo,{**values,"RN_PIN_METRIC_AUDIT":"0"*64})
   self.assertEqual(len(problems),2); self.assertTrue(all("RN_PIN_METRIC_AUDIT" in p for p in problems))

if __name__=="__main__":
 result=unittest.main(exit=False,verbosity=2).result
 sys.exit(0 if result.wasSuccessful() else 1)
