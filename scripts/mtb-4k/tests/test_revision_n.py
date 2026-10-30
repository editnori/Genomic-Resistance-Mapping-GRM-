#!/usr/bin/env python3
"""revision_n.py: every value is an unfilled placeholder at this commit and refuses; filled values
are checked by kind; a copy filled with the revision-9 facts runs; verify-checkout.

usage: python3 -I scripts/mtb-4k/tests/test_revision_n.py
"""
import hashlib,json,os,re,shutil,subprocess,sys,tempfile,unittest
from pathlib import Path
LANE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(LANE))
import revision_n as rn

R9={"REVISION_N":9,"RN_CHECKER_COMMIT":"b2fd4337ce88be64767a3a552e87df21e64803c4",
 "RN_PIN_CHECK_RESULTS":"b7e5a8a067179bc4578d6487aebe2baad51f4fe419a44aa2e3a942ae3d10a347",
 "RN_PIN_METRIC_AUDIT":"c1af2ead6fbedb869a067509efd2a431a218a9ddc520fdfb3625d2054e4af98b",
 "RN_PIN_SUMMARIZE":"fd5fdfab1bedb8118c95eaadef19452fc89786f371a5d1757bd784844eb3295f",
 "RN_PIN_RESULT_SCOPE":"3bcb2db34b23accd9186a595e6e0002cbd3bc292ab548442d8a53ba94b65ae2a",
 "RN_FIELD_CHANGES":{},"RN_CONVERGENCE_STATUS_RULE":"as-690815a"}

def run(script,*args,cwd):
 return subprocess.run([sys.executable,"-I","-B",str(script),*map(str,args)],cwd=cwd,capture_output=True,text=True)

def filled_copy(directory,values):
 """revision_n.py with each placeholder string replaced by repr(value), as a filling commit would."""
 text=(LANE/"revision_n.py").read_text()
 for name,value in values.items():
  needle=f'"{name}": "@@{name}@@",'
  assert text.count(needle)==1,name
  text=text.replace(needle,f'"{name}": {value!r},')
 path=Path(directory)/"revision_n.py"; path.write_text(text); return path

class Placeholders(unittest.TestCase):
 def test_every_committed_value_is_its_own_placeholder(self):
  # The commit that fills a value changes this assertion with it.
  self.assertEqual(set(rn.VALUES),set(rn.TABLE))
  self.assertEqual(rn.VALUES,{n:f"@@{n}@@" for n in rn.TABLE})
  self.assertEqual(rn.unfilled(),list(rn.TABLE))
  # The contract reads every value but the checker commit, which verify-checkout reads.
  self.assertEqual(set(rn.CONTRACT_NAMES)|{"RN_CHECKER_COMMIT"},set(rn.TABLE))
 def test_require_names_every_unfilled_placeholder(self):
  with self.assertRaises(rn.Unfilled) as stop: rn.require(rn.CONTRACT_NAMES)
  self.assertEqual(stop.exception.names,rn.CONTRACT_NAMES)
  with self.assertRaises(rn.Unfilled) as stop: rn.require(rn.CONTRACT_NAMES,{**R9,"RN_PIN_SUMMARIZE":"@@RN_PIN_SUMMARIZE@@"})
  self.assertEqual(stop.exception.names,["RN_PIN_SUMMARIZE"])
  with self.assertRaises(rn.Unfilled): rn.require(["REVISION_N"],{})
 def test_the_cli_refuses_while_unfilled(self):
  with tempfile.TemporaryDirectory() as cwd:
   p=run(LANE/"revision_n.py","status",cwd=cwd)
   self.assertEqual(p.returncode,2,p.stderr); self.assertIn("revision_n: STOP: unfilled placeholder(s)",p.stderr)
   self.assertEqual(json.loads(p.stdout),{n:False for n in rn.TABLE})
   for name in rn.TABLE:
    with self.subTest(name=name):
     p=run(LANE/"revision_n.py","get",name,cwd=cwd)
     self.assertEqual((p.returncode,p.stdout),(2,"")); self.assertIn(f"@@{name}@@",p.stderr)
   p=run(LANE/"revision_n.py","verify-checkout",cwd,cwd=cwd)
   self.assertEqual(p.returncode,2); self.assertIn("@@RN_CHECKER_COMMIT@@",p.stderr)

class Filled(unittest.TestCase):
 def test_the_revision_9_facts_are_well_formed(self):
  self.assertEqual(rn.require(list(rn.TABLE),R9),R9)
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
