#!/usr/bin/env python3
"""diff_reference_prefix.py: the rerun's prefix against revision N's MTB 2k candidate, every field
but wall-clock time, lasso included.

usage: FOLDGEN=<foldgen binary> [MTB4K_TEST_PROFILE=r8|rN] python3 -I scripts/mtb-4k/tests/test_diff_reference_prefix.py
"""
import contextlib,io,json,shutil,subprocess,sys,tempfile,unittest
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent)); sys.path.insert(0,str(HERE))
import profile_boot as pb  # noqa: E402  (MTB4K_TEST_PROFILE, before report_contract loads)
import fixtures as fx  # noqa: E402
import diff_reference_prefix as diff  # noqa: E402

SLUG="mycobacterium_tuberculosis"
REV=pb.vc.REVISION if pb.vc.REVISION>=9 else 9  # revision N is at least 9; the fixtures' reports are restamped to it

class Diff(unittest.TestCase):
 @classmethod
 def setUpClass(cls):
  cls.tmp=tempfile.TemporaryDirectory()
  rules=fx.make_rules({"prefix":("amikacin","isoniazid")})
  cls.prefix=fx.build_cohort(Path(cls.tmp.name)/"prefix","prefix",rules)/"data"
  for p in (cls.prefix/"species"/SLUG/"results").glob("*.json"):
   r=json.loads(p.read_text()); r["evaluation"]["revision"]=REV; p.write_text(json.dumps(r))
  cls.saved=dict(diff.rev_n.VALUES)
  diff.rev_n.VALUES["REVISION_N"]=REV
 @classmethod
 def tearDownClass(cls):
  diff.rev_n.VALUES.clear(); diff.rev_n.VALUES.update(cls.saved); cls.tmp.cleanup()

 def reference(self,label,edit=None,genomes=None):
  root=Path(self.tmp.name)/label/"data"; sp=root/"species"/SLUG
  (sp/"results").mkdir(parents=True)
  (sp/"genomes.txt").write_bytes(genomes if genomes is not None else (self.prefix/"species"/SLUG/"genomes.txt").read_bytes())
  for path in sorted((self.prefix/"species"/SLUG/"results").glob("*.json")):
   r=json.loads(path.read_text())
   r["seconds"]=r["seconds"]+50.0
   for split in r["splits"]:
    for row in split["learners"]: row["seconds"]=row["seconds"]+7.0
   if edit and path.stem.endswith("amikacin"): edit(r)
   (sp/"results"/path.name).write_text(json.dumps(r))
  return root

 def run_diff(self,ref):
  out=ref.parent/"diff.json"
  with contextlib.redirect_stdout(io.StringIO()),contextlib.redirect_stderr(io.StringIO()) as err:
   code=diff.main(["--prefix",str(self.prefix),"--reference",str(ref),"--out",str(out)])
  return code,(json.loads(out.read_text()) if out.exists() else None),err.getvalue()

 def test_only_wall_clock_seconds_may_differ(self):
  code,out,err=self.run_diff(self.reference("same"))
  self.assertEqual(code,0,err); self.assertEqual((out["datasets"],out["identical"],out["evaluation_revision"]),(2,2,REV))
  self.assertEqual(out["rows"][0]["not_compared"],["seconds","splits[*].learners[*].seconds"])

 def test_every_other_field_lasso_included_is_a_difference(self):
  lasso=lambda r:next(row for row in r["splits"][0]["learners"] if row["learner"]=="lasso")  # noqa: E731
  for label,edit,where in [
    ("lasso held-out score",lambda r:lasso(r)["scores"].__setitem__(0,0.123456),"splits[0].learners"),
    ("lasso fit status",lambda r:r["models"]["lasso_fit"].update(path_points=r["models"]["lasso_fit"]["path_points"]-1),"models.lasso_fit.path_points"),
    ("final lasso model",lambda r:r["models"].update(lasso={"intercept":0.0,"coef":[],"lambda":1.0}),"models.lasso"),
    ("feature cards",lambda r:r.update(features=[{"feature":1}]),"features"),
    ("cart score",lambda r:r["splits"][0]["learners"][1]["scores"].__setitem__(0,0.5),"splits[0].learners[1].scores"),
    ("fingerprint",lambda r:r.update(input_fingerprint="0"*16),"input_fingerprint"),
    ("evaluation reasons",lambda r:r["evaluation"].update(reasons=["x"]),"evaluation.reasons"),
    ("a field on one side only",lambda r:r["models"].update(extra=1),"models.extra"),
    ("a model's seconds",lambda r:r["models"].update(seconds=1.0),"models.seconds")]:
   with self.subTest(label):
    code,out,err=self.run_diff(self.reference("d-"+label.replace(" ","-").replace("'",""),edit))
    self.assertEqual(code,1,(label,err))
    row=next(r for r in out["rows"] if r["id"].endswith("amikacin"))
    self.assertTrue(any(d.startswith(where) for d in row["differences"]),(label,row["differences"]))
    self.assertTrue(next(r for r in out["rows"] if r["id"].endswith("isoniazid"))["identical"])

 def test_another_genome_list_revision_or_dataset_set_stops(self):
  genomes=(self.prefix/"species"/SLUG/"genomes.txt").read_text().splitlines()
  ref=self.reference("other-genomes",genomes=("\n".join(genomes[1:]+genomes[:1])+"\n").encode())
  code,_,err=self.run_diff(ref); self.assertEqual(code,2); self.assertIn("not the prefix's genome list",err)
  ref=self.reference("wrong-revision",lambda r:r["evaluation"].update(revision=REV-1))
  code,_,err=self.run_diff(ref); self.assertEqual(code,2); self.assertIn(f"not {REV}",err)
  ref=self.reference("missing-dataset"); next((ref/"species"/SLUG/"results").glob("*isoniazid.json")).unlink()
  code,_,err=self.run_diff(ref); self.assertEqual(code,2); self.assertIn("dataset sets differ",err)
  self.assertFalse((ref.parent/"diff.json").exists())

 def test_an_unfilled_revision_n_stops_the_cli(self):
  with tempfile.TemporaryDirectory() as tmp:
   lane=Path(tmp)/"mtb-4k"; lane.mkdir()
   for name in ("diff_reference_prefix.py","revision_n.py"): shutil.copyfile(HERE.parent/name,lane/name)
   text=(lane/"revision_n.py").read_text(); start=text.index("\nVALUES = {\n")+1; end=text.index("\n}\n",start)+3
   (lane/"revision_n.py").write_text(text[:start]+"VALUES = {\n"+"".join(f'    "{n}": "@@{n}@@",\n' for n in diff.rev_n.TABLE)+"}\n"+text[end:])
   ref=self.reference("cli")
   p=subprocess.run([sys.executable,"-I","-B",str(lane/"diff_reference_prefix.py"),"--prefix",str(self.prefix),"--reference",str(ref),
                     "--out",str(Path(tmp)/"diff.json")],cwd=tmp,capture_output=True,text=True)
   self.assertEqual(p.returncode,2,p.stderr); self.assertIn("diff_reference_prefix: STOP: unfilled placeholder(s): @@REVISION_N@@",p.stderr)
   self.assertFalse((Path(tmp)/"diff.json").exists())

if __name__=="__main__":
 result=unittest.main(exit=False,verbosity=2).result
 sys.exit(0 if result.wasSuccessful() else 1)
