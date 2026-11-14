# GRM's export for Givens

This directory holds GRM's entries for [Givens](https://github.com/editnori/givens), a shared browser of the definitions and assumptions of four projects. They are in one file, `grm.json`. `scripts/givens-export/export.py` writes it from `scripts/givens-export/spec.json`, from sentences of this repository's documents and from seven committed summary files. [docs/givens-export.md](../../docs/givens-export.md) describes each field, what is missing and what is left out.

| File | sha256 |
| --- | --- |
| `grm.json` | `f65dc03dffb066673ad22557d6338d3424f4f62146f62170d427d2eba2384f32` |

Write it again from the repository root, with Python 3.12 and the standard library only:

```sh
python3 -I scripts/givens-export/export.py > reference/givens-export/grm.json
python3 -I scripts/givens-export/export.py --check reference/givens-export/grm.json
```

Then put the new sha256 in the table above. The file names no commit, because a committed file cannot name its own. Givens records the GRM commit it copied the file from.
