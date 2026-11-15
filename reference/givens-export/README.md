# GRM's export for Givens

This directory holds GRM's entries for [Givens](https://github.com/editnori/givens), a shared browser of the definitions and assumptions of four projects. They are in one file, `grm.json`. `scripts/givens-export/export.py` writes it from `scripts/givens-export/spec.json`, from sentences of this repository's documents and from seven committed summary files. [docs/givens-export.md](../../docs/givens-export.md) describes each field, what is missing and what is left out. [docs/givens-terms.md](../../docs/givens-terms.md) shows the same entries as a page to read.

`SHA256SUMS` holds the file's sha256 on one line, as `sha256sum` prints it from the repository root. Givens can store that list beside its copy of the file.

Write both files and that page again from the repository root, with Python 3.12 and the standard library only:

```sh
python3 -I scripts/givens-export/export.py --write
python3 -I scripts/givens-export/export.py --check
sha256sum -c reference/givens-export/SHA256SUMS
```

The second command prints nothing and exits 1 unless both files and the page are what the script writes. The file names no commit, because a committed file cannot name its own. Givens records the GRM commit it copied the file from.
