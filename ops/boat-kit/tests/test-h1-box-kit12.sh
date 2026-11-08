#!/usr/bin/env bash
# KIT12: h1-box.sh's and r11-box.sh's add-on role (H1_ROLE=addon; G2(b), run.py --reference-from) against KIT10's fake
# run.py (no grm, no box data; every byte synthetic), under H1_HOME. An add-on's union, verify and run pass
# --reference-from the main part's pushed WORK after that set matches its pinned manifest and holds that one finished
# part; a run box passes none; r11-box.sh runs h1-box.sh run for an add-on and reads h1-run's guard file. One
# "RESULT id want got verdict detail" line per check, then "FAILS n"; exits 1 on any failure. The fake run.py and the
# helpers are test-h1-box.sh's, copied.
# usage: test-h1-box-kit12.sh WORK KIT_DIR
set -u
WK=${1:?work} KIT=${2:?kit}
rm -rf "$WK"; mkdir -p "$WK"
FAILS=0
res() { local v=as-predicted; [ "$2" = "$3" ] || { v=DIFFERS; FAILS=$((FAILS + 1)); }; printf 'RESULT %s want=%s got=%s %s %s\n' "$1" "$2" "$3" "$v" "${4:-}"; }
yn() { if "$@"; then echo yes; else echo no; fi; }
# the fake run.py: --estimate-only writes each root's genomes.txt (FAKE_UNION) and an estimate whose
# union_without_fasta counts the ids with no RAW/fasta/<id>.fna.gz (FAKE_UWF overrides) and whose ok is FAKE_OK;
# the stderr gate line is FAKE_GATE. A run writes roots with ksets, a matrix, a curated tree, a link and arm outputs.
mkfake() { # HOME
  mkdir -p "$1/grm/scripts/input-comparison/h1"
  cat > "$1/grm/scripts/input-comparison/run.py" << 'EOF'
import json, os, sys
a = sys.argv[1:]; g = lambda k: a[a.index(k) + 1]
raw, work = g("--raw"), g("--work")
open(os.environ["FAKE_ARGS"], "a").write(" ".join(a) + "\n")
os.makedirs(work, exist_ok=True)
ids = open(os.environ["FAKE_UNION"]).read().split()
if "--estimate-only" in a:
    for r, part in (("r1", ids[:1]), ("r2", ids)):
        d = os.path.join(work, "roots", r, "species", "s"); os.makedirs(d, exist_ok=True)
        open(os.path.join(d, "genomes.txt"), "w").write("\n".join(part) + "\n")
    miss = [i for i in ids if not os.path.isfile(os.path.join(raw, "fasta", i + ".fna.gz"))]
    uwf = int(os.environ.get("FAKE_UWF", len(miss)))
    print(json.dumps({"ok": os.environ.get("FAKE_OK", "1") == "1", "union_without_fasta": uwf}))
    if os.environ.get("FAKE_GATE", "1") == "1": print("revision-12 gate: PASSED (fake)", file=sys.stderr)
    sys.exit(0)
rc = int(os.environ.get("FAKE_RUN_RC", "0"))
for r in ("pool-k31", "default-k31"):
    top = os.path.join(work, "roots", r); sp = os.path.join(top, "species", "s")
    for x in ("ksets", "matrix", "curated"):
        os.makedirs(os.path.join(sp if x != "curated" else top, x), exist_ok=True)
    for i in ids: open(os.path.join(sp, "ksets", i + ".kset"), "wb").write(b"k" * 1000)
    open(os.path.join(sp, "matrix", "m.bin"), "wb").write(b"m" * 500)
    open(os.path.join(top, "curated", "provenance.jsonl"), "w").write("{}\n")
    open(os.path.join(sp, "genomes.txt"), "w").write("\n".join(ids) + "\n")
    open(os.path.join(top, "built.json"), "w").write("{}")
    if not os.path.islink(os.path.join(top, "fasta")): os.symlink(os.path.join(raw, "fasta"), os.path.join(top, "fasta"))
d = os.path.join(work, "arms", "default", "ds"); os.makedirs(d, exist_ok=True)
open(os.path.join(d, "predictions.jsonl"), "w").write("{}\n")
for f in ("summary.json", "spec.json"): open(os.path.join(work, f), "w").write("{}")
json.dump({"default_kset_bytes": 1000}, open(os.path.join(work, "disk-estimate.json"), "w"))
sys.exit(rc)
EOF
  echo '{"species": "s"}' > "$1/grm/scripts/input-comparison/h1/real.json"
  echo '{"species": "s", "synthetic": true}' > "$1/grm/scripts/input-comparison/h1/syn.json"
  printf '#!/bin/sh\nexit 0\n' > "$1/grm/grm"
}
# set_ NAME [SPEC]: a fresh home with the three pushed sets (metadata, vocabulary, the union's FASTAs) and stage.env
set_() {
  H=$WK/$1; chmod -R u+w "$H" 2> /dev/null; rm -rf "$H"; mkdir -p "$H/r11"; mkfake "$H"
  mkdir -p "$H/grm-meta/bvbrc" "$H/grm-meta/catalog" "$H/grm-meta/species-check" "$H/grm-vocab/vocabulary" "$H/grm-fasta/fasta"
  echo '{}' > "$H/grm-meta/bvbrc/genomes.jsonl"; echo v > "$H/grm-meta/catalog/version.txt"; echo '{}' > "$H/grm-meta/species-check/species_check.json"
  printf 'a\n' > "$H/grm-vocab/vocabulary/a.tsv"
  printf '1.1\n2.1\n' > "$H/union"; for i in 1.1 2.1; do echo ">$i" > "$H/grm-fasta/fasta/$i.fna.gz"; done
  export FAKE_UNION=$H/union FAKE_ARGS=$H/args H1_HOME=$H W=$H/r11
  unset FAKE_UWF FAKE_OK FAKE_GATE FAKE_RUN_RC
  SP=${2:-real}
}
man() { (cd "$1" && find . -type f ! -name .grm-input-manifest.sha256 | sed 's#^\./##' | LC_ALL=C sort | xargs sha256sum) > "$1/.grm-input-manifest.sha256"; sha256sum "$1/.grm-input-manifest.sha256" | cut -d' ' -f1; }
env_() { printf 'SPEC=scripts/input-comparison/h1/%s.json\nSOURCE_COMMIT=%s\nH1_ROLE=run\nH1_JOBS=2\nH1_TRACE=%s\nH1_EXTRA=%s\nH1_META_SHA256=%s\nH1_VOCAB_SHA256=%s\nH1_FASTA_SHA256=%s\n' \
  "$SP" "$(printf '0%.0s' {1..40})" "${TR:-0}" "${EX:-}" "$(man "$H/grm-meta")" "$(man "$H/grm-vocab")" "$(man "$H/grm-fasta")" > "$H/r11/stage.env"; }
hb() { bash "$KIT/h1-box.sh" "$@" > "$H/out" 2>&1; }
stop() { grep -o 'STOP .*' "$H/out" | head -n 1 | cut -c1-110; }
# addon_ [PIN]: the main part's set at grm-combine (a master manifest copy and parts/main/h1work), its manifest, and
# stage.env's add-on lines (H1_REFERENCE_SHA256 is the set's manifest sha256, or PIN)
addon_() {
  local m=$H/grm-combine/parts/main/h1work
  mkdir -p "$m/arms/default/ds" "$m/k"
  echo '{}' > "$m/summary.json"; echo '{}' > "$m/spec.json"; echo '{}' > "$m/k/manifest.json"; echo '{}' > "$m/arms/default/ds/runs.json"
  echo "$(printf 'c%.0s' {1..64})  bvbrc/genomes.jsonl" > "$H/grm-combine/MASTER-MANIFEST-test.sha256"
  local s; s=$(man "$H/grm-combine")
  sed -i 's/^H1_ROLE=.*/H1_ROLE=addon/' "$H/r11/stage.env"
  printf 'H1_REFERENCE=%s\nH1_REFERENCE_SHA256=%s\n' "${REFN-main}" "${1:-$s}" >> "$H/r11/stage.env"
}
REF="--reference-from $WK/%s/grm-combine/parts/main/h1work"
cnt() { local n=0; [ -s "$1" ] && n=$(grep -c . "$1"); echo "$n"; }   # one line, even for no file
refline() { grep -- "$1" "$H/args" | grep -c -- "$(printf -- "$REF" "$2")"; }

set_ b1; env_; addon_; hb union; rc=$?
res B1-addon-union-passes-the-reference "0,1" "$rc,$(refline "--work $H/h1union " b1)" "$(grep '^UNION' "$H/out")"
hb verify; rc=$?
res B2-addon-verify-passes-the-reference "0,1" "$rc,$(refline "--work $H/h1est " b1)" "$(grep '^verified' "$H/out" | cut -c1-80)"
hb run; rc=$?
res B3-addon-run-passes-the-reference "0,1" "$rc,$(grep -v -- '--estimate-only' "$H/args" | grep -c -- "$(printf -- "$REF" b1)")" "$(tail -n 1 "$H/out")"
hb ckpt c1; rc=$?
res B4-addon-checkpoint-is-a-run-s "0,yes,no" "$rc,$(yn test -f "$H/ckpt/c1/h1work/arms/default/ds/predictions.jsonl"),$(yn test -e "$H/ckpt/c1/h1combine")" ""
res B4-the-set-is-not-written "yes" "$(yn bash -c "cd '$H/grm-combine' && sha256sum -c --strict --quiet .grm-input-manifest.sha256")" ""
set_ b5; env_; addon_ 0000; hb union; rc=$?
res B5-reference-pin-mismatch-2-no-run "2,0" "$rc,$(cnt "$H/args")" "$(stop)"
set_ b6; env_; addon_; echo changed >> "$H/grm-combine/parts/main/h1work/summary.json"; hb union; rc=$?
res B6-reference-file-changed-2 "2,0" "$rc,$(cnt "$H/args")" "$(stop)"
set_ b7; mkdir -p "$H/grm-combine/parts/other/h1work"; echo '{}' > "$H/grm-combine/parts/other/h1work/summary.json"; env_; addon_; hb union; rc=$?
res B7-a-second-part-2 "2,1" "$rc,$(grep -c 'holds more or other than the main part main' "$H/out")" "$(stop)"
set_ b8; env_; REFN= addon_; hb union; rc=$?
res B8-no-reference-name-2 "2,1" "$rc,$(grep -c 'an add-on needs H1_REFERENCE' "$H/out")" "$(stop)"
set_ b9; env_; REFN=../main addon_; hb union; rc=$?
res B9-reference-name-a-path-2 "2,1" "$rc,$(grep -c 'an add-on needs H1_REFERENCE' "$H/out")" "$(stop)"
set_ b10; env_; addon_; rm "$H/grm-combine/parts/main/h1work/summary.json"; man "$H/grm-combine" > /dev/null
sed -i "s/^H1_REFERENCE_SHA256=.*/H1_REFERENCE_SHA256=$(sha256sum "$H/grm-combine/.grm-input-manifest.sha256" | cut -d' ' -f1)/" "$H/r11/stage.env"; hb union; rc=$?
res B10-unfinished-main-part-2 "2,1" "$rc,$(grep -c 'no finished WORK of the main part' "$H/out")" "$(stop)"
set_ b11; env_; hb union; hb verify; hb run; printf 'H1_REFERENCE=main\n' >> "$H/r11/stage.env"; hb union; rc=$?
res B11-a-run-box-passes-no-reference "0,0" "$rc,$(grep -c -- '--reference-from' "$H/args")" "(H1_REFERENCE is ignored for H1_ROLE=run)"
set_ b12; env_; addon_; hb union; hb verify; chmod -R u+w "$H/grm-combine"; rm -rf "$H/grm-combine/parts/main"; hb run; rc=$?
res B12-run-rechecks-the-set "2,0" "$rc,$(grep -v -- '--estimate-only' "$H/args" | grep -c .)" "$(stop)"

# r11-box.sh: an add-on's stage is h1-box.sh run, and its guard's pause file is h1-run's
rb() { bash "$KIT/r11-box.sh" "$@" > "$H/rb.out" 2>&1; }
poll() { local s; for _ in $(seq 60); do s=$(bash "$KIT/r11-box.sh" status 2> /dev/null); case "$s" in "STATE running"*) sleep 0.5;; *) break;; esac; done; echo "$s"; }
h1env() { env_; printf 'KIND=h1\nSTAGE=t1\nH1_PARTS=\nH1_COMBINE_SHA256=\n' >> "$H/r11/stage.env"; }
set_ d1; h1env; addon_; hb union; rb verify; rc=$?
res D1-addon-verify-goes-to-h1-box "0,1" "$rc,$(refline "--work $H/h1est " d1)" "$(grep '^verified' "$H/rb.out" | cut -c1-80)"
rb start; s=$(poll)
res D2-addon-start-runs-h1-box-run-to-done "done,1,yes" "$(awk '{print $2}' <<< "$s"),$(grep -v -- '--estimate-only' "$H/args" | grep -c -- "$(printf -- "$REF" d1)"),$(yn grep -q 'H-1 addon of scripts/input-comparison/h1/real.json' "$H/r11/status")" "$(cut -c1-90 <<< "$s")"
set_ d3; h1env; addon_; setsid sleep 60 & sp=$!; echo "$sp" > "$H/r11/stage.pid"; mkdir -p "$H/r11/gates"
echo "paused 2026-11-07T12:00:00Z disk_free=1073741824 eff_mem=99999999999" > "$H/r11/gates/h1-run.tripped"
s=$(bash "$KIT/r11-box.sh" status)
res D3-addon-reads-h1-run-s-pause paused "$(awk '{print $2}' <<< "$s")" "$(cut -c1-110 <<< "$s")"
kill "$sp" 2> /dev/null; wait "$sp" 2> /dev/null
set_ d4; h1env; addon_; sed -i 's/^H1_ROLE=.*/H1_ROLE=other/' "$H/r11/stage.env"; rb start; s=$(poll)
res D4-bad-role-names-addon "stopped,1" "$(awk '{print $2}' <<< "$s"),$(grep -c 'H1_ROLE is other, not run, addon or combine' "$H/r11/status")" "$(cut -c1-90 <<< "$s")"
chmod -R u+w "$WK" 2> /dev/null
echo "FAILS $FAILS"
[ "$FAILS" = 0 ]
