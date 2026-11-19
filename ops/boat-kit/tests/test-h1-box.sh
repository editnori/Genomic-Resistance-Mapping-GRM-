#!/usr/bin/env bash
# KIT10: h1-box.sh's own checks against a fake run.py (no grm, no box data; every byte synthetic), under H1_HOME.
# One "RESULT id want got verdict detail" line per check, then "FAILS n"; exits 1 on any failure. KIT10 round 2
# (REVIEW-KIT10): V7, V8 and R4 (run.py gets the free disk less guard.sh's pause floor and 1 GiB, L3), M3 (a marker
# that can't move stops with 2), and D1-D6, r11-box.sh's H-1 dispatch (verify, start to done, a bad role, paused from
# the guard's tripped file, stop, and the combine role).
# usage: test-h1-box.sh WORK KIT_DIR
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

set_ u1; env_; hb union; rc=$?
res U1-union "0,yes" "$rc,$(yn cmp -s "$H/union" "$H/r11/receipts/h1/union.txt")" "$(grep '^UNION' "$H/out")"
res U1-raw0-read-only no "$(yn test -w "$H/h1raw0/bvbrc")" ""
hb union; res U1-rerun 0 $? "(a second union removes the read-only RAW0 first)"
set_ u2; env_; export FAKE_UWF=5; hb union; rc=$?
res U2-union-count-differs-stops 4 "$rc" "$(stop)"
set_ u3; env_; sed -i 's/^H1_META_SHA256=.*/H1_META_SHA256=0000/' "$H/r11/stage.env"; hb union; rc=$?
res U3-metadata-pin-mismatch 2 "$rc" "$(stop)"
set_ u4; mkdir -p "$H/grm-vocab/other"; echo x > "$H/grm-vocab/other/x"; env_; hb union; rc=$?
res U4-vocabulary-set-holds-more 2 "$rc" "$(stop)"
set_ u5; env_; echo changed >> "$H/grm-meta/catalog/version.txt"; hb union; rc=$?
res U5-pushed-file-differs-from-its-manifest 2 "$rc" "$(stop)"
set_ v1; env_; hb union; hb verify; rc=$?
res V1-verify "0,no,no" "$rc,$(yn test -w "$H/h1raw/fasta"),$(yn test -w "$H/h1raw/fasta/1.1.fna.gz")" "$(grep '^RAW' "$H/out")"
res V1-estimate-at-measured-free yes "$(yn grep -q -- "--raw $H/h1raw --work $H/h1est --estimate-only --free-disk-bytes [0-9]" "$H/args")" ""
set_ v2; echo ">3.1" > "$H/grm-fasta/fasta/3.1.fna.gz"; env_; hb union; hb verify; rc=$?
res V2-fasta-set-not-the-union 2 "$rc" "$(stop)"
set_ v3; mkdir -p "$H/grm-fasta/catalog"; echo y > "$H/grm-fasta/catalog/y.txt"; env_; hb union; hb verify; rc=$?
res V3-fasta-set-holds-more 2 "$rc" "$(stop)"
set_ v4; env_; hb union; export FAKE_OK=0; hb verify; rc=$?
res V4-estimate-not-ok 4 "$rc" "$(stop)"
set_ v5; env_; hb union; export FAKE_GATE=0; hb verify; rc=$?
res V5-real-spec-without-the-gate-line 4 "$rc" "$(stop)"
set_ v6 syn; env_; hb union; export FAKE_GATE=0; hb verify; rc=$?
res V6-synthetic-spec-needs-no-gate-line 0 "$rc" ""
set_ m1; mkdir -p "$H/grm-meta/receipts"; echo '{}' > "$H/grm-meta/receipts/truth.json"; env_; hb union; rc=$?
res M1-marker-on-a-real-spec-refused 2 "$rc" "$(stop)"
set_ m2 syn; mkdir -p "$H/grm-meta/receipts"; echo '{}' > "$H/grm-meta/receipts/truth.json"; env_; hb union; hb verify; rc=$?
res M2-synthetic-marker-back-at-raw-top "0,yes,no" "$rc,$(yn test -f "$H/h1raw/truth.json"),$(yn test -e "$H/h1raw/receipts")" ""
set_ e1; EX=--foo env_; hb union; rc=$?
res E1-other-extra-refused 2 "$rc" "$(stop)"
set_ e1b syn; EX=--foo env_; hb union; rc=$?   # synthetic, so only the H1_EXTRA check can refuse it
res E1b-other-extra-on-a-synthetic-spec-refused 2 "$rc" "$(stop)"
set_ e2; EX=--below-panel-floor env_; hb union; rc=$?
res E2-below-floor-on-a-real-spec-refused 2 "$rc" "$(stop)"
set_ e3 syn; EX=--below-panel-floor env_; hb union; rc=$?
res E3-below-floor-on-a-synthetic-spec "0,1" "$rc,$(grep -c -- '--below-panel-floor' "$H/args")" ""
set_ r1; TR=1 env_; hb union; hb verify; hb run; rc=$?
res R1-run "0,1" "$rc,$(grep -c -- '--trace-training' "$H/args")" "$(tail -n 1 "$H/out")"
res R1-disk-measured "yes,1000" "$(yn test -s "$H/r11/receipts/h1/disk-measured.json"),$(python3 -I -c 'import json, sys; print(json.load(open(sys.argv[1]))["pool"]["mean_kset_apparent_bytes"])' "$H/r11/receipts/h1/disk-measured.json" 2> /dev/null)" ""
hb ckpt c1; rc=$?; C=$H/ckpt/c1
res R1-ckpt "0,yes" "$rc,$(yn bash -c "cd '$C' && sha256sum -c --strict --quiet SHA256SUMS")" "$(tail -n 1 "$H/out")"
res R1-ckpt-no-ksets-matrix-links-curated "0,0,0,0" "$(find "$C" -path '*ksets*' | wc -l),$(find "$C" -path '*matrix*' | wc -l),$(find "$C" -type l | wc -l),$(find "$C" -path '*curated*' | wc -l)" ""
res R1-ckpt-keeps-arms-lists-records "yes,yes,yes" "$(yn test -f "$C/h1work/arms/default/ds/predictions.jsonl"),$(yn test -f "$C/h1work/roots/pool-k31/species/s/genomes.txt"),$(yn test -f "$C/h1work/roots/pool-k31/built.json")" ""
set_ r2; env_; hb union; hb verify; export FAKE_RUN_RC=3; hb run; rc=$?
res R2-run-fails 1 "$rc" "$(stop)"
set_ r3; env_; hb union; hb verify; chmod u+w "$H/h1raw"; hb run; rc=$?
res R3-writable-raw-refused 2 "$rc" "$(stop)"

# REVIEW-KIT10 L3: a df shim reads FAKE_DF bytes free; with DISK_PAUSE_GIB 30, run.py gets FAKE_DF - 31 GiB
SH=$WK/shim; mkdir -p "$SH"
printf '#!/bin/sh\nprintf "    Avail\\n%%s\\n" "$FAKE_DF"\n' > "$SH/df"
printf '#!/bin/sh\ncase "$1" in */receipts/truth.json) echo "mv: refused (test)" >&2; exit 1;; esac\nexec %s "$@"\n' "$(command -v mv)" > "$SH/mv"
chmod +x "$SH/df" "$SH/mv"
G100=107374182400 U100=74088185856
set_ v7; env_; hb union; PATH="$SH:$PATH" FAKE_DF=$G100 hb verify; rc=$?
res V7-estimate-at-usable-free "0,1,$U100" "$rc,$(grep -- "--work $H/h1est " "$H/args" | grep -cE -- "--free-disk-bytes $U100( |\$)"),$(cat "$H/r11/receipts/h1/usable-at-estimate.txt" 2> /dev/null)" "$(grep -- "--work $H/h1est " "$H/args" | grep -o -- '--free-disk-bytes [0-9]*' | head -n 1)"
set_ v8; env_; hb union; PATH="$SH:$PATH" FAKE_DF=21474836480 hb verify; rc=$?
res V8-free-below-the-pause-floor-4 "4,0" "$rc,$(grep -c -- "--work $H/h1est " "$H/args")" "$(stop)"
set_ r4; env_; hb union; hb verify; PATH="$SH:$PATH" FAKE_DF=$G100 hb run; rc=$?
res R4-run-at-usable-free "0,1" "$rc,$(grep -v -- '--estimate-only' "$H/args" | grep -cE -- "--free-disk-bytes $U100( |\$)")" "$(grep -v -- '--estimate-only' "$H/args" | grep -o -- '--free-disk-bytes [0-9]*' | head -n 1)"
# REVIEW-KIT10 note: a synthetic marker that can't move back stops with 2 (round 1 went on with receipts/ in RAW)
set_ m3 syn; mkdir -p "$H/grm-meta/receipts"; echo '{}' > "$H/grm-meta/receipts/truth.json"; env_; PATH="$SH:$PATH" hb union; rc=$?
res M3-marker-that-cannot-move-2 2 "$rc" "$(stop)"

# r11-box.sh's H-1 dispatch (REVIEW-KIT10 gap): KIND=h1 in stage.env, as the drive writes it
rb() { bash "$KIT/r11-box.sh" "$@" > "$H/rb.out" 2>&1; }
poll() { local s; for _ in $(seq 60); do s=$(bash "$KIT/r11-box.sh" status 2> /dev/null); case "$s" in "STATE running"*) sleep 0.5;; *) break;; esac; done; echo "$s"; }
# the drive writes every H1_ variable, empty where the role has none (drive.sh step 5)
h1env() { env_; printf 'KIND=h1\nSTAGE=t1\nH1_PARTS=\nH1_COMBINE_SHA256=\n' >> "$H/r11/stage.env"; }
set_ d1; h1env; hb union; rb verify; rc=$?
res D1-verify-goes-to-h1-box "0,1,1" "$rc,$(grep -c '^verified: ' "$H/rb.out"),$(grep -c -- "--work $H/h1est " "$H/args")" "$(grep '^verified' "$H/rb.out" | cut -c1-100)"
rb start; s=$(poll)
res D2-start-runs-h1-box-run-to-done "done,1,yes" "$(awk '{print $2}' <<< "$s"),$(grep -v -- '--estimate-only' "$H/args" | grep -c -- '--work '),$(yn grep -q 'H-1 run of scripts/input-comparison/h1/real.json' "$H/r11/status")" "$(cut -c1-90 <<< "$s")"
set_ d3; h1env; sed -i 's/^H1_ROLE=.*/H1_ROLE=other/' "$H/r11/stage.env"; hb union; rb start; s=$(poll)
res D3-bad-role-stops "stopped,1,0" "$(awk '{print $2}' <<< "$s"),$(grep -c 'H1_ROLE is other, not run, addon or combine' "$H/r11/status"),$(grep -vc -- '--estimate-only' "$H/args")" "$(cut -c1-90 <<< "$s")"
set_ d4; h1env; setsid sleep 60 & sp=$!; echo "$sp" > "$H/r11/stage.pid"; mkdir -p "$H/r11/gates"
echo "paused 2026-11-07T12:00:00Z disk_free=1073741824 eff_mem=99999999999" > "$H/r11/gates/h1-run.tripped"
s=$(bash "$KIT/r11-box.sh" status)
res D4-guard-pause-reads-paused "paused,yes" "$(awk '{print $2}' <<< "$s"),$(yn grep -q '| guard disk: 1073741824 B free, below DISK_PAUSE_GIB 30 GiB (paused ' <<< "$s")" "$(cut -c1-110 <<< "$s")"
echo "resumed 2026-11-07T12:05:00Z disk_free=99999999999 eff_mem=99999999999" >> "$H/r11/gates/h1-run.tripped"
s=$(bash "$KIT/r11-box.sh" status)
res D4b-resumed-reads-running running "$(awk '{print $2}' <<< "$s")" ""
echo "$sp" > "$H/r11/gates/pgids"; rb stop; rc=$?
res D5-stop-ends-the-group "0,no,yes" "$rc,$(yn kill -0 "$sp" 2> /dev/null),$(yn test -e "$H/r11/stop")" "$(head -c 90 "$H/rb.out")"
wait "$sp" 2> /dev/null
set_ d6; h1env; sed -i 's/^H1_ROLE=.*/H1_ROLE=combine/; s/^H1_PARTS=.*/H1_PARTS=p1/' "$H/r11/stage.env"; rb start; s=$(poll)
res D6-combine-role-runs-h1-box-combine "stopped,1" "$(awk '{print $2}' <<< "$s"),$(grep -c 'H-1 combine failed: STOP ' "$H/r11/status")" "$(grep -o 'H-1 combine failed.*' "$H/r11/status" | cut -c1-90)"
chmod -R u+w "$WK" 2> /dev/null
echo "FAILS $FAILS"
[ "$FAILS" = 0 ]
