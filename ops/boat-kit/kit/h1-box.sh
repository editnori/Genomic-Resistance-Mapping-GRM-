#!/usr/bin/env bash
# H-1 on one Boat box (G5): the input comparison's run.py, driven by the laptop's drive.sh (KIND=h1) through
# grm-boat.sh run. r11-box.sh keeps host (build, speed check, CPU and binary sha256s), start, status and stop, and
# hands KIND=h1's verify and stage to this script. Settings come from W/stage.env: SPEC (the spec, relative to the
# pushed source), H1_ROLE (run|addon|combine), H1_TRACE (1: --trace-training), H1_JOBS, H1_PARTS (combine: the part
# folders under /home/user/grm-combine/parts, comma list, in --work-from order), H1_EXTRA (only --below-panel-floor,
# only for a spec that says "synthetic": true), SOURCE_COMMIT.
# KIT12, H1_ROLE=addon (G2(b)): a run box whose union, verify and run pass run.py --reference-from the main part's
# pulled WORK, /home/user/grm-combine/parts/$H1_REFERENCE/h1work; that set must match its manifest, pinned in stage.env
# as H1_REFERENCE_SHA256, and hold that one part. Its checkpoints are a run's.
#   prep      python3-numpy (J needs it), recorded with its version
#   union     the spec's genome union, from run.py --estimate-only on RAW0 (the pushed metadata, no FASTA): every
#             root's genomes.txt for the spec's species (a real root lists every species it holds), which must count
#             what the estimate reports as union_without_fasta; written to receipts/h1/union.txt for the laptop, which
#             derives the FASTA set from the master manifest
#   verify    every pushed set against its manifest, RAW = /home/user/h1raw hardlinked from the metadata and FASTA sets
#             (exactly the union's FASTAs), made read only, then run.py --estimate-only with the usable free disk
#             (/home/user's free disk minus guard.sh's pause floor DISK_PAUSE_GIB and 1 GiB, REVIEW-KIT10 L3) in a
#             scratch WORK: rc 0, ok, 0 union genomes without a FASTA, and the gate line. The union's and the
#             estimate's scratch WORKs are removed, so their disk goes to the run
#   run       (from r11-box.sh's detached stage) run.py on RAW into WORK under pg-run.sh's guard, with the same usable
#             free disk, then measure
#   measure   receipts/h1/disk-measured.json: disk blocks per root (ksets, matrix, the rest), all roots together
#             (hardlinked ksets once), the pool's mean kset, arm outputs, peak RSS, against disk-estimate.json
#   combine   (combine box) run.py combine over the pushed parts with the pushed master manifest; with one part, its
#             J files must equal the part's byte for byte and K's manifest and report after timings (acceptance 11)
#   ckpt L    /home/user/ckpt/L: the run's receipts (W) and WORK (or the combine's output) without ksets/, matrix/,
#             links or roots' curated trees, with SHA256SUMS over every file
#   prune DIR (KIT12, a box made from the data snapshot) after the FASTA set is staged: the box's own copy of the
#             snapshot's master root, DIR, once the sets match their pins and the metadata set is linked from DIR; the
#             sets keep their files, so the blocks only DIR held are freed (receipts/h1/prune.txt)
# Nothing here writes under the pushed sets; RAW is chmod a-w before any run.py call reads it.
# usage: h1-box.sh prep|union|verify|run|measure|combine|ckpt LABEL|prune DIR
set -u -o pipefail
HB=${H1_HOME:-/home/user}   # the box's home; a test points it elsewhere
W=${W:-$HB/r11}
set -a; . "$W/stage.env"; set +a
K=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
S=${SRC:-$HB/grm}
GRM=$S/target/release/grm RUNPY=$S/scripts/input-comparison/run.py
META=$HB/grm-meta FA=$HB/grm-fasta VOC=$HB/grm-vocab CI=$HB/grm-combine
RAW0=$HB/h1raw0 RAW=$HB/h1raw WORK=$HB/h1work OUT=$HB/h1combine
R=$W/receipts/h1
mkdir -p "$R" "$W/logs" "$W/gates"
log() { printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$W/logs/driver.log"; }
fail() { echo "STOP $*" >&2; log "h1 STOP $*"; exit "${RC:-1}"; }
[ -f "$S/$SPEC" ] || RC=2 fail "no spec $SPEC in the pushed source"
SPECF=$S/$SPEC
# H1_EXTRA: --below-panel-floor only, and only on a synthetic spec (the dry run)
EXTRA=()
if [ -n "${H1_EXTRA:-}" ]; then
  [ "$H1_EXTRA" = --below-panel-floor ] || RC=2 fail "H1_EXTRA may only be --below-panel-floor, not $H1_EXTRA"
  python3 -I -c 'import json, sys; sys.exit(0 if json.load(open(sys.argv[1])).get("synthetic") is True else 1)' "$SPECF" ||
    RC=2 fail "--below-panel-floor on a spec that isn't synthetic"
  EXTRA=(--below-panel-floor)
fi
TRACE=(); [ "${H1_TRACE:-0}" = 1 ] && TRACE=(--trace-training)
# KIT12: an add-on reads its main part's WORK (--reference-from); run.py checks it before anything is counted
REFARG=() REFW=""
if [ "${H1_ROLE:-run}" = addon ]; then
  [[ "${H1_REFERENCE:-}" =~ ^[a-z0-9][a-z0-9._-]*$ ]] || RC=2 fail "an add-on needs H1_REFERENCE, its main part's name: ${H1_REFERENCE:-unset}"
  REFW=$CI/parts/$H1_REFERENCE/h1work; REFARG=(--reference-from "$REFW")
fi
checkset() { # DIR: the pushed set matches its manifest, and its manifest is the one the laptop pinned in stage.env
  local d=$1 want=$2
  [ -f "$d/.grm-input-manifest.sha256" ] || RC=2 fail "$d has no manifest"
  [ "$(sha256sum "$d/.grm-input-manifest.sha256" | cut -d' ' -f1)" = "$want" ] || RC=2 fail "$d's manifest is not the pinned $want"
  (cd "$d" && sha256sum -c --strict --quiet .grm-input-manifest.sha256) || RC=2 fail "sha256sum -c failed in $d"
  echo "$(date -u +%FT%TZ) $d manifest $want lines=$(wc -l < "$d/.grm-input-manifest.sha256") rc=0" >> "$R/manifest-check.txt"
}
# refset: an add-on's pushed main-part set matches its pinned manifest and holds that one part's WORK
refset() {
  [ -z "$REFW" ] && return 0
  checkset "$CI" "${H1_REFERENCE_SHA256:-}"
  [ "$(find "$CI/parts" -mindepth 1 -maxdepth 1 -printf '%f\n' 2> /dev/null)" = "$H1_REFERENCE" ] || RC=2 fail "$CI/parts holds more or other than the main part $H1_REFERENCE"
  [ -d "$REFW" ] && [ -f "$REFW/summary.json" ] || RC=2 fail "no finished WORK of the main part at $REFW"
}
# lnset SRC DST: DST gets hard links of SRC's files (not its manifest), so RAW is a plain tree on the same disk
lnset() { (cd "$1" && find . -type f ! -name .grm-input-manifest.sha256 -printf '%P\n') | while IFS= read -r f; do
  mkdir -p "$2/$(dirname "$f")" && ln "$1/$f" "$2/$f" || exit 1; done; }
# marker DIR: a synthetic root's truth.json comes as receipts/truth.json (a receipt class); it goes back to DIR/truth.json
# for a synthetic spec only, and a real spec's metadata may hold no receipts/ at all
marker() {
  [ -e "$1/receipts" ] || return 0
  python3 -I -c 'import json, sys; sys.exit(0 if json.load(open(sys.argv[1])).get("synthetic") is True else 1)' "$SPECF" ||
    RC=2 fail "the metadata holds receipts/, which only a synthetic root's truth.json may"
  [ "$(find "$1/receipts" -mindepth 1)" = "$1/receipts/truth.json" ] || RC=2 fail "receipts/ holds more than truth.json"
  # REVIEW-KIT10 note: a failed move would leave receipts/ in RAW
  mv "$1/receipts/truth.json" "$1/truth.json" && rmdir "$1/receipts" || RC=2 fail "could not move receipts/truth.json to $1/truth.json"
}
# usable FREE: FREE minus guard.sh's pause floor (DISK_PAUSE_GIB, 30) and 1 GiB, as r11-box.sh holds back before a
# species (REVIEW-KIT10 L3): run.py's plan then can't pass while ending where guard.sh stops the group for good.
# Exit 4 when nothing is left.
usable() {
  local u
  [[ "${DISK_PAUSE_GIB:-30}" =~ ^[0-9]+(\.[0-9]+)?$ ]] || RC=2 fail "DISK_PAUSE_GIB is not a number: $DISK_PAUSE_GIB"
  [[ "$1" =~ ^[0-9]+$ ]] || RC=4 fail "no free disk reading: $1"
  u=$(awk -v f="$1" -v p="${DISK_PAUSE_GIB:-30}" 'BEGIN { printf "%.0f", f - (p + 1) * 1073741824 }')
  [[ "$u" =~ ^[0-9]+$ ]] && [ "$u" -gt 0 ] || RC=4 fail "$1 B free is not above guard.sh's pause floor of ${DISK_PAUSE_GIB:-30} GiB plus 1 GiB"
  echo "$u"
}
runpy() { python3 -B "$RUNPY" "$SPECF" --grm "$GRM" --source-commit "$SOURCE_COMMIT" --sigdrug-vocabulary "$VOC/vocabulary" "$@"; }

cmd_prep() {
  if ! python3 -I -c 'import numpy' 2> /dev/null; then
    sudo -n apt-get update -qq > "$W/logs/apt.log" 2>&1
    sudo -n env DEBIAN_FRONTEND=noninteractive apt-get install -y -qq --no-install-recommends python3-numpy >> "$W/logs/apt.log" 2>&1 ||
      RC=3 fail "apt-get install python3-numpy failed (logs/apt.log)"
  fi
  echo "numpy $(python3 -I -c 'import numpy; print(numpy.__version__)') python $(python3 -V 2>&1)" | tee "$R/python.txt"
}

cmd_union() {
  local n u slug
  checkset "$META" "$H1_META_SHA256"; checkset "$VOC" "$H1_VOCAB_SHA256"; refset
  [ "$(find "$VOC" -mindepth 1 -maxdepth 1 ! -name .grm-input-manifest.sha256 -printf '%f\n')" = vocabulary ] || RC=2 fail "$VOC holds more than vocabulary/"
  # run.py's union is the spec's own species (its genomes_txt: ROOT/species/<spec species>/genomes.txt), but a curated
  # root of real metadata lists every species it holds, so the union reads that one species' lists only
  slug=$(python3 -I -c 'import json, sys; s = json.load(open(sys.argv[1])).get("species"); print(s if isinstance(s, str) else "")' "$SPECF") ||
    RC=2 fail "could not read the spec's species"
  [[ "$slug" =~ ^[a-z0-9][a-z0-9_]*$ ]] || RC=2 fail "the spec names no species slug: ${slug:-none}"
  [ -e "$RAW0" ] && chmod -R u+w "$RAW0"
  rm -rf "$RAW0" "$HB/h1union"; mkdir -p "$RAW0"
  lnset "$META" "$RAW0" || RC=2 fail "could not link the metadata into $RAW0"
  marker "$RAW0"
  chmod -R a-w "$RAW0"
  runpy --raw "$RAW0" --work "$HB/h1union" "${EXTRA[@]}" "${REFARG[@]}" --estimate-only --free-disk-bytes 1000000000000000 \
    > "$R/union-estimate.json" 2> "$R/union-estimate.err" || RC=4 fail "run.py --estimate-only on the metadata failed: $(tail -n 2 "$R/union-estimate.err" | tr '\n' ' ')"
  find "$HB/h1union/roots" -mindepth 4 -maxdepth 4 -path "$HB/h1union/roots/*/species/$slug/genomes.txt" -type f -exec cat {} + |
    tr -s ' \t\r' '\n' | grep . | LC_ALL=C sort -u > "$R/union.txt"
  n=$(wc -l < "$R/union.txt")
  u=$(python3 -I -c 'import json, sys; print(json.load(open(sys.argv[1]))["union_without_fasta"])' "$R/union-estimate.json")
  [ "$n" -gt 0 ] && [ "$n" = "$u" ] || RC=4 fail "the roots hold $n genomes, the estimate counts $u union genomes without a FASTA"
  echo "UNION $n $(sha256sum "$R/union.txt" | cut -d' ' -f1)"
}

cmd_verify() {
  local free use n
  checkset "$META" "$H1_META_SHA256"; checkset "$VOC" "$H1_VOCAB_SHA256"; checkset "$FA" "$H1_FASTA_SHA256"; refset
  [ "$(cut -c67- "$FA/.grm-input-manifest.sha256" | sed -n 's#^fasta/\(.*\)\.fna\.gz$#\1#p' | LC_ALL=C sort)" = "$(cat "$R/union.txt")" ] ||
    RC=2 fail "the FASTA set is not exactly the union's FASTAs"
  [ -z "$(cut -c67- "$FA/.grm-input-manifest.sha256" | grep -v '^fasta/')" ] || RC=2 fail "the FASTA set holds more than fasta/"
  [ -e "$RAW" ] && { chmod -R u+w "$RAW"; rm -rf "$RAW"; }
  mkdir -p "$RAW"
  lnset "$META" "$RAW" && lnset "$FA" "$RAW" || RC=2 fail "could not link the sets into $RAW"
  marker "$RAW"
  chmod -R a-w "$RAW"
  [ -z "$(find "$RAW" ! -type f ! -type d | head -n 1)" ] || RC=2 fail "$RAW holds a link or special file"
  (cd "$RAW" && find . -type f -printf '%P %s\n' | LC_ALL=C sort) > "$R/raw-listing.txt"
  echo "RAW $(wc -l < "$R/raw-listing.txt") files, listing $(sha256sum "$R/raw-listing.txt" | cut -d' ' -f1)"
  # the union's scratch WORK (its curated roots) is spent once union.txt is in the receipts; its disk goes to the run
  rm -rf "$HB/h1union"
  free=$(df -B1 --output=avail "$HB" | tail -n 1 | tr -d ' ')
  use=$(usable "$free") || exit $?
  rm -rf "$HB/h1est"
  runpy --raw "$RAW" --work "$HB/h1est" "${EXTRA[@]}" "${REFARG[@]}" --estimate-only --free-disk-bytes "$use" \
    > "$R/estimate.json" 2> "$R/estimate.err" || RC=4 fail "run.py --estimate-only refused at $use B usable ($free B free): $(tail -n 2 "$R/estimate.err" | tr '\n' ' ')"
  python3 -I - "$R/estimate.json" << 'EOF' || RC=4 fail "the estimate does not pass"
import json, sys
d = json.load(open(sys.argv[1]))
sys.exit(0 if d.get("ok") is True and d.get("union_without_fasta") == 0 else 1)
EOF
  grep -m1 'revision-12 gate' "$R/estimate.err" | tee "$R/gate-line.txt"
  if ! python3 -I -c 'import json, sys; sys.exit(0 if json.load(open(sys.argv[1])).get("synthetic") is True else 1)' "$SPECF"; then
    grep -q 'revision-12 gate: PASSED' "$R/estimate.err" || RC=4 fail "no 'revision-12 gate: PASSED' line"
  fi
  echo "$free" > "$R/free-at-estimate.txt"; echo "$use" > "$R/usable-at-estimate.txt"
  # the estimate's scratch WORK, likewise: estimate.json and its stderr stay in the receipts
  rm -rf "$HB/h1est"
  echo "verified: sets match their manifests; RAW read only; estimate ok at $use B usable ($free B free, less ${DISK_PAUSE_GIB:-30} GiB + 1 GiB)"
}

# prune DIR (KIT12, a box made from the data snapshot): removes the box's own copy of the snapshot's master root once
# the metadata and FASTA sets are linked from it (grm-boat.sh stage-local) and match their pins. The sets' files keep
# their inodes, so only the blocks no set links are freed; the template's copy is not touched. DIR must be a real
# folder directly under the box's home and none of the sets or run folders. Only folders are made writable, since a
# file shares its mode with its set's hard link. verify checks every set again after this. Writes receipts/h1/prune.txt.
cmd_prune() {
  local d=${1:-} x before after
  [ -n "$d" ] && [ -d "$d" ] && [ ! -L "$d" ] || RC=2 fail "prune needs the snapshot's data folder, a real folder: ${d:-none}"
  x=$(basename -- "$d")
  [ "$d" = "$HB/$x" ] && [ "$x" != . ] && [ "$x" != .. ] || RC=2 fail "$d is not a folder directly under $HB, written plainly"
  for x in "$META" "$FA" "$VOC" "$CI" "$RAW0" "$RAW" "$WORK" "$OUT" "$S" "$W" "$HB/h1union" "$HB/h1est" "$HB/ckpt" "$HB/port"; do
    [ "$d" = "$x" ] && RC=2 fail "$d is a set or run folder, not the snapshot's data"
  done
  checkset "$META" "$H1_META_SHA256"; checkset "$FA" "$H1_FASTA_SHA256"; checkset "$VOC" "$H1_VOCAB_SHA256"
  # the metadata set must have been linked from DIR, so DIR is the root the sets came from
  [ "$d/bvbrc/genomes.jsonl" -ef "$META/bvbrc/genomes.jsonl" ] || RC=2 fail "$META was not linked from $d (bvbrc/genomes.jsonl is not one file)"
  before=$(df -B1 --output=avail "$HB" | tail -n 1 | tr -d ' ')
  find "$d" -type d ! -perm -u+w -exec chmod u+w {} + || RC=2 fail "could not open $d's folders for removal"
  rm -rf -- "$d" && [ ! -e "$d" ] || RC=2 fail "could not remove $d"
  after=$(df -B1 --output=avail "$HB" | tail -n 1 | tr -d ' ')
  printf 'pruned\t%s\nfree_before\t%s\nfree_after\t%s\n' "$d" "$before" "$after" > "$R/prune.txt"
  echo "PRUNED $d: $((after - before)) B freed ($before B free before, $after B after)"
}

cmd_run() {
  local free use
  [ -d "$RAW" ] && [ ! -w "$RAW" ] || RC=2 fail "RAW is missing or writable; run verify first"
  refset
  free=$(df -B1 --output=avail "$HB" | tail -n 1 | tr -d ' ')
  use=$(usable "$free") || exit $?
  echo "$free" > "$R/free-at-start.txt"; echo "$use" > "$R/usable-at-start.txt"
  log "h1 run: spec $SPEC jobs ${H1_JOBS:-2} trace ${H1_TRACE:-0} free $free usable $use${REFW:+ reference $REFW}"
  # pg-run.sh: its own process group under guard.sh (disk and RAM pauses), output to logs/h1-run.log
  DATA=$WORK bash "$K/pg-run.sh" h1-run -- python3 -B "$RUNPY" "$SPECF" --raw "$RAW" --work "$WORK" --grm "$GRM" \
    --source-commit "$SOURCE_COMMIT" --sigdrug-vocabulary "$VOC/vocabulary" --jobs "${H1_JOBS:-2}" --free-disk-bytes "$use" \
    "${TRACE[@]}" "${EXTRA[@]}" "${REFARG[@]}"
  local rc=$?
  [ $rc = 0 ] || RC=1 fail "run.py exited $rc (logs/h1-run.log: $(tail -n 2 "$W/logs/h1-run.log" | tr '\n' ' '))"
  cmd_measure
}

cmd_measure() {
  python3 -I - "$WORK" "$W/logs/h1-run.time" "$R/disk-measured.json" << 'EOF'
import json, os, re, subprocess, sys
work, tfile, out = sys.argv[1:]
def du(*paths):
    paths = [p for p in paths if os.path.exists(p)]
    if not paths: return 0
    r = subprocess.run(["du", "-scB1", *paths], capture_output=True, text=True, check=True)
    return int(r.stdout.strip().splitlines()[-1].split()[0])
roots = {}
rdir = os.path.join(work, "roots")
for r in sorted(os.listdir(rdir)) if os.path.isdir(rdir) else []:
    top = os.path.join(rdir, r); sp = os.path.join(top, "species")
    ks = [os.path.join(sp, s, "ksets") for s in (os.listdir(sp) if os.path.isdir(sp) else [])]
    mx = [os.path.join(sp, s, "matrix") for s in (os.listdir(sp) if os.path.isdir(sp) else [])]
    files = [os.path.join(k, f) for k in ks if os.path.isdir(k) for f in os.listdir(k)]
    kset_files = [f for f in files if os.path.isfile(f) and not f.endswith(".json")]
    roots[r] = {"blocks_bytes": du(top), "ksets_blocks_bytes": du(*ks), "matrix_blocks_bytes": du(*mx),
                "ksets_apparent_bytes": sum(os.path.getsize(f) for f in kset_files), "kset_files": len(kset_files)}
    roots[r]["other_blocks_bytes"] = roots[r]["blocks_bytes"] - roots[r]["ksets_blocks_bytes"] - roots[r]["matrix_blocks_bytes"]
peak = None
if os.path.exists(tfile):
    m = re.search(r"Maximum resident set size \(kbytes\): (\d+)", open(tfile).read())
    peak = int(m.group(1)) * 1024 if m else None
est = os.path.join(work, "disk-estimate.json")
rec = {"roots": roots,
       "all_roots_blocks_bytes": du(rdir),        # du counts a hardlinked kset once per call
       "arms_blocks_bytes": du(os.path.join(work, "arms")),
       "work_blocks_bytes": du(work),
       "pool": None, "peak_rss_bytes": peak,
       "disk_estimate": json.load(open(est)) if os.path.exists(est) else None}
pools = [r for r in roots if r.startswith("pool")] or [r for r, v in roots.items() if v["kset_files"] == max(x["kset_files"] for x in roots.values())][:1]
if pools and roots[pools[0]]["kset_files"]:
    p = roots[pools[0]]
    rec["pool"] = {"root": pools[0], "kset_files": p["kset_files"], "mean_kset_apparent_bytes": p["ksets_apparent_bytes"] // p["kset_files"]}
json.dump(rec, open(out, "w"), indent=1, sort_keys=True)
print("measured: work %d B in blocks, roots together %d B, arms %d B" % (rec["work_blocks_bytes"], rec["all_roots_blocks_bytes"], rec["arms_blocks_bytes"]))
EOF
}

cmd_combine() {
  local args=() p n=0
  checkset "$CI" "$H1_COMBINE_SHA256"
  M=$(find "$CI" -maxdepth 1 -name 'MASTER-MANIFEST*.sha256' -type f | head -n 1)
  [ -n "$M" ] || RC=2 fail "no master manifest in $CI"
  for p in ${H1_PARTS//,/ }; do [ -d "$CI/parts/$p/h1work" ] || RC=2 fail "no part $p"; args+=(--work-from "$CI/parts/$p/h1work"); n=$((n + 1)); done
  [ $n -gt 0 ] || RC=2 fail "H1_PARTS names no part"
  rm -rf "$OUT"
  DATA=$OUT bash "$K/pg-run.sh" h1-combine -- python3 -B "$RUNPY" combine "$SPECF" "${args[@]}" --master-manifest "$M" --out "$OUT"
  local rc=$?
  [ $rc = 0 ] || RC=1 fail "run.py combine exited $rc (logs/h1-combine.log: $(tail -n 2 "$W/logs/h1-combine.log" | tr '\n' ' '))"
  # acceptance 11: one part combined must give that part's own J and K, after timings (the tests' timeless)
  if [ $n = 1 ]; then
    (cd "$S" && python3 -B - "$CI/parts/${H1_PARTS}/h1work" "$OUT" "$R/combine-equal.json") << 'EOF' || RC=1 fail "combine of one part differs from the part (receipts/h1/combine-equal.json)"
import json, sys
from pathlib import Path
sys.path.insert(0, ".")
from tests.input_comparison import test_run as base   # the e2e tests' own timeless()
part, out, rec = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3]
res = {n: (out / "j" / f"{n}.json").read_bytes() == (part / "j" / f"{n}.json").read_bytes() for n in ("contrasts", "planned", "missing")}
res["timeless"] = base.timeless(out) == base.timeless(part)
json.dump(res, open(rec, "w"), indent=1, sort_keys=True); print(json.dumps(res, sort_keys=True))
sys.exit(0 if all(res.values()) else 1)
EOF
  fi
  echo "combined $n part(s) into $OUT"
}

cmd_ckpt() {
  local L=$1 C src
  [[ "$L" =~ ^[A-Za-z0-9._-]+$ ]] || RC=2 fail "bad label"
  C=$HB/ckpt/$L; rm -rf "$C"; mkdir -p "$C/run"
  rsync -a --no-links --exclude views --exclude prelaunch-empty --exclude ksets --exclude fasta --exclude features "$W/" "$C/run/"
  if [ "${H1_ROLE:-run}" = combine ]; then src=$OUT name=h1combine; else src=$WORK name=h1work; fi
  # WORK without ksets/, matrix/, links, or anything of a root but its records and genome lists
  [ -d "$src" ] && rsync -a --no-links --exclude 'ksets/' --exclude 'matrix/' \
    --include '/roots/' --include '/roots/*/' --include '/roots/*/*.json' --include '/roots/*/species/' \
    --include '/roots/*/species/*/' --include '/roots/*/species/*/genomes.txt' --include '/roots/*/species/*/structure.json' \
    --include '/roots/*/species/*/counted.json' --exclude '/roots/**' "$src/" "$C/$name/"
  (cd "$C" && find . -type f ! -name SHA256SUMS -print0 | sort -z | xargs -0 -r sha256sum) > "$C/SHA256SUMS"
  echo "CKPT $L files=$(wc -l < "$C/SHA256SUMS") bytes=$(du -sb "$C" | cut -f1) sums=$(sha256sum "$C/SHA256SUMS" | cut -d' ' -f1)"
}

case "${1:-}" in
  prep) cmd_prep;; union) cmd_union;; verify) cmd_verify;; run) cmd_run;; measure) cmd_measure;; combine) cmd_combine;;
  ckpt) cmd_ckpt "${2:?label}";; prune) cmd_prune "${2:-}";;
  *) sed -n '2,/^set -u/p' "$0" | sed '$d'; exit 2;;
esac
